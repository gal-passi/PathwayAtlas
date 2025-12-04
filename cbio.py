from utils import *
from definitions import *
import os
import pandas as pd
import pickle
from statistical_models import PermutationTest

"""
This module provides functions to interact with the cBioPortal API.
Run the `get_studies` function to download mutation data for all studies.
Run the `merge_studies_by_cancer` function to merge studies by cancer type.
Run the `expand_cancer_studies` function to add sequences to all mutations in 'filename' csv.
"""


def check_for_duplicates(dfs: dict) -> None:
    """Check for duplicate mutations in all the studies.
        @param dfs: dictionary of study id to df of that study, returned by `get_studies()`.
    """
    all_patients = []
    for study, df in dfs.items():
        patients = df[PATIENT_ID_COL].unique()
        all_patients.extend(patients)

    all_df = pd.DataFrame(all_patients, columns=[PATIENT_ID_COL, STUDY_ID_COL])

    overlaps = (all_df.groupby(PATIENT_ID_COL)[STUDY_ID_COL]
                .nunique()
                .reset_index()
                .query(STUDY_ID_COL + " > 1")
                )

    if overlaps.empty:
        print("No overlapping patients found.")
    else:
        print(f"found {len(overlaps)} overlapping patients.")


def get_tcga_studies(cbio: CbioApi) -> list:
    """Retrieve all TCGA studies from cBioPortal.
        @param cbio: cBioPortal API object
        @return: list of TCGA studies
    """
    all_studies = cbio.api.Studies.getAllStudiesUsingGET().result()
    tcga_studies = [study for study in all_studies if 'tcga' in study.studyId.lower()]
    return tcga_studies


def download_studies_dfs(cbio: CbioApi) -> dict:
    """Retrieve all cBioPortal studies.
        Downloads mutations for each study and saves them as CSV files.
        Returns a dictionary of study_id to DataFrame of mutations.
        @param cbio: cBioPortal API object
        @return: dictionary of study id to DataFrame of that study.
    """
    # ---------- Return existing dictionary if already downloaded ----------

    studies_dfs = open_df_pickle('studies_dfs.pkl')
    if studies_dfs != {}:
        return studies_dfs

    tcga_studies = get_tcga_studies(cbio) # get all TCGA studies

    for study in tcga_studies:
        study_id = study.studyId
        try:
            # ---------- Download mutations ----------

            results = cbio.download_study_mutations(study_id)
            outpath = STUDIES_PATH + "/" + study_id + MUTATIONS_CSV_SUFFIX

            if os.path.exists(outpath):
                print(f"{study_id} mutations already downloaded.")
                studies_dfs[study_id] = pd.read_csv(outpath)
            else:
                # ---------- Save DataFrame to CSV ----------

                df = cbio.study_to_csv(results, outpath=outpath)
                studies_dfs[study_id] = df
                print(f"Saved {study_id}_mutations.csv")

        except Exception as e:
            print(f"Skipping {study_id}: {e}")

    with open('studies_dfs.pkl', 'wb') as f:
        pickle.dump(studies_dfs, f)

    print("Downloaded all studies.")
    return studies_dfs

def merge_studies_by_cancer(cbio: CbioApi, dfs: dict, remove_duplicates: bool=True) -> dict:
    """
    merge the studies of the same cancer type into a single dataframe.
    @param cbio: cBioPortal API object
    @param dfs: dictionary of study id to DataFrame of that study, returned by 'get_studies()'.
    @param remove_duplicates: remove duplicate studies if True.
    """
    # ---------- Return existing dictionary if already merged ----------

    cancer_dfs = open_df_pickle('cancer_dfs.pkl')
    if cancer_dfs != {}:
        return cancer_dfs

    cancer_dfs = {}  # dictionary of cancer type to DataFrame
    cancer_types_dict = cbio.cancer_types_dict()

    for cancer_type, short_cancer_name in cancer_types_dict.items():

        # ---------- Get all studies for the given cancer type ----------

        study_ids, study_names = cbio.all_studies_by_keyword(short_cancer_name.lower())
        study_ids = [id for id in study_ids if id in dfs.keys()]

        if not study_ids:
            print(f"No studies found for cancer type {short_cancer_name.lower()}")
            continue  # skip to the next cancer type

        # ---------- Merge them all to one Dataframe ----------
        merged_df = pd.concat([dfs[study_id] for study_id in study_ids], ignore_index=True)

        if remove_duplicates:  # Optionally remove duplicates across studies
            merged_df.drop_duplicates(keep='first', inplace=True, ignore_index=True, subset=DUPLICATE_EXCLUSION_COLUMNS)

        outpath = CANCERS_PATH + "/" + short_cancer_name.lower() + MUTATIONS_CSV_SUFFIX

        if os.path.exists(outpath):
            print(f"{cancer_type} mutations already downloaded.")
            cancer_dfs[cancer_type] = pd.read_csv(outpath)
        else:
            # ---------- Save merged DataFrame to CSV ----------
            merged_df.to_csv(outpath, index=False)
            cancer_dfs[cancer_type] = merged_df
            print(f"Saved {short_cancer_name.lower() + MUTATIONS_CSV_SUFFIX}")

    with open('data/cancer_dfs.pkl', 'wb') as f:
        pickle.dump(cancer_dfs, f)

    print("Merged all cancer studies.")
    return cancer_dfs

def expand_cancer_studies(filename: str= '') -> None:
    if filename != '':
        expand_cancer_study(filename)
    else:
        for file in os.listdir(CANCERS_PATH):
            if file.endswith(MUTATIONS_CSV_SUFFIX):
                expand_cancer_study(file)

def expand_cancer_study(filename: str) -> None:
    """
    Add sequences to all mutations in the downloaded csvs.
    """
    print(f"Processing {filename}...")

    cancer_path = os.path.join(CANCERS_PATH, filename)

    try:
        df = pd.read_csv(cancer_path)
    except pd.errors.EmptyDataError:
        print(f"{filename} is empty. Skipping.")

    if REFERENCE_SEQ_COL and UNIPROT_ID_COL and KEGG_COL in df.columns:
        print(f"Sequences already added to {filename}. Skipping.")
        return

    # Initialize new columns
    df[REFERENCE_SEQ_COL] = None
    df[UNIPROT_ID_COL] = None
    df[KEGG_COL] = None

    df = add_seq_to_df(df)
    df.to_csv(cancer_path, index=False)

    print(f"All sequences added to {filename}.")

def add_seq_to_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Add sequences, Kegg ids and uniprot ids to all proteins in the given dataframe.
    @param df: DataFrame of proteins.
    @return: DataFrame with sequences added and updated protein_seq_dict.
    """
    uniprot = UniprotApi()

    for idx, row in df.iterrows():
        ref_name = row[PROTEIN_NAME_COL]
        ref_mut = row[VARIANT_COL]
        cache_file = os.path.join(SEQ_PATH, f"{ref_name}.pickle")

        uid, seq, kegg_id = None, None, None

        # Search in cache first
        if os.path.exists(cache_file):
            try:
                with open(cache_file, "rb") as f:
                    data = pickle.load(f)
                print(f"Loaded cached data for {ref_name}")
                uid = data.get("uniprot_id")
                seq = data.get("sequence")
                kegg_id = data.get("kegg_id")

            except Exception as e:
                pass
        else:
            # ---------- Fetch uid and sequence from Uniprot ----------

            iso_seq = None

            try:
                iso_seq = uniprot.expand_isoforms(ref_name, ref_mut)
            except Exception as e:
                print(e)

            if iso_seq is None or len(iso_seq) != 1:
                print(f"At row {idx}: no correct sequence found for {ref_name}")
            else:
                uid = list(iso_seq.keys())[0]
                seq = list(iso_seq.values())[0] # get the only sequence

            # ---------- Fetch Kegg ID from Kegg ----------

            try:
                kegg_id = ','.join(KeggApi.hugo_to_kegg_hsa(ref_name))
            except Exception as e:
                print(e)

            # ---------- Cache the result ----------

            try:
                data = {
                    "name": ref_name,
                    "uniprot_id": uid,
                    "kegg_id": kegg_id,
                    "sequence": seq,
                }
                with open(cache_file, "wb") as f:
                    pickle.dump(data, f)

            except Exception as e:
                print(f"Failed to cache data for {ref_name}: {e}")

        # ---------- Update DataFrame ----------

        print(f"Adding sequence for {ref_name}, at row {idx}")

        df.at[idx, UNIPROT_ID_COL] = uid
        df.at[idx, REFERENCE_SEQ_COL] = seq
        df.at[idx, KEGG_COL] = kegg_id

    return df
