from setuptools.dist import sequence

from utils import *
from Kegg import *
from definitions import *
import os

"""
This module provides functions to interact with the cBioPortal API.
Run the `get_studies` function to download mutation data for all studies.
"""


def check_for_duplicates(dfs):
    """Check for duplicate mutations in all the studies.
        @param dfs: dictionary of study id to df of that study, returned by `get_studies()`.
    """
    all_patients = []
    for study, df in dfs.items():
        patients = df['PatientId'].unique()
        all_patients.extend(patients)

    all_df = pd.DataFrame(all_patients, columns=['PatientId', 'StudyId'])

    overlaps = (all_df.groupby('PatientId')['StudyId']
                .nunique()
                .reset_index()
                .query("StudyId > 1")
                )

    if overlaps.empty:
        print("No overlapping patients found.")
    else:
        print(f"found {len(overlaps)} overlapping patients.")


def get_studies(cbio: CbioApi) -> dict:
    """Retrieve all cBioPortal studies.
        Downloads mutations for each study and saves them as CSV files.
        Returns a dictionary of study_id to DataFrame of mutations.
    """

    if os.path.exists('studies_dfs.pkl'):
        with open('studies_dfs.pkl', 'rb') as f:
            studies_dfs = pickle.load(f)
            print("Opened existing studies_dfs.pkl")
            return studies_dfs

    studies_dfs = {}
    all_studies = cbio.api.Studies.getAllStudiesUsingGET().result()

    for study in all_studies:
        study_id = study.studyId
        try:
            # Download mutations
            results = cbio.download_study_mutations(study_id)

            # Define file name for each study
            outpath = f"data/cbio/studies/{study_id}_mutations.csv"

            if os.path.exists(outpath):
                print(f"{study_id} mutations already downloaded.")
                studies_dfs[study_id] = pd.read_csv(outpath)
            else:
                # Save DataFrame to CSV
                df = cbio.study_to_csv(results, outpath=outpath)
                studies_dfs[study_id] = df
                print(f"Saved {study_id}_mutations.csv")
        except Exception as e:
            print(f"Skipping {study_id}: {e}")

    with open('studies_dfs.pkl', 'wb') as f:
        pickle.dump(studies_dfs, f)

    print("Downloaded all studies.")
    return studies_dfs


def merge_studies(cbio, dfs, remove_duplicates=True) -> dict:
    """
    merge the studies of the same cancer type into a single dataframe.
    @param cbio: cBioPortal API object
    @param dfs: dictionary of study id to DataFrame of that study, returned by 'get_studies()'.
    @param remove_duplicates: remove duplicate studies if True.
    """

    if os.path.exists('cancer_dfs.pkl'):
        with open('cancer_dfs.pkl', 'rb') as f:
            cancer_dfs = pickle.load(f)
            print("Opened existing cancer_dfs.pkl")
            return cancer_dfs

    cancer_dfs = {}
    cancer_types_dict = cbio.cancer_types_dict()

    for cancer_type, short_cancer_name in cancer_types_dict.items():
        # Get all studies for the given cancer type
        study_ids, study_names = cbio.all_studies_by_keyword(short_cancer_name.lower())

        study_ids = [id for id in study_ids if id in dfs.keys()]

        if not study_ids:
            print(f"No studies found for cancer type {short_cancer_name.lower()}")
            continue  # skip to the next cancer type

        # Merge them all to one Dataframe
        merged_df = pd.concat([dfs[study_id] for study_id in study_ids], ignore_index=True)

        # Optionally remove duplicates across studies
        if remove_duplicates:
            merged_df.drop_duplicates(keep='first', inplace=True, ignore_index=True, subset=DUPLICATE_EXCLUSION_COLUMNS)

        # Define file name for each cancer
        outpath = CANCERS_PATH + short_cancer_name.lower() + MUTATIONS_CSV_SUFFIX

        if os.path.exists(outpath):
            print(f"{cancer_type} mutations already downloaded.")
            cancer_dfs[cancer_type] = pd.read_csv(outpath)
        else:
            # Save merged DataFrame to CSV
            merged_df.to_csv(outpath, index=False)
            cancer_dfs[cancer_type] = merged_df
            print(f"Saved {short_cancer_name.lower() + MUTATIONS_CSV_SUFFIX}")

    with open('cancer_dfs.pkl', 'wb') as f:
        pickle.dump(cancer_dfs, f)

    print("Merged all cancer studies.")
    return cancer_dfs


def add_sequences_to_mutations():
    """
    Add sequences to all mutations in the downloaded csvs.
    """
    protein_seq_dict = {}

    if os.path.exists(PROTEIN_SEQUENCES_FILE):
        with open(PROTEIN_SEQUENCES_FILE, 'rb') as f:
            protein_seq_dict = pickle.load(f)
            print(f"Opened existing {PROTEIN_SEQUENCES_FILE}")

    for cancer_file in os.listdir(CANCERS_PATH):
        if cancer_file.endswith('_mutations.csv'):
            cancer_path = os.path.join(CANCERS_PATH, cancer_file)
            df = pd.read_csv(cancer_path)

            if REFERENCE_SEQ_COL and UNIPROT_ID_COL in df.columns:
                print(f"Sequences already added to {cancer_file}. Skipping.")
                continue

            df[REFERENCE_SEQ_COL] = None  # initialize the column
            df[UNIPROT_ID_COL] = None

            df, protein_seq_dict = add_seq_to_df(df, protein_seq_dict)
            df.to_csv(cancer_path, index=False)
            print(f"All sequences added to {cancer_file}.")

    # save updated cache
    with open(PROTEIN_SEQUENCES_FILE, 'wb') as f:
        pickle.dump(protein_seq_dict, f)
        print(f"Cache updated in {PROTEIN_SEQUENCES_FILE}")

def add_seq_to_df(df: pd.DataFrame, protein_seq_dict):
    """
    Add sequences to all proteins in the given dataframe.
    @param df: DataFrame of proteins.
    @param protein_seq_dict: dictionary of protein name to sequence.
    @return: DataFrame with sequences added and updated protein_seq_dict.
    """
    uniprot = UniprotApi()

    for idx, row in df.iterrows():
        ref_name = row[PROTEIN_NAME_COL]
        ref_mut = row[VARIANT_COL]

        # Search in cache first
        if ref_name in protein_seq_dict:
            uid, seq = protein_seq_dict[ref_name]
            df.at[idx, REFERENCE_SEQ_COL] = seq
            df.at[idx, UNIPROT_ID_COL] = uid
            continue

        # Fetch from Uniprot
        iso_seq = uniprot.expand_isoforms(ref_name, ref_mut)
        if iso_seq is None or len(iso_seq) != 1:
            print(f"Skipping row {idx}: no correct sequence found for {ref_name}")
            continue

        print(f"Adding sequence for {ref_name} at row {idx}")
        uid = list(iso_seq.keys())[0]
        seq = list(iso_seq.values())[0] # get the only sequence
        df.at[idx, UNIPROT_ID_COL] = uid
        df.at[idx, REFERENCE_SEQ_COL] = seq

        protein_seq_dict[ref_name] = uid, seq  # update cache

    return df, protein_seq_dict
