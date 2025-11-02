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
    studies_dfs = open_df_pickle('studies_dfs.pkl')

    if studies_dfs != {}:
        return studies_dfs

    tcga_studies = get_tcga_studies(cbio)

    for study in tcga_studies:
        study_id = study.studyId
        try:
            # Download mutations
            results = cbio.download_study_mutations(study_id)

            # Define file name for each study
            outpath = STUDIES_PATH + "/" + study_id + MUTATIONS_CSV_SUFFIX

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

def merge_studies_by_cancer(cbio, dfs, remove_duplicates=True) -> dict:
    """
    merge the studies of the same cancer type into a single dataframe.
    @param cbio: cBioPortal API object
    @param dfs: dictionary of study id to DataFrame of that study, returned by 'get_studies()'.
    @param remove_duplicates: remove duplicate studies if True.
    """
    cancer_dfs = open_df_pickle('cancer_dfs.pkl')
    if cancer_dfs != {}:
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
        outpath = CANCERS_PATH + "/" + short_cancer_name.lower() + MUTATIONS_CSV_SUFFIX

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

def expand_cancer_studies(cancer_file=''):
    protein_seq_dict = open_df_pickle(PROTEIN_SEQUENCES_FILE)

    if cancer_file != '':
        expand_cancer_study(protein_seq_dict, cancer_file)
    else:
        for file in os.listdir(CANCERS_PATH):
            if file.endswith(MUTATIONS_CSV_SUFFIX):
                protein_seq_dict = expand_cancer_study(protein_seq_dict, file)

def expand_cancer_study(protein_seq_dict, cancer_file):
    """
    Add sequences to all mutations in the downloaded csvs.
    """
    print(f"Processing {cancer_file}...")
    cancer_path = os.path.join(CANCERS_PATH, cancer_file)
    df = pd.read_csv(cancer_path)

    if REFERENCE_SEQ_COL and UNIPROT_ID_COL and KEGG_COL in df.columns:
        print(f"Sequences already added to {cancer_file}. Skipping.")
        return

    # Initialize new columns
    df[REFERENCE_SEQ_COL] = None
    df[UNIPROT_ID_COL] = None
    df[KEGG_COL] = None

    df, protein_seq_dict = add_seq_to_df(df, protein_seq_dict)
    df.to_csv(cancer_path, index=False)

    print(f"All sequences added to {cancer_file}.")
    return protein_seq_dict

def add_seq_to_df(df: pd.DataFrame, protein_seq_dict):
    """
    Add sequences and uniprot ids to all proteins in the given dataframe.
    @param df: DataFrame of proteins.
    @param protein_seq_dict: dictionary of protein name to sequence.
    @return: DataFrame with sequences added and updated protein_seq_dict.
    """
    uniprot = UniprotApi()

    for idx, row in df.iterrows():
        ref_name = row[PROTEIN_NAME_COL]
        ref_mut = row[VARIANT_COL]

        uid, seq, kegg_id = None, None, None

        # Search in cache first
        if ref_name in protein_seq_dict:
            uid, seq = protein_seq_dict[ref_name]
        else:
            # Fetch from Uniprot
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
 
        try:
            kegg_id = ','.join(KeggApi.hugo_to_kegg_hsa(ref_name))
        except Exception as e:
            print(e)

        print(f"Adding sequence for {ref_name}, Kegg ID: {kegg_id} at row {idx}")

        df.at[idx, UNIPROT_ID_COL] = uid
        df.at[idx, REFERENCE_SEQ_COL] = seq
        df.at[idx, KEGG_COL] = kegg_id

        protein_seq_dict[ref_name] = uid, seq  # update cache

    return df, protein_seq_dict



