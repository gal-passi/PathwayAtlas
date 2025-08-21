from utils import *
from Kegg import *


"""
This module provides functions to interact with the cBioPortal API.
Run the `get_studies` function to download mutation data for all studies.
"""


def check_for_duplicates(dfs):
    """Check for duplicate mutations in all the studies.
        @param dfs: dictionary of study id to df of that study
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


def get_studies():
    """Retrieve all cBioPortal studies.
        Downloads mutations for each study and saves them as CSV files.
        Returns a dictionary of study_id to DataFrame of mutations.
    """
    cbio = CbioApi()
    dfs = {}
    all_studies = cbio.api.Studies.getAllStudiesUsingGET().result()

    for study in all_studies:
        study_id = study.studyId
        try:
            # download mutations
            results = cbio.download_study_mutations(study_id)

            # define file name for each study
            outpath = f"data/cbio/studies/{study_id}_mutations.csv"

            if os.path.exists(outpath):
                print(f"{study_id} already downloaded.")
            else:
                # save
                df = cbio.study_to_csv(results, outpath=outpath)
                dfs[study_id] = df
                print(f"Saved {outpath}")
        except Exception as e:
            print(f"Skipping {study_id}: {e}")

    check_for_duplicates(dfs)
    return dfs
