import pandas as pd
import glob
from os.path import join as pjoin
import os
from cbio import CbioApi
from definitions import *
import openpyxl

class CancerResultsAnalyzer:

    def __init__(self, cbio: CbioApi, results_csv_path):
        self.cbio = cbio
        self.cancer_name = os.path.basename(results_csv_path).replace(MUTATIONS_CSV_SUFFIX, '')

        # Load cancer mutations DataFrame
        cancer_csv_path = pjoin(CANCER_CSVS_MUTATIONS, f"{self.cancer_name}{MUTATIONS_CSV_SUFFIX}")
        try:
            self.mutations_df = pd.read_csv(cancer_csv_path)
        except FileNotFoundError:
            print(f"Could not find cancer mutations CSV: {cancer_csv_path}")
            self.mutations_df = pd.DataFrame()

        # Load cancer pathways results DataFrame
        try:
            self.results_df = pd.read_csv(results_csv_path)
        except FileNotFoundError:
            print(f"Could not find cancer pathways results CSV: {results_csv_path}")
            self.results_df = pd.DataFrame()

    def get_sub_df_for_pathway(self, pathway_name: str) -> pd.DataFrame:
        """
    .
        """
        pathway_dict_path = pjoin("/cs/labs/dina/ophirmil12/PathwayAtlas/", KEGG_PATHWAY_OBJECTS_PATH)
        pathway_dict_path = pjoin(pathway_dict_path, f"{pathway_name}.pickle")

        try:
            with open(pathway_dict_path, 'rb') as f:
                pathway_dict = pickle.load(f)
        except (FileNotFoundError, EOFError) as e:
            print(f"Could not load pathway dictionary {pathway_dict_path}: {e}")
            return pd.DataFrame()

        pathway_gene_ids = set(pathway_dict.keys())

        if not pathway_gene_ids:
            print(f"No gene IDs found for pathway: {os.path.basename(pathway_dict_path)}")
            return pd.DataFrame()

        def is_in_pathway(kegg_id_cell: str) -> bool:
            return any(gid in pathway_gene_ids for gid in str(kegg_id_cell).split(','))

        mask = self.mutations_df[KEGG_COL].apply(is_in_pathway)
        return self.mutations_df[mask]

    def summarize_cancer_results(self):

        cancer_type = self.cbio.get_cancer_type(self.cancer_name)

        num_significant_01 = self.results_df['significant_0.01'].sum()
        num_significant_05 = self.results_df['significant_0.05'].sum()

        all_pathways = self.results_df['pathway_name'].tolist()

        merged_sub_df = pd.DataFrame()
        for pathway in all_pathways:
            pathway_sub_df = self.get_sub_df_for_pathway(pathway)
            if pathway_sub_df.empty:
                continue
            merged_sub_df = pd.concat([merged_sub_df, pathway_sub_df], ignore_index=True)
            merged_sub_df = merged_sub_df.drop_duplicates()

        num_unique_mutations = merged_sub_df[[VARIANT_COL, PROTEIN_NAME_COL]].drop_duplicates().shape[0]
        num_samples = merged_sub_df.shape[0]
        num_patients = merged_sub_df[PATIENT_ID_COL].nunique()

        return {
            'cancer_short_name': self.cancer_name,
            'cancer_type': cancer_type,
            'num_unique_mutations': num_unique_mutations,
            'num_samples_used': num_samples,
            'num_patients_used': num_patients,
            'num_proteins_used': merged_sub_df[PROTEIN_NAME_COL].nunique(),
            'significant_pathways_0.01': num_significant_01,
            'significant_pathways_0.05': num_significant_05
        }

def summarize_results_all_cancers(cbio: CbioApi) -> None:
    """
    Process all result CSV files in the CANCER_PATHWAY_RESULTS directory.
    Summarize the number of significant pathways for each cancer type and save to a new CSV.
    @param cbio: cBioPortal API object
    """
    all_files = glob.glob(pjoin(CANCER_SCORES_KL_PATH, "*.csv"))
    if not all_files:
        print(f"No result CSV files found in the specified directory: {CANCER_SCORES_KL_PATH}")
        return

    rows = []

    for filename in all_files:

        cancer_analyzer = CancerResultsAnalyzer(cbio, filename)
        summary = cancer_analyzer.summarize_cancer_results()
        rows.append(summary)

    summary_df = pd.DataFrame(rows)

    summary_df = summary_df.sort_values(by='significant_pathways_0.01', ascending=True)
    summary_df.to_csv(pjoin(LOTEM_RESULTS_PATH, "cancer_statistics_summary.csv"), index=False)
    print("Summary of significant pathways saved.")

def add_statistics_to_excel(cbio: CbioApi, xl_file: str):
    sheets = pd.read_excel(xl_file, sheet_name=None)
    all_files = glob.glob(pjoin(CANCER_SCORES_KL_PATH, "*.csv"))
    if not all_files:
        print(f"No result CSV files found in the specified directory: {CANCER_SCORES_KL_PATH}")
        return

    updated_sheets = {}
    cancer_types_dict = cbio.cancer_types_dict()
    print(cancer_types_dict)
    for sheet_name, df in sheets.items():
        print("Processing sheet:", sheet_name)
        if sheet_name == 'Uterine Carcinosarcoma-Uterine ':
            cancer_short_name = 'ucs'
        else:
            cancer_short_name = cbio.get_cancer_short_name(sheet_name)
        matching_file = next((f for f in all_files if cancer_short_name in os.path.basename(f)), None)

        df['num_unique_mutations'] = 0
        df['num_patients'] = 0
        df['num_proteins'] = 0
        df = df.rename(columns={"n": "num_samples"})

        if matching_file:
            cancer_analyzer = CancerResultsAnalyzer(cbio, matching_file)
            for pathway in df['pathway_name']:
                pathway_sub_df = cancer_analyzer.get_sub_df_for_pathway(pathway)
                if pathway_sub_df.empty:
                    continue
                num_unique_mutations = pathway_sub_df[[VARIANT_COL, PROTEIN_NAME_COL]].drop_duplicates().shape[0]
                num_patients = pathway_sub_df[PATIENT_ID_COL].nunique()
                df.loc[df['pathway_name'] == pathway, 'num_unique_mutations'] = num_unique_mutations
                df.loc[df['pathway_name'] == pathway, 'num_patients'] = num_patients
                df.loc[df['pathway_name'] == pathway, 'num_proteins'] = pathway_sub_df[PROTEIN_NAME_COL].nunique()

        updated_sheets[sheet_name] = df

    with pd.ExcelWriter(xl_file, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
        for sheet_name, df in updated_sheets.items():
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print("Updated Excel file with additional statistics.")

