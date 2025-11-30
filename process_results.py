import pandas as pd
import glob
from os.path import join as pjoin
import os
from cbio import CbioApi
from definitions import *

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
        num_samples = self.results_df['n'].sum()

        all_pathways = self.results_df['pathway_name'].tolist()

        num_unique_mutations, num_patients = 0, 0

        for pathway in all_pathways:
            pathway_sub_df = self.get_sub_df_for_pathway(pathway)
            if pathway_sub_df.empty:
                continue

            num_unique_mutations += pathway_sub_df[VARIANT_COL].nunique()
            num_patients += pathway_sub_df[PATIENT_ID_COL].nunique()

        return {
            'cancer_short_name': self.cancer_name,
            'cancer_type': cancer_type,
            'num_unique_mutations': num_unique_mutations,
            'num_samples_used': num_samples,
            'num_patients_used': num_patients,
            'significant_pathways_0.01': num_significant_01,
            'significant_pathways_0.05': num_significant_05
        }

def summarize_results(cbio: CbioApi) -> None:
    """
    Process all result CSV files in the CANCER_PATHWAY_RESULTS directory.
    Summarize the number of significant pathways for each cancer type and save to a new CSV.
    @param cbio: cBioPortal API object
    """
    all_files = glob.glob(pjoin(CANCER_PATHWAY_RESULTS, "*.csv"))
    if not all_files:
        print(f"No result CSV files found in the specified directory: {CANCER_PATHWAY_RESULTS}")
        return

    rows = []

    for filename in all_files:

        cancer_analyzer = CancerResultsAnalyzer(cbio, filename)
        summary = cancer_analyzer.summarize_cancer_results()
        rows.append(summary)

    summary_df = pd.DataFrame(rows)

    summary_df = summary_df.sort_values(by='significant_pathways_0.01', ascending=True)
    summary_df.to_csv(pjoin(LOTEM_RESULTS_PATH, "results_summary.csv"), index=False)
    print("Summary of significant pathways saved.")

