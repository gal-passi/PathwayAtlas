import sys

from disorder_predictor import DisorderPredict
from statistical_models import PermutationTest
from utils import *

import os
import pandas as pd
from tqdm import tqdm

from esm import pretrained      # for model choosing and alphabet
from cbio import *
from process_results import add_statistics_to_excel, summarize_results_all_cancers

def init_kegg_genome(recalc=False):
    """Initialize KEGG genome by creating KeggGene objects for all genes in KEGG database."""
    genes = KeggApi().get_all_genes().keys()
    if not recalc:
        genes = set(genes) - kegg_genes_in_dataset()
    genes = list(genes)
    def objects_creator(*gene_ids):
        kegg = KeggApi()
        genes_dict = kegg.genes_info(list(gene_ids))
        for kegg_id, data in genes_dict.items():
            gene = KeggGene(kegg_id, default_init=False)
            gene.create_from_dict(data)

    target = objects_creator
    #  at most 10 genes per request
    tasks = [genes[i: i + KEGG_MAX_IDS] for i in range(0, len(genes), KEGG_MAX_IDS)]
    multiprocess_task(tasks=tasks, target=target, workers=KEG_API_RECOMMENDED_WORKERS)


def genes_all_snvs(kegg: KeggApi, recalc=False):
    """Generate all SNVs for all KEGG genes and save them to files."""
    all_genes = list(kegg.get_all_genes().keys())
    if not recalc:
        all_genes = set(all_genes) - kegg_genes_in_dataset()
    all_genes = list(all_genes)
    print(f"Total genes to process SNVs for: {len(all_genes)}")

    for gene_id in tqdm(all_genes, desc="Generating gene SNVs", unit="gene"):
        try:
            gene = KeggGene(gene_id)
            out_file = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene.kegg_id}.csv")
            if not os.path.exists(out_file):
                gene.all_snvs(outpath=out_file, index=True)
        except Exception as e:
            print(f"[Error] {gene_id}: {e}")


def pathways_and_modules_dict(kegg: KeggApi):
    """Create a dictionary for each KEGG pathway and module mapping gene IDs to their SNV files."""

    def process_collection(collection, network_type):
        for kegg_id in tqdm(collection, desc=f"Mapping {network_type}s", unit=network_type):
            try:
                gene_list = kegg.get_gene_list(kegg_id)
                gene_map = {}

                for gene_id in gene_list:
                    snv_file = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene_id}.csv")
                    if os.path.exists(snv_file):
                        gene_map[gene_id] = snv_file
                    else:
                        print(f"[Missing] SNV file not found for gene: {gene_id} in p/m: {kegg_id}")
                        # If SNV file is not found, we still want to keep the gene in the map
                        # The main reason of not having SNV file is that the len(gene)%3!=0
                        gene_map[gene_id] = None

                out_file = os.path.join(KEGG_PATHWAY_OBJECTS_PATH, f"{kegg_id}.pickle")
                pd.to_pickle(gene_map, out_file)
            except Exception as e:
                print(f"[Error] {kegg_id}: {e}")

    # Process pathways
    pathways = kegg.get_all_pathways()
    process_collection(pathways, 'pathway')

    # Process modules
    modules = kegg.get_all_modules()
    process_collection(modules, 'module')


def disordered_aa_prediction(kegg: KeggApi, recalc=False):
    """
    Predicts disordered regions for all KEGG genes and ADDS the predictions as columns
    to the existing SNV files in KEGG_PATHWAY_MUTATIONS_PATH.

    This function caches the cleaned sequences to speed up subsequent runs.
    If `recalc=True`, it will re-process all genes. Otherwise, it will only process
    genes whose corresponding SNV file does not yet have a 'disorder_score' column.

    Args:
        kegg (KeggApi): An instance of the KeggApi class.
        recalc (bool): If True, re-calculate and overwrite for all genes.
    """
    os.makedirs(KEGG_DISORDERED_AA_PATH, exist_ok=True)
    sequences_pickle_file = os.path.join(KEGG_DISORDERED_AA_PATH, "sequences_to_predict.pkl")

    disorder_predict = DisorderPredict(kegg)

    # Step 1: Load or create the dictionary of cleaned sequences [V3].
    sequences_to_predict = disorder_predict.load_sequences_to_predict(sequences_pickle_file, recalc)
    if not sequences_to_predict:
        return

    # Step 2: Predict disorder for the sequences.
    disorder_predict.predict(sequences_to_predict)


def create_llr_and_dis_scoring(kegg: KeggApi):
    """Create LLR scoring and esm-disorder scoring
    for all KEGG proteins using ESM1b model"""
    # Load ESM1b model
    model, alphabet = pretrained.load_model_and_alphabet(ESM1B_MODEL)
    calculator = ScoringCalculator(model=model, alphabet=alphabet)

    # Get KEGG genes
    all_genes = kegg.get_all_genes().keys()
    print(f"Total genes to process: {len(all_genes)}")

    os.makedirs(KEGG_PATHWAY_MUTATIONS_PATH, exist_ok=True)

    for kegg_id in tqdm(all_genes, desc="Generating LLR scoring", unit="gene"):
        try:
            # get sequence for the gene
            gene = KeggGene(kegg_id)
            seq = gene.aa_seq

            if gene.coding_type != "CDS":
                print(f"Skipping {kegg_id}: not a CDS, a {gene.coding_type}")
                continue
            elif not seq:
                print(f"Skipping {kegg_id}: empty sequence")
                continue
            elif any(aa not in VALID_AA for aa in seq):
                print(f"Skipping {kegg_id}: invalid sequence")
                continue

            # Compute mutation scores [note that pathways and modules are already made]
            csv_path = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{kegg_id}.csv")
            calculator.save_mutation_scores_to_csv(seq, csv_path)

        except Exception as e:
            print(f"Error processing {kegg_id}: {e}")


# if __name__ == '__main__':
#     # Lab Notebook:
#     #    https://docs.google.com/document/d/1XR21LBpqW3q96BjExqsbH6JgEhV3Yc9sJuzG3abzUmY/edit?usp=sharing
#
#
#     kegg = KeggApi()
#
#     """
#     # Step 1: Initialize KEGG genome (download all genes) [./data/kegg/genes]
#     print("Downloading/Loading KEGG genome...")
#     init_kegg_genome()
#
#     # step 2: create all_snvs for all genes [./data/kegg/pathways/snvs]
#     print("Creating all SNVs in all genes...")
#     genes_all_snvs(kegg, recalc=True)
#
#     # Step 3: build dict objects for each pathway from gene_id to snvs file [./data/kegg/pathways/objects]
#     print("Mapping genes snvs to pathways and modules...")
#     pathways_and_modules_dict(kegg)
#     # Notice that some values in the dict may be None, which means that the SNV file is not available for that gene.
#     """
#     # Step 4: for each AA, determine if it`s in disordered region (with metapredict)
#     print("Disordered predictions...")
#     disordered_aa_prediction(kegg)  # , recalc=True
#
#     # Step 5: create LLR scoring and disorder scoring (clinvar_reg_dis_ordered_prob) from ClinVar regressor,
#     # And global regressor for all KEGG proteins
#     # using ESM1b model (GPU recommended) [./data/kegg/pathways/snvs]
#     print("Creating mutations scoring for all KEGG proteins...")
#     create_llr_and_dis_scoring(kegg)


if (__name__ == '__main__'):
    # Lab Notebook:
    # https://docs.google.com/document/d/1XR21LBpqW3q96BjExqsbH6JgEhV3Yc9sJuzG3abzUmY/edit?usp=sharing
    args = sys.argv[1:]
    if len(args) < 2:
        print("Usage: python main.py <cancer_name_mutations.csv>")
        sys.exit(1)

    cancer_file = args[0]
    distance_metric = args[1]

    cbio = CbioApi()
    # excel_path = os.path.join("results_and_graphs/pathway_analysis_combined_for_michal.xlsx")
    # add_statistics_to_excel(cbio, excel_path)
    # summarize_results_all_cancers(cbio)
    perm_tester = PermutationTest(cancer_file)
    print("Performing permutation test and FDR correction...")
    perm_tester.run_permutation_test(distance_metric)

    """
>>>>>>> origin/master


<<<<<<< HEAD
=======

    # Step 1: Initialize KEGG genome (download all genes) [./data/kegg/genes]
    print("Downloading/Loading KEGG genome...")
    init_kegg_genome()

    # step 2: create all_snvs for all genes [./data/kegg/pathways/snvs]
    print("Creating all SNVs in all genes...")
    genes_all_snvs(kegg, recalc=True)

    # Step 3: build dict objects for each pathway from gene_id to snvs file [./data/kegg/pathways/objects]
    print("Mapping genes snvs to pathways and modules...")
    pathways_and_modules_dict(kegg)
    # Notice that some values in the dict may be None, which means that the SNV file is not available for that gene.
    
    # Step 4: for each AA, determine if it`s in disordered region
    print("Disordered predictions...")
    disordered_aa_prediction(kegg, recalc=True)

    # Step 5: create LLR scoring and disorder scoring for all KEGG proteins
    # using ESM1b model (GPU recommended) and metapredict [./data/kegg/pathways/snvs]
    print("Creating LLR scoring for all KEGG proteins...")
    create_llr_scoring(kegg)
"""

