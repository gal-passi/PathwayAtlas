import itertools

from statistical_models import BackgroundModel
from utils import *

import os
import pandas as pd
from tqdm import tqdm

from esm import pretrained      # for model choosing and alphabet



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
#
#     # Step 6: Create GMMs distributions for all pathways
#     # TODO choose GMM number of component by BIC
#     # TODO

















# # TODO delete this part - gc content calc
# def update_gc_content_in_pickles(directory: str):
#     """
#     Iterates through all pickle files in a directory that contain dictionaries,
#     calculates the GC content from the 'na_seq' key, adds a 'gc_content' key,
#     and saves the updated dictionary back to the file.
#
#     Args:
#         directory (str): The path to the directory containing the pickle files.
#     """
#     if not os.path.isdir(directory):
#         print(f"Error: Directory not found at '{directory}'")
#         return
#
#     print(f"Starting GC content update process for files in '{directory}'...")
#     updated_count = 0
#     skipped_count = 0
#     error_count = 0
#
#     # Get a list of files to process to avoid issues with directory changes
#     files_to_process = [f for f in os.listdir(directory) if f.endswith(('.pickle', '.pkl'))]
#
#     for filename in files_to_process:
#         file_path = os.path.join(directory, filename)
#
#         try:
#             # --- 1. Safely Load the Pickle File (in 'read binary' mode) ---
#             with open(file_path, 'rb') as f:
#                 gene_dict = pickle.load(f)
#
#             # --- Safety Check: Ensure the loaded object is a dictionary ---
#             if not isinstance(gene_dict, dict):
#                 print(f"Warning: File '{filename}' does not contain a dictionary. Skipping.")
#                 error_count += 1
#                 continue
#
#             # --- 2. Check if Update is Needed ---
#             # Skips if 'gc_content' key already exists and has a valid value.
#             if 'gc_content' in gene_dict and gene_dict['gc_content'] is not None:
#                 skipped_count += 1
#                 continue
#
#             # --- 3. Robustly Calculate GC Content ---
#             # Safely get the sequence using .get(), which returns None if the key doesn't exist.
#             na_sequence = gene_dict.get('na_seq', None)
#
#             # Check if the sequence is a non-empty string before calculating.
#             if isinstance(na_sequence, str) and na_sequence:
#                 try:
#                     g_count = na_sequence.upper().count('G')
#                     c_count = na_sequence.upper().count('C')
#                     gc_fraction = (g_count + c_count) / len(na_sequence)
#                     # CHANGE: Use dictionary key assignment
#                     gene_dict['gc_content'] = gc_fraction
#                 except ZeroDivisionError:
#                     gene_dict['gc_content'] = 0.0
#             else:
#                 # Set to None if the sequence is missing, empty, or not a string.
#                 gene_dict['gc_content'] = None
#
#             # --- 4. Safely Save the Updated Dictionary Back to the File ---
#             # Open in 'write binary' mode to overwrite the original file.
#             with open(file_path, 'wb') as f:
#                 pickle.dump(gene_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
#
#             updated_count += 1
#             if updated_count % 100 == 0:
#                 print(f"  ... updated {updated_count} files.")
#
#         except pickle.UnpicklingError:
#             print(f"Warning: Could not unpickle '{filename}'. File might be corrupted. Skipping.")
#             error_count += 1
#             continue
#         except Exception as e:
#             # Catch any other unexpected errors (e.g., I/O issues, permissions).
#             print(f"An unexpected error occurred with '{filename}': {e}")
#             error_count += 1
#             continue
#
#     print("\n--- Update Complete ---")
#     print(f"Successfully updated: {updated_count} files.")
#     print(f"Skipped (already up-to-date): {skipped_count} files.")
#     print(f"Errors or invalid files: {error_count} files.")
#
# if __name__ == "__main__":
#     update_gc_content_in_pickles(KEGG_GENES_PATH)



# # TODO delete this part - it is the bic curves creator
# def create_GMM_debug():
#     """This function is only for choosing the number of components in the GMMs"""
#     print("BICing: collecting scores...")
#     model = BackgroundModel()
#     res = model.collect_scores()        # all pathways, all scores, from all genes
#     for pathway_name, scores_dict in itertools.islice(res.items(), 200):       # example for 5 pathways
#         joint_distribution, bins_edges = model.create_joint_distribution(res[pathway_name],
#                                                                          score_type="clinvar_reg_dis_ordered_prob")
#
#         # plot some gmm
#         #gmm, bic = model.GMM_the_distribution(joint_distribution, bins_edges)
#         #model.plot_1D_GMM(joint_distribution, bins_edges, gmm)
#
#         # plot gmm bic curve
#         print(f"BICing: {pathway_name}")
#         model.GMM_bic_curve(joint_distribution, bins_edges, filename=pathway_name + "_bic")
#
# if __name__ == "__main__":
#     create_GMM_debug()