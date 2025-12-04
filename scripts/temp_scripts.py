import os
import sys
try:
    back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(back_dir)
except NameError:
    print("Running in an environment where __file__ is not defined. Assuming project root is correctly configured.")














# # TODO  - removing columns that added by mistake to cancer CSVs files
# import os
# import pandas as pd
# import glob
#
# def clean_csv_files_in_folder(folder_path: str):
#     """
#     Scans a folder for CSV files and removes a predefined list of columns.
#
#     Args:
#         folder_path (str): The path to the folder containing the CSV files.
#     """
#     # 1. Define all the column names you want to delete
#     columns_to_delete = [
#         'esm_log_probs', 'disorder_score', 'is_disordered', 'esm_min_max_naive',
#         'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob',
#         'disorder_score.1', 'is_disordered.1', 'esm_log_probs.1',
#         'clinvar_reg_dis_ordered_prob.1', 'clinvar_reg_global_prob.1',
#         'disorder_score.2', 'is_disordered.2', 'esm_log_probs.2',
#         'clinvar_reg_dis_ordered_prob.2', 'clinvar_reg_global_prob.2'
#     ]
#
#     # 2. Find all CSV files in the target folder
#     csv_files = glob.glob(os.path.join(folder_path, '*.csv'))
#
#     if not csv_files:
#         print(f"No CSV files found in '{folder_path}'.")
#         return
#
#     print(f"Found {len(csv_files)} CSV files. Starting cleanup...")
#
#     # 3. Loop through each file
#     for file_path in csv_files:
#         try:
#             df = pd.read_csv(file_path)
#
#             # Find which of the columns to delete actually exist in this file
#             cols_in_file_to_drop = [col for col in columns_to_delete if col in df.columns]
#
#             if cols_in_file_to_drop:
#                 print(f"Processing '{os.path.basename(file_path)}'...")
#                 # Drop the columns
#                 df.drop(columns=cols_in_file_to_drop, inplace=True)
#
#                 # Save the modified DataFrame back to the original file path
#                 df.to_csv(file_path, index=False)
#                 print(f"  -> Removed columns: {', '.join(cols_in_file_to_drop)}")
#             else:
#                 print(f"Skipping '{os.path.basename(file_path)}' (no relevant columns found).")
#
#         except Exception as e:
#             print(f"Could not process file '{os.path.basename(file_path)}'. Error: {e}")
#
#     print("\nCleanup complete.")
#
# # --- MAIN EXECUTION ---
# if __name__ == "__main__":
#     # IMPORTANT: Change this path to your folder containing the CSVs
#     target_folder = '/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/cancers'
#
#     clean_csv_files_in_folder(target_folder)



# # TODO  - checking one-to-one match of UniprotId to ReferenceSeq in the cancer CSVs
# import os
# import csv
# from collections import defaultdict
#
# def check_uniprot_sequence_matches_across_all_files(folder_path):
#     """
#     Goes through all CSV files in a folder, aggregating UniprotId and ReferenceSeq
#     columns to check for a one-to-one match of ID to sequence across all files.
#     Files without both specified columns are skipped.
#
#     Args:
#         folder_path (str): The path to the folder containing the CSV files.
#     """
#
#     """
#         # --- Execution ---
#     # The specified folder path
#     folder_path = '/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/cancers/'
#     """
#
#     uniprot_to_sequences = defaultdict(set)
#     processed_files_count = 0
#     skipped_files_count = 0
#
#     if not os.path.isdir(folder_path):
#         print(f"Error: The specified folder does not exist: {folder_path}")
#         return
#
#     print(f"Starting to process files in: {folder_path}\n")
#
#     for filename in os.listdir(folder_path):
#         if filename.endswith('.csv'):
#             file_path = os.path.join(folder_path, filename)
#             try:
#                 with open(file_path, mode='r', encoding='utf-8') as infile:
#                     # First, check for the required headers
#                     reader = csv.reader(infile)
#                     try:
#                         header = next(reader)
#                     except StopIteration:
#                         # This handles empty files
#                         print(f"Skipping empty file: {filename}")
#                         skipped_files_count += 1
#                         continue
#
#                     if 'UniprotId' in header and 'ReferenceSeq' in header:
#                         # If headers are present, process the file
#                         infile.seek(0) # Reset file reader to the beginning
#                         dict_reader = csv.DictReader(infile)
#                         for row in dict_reader:
#                             uniprot_id = row.get('UniprotId')
#                             ref_seq = row.get('ReferenceSeq')
#                             if uniprot_id and ref_seq: # Ensure the values are not empty
#                                 uniprot_to_sequences[uniprot_id].add(ref_seq)
#                         processed_files_count += 1
#                     else:
#                         # Skip the file if the required columns are not in the header
#                         print(f"Skipping file: {filename} (missing 'UniprotId' or 'ReferenceSeq' column)")
#                         skipped_files_count += 1
#             except Exception as e:
#                 print(f"An error occurred while processing {filename}: {e}")
#                 skipped_files_count += 1
#
#
#     print("\n--- Processing Complete ---")
#     print(f"Total files processed: {processed_files_count}")
#     print(f"Total files skipped: {skipped_files_count}")
#     print("---------------------------\n")
#
#
#     # Identify UniprotIds with more than one unique sequence
#     mismatched_ids = {
#         uniprot_id: sequences
#         for uniprot_id, sequences in uniprot_to_sequences.items()
#         if len(sequences) > 1
#     }
#
#     # Report the findings
#     if mismatched_ids:
#         print("Mismatch Found: The following UniprotIds are linked to multiple unique ReferenceSeqs:")
#         for uniprot_id, sequences in mismatched_ids.items():
#             print(f"\n- UniprotId: {uniprot_id}")
#             for seq in sequences:
#                 print(f"  - Sequence: {seq}")
#     else:
#         print \
#             ("Success: All UniprotIds have a consistent one-to-one match with their ReferenceSeq across all processed files.")






# # TODO  - retry creating CSV files [all found to be non-CDS]
# problematic_ids = [
#     'hsa:693200',
#     'hsa:693210',
#     'hsa:7012',
#     'hsa:7296',
#     'hsa:768220',
#     'hsa:780851',
#     'hsa:780852',
#     'hsa:780853',
#     'hsa:780854',
#     'hsa:83642',
#     'hsa:85465',
#     'hsa:9403',
#     'hsa:94163',
#     ...
# ]
#
# print(f"Starting reprocessing for {len(problematic_ids)} problematic files.")
# kegg = KeggApi()
#
# # =================================================================================
# # Phase 1 & 2: Re-fetch KEGG data and generate base SNV CSVs
# # =================================================================================
# print("\n--- Phase 1/2: Re-fetching gene data and generating base SNV files ---")
# sequences_to_predict = {}
# for kegg_id in tqdm(problematic_ids, desc="Generating base SNVs"):
#     try:
#         # Step 1: Delete old pickle object to force a fresh download
#         gene_pickle_path = os.path.join(KEGG_GENES_PATH, f"{kegg_id.replace(':', '_')}.pickle")
#         if os.path.exists(gene_pickle_path):
#             os.remove(gene_pickle_path)
#
#         # This will now fetch from KEGG and save a new pickle
#         gene = KeggGene(kegg_id)
#
#         # Validate that the gene is suitable for SNV generation
#         if gene.coding_type != "CDS":
#             print(f"[Info] Skipping {kegg_id}: Not a CDS protein ({gene.coding_type}).")
#             continue
#         if not gene.na_seq or len(gene.na_seq) == 0:
#             print(f"[Info] Skipping {kegg_id}: Missing nucleic acid sequence.")
#             continue
#
#         # Step 2: Create the base SNV file
#         out_file = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{gene.kegg_id}.csv")
#         gene.all_snvs(outpath=out_file)
#
#         # Collect valid AA sequence for the next phase
#         if gene.aa_seq:
#             # Clean sequence based on metapredict V3 rules
#             processed_seq = gene.aa_seq.replace('B', 'N').replace('U', 'C').replace('X', 'G').replace('Z', 'Q')
#             processed_seq = processed_seq.replace(' ', '').replace('*', '').replace('-', '')
#             if processed_seq:
#                 sequences_to_predict[kegg_id] = processed_seq
#
#     except Exception as e:
#         print(f"[ERROR] Failed during base SNV generation for {kegg_id}: {e}")
#
# # =================================================================================
# # Phase 3: Add disorder predictions to the newly created CSVs
# # =================================================================================
# print("\n--- Phase 3: Predicting and merging disorder scores ---")
# if sequences_to_predict:
#     disorder_predictor = DisorderPredict(kegg)
#     disorder_predictor.predict(sequences_to_predict)
# else:
#     print("No valid sequences found to predict disorder.")
#
# # =================================================================================
# # Phase 4: Add ESM-based LLR and ClinVar regression scores
# # =================================================================================
# print("\n--- Phase 4: Calculating and adding final ESM-based scores (This may take a while) ---")
# try:
#     model, alphabet = pretrained.load_model_and_alphabet("esm1b_t33_650M_UR50S")
#     calculator = ScoringCalculator(model=model, alphabet=alphabet)
#
#     for kegg_id in tqdm(problematic_ids, desc="Applying ESM Scores"):
#         try:
#             # We need a fresh KeggGene instance to get the sequence
#             gene = KeggGene(kegg_id)
#             seq = gene.aa_seq
#
#             csv_path = os.path.join(KEGG_PATHWAY_MUTATIONS_PATH, f"{kegg_id}.csv")
#
#             # Check if the base file was created and the sequence is valid
#             if not os.path.exists(csv_path):
#                 print(f"[Warning] Skipping scoring for {kegg_id}: Base CSV not found.")
#                 continue
#             if not seq:
#                 print(f"[Warning] Skipping scoring for {kegg_id}: No amino acid sequence available.")
#                 continue
#
#             # The function reads the CSV, adds all score columns, and saves it in-place
#             calculator.save_mutation_scores_to_csv(seq, csv_path)
#
#         except Exception as e:
#             print(f"[ERROR] Failed during scoring for {kegg_id}: {e}")
#
# except Exception as e:
#     print(f"[FATAL ERROR] Could not initialize ESM model. Aborting scoring phase. Error: {e}")
#     print("Please ensure 'esm' and 'torch' are installed and a GPU is available if possible.")
#
# print("\nReprocessing complete.")






# # TODO  - check what csv files are missing
# import os
# import pickle
#
# # 1. Initialize an empty set to store the keys.
# #    A 'set' is a great choice here because it automatically prevents duplicates.
# #    If the same key is found in multiple files, it will only be stored once.
# missing_files = set()
#
# # Define the directory where the pickled object files are located.
# objects_dir = "./data/kegg/pathways/objects"
#
# # 2. Loop through every file in the specified directory.
# for path_file in os.listdir(objects_dir):
#     # Construct the full, correct path to the file.
#     file_path = os.path.join(objects_dir, path_file)
#
#     try:
#         # 3. Open each file in "read binary" ('rb') mode, which is required for pickle.
#         #    The 'with' statement ensures the file is closed automatically.
#         with open(file_path, "rb") as f:
#
#             # 4. Load the Python object (expected to be a dictionary) from the file.
#             #    Then, iterate through its key-value pairs.
#             for key, val in pickle.load(f).items():
#
#                 # 5. Check if the value 'val' is "falsy".
#                 #    In Python, "falsy" values include:
#                 #    - None
#                 #    - False
#                 #    - An empty list: []
#                 #    - An empty dictionary: {}
#                 #    - An empty string: ""
#                 #    - The number zero: 0
#                 if not val:
#                     # 6. If the value is falsy, add the corresponding 'key' to the set.
#                     missing_files.add(key)
#     except pickle.UnpicklingError:
#         print(f"Warning: Could not unpickle file '{path_file}'. It might be corrupted or not a pickle file.")
#     except Exception as e:
#         print(f"An error occurred with file '{path_file}': {e}")
#
#
# # 7. Convert the set of unique keys into a list.
# missing_list = list(missing_files)
#
# # 8. Sort the list alphabetically/numerically.
# missing_list.sort()
#
# # 9. Print the final, sorted list of keys.
# print(missing_list)
#
# # 10. check what files actually does not exists
# for file in missing_list:
#     if file in os.listdir("./data/kegg/pathways/snvs"):
#         print(file)





# # TODO  - check what csv files are missing columns
# import os
# import pandas as pd
# from tqdm import tqdm
#
# # Define the directory containing the CSV files and the reference file
# snvs_dir = "./data/kegg/pathways/snvs"
# reference_file_path = os.path.join(snvs_dir, "hsa:1.csv")
#
# # --- Step 1: Get the column names from the reference file ---
# try:
#     reference_df = pd.read_csv(reference_file_path)
#     reference_columns = set(reference_df.columns)
#     print(f"Successfully read {len(reference_columns)} columns from the reference file: {reference_file_path}")
# except FileNotFoundError:
#     print(f"Error: The reference file was not found at {reference_file_path}")
#     exit()
# except Exception as e:
#     print(f"An error occurred while reading the reference file: {e}")
#     exit()
#
# # --- Step 2: Get a list of all CSV files in the directory ---
# try:
#     all_csv_files = [f for f in os.listdir(snvs_dir) if f.endswith('.csv') and f != "hsa:1.csv"]
#     if not all_csv_files:
#         print("No other CSV files found in the directory to compare against.")
#         exit()
# except FileNotFoundError:
#     print(f"Error: The directory was not found: {snvs_dir}")
#     exit()
#
# # --- Step 3: Check each CSV file for the reference columns ---
# # Initialize a dictionary to store the missing columns and the files that are missing them
# missing_columns_report = {col: [] for col in reference_columns}
#
# print(f"\nChecking {len(all_csv_files)} CSV files for missing columns...")
#
# for csv_file in tqdm(all_csv_files, desc="Processing CSV files"):
#     file_path = os.path.join(snvs_dir, csv_file)
#     try:
#         df = pd.read_csv(file_path, nrows=0)  # Read only the header for efficiency
#         current_columns = set(df.columns)
#
#         # Find which of the reference columns are missing in the current file
#         missing_in_this_file = reference_columns - current_columns
#
#         # For each missing column, add the current file to its list
#         for col in missing_in_this_file:
#             missing_columns_report[col].append(csv_file)
#
#     except Exception as e:
#         print(f"\nCould not process file {csv_file}. Error: {e}")
#
# # --- Step 4: Print the final report ---
# print("\n--- Missing Columns Report ---")
# found_missing = False
# for column, files in missing_columns_report.items():
#     if files:
#         found_missing = True
#         print(f"\nColumn '{column}' is missing in the following {len(files)} files:")
#         for file_name in files:
#             print(f"  - {file_name}")
#
# if not found_missing:
#     print("\nAll checked CSV files contain all the columns from the reference file.")





# # TODO  - creat visual distributions for some pathways
# import os
# import random
# import pickle
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
#
# pathway_dir = "./data/kegg/pathways/objects"
# column = "clinvar_reg_global_prob"
#
# all_pathways = os.listdir(pathway_dir)
#
# # Set the number of pathways to sample
# num_samples = min(10, len(all_pathways))
# some_pathways = random.sample(all_pathways, num_samples)
#
# for path in some_pathways:
#     try:
#         with open(os.path.join(pathway_dir, path), "rb") as file:
#             # Assuming the pickled object is a dictionary where values are paths to CSV files
#             csv_files = pickle.load(file).values()
#             scores = []
#             for csv_file in csv_files:
#                 try:
#                     # Read each CSV file with pandas
#                     df = pd.read_csv(csv_file)
#                     # Check if the column exists
#                     if column in df.columns:
#                         # Extend the scores list with the column's values, dropping missing values
#                         scores.extend(df[column].fillna(1).tolist())
#                     else:
#                         print(f"Column '{column}' not found in {csv_file}")
#                 except Exception as e:
#                     print(f"Error processing file {csv_file}: {e}")
#
#             if scores:
#                 # Create a plot for the current pathway
#                 plt.figure(figsize=(10, 6))
#                 sns.histplot(scores, kde=True, bins=100)
#                 plt.title(f"Distribution of {column} for {os.path.splitext(path)[0]}")
#                 plt.xlabel(column)
#                 plt.ylabel("Frequency")
#                 plt.grid(True)
#                 plt.show()
#             else:
#                 print(f"No data to plot for pathway: {path}")
#
#     except FileNotFoundError:
#         print(f"File not found: {path}")
#     except Exception as e:
#         print(f"An error occurred with pathway {path}: {e}")





# # TODO  - gc content calc
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




# # TODO  - it is the bic curves creator
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