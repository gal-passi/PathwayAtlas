
# Goal:
# For each sequence from the cancer studies, save a file with the
# logits matrix from the esm embedding

# Core code:
# for file in cancer csvs
#   for seq in lines
#       embed (seq)
#       save logits to file

from utils import *
import os
import pandas as pd
from tqdm import tqdm
from esm import pretrained
import torch

def main():
    # --- Setup ---
    working_dir = "/"
    cancer_csv_dir = "/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers"
    embedding_dir = "/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/emb"

    os.makedirs(embedding_dir, exist_ok=True)

    # --- Efficiently get the list of already processed IDs ---
    print("Finding existing embeddings...")
    existing_ids = {f.replace('.pt', '') for f in os.listdir(embedding_dir) if f.endswith('.pt')}
    print(f"Found {len(existing_ids)} existing embeddings.")

    # --- Model Loading ---
    model, alphabet = pretrained.load_model_and_alphabet(ESM1B_MODEL)
    calculator = ScoringCalculator(model=model, alphabet=alphabet)

    cancer_csv_files = os.listdir(cancer_csv_dir)
    
    print(f"len(cancer_csv_files):  {len(cancer_csv_files)}\n\n\n")

    run_emb(cancer_csv_files, existing_ids, calculator, cancer_csv_dir, embedding_dir)

def run_emb(cancer_csv_files, existing_ids, calculator, cancer_csv_dir, embedding_dir):
    # Use tqdm on the outer loop to track progress through the list of files
    for file_name in tqdm(cancer_csv_files, desc="Processing Cancer Files"):
        if not file_name.endswith(".csv"):
            continue

        file_path = os.path.join(cancer_csv_dir, file_name)

        # --- Step 1: Load Data Safely ---
        try:
            # low_memory=False can help prevent errors with mixed data types in columns
            df = pd.read_csv(file_path, low_memory=False)
        except Exception as e:
            print(f"Could not read or parse {file_name}. Error: {e}")
            continue

        # --- Step 2: Validate DataFrame Structure ---
        if "ReferenceSeq" not in df.columns or "UniprotId" not in df.columns:
            # print(f"Skipping {file_name}: Missing required 'ReferenceSeq' or 'UniprotId' columns.")
            continue

        # --- Step 3: Efficiently Filter the DataFrame ---
        # This is the core optimization. Instead of looping and checking each row,
        # we create "boolean masks" to select all valid rows in a single operation.

        # Condition 1: The sequence must be valid (not null, not the string 'NA', and not an empty string).
        valid_seq_mask = df['ReferenceSeq'].notna() & (df['ReferenceSeq'] != 'NA') & (df['ReferenceSeq'] != '')

        # Condition 2: The Uniprot ID must not be null.
        valid_id_mask = df['UniprotId'].notna()

        # Apply the initial filters to the DataFrame.
        df_filtered = df[valid_seq_mask & valid_id_mask]

        # If no rows meet the basic validity criteria, we can skip this file.
        if df_filtered.empty:
            continue

        # Condition 3: The Uniprot ID must NOT already be in our set of processed IDs.
        # The .isin() method is highly optimized for this check against the existing_ids set.
        unprocessed_mask = ~df_filtered['UniprotId'].isin(existing_ids)
        df_to_process = df_filtered[unprocessed_mask]

        # --- Step 4: Process and Embed Only the Filtered Rows ---
        # If there are no new proteins in this file, this loop will be skipped.
        # We iterate over the small, clean 'df_to_process' DataFrame.
        # .itertuples() is significantly faster than the old .iterrows().
        for row in df_to_process.itertuples():
            uniprot_id = row.UniprotId
            sequence = row.ReferenceSeq

            try:
                logits = calculator.score_all_mutations(sequence)

                # --- Handle cases where no valid logits were produced ---
                if logits.shape[0] == 0:
                    print(
                        f"Skipping {uniprot_id} from {file_name}: No valid ESM logits produced (sequence likely too short/problematic for ESM).")
                    continue  # Skip saving and do not add to existing_ids

                # Save the resulting tensor to its own file
                output_path = os.path.join(embedding_dir, f"{uniprot_id}.pt")
                torch.save(logits, output_path)

                # IMPORTANT: Add the newly processed ID to our live set.
                # This prevents reprocessing it if it appears again in a later file
                # during the same run.
                existing_ids.add(uniprot_id)

            except Exception as e:
                # If a single protein fails, log the error and continue with the rest.
                print(f"Error processing {uniprot_id} from {file_name}: {e}")

if __name__ == "__main__":
    main()