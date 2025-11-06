import os
import pandas as pd
from tqdm import tqdm
from metapredict import meta
import numpy as np
from definitions import DISORDERED_THRESHOLD



def get_disorder_scores(sequence):
    """
    Calculates the disorder scores for a single protein sequence.
    
    Args:
        sequence (str): The amino acid sequence of the protein.

    Returns:
        numpy.ndarray: An array of disorder scores for each amino acid, or None if input is invalid.
    """
    # 1. Ensure sequence data is a valid string
    if not isinstance(sequence, str) or not sequence:
        return None

    # 2. Clean up the protein sequence
    # Replace non-standard amino acid characters and remove whitespace/symbols
    processed_seq = sequence.replace('B', 'N').replace('U', 'C').replace('X', 'G').replace('Z', 'Q')
    processed_seq = processed_seq.replace(' ', '').replace('*', '').replace('-', '')

    if not processed_seq:
        return None, None

    # 3. Predict disorder for the entire sequence
    # meta.predict_disorder returns a numpy array of scores for each amino acid
    disorder_scores = meta.predict_disorder(processed_seq)

    return processed_seq, disorder_scores


def main():
    """
    Main function to process cancer data files, calculate protein disorder scores,
    and save them to individual files.
    """
    # --- Setup ---
    disorder_scores_dir = "/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/disorder_scores"
    cancer_csv_dir = "/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers"

    # Create the output directory if it doesn't already exist
    os.makedirs(disorder_scores_dir, exist_ok=True)

    # --- Efficiently get the list of already processed IDs ---
    print("Finding existing disorder score files...")
    # We assume the file is named 'UniprotId.txt' and strip the extension
    existing_ids = {f.replace('.txt', '') for f in os.listdir(disorder_scores_dir) if f.endswith('.txt')}
    print(f"Found {len(existing_ids)} existing score files.")

    # --- Get and sort the list of input files ---
    try:
        cancer_csv_files = sorted([f for f in os.listdir(cancer_csv_dir) if f.endswith(".csv")])
    except FileNotFoundError:
        print(f"Error: Input directory not found at {cancer_csv_dir}")
        return

    # --- Main Processing Loop ---
    # Use tqdm on the outer loop to track progress through the cancer files
    for file_name in tqdm(cancer_csv_files, desc="Processing Cancer Files"):
        file_path = os.path.join(cancer_csv_dir, file_name)

        # 1. Load Data Safely
        try:
            df = pd.read_csv(file_path, low_memory=False)
        except Exception as e:
            print(f"Could not read or parse {file_name}. Error: {e}")
            continue

        # 2. Validate DataFrame Structure
        required_cols = ["ReferenceSeq", "UniprotId"]
        if not all(col in df.columns for col in required_cols):
            # print(f"Skipping {file_name}: Missing one or more required columns.")
            continue
        
        # 3. Efficiently Filter the DataFrame
        # Condition 1: The sequence must be a valid, non-empty string.
        valid_seq_mask = df['ReferenceSeq'].notna() & (df['ReferenceSeq'] != '')
        
        # Condition 2: The Uniprot ID must not be null.
        valid_id_mask = df['UniprotId'].notna()

        df_filtered = df[valid_seq_mask & valid_id_mask]

        if df_filtered.empty:
            continue

        # Condition 3: The Uniprot ID must NOT already be in our set of processed IDs.
        unprocessed_mask = ~df_filtered['UniprotId'].isin(existing_ids)
        df_to_process = df_filtered[unprocessed_mask].drop_duplicates(subset=['UniprotId'])


        # 4. Process and Save Scores for Only the Filtered Rows
        # .itertuples() is significantly faster than .iterrows()
        for row in df_to_process.itertuples(index=False):
            uniprot_id = row.UniprotId
            sequence = row.ReferenceSeq

            try:
                # Calculate the disorder scores for the sequence
                processed_seq, disorder_scores = get_disorder_scores(sequence)

                if disorder_scores is None:
                    # print(f"Skipping {uniprot_id} from {file_name}: Invalid sequence after cleaning.")
                    continue

                # Define the output path for the scores file
                output_path = os.path.join(disorder_scores_dir, f"{uniprot_id}.txt")

                # Save the scores to a 2-line text file
                with open(output_path, 'w') as f:
                    f.write(processed_seq + '\n')
                    # Convert numpy array to a space-separated string and write to file
                    f.write(' '.join(map(str, disorder_scores)) + '\n')

                # IMPORTANT: Add the newly processed ID to our live set.
                # This prevents reprocessing it if it appears again in a later file
                # during the same run.
                existing_ids.add(uniprot_id)

            except Exception as e:
                # If a single protein fails, log the error and continue with the rest.
                print(f"Error processing {uniprot_id} from {file_name}: {e}")

    print("\nProcessing complete.")


if __name__ == "__main__":
    main()