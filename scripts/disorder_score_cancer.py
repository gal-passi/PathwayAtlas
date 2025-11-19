




    # TODO I suspect that the second line is shortening the sequence, thus giving shorter prediction vector,
    #  and causing mismatch of the disorder scoring in the "print(f"\ndisorder: position {position} out of {len(disorder_scores)}")" error.
    #  This version is putting np.nan in places that dont have a score


import os
import pandas as pd
from tqdm import tqdm
from metapredict import meta
import numpy as np

from definitions import VALID_AA, V3_version_letters, CANCER_CSVS_MUTATIONS, CANCER_READY_DISORDER_PATH


def get_aligned_disorder_scores(sequence):
    """
    Calculates disorder scores and aligns them with the original protein sequence.

    Args:
        sequence (str): The amino acid sequence of the protein, which may contain
                        non-standard characters or gaps.

    Returns:
        tuple: A tuple containing:
               - processed_seq (str): The sequence after metapredict-like cleaning.
               - aligned_scores (numpy.ndarray): An array of disorder scores with the
                 same length as the original sequence, with np.nan at positions
                 corresponding to removed characters.
               Returns (None, None) if the input is invalid.
    """
    if not isinstance(sequence, str) or not sequence:
        return None, None

    # 1. Prepare the sequence for metapredict and create a map
    original_len = len(sequence)
    aligned_scores = np.full(original_len, np.nan)

    # These are the characters metapredict will process
    valid_chars = VALID_AA

    # This is the sequence metapredict will effectively see.
    # Non-standard amino acids are replaced as per metapredict's known behavior.
    processed_seq_list = []
    original_indices = []

    temp_seq = V3_version_letters(sequence)

    for i, char in enumerate(temp_seq):
        if char in valid_chars:
            processed_seq_list.append(char)
            original_indices.append(i)

    processed_seq = "".join(processed_seq_list)

    if not processed_seq:
        return sequence, aligned_scores  # Return the original sequence and an array of NaNs

    # 2. Predict disorder for the cleaned sequence
    disorder_scores = meta.predict_disorder(processed_seq)

    # 3. Use the map to align scores with the original sequence
    for idx, score in zip(original_indices, disorder_scores):
        aligned_scores[idx] = score

    return sequence, aligned_scores


def main():
    """
    Main function to process cancer data files, calculate protein disorder scores,
    and save them to individual files, ensuring alignment with the original sequence.
    """
    # --- Setup ---
    disorder_scores_dir = CANCER_READY_DISORDER_PATH
    cancer_csv_dir = CANCER_CSVS_MUTATIONS

    os.makedirs(disorder_scores_dir, exist_ok=True)

    # --- Efficiently get the list of already processed IDs ---
    print("Finding existing disorder score files...")
    existing_ids = {f.replace('.txt', '') for f in os.listdir(disorder_scores_dir) if f.endswith('.txt')}
    print(f"Found {len(existing_ids)} existing score files.")

    # --- Get and sort the list of input files ---
    try:
        cancer_csv_files = sorted([f for f in os.listdir(cancer_csv_dir) if f.endswith(".csv")])
    except FileNotFoundError:
        print(f"Error: Input directory not found at {cancer_csv_dir}")
        return

    # --- Main Processing Loop ---
    for file_name in tqdm(cancer_csv_files, desc="Processing Cancer Files"):
        file_path = os.path.join(cancer_csv_dir, file_name)

        try:
            df = pd.read_csv(file_path, low_memory=False)
        except Exception as e:
            print(f"Could not read or parse {file_name}. Error: {e}")
            continue

        required_cols = ["ReferenceSeq", "UniprotId"]
        if not all(col in df.columns for col in required_cols):
            continue

        valid_seq_mask = df['ReferenceSeq'].notna() & (df['ReferenceSeq'] != '')
        valid_id_mask = df['UniprotId'].notna()
        df_filtered = df[valid_seq_mask & valid_id_mask]

        if df_filtered.empty:
            continue

        unprocessed_mask = ~df_filtered['UniprotId'].isin(existing_ids)
        df_to_process = df_filtered[unprocessed_mask].drop_duplicates(subset=['UniprotId'])

        for row in df_to_process.itertuples(index=False):
            uniprot_id = row.UniprotId
            sequence = row.ReferenceSeq

            try:
                original_sequence, aligned_scores = get_aligned_disorder_scores(sequence)

                if aligned_scores is None:
                    continue

                output_path = os.path.join(disorder_scores_dir, f"{uniprot_id}.txt")

                with open(output_path, 'w') as f:
                    f.write(original_sequence + '\n')
                    # Convert numpy array to a space-separated string, representing NaNs as 'nan'
                    f.write(' '.join(map(str, aligned_scores)) + '\n')

                existing_ids.add(uniprot_id)

            except Exception as e:
                print(f"Error processing {uniprot_id} from {file_name}: {e}")

    print("\nProcessing complete.")


if __name__ == "__main__":
    main()