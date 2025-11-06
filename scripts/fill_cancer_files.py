from definitions import CANCER_READY_DISORDER_PATH, CANCER_READY_EMBEDDINGS_PATH
from utils import ScoringCalculator

import glob
import os
from tqdm import tqdm

#############
# going through the cancer CSV files and filling the scores columns
# for skipping files that are ready - see below the call for handle_cancer_csv, "recalc_scores"
#############

disorder_files = glob.glob(os.path.join(CANCER_READY_DISORDER_PATH, '*.txt'))
embedding_files = glob.glob(os.path.join(CANCER_READY_EMBEDDINGS_PATH, '*.pt'))

disorder_ids = {os.path.splitext(os.path.basename(f))[0] for f in disorder_files}
embedding_ids = {os.path.splitext(os.path.basename(f))[0] for f in embedding_files}

available_ids = disorder_ids.intersection(embedding_ids)

if not available_ids:
    print("[Error] No matching data files found. Aborting.")
    exit(1)

print(f"Found {len(available_ids)} proteins with complete data.")


# Initialize the calculator once
calculator = ScoringCalculator(None, None, testing=True)

# Define the folder path
cancers_folder = "/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers"

# Get a list of files to process
all_files = [f for f in os.listdir(cancers_folder) if f.endswith('.csv')]

if not all_files:
    print(f"No CSV files found in {cancers_folder}")
else:
    print(f"Found {len(all_files)} files to process...")

    file_iterator = tqdm(all_files, desc="Processing cancer files")

    for filename in file_iterator:
        # Construct the full path to the file
        full_path = os.path.join(cancers_folder, filename)

        # Optional: Update the progress bar with the current file name
        file_iterator.set_postfix_str(filename)

        # Process the file
        calculator.handle_cancer_csv(full_path, available_ids, recalc_scores=False)

print("\nAll files processed.")