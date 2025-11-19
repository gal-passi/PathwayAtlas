# Goal: collect all the scores for each pathway to a single new CSV file
# collecting from KEGG_PATHWAY_MUTATIONS_PATH using KEGG_PATHWAY_OBJECTS_PATH into PATHWAY_SCORES_PATH

import pandas as pd
import os, sys

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(os.path.dirname(script_dir))

from definitions import KEGG_PATHWAY_MUTATIONS_PATH, KEGG_PATHWAY_OBJECTS_PATH, PATHWAY_SCORES_PATH
import pickle

# Create the output directory if it doesn't exist
os.makedirs(PATHWAY_SCORES_PATH, exist_ok=True)

for pathway_obj_filename in sorted(os.listdir(KEGG_PATHWAY_OBJECTS_PATH)):
    # Construct the full path to the pickle file
    pathway_obj_filepath = os.path.join(KEGG_PATHWAY_OBJECTS_PATH, pathway_obj_filename)

    # Load the pathway dictionary from the pickle file
    try:
        with open(pathway_obj_filepath, 'rb') as f:
            pathway_dict = pickle.load(f)
    except EOFError:
        print(
            f"Warning: Could not load data from '{pathway_obj_filename}'. The file may be empty or corrupted. Skipping.")
        continue  # This command skips the rest of the loop and moves to the next file
    except Exception as e:
        print(f"An unexpected error occurred while reading '{pathway_obj_filename}': {e}. Skipping.")
        continue

    # List to hold DataFrames for the current pathway
    list_of_gene_scores_df = []

    print(f"Processing pathway: {pathway_obj_filename}")

    # Iterate through the gene score CSV filenames in the pathway dictionary
    for gene_scores_filename in pathway_dict.values():
        # Construct the full path to the gene scores CSV file
        if gene_scores_filename is None:
            continue

        try:
            # Open the gene_scores csv file
            gene_scores_df = pd.read_csv(gene_scores_filename)

            # Specify the columns to take
            columns_to_take = ['Ref', 'Alt', 'Protein', 'Variant', 'disorder_score', 'is_disordered', 'esm_log_probs', 'esm_min_max_naive',
                               'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob']

            # Ensure all required columns exist in the DataFrame
            existing_columns = [col for col in columns_to_take if col in gene_scores_df.columns]

            if not existing_columns:
                print(f"Warning: None of the required columns found in {gene_scores_filename}. Skipping this file.")
                continue

            # Take the specified columns
            subset_df = gene_scores_df[existing_columns]

            # Add them to our list of DataFrames
            list_of_gene_scores_df.append(subset_df)

        except FileNotFoundError:
            print(f"Warning: File not found at {gene_scores_filename}. Skipping.")
        except Exception as e:
            print(f"An error occurred while processing {gene_scores_filename}: {e}")

    # Concatenate all the dataframes for the pathway if any were collected
    if list_of_gene_scores_df:
        pathway_scores = pd.concat(list_of_gene_scores_df, ignore_index=True)

        # Create a new CSV filename for the combined scores
        output_csv_name = os.path.splitext(pathway_obj_filename)[0] + '.csv'
        output_csv_path = os.path.join(PATHWAY_SCORES_PATH, output_csv_name)

        # Save the combined scores to the new CSV file in PATHWAY_SCORES_PATH
        pathway_scores.to_csv(output_csv_path, index=False)
        print(f"Successfully saved scores to {output_csv_path}")
    else:
        print(f"No data collected for pathway {pathway_obj_filename}. No CSV file will be created.")