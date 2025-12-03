import os
import sys
import pandas as pd
import argparse

# --- Setup: Correctly set the working directory to the project root ---
# This ensures that relative paths in `definitions.py` will work correctly.
back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(back_dir)

# Now that the CWD is correct, we can import our project modules
from definitions import SCORES_RESULTS_PATH
from statistical_models import PathwayAtlasResults


def analyze_and_save_cancer_results(cancer_file_path: str, score_to_analyze: str, scoring_system: str):
    """
    Runs the full analysis for a single cancer file and saves the results
    as a DataFrame to the appropriate directory.
    """
    # 1. Validate that the cancer file exists
    if not os.path.isfile(cancer_file_path):
        print(f"Error: Cancer file not found at path: {cancer_file_path}", file=sys.stderr)
        sys.exit(1)

    print(f"--- Starting analysis for: {os.path.basename(cancer_file_path)} ---")
    print(f"Score Type: {score_to_analyze}")
    print(f"Distance Metric: {scoring_system}")

    # 2. Prepare the output directory structure
    # The sub-folder name is a combination of the score and system
    output_subdir_name = f"{score_to_analyze}-{scoring_system}"
    output_dir = os.path.join(SCORES_RESULTS_PATH, output_subdir_name)
    os.makedirs(output_dir, exist_ok=True)
    print(f"Results will be saved in: {output_dir}")

    # 3. Run the analysis
    resulter = PathwayAtlasResults()
    cancer_results_dict = resulter.calculate_distances_single_cancer(
        cancer_file_path=cancer_file_path,
        score_to_analyze=score_to_analyze,
        scoring_system=scoring_system
    )

    # 4. Convert the results to a DataFrame and save
    if not cancer_results_dict:
        print("Warning: No results were generated. This might happen if no relevant mutations were found.")
        return

    try:
        # Convert the nested dictionary to a DataFrame
        # The dictionary format is: {'hsa00010.pickle': {'distance': 0.1, 'n': 50}}
        df = pd.DataFrame.from_dict(cancer_results_dict, orient='index')
        df.reset_index(inplace=True)
        df.rename(columns={'index': 'pathway_name'}, inplace=True)

        # Clean up the pathway name by removing the file extension
        df['pathway_name'] = df['pathway_name'].str.replace('.pickle', '', regex=False)

        # 5. Construct the output file path and save
        cancer_filename = os.path.basename(cancer_file_path)
        output_csv_path = os.path.join(output_dir, cancer_filename)

        df.to_csv(output_csv_path, index=False)
        print(f"--- Successfully saved results to: {output_csv_path} ---")

    except Exception as e:
        print(f"An error occurred during DataFrame creation or saving: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """
    Parses command-line arguments to run the single-cancer analysis.
    """
    parser = argparse.ArgumentParser(description="Run pathway distance analysis for a single cancer type.")
    parser.add_argument('--cancer_file', type=str, required=True, help='Full path to the cancer mutations CSV file.')
    parser.add_argument('--score_type', type=str, required=True,
                        choices=['clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob'],
                        help='The score type to analyze.')
    parser.add_argument('--scoring_system', type=str, required=True, choices=['kl_divergence', 'wasserstein', 'dw_distance'],
                        help='The distance metric to use.')

    args = parser.parse_args()

    analyze_and_save_cancer_results(
        cancer_file_path=args.cancer_file,
        score_to_analyze=args.score_type,
        scoring_system=args.scoring_system
    )


if __name__ == "__main__":
    main()