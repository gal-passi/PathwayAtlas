# import os
# import sys
# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns
# import itertools
#
# # --- Setup: Correctly set the working directory to the project root ---
# # This ensures that relative paths in `definitions.py` will work correctly.
# # Note: __file__ is not defined in some interactive environments.
# # If running this cell by cell, you might need to set the path manually.
# try:
#     back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
#     sys.path.append(back_dir)
# except NameError:
#     print("Running in an environment where __file__ is not defined. Assuming project root is correctly configured.")
#
#
# def find_csv_files(directory):
#     """Find all CSV files in a directory and its subdirectories."""
#     csv_files = []
#     for root, _, files in os.walk(directory):
#         for file in files:
#             if file.endswith('.csv'):
#                 csv_files.append(os.path.join(root, file))
#     return csv_files
#
#
# def process_and_plot_data_efficiently():
#     """
#     Efficiently processes large CSV files by reading only necessary columns and
#     creates scatter plots of specified data combinations.
#     """
#     # Define the directory to search for CSV files and the columns of interest
#     search_directory = '/cs/labs/dina/ophirmil12
#     /PathwayAtlas/data/cbio/cancers'
#     columns_to_plot = [
#         'disorder_score', 'is_disordered', 'esm_log_probs',
#         'clinvar_reg_dis_ordered_prob', 'clinvar_reg_global_prob'
#     ]
#
#     # Define the output directory for the plots
#     output_directory = './results_and_graphs/plots_corr'
#
#     # Create the output directory if it doesn't exist
#     os.makedirs(output_directory, exist_ok=True)
#     print(f"Output directory '{output_directory}' is ready.")
#
#     # Find all CSV files
#     csv_files = find_csv_files(search_directory)
#     if not csv_files:
#         print(f"No CSV files found in '{search_directory}'.")
#         return
#
#     print(f"Found {len(csv_files)} CSV files. Processing files efficiently to conserve memory...")
#
#     # List to hold small, pre-processed DataFrames
#     all_data = []
#
#     for i, file_path in enumerate(csv_files):
#         try:
#             # OPTIMIZATION 1: Only read the required columns into memory using `usecols`.
#             # This is the most significant memory-saving step.
#             df = pd.read_csv(file_path, usecols=columns_to_plot)
#
#             # OPTIMIZATION 2: Drop rows with nulls immediately to reduce the DataFrame's size.
#             df.dropna(inplace=True)
#
#             # Add the cleaned, smaller DataFrame to our list
#             if not df.empty:
#                 all_data.append(df)
#
#             # Provide progress feedback
#             print(f"({i + 1}/{len(csv_files)}) Processed {os.path.basename(file_path)}: Found {len(df)} valid rows.")
#
#         except ValueError as ve:
#             # This error often occurs if a specified column in `usecols` is not in the CSV
#             print(
#                 f"Skipping {os.path.basename(file_path)}: It may be missing one of the required columns. Details: {ve}")
#         except Exception as e:
#             print(f"Error reading or processing {file_path}: {e}")
#
#     if not all_data:
#         print("No data to plot after processing all files.")
#         return
#
#     # OPTIMIZATION 3: Concatenate only once, after processing all files.
#     # The total memory is now determined by the final size of the cleaned data, not the raw data.
#     print("\nConcatenating all processed data...")
#     cleaned_df = pd.concat(all_data, ignore_index=True)
#     print(f"Final combined data has {len(cleaned_df)} rows.")
#
#     # Free up memory by deleting the list of dataframes
#     del all_data
#
#     if cleaned_df.empty:
#         print("No data available for plotting after cleaning.")
#         return
#
#     # Generate and save scatter plots for all combinations of columns
#     print("Generating scatter plots...")
#     for col1, col2 in itertools.permutations(columns_to_plot, 2):
#         plt.figure(figsize=(10, 6))
#         sns.scatterplot(data=cleaned_df, x=col1, y=col2, alpha=0.5,
#                         s=10)  # Added alpha and s for better visualization of dense data
#         plt.title(f'Scatter Plot of {col1} vs {col2}')
#
#         # Define a safe filename
#         filename = f"scatter_{col1}_vs_{col2}.png"
#         save_path = os.path.join(output_directory, filename)
#
#         plt.savefig(save_path, dpi=450)  # Set DPI for better resolution
#         plt.close()
#         print(f"Saved plot: {save_path}")
#
#     print("\nAll plots have been generated and saved successfully.")
#
#
# # Execute the main function
# process_and_plot_data_efficiently()

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from definitions import CANCERS_PATH

# --- Setup: Correctly set the working directory to the project root ---
# Note: __file__ is not defined in some interactive environments.
# If running this cell by cell, you might need to set the path manually.
try:
    back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(back_dir)
except NameError:
    print("Running in an environment where __file__ is not defined. Assuming project root is correctly configured.")


def find_csv_files(directory):
    """Find all CSV files in a directory and its subdirectories."""
    csv_files = []
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.csv'):
                csv_files.append(os.path.join(root, file))
    return csv_files


def create_colored_scatter_plot():
    """
    Efficiently processes large CSV files to create a single, specific scatter plot.
    """
    # Define the directory to search for CSV files and the columns of interest
    search_directory = CANCERS_PATH
    # We only need these three columns for the final plot
    required_columns = [
        'is_disordered', 'clinvar_reg_dis_ordered_prob', 'esm_log_probs'
    ]

    # Define the output directory for the plot
    output_directory = './results_and_graphs/plots_corr'

    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)
    print(f"Output directory '{output_directory}' is ready.")

    # Find all CSV files
    csv_files = find_csv_files(search_directory)
    if not csv_files:
        print(f"No CSV files found in '{search_directory}'.")
        return

    print(f"Found {len(csv_files)} CSV files. Processing files efficiently to conserve memory...")

    # List to hold small, pre-processed DataFrames
    all_data = []

    for i, file_path in enumerate(csv_files):
        try:
            # OPTIMIZATION: Only read the required columns into memory using `usecols`.
            df = pd.read_csv(file_path, usecols=required_columns)

            # Drop rows with nulls in these specific columns.
            df.dropna(inplace=True)

            if not df.empty:
                all_data.append(df)

            print(f"({i + 1}/{len(csv_files)}) Processed {os.path.basename(file_path)}: Found {len(df)} valid rows.")

        except ValueError as ve:
            print(
                f"Skipping {os.path.basename(file_path)}: It may be missing one of the required columns. Details: {ve}")
        except Exception as e:
            print(f"Error reading or processing {file_path}: {e}")

    if not all_data:
        print("No data to plot after processing all files.")
        return

    # Concatenate all the processed data into a single DataFrame.
    print("\nConcatenating all processed data...")
    cleaned_df = pd.concat(all_data, ignore_index=True)
    print(f"Final combined data has {len(cleaned_df)} rows.")

    # Free up memory
    del all_data

    if cleaned_df.empty:
        print("No data available for plotting after cleaning.")
        return

    # --- Generate and save the single colored scatter plot ---
    print("\nGenerating the colored scatter plot...")
    plt.figure(figsize=(12, 8))

    # Use the 'hue' parameter to color points by the 'is_disordered' column
    sns.scatterplot(
        data=cleaned_df,
        x='clinvar_reg_dis_ordered_prob',
        y='esm_log_probs',
        hue='is_disordered',
        palette='viridis',  # A color-blind friendly palette
        alpha=0.6,
        s=15
    )
    plt.title('esm_log_probs vs. Disordered Region Probabilities (Colored by Disorder Status)')
    plt.xlabel('ClinVar Disordered Region Probability')
    plt.ylabel('esm_log_probs')
    plt.legend(title='Is Disordered')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)  # Add a grid for readability

    # Define a safe filename
    filename = "scatter_clinvar_disorder_probs_esm_log_probs_colored_by_disorder.png"
    save_path = os.path.join(output_directory, filename)

    plt.savefig(save_path, dpi=450, bbox_inches='tight')  # Use bbox_inches to prevent labels from being cut off
    plt.close()

    print(f"Plot saved successfully to: {save_path}")


# Execute the main function
create_colored_scatter_plot()