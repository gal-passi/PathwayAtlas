
import glob
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Setup: Correctly set the working directory to the project root ---
# Note: __file__ is not defined in some interactive environments.
# If running this cell by cell, you might need to set the path manually.
try:
    back_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    sys.path.append(back_dir)
except NameError:
    print("Running in an environment where __file__ is not defined. Assuming project root is correctly configured.")


# Define paths
input_root_dir = "./results_and_graphs/scores"
output_dir = "./results_and_graphs/n_scores_threshold"

# Create the output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Get a list of all subdirectories in the input folder
subfolders = [f.path for f in os.scandir(input_root_dir) if f.is_dir()]

print(f"Found {len(subfolders)} folders to process.")

# Loop through each folder
for folder_path in subfolders:
    folder_name = os.path.basename(folder_path)
    print(f"Processing folder: {folder_name}...")

    # Find all CSV files in the current folder
    csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

    if not csv_files:
        print(f" - No CSV files found in {folder_name}. Skipping.")
        continue

    data_frames = []

    # Read each CSV and extract specific columns
    for file in csv_files:
        try:
            df = pd.read_csv(file)
            # Check if columns exist (case-sensitive)
            if 'distance' in df.columns and 'n' in df.columns:
                data_frames.append(df[['distance', 'n']])
            else:
                # Try case-insensitive check if needed, or warn
                print(f" - Warning: Columns 'distance' and 'n' not found in {os.path.basename(file)}")
        except Exception as e:
            print(f" - Error reading {file}: {e}")

    # If we have data, combine it and plot
    if data_frames:
        combined_df = pd.concat(data_frames, ignore_index=True)

        # Calculate N (Total number of samples/rows)
        N = len(combined_df)

        # Initialize the plot
        plt.figure(figsize=(10, 6))
        sns.set_theme(style="whitegrid")

        # Plot Scatter + Kernel Regression Line
        # lowess=True creates a locally weighted regression (smooth trend line)
        # scatter_kws sets the transparency and size of dots
        sns.regplot(
            data=combined_df,
            x='n',
            y='distance',
            lowess=True,
            line_kws={"color": "red", "linewidth": 2},
            scatter_kws={"alpha": 0.5, "s": 20}
        )

        # Add titles and labels
        plt.title(f"Folder: {folder_name} | Total Samples N={N}", fontsize=14)
        plt.xlabel("n")
        plt.ylabel("Distance")

        plt.xlim(0, 150)

        # Save the plot
        save_path = os.path.join(output_dir, f"{folder_name}_scatter_0_150.png")
        plt.savefig(save_path, dpi=300)
        plt.close()  # Close plot to free memory

        print(f" - Graph saved to: {save_path}")
    else:
        print(f" - No valid data extracted for {folder_name}.")

print("Processing complete.")