import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import math

# ==========================================
# 1. SETUP & STYLE
# ==========================================
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 7
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['legend.frameon'] = False

# ==========================================
# 2. DATA LOADING
# ==========================================
folder_path = "/results_and_graphs/scores/clinvar_reg_dis_ordered_prob-kl_divergence"
file_pattern = os.path.join(folder_path, "*.csv")
files = glob.glob(file_pattern)

print(f"Found {len(files)} files.")

all_data = []

# Columns to load based on your request
cols_to_use = [
    'pathway_name',
    'p_value',
    'q_value',
    'significant_0.05',
    'significant_0.01'
    # TODO , significant_0.001'
]

for file in files:
    try:
        # We try to load the specific columns.
        # Note: If significant_0.05 is missing in a specific CSV, pandas might raise an error.
        # We use a small check to ensure robustness.

        # Option A: Read just the header first to see what's available (safest)
        header = pd.read_csv(file, nrows=0).columns.tolist()
        # Intersection of what we want and what exists
        actual_cols = [c for c in cols_to_use if c in header]

        df_temp = pd.read_csv(file, usecols=actual_cols)
        df_temp['source_file'] = os.path.basename(file)
        all_data.append(df_temp)
    except Exception as e:
        print(f"Error reading {os.path.basename(file)}: {e}")

if not all_data:
    raise ValueError("No data found.")

full_df = pd.concat(all_data, ignore_index=True)
total_cancers = len(files)


# ==========================================
# 3. ANALYSIS & PLOTTING FUNCTION
# ==========================================

def plot_recurrence_faceted(dataframe, threshold, output_folder, items_per_row=90):
    print(f"\n--- Processing Threshold {threshold} ---")

    col_name = f"significant_{threshold}"

    # --- FILTERING LOGIC ---
    if col_name in dataframe.columns:
        # Use the pre-calculated column (handling string 'True' vs boolean True)
        # We convert to string first to handle both boolean and string types safely
        filtered = dataframe[dataframe[col_name].astype(str) == "True"].copy()
        method_used = f"Column '{col_name}'"
    else:
        # Fallback for 0.001 or if column is missing: Calculate manually
        print(f"Column '{col_name}' not found. Filtering by p_value < {threshold}...")
        filtered = dataframe[dataframe['p_value'] < threshold].copy()
        method_used = "p_value calculation"

    if filtered.empty:
        print(f"No pathways found for threshold {threshold}")
        return

    print(f"Using {method_used}. Found {len(filtered)} significant entries.")

    # Count distinct files per pathway
    counts = filtered.groupby('pathway_name')['source_file'].nunique()
    counts = counts.sort_values(ascending=False)

    num_pathways = len(counts)
    num_rows = math.ceil(num_pathways / items_per_row)

    print(f"Total recurrent pathways: {num_pathways}")
    print(f"Generating 1 figure with {num_rows} rows...")

    # Dynamic Height
    fig_height = 3.5 * num_rows
    fig, axes = plt.subplots(nrows=num_rows, ncols=1,
                             figsize=(12, fig_height),
                             sharey=True)

    if num_rows == 1:
        axes = [axes]

    for i, ax in enumerate(axes):
        start = i * items_per_row
        end = min((i + 1) * items_per_row, num_pathways)
        subset = counts.iloc[start:end]

        ax.bar(subset.index, subset.values, color='#4c72b0', width=0.7)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(0, total_cancers + 2)
        ax.set_xlim(-1, len(subset))
        ax.set_ylabel('Count')

        # X-axis
        ax.set_xticklabels(subset.index, rotation=90, ha='center', fontsize=6)
        ax.set_title(f"Rank {start + 1} - {end}", fontsize=9, pad=3)

    fig.suptitle(f"Recurrent Pathways (Threshold {threshold})", fontsize=12, y=1.00)
    plt.tight_layout()

    # Save PNG
    filename = f"recurrence_faceted_thresh_{threshold}.png"
    output_path = os.path.join(output_folder, filename)
    plt.savefig(output_path, format='png', dpi=300, bbox_inches='tight')

    plt.close(fig)
    print(f"Saved: {filename}")


# ==========================================
# 4. EXECUTION
# ==========================================

output_dir = "/results_and_graphs/plots_cancer_count_per_pathway_significant"
os.makedirs(output_dir, exist_ok=True)

# We include 0.05, 0.01 (using columns) and 0.001 (fallback to p-value)
thresholds_to_plot = [0.05, 0.01]       # TODO run also with 0.001 when have the q-val column

for thresh in thresholds_to_plot:
    plot_recurrence_faceted(full_df, thresh, output_dir, items_per_row=100)

print("\nAll plots generated successfully.")