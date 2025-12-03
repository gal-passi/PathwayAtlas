import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde

# ==========================================
# 1. SETUP & STYLE CONFIGURATION
# ==========================================
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.linewidth'] = 0.8
mpl.rcParams['legend.frameon'] = False

# ==========================================
# 2. PATH DEFINITIONS
# ==========================================

input_folder = "/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/scores/clinvar_reg_dis_ordered_prob-kl_divergence"
output_folder = "/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/aggregated_n_vs_pvalue_density_by_cancer"

# Create output directory if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

file_pattern = os.path.join(input_folder, "*.csv")
files = glob.glob(file_pattern)

print(f"Found {len(files)} files. Starting processing...")

# ==========================================
# 3. PROCESSING LOOP
# ==========================================

for file_path in files:
    try:
        # Get the cancer name from filename
        file_name = os.path.basename(file_path)
        cancer_name = os.path.splitext(file_name)[0]

        print(f"Processing: {cancer_name}...")

        # ---------------------------
        # Load Data
        # ---------------------------
        df = pd.read_csv(file_path, usecols=['n', 'q_value'])
        df = df.dropna(subset=['n', 'q_value'])

        df = df[df['n'] >= 20]

        if df.empty:
            print(f"  -> Skipping {cancer_name}: Empty file.")
            continue

        # ---------------------------
        # Filter: n <= 90th Percentile
        # ---------------------------
        # Calculate the 90th percentile for 'n' in this specific file
        n_limit = df['n'].quantile(0.90)

        # Apply the filter
        df = df[df['n'] <= n_limit]

        # Skip if not enough data points left for KDE
        if len(df) < 5:
            print(f"  -> Skipping {cancer_name}: Not enough data (n={len(df)}) after filtering.")
            continue

        x = df['n']
        y = df['q_value']

        # ---------------------------
        # Kernel Density Estimation
        # ---------------------------
        xy = np.vstack([x, y])

        try:
            z = gaussian_kde(xy)(xy)
            # Sort points by density (densest on top)
            idx = z.argsort()
            x_plot, y_plot, z_plot = x.iloc[idx], y.iloc[idx], z[idx]
        except Exception as kde_err:
            print(f"  -> KDE Warning for {cancer_name}: {kde_err}. Plotting without density.")
            x_plot, y_plot = x, y
            z_plot = None

        # ---------------------------
        # Plotting
        # ---------------------------
        fig, ax = plt.subplots(figsize=(3.5, 3.0))

        if z_plot is not None:
            scatter = ax.scatter(x_plot, y_plot, c=z_plot, s=10, cmap='viridis',
                                 edgecolor='none', alpha=0.8)
            # Add colorbar
            cbar = plt.colorbar(scatter, ax=ax)
            cbar.set_label('Kernel Density', rotation=270, labelpad=15)
            cbar.outline.set_visible(False)
            cbar.ax.tick_params(size=0)
        else:
            ax.scatter(x_plot, y_plot, c='blue', s=10, edgecolor='none', alpha=0.6)

        # ---------------------------
        # Beautification
        # ---------------------------
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_xlabel('Sample Size (n)')
        ax.set_ylabel('q-value')
        ax.set_title(cancer_name, pad=10)

        # Set X-axis limit to the calculated percentile + small buffer
        ax.set_xlim(left= 20 * 0.9, right=n_limit * 1.05)

        plt.tight_layout()

        # ---------------------------
        # Saving
        # ---------------------------
        save_path_png = os.path.join(output_folder, f"{cancer_name}.png")

        plt.savefig(save_path_png, format='png', dpi=300, bbox_inches='tight')

        plt.close(fig)  # Close figure to free memory

    except Exception as e:
        print(f"Error processing {file_path}: {e}")

print("All processing complete.")