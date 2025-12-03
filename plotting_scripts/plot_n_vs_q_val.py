import pandas as pd
import glob
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import gaussian_kde

# ==========================================
# 1. SETUP & STYLE CONFIGURATION ("NATURE" STYLE)
# ==========================================

# Nature guidelines:
# - Single column width: ~89mm (3.5 inches)
# - Double column width: ~183mm (7.2 inches)
# - Font: Arial or Helvetica
# - Font size: 7pt - 8pt for axis tick labels, slightly larger for titles

mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
mpl.rcParams['font.size'] = 8
mpl.rcParams['axes.labelsize'] = 9
mpl.rcParams['axes.titlesize'] = 10
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['axes.linewidth'] = 0.8 # Thinner spines
mpl.rcParams['xtick.major.width'] = 0.8
mpl.rcParams['ytick.major.width'] = 0.8
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['legend.frameon'] = False

# ==========================================
# 2. DATA LOADING
# ==========================================

folder_path = "/results_and_graphs/scores/clinvar_reg_dis_ordered_prob-kl_divergence"
file_pattern = os.path.join(folder_path, "*.csv")
files = glob.glob(file_pattern)

print(f"Found {len(files)} files. Loading data...")

data_frames = []

for file in files:
    try:
        # Read only necessary columns to save memory if files are huge
        print(f"Loading data from {file}...")
        df_temp = pd.read_csv(file, usecols=['n', 'q_value'])
        df_temp = df_temp.dropna(subset=['n', 'q_value'])
        data_frames.append(df_temp)
    except Exception as e:
        print(f"Error reading {os.path.basename(file)}: {e}")

if not data_frames:
    raise ValueError("No data found. Please check the path.")

# Concatenate all data
df = pd.concat(data_frames, ignore_index=True)

df = df[df['n'] <= 150]
#df = df[df['q_value'] <= 0.01]

x = df['n']
y = df['q_value']

print(f"Total data points: {len(df)}")

# ==========================================
# 3. KERNEL DENSITY ESTIMATION
# ==========================================
print("Calculating Gaussian Kernel Density...")

# Stack data for KDE
xy = np.vstack([x, y])

# Calculate density (z)
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
# (on top) for better visualization
idx = z.argsort()
x, y, z = x.iloc[idx], y.iloc[idx], z[idx]

# ==========================================
# 4. PLOTTING
# ==========================================

# Figure size: 3.5 inches is standard for single-column figures in journals
fig, ax = plt.subplots(figsize=(3.5, 3.0))

# Scatter plot with density coloring
# 'c' is the color array (density), 'cmap' is the colormap
scatter = ax.scatter(x, y, c=z, s=10, cmap='viridis',
                     edgecolor='none', alpha=0.8)

# Optional: Add a colorbar to indicate density (comment out if not needed)
cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('Kernel Density', rotation=270, labelpad=15)
cbar.outline.set_visible(False)
cbar.ax.tick_params(size=0) # remove ticks from colorbar

# ==========================================
# 5. BEAUTIFICATION
# ==========================================

# Remove top and right spines (classic scientific look)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Labels
ax.set_xlabel('Sample Size (n)')
ax.set_ylabel('q-value')

# Optional: Log scale?
# q-values are often clustered near zero. If your graph looks like an L-shape,
# uncomment the line below:
# ax.set_yscale('log')

# Title (Keep it simple or remove for publication as it usually goes in the caption)
ax.set_title('Pathway Analysis Scores', pad=10)

plt.tight_layout()

# ==========================================
# 6. SAVING
# ==========================================

output_path = os.path.join("/results_and_graphs", "aggregated_n_vs_pvalue_density_q_value.pdf")
# Save as PDF (vector graphics) is preferred for publication
# Also save a high-res PNG for quick viewing
plt.savefig(output_path.replace('.pdf', '.png'), format='png', dpi=300, bbox_inches='tight')

print(f"Plot saved to: {output_path}")
plt.show()