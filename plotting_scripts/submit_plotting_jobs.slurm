#!/bin/bash
#SBATCH --job-name=pathway_plotting
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --killable
#SBATCH --ntasks=1
#SBATCH --requeue
#SBATCH --output=slurm_out/plot_%x_%A_%a.out

# --- Arguments: score_type and scoring_system ---
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch submit_plotting_jobs.slurm <score_type> <scoring_system>"
    exit 1
fi

SCORE_TYPE=$1
SCORING_SYSTEM=$2

# Update job name for visibility
scontrol update job $SLURM_JOB_ID JobName=plot_${SCORE_TYPE}_${SCORING_SYSTEM}

# --- Set up file lists ---
# We use the same directory to ensure the array index matches the file list
CANCER_DATA_DIR="/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/cancers"

CANCER_FILES=($(find "$CANCER_DATA_DIR" -type f -name "*.csv"))

if [ ${#CANCER_FILES[@]} -eq 0 ]; then
    echo "Error: No cancer CSV files found in $CANCER_DATA_DIR"
    exit 1
fi

# Get specific file for this task
TASK_ID=$SLURM_ARRAY_TASK_ID
CANCER_FILE_PATH=${CANCER_FILES[$TASK_ID]}

echo "======================================================"
echo "Running PLOTTING for:"
echo "  - Cancer File: $CANCER_FILE_PATH"
echo "  - Score Type: $SCORE_TYPE"
echo "  - Scoring System: $SCORING_SYSTEM"
echo "======================================================"

cd /cs/labs/dina/ophirmil12/PathwayAtlas/
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

# --- Execute the Plotting Script ---
# Note: We are running a small python one-liner to call the specific function
# for JUST this one cancer file. This avoids the overhead of the script listing
# all files again, ensuring the array index controls exactly which file is processed.

python -u -c "
import sys
import os

# Add project root to path
sys.path.append(os.getcwd())

from statistical_models import PathwayAtlasResults

# Config
cancer_csv = '$CANCER_FILE_PATH'
score_type = '$SCORE_TYPE'
scoring_system = '$SCORING_SYSTEM'
pickle_dir = '/cs/labs/dina/ophirmil12/PathwayAtlas/data/kegg/pathways/objects'
scores_dir = '/cs/labs/dina/ophirmil12/PathwayAtlas/data/kegg/pathways/scores'
results_base_dir = '/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/scores'

# Construct path to the specific results CSV
cancer_filename = os.path.basename(cancer_csv)
results_dir = os.path.join(results_base_dir, f'{score_type}-{scoring_system}')
stats_csv = os.path.join(results_dir, cancer_filename)

if not os.path.exists(stats_csv):
    print(f'Skipping: Results file not found at {stats_csv}')
    sys.exit(0)

print(f'Generating plots for {cancer_filename}...')

atlas = PathwayAtlasResults(
    cancer_csvs_path=os.path.dirname(cancer_csv),
    pathways_dicts_path=pickle_dir,
    pathways_scores_path=scores_dir
)

atlas.generate_plots_for_cancer(
    cancer_file_path=cancer_csv,
    stats_csv_path=stats_csv,
    score_to_analyze=score_type
)
"

echo "Plotting job finished with exit code $?"