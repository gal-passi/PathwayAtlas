#!/bin/bash
#SBATCH --job-name=filling_scores
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=90:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --array=0-36
#SBATCH -o slurm.out.%A_%a.out
###########SBATCH --gres=gg:g4:1

DISTANCE_METRIC="kl_divergence"
CANCER_SCORES_DIR="/cs/labs/dina/ophirmil12/PathwayAtlas/results_and_graphs/scores/clinvar_reg_dis_ordered_prob-" + $DISTANCE_METRIC

# Get all CSV files into an array
mapfile -t FILES < <(find "$CANCER_SCORES_DIR" -type f -name "*.csv")

NUM_CANCER_FILES=${#FILES[@]}

if [ "$NUM_CANCER_FILES" -eq 0 ]; then
    echo "Error: No cancer files found. Cannot submit jobs."
    exit 1
fi

# FILE corresponding to this array task ID
CANCER_SCORES_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

echo "Processing file: $CANCER_SCORES_FILE"

# Run your Python script with the file as parameter
python main.py --cancer_file "CANCER_SCORES_FILE" "$DISTANCE_METRIC"

echo "Job completed for file: $CANCER_SCORES_FILE"