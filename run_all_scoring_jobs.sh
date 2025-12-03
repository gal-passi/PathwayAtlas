#!/bin/bash

# This script submits a Slurm job array for each combination of score type and scoring system.

echo "--- Starting Pathway Atlas Analysis Submission ---"

cd /cs/labs/dina/ophirmil12/PathwayAtlas/

# Create the directory for Slurm output files if it doesn't exist
mkdir -p slurm_out

# Define the combinations to run
SCORE_TYPES=("clinvar_reg_dis_ordered_prob")
#SCORE_TYPES=("clinvar_reg_dis_ordered_prob" "clinvar_reg_global_prob")
SCORING_SYSTEMS=("dw_distance")
#SCORING_SYSTEMS=("kl_divergence" "wasserstein")

# Get the total number of cancer files to set the array size
CANCER_DATA_DIR="/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers"
NUM_CANCER_FILES=$(find "$CANCER_DATA_DIR" -type f -name "*.csv" | wc -l)

if [ "$NUM_CANCER_FILES" -eq 0 ]; then
    echo "Error: No cancer files found. Cannot submit jobs."
    exit 1
fi

# Slurm arrays are 0-indexed
ARRAY_MAX=$((NUM_CANCER_FILES - 1))

# Loop through each combination and submit a job array
for score in "${SCORE_TYPES[@]}"; do
    for system in "${SCORING_SYSTEMS[@]}"; do
        echo "Submitting jobs for: Score=${score}, System=${system}"
        sbatch --array=0-${ARRAY_MAX} submit_scoring_job.slurm "$score" "$system"
    done
done

echo "--- All job arrays have been submitted. ---"