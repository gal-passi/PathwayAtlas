#!/bin/bash

# Submits Slurm job arrays for generating plots.

echo "--- Starting Plotting Submission ---"

cd /cs/labs/dina/ophirmil12/PathwayAtlas/
mkdir -p slurm_out

# Count files to determine array size
CANCER_DATA_DIR="/cs/labs/dina/ophirmil12/PathwayAtlas/data/cbio/cancers"
NUM_CANCER_FILES=$(find "$CANCER_DATA_DIR" -type f -name "*.csv" | wc -l)

if [ "$NUM_CANCER_FILES" -eq 0 ]; then
    echo "Error: No cancer files found."
    exit 1
fi

ARRAY_MAX=$((NUM_CANCER_FILES - 1))

# Submit jobs

#echo "Submitting PLOTTING jobs for: Score=clinvar_reg_dis_ordered_prob, System=kl_divergence"
#sbatch --array=0-${ARRAY_MAX} submit_plotting_jobs.slurm "clinvar_reg_dis_ordered_prob" "kl_divergence"

echo "Submitting PLOTTING jobs for: Score=clinvar_reg_dis_ordered_prob, System=dw_distance"
sbatch --array=0-${ARRAY_MAX} submit_plotting_jobs.slurm "clinvar_reg_dis_ordered_prob" "dw_distance"

echo "--- All plotting job arrays submitted. ---"