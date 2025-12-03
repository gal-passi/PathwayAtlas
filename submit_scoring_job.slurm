#!/bin/bash
#SBATCH --job-name=pathway_analysis
#SBATCH --time=90:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --requeue

# --- Dynamically set output file and job name based on arguments ---
#SBATCH --output=slurm_out/%x_%A_%a.out  # %x = job name, %A = job ID, %a = task ID

# --- This script expects two arguments: score_type and scoring_system ---
if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch submit_analysis_jobs.slurm <score_type> <scoring_system>"
    exit 1
fi

SCORE_TYPE=$1
SCORING_SYSTEM=$2

# Update the job name to be descriptive for the queue (squeue)
# This command must be run before the job logic starts
scontrol update job $SLURM_JOB_ID JobName=${SCORE_TYPE}_${SCORING_SYSTEM}

# --- Set up the environment and file lists ---
CANCER_DATA_DIR="/cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/data/cbio/cancers"

echo "Listing cancer files from: $CANCER_DATA_DIR"
# Create a bash array of all cancer CSV files
CANCER_FILES=($(find "$CANCER_DATA_DIR" -type f -name "*.csv"))

# Check if files were found
if [ ${#CANCER_FILES[@]} -eq 0 ]; then
    echo "Error: No cancer CSV files found in $CANCER_DATA_DIR"
    exit 1
fi

# Get the specific cancer file for this job array task
TASK_ID=$SLURM_ARRAY_TASK_ID
CANCER_FILE_PATH=${CANCER_FILES[$TASK_ID]}

# --- Print job information to the output file for debugging ---
echo "======================================================"
echo "Slurm Job ID: $SLURM_JOB_ID"
echo "Slurm Array Task ID: $TASK_ID"
echo "Running analysis for:"
echo "  - Cancer File: $CANCER_FILE_PATH"
echo "  - Score Type: $SCORE_TYPE"
echo "  - Scoring System: $SCORING_SYSTEM"
echo "======================================================"

# --- Go to the project directory and activate environment ---
cd /cs/labs/dina/ophirmil12/PathwayAtlas/
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh
conda activate project_env

# --- Execute the Python worker script ---
python -u scripts/run_single_cancer_analysis.py \
    --cancer_file "$CANCER_FILE_PATH" \
    --score_type "$SCORE_TYPE" \
    --scoring_system "$SCORING_SYSTEM"

echo "Job finished with exit code $?"