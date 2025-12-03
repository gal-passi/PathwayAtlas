#!/bin/bash
#SBATCH --job-name=plots
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH -o slurm.out.%A_%a.out
###########SBATCH --gres=gg:g4:1

echo "Script Starting..."

# IMPORTANT: Replace /cs/labs/dina/ophirmil12/miniforge3 with the actual path to your miniforge3 installation
# This initializes Conda in the current shell session.
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh

# Activate your specific Conda environment
conda activate project_env

# go to PathwayAtlas folder [local git folder - cd .]
cd /cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/

## startup venv       NO NEED WHEN HAVE CONDA
#source ./venv/bin/activate

# install all requirements
#pip install --no-cache-dir -r requirements.txt

# export esm model path to local venv
export TORCH_HOME=./torch_cache




echo "Starting python code..."

# run code (output in slurm .out file)                              ####### PY FILE NAME #######
#python -u ./scripts/embed_cancer_seqs.py       ## need gpu
#python -u ./scripts/disorder_score_cancer.py
#python -u ./scripts/fill_cancer_files.py
#python -u ./scripts/n_scores_threshold.py
python -u ./scripts/plot_n_vs_q_val_by_cancer.py
#python -u ./scripts/plot_cancer_count_per_pathway_significant.py

# give access to bew files added to data
############################################################## chmod -R g+w ./data

echo "Script Finished."
