#!/bin/bash
#SBATCH --job-name=main
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=18:00:00
#SBATCH --gres=gg:g4:1
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH -o slurm.out.%A_%a.out



# IMPORTANT: Replace /cs/labs/dina/ophirmil12/miniforge3 with the actual path to your miniforge3 installation
# This initializes Conda in the current shell session.
source /cs/labs/dina/ophirmil12/miniforge3/etc/profile.d/conda.sh

# Activate your specific Conda environment
conda activate project_env



# go to PathwayAtlas folder [local git folder - cd .]
cd /cs/labs/dina/ophirmil12/PathwayAtlas/

## startup venv       NO NEED WHEN HAVE CONDA
#source ./venv/bin/activate

# install all requirements
# pip install --no-cache-dir -r requirements.txt

# export esm model path to local venv
export TORCH_HOME=./torch_cache

# run code (output in slurm .out file)
python -u main.py

### #SBATCH --gres=gg:g4:1      # for GPU [g0 or g4, g10... we have access to g4 only]