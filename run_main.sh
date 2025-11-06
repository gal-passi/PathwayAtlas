#!/bin/bash
#SBATCH --job-name=emb_cancer
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH -o slurm.out.%A_%a.out
#SBATCH --gres=gg:g4:1

echo "Script Starting..."

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
#python -u embed_cancer_seqs.py       ## need gpu
#python -u disorder_score_cancer.py                              ####### PY FILE NAME #######
#





# give access to bew files added to data
chmod -R g+w ./data

echo "Script Finished."

### #SBATCH --gres=gg:g4:1      # for GPU [g0 or g4, g10... we have access to g4 only]