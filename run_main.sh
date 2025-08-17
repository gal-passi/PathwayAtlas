#!/bin/bash
#SBATCH --job-name=main
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --gres=gg:g4:1      # for GPU [g0 or g4...]
#SBATCH --ntasks=1
#SBATCH -o slurm.out.%A_%a.out


# go to PathwayAtlas folder [local git folder - cd .]
cd /cs/labs/dina/ophirmil12/PathwayAtlas/

# startup venv
source ./venv/bin/activate

# install all requirements
# pip install --no-cache-dir -r requirements.txt

# export esm model path to local venv
export TORCH_HOME=./venv/torch_cache

# run code (output in slurm .out file)
python -u main.py

