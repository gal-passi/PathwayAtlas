#!/bin/bash
#SBATCH --job-name=main
#SBATCH --killable
#SBATCH --requeue
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH -o slurm.out.%A_%a.out


# go to PathwayAtlas folder [local git folder - cd .]
cd /cs/labs/dina/lotem.senderov/PycharmProjects/PathwayAtlas/

# startup venv
source ./venv/bin/activate

# install all requirements
#pip install --no-cache-dir -r requirements.txt

# export esm model path to local venv
#export TORCH_HOME=./venv/torch_cache

# run code (output in slurm .out file)
python -u main.py

### #SBATCH --gres=gg:g4:1      # for GPU [g0 or g4, g10... we have access to g4 only]