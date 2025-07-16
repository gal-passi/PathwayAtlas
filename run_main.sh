#!/bin/bash
#SBATCH --job-name=main
#SBATCH --killable
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH -o slurm.out.%A_%a.out


# go to PathwayAtlas folder
cd /sci/labs/itamarsi/ophirmil12/dina_temp/PathwayAtlas

# startup venv
source ./venv/bin/activate

# pip install --no-cache-dir -r requirements.txt

# run code (output in slurm .out file)
python -u main.py

