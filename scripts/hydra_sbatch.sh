#!/bin/bash
#SBATCH --job-name=hydra
#SBATCH --partition=luna-cpu-long
#SBATCH --qos=anw-cpu
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=16
#SBATCH --time=5-00:00:00
#SBATCH --nice=2000

# specify environmental variables
module load Anaconda3
conda activate /home/radv/gpontillo/python-env/mlni

python hydra.py
python hydra414.py