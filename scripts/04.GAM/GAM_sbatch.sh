#!/bin/bash
#SBATCH --job-name=GAM
#SBATCH --partition=luna-cpu-long
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=4-00:00:00

module load R

Rscript GAM_modelling.R
