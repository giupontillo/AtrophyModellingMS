#!/bin/bash

#SBATCH --job-name=nodalstress             #a name for your job
#SBATCH --mem=10G                     	 #max memory per node
#SBATCH --partition=luna-cpu-long        #using luna short queue
#SBATCH --cpus-per-task=8                #max CPU cores per process
#SBATCH --time=07-00:00:00               #time limit (DD-HH:MM)
#SBATCH --nice=2000                  
#SBATCH --qos=anw                     
#SBATCH --output=slurm_logs/slurm-%x.%j_%A_%a.out

#----------------------------------------------------------------------
#               run epicenter mapping server
#----------------------------------------------------------------------

#author: m.barrantescepas@amsterdamumc.nl
#last update: 30 January 2025


ml Anaconda3 
conda activate /data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling

echo "Starting processing..."

python /data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/scripts/02.indv/NodalStress_Indv.py

conda deactivate
#----------------------------------------------------------------------

