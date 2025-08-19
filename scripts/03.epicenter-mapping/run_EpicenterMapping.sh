#!/bin/bash

#SBATCH --job-name=epicenter             #a name for your job
#SBATCH --mem=4G                     	 #max memory per node
#SBATCH --partition=luna-cpu-long        #using luna short queue
#SBATCH --cpus-per-task=4                #max CPU cores per process
#SBATCH --time=03-10:00:00               #time limit (DD-HH:MM)
#SBATCH --nice=4000                  
#SBATCH --qos=anw-cpu                     
#SBATCH --output=slurm_logs/slurm-%x.%j_%A_%a.out

#----------------------------------------------------------------------
#               run epicenter mapping server
#----------------------------------------------------------------------

#author: m.barrantescepas@amsterdamumc.nl
#last update: 30 January 2025


ml Anaconda3 
conda activate /data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling

echo "Starting processing..."

python /data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/scripts/epicenter-mapping/EpicenterMapping_long.py

conda deactivate
#----------------------------------------------------------------------

