#!/bin/bash

#SBATCH --job-name=epicenter             #a name for your job
#SBATCH --mem=1G                     	 #max memory per node
#SBATCH --partition=luna-cpu-long        #using luna short queue
#SBATCH --cpus-per-task=2                #max CPU cores per process
#SBATCH --time=07-00:00:00               #time limit (DD-HH:MM)
#SBATCH --nice=4000                  
#SBATCH --qos=anw-cpu                     
#SBATCH --output=slurm_logs_indv/slurm-%x.%j_%A_%a.out
#SBATCH --array=1-1964%50

#----------------------------------------------------------------------
#               run epicenter mapping server
#----------------------------------------------------------------------

#author: m.barrantescepas@amsterdamumc.nl
#last update: 30 January 2025

ml Anaconda3 
conda activate /data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling

curdir=`pwd`
parent_dir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
atrophydir="${parent_dir}/data/atrophy_cross"
csv_file="${atrophydir}/atrophy_baseline_SchaeferSubcortex114Parcels_individuals_wscore.csv"

full_row=$(tail -n +2 "$csv_file" | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
subject_id=$(echo $full_row | cut -d',' -f1)
subject_id=${subject_id//\"/}

echo $subject_id
touch_file="${parent_dir}/output/june/cross_indv_epi_results/${subject_id}.touch"
echo $touch_file			

echo "Processing row: $full_row"

if [ ! -e ${touch_file} ]; then
        touch "$touch_file"
        
	python ${parent_dir}/scripts/03.epicenter-mapping/epicenter-mapping_indv.py --row $full_row

else echo "${subject_id} already exists!; skipping..."
fi

conda deactivate
#----------------------------------------------------------------------

