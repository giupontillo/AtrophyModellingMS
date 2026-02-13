#!/bin/bash

#SBATCH --job-name=indv	                 #a name for your job
#SBATCH --mem=1G                     	 #max memory per node
#SBATCH --partition=luna-cpu-long        #using luna short queue
#SBATCH --cpus-per-task=2                #max CPU cores per process
#SBATCH --time=07-00:00:00               #time limit (DD-HH:MM)
#SBATCH --nice=4000                  
#SBATCH --qos=anw-cpu                     
#SBATCH --output=slurm_logs_indv/slurm-%x.%j_%A_%a.out
#SBATCH --array=2-1965%50

# Extract the current line (skip header)
CSV="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data/atrophy_cross/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv"
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CSV")

# Parse subject and atrophy values
subject=$(echo "$LINE" | cut -d',' -f1)
atrophy=$(echo "$LINE" | cut -d',' -f2-)

# Create a file named after the subject
#touch "${subject}.done"

# Print info
echo "Processing subject: $subject"
echo "Atrophy values: $atrophy"

# Run the Python script, passing subject and atrophy as arguments
#python3 process_atrophy.py "$subject" "$atrophy"






