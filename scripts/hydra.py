#!/home/radv/gpontillo/python-env/mlni/bin/python

"""
Clustering based on the HYDRA method
=============================================================================================
"""

# import all relevant modules
from mlni.hydra_clustering import clustering

feature_tsv="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data/feature_SchaeferSubcortex114Parcels.tsv"
output_dir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/hydra_SchaeferSubcortex114Parcels"
k_min=1
k_max=5
cv_repetition=100
clustering(feature_tsv, output_dir, k_min, k_max, cv_repetition, n_threads=16)