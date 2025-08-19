#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from enigmatoolbox.plotting import plot_subcortical
import numpy as np

# import values and separate cortical and subcortical 
atrophy_map = np.loadtxt("/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data/atrophy_cross/atrophy_baseline_SchaeferSubcortex114Parcels_group1_cohend.csv", delimiter=',', skiprows=1)
atrophy_cortical = atrophy_map[0:100] # values for cortical regions ordered according to the 100 parcels 7 Networks version of the Schaefer atlas  
atrophy_subcortical = atrophy_map[100:114] # values for subcortical regions ordered alphabetically, with all left hemisphere first followed by right hemisphere 

# map parcels to vertices 
atrophy_fsa5 = parcel_to_surface(atrophy_cortical, 'schaefer_100_fsa5')

# plot 
plot_cortical(array_name=atrophy_fsa5, surface_name="fsa5", size=(800, 400), cmap='RdBu_r', color_bar=True, color_range=(-1.5, 1.5))
plot_subcortical(array_name=atrophy_subcortical, ventricles=False, size=(800, 400), cmap='RdBu_r', color_bar=True, color_range=(-1.5, 1.5))