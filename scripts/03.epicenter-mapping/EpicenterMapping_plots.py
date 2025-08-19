#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 22 Jan 2024

#%%

import os 
import numpy as np 
import pandas as pd
from enigmatoolbox.plotting import plot_cortical, plot_subcortical
from enigmatoolbox.utils.parcellation import parcel_to_surface

#%%

def load_csv_to_numpy(file_mapping, base_path):
    "reads multiple csv files into numpy arrays"
    arrays = {}
    for key, filename in file_mapping.items():
        arrays[key] = pd.read_csv(f"{base_path}/{filename}.csv", header=None).to_numpy()
    return arrays

def plot_cortical_epi_func(sfc_ctx_epi, sfc_ctx_epi_p, atlas, output_folder, filename, color_range=(-0.5,0.5), cmap="RdBu"):
    sfc_ctx_epi_p_sig = np.zeros_like(sfc_ctx_epi_p)
    sfc_ctx_epi_p_sig[np.argwhere(sfc_ctx_epi_p < 0.1)] = sfc_ctx_epi[np.argwhere(sfc_ctx_epi_p < 0.1)]
    surf = parcel_to_surface(sfc_ctx_epi_p_sig, atlas) 
    plot_cortical(array_name=surf, surface_name="fsa5", size=(2000,1500), 
                  cmap=cmap, color_bar=True, color_range=color_range, 
                  screenshot=True, transparent_bg=False, 
                  filename=f"{output_folder}/{filename}.png")
    
def plot_subcortical_epi_func(sfc_sctx_epi, sfc_sctx_epi_p, output_folder, filename, color_range=(-0.5,0.5), cmap="RdBu"):
    sfc_sctx_epi=sfc_sctx_epi.squeeze()
    sfc_sctx_epi_p=sfc_sctx_epi_p.squeeze()
    sfc_sctx_epi_p_sig = np.zeros_like(sfc_sctx_epi_p)
    sfc_sctx_epi_p_sig[np.argwhere(sfc_sctx_epi_p < 0.1)] = sfc_sctx_epi[np.argwhere(sfc_sctx_epi_p < 0.1)]    
    plot_subcortical(array_name=sfc_sctx_epi_p_sig, ventricles=False, 
                     size=(2000,1500), cmap=cmap, color_bar=True, color_range=color_range, 
                     screenshot=True, transparent_bg=False,
                     filename=f"{output_folder}/{filename}.png")

#%%
#----------------------------------------------------------------------
# LOAD DATA 
#----------------------------------------------------------------------

cortical_epicenters_data = {
    "fc_ctx_epi_d100": "fc_ctx_epi_d100",
    "fc_ctx_epi_p_d100": "fc_ctx_epi_p_d100",
    "sc_ctx_epi_d100": "sc_ctx_epi_d100",
    "sc_ctx_epi_p_d100": "sc_ctx_epi_p_d100",
    "fc_ctx_epi_d400": "fc_ctx_epi_d400",
    "fc_ctx_epi_p_d400": "fc_ctx_epi_p_d400",
    "sc_ctx_epi_d400": "sc_ctx_epi_d400",
    "sc_ctx_epi_p_d400": "sc_ctx_epi_p_d400"}

subcortical_epicenters_data = {
    "fc_sctx_epi_d100": "fc_sctx_epi_d100",
    "fc_sctx_epi_p_d100": "fc_sctx_epi_p_d100",
    "sc_sctx_epi_d100": "sc_sctx_epi_d100",
    "sc_sctx_epi_p_d100": "sc_sctx_epi_p_d100",
    "fc_sctx_epi_d400": "fc_sctx_epi_d400",
    "fc_sctx_epi_p_d400": "fc_sctx_epi_p_d400",
     "sc_sctx_epi_d400": "sc_sctx_epi_d400",
     "sc_sctx_epi_p_d400": "sc_sctx_epi_p_d400"}

dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
output_folder = f"{dir}/output/june/cross_group_epi_results"

epi_ctx_df = load_csv_to_numpy(cortical_epicenters_data, output_folder)
epi_sctx_df = load_csv_to_numpy(subcortical_epicenters_data, output_folder)

#%%
#----------------------------------------------------------------------
# PLOT RESULTS ON BRAIN SURFACE
#----------------------------------------------------------------------
# Project the results on the surface brain
# Selecting only regions with p < 0.1 (functional epicenters)
plot_folder=f"{dir}/output/june/cross_group_epi_results"
os.makedirs(plot_folder, exist_ok=True)

plot_cortical_epi_func(epi_ctx_df["fc_ctx_epi_d100"], epi_ctx_df["fc_ctx_epi_p_d100"], 'schaefer_100_fsa5', plot_folder, 'plot_fc_ctx_epi_d100', cmap="GyRd_r")
plot_cortical_epi_func(epi_ctx_df["sc_ctx_epi_d100"], epi_ctx_df["sc_ctx_epi_p_d100"], 'schaefer_100_fsa5', plot_folder, 'plot_sc_ctx_epi_d100', cmap="GyBu_r")
plot_cortical_epi_func(epi_ctx_df["fc_ctx_epi_d400"], epi_ctx_df["fc_ctx_epi_p_d400"], 'schaefer_400_fsa5', plot_folder, 'plot_fc_ctx_epi_d400', cmap="GyRd_r")
plot_cortical_epi_func(epi_ctx_df["sc_ctx_epi_d400"], epi_ctx_df["sc_ctx_epi_p_d400"], 'schaefer_400_fsa5', plot_folder, 'plot_sc_ctx_epi_d400', cmap="GyBu_r")

#%%
plot_subcortical_epi_func(epi_sctx_df["fc_sctx_epi_d100"], epi_sctx_df["fc_sctx_epi_p_d100"], plot_folder, 'plot_fc_sctx_epi_d100', cmap="GyRd_r")
plot_subcortical_epi_func(epi_sctx_df["sc_sctx_epi_d100"], epi_sctx_df["sc_sctx_epi_p_d100"], plot_folder, 'plot_sc_sctx_epi_d100', cmap="GyBu_r")
plot_subcortical_epi_func(epi_sctx_df["fc_sctx_epi_d400"], epi_sctx_df["fc_sctx_epi_p_d400"], plot_folder, 'plot_fc_sctx_epi_d400', cmap="GyRd_r")
plot_subcortical_epi_func(epi_sctx_df["sc_sctx_epi_d400"], epi_sctx_df["sc_sctx_epi_p_d400"], plot_folder, 'plot_sc_sctx_epi_d400', cmap="GyBu_r")
