#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 16:57:24 2024

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 09 Dec 2024
"""

#%%
from enigmatoolbox.datasets import load_sc, load_fc
from nilearn import plotting 
import numpy as np 
from enigmatoolbox.plotting import plot_cortical
from enigmatoolbox.plotting import plot_subcortical
from enigmatoolbox.utils.parcellation import parcel_to_surface 
from enigmatoolbox.datasets import load_summary_stats, load_example_data
from enigmatoolbox.utils.useful import reorder_sctx, zscore_matrix #mega version
from enigmatoolbox.permutation_testing import spin_test, shuf_test
import matplotlib.pyplot as plt

#%%
# load cortical and subcortical connectivity matrices
sc_ctx, sc_ctx_labels, sc_sctx, sc_sctx_labels = load_sc() 
fc_ctx, fc_ctx_labels, fc_sctx, fc_sctx_labels = load_fc() 

#sc_plot = plotting.plot_matrix(sc_ctx, figure = (9,9), labels=sc_ctx_labels, vmax = 10, vmin = 0, cmap = 'Blues')
#fc_plot = plotting.plot_matrix(fc_ctx, figure = (9,9), labels=fc_ctx_labels, vmax = 0.8, vmin = 0, cmap = 'Reds')

#%%
# compute weighted degree centrality measures from the connectivity data 
fc_ctx_dc = np.sum(fc_ctx, axis = 0)
sc_ctx_dc = np.sum(sc_ctx, axis = 0)

#map parcelled data to the surface 
fc_ctx_dc_fsa5 = parcel_to_surface(fc_ctx_dc, "aparc_fsa5") #try with schaefer 
sc_ctx_dc_fsa5 = parcel_to_surface(sc_ctx_dc, "aparc_fsa5")

#plot results on brain surface 
#plot_cortical(array_name=fc_ctx_dc_fsa5, surface_name="fsa5", size=(800,400), cmap="Reds", color_bar=True, color_range=(20,30))
#plot_cortical(array_name=sc_ctx_dc_fsa5, surface_name="fsa5", size=(800,400), cmap="Blues", color_bar=True, color_range=(100,300))

#%%
# same as before for subcortical 

# plot subcortical connectivity matrices (data loaded previously)
ssc_plot = plotting.plot_matrix(sc_sctx, figure = (9,9), labels=sc_sctx_labels, vmax = 10, vmin = 0, cmap = 'Blues')
sfc_plot = plotting.plot_matrix(fc_sctx, figure = (9,9), labels=fc_sctx_labels, vmax = 0.8, vmin = 0, cmap = 'Reds')

#compute weighted degree centrality measures from the connectivity data 
fc_sctx_dc = np.sum(fc_sctx, axis =1)
sc_sctx_dc = np.sum(sc_sctx, axis =1)

# project results on the surface brain (same error as below: array shape is not valid!) 
#plot_subcortical(array_name=fc_sctx_dc, ventricles=False, size=(800,400), cmap='Reds', color_bar=True, color_range=(5,10)) 
#plot_subcortical(array_name=sc_sctx_dc, ventricles=False, size=(800,400), cmap='Blues', color_bar=True, color_range=(100,300)) 
#%%

#load summary stats data 
cov, metr1_SubVol, metr2_CortThick, metr3_CortSurf = load_example_data()

# re-order subcortical data alphabetically and by hemisphere 
metr1_SubVol_r = reorder_sctx(metr1_SubVol)

#z-score patients' data relative to controls (lower z-score = more atrophy)
group = cov['Dx'].to_list()
controlCode = 0 
SV_z = zscore_matrix(metr1_SubVol_r.iloc[:, 1:-1], group, controlCode )
CT_z = zscore_matrix(metr2_CortThick.iloc[:, 1:-5], group, controlCode )
SA_z = zscore_matrix(metr3_CortSurf.iloc[:, 1:-5], group, controlCode )

# mean z-score values across individuals with from a specific group 
SV_z_mean = SV_z.iloc[cov[cov["SDx"] == 3].index, :].mean(axis=0) 
CT_z_mean = CT_z.iloc[cov[cov["SDx"] == 3].index, :].mean(axis=0) 
SA_z_mean = SA_z.iloc[cov[cov["SDx"] == 3].index, :].mean(axis=0) 

#remove subcortical values corresponding to the ventricles 
SV_z_mean_noVent = SV_z_mean.drop(['LLatVent', 'RLatVent']).reset_index(drop=True)

#%%

# perform spatial correlation between funtional/structural hubs and z-scores 
fc_ctx_r = np.corrcoef(fc_ctx_dc, CT_z_mean)[0,1]  
fc_sctx_r = np.corrcoef(fc_sctx_dc, SV_z_mean_noVent)[0,1]  

sc_ctx_r = np.corrcoef(sc_ctx_dc, CT_z_mean)[0,1]  
sc_sctx_r = np.corrcoef(sc_sctx_dc, SV_z_mean_noVent)[0,1]  

#store correlation coefficients 
rvals = {'functional cortical hubs': fc_ctx_r, 'functional subcortical hubs': fc_sctx_r,
         'structural cortical hubs': sc_ctx_r, 'structural subcortical hubs': sc_sctx_r}

#%%

# Spin permutation testing for two cortical maps
fc_ctx_p, fc_ctx_d = spin_test(fc_ctx_dc, CT_z_mean, surface_name='fsa5', parcellation_name='aparc',
                               type='pearson', n_rot=1000, null_dist=True)
sc_ctx_p, sc_ctx_d = spin_test(sc_ctx_dc, CT_z_mean, surface_name='fsa5', parcellation_name='aparc',
                               type='pearson', n_rot=1000, null_dist=True)

# Shuf permutation testing for two subcortical maps
fc_sctx_p, fc_sctx_d = shuf_test(fc_sctx_dc, SV_z_mean_noVent, n_rot=1000,
                                 type='pearson', null_dist=True)
sc_sctx_p, sc_sctx_d = shuf_test(sc_sctx_dc, SV_z_mean_noVent, n_rot=1000,
                                 type='pearson', null_dist=True)

# Store p-values and null distributions
p_and_d = {'functional cortical hubs': [fc_ctx_p, fc_ctx_d], 'functional subcortical hubs': [fc_sctx_p, fc_sctx_d],
           'structural cortical hubs': [sc_ctx_p, sc_ctx_d], 'structural subcortical hubs': [sc_sctx_p, sc_sctx_d]}

#%%

# we can also plot the null distributions of generated correlations
fig, axs = plt.subplots(1, 4, figsize=(15, 3))

for k, (fn, dd) in enumerate(p_and_d.items()):
    # Define plot colors
    if k <= 1:
        col = '#A8221C'     # red for functional hubs
    else:
        col = '#324F7D'     # blue for structural hubs

    # Plot null distributions
    axs[k].hist(dd[1], bins=50, density=True, color=col, edgecolor='white', lw=0.5)
    axs[k].axvline(rvals[fn], lw=1.5, ls='--', color='k', dashes=(2, 3),
                   label='$r$={:.2f}'.format(rvals[fn]) + '\n$p$={:.3f}'.format(dd[0]))
    axs[k].set_xlabel('Null correlations \n ({})'.format(fn))
    axs[k].set_ylabel('Density')
    axs[k].spines['top'].set_visible(False)
    axs[k].spines['right'].set_visible(False)
    axs[k].legend(loc=1, frameon=False)

fig.tight_layout()
plt.show()

#%%

# Store degree centrality and atrophy measures
meas = {('functional cortical hubs', 'cortical thickness'): [fc_ctx_dc, CT_z_mean],
        ('functional subcortical hubs', 'subcortical volume'): [fc_sctx_dc, SV_z_mean_noVent],
        ('structural cortical hubs', 'cortical thickness'): [sc_ctx_dc, CT_z_mean],
        ('structural subcortical hubs', 'subcortical volume'): [sc_sctx_dc, SV_z_mean_noVent]}

fig, axs = plt.subplots(1, 4, figsize=(15, 3))

for k, (fn, dd) in enumerate(meas.items()):
    # Define scatter colors
    if k <= 1:
        col = '#A8221C'
    else:
        col = '#324F7D'

    # Plot relationships between hubs and atrophy
    axs[k].scatter(meas[fn][0], meas[fn][1], color=col,
                   label='$r$={:.2f}'.format(rvals[fn[0]]) + '\n$p$={:.3f}'.format(p_and_d[fn[0]][0]))
    m, b = np.polyfit(meas[fn][0], meas[fn][1], 1)
    axs[k].plot(meas[fn][0], m * meas[fn][0] + b, color=col)
    axs[k].set_ylim((-3.5, 1.5))
    axs[k].set_xlabel('{}'.format(fn[0].capitalize()))
    axs[k].set_ylabel('{}'.format(fn[1].capitalize()))
    axs[k].spines['top'].set_visible(False)
    axs[k].spines['right'].set_visible(False)
    axs[k].legend(loc=1, frameon=False, markerscale=0)

fig.tight_layout()
plt.show()

#%%

#==============================================================================
# EPICENTER MAPPING 
#==============================================================================

# Identify cortical epicenters (from functional connectivity)
fc_ctx_epi = []
fc_ctx_epi_p = []
for seed in range(fc_ctx.shape[0]):
    seed_con = fc_ctx[:, seed]
    fc_ctx_epi = np.append(fc_ctx_epi, np.corrcoef(seed_con, CT_z_mean)[0, 1])
    fc_ctx_epi_p = np.append(fc_ctx_epi_p,
                             spin_test(seed_con, CT_z_mean, surface_name='fsa5', parcellation_name='aparc',
                                       type='pearson', n_rot=10, null_dist=False))

# Identify cortical epicenters (from structural connectivity)
sc_ctx_epi = []
sc_ctx_epi_p = []
for seed in range(sc_ctx.shape[0]):
    seed_con = sc_ctx[:, seed]
    sc_ctx_epi = np.append(sc_ctx_epi, np.corrcoef(seed_con, CT_z_mean)[0, 1])
    sc_ctx_epi_p = np.append(sc_ctx_epi_p,
                             spin_test(seed_con, CT_z_mean, surface_name='fsa5', parcellation_name='aparc',
                                       type='pearson', n_rot=10, null_dist=False))
    
#%%
# Identify subcortical epicenters (from functional connectivity)
fc_sctx_epi = []
fc_sctx_epi_p = []
for seed in range(fc_sctx.shape[0]):
    seed_con = fc_sctx[seed, :]
    fc_sctx_epi = np.append(fc_sctx_epi, np.corrcoef(seed_con, CT_z_mean)[0, 1])
    fc_sctx_epi_p = np.append(fc_sctx_epi_p,
                              spin_test(seed_con, CT_z_mean, surface_name='fsa5', n_rot=100))

# Identify subcortical epicenters (from structural connectivity)
sc_sctx_epi = []
sc_sctx_epi_p = []
for seed in range(sc_sctx.shape[0]):
    seed_con = sc_sctx[seed, :]
    sc_sctx_epi = np.append(sc_sctx_epi, np.corrcoef(seed_con, CT_z_mean)[0, 1])
    sc_sctx_epi_p = np.append(sc_sctx_epi_p,
                              spin_test(seed_con, CT_z_mean, surface_name='fsa5', n_rot=100))


#%% 
dir ="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
output_folder= f"{dir}/output/epicentermapping"

np.savetxt(f"{output_folder}/fc_ctx_epi_example.csv", fc_ctx_epi, delimiter = ",")
np.savetxt(f"{output_folder}/fc_ctx_epi_p_example.csv", fc_ctx_epi_p, delimiter = ",")
np.savetxt(f"{output_folder}/sc_ctx_epi_example.csv", sc_ctx_epi, delimiter = ",")
np.savetxt(f"{output_folder}/sc_ctx_epi_p_example.csv", sc_ctx_epi_p, delimiter = ",")

np.savetxt(f"{output_folder}/fc_sctx_epi_example.csv", fc_sctx_epi, delimiter = ",")
np.savetxt(f"{output_folder}/fc_sctx_epi_p_example.csv", fc_sctx_epi_p, delimiter = ",")
np.savetxt(f"{output_folder}/sc_sctx_epi_example.csv", sc_sctx_epi, delimiter = ",")
np.savetxt(f"{output_folder}/sc_sctx_epi_p_example.csv", sc_sctx_epi_p, delimiter = ",")
#%%
# Project the results on the surface brain
# Selecting only regions with p < 0.1 (functional epicenters)
fc_ctx_epi_p_sig = np.zeros_like(fc_ctx_epi_p)
fc_ctx_epi_p_sig[np.argwhere(fc_ctx_epi_p < 0.1)] = fc_ctx_epi[np.argwhere(fc_ctx_epi_p < 0.1)]
plot_cortical(array_name=parcel_to_surface(fc_ctx_epi_p_sig, 'aparc_fsa5'), surface_name="fsa5", size=(800, 400),
              cmap='GyRd_r', color_bar=True, color_range=(-0.5, 0.5))

# Selecting only regions with p < 0.1 (structural epicenters)
sc_ctx_epi_p_sig = np.zeros_like(sc_ctx_epi_p)
sc_ctx_epi_p_sig[np.argwhere(sc_ctx_epi_p < 0.1)] = sc_ctx_epi[np.argwhere(sc_ctx_epi_p < 0.1)]
plot_cortical(array_name=parcel_to_surface(sc_ctx_epi_p_sig, 'aparc_fsa5'), surface_name="fsa5", size=(800, 400),
              cmap='GyBu_r', color_bar=True, color_range=(-0.5, 0.5))

# Project the results on the surface brain
# Selecting only regions with p < 0.1 (functional epicenters)
fc_sctx_epi_p_sig = np.zeros_like(fc_sctx_epi_p)
fc_sctx_epi_p_sig[np.argwhere(fc_sctx_epi_p < 0.1)] = fc_sctx_epi[np.argwhere(fc_sctx_epi_p < 0.1)]
plot_subcortical(fc_sctx_epi_p_sig, ventricles=False, size=(800, 400),
                 cmap='GyRd_r', color_bar=True, color_range=(-0.5, 0.5))

# Selecting only regions with p < 0.1 (functional epicenters)
sc_sctx_epi_p_sig = np.zeros_like(sc_sctx_epi_p)
sc_sctx_epi_p_sig[np.argwhere(sc_sctx_epi_p < 0.1)] = sc_sctx_epi[np.argwhere(sc_sctx_epi_p < 0.1)]
plot_subcortical(sc_sctx_epi_p_sig, ventricles=False, size=(800, 400),
                 cmap='GyBu_r', color_bar=True, color_range=(-0.5, 0.5))