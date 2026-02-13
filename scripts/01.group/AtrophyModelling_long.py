#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 25 May 2025

#%%

import os 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from enigmatoolbox.plotting import plot_subcortical
from enigmatoolbox.permutation_testing import spin_test, shuf_test
#%%
#==============================================================================
# FUNCTIONS 
#==============================================================================

def shuf_test_1(map1, map2, n_rot=10000, type='pearson', null_dist=False):
    """Shuf permuation (author: @saratheriver)

    Parameters
    ----------
    map1 : narray, ndarray, or pandas.Series
        One of two map to be correlated
    map2 : narray, ndarray, or pandas.Series
        The other map to be correlated
    n_rot : int, optional
        Number of shuffles. Default is 1000.
    type : string, optional
        Correlation type {'pearson', 'spearman'}. Default is 'pearson'.
    null_dist : bool, optional
        Output null correlations. Default is False.

    Returns
    -------
    p_shuf : float
        Permutation p-value
    r_dist : 1D ndarray
        Null correlations, shape = (n_rot*2,). Only if ``null_dist is True``.
    """
    r = 0  # count successful (r) iterations unsuccessful (c) iterations
    c = 0  # count unsuccessful (c) iterations
    nroi = map1.shape[0]  # number of regions

    # generate random permutations
    perm_id = np.zeros((nroi, n_rot))
    while (r < n_rot):
        rot_lr_sort = np.random.permutation(nroi)

        # verify that permutation does not map to itself
        if np.all(rot_lr_sort == range(nroi)) is not True:
            perm_id[:, r] = rot_lr_sort
            r = r + 1
        elif np.all(rot_lr_sort == range(nroi)) is True:
            c = c + 1
            print('map to itself n.' + str(c))

        # track progress
        if np.mod(r, 100) == 0:
            print('permutation ' + str(r) + ' of ' + str(n_rot))

    # empirical correlation
    #rho_emp = pd.Series(map1).corr(pd.Series(map2), method=type)
    rho_emp =  np.corrcoef(map1, map2)[0,1]  

    # permutation of measures
    x_perm = np.empty((0, nroi))
    y_perm = np.empty((0, nroi))
    for rr in range(n_rot):
        x_perm2 = []
        y_perm2 = []
        for ii in range(nroi):
            x_perm2 = np.append(x_perm2, map1[int(perm_id[ii, rr])])
            y_perm2 = np.append(y_perm2, map2[int(perm_id[ii, rr])])
        x_perm = np.vstack((x_perm, x_perm2))
        y_perm = np.vstack((y_perm, y_perm2))

    x_perm = np.transpose(x_perm)
    y_perm = np.transpose(y_perm)

    # correlation to unpermuted measures
    rho_null_xy = []
    rho_null_yx = []
    for rr in range(n_rot):
        #rho_null_xy = np.append(rho_null_xy, pd.Series(x_perm[:, rr]).corr(map2, method=type))
        #rho_null_yx = np.append(rho_null_yx, pd.Series(map1).corr(pd.Series(y_perm[:, rr]), method=type))
        rho_null_xy = np.append(rho_null_xy, np.corrcoef(x_perm[:, rr], map2)[0,1])
        rho_null_yx = np.append(rho_null_yx, np.corrcoef(map1, y_perm[:, rr])[0,1])

    # p-value definition depends on the sign of the empirical correlation
    if rho_emp >= 0:
        p_perm_xy = np.sum((rho_null_xy > rho_emp).astype(int)) / n_rot
        p_perm_yx = np.sum((rho_null_yx > rho_emp).astype(int)) / n_rot
    else:
        p_perm_xy = np.sum((rho_null_xy < rho_emp).astype(int)) / n_rot
        p_perm_yx = np.sum((rho_null_yx < rho_emp).astype(int)) / n_rot

    # average p-values
    p_shuf = (p_perm_xy + p_perm_yx) / 2

    # null distribution
    r_dist = np.append(rho_null_xy, rho_null_yx)


    if null_dist is True:
        return p_shuf, r_dist
    elif null_dist is not True:
        return p_shuf

def load_csv_to_numpy(file_mapping, base_path):
    "reads multiple csv files into numpy arrays"
    arrays = {}
    for key, filename in file_mapping.items():
        arrays[key] = pd.read_csv(f"{base_path}/{filename}", header=None).to_numpy()
    return arrays

def compute_dc(fc_ctx, sc_ctx, fc_sctx, sc_sctx): 
    "computes degree centrality for fmri & dmri"
    fc_ctx_dc = np.sum(fc_ctx, axis = 0)
    sc_ctx_dc = np.sum(sc_ctx, axis = 0)
    fc_sctx_dc = np.sum(fc_sctx, axis =1)
    sc_sctx_dc = np.sum(sc_sctx, axis =1)
    
    return fc_ctx_dc, sc_ctx_dc, fc_sctx_dc, sc_sctx_dc

def plot_cortical_atrophy(atrophy, atlas, output_folder, filename, color_range=(0,20), cmap="RdBu_r"):
    #atrophy_ctx = atrophy.iloc[:,:-14] #remove subcortical 
    #atrophy_ctx.to_numpy()   
    atrophy_surf = parcel_to_surface(atrophy, atlas) 
    plot_cortical(array_name=atrophy_surf, surface_name="fsa5", size=(2000,1500), 
                  cmap=cmap, color_bar=True, color_range=color_range, 
                  #label_text={'left': ['left']},
                  screenshot=True, transparent_bg=False, 
                  filename=f"{output_folder}/{filename}.png")
    
def plot_subcortical_atrophy(atrophy, output_folder, filename, color_range=(0,10), cmap="RdBu_r"):
    #atrophy_sctx = atrophy.iloc[:,-14:] #remove cortical 
    #atrophy_sctx = atrophy_sctx.to_numpy().T.squeeze()     
    plot_subcortical(array_name=atrophy, ventricles=False, 
                     size=(2000,1500), cmap=cmap, color_bar=True, color_range=color_range, 
                     screenshot=True, transparent_bg=False,
                     filename=f"{output_folder}/{filename}.png")
    
def nodal_stress(fc_ctx_dc, fc_sctx_dc, sc_ctx_dc, sc_sctx_dc, atrophy, parcellation): 
    " performs spatial correlation between funtional/structural hubs and atrophy z-scores "
    
    ################# spatial correlation ################# 
    atrophy_ctx = atrophy.iloc[:,:-14].reset_index(drop=True).squeeze()
    atrophy_sctx = atrophy.iloc[:,-14:].reset_index(drop=True).squeeze()
    
    fc_ctx_r = np.corrcoef(fc_ctx_dc, atrophy_ctx)[0,1]  
    fc_sctx_r = np.corrcoef(fc_sctx_dc, atrophy_sctx)[0,1]  
    sc_ctx_r = np.corrcoef(sc_ctx_dc, atrophy_ctx)[0,1]  
    sc_sctx_r = np.corrcoef(sc_sctx_dc, atrophy_sctx)[0,1]  

    #store correlation coefficients 
    rvals = {'functional cortical hubs': fc_ctx_r, 
             'functional subcortical hubs': fc_sctx_r,
             'structural cortical hubs': sc_ctx_r, 
             'structural subcortical hubs': sc_sctx_r}
    
    ################# premutation testing ################# 
    fc_ctx_p, fc_ctx_d = spin_test(fc_ctx_dc, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=10000, null_dist=True)
    sc_ctx_p, sc_ctx_d = spin_test(sc_ctx_dc, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=10000, null_dist=True)
       
    # Shuf permutation testing for two subcortical maps
    fc_sctx_p, fc_sctx_d = shuf_test_1(fc_sctx_dc, atrophy_sctx, n_rot=10000,
                                     type='pearson', null_dist=True)
    sc_sctx_p, sc_sctx_d = shuf_test_1(sc_sctx_dc, atrophy_sctx, n_rot=10000,
                                     type='pearson', null_dist=True)

    # Store p-values and null distributions
    p_and_d = {'functional cortical hubs': [fc_ctx_p, fc_ctx_d], 
               'functional subcortical hubs': [fc_sctx_p, fc_sctx_d],
               'structural cortical hubs': [sc_ctx_p, sc_ctx_d], 
               'structural subcortical hubs': [sc_sctx_p, sc_sctx_d]}
    
    # Store meas
    meas = {('functional cortical hubs', 'cortical atrophy (t-values)'): [fc_ctx_dc, atrophy_ctx],
            ('functional subcortical hubs', 'subcortical atrophy (t-values)'): [fc_sctx_dc, atrophy_sctx],
            ('structural cortical hubs', 'cortical atrophy (t-values)'): [sc_ctx_dc, atrophy_ctx],
            ('structural subcortical hubs', 'subcortical atrophy (t-values)'): [sc_sctx_dc, atrophy_sctx]}
    
    return rvals, p_and_d, meas

def neighbour_atrophy(fc_ctx_ngh, fc_sctx_ngh, sc_ctx_ngh, sc_sctx_ngh, atrophy, parcellation): 
    """
    performs spatial correlation between atrophy z-scores and neighbour atrophy 
    weighted by funtional/structural connection strenght  
    """
    
    ################# spatial correlation ################# 
    atrophy_ctx = atrophy.iloc[:,:-14].reset_index(drop=True).squeeze()
    atrophy_sctx = atrophy.iloc[:,-14:].reset_index(drop=True).squeeze()
    
    fc_ctx_r = np.corrcoef(fc_ctx_ngh, atrophy_ctx)[0,1]  
    fc_sctx_r = np.corrcoef(fc_sctx_ngh, atrophy_sctx)[0,1]  
    sc_ctx_r = np.corrcoef(sc_ctx_ngh, atrophy_ctx)[0,1]  
    sc_sctx_r = np.corrcoef(sc_sctx_ngh, atrophy_sctx)[0,1]  

    #store correlation coefficients 
    rvals = {'functional cortical neighbour atrophy': fc_ctx_r, 
             'functional subcortical neighbour atrophy': fc_sctx_r,
             'structural cortical neighbour atrophy': sc_ctx_r, 
             'structural subcortical neighbour atrophy': sc_sctx_r}
    
    ################# premutation testing ################# 
    fc_ctx_p, fc_ctx_d = spin_test(fc_ctx_ngh, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=10000, null_dist=True)
    sc_ctx_p, sc_ctx_d = spin_test(sc_ctx_ngh, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=10000, null_dist=True)
       
    # Shuf permutation testing for two subcortical maps
    fc_sctx_p, fc_sctx_d = shuf_test_1(fc_sctx_ngh, atrophy_sctx, n_rot=10000,
                                     type='pearson', null_dist=True)
    sc_sctx_p, sc_sctx_d = shuf_test_1(sc_sctx_ngh, atrophy_sctx, n_rot=10000,
                                     type='pearson', null_dist=True)

    # Store p-values and null distributions
    p_and_d = {'functional cortical neighbour atrophy': [fc_ctx_p, fc_ctx_d], 
               'functional subcortical neighbour atrophy': [fc_sctx_p, fc_sctx_d],
               'structural cortical neighbour atrophy': [sc_ctx_p, sc_ctx_d], 
               'structural subcortical neighbour atrophy': [sc_sctx_p, sc_sctx_d]}
    
    meas = {('functional cortical neighbour atrophy', 'cortical atrophy (t-values)'): [fc_ctx_ngh, atrophy_ctx],
            ('functional subcortical neighbour atrophy', 'subcortical atrophy (t-values)'): [fc_sctx_ngh, atrophy_sctx],
            ('structural cortical neighbour atrophy', 'cortical atrophy (t-values)'): [sc_ctx_ngh, atrophy_ctx],
            ('structural subcortical neighbour atrophy', 'subcortical atrophy (t-values)'): [sc_sctx_ngh, atrophy_sctx]}
    
    return rvals, p_and_d, meas

def axonal_disc(sc_ctx_dc, sc_sctx_dc, atrophy, parcellation): 
    " performs spatial correlation between funtional/structural hubs and atrophy z-scores "
    
    ################# spatial correlation ################# 
    atrophy_ctx = atrophy.iloc[:,:-14].reset_index(drop=True).squeeze()
    atrophy_sctx = atrophy.iloc[:,-14:].reset_index(drop=True).squeeze()
    
    sc_ctx_r = np.corrcoef(sc_ctx_dc, atrophy_ctx)[0,1]  
    sc_sctx_r = np.corrcoef(sc_sctx_dc, atrophy_sctx)[0,1]  

    #store correlation coefficients 
    rvals = {'structural cortical disconnectome': sc_ctx_r, 
             'structural subcortical disconnectome': sc_sctx_r}
    
    ################# premutation testing ################# 
    sc_ctx_p, sc_ctx_d = spin_test(sc_ctx_dc, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=10000, null_dist=True)
       
    # Shuf permutation testing for two subcortical maps
    sc_sctx_p, sc_sctx_d = shuf_test_1(sc_sctx_dc, atrophy_sctx, n_rot=10000,
                                     type='pearson', null_dist=True)

    # Store p-values and null distributions
    p_and_d = {'structural cortical disconnectome': [sc_ctx_p, sc_ctx_d], 
               'structural subcortical disconnectome': [sc_sctx_p, sc_sctx_d]}
    
    meas = {('structural cortical disconnectome', 'cortical atrophy (t-values)'): [sc_ctx_dc, atrophy_ctx],
            ('structural subcortical disconnectome', 'subcortical atrophy (t-values)'): [sc_sctx_dc, atrophy_sctx]}
    
    
    return rvals, p_and_d, meas

def transcriptomic_vulnerablity(cge_ctx, cge_sctx, atrophy, parcellation): 
    """
    performs spatial correlation between atrophy z-scores and neighbour atrophy 
    weighted by funtional/structural connection strenght  
    """
    
    ################# spatial correlation ################# 
    atrophy_ctx = atrophy.iloc[:,:-14].reset_index(drop=True).squeeze()
    atrophy_sctx = atrophy.iloc[:,-14:].reset_index(drop=True).squeeze()
    
    cge_ctx_r = np.corrcoef(cge_ctx, atrophy_ctx)[0,1]  
    cge_sctx_r = np.corrcoef(cge_sctx, atrophy_sctx)[0,1]   

    #store correlation coefficients 
    rvals = {'cortical gene expression': cge_ctx_r, 
             'subcortical gene expression': cge_sctx_r}
    
    ################# premutation testing ################# 
    fc_ctx_p, fc_ctx_d = spin_test(cge_ctx, atrophy_ctx, 
                                   surface_name='fsa5', parcellation_name=parcellation,
                                   type='pearson', n_rot=100, null_dist=True)
       
    # Shuf permutation testing for two subcortical maps
    fc_sctx_p, fc_sctx_d = shuf_test_1(cge_sctx, atrophy_sctx, n_rot=100,
                                     type='pearson', null_dist=True)

    # Store p-values and null distributions
    p_and_d = {'cortical gene expression': [fc_ctx_p, fc_ctx_d], 
               'subcortical gene expression': [fc_sctx_p, fc_sctx_d]}
    
    meas = {('cortical gene expression', 'cortical atrophy (t-values)'): [cge_ctx, atrophy_ctx],
            ('subcortical gene expression', 'subcortical atrophy (t-values)'): [cge_sctx, atrophy_sctx]}
        
    return rvals, p_and_d, meas

def plot_null_distributions(p_and_d, rvals, filename): 
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
    fig.savefig(f"{output_folder}/{filename}.png", dpi=300, bbox_inches='tight')
   # plt.show()

def plot_hubs_vs_atrophy(p_and_d, rvals, meas, filename):
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
        axs[k].set_ylim((-5, 17.5))
        axs[k].set_xlabel('{}'.format(fn[0].capitalize()))
        axs[k].set_ylabel('{}'.format(fn[1].capitalize()))
        axs[k].spines['top'].set_visible(False)
        axs[k].spines['right'].set_visible(False)
        axs[k].legend(loc=1, frameon=False, markerscale=0)
    
    fig.tight_layout()
    fig.savefig(f"{output_folder}/{filename}.png", dpi=300, bbox_inches='tight')
    #plt.show()
    
#%%
#==============================================================================
# IMPORT DATA 
#==============================================================================

# ################# NORMATIVE DATA HUMAN CONNECTOME PROJECT ################# 
sc_ctx_100, sc_ctx_labels_100, sc_sctx_100, sc_sctx_labels_100 = load_sc(parcellation="schaefer_100") 
fc_ctx_100, fc_ctx_labels_100, fc_sctx_100, fc_sctx_labels_100 = load_fc(parcellation="schaefer_100") 
sc_ctx_400, sc_ctx_labels_400, sc_sctx_400, sc_sctx_labels_400 = load_sc(parcellation="schaefer_400") 
fc_ctx_400, fc_ctx_labels_400, fc_sctx_400, fc_sctx_labels_400 = load_fc(parcellation="schaefer_400") 

# structural subcortical are [14,X14], remove 14 subcortical regions 
sc_sctx_100 = sc_sctx_100[:,:-14]
sc_sctx_400 = sc_sctx_400[:,:-14]

#normative enigma data
fc_ctx_dc_100, sc_ctx_dc_100, fc_sctx_dc_100, sc_sctx_dc_100 = compute_dc(fc_ctx_100, sc_ctx_100, fc_sctx_100, sc_sctx_100)
fc_ctx_dc_400, sc_ctx_dc_400, fc_sctx_dc_400, sc_sctx_dc_400 = compute_dc(fc_ctx_400, sc_ctx_400, fc_sctx_400, sc_sctx_400)

#%%
# ################# ATROPHY LONG DATA FROM AMSTERDAM ################# 
dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
atrophydir = f'{dir}/data/atrophy_long'

output_folder = f"{dir}/output/oct/long_group_results"
os.makedirs(output_folder, exist_ok=True)

atrophy_long_group_100 = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex114Parcels_group_t-values.csv")
atrophy_long_group_400 = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex414Parcels_group_t-values.csv")

atrophy_long_group_100_fdr = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex114Parcels_group_pfdr.csv")
atrophy_long_group_400_fdr = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex414Parcels_group_pfdr.csv")

# ################# ATROPHY LONG WSCORES DATA FROM AMSTERDAM ################# 

atrophy_wscores_group_100 = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex114Parcels_group_t-values_wscore.csv")
atrophy_wscores_group_400 = pd.read_csv(f"{atrophydir}/atrophy_long_SchaeferSubcortex414Parcels_group_t-values_wscore.csv")

#%%
# ################# NEIGHBOUR DATA FROM AMSTERDAM ################# 
neighdir = f'{dir}/data/neighbourhood_atrophy/'

neighcross_group_100_fc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_fc_ctx.csv")
neighcross_group_100_fc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_fc_sctx.csv")
neighcross_group_100_sc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_sc_ctx.csv")
neighcross_group_100_sc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_sc_sctx.csv")

neighcross_group_400_fc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_fc_ctx.csv")
neighcross_group_400_fc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_fc_sctx.csv")
neighcross_group_400_sc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_sc_ctx.csv")
neighcross_group_400_sc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_sc_sctx.csv")

#%%
# ################# DISCONNECTOME DATA FROM AMSTERDAM ################# 
discdir = f'{dir}/data/disc/'

#disc_group_100_ctx = pd.read_csv(f"{discdir}/disc_SchaeferSubcortex114Parcels_group_long_ctx.csv", header=None)
#disc_group_100_sctx = pd.read_csv(f"{discdir}/disc_SchaeferSubcortex114Parcels_group_long_sctx.csv", header=None)

#disc_group_400_ctx = pd.read_csv(f"{discdir}/disc_SchaeferSubcortex414Parcels_group_long_ctx.csv", header=None)
#disc_group_400_sctx = pd.read_csv(f"{discdir}/disc_SchaeferSubcortex414Parcels_group_long_sctx.csv", header=None)

#%%
# ################# ALTERNATIVE DISCONNECTOME DATA #################
parcel_group_100_ctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex114Parcels_group_long_ctx.csv", header=None)
parcel_group_100_sctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex114Parcels_group_long_sctx.csv", header=None)

parcel_group_400_ctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex414Parcels_group_long_ctx.csv", header=None)
parcel_group_400_sctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex414Parcels_group_long_sctx.csv", header=None)

#%%
# ################# NEIGHBOUR DATA WEIGHTED BY GENE #################

cge_group_100_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_cge_ctx.csv")
cge_group_100_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_group_long_cge_sctx.csv")
cge_group_400_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_cge_ctx.csv")
cge_group_400_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_group_long_cge_sctx.csv")

#%%
# =============================================================================
# #Load functional data from programs
# dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
# fc = f'{dir}/data/fc'
# 
# connectivities = ["fc"] 
# regions = ["ctx", "sctx"] 
# atlas = ["114", "414"] 
# groups = ["HC", "MS"] 
# 
# file_mapping = {
#     f"{conn}_{region}_{res}_{group}": f"{conn}_SchaeferSubcortex{res}Parcels_{group}_{region}.csv"   
#     for conn in connectivities
#     for region in regions 
#     for res in atlas 
#     for group in groups
# } 
# 
# fcp = load_csv_to_numpy(file_mapping, fc)
#
# #programs data
# fc_ctx_dc_MS100,  sc_ctx_dc_100, fc_sctx_dc_MS100, sc_sctx_dc_100 = compute_dc(fcp["fc_ctx_114_MS"], sc_ctx_100, fcp["fc_sctx_114_MS"], sc_sctx_100)
# fc_ctx_dc_MS400,  sc_ctx_dc_400, fc_sctx_dc_MS400, sc_sctx_dc_400 = compute_dc(fcp["fc_ctx_414_MS"], sc_ctx_400, fcp["fc_sctx_414_MS"], sc_sctx_400)
# 
# fc_ctx_dc_HC100,  sc_ctx_dc_100, fc_sctx_dc_HC100, sc_sctx_dc_100 = compute_dc(fcp["fc_ctx_114_HC"], sc_ctx_100, fcp["fc_sctx_114_HC"], sc_sctx_100)
# fc_ctx_dc_HC400,  sc_ctx_dc_400, fc_sctx_dc_HC400, sc_sctx_dc_400 = compute_dc(fcp["fc_ctx_414_HC"], sc_ctx_400, fcp["fc_sctx_414_HC"], sc_sctx_400)
# 
# =============================================================================
#%%
#==============================================================================
# PLOT DATA 
#==============================================================================

# ################# ATROPHY LONG DATA FROM AMSTERDAM ################# 
atrophy_d100_ctx = atrophy_long_group_100.iloc[:,:-14].to_numpy() #remove sctx
atrophy_d400_ctx = atrophy_long_group_400.iloc[:,:-14].to_numpy() #remove sctx

plot_cortical_atrophy(atrophy_d100_ctx, "schaefer_100_fsa5", output_folder, "atrophy_long_group_100_ctx_map", color_range=(-15,15), cmap="RdBu_r")
plot_cortical_atrophy(atrophy_d400_ctx, "schaefer_400_fsa5", output_folder, "atrophy_long_group_400_ctx_map", color_range=(-15,15), cmap="RdBu_r")

atrophy_d100_sctx = atrophy_long_group_100.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 
atrophy_d400_sctx = atrophy_long_group_100.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 
 
plot_subcortical_atrophy(atrophy_d100_sctx, output_folder, "atrophy_long_group_100_sctx_map", color_range=(-15, 15), cmap="RdBu_r")
plot_subcortical_atrophy(atrophy_d400_sctx, output_folder, "atrophy_long_group_400_sctx_map", color_range=(-15, 15), cmap="RdBu_r")

#%%

atrophy_pfdr100_ctx = atrophy_long_group_100_fdr.iloc[:,:-14].to_numpy() #remove sctx
atrophy_pfdr400_ctx = atrophy_long_group_400_fdr.iloc[:,:-14].to_numpy() #remove sctx

atrophy_pfdr100_sctx = atrophy_long_group_100_fdr.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 
atrophy_pfdr400_sctx = atrophy_long_group_400_fdr.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 

log_atrophy_pfdr100_ctx = -np.log(atrophy_pfdr100_ctx)
log_atrophy_pfdr400_ctx = -np.log(atrophy_pfdr400_ctx)

log_atrophy_pfdr100_sctx = -np.log(atrophy_pfdr100_sctx)
log_atrophy_pfdr400_sctx = -np.log(atrophy_pfdr400_sctx)

plot_cortical_atrophy(log_atrophy_pfdr100_ctx, "schaefer_100_fsa5", output_folder, "atrophy_long_group_100_ctx_pfdr_map", color_range=(1.3,20), cmap="Reds")
plot_cortical_atrophy(log_atrophy_pfdr400_ctx, "schaefer_400_fsa5", output_folder, "atrophy_long_group_400_ctx_pfdr_map", color_range=(1.3,20), cmap="Reds")
 
plot_subcortical_atrophy(log_atrophy_pfdr100_sctx, output_folder, "atrophy_long_group_100_sctx_pfdr_map", color_range=(1.3,20), cmap="Reds")
plot_subcortical_atrophy(log_atrophy_pfdr400_sctx, output_folder, "atrophy_long_group_400_sctx_pfdr_map", color_range=(1.3,20), cmap="Reds")

#%%

# ################# ATROPHY LONG DATA FROM AMSTERDAM ################# 
# atrophy_d100_ctx = atrophy_wscores_group_100.iloc[:,:-14].to_numpy() #remove sctx
# atrophy_d400_ctx = atrophy_wscores_group_100.iloc[:,:-14].to_numpy() #remove sctx

# plot_cortical_atrophy(atrophy_d100_ctx, "schaefer_100_fsa5", output_folder, "atrophy_long_wscore_group_100_ctx_map", color_range=(-10,10), cmap="RdBu_r")
# plot_cortical_atrophy(atrophy_d400_ctx, "schaefer_400_fsa5", output_folder, "atrophy_long_wscore_group_400_ctx_map", color_range=(-10,10), cmap="RdBu_r")

# atrophy_d100_sctx = atrophy_wscores_group_100.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 
# atrophy_d400_sctx = atrophy_wscores_group_100.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 
 
# plot_subcortical_atrophy(atrophy_d100_sctx, output_folder, "atrophy_long_wscore_group_100_sctx_map", color_range=(-15, 15), cmap="RdBu_r")
# plot_subcortical_atrophy(atrophy_d400_sctx, output_folder, "atrophy_long_wscore_group_400_sctx_map", color_range=(-15, 15), cmap="RdBu_r")

#%%
# ################# NORMATIVE DATA HUMAN CONNECTOME PROJECT ################# 

plot_cortical_atrophy(fc_ctx_dc_100, "schaefer_100_fsa5", output_folder, "fc_ctx_dc_100_cortical_map", color_range=(-40,40), cmap="RdBu_r")
plot_cortical_atrophy(sc_ctx_dc_100, "schaefer_100_fsa5", output_folder, "sc_ctx_dc_100_cortical_map", color_range=(-340,340), cmap="RdBu")

plot_subcortical_atrophy(fc_sctx_dc_100, output_folder, "fc_sctx_dc_100_subcortical_map", color_range=(-20,20), cmap="RdBu_r")
plot_subcortical_atrophy(sc_sctx_dc_100, output_folder, "sc_sctx_dc_100_subcortical_map", color_range=(-410,410), cmap="RdBu")

plot_cortical_atrophy(fc_ctx_dc_400, "schaefer_400_fsa5", output_folder, "fc_ctx_dc_400_cortical_map", color_range=(-40,40), cmap="RdBu_r")
plot_cortical_atrophy(sc_ctx_dc_400, "schaefer_400_fsa5", output_folder, "sc_ctx_dc_400_cortical_map", color_range=(-340,340), cmap="RdBu")

plot_subcortical_atrophy(fc_sctx_dc_400, output_folder, "fc_sctx_dc_400_subcortical_map", color_range=(-20,20), cmap="RdBu_r")
plot_subcortical_atrophy(sc_sctx_dc_400, output_folder, "sc_sctx_dc_400_subcortical_map", color_range=(-410,410), cmap="RdBu")

#%%
# ################# NEIGHBOURHOOD DATA ################# 
neighcross_group_100_fc_ctx=neighcross_group_100_fc_ctx.T.squeeze().astype(np.float64)
neighcross_group_100_sc_ctx=neighcross_group_100_sc_ctx.T.squeeze().astype(np.float64)
neighcross_group_100_fc_sctx=neighcross_group_100_fc_sctx.T.squeeze().astype(np.float64)
neighcross_group_100_sc_sctx=neighcross_group_100_sc_sctx.T.squeeze().astype(np.float64)

neighcross_group_400_fc_ctx=neighcross_group_400_fc_ctx.T.squeeze().astype(np.float64)
neighcross_group_400_sc_ctx=neighcross_group_400_sc_ctx.T.squeeze().astype(np.float64)
neighcross_group_400_fc_sctx=neighcross_group_400_fc_sctx.T.squeeze().astype(np.float64)
neighcross_group_400_sc_sctx=neighcross_group_400_sc_sctx.T.squeeze().astype(np.float64)

#%%
plot_cortical_atrophy(neighcross_group_100_fc_ctx, "schaefer_100_fsa5", output_folder, "neighcross_group_100_fc_ctx_map", color_range=(-0.2,0.2), cmap="RdBu_r")
plot_cortical_atrophy(neighcross_group_100_sc_ctx, "schaefer_100_fsa5", output_folder, "neighcross_group_100_sc_ctx_map", color_range=(-4,4), cmap="RdBu")

plot_cortical_atrophy(neighcross_group_400_fc_ctx, "schaefer_400_fsa5", output_folder, "neighcross_group_400_fc_ctx_map", color_range=(-0.2,0.2), cmap="RdBu_r")
plot_cortical_atrophy(neighcross_group_400_sc_ctx, "schaefer_400_fsa5", output_folder, "neighcross_group_400_sc_ctx_map", color_range=(-4,4), cmap="RdBu")

plot_subcortical_atrophy(neighcross_group_100_fc_sctx, output_folder, "neighcross_group_100_fc_sctx_map", color_range=(-0.1,0.1), cmap="RdBu_r")
plot_subcortical_atrophy(neighcross_group_100_sc_sctx, output_folder, "neighcross_group_100_sc_sctx_map", color_range=(-3,3), cmap="RdBu")

plot_subcortical_atrophy(neighcross_group_400_fc_sctx, output_folder, "neighcross_group_400_fc_sctx_map", color_range=(-0.1,0.1), cmap="RdBu_r")
plot_subcortical_atrophy(neighcross_group_400_sc_sctx, output_folder, "neighcross_group_400_sc_sctx_map", color_range=(-3,3), cmap="RdBu")

#%%
# ################# DISCONNECTOME DATA #################

#TODO: warning change 400 in the future
# disc_ctx_dc_100, disc_ctx_dc_400, disc_sctx_dc_100, disc_sctx_dc_400 = compute_dc(disc_group_100_ctx, disc_group_100_ctx, disc_group_100_sctx, disc_group_100_sctx)

# plot_cortical_atrophy(disc_ctx_dc_100, "schaefer_100_fsa5", output_folder, "disc_group_100_ctx_map", color_range=(0,2600), cmap="Blues")
# #plot_cortical_atrophy(disc_ctx_dc_400, "schaefer_400_fsa5", output_folder, "disc_group_400_ctx_map", color_range=(0,2600), cmap="Blues")

# plot_subcortical_atrophy(disc_sctx_dc_100, output_folder, "disc_group_100_sctx_map", color_range=(0,1200), cmap="Blues")
# #plot_subcortical_atrophy(disc_sctx_dc_400, output_folder, "disc_group_400_sctx_map", color_range=(0,1200), cmap="Blues")

#%%
# ################# ALTERNATIVE DISCONNECTOME DATA #################
plot_cortical_atrophy(parcel_group_100_ctx.squeeze(), "schaefer_100_fsa5", output_folder, "parcel_group_100_ctx_map", color_range=(0,0.3), cmap="Blues")
plot_cortical_atrophy(parcel_group_400_ctx.squeeze(), "schaefer_400_fsa5", output_folder, "parcel_group_400_ctx_map", color_range=(0,0.3), cmap="Blues")

plot_subcortical_atrophy(parcel_group_100_sctx.squeeze(), output_folder, "parcel_group_100_sctx_map", color_range=(0,0.3), cmap="Blues")
plot_subcortical_atrophy(parcel_group_400_sctx.squeeze(), output_folder, "parcel_group_400_sctx_map", color_range=(0,0.3), cmap="Blues")

#%%

# ################# TRANSCRIPTOMIC DATA #################
cge_group_100_ctx.fillna(0, inplace=True) 
cge_group_100_sctx.fillna(0, inplace=True) 
cge_group_400_ctx.fillna(0, inplace=True) 
cge_group_400_sctx.fillna(0, inplace=True) 

cge_group_100_ctx=cge_group_100_ctx.T.squeeze().astype(np.float64)
cge_group_400_ctx=cge_group_400_ctx.T.squeeze().astype(np.float64)

cge_group_100_sctx=cge_group_100_sctx.T.squeeze().astype(np.float64)
cge_group_400_sctx=cge_group_400_sctx.T.squeeze().astype(np.float64)

plot_cortical_atrophy(cge_group_100_ctx, "schaefer_100_fsa5", output_folder, "cge_group_100_ctx_map", color_range=(0,0.2), cmap="Blues")
plot_cortical_atrophy(cge_group_400_ctx, "schaefer_400_fsa5", output_folder, "cge_group_400_ctx_map", color_range=(0,0.2), cmap="Blues")

plot_subcortical_atrophy(cge_group_100_sctx, output_folder, "cge_group_100_sctx_map", color_range=(-0.2,0), cmap="Blues")
plot_subcortical_atrophy(cge_group_400_sctx, output_folder, "cge_group_400_sctx_map", color_range=(-0.2,0), cmap="Blues")

#%%
#==============================================================================
# NODAL STRESS 
#==============================================================================

rvals_longgroup_100, p_and_d_longgroup_100, meas_longgroup_100 = nodal_stress(fc_ctx_dc_100, fc_sctx_dc_100, sc_ctx_dc_100, sc_sctx_dc_100, atrophy_long_group_100, "schaefer_100")
rvals_longgroup_400, p_and_d_longgroup_400, meas_longgroup_400 = nodal_stress(fc_ctx_dc_400, fc_sctx_dc_400, sc_ctx_dc_400, sc_sctx_dc_400, atrophy_long_group_400, "schaefer_400")
    
plot_null_distributions(p_and_d_longgroup_100, rvals_longgroup_100, "plot_nodalstress_group_long_100_null")
plot_null_distributions(p_and_d_longgroup_400, rvals_longgroup_400, "plot_nodalstress_group_long_400_null")

plot_hubs_vs_atrophy(p_and_d_longgroup_100, rvals_longgroup_100, meas_longgroup_100, "plot_nodalstress_group_long_100")
plot_hubs_vs_atrophy(p_and_d_longgroup_400, rvals_longgroup_400, meas_longgroup_400, "plot_nodalstress_group_long_400")

#%%
#==============================================================================
# TRANSNEURONAL DEGENERATION (neighbour atrophy)
#==============================================================================
#neighcross_group_100_fc_ctx_dc,  neighcross_group_100_sc_ctx_dc, neighcross_group_100_fc_sctx_dc, neighcross_group_100_sc_sctx_dc = compute_dc(neighcross_group_100_fc_ctx,  neighcross_group_100_sc_ctx, neighcross_group_100_fc_sctx, neighcross_group_100_sc_sctx)

rvals_neighcross_100, p_and_d_neighcross_100, meas_neighcross_100 = neighbour_atrophy(neighcross_group_100_fc_ctx, neighcross_group_100_fc_sctx, neighcross_group_100_sc_ctx, neighcross_group_100_sc_sctx, atrophy_long_group_100, "schaefer_100")
rvals_neighcross_400, p_and_d_neighcross_400, meas_neighcross_400 = neighbour_atrophy(neighcross_group_400_fc_ctx, neighcross_group_400_fc_sctx, neighcross_group_400_sc_ctx, neighcross_group_400_sc_sctx, atrophy_long_group_400, "schaefer_400")

plot_null_distributions(p_and_d_neighcross_100, rvals_neighcross_100, "plot_neigh_group_cross_100_null")
plot_null_distributions(p_and_d_neighcross_400, rvals_neighcross_400, "plot_neigh_group_cross_400_null")

plot_hubs_vs_atrophy(p_and_d_neighcross_100, rvals_neighcross_100, meas_neighcross_100, "plot_neigh_group_cross_100")
plot_hubs_vs_atrophy(p_and_d_neighcross_400, rvals_neighcross_400, meas_neighcross_400, "plot_neigh_group_cross_400")

#%%
#==============================================================================
# AXONAL TRANSECTION
#==============================================================================

# rvals_disc_100, p_and_d_disc_100, meas_disc_100  = axonal_disc(disc_ctx_dc_100, disc_sctx_dc_100, atrophy_long_group_100, "schaefer_100")
# #rvals_disc_400, p_and_d_disc_400, meas_disc_400  = axonal_disc(disc_ctx_dc_400, disc_sctx_dc_100, atrophy_long_group_400, "schaefer_400")

# plot_null_distributions(p_and_d_disc_100, rvals_disc_100)
# #plot_null_distributions(p_and_d_disc_400, rvals_disc_400)

# plot_hubs_vs_atrophy(p_and_d_disc_100, rvals_disc_100, meas_disc_100)
# #plot_hubs_vs_atrophy(p_and_d_disc_400, rvals_disc_400, meas_disc_400)

#%%

rvals_parcel_100, p_and_d_parcel_100, meas_parcel_100  = axonal_disc(parcel_group_100_ctx.squeeze(), parcel_group_100_sctx.squeeze(), atrophy_long_group_100, "schaefer_100")
rvals_parcel_400, p_and_d_parcel_400, meas_parcel_400  = axonal_disc(parcel_group_400_ctx.squeeze(), parcel_group_400_sctx.squeeze(), atrophy_long_group_400, "schaefer_400")

plot_null_distributions(p_and_d_parcel_100, rvals_parcel_100, "plot_parcel_group_cross_100_null")
plot_null_distributions(p_and_d_parcel_400, rvals_parcel_400, "plot_parcel_group_cross_400_null")

plot_hubs_vs_atrophy(p_and_d_parcel_100, rvals_parcel_100, meas_parcel_100, "plot_parcel_group_cross_100")
plot_hubs_vs_atrophy(p_and_d_parcel_400, rvals_parcel_400, meas_parcel_400, "plot_parcel_group_cross_400")


#%%
#==============================================================================
# TRANSCRIPTOMIC VULNERABILITY
#==============================================================================

rvals_cge_100, p_and_d_cge_100, meas_cge_100 = transcriptomic_vulnerablity(cge_group_100_ctx, cge_group_100_sctx, atrophy_long_group_100, "schaefer_100")
rvals_cge_400, p_and_d_cge_400, meas_cge_400 = transcriptomic_vulnerablity(cge_group_400_ctx, cge_group_400_sctx, atrophy_long_group_400, "schaefer_400")

plot_null_distributions(p_and_d_cge_100, rvals_cge_100, "plot_cge_group_cross_100_null")
plot_null_distributions(p_and_d_cge_400, rvals_cge_400, "plot_cge_group_cross_400_null")

plot_hubs_vs_atrophy(p_and_d_cge_100, rvals_cge_100, meas_cge_100, "plot_cge_group_cross_100")
plot_hubs_vs_atrophy(p_and_d_cge_400, rvals_cge_400, meas_cge_400, "plot_cge_group_cross_400")
