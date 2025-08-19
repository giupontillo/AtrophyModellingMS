#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 22 Jan 2024

import numpy as np
import pandas as pd
#from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.permutation_testing import spin_test, shuf_test

#TODO: add repeated functions separated code
#from helpers_functions import *
#%%
#==============================================================================
# FUNCTIONS 
#==============================================================================

def shuf_test_1(map1, map2, n_rot=1000, type='pearson', null_dist=False):
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
    
    return rvals, p_and_d

#%%

#==============================================================================
# IMPORT DATA 
#==============================================================================
# ################# ATROPHY DATA FROM AMSTERDAM ################# 
dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
atrophydir = f'{dir}/data/atrophy_cross'

atrophy_long_indv_100 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex114Parcels_individuals_wscore.csv", index_col="id")
atrophy_long_indv_400 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv", index_col="id")

#atrophy_long_group_100 = pd.read_csv(f"{atrophydir}/PrograMS_atrophy_long_SchaeferSubcortex114Parcels_group_t-values.csv")
#atrophy_long_group_400 = pd.read_csv(f"{atrophydir}/PrograMS_atrophy_long_SchaeferSubcortex414Parcels_group_t-values.csv")

# ################# NEIGHBOUR DATA FROM AMSTERDAM ################# 
neighdir = f'{dir}/data/neighbourhood_atrophy/'

neighcross_indv_100_fc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_fc_ctx.csv", index_col="id")
neighcross_indv_100_fc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_fc_sctx.csv", index_col="id")
neighcross_indv_100_sc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_sc_ctx.csv", index_col="id")
neighcross_indv_100_sc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex114Parcels_individuals_sc_sctx.csv", index_col="id")

neighcross_indv_400_fc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_fc_ctx.csv", index_col="id")
neighcross_indv_400_fc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_fc_sctx.csv", index_col="id")
neighcross_indv_400_sc_ctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_sc_ctx.csv", index_col="id")
neighcross_indv_400_sc_sctx = pd.read_csv(f"{neighdir}/neighbourhood_atrophy_cross_SchaeferSubcortex414Parcels_individuals_sc_sctx.csv", index_col="id")

#%%
# ################# OUTPUT DIRECTORY ################# 

outdir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/june/cross_indv_results"

#%%
#==============================================================================
# COMPUTE TRANSNEURNAL DEGENERATION
#==============================================================================
neighbour_atrophy_d100={}

for indv, row in atrophy_long_indv_100.iterrows():
    
    print(indv)
    
    #select subject from neighbour data
    neighcross_100_fc_ctx = neighcross_indv_100_fc_ctx.loc[neighcross_indv_100_fc_ctx.index.astype(str) == str(indv)]
    neighcross_100_sc_ctx = neighcross_indv_100_sc_ctx.loc[neighcross_indv_100_sc_ctx.index.astype(str) == str(indv)]
    neighcross_100_fc_sctx = neighcross_indv_100_fc_sctx.loc[neighcross_indv_100_fc_sctx.index.astype(str) == str(indv)]
    neighcross_100_sc_sctx = neighcross_indv_100_sc_sctx.loc[neighcross_indv_100_fc_sctx.index.astype(str) == str(indv)]
    #transform neighbour data to the needed format 
    neighcross_100_fc_ctx=neighcross_100_fc_ctx.T.squeeze().astype(np.float64)
    neighcross_100_sc_ctx=neighcross_100_sc_ctx.T.squeeze().astype(np.float64)
    neighcross_100_fc_sctx=neighcross_100_fc_sctx.T.squeeze().astype(np.float64)
    neighcross_100_sc_sctx=neighcross_100_sc_sctx.T.squeeze().astype(np.float64)
    
    row = row.to_frame().T 
    rvals, p_and_d = neighbour_atrophy(neighcross_100_fc_ctx, neighcross_100_fc_sctx, 
                                       neighcross_100_sc_ctx, neighcross_100_sc_sctx,
                                       row, "schaefer_100")
    neighbour_atrophy_d100[indv] = {'rvals': rvals, 'p_and_d': p_and_d}

    print(f"{indv} done!")
    
#%%
neighbour_atrophy_d100 =  pd.DataFrame(neighbour_atrophy_d100.items()) 
np.save(f'{outdir}/2025-06-18_neighbour-atrophy-indv_Schaefer100_rot-10000.npy', neighbour_atrophy_d100) 
    
#%%
neighbour_atrophy_d400={}

for indv, row in atrophy_long_indv_400.iterrows():
    
    print(indv)
    
    #select subject from neighbour data
    neighcross_400_fc_ctx = neighcross_indv_400_fc_ctx.loc[neighcross_indv_400_fc_ctx.index.astype(str) == str(indv)]
    neighcross_400_sc_ctx = neighcross_indv_400_sc_ctx.loc[neighcross_indv_400_sc_ctx.index.astype(str) == str(indv)]
    neighcross_400_fc_sctx = neighcross_indv_400_fc_sctx.loc[neighcross_indv_400_fc_sctx.index.astype(str) == str(indv)]
    neighcross_400_sc_sctx = neighcross_indv_400_sc_sctx.loc[neighcross_indv_400_fc_sctx.index.astype(str) == str(indv)]
    
    #transform neighbour data to the needed format 
    neighcross_400_fc_ctx=neighcross_400_fc_ctx.T.squeeze().astype(np.float64)
    neighcross_400_sc_ctx=neighcross_400_sc_ctx.T.squeeze().astype(np.float64)
    neighcross_400_fc_sctx=neighcross_400_fc_sctx.T.squeeze().astype(np.float64)
    neighcross_400_sc_sctx=neighcross_400_sc_sctx.T.squeeze().astype(np.float64)
    
    row = row.to_frame().T 
    rvals, p_and_d = neighbour_atrophy(neighcross_400_fc_ctx, neighcross_400_fc_sctx, 
                                       neighcross_400_sc_ctx, neighcross_400_sc_sctx,
                                       row, "schaefer_400")
    neighbour_atrophy_d400[indv] = {'rvals': rvals, 'p_and_d': p_and_d}

    print(f"{indv} done!")
    
#%%
neighbour_atrophy_d400 =  pd.DataFrame(neighbour_atrophy_d400.items()) 
np.save(f'{outdir}/2025-06-18_neighbour-atrophy-indv_Schaefer400_rot-10000.npy', neighbour_atrophy_d400) 