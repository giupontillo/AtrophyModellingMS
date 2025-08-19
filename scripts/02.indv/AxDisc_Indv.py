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

def split_and_clean_header(df):
    df[['id', 'session']] = df['ID'].str.split('_', expand=True)
    df = df.drop(columns=['ID', 'session'])
    df.set_index('id', inplace=True)
    return df

def compute_dc(disc_ctx, disc_sctx): 
    "computes degree centrality for fmri & dmri"
    disc_ctx_dc = np.sum(disc_ctx, axis = 1)
    disc_sctx_dc = np.sum(disc_sctx, axis=1)
    
    return disc_ctx_dc, disc_sctx_dc

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
    
    return rvals, p_and_d

#%%
# #==============================================================================
# # IMPORT DATA 
# #==============================================================================
# # ################# ATROPHY DATA FROM AMSTERDAM ################# 
# dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
# atrophydir = f'{dir}/data/atrophy_long'
# discdir = f'{dir}/data/disc'

# atrophy_long_indv_100 = pd.read_csv(f"{atrophydir}/PrograMS_atrophy_long_SchaeferSubcortex114Parcels_individual_slopes.csv", index_col="id")
# atrophy_long_indv_400 = pd.read_csv(f"{atrophydir}/PrograMS_atrophy_long_SchaeferSubcortex414Parcels_individual_slopes.csv", index_col="id")

# # ################# DISCONNECTOME DATA FROM AMSTERDAM ################# 

# disc_indv_100_ctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex114Parcels_individuals_ctx.csv")
# disc_indv_400_ctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex414Parcels_individuals_ctx.csv")

# disc_indv_100_sctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex114Parcels_individuals_sctx.csv")
# disc_indv_400_sctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex414Parcels_individuals_sctx.csv")

# # ################# clean disc data ################# 

# disc_indv_100_ctx = split_and_clean_header(disc_indv_100_ctx)
# disc_indv_100_sctx = split_and_clean_header(disc_indv_100_sctx)

# disc_indv_400_ctx = split_and_clean_header(disc_indv_400_ctx)
# disc_indv_400_sctx = split_and_clean_header(disc_indv_400_sctx)

# #%%
# # ################# OUTPUT DIRECTORY ################# 
# outdir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/atrophy_cross_results_v042025"

#%%

#==============================================================================
# IMPORT DATA 
#==============================================================================
# ################# ATROPHY DATA FROM AMSTERDAM ################# 
dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
atrophydir = f'{dir}/data/atrophy_cross'
discdir = f'{dir}/data/disc'

atrophy_long_indv_100 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex114Parcels_individuals_wscore.csv", index_col="id")
atrophy_long_indv_400 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv", index_col="id")

# ################# DISCONNECTOME DATA FROM AMSTERDAM ################# 

parcel_disc_indv_100_ctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex114Parcels_individuals_ctx.csv")
parcel_disc_indv_400_ctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex414Parcels_individuals_ctx.csv")

parcel_disc_indv_100_sctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex114Parcels_individuals_sctx.csv")
parcel_disc_indv_400_sctx = pd.read_csv(f"{discdir}/parcel_disc_SchaeferSubcortex414Parcels_individuals_sctx.csv")

#disc_indv_400_ctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex414Parcels_individuals_ctx.csv")
#disc_indv_400_sctx = pd.read_csv(f"{discdir}/PrograMS_disc_SchaeferSubcortex414Parcels_individuals_sctx.csv")

# ################# clean disc data ################# 

disc_indv_100_ctx = split_and_clean_header(parcel_disc_indv_100_ctx)
disc_indv_100_sctx = split_and_clean_header(parcel_disc_indv_100_sctx)

disc_indv_400_ctx = split_and_clean_header(parcel_disc_indv_400_ctx)
disc_indv_400_sctx = split_and_clean_header(parcel_disc_indv_400_sctx)

#disc_indv_400_ctx = split_and_clean_header(disc_indv_400_ctx)
#disc_indv_400_sctx = split_and_clean_header(disc_indv_400_sctx)

#%%
# ################# OUTPUT DIRECTORY ################# 
outdir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/june/cross_indv_results"

#%%
#==============================================================================
# AXONAL DISCONNECTOME 
#==============================================================================

axonal_disc_d100={}

for indv, row in atrophy_long_indv_100.iterrows():
    
    print(indv)
    #select subject from neighbour data 
    disc_subj_100_ctx = disc_indv_100_ctx.loc[disc_indv_100_ctx.index.astype(str) == str(indv)]
    disc_subj_100_sctx = disc_indv_100_sctx.loc[disc_indv_100_sctx.index.astype(str) == str(indv)]
    
    ##convert data to matrix format 
    #values = disc_subj_100_ctx.values  
    #matrix = np.full((100, 100), np.nan)
    #triu_indices = np.triu_indices(100, k=1)  # k=1 excludes diagonal
    #matrix[triu_indices] = values
    #matrix[(triu_indices[1], triu_indices[0])] = values # Mirror to lower-triangle for symmetry
    #np.fill_diagonal(matrix, 0)
    #matrix_df = pd.DataFrame(matrix)
    
    ##convert data to matrix format 
    #values_sctx  = disc_subj_100_sctx.values 
    #matrix_sctx = values_sctx.reshape(14, 100)  
    #matrix_df_sctx  = pd.DataFrame(matrix_sctx)
    
    #compute strenght
    #disc_subj_100_ctx_dc, disc_subj_100_sctx_dc = compute_dc(matrix_df, matrix_df_sctx)
    
    row = row.to_frame().T 
    rvals, p_and_d = axonal_disc(disc_subj_100_ctx.T.squeeze(), disc_subj_100_sctx.T.squeeze(), row, "schaefer_100")
    axonal_disc_d100[indv] = {'rvals': rvals, 'p_and_d': p_and_d}

    print(f"{indv} done!")
    
#%%

axonal_disc_d100_df =  pd.DataFrame(axonal_disc_d100.items()) 
np.save(f'{outdir}/2025-06-18_axonal-disc-parcel-indv_Schaefer100_rot-10000.npy', axonal_disc_d100) 
    
#%%

axonal_disc_d400={}

for indv, row in atrophy_long_indv_400.iterrows():
    
    print(indv)
    #select subject from neighbour data 
    disc_subj_400_ctx = disc_indv_400_ctx.loc[disc_indv_400_ctx.index.astype(str) == str(indv)]
    disc_subj_400_sctx = disc_indv_400_sctx.loc[disc_indv_400_sctx.index.astype(str) == str(indv)]
    
    ##convert data to matrix format 
    #values = disc_subj_400_ctx.values  
    #matrix = np.full((400, 400), np.nan)
    #triu_indices = np.triu_indices(400, k=1)  # k=1 excludes diagonal
    #matrix[triu_indices] = values
    #matrix[(triu_indices[1], triu_indices[0])] = values # Mirror to lower-triangle for symmetry
    #np.fill_diagonal(matrix, 0)
    #matrix_df = pd.DataFrame(matrix)
    
    ##convert data to matrix format 
    #values_sctx  = disc_subj_400_sctx.values 
    #matrix_sctx = values_sctx.reshape(14, 400)  
    #matrix_df_sctx  = pd.DataFrame(matrix_sctx)
    
    #compute strenght
    #disc_subj_400_ctx_dc, disc_subj_400_sctx_dc = compute_dc(matrix_df, matrix_df_sctx)
    
    row = row.to_frame().T 
    rvals, p_and_d = axonal_disc(disc_subj_400_ctx.T.squeeze(), disc_subj_400_sctx.T.squeeze(), row, "schaefer_400")
    axonal_disc_d400[indv] = {'rvals': rvals, 'p_and_d': p_and_d}

    print(f"{indv} done!")
    
#%%
axonal_disc_d400_df =  pd.DataFrame(axonal_disc_d400.items())    
np.save(f'{outdir}/2025-06-18_axonal-disc-parcel-indv_Schaefer400_rot-10000.npy', axonal_disc_d400_df) 
