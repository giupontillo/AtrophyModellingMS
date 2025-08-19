#!/data/anw/knw-work/g.pontillo/predictive_modelling/python-env/atrophymodelling/bin/python

#authors: m.barrantescepas@amsterdamumc.nl
#updates: 18 june 2025

#%%
import os
import argparse 
import numpy as np
import pandas as pd 
from io import StringIO
from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.permutation_testing import spin_test#, shuf_test

#%%
def load_csv_to_numpy(file_mapping, base_path):
    "reads multiple csv files into numpy arrays"
    arrays = {}
    for key, filename in file_mapping.items():
        arrays[key] = pd.read_csv(f"{base_path}/{filename}", header=None).to_numpy()
    return arrays

# Identify cortical epicenters (from structural/functional connectivity)
def cortical_epicenters(fsc_ctx, atrophy_ctx): 
    'function to identify cortical epicenters from structural/functional connectivity'
    fsc_ctx_epi = []
    fsc_ctx_epi_p = []
    for seed in range(fsc_ctx.shape[0]):
        seed_con = fsc_ctx[:, seed]
        fsc_ctx_epi = np.append(fsc_ctx_epi, np.corrcoef(seed_con, atrophy_ctx)[0, 1])
        fsc_ctx_epi_p = np.append(fsc_ctx_epi_p,
                                 spin_test(seed_con, atrophy_ctx, 
                                           surface_name='fsa5', parcellation_name='schaefer_100',
                                           type='pearson', n_rot=10000, null_dist=False))
    
    return fsc_ctx_epi, fsc_ctx_epi_p 

# Identify subcortical epicenters (from structural/functional connectivity)
def subcortical_epicenters(fsc_sctx, atrophy_ctx): 
    'function to identify subcortical epicenters from structural/functional connectivity'
    fsc_sctx_epi = []
    fsc_sctx_epi_p = []
    for seed in range(fsc_sctx.shape[0]):
        seed_con = fsc_sctx[seed, :]
        fsc_sctx_epi = np.append(fsc_sctx_epi, np.corrcoef(seed_con, atrophy_ctx)[0, 1])
        fsc_sctx_epi_p = np.append(fsc_sctx_epi_p,
                                  spin_test(seed_con, atrophy_ctx,
                                            surface_name='fsa5', n_rot=10000))
    
    return fsc_sctx_epi, fsc_sctx_epi_p

def save_epicenters(data, output_folder, subjectid):
    os.makedirs(output_folder, exist_ok=True)
    for var_name, array in data.items():
        output_path = os.path.join(output_folder, f"{subjectid}_{var_name}.csv")
        np.savetxt(output_path, array, delimiter=",")
        print(f"Saved {var_name} to {output_path}")

#%%
#----------------------------------------------------------------------
# IMPORT PROGRAMS + NORMATIVE DATA
#----------------------------------------------------------------------
sc_ctx_100, sc_ctx_labels_100, sc_sctx_100, sc_sctx_labels_100 = load_sc(parcellation="schaefer_100") 
fc_ctx_100, fc_ctx_labels_100, fc_sctx_100, fc_sctx_labels_100 = load_fc(parcellation="schaefer_100") 
sc_ctx_400, sc_ctx_labels_400, sc_sctx_400, sc_sctx_labels_400 = load_sc(parcellation="schaefer_400") 
fc_ctx_400, fc_ctx_labels_400, fc_sctx_400, fc_sctx_labels_400 = load_fc(parcellation="schaefer_400") 

# structural subcortical are [14,X14], remove 14 subcortical regions 
sc_sctx_100 = sc_sctx_100[:,:-14]
sc_sctx_400 = sc_sctx_400[:,:-14]

#%%
#Load atrophy data from programs
dir = "/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling"
atrophydir = f'{dir}/data/atrophy_cross'

parser = argparse.ArgumentParser()
parser.add_argument('--row', required=True, help='Full CSV row')
args = parser.parse_args()

#%%
#atrophy_d100 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex114Parcels_individuals_wscore.csv")
#atrophy_d400 = pd.read_csv(f"{atrophydir}/atrophy_baseline_SchaeferSubcortex414Parcels_individuals_wscore.csv")
#atrophy_d100.columns = range(atrophy_d100.shape[1])
#atrophy_d100 = atrophy_d100.iloc[[0]]

#%%

atrophy_d100=pd.read_csv(StringIO(args.row), header=None)

atrophy_d100.set_index(atrophy_d100.columns[0], inplace= True)
print(atrophy_d100)
subject_id = str(atrophy_d100.index[0])

#%%
#split atrophy maps in cortical and subcortical regions 
atrophy_d100_ctx = atrophy_d100.iloc[:,:-14].to_numpy().T.squeeze()  #remove sctx
atrophy_d100_sctx = atrophy_d100.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 

# cortical epicenters from normative data 
fc_ctx_epi_d100,  fc_ctx_epi_p_d100 = cortical_epicenters(fc_ctx_100, atrophy_d100_ctx)
sc_ctx_epi_d100,  sc_ctx_epi_p_d100 = cortical_epicenters(sc_ctx_100, atrophy_d100_ctx)

# subcortical epicenters from normative data 
fc_sctx_epi_d100,  fc_sctx_epi_p_d100 = subcortical_epicenters(fc_sctx_100, atrophy_d100_ctx)
sc_sctx_epi_d100,  sc_sctx_epi_p_d100 = subcortical_epicenters(sc_sctx_100, atrophy_d100_ctx)

# save data into csv files
cortical_epicenters_data = {
    "fc_ctx_epi_d100": fc_ctx_epi_d100,
    "fc_ctx_epi_p_d100": fc_ctx_epi_p_d100,
    "sc_ctx_epi_d100": sc_ctx_epi_d100,
    "sc_ctx_epi_p_d100": sc_ctx_epi_p_d100}

subcortical_epicenters_data = {
    "fc_sctx_epi_d100": fc_sctx_epi_d100,
    "fc_sctx_epi_p_d100": fc_sctx_epi_p_d100,
    "sc_sctx_epi_d100": sc_sctx_epi_d100,
    "sc_sctx_epi_p_d100": sc_sctx_epi_p_d100}

output_folder = f"{dir}/output/june/cross_indv_epi_results/"
save_epicenters(cortical_epicenters_data, output_folder, subject_id) 
save_epicenters(subcortical_epicenters_data, output_folder, subject_id) 

#%%

# atrophy_d400_ctx = atrophy_d400.iloc[:,:-14].to_numpy().T.squeeze()  #remove sctx
# atrophy_d400_sctx = atrophy_d400.iloc[:,-14:].to_numpy().T.squeeze()  #remove cortical 

# fc_ctx_epi_d400,  fc_ctx_epi_p_d400 = cortical_epicenters(fc_ctx_400, atrophy_d400_ctx)
# sc_ctx_epi_d400,  sc_ctx_epi_p_d400 = cortical_epicenters(sc_ctx_400, atrophy_d400_ctx)

# fc_sctx_epi_d400,  fc_sctx_epi_p_d400 = subcortical_epicenters(fc_sctx_400, atrophy_d400_ctx)
# sc_sctx_epi_d400,  sc_sctx_epi_p_d400 = subcortical_epicenters(sc_sctx_400, atrophy_d400_ctx)