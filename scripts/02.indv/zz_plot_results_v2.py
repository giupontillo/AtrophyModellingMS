#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 02:40:12 2025

@author: mbarrantescepas
"""

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_1samp, ttest_ind, f_oneway

#%%
#==============================================================================
# LOAD DEMOGRAPHIC DATA  
#==============================================================================

datadir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data"

stages = pd.read_csv(f'{datadir}/stages.csv')
demographics = pd.read_csv(f'{datadir}/df_baseline.csv')
global_atrophy = pd.read_csv(f'{datadir}/globalatrophy.csv')
hydra_114 = pd.read_csv(f'{datadir}/hydra_SchaeferSubcortex114Parcels_2subtypes.csv')
#hydra_414 = pd.read_csv(f'{datadir}/hydra_SchaeferSubcortex414Parcels_2subtypes.csv')

#%%

hydra_114 = hydra_114.rename(columns={'participant_id': 'id', 'session_id': 'session' })
#hydra_414 = hydra_414.rename(columns={'participant_id': 'id', 'session_id': 'session' })

#%%
#==============================================================================
# LOAD INDIVIDUAL DATA 
#==============================================================================

indvdir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output/june/cross_indv_results"

#nodalstress
df_nodalstr_100 = np.load(f'{indvdir}/2025-05-24_nodal-stress-indv_Schaefer100.npy',allow_pickle='TRUE').item()

#neighbour atrophy
df_neigh_100 = np.load(f'{indvdir}/2025-05-24_neighbour-atrophy-indv_Schaefer100.npy', allow_pickle='TRUE')
df_neigh_100 = pd.DataFrame(df_neigh_100, columns=['id', 'values'])
df_neigh_100 = pd.Series(df_neigh_100['values'].values, index=df_neigh_100['id']).to_dict()

#axonal disconnection 
df_axdisc_100 = np.load(f'{indvdir}/2025-06-18_axonal-disc-parcel-indv_Schaefer100_rot-10000.npy',allow_pickle='TRUE').item() #new

#trasncriptomic vulnerability
df_trasnvul_100 = np.load(f'{indvdir}/2025-06-18_transvulnerability-indv_Schaefer100_rot-10000.npy', allow_pickle='TRUE') #new
df_trasnvul_100 = pd.DataFrame(df_trasnvul_100, columns=['id', 'values'])
df_trasnvul_100 = pd.Series(df_trasnvul_100['values'].values, index=df_trasnvul_100['id']).to_dict()

#%%
#==============================================================================
# combine all data same dataset and plot
#==============================================================================

datasets=[df_nodalstr_100, df_neigh_100, df_axdisc_100, df_trasnvul_100]

rvals_list = []

for data in datasets: 
    subjects = list(data.keys())
    categories = list(data[subjects[0]]['rvals'].keys())

    for subject in subjects:
        for category in categories:
            rvals_list.append({
                'id': subject,
                'category': category,
                'r_value': data[subject]['rvals'][category], 
                'p_value': data[subject]['p_and_d'][category][0]  # p-values first element
            })
        

df_rvals = pd.DataFrame(rvals_list)
df_rvals['z_value'] = np.arctanh(df_rvals['r_value'])

#%%

#merge global atrophy & stages  data 
df_rvals = df_rvals.merge(global_atrophy.drop(columns=["session"]), on='id', how='left')
df_rvals = df_rvals.merge(stages.drop(columns=["session"]), on='id', how='left')

#merge hydra data 
df_rvals = df_rvals.merge(hydra_114, on='id', how='left')
df_rvals = df_rvals.sample(frac=1).reset_index(drop=True)
df_rvals = df_rvals.dropna(subset=['assignment_2'])

#%%
#save dataframe
df_rvals.to_csv(f"{indvdir}/2025.06.23_all_indv_rvals.csv", index=False)

#%%

#categories = df_rvals["category"].unique()
categories = [ "functional cortical hubs", "structural cortical hubs", 
              "functional subcortical hubs", "structural subcortical hubs", 
              "functional cortical neighbour atrophy", "structural cortical neighbour atrophy", 
              "functional subcortical neighbour atrophy", "structural subcortical neighbour atrophy", 
              "structural cortical disconnectome", "structural subcortical disconnectome", 
              "cortical gene expression", "subcortical gene expression"]

#%%
results = []
groups = df_rvals['assignment_2'].dropna().unique()

for category in categories:
    df_cat = df_rvals[df_rvals['category'] == category]
    
    # One-sample t-test against 0
    z_all = df_cat['z_value'].dropna()
    t1, p1 = ttest_1samp(z_all, popmean=0)

    # Two-sample t-test (only if exactly 2 groups)
    if len(groups) == 2:
        z1 = df_cat[df_cat['assignment_2'] == groups[0]]['z_value'].dropna()
        z2 = df_cat[df_cat['assignment_2'] == groups[1]]['z_value'].dropna()
        t2, p2 = ttest_ind(z1, z2, equal_var=False)  # Welch’s test
    else:
        t2, p2 = None, None  # Skip if not two groups

    # Save results
    results.append({
        'Category': category,
        'OneSample_t': t1,
        'OneSample_p': p1,
        'TwoSample_t': t2,
        'TwoSample_p': p2,
        'N_total': len(z_all),
        f'N_{groups[0]}': len(z1) if len(groups) >= 1 else None,
        f'N_{groups[1]}': len(z2) if len(groups) >= 2 else None,
    })

df_stats = pd.DataFrame(results)

# Optional: adjust p-values for multiple comparisons (e.g., FDR)
from statsmodels.stats.multitest import multipletests
df_stats['OneSample_p_FDR'] = multipletests(df_stats['OneSample_p'], method='fdr_bh')[1]
df_stats['TwoSample_p_FDR'] = multipletests(df_stats['TwoSample_p'].dropna(), method='fdr_bh')[1].tolist() + [None] * (len(df_stats) - len(df_stats['TwoSample_p'].dropna()))

print(df_stats)

#%%

results = []
groups = df_rvals['globalatrophy'].dropna().unique()

for category in categories:
    df_cat = df_rvals[df_rvals['category'] == category]
    
    # One-sample t-test against 0
    z_all = df_cat['z_value'].dropna()
    t1, p1 = ttest_1samp(z_all, popmean=0)

    # Two-sample t-test (only if exactly 2 groups)
    if len(groups) == 2:
        z1 = df_cat[df_cat['globalatrophy'] == groups[0]]['z_value'].dropna()
        z2 = df_cat[df_cat['globalatrophy'] == groups[1]]['z_value'].dropna()
        t2, p2 = ttest_ind(z1, z2, equal_var=False)  # Welch’s test
    else:
        t2, p2 = None, None  # Skip if not two groups

    # Save results
    results.append({
        'Category': category,
        'OneSample_t': t1,
        'OneSample_p': p1,
        'TwoSample_t': t2,
        'TwoSample_p': p2,
        'N_total': len(z_all),
        f'N_{groups[0]}': len(z1) if len(groups) >= 1 else None,
        f'N_{groups[1]}': len(z2) if len(groups) >= 2 else None,
    })

df_stats = pd.DataFrame(results)

# Optional: adjust p-values for multiple comparisons (e.g., FDR)
from statsmodels.stats.multitest import multipletests
df_stats['OneSample_p_FDR'] = multipletests(df_stats['OneSample_p'], method='fdr_bh')[1]
df_stats['TwoSample_p_FDR'] = multipletests(df_stats['TwoSample_p'].dropna(), method='fdr_bh')[1].tolist() + [None] * (len(df_stats) - len(df_stats['TwoSample_p'].dropna()))

print(df_stats)
    
#%%
sites = ['prague', 'amsterdam', 'graz']
palette_colors = sns.color_palette("rocket_r", n_colors=len(sites))
site_palette = {
    "amsterdam": palette_colors[2],  
    "prague": palette_colors[0],      
    "graz": palette_colors[1]
}
#%%
unique_sites = df_rvals['assignment_2'].unique()
base_palette = sns.color_palette("rocket_r", n_colors=len(unique_sites))
site_palette = dict(zip(unique_sites, base_palette))


n_categories = len(categories)
n_cols = 4
n_rows = int(np.ceil(n_categories / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 5 * n_rows), sharey=True)
axes=axes.flatten()

for i, category in enumerate(categories):
    print(category)
    ax = axes[i]
    sns.violinplot(data=df_rvals[df_rvals['category'] == category],
        x='assignment_2', y='z_value', hue='assignment_2', palette=site_palette, ax=ax)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.5)
    ax.set_ylim(-2, 2)
    ax.set_xticks([])  
    ax.set_title(category, fontsize=10)
    ax.legend(loc='upper right', fontsize='small', frameon=True)
    
    if i % n_cols == 0:
        ax.set_ylabel('Fisher z')
    else:
        ax.set_ylabel('hydra cluster')
        
        
for j in range(n_categories, n_rows * n_cols):
    fig.delaxes(axes[j])


plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

#%%

n_categories = len(categories)

n_cols = 4
n_rows = int(np.ceil(n_categories / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 5 * n_rows), sharey=True)
axes=axes.flatten()

for i, category in enumerate(categories):
    print(category)
    ax = axes[i]
    sns.violinplot(data=df_rvals[df_rvals['category'] == category],
        x='globalatrophy', y='z_value', hue='globalatrophy', palette=site_palette, ax=ax)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.5)
    ax.set_ylim(-2, 2)
    ax.set_xticks([])  
    ax.set_title(category, fontsize=10)
    ax.legend(loc='upper right', fontsize='small', frameon=True)
    
    if i % n_cols == 0:
        ax.set_ylabel('Fisher z')
    else:
        ax.set_ylabel('hydra cluster')
        
        
for j in range(n_categories, n_rows * n_cols):
    fig.delaxes(axes[j])


plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

#%%

sites = ['LR', 'SP', 'PP', "ER"]
palette_colors = sns.color_palette("rocket_r", n_colors=len(sites))
site_palette = {
    "LR": palette_colors[2],  
    "SP": palette_colors[0],      
    "PP": palette_colors[1],
    "ER": palette_colors[3]
}

n_categories = len(categories)

n_cols = 4
n_rows = int(np.ceil(n_categories / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 5 * n_rows), sharey=True)
axes=axes.flatten()

for i, category in enumerate(categories):
    print(category)
    ax = axes[i]
    sns.boxplot(data=df_rvals[df_rvals['category'] == category],
        x='stage', y='z_value', hue='stage', palette=site_palette, ax=ax)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.5)
    ax.set_ylim(-2, 2)
    ax.set_xticks([])  
    ax.set_title(category, fontsize=10)
    ax.legend(loc='upper right', fontsize='small', frameon=True)
    
    if i % n_cols == 0:
        ax.set_ylabel('Fisher z')
    else:
        ax.set_ylabel('hydra cluster')
        
        
#for j in range(n_categories, n_rows * n_cols):
 #   fig.delaxes(axes[j])


plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()