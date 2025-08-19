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

#%%
datadir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data"
demographics = pd.read_csv(f'{datadir}/df_baseline.csv')
subset = demographics[['site', 'id']]

subset = subset.rename(columns={'id': 'Subject'})


#%%
#read data to plot it 
outdir="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/output_may/indv_cross_results_v2025.05"

nodalstress_100 = np.load(f'{outdir}/2025-05-24_nodal-stress-indv_Schaefer100.npy',allow_pickle='TRUE').item()

#transform npy to same structure as nodal stress 
neigh_100 = np.load(f'{outdir}/2025-05-24_neighbour-atrophy-indv_Schaefer100.npy', allow_pickle='TRUE')
neigh_100 = pd.DataFrame(neigh_100, columns=['id', 'values'])
neigh_100 = pd.Series(neigh_100['values'].values, index=neigh_100['id']).to_dict()

trasnvul_100 = np.load(f'{outdir}/2025-05-24_transvulnerability-indv_Schaefer100.npy', allow_pickle='TRUE')
trasnvul_100 = pd.DataFrame(trasnvul_100, columns=['id', 'values'])
trasnvul_100 = pd.Series(trasnvul_100['values'].values, index=trasnvul_100['id']).to_dict()

#nodalstress_100 = np.load(f'{outdir}/2025-03-13_axonal-disc-indv_Schaefer100.npy',allow_pickle='TRUE').item()
#nodalstress_100 = np.load(f'{outdir}/2025-03-13_transvulnerability-indv_Schaefer100.npy', allow_pickle='TRUE')

#%%
data=trasnvul_100

subjects = list(data.keys())
categories = list(data[subjects[0]]['rvals'].keys())

# Prepare a DataFrame for easier seaborn plotting
# We'll create one DataFrame for r-values and one for p-values

rvals_list = []
pvals_list = []

for subject in subjects:
    for category in categories:
        rvals_list.append({
            'Subject': subject,
            'Category': category,
            'R_value': data[subject]['rvals'][category]
        })
        pvals_list.append({
            'Subject': subject,
            'Category': category,
            'P_value': data[subject]['p_and_d'][category][0]  # p-values first element
        })

df_rvals = pd.DataFrame(rvals_list)
df_pvals = pd.DataFrame(pvals_list)


df_rvals = df_rvals.merge(subset, on='Subject', how='left')
df_pvals = df_pvals.merge(subset, on='Subject', how='left')

df_rvals = df_rvals.sample(frac=1).reset_index(drop=True)
df_pvals = df_pvals.sample(frac=1).reset_index(drop=True)

#%%

unique_sites = df_rvals['site'].unique()
base_palette = sns.color_palette("rocket_r", n_colors=len(unique_sites))
site_palette = dict(zip(unique_sites, base_palette))

n_categories = len(categories)
fig, axes = plt.subplots(1, n_categories, figsize=(4 * n_categories, 5), sharey=True)

for i, category in enumerate(categories):
    ax = axes[i]
    sns.scatterplot(data=df_rvals[df_rvals['Category'] == category],
        x='Subject', y='R_value', hue='site', palette=site_palette, s=15, ax=ax)
    ax.axhline(0.0, color='black', linestyle='--', linewidth=0.5)
    ax.set_ylim(-1, 1)
    ax.set_xticks([])  
    ax.set_title(category, fontsize=10)
    ax.legend(loc='upper right', fontsize='small', frameon=True)
    
    if i == 0:
        ax.set_ylabel('r-values')
    else:
        ax.set_ylabel('')

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()

#%%
i=0
n_categories = len(categories)
fig, axes = plt.subplots(1, n_categories, figsize=(4 * n_categories, 5), sharey=True)

for i, category in enumerate(categories):
    ax = axes[i]
    sns.scatterplot(data=df_pvals[df_pvals['Category'] == category],
        x='Subject', y='P_value', hue='site', palette=site_palette, s=15, ax=ax)
    ax.axhline(0.05, color='black', linestyle='--', linewidth=0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([])  
    ax.set_title(category, fontsize=10)
    ax.legend(loc='upper right', fontsize='small', frameon=True)
    
    if i == 0:
        ax.set_ylabel('p-values')
    else:
        ax.set_ylabel('')

plt.tight_layout(rect=[0, 0, 0.85, 1])
plt.show()