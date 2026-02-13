#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  1 11:18:16 2025

@author: mbarrantescepas
"""


import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from enigmatoolbox.datasets import load_sc, load_fc
from enigmatoolbox.utils.parcellation import parcel_to_surface
from enigmatoolbox.plotting import plot_cortical
from enigmatoolbox.plotting import plot_subcortical

#%%

data="/data/anw/knw-work/g.pontillo/predictive_modelling/code/stats/AtrophyModelling/data"

sc_ctx_100 = pd.read_csv(f"{data}/hcp_connectivity/strucMatrix_ctx_schaefer_100.csv", header=None)
fc_ctx_100 = pd.read_csv(f"{data}/hcp_connectivity/funcMatrix_ctx_schaefer_100.csv", header=None)
ge_100 = pd.read_csv(f"{data}/ge/allgenes_stable_r0.2_schaefer_100.csv", index_col="label")



#%%

# plot heat maps
plt.figure(figsize=(8,6))
sns.heatmap(sc_ctx_100, annot=False, cmap="Blues", xticklabels=False, yticklabels=False)
plt.show()

plt.figure(figsize=(8,6))
sns.heatmap(fc_ctx_100, annot=False, cmap="Reds", xticklabels=False, yticklabels=False)
plt.show()

#%%
plt.figure(figsize=(20,6))
sns.heatmap(ge_100, annot=False, cmap="viridis", xticklabels=False, yticklabels=False)
plt.show()