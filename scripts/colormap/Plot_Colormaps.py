#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 12:53:16 2025

@author: mbarrantescepas
"""

import matplotlib.pyplot as plt 

#plot color map 
def plot_cmap(title):
    fig, ax = plt.subplots(figsize=(0.25, 4))  # Wider figure for better visibility
    fig.subplots_adjust(left=0.2)  # Reduce left margin
    
    cmap = plt.cm.Blues
    norm = plt.Normalize(vmin=0, vmax=0.3)

    colorbar = plt.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=ax, orientation='vertical'
    )
    
    colorbar.ax.set_title(title, fontsize=12, pad=20)
    
    # Show numbers (ticks) on color bar
    colorbar.ax.tick_params(labelsize=12)
    colorbar.set_ticks([norm.vmin, norm.vmax])
    colorbar.ax.set_yticklabels([f"{norm.vmin}", f"{norm.vmax}"])
    
    # Place vmax at the top
    colorbar.ax.yaxis.set_ticks_position('right')
    
    # Show plot
    plt.show()
    
plot_cmap("")

#%%

fig, ax = plt.subplots(figsize=(0.25, 4))

# Plot the data on the primary axis
cmap = plt.cm.RdBu_r
right_norm2 = plt.Normalize(vmin=-10, vmax=10)
right_norm = plt.Normalize(vmin=-15, vmax=15)
ax2 = ax.twinx()

# Add colorbar for the primary axis
colorbar1 = plt.colorbar(plt.cm.ScalarMappable(norm=right_norm, cmap=cmap), cax=ax, orientation='vertical')
colorbar2 = plt.colorbar(plt.cm.ScalarMappable(norm=right_norm2, cmap=cmap), cax=ax2, orientation='vertical')

colorbar1.ax.set_title('Cortical | Subcortical \n (t-values)', fontsize=12, pad=20)
colorbar1.ax.tick_params(labelsize=12)
colorbar1.set_ticks([right_norm.vmin, right_norm.vmax])
colorbar2.ax.tick_params(labelsize=12)
colorbar2.ax.yaxis.set_ticks_position('left')
colorbar2.set_ticks([right_norm2.vmin, right_norm2.vmax])


plt.show()