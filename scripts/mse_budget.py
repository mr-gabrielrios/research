import numpy as np, pandas as pd, xarray as xr

import os

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

def budget_planar_visualization():
    
    nrows, ncols = 2, 4
    fig = plt.figure(figsize=(2*ncols, 3))
    gs = matplotlib.gridspec.GridSpec(ncols=ncols, nrows=nrows, height_ratios=(1, 0.1))
    
    fields = {'thf': None, 'net_lw': None, 'net_sw': None, 'conv': None}
    axes = {}
    
    for i, item in enumerate(fields.items()):
        key, value = item
        axes[key] = fig.add_subplot(gs[0, i])
        axes[key].set_title(key)
        
        ''' Colorbar. '''
        # Assign axis
        axes['colorbar_{0}'.format(key)] = fig.add_subplot(gs[1, i])
        # Create shorthand for the axis name
        cax = axes['colorbar_{0}'.format(key)]
        
if __name__ == '__main__':
    budget_planar_visualization()