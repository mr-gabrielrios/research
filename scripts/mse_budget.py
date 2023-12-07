''' Import packages. '''
# Analytical packages
import numpy as np, pandas as pd, xarray as xr
# Utilities
import os
# Visualization
import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt
# Internal methods
import colormap_normalization as cn
import tc_processing
import utilities

def budget_timeseries(data, budget_type='mse', level=500):
    
    ''' Get timeseries of a chosen MSE-based (let MSE = h) budget. 
    1. Option 'h': h = c_p*T + g*z + L_v*q --> this is given at a specific level
    2. Option 'dh_dt': dh_dt = (shflx + lhflx) + net_lw + net_sw + div(uh)
    '''
    
    # Get relevant constants
    c_p = utilities.get_constants('c_p')
    L_v = utilities.get_constants('L_v')
    R_d = utilities.get_constants('R_d')
    eps = utilities.get_constants('eps')
    g = utilities.get_constants('g')
    
    if budget_type == 'mse':
        # Initialize budget dictionary
        budget = {'type': None, 'data': {}, 'attrs': None}
        # Pre-populate data subdictionary with timestamp keys for proper time-data pairing
        budget['data'] = {'c_p T': {}, 'gz': {}, 'L_v q': {}}
        # Iterate over each timestamp
        for timestamp in data['tc_vertical_output'].time.values:
            # Shorthand for tc_vertical_output
            iterdata = data['tc_vertical_output'].sel(time=timestamp).sel(pfull=level, method='nearest').dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            budget['data']['c_p T'][timestamp] = c_p*iterdata['temp']
            budget['data']['L_v q'][timestamp] = L_v*iterdata['sphum']
            budget['data']['gz'][timestamp] = (100*iterdata.pfull)/iterdata['rho']
        # Append budget metadata
        budget['type'] = budget_type
        budget['attrs'] = data['tc_vertical_output'][budget_type].attrs
        
    return budget

def budget_timeseries_visualization(budget):
    
    # Select point for analysis - defaults to domain center
    util_key = [k for k in budget['data'].keys()][0]
    util_value = [v for v in budget['data'][util_key].values()][0] # get term to obtain domain dimensions
    # Pick domain center
    x, y = [len(util_value['grid_xt']) // 2, len(util_value['grid_yt']) // 2]

    fig, ax = plt.subplots(figsize=(5, 2))
    # Iterate over each term
    for term in budget['data'].keys():
        # Plot date on x-axis, term value on y-axis  
        ax.plot([k.isoformat() for k in budget['data'][term].keys()],
                [v.isel(grid_xt=x, grid_yt=y) for v in budget['data']term].values()], label='{0}'.format(term))
    
    # Format date
    ax.xaxis.set_major_formatter(matplotlib.dates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    # Format legend
    fig.legend(ncols=len(budget.keys()), frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.1))
    
def budget_planar_visualization(data):
    
    fields = {'wind': None, 'thflx': None, 'net_lw': None, 'net_sw': None}
    for field in fields.keys():
        fields[field] = data[field].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    
    nrows, ncols = 2, len(fields.keys())
    fig = plt.figure(figsize=(3*ncols, 3))
    gs = matplotlib.gridspec.GridSpec(ncols=ncols, nrows=nrows, height_ratios=(1, 0.1))
    
    axes = {}
    
    for i, item in enumerate(fields.items()):
        field, values = item
        axes[field] = fig.add_subplot(gs[0, i])
        axes[field].set_title(fields[field].attrs['long_name'])
        
        # Get normalization and colormap for the field of interest
        norm, cmap = cn.norm_cmap(data=fields, param=field)
        # Plot field of interest
        im = axes[field].pcolormesh(fields[field].grid_xt, fields[field].grid_yt, fields[field],
                                    norm=norm, cmap=cmap)
        # Plot storm center
        axes[field].scatter(data['center_lon'], data['center_lat'], s=50, c='r', ec='k')
        
        ''' Colorbar. '''
        # Assign axis
        axes['colorbar_{0}'.format(field)] = fig.add_subplot(gs[1, i])
        axes['colorbar_{0}'.format(field)].set_axis_off()
        # Create shorthand for the axis name
        cax = axes['colorbar_{0}'.format(field)].inset_axes([0, 0, 1, 0.5])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                orientation='horizontal', cax=cax)
        
    # Plot storm metadata at given time
    fig.suptitle(data.time.values, y=1.1)
        
if __name__ == '__main__':
    data = tc_processing.main('HIRAM', 'C15w', storm_id='2052-0030')
    
    budget = budget_timeseries(data)
    
    # for time in range(0, len(data['tc_model_output'].time)):    
    #     budget_planar_visualization(data['tc_model_output'].isel(time=time))