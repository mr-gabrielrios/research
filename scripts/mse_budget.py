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

def budget_timeseries(data, budget_type='h', level=500):
    
    ''' Get timeseries of a chosen MSE-based (let MSE = h) budget. 
    1. Option 'h': h = c_p*T + g*z + L_v*q --> this is given at a specific level
    2. Option 'dh_dt': dh_dt = (shflx + lhflx) + net_lw + net_sw + div(uh)
    '''
    
    ''' NEED TO INTEGRATE INTO XARRAY DATASET FOR EASIER ACCESSING'''
    
    # Get relevant constants
    c_p = utilities.get_constants('c_p')
    L_v = utilities.get_constants('L_v')
    R_d = utilities.get_constants('R_d')
    eps = utilities.get_constants('eps')
    g = utilities.get_constants('g')
    
    ''' Moist static energy (h). '''
    if budget_type == 'h':
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
        
    ''' Time tendency of moist static energy (dh_dt). '''
    # Define residual term - this determines which term to calculate as a residua
    residual_term = 'vi_flux_h' 
    if budget_type == 'dh_dt':
        # Initialize budget dictionary
        budget = {'type': None, 'data': {}, 'attrs': None}
        # Pre-populate data subdictionary with timestamp keys for proper time-data pairing
        if residual_term == 'vi_flux_h':
            budget['data'] = {'thflx': {}, 'net_lw': {}, 'net_sw': {}, 'dh_dt': {}}
        else:
            budget['data'] = {'thflx': {}, 'net_lw': {}, 'net_sw': {}, 'vi_flux_h': {}}
        # Initialize container array for the residual term
        residual_budget = {}
        # Iterate over each timestamp
        # If the residual term is not dh_dt, then pre-calculate it [(h[t] - h[t-1])/dt]
        if residual_term != 'dh_dt':
            for t, timestamp in enumerate(data['tc_model_output'].time.values):
                # Get h[t] (curr)
                curr = data['tc_model_output']['vi_h'].isel(time=t)
                # Get h[t-1] (prev) and dh_dt
                if t > 0:
                    prev = data['tc_model_output']['vi_h'].isel(time=(t-1))
                    num = curr.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all') - prev.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
                    den = (curr.time.values - prev.time.values).total_seconds()
                    budget['data']['dh_dt'][timestamp] = num/den
                    # fig, ax = plt.subplots()
                    # budget['data']['dh_dt'][timestamp].plot(ax=ax)
                # Fill the first timestep with nans
                else:
                    temp = xr.where(curr != -9999, 0, curr)
                    budget['data']['dh_dt'][timestamp] = temp
        # Iterate over each term
        for t, timestamp in enumerate(data['tc_model_output'].time.values):
            # Shorthand for tc_vertical_output
            iterdata = data['tc_model_output'].sel(time=timestamp)
            for k in budget['data'].keys():
                # Skip dh_dt, which is already in the budget
                if k == 'dh_dt':
                    continue
                else:
                    budget['data'][k][timestamp] = iterdata[k].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            # Calculate residual term
            if residual_term == 'vi_flux_h':
                residual_budget[timestamp] = budget['data']['net_lw'][timestamp] + budget['data']['net_sw'][timestamp] + \
                    budget['data']['thflx'][timestamp] - budget['data']['dh_dt'][timestamp]
            elif residual_term == 'dh_dt':
                residual_budget[timestamp] = budget['data']['net_lw'][timestamp] + budget['data']['net_sw'][timestamp] + \
                    budget['data']['thflx'][timestamp] - budget['data']['vi_flux_h'][timestamp]
        # Integrate the residual term into the budget
        budget['data'][residual_term] = residual_budget
        # Append budget metadata
        budget['type'] = budget_type
        
    return budget

def budget_timeseries_visualization(budget):
    
    # Select point for analysis - defaults to domain center
    util_key = [k for k in budget['data'].keys()][0] # get sample key
    util_value = [v for v in budget['data'][util_key].values()][0] # get sample data
    # Pick domain center
    center_x, center_y = [len(util_value['grid_xt']) // 2, len(util_value['grid_yt']) // 2]
    # Define extent for averaging (note: average grid cell is ~0.5 degrees, "5" chosen to target +/- 2.5 degrees from center)
    core_extent = 5
    # Perform averaging
    range_x, range_y = [slice(center_x - core_extent, center_x + core_extent), 
                        slice(center_y - core_extent, center_y + core_extent)]

    # Initialize figure
    fig, ax = plt.subplots(figsize=(6, 3))
    # Plot zero line
    ax.axhline(0, c='k')
    # Iterate over each term
    for term in budget['data'].keys():
        # Plot date on x-axis, term value on y-axis  
        ax.plot([k.isoformat() for k in budget['data'][term].keys()],
                [v.isel(grid_xt=range_x, grid_yt=range_y).mean().values for v in budget['data'][term].values()], 
                label='{0}'.format(term), lw=2, marker='o')
    
    # Format date
    ax.xaxis.set_major_formatter(matplotlib.dates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
    # Format legend
    fig.legend(ncols=len(budget.keys()), frameon=False, loc='upper center', bbox_to_anchor=(0.5, 1.1))
    
def budget_planar_visualization(budget, time_index=0):
    
    ''' REDO FROM THE BUDGET DECOMP METHOD AVOBE'''
    
    print(budget.keys())
    
    fields = {'dh_dt': None, 'thflx': None, 'net_lw': None, 'net_sw': None, 'vi_flux_h': None}
    for field in fields.keys():
        print(field)
        for index, timestamp in enumerate(budget['data'][field].keys()):
            if index == time_index:
                fields[field] = budget['data'][field][timestamp].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    
    nrows, ncols = 2, len(fields.keys())
    fig = plt.figure(figsize=(3*ncols, 3))
    gs = matplotlib.gridspec.GridSpec(ncols=ncols, nrows=nrows, height_ratios=(1, 0.1))
    
    axes = {}
    
    for i, item in enumerate(fields.items()):
        field, values = item
        axes[field] = fig.add_subplot(gs[0, i])
        axes[field].set_title(field)
        
        # Get normalization and colormap for the field of interest
        norm, cmap = cn.norm_cmap(data=fields, param=field)
        # Plot field of interest
        im = axes[field].pcolormesh(fields[field].grid_xt, fields[field].grid_yt, fields[field],
                                    norm=norm, cmap=cmap)
        # Plot storm center
        # axes[field].scatter(data['center_lon'], data['center_lat'], s=50, c='r', ec='k')
        
        ''' Colorbar. '''
        # Assign axis
        axes['colorbar_{0}'.format(field)] = fig.add_subplot(gs[1, i])
        axes['colorbar_{0}'.format(field)].set_axis_off()
        # Create shorthand for the axis name
        cax = axes['colorbar_{0}'.format(field)].inset_axes([0, 0, 1, 0.5])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                orientation='horizontal', cax=cax)
        
    # Plot storm metadata at given time
    # fig.suptitle(data.time.values, y=1.1)
        
if __name__ == '__main__':
    print('')
    data = tc_processing.main('HIRAM', 'C15w', storm_id='2056-0039')
    
    budget = budget_timeseries(data, budget_type='dh_dt')
    
    