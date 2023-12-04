# %%

import numpy as np
import os
import pandas as pd
import scipy
import time
import xarray as xr

import cartopy, cartopy.crs as ccrs
from cartopy.util import add_cyclic_point
import matplotlib, matplotlib.pyplot as plt

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

from colormap_normalization import get_cmap, norm_cmap

def coords_to_dist(a, b):
    ''' Convert coordinates to distance in meters. '''
    
    R = 6371e3
    
    lon_a, lat_a = np.array(a)*np.pi/180
    lon_b, lat_b = np.array(b)*np.pi/180
    
    dlon, dlat = lon_b - lon_a, lat_b - lat_a
    
    a = np.sin(dlat/2)**2 + np.cos(lat_a)*np.cos(lat_b)*np.sin(dlon/2)**2    
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    distance = R*c
    
    return distance

def get_data(model, storm_type, experiments, years=None, nums=None):

    ''' Load TC output data from 'tc_storage_algorithm.ipynb'. If intensity-binned data exists, find and load it. Else, generate it. '''

    # Location where TC data is stored
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage'
    # Initialize container dictionaries
    data = {}
    # Define a storage boolean - True if binned data needs to be stored
    storage_bool = False
    # Iterate through experiments
    for experiment in experiments:
        num = nums[experiment]
        # Load data
        try:
            data[experiment] = np.load('{0}/tc_output.model-{1}.exp-{2}.storm_type-{3}.years-{4}_{5}.num-{6}.npy'.format(storage_dirname, model, experiment, storm_type, min(years), max(years), num),
                                      allow_pickle='TRUE').item()
        except:
            data[experiment] = np.load('{0}/tc_output.model-{1}.exp-{2}.storm_type-{3}.years-{4}_{5}.num-{6}.npy'.format(storage_dirname, model, experiment, storm_type, min(years), max(years), num),
                                      allow_pickle='TRUE')
        # Iterate through all files in the storage directory to see if pre-binned data already exists. Outputs a boolean.
        storage_check = [(('intensity_bins' in f) and (model in f)) for f in os.listdir('/projects/GEOCLIM/gr7610/analysis/tc_storage')]
        # If the boolean array sum is nonzero, data exists, so load it.
        # To-do: add method to stitch all files that meet this criteria to accommodate different years
        if sum(storage_check) > 0:
            data = np.load('{0}/tc_output.model-{1}.intensity_bins.storm_type-{2}.years-{3}_{4}.num-{5}.npy'.format(storage_dirname, model, storm_type, min(years), max(years), num), 
                           allow_pickle='TRUE').item()
            # Escape the loop
            break
        # Else, load
        else:
            data[experiment] = intensity_binning(data[experiment], benchmarking=True)
            storage_bool = True
    # Store binned data
    if storage_bool:
        output_path = '{0}/tc_output.model-{1}.intensity_bins.storm_type-{2}.years-{3}_{4}.num-{5}.npy'.format(storage_dirname, model, storm_type, min(years), max(years), num)
        try:
            np.save(output_path, data)
        except:
            import pickle
            with open(output_path, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
    
    return data

def intensity_binning(data, bins=None, benchmarking=False):

    lap = time.time()
    if benchmarking:
        checkpoint = 0
        print('\t Binning time, checkpoint {0}: {1:.3f} s'.format(checkpoint, (time.time() - lap)/1e9))
        checkpoint += 1
        lap = time.time()
        
    # Define surface-level TC-centric data and vertical Datasets with shorthand variable names
    planar, vertical = [data['tc_model_output'], data['tc_vertical_output']]
    # Gather the storm IDs from the planar Datasets
    storm_ids = {}
    # Iterate over unique storm IDs in the planar Dataset
    for storm_id in planar['storm_id'].values:
        # Drop all data not corresponding to the iterand storm ID
        storm = planar.where(planar['storm_id'] == storm_id).dropna(dim='time', how='all')
        # Add storm intensity data (use maximum winds as a proxy for intensity)
        storm_ids[storm_id] = storm['max_wind'].values
    # Rewrite storm IDs into a DataFrame for easier grouping
    storm_ids = pd.DataFrame.from_dict(storm_ids, orient='index', 
                                       columns=['max_wind']).reset_index(names=['storm_id'])
    
    if benchmarking:
        print('\t Binning time, checkpoint {0}: {1:.3f} s'.format(checkpoint, (time.time() - lap)/1e9))
        checkpoint += 1
        lap = time.time()
        
    # Define bins to follow those defined in Kim et al. (2019)
    bins = np.array([6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 35, 38, 41, 44, 47, 50])
    bins = np.array([6, 18, 27, 33, np.inf])
    # Initialize dictionary with binned data
    binned_storms = {'b{0}'.format(i): 
                     {'max_wind': None, 'planar': [], 'vertical': []}
                     for i, bin in enumerate(bins[:-1])}
    # Group the input data by intensity bin and iterate over each bins to find corresponding storms
    for i, grouper in enumerate(storm_ids.groupby(pd.cut(storm_ids['max_wind'], bins))):
        # Unpack Pandas grouper data
        g, v = grouper
        # Add wind range for the given bin
        binned_storms['b{0}'.format(i)]['max_wind'] = [g.left, g.right]
        # Within each bin, iterate over all matching storm IDs and extract relevant planar and vertical data
        for storm_id in v['storm_id'].values:
            # Filter out planar and vertical TC data corresponding to storm ID
            iter_planar = planar.where(planar['storm_id'] == storm_id).dropna(dim='time', how='all')
            iter_vertical = vertical.where(vertical['storm_id'] == storm_id).dropna(dim='time', how='all')
            # Add data to the dictionary
            binned_storms['b{0}'.format(i)]['planar'].append(iter_planar)
            binned_storms['b{0}'.format(i)]['vertical'].append(iter_vertical)
        
        if benchmarking:
            print('\t \t Binning time, checkpoint {0}: {1:.3f} s'.format(checkpoint, (time.time() - lap)/1e9))
            checkpoint += 1
            lap = time.time()

    return binned_storms

def planar_composite(ctrl, ktc, param='U', comparison_type='simple_difference'):
    
    # Visualization tools
    import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt
    from cartopy.util import add_cyclic_point
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    
    '''
    Algorithm to compare KillTC and control normalized datasets 
    for a given parameter using a selected comparison method.
    
    Inputs:
        - ctrl (dict):           Dictionary containing data for the control run.
        - ktc (dict):            Dictionary containing data for the KillTC run.
        - param (str):           String denoting the desired parameter for comparison.
        - comparison_type (str): Method for comparison. Default is 'simple_difference.'
    Outputs:
        None.
    '''
    
    # Get parameter long name for labeling
    try:
        long_name, units = [ctrl['tc_model_output'][param].attrs['long_name'], 
                            ctrl['tc_model_output'][param].attrs['units']]
    except:
        long_name, units = '', ''
    
    # Normalize datasets
    norm_ktc, norm_ctrl = [tc_normalization(ktc['tc_model_output'], var=param, test_num=None), 
                           tc_normalization(ctrl['tc_model_output'], var=param, test_num=None)]
    
    ''' Compare datasets. Defaults to a simple differencing method. '''
    # Relative difference (in other words, percent differene)
    if comparison_type == 'rel_diff':
        diff = 100*(norm_ktc.mean('time') - norm_ctrl.mean('time'))/(norm_ktc.mean('time') + norm_ctrl.mean('time'))
    elif comparison_type == 'std_anom':
    # Standardized anomaly differencing
        diff = (norm_ktc.mean('time') - norm_ctrl.mean('time'))/norm_ctrl.std('time')
    else:
        diff = norm_ktc.mean('time') - norm_ctrl.mean('time')
    
    # Define plot data
    plot_data = {'Control': norm_ctrl.mean('time'), 
                 'SWISHE': norm_ktc.mean('time'), 
                 '(SWISHE - control)': diff}
    
    # Restore parameter definition if divergence field calculated
    if 'div' in param:
        param = '{0}{1}'.format(div_param, level)
    
    ''' Define normalization (normzn); absolute and difference '''
    ### Absolute
    vmin_abs, vmax_abs = [np.nanmin([np.nanmin(norm_ctrl.mean('time')), np.nanmin(norm_ktc.mean('time'))]),
                          np.nanmax([np.nanmax(norm_ctrl.mean('time')), np.nanmax(norm_ktc.mean('time'))])]
    cmap = prop_dict(ctrl)[param]['cmap']
    # Choose colormap based on diverging or sequential values
    if ((vmin_abs < 0) & (vmax_abs > 0)):
        normzn_abs = matplotlib.colors.CenteredNorm()
        cmap_abs = cmap['div']
    else:
        cmap_abs = cmap['seq']
        normzn_abs = matplotlib.colors.Normalize(vmin=vmin_abs, vmax=vmax_abs) 
    
    ### Difference
    # Get extrema
    vmin_diff, vmax_diff = [np.nanmin(diff), np.nanmax(diff)]
    # Choose colormap based on diverging or sequential values
    if ((vmin_diff < 0) & (vmax_diff > 0)):
        if comparison_type == 'rel_diff':
            if vmin_diff < -100 or vmax_diff > 100:
                normzn_diff = matplotlib.colors.TwoSlopeNorm(vmin=-100, vmax=100, vcenter=0)
            else:
                normzn_diff = matplotlib.colors.CenteredNorm(vcenter=0)
            cmap_diff = cmap['div']
        elif comparison_type == 'std_anom':
            if vmin_diff < -3 or vmax_diff > 3:
                normzn_diff = matplotlib.colors.TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
            else:
                normzn_diff = matplotlib.colors.CenteredNorm(vcenter=0)
            cmap_diff = cmap['div']
        else:
            normzn_diff = matplotlib.colors.CenteredNorm()
            cmap_diff = cmap['div']
    else:
        cmap_diff = cmap_abs
        normzn_diff = matplotlib.colors.Normalize(vmin=vmin_diff, vmax=vmax_diff) 
    
    ''' Visualization. '''
    fig = plt.figure(figsize=(7, 3), layout='constrained')
    gs = matplotlib.gridspec.GridSpec(1, len(plot_data), figure=fig)
    
    # Define subplots
    axs = [fig.add_subplot(gs[0, i]) for i in range(0, len(plot_data))]
    # Iterate over axes
    for i, ax in enumerate(fig.axes):
        # Get name and data for dataset corresponding to iterand axis
        ax_name, ax_data = list(plot_data.items())[i]
        
        # Choose normalization corresponding to axis
        norm = normzn_abs if i != (len(plot_data)-1) else normzn_diff
        # Choose normalization corresponding to axis
        cmap = cmap_abs if i != (len(plot_data)-1) else cmap_diff
        
        # Plot dataset corresponding to axis
        im = ax.pcolormesh(ax_data, norm=norm, cmap=cmap)
        # Hide y-axis tick labels
        if i > 0:
            ax.set_yticklabels([])
            
        # Define colorbar
        cax = inset_axes(ax, width='100%', height='3%', 
                         loc='lower left', bbox_to_anchor=(0, 1.025, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0)
        
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), 
                                cax=cax, orientation='horizontal', format=lambda x, _: f'{x:.0f}')
        cax.xaxis.set_ticks_position('top')
        colorbar.ax.tick_params(labelsize=9)
        
        # Axis metadata
        if i == (len(plot_data)-1):
            if comparison_type == 'std_anom':
                ax.set_title('{0} [$\sigma$]'.format(ax_name), y=1.2, fontsize=10)
            elif comparison_type == 'rel_diff':
                ax.set_title('{0} [%]'.format(ax_name), y=1.2, fontsize=10)
            else:
                ax.set_title('{0} [{1}]'.format(ax_name, units), y=1.2, fontsize=10)
        else:
            ax.set_title('{0} [{1}]'.format(ax_name, units), y=1.2, fontsize=10)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        fig.suptitle(long_name, y=1.125)
        
def radial_tangential_velocity(data, x, y, R):
    '''
    Calculate radial and tangential velocity components from zonal and meridional velocities.
    '''

    u, v = data['ucomp'], data['vcomp']

    if np.nanmin(data.grid_yt) < 0:
        u, v = -u, -v

    # If sizes are different, reconcile by trimming the edges
    if u.shape != x.shape:
        min_dim = np.nanmin([np.nanmin([u.shape[0], x.shape[0]]), np.nanmin([u.shape[1], x.shape[1]])])
        u = u[0:min_dim, 0:min_dim]
        x = x[0:min_dim, 0:min_dim]
    if v.shape != y.shape:
        min_dim = np.nanmin([np.nanmin([v.shape[0], y.shape[0]]), np.nanmin([v.shape[1], y.shape[1]])])
        v = v[0:min_dim, 0:min_dim]
        y = y[0:min_dim, 0:min_dim]
    
    data['v_radial'] = (u*x + v*y)/R
    data['v_tangential'] = (v*x - u*y)/R
    
    data['v_radial'].attrs = {'long_name': 'radial velocity', 'units': 'm s$^{-1}$'}
    data['v_tangential'].attrs = {'v_tangential': 'tangential velocity', 'units': 'm s$^{-1}$'}
    
    return data

def add_radius(storm, test_num=None, debug=False, benchmarking=False):
    
    start = time.time()
        
    storm = storm.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    storm = radius_estimate(storm, box=10, overlay_check=False, benchmarking=False, debug=debug)
    
    if benchmarking:
        print('\t \t Dropping method: {0:.3f} s'.format(time.time() - lap))
    
    return storm

def radius_estimate(storm, box, overlay_check=True, benchmarking=False, debug=False):
    
    '''
    Algorithm to estimate the radius of a TC based on filters set below.
    Returns a radius in meters.
    
    Note 1: the idea here is: 
            (1) regrid data to higher resolution by interpolation, 
            (2) smooth data and use gradient-based filtering to prevent dual-vortex identification,
            (3) identify thresholds for value-based filtering,
            (4) use gradient- and value-based filtering to identify data points that match criteria,
            (5) clip data and estimate radius from filtered data
    '''
    
    if benchmarking:
        start = time.time()
    
    # Ensure that only one timestamp is in the Dataset
    try:
        storm = storm.isel(time=0)
    except:
        storm = storm

    # Derives horizontal wind speed (proxy for azimuthal wind)
    if 'U' not in storm.data_vars.keys():
        # Handle data from 'atmos_4xdaily' files
        if storm.attrs['filename'].split('.')[0] == 'atmos_4xdaily':
            storm['U'] = np.sqrt(storm['u_ref']**2 + storm['v_ref']**2)
        # Handle data from 'atmos_daily' files
        elif storm.attrs['filename'].split('.')[0] == 'atmos_daily':
            # Get wind data from the lowest vertical level
            storm['U'] = np.sqrt(storm['ucomp'].isel(pfull=-1)**2 + storm['vcomp'].isel(pfull=-1)**2)
        
    if benchmarking:
        print('\t \t Dataset checks: {0:.3f} s'.format(time.time() - start))
        lap = time.time()

    ''' Perform linear interpolation to allow for better gradient estimation for future filtering. '''
    # Pull numerical data from parameters relevant to radius estimation
    params = ['grid_xt', 'grid_yt', 'vort850', 'tm', 'slp', 'U']
    # Define the interpolation resolution (in degrees) and spatial extent of clipping
    resolution, extent = 0.5, 15
    # Define dictionary to store data in
    params = {param: [] for param in params}
    # Iterate through parameters and interpolate
    for param in params.keys():
        if param == 'grid_xt':
            params[param] =  np.arange(storm['center_lon'] - extent, storm['center_lon'] + extent, resolution)
        elif param == 'grid_yt':
            params[param] =  np.arange(storm['center_lat'] - extent, storm['center_lat'] + extent, resolution)
        else:
            params[param] = storm[param].interp(grid_xt=np.arange(storm['center_lon'] - extent, storm['center_lon'] + extent, resolution), 
                                                grid_yt=np.arange(storm['center_lat'] - extent, storm['center_lat'] + extent, resolution),
                                                method='nearest').values

    if debug:
        # Print center and storm bounds to ensure center is within bounds
        print('x: ', params['grid_xt'].min() < storm['center_lon'] < params['grid_xt'].max())
        print('y: ', params['grid_yt'].min() < storm['center_lat'] < params['grid_yt'].max())
    
    if benchmarking:
        print('\t \t Interpolation elapsed time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
        
    ''' 
    Data smoothing and gradient filtering algorithm. 
    The idea here is to use gradients for a chosen field to isolate TC extent and prevent dual-vortex pickup for a given storm, 
    which distorts radius calculation. 
    '''
    # 1 hPa/deg pressure gradient (attempt at similarity to Harris)
    diff_var, diff_val = 'slp', 1*resolution
    # Use a 1-sigma Gaussian smoothing filter
    smoothed = scipy.ndimage.gaussian_filter(np.abs(np.diff(np.diff(params[diff_var], axis=0), axis=1)), sigma=1)
    # Apply the filter and resize such that filter boolean array shape matches the data array shape
    diff_filter = smoothed > diff_val
    diff_filter = np.hstack((diff_filter, np.full((diff_filter.shape[0], 1), False)))
    diff_filter = np.vstack((diff_filter, np.full((1, diff_filter.shape[1]), False)))
    
    if benchmarking:
        print('\t \t Data smoothing and filter definition: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    ''' Apply thresholds to absolute values of identified parameters. '''
    # Define number of standard deviations to analyze
    sigma = 1
    # Exact magnitude thresholds
    filters = {'vort850': np.abs(params['vort850']) > 1e-4,
               'tm': params['tm'] > (np.nanmean(params['tm']) + sigma*np.nanstd(params['tm'].std())),
               'slp': params['slp'] < 1005,
               'U': params['U'] >= 5}

    
    ''' Perform the filtering and associated array clipping. '''
    # Define the conditional based on the threshold and gradient filters
    conditional = (filters['slp'] & filters['U'])
    # Define variable for filtering on
    filter_var = 'slp'
    # Perform filtering based on chosen variable
    filtered = np.where(conditional, params[filter_var], np.nan)
    
    # Crop all-nan rows in the zonal and meridional (x- and y-) array axes
    crop_x, crop_x_idx = filtered[~np.all(np.isnan(filtered), axis=1), :], ~np.all(np.isnan(filtered), axis=1)
    crop_y, crop_y_idx = crop_x[:, ~np.all(np.isnan(crop_x), axis=0)], ~np.all(np.isnan(crop_x), axis=0)
    # Output masked array for visualization of algorithm output
    arr = np.ma.masked_values(filtered, np.nan)
    
    if benchmarking:
        print('\t \t Filtering: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    ''' Derive radius, if the filtered array is not empty. '''
    # If filtering results in populated output array, get a radius
    if crop_y.shape != (0, 0):
        # If there's a mismatch in grid sizes, crop the larger one
        if params['grid_xt'].shape != params['grid_yt'].shape:
            if params['grid_xt'].shape[0] > params['grid_yt'].shape[0]:
                cut = params['grid_yt'].shape[0]
                # Perform the slicing
                params['grid_xt'] = params['grid_xt'][:cut]
                crop_y_idx = crop_y_idx[:cut]
            else:
                cut = params['grid_xt'].shape[0]
                # Perform the slicing
                params['grid_yt'] = params['grid_yt'][:cut]
                crop_x_idx = crop_x_idx[:cut]
                
        # Get the longitude and latitude extrema corresponding to the filtered array
        lons = np.min(params['grid_xt'][crop_x_idx]), np.max(params['grid_xt'][crop_x_idx])
        lats = np.min(params['grid_yt'][crop_y_idx]), np.max(params['grid_yt'][crop_y_idx])
        # Get the coordinate extrema for radius derivation
        coords = [lons[0], lats[0]], [lons[1], lats[1]]
        # Derive radius from coordinate pairs (divide by 2 and divide by 1000 to go from diameter to radius and m to km)
        radius = coords_to_dist(coords[0], coords[1])/2000
        # Add radius to the storm Dataset
        storm['radius'] = radius
        
        if benchmarking:
            print('\t \t Radius derivation: {0:.3f} s'.format(time.time() - lap))
            lap = time.time()
        
        # Overlay the storm size algorithm output on maps of the storms 
        if overlay_check:
            # Define the plot basics
            fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
            ax.coastlines()
            ax.set_extent([storm['center_lon']-extent, storm['center_lon']+extent, storm['center_lat']-extent, storm['center_lat']+extent])
            # Plot data
            im = ax.contourf(params['grid_xt'][:-1], params['grid_yt'][:-1], smoothed, levels=16)
            ax.pcolormesh(params['grid_xt'], params['grid_yt'], arr[:-1, :-1], 
                          zorder=9, cmap='Reds', transform=ccrs.PlateCarree())
            # Plot metadata
            ax.set_title('radius: {0:.2f} km'.format(radius))
            fig.colorbar(im)
            
        return storm
    # Else, return nan
    else: 
        print('Unable to get radius for {0}'.format(storm['storm_id'].values))
        if debug:
            print('------------------------------------------------')
            print(storm[['min_slp', 'max_wind']])
            print('------------------------------------------------')
        return np.nan
    
def tc_normalization(storm, param='U', test_num=None):
    
    start = time.time()
    
    num_pts = 40 # chosen because it's twice the resolution of incoming plots
    
    # Ensure that only one timestamp is in the Dataset
    try:
        storm = storm.isel(time=0)
    except:
        storm = storm
    # Get grid data and field of interest
    x, y, v = storm.grid_xt.values, storm.grid_yt.values, storm[param]

    # Only allow nonzero radii
    if (np.isfinite(storm['radius'].values)) and (storm['radius'].values > 0):
        # Normalize grid values to 1
        x = (x - x.min())*storm['radius'].values
        x = x/x.max()
        y = (y - y.min())*storm['radius'].values
        y = y/y.max()
        
        # If there's a mismatch in grid sizes, crop the larger one
        if len(x) != len(y):
            if len(x) > len(y):
                cut = len(y)
                # Perform the slicing
                x = x[:cut]
                v = v[:cut, :]
            else:
                cut = len(x)
                # Perform the slicing
                y = y[:cut]
                v = v[:, :cut]
                
        if v.shape[0] != len(y):
            v = v[:len(y), :]
        
        if v.shape[1] != len(x):
            v = v[:, :len(x)]
        
        # Create nominal Dataset
        ds = xr.DataArray(data=v, dims=['x', 'y'], coords={'x': x, 'y': y})
        # Regrid
        ds = ds.interp(x=np.linspace(0, 1, num_pts), y=np.linspace(0, 1, num_pts))
        # Add storm ID
        ds['storm_id'] = storm['storm_id']
    
        return ds
    else:
        return None
    
def tc_azimuthal_mean(data, param, debug=False, benchmarking=False):

    '''
    Method to generate azimuthal means for identified TCs.
    '''

    start = time.time()

    # Define dataset time values
    times = []
    # Container list to hold domain widths for each identified storm
    N_arr = []
    # First-pass to obtain domain dimensions across all storms
    for t, timestamp in enumerate(data.time.values):
        # Grab storm-specific data
        temp = data.isel(time=t).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        if debug:
            print('Domain dimensions: {0}, {1}'.format(len(temp.grid_xt), len(temp.grid_yt)))
        times.append(timestamp)
        # Determine minimum domain dimension (in other words, is grid_xt or grid_yt the constraint)
        N = min(len(temp.grid_xt), len(temp.grid_yt))
        # Add to container list
        N_arr.append(N)
        
    if debug:
        print('Domain constraining size: {0}'.format(min(N_arr)))

    del temp
        
    # Identify minimum domain dimension across all storms
    N = min(N_arr)
    # Create container array width dimensions (time, domain dimension, vertical level)
    container = np.full(shape=(len(times), N, len(data.pfull.values)), fill_value=np.nan)
    # Create empty distance value
    dist = np.nan
    # Create empty distance list to support normalization
    distances = []

    if benchmarking:
        print('\t Domain size retrieval time: {0:.3f} s'.format(time.time() - start))

    # Catalog all storm_ids for output wind field analysis - used for troubleshooting on composite means
    storm_times = []
    
    # Iterate over each storm
    for t, timestamp in enumerate(times):
        # Pull data for the selected parameter
        if param not in ['v_radial', 'v_tangential']:
            temp = data.sel(time=timestamp)[param].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        else:
            # Use omega as a "dynamical filler" to allow for domain computation first, 
            # then feed that domain data for radial and tangential velocity calculation
            temp = data.isel(time=t)[['ucomp', 'vcomp']].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        storm_times.append(data.sel(time=timestamp)['time'].values)
        # Calculate domain half-width, assuming quasi-axisymmetry (in other words, center to edge of domain)
        domain_hw = N // 2
        # Generate basis vectors for the domain
        basis_x, basis_y = np.linspace(0, N, N), np.linspace(0, N, N)
        # Generate domain coordinate system
        x, y = np.meshgrid(basis_x, basis_y)
        # Adjust to center coordinate system to domain center (which hopefully aligns with storm center)
        x, y = x - domain_hw, y - domain_hw
        # Calculate radii for each meshgrid grid cell
        R = np.sqrt(x**2 + y**2)

        # Handle domain trimming based on whether N is even or odd
        index_adjust = N % 2
        if debug:
            print(timestamp, N, x.shape, y.shape)
        
        # Calculate the azimuthal mean over each vertical level
        for level in range(0, len(data.pfull)):
            
            if benchmarking:
                print('\t pfull retrieval time: {0:.3f} s'.format(time.time() - start)) 
                
            if debug:
                print('Vertical level: {0:.2f} hPa'.format(level))
        
            # Clip the domain to the constraining dimensions
            if param not in ['v_radial', 'v_tangential']: 
                image = temp.isel(pfull=level)
                image_h, image_w = image.shape[0], image.shape[1]
                image = image[(image_h//2 - domain_hw):(image_h//2 + domain_hw + index_adjust), 
                              (image_w//2 - domain_hw):(image_w//2 + domain_hw + index_adjust)]
            else:
                image = temp[['ucomp', 'vcomp']].isel(pfull=level)
                image_h, image_w = image['ucomp'].shape[0], image['ucomp'].shape[1]
                image = image.isel(grid_xt=range(image_w//2 - domain_hw, image_w//2 + domain_hw + index_adjust),
                                   grid_yt=range(image_h//2 - domain_hw, image_h//2 + domain_hw + index_adjust))
                if debug:
                    print('\t', image_h, image_w, range(image_w//2 - domain_hw, image_w//2 + domain_hw + 1))
                    print('ucomp size: ({0}); x-basis vector size: ({1})'.format(image['ucomp'].shape, x.shape))
                    
                # Calculate radial and tangential velocities
                image = radial_tangential_velocity(image, x, y, R)[param]

            if image_w < 10 or image_h < 10:
                print(data.sel(time=timestamp)['time'].values)
        
            if debug:
                print('Original domain dimensions: {0}, {1}'.format(image_w, image_h))
                print('Clipped domain dimensions: {0}, {1}'.format(image.shape[1], image.shape[0]))

            # Calculate coordinates for relative domain points for distance calculation
            center_index_x, center_index_y = image.shape[1]//2, image.shape[0]//2
            center_coord = (image.grid_xt[center_index_x].values, image.grid_yt[center_index_y].values)
            outer_coord = (image.grid_xt[-1].values, image.grid_yt[center_index_y].values)
            # Compute center-outer distance in km
            dist = coords_to_dist(center_coord, outer_coord) / 1000
            
            if debug:
                print('Center coordinate: {0}, outer coordinate: {1}'.format(center_coord, outer_coord))
                print('Center-outer distance: {0:.2f} km'.format(dist))
                
                # fig, ax = plt.subplots()
                # image.plot(ax=ax)
                # ax.scatter(center_coord[0], center_coord[1], c='g', s=10)
                # plt.gca()

            # Get azimuthal mean
            f = lambda r : image.values[(R >= r-.5) & (R < r+.5)].mean()
            r = np.linspace(1, N, N)
            mean = np.vectorize(f)(r)
            
            if debug:
                print(mean)

            # Assign mean to the container array
            container[t, :, level] = mean

        # Normalize horizontal distance by derived radius
        dist = dist/data.sel(time=timestamp)['radius']
        distances.append(dist.values.item())

    if benchmarking:
        print('\t Composite retrieval time: {0:.3f} s'.format(time.time() - start))

    horiz_ax = np.linspace(0, 20, N)
    for i in range(0, container.shape[0]):
        for level in range(0, container[i, :, :].shape[1]):
            container[i, :, level] = np.interp(horiz_ax, range(0, len(container[i, :, level])), container[i, :, level])
    
    # Generate xArray Dataset for data
    # output = xr.Dataset(data_vars={param: (['level', 'distance'], np.nanmean(container, axis=0).T)},
    #                     coords={'distance': (['distance'], np.linspace(0, dist, N)),
    #                             'level': (['level'], data.pfull.values)})
    output = xr.Dataset(data_vars={param: (['level', 'distance'], np.nanmean(container, axis=0).T)},
                        coords={'distance': (['distance'], horiz_ax),
                                'level': (['level'], data.pfull.values)})
    output = output.dropna(dim='distance')

    if param not in ['v_radial', 'v_tangential']:
        output[param].attrs = {'long_name': data[param].attrs['long_name'],
                               'units': data[param].attrs['units']}

    storm_times = set([t.item() for t in storm_times])
    
    return output, storm_times

def azimuthal_composite(data, params=None, intensity_bin=None, visualization=False, individual_plots=False, debug=False):

    # Append radii
    planar_profiles = [add_radius(data[intensity_bin]['planar'][i]) 
                       for i in range(0, len(data[intensity_bin]['planar']))]
    planar_profiles = [planar_profiles[i] for i in range(0, len(planar_profiles)) 
                       if 'Dataset' in str(type(planar_profiles[i]))]
    planar_profiles = xr.concat(planar_profiles, dim='time').drop_duplicates(dim='time')
    
    if intensity_bin:
        vertical_profiles = xr.concat(data[intensity_bin]['vertical'], dim='time').drop_duplicates(dim='time')
    else:
        vertical_profiles = data['vertical'].drop_duplicates(dim='time')
        
    timestamps = []
    
    radii = []
    for storm_id in planar_profiles['storm_id'].values:
        storm_timestamp = xr.where(planar_profiles['storm_id'] == storm_id, planar_profiles['time'], np.nan).dropna(dim='time')
        temp = vertical_profiles.sel(time=storm_timestamp)
        temp['radius'] = planar_profiles.sel(time=storm_timestamp)['radius']
        radii.append(temp)
        timestamps.append(storm_timestamp.item())

    vertical_profiles = xr.concat(radii, dim='time')
    vertical_profiles = vertical_profiles.sel(time=timestamps)

    if not params:
        params = [('v_tangential', 'v_radial'),
                  ('rh', 'omega'),
                  ('temp_anom', 'omega')]
    
    if visualization:
        fig, axes = plt.subplots(figsize=(4*len(params), 4), ncols=len(params))
    
    for i, param in enumerate(params):
        fill_param, contour_param = param
        fill, storm_times = tc_azimuthal_mean(vertical_profiles, fill_param, debug=debug, benchmarking=False)

        if visualization:
            if contour_param:
                contour, _ = tc_azimuthal_mean(vertical_profiles, contour_param)
            
            levels = 16
            norm, cmap = norm_cmap(fill, fill_param, mode='raw', bounds=True, n_bounds=levels+1, difference=False)
            print(param, cmap)
            im_fill = axes[i].contourf(fill.distance, fill.level, fill[fill_param], 
                                       norm=norm, cmap=cmap, levels=levels)
            
            if contour_param:
                if contour_param == 'omega':
                    level_step = 0.25
                    contour_levels = np.arange(-1, 1+level_step, level_step)
                    contour_levels = contour_levels[contour_levels != 0]
                else:
                    extremum = round(abs(np.nanmax(contour[contour_param])))
                    contour_levels = np.linspace(-extremum, extremum, 13)
                    contour_levels = contour_levels[contour_levels != 0]
                im_contour = axes[i].contour(contour.distance, contour.level, contour[contour_param], norm=norm, colors=['k'], 
                                        alpha=0.5, linewidths=1, levels=contour_levels)
                axes[i].clabel(im_contour, im_contour.levels[::2], inline=True, fontsize=9)
            
            
            axes[i].set_ylim(axes[i].get_ylim()[::-1])
            axes[i].set_title('{0}'.format(fill_param))
            axes[i].set_xlabel('distance from TC center [km]', labelpad=10)
            if i == 0:
                axes[i].set_ylabel('pressure [hPa]', labelpad=10)
            
            if fill_param in ['v_radial', 'v_tangential']:
                long_names = {'v_radial': 'radial velocity', 'v_tangential': 'tangential velocity'}
                long_name = long_names[fill_param]
                units = 'm s$^{-1}$'
            else:
                long_name, units = vertical_profiles[fill_param].attrs['long_name'], vertical_profiles[fill_param].attrs['units']
        
            cax = axes[i].inset_axes([0, -0.25, 1, 0.03])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
            colorbar.set_label('{0} [{1}]'.format(long_name, units), labelpad=10)
            [l.set_visible(False) for i, l in enumerate(cax.get_xticklabels()) if i % 2 != 0]

    if visualization:
        fig.suptitle('Azimuthal composites for bin {0}, v = {1}'.format(intensity_bin, data[intensity_bin]['max_wind']))
        plt.gca()

    if individual_plots:
        individual_plot_param = 'U'
        for timestamp in storm_times:
            fig, ax = plt.subplots(figsize=(3, 3))
            temp = xr.concat(data[intensity_bin]['planar'], dim='time').sel(time=timestamp)
            temp[individual_plot_param].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all').plot(ax=ax)
            ax.set_title(temp['storm_id'])
            plt.gca()
            del temp

    return fill[fill_param]

def azimuthal_compositor(data, experiments, params=None, intensity_bins=None):
    
    azimuthal_composites = {}
    
    if not params:
        params = [('rh', 'omega'),
                  ('mse_anom', 'omega'),
                  ('omega', 'omega'),
                  ('temp_anom', 'omega')]

    for experiment in experiments:
        azimuthal_composites[experiment] = {}
        for intensity_bin in intensity_bins:
            azimuthal_composites[experiment][intensity_bin] = {}
            for param in params:
                azimuthal_composites[experiment][intensity_bin][param[0]] = azimuthal_composite(data[experiment], params=[param], intensity_bin=intensity_bin, 
                                                                                debug=False, individual_plots=False)
                
    return azimuthal_composites

def mean_profile_vertical(data, experiments, params, intensity_bin, ax):
    
    ''' Method to obtain storm-centered vertical profiles. '''
    
    if type(params) is not list:
        params = [params]
    if type(experiments) is not list:
        experiments = [experiments]
    
    resolution = 0.5
    degrees = 5
    extent = round(degrees/resolution/2)    
    
    plot_data = {}
    for experiment in experiments:
        plot_data[experiment] = {param: np.array([]) for param in params}
        times = len(data[experiment][intensity_bin]['vertical'])
        for time_index in range(0, times):
            storm = data[experiment][intensity_bin]['vertical'][time_index]
            storm = storm.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            # storm = add_radius(storm)
            center_x, center_y = len(storm['grid_xt'].values)//2, len(storm['grid_yt'].values)//2
            for param in params:
                if param == 'rh':
                    values = storm[param].isel(grid_xt=slice(center_x-extent, center_x+extent), 
                                               grid_yt=slice(center_y-extent, center_y+extent)).mean(dim='grid_xt').mean(dim='grid_yt').values.T
                else:
                    values = storm[param].isel(grid_xt=slice(center_x-extent, center_x+extent), 
                                               grid_yt=slice(center_y-extent, center_y+extent)).mean(dim='grid_xt').mean(dim='grid_yt').values
                
                if len(plot_data[experiment][param]) == 0:
                    plot_data[experiment][param] = values
                else:                        
                    plot_data[experiment][param] = np.append(plot_data[experiment][param], values, axis=0)
    
        for param in params:
            plot_data[experiment][param] = {'mean': np.nanmean(plot_data[experiment][param], axis=0),
                                            'std': np.nanstd(plot_data[experiment][param], axis=0)}
    
    # fig, ax = plt.subplots(figsize=(3, 4))
    
    colors = ['b', 'r']
    linestyles = ['-', ':', '--']
    
    ax.axvline(0, c='k', alpha=0.5)
    
    for i, experiment in enumerate(experiments):
        for j, param in enumerate(params):
            ax.fill_betweenx(data[experiment][intensity_bin]['vertical'][0].pfull,
                             plot_data[experiment][param]['mean'] - plot_data[experiment][param]['std'], 
                             plot_data[experiment][param]['mean'] + plot_data[experiment][param]['std'],
                             fc=colors[i], alpha=0.1)
            ax.plot(plot_data[experiment][param]['mean'], data[experiment][intensity_bin]['vertical'][0].pfull,
                    lw=3, c=colors[i], ls=linestyles[j], label=experiment)

    xlabel = '{0} [{1}]'.format(data[experiment][intensity_bin]['vertical'][0][param].attrs['long_name'],
                              data[experiment][intensity_bin]['vertical'][0][param].attrs['units'])
    ax.set_xlabel(xlabel)
    
    ytick_step = 100

    pressure_bounds = [100, 1000] if params == ['omega'] else [200, 1000]
    yticks = np.arange(min(pressure_bounds), max(pressure_bounds)+ytick_step, ytick_step)

    # Note: logarithmic scaling omitted to emphasize boundary layer processes
    # ax.set_yscale('log')
    
    ax.set_title('bin {0}, v = {1}'.format(intensity_bin, data[experiment][intensity_bin]['max_wind']), y=1.05)
    ax.set_yticks(yticks, ['{0:0d}'.format(ytick) for ytick in yticks])
    ax.set_ylim(pressure_bounds)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.legend(frameon=False)
    
def planar_composite(data, model_name, experiments=['control', 'swishe'], param='U', vertical_level=500, intensity_bin='b5'):

    fig = plt.figure(figsize=(9, 4))
    gs = matplotlib.gridspec.GridSpec(nrows=1, ncols=3, wspace=0.25)
    
    # Correction factor for precip from kg m^-2 s^-1 to mm d^-1
    factor = 86400 if param in ['precip', 'evap'] else 1
    # Rig the script to get latent heat flux from evaporation
    lhflx = False
    if param == 'lhflx':
        lhflx = True
        param = 'evap'
        factor = 2.5e6
    
    planar_params = data[experiments[0]][intensity_bin]['planar'][1].copy().data_vars
    vertical_params = data[experiments[0]][intensity_bin]['vertical'][1].copy().data_vars
    
    if param in planar_params:
        dataset = 'planar'
    elif param in vertical_params:
        dataset = 'vertical'
    
    plot_data = {}
    
    for experiment in experiments:
        plot_data[experiment] = {param: np.array([])}
        times = len(data[experiment][intensity_bin][dataset])
        for time_index in range(0, times):
            storm = data[experiment][intensity_bin][dataset][time_index]
            storm = storm.sel(pfull=vertical_level, method='nearest') if dataset == 'vertical' else storm
    
            # Needs work to map vertical profiles onto planes for compositing
            if dataset == 'vertical':
                temp = data[experiment][intensity_bin]['planar'][time_index].copy().dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
                storm_slice = storm[param].interp(grid_xt=temp.grid_xt.values, grid_yt=temp.grid_yt.values)
                # print('Planar: ', temp.grid_xt.min().item(), temp.grid_xt.max().item(), temp.grid_yt.min().item(), temp.grid_yt.max().item())
                # print('Storm: ', storm.grid_xt.min().item(), storm.grid_xt.max().item(), storm.grid_yt.min().item(), storm.grid_yt.max().item())
                # print('---')
                temp[param] = storm_slice
                storm = temp
                del temp
                
            temp = storm.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            temp = add_radius(temp, debug=False)
            try:
                temp = tc_normalization(temp, param=param)
            except:
                continue
            
            if len(plot_data[experiment][param]) == 0:
                if temp is not None:
                    plot_data[experiment][param] = [temp]
            else:       
                if temp is not None:
                    plot_data[experiment][param] = np.append(plot_data[experiment][param], 
                                                             [temp], axis=0)
        
        plot_data[experiment][param] = np.nanmean(plot_data[experiment][param], axis=0)   
        plot_data[experiment][param] = (plot_data[experiment][param])[:, ~np.isnan(plot_data[experiment][param]).all(axis=0)]
        plot_data[experiment][param] = (plot_data[experiment][param])[~np.isnan(plot_data[experiment][param]).all(axis=1), :]
        plot_data[experiment]['shape'] = plot_data[experiment][param].shape
    
    ''' Domain reconciliation algorithm. '''
    # 1. Get shapes of each composite domain. (done above)
    # 2. Get minimum dimensions in x and y.
    # 3. Get the minimum of these dimensions.
    # 4. Trim each domain to the minimum obtained in Step 3.
    
    # Step 2. Get minimum dimensions in x and y.
    min_x, min_y = np.inf, np.inf
    for experiment in experiments:
        if plot_data[experiment]['shape'][0] < min_x:
            min_x = plot_data[experiment]['shape'][0]
        if plot_data[experiment]['shape'][1] < min_y:
            min_y = plot_data[experiment]['shape'][1]
    # Step 3. Get the minimum of these dimensions.
    min_dim = min(min_x, min_y) if (np.isfinite(min_x) and np.isfinite(min_y)) else np.nan
    # Step 4. Trim each domain to the minimum obtained in Step 3.
    for experiment in experiments:
        center_x, center_y = plot_data[experiment]['shape'][0] // 2, plot_data[experiment]['shape'][1] // 2
        plot_data[experiment][param] = plot_data[experiment][param][(center_x-min_dim//2):(center_x+min_dim//2), (center_y-min_dim//2):(center_y+min_dim//2)]
        plot_data[experiment][param] = plot_data[experiment][param][3:-3, 3:-3]
    
    param = 'lhflx' if lhflx else param
    if param == 'lhflx':
        for experiment in experiments:
            plot_data[experiment]['lhflx'] = plot_data[experiment]['evap']
        
    norm, cmap = norm_cmap([plot_data[experiment] for experiment in experiments], 
                           param, mode='raw', n_bounds=25)
    
    for i, experiment in enumerate(experiments):
        ax = fig.add_subplot(gs[0, i])
        im = ax.pcolormesh(factor*plot_data[experiment][param], cmap=cmap, norm=norm)
        ax.set_aspect('equal')
    
        cax = ax.inset_axes([0, 1.05, 1, 0.05])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                cax=cax, orientation='horizontal')
        cax.xaxis.set_ticks_position('top')
        cax.tick_params(labelsize=9)
    
        ax.set_title('{0} (N={1})'.format(experiment, len(data[experiment][intensity_bin][dataset])), fontsize=10, y=1.25)
    
        if (abs(norm.vmin) < 0.1) and (abs(norm.vmax) < 0.1):
            cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
        elif (abs(norm.vmin) >= 1e3) and (abs(norm.vmax) >= 1e3):
            cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
        else:
            cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))
        [c.set_visible(False) for i, c in enumerate(cax.get_xticklabels()) if i % 3 != 0 ]
        cax.tick_params(labelsize=8)
    
    if param == 'U':
        long_name, units = 'Horizontal wind speed', 'm s^-1'
    elif param == 'lhflx':
        long_name, units = 'latent heat flux', 'W m$^{-2}$'
    elif param in ['precip', 'evap']:
        long_name, units = storm[param].attrs['long_name'], 'mm d^-1'
    else:
        long_name, units = storm[param].attrs['long_name'], storm[param].attrs['units']
    
    difference_mode = 'simple'
    if difference_mode == 'relative':
        difference = (plot_data['swishe'][param] - plot_data['control'][param])/(plot_data['swishe'][param] + plot_data['control'][param])
        difference = np.where(np.abs(difference) < 100, difference, np.nan)
    else:
        difference = plot_data['swishe'][param] - plot_data['control'][param]
    
    if np.sign(np.nanmin(difference)) == np.sign(np.nanmax(difference)):
        norm, cmap = norm_cmap(difference, param, mode='raw', n_bounds=17, difference=False)
    else:
        norm, cmap = norm_cmap(difference, param, mode='raw', n_bounds=17, difference=True)
    
    ax = fig.add_subplot(gs[0, -1])
    # cmap = 'Blues_r'
    im = ax.pcolormesh(factor*difference, cmap=cmap, norm=norm)
    ax.set_aspect('equal')
    
    ax.set_title('difference', fontsize=10, y=1.25)
    
    cax = ax.inset_axes([0, 1.05, 1, 0.05])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                            cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    
    if (abs(norm.vmin) < 0.1) and (abs(norm.vmax) < 0.1):
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
    elif (abs(norm.vmin) >= 1e3) and (abs(norm.vmax) >= 1e3):
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
    else:
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))
    [c.set_visible(False) for i, c in enumerate(cax.get_xticklabels()) if i % 3 != 0 ]
    cax.tick_params(labelsize=8)
    
    # Figure title definition
    if dataset == 'vertical':
        figure_title = '{0} [{1}] at {2} hPa\n{3}, bin {4}'.format(long_name.capitalize(), units, vertical_level, 
                                                                model_name, intensity_bin)  
    else:
        figure_title = '{0} [{1}]\n{2}, bin {3}'.format(long_name.capitalize(), units, 
                                                               model_name, intensity_bin) 
    
    fig.suptitle(figure_title, y=1.025, fontsize=12)
    
def mse(data, experiments, intensity_bin):

    # Reference: Emanuel (1994) 
    c_p = 1005.7
    L_v = 2.5e6
    R_d = 287.04
    eps = 0.622
    g = 9.81
    # Reference: Huang (2018), doi:10.1175/JAMC-D-17-0334.1
    e_s_huang_water = lambda t: np.exp(34.494 - 4924.99/((t-273.15)+237.1))/((t-273.15) + 105)**(1.57)
    e_s_huang_ice = lambda t: np.exp(43.494 - 6545.8/((t-273.15)+278))/((t-273.15) + 868)**2

    for experiment in experiments:
        for i in range(0, len(data[experiment][intensity_bin]['vertical'])):
    
            iterdata = data[experiment][intensity_bin]['vertical'][i]
            # Modify e_s approximation based on temperature
            e_s = xr.where(iterdata['temp'] >= 273.15, e_s_huang_water(iterdata['temp']), 
                           e_s_huang_ice(iterdata['temp']))
            e = iterdata['rh']*e_s/100
            iterdata['e'] = e/100
            iterdata['e'].attrs = {'long_name': 'vapor pressure', 'units': 'hPa'}
        
            iterdata['rho'] = 100*(iterdata.pfull - iterdata['e'])/(R_d*iterdata['temp']) + 100*iterdata['e']/(461.5*iterdata['temp'])
            iterdata['mse'] = c_p*iterdata['temp'] + L_v*iterdata['sphum'] + (100*iterdata.pfull)/iterdata['rho']
            iterdata['mse'].attrs = {'long_name': 'moist static energy', 'units': 'J kg^-1'}

            # Beware of indexing here - indexing errors may arise from implictly-defined indices
            # for time, pfull, grid_yt, and grid_xt.
            cmse = np.full(iterdata['mse'].values.shape, np.nan)
            for j in range(1, len(iterdata.pfull.values)):
                cmse[:, j, :, :] = iterdata['mse'].isel(pfull=j).values * 100*(iterdata.pfull.values[j] - iterdata.pfull.values[j-1])/g
            cmse = np.where(np.nansum(cmse, axis=1) > 0, np.nansum(cmse, axis=1), np.nan)

            iterdata['cmse'] = xr.DataArray(dims=['time', 'grid_yt', 'grid_xt'], data=cmse)
            iterdata['cmse'].attrs = {'long_name': 'column-integrated moist static energy', 'units': 'J m^-2'}

            data[experiment][intensity_bin]['vertical'][i]['mse'] = iterdata['mse']
            data[experiment][intensity_bin]['vertical'][i]['cmse'] = iterdata['cmse']
            
            data[experiment][intensity_bin]['vertical'][i]['mse_anom'] = iterdata['mse'] - iterdata['mse'].mean(dim='grid_xt').mean(dim='grid_yt')
            data[experiment][intensity_bin]['vertical'][i]['mse_anom'].attrs = {'long_name': 'domainwise moist static energy anomaly', 'units': 'J kg^-1'}
            
            del iterdata

    return data

def histogram(data):
    
    """
    Method to create histogram of storms per intensity bin.
    
    Arguments:
        data {xArray}
    """    
    
    fig, ax = plt.subplots(figsize=(6, 3))

    colors = ['b', 'r']

    for i, k in enumerate(data.keys()):
        for j, intensity_bin in enumerate(data[k].keys()):
            width = 0.4
            label = k if j == 0 else '_nolegend_'
            ax.bar(j+width*i, len(data[k][intensity_bin]['planar']), fc=colors[i], width=width, label=label)

    ax.legend(frameon=False)
    ax.set_xticks(range(0, len(data['control'].keys())), data['control'].keys());

def azimuthal_plot(azimuthal_composites, param, intensity_bins, model_name, levels=16):

    ''' Visualization. '''
    fig = plt.figure(figsize=(3*len(intensity_bins), 4))
    gs = matplotlib.gridspec.GridSpec(nrows=2, ncols=len(intensity_bins), height_ratios=(1, 0.1))

    differences = {}
    for intensity_bin in intensity_bins:
        differences[intensity_bin] = azimuthal_composites['swishe'][intensity_bin][param] - azimuthal_composites['control'][intensity_bin][param]

    norm, cmap = norm_cmap([differences[i] for i in intensity_bins], 
                        param, mode='raw', bounds=True, n_bounds=levels+1, difference=True)

    for i in range(0, len(intensity_bins)):

        ax = fig.add_subplot(gs[0, i])
        im = ax.contourf(differences[intensity_bins[i]].distance, differences[intensity_bins[i]].level, differences[intensity_bins[i]], norm=norm, cmap=cmap, levels=levels)
                    
        ax.set_xlabel('normalized distance\nfrom TC center', labelpad=10)
        ax.set_title('{0}, bin {1}'.format(model_name, intensity_bins[i]))
        if i == 0:
            ax.set_ylabel('pressure [hPa]', labelpad=10)
        else:
            ax.set_yticklabels([])

        ax.set_ylim([100, 1000])
        ax.set_ylim(ax.get_ylim()[::-1])
        
        if param in ['v_radial', 'v_tangential']:
            long_names = {'v_radial': 'radial velocity', 'v_tangential': 'tangential velocity'}
            long_name = long_names[param]
            units = 'm s$^{-1}$'
        else:
            long_name, units = [data['control'][intensity_bins[i]]['vertical'][0][param].attrs['long_name'], 
                                data['control'][intensity_bins[i]]['vertical'][0][param].attrs['units']]

    colorbar_ax = fig.add_subplot(gs[1, :])
    colorbar_ax.set_axis_off()
    cax_width = 0.5
    cax = colorbar_ax.inset_axes([0.5-cax_width/2, -2, cax_width, 0.5])
            
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
    colorbar.set_label('{0} [{1}]'.format(long_name, units), labelpad=10)

    if (abs(norm.vmin) < 0.1) and (abs(norm.vmax) < 0.1):
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
    elif (abs(norm.vmin) >= 1e3) and (abs(norm.vmax) >= 1e3):
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.1e}'))
    else:
        cax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:.2f}'))

def main(model_name, storm_category, experiments, years):
    data = get_data(model_name, storm_category, experiments, years=years, nums={'control': 101, 'swishe': 101})

    for intensity_bin in ['b0', 'b1', 'b2', 'b3']:
        data = mse(data, ['control', 'swishe'], intensity_bin)
        
    return data

if __name__ == '__main__':
    override = False
    if override or 'data' not in globals().keys():
        data = main('FLOR', 'TS', ['control', 'swishe'], (2001, 2050))