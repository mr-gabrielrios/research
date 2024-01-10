import cftime, datetime
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle, random

import utilities

def planar_compositor(model, datasets, intensity_bin, field, pressure_level=None):
    
    """
    Method to generate planar composites for a list of given TCs for a field, 
    pressure level (if not a surface or integrated quantity), and an intensity bin.

    Args:
        model (str): name of model
        datasets (list): list of data dictionaries (.pkl file)
        intensity_bin (str): string
        field (str): field to be evaluated
        pressure_level (numeric): pressure level at which to evaluate the data, if applicable. Default is None.
    Returns:
        data (dictionary): dictionary with list of processed data. Dictionary is 2-tiered: (1) storm ID --> (2) field
        composite_mean (dictionary): dictionary with center-relative distances, in degrees, and composite mean data.  
    """
    
    # Define key with format {FIELD}{HEIGHT} to match GFDL diagnostics names.
    # For example, specific humidity at 700 hPa would be 'q700'.
    key = '{0}{1}'.format(field, pressure_level) if pressure_level else field
    # Grab minimum extents for the datasets for future trimming based on longitude (grid_xt) and latitude (grid_yt) extents
    min_x, min_y = None, None
    # Initialize dictionary to hold data for each storm corresponding to the intensity bin and field of choice
    data = {}
    # Iterate over datasets
    for ds in datasets:
        # Get storm ID
        storm_id = ds['track_output']['storm_id'].unique().item()
        print('\t Processing {0}...'.format(storm_id))
        # Get timestamps relevant to the intensity bin at hand
        times = ds['track_output']['time'].loc[ds['track_output']['intensity_bin'] == intensity_bin]
        # Change from Pandas to cftime formats
        times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times]
        # Determine which dictionary the field data is in
        subdict = None
        for sd in ['tc_model_output', 'tc_vertical_output']:
            if field in ds[sd].data_vars:
                subdict = sd
        # Print error message
        if subdict is None:
            print('Field not found in the data, please try again...')
            return None, None
        # Get the matching timestamps between the track output and the corresponding dictionary with field data
        index_timestamps = [t for t in times if t in ds[subdict][field].time.values]
        # and field + vertical level from the vertical data
        ds = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[subdict][field].sel(time=index_timestamps)
        # If timestamps match, proceed. Else, continue.
        if len(ds.time.values) > 0:
            # Populate dictionary with data. Limit to one entry per storm by choosing the first timestamp.
            data[storm_id] = {key: ds.isel(time=0).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
            # Check to see if the extents are less than the minima. If so, make the new minima
            if (min_x == None) or (len(data[storm_id][key].grid_xt) < min_x):
                min_x = len(data[storm_id][key].grid_xt)
            if (min_y == None) or (len(data[storm_id][key].grid_yt) < min_y):
                min_y = len(data[storm_id][key].grid_yt)
        else:
            continue
    
    # Define buffer to prevent edge effects of domain trimming. 
    # Will be applied to each edge, such that the final domain will be 2*buffer less after trimming.
    buffer = 1
    
    # Iterate over each storm and trim
    for storm_id in data.keys():
        # Make minimum domain lengths even by reducing frame size, if not even
        min_x = min_x - 1 if min_x % 2 == 1 else min_x
        min_y = min_y - 1 if min_y % 2 == 1 else min_y
        # Get midpoints
        center_x, center_y = len(data[storm_id][key].grid_xt) // 2, len(data[storm_id][key].grid_yt) // 2
        # Define slices
        slice_x = slice(center_x - min_x//2 + buffer, center_x + min_x//2 - buffer)
        slice_y = slice(center_y - min_y//2 + buffer, center_y + min_y//2 - buffer)
        # Slice the domain
        data[storm_id][key] = data[storm_id][key].isel(grid_xt=slice_x, grid_yt=slice_y)
        # Flip the data about the x-axis if the center latitude is < 0 (for Southern Hemisphere storms).
        if data[storm_id][key].isel(grid_yt=center_y)['grid_yt'] < 0:
            data[storm_id][key] = np.flip(data[storm_id][key], axis=0)
    
    # Get composite mean over all storms by averaging along the 0th axis (storm ID axis)
    composite_mean = np.nanmean(np.stack([v[key] for v in data.values()]), axis=0)
    # Get degree spacing for the data. Assumes all spacings are equal for each dimension
    dx, dy = [[entry[key] for entry in data.values()][0].grid_xt.diff(dim='grid_xt').values[0],
              [entry[key] for entry in data.values()][0].grid_yt.diff(dim='grid_yt').values[0]]
    # Define dictionary to hold distance data (x and y) from composite_mean shape
    # The domains are define as +/- 1/2 the width (centered on 0), stepping by 1 and scaling by dx and dy
    composite_mean_x = dx*np.arange(-composite_mean.shape[1]//2, composite_mean.shape[1]//2, 1)
    composite_mean_y = dy*np.arange(-composite_mean.shape[0]//2, composite_mean.shape[0]//2, 1)
    # Populate with the dimensional data
    composite_mean = {'grid_xt': composite_mean_x, 'grid_yt': composite_mean_y, 'data': composite_mean}    
    # Note: composite_mean can be made into an xArray DataArray directly from a dictionary.
    
    return data, composite_mean

def azimuthal_compositor(model, datasets, intensity_bin, field):
    
    # Define key with format {FIELD}{HEIGHT} to match GFDL diagnostics names.
    # For example, specific humidity at 700 hPa would be 'q700'.
    key = field
    # Determine which dictionary the field is in. Used externally as a placeholder for future averaging.
    subdict = None
    # Grab minimum extents for the datasets for future trimming based on longitude (grid_xt) and latitude (grid_yt) extents
    min_x, min_y = None, None
    # Initialize dictionary to hold data for each storm corresponding to the intensity bin and field of choice
    data = {}
    # Iterate over datasets
    for ds in datasets:
        # Get storm ID
        storm_id = ds['track_output']['storm_id'].unique().item()
        print('\t Processing {0}...'.format(storm_id))
        # Get timestamps relevant to the intensity bin at hand
        times = ds['track_output']['time'].loc[ds['track_output']['intensity_bin'] == intensity_bin]
        # Change from Pandas to cftime formats
        times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times]
        # Determine which dictionary the field data is in
        for sd in ['tc_model_output', 'tc_vertical_output']:
            if field in ds[sd].data_vars:
                subdict = sd
        # Print error message
        if subdict is None:
            print('Field not found in the data, please try again...')
            return None, None
        # Get the matching timestamps between the track output and the corresponding dictionary with field data
        index_timestamps = [t for t in times if t in ds[subdict][field].time.values]
        # and field + vertical level from the vertical data
        ds = ds[subdict][field].sel(time=index_timestamps)
        if len(ds.time.values) > 0:
            # Populate dictionary with data. Limit to one entry per storm by choosing the first timestamp.
            data[storm_id] = {key: ds.isel(time=0).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
            # Check to see if the extents are less than the minima. If so, make the new minima
            if (min_x == None) or (len(data[storm_id][key].grid_xt) < min_x):
                min_x = len(data[storm_id][key].grid_xt)
            if (min_y == None) or (len(data[storm_id][key].grid_yt) < min_y):
                min_y = len(data[storm_id][key].grid_yt)
        else:
            continue
    
    # Define buffer to prevent edge effects of domain trimming. 
    # Will be applied to each edge, such that the final domain will be 2*buffer less after trimming.
    buffer = 1
    
    # Iterate over each storm and trim
    for storm_id in data.keys():
        # Make minimum domain lengths even by reducing frame size, if not even
        min_x = min_x - 1 if min_x % 2 == 1 else min_x
        min_y = min_y - 1 if min_y % 2 == 1 else min_y
        # Get midpoints
        center_x, center_y = len(data[storm_id][key].grid_xt) // 2, len(data[storm_id][key].grid_yt) // 2
        # Define slices
        slice_x = slice(center_x - min_x//2 + buffer, center_x + min_x//2 - buffer)
        slice_y = slice(center_y - min_y//2 + buffer, center_y + min_y//2 - buffer)
        # Slice the domain
        data[storm_id][key] = data[storm_id][key].isel(grid_xt=slice_x, grid_yt=slice_y)
        # Flip the data about the x-axis if the center latitude is < 0 (for Southern Hemisphere storms).
        if data[storm_id][key].isel(grid_yt=center_y)['grid_yt'] < 0:
            data[storm_id][key] = np.flip(data[storm_id][key], axis=0)
            
    # Now, all data should be trimmed to equal domain sizes.
    # Iterate over each storm and get radii
    for storm_id in data.keys():
        # Get new midpoints for trimmed data
        center_x, center_y = len(data[storm_id][key].grid_xt) // 2, len(data[storm_id][key].grid_yt) // 2
        # Create coordinates for a TC-centric coordinate system (TC_xt, TC_yt)
        data[storm_id][key] = data[storm_id][key].assign_coords({'TC_xt': data[storm_id][key]['grid_xt'] - data[storm_id][key]['grid_xt'].isel(grid_xt=center_x),
                                                                 'TC_yt': data[storm_id][key]['grid_yt'] - data[storm_id][key]['grid_yt'].isel(grid_yt=center_y)})
        # Calculate radial distance in the TC-centric coordinate system from the domain center (center_x, center_y)
        X, Y = np.meshgrid(data[storm_id][key]['TC_xt'].values, data[storm_id][key]['TC_yt'].values) # radius substep 1: create meshgrid
        data[storm_id][key]['radius'] = xr.DataArray(data=np.sqrt(X**2 + Y**2), dims=['grid_yt', 'grid_xt'])
        # Initiate radius vector beginning from TC center to domain edge (this is selected using -1)
        # Note: assume minimum domain dimension is 'grid_xt'
        radius_vector = data[storm_id][key]['radius'].isel(grid_xt=slice(center_x, -1), grid_yt=center_y)
        
        # Begin the azimuthal mean sequence here.
        # This sequence will essentially step out radially from the center by establishing rings of tolerance at half the grid resolution
        # Create a tolerance of half the grid resolution
        tolerance =  data[storm_id][key]['TC_xt'].diff(dim='grid_xt').values[0] / 2
        
        # If data has a vertical dimension, process along the vertical axis
        if 'pfull' in data[storm_id][key].dims:
            # Initialize the output NumPy array
            out = np.full(shape=(len(data[storm_id][key].pfull.values), len(radius_vector)), fill_value=np.nan)
            # Iterate over each pressure level found in the data
            for j, pressure_level in enumerate(data[storm_id][key].pfull.values):
                # Iterate over all radii in radius_vector
                for i, r in enumerate(radius_vector.values):
                    # Define the radial limits for the 'i' iterand
                    inner_radius, outer_radius = r - tolerance, r + tolerance
                    # Extract data from the radial ring at the given pressure level
                    temp = data[storm_id][key].isel(pfull=j).where((data[storm_id][key].isel(pfull=j)['radius'] >= inner_radius) & 
                                                                   (data[storm_id][key].isel(pfull=j)['radius'] < outer_radius)).values
                    # Filter out missing data and assign to the corresponding array location
                    out[j, i] = np.nanmean(temp[~np.isnan(temp)])
            # Generate the output xArray DataArray
            data[storm_id][key] = xr.DataArray(data=out, dims=['pfull', 'radius'], 
                                  coords={'pfull': (['pfull'], data[storm_id][key].pfull.values),
                                  'radius': (['radius'], radius_vector.values)})
        # Else, process on the single available level
        else:
            # Initialize the output NumPy array
            out = np.full(shape=(len(radius_vector)), fill_value=np.nan)
            # Iterate over all radii in radius_vector
            for i, r in enumerate(radius_vector.values):
                # Define the radial limits for the 'i' iterand
                inner_radius, outer_radius = r - tolerance, r + tolerance
                # Extract data from the radial ring at the given pressure level
                temp = data[storm_id][key].where((data[storm_id][key]['radius'] >= inner_radius) & 
                                                                (data[storm_id][key]['radius'] < outer_radius)).values
                # Filter out missing data and assign to the corresponding array location
                out[i] = np.nanmean(temp[~np.isnan(temp)])
            # Generate the output xArray DataArray
            data[storm_id][key] = xr.DataArray(data=out, dims=['radius'], coords={'radius': (['radius'], radius_vector.values)})
    
    # Get sample dataset for axis building in the composite array.
    # Note: this assumes all radii are equal.
    sample_dataset = [v[key] for v in data.values()][0]
    # Build xArray with dimensions corresponding to the requested field
    if subdict == 'tc_vertical_output':
        # Get composite mean over all storms by averaging along the 0th axis (storm ID axis)
        # Generate the output xArray DataArray
        composite_mean = xr.DataArray(data=np.nanmean(np.stack([v[key] for v in data.values()]), axis=0), 
                                      dims=['pfull', 'radius'], 
                                      coords={'pfull': (['pfull'],sample_dataset['pfull'].values),
                                              'radius': (['radius'], sample_dataset['radius'].values)})
    elif subdict == 'tc_model_output':
        # Get composite mean over all storms by averaging along the 0th axis (storm ID axis)
        # Generate the output xArray DataArray
        composite_mean = xr.DataArray(data=np.nanmean(np.stack([v[key] for v in data.values()]), axis=0), 
                                      dims=['radius'], coords={'radius': (['radius'], sample_dataset['radius'].values)})
    
    return data, composite_mean

if __name__ == '__main__':
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed'
    model, experiment = 'HIRAM', 'swishe'
    num_storms = 5 # enter -1 to get all
    storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] 
                 for f in os.listdir(dirname)
                 if (model in f) and (experiment in f)][:num_storms]
    print('Storms to be processed: {0}'.format(storm_ids))
    datasets = []
    for storm_id in storm_ids:
        _, dataset = utilities.access(model, experiment, storm_type='C15w', storm_id=storm_id, processed=True)
        datasets.append(dataset)
