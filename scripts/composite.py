import cftime, datetime
import multiprocessing
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
        times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times if not ((t.month == 2) & (t.day == 29))]
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
        # This try/except is supposed to catch duplicate time indices in ds[subdict]
        try:
            dataset = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[subdict][field].sel(time=index_timestamps)
        except:
            # Check for duplicate time values
            ds[subdict] = ds[subdict].drop_duplicates('time')
            # Now try reloading
            dataset = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[subdict][field].sel(time=index_timestamps)
        # If timestamps match, proceed. Else, continue.
        if len(dataset.time.values) > 0:
            print('\t \t {0} has {1} matching timestamps! Further processing them now...'.format(storm_id, len(dataset.time.values)))
            # Methodology: extract all timestamps with matching intensity for the given storm and assign them a storm 'sub-ID', which preserves the storm ID but gives it a unique identifier.
            storm_subids = [] # create storage list for all storm sub-IDs to be generated
            ''' Populate dictionary with data. '''
            # If one intensity bin is found in the dataset, limit to N strongest entries per storm by filtering the storm according to the matching timestamps.
            num_timestamps = 5 if len(dataset.time.values) > 5 else len(dataset.time.values) # limit number of timestamp indices per storm to this number
            timestamps = [utilities.time_adjust(model, t, method='cftime_to_pandas') for t in dataset.time.values if not ((t.month == 2) & (t.day == 29))] # Get Pandas-friendly timestamps for DataFrame indexing
            timestamps = ds['track_output'].loc[ds['track_output']['time'].isin(timestamps)].sort_values('max_wind', ascending=False).iloc[0:num_timestamps]['time'] # Get N strongest timestamps
            timestamps = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in timestamps] # Revert to cftime for DataArray indexing
            # Process data for each filtered timestamp
            for i, timestamp in enumerate(timestamps):   
                # Assign the storm sub-ID
                storm_subid = '{0}_{1:03d}'.format(storm_id, i)
                # Extract all nans
                data[storm_subid] = {key: dataset.sel(time=timestamp).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
                # Append to storage list
                storm_subids.append(storm_subid)
                print('\t \t \t {0} created, continuing...'.format(storm_subid))
            for storm_subid in storm_subids:
                # Check to see if the extents are less than the minima. If so, make the new minima
                if (min_x == None) or (len(data[storm_subid][key].grid_xt) < min_x):
                    min_x = len(data[storm_subid][key].grid_xt)
                if (min_y == None) or (len(data[storm_subid][key].grid_yt) < min_y):
                    min_y = len(data[storm_subid][key].grid_yt)
        else:
            # print('\t \t {0} was not processed since no timestamps were found...'.format(storm_id))
            continue
    
    # Define buffer to prevent edge effects of domain trimming. 
    # Will be applied to each edge, such that the final domain will be 2*buffer less after trimming.
    buffer = 1
    
    # Iterate over each storm and its sub-IDs and trim
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
        # if data[storm_id][key].isel(grid_yt=center_y)['grid_yt'] < 0:
        #     data[storm_id][key] = np.flip(data[storm_id][key], axis=0)
    
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
    composite_mean = xr.DataArray(data=composite_mean['data'], dims=['grid_yt', 'grid_xt'], 
                                  coords={'grid_yt': (['grid_yt'], composite_mean['grid_yt']), 
                                          'grid_xt': (['grid_xt'], composite_mean['grid_xt'])})
    
    return key, data, composite_mean

''' Begin azimuthal compositing support methods. '''
# These methods enable the radial variation of azimuthal means to be parallelized. Speedups >4x, especially with larger numbers of storms.
def azimuthal_mean(data, storm_id, field):
    """
    Core method for azimuthal compositing. This method is separate to enable faster parallelization.

    Args:
        data (dict): dictionary with key of storm ID and values of an xArray Dataset
        storm_id (str): storm ID
        field (str): field in Dataset

    Returns:
        storm_id (str): storm ID
        data (dict): dictionary with key of storm ID and values of an xArray Dataset
    """
    
    ''' Methodology.
    1. Create radius field from storm-centered coordinates.
    2. Define grid resolution to define annular thickness.
    3. Define the starting radius to initialize annulus
    4. Create annular basis vector to determine where the annular radii will be defined.
    5. Initialize output array, with dimensions of (pfull, basis_vector).
    6. Create annulus, whose area will be averaged over.
    7. Average over annular area.
    8. Expand annulus by a grid cell and repeat 6-8 until domain edge is reached.
    '''
    
    # Define shorthand for data
    dataset = data[storm_id][field]
    # Initialize list to hold figures for future plotting for debugging
    visuals = {}
    
    # 1. Create radial field relative to storm center.
    X, Y = np.meshgrid(dataset['TC_xt'], dataset['TC_yt'])
    dataset['radius'] = xr.DataArray(data=np.sqrt(X**2 + Y**2), dims=('grid_yt', 'grid_xt'))
    # 2. Define grid resolution from coarser of the two domain dimensions ('TC_xt' or 'TC_yt')
    # Note: because the spatial extent is square, the dimension with fewer elements is coarser by definition
    limiting_dimension = 'TC_xt' if len(dataset['TC_xt']) < len(dataset['TC_yt']) else 'TC_yt'
    resolution = 1.5*np.diff(dataset[limiting_dimension])[0] # assume all elements are equal
    # 3. Create the first annulus inner radius and the radial limit
    r_i = 0
    r_limit = dataset[limiting_dimension].max()
    # 4. Create a basis vector over the radius to guide the annular expansion
    rs = np.arange(r_i, r_limit + resolution, resolution)
    # 5. Create output array (conditional on pressure levels being present)
    out = np.full(shape=(len(dataset.pfull.values), len(rs)-1), fill_value=np.nan) if 'pfull' in dataset.dims else np.full(shape=(len(rs)-1), fill_value=np.nan) 
    # 6. Create the annulus using a loop based on the limiting dimension and resolution
    # Note: use a basis vector from the center (0) to the domain edge
    # Note: begin from index 1 to establish the outer annular radius
    for index, r_o in enumerate(rs[1:]):
        # Uncomment to get metrics
        # print('Annulus #{0}: inner radius = {1:.2f}, outer radius = {2:.2f}'.format(index, r_i, r_o))

        # 7. Average over annular area.
        average = dataset.where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o)).mean(dim='grid_xt').mean(dim='grid_yt')
        out[:, index] = average.sortby('pfull', ascending=True)
        # Append this output to a storage list for easy debugging, if needed.
        if 'pfull' in dataset.dims: 
            visuals['Annulus #{0}:\ninner rad. = {1:.2f}, out. rad. = {2:.2f}'.format(index, r_i, r_o)] = \
            dataset.where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o)).sel(pfull=850, method='nearest')   
        else:
            visuals['Annulus #{0}:\ninner rad. = {1:.2f}, out. rad. = {2:.2f}'.format(index, r_i, r_o)] = \
            dataset.where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o))
        
        # 8. Expand annulus by a grid cell and repeat steps 6-8 until domain edge is reached.
        r_i += resolution
    # Populate the Dataset with the composite output
    data[storm_id][field] = xr.DataArray(data=out, dims=('pfull', 'radius'), coords={'pfull': (['pfull'], dataset.pfull.values), 
                                                                                    'radius': (['radius'], rs[:-1])},
                                        attrs=dataset.attrs).sortby('pfull', ascending=True)

    return storm_id, data[storm_id]

def pproc(data, field):
    """
    [P]arallel [proc]essing ([P] + [proc] = pproc) method for azimuthal compositing.

    Args:
        data (dict): dictionary with keys of storm IDs and values of an xArray Dataset
        field (str): field to composite within the xArray Dataset

    Returns:
        data (dict): dictionary with keys of storm IDs and values of an xArray Dataset
    """
    
    # Assemble the input lists for Pool.starmap
    inputs = [(data, id_, field) for id_ in list(data.keys())]
    # Use 16 processors to calculate azimuthal means in parallel
    with multiprocessing.get_context("spawn").Pool(16) as p:
        results = [result for result in p.starmap(azimuthal_mean, inputs)]
    # Rebuild dictionary with composites
    for result in results:
        data[result[0]] = result[1]
        
    return data
''' End azimuthal compositing support methods. '''

def azimuthal_compositor(model, datasets, intensity_bin, field, parallel = True):
    
    """
    Method to generate azimuthal composites for a list of given TCs for a field, 
    pressure level (if not a surface or integrated quantity), and an intensity bin.

    Args:
        model (str): name of model
        datasets (list): list of data dictionaries (.pkl file)
        intensity_bin (str): string
        field (str): field to be evaluated
    Returns:
        data (dictionary): dictionary with list of processed data. Dictionary is 2-tiered: (1) storm ID --> (2) field
        composite_mean (xArray DataArray): composite mean in radial (and pressure, if applicable) coordinates 
    """
    
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
        times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times if not ((t.month == 2) & (t.day == 29))]
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
        # This try/except is supposed to catch duplicate time indices in ds[subdict]
        try:
            dataset = ds[subdict][field].sel(time=index_timestamps)
        except:
            # Check for duplicate time values
            ds[subdict] = ds[subdict].drop_duplicates('time')
            # Now try reloading
            dataset = ds[subdict][field].sel(time=index_timestamps)
            
         # If timestamps match, proceed. Else, continue.
        if len(dataset.time.values) > 0:
            print('\t \t {0} has {1} matching timestamps! Further processing them now...'.format(storm_id, len(dataset.time.values)))
            # Methodology: extract all timestamps with matching intensity for the given storm and assign them a storm 'sub-ID', 
            # which preserves the storm ID but gives it a unique identifier.
            storm_subids = [] # create storage list for all storm sub-IDs to be generated
            
            ''' Populate dictionary with data. '''
            # Methodology:  (1) get matching timestamps between GCM output and TC tracker output. 
            #               (2) then, get the timestamp corresponding to LMI from the TC tracker output.
            #               (3) finally, get timestamps from +/- N days from the LMI timestamp, where N is defined in the script as 'n_days'
            n_days = 1
            # Get Pandas-friendly timestamps for DataFrame indexing
            timestamps = [utilities.time_adjust(model, t, method='cftime_to_pandas') 
                          for t in dataset.time.values if not ((t.month == 2) & (t.day == 29))] 
            # Get strongest timestamp
            timestamp_lmi = ds['track_output'].loc[ds['track_output']['time'].isin(timestamps)].sort_values('max_wind', ascending=False).iloc[0]['time'] 
            # Get the timestamps found from +/- n_days and revert to cftime for DataArray indexing
            timestamps = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in timestamps 
                          if abs(t - timestamp_lmi) <= pd.Timedelta(n_days, "d")]
            
            # Process data for each filtered timestamp
            for i, timestamp in enumerate(timestamps):   
                # Assign the storm sub-ID
                storm_subid = '{0}_{1:03d}'.format(storm_id, i)
                # Extract all nans
                data[storm_subid] = {key: dataset.sel(time=timestamp).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
                # Append to storage list
                storm_subids.append(storm_subid)
                print('\t \t \t {0} created for timestamp {1}, continuing...'.format(storm_subid, timestamp))
            for storm_subid in storm_subids:
                # Check to see if the extents are less than the minima. If so, make the new minima
                if (min_x == None) or (len(data[storm_subid][key].grid_xt) < min_x):
                    min_x = len(data[storm_subid][key].grid_xt)
                if (min_y == None) or (len(data[storm_subid][key].grid_yt) < min_y):
                    min_y = len(data[storm_subid][key].grid_yt)
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
            
    # Apply a parallel implementation to perform the azimuthal composite
    data = pproc(data, key)
    
    # Iterate through all datasets and obtain the limiting dimensions (height == pfull, width == radius)
    min_height, min_width = min([v[key].shape[0] for v in data.values()]), min([v[key].shape[1] for v in data.values()])
    # Trim each dataset to the limiting dimension
    for k in data.keys():
        data[k] = {key: data[k][key][0:min_height, 0:min_width].sortby('pfull', ascending=True)}
    # Get sample dataset for axis building in the composite array.
    sample_dataset = [v[key] for v in data.values()][0]
    
    # Build xArray with dimensions corresponding to the requested field
    if subdict == 'tc_vertical_output':
        # Get composite mean over all storms by averaging along the 0th axis (storm ID axis)
        # Generate the output xArray DataArray
        composite_mean = xr.DataArray(data=np.nanmean(np.stack([v[key] for v in data.values()]), axis=0), 
                                      dims=['pfull', 'radius'], 
                                      coords={'pfull': (['pfull'],sample_dataset['pfull'].values),
                                              'radius': (['radius'], sample_dataset['radius'].values)}).sortby('pfull', ascending=True)
    elif subdict == 'tc_model_output':
        # Get composite mean over all storms by averaging along the 0th axis (storm ID axis)
        # Generate the output xArray DataArray
        composite_mean = xr.DataArray(data=np.nanmean(np.stack([v[key] for v in data.values()]), axis=0), 
                                      dims=['radius'], coords={'radius': (['radius'], sample_dataset['radius'].values)})
    
    return key, data, composite_mean
