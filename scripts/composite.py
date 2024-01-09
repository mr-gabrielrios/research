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
        pressure_level: ()
    Returns:
        _type_: _description_
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
        # Get timestamps relevant to the intensity bin at hand
        times = ds['track_output']['time'].loc[ds['track_output']['intensity_bin'] == intensity_bin]
        # Change from Pandas to cftime formats
        times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times]
        # Determine which dictionary the field data is in
        subdict = None
        for sd in ['tc_model_output', 'tc_vertical_output']:
            subdict = sd if field in ds[sd].data_vars else None
        # Get the matching timestamps between the track output and the corresponding dictionary with field data
        index_timestamps = [t for t in times if t in ds[subdict][field].time.values]
        # and field + vertical level from the vertical data
        ds = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[field].sel(time=times)
        # If timestamps match, proceed. Else, continue.
        if len(ds.time.values) > 0:
            # Populate dictionary with data. Limit to one entry per storm by choosing the first timestamp.
            data[storm_id] = {key: ds.isel(time=0).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
            # Check to see if the extents are less than the minima. If so, make the new minima
            if (min_x == None) or (len(data[storm_id][key].grid_xt.values) < min_x):
                min_x = len(data[storm_id][key].grid_xt.values)
            if (min_y == None) or (len(data[storm_id][key].grid_yt.values) < min_y):
                min_y = len(data[storm_id][key].grid_yt.values)
        else:
            continue
    
    # Define buffer to prevent edge effects of domain trimming. 
    # Will be applied to each edge, such that the final domain will be 2*buffer less after trimming.
    buffer = 1
    
    for storm_id in data.keys():
        domain_x, domain_y = len(data[storm_id][key].grid_xt.values), len(data[storm_id][key].grid_yt.values)
        center_x, center_y = domain_x // 2, domain_y // 2
        data[storm_id][key] = data[storm_id][key].iloc(grid_xt=)
        
    return data

# if __name__ == '__main__':
#     dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed'
#     model, experiment = 'HIRAM', 'control'
#     storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] 
#                  for f in os.listdir(dirname)
#                  if (model in f) and (experiment in f)]
#     datasets = []
#     for storm_id in storm_ids:
#         _, dataset = utilities.access(model, experiment, storm_type='C15w', storm_id=storm_id, processed=True)
#         datasets.append(dataset)