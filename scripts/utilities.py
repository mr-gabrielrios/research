import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle

def access(model, experiment, storm_type, storm_id=None):
    
    """Access a random TC pickle file.
    
    Arguments:
        model (str): model name (AM2.5, HIRAM, FLOR)
        storm_type (str): TS or C15w
        storm_id (str, default: None): ID of storm of interest

    Returns:
        data (dict): 3-element dictionary with track outputs from Lucas Harris' TC tracker, planar model outputs from the chosen GCM, and vertical outputs from the chosen GCM
    """
    
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    files = [os.path.join(dirname, filename) for filename in os.listdir(dirname)
            if model in filename and experiment in filename and storm_type in filename]
    
    # If a specific storm ID is given, check for it
    storm_exists = False # check for storm existence - True if exists
    if storm_id:
        # Get list of storms that match the storm ID
        storm_check_list = [file for file in files if storm_id in file]
        # If the storm exists and the list length is 1, get the filename
        if len(storm_check_list) == 1:
            storm_exists = True
            filename = storm_check_list[0]
        elif len(storm_check_list) > 1:
            print(storm_check_list)
        
    # Load the storm - random if no storm ID is given or found, storm ID if given and found
    if storm_exists:
        with open(filename, 'rb') as f:
            data = pickle.load(f)
    else:
        with open(random.choice(files), 'rb') as f:
            data = pickle.load(f)
        
    return data

def get_constants(name):
    
    # Reference: Emanuel (1994)
    constants = {'c_p': 1005.7,
                 'L_v': 2.5e6, 
                 'R_d': 287.04,
                 'R_v': 461.5,
                 'eps': 0.622,
                 'g': 9.81}
    
    return constants[name]

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

def distance_grid(data):

    ''' Method to approximate the distance between grid points on the GFDL cubed-sphere grid. '''
    
    distances = np.full(shape=(len(data.grid_yt), len(data.grid_xt)), fill_value=np.nan)
    for i, lat in enumerate(data.grid_yt.values):
        for j, lon in enumerate(data.grid_xt.values):
            if i < (len(data.grid_yt.values) - 1) and j < (len(data.grid_xt.values) - 1):
                lon_next, lon_curr = data.grid_xt.values[j+1], data.grid_xt.values[j]
            # Handle the boundary condition by assuming periodicity in x
            else:
                lon_next, lon_curr = data.grid_xt.values[0], data.grid_xt.values[j]
            distances[i, j] = coords_to_dist((lon_curr, lat), (lon_next, lat))
    out = xr.DataArray(data=distances,
                       dims=('grid_yt', 'grid_xt'),
                       coords=[data.grid_yt, data.grid_xt])
    return out

def domain_differentiation(data, distance, field, dim):

    ''' 
    Differentiate an xarray DataArray field along a given dimension to preserve shape. 
    Last row or column in the bulk differentiation array will be appended, effectively repeating that last row or column.
    Data must be 3-dimensional (time, grid_xt, grid_yt).
    '''

    # Trim unused dimensions
    for drop_dim in ['bnds', 'phalf']:
        data = data.drop_dims(drop_dim) if drop_dim in data.dims else data
    
    # Ensure proper dimensional order
    data = data.transpose('time', 'grid_xt', 'grid_yt', 'pfull')
    
    # Get the bulk differentiation (will result in an output with the shape of data, minus one entry in the differentiation dimension)
    a = data[field].diff(dim=dim)/distance
    # Get the last row/columns of the differentiation array and repeat to append to bulk array and preserve original dimensions
    b = a[{dim: -1}]
    # Concatenate along respective axes
    if dim == 'grid_xt':
        b_ = b.values[:, np.newaxis, :, :]
        c = xr.DataArray(data=np.concatenate((a.values, b_), axis=1), dims=a.dims)
    elif dim == 'grid_yt':
        b_ = b.values[:, :, np.newaxis, :]
        c = xr.DataArray(data=np.concatenate((a.values, b_), axis=2), dims=a.dims)
    return c

def compositor(model, experiment, intensity_range):
    """
    Method to generate composite for a given set of parameters from offloaded model data (.pkl files)
    """
    
    