import numpy as np, pandas as pd, scipy as sp, xarray as xr
import utilities

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
            if i < len(data.grid_yt.values) - 1:
                lon_next, lon_curr = data.grid_xt.values[i+1], data.grid_xt.values[i]
            # Handle the boundary condition by assuming periodicity in x
            else:
                lon_next, lon_curr = data.grid_xt.values[0], data.grid_xt.values[i]
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