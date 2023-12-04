import numpy as np, pandas as pd
import matplotlib

def get_cmap(param, difference=False):
    
    colormaps = pd.read_csv('/projects/GEOCLIM/gr7610/reference/param_colormaps.csv')

    if param not in colormaps['param'].values:
        param = 'other'

    if difference:
        return colormaps.loc[colormaps['param'] == param]['difference'].values[0]
    else:
        return colormaps.loc[colormaps['param'] == param]['normal'].values[0]

def norm_cmap(data, param, mode='raw', bounds=True, n_bounds=None, difference=False):
    ''' 
    This function takes input data and return a corresponding normalization
    and colormap based on property and mode. This function works for both
    single-entry and multi-entry datasets. Note that multi-entry datasets
    must be of the same property, since this function makes a common normalization
    across all entry datasets.
    '''
    
    # Convert data input to list if single-entry
    if type(data) is not list:
        data = [data]
    # Extract parameter
    if not difference:
        if 'numpy.ndarray' in str(type(data[0])):
            data = [data]
        else:
            print('update')
            data = [d[param] for d in data]

    ''' Normalization definition. '''
    # Get extrema from input dataset(s)
    vmin, vmax = np.nan, np.nan
    for i, dataset in enumerate(data): 
        # Correction factor for different parameters
        if (param in ['precip', 'evap']):
            factor = 86400 
        elif param == 'lhflx':
            factor = 2.5e6
        else:
            factor = 1
        dataset_min, dataset_max = np.nanmin(dataset)*factor, np.nanmax(dataset)*factor
        # Overwrite extremum if extrema surpassed
        if i == 0:
            vmin = dataset_min
            vmax = dataset_max
        else:
            vmin = dataset_min if dataset_min < vmin else vmin
            vmax = dataset_max if dataset_max > vmax else vmax

    if param == 'precip':
        vmax = vmax / 2
    if param == 'evap' and not difference:
        vmin = 0
        
    # Define extremum for CenteredNorm and bounded normalizations
    extremum = np.nanmax([abs(vmin), abs(vmax)])
    # Number of bounded levels
    n_bounds = 11 if not n_bounds else n_bounds

    # Catch floating point errors to accuracy of 1e-9
    vmin = np.around(vmin, decimals=9) + 0.0
    vmax = np.around(vmax, decimals=9) + 0.0
    
    print('Colormap extrema: ', vmin, vmax)
    
    # Define the type of normalization based on generated extrema
    if (np.sign(vmin) == np.sign(vmax)) or vmin == 0 or vmax == 0:
        if bounds:
            norm = matplotlib.colors.BoundaryNorm(np.linspace(vmin, vmax, n_bounds), 256)
        else:
            norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        # Get colormap
        cmap = get_cmap(param, difference=False)
    else:
        if bounds:
            norm = matplotlib.colors.BoundaryNorm(np.linspace(-extremum, extremum, n_bounds), 256)
        else:
            norm = matplotlib.colors.Normalize(vmin=-extremum, vmax=extremum)
        # Get colormap
        cmap = get_cmap(param, difference=True)
    
    return norm, cmap