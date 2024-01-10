import cftime, datetime
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle, random

def directories(model, experiment, data_type='model_output'):
    """
    Method to log the directories of all raw model runs and return a corresponding path.
    """
    
    dirnames = {'AM2.5': {'control': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_tigercpu_intelmpi_18_540PE/analysis_lmh',
                                      'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_tigercpu_intelmpi_18_540PE/POSTP'},
                          'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh',
                                      'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'}},
                'AM2.5C360': {'control': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                                          'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/POSTP'},
                              'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                                         'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/POSTP'}},
                'FLOR': {'control': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_tigercpu_intelmpi_18_576PE/analysis_lmh',
                                     'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_tigercpu_intelmpi_18_576PE/POSTP'},
                          'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_swishe_tigercpu_intelmpi_18_576PE/analysis_lmh',
                                     'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_swishe_tigercpu_intelmpi_18_576PE/POSTP'}},
                'HIRAM': {'control': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh',
                                      'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_tigercpu_intelmpi_18_540PE/POSTP'},
                          'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh',
                                     'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'}}}
    
    return dirnames[model][experiment][data_type]

def access(model, experiment, storm_type, storm_id=None, processed=False):
    
    """Access a random TC pickle file.
    
    Arguments:
        model (str): model name (AM2.5, HIRAM, FLOR)
        storm_type (str): TS or C15w
        storm_id (str, default: None): ID of storm of interest

    Returns:
        filename (str): filename string
        data (dict): 3-element dictionary with track outputs from Lucas Harris' TC tracker, planar model outputs from the chosen GCM, and vertical outputs from the chosen GCM
    """
    
    # Define addendum if proessed data is being accessed
    addendum = '/processed' if processed else ''
    # Retrieve filenames
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs{0}'.format(addendum)
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
        
    return filename, data

def time_adjust(model, timestamp, method='pandas_to_cftime'):
    """Method to adjust datetime conventions.

    Arguments:
        model (str): model name
        timestamp (multiple): Pandas, datetime, or cftime object.
        method (str, optional): conversion type. Defaults to 'pandas_to_cftime'.

    Returns:
        out (multiple): Pandas, datetime, or cftime object.
    """
    
    if method == 'pandas_to_cftime':
        year_adjust = 0 if model == 'FLOR' else 1900
        out = cftime.DatetimeNoLeap(year=timestamp.year-year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        return out
    else:
        return None
    
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

def intensity_binning(data, intensity_metric='max_wind'):
    """
    Method to generate bins for TC intensities to support compositing based on a given intensity metric.
    Note: this is chosen to be performed post-storage to keep raw TC data as clean as possible from metadata additions, in case of future definition changes.
    Note: the intensity metric will either be by maximum winds (max_wind) or minimum sea-level pressure (min_slp).
    """

    # Define the intensity bins
    if intensity_metric == 'max_wind':
        intensity_bin_limits = [0, 15, 20, 25, 30, 35, 40, np.inf]
    elif intensity_metric == 'min_slp':
        intensity_bin_limits = [np.inf, 1000, 990, 980, 960, 940, 920, 0]
    # Create the bin data structure, with bin numbers as keys and bin bounds and data as subdictionaries
    intensity_bins = {'b{0}'.format(i): {'bounds': (intensity_bin_limits[i], intensity_bin_limits[i+1]), 'storms': []} 
                    for i in range(0, len(intensity_bin_limits)-1)} 
    # Initialize data intensity bin column
    data['track_output']['intensity_bin'] = np.nan
    # Iterate over all bin and assign to each timestamp
    for bin, bin_data in intensity_bins.items():
        # Get intensity bounds
        min_bound, max_bound = min(intensity_bins[bin]['bounds']), max(intensity_bins[bin]['bounds'])
        # Assign bin name to timestamp
        data['track_output'].loc[(min_bound <= data['track_output'][intensity_metric]) & (data['track_output'][intensity_metric] < max_bound), 'intensity_bin'] = bin
            
    return data

def lmh_parser(path):
    
    ''' 
    This method parses through text files from Lucas Harris' run outputs (held in directories titled 'analysis_lmh') 
    and produces an output DataFrame. 
    
    Input(s):
    - path (str):            path containing raw tracker data from Lucas Harris' runs.
    Output(s):
    - df (Pandas DataFrame): Pandas DataFrame containing tracked TC data
    '''
    
    # Create file object instance
    fobj = open(path, 'r').readlines()
    # Initialize dictionary to hold data
    data = {'storm_num': {}}
    # Initialize storm counter
    count = 1
    # Iterate through text file
    for line in fobj:
        # Extract information from the line
        content = line.strip()
        # Creates new storm-specific dict in the parent dict. The '+++' demarcates a new storm.
        if '+++' in line:
            storm_num = '{0:04d}'.format(count)
            data['storm_num'][storm_num] = {'storm_id': [], 'time': [], 'lon': [], 'lat': [], 'slp': [], 'max_wnd': [], 'flag': []}
            count += 1
        # Populates the storm-specific dict
        else:
            storm_num = '{0:04d}'.format(count-1) 
            tc_info = [x for x in content.split(' ') if x]
            year = tc_info[0][0:4] # get 4-digit year
            data['storm_num'][storm_num]['storm_id'].append('{0}-{1:04d}'.format(year, count-1))
            data['storm_num'][storm_num]['time'].append(tc_info[0])
            data['storm_num'][storm_num]['lon'].append(tc_info[1])
            data['storm_num'][storm_num]['lat'].append(tc_info[2])
            data['storm_num'][storm_num]['slp'].append(tc_info[3])
            data['storm_num'][storm_num]['max_wnd'].append(tc_info[4])
            data['storm_num'][storm_num]['flag'].append(tc_info[5])
    
    try:
        # Converts the dictionary into a DataFrame
        df = pd.concat({k: pd.DataFrame(v).T for k, v in data.items()}, axis=1)['storm_num']
        df = df.explode(df.columns.to_list()).reset_index().rename(columns={'index': 'storm_num'})
        # Re-cast column data types
        df = df.astype({'lon': 'float', 'lat': 'float', 'slp': 'float', 'max_wnd': 'float', 'flag': 'float'})
    except:
        df = pd.DataFrame(columns=['storm_id', 'time', 'lon', 'lat', 'slp', 'max_wnd', 'flag'])
    
    ''' DataFrame refinement. '''
    # Remove cold-core data points (flag == -1)
    df = df.loc[df['flag'] != -1].reset_index(drop=True)
    # Convert timestamps to datetime objects
    df['time'] = pd.to_datetime(df['time'], format='%Y%m%d%H')
    
    return df

def retrieve_tracked_TCs(model, experiment, storm_type, year_range):
    
    '''
    Function to collect tracked TC data and add derived data, such as duration and storm speed.
    
    Input(s):
    - model (str):                name of model
    - experiment (str):           name of experiment
    - storm_type (str):           type of storm to evaluate from TC tracks data ("TS" for all storms or "C15w" for hurricanes)
    - year_range (tuple of ints): 2-element tuple with a start and end year
    Output(s):
    - data (Pandas DataFrame):    Pandas DataFrame with tracked TC data
    '''
    
    # Retrieve directory containing parent directory for track data
    dirname = directories(model, experiment, data_type='track_data')
    
    ''' File collection. '''
    # Get filenames for all files within the specified directory 
    # Filenames will correspond to the determined storm type    
    fnames = [[os.path.join(dirname, file, 'Harris.TC', f) for f in os.listdir(os.path.join(dirname, file, 'Harris.TC')) 
               if '{0}.world'.format(storm_type) in f]
               for file in sorted(os.listdir(dirname))]
    # Compress 2D list to 1D list
    fnames = [item for sublist in fnames for item in sublist]

    # Select files with dates within 'year_range'
    # Note: the '+ 1900' is added because tracked TCs are on the 2000 year range, whereas model output is on the 100 year range
    # Note: conditional added as exception for provisional FLOR data ahead of 10th NE Tropical Workshop
    if dirname == '/tigress/GEOCLIM/grios/HIRAM/exp/CTL1990_v201905/analysis_lmh/cyclones_gav_ro110_1C_330k':
        fnames = [f for f in fnames]
    else:
        year_adjust = 0 if 'FLOR' in dirname else 1900
        fnames = [f for f in fnames 
                  if min(year_range) + year_adjust <= pd.to_datetime(f.split('.')[-2].split('-')[0]).year <= max(year_range) + year_adjust]
    
    # Concatenate all tracked TC data from the filename list
    data = pd.concat([lmh_parser(os.path.join(dirname, fname)) for fname in fnames])
    
    ''' Derived track-based data algorithm. Storm-specific derived properties will be generated in here. '''
    
    # Initialize empty duration column to populate iteratively
    data[['duration', 'speed', 'direction']] = np.nan
    # Initialize list to populate iteratively for each storm, then concatenate
    storms = []
    # Iterate through each unique storm (identify by 'storm_id') and get duration
    for storm_id in data['storm_id'].unique():
        # Define iterand storm
        storm = data.loc[data['storm_id'] == storm_id].copy().reset_index(drop=True)
        
        ''' Duration derivation. '''
        # Get difference between minimum and maximum timestamps
        dt = (storm['time'].max() - storm['time'].min())
        # Convert difference timedelta into hours
        dt = dt.days + dt.seconds/86400
        # Add duration to the outer DataFrame for the corresponding storm
        data.loc[data['storm_id'] == storm_id, 'duration'] = dt
        # Re-define iterand storm to incorporate duration
        storm = data.loc[data['storm_id'] == storm_id].copy().reset_index(drop=True)
        
        ''' Velocity (speed, direction) derivation. '''
        # Initialize dictionary for preliminary storage. Will be reassigned into the DataFrame by the join() method using time as the matching criterion.
        velocity = {'time': [storm.iloc[0]['time']], 'speed': [np.nan], 'direction': [np.nan]}
        # Iterate over all of the iterand storm timestamps
        for i in range(1, len(storm)):
            # Define coordinates for two points considered (i, i-1)
            lon_a, lat_a = [storm.iloc[i-1]['lon'], storm.iloc[i-1]['lat']]
            lon_b, lat_b = [storm.iloc[i]['lon'], storm.iloc[i]['lat']]
            # Determine timedelta between points (i, i-1)
            dt = storm.iloc[i]['time'] - storm.iloc[i-1]['time']
            # Derive speed (distance / time in m s^-1)
            speed = coords_to_dist((lon_b, lat_b), (lon_a, lat_a))/dt.seconds
            # Get changes in longtiude and latitude
            dlon, dlat = lon_b - lon_a, lat_b - lat_a
            # Derive direction relative to north (range of 0 to 360)
            direction = 180*np.arctan(dlon/dlat)/np.pi % 360
            # Append quantities to the 'velocity' dictionary
            velocity['time'].append(storm.iloc[i]['time'])    
            velocity['speed'].append(speed)    
            velocity['direction'].append(direction)
        # Build DataFrame
        velocity = pd.DataFrame(velocity)
        # Re-cast time column as a datetime object
        velocity['time'] = pd.to_datetime(velocity['time'])
        # Merge the storm and velocity DataFrames
        storm = storm.merge(velocity, how='left', on='time', suffixes=['_x', None]).drop(columns={'speed_x', 'direction_x'}).reset_index(drop=True)
        # Append to the list for future concatenation
        storms.append(storm)
        
    # Concatenate DataFrames
    data = pd.concat(storms)   
    # Rename columns for future addition into xArray Dataset, and reset index
    data = data.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp'}).reset_index(drop=True)
        
    return data