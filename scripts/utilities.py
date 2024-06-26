import cftime, datetime
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle, random

import warnings
warnings.filterwarnings("ignore")

def directories(model, experiment, data_type='model_output'):
    """
    Method to log the directories of all raw model runs and return a corresponding path.
    """
    
    dirnames = {'AM2.5': {'control': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s/POSTP'},
                          'swishe': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_swishe/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_swishe/POSTP'}},
                'AM2.5C360': {'control': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                                          'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/POSTP'},
                              'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                                         'model_output': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/POSTP'}},
                'FLOR': {'control': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s/POSTP'},
                          'swishe': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s/POSTP'}},
                'HIRAM': {'control': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s/POSTP'},
                          'swishe': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s_swishe/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                     'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s_swishe/POSTP'}},
                'HIRAM-8xdaily': {'control': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s-8xdaily_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                              'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s-8xdaily_tigercpu_intelmpi_18_540PE/POSTP'},
                                  'swishe': {'track_data': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                             'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_540PE/POSTP'}}}
    
    return dirnames[model][experiment][data_type]

def month_letter(month):
    month_letters = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    return month_letters[month-1]

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
    
    # Define addendum if processed data is being accessed
    addendum = '/processed' if processed else ''
    # Retrieve filenames
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs{0}'.format(addendum)
    files = [os.path.join(dirname, filename) for filename in os.listdir(dirname)
            if model in filename and experiment in filename and storm_type in filename]
    # Initialize filename container
    filename = None
    # If a specific storm ID is given, check for it
    storm_exists = False # check for storm existence - True if exists
    
    print(files, storm_id)
    
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

def time_adjust(model=None, timestamp=None, method='pandas_to_cftime'):
    """Method to adjust datetime conventions.

    Arguments:
        model (str): model name
        timestamp (multiple): Pandas, datetime, or cftime object.
        method (str, optional): conversion type. Defaults to 'pandas_to_cftime'.

    Returns:
        out (multiple): Pandas, datetime, or cftime object.
    """
    # Convert formats.
    # Note: converting from 'cftime' to 'datetime' results in an addition of 1900 for models not affected by 'year_adjust' years 
    if method == 'pandas_to_cftime':
        # Define adjustment interval for select models
        year_adjust = 1900 if ((timestamp.year >= 1900) & (model != 'FLOR')) else 0
        out = cftime.DatetimeNoLeap(year=timestamp.year-year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        return out
    elif 'cftime_to_pandas':
        # Define adjustment interval for select models
        year_adjust = 1900 if timestamp.year >= 1900 else 0
        out = datetime.datetime(year=timestamp.year + 1900 - year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        return pd.to_datetime(out)
    
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
    ''' Convert coordinates to distance in meters using the Haversine formula. '''
    
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
    
    distance_lon = np.full(shape=(len(data.grid_yt), len(data.grid_xt)), fill_value=np.nan)
    for i, lat in enumerate(data.grid_yt.values):
        for j, lon in enumerate(data.grid_xt.values):
            if i < (len(data.grid_yt.values) - 1) and j < (len(data.grid_xt.values) - 1):
                lon_next, lon_curr = data.grid_xt.values[j+1], data.grid_xt.values[j]
            # Handle the boundary condition by assuming periodicity in x
            else:
                lon_next, lon_curr = data.grid_xt.values[0], data.grid_xt.values[j]
            distance_lon[i, j] = coords_to_dist((lon_curr, lat), (lon_next, lat))
    distance_lon = xr.DataArray(data=distance_lon, dims=('grid_yt', 'grid_xt'), coords=[data.grid_yt, data.grid_xt])
    
    sample_index = 6
    distance_lat = np.full(shape=(len(data.grid_yt), len(data.grid_xt)), 
                   fill_value=coords_to_dist((data.grid_xt.values[sample_index], data.grid_yt.values[sample_index]), 
                                             (data.grid_xt.values[sample_index], data.grid_yt.values[sample_index+1])))
    distance_lat = xr.DataArray(data=distance_lat, dims=('grid_yt', 'grid_xt'), coords=[data.grid_yt, data.grid_xt])
    
    return distance_lon, distance_lat

def area_weighted(data, function='average', extent=[0, 360, -60, 60]):
    """
    Helper function to generate area-weighted averages for a given field and domain (atmosphere or ocean).
    """
    
    # Check if 'ocean' is in the data dimensions. Rename dimensions if so
    if 'ocean' in ''.join(data.dims):
        data = data.rename({'xt_ocean': 'grid_xt', 'yt_ocean': 'grid_yt', 'st_ocean': 'pfull'})
        
    if function == 'average':
        mean = data.weighted(np.cos(np.deg2rad(data.grid_yt))).mean(['grid_xt', 'grid_yt'])
        std = data.weighted(np.cos(np.deg2rad(data.grid_yt))).std(['grid_xt', 'grid_yt'])
        return mean
    elif function == 'sum':
        sum = data.weighted(np.cos(np.deg2rad(data.grid_yt))).sum(['grid_xt', 'grid_yt'])
        return sum

def month_selector(month, start, end=None):
    ''' Helper function that takes a DataArray time.month index and integer arguments for start and end months for xArray time indexing, end-inclusive. '''
    end = start if not end else end
    return (month >= start) & (month <= end)

def domain_differentiation(data, dim):

    ''' 
    Differentiate an xarray DataArray field along a given dimension to preserve shape. 
    Last row or column in the bulk differentiation array will be appended, effectively repeating that last row or column.
    Data must be 3-dimensional (time, grid_xt, grid_yt).
    '''

    # Trim unused dimensions
    for drop_dim in ['bnds', 'phalf']:
        data = data.drop_dims(drop_dim) if drop_dim in data.dims else data
    
    # Ensure proper dimensional order
    if 'pfull' in data.dims:
        data = data.transpose('time', 'grid_xt', 'grid_yt', 'pfull')
    else:
        data = data.transpose('time', 'grid_xt', 'grid_yt')
    
    # Get distance between grid points (essentially, distance_lon = dx, distance_lat = dy)
    distance_lon, distance_lat = distance_grid(data)
    distance = distance_lon if dim == 'grid_xt' else distance_lat
    
    # Get the bulk differentiation (will result in an output with the shape of data, minus one entry in the differentiation dimension)
    a = data.diff(dim=dim)/distance
    # Get the last row/columns of the differentiation array and repeat to append to bulk array and preserve original dimensions
    b = a[{dim: -1}]
    
    # Concatenate along respective axes
    if dim == 'grid_xt':
        b_ = b.values[:, np.newaxis, :, :] if 'pfull' in data.dims else b.values[:, np.newaxis, :]
        c = xr.DataArray(data=np.concatenate((a.values, b_), axis=1), dims=a.dims)
    elif dim == 'grid_yt':
        b_ = b.values[:, :, np.newaxis, :] if 'pfull' in data.dims else b.values[:, :, np.newaxis]
        c = xr.DataArray(data=np.concatenate((a.values, b_), axis=2), dims=a.dims)
    return c

def intensity_binning(mode='track_output', data=None, intensity_metric='max_wind'):
    """
    Method to generate bins for TC intensities to support compositing based on a given intensity metric.
    Note: this is chosen to be performed post-storage to keep raw TC data as clean as possible from metadata additions, in case of future definition changes.
    Note: the intensity metric will either be by maximum winds (max_wind) or minimum sea-level pressure (min_slp).
    
    Args:
        mode (str):  analysis mode for the counting. Can either be (1) 'track_output' or (2) 'model_output'.
                     (1) 'track_output' refers to output from tc_analysis.tc_track_data(). 
                         This is meant to catalog all storms detected in the model runs, but not necessarily all analyzed for planar/azimuthal fields.
                     (2) 'model_output' refers to the 'track_output' from tc_analysis.tc_model_data().
                         This is meant to catalog all storms used for analysis for planar/azimuthal fields.
                     (3) all other entries will return a list of the intensity bins
        data (dict): dictionary output to match data accepted by 'track_output' or 'model_output'. See above for descriptions.
    """

    # Define the intensity bins
    if intensity_metric == 'max_wind':
        intensity_bin_limits = [0, 15, 32.5, np.inf] 
    elif intensity_metric == 'min_slp':
        intensity_bin_limits = [np.inf, 1000, 970, 940, 0]
    # Create the bin data structure, with bin numbers as keys and bin bounds and data as subdictionaries
    intensity_bins = {'b{0}'.format(i): {'bounds': (intensity_bin_limits[i], intensity_bin_limits[i+1]), 'storms': []} 
                    for i in range(0, len(intensity_bin_limits)-1)} 
     
    # Create copy of input data
    if mode == 'track_output':
        binned_data = data.copy()
    elif mode == 'model_output':
        binned_data = data['track_output'].copy()
    
    # Return intensity bins if a mode isn't specified, otherwise process as normal
    if mode in ['track_output', 'model_output']:
        # Initialize data intensity bin column
        binned_data['intensity_bin'] = np.nan
        
        # Iterate over all bin and assign to each timestamp
        for bin, bin_data in intensity_bins.items():
            # Get intensity bounds
            min_bound, max_bound = min(intensity_bins[bin]['bounds']), max(intensity_bins[bin]['bounds'])
            # Assign bin name to timestamp
            binned_data.loc[(min_bound <= binned_data[intensity_metric]) & (binned_data[intensity_metric] < max_bound), 'intensity_bin'] = bin
        
        # Append Series to the input dataset to return data in the same imput format
        if mode == 'track_output':
            data['intensity_bin'] = binned_data['intensity_bin']
        elif mode == 'model_output':
            data['track_output']['intensity_bin'] = binned_data['intensity_bin']
            
        return data
    else:
        return intensity_bins

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

def storm_snapshot(storm, mode='lmi'):
    ''' Grab single instance of storm from Pandas DataFrame (typically the 'track_output' key of 3-tiered data dictionaries). Could be genesis, LMI, or other. '''

    # Pick selected mode
    if mode == 'genesis': 
        # Sort storm by time
        storm = storm.sort_values('time', ascending=True)
    elif mode == 'lmi': 
        # Sort storm by maximum wind
        storm = storm.sort_values('max_wind', ascending=False)

    return storm.iloc[[0]]

def retrieve_tracked_TCs(model, experiment, storm_type, year_range, config=None):
    
    '''
    Function to collect tracked TC data and add derived data, such as duration and storm speed.
    
    Input(s):
    - model (str):                name of model
    - experiment (str):           name of experiment
    - storm_type (str):           type of storm to evaluate from TC tracks data ("TS" for all storms or "C15w" for hurricanes)
    - year_range (tuple of ints): 2-element tuple with a start and end year
    - config (str or None):       string to indicate functionality depending on which script calls it
    Output(s):
    - data (Pandas DataFrame):    Pandas DataFrame with tracked TC data
    '''
    
    # Retrieve directory containing parent directory for track data
    dirname = directories(model, experiment, data_type='track_data')
    
    print(year_range, dirname)
    
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
        year_adjust = 0 if 'FLOR' in dirname else 0
        print('Years evaluated from LMH output: {0} to {1}'.format(min(year_range) + year_adjust, max(year_range) + year_adjust))
        # fnames = [f for f in fnames 
        #           if min(year_range) + year_adjust <= pd.to_datetime(f.split('.')[-2].split('-')[0]).year <= max(year_range) + year_adjust]
        fnames = [f for f in fnames 
                  if min(year_range) + year_adjust <= datetime.datetime(year=int(f.split('/')[-3].split('_')[-1]), month=1, day=1).year <= max(year_range) + year_adjust]
    
    
    [print(fn, pd.to_datetime(f.split('.')[-2].split('-')[0]), f) for fn, f in enumerate(sorted(fnames))]
    
    # Concatenate all tracked TC data from the filename list
    data = pd.concat([lmh_parser(os.path.join(dirname, fname)) for fname in fnames])
    
    ''' Derived track-based data algorithm. Storm-specific derived properties will be generated in here. '''
    
    # Initialize empty duration column to populate iteratively
    data[['duration', 'speed', 'direction']] = np.nan
    # Initialize list to populate iteratively for each storm, then concatenate
    storms = {} if config == 'individual_tc_storage' else []
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
        velocity = {'time': [storm.iloc[0]['time']], 'speed': [np.nan], 'direction': [np.nan], 'ucomp': [np.nan], 'vcomp': [np.nan]}
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
            direction = np.mod(np.arctan2(dlon, dlat)*180/np.pi, 360)
            # Derive storm motion zonal and meridional
            ucomp = speed*np.sin(np.deg2rad(360-direction))
            vcomp = speed*np.cos(np.deg2rad(360-direction))
            # Append quantities to the 'velocity' dictionary
            velocity['time'].append(storm.iloc[i]['time'])    
            velocity['speed'].append(speed)    
            velocity['direction'].append(direction)
            velocity['ucomp'].append(ucomp)
            velocity['vcomp'].append(vcomp)
        # Build DataFrame
        velocity = pd.DataFrame(velocity)
        # Re-cast time column as a datetime object
        velocity['time'] = pd.to_datetime(velocity['time'])
        # Merge the storm and velocity DataFrames
        storm = storm.merge(velocity, how='left', on='time', suffixes=['_x', None]).drop(columns={'speed_x', 'direction_x'}).reset_index(drop=True)
        # Allow functionality for individual TC storage (see individual_tc_storage.py)
        if config == 'individual_tc_storage':   # Rename columns for future addition into xArray Dataset, and reset index
            storm = storm.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp', 'max_wnd': 'max_wind'}).reset_index(drop=True)
            # Append to the list for future concatenation
            storms[storm_id] = storm
        else:
            # Append to the list for future concatenation
            storms.append(storm)
       
    
       
    # Allow functionality for individual TC storage (see individual_tc_storage.py)
    if config == 'individual_tc_storage':
        # Concatenate DataFrames if storms are found. If not, return None.
        if len(storms.values()) > 0:
            data = pd.concat(storms.values())   
        else:
            data, storms = None, None
        return data, storms
    else:
         # Concatenate DataFrames
        data = pd.concat(storms)   
        # Rename columns for future addition into xArray Dataset, and reset index
        data = data.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp'}).reset_index(drop=True)
        return data

def file_counter(model_names=['AM2.5', 'FLOR', 'HIRAM'], num_files=10):
    """
    This method takes all tracked TCs from the Harris TC tracker and gets the ratio of SWISHE to control storms to ensure representative sampling for analytical methods.

    Args:
        num_files (int, optional): number of files desired. Defaults to 10.

    Returns:
        file_counts (dict): dictionary with the number of files to process per intensity bin.
    """
    # Load track data from a saved file to obtain the distribution of storm intensities
    filename = '/projects/GEOCLIM/gr7610/analysis/tc_storage/track_data-TS-2001_2050.pkl'
    with open(filename, 'rb') as f:
        data = pickle.load(f)
        
    # Collect data for output counts per model
    counts = {model: {} for model in model_names}
    # Parameter used for intensity binning
    param = 'max_wind'
    # Determine number of bins for wind speed
    wind_bins = np.arange(10, 50, 4)
    # Get the approximate number of storms per bin
    binwise_multiplier = int(np.ceil(num_files/len(wind_bins)))

    # Iterate over all models and experiments
    for model in model_names:
        counts[model] = {experiment: {} for experiment in data[model].keys()}
        for experiment in data[model].keys():
            # Get unique values for the parameter passed
            out = data[model][experiment]['unique'][param].dropna()
            # Filter out nans and infs
            out = out.loc[np.isfinite(out.values)]
            # Get bin values and bin edges
            n, bins = np.histogram(out, bins=wind_bins) 
            # Add approximate data count per intensity bin
            counts[model][experiment] = {bins[i]: n[i] for i in range(0, len(bins)-1)}
    
    # Build DataFrame to facilitate relative file number calculation per model/experiment configuration
    counts = {(model, experiment): values for model, model_data in counts.items() for experiment, values in model_data.items()}
    counts = pd.DataFrame.from_dict(counts, orient='columns')
    # Get the ratios
    for model in model_names:
        sample_model = model_names[0]
        counts[model, 'file_count_control'] = ((1-counts[sample_model]['swishe']/(counts[sample_model]['control'] + counts[sample_model]['swishe'])).fillna(0)*binwise_multiplier).round(0).astype(int)
        counts[model, 'file_count_swishe'] = ((counts[sample_model]['swishe']/(counts[sample_model]['control'] + counts[sample_model]['swishe'])).fillna(0)*binwise_multiplier).round(0).astype(int)
    # Get file counts and output as a dictionary
    file_counts = {model: {} for model in model_names}
    for model in model_names:
        for experiment in data[model].keys():
            file_counts[model][experiment] = counts[model, 'file_count_{0}'.format(experiment)].to_dict()
        
    return file_counts