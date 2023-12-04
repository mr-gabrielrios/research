''' Import packages. '''
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy
# Local utility packages
import dill, multiprocessing, os, pickle
# Data structure packages
import pandas as pd, xarray as xr
# Visualization tools
import cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")
def tc_storage(model, experiment, year_range, storm_type, random_num=None, benchmarking=True, reference_time='lmi', offset=0):
    
    ''' 
    Function to build unified xArray Dataset to hold model and track data for tracked TCs in a given model and experiment.
    
    Input(s):
    - model (str):                name of model to access (typically "AM2.5" or "HIRAM")
    - experiment (str):           name of experiment to access (typically "control", "ktc", "plus2K", etc.)
    - year_range (tuple of ints): 2-element tuple with a start and end year
    - storm_type (str):           type of storm to evaluate from TC tracks data ("TS" for all storms or "C15w" for hurricanes)
    - benchmarking (bool):        boolean to enable time benchmarking
    - reference_time (str):       time at which to get data 
    - offset (int):               number of 6-hour periods by which to offset time selection (for examples, 2 = 12 hours ahead)
    Output(s):
    - dirs (tuple):               2-element tuple with pathnames for the model- and experiment-specific GCM and track directories.
    '''
    
    if benchmarking:
        start = time.time()
    
    # Define paths to model- and experiment-specific data.
    model_dir, track_dir = dirnames(model, experiment)
    
    if benchmarking:
        print('Directory access time: {0:.3f} s'.format(time.time() - start))
        lap = time.time()
    
    # Retrieve tracked TCs over the specified year range. Note: year_range re-specified to ensure increasing order in tuple.
    storm_track_output = retrieve_tracked_TCs(track_dir, storm_type, year_range=(min(year_range), max(year_range)),
                                              basins=None)
    
    if benchmarking:
        print('Track access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    # Retrieve model data over specified year range. 
    # This is here so that both TC-specific and associated global model data can be analyzed together.
    model_output = retrieve_model_data(model, model_dir, year_range=(min(year_range), max(year_range)), output_type='atmos_4xdaily', benchmarking=benchmarking)
    
    if benchmarking:
        print('Model output access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    # Retrieve model data specified to tracked TCs.
    # Note: default model output type is 'atmos_4xdaily'. If 'output_type' matches 'model_output', then access pre-loaded data from 'model_output'.
    storm_model_output, storm_ids = retrieve_model_TCs(model_dir, (min(year_range), max(year_range)), storm_track_output, 
                                                       model_output, output_type='atmos_4xdaily', random_num=random_num)
    
    if benchmarking:
        print('TC accessing output access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    # Get radius for each storm to use for future normalization
    storm_model_output = add_radius(storm_model_output, benchmarking=benchmarking)
    
    if benchmarking:
        print('Radius derivation access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()

    # Get full vertical data for the corresponding storm dates and locations
    vertical_storm_output = vertical_profile_selection(storm_model_output, model, experiment)
    
    if benchmarking:
        print('Vertical profile/daily data access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    if benchmarking:
        print('Total runtime: {0:.3f} s'.format(time.time() - start))
    
    return storm_track_output, storm_model_output, vertical_storm_output

def dirnames(model, experiment):

    ''' 
    Function to store pathnames for selected models and experiments. 
    
    Input(s):
    - model (str):  name of model to access (typically "AM2.5" or "HIRAM").
    Output(s):
    - dirs (tuple): 2-element tuple with pathnames for the model- and experiment-specific GCM and track directories.
    '''
    # Model directories 
    model_dirs = {'AM2.5': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_tigercpu_intelmpi_18_540PE/POSTP',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'},
                  'AM2.5C360': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/POSTP',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/POSTP'},
                  'FLOR': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_tigercpu_intelmpi_18_576PE/POSTP',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_swishe_tigercpu_intelmpi_18_576PE/POSTP'},
                  'HIRAM': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_tigercpu_intelmpi_18_540PE/POSTP',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'}}
    # Track directories
    track_dirs = {'AM2.5': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k'},
                  'AM2.5C360': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k'},
                  'FLOR': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_v201905_swishe_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k'},
                  'HIRAM': {'control': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                            'swishe': '/tiger/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k'}}
        
    dirs = (model_dirs[model][experiment], track_dirs[model][experiment])
    
    return dirs

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
            data['storm_num'][storm_num] = {'storm_id': [], 'time': [], 'lon': [], 'lat': [], 'slp': [], 'max_wind': [], 'flag': []}
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
            data['storm_num'][storm_num]['max_wind'].append(tc_info[4])
            data['storm_num'][storm_num]['flag'].append(tc_info[5])
    
    try:
        # Converts the dictionary into a DataFrame
        df = pd.concat({k: pd.DataFrame(v).T for k, v in data.items()}, axis=1)['storm_num']
        df = df.explode(df.columns.to_list()).reset_index().rename(columns={'index': 'storm_num'})
        # Re-cast column data types
        df = df.astype({'lon': 'float', 'lat': 'float', 'slp': 'float', 'max_wind': 'float', 'flag': 'float'})
    except:
        df = pd.DataFrame(columns=['storm_id', 'time', 'lon', 'lat', 'slp', 'max_wind', 'flag'])
    
    ''' DataFrame refinement. '''
    # Remove cold-core data points (flag == -1)
    df = df.loc[df['flag'] != -1].reset_index(drop=True)
    # Convert timestamps to datetime objects
    df['time'] = pd.to_datetime(df['time'], format='%Y%m%d%H')
    
    return df

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

def retrieve_tracked_TCs(dirname, storm_type, year_range, basins=None):

    '''
    Function to collect tracked TC data and add derived data, such as duration and storm speed.
    
    Input(s):
    - dirname (str):              name of directory containing files of interest
    - storm_type (str):           type of storm to evaluate from TC tracks data ("TS" for all storms or "C15w" for hurricanes)
    - year_range (tuple of ints): 2-element tuple with a start and end year
    Output(s):
    - data (Pandas DataFrame):    Pandas DataFrame with tracked TC data
    '''
    
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
                  if min(year_range) + year_adjust <= pd.to_datetime(f.split('.')[-2].split('-')[0]).year < max(year_range) + year_adjust]
    
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
            # Define coordinates fofr two points considered (i, i-1)
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
    # Drop high-latitude storms
    data = data.loc[abs(data['center_lat']) >= 30]

    # Select specific basins, if option is chosen
    if basins:
        basin_dict = {}
        for basin in basins:
            if basin == 'EP':
                # Only get EPAC storms
                basin_dict[basin] = data.where((data['center_lat'] > 0) & 
                                               (data['center_lat'] < 40) &
                                               (data['center_lon'] > 220) &
                                               (data['center_lon'] <= 285)).dropna()
    
    return data

def pull_gcm_data(dirname, year, output_type):
    
    ''' Method to read data for given parameters. '''
    
    # Get filenames for corresponding files
    filename = '{0:04d}0101.{1}.nc'.format(year, output_type)
    try:
        fname = os.path.join(dirname, filename)
        # Retrieve data
        data = xr.open_dataset(fname)
    except:
        return None
    
    return data

def retrieve_model_data(model, dirname, year_range, output_type='atmos_month', benchmarking=False):
    
    '''
    Method to open experiment-specific GCM output data for a specified year range and model output type.
    
    Input(s):
    - dirname (str):            directory name for the given model and experiment type.
    - year_range (tuple, list): tuple or list (minimum of 2 items) containing year range for desired data.
    - output_type (str):        string denoting GCM output desired.
    Output(s):
    - data (xArray Dataset):    xArray Dataset containing concatenated data for the year range selected.
    '''
    
    if benchmarking:
        start = time.time()
    
    # Identifier substring definition. This will be used for splitting the filename for identification purposes.
    substring = '0101.' if model in ['AM2.5', 'HIRAM'] else '0001.'
    # Access parent directory with experiment-specific model data and list directories corresponding to year range
    files = [os.path.join(dirname, file) for file in os.listdir(dirname) for year in range(min(year_range), max(year_range)) 
             if '{0:04d}'.format(year) in file.split(substring)[0] and output_type in file]
    # Store file data into an xArray Dataset.
    # Note: benchmarking showed ~1.5 s for opening 4 files using 'open_mfdataset', and at least 10x longer using 'open_dataset' + 'xr.concat'
    data = xr.open_mfdataset(files)
    
    if benchmarking:
        print('\t Model data retrieval time: {0:.3f} s'.format(time.time() - start))
    
    return data

def retrieve_daily_data(dirname, model_output, storm_ids, output_type='atmos_daily', reference_time='lmi', offset=0, benchmarking=False):
    
    '''
    Method to open experiment-specific GCM daily output data specific to TC dates and locations.
    Time selection dependent on chosen reference time (genesis or LMI) and corresponding offset from the reference time.
    
    Input(s):
    - dirname (str):                   directory name for the given model and experiment type.
    - model_output (Pandas DataFrame): Pandas DataFrame with TC track data.
    - output_type (str):               string denoting GCM output desired.
    - reference_time (str):                      string denoting reference time desired ('genesis' or 'LMI').
    - offset (int):                    number of 6-hour periods offset from the reference time.
    Output(s):
    - ds (xArray Dataset):             xArray Dataset containing concatenated data, TC-specific.
    - ds_clim (xArray Dataset):        xArray Dataset containing concatenated data, global for specific variables.
    '''
    
    if benchmarking:
        start = time.time()
    
    # Spatial extent of clipping around TCs
    extent = 15
    # Initialize lists to hold dates
    ds, ds_clim = [], []
    # Define range of storm IDs to iterate over. 
    # If 'random_num' is defined, get that many randomized storms. Otherwise, get all.
    # storm_ids = random.sample(list(model_output['storm_id'].unique()), random_num) if random_num else model_output['storm_id'].unique()
    # for each storm in the DataFrame
    for storm_id in storm_ids:
        # Get data for the iterand storm
        storm = model_output.loc[model_output['storm_id'] == storm_id]
        ### reference time conditional
        ### if genesis, get index of the first time
        if reference_time == 'genesis':
            index = storm.drop_duplicates(subset=['storm_id'], keep='first')
        ### if lmi, get index of the maximum wind and, subsequently, minimum pressure
        else:
            index = storm.loc[storm['max_wind'] == storm['max_wind'].max()]
            if len(index) > 1:
                index = index.loc[index['min_slp'] == index['min_slp'].min()]
        # get date corresponding to the index chosen in the conditional and apply offset
        date = pd.to_datetime(index.time.values[0]) + datetime.timedelta(hours=offset*6)
        
        # Define cftime date equivalent (hour 12 used to match model output data convention)
        if date.year > 2050:
            cfdate = cftime.DatetimeNoLeap(year=date.year-500, month=date.month, day=date.day, hour=12)
        else:
            cfdate = cftime.DatetimeNoLeap(year=date.year-1900, month=date.month, day=date.day, hour=12)
        # Identifier substring definition. This will be used for splitting the filename for identification purposes.
        substring = '0101.'
        # Define file name
        fname = [os.path.join(dirname, file) for file in os.listdir(dirname) 
                 if str(cfdate.year) in file.split(substring)[0] and output_type in file][0]
        # Pull data
        try:
            temp = xr.open_dataset(fname).sel(time=cfdate)
        except:
            pass
        # clip spatially with extent larger than storm radius to understand surrounding environment
        try:
            temp = temp.sel(grid_xt=np.arange(index['center_lon'].values-extent, index['center_lon'].values+extent), 
                            grid_yt=np.arange(index['center_lat'].values-extent, index['center_lat'].values+extent), method='nearest')
            # Append to list for future concatenation
            ds.append(temp)
        except:
            pass
        # For large-scale, pare down the data to specific variables
        clim_vars = ['precip']
        temp_clim = xr.open_dataset(fname)[clim_vars]
        ds_clim.append(temp_clim)
    # concatenate
    ds = xr.concat(ds, dim='time')
    ds_clim = xr.concat(ds_clim, dim='time')
    
    if benchmarking:
        print('\t Model data retrieval time: {0:.3f} s'.format(time.time() - start))
        
    return ds, ds_clim

def retrieve_model_TCs(dirname, year_range, storms, model_output=None, output_type='atmos_daily', extent=10, random_num=None, storm_id_list=None, max_lat=40):

    # This boolean controls what frequency data is pulled at. Default to daily to allow for azimuthal compositing of vertical fields.
    timestamp_daily_freq = True
    
    # Check to see if model_output (previously-access model data) is the same output type as desired for TCs. If so, pull from that Dataset.
    if model_output and output_type == model_output.attrs['filename'].split('.')[0]:
        data = model_output
    else:
        # Identifier substring definition. This will be used for splitting the filename for identification purposes.
        substring = '0101.'
        # Access parent directory with experiment-specific model data and list directories corresponding to year range
        files = [os.path.join(dirname, file) for file in os.listdir(dirname) for year in range(min(year_range), max(year_range)) 
                 if '{0:04d}'.format(year) in file.split(substring)[0] and output_type in file]
        # Store file data into an xArray Dataset.
        data = xr.open_mfdataset(files)

    # Initialize list to hold each storm Dataset for future concatenation
    storms_xr = []
    # Define range of storm IDs to iterate over. 
    # If 'random_num' is defined, get that many randomized storms. Otherwise, get all.
    ''' Previous attempt. '''
    if storm_id_list is not None:
        storm_ids = storm_id_list
    else:
        # storm_ids = random.sample(list(storms['storm_id'].unique()), random_num) if random_num else storms['storm_id'].unique()
        storm_ids = random.sample(list(storms['storm_id'].unique()), 
                                len(list(storms['storm_id'].unique())))

    print('\t ', storm_ids)
    print('\t Number of storms to be processed: {0}'.format(len(storm_ids)))
    
    storm_counter = 0
    # Access model data that are specific to each tracked TC.
    for storm_id in storm_ids:
        print('\t \t', storm_id, storm_counter, random_num)
        if storm_counter > random_num:
            break
        # Pull tracked TC data relative to storm
        storm = storms.loc[storms['storm_id'] == storm_id]
        # Initialize list to hold Dataset entries for each storm timestamp
        storm_xr = []
        # Iterate over each storm Series
        for i in range(0, len(storm)):        
            # Convert from tracked TC timestamp convention (datetime) to model timestamp convention (cftime DatetimeNoLeap)
            try:
                if storm.iloc[i]['time'].year < 2050:
                    year_adjust = 0 if 'FLOR' in dirname else 1900
                    cf_timestamp = cftime.DatetimeNoLeap(year=storm.iloc[i]['time'].year-year_adjust, month=storm.iloc[i]['time'].month, 
                                                         day=storm.iloc[i]['time'].day, hour=storm.iloc[i]['time'].hour)
                else:
                    cf_timestamp = cftime.DatetimeNoLeap(year=storm.iloc[i]['time'].year-500, month=storm.iloc[i]['time'].month, 
                                                         day=storm.iloc[i]['time'].day, hour=storm.iloc[i]['time'].hour)
                # Take snapshot of data at this timestamp
                # Note 1: this is where the connection between track and model output data happens
                # Note 2: any operations beyond selection result in active computation, which disrupts the 'lazy' approach
                if cf_timestamp in data.time.values:
                    snapshot = data.sel(time=cf_timestamp)
                    # Get iterand information to append to Dataset for the given timestamp
                    snapshot[['center_lon', 'center_lat', 'min_slp', 
                              'max_wind', 'core_temp', 'speed', 'heading']] = [storm.iloc[i]['center_lon'], storm.iloc[i]['center_lat'], 
                                                                              storm.iloc[i]['min_slp'], storm.iloc[i]['max_wind'], 
                                                                              storm.iloc[i]['core_temp'], storm.iloc[i]['speed'], storm.iloc[i]['direction']]  
                    snapshot = snapshot.drop_duplicates(['grid_xt'])
                    # Append to list
                    storm_xr.append(snapshot)
            except:
                print('Unable to process:', storm_id)
                pass
            
        # Concatenate into unified Dataset for a storm   
        try:
            storm_xr = xr.concat(storm_xr, dim='time')
            storm_xr['storm_id'] = storm_id
        except:
            for i in range(0, len(storm_xr)):
                print(storm_xr[i].grid_xt.values)
        # Assign storm identifier (storm_id)
        # Note: to select specific storm from concatenated Dataset, use command as follows:
        #       storms.where(storms['storm_id'] == <storm_id>, drop=True)
        
        ''' Spatial clipping - only at maximum intensity to save computation time. '''
        # Grab storm at timestamp with maximum wind speed
        # Calling this 'lifetime maximum intensity', or 'LMI', to use terminology from Wing et al, 2016 (10.1175/JCLI-D-18-0599.1)
        # Note: this will be revisited if a non-LMI timestamp is requested
        if timestamp_daily_freq:
            try:
                # Alternate for azimuthal compositing - try finding storm point of interest, then get nearest day
                pull_timestamp = storm_xr.where(storm_xr['max_wind'] == storm_xr['max_wind'].max(), drop=True).time.values[0]
                # print(storm_xr.time.values, pull_timestamp)
                pull_timestamp = datetime.datetime(pull_timestamp.year, pull_timestamp.month, pull_timestamp.day, pull_timestamp.hour)
                # Rewrite timestamp to the corresponding day as a cftime object. Pick the nearest full day
                if pull_timestamp.hour <= 12:
                    pull_timestamp = cftime.datetime(pull_timestamp.year, pull_timestamp.month, pull_timestamp.day, calendar='noleap')
                else:
                    pull_timestamp = pull_timestamp + datetime.timedelta(days=1)
                    try:
                        pull_timestamp = cftime.datetime(pull_timestamp.year, pull_timestamp.month, pull_timestamp.day, calendar='noleap')
                    except:
                        print('Unable to process:', storm_id)
                        pass
                snapshot = storm_xr.sel(time=pull_timestamp, method='nearest')
            except:
                print('Unable to process:', storm_id)
                pass
        else:
            snapshot = storm_xr.where(storm_xr['max_wind'] == storm_xr['max_wind'].max(), drop=True)
            # If multiple wind maxima are detected, get the storm with lower pressure
            snapshot = lmi.where(snapshot['min_slp'] == snapshot['min_slp'].min(), drop=True) if len(snapshot.time.values) > 1 else snapshot

        # Clip storm  
        # Identify center as location of maximum absolute value of vorticity
        center_param = 'vort850'
        try:
            # Restrict search area in model data
            temp_center_lon, temp_center_lat = [storm.loc[storm['max_wind'] == storm['max_wind'].max()]['center_lon'].values,
                                                storm.loc[storm['max_wind'] == storm['max_wind'].max()]['center_lat'].values]
            # print(temp_center_lon, temp_center_lat)
            search_extent = 10
            temp_snapshot = snapshot.sel(grid_xt=np.arange(temp_center_lon-search_extent, temp_center_lon+search_extent), 
                                         grid_yt=np.arange(temp_center_lat-search_extent, temp_center_lat+search_extent), method='nearest').drop_duplicates(['grid_xt'])
            if center_param == 'vort850':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).max().values, temp_snapshot[center_param], np.nan)
            elif center_param == 'slp':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).min().values, temp_snapshot[center_param], np.nan)
            center_detector = center_detector.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            center_lon, center_lat = center_detector.grid_xt.values[0], center_detector.grid_yt.values[0]
            # print('\t', center_lon, center_lat)
            snapshot = snapshot.sel(grid_xt=np.arange(center_lon-extent, center_lon+extent), 
                                    grid_yt=np.arange(center_lat-extent, center_lat+extent), method='nearest').drop_duplicates(['grid_xt'])
            snapshot['center_lon'] = center_lon
            snapshot['center_lat'] = center_lat
            # Append to list for future concatenation assuming latitude of center is in tropics
            if abs(center_lat) <= 32:
                print(storm_counter)
                storms_xr.append(snapshot)
                storm_counter += 1
        except:
            print('Unable to process:', storm_id)
            pass
        
    # Concatenate into single xArray Dataset
    storms_xr = xr.concat(storms_xr, dim='time')
    
    return storms_xr, storm_ids

def radius_estimate(storm, box, overlay_check=True, benchmarking=False, troubleshooting=True):
    
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

    # Troubleshooting boolean added to prevent radius addition to see if it allows more weak TCs to be identified
    if troubleshooting: 
        return storm
    else:
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
            print('\t Unable to get radius for {0}'.format(storm['storm_id'].values))
            print('------------------------------------------------')
            print(storm)
            print('------------------------------------------------')
            return np.nan
        
def add_radius(storms, test_num=None, benchmarking=False):
    
    start = time.time()
    
    samples = random.sample(list(set(storms['storm_id'].values)), test_num) if test_num else list(set(storms['storm_id'].values))

    container = []
    # Embarrassingly parallel
    for storm_id in samples:
        
        if benchmarking:
            lap = time.time()
        
        storm = storms.where(storms['storm_id'] == storm_id, drop=True).load()
        
        storm = storm.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        storm = radius_estimate(storm, box=10, overlay_check=False, benchmarking=False)
        
        if benchmarking:
            print('\t \t Dropping method: {0:.3f} s'.format(time.time() - lap))

        if 'xarray' in str(type(storm)):
            container.append(storm)
    storms = xr.concat(container, dim='time')
    
    print('Radius estimation runtime per storm: {0:.3f} s'.format((time.time() - start)/len(samples)))
    
    return storms

def vertical_profile_selection(storms, model, experiment):

    start = time.time()

    container = []

    for storm_id in list(set(storms['storm_id'].values)):
        # Get reference storm from atmos_4xdaily output
        storm_reference = storms.where(storms['storm_id'] == storm_id, drop=True)
        # Remove extraneous data, spatially, to obtain longitude (grid_xt) and latitude (grid_yt) bounds
        storm_reference = storm_reference.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        # Get timestamp
        timestamp = storm_reference.time
        # Define bounds from xt/yt maxima
        grid_xt_bounds, grid_yt_bounds = [(storm_reference.grid_xt.min(), storm_reference.grid_xt.max()),
                                          (storm_reference.grid_yt.min(), storm_reference.grid_yt.max())]
        
        # Rewrite timestamp to the corresponding day as a cftime object. Pick the nearest full day
        if timestamp.dt.hour <= 12:
            timestamp = cftime.datetime(timestamp.dt.year, timestamp.dt.month, timestamp.dt.day, calendar='noleap')
        else:
            timestamp = timestamp + datetime.timedelta(days=1)
            timestamp = cftime.datetime(timestamp.dt.year, timestamp.dt.month, timestamp.dt.day, calendar='noleap')
        # Get path name corresponding to model and experiment
        model_dir, _ = dirnames(model, experiment)
        # Pull full vertical output and select corresponding day
        storm_full = xr.open_dataset('{0}/{1:04d}0101.atmos_daily.nc'.format(model_dir, timestamp.year))
        # Filter by time and spatial domain
        try:
            storm_full = storm_full.sel(time=timestamp, 
                                        grid_xt=slice(min(grid_xt_bounds), max(grid_xt_bounds)), 
                                        grid_yt=slice(min(grid_yt_bounds), max(grid_yt_bounds))).load()
        except:
            print('\t Unable to retrieve vertical data for {0}'.format(storm_id))
            pass
        # Assign storm-specific snapshot values
        storm_params = ['center_lon', 'center_lat', 'max_wind', 'min_slp', 'core_temp', 'speed', 'heading', 'storm_id', 'U']

        for storm_param in storm_params:
            storm_param_name = storm_param if storm_param != 'U' else 'U_ref'
            storm_full[storm_param] = storm_reference[storm_param]

        container.append(storm_full)

    output = xr.concat(container, dim='time')

    return output

def rel_hum(data):
    ''' 
    Calculate relative humidity (H) as a function of atmospheric pressure, specific humidity, and temperature.

    H = e(p, q)/e_s(T)
    e = p/(eps/q - e + 1) # see Emanuel (1994), Eq. 4.1.4
    e_s = 6.112*exp(17.67*t/(t+243.5)) # see Bolton (1980), Eq. 10. 
    
    Bolton (1980) doi: https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    ''' 

    p, q, t = data.pfull, data['sphum'], data['temp']

    eps = 0.622 # ratio of R_d to R_v
    e = p/(eps/q - eps + 1)

    tc = t - 273.16
    e_s = 6.112*np.exp(17.67*tc/(tc + 243.5))

    data['rh'] = 100*e/e_s
    print((100*e/e_s).shape)
    data['rh'].attrs['long_name'] = 'relative humidity'
    data['rh'].attrs['units'] = '%'

    return data

def moisture_flux(data):
    ''' 
    Calculate moisture convergence.
    ''' 

    du_dx = data['ucomp'].diff(dim='grid_xt')
    dv_dy = data['vcomp'].diff(dim='grid_yt')
    dq_dx = data['sphum'].diff(dim='grid_xt')
    dq_dy = data['sphum'].diff(dim='grid_yt')
    
    data['moisture_adv'] = (du_dx + dv_dy)*data['sphum']
    data['moisture_adv'].attrs['long_name'] = 'horizontal moisture advection'
    data['moisture_adv'].attrs['units'] = 's^-1'
    
    data['moisture_div'] = dq_dx*data['ucomp'] + dq_dy*data['vcomp']
    data['moisture_div'].attrs['long_name'] = 'horizontal moisture divergence'
    data['moisture_div'].attrs['units'] = 's^-1'
    
    data['moisture_flux'] = data['moisture_div'] + data['moisture_adv']
    data['moisture_flux'].attrs['long_name'] = 'horizontal moisture flux'
    data['moisture_flux'].attrs['units'] = 's^-1'
    
    return data

def radial_tangential_velocity(data, x, y, R):
    '''
    Calculate radial and tangential velocity components from zonal and meridional velocities.
    '''

    u, v = data['ucomp'], data['vcomp']

    if np.nanmin(data.grid_yt) < 0:
        u, v = -u, -v
    
    theta = np.arctan(v/u)
    data['v_radial'] = (u*x + v*y)/R
    data['v_tangential'] = (v*x - u*y)/R
    
    data['v_radial'].attrs = {'long_name': 'radial velocity', 'units': 'm s$^{-1}$'}
    data['v_tangential'].attrs = {'v_tangential': 'tangential velocity', 'units': 'm s$^{-1}$'}
    
    return data

def domain_temp_anomaly(data):

    data['temp_anom'] = data['temp'] - data['temp'].isel(time=0).mean(dim='grid_xt').mean(dim='grid_yt')
    
    data['temp_anom'].attrs = {'long_name': 'domainwise temperature anomaly', 'units': 'K'}

    return data    

def differences(data, param, mode='relative', data_type='atmos_month', months=[1, 12]):

    ctl, exp = [data['control'][param].sel(time=((data['control'].time.dt.month >= min(months)) & 
                                                 (data['control'].time.dt.month < max(months)))).mean(dim='time').load(),
                data['swishe'][param].sel(time=((data['swishe'].time.dt.month >= min(months)) &
                                                (data['swishe'].time.dt.month < max(months)))).mean(dim='time').load()]
    
    diff = (exp - ctl) if mode == 'absolute' else (100*(exp - ctl)/(exp + ctl))

    # Custom parameter bounds
    if mode == 'relative':
        if param in ['precip', 'evap', 'netrad_toa']:
            diff = diff.where(abs(diff) <= 25)
        else:
            diff = diff.where(abs(diff) <= 100)

    return diff

def output_storage(data, models, experiments, storm_type, year_range):
    ''' Method to store processed data, TC-specific. '''

    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage'
    for model in models:
        for experiment in experiments:
            num = len(data[model][experiment]['tc_model_output'].time.values)
            filename = 'tc_output.model-{0}.exp-{1}.storm_type-{2}.years-{3}_{4}.num-{5}.npy'.format(model, experiment, storm_type, min(year_range), max(year_range), num)
            output_path = os.path.join(dirname, filename)

            import sys
            print('Object size: {0:.2f} B'.format(sys.getsizeof(data[model][experiment])))

            try:
                np.save(output_path, data[model][experiment])
            except:
                import pickle
                with open(output_path, 'wb') as f:
                    pickle.dump(data[model][experiment], f, protocol=pickle.HIGHEST_PROTOCOL)

def main(model, experiments, storm_type, year_range, num_storms, storage=False):

    data = {model: {}}
    
    num_procs = 6
    # Create chunks from year range based on number of processors
    chunks = np.linspace(min(year_range), max(year_range), num_procs+1, dtype=int)
    # Create chunks
    chunks = [(chunks[i], chunks[i+1]) for i, _ in enumerate(chunks[:-1])]
    
    # Iterate over experiments to load and process data
    for experiment_num, experiment in enumerate(experiments):
        # Generate inputs for parallel procesing
        inputs = []
        # Initialize expeirment-specific dictionary entry 
        data[model][experiment] = {}
        # Create inputs by iterating over chunks
        for chunk_num, year_chunk in enumerate(chunks):
            input_list = [model, experiment, year_chunk, storm_type, num_storms, True, 'lmi', 0]
            inputs.append(input_list)
        [print('{0}: {1}\n'.format(experiment, input)) for input in inputs]
        # Process data in parallel (each processor corresponds to a year chunk and field)
        with multiprocessing.get_context("spawn").Pool(num_procs) as p:
            container = p.starmap(tc_storage, inputs)
        p.close()
        
        # Load starmap data into dictionary by concatenating like items 
        # (index 0 = track output, index 1 = planar model output at 6 h frequency, index 2 = full vertical output at daily frequency)
        data[model][experiment]['track_output'] = pd.concat([container[i][0] for i in range(0, len(container))])
        data[model][experiment]['tc_model_output'] = xr.concat([container[i][1] for i in range(0, len(container))], dim='time')
        data[model][experiment]['tc_vertical_output'] = xr.concat([container[i][2] for i in range(0, len(container))], dim='time')
    
    #         data[model][experiment]['tc_vertical_output'] = moisture_flux(data[model][experiment]['tc_vertical_output'])
    #         data[model][experiment]['tc_vertical_output'] = rel_hum(data[model][experiment]['tc_vertical_output'])
    #         data[model][experiment]['tc_vertical_output'] = domain_temp_anomaly(data[model][experiment]['tc_vertical_output'])
       
    if storage:     
        output_storage(data, model, experiments, storm_type, year_range)
        
    return data

if __name__ == '__main__':
    start = time.time()
    data = main('HIRAM', experiments=['control', 'swishe'], storm_type='TS', year_range=[101, 150], num_storms=100, storage=False)
    print('Elapsed total runtime: {0:.3f}s'.format(time.time() - start))

# Note to self
### Serial test run: 1 proc, 10 years, 5 storms --> 307 s
### Serial test run: 4 proc, 10 years, 5 storms --> 117 s
### Serial test run: 8 proc, 10 years, 5 storms --> 140 s

# fig, axes = plt.subplots(ncols=2, figsize=(9, 3.5))

# model_name = 'AM2.5C360'
# experiment = 'control'
# vertical_param = 'omega'

# index = 0
# data[model_name][experiment]['tc_model_output'].isel(time=index).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')['U'].plot(ax=axes[0])
# vertical_temp = data[model_name][experiment]['tc_vertical_output'].sel(time=data[model_name][experiment]['tc_model_output'].isel(time=index).time.values).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
# center_x, center_y = len(vertical_temp.grid_xt)//2, len(vertical_temp.grid_yt)//2
# vertical_temp.isel(grid_yt=slice(center_y - 2, center_y + 2))[vertical_param].mean(dim='grid_yt').plot(ax=axes[1])

# storm_id = data[model_name][experiment]['tc_model_output'].isel(time=index)['storm_id']

# # Check to see center alignment
# x, y = [data[model_name][experiment]['tc_model_output'].isel(time=index).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all').grid_xt.values,
#         data[model_name][experiment]['tc_model_output'].isel(time=index).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all').grid_yt.values]

# axes[0].scatter(x[len(x)//2], y[len(y)//2], c='r', s=100)

# fig.suptitle(data[model_name][experiment]['tc_model_output'].isel(time=index)['storm_id'].values)
# fig.tight_layout()