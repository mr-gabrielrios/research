''' Import packages. '''
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy
# Local utility packages
import dill, multiprocessing, os, pickle, utilities
# Data structure packages
import pandas as pd, xarray as xr
# Visualization tools
import cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

def tc_storage(model, experiment, year_range, storm_type, random_num=None, benchmarking=True, reference_time=None, offset=0):
    
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
    
    # Define output type based on storage conventions
    output_type = 'atmos_4xdaily'
    
    # Retrieve model data over specified year range. 
    # This is here so that both TC-specific and associated global model data can be analyzed together.
    model_output = retrieve_model_data(model, model_dir, year_range=(min(year_range), max(year_range)), output_type=output_type, benchmarking=benchmarking)
    
    if benchmarking:
        print('Model output access time: {0:.3f} s'.format(time.time() - lap))
        lap = time.time()
    
    # Retrieve model data specified to tracked TCs.
    # Note: default model output type is 'atmos_4xdaily'. If 'output_type' matches 'model_output', then access pre-loaded data from 'model_output'.
    storm_model_output, storm_ids = retrieve_model_TCs(model_dir, (min(year_range), max(year_range)), storm_track_output, 
                                                       model_output, output_type=output_type, random_num=random_num)
    
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
    storms = {}
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
        # Rename columns for future addition into xArray Dataset, and reset index
        storm = storm.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp'}).reset_index(drop=True)
        # Append to the list for future concatenation
        storms[storm_id] = storm
    
    # Concatenate DataFrames
    data = pd.concat(storms.values())   
    
    return data, storms

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

def retrieve_model_TCs(dirname, year_range, storms, model_output=None, output_type='atmos_daily', extent=15, num_storms=None, reference_time=None):

    # This boolean controls what frequency data is pulled at. Default to daily to allow for azimuthal compositing of vertical fields.
    timestamp_daily_freq = False if min(year_range) >= 157 else True
    
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
    storms_xr = {}
    # Define range of storm IDs to iterate over. 
    # If 'random_num' is defined, get that many randomized storms. Otherwise, get all.
    storm_ids = random.sample(list(storms['storm_id'].unique()), len(list(storms['storm_id'].unique())))
    # Container for list of storm IDs for processed storms
    output_ids = []

    print('\t ', storm_ids)
    print('\t Number of storms to be processed: {0}'.format(len(storm_ids)))
    
    storm_counter = 0
    # Access model data that are specific to each tracked TC.
    for storm_index, storm_id in enumerate(storm_ids):
        # Pull tracked TC data relative to storm
        storm = storms.loc[storms['storm_id'] == storm_id]
        # Find timestamp corresponding to storm LMI. LMI is defined as timestamp with maximum wind.
        lmi = storm.loc[storm['max_wind'] == storm['max_wind'].max()] 
        # Check for multiple instances of LMI - if multiple exist, get instance with lowest pressure
        if len(lmi) > 1:
            lmi = lmi.loc[lmi['min_slp'] == lmi['min_slp'].min()]
        
        ''' Pre-process filters. '''
        # 1. Ensure LMI occurs below a prescribed latitude
        if abs(lmi['center_lat'].values) >= 30:
            print('\t \t \t LMI occurs too far poleward for {0}, skipping...'.format(storm_id))
            continue
        # 2. Ensure LMI is not the first timestamp - this implies that it's a baroclinic storm.
        if lmi['time'].values == storm.sort_values('time').iloc[0]['time']:
            print('\t \t \t LMI occurs at first timestamp, skipping...'.format(storm_id))
            continue
        # 3: Ensure last TC timestamp is before 12/31 to avoid year-to-year overruns
        if storm['time'].max().month == 12 and storm['time'].max().day > 30:
            print('\t \t \t TC overruns into next year, skipping...'.format(storm_id))
            continue
        
        # Initialize list to hold Dataset entries for each storm timestamp
        storm_xr = []
        # Iterate over each storm Series
        for i in range(0, len(storm)):        
            # Convert from tracked TC timestamp convention (datetime) to model timestamp convention (cftime DatetimeNoLeap)
            try:
                if storm.iloc[i]['time'].year < 2101:
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
            # Assign storm identifier (storm_id)
            storm_xr = xr.concat(storm_xr, dim='time')
            storm_xr['storm_id'] = storm_id
        except:
            for i in range(0, len(storm_xr)):
                print(storm_xr[i].grid_xt.values)
            pass
        
        ''' Spatial clipping. '''
        # Part 1: compile model output from each timestamp in the storm dataset
        if timestamp_daily_freq:
            try:
                # Only grab single timestamp at LMI
                if reference_time == 'lmi':
                    # Alternate for azimuthal compositing - try finding storm point of interest, then get nearest day
                    pull_timestamp = storm_xr.where(storm_xr['max_wind'] == storm_xr['max_wind'].max(), drop=True).time.values[0]
                    pull_timestamp = datetime.datetime(pull_timestamp.year, pull_timestamp.month, pull_timestamp.day, pull_timestamp.hour)
                    # Rewrite timestamp to the corresponding day as a cftime object. Pick the nearest full day
                    if pull_timestamp.hour <= 12:
                        pull_timestamp = cftime.datetime(year=pull_timestamp.year, 
                                                         month=pull_timestamp.month, 
                                                         day=pull_timestamp.day, hour=0, calendar='noleap')
                    else:
                        try:
                            pull_timestamp = cftime.datetime(year=pull_timestamp.year, 
                                                             month=pull_timestamp.month, 
                                                             day=pull_timestamp.day+1, hour=0, calendar='noleap')
                        except:
                            print('Unable to process, part 1a:', storm_id)
                            pass
                    snapshot = storm_xr.sel(time=pull_timestamp, method='nearest')
                # If no reference time is specified, grab all timestamps for a given storm
                else:
                    # Initialize container list of storm snapshots
                    snapshots = []
                    # Iterate over each timestamp
                    for pull_timestamp in storm_xr.time.values:
                        # Rewrite timestamp to the corresponding day as a cftime object. Pick the nearest full day
                        if pull_timestamp.hour <= 12:
                            pull_timestamp = cftime.datetime(year=pull_timestamp.year, month=pull_timestamp.month, 
                                                             day=pull_timestamp.day, hour=0, calendar='noleap')
                        else:
                            try:
                                pull_timestamp = cftime.datetime(year=pull_timestamp.year, month=pull_timestamp.month,
                                                                 day=pull_timestamp.day+1, hour=0, calendar='noleap')
                            except:
                                print('Unable to process, part 1b:', storm_id)
                                pass
                        snapshot = storm_xr.sel(time=pull_timestamp, method='nearest')
                        snapshots.append(snapshot)
                    snapshot = xr.concat(snapshots, dim='time').drop_duplicates('time')
            except:
                print('Unable to process, part 1c:', storm_id)
                pass
        else:
            print('\t Using the new method!')
            # Initialize container list of storm snapshots
            snapshots = []
            # Iterate over each timestamp
            for pull_timestamp in storm_xr.time.values:
                snapshot = storm_xr.sel(time=pull_timestamp, method='nearest')
                snapshots.append(snapshot)
            snapshot = xr.concat(snapshots, dim='time').drop_duplicates('time')

        # Part 2: clip model data for each storm timestamp
        # Identify center as location of maximum absolute value of vorticity or minimum pressure
        center_param = 'slp'
        # try:
        clip_container = []
        print('Snapshot timestamps: {0}'.format(snapshot.time.values))
        # Iterate over each timestep
        for index, timestamp in enumerate(storm_xr.time.values):
            # Create year-adjusted timestamp
            year_adjust = 1900 if timestamp.year < 1900 else 0 # adjust year for AM2.5, HIRAM
            dt_timestamp = datetime.datetime(year=timestamp.year + year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
            # Match timestamps
            print('\t Storm ID: {0}, iterand timestamp: {1}'.format(storm_id, timestamp))
            # Get storm track center to locate the tracker-driven TC center
            temp_center_lon, temp_center_lat = [storm.loc[storm['time'] == dt_timestamp]['center_lon'].values,
                                                storm.loc[storm['time'] == dt_timestamp]['center_lat'].values]
            # Get temporary storm snapshot to locate the GCM output-driven TC center
            temp_snapshot = snapshot.sel(time=timestamp, 
                                        grid_xt=np.arange(temp_center_lon-extent, temp_center_lon+extent),
                                        grid_yt=np.arange(temp_center_lat-extent, temp_center_lat+extent), method='nearest').drop_duplicates(['grid_xt'])
            # Center-detection using chosen method (e.g., look for max(vort850), look for min(slp)
            if center_param == 'vort850':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).max().values, temp_snapshot[center_param], np.nan)
            elif center_param == 'slp':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).min().values, temp_snapshot[center_param], np.nan)
            # Get extrema for chosen method from GCM output
            center_detector = center_detector.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            # Distill to single float values
            center_lon, center_lat = center_detector.grid_xt.values[0], center_detector.grid_yt.values[0]
            print('\t Storm ID: {0}, coordinates: {1}, {2}'.format(storm_id, center_lat, center_lon))
            # Get snapshot
            temp_snapshot = snapshot.sel(time=timestamp, 
                                         grid_xt=np.arange(center_lon-extent, center_lon+extent), 
                                         grid_yt=np.arange(center_lat-extent, center_lat+extent), 
                                         method='nearest').drop_duplicates(['grid_xt'])
            temp_snapshot['center_lon'], temp_snapshot['center_lat'] = center_lon, center_lat
            clip_container.append(temp_snapshot)
            print('\t ', temp_snapshot.time.values)
        clip_container = xr.concat(clip_container, dim='time').drop_duplicates(dim='time')
        storms_xr[storm_id] = clip_container
        output_ids.append(storm_id)
        print('\t Append success for ', storm_id, clip_container.time.values, '\n')
        storm_counter += 1
        
        # Limit number of storms read if a limit is given
        if num_storms and (storm_counter >= num_storms):
            return storms_xr, output_ids

def vertical_profile_selection(storms, model, experiment):

    '''
    Arguments:
    - storms (list)
    '''
    start = time.time()
    
    # Vertical data before year 157 of HIRAM runs is output daily, 8xdaily after
    output_type = 'atmos_8xdaily' if min(year_range) >= 157 else 'atmos_daily'
    # Determine output frequency
    output_daily = True if output_type == 'atmos_daily' else False
    # Initialize container dictionary
    container = {}

    # Iterate over each storm xarray Dataset
    for storm in storms.values():
        print('Storm ID: {0}'.format(storm['storm_id']))
        # Get unique storm ID. Skip if more than one is found.
        storm_ids = list(set(storm['storm_id'].values))
        if len(storm_ids) == 1:
            storm_id = storm_ids[0]
        else:
            continue
        
        # Initialize storm-specific data container
        output = []
        # Get timestamp
        timestamps = storm.time.values
        # Iterate over timestamps
        for timestamp in timestamps:
            print('\t Processing timestamp: {0} at {1} frequency'.format(timestamp, output_type))
            snapshot = storm.sel(time=timestamp)
            # Drop null values
            snapshot = snapshot.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            # Define bounds from xt/yt maxima
            grid_xt_bounds, grid_yt_bounds = [(snapshot.grid_xt.min(), snapshot.grid_xt.max()), (snapshot.grid_yt.min(), snapshot.grid_yt.max())]
            # Rewrite timestamp to the corresponding day as a cftime object. 
            # Pick the nearest full day if daily data is used, otherwise use the raw time
            if output_daily:
                if timestamp.hour <= 12:
                    timestamp = cftime.datetime(timestamp.year, timestamp.month, timestamp.day, hour=0, calendar='noleap')
                else:
                    timestamp = timestamp + datetime.timedelta(days=1)
                    timestamp = cftime.datetime(timestamp.year, timestamp.month, timestamp.day, hour=0, calendar='noleap')
            print('\t ---> Vertical timestamp out: {0}'.format(timestamp))
            # Get path name corresponding to model and experiment
            model_dir, _ = dirnames(model, experiment)
            # Pull full vertical output and select corresponding day
            vertical_profile = xr.open_dataset('{0}/{1:04d}0101.{2}.nc'.format(model_dir, timestamp.year, output_type))
            # Filter by time and spatial domain
            try:
                vertical_profile = vertical_profile.sel(time=timestamp, 
                                                        grid_xt=slice(min(grid_xt_bounds), max(grid_xt_bounds)), 
                                                        grid_yt=slice(min(grid_yt_bounds), max(grid_yt_bounds))).load()
            except:
                print('\t Unable to retrieve vertical data for {0}'.format(snapshot['storm_id']))
                pass
            # Assign storm-specific snapshot values
            storm_params = ['center_lon', 'center_lat', 'max_wind', 'min_slp', 'core_temp', 'speed', 'heading', 'storm_id']
            # Assign fields to output
            for storm_param in storm_params:
                vertical_profile[storm_param] = snapshot[storm_param]

            output.append(vertical_profile)
            print('\t Successful processing of vertical timestamp.')
        # Concatenate all timestamps
        container[storm_id] = xr.concat(output, dim='time')

    return container

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

def main(model, experiments, storm_type, year_range, num_storms, storage=False, override=False):

    data = {model: {}}
    
    # Iterate over experiments to load and process data
    for experiment_num, experiment in enumerate(experiments):
        # data[model][experiment] = {}
        # Define paths to model- and experiment-specific data.
        model_dir, track_dir = dirnames(model, experiment)
    
        # For a given year:
        for year in year_range:
            data[model][experiment] = {}
            years = [year, year+1]
            print('\n-----------------------------------------------------------------------')
            print('Processing year: {0}'.format(year))
            # Get model data for future use in TC processing
            if model == 'HIRAM' and year >= 157:
                output_type = 'atmos_8xdaily'
            elif model == 'HIRAM' and year >= 152 and year < 157:
                output_type = 'atmos_4xdaily'
            else:
                output_type = 'atmos_daily'
            model_output = retrieve_model_data(model, model_dir, year_range=years, output_type=output_type)
            # Get tracked TCs from Lucas Harris' TC tracker
            track_output, storm_track_output  = retrieve_tracked_TCs(track_dir, storm_type, year_range=years, basins=None)
            # Pair model data to tracked TC for entire storm lifetime
            storm_model_output, storm_ids = retrieve_model_TCs(model_dir, (min(year_range), max(year_range)), 
                                                       track_output, model_output, output_type=output_type, num_storms=num_storms)
            # Pair vertical model data to tracked TC for entire storm lifetime
            vertical_storm_output = vertical_profile_selection(storm_model_output, model, experiment)
            
            for storm_id in storm_ids:
                data[model][experiment][storm_id] = {}
                data[model][experiment][storm_id]['track_output'] = storm_track_output[storm_id]
                data[model][experiment][storm_id]['tc_model_output'] = storm_model_output[storm_id]
                data[model][experiment][storm_id]['tc_vertical_output'] = vertical_storm_output[storm_id]
                data[model][experiment][storm_id]['run_name'] = model_dir

                # Storage for each individual storm
                if storage:
                    # Define directory name and file name
                    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
                    # Define filename using max wind and min SLP for future binning
                    max_wind, min_slp = storm_track_output[storm_id]['max_wind'].max(), storm_track_output[storm_id]['min_slp'].min()
                    storage_filename = 'TC-{0}-{1}-{2}-{3}-max_wind-{4:0.0f}-min_slp-{5:0.0f}.pkl'.format(model, experiment, storm_type, storm_id, max_wind, min_slp)
                    storage_path = os.path.join(storage_dirname, storage_filename)
                    # If file doesn't exist, save
                    if not os.path.isfile(os.path.join(storage_dirname, storage_filename)) and override:
                        with open(os.path.join(storage_dirname, storage_filename), 'wb') as f:
                            print('Saving to: ', os.path.join(storage_dirname, storage_filename))
                            pickle.dump(data[model][experiment][storm_id], f)

    return data

if __name__ == '__main__':
    start = time.time()
    # Use range(start, stop) for a range of years between 'start' and 'stop', and a list [start, stop] for specific years.
    year_range = range(152, 157)
    data = main('HIRAM', experiments=['control', 'swishe'], storm_type='C15w', year_range=year_range, 
                num_storms=10, storage=True, override=True)
    print('Elapsed total runtime: {0:.3f}s'.format(time.time() - start))