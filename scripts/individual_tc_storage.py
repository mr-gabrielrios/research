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

# Import custom scripts
import derived, tc_processing, importlib, utilities
importlib.reload(utilities)
importlib.reload(tc_processing)
# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

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
    substring = '0101.' if model in ['AM2.5', 'HIRAM'] else '0101.'
    # Access parent directory with experiment-specific model data and list directories corresponding to year range
    files = [os.path.join(dirname, file) for file in os.listdir(dirname) for year in range(min(year_range), max(year_range)) 
             if '{0:04d}'.format(year) in file.split(substring)[0] and output_type in file]
    # Store file data into an xArray Dataset.
    # Note: benchmarking showed ~1.5 s for opening 4 files using 'open_mfdataset', and at least 10x longer using 'open_dataset' + 'xr.concat'
    data = xr.open_mfdataset(files)
    
    if benchmarking:
        print('\t Model data retrieval time: {0:.3f} s'.format(time.time() - start))
    
    return data

def retrieve_model_TCs(dirname, experiment, year_range, storms, model_output=None, output_type='atmos_4xdaily', extent=15, 
                       num_storms=None, maximum_latitude=35, reference_time=None, override=False, single_storm_id=None, diagnostic=False):

    # This boolean controls what frequency data is pulled at. Default to daily to allow for azimuthal compositing of vertical fields.
    timestamp_daily_freq = False if '8xdaily' in dirname else True
    output_type = 'atmos_4xdaily' # temporary override, 2024-05-16
    
    # Check to see if model_output (previously-access model data) is the same output type as desired for TCs. If so, pull from that Dataset.
    if model_output and output_type == model_output.attrs['filename'].split('.')[0]:
        data = model_output
    else:
        # Identifier substring definition. This will be used for splitting the filename for identification purposes.
        substring = '0101.'
        print(min(year_range))
        # Access parent directory with experiment-specific model data and list directories corresponding to year range
        files = [os.path.join(dirname, file) for file in os.listdir(dirname) for year in range(min(year_range), max(year_range)) 
                 if '{0:04d}'.format(year) in file.split(substring)[0] and output_type in file]
        print(dirname, output_type)
        # Store file data into an xArray Dataset.
        data = xr.open_mfdataset(files)

    # Initialize list to hold each storm Dataset for future concatenation
    storms_xr = {}
    # Define range of storm IDs to iterate over.
    # If 'random_num' is defined, get that many randomized storms. Otherwise, get all.
    storm_ids = random.sample(list(storms['storm_id'].unique()), len(list(storms['storm_id'].unique())))
    # Container for list of storm IDs for processed storms
    output_ids = []

    if diagnostic:
        print('[individual_tc_storage.py, retrieve_model_TCs()] pulling up storm IDs: {0}'.formt(storm_ids))
        print('[individual_tc_storage.py, retrieve_model_TCs()] number of storms to be processed: {0}'.format(len(storm_ids)))
    
    # Truncate number of storms if less than requested are found
    if len(storm_ids) < num_storms:
        num_storms = len(storm_ids)
    
    storm_counter = 0
    # Access model data that are specific to each tracked TC.
    for storm_index, storm_id in enumerate(storm_ids):
        # Use this for debugging a single storm ID
        if single_storm_id and storm_id != single_storm_id:
            continue
        # Check if already processed. If so, skip it assuming override is False
        if not override:
            override_flag = False
            # Check for matching storm ID and experiment
            for f in os.listdir('/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'):
                if storm_id in f and '{0}-'.format(experiment) in f:
                    override_flag = True
                    print('[individual_tc_storage.py, retrieve_model_TCs()] file {0} already exists.'.format(f))
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] override status: {0}'.format(override))
            if override_flag:
                print('[individual_tc_storage.py, retrieve_model_TCs()] {0} has already been processed and overrides are not being processed. Continuing...'.format(storm_id))
                continue
        # Pull tracked TC data relative to storm
        storm = storms.loc[storms['storm_id'] == storm_id]
        # Find timestamp corresponding to storm LMI. LMI is defined as timestamp with maximum wind.
        lmi = storm.loc[storm['min_slp'] == storm['min_slp'].min()] 
        # Check for multiple instances of LMI - if multiple exist, get instance with lowest pressure
        if len(lmi) > 1:
            lmi = lmi.loc[lmi['max_wind'] == lmi['max_wind'].max()]
        
        ''' Pre-process filters. '''
        # 1. Ensure LMI occurs below a prescribed latitude
        if abs(lmi['center_lat'].values) >= maximum_latitude:
            print(abs(lmi['center_lat'].values), maximum_latitude)
            print('[individual_tc_storage.py, retrieve_model_TCs()] LMI occurs too far poleward for {0}, skipping...'.format(storm_id))
            continue        
        # 2. Ensure LMI is not the first timestamp - this implies that it's a baroclinic storm.
        if lmi['time'].values == storm.sort_values('time').iloc[0]['time']:
            print('[individual_tc_storage.py, retrieve_model_TCs()] LMI occurs at first timestamp, skipping...'.format(storm_id))
            continue
        # 3: Ensure last TC timestamp is before 12/31 to avoid year-to-year overruns
        if storm['time'].max().month == 12 and storm['time'].max().day > 30:
            print('[individual_tc_storage.py, retrieve_model_TCs()] TC overruns into next year, skipping...'.format(storm_id))
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
                print('[individual_tc_storage.py, retrieve_model_TCs()] unable to process {0}, skipping the timestamp...'.format(storm_id))
                pass
        
        # Concatenate into unified Dataset for a storm. 
        # If it doesn't work, continue to the next loop iteration. 
        try:
            # Assign storm identifier (storm_id)
            storm_xr = xr.concat(storm_xr, dim='time')
            storm_xr['storm_id'] = storm_id
        except:
            continue
        
        if diagnostic:
            print('[individual_tc_storage.py, retrieve_model_TCs()] storm GCM output timestamps: {0}'.format(storm_xr.time.values))
        
        ''' Spatial clipping. '''
        # Part 1: compile model output from each timestamp in the storm dataset
        if timestamp_daily_freq:
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Getting daily timestamps...')
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
                            print('[individual_tc_storage.py, retrieve_model_TCs()] Unable to process, part 1a:', storm_id)
                            pass
                    if diagnostic:
                        print('[individual_tc_storage.py, retrieve_model_TCs()] Pulling data snapshot at {0}'.format(pull_timestamp))
                    snapshot = storm_xr.sel(time=pull_timestamp, method='nearest')
                # If no reference time is specified, grab all timestamps for a given storm
                else:
                    # Initialize container list of storm snapshots
                    snapshots = []
                    # Iterate over each timestamp
                    for pull_timestamp in storm_xr.time.values:
                        # # Rewrite timestamp to the corresponding day as a cftime object. Pick the nearest full day
                        # if pull_timestamp.hour <= 12:
                        #     pull_timestamp = cftime.datetime(year=pull_timestamp.year, month=pull_timestamp.month, 
                        #                                      day=pull_timestamp.day, hour=0, calendar='noleap')
                        # else:
                        #     try:
                        #         pull_timestamp = cftime.datetime(year=pull_timestamp.year, month=pull_timestamp.month,
                        #                                          day=pull_timestamp.day+1, hour=0, calendar='noleap')
                        #     except:
                        #         print('Unable to process for timestamp {0}, part 1b: {1}'.format(pull_timestamp, storm_id))
                        #         pass
                        if pull_timestamp.hour == 0:
                            pull_timestamp = cftime.datetime(year=pull_timestamp.year, month=pull_timestamp.month, 
                                                             day=pull_timestamp.day, hour=0, calendar='noleap')
                        else:
                            if diagnostic:
                                print('[individual_tc_storage.py, retrieve_model_TCs()] Unable to process for timestamp {0} due to mismatching hour, part 1b: {1}'.format(pull_timestamp, storm_id))
                            continue
                        if diagnostic:
                            print('[individual_tc_storage.py, retrieve_model_TCs()] Pulling data snapshot at {0}'.format(pull_timestamp))
                        snapshot = storm_xr.sel(time=pull_timestamp)
                        snapshots.append(snapshot)
                    snapshot = xr.concat(snapshots, dim='time').drop_duplicates('time')
            except:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Unable to process, part 1c:', storm_id)
                pass
        else:
            print('[individual_tc_storage.py, retrieve_model_TCs()] Using all timestamps!')
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
        if diagnostic:
            print('[individual_tc_storage.py, retrieve_model_TCs()] Snapshot timestamps: {0}'.format(snapshot.time.values))
            print('[individual_tc_storage.py, retrieve_model_TCs()] Storm_xr timestamps: {0}'.format(storm_xr.time.values))
        # Get common snapshot and storm timestamps
        common_timestamps = sorted(list(set(snapshot.time.values) & set(storm_xr.time.values)))
        # Iterate over each timestep
        for index, timestamp in enumerate(common_timestamps):
            # Create year-adjusted timestamp
            year_adjust = 1900 if timestamp.year < 1900 else 0 # adjust year for AM2.5, HIRAM
            dt_timestamp = datetime.datetime(year=timestamp.year + year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
            # Match timestamps
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Storm ID: {0}, iterand timestamp: {1}'.format(storm_id, timestamp))
            # Get storm track center to locate the tracker-driven TC center
            temp_center_lon, temp_center_lat = [storm.loc[storm['time'] == dt_timestamp]['center_lon'].values,
                                                storm.loc[storm['time'] == dt_timestamp]['center_lat'].values]
            temp_center_pressure = storm.loc[storm['time'] == dt_timestamp]['min_slp'].values[0]
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Tracker center: ({0:.2f}, {1:.2f}) with tracker min SLP {2:.2f} hPa'.format(temp_center_lon[0], temp_center_lat[0], temp_center_pressure))
            # Get temporary storm snapshot to locate the GCM output-driven TC center
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Looking at: {0} to {1}, {2} to {3}'.format(temp_center_lon-(extent), temp_center_lon+(extent),
                                                                temp_center_lat-(extent), temp_center_lat+(extent)))
            temp_snapshot = snapshot.sel(time=timestamp).sel(grid_xt=np.arange(temp_center_lon-(extent), temp_center_lon+(extent)),
                                                             grid_yt=np.arange(temp_center_lat-(extent), temp_center_lat+(extent)), method='nearest')
            # temp_snapshot = temp_snapshot.where(temp_snapshot['slp'] >= (temp_center_pressure))
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Minimum GCM pressure: {0:.2f} hPa over longitudes {1} to {2} and latitudes {3} to {4}'.format(temp_snapshot['slp'].min().values, 
                                                                                                                 temp_snapshot['grid_xt'].min().values,
                                                                                                                 temp_snapshot['grid_xt'].max().values,
                                                                                                                 temp_snapshot['grid_yt'].min().values,
                                                                                                                 temp_snapshot['grid_yt'].max().values))
            if diagnostic:
                print(temp_snapshot.data_vars)
            # Center-detection using chosen method (e.g., look for max(vort850), look for min(slp)
            if center_param == 'vort850':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).max().values, temp_snapshot[center_param], np.nan)
            elif center_param == 'slp':
                center_detector = xr.where(abs(temp_snapshot[center_param]) == abs(temp_snapshot[center_param]).min().values, temp_snapshot[center_param], np.nan)

            # Get extrema for chosen method from GCM output
            center_detector = center_detector.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] GCM-derived center: {0}, {1}'.format(center_detector['grid_xt'].values, center_detector['grid_yt'].values))
            # Distill to single float values
            center_lon, center_lat = center_detector.grid_xt.values[0], center_detector.grid_yt.values[0]
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] Storm ID: {0}, coordinates: {1}, {2}'.format(storm_id, center_lat, center_lon))
            # Get snapshot from original data
            temp_snapshot = snapshot.sel(time=timestamp, 
                                         grid_xt=np.arange(center_lon-extent, center_lon+extent), 
                                         grid_yt=np.arange(center_lat-extent, center_lat+extent), 
                                         method='nearest').drop_duplicates(['grid_xt'])
            temp_snapshot['center_lon'], temp_snapshot['center_lat'] = center_lon, center_lat
            clip_container.append(temp_snapshot)
            if diagnostic:
                print('[individual_tc_storage.py, retrieve_model_TCs()] clipped snapshot timestamps: {0}'.format(temp_snapshot.time.values))
        clip_container = xr.concat(clip_container, dim='time').drop_duplicates(dim='time')
        storms_xr[storm_id] = clip_container
        output_ids.append(storm_id)
        if diagnostic:
            print('[individual_tc_storage.py, retrieve_model_TCs()] append success for storm ID {0} at times {1}'.format(storm_id, clip_container.time.values))
        storm_counter += 1
        
        print('[individual_tc_storage.py, retrieve_model_TCs()] Number of storms to process: {0}, currently on {1}'.format(num_storms, storm_counter))
        
        # Limit number of storms read if a limit is given
        if num_storms and (storm_counter >= num_storms):
            return storms_xr, output_ids
        else:
            continue
        
    return storms_xr, output_ids

def vertical_profile_selection(storms, model, experiment, diagnostic=False):

    '''
    Arguments:
    - storms (list)
    '''
    start = time.time()
    
    # Vertical data before year 157 of HIRAM runs is output daily, 8xdaily after
    output_type = 'atmos_8xdaily' if '8xdaily' in model else 'atmos_daily'
    # Determine output frequency
    output_daily = True if output_type == 'atmos_daily' else False
    # Initialize container dictionary
    container = {}
    
    if diagnostic:
        print('[individual_tc_storage.py, vertical_profile_selection()] GCM output frequency: {0}'.format(output_type))

    # Iterate over each storm xarray Dataset
    for storm in storms.values():
        print('[individual_tc_storage.py, vertical_profile_selection()] Processing vertical data for storm ID: {0}'.format(storm['storm_id']))
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
            if diagnostic:
                print('[individual_tc_storage.py, vertical_profile_selection()] Processing timestamp: {0} at {1} frequency'.format(timestamp, output_type))
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
            if diagnostic:
                print('[individual_tc_storage.py, vertical_profile_selection()] Vertical timestamp out: {0}'.format(timestamp))
            # Get path name corresponding to model and experiment
            model_dir = utilities.directories(model, experiment, data_type='model_output')
            # Pull full vertical output and select corresponding day
            vertical_profile = xr.open_dataset('{0}/{1:04d}0101.{2}.nc'.format(model_dir, timestamp.year, output_type))
            # Filter by time and spatial domain
            try:
                if diagnostic:
                    print('[individual_tc_storage.py, vertical_profile_selection()] Vertical extraction at {0} for filename {1}/{2:04d}0101.{3}.nc'.format(timestamp, model_dir, timestamp.year, output_type))
                vertical_profile = vertical_profile.sel(time=timestamp, 
                                                        grid_xt=slice(min(grid_xt_bounds), max(grid_xt_bounds)), 
                                                        grid_yt=slice(min(grid_yt_bounds), max(grid_yt_bounds))).load()
            except:
                print('[individual_tc_storage.py, vertical_profile_selection()] Unable to retrieve vertical data for {0}'.format(snapshot['storm_id']))
                continue
            # Assign storm-specific snapshot values
            storm_params = ['center_lon', 'center_lat', 'max_wind', 'min_slp', 'core_temp', 'speed', 'heading', 'storm_id']
            # Assign fields to output
            for storm_param in storm_params:
                vertical_profile[storm_param] = snapshot[storm_param]

            output.append(vertical_profile)
            print('[individual_tc_storage.py, vertical_profile_selection()] Successful processing of vertical timestamp at {0}.'.format(timestamp))
        # The iteration is done to combat datasets with month 1 and day 1, which grab the whole year.
        output = [out.isel(time=0) if 'time' in out.dims else out for out in output]
        # Concatenate all timestamps and ensure each output only has one time index. 
        print('[individual_tc_storage.py, vertical_profile_selection()] Pre-concatenation data check for storm {0}: {1}'.format(storm_id, output))
        container[storm_id] = xr.concat(output, dim='time')
    return container

def main(model, experiments, storm_type, year_range, num_storms, storage=False, override=False, domain='atmos'):

    data = {model: {}}
    
    # Iterate over experiments to load and process data
    for experiment_num, experiment in enumerate(experiments):
        # data[model][experiment] = {}
        # Define paths to model- and experiment-specific data.
        model_dir, track_dir = utilities.directories(model, experiment, data_type='model_output'), utilities.directories(model, experiment, data_type='track_data')
    
        # For a given year:
        for year in year_range:
            data[model][experiment] = {}
            years = [year, year+1]
            print('\n-----------------------------------------------------------------------')
            print('Processing year: {0}'.format(year))
            # Get model data for future use in TC processing
            if domain == 'atmos':
                if model == 'HIRAM-8xdaily':
                    output_type = 'atmos_8xdaily'
                else:
                    output_type = 'atmos_daily'
            else:
                output_type = 'ocean_daily'
                
            print('serial')
            print(years)
            
            model_output = retrieve_model_data(model, model_dir, year_range=years, output_type=output_type)
            # Get tracked TCs from Lucas Harris' TC tracker
            track_output, storm_track_output = utilities.retrieve_tracked_TCs(model, experiment, storm_type, years, config='individual_tc_storage')
            # Check if any storms are found - if not, go to next loop
            if track_output is None or storm_track_output is None:
                continue
            # Pair model data to tracked TC for entire storm lifetimeprint
            storm_model_output, storm_ids = retrieve_model_TCs(model_dir, experiment, (min(years), max(years)), 
                                                               track_output, model_output, output_type=output_type, num_storms=num_storms, override=override)
            # Pair vertical model data to tracked TC for entire storm lifetime
            vertical_storm_output = vertical_profile_selection(storm_model_output, model, experiment)
            
            for storm_id in storm_ids:
                data[model][experiment][storm_id] = {}
                data[model][experiment][storm_id]['track_output'] = storm_track_output[storm_id]
                data[model][experiment][storm_id]['tc_model_output'] = storm_model_output[storm_id]
                data[model][experiment][storm_id]['tc_vertical_output'] = vertical_storm_output[storm_id]
                data[model][experiment][storm_id]['run_name'] = model_dir

                # Assign intensity bins to each timestamp
                data[model][experiment][storm_id] = utilities.intensity_binning(mode='model_output', data=data[model][experiment][storm_id], intensity_metric='max_wind')
                # # Process data and obtain derived fields
                # data[model][experiment][storm_id] = tc_processing.derived_fields(data[model][experiment][storm_id])
        
                # Storage for each individual storm
                if storage:
                    # Define directory name and file name
                    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
                    # Define filename using max wind and min SLP for future binning
                    max_wind, min_slp = storm_track_output[storm_id]['max_wind'].max(), storm_track_output[storm_id]['min_slp'].min()
                    storage_filename = 'TC-{0}-{1}-{2}-{3}-max_wind-{4:0.0f}-min_slp-{5:0.0f}.pkl'.format(model, experiment, storm_type, storm_id, max_wind, min_slp)
                    storage_path = os.path.join(storage_dirname, storage_filename)
                    # If file doesn't exist, save
                    if not os.path.isfile(os.path.join(storage_dirname, storage_filename)):
                        with open(os.path.join(storage_dirname, storage_filename), 'wb') as f:
                            print('Saving to: ', os.path.join(storage_dirname, storage_filename))
                            pickle.dump(data[model][experiment][storm_id], f)

    return data

def main_parallel(model, experiment, storm_type, year, num_storms, storage=False, override=False, domain='atmos', single_storm_id=None):

    data = {model: {}}
    
    # Iterate over experiments to load and process data
    # Define paths to model- and experiment-specific data.
    model_dir, track_dir = utilities.directories(model, experiment, data_type='model_output'), utilities.directories(model, experiment, data_type='track_data')

    # For a given year:
    data[model][experiment] = {}
    years = [year, year+1]
    print('\n-----------------------------------------------------------------------')
    print('Processing year: {0}, override status: {1}'.format(year, override))
    # Get model data for future use in TC processing
    if domain == 'atmos':
        if model == 'HIRAM-8xdaily':
            output_type = 'atmos_8xdaily'
        else:
            output_type = 'atmos_daily'
    else:
        output_type = 'ocean_daily'
        
    model_output = retrieve_model_data(model, model_dir, year_range=years, output_type=output_type)
    # Get tracked TCs from Lucas Harris' TC tracker
    intensity_metric, intensity_threshold = ('min_slp', 990) if storm_type == 'TS' else (None, None)
    track_output, storm_track_output = utilities.retrieve_tracked_TCs(model, experiment, storm_type, years, 
                                                                      intensity_metric=intensity_metric, intensity_threshold=intensity_threshold, config='individual_tc_storage')

    print(track_output)

    # Check if any storms are found - if not, go to next loop
    if track_output is None or storm_track_output is None:
        return None
    
    # Pair model data to tracked TC for entire storm lifetimeprint
    storm_model_output, storm_ids = retrieve_model_TCs(model_dir, experiment, (min(years), max(years)),
                                                       track_output, model_output, output_type='atmos_4xdaily', maximum_latitude=35, num_storms=num_storms, override=override, single_storm_id=single_storm_id)
    
    # Pair vertical model data to tracked TC for entire storm lifetime
    vertical_storm_output = vertical_profile_selection(storm_model_output, model, experiment)
     
    for storm_id in storm_ids:
        data[model][experiment][storm_id] = {}
        data[model][experiment][storm_id]['track_output'] = storm_track_output[storm_id]
        data[model][experiment][storm_id]['tc_model_output'] = storm_model_output[storm_id].load()
        data[model][experiment][storm_id]['tc_vertical_output'] = vertical_storm_output[storm_id].load()
        data[model][experiment][storm_id]['run_name'] = model_dir
        
        # Assign intensity bins to each timestamp
        data[model][experiment][storm_id] = utilities.intensity_binning(mode='model_output', data=data[model][experiment][storm_id], intensity_metric='max_wind')
        # Process data and obtain derived fields
        data[model][experiment][storm_id] = tc_processing.derived_fields(data[model][experiment][storm_id])

        # Storage for each individual storm
        if storage:
            # Define directory name and file name
            storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
            # Define filename using max wind and min SLP for future binning
            max_wind, min_slp = storm_track_output[storm_id]['max_wind'].max(), storm_track_output[storm_id]['min_slp'].min()
            storage_filename = 'TC-{0}-{1}-{2}-{3}-max_wind-{4:0.0f}-min_slp-{5:0.0f}.pkl'.format(model, experiment, storm_type, storm_id, max_wind, min_slp)
            storage_path = os.path.join(storage_dirname, storage_filename)
            # If file doesn't exist, save
            if not os.path.isfile(os.path.join(storage_dirname, storage_filename)) or override:
                with open(os.path.join(storage_dirname, storage_filename), 'wb') as f:
                    print('Saving to: ', os.path.join(storage_dirname, storage_filename))
                    pickle.dump(data[model][experiment][storm_id], f)

    # return track_output

if __name__ == '__main__':
    start = time.time()
    # Use range(start, stop) for a range of years between 'start' and 'stop', and a list [start, stop] for specific years.
    year_range = range(101, 102)
    # Enable parallel processing
    parallel_implementation = True

    models = ['AM2.5']
    experiments = ['CTL1990s_tiger3']
    storm_type = 'C15w'
    num_storms = 5
    domain = 'atmos'
    single_storm_id = None
        
    if parallel_implementation:
        # Prepare list of model inputs for parallel processing
        inputs = []
        # Assemble the input lists for Pool.starmap
        for model in models:
            print(model)
            for experiment in experiments:
                for year in year_range:
                    inputs.append([model, experiment, storm_type, year, num_storms, True, True, domain, single_storm_id])
        [print(input_, '\n') for input_ in inputs]
        # Use N processors to process storms in parallel
        max_procs = 16
        num_procs = num_storms if num_storms < max_procs else max_procs
        with multiprocessing.get_context("spawn").Pool(num_procs) as p:
            track_data = p.starmap(main_parallel, inputs)
            
    else:
        for model in models:
            data = main(model, experiments=experiments, storm_type=storm_type, year_range=year_range, 
                        num_storms=num_storms, storage=True, override=True, domain=domain)
            print('Elapsed total runtime: {0:.3f}s'.format(time.time() - start))
            
    print('Runtime: {0:.2f} s'.format(time.time() - start))