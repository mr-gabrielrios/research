'''

Storm tracking methodology

__Storm datasets__: 
- `GCM output`: output from global climate model
- `track data`: output from GFDL QuickTracks model that is run on `GCM output`

__For each storm__:
1. Find a candidate storm from `track data`
2. Get candidate storm timestamps
3. Get candidate storm coordinates
4. For each candidate storm timestamp, find corresponding `GCM output` file
5. For each candidate storm timestamp, use storm coordinates to trim time of `GCM output` file in 
6. For each candidate storm timestamp, use storm coordinates to trim spatial extent of `GCM output` file
7. Append information from `track data` to netCDF object containing GCM output
8. Save xArray Dataset to netCDF file

'''

# User input package
import argparse
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import functools, os, pickle, collections, sys
# Data structure packages
import pandas as pd, xarray as xr, nc_time_axis
xr.set_options(keep_attrs=True)
# Visualization tools
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt
# Local imports
import accessor, composite, composite_snapshots, derived, utilities, socket, visualization, tc_analysis, tc_processing, track_TCs

import multiprocessing
multiprocessing.set_start_method("spawn", force=True)
from multiprocessing import Pool

def access_storm_tracks(model_name: str,
                        experiment_name: str,
                        year_range: tuple[int, int]) -> pd.DataFrame:

    ''' Obtain track data for a given model, experiment, and year range. '''

    track_data = tc_analysis.load_TC_tracks(model_name, experiment_name, year_range)[model_name][experiment_name]['raw']

    return track_data

def interpolate_storm_tracks(track_data: pd.DataFrame,
                             frequency: str='6H'):
    
    ''' Interpolate TC tracks to a given temporal frequency. '''
    
    # Initialize container to hold all interpolate DataFrames
    container = [] 
    # Iterate over all TCs in the input track data dataset
    for TC_ID, TC in track_data.groupby('storm_id'):
        # Convert the index to a DatetimeIndex type, resample to the given frequency, interpolate numeric values, and forward fill strings
        temp = TC.set_index('time').resample(frequency).interpolate(method='linear', limit_direction='forward', axis=0).ffill()
        # Append the iterand interpolated data to the container
        container.append(temp.reset_index())

    interpolated_track_data = pd.concat(container).sort_values('time')
    
    return interpolated_track_data

def test_single_storm(storm_track_data: pd.DataFrame):

    ''' Provide tests to ensure tracked storm is legitimate. '''
    criterion_flag = True

    # Look for duplicate storms
    # assert len(storm_track_data['duration'].unique()) == 1, 'It looks like two tracked storms are using the same storm ID.'
    if len(storm_track_data['duration'].unique()) != 0:
        criterion_flag = False

    # Constrain the distance the storm can move between two timestamps
    threshold_delta_longitude = 10 # units of degrees
    threshold_delta_latitude = 10 # units of degrees
    
    delta_longitude = storm_track_data['center_lon'].diff().dropna().abs()
    # assert max(delta_longitude) < threshold_delta_longitude, f'Longitude threshold exceeded, value = {max(delta_longitude):.2f} degrees.'
    if max(delta_longitude) >= threshold_delta_longitude:
        criterion_flag = False
        
    delta_latitude = storm_track_data['center_lat'].diff().dropna().abs()
    # assert max(delta_latitude) < threshold_delta_latitude, f'Longitude threshold exceeded, value = {max(delta_latitude):.2f} degrees.'
    if max(delta_latitude) >= threshold_delta_latitude:
        criterion_flag = False
        
    # Ensure the necessary columns are in the storm DataFrame
    # assert ('cftime' in storm_track_data.columns), f'cftime data is not in the track DataFrame.'
    if 'cftime' not in storm_track_data.columns:
        criterion_flag = False

    return 

def intensity_filter(track_data: pd.DataFrame,
                     intensity_parameter: str='min_slp',
                     intensity_range: tuple[int|float, int|float]=(0, np.inf)) -> pd.DataFrame:

    ''' 
    Method to filter QuickTracks output track data to find storms that fit a given intensity parameter and range.
    Note: storms are only filtered by maximum winds ('max_wind') and minimum sea-level pressure ('min_slp').
    '''

    assert intensity_parameter in ['min_slp', 'max_wind'], f'Parameter {intensity_parameter} not recognized, please use `max_wind` or `min_slp`.'

    # Obtain maximum intensities for each storm
    maximum_intensity_by_storm = track_data.groupby('storm_id')[intensity_parameter].min() if intensity_parameter == 'min_slp' else track_data.groupby('storm_id')[intensity_parameter].max()
    # Obtain the storm IDs that contain maximum intensities within the given intensity band
    threshold_storm_IDs = maximum_intensity_by_storm[((maximum_intensity_by_storm >= min(intensity_range)) & 
                                                      (maximum_intensity_by_storm < max(intensity_range)))].index.values
    # Obtain track entries matching the storms tht meet the threshold
    threshold_storms = track_data.loc[track_data['storm_id'].isin(threshold_storm_IDs)]

    return threshold_storms

def latitude_filter(track_data: pd.DataFrame,
                    intensity_parameter: str='min_slp',
                    latitude_range: tuple[int|float, int|float]=(-40, 40)) -> pd.DataFrame:

    ''' Filter out TCs with LMIs occurring outside a given latitude band to prevent capturing ETCs. '''

    # Initialize a container list for storm IDs that meet the criteria
    storm_IDs = []

    # Iterate over each storm to determine if its LMI occurs within the given latitude band
    for storm_ID, storm_track_data in track_data.groupby('storm_id'):
        # Get maximum storm intensity
        maximum_intensity = storm_track_data[intensity_parameter].max() if intensity_parameter == 'max_wind' else storm_track_data[intensity_parameter].min()
        # Get latitude corresponding to LMI occurrence.
        # Note the absolute and maximum functions. This prevents duplicate timestamps with identical maximum intensities of being loaded, gets the more poleward occurrence.
        latitude_of_maximum_intensity = abs(storm_track_data.loc[storm_track_data[intensity_parameter] == maximum_intensity]['center_lat']).max()
        # Append to list if the latitude of LMI is in the given range
        if min(latitude_range) <= latitude_of_maximum_intensity <= max(latitude_range):
            storm_IDs.append(storm_ID)

    threshold_track_data = track_data.loc[track_data['storm_id'].isin(storm_IDs)]

    return threshold_track_data

def pick_storm_IDs(track_data: pd.DataFrame,
                   number_of_storms: int,
                   TC_storage_dirname: str|None = None) -> list:

    ''' Method to obtain `number_of_storms` random storm IDs from a provided track dataset. '''
    
    # Define the directory where TCs are stored
    TC_storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs' if not TC_storage_dirname else TC_storage_dirname
    # Load filenames for future searchinf
    TC_filenames = [filename for filename in os.listdir(TC_storage_dirname) if filename.endswith('.nc')]
    
    # Get list of unique storm IDs
    unique_storm_IDs = list(track_data['storm_id'].unique())
    # Make sure the number of requested storms is less than the number of unique IDs; 
    number_of_storms = len(unique_storm_IDs) if (number_of_storms > len(unique_storm_IDs)) or (number_of_storms < 0) else number_of_storms
    
    # Remove storm IDs that already have a corresponding TC saved
    # Search the storage directory to see if a TC with the iterand storm ID has already been saved
    existing_storm_IDs = [existing_storm_ID for existing_storm_ID in TC_filenames if
                          existing_storm_ID in TC_filenames]
    # Remove the storm IDs that already exist from the list of unique storm IDs
    [unique_storm_IDs.remove(existing_storm_ID) for existing_storm_ID in existing_storm_IDs if
     existing_storm_ID in unique_storm_IDs]
    
    # If not, make the number of storms equal to `unique_storm_IDs`
    # Get indices for `number_of_storms` random storm IDs
    storm_ID_indices = np.random.choice(range(len(unique_storm_IDs)), size=number_of_storms, replace=False)
    # Get randomized storm IDs
    storm_IDs = track_data['storm_id'].unique()[storm_ID_indices]

    return storm_IDs

def pick_storm(track_data: pd.DataFrame,
               selection_method: str='random',
               storm_ID: str|None=None):

    ''' Pick a single storm from the track data. '''

    if selection_method == 'random':
        storm_index = np.random.randint(0, len(track_data['storm_id'].unique())) # choose random storm
        storm_ID = track_data['storm_id'].unique()[storm_index]
    elif selection_method == 'storm_number' and storm_ID:
        assert isinstance(storm_ID, str), 'Storm ID must be a string.'
        assert storm_ID in track_data['storm_id'].values, f'Storm ID {storm_ID} not found in track dataset.'
    else:
        print('Please provide a storm ID or set `selection_method` to `random`. Exiting.')
        sys.exit()

    # Pull data for a specific storm and sort values by time
    storm_track_data = track_data.loc[track_data['storm_id'] == storm_ID].sort_values('cftime')
    # Check for storm quality
    test_single_storm(storm_track_data)
    
    return storm_track_data

def storm_GCM_calendar_alignment(storm_timestamps: list[cftime.datetime], 
                                 GCM_timestamps: list[cftime.datetime]) -> list[cftime.datetime]:

    ''' 
    Ensure GFDL QuickTracks track data and GCM output data have same calendars. 
    GCM data takes precedence since GCM timestamps are less mutable than QuickTracks timestamps. 
    '''
    
    # Scrape calendar type from variable type
    get_calendar_type = lambda timestamp: (utilities.cftime_calendar_type(timestamp), timestamp.has_year_zero)
    # Function to convert timestamp formats for a given timestamp and calendar
    timestamp_conversion = lambda t, calendar, has_year_zero: cftime.datetime(year=t.year,
                                                                              month=t.month,
                                                                              day=t.day,
                                                                              hour=t.hour,
                                                                              calendar=calendar,
                                                                              has_year_zero=has_year_zero)
    
    # Iterate through timestamps to get types for storm and GCM data.
    # Assume all entries have the same data type.
    storm_timestamp_type = str(type(storm_timestamps[0]))
    GCM_timestamp_type = str(type(GCM_timestamps[0]))
    # print(f'Track data timestamp type: {storm_timestamp_type}; GCM timestamp type: {GCM_timestamp_type}')
    
    # Reformat storm timestamps to the GCM timestamp format, if different
    # Assume all entries have the same data type.
    GCM_calendar_type, GCM_has_year_zero = get_calendar_type(GCM_timestamps[0])
    # print(GCM_calendar_type, GCM_has_year_zero)
    storm_timestamps_reformatted = [timestamp_conversion(t, GCM_calendar_type, GCM_has_year_zero) 
                                    for index, t in enumerate(storm_timestamps)]
    
    return storm_timestamps_reformatted

def get_storm_GCM_data(model_name: str,
                       experiment_name: str,
                       storm_track_timestamps,
                       GCM_data_type: str='atmos_4xdaily') -> list:
    
    ''' For each candidate storm timestamp, find corresponding `GCM output` file. '''

    # Define top-level directory where GCM output is kept
    gcm_container_dirname = '/tigress/GEOCLIM/gr7610/MODEL_OUT' 
    # Define the configuration-specific directory
    gcm_dirname = os.path.join(gcm_container_dirname, model_name, experiment_name, 'POSTP') 
    # Obtain filenames in the configuration-specific directory for the chosen data type
    gcm_pathnames = [os.path.join(gcm_dirname, gcm_filename) for gcm_filename in os.listdir(gcm_dirname)
                     if gcm_filename.endswith('.nc') and 
                     GCM_data_type in gcm_filename] 
    # Ensure files are found in the directory. If not, exit.
    if len(gcm_pathnames) > 0:
        # Filter pathname list for paths containing years relevant to storm
        # Note that a minimum and maximum year is obtained to handle storms that run into the following year
        storm_track_year_min, storm_track_year_max = min(storm_track_timestamps).year, max(storm_track_timestamps).year
        # Filter by year, ends-inclusive. 
        # Note that the year is obtained crudely, assuming that GCM model output is in YYYYMMDD.{gcm_data_type}.nc format.
        storm_gcm_pathnames = [pathname for pathname in gcm_pathnames 
                               if int(os.path.basename(pathname).split('.')[0][0:4]) >= storm_track_year_min 
                               and int(os.path.basename(pathname).split('.')[0][0:4]) <= storm_track_year_max] 
        # Check on pathnames length
        assert (len(storm_gcm_pathnames) > 0), f'No files found in {gcm_dirname} for model {model_name} and experiment {experiment_name} for the years {storm_track_year_min}-{storm_track_year_max}. Please check the directory and retry.'
    
    else:
        print(f'No files found in {gcm_dirname} for model {model_name} and experiment {experiment_name}. Please check the directory and retry.')
        sys.exit()

    ''' Timestamp alignment. '''
    # Get storm GCM timestamps for dataset calendar alignment.
    storm_gcm_data_timestamps = xr.open_mfdataset(storm_gcm_pathnames).time.values
    # Check if GCM data timestamps are non-noleap. 
    # If so, assume Julian and convert to match the QuickTracks data with GCM data calendar conventions.
    storm_track_timestamps_reformatted = storm_GCM_calendar_alignment(storm_track_timestamps.values,
                                                                      storm_gcm_data_timestamps)

    return storm_gcm_pathnames, storm_track_timestamps_reformatted

def get_storm_coordinates(storm_track_data: pd.DataFrame,
                          storm_track_timestamps) -> dict:

    ''' Get candidate storm coordinates. '''

    # Note the variable name convention ('storm_track' instead of just 'storm').
    # This is intended to distinguish tracker coordinates from actual storm centered coordinates, which may be different in GCM output.
    
    # Initialize coordinates dictionary, which will save a longitude and latitude as values for a timestamp key
    storm_track_coordinates = {}
    # Iterate over timestamps to pair a timestamp to corresponding coordinates
    for storm_track_timestamp in storm_track_timestamps:
        # Obtain longitude and latitudes for each timestamp
        storm_track_longitude = storm_track_data.loc[storm_track_data['cftime'] == storm_track_timestamp]['center_lon'].item()
        storm_track_latitude = storm_track_data.loc[storm_track_data['cftime'] == storm_track_timestamp]['center_lat'].item()
        storm_track_coordinates[storm_track_timestamp] = {'lon': storm_track_longitude, 'lat': storm_track_latitude}

    return storm_track_coordinates

def load_GCM_data(storm_gcm_pathnames, 
                  storm_track_timestamps,
                  storm_track_coordinates) -> (xr.Dataset, ):

    ''' Load and trim the data in time and space. '''
    
    # Load the data.
    storm_gcm_data = xr.open_mfdataset(storm_gcm_pathnames)
    # # Check if GCM data timestamps are non-noleap. If so, assume Julian and convert.
    # storm_track_timestamps_reformatted = storm_GCM_calendar_alignment(storm_track_timestamps.values,
    #                                                                   storm_gcm_data.time.values)
    # Obtain timestamps shared by GCM data and track timestamps.
    storm_gcm_timestamps = list(set(storm_track_timestamps) & set(storm_gcm_data.time.values))
    # Trim storm track coordinates to match the shared timestamps
    storm_track_coordinates = {storm_gcm_timestamp: storm_track_coordinates[storm_gcm_timestamp] for storm_gcm_timestamp in storm_gcm_timestamps}
    # Trim GCM output data in time.
    storm_gcm_data = storm_gcm_data.sel(time=storm_gcm_timestamps)

    return storm_gcm_data, storm_gcm_timestamps

def trim_GCM_data(storm_gcm_data: xr.Dataset,
                  storm_gcm_timestamps,
                  storm_track_coordinates, 
                  storm_gcm_window_size: int | float=12,
                  diagnostic: bool=False) -> xr.Dataset:
    
    ''' For each candidate storm timestamp, use storm coordinates to trim spatial extent of `GCM output` file. '''

    # Start a timer for performance profiling
    start_time = time.time()

    # Initialize a container to hold GCM output connected to each storm timestamp and the corresponding spatial extent
    storm_gcm_container = {}
    # Generate trimming window extents for each timestamp.
    # Window extents are defined as: 
    # 'grid_xt' = (longitude - window_extent, longitude + window_extent), 
    # 'grid_yt' = (latitude - window_extent, latitude + window_extent)
    storm_gcm_window_extent = {}
    for storm_gcm_timestamp in storm_gcm_timestamps:
        storm_gcm_window_extent[storm_gcm_timestamp] = {}
        # Assign zonal window
        storm_gcm_window_extent[storm_gcm_timestamp]['grid_xt'] = slice(storm_track_coordinates[storm_gcm_timestamp]['lon'] - storm_gcm_window_size,
                                                                        storm_track_coordinates[storm_gcm_timestamp]['lon'] + storm_gcm_window_size)
        # Assign meridional window
        storm_gcm_window_extent[storm_gcm_timestamp]['grid_yt'] = slice(storm_track_coordinates[storm_gcm_timestamp]['lat'] - storm_gcm_window_size,
                                                                        storm_track_coordinates[storm_gcm_timestamp]['lat'] + storm_gcm_window_size)
        # Extract GCM data for the given timestamp and spatial extent
        storm_gcm_container[storm_gcm_timestamp] = storm_gcm_data.sel(time=storm_gcm_timestamp,
                                                                      grid_xt=storm_gcm_window_extent[storm_gcm_timestamp]['grid_xt'],
                                                                      grid_yt=storm_gcm_window_extent[storm_gcm_timestamp]['grid_yt'])
    
    # Concatenate all GCM output data corresponding to storm into a single xArray Dataset
    storm_gcm_data = xr.concat(storm_gcm_container.values(), dim='time').sortby('time')
    
    if diagnostic:
        print(f'[trim_GCM_data()] Elapsed time for storm: {(time.time() - start_time):.2f} s.')

    return storm_gcm_data

def join_track_GCM_data(storm_track_data: pd.DataFrame,
                        storm_gcm_data: xr.Dataset,
                        storm_time_variable: str='cftime'):

    ''' Append information from `track data` to netCDF object containing GCM output. '''
    
    # Filter storm track data that has matching timestamps with xArray Dataset `storm_gcm_data`
    # These should already match, but this is for posterity
    storm_track_data_gcm = storm_track_data.loc[storm_track_data[storm_time_variable].isin(storm_gcm_data.time.values)]
    # Define variables to append to netCDF
    storm_track_gcm_vars = ['center_lon', 'center_lat', 'min_slp', 'max_wind', 'storm_id']
    # Test to make sure that the variables requested are all in the track data DataFrame
    assert set(storm_track_gcm_vars) <= set(storm_track_data_gcm.columns), 'Not all requested track data columns are available in the given track dataset.'
    
    # Add the data to the xArray Dataset `storm_gcm_data`
    for storm_track_gcm_var in storm_track_gcm_vars:
        # Handle storm ID as an attribute since it's time-invariant
        if storm_track_gcm_var == 'storm_id':
            # Get unique storm ID value
            storm_id = storm_track_data_gcm['storm_id'].unique().item()
            # Add to xArray Dataset attributes
            storm_gcm_data.attrs['storm_id'] = storm_track_data_gcm['storm_id'].unique().item()
        # Otherwise, append track data along the time axis
        else:
            storm_gcm_data[storm_track_gcm_var] = xr.DataArray(data=storm_track_data_gcm[storm_track_gcm_var].values,
                                                               dims=['time'],
                                                               coords={'time': ('time', storm_gcm_data.time.values)})

    return storm_gcm_data

def get_storm_basin_name(storm: xr.Dataset) -> str:

    ''' 
    Method to obtain the basin a TC belongs to. Assume the TC's first timestamp is representative of its basin. 
    
    Algorithm: for a given TC, use the coordinates at the first timestamp to determine which basin the TC is in.
    '''

    # Ensure storm is sorted by time
    storm = storm.sortby('time')
    # Get coordinates at first timestamp
    lon, lat = storm.isel(time=0)['center_lon'].item(), storm.isel(time=0)['center_lat'].item()
    # Pull archived basin masks for GCM output dats with 0.5 degree spatial reslution
    basin_masks = xr.open_dataset('/projects/GEOCLIM/gr7610/tools/basin_mask.nc')
    # Iterate through basin names until a match is found
    basin_name = [basin_name for basin_name in basin_masks.keys() if
                  basin_name not in ['global', 'IPWP', 'ENSO'] and 
                  basin_masks[basin_name].sel(grid_xt=lon, method='nearest').sel(grid_yt=lat, method='nearest').item() == 1]
    # Ensure the length of the resulting list only has 0 or 1 elements
    assert len(basin_name) < 2, f'Storm {storm.attrs['storm_id']} found in basins {basin_name}, investigate why this is the case.'
    # Extract the basin name, or if none is found, assign as extratropical
    basin_name = 'ET' if len(basin_name) == 0 else basin_name[0] # handle storms that occur outside conventional basis as extratropical (ET)

    return basin_name

def save_storm_netcdf(storm_gcm_data: xr.Dataset,
                      model_name: str,
                      experiment_name: str,
                      GCM_data_type: str,
                      storage_dirname: str|None=None,
                      number_of_lag_days: int=0,
                      overwrite: bool=False):
    
    ''' Save xArray Dataset to netCDF file. '''

    # Define reanalysis dataset names for custom handling
    reanalysis_dataset_names = ['ERA5']
    # Define storage directory
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs' if storage_dirname is None else storage_dirname
    
    # Obtain parameters for filename construction
    storm_ID = storm_gcm_data.attrs['storm_id']
    max_wind = f'{np.round(storm_gcm_data['max_wind'].max()):.0f}' # round to nearest integer for brevity
    min_slp = f'{np.round(storm_gcm_data['min_slp'].min()):.0f}' # round to nearest integer for brevity

    # Obtain the basin name for the storm filename
    storm_basin = get_storm_basin_name(storm_gcm_data)
    # Construct lag day substring
    number_of_lag_days_str = f'p{abs(number_of_lag_days)}' if number_of_lag_days > 0 else f'm{abs(number_of_lag_days)}'
    lag_day_substring = f'.lag_days-{number_of_lag_days_str}' if number_of_lag_days != 0 else f'.lag_days-0'
    # Build storm filename
    storm_filename = f'TC.model-{model_name}.experiment-{experiment_name}.storm_ID-{storm_ID}.max_wind-{max_wind}.min_slp-{min_slp}.data_type-{GCM_data_type}.basin-{storm_basin}{lag_day_substring}.nc'
    storm_pathname = os.path.join(storage_dirname, storm_filename)
    
    # Load the data into memory before saving to ensure output is fully there
    print(f'[save_storm_netcdf] Loading data for {storm_filename}...')

    # Profile loading time
    start_time = time.time()
    storm_gcm_data.load()
    print(f'Elapsed loading time for {storm_filename}: {(time.time() - start_time):.2f} s')
    
    # Print output file size as a diagnostic
    print(f'File size for {storm_filename}: {(storm_gcm_data.nbytes / 1e6):.2f} MB\n')
    
    # Check if file exists
    does_file_exist = os.path.isfile(storm_pathname)
    if (does_file_exist and overwrite) or not does_file_exist:
        # Save the data
        storm_gcm_data.to_netcdf(storm_pathname)
    
def storm_generator(model_name: str,
                    experiment_name: str,
                    track_data: pd.DataFrame,
                    storage_dirname: str|None,
                    GCM_data_type: str,
                    number_of_lag_days: int,
                    storm_ID: str):

    ''' Method to perform all steps related to binding corresponding GFDL QuickTracks and GCM model output together for a given TC. '''

    print(f'[storm_generator] Processing storm ID {storm_ID}...')
    
    # 3. Find a candidate storm from the track data
    storm_track_data = pick_storm(track_data, selection_method='storm_number', storm_ID=storm_ID)
    if storm_track_data is not None:
        # 4. Get candidate storm timestamps
        storm_track_timestamps = storm_track_data['cftime']
        # 5. For each candidate storm timestamp, find corresponding `GCM output` file
        storm_track_timestamps = storm_track_timestamps + pd.Timedelta(number_of_lag_days, 'D') # adjust days based on lag input
        storm_gcm_pathnames, storm_track_timestamps = get_storm_GCM_data(model_name, experiment_name, storm_track_timestamps, GCM_data_type=GCM_data_type)
        # 5a. Correct GFDL QuickTracks cftime timestamp format to match GCM output format
        storm_track_data['cftime'] = storm_track_timestamps
        # 6. Get candidate storm coordinates
        storm_track_coordinates = get_storm_coordinates(storm_track_data, storm_track_timestamps)
        # 7. For each candidate storm timestamp, load GCM data and use storm timestamps to trim time of `GCM output` file
        storm_gcm_data, storm_gcm_timestamps = load_GCM_data(storm_gcm_pathnames, storm_track_timestamps, storm_track_coordinates)
        # 8. For each candidate storm timestamp, use storm coordinates to trim spatial extent of `GCM output` file
        storm_gcm_data = trim_GCM_data(storm_gcm_data, storm_gcm_timestamps, storm_track_coordinates)
        # 9. Append information from `track data` to netCDF object containing GCM output
        storm_gcm_data = join_track_GCM_data(storm_track_data, storm_gcm_data)
        # 10. Save xArray Dataset to netCDF file
        save_storm_netcdf(storm_gcm_data, model_name, experiment_name, GCM_data_type=GCM_data_type, storage_dirname=storage_dirname, number_of_lag_days=number_of_lag_days)
    
def main(model_name: str, 
         experiment_name: str, 
         year_range: tuple[int, int],
         intensity_parameter: str,
         intensity_range: tuple[int|float, int|float]=(0, np.inf),
         latitude_range: tuple[int|float, int|float]=(-40, 40),
         number_of_storms: int=1,
         storage_dirname: str|None=None,
         GCM_data_type: str='atmos_4xdaily',
         number_of_lag_days: int|list[int]=0,
         storm_IDs: list[str] | None=None):

    # 0. Obtain track data for a given model, experiment, and year range
    track_data = access_storm_tracks(model_name, experiment_name, year_range)
    # 1a. Filter storms by intensity
    track_data = intensity_filter(track_data, intensity_parameter, intensity_range)
    # 1b. Filter storms by latitude
    track_data = latitude_filter(track_data, intensity_parameter, latitude_range)
    # 2. If no storm IDs are provided, obtain N randomized storm IDs from the filtered data, where 'N' is `number_of_storms`
    storm_IDs = pick_storm_IDs(track_data, number_of_storms) if storm_IDs is None else storm_IDs

    # 3. If the number of lag days provided is an int, convert to a string
    if isinstance(number_of_lag_days, int):
        number_of_lag_days = [number_of_lag_days]
    else:
        assert isinstance(number_of_lag_days, list), 'number_of_lag_days must be an integer or list of integers.'

    ''' Offload TC-specific data generation onto parallel processes. '''
    # Maximum number of processors for computation
    max_number_procs = 8
    # Specify number of processors to use
    number_procs = len(storm_IDs) if len(storm_IDs) < max_number_procs else max_number_procs
    
    for lag_day in number_of_lag_days:
        print(f'[main] Working on lag day {lag_day}...')
        # Define partial function to allow for using Pool.map since `track_data` is equivalent for all subprocesses
        preloaded_storm_generator = functools.partial(storm_generator, model_name, experiment_name, track_data, storage_dirname, GCM_data_type, lag_day)

        with Pool(processes=number_procs) as pool:
            pool.map(preloaded_storm_generator, storm_IDs)
            pool.close()
            
        del preloaded_storm_generator
        
    # Print completion message
    print('Processing of TCs completed.')

if __name__ == '__main__':
    
    ''' Initialize user input/output section. '''
    # Catch user inputs
    parser = argparse.ArgumentParser()
    
    # Required inputs
    parser.add_argument('--model_name', type=str, help="Name of the model to extract data from.", required=True)
    parser.add_argument('--experiment_name', type=str, help="Name of the experiment to extract data from.", required=True)
    parser.add_argument('--year_range', type=str, help="Range of years to extract data from. Provide input in `YYYY:YYYY` format.", required=True)
    # Optional inputs
    
    # Required inputs
    parser.add_argument('--intensity_parameter', nargs='?', default='min_slp', type=str, help="Intensity parameter for TC filtering. Either `min_slp` or `max_wind`.")
    parser.add_argument('--intensity_range', nargs='?', default='0:1020', type=str, help="Range of values over which to search for TCs for the given intensity range. Provide input in `X:Y` format.")
    parser.add_argument('--number_of_storms', nargs='?', default=1, type=int, help="Number of TCs to generate data for.")
    parser.add_argument('--GCM_data_type', nargs='?', default='atmos_4xdaily', type=str, help="Data type to extract files from. Usually `atmos_4xdaily` or `atmos_month`.")
    parser.add_argument('--number_of_lag_days', nargs='?', default=0, type=str, help="Number of days from TC occurrence to sample GCM data from. Can be positive or negative, where a positive value is days after TC passage.")
    parser.add_argument('--storage_dirname', nargs='?', default='/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs', type=str, help="Directory on the Tiger cluster to save data to.")
    parser.add_argument('--storm_IDs', nargs='?', default='None', type=str, help="List of storm IDs saved as strings, delimited by a colon. Defaults to None.")
    
    
    # Process inputs
    args = parser.parse_args()
    
    # Process format-dependent arguments
    args.year_range = tuple([int(year) for year in args.year_range.split(':')])
    args.intensity_range = tuple([int(bound) for bound in args.intensity_range.split(':')])
    args.storm_IDs = None if args.storm_IDs.lower() == 'none' else [storm_ID for storm_ID in args.storm_IDs.split(':')]
    args.number_of_lag_days = [int(year) for year in args.number_of_lag_days.split(':')] if ':' in args.number_of_lag_days else int(args.number_of_lag_days)
    # Print arguments for diagnostic
    print(f'------------------------------------------------------------------\n{args}\n------------------------------------------------------------------')
    
    main(model_name=args.model_name, 
         experiment_name=args.experiment_name,
         year_range=args.year_range, 
         intensity_parameter=args.intensity_parameter, 
         intensity_range=args.intensity_range, 
         number_of_storms=args.number_of_storms,
         number_of_lag_days=args.number_of_lag_days,
         GCM_data_type=args.GCM_data_type,
         storage_dirname=args.storage_dirname,
         storm_IDs=args.storm_IDs)
    ''' End user input/output section. '''
    
''' 
Example snippet

Parameters: model FLOR, experiment CTL1990s_FA_tiger3, years 1901 to 1905 with an intensity range of 980 to 1000 hPa:

### 

model_name = 'FLOR'
experiment_name = 'CTL1990s_FA_tiger3'
year_range = (1901, 1905)

intensity_parameter = 'min_slp'
intensity_range = (980, 1000)

main(model_name, experiment_name, year_range, intensity_parameter, intensity_range, number_of_storms=5)

###
'''