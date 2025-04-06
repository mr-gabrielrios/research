import pandas as pd
import os
import importlib
import functools
import time
from multiprocessing import Pool
import timezonefinder
import datetime
import pytz

import cftime
import nc_time_axis
import xarray as xr
import numpy as np
import matplotlib, matplotlib.pyplot as plt, matplotlib.patheffects as pe

import derived
import tc_analysis
import utilities
import visualization

importlib.reload(derived);
importlib.reload(tc_analysis);
importlib.reload(utilities);
importlib.reload(visualization);


def load_track_data(model_name: str,
                    experiment_name: str,
                    year_range: tuple[int, int]):

    # Handle ERA5 case for IBTrACS data
    model_name = 'IBTrACS' if model_name in ['CERES', 'ERA5', 'MERRA'] else model_name
    experiment_name = '' if model_name == 'IBTrACS' else experiment_name
    
    # Load track data for the model-experiment configuration
    track_dataset = tc_analysis.load_TC_tracks(model_name, 
                                               experiment_name, 
                                               year_range, 
                                               diagnostic=False)[model_name][experiment_name]['raw']

    print(f"[load_track_data()] Number of unique storms in track dataset: {track_dataset['storm_id'].nunique()}")

    return track_dataset

def filter_track_data(model_name: str,
                      experiment_name: str,
                      year_range: tuple[int, int],
                      basin_name: str,
                      track_dataset: pd.DataFrame) -> (pd.DataFrame, list[str]):
    
    # Load storm IDs from generated TC files
    TC_storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    # Get filenames for TCs with generated data that correspond to the model and experiment name
    pathnames = [os.path.join(TC_storage_dirname, filename) for filename in os.listdir(TC_storage_dirname) if
                 filename.endswith('.nc') and
                 f'model-{model_name}' in filename and
                 f'experiment-{experiment_name}' in filename]
    # Filter pathnames by basin, if basin_name is not global
    if basin_name != 'global':
        pathnames = [pathname for pathname in pathnames if basin_name in pathname]
    # Stop the script if no files matching the input criteria are found.
    assert len(pathnames) > 0, f'No files found matching the criteria provided.'
    
    # Scrape storm IDs from the filenames
    storm_IDs = [pathname.split('storm_ID-')[-1].split('.')[0] for pathname in pathnames]
    # Filter track dataset by storm IDs corresponding to TCs with generated data
    configuration_track_dataset = track_dataset.loc[track_dataset['storm_id'].isin(storm_IDs)]
    
    print(f"[filter_track_data()] Number of unique storms in filtered track dataset: {configuration_track_dataset['storm_id'].nunique()}")
    
    return configuration_track_dataset, pathnames

def get_snapshots_from_track_data(track_dataset: pd.DataFrame,
                                  pathnames: list[str],
                                  intensity_parameter: str,
                                  intensity_range: tuple,
                                  number_of_snapshots: int) -> pd.DataFrame:

    ''' This method filters through track data to identify snapshots to be used. '''
    
    assert intensity_parameter in track_dataset.columns, f'Intensity parameter not found in track data. Parameter must be one of: {track_dataset.columns}.'

    # 1. Filter track dataset by intensity
    track_dataset_bin = track_dataset.loc[((track_dataset[intensity_parameter] >= min(intensity_range)) &
                                           (track_dataset[intensity_parameter] <= max(intensity_range)))]
    # 2. Determine some number (`number_of_snapshots`) of random indices of storms in filtered track dataset
    number_of_snapshots = len(track_dataset_bin) if number_of_snapshots < 0 else number_of_snapshots # allow for all storms to be chosen if `number_of_storms` < 0
    print(f"[get_snapshots_from_track_data()] Number of snapshots: {number_of_snapshots}; number of track dataset entries: {len(track_dataset_bin)}")
    # 3. Use randomized indices to extract `number_of_snapshots` storms from the track dataset
    snapshots = track_dataset_bin.sample(number_of_snapshots)
    # 4. Get "uniqueness" value of the snapshot dataset. 
    #    0 means all snapshots are from the same TC, 1 means all snapshots are from different TCs.
    #    This measures the diversity of the TCs in the snapshot dataset. 
    snapshot_diversity = snapshots['storm_id'].nunique() / len(snapshots)
    if snapshot_diversity < 0.2:
        print(f'Warning: this snapshot group risks having too little diversity in unique TCs with {snapshots['storm_id'].nunique()} unique TCs in {len(snapshots)} snapshots.')

    # Define lambda function to apply to the `snapshots` DataFrame
    retrieve_filename_from_storm_ID = lambda f: [pathname for pathname in pathnames if f['storm_id'] in pathname][0]
    # Pull filenames for the corresponding storm IDs
    snapshots['pathname'] = snapshots.apply(retrieve_filename_from_storm_ID, axis=1)

    return snapshots

def get_snapshot_timestamp(snapshot_entry: dict,
                           snapshot: xr.Dataset):

    ''' Ensure calendar type of timestamp is aligned with calendar type of model data. '''

    # Get iterand timestamp
    snapshot_entry_timestamp = snapshot_entry['cftime']
    # Get timestamp format of model data time axis as a cftime object. 
    # Assumes all timestamps are of the same format.
    snapshot_timestamp = snapshot.time.values[0]
    # Get calendar format. Assumes calendar types are equivalent for all timestamps.
    timestamp_calendar_format = utilities.cftime_calendar_type(snapshot_timestamp)
    # Define the identifier used as a reference time for cftime timedelta calculation
    reference_identifier = f'seconds since {snapshot_timestamp}'

    # Perform calendar type conversion
    snapshot_entry_timestamp_seconds = cftime.date2num(snapshot_entry_timestamp, 
                                                       reference_identifier, 
                                                       calendar=timestamp_calendar_format)                         
    snapshot_entry_timestamp = cftime.num2date(snapshot_entry_timestamp_seconds, 
                                               reference_identifier,
                                               calendar=timestamp_calendar_format)

    assert (type(snapshot_timestamp) == type(snapshot_entry_timestamp)), f'Timestamps are of dissimilar types.'
    
    return snapshot_entry_timestamp

def storm_interpolation_grid_basis_vectors(storm: xr.Dataset,
                                           window_size: int,
                                           equal_number_of_points: bool=True,
                                           diagnostic: bool=False) -> tuple[np.array, np.array]:

    ''' 
    Method to generate uniform interpolation basis vectors for a TC-centered grid. 
    Boolean `equal_number_of_points` is used to ensure equal grid point numbers, and is optional.
    '''

    # Ensure necessary basis vector dimensions are available
    assert ('grid_xt' in storm.dims) and ('grid_yt' in storm.dims)
    # Get differences in grid spacing along each vector
    d_grid_yt = storm['grid_yt'].diff(dim='grid_yt')
    d_grid_xt = storm['grid_xt'].diff(dim='grid_xt')

    # Ensure that the differences are equivalent for all indices to ensure equal spacing
    grid_tolerance = 1e-6

    assert sum(d_grid_xt.diff(dim='grid_xt') < grid_tolerance) == len(d_grid_xt) - 1, 'Grid is irregular along the `grid_xt` axis.'
    assert sum(d_grid_yt.diff(dim='grid_yt') < grid_tolerance) == len(d_grid_yt) - 1, 'Grid is irregular along the `grid_yt` axis.'
    # Get grid spacing along each direction
    dx, dy = d_grid_xt.isel(grid_xt=0).item(), d_grid_yt.isel(grid_yt=0).item()
    # Get the number of grid points in each direction for vector construction
    # Use the number of y-grid points if an equal number of points is desired in each direction, since y-resolution is typically smaller
    number_y_grid_points = int(np.around((window_size * 2) / dy))
    number_x_grid_points = number_y_grid_points if equal_number_of_points else int(np.around((window_size * 2) / dx))
    # Construct the basis vectors for Dataset interpolation
    arr_x_interp = np.linspace(-window_size, window_size, number_x_grid_points)
    arr_y_interp = np.linspace(-window_size, window_size, number_y_grid_points)

    if diagnostic:
        print(f'[storm_interpolation_gridpoints]: Number of grid points: x = {number_x_grid_points}, y = {number_y_grid_points}')

    return arr_x_interp, arr_y_interp

def storm_centered_interpolation(snapshot: xr.Dataset,
                                 field_name: str,
                                 window_size: int=10) -> xr.DataArray:

    # Generate storm-centered coordinates.
    # This will remove dependence on global coordinates and allow for storm-centered compositing.
    arr_x = snapshot.grid_xt - snapshot['center_lon']
    arr_y = snapshot.grid_yt - snapshot['center_lat']

    # Generate storm ID list to serve as an axis for the storm identifier that will enable xArray-based compositing
    storm_ID = [snapshot.attrs['storm_id']]
    # Expand dimensions of 2D data for the storm ID axis
    snapshot_dataset = np.expand_dims(snapshot[field_name].data, axis=0)

    # Generate the xArray DataArray
    snapshot_dataset = xr.DataArray(data=snapshot_dataset,
                          dims=['storm_id', 'grid_yt_TC', 'grid_xt_TC'],
                          coords={'grid_xt_TC': (['grid_xt_TC'], arr_x.data),
                                  'grid_yt_TC': (['grid_yt_TC'], arr_y.data),
                                  'storm_id': (['storm_id'], storm_ID)})

    # To allow for the combination of different TCs, interpolate the storm-centered coordinates to common axis values based on the window size.
    arr_x_interp, arr_y_interp = storm_interpolation_grid_basis_vectors(snapshot, window_size)

    # Perform an interpolation from the storm-centered coordinates to the interpolated coordinates
    interpolated_snapshot_dataset = snapshot_dataset.interp(grid_xt_TC=arr_x_interp).interp(grid_yt_TC=arr_y_interp)

    return interpolated_snapshot_dataset

def get_local_time(timestamp: cftime.datetime, 
                   latitude: float, 
                   longitude: float,
                   diagnostic: bool=False):

    if diagnostic:
        print(f'[get_local_time()] Timestamp value input: {timestamp}.')

    # Account for International Date Line
    longitude = longitude if longitude <= 180 else longitude - 360
    # Get timezone designator
    timezone_str = timezonefinder.TimezoneFinder().timezone_at(lng=longitude, lat=latitude)

    if timezone_str:
        timezone = pytz.timezone(timezone_str)
        timestamp_datetime_GMT = datetime.datetime(year=timestamp.year, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        timestamp_datetime_local = timestamp_datetime_GMT + timezone.utcoffset(timestamp_datetime_GMT)
        timestamp_cftime_local = cftime.datetime(year=timestamp_datetime_local.year,
                                                 month=timestamp_datetime_local.month,
                                                 day=timestamp_datetime_local.day,
                                                 hour=timestamp_datetime_local.hour)
        return timestamp_cftime_local
    else:
        return "Timezone not found."

def hour_filter(composite_snapshots: dict,
                start_hour: int=6,
                end_hour: int=18) -> dict:

    ''' Method to filter TC samples by the hour at which they were sampled. '''

    assert start_hour != end_hour, 'Hours must be different.'
    assert isinstance(start_hour, int) & isinstance(end_hour, int), 'Hours must be integers.'
    
    filtered_composite_snapshots = {}
    for configuration in composite_snapshots.keys():
        if start_hour < end_hour:
            filtered_composite_snapshots[configuration] = composite_snapshots[configuration].where((composite_snapshots[configuration]['local_time'].dt.hour >= start_hour) & 
                                                                                                   (composite_snapshots[configuration]['local_time'].dt.hour <= end_hour)).dropna('storm_id', how='all')
        else:
            
            filtered_composite_snapshots[configuration] = composite_snapshots[configuration].where((composite_snapshots[configuration]['local_time'].dt.hour >= start_hour) | 
                                                                                                   (composite_snapshots[configuration]['local_time'].dt.hour <= end_hour)).dropna('storm_id', how='all')
            
    return filtered_composite_snapshots

def sampling_hour_histogram(composite_snapshots):

    ''' This method plots a histogram of the hours at which samples were taken, in each samples' local time. '''

    storm_sampling_hour = {configuration: [] for configuration in composite_snapshots.keys()}
    for configuration in composite_snapshots.keys():
        for storm_ID in range(len(composite_snapshots[configuration]['storm_id'].values)):
            storm_sampling_hour[configuration].append(composite_snapshots[configuration].isel(storm_id=storm_ID)['local_time'].dt.hour)
    
    hour_bins = range(24)
    fig, ax = plt.subplots(figsize=(4, 2))
    for configuration in composite_snapshots.keys():
        number_of_samples = len(storm_sampling_hour[configuration])
        ax.hist(storm_sampling_hour[configuration], bins=hour_bins, label=f'{configuration}: N = {number_of_samples}', histtype='step', lw=2)
    fig.tight_layout()
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1.1))
    ax.set_xlabel('Local hour of day');

def derived_quantities(model_name: str,
                       dataset: xr.Dataset):

    # Correct sign convention from positive down to positive up for upwards-directed quantities
    dataset['lwup_sfc'] = dataset['lwup_sfc'] * -1 if model_name in ['ERA5'] else dataset['lwup_sfc']
    dataset['swup_sfc'] = dataset['swup_sfc'] * -1 if model_name in ['ERA5'] else dataset['swup_sfc']

    # Correct for units from kg/m^2/s to mm/d
    if model_name in ['HIRAM', 'AM2.5', 'FLOR']:
        for field_name in dataset.data_vars:
            factor = 86400 if field_name in ['precip', 'evap', 'p-e'] else 1
            dataset[field_name] = dataset[field_name] * factor

    # Generate derived fields
    dataset = derived.net_lw(dataset)
    dataset = derived.net_sw(dataset)
    if model_name in ['HIRAM', 'AM2.5', 'FLOR']:
        dataset = derived.TC_lhflx(dataset)
    if model_name not in ['CERES', 'MERRA']:
        dataset = derived.TC_surface_wind_speed(dataset)
        # dataset = derived.atmospheric_heating(dataset)

    return dataset

def generate_TC_snapshot(model_name: str,
                         field_name: str,
                         window_size: int,
                         snapshot_entry: dict) -> pd.DataFrame:

    ''' Filter and select snapshots from track data. '''
    
    # Get iterand pathname
    snapshot_pathname = snapshot_entry['pathname']
    
    # Pull data for the TC at the iterand timestamp
    snapshot = xr.open_dataset(snapshot_pathname, use_cftime=True)

    # Acquire and process the snapshot timestamp
    snapshot_timestamp = get_snapshot_timestamp(snapshot_entry, snapshot)
    # Get local time corresponding to timestamp
    snapshot_local_time = get_local_time(snapshot_timestamp, 
                                         snapshot.sel(time=snapshot_timestamp, method='nearest')['center_lat'].item(), 
                                         snapshot.sel(time=snapshot_timestamp, method='nearest')['center_lon'].item())
    
    # Select the converted timestamp and remove null data
    snapshot = snapshot.sel(time=snapshot_timestamp, method='nearest').dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    
    # Get derived fields
    snapshot = derived_quantities(model_name, snapshot)
    # Perform grid redefinition and associated spatial interpolation
    interpolated_snapshot = storm_centered_interpolation(snapshot, 
                                                         field_name=field_name, 
                                                         window_size=window_size)
    interpolated_snapshot['local_time'] = snapshot_local_time
    
    return interpolated_snapshot

def generate_TC_composite(model_name: str,
                          experiment_name: str,
                          year_range: tuple[int, int],
                          field_name: str,
                          basin_name: str='global',
                          intensity_parameter: str='min_slp',
                          intensity_range: tuple[int|float, int|float]=(0, np.inf),
                          number_of_snapshots: int=10,
                          composite_window_size: int=10,
                          parallel: bool=False) -> tuple[xr.Dataset, dict]:

    ''' 1, Pull track data, filter by inputs, and generate QuickTracks-like entries for snapshots. '''

    # Load track data
    # Note: loading track is handled separately to handle IBTrACS
    track_dataset = load_track_data(model_name,
                                    experiment_name,
                                    year_range)
    
    # Filter GFDL Quick Tracks data based on the provided inputs
    track_dataset, pathnames = filter_track_data(model_name,
                                                 experiment_name,
                                                 year_range,
                                                 basin_name,
                                                 track_dataset)
    
    # Return track dataset for snapshots
    snapshots = get_snapshots_from_track_data(track_dataset,
                                              pathnames,
                                              intensity_parameter,
                                              intensity_range,
                                              number_of_snapshots)
    # Generate dictionary for iterating through each entry, corresponding to a snapshot
    snapshot_entries = snapshots.to_dict('records')

    ''' 2. Use snapshot data to pull TC data at individual timestamps and concatenate them into single data object. '''
    
    # Initialize list to contain all individual TCs samples
    snapshot_container = []
    # Define partial function to allow for using Pool.map since all inputs except `pathname` are the same for each storm
    # Generate a partial function for easier parallelization
    partial_generate_TC_snapshot = functools.partial(generate_TC_snapshot,
                                                     model_name,
                                                     field_name,
                                                     composite_window_size)
    
    # Gather TC samples for compositing
    if parallel:
        # Maximum number of processors for computation
        max_number_procs = 20
        # Specify number of processors to use
        number_procs = len(snapshot_entries) if len(snapshot_entries) < max_number_procs else max_number_procs
        
        with Pool(processes=number_procs) as pool:
            snapshot_container = pool.map(partial_generate_TC_snapshot, snapshot_entries)
            print(f'Number of snapshots in container: {len(snapshot_container)}')
            pool.close()
    # Perform it serially
    else:
        for snapshot_entry in snapshot_entries:
            # Perform grid redefinition and associated spatial interpolation
            TC_snapshot = partial_generate_TC_snapshot(snapshot_entry)
            # Append to a container list for future concatenation
            snapshot_container.append(TC_snapshot)
    
    # Merge the data
    snapshot_composite = xr.concat(snapshot_container, dim='storm_id')
    # Add model and experiment designators
    snapshot_composite.attrs['model_name'] = model_name
    snapshot_composite.attrs['experiment_name'] = experiment_name

    return snapshot_composite, snapshot_entries

def get_time_window(dataset: xr.Dataset,
                    timestamp,
                    window_day_size: int):

    ''' Filter a dataset by a given timestamp +/- a specific number of days. '''

    # Obtain day of year for the timestamp
    timestamp_day_of_year = timestamp.dayofyr if 'cftime' in str(type(timestamp)) else timestamp.dt.dayofyear
    # Dataset time array days of year (handled differently by time object type)
    dataset_day_of_year = dataset.time.dt.dayofyear if 'cftime' in str(type(timestamp)) else dataset.time.dt.dayofyear
    # Get start and end days of year
    start_day_of_year, end_day_of_year = timestamp_day_of_year - window_day_size, timestamp_day_of_year + window_day_size
    # Mask by window from start_day_of_year to end_day_of_year
    window = (dataset_day_of_year >= start_day_of_year) & (dataset_day_of_year <= end_day_of_year)
    # Mask data by the window
    dataset_window = dataset.sel(time=window)

    return dataset_window

def sample_GCM_constructor(snapshot: xr.Dataset,
                           GCM_snapshot: xr.Dataset,
                           field_name: str,
                           window_size: int=10):

    ''' Modify GCM data into a TC-centered dataset. ''' 

    # Add center coordinates to the GCM dataset
    GCM_snapshot['center_lon'] = snapshot['center_lon']
    GCM_snapshot['center_lat'] = snapshot['center_lat']
    GCM_snapshot.attrs = snapshot.attrs
    # Perform storm-centered interpolation for GCM data
    interpolated_GCM_snapshot = storm_centered_interpolation(snapshot=GCM_snapshot,
                                                             field_name=field_name,
                                                             window_size=window_size)    

    return interpolated_GCM_snapshot

def get_sample_GCM_data(model_name: str,
                        experiment_name: str,
                        field_name: str,
                        year_range: tuple[int, int],
                        snapshot_timestamp: pd.Timestamp,
                        longitude: int|float,
                        latitude: int|float,
                        window_size: int,
                        sampling_day_window: int=5):

    ''' Method to pull GCM data corresponding to a given TC snapshot. '''

    # Construct field dictionary for postprocessed data loading
    # See `utilities.postprocessed_data_load` for details.
    # Note: this currently only supports single-surface atmospheric data
    field_dictionary = {field_name: {'domain': 'atmos', 'level': None}}
    # Extract month from the iterand timestamp to perform initial climatology filtering
    snapshot_year, snapshot_month, snapshot_day = [snapshot_timestamp.year,
                                                   snapshot_timestamp.month,
                                                   snapshot_timestamp.day,]
    # Load the data
    sample_GCM_data = utilities.postprocessed_data_load(model_name,
                                                        experiment_name,
                                                        field_dictionary,
                                                        year_range,
                                                        data_type='mean_daily',
                                                        month_range=(snapshot_month, snapshot_month),
                                                        load_full_time=True)[model_name][experiment_name]
    
    # Define spatial extent for sample clipping
    grid_xt_extent = slice(longitude - window_size, longitude + window_size)
    grid_yt_extent = slice(latitude - window_size, latitude + window_size)
    # Trim the data spatially
    sample_GCM_data_filtered_space = sample_GCM_data.sortby('grid_yt').sel(grid_xt=grid_xt_extent).sel(grid_yt=grid_yt_extent)
    # Subsample over the time window specified: (iterand timestamp - sampling_day_window) to (iterand_timestamp + sampling_day_window)
    sample_GCM_data_filtered_time = get_time_window(sample_GCM_data_filtered_space, snapshot_timestamp, sampling_day_window)
    
    # Average in time
    sample_GCM_data_filtered = sample_GCM_data_filtered_time.mean(dim='time')

    return sample_GCM_data_filtered

def generate_composite_climatology_sample(model_name: str,
                                          experiment_name: str,
                                          field_name: str,
                                          year_range: tuple[int, int],
                                          window_size: int, 
                                          troubleshooting: bool,
                                          snapshot_entry: str) -> xr.Dataset:    
    ''' 
    Generates a climatological snapshot for a given TC to be used as a sample for composite analysis. 
    Snapshot is based on provided conditions. 
    '''

    # Get iterand pathname
    snapshot_pathname = snapshot_entry['pathname']
    
    # Pull data for the TC at the iterand timestamp
    snapshot = xr.open_dataset(snapshot_pathname, use_cftime=True)
    
    # Acquire and process the snapshot timestamp
    snapshot_timestamp = get_snapshot_timestamp(snapshot_entry, snapshot)    
    # Select the converted timestamp and remove null data
    snapshot = snapshot.sel(time=snapshot_timestamp, method='nearest').dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    
    # Get sample TC center coordinates
    snapshot_center_longitude = snapshot['center_lon'].item()
    snapshot_center_latitude = snapshot['center_lat'].item()

    # Load GCM data corresponding to the given TC sample
    sample_GCM_data = get_sample_GCM_data(model_name, 
                                          experiment_name,
                                          field_name,
                                          year_range,
                                          snapshot_timestamp,
                                          snapshot_center_longitude,
                                          snapshot_center_latitude,
                                          window_size)
    
    # Construct TC-centered GCM xarray object
    sample_GCM_data = sample_GCM_constructor(snapshot, sample_GCM_data, field_name)

    return sample_GCM_data

def generate_composite_climatology(model_name: str,
                                   experiment_name: str,
                                   year_range: tuple[int, int],
                                   field_name: str,
                                   snapshot_entries: dict,
                                   window_size: int=12,
                                   troubleshooting: bool=False,
                                   parallel: bool=False):

    ''' Gather climatological samples. '''
    # Initialize list to contain all individual TCs samples
    GCM_sample_container = []
    # Generate partial function
    partial_generate_composite_climatology_sample = functools.partial(generate_composite_climatology_sample, 
                                                                      model_name,
                                                                      experiment_name,
                                                                      field_name,
                                                                      year_range,
                                                                      window_size,
                                                                      troubleshooting)
    # Perform loading in parallel
    if parallel:
        ''' Offload TC-specific data generation onto parallel processes. '''
        # Maximum number of processors for computation
        max_number_procs = 20
        # Specify number of processors to use
        number_procs = len(snapshot_entries) if len(snapshot_entries) < max_number_procs else max_number_procs
        
        with Pool(processes=number_procs) as pool:
            GCM_sample_container = pool.map(partial_generate_composite_climatology_sample, snapshot_entries)
            pool.close()
    # Perform it serially
    else:
        for snapshot_entry in snapshot_entries:
            # Construct TC-centered GCM xarray object
            sample_GCM_data = partial_generate_composite_climatology_sample(snapshot_entry)
            # Append to container list
            GCM_sample_container.append(sample_GCM_data)
        
    # Concatenate all samples 
    composite_GCM_samples = xr.concat(GCM_sample_container, dim='storm_id')
    
    return composite_GCM_samples

def add_colorbar(fig, 
                 ax,
                 field_name: str,
                 norm, 
                 cmap,
                 orientation: str='vertical'):

    # Define colorbar axis
    colorbar_axis_thickness = 0.05 # in units of axes fraction
    colorbar_axis_padding = 0.025 # in units of axes fraction
    # Generate the colorbar axis size and position
    colorbar_axis_extent = [1 + colorbar_axis_padding, 0, colorbar_axis_thickness, 1] if orientation == 'vertical' else [0, 0 - colorbar_axis_padding, 1, colorbar_axis_thickness]
    # Initialize colorbar axis
    colorbar_axis = ax.inset_axes(colorbar_axis_extent)
    # Generate the colorbar
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=colorbar_axis)

    return colorbar_axis

def domain_masking(composite_samples: xr.Dataset,
                   outer_radius: int,
                   inner_radius: int=0):

    # Get composite means, but handle averaging depending on present dimensions
    composite_samples_mean = composite_samples.mean('storm_id') if 'storm_id' in composite_samples.dims else composite_samples
    # Get mask for composite TC domain center
    domain_mask = utilities.circular_mask(X=composite_samples_mean.grid_xt_TC,
                                          Y=composite_samples_mean.grid_yt_TC,
                                          x=0,
                                          y=0,
                                          inner_radius=inner_radius,
                                          outer_radius=outer_radius)

    return composite_samples.where(domain_mask)

def plot_metadata(ax,
                  composite_samples: xr.Dataset | xr.DataArray,
                  composite_mean: xr.DataArray,
                  field_name: str,
                  configuration_name: str,
                  window_radius: int=10,
                  center_window_radius: int=5):

    # Initialize list of annotations for common formatting functions to be applied at end of function
    annotations = []
    # Get field-specific metadata
    long_name, units = visualization.field_properties(field_name)

    ''' Insert annotation for TC center mean, with radial extent from center set by `center_window_radius`. '''
    # Get mask for composite TC domain center
    center_mask = domain_masking(composite_samples, inner_radius=0, outer_radius=center_window_radius)
    # Apply mask to data, get mean and standard deviation with respect to TC snapshots
    composite_mean_center = center_mask.mean(['grid_xt_TC', 'grid_yt_TC'])
    composite_std_center = center_mask.mean(['grid_xt_TC', 'grid_yt_TC'])
    # Get mean and standard deviation values over central window
    center_field_mean_value = composite_mean_center.mean(dim='storm_id') if 'storm_id' in composite_mean_center.dims else composite_mean_center
    center_field_std_value = composite_std_center.std(dim='storm_id') if 'storm_id' in composite_std_center.dims else composite_std_center
    # Generte annotation string
    center_value_annotation_string = f'{center_window_radius} deg. avg.: {center_field_mean_value:.1f} {units}'
    # Generate annotation
    center_value_annotation_string = ax.annotate(center_value_annotation_string, 
                                                 xy=(0.03, 0.97), 
                                                 xycoords='axes fraction', 
                                                 fontsize=10, 
                                                 ha='left', 
                                                 va='top')
    annotations.append(center_value_annotation_string)
    
    ''' Insert annotation for TC domain, with radial extent from center set by `window_radius`. '''
    # Get mask for composite TC domain center
    outer_mask = domain_masking(composite_samples, 
                                inner_radius=center_window_radius, 
                                outer_radius=window_radius)
    # Apply mask to data
    composite_mean_domain = outer_mask.mean('storm_id') if 'storm_id' in composite_std_center.dims else outer_mask
    # Get mean value over central window
    domain_field_mean_value = composite_mean_domain.mean()
    # Generte annotation string
    domain_value_annotation_string = f'{window_radius} deg. avg.: {domain_field_mean_value:.1f} {units}'
    # Generate annotation
    domain_value_annotation_string = ax.annotate(domain_value_annotation_string, 
                                                 xy=(0.03, 0.90), 
                                                 xycoords='axes fraction', 
                                                 fontsize=10, 
                                                 ha='left', 
                                                 va='top')
    annotations.append(domain_value_annotation_string)
    
    ''' Insert annotation for sample count. '''
    # Get sample count
    number_of_samples = f'N = {len(composite_samples['storm_id'].values)}' if 'storm_id' in composite_samples.dims else ''
    # Generate annotation
    sample_number_annotation = ax.annotate(number_of_samples, 
                                           xy=(0.03, 0.03), 
                                           xycoords='axes fraction', 
                                           fontsize=10, 
                                           ha='left', 
                                           va='bottom')
    annotations.append(sample_number_annotation)
    
    # Provide path effect for all annotations
    for annotation in annotations:
        annotation.set_path_effects([pe.Stroke(linewidth=1.5, foreground='white'), pe.Normal()])
    
    title_string = f'Composite mean {long_name} [{units}]\n{configuration_name}'
    ax.set_title(title_string, loc='left', ha='left', fontsize=10);

def plot_guidelines(ax,
                    domain_radius: int,
                    center_radius: int,):

    ''' Helper function to plot guidelines for TC composite plots .'''

    ax.axhline(0, c='k', lw=0.5, ls='--', alpha=0.25)
    ax.axvline(0, c='k', lw=0.5, ls='--', alpha=0.25)
    
    # Plot TC center extent for analysis purposes
    TC_center_guideline = matplotlib.patches.Circle((0, 0), 
                                                    radius=domain_radius, 
                                                    fc='none', 
                                                    lw=0.5, 
                                                    ls='--', 
                                                    ec='k',
                                                    alpha=0.55)
    ax.add_patch(TC_center_guideline)
    # Plot TC domain extent for analysis purposes
    TC_domain_guideline = matplotlib.patches.Circle((0, 0), 
                                                    radius=center_radius, 
                                                    fc='none', 
                                                    lw=0.5, 
                                                    ls='--', 
                                                    ec='k',
                                                    alpha=0.55)
    ax.add_patch(TC_domain_guideline)    

def plot_composite_single(composite_samples: xr.Dataset | xr.DataArray,
                          field_name: str,
                          configuration_name: str,
                          plotting_method: str='contourf',
                          number_of_normalization_bins: int=16,
                          fig=None,
                          ax=None,
                          norm=None,
                          cmap=None):
    
    ''' Plot composite values for a given model and experiment. '''

    # Define radial extents for analytical outputs and guidelines
    composite_domain_radius = 10
    composite_center_radius = 2.5

    # Get mean and standard deviation over all TC samples
    composite_mean = composite_samples.mean(dim='storm_id') if 'storm_id' in composite_samples.dims else composite_samples

    # Get normalization and colormap for the composite plot if not provided
    if not norm and not cmap:
        norm, cmap = visualization.norm_cmap(composite_mean, 
                                             field=field_name, 
                                             num_bounds=number_of_normalization_bins)

    ''' Plotting. '''
    
    # Initialize figure if an initial figure not provided
    if not fig and not ax:
        fig, ax = plt.subplots(figsize=(4, 4))
    # Define plotting method
    if plotting_method == 'contourf':
        im = ax.contourf(composite_mean.grid_xt_TC, composite_mean.grid_yt_TC, composite_mean, 
                         norm=norm, cmap=cmap, levels=len(norm.boundaries))
    else:
        im = ax.pcolormesh(composite_mean.grid_xt_TC, composite_mean.grid_yt_TC, composite_mean, 
                           norm=norm, cmap=cmap)
        
    # Plot guidelines
    plot_guidelines(ax, composite_domain_radius, composite_center_radius)
    
    # Ensure subplot is square
    ax.set_aspect('equal')

    # Append the colorbar
    colorbar_axis = add_colorbar(fig, ax, field_name, norm, cmap)

    # Append plot metadata
    plot_metadata(ax, composite_samples, composite_mean, field_name, configuration_name, composite_domain_radius, composite_center_radius)

    return fig, ax

def plot_composite(composite_TC_samples: dict,
                   configurations: str | list[str] | None=None,
                   plotting_method: str='contourf',
                   number_of_normalization_bins: int=16,
                   get_difference: bool=False,
                   difference_order: str='ltr',
                   extrema: tuple[int|float, int|float] | None=None,
                   plot_orientation: str='horizontal',
                   dpi: int=144):

    ''' Composite plotting method. '''

    # Assign configurations and ensure proper variable typing
    if isinstance(configurations, str):
        configurations = [configurations]
    elif not configurations:
        configurations = list(composite_TC_samples.keys())

    # Get number of configurations for determining figure size
    number_of_configurations = len(configurations)

    # Ensure that a difference plot is generate only if there are two experiments
    if get_difference:
        assert number_of_configurations == 2, 'Cannot get difference with fewer or more than 2 experiments.'
        number_of_configurations = number_of_configurations + 1
    
    # Combine all composite samples to generate uniform extrema
    combined_composite_samples = xr.concat([sample.mean(dim='storm_id') for sample in composite_TC_samples.values()], 
                                           dim='storm_id')
    # Get common normalization and colormap
    norm, cmap = visualization.norm_cmap(combined_composite_samples, 
                                         field=field_name, 
                                         num_bounds=number_of_normalization_bins,
                                         extrema=extrema)
    
    ''' Begin plotting. '''
    
    # Define subplot size (width = height)
    axis_size = 3.5
    # Initialize figure and grid
    if plot_orientation == 'horizontal':
        fig = plt.figure(figsize=(axis_size * number_of_configurations, axis_size), dpi=dpi)
        gs = matplotlib.gridspec.GridSpec(nrows=1, ncols=number_of_configurations)
    else:
        fig = plt.figure(figsize=(axis_size, axis_size * number_of_configurations), dpi=dpi)
        gs = matplotlib.gridspec.GridSpec(nrows=number_of_configurations, ncols=1)
    
    # Plot composite mean for each model-experiment configuration
    for axis_index, configuration_name in enumerate(configurations):
        # Initialize axis
        configuration_axis = fig.add_subplot(gs[0, axis_index]) if plot_orientation == 'horizontal' else fig.add_subplot(gs[axis_index, 0])
        # Plot the composite for the corresponding configuration
        plot_composite_single(composite_TC_samples[configuration_name], 
                              field_name, 
                              plotting_method=plotting_method,
                              configuration_name=configuration_name,
                              fig=fig,
                              ax=configuration_axis,
                              norm=norm,
                              cmap=cmap)

    # If the difference boolean is enabled, get the difference
    if get_difference:
        # Initialize the corresponding subplot axis
        configuration_axis = fig.add_subplot(gs[0, -1]) if plot_orientation == 'horizontal' else fig.add_subplot(gs[-1, 0])
        # Determine control and experiment configurations
        if difference_order == 'ltr':
            configuration_CTL, configuration_EXP = configurations[0], configurations[1]
        else:
            configuration_CTL, configuration_EXP = configurations[1], configurations[0] 
        experiment_difference_configuration = f'{configuration_CTL} - {configuration_EXP}'
        # Calculate difference
        composite_difference = composite_TC_samples[configuration_CTL].mean('storm_id') - composite_TC_samples[configuration_EXP].mean('storm_id')
        # Get difference normalization and colormap
        difference_norm, difference_cmap = visualization.norm_cmap(composite_difference, 
                                                                   field=field_name, 
                                                                   num_bounds=number_of_normalization_bins,
                                                                   white_adjust=True)
        # Plot the composite for the corresponding configuration
        plot_composite_single(composite_difference, 
                              field_name, 
                              plotting_method=plotting_method,
                              configuration_name=experiment_difference_configuration,
                              fig=fig,
                              ax=configuration_axis,
                              norm=difference_norm,
                              cmap=difference_cmap)
        
    fig.tight_layout()

def save_composites(composite_TC_samples: dict,
                    configuration_name: str,
                    field_name: str,
                    year_range: list[int, int] | tuple[int, int],
                    intensity_range: list[int, int] | tuple[int, int],
                    basin_name: str,
                    dirname: str):
    
    filename = f'configuration.{configuration_name}.field_name.{field_name}.year_range.{min(year_range)}_{max(year_range)}.intensity_range.{min(intensity_range)}_{max(intensity_range)}.basin.{basin_name}.nc'
    pathname = os.path.join(dirname, filename)

    print(f'Saving composite data to {pathname}...')

    if os.path.exists(pathname):
        override = input('File already exists - do you wish to override? [y/n]\n')
        if 'y' in override:
            os.remove(pathname)
            composite_TC_samples[configuration_name].to_dataset(name=field_name).to_netcdf(pathname)
        else:
            import sys
            print(f'File at {pathname} not overridden. Exiting...')
            sys.exit()
    else:
        composite_TC_samples[configuration_name].to_dataset(name=field_name).to_netcdf(pathname)

def composite_convergence_check(composite_samples: xr.Dataset,
                                configuration_name: str,
                                field_name: str,
                                inner_radius: int|float=0,
                                outer_radius: int|float=5,
                                visualize: bool=True):
    
    ''' 
    Method to check if the composite reaches an asymptote using a running mean over some spatial window. 
    If yes, that means an adequate number of samples have been used for compositing.
    '''

    # Initialize dictionary to hold sample number as a key, with running mean as a value
    running_mean = {}
    # Iterate along the storm ID axis
    for index in range(len(composite_samples.storm_id)):
        # Get the running mean over the selected region
        window_mask = domain_masking(composite_samples.isel(storm_id=range(index)),
                                     outer_radius=outer_radius,
                                     inner_radius=inner_radius)
        window_mean = window_mask.mean().item()
        running_mean[index] = window_mean

    # Get numeric array from samples
    running_mean_samples = np.fromiter(running_mean.values(), dtype=float)
    # Check for steady state
    steady_state_tolerance = 0.1
    correlations, lag = check_steady_state(running_mean_samples,
                                           steady_state_tolerance=steady_state_tolerance)

    ''' Visualization. '''
    if visualization:
        long_name, units = visualization.field_properties(field_name)
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.plot(running_mean_samples)
        ax.set_xlabel('Sample number')    
        ax.set_ylabel(units)
    
        autocorrelation_ax = ax.twinx()
        autocorrelation_ax.axhline(1 - steady_state_tolerance, c='k', lw=0.5, ls='--', alpha=0.5)
        autocorrelation_ax.axhline(1 + steady_state_tolerance, c='k', lw=0.5, ls='--', alpha=0.5)
        autocorrelation_ax.plot(correlations.keys(), correlations.values(), c='tab:green')
        autocorrelation_ax.set_ylim([0, 2])
        autocorrelation_ax.set_ylabel('Autocorrelation', color='tab:green', rotation=270, labelpad=15)
        
        ax.set_title(f'Composite convergence of {field_name}, {configuration_name}\nAutocorrelation lag size: {lag} samples', loc='left', ha='left', fontsize=10)
    
        fig.tight_layout()

def check_steady_state(sample: np.array,
                       lag_window_percentile: int=10,
                       steady_state_tolerance: float=0.001):

    '''
    lag_window_percentile: percent of array length that is used as a lag window.
    steady_state_tolerance:  = 0.01
    '''

    # Define shorthand for number of samples
    sample_length = sample.shape[0]
    # Number of indices used for autocorrelation
    lag = int(np.round(sample_length * lag_window_percentile / 100)) 
    
    # Initialize dictionary to hold correlation coefficients as a function of sample index
    correlations = {}
    # Iterate over indices to obtain lag autocorrelations
    for index in range(lag + 1, sample_length):
        # Get correlation values
        corr = np.correlate(sample[index - lag:index], sample[index - lag:index], mode='full')
        result = corr[corr.size // 2:]
        # Obtain last value in correlation window, multiply by lag window sample size to normalize to 1
        correlations[index] = (result / result.max())[-1] * lag

    # Convert correlation values to a numeric array
    autocorrelation_values = np.fromiter(correlations.values(), dtype=float)
    print(f'Number of autocorrelation values: {len(autocorrelation_values)}')
    
    # Iterate through autocorrelation to determine if and when steady-state is reached
    # Steady state is defined as an autocorrelation 1 +/- some defined tolerance
    steady_state_window = (1 - steady_state_tolerance, 1 + steady_state_tolerance)
    
    # Use while loop, which is broken when the running mean autocorrelation value falls within the threshold window
    steady_state = False
    autocorrelation_index = lag
    while not steady_state and autocorrelation_index < len(autocorrelation_values):
        # Define start and end indices
        start_index, end_index = autocorrelation_index - lag, autocorrelation_index
        # Obtain mean of autocorrelation values over the lag window
        autocorrelation_window = autocorrelation_values[start_index:end_index]

        autocorrelation_window_check = np.where((min(steady_state_window) < autocorrelation_window) &
                                                (autocorrelation_window < max(steady_state_window)), True, False)
        
        # if autocorrelation_index % 50 == 0:
        #     print(start_index, end_index)
            
        # If a steady-state is reached, end the while loop
        if sum(autocorrelation_window_check) == len(autocorrelation_window):
            steady_state = True
            # print(start_index, end_index, autocorrelation_window, sum(autocorrelation_window_check) == len(autocorrelation_window))
            print(f'Steady-state achieved by sample {autocorrelation_index}.')
        else:
            autocorrelation_index += 1
    
    if not steady_state:
        print(f'Steady-state not achieved.')

    return correlations, lag

def main(compositing_mode: str,
         model_name: str,
         experiment_name: str,
         year_range: tuple[int|float, int|float],
         field_name: str,
         basin_name: str='global',
         intensity_parameter: str='min_slp',
         intensity_range: tuple[int|float, int|float]=(0, np.inf),
         number_of_snapshots: int=-1,
         parallel: bool=False) -> xr.Dataset | tuple[xr.Dataset, xr.Dataset]:

    # Constrain analysis modes to TC (tropical cyclone composites only) or anomaly (tropical cyclone and climatological composites)
    compositing_modes = ['TC', 'anomaly']
    assert compositing_mode in compositing_modes, f'Compositing mode must be one of {compositing_modes}. Please retry by specifying argument `compositing_mode` as one of these.' 

    start_time = time.time()
    
    # Obtain composites for TCs only
    if compositing_mode == 'TC':
        composite_TC_snapshots, _ = generate_TC_composite(model_name,
                                                       experiment_name,
                                                       year_range,
                                                       field_name,
                                                       intensity_parameter=intensity_parameter,
                                                       intensity_range=intensity_range,
                                                       number_of_snapshots=number_of_snapshots,
                                                       parallel=parallel)
        
        print(f'Time elapsed for TCs: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(composite_TC_snapshots.storm_id.values)):.2f} s.')
        
        return composite_TC_snapshots
        
    # Obtain composites for TCs and GCM data only
    else:
        # Obtain TC composites
        composite_TC_snapshots, TC_snapshot_entries = generate_TC_composite(model_name,
                                                                            experiment_name,
                                                                            year_range,
                                                                            field_name,
                                                                            intensity_parameter=intensity_parameter,
                                                                            intensity_range=intensity_range,
                                                                            number_of_snapshots=number_of_snapshots,
                                                                            parallel=parallel)
        
        print(f'Time elapsed for TCs: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(TC_snapshot_entries)):.2f} s.')
        start_time = time.time()
        
        # Obtain GCM composites for times and locations corresponding to provided TC storm IDs
        composite_GCM_snapshots = generate_composite_climatology(model_name,
                                                                 experiment_name,
                                                                 year_range, 
                                                                 field_name,
                                                                 TC_snapshot_entries,
                                                                 parallel=parallel)

        print(f'Time elapsed for GCM data: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(TC_snapshot_entries)):.2f} s.')
        

        return (composite_TC_snapshots, composite_GCM_snapshots)

