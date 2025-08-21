# Functional and system imports
import os
import importlib
import functools
import time
from multiprocess import Pool
import timezonefinder
import datetime
import pytz
import random
# Analytical imports
import cftime
import nc_time_axis
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import xarray as xr
# Local imports
import derived
import tc_analysis
import utilities
import visualization


def derived_quantities(model_name: str,
                       dataset: xr.Dataset):

    # Correct sign convention from positive down to positive up for upwards-directed quantities
    if model_name in ['CERES']:
        if 'lwup_sfc' in dataset.data_vars:
            dataset['lwup_sfc'] = dataset['lwup_sfc'] * - \
                1 if model_name in ['ERA5'] else dataset['lwup_sfc']
            dataset = derived.net_lw(dataset)
        if 'swup_sfc' in dataset.data_vars:
            dataset['swup_sfc'] = dataset['swup_sfc'] * - \
                1 if model_name in ['ERA5'] else dataset['swup_sfc']
            dataset = derived.net_sw(dataset)

    # Correct for units from kg/m^2/s to mm/d
    if model_name in ['HIRAM', 'AM2.5', 'FLOR', 'AM4']:
        for field_name in dataset.data_vars:
            factor = 86400 if field_name in ['precip', 'evap', 'p-e'] else 1
            dataset[field_name] = dataset[field_name] * factor

    # Generate derived fields
    if model_name in ['HIRAM']:
        dataset = derived.TC_lhflx(dataset)
    if model_name not in ['CERES', 'MERRA']:
        dataset = derived.TC_surface_wind_speed(dataset)
        # dataset = derived.atmospheric_heating(dataset)

    return dataset


def field_name_data_types(field_name: str) -> str:
    ''' Retrieve data type associated with the input field. Assists with selection of file corresponding to TC. '''

    field_names = {'netrad_toa': 'atmos_4xdaily',
                   'olr': 'atmos_4xdaily',
                   'wind': 'atmos_4xdaily',
                   'slp': 'atmos_4xdaily',
                   'WVP': 'atmos_4xdaily',
                   'rh200': 'atmos_4xdaily',
                   'rh500': 'atmos_4xdaily',
                   'rh700': 'atmos_4xdaily',
                   'rh850': 'atmos_4xdaily',
                   'vort850': 'atmos_4xdaily',
                   't_surf': 'atmos_4xdaily',
                   'precip': 'atmos_4xdaily',
                   'evap': 'atmos_4xdaily',
                   'shflx': 'atmos_4xdaily',
                   'swup_toa': 'atmos_4xdaily',
                   'ucomp': 'atmos_daily',
                   'vcomp': 'atmos_daily',
                   'temp': 'atmos_daily',
                   'sphum': 'atmos_daily',
                   'cld_amt': 'atmos_daily', }

    assert field_name in field_names.keys(
    ), f'Field name {field_name} not registered. Please try again.'

    return field_names[field_name]


def load_track_data(model_name: str,
                    experiment_name: str,
                    year_range: tuple[int, int]):

    # Handle ERA5 case for IBTrACS data
    model_name = 'IBTrACS' if model_name in [
        'CERES', 'ERA5', 'MERRA'] else model_name
    experiment_name = '' if model_name == 'IBTrACS' else experiment_name

    # Load track data for the model-experiment configuration
    track_dataset = tc_analysis.load_TC_tracks(model_name,
                                               experiment_name,
                                               year_range,
                                               diagnostic=False)[model_name][experiment_name]['raw']

    print(
        f"[load_track_data()] Number of unique storms in track dataset: {track_dataset['storm_id'].nunique()}")

    return track_dataset


def filter_track_data(model_name: str,
                      experiment_name: str,
                      year_range: tuple[int, int],
                      field_name: str,
                      basin_name: str,
                      pressure_level: int | None,
                      track_dataset: pd.DataFrame,
                      TC_source_dirname: str = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs_3D') -> (pd.DataFrame, list[str]):

    # Get data type corresponding to field name
    data_type = field_name_data_types(field_name)
    # Get filenames for TCs with generated data that correspond to the model and experiment name
    pathnames = [os.path.join(TC_source_dirname, filename) for filename in os.listdir(TC_source_dirname) if
                 filename.endswith('.nc') and
                 f'model-{model_name}' in filename and
                 f'experiment-{experiment_name}' in filename]
    # Filter pathnames by basin, if basin_name is not global
    if basin_name != 'global':
        pathnames = [
            pathname for pathname in pathnames if basin_name in pathname]
    # Stop the script if no files matching the input criteria are found.
    assert len(pathnames) > 0, f'No files found matching the criteria provided.'

    # Scrape storm IDs from the filenames
    storm_IDs = [pathname.split('storm_ID-')[-1].split('.')[0]
                 for pathname in pathnames]
    # Filter track dataset by storm IDs corresponding to TCs with generated data
    configuration_track_dataset = track_dataset.loc[track_dataset['storm_id'].isin(
        storm_IDs)]

    print(
        f"[filter_track_data()] Number of unique storms in filtered track dataset: {configuration_track_dataset['storm_id'].nunique()}")

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
    # allow for all storms to be chosen if `number_of_storms` < 0
    number_of_snapshots = len(
        track_dataset_bin) if number_of_snapshots < 0 else number_of_snapshots
    print(
        f"[get_snapshots_from_track_data()] Number of snapshots: {number_of_snapshots}; number of track dataset entries: {len(track_dataset_bin)}")
    # 3. Use randomized indices to extract `number_of_snapshots` storms from the track dataset
    snapshots = track_dataset_bin.sample(number_of_snapshots)
    # 4. Get "uniqueness" value of the snapshot dataset.
    #    0 means all snapshots are from the same TC, 1 means all snapshots are from different TCs.
    #    This measures the diversity of the TCs in the snapshot dataset.
    snapshot_diversity = snapshots['storm_id'].nunique() / len(snapshots)
    if snapshot_diversity < 0.2:
        print(
            f'Warning: this snapshot group risks having too little diversity in unique TCs with {snapshots['storm_id'].nunique()} unique TCs in {len(snapshots)} snapshots.')

    # Define lambda function to apply to the `snapshots` DataFrame
    def retrieve_filename_from_storm_ID(
        f): return [pathname for pathname in pathnames if f['storm_id'] in pathname][0]
    # Pull filenames for the corresponding storm IDs
    snapshots['pathname'] = snapshots.apply(
        retrieve_filename_from_storm_ID, axis=1)

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
    timestamp_calendar_format = utilities.cftime_calendar_type(
        snapshot_timestamp)
    # Define the identifier used as a reference time for cftime timedelta calculation
    reference_identifier = f'seconds since {snapshot_timestamp}'

    # Perform calendar type conversion
    snapshot_entry_timestamp_seconds = cftime.date2num(snapshot_entry_timestamp,
                                                       reference_identifier,
                                                       calendar=timestamp_calendar_format)
    snapshot_entry_timestamp = cftime.num2date(snapshot_entry_timestamp_seconds,
                                               reference_identifier,
                                               calendar=timestamp_calendar_format)

    assert (type(snapshot_timestamp) == type(snapshot_entry_timestamp)
            ), f'Timestamps are of dissimilar types.'

    return snapshot_entry_timestamp


def storm_interpolation_grid_basis_vectors(storm: xr.Dataset,
                                           window_size: int,
                                           equal_number_of_points: bool = True,
                                           diagnostic: bool = False) -> tuple[np.array, np.array]:
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

    assert sum(d_grid_xt.diff(dim='grid_xt') < grid_tolerance) == len(
        d_grid_xt) - 1, 'Grid is irregular along the `grid_xt` axis.'
    assert sum(d_grid_yt.diff(dim='grid_yt') < grid_tolerance) == len(
        d_grid_yt) - 1, 'Grid is irregular along the `grid_yt` axis.'
    # Get grid spacing along each direction
    dx, dy = d_grid_xt.isel(grid_xt=0).item(), d_grid_yt.isel(grid_yt=0).item()
    # Get the number of grid points in each direction for vector construction
    # Use the number of y-grid points if an equal number of points is desired in each direction, since y-resolution is typically smaller
    number_y_grid_points = int(np.around((window_size * 2) / dy))
    number_x_grid_points = number_y_grid_points if equal_number_of_points else int(
        np.around((window_size * 2) / dx))
    # Construct the basis vectors for Dataset interpolation
    arr_x_interp = np.linspace(-window_size, window_size, number_x_grid_points)
    arr_y_interp = np.linspace(-window_size, window_size, number_y_grid_points)

    if diagnostic:
        print(
            f'[storm_interpolation_gridpoints]: Number of grid points: x = {number_x_grid_points}, y = {number_y_grid_points}')

    return arr_x_interp, arr_y_interp


def storm_centered_interpolation(snapshot: xr.Dataset,
                                 field_name: str,
                                 window_size: int = 10,
                                 model_name: str = None,
                                 diagnostic: bool = False) -> xr.DataArray:

    if diagnostic:
        print(
            f'[storm_centered_interpolation] available snapshot data variables: {snapshot.data_vars}')

    # Define dimension names based on the model
    dimension_name_x = 'longitude' if model_name == 'CERES' else 'grid_xt'
    dimension_name_y = 'latitude' if model_name == 'CERES' else 'grid_yt'

    assert dimension_name_x in snapshot.dims
    assert dimension_name_y in snapshot.dims

    # Generate storm-centered coordinates.
    # This will remove dependence on global coordinates and allow for storm-centered compositing.
    arr_x = snapshot[dimension_name_x] - snapshot['center_lon']
    arr_y = snapshot[dimension_name_y] - snapshot['center_lat']

    # Generate storm ID list to serve as an axis for the storm identifier that will enable xArray-based compositing
    storm_ID = [snapshot.attrs['storm_id']]
    center_lat = [snapshot['center_lat'].item()]
    # Expand dimensions of 2D data for the storm ID axis
    snapshot_dataset = np.expand_dims(np.expand_dims(
        snapshot[field_name].data, axis=0), axis=0)

    # Generate the xArray DataArray
    snapshot_dataset = xr.DataArray(data=snapshot_dataset,
                                    dims=['storm_id', 'grid_yt',
                                          'grid_yt_TC', 'grid_xt_TC'],
                                    coords={'grid_xt_TC': (['grid_xt_TC'], arr_x.data),
                                            'grid_yt_TC': (['grid_yt_TC'], arr_y.data),
                                            'storm_id': (['storm_id'], storm_ID),
                                            'center_lat': (['grid_yt'], center_lat)})

    if model_name == 'CERES':
        snapshot = snapshot.rename_dims(
            {dimension_name_x: 'grid_xt', dimension_name_y: 'grid_yt'})

    # To allow for the combination of different TCs, interpolate the storm-centered coordinates to common axis values based on the window size.
    arr_x_interp, arr_y_interp = storm_interpolation_grid_basis_vectors(
        snapshot, window_size)

    # Perform an interpolation from the storm-centered coordinates to the interpolated coordinates
    interpolated_snapshot_dataset = snapshot_dataset.interp(
        grid_xt_TC=arr_x_interp).interp(grid_yt_TC=arr_y_interp)

    return interpolated_snapshot_dataset


def get_local_time(timestamp: cftime.datetime,
                   latitude: float,
                   longitude: float,
                   method: str = 'simple',
                   diagnostic: bool = False):
    '''
    Methods to get local time.
    1. 'timezonefinder': use political timezone boundaries.
    2. 'simple' (default): divide the globe into equidistant 15-degree bins, one for each hour of offset.
    '''

    if diagnostic:
        print(f'[get_local_time()] Timestamp value input: {timestamp}.')

    # Account for International Date Line
    longitude = longitude if longitude <= 180 else longitude - 360
    # Get timestamp
    timestamp_datetime_GMT = datetime.datetime(
        year=timestamp.year, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)

    if method == 'timezonefinder':
        # Get timezone designator
        timezone_str = timezonefinder.TimezoneFinder(
        ).timezone_at(lng=longitude, lat=latitude)
        timezone = pytz.timezone(timezone_str)
        timestamp_datetime_GMT = datetime.datetime(
            year=timestamp.year, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        timestamp_datetime_local = timestamp_datetime_GMT + \
            timezone.utcoffset(timestamp_datetime_GMT)
        timestamp_cftime_local = cftime.datetime(year=timestamp_datetime_local.year,
                                                 month=timestamp_datetime_local.month,
                                                 day=timestamp_datetime_local.day,
                                                 hour=timestamp_datetime_local.hour)
        return timestamp_cftime_local

    elif method == 'simple':
        number_of_degrees = 360
        number_of_hours = 24
        degrees_per_hour = number_of_degrees / number_of_hours
        timezone_offset_hours = np.floor(longitude / degrees_per_hour)

        timestamp_datetime_local = timestamp_datetime_GMT + \
            datetime.timedelta(hours=timezone_offset_hours)
        timestamp_cftime_local = cftime.datetime(year=timestamp_datetime_local.year,
                                                 month=timestamp_datetime_local.month,
                                                 day=timestamp_datetime_local.day,
                                                 hour=timestamp_datetime_local.hour)
        return timestamp_cftime_local

    else:
        return "Timezone not found."


def hour_filter(composite_snapshots: dict,
                start_hour: int = 6,
                end_hour: int = 18) -> dict:
    ''' Method to filter TC samples by the hour at which they were sampled. '''

    assert start_hour != end_hour, 'Hours must be different.'
    assert isinstance(start_hour, int) & isinstance(
        end_hour, int), 'Hours must be integers.'

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

    storm_sampling_hour = {configuration: []
                           for configuration in composite_snapshots.keys()}
    for configuration in composite_snapshots.keys():
        for storm_ID in range(len(composite_snapshots[configuration]['storm_id'].values)):
            storm_sampling_hour[configuration].append(
                composite_snapshots[configuration].isel(storm_id=storm_ID)['local_time'].dt.hour)

    hour_bins = range(24)
    fig, ax = plt.subplots(figsize=(4, 2))
    for configuration in composite_snapshots.keys():
        number_of_samples = len(storm_sampling_hour[configuration])
        ax.hist(storm_sampling_hour[configuration], bins=hour_bins,
                label=f'{configuration}: N = {number_of_samples}', histtype='step', lw=2)
    fig.tight_layout()
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1, 1.1))
    ax.set_xlabel('Local hour of day')


def get_intensity_rate_of_change(snapshot: xr.Dataset,
                                 snapshot_timestamp: cftime,
                                 intensity_parameter: str = 'min_slp') -> float:
    ''' 
    Get rate of intensity change with time.
    1. Get the index of the snapshot timestamp nearest to the requested timestamp along the time axis
    3. Pull the indices immediately preceding and succeeding the snapshot timestamp index
    4. Sample the intensity parameter (I) and timestamp (t) at each point
    5. Get the rate of intensity change (dI) with respect to time (dt), dI/dt.
    '''

    assert intensity_parameter in [
        'min_slp', 'max_wind'], f"Intensity parameter provided is {intensity_parameter}, but must be 'min_slp' or 'max_wind'."

    # 1. Get the index of the snapshot timestamp nearest to the requested timestamp along the time axis
    snapshot_time_index = np.where(snapshot.sel(
        time=snapshot_timestamp, method='nearest').time.values == snapshot.indexes['time'])[0].item()

    # 2. Check if time index is at axis extremum (first or last element).
    # If yes, use 1-sided rate of change calculation. Else, use centered difference.
    is_index_first, is_index_last = False, False
    if snapshot_time_index == 0:
        is_index_first = True
    elif snapshot_time_index == len(snapshot.time.values) - 1:
        is_index_last = True

    # 3. Pull the indices immediately preceding and succeeding the snapshot timestamp index
    if is_index_first:
        start_index = snapshot_time_index
        end_index = snapshot_time_index + 1
    elif is_index_last:
        start_index = snapshot_time_index - 1
        end_index = snapshot_time_index
    else:
        start_index = snapshot_time_index - 1
        end_index = snapshot_time_index + 1

    # 4. Sample the intensity parameter (I) and timestamp (t) at each point.
    # Note: dt is originally computed in units of nanoseconds and is converted to hours.
    dI = (snapshot[intensity_parameter].isel(time=end_index) -
          snapshot[intensity_parameter].isel(time=start_index)).item()
    dt = (snapshot['time'].isel(time=end_index) -
          snapshot['time'].isel(time=start_index)).item() / 1e9 / 3600

    # 5. Get the rate of intensity change (dI) with respect to time (dt), dI/dt.
    dI_dt = dI/dt

    return dI_dt


def generate_TC_snapshot(model_name: str,
                         field_name: str,
                         window_size: int,
                         pressure_level: int | float,
                         snapshot_entry: dict,) -> pd.DataFrame:
    ''' Filter and select snapshots from track data. '''

    # Get iterand pathname
    snapshot_pathname = snapshot_entry['pathname']

    # Pull data for the TC at the iterand timestamp
    snapshot = xr.open_dataset(snapshot_pathname, use_cftime=True)

    # Get specific pressure level for snapshot if a vertical dimension (`pfull`) exists
    snapshot = snapshot.sel(
        pfull=pressure_level, method='nearest') if 'pfull' in snapshot.dims else snapshot

    # Acquire and process the snapshot timestamp
    snapshot_timestamp = get_snapshot_timestamp(snapshot_entry, snapshot)
    # Get local time corresponding to timestamp
    snapshot_local_time = get_local_time(snapshot_timestamp,
                                         snapshot.sel(time=snapshot_timestamp, method='nearest')[
                                             'center_lat'].item(),
                                         snapshot.sel(time=snapshot_timestamp, method='nearest')['center_lon'].item())

    # Get rate of intensity change with time
    intensity_rate_parameter = 'min_slp'
    intensity_rate_of_change = get_intensity_rate_of_change(snapshot=snapshot,
                                                            snapshot_timestamp=snapshot_timestamp,
                                                            intensity_parameter=intensity_rate_parameter)

    # Select the converted timestamp and remove null data
    snapshot = snapshot.sel(time=snapshot_timestamp, method='nearest').dropna(
        dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')

    # Get derived fields
    snapshot = derived_quantities(model_name, snapshot)
    # Perform grid redefinition and associated spatial interpolation
    interpolated_snapshot = storm_centered_interpolation(snapshot,
                                                         field_name=field_name,
                                                         window_size=window_size)

    interpolated_snapshot['local_time'] = snapshot_local_time
    interpolated_snapshot['min_slp'] = snapshot['min_slp'].item()
    interpolated_snapshot['dI_dt'] = intensity_rate_of_change

    return interpolated_snapshot


def generate_TC_composite(model_name: str,
                          experiment_name: str,
                          year_range: tuple[int, int],
                          field_name: str,
                          basin_name: str = 'global',
                          intensity_parameter: str = 'min_slp',
                          intensity_range: tuple[int | float, int | float] = (
                              0, np.inf),
                          number_of_snapshots: int = 10,
                          composite_window_size: int = 10,
                          pressure_level: int | float = 850,
                          TC_source_dirname: str='/projects/GEOCLIM/gr7610/analysis/TC-AQP/data/individual_TCs',
                          parallel: bool = False) -> tuple[xr.Dataset, dict]:
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
                                                 field_name,
                                                 basin_name,
                                                 pressure_level,
                                                 track_dataset,
                                                 TC_source_dirname=TC_source_dirname)

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
                                                     composite_window_size,
                                                     pressure_level)

    # Gather TC samples for compositing
    if parallel:
        # Maximum number of processors for computation
        max_number_procs = 24
        # Specify number of processors to use
        number_procs = len(snapshot_entries) if len(
            snapshot_entries) < max_number_procs else max_number_procs

        with Pool(processes=number_procs) as pool:
            snapshot_container = pool.map(
                partial_generate_TC_snapshot, snapshot_entries)
            print(
                f'Number of snapshots in container: {len(snapshot_container)}')
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
    timestamp_day_of_year = timestamp.dayofyr if 'cftime' in str(
        type(timestamp)) else timestamp.dt.dayofyear
    # Dataset time array days of year (handled differently by time object type)
    dataset_day_of_year = dataset.time.dt.dayofyear if 'cftime' in str(
        type(timestamp)) else dataset.time.dt.dayofyear
    # Get start and end days of year
    start_day_of_year, end_day_of_year = timestamp_day_of_year - \
        window_day_size, timestamp_day_of_year + window_day_size
    # Mask by window from start_day_of_year to end_day_of_year
    window = (dataset_day_of_year >= start_day_of_year) & (
        dataset_day_of_year <= end_day_of_year)
    # Mask data by the window
    dataset_window = dataset.sel(time=window)

    return dataset_window


def sample_GCM_constructor(snapshot: xr.Dataset,
                           GCM_snapshot: xr.Dataset,
                           field_name: str,
                           model_name: str = None,
                           window_size: int = 10):
    ''' Modify GCM data into a TC-centered dataset. '''

    # Add center coordinates to the GCM dataset
    GCM_snapshot['center_lon'] = snapshot['center_lon']
    GCM_snapshot['center_lat'] = snapshot['center_lat']
    GCM_snapshot.attrs = snapshot.attrs
    # Perform storm-centered interpolation for GCM data
    interpolated_GCM_snapshot = storm_centered_interpolation(snapshot=GCM_snapshot,
                                                             field_name=field_name,
                                                             window_size=window_size,
                                                             model_name=model_name)

    return interpolated_GCM_snapshot


def get_sample_GCM_data(model_name: str,
                        experiment_name: str,
                        field_name: str,
                        year_range: tuple[int, int],
                        pressure_level: int | float | None,
                        snapshot_timestamp: pd.Timestamp,
                        longitude: int | float,
                        latitude: int | float,
                        window_size: int,
                        sampling_day_window: int = 5):
    ''' Method to pull GCM data corresponding to a given TC snapshot. '''
    # Define dimension names based on the model
    dimension_name_x = 'longitude' if model_name == 'CERES' else 'grid_xt'
    dimension_name_y = 'latitude' if model_name == 'CERES' else 'grid_yt'

    # Construct field dictionary for postprocessed data loading
    # See `utilities.postprocessed_data_load` for details.
    # Note: this currently only supports single-surface atmospheric data
    pressure_level = str(
        pressure_level) if pressure_level is not None else pressure_level
    field_dictionary = {field_name: {
        'domain': 'atmos', 'level': pressure_level}}
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
                                                        month_range=(
                                                            snapshot_month - 1, snapshot_month + 1),
                                                        load_full_time=True, diagnostic=False)[model_name][experiment_name]

    # Define spatial extent for sample clipping
    grid_xt_extent = slice(longitude - window_size, longitude + window_size)
    grid_yt_extent = slice(latitude - window_size, latitude + window_size)
    # Trim the data spatially
    sample_GCM_data_filtered_space = sample_GCM_data.sortby(dimension_name_y).sel(
        {dimension_name_x: grid_xt_extent}).sel({dimension_name_y: grid_yt_extent})
    # Subsample over the time window specified: (iterand timestamp - sampling_day_window) to (iterand_timestamp + sampling_day_window)
    sample_GCM_data_filtered_time = get_time_window(
        sample_GCM_data_filtered_space, snapshot_timestamp, sampling_day_window)

    # Average in time
    sample_GCM_data_filtered = sample_GCM_data_filtered_time.mean(dim='time')

    return sample_GCM_data_filtered


def generate_composite_climatology_sample(model_name: str,
                                          experiment_name: str,
                                          field_name: str,
                                          year_range: tuple[int, int],
                                          pressure_level: int | float,
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

    # Get specific pressure level for snapshot if a vertical dimension (`pfull`) exists
    if 'pfull' in snapshot.dims:
        snapshot = snapshot.sel(pfull=pressure_level, method='nearest')
    else:
        pressure_level = None

    # Acquire and process the snapshot timestamp
    snapshot_timestamp = get_snapshot_timestamp(snapshot_entry, snapshot)
    # Select the converted timestamp and remove null data
    snapshot = snapshot.sel(time=snapshot_timestamp, method='nearest').dropna(
        dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')

    # Get sample TC center coordinates
    snapshot_center_longitude = snapshot['center_lon'].item()
    snapshot_center_latitude = snapshot['center_lat'].item()

    # Load GCM data corresponding to the given TC sample
    sample_GCM_data = get_sample_GCM_data(model_name,
                                          experiment_name,
                                          field_name,
                                          year_range,
                                          pressure_level,
                                          snapshot_timestamp,
                                          snapshot_center_longitude,
                                          snapshot_center_latitude,
                                          window_size)

    # Construct TC-centered GCM xarray object
    sample_GCM_data = sample_GCM_constructor(
        snapshot, sample_GCM_data, field_name, model_name=model_name)

    return sample_GCM_data


def generate_composite_climatology(model_name: str,
                                   experiment_name: str,
                                   year_range: tuple[int, int],
                                   field_name: str,
                                   pressure_level: int | float,
                                   snapshot_entries: dict,
                                   window_size: int = 12,
                                   troubleshooting: bool = False,
                                   parallel: bool = False):
    ''' Gather climatological samples. '''
    # Initialize list to contain all individual TCs samples
    GCM_sample_container = []
    # Generate partial function
    partial_generate_composite_climatology_sample = functools.partial(generate_composite_climatology_sample,
                                                                      model_name,
                                                                      experiment_name,
                                                                      field_name,
                                                                      year_range,
                                                                      pressure_level,
                                                                      window_size,
                                                                      troubleshooting)
    # Perform loading in parallel
    if parallel:
        ''' Offload TC-specific data generation onto parallel processes. '''
        # Maximum number of processors for computation
        max_number_procs = 24
        # Specify number of processors to use
        number_procs = len(snapshot_entries) if len(
            snapshot_entries) < max_number_procs else max_number_procs

        with Pool(processes=number_procs) as pool:
            GCM_sample_container = pool.map(
                partial_generate_composite_climatology_sample, snapshot_entries)
            pool.close()
    # Perform it serially
    else:
        for snapshot_entry in snapshot_entries:
            # Construct TC-centered GCM xarray object
            sample_GCM_data = partial_generate_composite_climatology_sample(
                snapshot_entry)
            # Append to container list
            GCM_sample_container.append(sample_GCM_data)

    # Concatenate all samples
    composite_GCM_samples = xr.concat(GCM_sample_container, dim='storm_id')

    return composite_GCM_samples


def get_composites(compositing_mode: str,
                   model_name: str,
                   experiment_name: str,
                   year_range: tuple[int | float, int | float],
                   field_name: str,
                   basin_name: str = 'global',
                   intensity_parameter: str = 'min_slp',
                   intensity_range: tuple[int | float,
                                          int | float] = (0, np.inf),
                   pressure_level: int | float | None = 850,
                   number_of_snapshots: int = -1,
                   TC_source_dirname: str='/projects/GEOCLIM/gr7610/analysis/TC-AQP/data/individual_TCs',
                   parallel: bool = False) -> xr.Dataset | tuple[xr.Dataset, xr.Dataset]:

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
                                                          pressure_level=pressure_level,
                                                          TC_source_dirname=TC_source_dirname,
                                                          parallel=parallel)

        print(
            f'Time elapsed for TCs: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(composite_TC_snapshots.storm_id.values)):.2f} s.')

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
                                                                            pressure_level=pressure_level,
                                                                            TC_source_dirname=TC_source_dirname,
                                                                            parallel=parallel)

        print(
            f'Time elapsed for TCs: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(TC_snapshot_entries)):.2f} s.')
        start_time = time.time()

        # Obtain GCM composites for times and locations corresponding to provided TC storm IDs
        composite_GCM_snapshots = generate_composite_climatology(model_name,
                                                                 experiment_name,
                                                                 year_range,
                                                                 field_name,
                                                                 pressure_level,
                                                                 TC_snapshot_entries,
                                                                 parallel=parallel)

        print(
            f'Time elapsed for GCM data: {(time.time() - start_time):.2f} s; per snapshot: {((time.time() - start_time)/len(TC_snapshot_entries)):.2f} s.')

        return (composite_TC_snapshots, composite_GCM_snapshots)


def save_composites(composite_TC_samples: dict,
                    configuration_name: str,
                    composite_type: str,
                    field_name: str,
                    year_range: list[int, int] | tuple[int, int],
                    intensity_range: list[int, int] | tuple[int, int],
                    basin_name: str,
                    pressure_level: int | float | None,
                    dirname: str):

    assert composite_type in ['TC', 'GCM'], 'Composite mode must be specified as `TC` or `GCM`.'

    pressure_level_substr = f'.pressure_level.{pressure_level:0d}' if pressure_level is not None else ''
    filename = f'{composite_type}.configuration.{configuration_name}.field_name.{field_name}{pressure_level_substr}.year_range.{min(year_range)}_{max(year_range)}.intensity_range.{min(intensity_range)}_{max(intensity_range)}.basin.{basin_name}.nc'
    pathname = os.path.join(dirname, filename)

    print(f'Saving composite data to {pathname}...')

    composite_TC_samples.to_dataset(name=field_name).to_netcdf(pathname)
    # if os.path.exists(pathname):
    #     override = input(
    #         'File already exists - do you wish to override? [y/n]\n')
    #     if 'y' in override:
    #         os.remove(pathname)
    #         composite_TC_samples.to_dataset(
    #             name=field_name).to_netcdf(pathname)
    #     else:
    #         import sys
    #         print(f'File at {pathname} not overridden. Exiting...')
    #         sys.exit()
    # else:
    #     composite_TC_samples.to_dataset(name=field_name).to_netcdf(pathname)


def main(configuration_name: str,
         field_name: str,
         year_range: list[int | float, int | float] | tuple[int | float, int | float],
         intensity_parameter: str = 'min_slp',
         intensity_range: list[int | float, int | float] | tuple[int | float, int | float] = (0, 1000),
         number_of_snapshots: int = 500,
         compositing_mode: str = 'TC',
         parallel: bool = False,
         save_data: bool = False,
         TC_source_dirname: str='/projects/GEOCLIM/gr7610/analysis/TC-AQP/data/individual_TCs',
         storage_dirname: str = '/projects/GEOCLIM/gr7610/analysis/TC-RAD/data'):

    basin_name = 'global'
    pressure_level = None

    model_name = configuration_name.split('-')[0]
    experiment_name = configuration_name.split('-')[1]

    composite_TC_samples, composite_GCM_samples = {}, {}

    print(
        f'Processing composites for {model_name}, {experiment_name} over years {year_range}...')

    if compositing_mode == 'TC':
        composite_TC_samples = get_composites(compositing_mode,
                                    model_name,
                                    experiment_name,
                                    year_range,
                                    field_name,
                                    basin_name=basin_name,
                                    intensity_parameter=intensity_parameter,
                                    intensity_range=intensity_range,
                                    number_of_snapshots=number_of_snapshots,
                                    pressure_level=pressure_level,
                                    TC_source_dirname=TC_source_dirname,
                                    parallel=parallel)

        if savefig:
            save_composites(composite_TC_samples,
                            configuration_name,
                            field_name=field_name,
                            year_range=year_range,
                            intensity_range=intensity_range,
                            basin_name=basin_name,
                            pressure_level=pressure_level,
                            composite_type=compositing_mode,
                            dirname=storage_dirname)
        else:
            return composite_TC_samples

    else:
        TMP = get_composites(compositing_mode,
                   model_name,
                   experiment_name,
                   year_range,
                   field_name,
                   basin_name=basin_name,
                   intensity_parameter=intensity_parameter,
                   intensity_range=intensity_range,
                   number_of_snapshots=number_of_snapshots,
                   pressure_level=pressure_level,
                   TC_source_dirname=TC_source_dirname,
                   parallel=parallel)

        composite_TC_samples, composite_GCM_samples = TMP

        # Interpolate the climate field data for non-GCM data, which is typically coarser than GCM data
        # For example, CERES data is available at ~1 deg resolution
        if model_name in ['CERES']:
            composite_GCM_samples = composite_GCM_samples.interp(grid_xt_TC=composite_TC_samples.grid_xt_TC,
                                                                 grid_yt_TC=composite_TC_samples.grid_xt_TC)

        if save_data:
            save_composites(composite_TC_samples,
                            configuration_name,
                            composite_type='TC',
                            field_name=field_name,
                            year_range=year_range,
                            intensity_range=intensity_range,
                            basin_name=basin_name,
                            pressure_level=pressure_level,
                            dirname=storage_dirname)
            save_composites(composite_GCM_samples,
                            configuration_name,
                            composite_type='GCM',
                            field_name=field_name,
                            year_range=year_range,
                            intensity_range=intensity_range,
                            basin_name=basin_name,
                            pressure_level=pressure_level,
                            dirname=storage_dirname)


''' 
Sample run script.

configuration_name = 'HIRAM-CTL1990.0N'
number_of_snapshots = -1
year_range = (2, 4)

for field_name in ['vort850']:
    for intensity_range in [(0, 1000)]:
        
        print(f'Processing field name {field_name} over intensity range {intensity_range}...')
        
        run(configuration_name=configuration_name,
            field_name=field_name,
            year_range=year_range,
            intensity_range=intensity_range,
            number_of_snapshots=number_of_snapshots,
            compositing_mode='TC',
            parallel=True,
            save_data=True,
            storage_dirname='/projects/GEOCLIM/gr7610/analysis/TC-AQP/data')
'''
