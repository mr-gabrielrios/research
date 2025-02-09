''' 
Script name:
TC anomaly generation script
---
Objective: 
Generate anomalies from the climatology pertaining to a given TC for a given field name. 
Outputs from this script are netCDF files containing xArray Dataset with TC anomalies
from the mean. These data can be paired with the individual TC netCDFs from which
they were derived.
---
Definitions:

A TC is defined as a storm whose coordinates are given by GFDL QuickTracks,
and paired in space and time with the global climate model data QuickTracks 
is run on.

The climatology is defined as the mean over a given time period (typically 25 years).

An anomaly is defined as the difference between a TC at a given position and time
from the climatology at the same position and time. A window in both space and time
is defined. In space, some distance `window_size` is defined in units of degrees 
to provide a spatial extent outside (but including) the TC. In time, some number of
days is defined such that the days before and after a given timestamp are included in the
climatology.
'''

import xarray as xr
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt
import importlib
import cftime
import functools
import numpy as np
import os
import pandas as pd
import random
import time
from multiprocessing import Pool

import derived
import utilities
import visualization

def get_filename_year(filename: str) -> int:
    delimiter_string = 'storm_ID'
    filename_year = int(filename.split(f'{delimiter_string}-')[-1].split('-')[0])
    return filename_year

def get_filename_intensity(filename: str,
                           delimiter_string: str) -> int:
    filename_intensity = int(filename.split(f'{delimiter_string}-')[-1].split('.')[0])
    return filename_intensity

def field_correction(model_name: str, 
                     dataset: xr.DataArray,
                     field_name: str):

    if field_name in ['precip', 'evap', 'p-e']:
        # Conversion of total precipitation per hour to daily precipitation rate for ERA5 data, instantaneous rate to daily for GFDL GCM data
        factor = (1 / 3600) * 1000 * 86400 if model_name == 'ERA5' else 86400
    else:
        factor = 1
    dataset[field_name] = dataset[field_name] * factor

    return dataset

def derived_quantities(model_name: str,
                       dataset: xr.Dataset):

    # Correct sign convention from positive down to positive up for upwards-directed quantities
    dataset['lwup_sfc'] = dataset['lwup_sfc'] * -1 if model_name in ['ERA5'] else dataset['lwup_sfc']
    dataset['swup_sfc'] = dataset['swup_sfc'] * -1 if model_name in ['ERA5'] else dataset['swup_sfc']

    # Generate derived fields
    dataset = derived.TC_surface_wind_speed(dataset)
    dataset = derived.net_lw(dataset)
    dataset = derived.net_sw(dataset)
    dataset = derived.atmospheric_heating(dataset)

    return dataset

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

def grid_check(dataset: xr.Dataset) -> tuple[float, float]:

    ''' Perform grid checks and get grid spacing for a generic xArray coordinate system. '''

    # Ensure necessary basis vector dimensions are available
    assert ('grid_xt' in dataset.dims) and ('grid_yt' in dataset.dims)
    # Get differences in grid spacing along each vector
    d_grid_yt = dataset['grid_yt'].diff(dim='grid_yt')
    d_grid_xt = dataset['grid_xt'].diff(dim='grid_xt')
    # Ensure that the differences are equivalent for all indices to ensure equal spacing
    grid_tolerance = 1e-6
    assert sum(d_grid_xt.diff(dim='grid_xt') < grid_tolerance) == len(d_grid_xt) - 1, 'Grid is irregular along the `grid_xt` axis.'
    assert sum(d_grid_yt.diff(dim='grid_yt') < grid_tolerance) == len(d_grid_yt) - 1, 'Grid is irregular along the `grid_yt` axis.'
    # Get grid spacing along each direction
    dx, dy = d_grid_xt.isel(grid_xt=0).item(), d_grid_yt.isel(grid_yt=0).item()

    return abs(dx), abs(dy)

def grid_interpolation(working_grid: xr.Dataset,
                       reference_grid: xr.DataArray,
                       diagnostic: bool=False) -> xr.DataArray:

    ''' 
    Method to generate uniform interpolation basis vectors for a TC-centered grid. 
    Boolean `equal_number_of_points` is used to ensure equal grid point numbers, and is optional.
    '''

    # Ensure necessary basis vector dimensions are available
    assert ('grid_xt' in working_grid.dims) and ('grid_yt' in working_grid.dims)
    # Get differences in grid spacing along each vector
    d_grid_yt = reference_grid['grid_yt'].diff(dim='grid_yt')
    d_grid_xt = reference_grid['grid_xt'].diff(dim='grid_xt')
    # Ensure that the differences are equivalent for all indices to ensure equal spacing
    grid_tolerance = 1e-6
    assert sum(d_grid_xt.diff(dim='grid_xt') < grid_tolerance) == len(d_grid_xt) - 1, 'Grid is irregular along the `grid_xt` axis.'
    assert sum(d_grid_yt.diff(dim='grid_yt') < grid_tolerance) == len(d_grid_yt) - 1, 'Grid is irregular along the `grid_yt` axis.'
    # Get grid spacing along each direction
    dx, dy = d_grid_xt.isel(grid_xt=0).item(), d_grid_yt.isel(grid_yt=0).item()
    
    # Padding on window search for storm coordinate grid
    # This expands the window `padding_factor` grid cells in each direction of the window extent for a given storm timestamp
    padding_factor = 0
    
    # Get extent of longitudes
    minimum_longitude, maximum_longitude = [reference_grid['grid_xt'].min().item() - dx * padding_factor, 
                                            reference_grid['grid_xt'].max().item() + dx * padding_factor]
    # Get extent of latitudes
    minimum_latitude, maximum_latitude = [reference_grid['grid_yt'].min().item() - dy * padding_factor, 
                                          reference_grid['grid_yt'].max().item() + dy * padding_factor]

    # Round values due to weird GFDL GCM output behavior
    minimum_longitude, maximum_longitude = [np.round(minimum_longitude, decimals=4),
                                            np.round(maximum_longitude, decimals=4)]
    minimum_latitude, maximum_latitude = [np.round(minimum_latitude, decimals=4),
                                          np.round(maximum_latitude, decimals=4)]
    dx, dy = [np.round(dx, decimals=4),
              np.round(dy, decimals=4)]

    if diagnostic:
        print(f'[grid_interpolation()] Window extent: longitudes = {(minimum_longitude, maximum_longitude)} and latitudes= {(minimum_latitude, maximum_latitude)}')
        print(f'[grid_interpolation()] dx = {dx}; dy = {dy}')

    # Construct interpolation arrays for interpolating the area grid onto the storm grid
    interpolation_array_x = np.arange(minimum_longitude, maximum_longitude + dx, dx)
    interpolation_array_y = np.arange(minimum_latitude, maximum_latitude + dy, dy)
    
    if diagnostic:
        print(f'[grid_interpolation()] interpolated latitudes: {interpolation_array_y}')
    
    # Perform interpolation of area grid onto the storm grid
    interpolated_working_grid = working_grid.interp(grid_xt=interpolation_array_x).interp(grid_yt=interpolation_array_y)

    return interpolated_working_grid

def get_surface_area(storm_dataset_timestamp: xr.DataArray,
                     diagnostic: bool=False) -> xr.DataArray:

    # Load surface area data from GFDL GCM output
    surface_area = xr.open_dataset('/projects/GEOCLIM/gr7610/tools/AM2.5_atmos_area.nc')['__xarray_dataarray_variable__']

    # Shed null values
    storm_dataset_timestamp = storm_dataset_timestamp.dropna('grid_xt', how='all').dropna('grid_yt', how='all')
    # Get minimum and maximum spatial extent values
    minimum_longitude, maximum_longitude = [storm_dataset_timestamp['grid_xt'].min().item(), 
                                            storm_dataset_timestamp['grid_xt'].max().item()]
    minimum_latitude, maximum_latitude = [storm_dataset_timestamp['grid_yt'].min().item(), 
                                          storm_dataset_timestamp['grid_yt'].max().item()]
    # Round values due to weird GFDL GCM output behavior
    minimum_longitude, maximum_longitude = [np.round(minimum_longitude, decimals=4),
                                            np.round(maximum_longitude, decimals=4)]
    minimum_latitude, maximum_latitude = [np.round(minimum_latitude, decimals=4),
                                          np.round(maximum_latitude, decimals=4)]
    # Get surface area at iterand timestamp
    surface_area_timestamp = surface_area.sel(grid_xt=slice(minimum_longitude, maximum_longitude),
                                              grid_yt=slice(minimum_latitude, maximum_latitude))

    # Interpolate area onto storm coordinates
    interpolated_surface_area_timestamp = grid_interpolation(surface_area_timestamp, storm_dataset_timestamp)

    if diagnostic:
        print('------------------------------------')
        print(f'Window extent: longitudes = {(minimum_longitude, maximum_longitude)} and latitudes= {(minimum_latitude, maximum_latitude)}')
        print(f'[get_surface_area()]: Storm dataset latitudes:\n{storm_dataset_timestamp.grid_yt.values}')
        print(f'[get_surface_area()]: Surface area dataset latitudes:\n{interpolated_surface_area_timestamp.grid_yt.values}')

    # Make sure all longitudes and latitudes are within some tolerance of each other
    assert np.allclose(storm_dataset_timestamp.grid_xt, interpolated_surface_area_timestamp.grid_xt), f'\nData longitudes:\n{storm_dataset_timestamp.grid_xt.values}; \n Surface area longitudes:\n{surface_area_timestamp.grid_xt.values}'
    assert np.allclose(storm_dataset_timestamp.grid_yt, interpolated_surface_area_timestamp.grid_yt), f'\nData latitude:\n{storm_dataset_timestamp.grid_yt.values}; \n Surface area latitudes:\n{surface_area_timestamp.grid_yt.values}'

    return interpolated_surface_area_timestamp

def get_surface_area(storm_dataset_timestamp: xr.DataArray,
                     diagnostic: bool=False) -> xr.DataArray:

    # Load surface area data from GFDL GCM output
    surface_area = xr.open_dataset('/projects/GEOCLIM/gr7610/tools/AM2.5_atmos_area.nc')['__xarray_dataarray_variable__']

    # Shed null values
    storm_dataset_timestamp = storm_dataset_timestamp.dropna('grid_xt', how='all').dropna('grid_yt', how='all')
    # Get minimum and maximum spatial extent values
    minimum_longitude, maximum_longitude = [storm_dataset_timestamp['grid_xt'].min().item(), 
                                            storm_dataset_timestamp['grid_xt'].max().item()]
    minimum_latitude, maximum_latitude = [storm_dataset_timestamp['grid_yt'].min().item(), 
                                          storm_dataset_timestamp['grid_yt'].max().item()]
    # Round values due to weird GFDL GCM output behavior
    minimum_longitude, maximum_longitude = [np.round(minimum_longitude, decimals=4),
                                            np.round(maximum_longitude, decimals=4)]
    minimum_latitude, maximum_latitude = [np.round(minimum_latitude, decimals=4),
                                          np.round(maximum_latitude, decimals=4)]
    # Get surface area at iterand timestamp
    surface_area_timestamp = surface_area.sel(grid_xt=slice(minimum_longitude, maximum_longitude),
                                              grid_yt=slice(minimum_latitude, maximum_latitude))

    # Interpolate area onto storm coordinates
    interpolated_surface_area_timestamp = grid_interpolation(surface_area_timestamp, storm_dataset_timestamp)

    if diagnostic:
        print('------------------------------------')
        print(f'Window extent: longitudes = {(minimum_longitude, maximum_longitude)} and latitudes= {(minimum_latitude, maximum_latitude)}')
        print(f'[get_surface_area()]: Storm dataset latitudes:\n{storm_dataset_timestamp.grid_yt.values}')
        print(f'[get_surface_area()]: Surface area dataset latitudes:\n{interpolated_surface_area_timestamp.grid_yt.values}')

    # Make sure all longitudes and latitudes are within some tolerance of each other
    assert np.allclose(storm_dataset_timestamp.grid_xt, interpolated_surface_area_timestamp.grid_xt), f'\nData longitudes:\n{storm_dataset_timestamp.grid_xt.values}; \n Surface area longitudes:\n{surface_area_timestamp.grid_xt.values}'
    assert np.allclose(storm_dataset_timestamp.grid_yt, interpolated_surface_area_timestamp.grid_yt), f'\nData latitude:\n{storm_dataset_timestamp.grid_yt.values}; \n Surface area latitudes:\n{surface_area_timestamp.grid_yt.values}'

    return interpolated_surface_area_timestamp

def get_sample_GCM_data(model_name: str,
                        experiment_name: str,
                        field_name: str,
                        year_range: tuple[int, int],
                        sampling_timestamp: pd.Timestamp,
                        longitude: int|float,
                        latitude: int|float,
                        time_window_size: int|float,
                        space_window_size: int|float,
                        diagnostic: bool=False):

    ''' Method to pull GCM data corresponding to a given TC snapshot. '''

    # Construct field dictionary for postprocessed data loading
    # See `utilities.postprocessed_data_load` for details.
    # Note: this currently only supports single-surface atmospheric data
    field_dictionary = {field_name: {'domain': 'atmos', 'level': None}}
    # Extract month from the iterand timestamp to perform initial climatology filtering
    sampling_year, sampling_month, sampling_day = [sampling_timestamp.year,
                                                   sampling_timestamp.month,
                                                   sampling_timestamp.day,]
    # Load the data and pathnames
    sample_GCM_data, sample_GCM_path = utilities.postprocessed_data_load(model_name,
                                                                        experiment_name,
                                                                        field_dictionary,
                                                                        year_range,
                                                                        data_type='mean_daily',
                                                                        month_range=(sampling_month, sampling_month),
                                                                        load_full_time=True,
                                                                        return_paths=True)
    # Subfilter by model name and experiment
    sample_GCM_data = sample_GCM_data[model_name][experiment_name]
    sample_GCM_path = sample_GCM_path[model_name][experiment_name][field_name]
    # Get GCM grid spacing
    GCM_dx, GCM_dy = grid_check(sample_GCM_data)
    # Define spatial extent for sample clipping
    grid_xt_extent = slice(longitude - space_window_size, longitude + space_window_size)
    grid_yt_extent = slice(latitude - space_window_size, latitude + space_window_size)
    # Trim the data spatially
    sample_GCM_data_filtered_space = sample_GCM_data.sortby('grid_yt').sel(grid_xt=grid_xt_extent).sel(grid_yt=grid_yt_extent)
    # Subsample over the time window specified: (iterand timestamp - sampling_day_window) to (iterand_timestamp + sampling_day_window)
    sample_GCM_data_filtered_time = get_time_window(sample_GCM_data_filtered_space, sampling_timestamp, time_window_size)
    # Average in time
    sample_GCM_data_filtered = sample_GCM_data_filtered_time.mean(dim='time')
    # Append time coordinate to GCM sample
    sample_GCM_data_filtered['time'] = sampling_timestamp
    # Add reference climatology filename to attributes
    sample_GCM_data_filtered.attrs['climatology_file'] = sample_GCM_path
    
    if diagnostic:
        print(f'Storm timestamp center: {longitude}, {latitude}.')
        print(f'GCM grid spacing: dx = {GCM_dx}, dy = {GCM_dy}.')
        print(f'Storm timestamp extent = longitude: {grid_xt_extent}, latitude: {grid_yt_extent}.')
        print(f'GCM extent = longitude: {sample_GCM_data.grid_xt.values}, latitude: {sample_GCM_data.grid_yt.values}.')
        print(f'Filtered GCM extent = longitude: {sample_GCM_data_filtered.grid_xt.values}, latitude: {sample_GCM_data_filtered.grid_yt.values}.')

    return sample_GCM_data_filtered

def get_TC_anomaly_timestamp(model_name: str,
                             experiment_name: str,
                             year_range: tuple[int, int],
                             storm_reanalysis_data: xr.Dataset,
                             field_name: str,
                             time_window_size: int | float,
                             space_window_size: int | float,
                             sampling_timestamp: cftime.datetime,
                             diagnostic: bool=False):
    
    # Get timestamp
    storm_sample = storm_reanalysis_data.sel(time=sampling_timestamp)
    sampling_timestamp = storm_sample.time.item()

    # Get sample TC month and day
    sample_month = sampling_timestamp.month
    sample_day = sampling_timestamp.day

    # Get sample TC center coordinates
    sample_center_longitude = storm_sample['center_lon'].item()
    sample_center_latitude = storm_sample['center_lat'].item()

    # Load GCM data according to the given sample
    sample_GCM_data = get_sample_GCM_data(model_name, 
                                          experiment_name,
                                          field_name,
                                          year_range,
                                          sampling_timestamp,
                                          sample_center_longitude,
                                          sample_center_latitude,
                                          time_window_size=time_window_size,
                                          space_window_size=space_window_size)

    # Interpolate the GCM data to the storm data
    sample_GCM_data = grid_interpolation(sample_GCM_data, storm_sample)

    # Get simple anomaly
    TC_climatological_anomaly_timestamp = storm_sample[field_name] - sample_GCM_data[field_name]

    storm_ID = storm_sample.attrs['storm_id']
    TC_climatological_anomaly_timestamp = xr.Dataset(data_vars={field_name: (['storm_id', 'grid_yt', 'grid_xt'], 
                                                                             np.expand_dims(TC_climatological_anomaly_timestamp.data, axis=0))},
                                                     coords={'storm_id': ('storm_id', [storm_ID]),
                                                             'grid_yt': ('grid_yt', TC_climatological_anomaly_timestamp.grid_yt.values),
                                                             'grid_xt': ('grid_xt', TC_climatological_anomaly_timestamp.grid_xt.values)},
                                                     attrs=sample_GCM_data.attrs)
    
    return TC_climatological_anomaly_timestamp

def get_TC_anomaly(model_name: str, 
                   experiment_name: str,
                   field_name: str,
                   year_range: tuple,
                   time_window_size: int|float,
                   space_window_size: int|float,
                   pathname: str,
                   parallel: bool):

    ''' Method to generate an anomaly from the climatology for a given TC and field name. '''

    # Load data and perform pre-processing
    storm_reanalysis_data = xr.open_dataset(pathname, use_cftime=True)
    storm_reanalysis_data = utilities.field_correction(model_name, storm_reanalysis_data)
    storm_reanalysis_data = derived_quantities(model_name, storm_reanalysis_data)
    
    # Log runtime for performance profiling
    start_time = time.time()

    # Define partial function for acquiring TC anomaly
    # This standardizes inputs and only leaves the storm timestamp as the free variable
    partial_TC_anomaly_timestamp = functools.partial(get_TC_anomaly_timestamp,
                                                     model_name,
                                                     experiment_name,
                                                     year_range,
                                                     storm_reanalysis_data,
                                                     field_name,
                                                     time_window_size,
                                                     space_window_size)
    # Pull storm timestamps to iterate over
    sampling_timestamps = storm_reanalysis_data.time.values
    # If chosen, run in parallel using the partial function
    if parallel:
        with Pool() as pool:
            TC_anomaly_timestamps = pool.map(partial_TC_anomaly_timestamp, sampling_timestamps)
            pool.close()
    # Else, run serial. Serial is usually better for troubleshooting and debugging.
    else:
        # Initialize container lists for appending outputs for each timestamp
        TC_anomaly_timestamps = []
        
        for sampling_timestamp in sampling_timestamps:
            TC_anomaly_timestamp = partial_TC_anomaly_timestamp(sampling_timestamp)
            TC_anomaly_timestamps.append(TC_anomaly_timestamp)
    # Concatenate all timestamps and ensure sorting by time
    TC_anomaly_dataset = xr.concat(TC_anomaly_timestamps, dim='time').sortby('time')
    
    # Append metadata corresponding to the TC
    TC_anomaly_dataset['center_lon'] = storm_reanalysis_data['center_lon']
    TC_anomaly_dataset['center_lat'] = storm_reanalysis_data['center_lat']
    TC_anomaly_dataset['max_wind'] = storm_reanalysis_data['max_wind']
    TC_anomaly_dataset['min_slp'] = storm_reanalysis_data['min_slp']
    TC_anomaly_dataset.attrs['storm_id'] = storm_reanalysis_data.attrs['storm_id']

    print(f"Elapsed time for {storm_reanalysis_data.attrs['storm_id']}: {(time.time() - start_time):.2f} s")
    
    return TC_anomaly_dataset

def save_TC_anomaly(pathname: str,
                    anomaly_dataset: xr.Dataset,
                    field_name: str,
                    time_window_size: int|float,
                    space_window_size: int|float):
    
    ''' Save TC anomaly as an xarray Dataset object encoded in a netCDF file. '''
    
    # Define directory to where netCDF files will be saved
    dirname = '/tigress/GEOCLIM/gr7610/analysis/tc_storage/individual_TC_anomalies'
    # Get filename from pathname
    TC_filename = os.path.basename(pathname)
    
    # Append anomaly-specific information to the pathname
    anomaly_filename = f'{TC_filename.split('.nc')[0]}.anomaly-{field_name}.anomaly_window_size-space_{space_window_size}-time_{time_window_size}.nc'
    # Create pathname
    anomaly_pathname = os.path.join(dirname, anomaly_filename)
    
    # Save Dataset to file
    print(f'Saving anomaly dataset for {field_name} with spatial window size {space_window_size} deg and temporal window size {time_window_size} days to file {anomaly_pathname}...')
    anomaly_dataset.to_netcdf(anomaly_pathname)

def main(model_name: str,
         experiment_name: str,
         year_range: tuple[int, int],
         field_name: str,
         basin_name: str='global',
         intensity_parameter: str='min_slp',
         intensity_range: tuple[int|float, int|float]=(0, np.inf),
         number_of_storms: int=1,
         time_window_size: int|float=3,
         space_window_size: int|float=10,
         parallel: bool=True):

    # Use TC characteristics to find corresponding saved TCs in defined directory
    dirname = '/tigress/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    # Generate list of filenames
    filenames = [filename for filename in os.listdir(dirname) if
                 model_name in filename and
                 experiment_name in filename and
                 min(year_range) <= get_filename_year(filename) < max(year_range) and
                 min(intensity_range) <= get_filename_intensity(filename, intensity_parameter) < max(intensity_range)]
    # Generate list of pathnames
    pathnames = [os.path.join(dirname, filename) for filename in filenames]
    # If a number of storms provided, sample at random
    pathnames = random.sample(pathnames, number_of_storms) if number_of_storms > 0 else pathnames
    # Filter pathnames by basin, if basin_name is not global
    if basin_name != 'global':
        pathnames = [pathname for pathname in pathnames if basin_name in pathname]
    # Stop the script if no files matching the input criteria are found.
    assert len(pathnames) > 0, f'No files found matching the criteria provided.'

    # Process each path
    for pathname in pathnames:
        print(f'Pulling data from filename {pathname}.')
        # Get anomaly for an individual TC
        TC_anomaly_dataset = get_TC_anomaly(model_name, 
                                            experiment_name,
                                            field_name,
                                            year_range,
                                            time_window_size,
                                            space_window_size,
                                            pathname,
                                            parallel=parallel)
        # Save the data to a custom location
        save_TC_anomaly(pathname=pathname,
                        anomaly_dataset=TC_anomaly_dataset,
                        field_name=field_name,
                        time_window_size=time_window_size,
                        space_window_size=space_window_size)

if __name__ == '__main__':
    model_name = 'FLOR'
    experiment_name = 'CTL1990s_FA_tiger3'
    year_range = (1901, 1902)
    field_name = 'WVP'
    intensity_parameter = 'min_slp'
    intensity_range = (970, 980)
        
    main(model_name=model_name,
         experiment_name=experiment_name,
         year_range=year_range,
         field_name=field_name,
         intensity_parameter=intensity_parameter,
         intensity_range=intensity_range)