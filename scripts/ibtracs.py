import numpy as np
import xarray as xr
import matplotlib, matplotlib.pyplot as plt
import cartopy, cartopy.crs as ccrs
import pandas as pd
import functools
import os
import cftime
import time

import TC_tracker
import visualization

import warnings
warnings.filterwarnings('ignore')

def load(dirname: str | None=None,
         filename: str | None=None) -> xr.Dataset:
    
    ''' Load an xArray Dataset from IBTrACS repository. '''
    
    # Generate path (use provided defaults if none are given)
    dirname = '/scratch/gpfs/GEOCLIM/gr7610/tiger3/reference/datasets/IBTrACS' if not dirname else dirname
    # Use data post-1980
    filename = 'IBTrACS.since1980.v04r01.nc' if not filename else filename
    pathname = os.path.join(dirname, filename)
    # Load the data
    dataset = xr.open_dataset(pathname)

    return dataset

def format_adjustment(dataset: xr.Dataset) -> xr.Dataset:
    
    ''' Adjust IBTrACS data to GFDL QuickTracks naming conventions. '''
    
    dataset = intake_filtering(dataset)
    dataset = unit_conversion(dataset).load()

    return dataset

def intake_filtering(data: xr.Dataset) -> xr.Dataset:
    ''' Filtering to perform on reading of initial loadingo of IBTrACS data. '''    
    
    # Keep useful variables
    data_variables = ['numobs', 'sid', 'season', 'number', 'basin', 'subbasin', 'name', 'iso_time', 'nature', 'dist2land', 'landfall', 'iflag', 'storm_speed', 'storm_dir', 'wmo_wind', 'wmo_pres']
    data = data[data_variables]
    
    # Remove quadrant data
    data_variables = [var for var in data.data_vars 
                      if 'quadrant' not in data[var].dims]
    data = data[data_variables]
    
    return data

def valid_entry(storm, check_variable: str='wmo_pres') -> bool:
    # Variable used to inspect for null values. 
    # If the number of nulls for the variable equals the number of entries, categorize storm as a bad entry.
    check_variable = 'wmo_pres'
    if sum(storm[check_variable].isnull()) == len(storm[check_variable]):
        return False
    else:
        return True
    
def unit_conversion(storm: xr.Dataset) -> xr.Dataset:
    storm['wmo_wind-ms'] = storm['wmo_wind'] / 1.944
    return storm

def clean_dataframe_IBTrACS(data_variables: dict,
                            dataset: pd.DataFrame) -> pd.DataFrame:

    ''' Method to perform data cleaning operations to IBTrACS data. '''
    
    # Round to nearest hour
    dataset['time'] = pd.to_datetime(dataset['time'].dt.floor('H'))
    # Rename columns
    dataset = dataset.rename(columns=data_variables)
    # Generate cftime object
    dataset['cftime'] = dataset['time'].apply(lambda x: cftime.datetime(year=x.year, month=x.month, day=x.day, hour=x.hour))
    # Change longitude range from (-180, 180) to (0, 360)
    dataset['center_lon'].loc[dataset['center_lon'] < 0] = dataset['center_lon'].loc[dataset['center_lon'] < 0] + 360

    return dataset

def intensity_categories(intensity_scale: str='SSHWS') -> pd.DataFrame:

    intensity_scales = {'SSHWS': {'TD': (0, 17),
                                  'TS': (18, 32),
                                  'C1': (33, 42),
                                  'C2': (43, 49),
                                  'C3': (50, 58),
                                  'C4': (59, 70),
                                  'C5': (71, np.inf)}}
    
    return intensity_scales[intensity_scale]

def generate_dataframe(dataset_xr: xr.core.dataset.Dataset) -> pd.DataFrame:

    # Define xArray data variable names to select, and corresponding DataFrame names for renaming
    data_variables = {'time': 'time', 
                      'lat': 'center_lat', 
                      'lon': 'center_lon', 
                      'wmo_wind-ms': 'max_wind', 
                      'wmo_pres': 'min_slp'}
    # Create partial function due to fixed input parameter
    dataframe_cleaner = functools.partial(clean_dataframe_IBTrACS, data_variables) 
    
    # Initialize empty container to store generated DataFrames for future concatenation
    storm_dataframes = []
    
    # Iterate over each strm
    for storm_ID in dataset_xr.storm.values:
        # Define shorthand for the iterand storm
        storm_dataset = dataset_xr.sel(storm=storm_ID)
        # Gather number of observations for index-based selection
        # In other words, select indices 0 to `number_of_observations` along the `date_time` axis
        number_of_observations = int(np.round(storm_dataset['numobs']))
        # Subselect data to place in DataFrame
        storm_dataframe = storm_dataset.isel(date_time=range(number_of_observations))[data_variables.keys()]
        # Convert object to DataFrame
        storm_dataframe = storm_dataframe.to_dataframe().reset_index(drop=True)
        # Append storm ID to DataFrame (assumes string is UTF-8 byte object)
        # Add a dash after the 4-digit year to conform with GFDL QuickTracks storm naming convention
        storm_ID = storm_dataset['sid'].item().decode('utf-8')
        storm_dataframe['storm_id'] = f'{storm_ID[:4]}-{storm_ID[4:]}'
        # Derive storm_duration
        storm_dataframe['duration'] = (storm_dataframe['time'].max() - storm_dataframe['time'].min()).total_seconds() / 86400
        # Perform data cleaning
        storm_dataframe = dataframe_cleaner(storm_dataframe)
        # Append to container
        storm_dataframes.append(storm_dataframe)
    # Concatenate into single DataFrame and drop repeat occurrences for each given TC
    storm_dataframes = pd.concat(storm_dataframes).drop_duplicates(subset=['storm_id', 'time'])

    return storm_dataframes

def intensity_categorization(dataset: pd.DataFrame,
                             intensity_scale: str='SSHWS'):

    container = []
    intensity_category_bins = intensity_categories(intensity_scale)
    
    # Get bins for grouping
    cut_bins = np.array([min(intensity_bin) for intensity_bin in sorted(intensity_category_bins.values())])
    # Append maximum value to ensure most intense bin is captured
    maximum_bin_value = max([max(intensity_bin) for intensity_bin in sorted(intensity_category_bins.values())])
    cut_bins = np.append(cut_bins, maximum_bin_value)
    
    for intensity_category, intensity_category_entries in dataset.groupby(pd.cut(dataset['max_wind'], cut_bins)):
        
        intensity_category_name = None
        
        for category_name, category_bin in intensity_category_bins.items():
            category_tuple = (intensity_category.left, intensity_category.right - 1)
            if category_bin == category_tuple:
                intensity_category_name = category_name

        intensity_category_check_string = f"{intensity_scale}; minimum winds: {intensity_category_entries['max_wind'].min():.2f}"
        assert (intensity_category_name is not None), f'Intensity category name not found in scale: {intensity_category_check_string}'
        
        intensity_category_entries[f'category_{intensity_scale}'] = intensity_category_name
    
        container.append(intensity_category_entries)
    
    dataset = pd.concat(container).sort_values(['storm_id', 'time'])

    return dataset

def filter_storms(dataset, 
                  filter_values: dict,
                  diagnostic: bool=True):

    ''' 
    Filter storms by a certain variable. 
    Intensity filtering is not handled here.
    Input dictionary must have keys pertaining to dataset variables and values pertaining to filter quantities.
    '''

    # Iterate through all parameters and perform filtering
    for filter_variable, filter_range in filter_values.items():

        if diagnostic:
            start_time = time.time()
            print(f'Checkpoint 0: {time.time() - start_time}')
            print(f'[filter_storms] Filtering for variable {filter_variable} over range {filter_range}.')

        if isinstance(filter_range, str):
            # Adjust encoding to match IBTrACS data encoding format
            filter_range_binary = filter_range.encode('ascii')
            tmp = dataset.where(dataset[filter_variable] == filter_range_binary)
        elif isinstance(filter_range, tuple):
            tmp = dataset.where((dataset[filter_variable] >= min(filter_range)) &
                                (dataset[filter_variable] < max(filter_range)))
        elif isinstance(filter_range, list):
            tmp = dataset.where(dataset[filter_variable].isin(filter_range))

        if diagnostic:
            print(f'Checkpoint 1: {time.time() - start_time}')
            
        # Generate boolean mask based on filter matches
        mask = np.full(shape=(len(dataset['storm'])), fill_value=True)
        # Use a boolean mask to filter null entries (corresponding to False values from filter not being satisftied) per variable
        for data_var in tmp.data_vars:
            storm_indices = np.where(tmp[data_var].isnull().sum(dim='date_time') == len(tmp['date_time']), False, True)
            mask = np.multiply(storm_indices, mask)
        
        if diagnostic:
            print(f'Checkpoint 2: {time.time() - start_time}')
            
        dataset = dataset.isel(storm=mask)
        
        if diagnostic:
            print(f'Checkpoint 3: {time.time() - start_time}')
    
    return dataset

def pick_storm(track_data: pd.DataFrame,
               visualization: bool=False) -> pd.DataFrame:
    
    # Pick random storm
    storm_IDs = track_data['storm_id'].unique()
    random_storm_ID_index = np.random.randint(0, high=len(storm_IDs), size=1, dtype=int)
    random_storm_ID = storm_IDs[random_storm_ID_index].item()
    random_storm = track_data.loc[track_data['storm_id'] == random_storm_ID].sort_values('time')

    if visualization:
        fig, ax = plt.subplots(figsize=(3, 3))
        ax.plot(random_storm['center_lon'], random_storm['center_lat'], c='k', zorder=1)
        ax.scatter(random_storm['center_lon'], random_storm['center_lat'], c=random_storm['max_wind'], zorder=2)
        
        maximum_wind = random_storm['max_wind'].max()
        minimum_pressure = random_storm['min_slp'].min()
        ax.set_title(f'Storm ID: {random_storm['storm_id'].unique().item()}\nMaximum wind: {maximum_wind:.2f} m/s\nMinimum pressure: {minimum_pressure:.2f} hPa', loc='left', ha='left', fontsize=10)

    return random_storm

def main(basin_name: str='global',
         date_range: tuple[str, str]=('2010-01-01', '2011-01-01'),
         intensity_parameter: str='max_wind',
         intensity_range: tuple[int, int]=(0, np.inf)):
    
    ''' Data loading. '''
    # Load and reformat data to GFDL conventions
    dataset = load()
    dataset = format_adjustment(dataset)
    
    ''' Data filtering. '''
    # Ensure valid basin name is provided
    basin_data, _ = visualization.basins()
    assert basin_name in basin_data.keys(), f'Basin name {basin_name} not found, please try again.'
    # Ensure date range has appropriate conditions
    assert len(date_range) == 2, f'Date range has inappropriate formatting, must be 2 items long'
    
    # Modify date to approriate data type
    date_range = tuple([pd.to_datetime(date_range_entry) for date_range_entry in date_range])
    
    # Fill filtering dictionary accordingly
    filter_values = {'basin': basin_name, 
                     'time': date_range}
    # Remove basin entry if 'global' is the provided value to avoid filtering
    if basin_name == 'global':
        filter_values.pop('basin')

    # Perform non-intensity filtering and metadata generation
    filtered_dataset = filter_storms(dataset, filter_values=filter_values, diagnostic=False)
    
    # Generate GFDL QuickTracks-style DataFrame
    track_data = generate_dataframe(filtered_dataset)

    # Perform intensity filtering
    # Note: intensity filtering is easier after initial filtering and DataFrame conversion
    track_data = TC_tracker.intensity_filter(track_data, intensity_parameter, intensity_range)
    
    print(track_data['max_wind'].describe())

    # Assign intensity category
    track_data = intensity_categorization(track_data)
    
    return track_data