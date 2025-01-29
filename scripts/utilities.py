import cftime
import calendar
import datetime
import numba
import numpy as np
import scipy
import pandas as pd
import xarray as xr
import os
import pickle
import visualization
import random
import time
from multiprocess import Pool
import collections

import warnings
warnings.filterwarnings("ignore")

import track_TCs

def cftime_calendar_type(timestamp: cftime.datetime) -> str:
    
    ''' 
    Method to store calendar definitions for `cftime` to allow for calendar type assignment for timestamps. 
    Reference: https://unidata.github.io/cftime/api.html#cftime.datetime
    '''
    
    # Ensure timestamp is a cftime object
    assert 'cftime' in str(type(timestamp)), f'[utilities.cftime_calendar_reference()] Timestamp is not a cftime object.'
    # Get timestamp calendar type
    timestamp_calendar = str(type(timestamp))
    # Define calendar types, with keys being the string outputs of type(`timestamp``), where timestamp is a cftime timestamp value.
    #                        and values being the calendar type to be assigned to a new cftime object
    calendar_types = {"<class 'cftime._cftime.DatetimeJulian'>": 'julian',
                      "<class 'cftime._cftime.datetime'>": 'noleap',
                      "<class 'cftime._cftime.DatetimeNoLeap'>": 'noleap'}
    
    assert timestamp_calendar in calendar_types.keys(), f'[utilities.cftime_calendar_reference()] Calendar type {timestamp_calendar} not found in reference dictionary.'
    
    return calendar_types[timestamp_calendar]

def directories(model, experiment, data_type='model_output'):
    
    """
    Method to log the directories of all raw model runs and return a corresponding path.
    """
    
    dirnames = {'AM2.5': {'CTL1990s': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                       'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s/POSTP'},
                          'CTL1990s_swishe': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'},
                          'CTL1990s_ewishe': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_ewishe_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_ewishe_tigercpu_intelmpi_18_540PE/POSTP'},
                          'CTL1990s_ewishe_5x': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_ewishe_5x_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_ewishe_5x_tigercpu_intelmpi_18_540PE/POSTP'},
                          'CTL1990s_plus2K': {'track_data': '/projects2/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_plus2K/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/projects2/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_plus2K/POSTP'},
                          'CTL1990s_tiger3': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/AM2.5/work/CTL1990s_tiger3_tiger3_intelmpi_24_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/AM2.5/work/CTL1990s_tiger3_tiger3_intelmpi_24_540PE/POSTP'},
                          'CTL1990s_swishe_tiger3': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/AM2.5/work/CTL1990s_swishe_tiger3_tiger3_intelmpi_24_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/AM2.5/work/CTL1990s_swishe_tiger3_tiger3_intelmpi_24_540PE/POSTP'},
                          'CTL1990s_ewishe2X_tiger3': {'track_data': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe2X_tiger3/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe2X_tiger3/POSTP'},
                          'CTL1990s_ewishe4X_tiger3': {'track_data': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe4X_tiger3/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe4X_tiger3/POSTP'},
                          'CTL1990s_ewishe8X_tiger3': {'track_data': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe8X_tiger3/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                              'model_output': '/projects/GEOCLIM/gr7610/MODEL_OUT/AM2.5/CTL1990s_ewishe8X_tiger3/POSTP'},
                          'CTL1990s_swishe_plus2K': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_plus2K_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                                     'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5/work/CTL1990s_swishe_plus2K_tigercpu_intelmpi_18_540PE/POSTP'}},
                'AM2.5C360': {'CTL1990s': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                           'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_tigercpu_intelmpi_18_1080PE/POSTP'},
                              'CTL1990s_swishe': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/analysis_lmh/cyclones_gav_ro110_330k',
                                                  'model_output': '/scratch/gpfs/GEOCLIM/gr7610/AM2.5C360/work/CTL1990s_swishe_tigercpu_intelmpi_18_1080PE/POSTP'}},
                'FLOR': {'CTL1990s': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s/POSTP'},
                         'CTL1990s_FA': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_FA/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_FA/POSTP'},
                         'CTL1990s_FA_tiger3': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/FLOR/exp/CTL1990s_FA_tiger3/work/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/FLOR/exp/CTL1990s_FA_tiger3/work/POSTP'},
                         'CTL1990s-8xdaily': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s-8xdaily_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s-8xdaily_tigercpu_intelmpi_18_576PE/POSTP'},
                         'CTL1990s_swishe': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe/POSTP'},
                         'CTL1990s_swishe-ens01': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe-ens01/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                                   'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe-ens01/POSTP'},
                         'CTL1990s_swishe-ens02': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe-ens02/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                                   'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe-ens02/POSTP'},
                         'CTL1990s_swishe-8xdaily': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_576PE/POSTP'},
                         'CTL1990s_swishe_FA': {'track_data': '/projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe_FA/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_swishe_FA/POSTP'},
                         'CTL1990s_swishe_FA_tiger3': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/FLOR/exp/CTL1990s_swishe_FA_tiger3/work/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                      'model_output': '/scratch/gpfs/GEOCLIM/gr7610/tiger3/FLOR/exp/CTL1990s_swishe_FA_tiger3/work/POSTP'},
                         'CTL1990s_ewishe_5x': {'track_data': '/projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_ewishe_5x/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/projects/GEOCLIM/gr7610/MODEL_OUT/FLOR/CTL1990s_ewishe_5x/POSTP'},
                         'CTL1990s_2xCO2': {'track_data': '/projects/w/wenchang/MODEL_OUT/FLOR/CTL1990_v201905_2xCO2_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/scratch/gpfs/wenchang/FLOR/work/CTL1990_v201905_ic2001_2xCO2_tigercpu_intelmpi_18_576PE/POSTP'},
                         'CTL1990s_swishe_2xCO2': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_swishe_2xCO2_tigercpu_intelmpi_18_576PE/analysis_lmh/cyclones_gav_ro110_1C_330k',
                                     'model_output': '/scratch/gpfs/GEOCLIM/gr7610/FLOR/work/CTL1990s_swishe_2xCO2_tigercpu_intelmpi_18_576PE/POSTP'}},
                'HIRAM': {'CTL1990s': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                      'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s/POSTP'},
                          'CTL1990s_swishe': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s_swishe/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                     'model_output': '/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/tmp/CTL1990s_swishe_tigercpu_intelmpi_18_540PE/POSTP'}},
                'HIRAM-8xdaily': {'control': {'track_data': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s-8xdaily_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                              'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s-8xdaily_tigercpu_intelmpi_18_540PE/POSTP'},
                                  'swishe': {'track_data': '/scratch/gpfs/GEOCLIM/gr7610/HIRAM/work/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_540PE/analysis_lmh/cyclones_gav_ro110_2p5C_330k',
                                             'model_output': '/tigress/GEOCLIM/gr7610/MODEL_OUT/HIRAM/CTL1990s_swishe-8xdaily_tigercpu_intelmpi_18_540PE/POSTP'}}}
    
    return dirnames[model][experiment][data_type]

def month_letter(month):
    month_letters = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    return month_letters[month-1]

def access(model, experiment, storm_type, storm_id=None, processed=False, diagnostic=False):
    
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
            if model in filename and '{0}-'.format(experiment) in filename and storm_type in filename]
    # Initialize filename container
    filename = None
    # If a specific storm ID is given, check for it
    storm_exists = False # check for storm existence - True if exists
    
    if diagnostic:
        print(storm_id)
        [print(file) for file in files]
    
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
            try:
                data = pickle.load(f)
            except:
                print('File could not be accessed: {0}'.format(filename))
                pass
    else:
        with open(random.choice(files), 'rb') as f:
            data = pickle.load(f)
    
    if filename:
        print('[utilities.py, access()] Loaded filename {0} for model {1} and experiment {2}.'.format(filename, model, experiment))  
    return filename, data

def postprocess_access(models: list, 
                       experiments: list, 
                       fields: dict, 
                       year_range: tuple[int, int], 
                       FLOR_year_adjustment: int=2050,
                       diagnostic: bool=False) -> dict:

    '''
    Determine path names for files from which climate model output data will be obtained.

    Args:
    - models (list): list of model names 
    - experiments (list): list of experiment names corresponding to each model
    - fields (list): list of fields from which to obtain data
    - year_range (list or tuple, 2-item): list of years over which to obtain model data
    - FLOR_year_adjustment (int, optional): number of years to adjust data for in FLOR
    Returns:
    - paths (dict): dictionary holding path names for model data, as well as model year ranges. 
                    Dictionary structure is: paths{model: {year_range, path}}
    '''

    # Directory holding post-processed climate model data
    dirname = '/tigress/GEOCLIM/gr7610/analysis/model_out'
    # Initialize data structure to hold pathnames and model year ranges for each model and experiment configuration
    paths = {model: {} for model in models}
    # Find matching data for each model and experiment
    for model in models:
        # Adjust the year range for FLOR runs
        model_year_range = (min(year_range) + FLOR_year_adjustment, max(year_range) + FLOR_year_adjustment) if 'FLOR' in model else year_range
        paths[model]['year_range'] = model_year_range
        # Iterate by experiment
        for experiment in experiments:
            paths[model][experiment] = {}
            for field, field_information in fields.items():
                # Fetch metadata for field
                domain, level = field_information['domain'], field_information['level']
                # Redefine vertical level string
                level = '-full-' if not level else level
                # Filter files by model, experiment, field, and vertical level
                filenames = [f for f in os.listdir(dirname) if
                             model in f and 
                             '{0}-'.format(experiment) in f and 
                             field in f and
                             domain in f and
                             level in f and
                             f.endswith('nc')]
                # Filter files to ensure the year range requested is covered by the found data
                filenames = [f for f in filenames if
                             int(f.split('.nc')[0].split('-')[-1].split('_')[0]) <= min(model_year_range) and
                             int(f.split('.nc')[0].split('-')[-1].split('_')[1]) >= max(model_year_range)]
                # Report how many files found for the requested field
                if diagnostic:
                    print('{0} files found for the {1} model in the {2} experiment configuration for field {3}...'.format(len(filenames), model, experiment, field))
                paths[model][experiment][field] = os.path.join(dirname, filenames[0]) if len(filenames) > 0 else None

    return paths

def postprocessed_data_load(models: str | list[str], 
                           experiments: str | list[str],
                           fields: dict, 
                           year_range: tuple[int, int], 
                           month_range: tuple[int, int] = (1, 12),
                           difference_experiment: tuple[str, str] | None = False,
                           load_full_time: bool=False):

    '''
    Loads data for a given list of models, experiments, and fields.

    Args:
    - models (list): string or list of model names 
    - experiments (list): string or list of experiment names corresponding to each model
    - fields (dict): dictionary consisting of field name as the top-level key, with the climate domain (atmos or ocean) and vertical level as subkeys. See below for an example.
    - year_range (tuple of int, 2 items): minimum and maximum years for data processing
    - difference_experiment (tuple of str, 2 items; default: False)
    - load_full_time (bool): boolean to dictate whether the full postprocessed dataset will be loaded and ignore the given year range
    Returns:
    - data (dict):  dictionary with keys corresponding to model names, subkeys corresponding to experiment names, 
                    and subkey values corresponding to xArray Datasets containing the requested fields

    Example: `fields` dictionary for ocean temperature, precipitation, and TOA net radiation
    - fields = {'temp': {'domain': 'ocean',
                         'level': '5m'},
                'precip': {'domain': 'atmos',
                           'level': None},
                'netrad_toa': {'domain': 'atmos',
                               'level': None}}

    '''

    # Ensure year range has 2 items
    assert len(year_range) == 2
    # Intake inputs based on data type
    models = [models] if isinstance(models, str) else models if isinstance(models, list) else print('Model input type not recognized. Must be a string or iterable of strings.')
    experiments = [experiments] if isinstance(experiments, str) else experiments if isinstance(experiments, list) else print('Model input type not recognized. Must be a string or iterable of strings.')
    # Ensure difference_experiment has 2 items
    if difference_experiment:
        assert len(difference_experiment) == 2
    
    # Initialize container dictionary
    data = {model: {} for model in models}
    # Obtain path names and model years for the requested models, experiments, and fields
    paths = postprocess_access(models, experiments, fields, year_range, FLOR_year_adjustment=0)
    # Iterate over each model to load data
    for model in models:
        start_year, end_year = ['{0:04d}-01-01'.format(min(paths[model]['year_range'])), 
                                '{0:04d}-01-01'.format(max(paths[model]['year_range'])+1)]
        for experiment in experiments:
            data[model][experiment] = {}
            # Load data for the given configuration, and slice by time to minimize data loading into memorry
            for field in fields.keys():
                data[model][experiment][field] = xr.open_dataset(paths[model][experiment][field])[field] if load_full_time else xr.open_dataset(paths[model][experiment][field])[field].sel(time=slice(start_year, end_year))
                data[model][experiment][field] = ocean_grid_alignment(data[model][experiment][field])
            # Merge dictionary values for each model and experiment pair into an xArray Dataset for a concise data structure
            data[model][experiment] = xr.merge(data[model][experiment].values())
            # Filter by month
            month_range_filter = month_selector(data[model][experiment]['time.month'], min(month_range), max(month_range))
            data[model][experiment] = data[model][experiment].sel(time=month_range_filter)
            # Correct field data for units
            data[model][experiment] = field_correction(data[model][experiment])

    # Generate the difference experiment dataset based on inputs, iterating over each model
    if difference_experiment:
        for model in models:
            # Assign variable names to difference experiment names
            experiment_CTL, experiment_EXP = difference_experiment
            # Initialize a counter to determine if experiment needed to compute the difference are present
            experiment_name_counter = 0
            for experiment_name in experiments:
                # Increment the field counter for each match found
                experiment_name_counter += 1 if experiment_name in experiments else 0
            if experiment_name_counter == len(experiments):
                experiment_DIFF = f'{experiment_CTL}-{experiment_EXP}'
                data[model][experiment_DIFF] = data[model][experiment_CTL] - data[model][experiment_EXP]
            else:
                print('Not all experiments loaded to compute difference between experiments {0} and {1}. The only experiments loaded are: {2}'.format(experiment_CTL, experiment_EXP, experiments))
    
    return data

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
        try:
            year_adjust = 1900 if ((timestamp.year >= 1900) & (model != 'FLOR')) else 0
            out = cftime.DatetimeNoLeap(year=timestamp.year-year_adjust, month=timestamp.month, day=timestamp.day, hour=timestamp.hour)
        except:
            year_adjust = 1900 if ((timestamp.dt.year >= 1900) & (model != 'FLOR')) else 0
            out = cftime.DatetimeNoLeap(year=timestamp.dt.year-year_adjust, month=timestamp.dt.month, day=timestamp.dt.day, hour=timestamp.dt.hour)
            
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
    else:
        output = data.weighted(np.cos(np.deg2rad(data.grid_yt)))
        return output

def latitude_weighting(data):

    ''' Check to ensure data has the right fields. '''
    coordinates = ['grid_yt', 'grid_xt']
    if contains(['grid_yt'], data.coords):
        return data.weighted(np.cos(np.deg2rad(data.grid_yt)))
    else:
        print('[latitude_weighting()] Input data does not contain requisite data variables and/or coordinates.')
        return None
        
def month_selector(month, start, end=None):
    ''' Helper function that takes a DataArray time.month index and integer arguments for start and end months for xArray time indexing, end-inclusive. '''
    end = start if not end else end
    return (month >= start) & (month <= end)

def domain_differentiation(data, dim, diff_method='2nd-order'):

    ''' 
    Differentiate an xarray DataArray field along a given dimension to preserve shape. 
    Last row or column in the bulk differentiation array will be appended, effectively repeating that last row or column.
    Data must be at least 3-dimensional (time, grid_xt, grid_yt OR time, grid_xt, grid_yt, pfull OR grid_xt, grid_yt, pfull).
    '''

    # Trim unused dimensions
    for drop_dim in ['bnds', 'phalf']:
        data = data.drop_dims(drop_dim) if drop_dim in data.dims else data
    
    # Ensure proper dimensional order
    if 'pfull' in data.dims and 'time' in data.dims:
        data = data.transpose('time', 'grid_xt', 'grid_yt', 'pfull')
    elif 'pfull' in data.dims and 'time' not in data.dims:
        data = data.transpose('grid_xt', 'grid_yt', 'pfull')
    else:
        data = data.transpose('time', 'grid_xt', 'grid_yt')
    
    # Get distance between grid points (essentially, distance_lon = dx, distance_lat = dy)
    distance_lon, distance_lat = distance_grid(data)
    distance = distance_lon if dim == 'grid_xt' else distance_lat
    
    # Get the bulk differentiation (will result in an output with the shape of data, minus one entry in the differentiation dimension)
    if diff_method == 'first-order':
        a = data.diff(dim=dim)/distance
    else:
        a = data.differentiate(dim)/(2*distance)
    
    # Get the last row/columns of the differentiation array and repeat to append to bulk array and preserve original dimensions
    b = a[{dim: -1}]

    # print(a, b)
    
    # Concatenate along respective axes
    if dim == 'grid_xt':
        if 'time' in data.dims:
            b_ = b.values[:, np.newaxis, :, :] if 'pfull' in data.dims else b.values[:, np.newaxis, :]
            c = xr.DataArray(data=np.concatenate((a.values, b_), axis=1), dims=a.dims).dropna(dim=dim).isel(grid_xt=range(0, len(data.grid_xt)))
        else: 
            b_ = b.values[np.newaxis, :, :]
            c = xr.DataArray(data=np.concatenate((a.values, b_), axis=0), dims=a.dims).dropna(dim=dim).isel(grid_xt=range(0, len(data.grid_xt)))
    elif dim == 'grid_yt':
        if 'time' in data.dims:
            b_ = b.values[:, :, np.newaxis, :] if 'pfull' in data.dims else b.values[:, :, np.newaxis]
            c = xr.DataArray(data=np.concatenate((a.values, b_), axis=2), dims=a.dims).dropna(dim=dim).isel(grid_yt=range(0, len(data.grid_yt)))
        else:
            b_ = b.values[:, :, np.newaxis]
            c = xr.DataArray(data=np.concatenate((a.values, b_), axis=2), dims=a.dims).dropna(dim=dim).isel(grid_yt=range(0, len(data.grid_yt)))

    # print('----------------------------------------------------------------------------------------')
    # print(c)
    # print('\n')
    
    return c

def intensity_binning(mode='track_output', data=None, intensity_metric='max_wind'):
    """
    Method to generate bins for TC intensities to support compositing based on a given intensity metric.
    Note: this is chosen to be performed post-storage to keep raw TC data as clean as possible from metadata additions, in case of future definition changes.
    Note: the intensity metric will either be by maximum winds (max_wind) or minimum sea-level pressure (min_slp).
    References: levels picked as in Harris et al (2016), 10.1175/JCLI-D-15-0389.1
    
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
        intensity_bin_limits = [0, 17.5, 30, np.inf] 
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

def retrieve_tracked_TCs_OLD(model, experiment, storm_type, year_range, intensity_metric=None, intensity_threshold=None, config=None, diagnostic=False):
    
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

    if not intensity_metric:
        intensity_metric = 'min_slp'
    if not intensity_threshold:
        intensity_threshold = 1000
    
    # Retrieve directory containing parent directory for track data
    dirname = directories(model, experiment, data_type='track_data')
    if diagnostic:
        print(f"[utilities.py, retrieve_tracked_TCs()] retrieving data from {dirname} for years {year_range}")
    
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
        if diagnostic:
            print('Years evaluated from LMH output: {0} to {1}'.format(min(year_range) + year_adjust, max(year_range) + year_adjust))
        fnames = [f for f in fnames 
                  if min(year_range) + year_adjust <= datetime.datetime(year=int(f.split('/')[-3].split('_')[-1]), month=1, day=1).year <= max(year_range) + year_adjust]
    
    if diagnostic:
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
        # storm['time'] = storm['time'] + pd.offsets.DateOffset(years=start_year_adjust)
        ''' Velocity (speed, direction) derivation. '''
        # Initialize dictionary for preliminary storage. Will be reassigned into the DataFrame by the join() method using time as the matching criterion.
        velocity = {'time': [storm.iloc[0]['time']], 'speed': [np.nan], 'direction': [np.nan], 'storm_motion_x': [np.nan], 'storm_motion_y': [np.nan]}
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
            storm_motion_x = speed*np.sin(np.deg2rad(360-direction))
            storm_motion_y = speed*np.cos(np.deg2rad(360-direction))
            # Append quantities to the 'velocity' dictionary
            velocity['time'].append(storm.iloc[i]['time'])    
            velocity['speed'].append(speed)    
            velocity['direction'].append(direction)
            velocity['storm_motion_x'].append(storm_motion_x)
            velocity['storm_motion_y'].append(storm_motion_y)
        # Build DataFrame
        velocity = pd.DataFrame(velocity)
        # Re-cast time column as a datetime object
        velocity['time'] = pd.to_datetime(velocity['time'])
        # Merge the storm and velocity DataFrames
        storm = storm.merge(velocity, how='left', on='time', suffixes=['_x', None]).drop(columns={'speed_x', 'direction_x'}).reset_index(drop=True)
        # Allow functionality for individual TC storage (see individual_tc_storage.py)
        if config == 'individual_tc_storage' and intensity_metric and intensity_threshold:   # Rename columns for future addition into xArray Dataset, and reset index
            storm = storm.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp', 'max_wnd': 'max_wind'}).reset_index(drop=True)
            # Filter by intensity bin. If the filtered storm entry is longer than 0, then the storm can be used.
            if intensity_metric and intensity_threshold:
                filtered_storm = []
                # If the intensity metric is minimum sea-level pressure, find all entries that are stronger than the threshold
                if intensity_metric == 'min_slp':
                    if diagnostic:
                        print('[utilities.py, retrieve_tracked_TCs()]: storm ID = {0}; minimum SLP = {1}'.format(storm_id, storm[intensity_metric].min()))
                    filtered_storm = storm.loc[storm[intensity_metric] <= intensity_threshold]
                else:
                    if diagnostic:
                        print('[utilities.py, retrieve_tracked_TCs()]: storm ID = {0}; maximum wind = {1}'.format(storm_id, storm[intensity_metric].max()))
                    filtered_storm = storm.loc[storm[intensity_metric] >= intensity_threshold]
                if len(filtered_storm) > 0:
                    # Append to the list for future concatenation
                    storms[storm_id] = storm
                    if diagnostic:
                        print('---> [utilities.py, retrieve_tracked_TCs()]: appending storm ID = {0}; {1}'.format(storm_id, filtered_storm['max_wind']))
                else: 
                    continue
            else:
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
        
        if diagnostic:
            print(f'[utilities.py, retrieve_tracked_TCs()] length of dataset: {len(data)}')
        
        return data

def retrieve_tracked_TCs(model, experiment, storm_type, year_range, diagnostic=False):
    
    data = track_TCs.main(model_name=model, 
                          experiment_name=experiment, 
                          year_range=year_range)
    
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

def land_mask(data, mask_type='land'):
    
    mask = xr.open_dataset('/projects2/GEOCLIM/gr7610/tools/land_mask.nc')['land_mask']
    mask_value = 1 if mask_type == 'ocean' else 0
    data = np.where(mask == mask_value, data, np.nan)

    return data

def load_parallel(filenames, field, num_cores=1, resampling_dims=None):
    # Initialize container list. Data will be loaded here for concatenation once all data is loaded.
    container = []
    # Define native loading function for xArray datasets
    def xr_read(filename):
        # Load the dataset for the prescribed field
        temp = xr.open_dataset(filename)[field]
        # If we're resampling, then iterate over the dictionary and resample
        if resampling_dims:
            for resample_dim, resample_size in resampling_dims.items():
                if resample_dim in temp.dims:
                    if resample_dim != 'time': 
                        temp = temp.coarsen({resample_dim: resample_size}).mean(resample_dim)
                    else:
                        temp = temp.resample({resample_dim: resample_size}).mean(resample_dim)
        return temp.load()
     
    start_time = time.time()
    
    # Create num_cores threads, where num_cores is passed into the method
    pool = Pool(num_cores)
    # Populate the container list with the filenames passed in
    container = pool.map(xr_read, filenames)
    pool.close()
    pool.join()
    
    # Concatenate the data over time
    output = xr.concat(container, dim='time').sortby('time')
    
    print('Total elapsed time: {0:.2f}s, elapsed time per file: {1:.2f}s'.format((time.time() - start_time), (time.time() - start_time)/len(filenames)))
    
    # Offload memory
    del pool, container

    return output

def pull_resampled_data(filename, field, month=None, troubleshooting=False):

    # Get model and experiment name
    model_name, experiment_name = filename.split('/')[5], filename.split('/')[6]
    # Get iterand year
    model_year = int(filename.split('/')[-1].split('.')[0][0:4])
    # If not found in the data, try finding in the resampled data directory
    resample_dirname = '/projects2/GEOCLIM/gr7610/analysis/model_out'

    if troubleshooting:
        print(filename, model_name, experiment_name, model_year)
    
    # Look for previously-resampled files and see if matching data exists
    resample_filenames = [os.path.join(resample_dirname, f) for f in os.listdir(resample_dirname) 
                          if f.endswith('nc') and (model_name in f) and ('-exp_{0}-'.format(experiment_name) in f) and 
                          (field in f) and ('atmos' in f) and ('mean_month' in f)]
    
    if troubleshooting:  
        print('Resample filenames pre-preselect: {0}'.format(resample_filenames))
    # Check for vertical dimension if the iterand field is supposed to be vertical
    if field in ['temp', 'sphum', 'slp']:
        resample_filenames = [f for f in resample_filenames if ('full' in f)]
        if troubleshooting:  
            print('Resample filenames preselect: {0}'.format(resample_filenames))
    # If files are found, then find files with matching year data.
    resample_filenames = [f for f in resample_filenames 
                          if (int(f.split('.')[-2].split('-')[-1].split('_')[0]) <= model_year) 
                          and (int(f.split('.')[-2].split('-')[-1].split('_')[1]) >= model_year)]
    
    if troubleshooting:  
        print('Resample filenames selected: {0}'.format(resample_filenames))
        
    # If files exist that satisfy the conditions, use that file. Else, end script.
    if len(resample_filenames) > 0:
        resample_filename = resample_filenames[0]
        if troubleshooting:  
            print('Using file {0} to compensate for fied {1}...'.format(resample_filename, field))

        if month:
            end_day = calendar.monthrange(1990, month+1)[1]
            output = xr.open_dataset(resample_filename)[field].sel(time=slice('{0:04d}-{1:02d}-01'.format(model_year, month+1), '{0:04d}-{1:02d}-{2:02d}'.format(model_year, month+1, end_day)))
        else:
            output = xr.open_dataset(resample_filename)[field]
            
        return output
    else:
        print('Data not available for {0} {1} {2}, consider pulling it using gcm_scraper_jupyter.ipynb...'.format(model_name, experiment_name, field))

def field_correction(data):
    '''
    Function to correct sign and apply multiplicative factor for specific fields.
    '''

    def assignment(d):
        output = d.copy()
        for field in d.data_vars:
            factor = 86400 if field in ['precip', 'evap', 'p-e'] else 1 # kg m^-2 s^-1 to mm d^-1
            output[field] = d[field] * factor
            _, output[field].attrs['units'] = visualization.field_properties(field)
        return output

    if isinstance(data, dict):
        for k in data.keys():
            if isinstance(data[k], dict):
                for sk in data[k].keys():
                    data[k][sk] = assignment(data[k][sk])
            elif isinstance(data[k], xr.Dataset):
                data[k] = assignment(data[k])
    
    else:
        data = assignment(data)

    return data

def contains(items: list, references: list):
    ''' 
    Check to ensure that a given list is a sublist of another list. 
    Useful for ensuring that the requisite data is in a given dataset.
    '''
    return set([item for item in items if item in references]) == set(items)

def in_basin(basin_masks: dict,
             basin_name: str,
             longitude: int|float|list, 
             latitude: int|float|list):
    
    ''' Method to determine if coordinate is in a pre-defined basin. '''

    # If input coordinates are not iterables, make them so. Else, make sure they're numeric arrays.
    longitude = np.array(longitude) if isinstance(longitude, collections.abc.Iterable) else np.array([longitude])
    latitude = np.array(latitude) if isinstance(latitude, collections.abc.Iterable) else np.array([latitude])
    # longitude = [longitude] if not iter(longitude) else np.array(longitude)
    # latitude = [latitude] if not iter(latitude) else np.array(latitude)
    # Make sure the input arguments are equal length.
    assert len(longitude) == len(latitude)
    input_length = len(longitude)

    # Load basin data if not provided. This may take a while to compute, not recommended for iterative solutions.
    if not basin_masks:
        _, basin_masks = visualization.basins()

    # Generate a boolean mask based on the coordinates provided by the basins dictionary input
    basin_mask = np.where(basin_masks[basin_name] > 0, True, False).ravel()
    # Obtain longitude and latitude value arrays
    x, y = basin_masks[basin_name].grid_xt.values, basin_masks[basin_name].grid_yt.values
    # Define a meshgrid, build a coordinate pair array, and filter based on the raveled mask
    X, Y = np.meshgrid(x, y)
    coordinates = np.c_[X.ravel(), Y.ravel()][basin_mask]
    # Derive the minimum grid spacing in each direction
    minimum_dx, minimum_dy = np.min(np.diff(x)), np.min(np.diff(y))
    # Calculate the minimum possible distance a point can be from a neighbor to be classified as "in the basin"
    minimum_grid_spacing = np.linalg.norm([minimum_dx, minimum_dy])
    # Use a nearest-neighbor lookup to return the distance of the point to the nearest-neighbor in `coordinates`
    tree = scipy.spatial.KDTree(coordinates)

    condition = []
    for index in range(input_length):
        distance, index = tree.query([longitude[index], latitude[index]], k=1)
        condition.append(True) if distance <= minimum_grid_spacing else condition.append(False)

    # If the distance is smaller than the minimum grid spacing, it is considered in the basin. Else, it's not.
    return condition

def ocean_grid_alignment(dataset):
    '''
    Aligns the FLOR ocean model coordinate system with the atmosphere model coordinate system.

    Args:
    - dataset (xArray DataArray or Dataset): data structure containing ocean data
    Returns:
    - dataset (xArray DataArray or Dataset): data structure containing ocean data with modified coordinate system
    '''

    ocean_coord_names = {'xt_ocean': 'grid_xt', 'xu_ocean': 'grid_xt', 'yt_ocean': 'grid_yt', 'yu_ocean': 'grid_yt'}

    # Gets modulo of longitudes and corrects it
    for coord_name in ['xt_ocean', 'xu_ocean']:
        if coord_name in dataset.coords:
            dataset[coord_name] = np.where(dataset[coord_name] < 0, dataset[coord_name] + 360, dataset[coord_name])
            # Sorts longitudes to properly arrange values in ascending order
            dataset = dataset.sortby(coord_name)
    # Modify coordinate names to unify with atmospheric grid
    if 'st_ocean' in dataset.dims:
        dataset = dataset.rename({'xt_ocean': 'grid_xt', 'yt_ocean': 'grid_yt', 'st_ocean': 'pfull'})
    else:
        for dim in dataset.dims:
            if dim in ocean_coord_names.keys():
                dataset = dataset.rename({dim: ocean_coord_names[dim]})

    return dataset

############################################################
# Begin numba-specific parallelization methods.

@numba.njit
def apply_along_axis_0(func1d, arr):
    """Like calling func1d(arr, axis=0)"""
    if arr.size == 0:
        raise RuntimeError("Must have arr.size > 0")
    ndim = arr.ndim
    if ndim == 0:
        raise RuntimeError("Must have ndim > 0")
    elif 1 == ndim:
        return func1d(arr)
    else:
        result_shape = arr.shape[1:]
        out = np.empty(result_shape, arr.dtype)
        _apply_along_axis_0(func1d, arr, out)
        return out

@numba.njit
def _apply_along_axis_0(func1d, arr, out):
    """Like calling func1d(arr, axis=0, out=out). Require arr to be 2d or bigger."""
    ndim = arr.ndim
    if ndim < 2:
        raise RuntimeError("_apply_along_axis_0 requires 2d array or bigger")
    elif ndim == 2:  # 2-dimensional case
        for i in range(len(out)):
            out[i] = func1d(arr[:, i])
    else:  # higher dimensional case
        for i, out_slice in enumerate(out):
            _apply_along_axis_0(func1d, arr[:, i], out_slice)

@numba.njit
def nb_mean_axis_0(arr):
    return apply_along_axis_0(np.nanmean, arr)

def bootstrap_3d(a, b, plane=None, N=1000, level=0.95):
    '''
    Method to run bootstrap statistical testing on 3-dimensional model output. Dimensions include time, grid_xt (longitude), grid_yt (latitude).
    Current method diagnoses difference in means and establishes statistical significance for a given level from 0 to 1.
    Uses a 2-tailed approach.
    Parallel approach produces a speedup of ~4x.
    See Delsole and Tippett (2013), Chapter 3.6 for details.
    '''

    x_axis = 2
    y_axis = 1
    time_axis = 0
            
    # Extract numeric values from the xArray Datasets
    x, y = a, b
    # Initialize empty arrays to contain output data
    out_x, out_y = [np.full((N, x.shape[y_axis], x.shape[x_axis]), np.nan), 
                    np.full((N, y.shape[y_axis], y.shape[x_axis]), np.nan)]

    @numba.njit(parallel=True)
    def bootstrap_resample(in_x, in_y, out_x, out_y, N):
        time_axis_length = min(in_x.shape[time_axis], in_y.shape[time_axis])
        for repetition in numba.prange(0, N):
            out_x_temp, out_y_temp = np.full(in_x.shape, np.nan, dtype=np.float32), np.full(in_y.shape, np.nan, dtype=np.float32)
            for k in range(0, time_axis_length):
                out_x_temp[k] = in_x[np.random.randint(0, in_x.shape[time_axis]), :, :]
                out_y_temp[k] = in_y[np.random.randint(0, in_y.shape[time_axis]), :, :]
            out_x[repetition, :, :] = nb_mean_axis_0(out_x_temp)
            out_y[repetition, :, :] = nb_mean_axis_0(out_y_temp)
        return out_x, out_y

    out_x, out_y = bootstrap_resample(x, y, out_x, out_y, N)
    
    # Get difference between datasets
    delta = out_x - out_y # Control - SWISHE
    # Get values at each respective tail 
    ci_min, ci_max = np.quantile(delta, (1 - level)/2, axis=0), np.quantile(delta, (1 + level)/2, axis=0)
    # Wherever the signs are equal, output the mean. 
    # This indicates that the confidence interval does not intersect 0, such that the null hypothesis is rejected.
    out_binary = np.where(np.sign(ci_min) == np.sign(ci_max), 1, np.nan).T
    out_full = np.where(np.sign(ci_min) == np.sign(ci_max), np.nanmedian(delta, axis=0), np.nan).T
    median = np.nanmedian(delta, axis=0).T

    return out_binary, out_full, median, delta, ci_min, ci_max
    
# End numba-specific parallelization methods.
############################################################