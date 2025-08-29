''' Import packages. '''
# Time packages
import calendar, cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import functools, importlib, os, pickle, collections, sys

import pandas as pd, xarray as xr, nc_time_axis
xr.set_options(keep_attrs=True)
# Visualization tools
import cartopy, cartopy.util as cutil, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import derived, tc_analysis, utilities, visualization, track_TCs
importlib.reload(utilities) 
importlib.reload(visualization)
importlib.reload(tc_analysis)
importlib.reload(derived)
importlib.reload(track_TCs)

def define_GCM_fields(fields: list[str]) -> dict:

    ''' Helper function to construct an input dictionary to pull post-processed GCM data. '''

    fields_reference = {'precip': {'domain': 'atmos', 'level': None},
                        'evap': {'domain': 'atmos', 'level': None},
                        'shflx': {'domain': 'atmos', 'level': None},
                        'WVP': {'domain': 'atmos', 'level': None},
                        'olr': {'domain': 'atmos', 'level': None},
                        'netrad_toa': {'domain': 'atmos', 'level': None},
                        'rh': {'domain': 'atmos', 'level': '500'},
                        'sphum': {'domain': 'atmos', 'level': None},
                        't_surf': {'domain': 'atmos', 'level': None},
                        'cld_amt': {'domain': 'atmos', 'level': '200'},
                        'rh200': {'domain': 'atmos', 'level': None},
                        'rh500': {'domain': 'atmos', 'level': None},
                        'rh850': {'domain': 'atmos', 'level': None},
                        'vort850': {'domain': 'atmos', 'level': None},
                        'omega': {'domain': 'atmos', 'level': '500'},
                        'swfq': {'domain': 'atmos', 'level': None},
                        'tau_x': {'domain': 'atmos', 'level': None},
                        'tau_y': {'domain': 'atmos', 'level': None},
                        'ucomp': {'domain': 'atmos', 'level': 1000},
                        'vcomp': {'domain': 'atmos', 'level': 1000},
                        'wind': {'domain': 'atmos', 'level': None},
                        'u_ref': {'domain': 'atmos', 'level': None},
                        'v_ref': {'domain': 'atmos', 'level': None},}

    field_dictionaries = {}
    for field_name in fields:
        if field_name in fields_reference.keys():
            field_dictionaries[field_name] = fields_reference[field_name]
        else:
            field_dictionaries[field_name] = {'domain': 'atmos', 'level': None}
    
    return field_dictionaries

def load_TC_track_data(experiment_configurations: dict[str: tuple[int, int]],
                       diagnostic: bool=False) -> dict:

    ''' 
    Loads TC track data for GCM experiments for a given model.
    
    Prerequisites:
    - TC tracker must have already been run. This is confirmed with an "analysis_lmh" directory corresponding to the experiment.
    - Experiment name must be logged in utilities.py.

    Workflow
    1. For each experiment, find the corresponding track directory
    2. Scrape track data for the experiment
    2a. If the track data exists on /projects already, pull from there
    2b. Else, manually load the track data 
    '''

    diagnostic_tag = f'[load_TC_track_data()]'

    # Ensure data types are appropriate. Convert to list if a string.
    assert isinstance(experiment_configurations, dict), "`experiment_names` must be a dictionary."
        
    # 1. Load paths to directories 
    directories = {}
    for experiment_configuration in experiment_configurations.keys():
        # Get configuration details
        model_name, experiment_name = experiment_configuration.split('-') if '-' in experiment_configuration else (experiment_configuration, '')
        # Get configuration-specific directory
        experiment_dirname = utilities.directories(model=model_name, experiment=experiment_name, data_type='track_data') 
        if experiment_dirname is not None:
            directories[experiment_name] = experiment_dirname

    # 2. Scrape track data for the experiment
    track_data = {}
    for experiment_configuration, year_range in experiment_configurations.items():
        # Get configuration details
        model_name, experiment_name = experiment_configuration.split('-') if '-' in experiment_configuration else (experiment_configuration, '')
        # Check to see if data is pre-loaded
        temporary_track_data = tc_analysis.load_TC_tracks(model=model_name,
                                                          experiment=experiment_name,
                                                          year_range=year_range,
                                                          print_statistics=False)
        # If so, append to the dictionary
        if temporary_track_data is not None:
            if diagnostic: print(f'{diagnostic_tag} Preloaded data being used...')
            track_data[experiment_configuration] = temporary_track_data[model_name][experiment_name]
        # Else, load manually
        else:
            if diagnostic:  print('Data must be generated to be used...')
            # Get year range available (must be continuous)
            available_track_years = [int(subdir.split('_')[-1]) for subdir in os.listdir(directories[experiment_name])
                                     if 'atmos_' in subdir]
            # year_range = (min(available_track_years), max(available_track_years))
            # Determine if data can be loaded in parallel based on number of year available
            parallel_load = True if max(year_range) - min(year_range) > 2 else False
            # Load data manually
            track_data[experiment_configuration] = tc_analysis.tc_track_data(models=[model_name], 
                                                                             experiments=[experiment_name], 
                                                                             year_range=year_range,
                                                                             parallel_load=parallel_load,
                                                                             diagnostic=False)[model_name][experiment_name]

        TMP = track_data[experiment_configuration]['unique']
        num_years = np.round((TMP['time'].max() -TMP['time'].min()).days / 365)
        if diagnostic:
            print(f'{diagnostic_tag} Number of storms per year, experiment {experiment_name}: {(len(TMP) / num_years):.2f}')
            print('\t#######################################################')
        
    return track_data

def load_GCM_data(model_name: str,
                  experiment_names: str|list[str],
                  field_names: list[str],
                  year_range: tuple[int, int] | None = None,
                  data_type: str='atmos_month') -> dict:

    ''' 
    Loads GCM model data for a given list of model and experiment configurations.
    
    Prerequisites:
    - Experiment name must be logged in utilities.py.
    '''

    # Ensure data types are appropriate. Convert to list if a string.
    assert isinstance(experiment_names, str) or isinstance(experiment_names, list), "`experiment_names` must be a string or list of strings."
    if isinstance(experiment_names, str):
        experiment_names = [experiment_names]

    # 1. Construct field dictionary for function input
    field_dictionary = define_GCM_fields(field_names)
    
    # 2. Load paths to directories 
    directories = {}
    for experiment_name in experiment_names:
        experiment_dirname = utilities.directories(model=model_name, experiment=experiment_name, data_type='model_output') 
        if experiment_dirname is not None:
            directories[experiment_name] = experiment_dirname

    # 3. Scrape GCM data for the experiment
    GCM_data = {model_name: {}}
    for experiment_name in experiment_names:
        # Check to see if data is pre-loaded
        # If so, append to the dictionary
        try:
            # Get year range for iterand experiment
            filenames = os.listdir(directories[experiment_name])
            filename_years = [int(f[:4]) for f in filenames if
                            f.endswith('.nc') and
                            data_type in f]
            year_range = (min(filename_years), max(filename_years)) if year_range is None else year_range
            # Note: year range temporarily suppressed because it is not relevant to preliminary analysis
            print(model_name, experiment_name, year_range)
            temporary_GCM_data = utilities.postprocessed_data_load(model_name, 
                                                                   experiment_name, 
                                                                   field_dictionary, 
                                                                   year_range=year_range,
                                                                   data_type=data_type,
                                                                   diagnostic=False)
            GCM_data[model_name][experiment_name] = temporary_GCM_data[model_name][experiment_name]
            print('Data loaded from refined postprocessed data.')
        # Else, load manually
        except:
            print('Data loaded from unrefined postprocessed data.')
            # Pull pathnames to raw postprocessed data
            pathnames = [os.path.join(directories[experiment_name], f) for f in os.listdir(directories[experiment_name]) if
                         f.endswith('.nc') and
                         data_type in f]
            # Load the data
            temporary_GCM_data = xr.open_mfdataset(pathnames).load()
            # Obtain fields common to the requested field set and the data variables available
            temporary_GCM_fields = set(field_names) & set(temporary_GCM_data.data_vars)
            temporary_GCM_data = temporary_GCM_data[temporary_GCM_fields]
            # Iterate over all fields to ensure any data with vertical levels are subselected to a single level
            processed_field_names = []
            for field_name, field_properties in field_dictionary.items():
                if field_name in temporary_GCM_data.data_vars:
                    if 'pfull' in temporary_GCM_data[field_name].dims and 'level' in field_properties.keys() and field_properties['level'] is not None:
                        processed_field_name = f'{field_name}{field_properties['level']}'
                        temporary_GCM_data[processed_field_name] = temporary_GCM_data[field_name].sel(pfull=field_properties['level'], method='nearest')
                    else:
                        processed_field_name = f'{field_name}'
                    processed_field_names.append(processed_field_name)
            GCM_data[model_name][experiment_name] = temporary_GCM_data[processed_field_names]
                
    return GCM_data

def derived_fields(data: dict):
    ''' Derive fields that are not automatically available in post-processed data. '''
    for model_name in data.keys():
        for experiment_name in data[model_name].keys():
            data[model_name][experiment_name] = derived.TC_surface_wind_speed(data[model_name][experiment_name])
            data[model_name][experiment_name] = derived.TC_surface_moisture_flux(data[model_name][experiment_name])
    return data

def unit_conversion(data: dict):
    ''' Convert units on post-processed data. '''
    for model_name in data.keys():
        for experiment_name in data[model_name].keys():
            for field_name in data[model_name][experiment_name].data_vars:
                if field_name in ['precip', 'evap']:
                    data[model_name][experiment_name][field_name] = data[model_name][experiment_name][field_name] * 86400
    return data

def get_TC_statistics(configurations: dict, 
                      diagnostic: bool=False) -> dict:

    diagnostic_tag = f'[get_TC_statistics()]'
    # Configuration name format is '{MODEL_NAME}-{EXPERIMENT_NAME}.{CONFIGURATION_TYPE}'

    # Working and output dictionary for statistics
    statistics_TMP, statistics_TC = {}, {}
    
    # Load track data for all configurations
    config_track_data = load_TC_track_data(experiment_configurations=configurations)

    # First pass: construct the custom dictionary
    for config_name in config_track_data.keys():
        if diagnostic: print(f'{diagnostic_tag} working on configuration {config_name}...')
        if 'IBTrACS' in config_name: continue # we don't need to process IBTrACS data in this function
        # Pull configuration-specific metadata
        model_name = config_name.split('-')[0] # get model name
        experiment_type, config_type = config_name.split('-')[1].split('.') # get experiment and configuration names
        statistics_TMP[f'{model_name}-{config_type}'] = {'CTL': None, 'EXP': None}
        statistics_TC[f'{model_name}-{config_type}'] = {'CTL': {}, 'EXP': {}}

    # Second pass: scrape unique TC track data for all configurations and place into a custom data structure for control-perturbation comparison
    for config_name in config_track_data.keys():
        if 'IBTrACS' in config_name: continue # we don't need to process IBTrACS data in this function
        # Pull configuration-specific metadata
        model_name = config_name.split('-')[0] # get model name
        experiment_type, config_type = config_name.split('-')[1].split('.') # get experiment and configuration names

        # Assign configuration type as a binary of either control or SWISHE
        assert experiment_type in ['CTL1990', 'CTL1990_SWISHE'], f'Experiment type is currently restricted to two entries: `CTL1990` and `CTL1990_SWISHE`. Please check this configuration, {experiment_type}.'
        experiment_type = 'CTL' if experiment_type == 'CTL1990' else 'EXP'

        # Assign dictionary 
        statistics_TMP[f'{model_name}-{config_type}'][experiment_type] = config_track_data[config_name]['unique']

    # Third pass: get statistics for TC activity for all configurations
    for config_name in statistics_TMP.keys():
        if 'IBTrACS' in config_name: continue # we don't need to process IBTrACS data in this function
        # Get model and configuration names
        model_name, config_type = config_name.split('-')
        # Iterate over each expeirment type to get statistics
        for experiment_type in statistics_TMP[config_name]:
            if diagnostic: print(f'{diagnostic_tag} Parsing statistics on the {experiment_type} experiment for configuration {config_name}...')

            TMP = statistics_TMP[config_name][experiment_type] # define a shorthand for the iterand configuration data
            TMP_HU = TMP.loc[TMP['max_wind'] >= 33] # filter for hurricane-strength TCs
            number_of_years = np.round((TMP['cftime'].max() - TMP['cftime'].min()).days / 365).astype(int) # number of years

            # Load statistics of interest into the dictionary
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['NUM_YR'] = number_of_years
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['COUNT_TS'] = TMP['storm_id'].count()
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['COUNT_PER_YR_TS'] = TMP['storm_id'].count() / number_of_years
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['MEDIAN_SLP_TS'] = TMP['min_slp'].median()
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['COUNT_HU'] = TMP_HU['storm_id'].count()
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['COUNT_PER_YR_HU'] = TMP_HU['storm_id'].count() / number_of_years
            statistics_TC[f'{model_name}-{config_type}'][experiment_type]['MEDIAN_SLP_HU'] = TMP_HU['min_slp'].median()
            
    # Fourth pass: get differences between control and perturbations
    for config_name in statistics_TC.keys():
        # Get model and configuration names
        model_name, config_type = config_name.split('-')
        # Get experiment types
        experiment_types = list(statistics_TC[f'{model_name}-{config_type}'].keys())
        assert len(experiment_types) == 2

        # Generate the difference experiment type
        statistics_TC[f'{model_name}-{config_type}']['CTL-EXP'] = {}
        # Get statistics differences between control and perturbations
        for statistic in statistics_TC[f'{model_name}-{config_type}'][experiment_types[0]]:
            get_percent_difference = lambda x, y: 100 * ((x - y) / y)
            if statistic in ['NUM_YR']: continue # skip over useless statistics to compare
            statistics_TC[f'{model_name}-{config_type}']['CTL-EXP'][statistic] = get_percent_difference(statistics_TC[f'{model_name}-{config_type}']['EXP'][statistic], statistics_TC[f'{model_name}-{config_type}']['CTL'][statistic])

    # Re-organize dictionary for re-assignment as a DataFrame
    statistics_TC_nested = {(config_name, experiment_name): TMP_EXP 
                            for config_name, TMP_CONFIG in statistics_TC.items() # outer loop
                            for experiment_name, TMP_EXP in TMP_CONFIG.items()} # inner loop
    # Generate a DataFrame, sort indices,a nd round to two decimal places
    statistics_TC = pd.DataFrame.from_dict(statistics_TC_nested, orient='index').sort_index().round(2)

    # Perform reindexing for legibility
    # Step 1: get level 0 (configuration type) order to match input experiment order
    get_config_name = lambda x: '-'.join([x.split('-')[0], x.split('-')[1].split('.')[1]])
    configuration_type_list = [get_config_name(k) for k in configurations.keys() if 'IBTrACS' not in k]
    configuration_type_dict = {k: index for index, k in enumerate(configuration_type_list)}
    configuration_type_list = list(configuration_type_dict.keys())
    # Step 2: get level 1 (experiment type) order
    experiment_type_list = ['CTL', 'EXP', 'CTL-EXP']
    # Step 3: reorder accordingly
    statistics_TC = statistics_TC.reindex(index=configuration_type_list, level=0).reindex(index=experiment_type_list, level=1)

    return statistics_TC