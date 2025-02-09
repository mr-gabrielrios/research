import datetime
import numpy as np
import pandas as pd
import scipy as sp
import os
import pickle
import random
import time
import cftime

from concurrent.futures import ProcessPoolExecutor as Pool

import utilities


def tc_model_data(models, experiments, storm_type, num_storms=-1, single_storm_id=None, diagnostic=False):
    """
    Method to load data into a 3-tiered dictionary:
    --> (1) model name, (2) experiment type, (3) TC data with [a] track output, [b] model output (planar), [c] model output (vertical)

    Args:
        models (list): list of strings with model names
        experiments (list): list of strings with model names
        storm_type (str): type of storm (TS = tropical storm or C15w = hurricane-strength storm)
        num_storms (int): number of storms to process. -1 is default and processes all storms found.
    Returns:
        data (dict): 3-tiered dictionary.
    """

    # Define directory where processed, single-storm data is kept
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    # Initialize a storage dictionary with experiment names as the top-level keys
    data = {model: {} for model in models}
    # Iterate over requested models
    for model in models:
        # Initialize a storage dictionary for the model experiments
        data[model] = {}
        # Iterate over requested experiments
        for experiment in experiments:
            # Pull storm IDs from matching file names, with storm ID of the format {YEAR}-{NUMBER}
            # Assume filename is of form: TC-{MODEL}-{EXPERIMENT}-{CATEGORY}-{YEAR}-{NUMBER}-max_wind-{MAX_WIND}-min_slp-{MIN_SLP}.pkl
            # If the storm_type is TS, limit storm LMI winds to 30 m s^-1
            if storm_type == 'TS':
                storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(dirname)
                             if (model in f) and (experiment in f) and (storm_type in f) and (int(f.split('-')[7]) < 35)]
            else:
                storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(dirname)
                             if (model in f) and (experiment in f) and (storm_type in f)]
            # and randomly select 'num_storms' number of storms
            if single_storm_id:
                storm_id_subset = storm_ids
            else:
                storm_id_subset = random.sample(storm_ids, num_storms) if num_storms <= len(
                    storm_ids) and num_storms != -1 else storm_ids
            if diagnostic:
                print('\t Storms to be processed: {0}'.format(storm_id_subset))
            # Initialize a storage list for the storms
            data[model][experiment] = {'data': []}
            # Iterate over all storms found
            for storm_id in storm_id_subset:
                # Access the processed data
                if single_storm_id and single_storm_id == storm_id:
                    if diagnostic:
                        print('Only picking up: {0}'.format(single_storm_id))
                    storm_filename, storm = utilities.access(
                        model, experiment, storm_type=storm_type, storm_id=storm_id, processed=False)
                elif single_storm_id and single_storm_id != storm_id:
                    continue
                else:
                    storm_filename, storm = utilities.access(
                        model, experiment, storm_type=storm_type, storm_id=storm_id, processed=False)
                if diagnostic:
                    print('\t \t Processing {0} from {1}'.format(
                        storm_id, storm_filename))
                # Append to storage list
                data[model][experiment]['data'].append(storm)

    return data

def tc_track_data(models: list[str], 
                  experiments: list[str], 
                  storm_type: str='TS', 
                  snapshot_type: str='lmi', 
                  year_range=None, 
                  parallel_load: bool=True,
                  diagnostic: bool=False, 
                  num_procs: int=16):
    """
    Method to extract track data (output from Lucas Harris' TC tracker).
    Sample uses include aggregate statistics, track plotting, etc.

    Args:
        models (list): list of models to get data from
        experiments (list): list of experiment types to evaluate
        storm_type (str, optional): type of storm to evaluate, either 'TS' or 'C15w'. Defaults to 'C15w'.
        snapshot_type (str, optional): type at which to get track data for storm. Typically 'lmi' or 'genesis'. Defaults to 'lmi'.
        year_range (tuple or list, optional): range of years (min, max) to evaluate over. Defaults to None.

    Returns:
        data: dictionary with the following hierarchy: 
              (1) model --> (2) experiment --> (3) [a] raw data (all data for all TCs), [b] unique TC data (single snapshot for each TC)
    """

    start_time = time.time()
    checkpoint_counter = 0
    print('[tc_track_data] Checkpoint {0}: {1:.3f} s elapsed'.format(
        checkpoint_counter, time.time() - start_time))

    # Initialize dictionary for gridded data
    data = {}
    # Iterate over each model
    for model in models:
        # Adjust years for FLOR, if the iterand model
        model_year_range = year_range
        # Initialize dictionary for the model
        data[model] = {}
        # Iterate over each experiment
        for experiment in experiments:
            print(f'[tc_track_data] working on tracks for configuration: {
                  model}, {experiment}')
            print(f'[tc_track_data] Checkpoint 1: {
                  (time.time() - start_time):.3f} s elapsed')
            if diagnostic:
                print('[tc_track_data] Processing {0} model data from the {1} experiment for {2} over years {3}'.format(
                    model, experiment, storm_type, year_range))
            # Initialize dictionary for the model + experiment configuration
            # Note: two keys are added - raw (for raw data, includes all data for each TC) and unique (includes single data point for each TC)
            if parallel_load:

                # Create num_procs threads, where num_procs controls the outer loop processor count
                # get number of years
                num_years = (max(year_range) - min(year_range))
                # prevent there from being more processors than year chunks
                num_procs = num_years // 2 if num_years // 2 <= num_procs else num_procs

                # Define year chunks
                year_chunks = np.linspace(min(model_year_range), max(model_year_range), num_procs + 1, dtype=int)

                # Create non-coterminous chunks to avoid duplicates
                processed_chunks = []
                for chunk_index, chunk in enumerate(year_chunks[:-1]):
                    if chunk_index == 0:
                        processed_chunks.append(
                            (year_chunks[chunk_index], year_chunks[chunk_index+1]))
                    else:
                        processed_chunks.append(
                            (year_chunks[chunk_index]+1, year_chunks[chunk_index+1]))

                if diagnostic:
                    print('[tc_track_data] Number of processors used: {0} for year chunks {1}'.format(num_procs, processed_chunks))

                # Create the outer loop pool
                pool = Pool(num_procs)
                # Populate the container list with the filenames passed in
                container = pool.map(utilities.retrieve_tracked_TCs, [model]*num_procs,
                                     [experiment]*num_procs, [storm_type]*num_procs, processed_chunks)
                # Concatenate all DataFrames
                # print(f'[tc_track_data] length of container list: {len(list(container))}')
                raw_tracks = pd.concat(
                    container, ignore_index=True).drop_duplicates()

                data[model][experiment] = {'raw': raw_tracks, 'unique': None}
            else:
                data[model][experiment] = {'raw': utilities.retrieve_tracked_TCs(model, experiment, storm_type, model_year_range, diagnostic=diagnostic),
                                           'unique': None}

            print('[tc_track_data] Checkpoint 2: {0:.3f} s elapsed'.format(
                time.time() - start_time))
            # Rename maximum wind columns
            data[model][experiment]['raw'] = data[model][experiment]['raw'].rename(
                columns={'max_wnd': 'max_wind'})
            # Get intensity bins for all storms
            data[model][experiment]['raw'] = utilities.intensity_binning(
                mode='track_output', data=data[model][experiment]['raw'])
            # Get unique TC storm IDs from the raw data
            storm_ids = data[model][experiment]['raw']['storm_id'].unique()
            # Initialize list to hold unique TC data
            storms = []
            # Iterate over individual TCs to extract data for the unique key
            for storm_id in storm_ids:
                if diagnostic:
                    print('\t Working on {0}...'.format(storm_id))
                # Select data for the iterand storm ID
                storm = data[model][experiment]['raw'].loc[data[model][experiment]['raw']['storm_id'] == storm_id]
                # Get snapshot of TC based on argument
                storm = utilities.storm_snapshot(storm, mode=snapshot_type)
                ###
                # Storm specific filtering can be added here!
                ###
                # Append to holding list
                storms.append(storm)

            print('[tc_track_data] Checkpoint 3: {0:.3f} s elapsed'.format(
                time.time() - start_time))
            # Concatenate storms into single DataFrame
            data[model][experiment]['unique'] = pd.concat(
                storms).sort_values('storm_id')
            # Read-across storm filtering
            # Filter out storms lasting more than 30 days
            data[model][experiment]['unique'] = data[model][experiment]['unique'].loc[data[model]
                                                                                      [experiment]['unique']['duration'] <= 30]
            # Apply additional storm_data
            # data[model][experiment]['raw'] = TC_statistics(data[model][experiment]['raw'])

    return data

def counts(mode='track_output', data=None):
    """
    Method to gather number of storms in each experiment.

    Args:
        mode (str):  analysis mode for the counting. Can either be (1) 'track_output' or (2) 'model_output'.
                     (1) 'track_output' refers to output from tc_analysis.tc_track_data(). 
                         This is meant to catalog all storms detected in the model runs, but not necessarily all analyzed for planar/azimuthal fields.
                     (2) 'model_output' refers to the 'track_output' from tc_analysis.tc_model_data().
                         This is meant to catalog all storms used for analysis for planar/azimuthal fields.
        data (dict): dictionary output to match data accepted by 'track_output' or 'model_output'. See above for descriptions.

    Returns:
        counts (Pandas DataFrame): table with data of interest
    """

    # Define intensity bins. This is subject to change every time the binning algorithm is updated.
    # Pre-definition is used to catalog instances where 0 records are found.
    intensity_bins = [k for k in utilities.intensity_binning(
        mode='bin_data', data=None, intensity_metric='max_wind').keys()]

    # Initialize dictionaries: one for TC counts, one for snapshot counts by intensity bin
    storm_counts, bin_counts = {}, {}
    # Iterate over models
    for model in data.keys():
        storm_counts[model], bin_counts[model] = {}, {}
        # Iterate over experiments
        for experiment in data[model].keys():
            # Build an aggregate DataFrame based on data input mode
            if mode == 'model_output':
                # Iterate through all track output data to obtain an aggregate DataFrame
                aggregate = pd.concat([data[model][experiment]['data'][i]['track_output']
                                      for i in range(0, len(data[model][experiment]))])
                print(aggregate['storm_id'].unique())
                # Log how many individual TCs exist in this configuration
                storm_count = len(aggregate['storm_id'].unique())
            elif mode == 'track_output':
                # Iterate through all track output data to obtain an aggregate DataFrame
                aggregate = data[model][experiment]['raw']
                # Log how many individual TCs exist in this configuration
                storm_count = len(aggregate['storm_id'].unique())
            # Initialize dictionary exclusive to the model + experiment configuration
            storm_counts[model][experiment], bin_counts[model][experiment] = storm_count, {
            }
            # Get number of timestamps per intensity bin and append to dictionary
            for intensity_bin in intensity_bins:
                # Get instances in this intensity bin
                instances = aggregate.loc[aggregate['intensity_bin']
                                          == intensity_bin]
                # Add number of instances to dictionary
                bin_counts[model][experiment][intensity_bin] = len(instances)

    # For storm bin counts
    # Re-arrange dictionary to allow for MultiIndex DataFrame construction (MultiIndex is the model name, column is the experiment name, index is the intensity bin)
    storm_counts = {(model_key, experiment_key): {'count': experiment_values} for model_key, model_values in storm_counts.items()
                    for experiment_key, experiment_values in model_values.items()}
    # Build the DataFrame
    storm_counts = pd.DataFrame.from_dict(storm_counts, orient='columns')

    # For intensity bin counts
    # Re-arrange dictionary to allow for MultiIndex DataFrame construction (MultiIndex is the model name, column is the experiment name, index is the intensity bin)
    bin_counts = {(model_key, experiment_key): experiment_values for model_key, model_values in bin_counts.items()
                  for experiment_key, experiment_values in model_values.items()}
    # Build the DataFrame
    bin_counts = pd.DataFrame.from_dict(bin_counts, orient='columns')

    return storm_counts, bin_counts

def tc_activity_overlay(model_name, experiment, year_range, bin_resolution=2.5, FLOR_year_adjust=True,
                        FLOR_year_adjustment=1950):
    '''
    This method generates contour overlays on an existing figure to delineate TC activity spatial extent for the model and experiment over a given year range.
    This method doesn't rely on data from the script being called, because it pulls TC activity from a separate dataset.

    Arguments:
    - model_name (str):          name of the model
    - experiment (str):          name of the experiment
    - year_range (tuple of int): list of years over which TC activity is plotted
    - bin_resolution (numeric):  interval (in degrees) over which TC activity should be subdivided for density estimates

    Returns:
    - x (numeric array): list of unique longitudes corresponding to the track density output in densities
    - y (numeric array): list of unique longitudes corresponding to the track density output in densities
    - densities (dict):  dictionary containing track density in NumPy arrays with subdictionaries for model name and experiment

    '''

    start_time = time.time()
    checkpoint_counter = 0
    print('[tc_activity_overlay] Checkpoint {0}: {1:.3f} s elapsed'.format(
        checkpoint_counter, time.time() - start_time))
    # Get the start and end years from the year_range argument
    start_year, end_year = min(year_range), max(year_range)
    # Get track data for all storms (storm type 'TS') for the given model name and experiment
    track_data = tc_track_data([model_name], [experiment], year_range=(start_year, end_year), storm_type='TS',
                               FLOR_year_adjust=FLOR_year_adjust, FLOR_year_adjustment=FLOR_year_adjustment, diagnostic=False)

    ''' Data processing. '''
    # Initialize dictionary for spatial density
    density = {}
    # Iterate over each model provided
    for model in track_data.keys():
        print('[tc_activity_overlay] Checkpoint 1: {0:.3f} s elapsed'.format(
            time.time() - start_time))
        # Initialize model-specific subdictionary
        density[model] = {}
        # Iterate over each experiment provided
        for experiment in track_data[model].keys():
            # Define storage list for TC data
            dataset = []
            # Iterate through each unique TC and resample by day to get number of TC days per grid point
            for storm_id in track_data[model][experiment]['raw']['storm_id'].unique():
                dataset.append(track_data[model][experiment]['raw'].loc[track_data[model][experiment]
                               ['raw']['storm_id'] == storm_id].resample('D', on='time').first().reset_index())
            dataset = pd.concat(dataset)
            # Define dictionary to hold data relevant to track density heatmap
            counts = {'lon': [], 'lat': [], 'count': [], 'num_storms': []}
            # Group DataFrame into defined bins and add storm counts to the dictionary
            for lon_g, lon_v in dataset.groupby(pd.cut(dataset['center_lon'], np.arange(0, 360+bin_resolution, bin_resolution))):
                for lat_g, lat_v in lon_v.groupby(pd.cut(lon_v['center_lat'], np.arange(-90, 90+bin_resolution, bin_resolution))):
                    counts['lon'].append(lon_g.left)
                    counts['lat'].append(lat_g.left)
                    counts['count'].append(len(lat_v))
                    counts['num_storms'].append(len(dataset))
            # Create DataFrame from the dictionary
            counts = pd.DataFrame(counts)
            # Add time metadata to the DataFrame for bin per year estimate
            counts['time_min'], counts['time_max'] = dataset['time'].min(
            ), dataset['time'].max()
            # Concatenate to get comprehensive DataFrame
            density[model][experiment] = counts

    # Collect processed density maps. Note: the difference will be calculated for the last subplot.
    densities = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for model in [model_name]:
        print('[tc_activity_overlay] Checkpoint 2: {0:.3f} s elapsed'.format(
            time.time() - start_time))
        # Initialize subdictionary for the iterand model
        densities[model] = {}
        # Get longitude and latitude bins
        x, y = density[model][experiment].lon.unique(
        ), density[model][experiment].lat.unique()
        # Get density array
        v = np.reshape(density[model][experiment]
                       ['count'].values, (len(x), len(y)))
        # Assign to dictionary for the given model/experiment configuration.
        # Normalize by number of years and by day (assume 6-hourly data, so 6 x 4 = 1 day)
        densities[model][experiment] = v.T/(end_year - start_year)

    return x, y, densities

def track_data_loading(models=['HIRAM', 'AM2.5', 'FLOR'], experiments=['CTL1990s', 'CTL1990s_swishe'],
                       storm_type='TS', year_range=(101, 151), save_data=False, diagnostic=True):

    dirname = '/projects2/GEOCLIM/gr7610/analysis/tc_storage/track_data'
    # Generate strings and filename for the potential output file
    models_str = '-'.join(models)
    experiments_str = '-'.join(experiments)
    date_str = datetime.datetime.today().strftime('%Y%m%d')
    filename = 'TC_track_data.s{0}_e{1}.models-{2}.experiments-{3}.storm_type-{4}.generated-{5}.pkl'.format(
        min(year_range), max(year_range), models_str, experiments_str, storm_type, date_str)
    # Generate pathname
    pathname = os.path.join(dirname, filename)
    # Get list of all files in the target directory to search for pre-generated data
    existing_files = [f for f in os.listdir(dirname) if f.endswith('pkl')]

    if diagnostic:
        print('[tc_analysis.py, track_data_loading()] The requested pathname would be: {0}'.format(
            pathname))
    # This conditions looks to see if the requested track data exists as-is
    if os.path.isfile(pathname):
        print('[tc_analysis.py, track_data_loading()] File exists, opening it now...'.format(
            pathname))
        with open(pathname, 'rb') as handle:
            track_data = pickle.load(handle)

        return track_data
    # This conditional looks to see if the requested track data exists as a subset of existing data to avoid regeneration
    elif len(existing_files) > 0:
        for existing_file in existing_files:
            if diagnostic:
                print('\t[tc_analysis.py, track_data_loading()] Existing file name: {0}'.format(
                    existing_file))
            # Get the model string for the iterand file
            existing_model_str = sorted(existing_file.split(
                '.models-')[-1].split('.experiments')[0].split('-'))
            # Get the beginning and end years
            existing_start_year, existing_end_year = [int(existing_file.split('.')[1].split('_')[0][1:]),
                                                      int(existing_file.split('.')[1].split('_')[1][1:])]
            # Get the list of experiments
            existing_experiments = sorted(existing_file.split(
                'experiments-')[1].split('.')[0].split('-'))
            # Use the conditions
            model_name_filter = set(models).issubset(existing_model_str)
            time_filter = min(year_range) >= existing_start_year and max(
                year_range) <= existing_end_year
            experiment_filter = set(experiments).issubset(existing_experiments)

            if diagnostic:
                print('\t[tc_analysis.py, track_data_loading()] Model filter: {1}; time filter: {2}; experiment filter: {3}'.format(existing_file, model_name_filter,
                                                                                                                                    time_filter, experiment_filter))
            # If all criteria are met, load the data.
            if model_name_filter and time_filter and experiment_filter:
                existing_path = os.path.join(dirname, existing_file)
                print('[tc_analysis.py, track_data_loading()] File does exist, opening it now...'.format(
                    existing_path))
                with open(existing_path, 'rb') as handle:
                    track_data = pickle.load(handle)
                break
        return track_data
    # Else, generate
    else:
        print('[tc_analysis.py, track_data_loading()] File does not exist, making it now...'.format(
            pathname))
        track_data = tc_track_data(models, experiments, storm_type=storm_type, year_range=year_range,
                                   FLOR_year_adjustment=1900, parallel_load=True, diagnostic=False)
        if save_data:
            with open(pathname, 'wb') as handle:
                pickle.dump(track_data, handle,
                            protocol=pickle.HIGHEST_PROTOCOL)

def TC_statistics(track_data: pd.DataFrame):
    
    ''' Helper function to append TC statistics (e.g., ACE, PDI, etc) to each storm. '''
    
    # Append accumulated cyclone energy [m^2 s^-1]
    # Method: get squared sum of each storm's maximum winds and multiply it by the duration (in days) and correct for units
    track_data['ACE'] = track_data.groupby('storm_id').apply(lambda x: ((x['max_wind']**2).sum())).rename('ACE').reset_index(drop=True)
    # Append power dissipation index [m^3 s^-2]
    # Method: get cubed sum of each storm's maximum winds and multiply it by the duration (in days) and correct for units
    track_data['PDI'] = track_data.groupby('storm_id').apply(lambda x: ((x['max_wind']**3).sum())).rename('PDI').reset_index(drop=True)
    
    return track_data

def describe(model_name: str, 
             experiment_name: str, 
             track_data: pd.DataFrame):
    
    ''' Information function used to provide statistics for TCs in an experiment. '''

    # Get number of TCs
    storm_count = len(track_data['raw']['storm_id'].unique())
    # Get approximate number of storms per year
    number_of_years = np.round((track_data['raw']['cftime'].max() - track_data['raw']['cftime'].min()).total_seconds() / 86400 / 365) # get approximate number of years
    storm_count_per_year = np.round(storm_count / number_of_years)
    # Get duration of TCs
    storm_duration_raw = track_data['raw'].groupby('storm_id').first()['duration']
    storm_duration = storm_duration_raw[np.abs(sp.stats.zscore(storm_duration_raw) < 5)] # remove extreme outliers that may be a data processing anomaly
    # Get intensity information
    storm_max_wind = track_data['raw'].groupby('storm_id')['max_wind'].max()
    storm_min_pressure = track_data['raw'].groupby('storm_id')['min_slp'].min()
    # Get total run energy information
    # storm_ACE = track_data['raw'].groupby('storm_id')['ACE'].first()
    # storm_PDI = track_data['raw'].groupby('storm_id')['PDI'].first()

    print('-------------------------------------------------------------')
    print(f'Statistics for TCs in model: {
          model_name}; experiment: {experiment_name}')
    print(f'Number of storms: total = {storm_count}; per year = {storm_count_per_year}')
    print(f'Storm duration: mean = {storm_duration.mean():.2f} +/- {storm_duration.std():.2f} days')
    print(f'Storm maximum winds: mean = {storm_max_wind.mean():.2f} +/- {storm_max_wind.std():.2f} m/s')
    print(f'Storm minimum pressure: mean = {storm_min_pressure.mean():.2f} +/- {storm_min_pressure.std():.2f} hPa')
    # print(f'Run-integrated global energy statistics: ACE = {storm_ACE.sum():.2e} m^2 s^-1; PDI = {storm_PDI.sum():.2e} m^3 s^-2')
    print('-------------------------------------------------------------\n')

def load_TC_tracks(model: str,
                   experiment: str,
                   year_range: tuple[int, int],
                   month_range: tuple[int, int]=(1, 13),
                   storm_type: str = 'TS',
                   print_statistics: bool=False,
                   diagnostic: bool = False):

    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/track_data'

    model_identifier = f'model_{model}.'
    experiment_identifier = f'experiment_{experiment}.'

    if diagnostic:
        print(f'[load_TC_tracks] model name: {model}, experiment name: {experiment}, year_range: {year_range}')

    for filename in os.listdir(storage_dirname):
        filename_min_year_str, filename_max_year_str = filename.split('.')[1].split('_')
        filename_min_year, filename_max_year = filename_min_year_str.strip('s'), filename_max_year_str.strip('e')

        # Define conditions for data to meet
        model_flag = model_identifier in filename
        experiment_flag = experiment_identifier in filename
        year_range_flag = min(year_range) >= int(filename_min_year) and max(year_range) <= int(filename_max_year)

        # If all conditions are met, load the data
        if model_flag and experiment_flag and year_range_flag:
            pathname = os.path.join(storage_dirname, filename)
            print(f'Loading data for {model}, experiment {experiment} from {pathname}...')
            with open(pathname, 'rb') as f:
                track_data = pickle.load(f)
            
            for key, value in track_data[model][experiment].items():
                # Get cftime calendar time from iterand dataset, assuming all timestamp types are equivalent
                calendar_type = value['cftime'].iloc[0].calendar

                # Filter the data by year range by creating cftime objects
                track_year_min_cftime, track_year_max_cftime = [cftime.datetime(year=min(year_range), month=1, day=1, calendar=calendar_type),
                                                                cftime.datetime(year=max(year_range), month=1, day=1, calendar=calendar_type)]
                
                if diagnostic:
                    print(f'Obtaining data from potential years {track_year_min_cftime} to {track_year_max_cftime}.')
            
                value = value.loc[(value['cftime'] >= track_year_min_cftime) &
                                  (value['cftime'] <= track_year_max_cftime)]
                
                # Filter by month
                # Note the inconsistency in time variables. Pandas datetime objects are easier for month filtering and month should be agnostic between time and cftime.
                track_data[model][experiment][key] = value.loc[(value['time'].dt.month >= min(month_range)) & 
                                                               (value['time'].dt.month < max(month_range))]
                
            # Filter all storms for hurricanes, if the corresponding storm type (`C15w`) is selected
            if storm_type == 'C15w':
                # per Harris et al. (2016), doi.org/10.1175/JCLI-D-15-0389.1
                hurricane_threshold = {'field': 'max_wind', 'value': 32}
                for dataset_key, dataset_value in track_data[model][experiment].items():
                    hurricane_storm_IDs = dataset_value.loc[dataset_value[hurricane_threshold['field']] > hurricane_threshold['value']]['storm_id'].unique()
                    track_data[model][experiment][dataset_key] = dataset_value.loc[dataset_value['storm_id'].isin(hurricane_storm_IDs)]

            # Print statistics
            if print_statistics:
                describe(model, experiment, track_data[model][experiment])

    
    return track_data

def TC_density(model_names: str | list,
               experiment_names: str | list,
               year_range: tuple[int, int],
               month_range: tuple[int, int]=(1, 13),
               year_adjustment: int=0,
               bin_resolution: int=5,
               storm_type: str = 'TS'):
    """
    Method to plot the daily spatial density of all TC occurrences given a track data dictionary. 
    See tc_analysis.py --> tc_track_data() for more information.

    Args:
        data (dictionary): 3-tiered dictionary.
        model_name (list): list of strings with names of model (usually 'AM2.5', 'HIRAM', 'FLOR'.)
        bin_resolution (int, optional): resolution at which to generate spatial density, in degrees. Defaults to 5.
    """

    # Ensure input data is of proper form
    if isinstance(model_names, str):
        model_names = [model_names]
    elif isinstance(model_names, list) or isinstance(model_names, tuple):
        model_names = model_names

    if isinstance(experiment_names, str):
        experiment_names = [experiment_names]
    elif isinstance(experiment_names, list) or isinstance(experiment_names, tuple):
        experiment_names = experiment_names

    # Load track data
    data = {}
    for model in model_names:
        data[model] = {}
        for experiment in experiment_names:
            model_year_range = tuple([year + year_adjustment for year in year_range]) if model == 'IBTrACS' else year_range
            experiment_name = '' if model == 'IBTrACS' else experiment
            print(model, experiment_name, model_year_range)
            data[model][experiment_name] = load_TC_tracks(model, 
                                                          experiment_name, 
                                                          model_year_range, 
                                                          month_range=month_range, 
                                                          storm_type=storm_type)[model][experiment_name]

    # Define binning windows
    longitude_bins = np.arange(0, 360 + bin_resolution, bin_resolution)
    latitude_bins = np.arange(-90, 90 + bin_resolution, bin_resolution)

    ''' Data processing. '''
    # Initialize dictionary for spatial density
    density = {}
    # Iterate over each model provided
    for model in data.keys():
        # Initialize model-specific subdictionary
        density[model] = {}
        # Iterate over each experiment provided
        for experiment in data[model].keys():
            dataset = data[model][experiment]['raw']
            # Cut DataFrame by latitude and longitude bins
            densities = dataset.groupby([pd.cut(dataset['center_lat'], latitude_bins), pd.cut(
                dataset['center_lon'], longitude_bins)]).count()['storm_id'].unstack()
            densities.index, densities.columns = [
                densities.index.categories.left.values, densities.columns.categories.left.values]
            # Concatenate to get comprehensive DataFrame
            # Divide by 4, given that the output unit is in days per year and data is provided at a 4x daily frequency
            # Divide by number of years, given that the output unit is in days per year and data is provided over all years
            density[model][experiment] = densities / 4 / (max(year_range) - min(year_range))

    return density
