import cftime, datetime
import cartopy, matplotlib, matplotlib.pyplot as plt
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle, random

import composite, utilities

def load(models, experiments, num_storms=-1):
    
    """
    Method to load data into a 3-tiered dictionary:
    --> (1) model name, (2) experiment type, (3) TC data with [a] track output, [b] model output (planar), [c] model output (vertical)

    Args:
        models (list): list of strings with model names
        experiments (list): list of strings with model names
        num_storms (int): number of storms to process. -1 is default and processes all storms found.
    Returns:
        data (dict): 3-tiered dictionary.
    """
    
    # Define directory where processed, single-storm data is kept
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed'
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
            storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(dirname)
                        if (model in f) and (experiment in f)][:num_storms]
            print('\t Storms to be processed: {0}'.format(storm_ids))
            # Initialize a storage list for the storms 
            data[model][experiment] = {'data': []}
            # Iterate over all storms found
            for storm_id in storm_ids:
                # Access the processed data
                storm_filename, storm = utilities.access(model, experiment, storm_type='C15w', storm_id=storm_id, processed=True)
                # Obtain the radial and tangential wind components
                # Note: I don't like having two methods lumped into one (in this case, loading and obtaining wind components),
                #       but it's the easiest way to obtain radius-dependent data
                print('\t \t Processing {0} from {1}'.format(storm_id, storm_filename))
                # Append to storage list
                data[model][experiment]['data'].append(storm)
        
    return data

def counts(data):
    
    """
    Method to gather number of storms in each experiment.

    Args:
        data (dict): 3-tiered dictionary: (1) model name, 
                                          (2) experiment type,
                                          (3) TC data with [a] track output, [b] model output (planar), [c] model output (vertical)

    Returns:
        counts (Pandas DataFrame): table with data of interest
    """
    
    # Define intensity bins. This is subject to change every time the binning algorithm is updated.
    # Pre-definition is used to catalog instances where 0 records are found.
    intensity_bins = ['b{0}'.format(i) for i in range(0, 7)]
    
    # Initialize dictionaries: one for TC counts, one for snapshot counts by intensity bin
    storm_counts, bin_counts = {}, {}
    # Iterate over models
    for model in data.keys():
        storm_counts[model], bin_counts[model] = {}, {}
        # Iterate over experiments
        for experiment in data[model].keys():
            # Iterate through all track output data to obtain an aggregate DataFrame
            aggregate = pd.concat([data[model][experiment][i]['track_output'] for i in range(0, len(data[model][experiment]))])
            # Log how many individual TCs exist in this configuration
            storm_count = len(aggregate['storm_id'].unique())
            # Initialize dictionary exclusive to the model + experiment configuration
            storm_counts[model][experiment], bin_counts[model][experiment] = storm_count, {}
            # Get number of timestamps per intensity bin and append to dictionary
            for intensity_bin in intensity_bins:
                # Get instances in this intensity bin
                instances = aggregate.loc[aggregate['intensity_bin'] == intensity_bin]
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
    
def tc_track_data(models, experiments, storm_type='C15w', snapshot_type='lmi', year_range=None):
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
    
    # Initialize dictionary for gridded data
    data = {}
    # Iterate over each model
    for model in models:
        # Adjust years for FLOR, if the iterand model
        model_year_range = [year + 1900 for year in year_range] if model == 'FLOR' else year_range
        # Initialize dictionary for the model
        data[model] = {}
        # Iterate over each experiment
        for experiment in experiments:
            print('Processing {0} model data from the {1} experiment for {2} over years {3}'.format(model, experiment, storm_type, year_range))
            # Initialize dictionary for the model + experiment configuration
            # Note: two keys are added - raw (for raw data, includes all data for each TC) and unique (includes single data point for each TC)
            data[model][experiment] = {'raw': utilities.retrieve_tracked_TCs(model, experiment, storm_type, model_year_range), 'unique': None}
            # Get unique TC storm IDs from the raw data
            storm_ids = data[model][experiment]['raw']['storm_id'].unique()
            # Initialize list to hold unique TC data
            storms = []
            # Iterate over individual TCs to extract data for the unique key
            for storm_id in storm_ids:
                # Select data for the iterand storm ID
                storm = data[model][experiment]['raw'].loc[data[model][experiment]['raw']['storm_id'] == storm_id]
                # Get snapshot of TC based on argument
                storm = utilities.storm_snapshot(storm, mode=snapshot_type)
                ###
                # Storm specific filtering can be added here!
                ###
                # Append to holding list
                storms.append(storm)
            # Concatenate storms into single DataFrame
            data[model][experiment]['unique'] = pd.concat(storms).sort_values('storm_id')
            ### Read-across storm filtering
            # Filter out storms lasting more than 30 days
            data[model][experiment]['unique'] = data[model][experiment]['unique'].loc[data[model][experiment]['unique']['duration'] <= 30]
            
    return data

if __name__ == '__main__':
    models, experiments = ['HIRAM'], ['control', 'swishe']
    data = load(models, experiments, num_storms=5)