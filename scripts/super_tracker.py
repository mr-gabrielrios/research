import pickle
import numpy as np
import os
import xarray as xr
import matplotlib.pyplot as plt
import importlib

import tc_analysis, visualization, utilities, track_TCs

importlib.reload(tc_analysis);
importlib.reload(utilities);
importlib.reload(visualization);

def find_experiments(root_dirname: str='/projects/GEOCLIM/gr7610/MODEL_OUT',
                     maximum_year_length: int=200) -> dict:

    ''' 
    Crawler method that goes through the root directory `root_dirname` and finds experiments with TC tracks.
    The output will be inputs for the tc_analysis.tc_track_data method.

    Arguments:
    - root_dirname: parent directory from which experiments in subdirectories will be found.
    - maximum_year_length: maximum number of years that track data can be loaded for.
    Returns:
    - track_names: multi-level dictionary that provides inputs for the track data parser method. Dictionary structure below:

    Example:
    - track_names = {'model_name': {'experiment_name': (minimum_year, maximum_year)}}
    '''

    # Define strings corresponding to TC tracking algorithm subdirectories
    track_directory_identifier = 'analysis_lmh'
    track_subdirectory_identifier = 'cyclones_gav'
    
    # Specify substrings that belong to experiments that TC tracks shouldn't be processed for
    ignore_identifiers = ['8xdaily']

    # Initialize the dictionary
    track_names = {}
    # Search through the root directory provided
    for model_name in os.listdir(root_dirname):
        # Get pathname, save to check if it's a directory. Assume the pathname is a model.
        model_pathname = os.path.join(root_dirname, model_name)
        # If the path is a directory, assume it's a model. 
        # Step into it and search through model subdirectories
        if os.path.isdir(model_pathname):
            track_names[model_name] = {}
            # Iterate over assumed model directories
            for experiment_name in os.listdir(model_pathname):
                # Ignore files with specific identifiers in the name
                ignore_boolean = sum([identifier in experiment_name for identifier in ignore_identifiers])
                
                # Get pathname, save to check if it's a directory
                experiment_pathname = os.path.join(model_pathname, experiment_name)
                # If the path is a directory, assume it's an experiment. 
                # Step into it and search through experiment subdirectories
                if os.path.isdir(experiment_pathname) and ignore_boolean == 0:
                    # Search through experiment subdirectories
                    for experiment_subdirectory_name in os.listdir(experiment_pathname):
                        track_directory = os.path.join(experiment_pathname, experiment_subdirectory_name)
                        # The track directory should be an experiment subdirectory. 
                        # Check here for a subdirectory or symbolic link that matches the identifier.
                        if (os.path.isdir(track_directory) or (os.path.islink(track_directory))) and experiment_subdirectory_name == track_directory_identifier:
                            for track_subdirectory_name in os.listdir(track_directory):
                                # Make sure identified subdirectory is a parent directory for tracked TC data.
                                if track_subdirectory_identifier in track_subdirectory_name:
                                    track_subdirectory = os.path.join(track_directory, track_subdirectory_name)
                                    # Get minimum and maximum years in the listed track data
                                    track_year_range = [int(f.split('_')[-1]) for f in os.listdir(track_subdirectory)]
                                    # Make sure year range isn't too long, redefine lengths based on model type
                                    if (max(track_year_range) - min(track_year_range)) > maximum_year_length:
                                        track_year_range = (2001, 2200) if model_name == 'FLOR' else (101, 300)
                                    # Assign year range to the model and experiment configuration entry
                                    track_names[model_name][experiment_name] = (min(track_year_range), max(track_year_range))

    return track_names

def save_tracks(model_name: str,
                experiment_name: str,
                year_range: tuple[int, int],
                track_dataset: dict,
                pathname: str,
                overwrite: bool=False):

    with open(pathname, 'wb') as f:
        pickle.dump(track_dataset, f)

def load(track_names: dict, 
         overwrite: bool=False) -> dict:

    '''
    Method that loads TC track data given a dictionary with track names output by `find_experiments()`.
    '''
    
    track_datasets = {}
    for model in track_names.keys():
        track_datasets[model] = {}
        for experiment in track_names[model].keys():
            year_range = track_names[model][experiment]
            print(f'[load()] Model {model}; experiment {experiment}; year range = {year_range}')

            root_dirname: str='/projects/GEOCLIM/gr7610/analysis/tc_storage/track_data'
            filename = f'TC_track_data.s{min(year_range)}_e{max(year_range)}.model_{model}.experiment_{experiment}.pkl'
            pathname = os.path.join(root_dirname, filename)
        
            if os.path.isfile(pathname) and not overwrite:
                print(f'Data already exists for {pathname}. Please set the `overwrite` argument to `True` if you wish to overwrite this file.')
            else:
                track_datasets[model][experiment] = tc_analysis.tc_track_data(models=[model], 
                                                                              experiments=[experiment], 
                                                                              year_range=year_range, 
                                                                              storm_type='TS',
                                                                              parallel_load=False)
    
                save_tracks(model_name=model,
                            experiment_name=experiment,
                            year_range=year_range,
                            track_dataset=track_datasets[model][experiment],
                            pathname=pathname)            

    return track_datasets

def main():
    track_names = find_experiments()
    track_datasets = load(track_names)

if __name__ == '__main__':
    main()