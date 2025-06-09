import itertools
import calendar
import pandas as pd
import numpy as np
import os
import xarray
import cftime
import sys
import time
from multiprocessing import Pool

import ibtracs
import utilities

def table_constructor(storm_number: str,
                      storm_track: list,
                      track_year: int,
                      diagnostic: bool=False) -> pd.DataFrame:

    ''' Constructs a DataFrame. '''
    
    # Define column names for DataFrame construction
    columns = ['time', 'center_lon', 'center_lat', 'min_slp', 'max_wind', 'warm_core']
    # Construct the DataFrame
    storm_df = pd.DataFrame(data=storm_track, columns=columns)
    # Correct data types
    storm_df = storm_df.astype({column: 'float' for column in columns if column != 'time'})
    # Assign storm ID
    storm_df['storm_id'] = f'{track_year:04d}-{storm_number}'
    # Convert Pandas index to Pandas datetime for time element extraction
    storm_df['time'] = pd.to_datetime(storm_df['time'], format='%Y%m%d%H')
    # Derive storm_duration
    storm_df['duration'] = (storm_df['time'].max() - storm_df['time'].min()).total_seconds() / 86400
    # Check for leap dates - cftime is not handling them for compatibility with GCM outputs
    # Note 1: if a storm timestamp falls on a leap year, skip it. A better algorithm has to be made for this, but is okay for the time being.
    has_leap_year = storm_df['time'].apply(lambda x: (x.month == 2) & calendar.isleap(x.year))
    # If a leap year is in the storm data, return null
    # Else, create column with cftime time data type to allow for Pandas operations while retaining actual model date
    if has_leap_year.sum() == 0:
        
        storm_df['cftime'] = storm_df['time'].apply(lambda x: cftime.datetime(year=track_year,  
                                                                              month=x.month,
                                                                              day=x.day,
                                                                              hour=x.hour,
                                                                              calendar='noleap',
                                                                              has_year_zero=True))
        return storm_df
    else:
        if diagnostic:
            print(f'[tc_tracks.table_constructor()] Leap day found in storm ID {track_year:04d}-{storm_number} - skipping this storm.')
        return None

def storm_criteria(storm: pd.DataFrame) -> pd.DataFrame:

    ''' 
    Filter tropical cyclone track data using criteria from Harris et al. (2016).
    doi.org/10.1175/JCLI-D-15-0389.1
    '''

    storm_duration_threshold = 72 # units of hours
    storm_track_frequency = 6 # units of hours
    tropical_storm_wind_threshold = 17.5 # units of m s^-1
    tropical_storm_duration_threshold = 36 # units of hours
    warm_core_duration_threshold = 48 # units of hours

    # Criterion 1. storm duration is at least 72 hours.
    duration_criterion = (storm.time.max() - storm.time.min()).total_seconds() / 3600 >= storm_duration_threshold # units of hours
    
    # Criterion 2. 36 hours of consecutive tropical-storm strength winds with a warm core
    # Get subselection of timestamps with tropical-storm strength winds
    tropical_storm_filter = storm.loc[(storm['max_wind'] >= tropical_storm_wind_threshold) & (storm['warm_core'] > 0)]['time']
    # Find where filter has consecutive values meeting the filter (consecutive meaning the timedelta is 6 hours)
    consecutive_threshold_timestamps = tropical_storm_filter.diff().values == np.timedelta64(storm_track_frequency, 'h')
    # Get number of consecutive entries meeting the filter
    consecutive_threshold_timestamp_numbers = [sum(1 for _ in group) for key, group in itertools.groupby(consecutive_threshold_timestamps) if key]
    consecutive_threshold_timestamp_numbers = [0] if len(consecutive_threshold_timestamp_numbers) == 0 else consecutive_threshold_timestamp_numbers
    # If the longest span of consecutive entries is greater than 36 hours
    tropical_storm_criterion = max(consecutive_threshold_timestamp_numbers) > (tropical_storm_duration_threshold / storm_track_frequency)

    # Criterion 3. Cumulative warm core time of 48 hours
    warm_core_filter = len(storm.loc[storm['warm_core'] > 0]) > (warm_core_duration_threshold / storm_track_frequency)

    if duration_criterion and tropical_storm_criterion and warm_core_filter:
        return storm.loc[storm['warm_core'] >= 0]
    else:
        return None
    
def access(model_name: str,
           experiment_name: str,
           year_range: tuple[int, int],
           diagnostic: bool=False) -> dict:

    ''' Obtains pathnames for GFDL QuickTracks outputs as a function of model, experiment, and year range. '''

    # Container dictionary with keys of year names and values of corresponding GFDL QuickTracks pathnames
    track_pathnames = {}

    # Retrieve directory containing parent directory for track data
    root_dirname = utilities.directories(model=model_name, 
                                         experiment=experiment_name, 
                                         data_type='track_data')
    
    if diagnostic:
        print(f'[track_TCs.access] pulling tracks from root directory: {root_dirname}...')
        
    # Retrieve subdirectories containing GFDL QuickTracks model outputs corresponding to a given year
    # Conditions are that (1) each subdirectory must contain the 'Harris.TC' subsubdirectory and (2) year is within `year_range` bounds
    # Assumes that year subdirectories are named with format `atmos_YYYY_YYYY`, where YYYY is a year.
    track_dirnames = [os.path.join(root_dirname, year_dirname, 'Harris.TC') for year_dirname in os.listdir(root_dirname)
                   if 'Harris.TC' in os.listdir(os.path.join(root_dirname, year_dirname))
                   and int(year_dirname.split('_')[-1]) >= min(year_range)
                   and int(year_dirname.split('_')[-1]) <= max(year_range)]
    
    # Iterate over each annual subdirectory
    for track_dirname in track_dirnames:
        if diagnostic:
            print(f'[track_TCs.access] iterating over track directory: {track_dirname}...')
            
        # Retrieve names of data.
        track_filenames = [track_filename for track_filename in os.listdir(track_dirname)
                           if 'allstorms.world' in track_filename
                           and track_filename.endswith('.txt')]
        # Check on file length
        assert len(track_filenames) == 1, f'Number of TCs files matching the criterion is {len(track_filenames)}, must be 1.'
        # Build track path
        track_pathname = os.path.join(track_dirname, track_filenames[0])
        # Obtain year for which the path is constructed
        track_year = int(track_pathname.split('/atmos_')[1].split('/')[0].split('_')[0])
        # Add to container dictionary
        track_pathnames[track_year] = track_pathname

    return track_pathnames

def access_global(root_dirname: str,
           year_range: tuple[int, int],
           diagnostic: bool=False) -> dict:
    
    ''' 
    Obtains pathnames for GFDL QuickTracks outputs as a function of directory name and year range. 
    
    root_dirname: path to parent directory of Lucas Harris' TC tracker
    --- path must be one level above the directory containing `atmos_XXX` subdirectories
    --- this path will usually start with 'cyclones_gav'
    '''

    # Container dictionary with keys of year names and values of corresponding GFDL QuickTracks pathnames
    track_pathnames = {}
    
    if diagnostic:
        print(f'[track_TCs.access] pulling tracks from root directory: {root_dirname}...')
        
    # Retrieve subdirectories containing GFDL QuickTracks model outputs corresponding to a given year
    # Conditions are that (1) each subdirectory must contain the 'Harris.TC' subsubdirectory and (2) year is within `year_range` bounds
    # Assumes that year subdirectories are named with format `atmos_YYYY_YYYY`, where YYYY is a year.
    track_dirnames = [os.path.join(root_dirname, year_dirname, 'Harris.TC') for year_dirname in os.listdir(root_dirname)
                    if 'Harris.TC' in os.listdir(os.path.join(root_dirname, year_dirname))
                    and int(year_dirname.split('_')[-1]) >= min(year_range)
                    and int(year_dirname.split('_')[-1]) <= max(year_range)]
    
    if diagnostic:
        print(f'[track_TCs.access] years found in {root_dirname} for requested range {year_range}: {track_dirnames}...')
    
    # Iterate over each annual subdirectory
    for track_dirname in track_dirnames:
        if diagnostic:
            print(f'[track_TCs.access] iterating over track directory: {track_dirname}...')
            
        # Retrieve names of data.
        track_filenames = [track_filename for track_filename in os.listdir(track_dirname)
                           if 'allstorms.world' in track_filename
                           and track_filename.endswith('.txt')]
        # Check on file length
        assert len(track_filenames) == 1, f'Number of TCs files matching the criterion is {len(track_filenames)}, must be 1.'
        # Build track path
        track_pathname = os.path.join(track_dirname, track_filenames[0])
        # Obtain year for which the path is constructed
        track_year = int(track_pathname.split('/atmos_')[1].split('/')[0].split('_')[0])
        # Add to container dictionary
        track_pathnames[track_year] = track_pathname

    return track_pathnames

def storm_parser(track_year, track_pathname):
    
    annual_storm_tracks = {} # initialize container to hold year-specific storm tracks
    
    storm_delimiter = '+++' # string used to demarcate new storm
    annual_track_number = 0 # initialize track number
    storm_track = [] # initialize temporary list that will hold data for an individual storm and be reset
    
    # Open the text file for parsing
    with open(track_pathname, 'r') as file:
        for line in file:
            if storm_delimiter in line:
                # If there's data loaded and a delimiter is met, save the track data to its respective storm number
                if storm_track:
                    annual_storm_number = f'{annual_track_number:04d}'
                    # Construct the DataFrame
                    storm_table = table_constructor(storm_number=annual_storm_number,
                                                    storm_track=storm_track,
                                                    track_year=track_year)
                    # Filter the DataFrame per criteria defined in Lucas Harris' 2016 paper
                    if storm_table is not None:
                        storm_table = storm_criteria(storm_table)
                        annual_storm_tracks[annual_storm_number] = storm_table
                        annual_track_number += 1
                    # Reset storm track list
                    storm_track = []
            else:
                # Append item to list, separate each column (separated by strings) into list elements
                storm_track.append(line.strip().split())
    
    # Concatenate storm DataFrames if non-empty, else return empty list
    if len(annual_storm_tracks.values()) > 0:
        annual_storm_tracks = pd.concat(annual_storm_tracks.values())
        return annual_storm_tracks
    else:
        return []

def constructor(track_pathnames: dict,
                parallel: bool = False,
                diagnostic: bool = True) -> pd.DataFrame:

    ''' Builds a DataFrame concatenating filtered track data from all pathnames provided. '''
    
    # Begin profiling
    start_time = time.time()
    
    # Parse and filter data for each storm and concatenate into singular DataFrame
    # If the parallel choice is enabled, use multiprocessing.starmap
    if parallel:
        # Initialize root container for all track data
        storm_tracks = []
        # Profile time
        start_time = time.time()
        with Pool(processes=20) as pool:
            # Get track years and pathnames for zipping
            func_input = [(k, v) for k, v in track_pathnames.items()]
            # Parallelize the call to `storm_parser`
            annual_storm_tracks = pool.starmap(storm_parser, func_input)
            # Concatenate this Pool's DataFrames
            annual_storm_track_entries = [annual_storm_track_entry for annual_storm_track_entry in annual_storm_tracks if
                                          len(annual_storm_track_entry) > 0]
            annual_storm_tracks = pd.concat(annual_storm_track_entries)
            # Append to container list for ultimate concatenation
            storm_tracks.append(annual_storm_tracks)
        print(f'Parallelized elapsed time for data construction: {(time.time() - start_time):.2f}s')
        # Close the pool
        pool.close()
        # Concatenate the DataFrames
        storm_tracks = pd.concat(storm_tracks).sort_values('time')
    else:
        # Initialize root container for all track data
        storm_tracks = {}
        # Iterate over every provided year, parse data, format it, and filter it
        for track_year, track_pathname in track_pathnames.items():

            annual_storm_tracks = {} # initialize container to hold year-specific storm tracks
            storm_delimiter = '+++' # string used to demarcate new storm
            annual_track_number = 0 # initialize track number
            storm_track = [] # initialize temporary list that will hold data for an individual storm and be reset
            
            # Open the text file for parsing
            with open(track_pathname, 'r') as file:
                for line in file:
                    if storm_delimiter in line:
                        # If there's data loaded and a delimiter is met, save the track data to its respective storm number
                        if storm_track:
                            annual_storm_number = f'{annual_track_number:04d}'
                            # Construct the DataFrame
                            storm_table = table_constructor(storm_number=annual_storm_number,
                                                            storm_track=storm_track,
                                                            track_year=track_year)
                            # Filter the DataFrame per criteria defined in Lucas Harris' 2016 paper
                            if storm_table is not None:
                                storm_table = storm_criteria(storm_table)
                                if storm_table is not None:
                                    annual_storm_tracks[annual_storm_number] = storm_table
                                    annual_track_number += 1
                            # Reset storm track list
                            storm_track = []
                    else:
                        # Append item to list, separate each column (separated by strings) into list elements
                        storm_track.append(line.strip().split())
            
            # Concatenate storm DataFrames
            if len(annual_storm_tracks) > 0:
                storm_tracks[track_year] = pd.concat(annual_storm_tracks.values())
            else:
                print(f'No storms found that match criteria for year {track_year}.')
                storm_tracks[track_year] = None
            
            if diagnostic:
                print(f'[tc_tracks.constructor()] Elapsed time for year {track_year}: {(time.time() - start_time):.2f}s')

        # Concatenate annual DataFrames
        storm_tracks = pd.concat(storm_tracks.values())

    return storm_tracks

def load_IBTrACS(year_range: tuple[int, int]) -> pd.DataFrame:
    
    ''' Parse IBTrACS track data. '''
    
    # Ensure inputs are properly formatted. 
    # Note the use of string comparison insteado f `isinstance` for type checking - this is meant to catch general integer data types
    assert len(year_range) == 2, 'Year range must have 2 elements.'
    assert ('int' in str(type(year_range[0]))) & ('int' in str(type(year_range[1]))), 'Year range must have 2 integer elements.'
    # Construct date range string for IBTrACS function input
    date_range = tuple([f'{year}-01-01' for year in sorted(year_range)])
    
    storm_tracks = ibtracs.main(date_range=date_range)
    
    return storm_tracks

def main(year_range: tuple[int, int],
         model_name: str=None,
         experiment_name: str=None,
         pathname: str=None,
         run_parallel: bool=True,
         diagnostic: bool=False) -> pd.DataFrame:

    # Get storm track data for IBTrACS data
    if model_name == 'IBTrACS':
        storm_tracks = load_IBTrACS(year_range)
    # Get pathnames for GFDL QuickTracks data
    else:
        if (model_name is not None) and (experiment_name is not None):
            track_pathnames = access(model_name, experiment_name, year_range, diagnostic=diagnostic)
        elif pathname is not None:
            track_pathnames = access_global(pathname, year_range, diagnostic=diagnostic)
        else:
            print('A model and experiment name, or a valid directory path, must be provided to track TCs. Please retry.')
            sys.exit()
        # Get storm track data
        storm_tracks = constructor(track_pathnames, parallel=run_parallel, diagnostic=diagnostic)
    
    return storm_tracks

if __name__ == '__main__':
    model_name = 'AM2.5'
    experiment_name = 'CTL1990s'
    year_range = (101, 105)
    run_parallel = False
    main(year_range, model_name, experiment_name, run_parallel)