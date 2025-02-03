import numpy as np
import multiprocessing, os
import pandas as pd
import pickle, random
import xarray as xr

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

import derived, utilities
    
def derived_fields(data):
    
    ''' Control which derived fields are processed for individual storms. '''
    
    # Add horizontal wind
    data = derived.wind(data)
    # Add latent heat flux
    data = derived.lhflx(data)
    # Add net surface heat flux
    data = derived.thflx(data)
    # Add relative humidity
    data = derived.rh(data)
    # Add domainwise temperature anomaly
    data = derived.domainwise_anomaly(data, 'temp')
    # Add domainwise specific humidity anomaly
    data = derived.domainwise_anomaly(data, 'sphum')
    # Add moist static energy
    data = derived.mse(data)
    # Add radial and tangential velocity components
    data = derived.radial_tangential_velocities(data)
    # Add column-integrated moist static energy
    data = derived.vertical_integral(data, field='h')
        
    return data
    
def storage(data, filename=None, override=True):
    
    ''' 
    Store data for a given TC. 
    
    Note: this isn't in utilities.py since this method is script-specific for processed TCs and should only be used locally.
    '''
    
    # Define the storage filename if not already defined
    if not filename:
        # Define directory name and file name
        # Define filename using max wind and min SLP for future binning
        max_wind, min_slp = data['track_output']['max_wind'].max(), data['track_output']['min_slp'].min()
        # Extract storm ID
        storm_id = data['track_output']['storm_id'].unique().item()
        # Define storage filename and full path
        storage_filename = 'TC-{0}-{1}-{2}-{3}-max_wind-{4:0.0f}-min_slp-{5:0.0f}.pkl'.format(model, experiment, storm_type, storm_id, max_wind, min_slp)
    else:
        storage_filename = filename.split('/')[-1]
    
    # Store in subdirectory with processed data.
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    storage_path = os.path.join(storage_dirname, storage_filename)
        
    # If file doesn't exist, save
    if not os.path.isfile(storage_path) or override:
        with open(storage_path, 'wb') as f:
            print('Storing data at {0}'.format(storage_path))
            pickle.dump(data, f)

def main(model, experiment, storm_type, storm_id=None):
    print('Processing {0}'.format(storm_id))
    # Load data for a given storm
    filename, data = utilities.access(model, experiment, storm_type, storm_id)
    print(data['tc_vertical_output'].dims)
    # Assign intensity bins to each timestamp
    data = utilities.intensity_binning(mode='model_output', data=data, intensity_metric='max_wind')
    # Process data and obtain derived fields
    data = derived_fields(data)
    # Store processed data
    storage(data, filename=filename, override=False)
    
    return data
    
if __name__ == '__main__':
    # Define parallel implementation selection
    parallel_implementation = False
    # Define loading parameters
    models, experiments, storm_type = ['HIRAM-8xdaily'], ['swishe'], 'C15w'
    # Single storm load    
    storm_ids = None
    # Multi-storm load
    num_storms = -1 # enter -1 to get all
    # Create data load ratios to reflect total storm counts
    # file_counts = utilities.file_counter(models, num_storms)
    # print(file_counts)
    data_dirname, processed_dirname = ['/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs', 
                                       '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed']
    
    import importlib
    importlib.reload(derived)
    importlib.reload(utilities)
    
    # Load storms
    for model in models:
        for experiment in experiments:
            # initialize dictionary to hold binned IDs
            # storm_ids = {}
            # # Iterate over bins derived from 'file_counts'
            # for bin_index in range(0, len(file_counts[model][experiment])-1):
            #     # Define intensity bounds
            #     lower_bound, upper_bound = list(file_counts[model][experiment].keys())[bin_index], list(file_counts[model][experiment].keys())[bin_index + 1]
            #     # Filter storms by storm type
            #     if storm_type == 'TS':
            #         bin_ids_all = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(data_dirname) if (model in f) and (experiment in f) 
            #                     and (storm_type in f) and (int(f.split('-')[7]) >= lower_bound) and (int(f.split('-')[7]) <= upper_bound) and (int(f.split('-')[7]) < 34)]
            #     else:
            #         bin_ids_all = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(data_dirname) if (model in f) and (experiment in f) 
            #                     and (storm_type in f) and (int(f.split('-')[7]) >= lower_bound) and (int(f.split('-')[7]) <= upper_bound)]
            #     # Get random sample 
            #     bin_ids = random.sample(bin_ids_all, file_counts[model][experiment][lower_bound]) if file_counts[model][experiment][lower_bound] < len(bin_ids_all) else bin_ids_all
                
            #     # print('{0}, {1}: lower bound: {2:.1f} m/s, {3} files selected'.format(model, experiment, lower_bound, len(bin_ids)))
            #     storm_ids[(lower_bound, upper_bound)] = bin_ids
            
            
            # [print('{0}, {1} | bin {2} - {3}/{4} possible storms retrieved: {5}'.format(model, experiment, k, len(v), file_counts[model][experiment][k[0]], v)) for k, v in storm_ids.items()]
            
            # storm_ids = [sv for k, v in storm_ids.items() for sv in v]
            
            # # Check if storm_ids are already processed. If so, remove from list
            # processed_storms = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(processed_dirname)
            #                     if (model in f) and (experiment in f)]
            # unprocessed_ids = list(set(storm_ids) - set(processed_storms))
            # # print('{0}, {1} - storms to process: {2}'.format(model, experiment, sorted(unprocessed_ids)))  
            # # Get random shuffle of IDs
            # unprocessed_ids = random.sample(unprocessed_ids, num_storms) if num_storms < len(unprocessed_ids) and num_storms != -1 else unprocessed_ids
            
            if parallel_implementation:
                # Assemble the input lists for Pool.starmap
                inputs = [(model, experiment, storm_type, id_) for id_ in storm_ids]
                # Use 16 processors to process storms in parallel
                with multiprocessing.get_context("spawn").Pool(8) as p:
                    p.starmap(main, inputs)
            
            else:
                print(storm_ids)
                if storm_ids:
                    for storm_id in storm_ids:
                        print(storm_id)
                        if storm_id != ['2063-0002', '2072-0002']:
                            _ = main(model, experiment, storm_type, storm_id=storm_id)
                else:
                    storm_ids = [os.path.join(data_dirname, filename) for filename in os.listdir(data_dirname)
                                 if experiment in filename and model in filename]
                    storm_ids = ['-'.join(file.split('/')[-1].split('-')[5:7]) for file in storm_ids]
                    for storm_id in storm_ids:
                        print(storm_id)
                        if storm_id != ['2063-0002', '2072-0002']:
                            _ = main(model, experiment, storm_type, storm_id=storm_id)