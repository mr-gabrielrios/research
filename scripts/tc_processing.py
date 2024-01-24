import numpy as np
import os
import pandas as pd
import pickle, random
import xarray as xr

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

import derived, utilities
    
def derived_fields(model_name, data):
    
    ''' Control which derived fields are processed for individual storms. '''
    # Get storm year to determine which fields can be calculated.
    storm_year = int(data['track_output']['storm_id'].unique().item().split('-')[0])
    
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
    # Horizontal moisture divergence
    data = derived.scalar_divergence(data, field='h')
    # Add column-integrated moist static energy
    data = derived.vertical_integral(data, field='h')
    # Column-integrated horizontal MSE divergence
    data = derived.vertical_integral(data, field='flux_h')
    
    # If storm year is greater than or equal to 2052 for HIRAM runs, process all the fields
    if storm_year >= 2052 and model_name == 'HIRAM':
        # Columnwise net longwave radiation
        data = derived.net_lw(data)
        # Columnwise net shortwave radiation
        data = derived.net_sw(data)
        
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
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed'
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
    data = derived_fields(model, data)
    # Store processed data
    storage(data, filename=filename, override=False)
    
    return data
    
if __name__ == '__main__':
    # Define loading parameters
    models, experiments, storm_type = ['AM2.5'], ['swishe'], 'TS'
    # Single storm load    
    storm_ids = ['2052-0034']
    # Multi-storm load
    num_storms = 100 # enter -1 to get all
    # Load storms
    for model in models:
        for experiment in experiments:
            data_dirname, processed_dirname = ['/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs', 
                                               '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed']
            if storm_type == 'TS':
                storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(data_dirname)
                            if (model in f) and (experiment in f) and (storm_type in f) and (int(f.split('-')[7]) <= 32)]
            else:
                storm_ids = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(data_dirname)
                            if (model in f) and (experiment in f) and (storm_type in f)]
            
            print('{0}, {1} - all storms: {2}'.format(model, experiment, sorted(storm_ids)))
            
            # Check if storm_ids are already processed. If so, remove from list
            processed_storms = [f.split('-')[4] + '-' + f.split('-')[5] for f in os.listdir(processed_dirname)
                                if (model in f) and (experiment in f)]
            unprocessed_ids = list(set(storm_ids) - set(processed_storms))
            print('{0}, {1} - storms to process: {2}'.format(model, experiment, sorted(unprocessed_ids)))  
            # Get random shuffle of IDs
            unprocessed_ids = random.sample(unprocessed_ids, num_storms) if num_storms < len(unprocessed_ids) and num_storms != -1 else unprocessed_ids
            
            for storm_id in unprocessed_ids:
                if storm_id != ['2063-0002', '2072-0002']:
                    data = main(model, experiment, storm_type, storm_id=storm_id)