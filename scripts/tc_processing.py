import numpy as np
import os
import pandas as pd
import pickle, random
import xarray as xr

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

import derived

# Make any mathematical text render in LaTeX
matplotlib.rcParams['mathtext.fontset'] = 'cm'

def access(model, storm_type, storm_id=None):
    
    """Access a random TC pickle file.
    
    Arguments:
        model (str): model name (AM2.5, HIRAM, FLOR)
        storm_type (str): TS or C15w
        storm_id (str, default: None): ID of storm of interest

    Returns:
        data (dict): 3-element dictionary with track outputs from Lucas Harris' TC tracker, planar model outputs from the chosen GCM, and vertical outputs from the chosen GCM
    """
    
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    model = 'HIRAM'
    storm_type = 'C15w'
    files = [os.path.join(dirname, filename) for filename in os.listdir(dirname)
            if model in filename and storm_type in filename]
    
    # If a specific storm ID is given, check for it
    storm_exists = False # check for storm existence - True if exists
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
            data = pickle.load(f)
    else:
        with open(random.choice(files), 'rb') as f:
            data = pickle.load(f)
        
    return data
    
def derived_fields(data):
    
    # Add horizontal wind
    data = derived.wind(data)
    # Add latent heat flux
    data = derived.lhflx(data)
    # Add net surface heat flux
    data = derived.thflx(data)
    # Add relative humidity
    data = derived.rh(data)
    # Add moist static energy
    data = derived.mse(data)
    # Columnwise net longwave radiation
    data = derived.net_lw(data)
    # Columnwise net shortwave radiation
    data = derived.net_sw(data)
    # Horizontal moisture divergence
    data = derived.scalar_divergence(data, field='h')
    # Column-integrated horizontal MSE divergence
    data = derived.vertical_integral(data, field='flux_h')
    # Add column-integrated moist static energy
    data = derived.vertical_integral(data, field='h')
        
    return data
    
def main(model, storm_type, storm_id=None):
    # Load data for a given storm
    data = access(model, storm_type, storm_id)
    # Process data and obtain derived fields
    data = derived_fields(data)
    return data
    
if __name__ == '__main__':
    model, storm_type = 'HIRAM', 'C15w'
    data = main(model, storm_type, storm_id='2056-0039')