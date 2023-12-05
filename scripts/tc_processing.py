import numpy as np
import os
import pandas as pd
import pickle, random
import xarray as xr

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

import multiprocess

# Make any mathematical text render in LaTeX
matplotlib.rcParams['mathtext.fontset'] = 'cm'

def access(model, storm_type):
    
    """Access a random TC pickle file.
    
    Arguments:
        model (str): model name (AM2.5, HIRAM, FLOR)
        storm_type (str): TS or C15w

    Returns:
        data (dict): 3-element dictionary with track outputs from Lucas Harris' TC tracker, planar model outputs from the chosen GCM, and vertical outputs from the chosen GCM
    """
    
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs'
    model = 'HIRAM'
    storm_type = 'C15w'
    files = [os.path.join(dirname, filename) for filename in os.listdir(dirname)
            if model in filename and storm_type in filename]

    with open(random.choice(files), 'rb') as f:
        data = pickle.load(f)
        
    return data
    
def derived_fields(data):
    
    if 'wind' not in data['tc_model_output'].data_vars:
        data['tc_model_output']['wind'] = np.sqrt(data['tc_model_output']['u_ref']**2 + data['tc_model_output']['v_ref']**2)
        data['tc_model_output']['wind'].attrs = {'long_name': 'horizontal wind at 10 m', 'units': 'm/s'}
    
    return data
    
def main(model, storm_type):
    # Load data for a given storm
    data = access(model, storm_type)
    return data
    
if __name__ == '__main__':
    model, storm_type = 'HIRAM', 'C15w'
    data = main(model, storm_type)