import numpy as np
import os
import pandas as pd
import pickle, random
import xarray as xr

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt

import derived, utilities

# Make any mathematical text render in LaTeX
matplotlib.rcParams['mathtext.fontset'] = 'cm'
    
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
    # Add column-integrated moist static energy
    data = derived.vertical_integral(data, field='h')
    # Column-integrated horizontal MSE divergence
    data = derived.vertical_integral(data, field='flux_h')
        
    return data
    
def main(model, experiment, storm_type, storm_id=None):
    # Load data for a given storm
    data = utilities.access(model, experiment, storm_type, storm_id)
    # Process data and obtain derived fields
    data = derived_fields(data)
    return data
    
if __name__ == '__main__':
    model, experiment, storm_type = 'HIRAM', 'control', 'C15w'
    data = main(model, experiment, storm_type, storm_id='2056-0039')