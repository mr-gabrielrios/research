''' Import packages. '''
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy
# Local utility packages
import dill, multiprocessing, os, pickle, utilities
# Data structure packages
import pandas as pd, xarray as xr
# Visualization tools
import cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

def main(model, storm_id, storage=False, override=False):

    # For a given TC, get the ocean data corresponding to the TC's location
    # This is mostly meant for FLOR. If run for AM2.5 and HIRAM, extract SST data from the Hadley SSTs

    # Get storm_id
    # Load track data

    # Iterate over experiments to load and process data