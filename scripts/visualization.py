import cartopy, cartopy.crs as ccrs
import numpy as np
import matplotlib, matplotlib.pyplot as plt

def cmap_white_adjust(cmap):
    ''' Adjust a given sequential monochromatic colormap to go from white to the monochrome. 
    
    Arguments:
    - cmap (str): colormap of choice
    Returns:
    - cm (ListedColormap)
    '''
    # Get the data for the given colormap
    cm = matplotlib.colormaps[cmap].resampled(256)
    # Get numerical entries at the resample points
    cm = cm(np.arange(cm.N))
    # Replace the first instance with a pure white entry
    cm = cm[:, :] + np.array([0, 1-cm[0, 1], 1-cm[0, 2], 0])
    # Cast the adjusted colormap as a new map
    cm = matplotlib.colors.ListedColormap(cm)
    
    return cm