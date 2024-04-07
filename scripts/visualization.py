import cartopy, cartopy.crs as ccrs
import numpy as np, os, pandas as pd, scipy as sp, xarray as xr
import cmap as cmap_pkg, matplotlib, matplotlib.pyplot as plt

import accessor, composite, utilities, tc_analysis

class ScalarFormatterForceFormat(matplotlib.ticker.ScalarFormatter):
    ''' This method is an override to the ScalarFormatter used to create a common scientific exponent for an axis. '''
    def _set_format(self):  # Override function that finds format to use.
        self.format = '%.3g'

def cycler(i):
    
    colors = ['b', 'g', 'r', 'm']
    linestyles = ['-', ':', '--', '-.']
    property_dict = {'c': colors[i], 'ls': linestyles[i]}
    
    return property_dict

def subplot_title(ax, model, start_year, end_year, metadata_str=None):
    # Subplot labeling
    title_y = 1.05
    # Define year range string plotting
    year_range_str = ', {0} to {1}'.format(start_year, end_year)
    # Set left-hand side to be {model name}, {min year} to {max year}
    title = ax.annotate('{0}{1}'.format(model, year_range_str), (0, title_y), 
                        va='baseline', ha='left', xycoords='axes fraction', fontsize=10)
    # Set right-hand side to be {metadata}
    if metadata_str:
        metadata = ax.annotate(metadata_str, (1, title_y), va='baseline', ha='right', xycoords='axes fraction', fontsize=10)

    return ax

def cmap_white_adjust(cmap, levels=16):
    ''' Adjust a given sequential monochromatic colormap to go from white to the monochrome. 
    
    Arguments:
    - cmap (str): colormap of choice
    Returns:
    - cm (ListedColormap)
    '''
    # Total number of sampling points
    samples = 256
    # Get the data for the given colormap
    cm = matplotlib.colormaps[cmap].resampled(samples)
    # Get numerical entries at the resample points
    cm = cm(np.arange(cm.N))
    # Replace the first instance with a pure white entry
    # cm = cm[:, :] + np.array([1-cm[:10, 0], 1-cm[:10, 1], 1-cm[:10, 2], 1-cm[:10, -1]])
    cm[:samples//(levels-2), :] = np.array([1, 1, 1, 1])
    # cm[:samples//(levels//2), :] = np.array([1, 1, 1, 1])
    # Cast the adjusted colormap as a new map
    cm = matplotlib.colors.ListedColormap(cm)
    
    return cm

def get_cmap(field, norm=None):
    """
    Method to pull colormap for a given field from a .csv with colormap info.

    Args:
        field (str): field of interest
        norm (matplotlib.colors.Normalize): normalization with data to extract values from

    Returns:
        colormap: matplotlib.colors.colormap
    """
    
    # print('Colormap field: {0}'.format(field))
    # Get colormap DataFrame from .csv
    colormaps = pd.read_csv('/projects/GEOCLIM/gr7610/reference/param_colormaps.csv')
    # Get field-specific colormap. If field is not in colormaps, default to 'other'.
    if field in colormaps['param'].unique():
        colormap = colormaps.loc[colormaps['param'] == field]
    else:
        colormap = colormaps.loc[colormaps['param'] == 'other']
    
    if norm:
        return colormap['difference'].item() if (norm.vmin < 0 and norm.vmax > 0) else colormap['normal'].item()
    else:
        return colormap['normal'].item()

def diverging_colormap_adjust(cmap, fraction=0.1):
    """
    Insert white into the center of a diverging colormap using an arctan normalization between 0 and 1.
    The arctan function was chosen to add more anchor points, or nodes, near the edges to preserve the colormap.

    Args:
        cmap (matplotlib.colors.LinearSegmentedColormap): matplotlib-registered colormap
        fraction (float, optional): width of the color scale to make white, centered at halfway. Defaults to 0.1.

    Returns:
        cmap (matplotlib.colors.LinearSegmentedColormap): matplotlib-registered colormap: white-adjusted colormap
    """
    # Define the endpoints for the arctan curve and number of sampling points
    endpoint, num_points = 10, 100
    # Initialize list of nodes
    nodes = [(np.arctan(i)+np.pi/2)/np.pi for i in np.linspace(-endpoint, endpoint, num_points)]
    # Initialize list of colors
    colors = [(1, 1, 1) if (node > 0.5 - fraction/2) and (node < 0.5 + fraction/2) 
              else matplotlib.colormaps[cmap](node) for i, node in enumerate(nodes)]
    # Re-anchor points at extrema to satisfy matplotlib requirements for colormaps
    nodes[0], nodes[-1] = 0, 1
    
    return matplotlib.colors.LinearSegmentedColormap.from_list(cmap, list(zip(nodes, colors)))

def norm_cmap(data, field, num_bounds=16, extrema=None, white_adjust=False, cmap_white_fraction=0.1):
    """
    Method to derive normalization and colormap for a given list of numeric arrays and field.

    Args:
        data (_type_): _description_
        field (_type_): _description_
    """
    
    # Ensure data is a list 
    if type(data) is not list:
        # Handle dictionary case by extracting values
        if type(data) is dict:
            data = [v for v in data.values()]
        else:
            data = [data]
    
    # Initialize minimum and maximum values
    vmin, vmax = np.nan, np.nan
    # Generate extrema if not prescribed
    if not extrema:
        # Iterate over each data item
        for item in data:
            # Find extrema for this item
            vmin_, vmax_ = np.nanmin(item), np.nanmax(item)
            # If the extrema are larger than their corresponding extrema, update
            vmin, vmax = [vmin_ if (vmin_ <= np.nanmin([vmin, vmin_])) else vmin, 
                          vmax_ if (vmax_ > np.nanmax([vmax, 0])) else vmax]
            vmin = vmin_ if np.isnan(vmin) else vmin
            vmax = vmax_ if np.isnan(vmax) else vmax
    else:
        # Define minimum and maximum values
        vmin, vmax = min(extrema), max(extrema)
        bins = np.linspace(min(extrema), max(extrema), num_bounds+1)
        
    # Subdivide normalization into sequential (all one sign) or diverging (two signs, including 0) bins
    if vmin < 0 and vmax > 0:
        # Get larger of the two bounds
        extremum = max([abs(vmin), abs(vmax)])
        bins = np.concatenate([np.linspace(-extremum, 0, int(num_bounds/2))[:-1], [0], np.linspace(0, extremum, int(num_bounds/2))[1:]])
    else:
        bins = np.linspace(vmin, vmax, num_bounds+1)
    
    # Assign to normalization dictionary
    norm = matplotlib.colors.BoundaryNorm(bins, 256)
    # Get colormap
    cmap = get_cmap(field, norm)
    # If the colormap is a "white-to-color" sequential colormap, apply an adjustment for a white basis
    sequential_cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds','YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 
                        'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    # Perform white adjustment for colormaps (set lowest color to white for sequential, center to white for diverging)
    if white_adjust:
        cmap = cmap_white_adjust(cmap, levels=num_bounds) if cmap in sequential_cmaps else diverging_colormap_adjust(cmap, fraction=cmap_white_fraction)
        
    return norm, cmap

def field_properties(field):
    """
    Method to generate long name and units for a given field. Returns the field name and no units if field doesn't have a property entry.

    Args:
        field (str): GCM field name.

    Returns:
        long_name, units: str, str
    """
    
    properties = {'lhflx': {'long_name': 'latent heat flux', 'units': 'W m$^{-2}$'},
                  'shflx': {'long_name': 'sensible heat flux', 'units': 'W m$^{-2}$'},
                  'olr': {'long_name': 'outgoing longwave radiation', 'units': 'W m$^{-2}$'},
                  'netrad_toa': {'long_name': 'net radiation, TOA', 'units': 'W m$^{-2}$'},
                  'wind': {'long_name': '10m horizontal wind speed', 'units': 'm s$^{-1}$'},
                  'WVP': {'long_name': 'column-integrated water vapor\n', 'units': 'kg m$^{-2}$'},
                  'h': {'long_name': 'moist static energy', 'units': 'J kg$^{-1}$'},
                  'slp': {'long_name': 'sea level pressure', 'units': 'hPa'},
                  'omega': {'long_name': 'pressure velocity', 'units': 'Pa s$^{-1}$'},
                  'rh': {'long_name': 'relative humidity', 'units': '%'},
                  'vi_h': {'long_name': 'vertically-integrated MSE', 'units': 'J m$^{-2}$'},
                  'sphum': {'long_name': 'specific humidity', 'units': 'kg kg$^{-1}$'},
                  'precip': {'long_name': 'inst. precip. rate', 'units': 'mm h$^{-1}$'},
                  'h_anom': {'long_name': 'domainwise MSE anomaly\n', 'units': 'J kg$^{-1}$'},
                  'temp': {'long_name': 'temperature', 'units': 'K'},
                  'temp_anom': {'long_name': 'domainwise temperature anomaly\n', 'units': 'K'},
                  'sphum_anom': {'long_name': 'domainwise specific humidity anomaly\n', 'units': 'kg kg$^{-1}$'},
                  'wind_tangential': {'long_name': 'tangential velocity', 'units': 'm s$^{-1}$'},
                  'wind_radial': {'long_name': 'radial velocity', 'units': 'm s$^{-1}$'},
                  'tm': {'long_name': '300-500 hPa mean temperature\n', 'units': 'K'},
                  't_surf': {'long_name': 'surface temperature', 'units': 'K'},
                  'swfq': {'long_name': 'SWISHE frequency', 'units': '% of timestamps'},
                  'p-e': {'long_name': 'surface moisture flux', 'units': 'mm d$^{-1}$'}}
    
    if field in properties.keys():
        return properties[field]['long_name'], properties[field]['units']
    else:
        return field, '-'

def planar_composite(data, model_names, intensity_bins, field, pressure_level, 
                     experiment_plots=['control', 'swishe', 'swishe-control'], dpi=96, inline_statistics=False):
    """
    Method to obtain planar composite for prescribed data for a given field.
    Note: input must be in the form of a 3-tiered dictionary that is the output of tc_analyis.tc_model_data().

    Args:
        data (dict): 3-tiered dictionary
        model_names (list of str): list of strs with model names
        intensity_bins (list of str): list of strs with intensity bin listings
        field (str): name of field to evaluate
        pressure_level (numeric): level at which to evaluate field
        experiment_plot (str): experiment to plot (typically control, swishe, or swishe-control)
        dpi (int): dots per inch - 96 is default, 300 is recommended for publication-quality
        inline_staitstics (bool): boolean to control wheether or not statistics are published in the figures
    """
    
    # Define control and experiment names (assume both are in the density dictionary)
    run_control, run_experiment = ['control', 'swishe']
    # Define name for difference subdictionary
    run_difference = '{1}-{0}'.format(run_control, run_experiment)
    # Determine which experiments to evaluate
    experiments = [run_control, run_experiment, run_difference] if '-' in ''.join(experiment_plots) else experiment_plots
    
    ''' Generate composites. '''
    # Iterate over the models provided in the dictionary
    for model in model_names:
        # Iterate over the experiments provided in the dictionary
        for experiment in data[model].keys():
            print('\t Processing experiment: {0}'.format(experiment))
            # Create subdictionaries for prescribed intensity bins
            data[model][experiment]['composite'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Create dictionary to hold number of records found for each intensity bin per experiment per model
            data[model][experiment]['composite_records'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Iterate over intensity bins
            for intensity_bin in intensity_bins:
                # Pull composite and store
                key, N, composite_mean = composite.planar_compositor(model, data[model][experiment]['data'], intensity_bin, field, pressure_level=pressure_level)
                print('\t {0} records found for {1}, {2}, {3}'.format(len(N), model, experiment, intensity_bin))
                print(composite_mean.shape)
                # Populate the dictionaries for the given intensity bin
                data[model][experiment]['composite'][intensity_bin] = composite_mean
                data[model][experiment]['composite_records'][intensity_bin] = len(N)

    # Collect composites. Note: the differences will be calculated after this loop series.
    # Note: composite dictionary is structured top-down as: (1) intensity bin, (2) model, (3) experiment
    composites = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for intensity_bin in intensity_bins:
        # Initialize subdictionary for intensity bins
        composites[intensity_bin] = {}
        # Iterate over all models
        for model in model_names:
            # Initialize subdictionary for the iterand model
            composites[intensity_bin][model] = {}
            # Iterate over each given experiment
            for experiment in data[model].keys():
                composites[intensity_bin][model][experiment] = data[model][experiment]['composite'][intensity_bin]
                
    # Pass 2: Iterate through experiments to get the differences if a difference plot type (anything with a subtraction, or '-', is detected in the "experiment_plots" argument)
    if '-' in ''.join(experiment_plots):
        for intensity_bin in intensity_bins:
            for model in composites[intensity_bin].keys():
                try:
                    # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
                    delta = composites[intensity_bin][model][run_experiment].values - composites[intensity_bin][model][run_control].values
                    # Populate the dictionary with reconstructed xArray
                    composites[intensity_bin][model][run_difference] = xr.DataArray(data=delta, dims=['grid_yt', 'grid_xt'],
                                                                                coords={'grid_yt': (['grid_yt'], composites[intensity_bin][model][run_control]['grid_yt'].values), 
                                                                                        'grid_xt': (['grid_xt'], composites[intensity_bin][model][run_control]['grid_xt'].values)})
                except:
                    # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
                    # Methodology: trim to common centers of data. In other words, find the smaller composite and trim the bigger one to their dimensions.
                    # Note: I wrote this stoned, so it can likely be optimized
                    
                    # Get shape
                    shape_control, shape_exp = composites[intensity_bin][model][run_control].shape, composites[intensity_bin][model][run_experiment].shape
                    # Round decimal places for the coordinates (floating point precision errors lead to weird indexing)
                    composites[intensity_bin][model][run_control]['grid_xt'] = composites[intensity_bin][model][run_control]['grid_xt'].round(3)
                    composites[intensity_bin][model][run_control]['grid_yt'] = composites[intensity_bin][model][run_control]['grid_yt'].round(3)
                    composites[intensity_bin][model][run_experiment]['grid_xt'] = composites[intensity_bin][model][run_experiment]['grid_xt'].round(3)
                    composites[intensity_bin][model][run_experiment]['grid_yt'] = composites[intensity_bin][model][run_experiment]['grid_yt'].round(3)
                    
                    # Compare latitudes:
                    if shape_control[0] < shape_exp[0]:
                        temp = composites[intensity_bin][model][run_experiment].where((composites[intensity_bin][model][run_experiment]['grid_yt'] >= composites[intensity_bin][model][run_control]['grid_yt'].min())
                                                                                    & (composites[intensity_bin][model][run_experiment]['grid_yt'] <= composites[intensity_bin][model][run_control]['grid_yt'].max()), drop=True)
                        del composites[intensity_bin][model][run_experiment]
                        composites[intensity_bin][model][run_experiment] = temp
                    else:
                        temp = composites[intensity_bin][model][run_control].where((composites[intensity_bin][model][run_control]['grid_yt'] >= composites[intensity_bin][model][run_experiment]['grid_yt'].min())
                                                                                 & (composites[intensity_bin][model][run_control]['grid_yt'] <= composites[intensity_bin][model][run_experiment]['grid_yt'].max()), drop=True)
                        del composites[intensity_bin][model][run_control]
                        composites[intensity_bin][model][run_control] = temp
                    # now longitudes
                    if shape_control[1] < shape_exp[1]:
                        temp = composites[intensity_bin][model][run_experiment].where((composites[intensity_bin][model][run_experiment]['grid_xt'] >= composites[intensity_bin][model][run_control]['grid_xt'].min())
                                                                                    & (composites[intensity_bin][model][run_experiment]['grid_xt'] <= composites[intensity_bin][model][run_control]['grid_xt'].max()), drop=True)
                        del composites[intensity_bin][model][run_experiment]
                        composites[intensity_bin][model][run_experiment] = temp
                    else:
                        temp = composites[intensity_bin][model][run_control].where((composites[intensity_bin][model][run_control]['grid_xt'] >= composites[intensity_bin][model][run_experiment]['grid_xt'].min())
                                                                                 & (composites[intensity_bin][model][run_control]['grid_xt'] <= composites[intensity_bin][model][run_experiment]['grid_xt'].max()), drop=True)
                        del composites[intensity_bin][model][run_control]
                        composites[intensity_bin][model][run_control] = temp
                    # Crude fix written when I was stoned, but it worked
                    x_, y_ = composites[intensity_bin][model][run_control].values, composites[intensity_bin][model][run_experiment].values
                    if x_.shape == y_.shape:
                        delta = y_ - x_
                        grid_xt = composites[intensity_bin][model][run_control]['grid_xt'].values
                        grid_yt = composites[intensity_bin][model][run_control]['grid_yt'].values
                    else:
                        print(x_.shape, y_.shape)
                        if x_.shape[0] > y_.shape[0]:
                            # if x is larger
                            x_ = np.vstack((np.full((1, x_.shape[1]), np.nan), x_[:y_.shape[0], :], np.full((1, x_.shape[1]), np.nan)))
                            grid_yt = np.concatenate(([composites[intensity_bin][model][run_experiment].grid_yt.values[0]], composites[intensity_bin][model][run_control].grid_yt.values, [composites[intensity_bin][model][run_experiment].grid_yt.values[-1]]))
                        elif y_.shape[0] > x_.shape[0]: 
                            y_ = np.vstack((np.full((1, y_.shape[1]), np.nan), y_[:x_.shape[0], :], np.full((1, y_.shape[1]), np.nan)))
                            grid_yt = np.concatenate(([composites[intensity_bin][model][run_control].grid_yt.values[0]], composites[intensity_bin][model][run_experiment].grid_yt.values, [composites[intensity_bin][model][run_control].grid_yt.values[-1]]))
                        else:
                            grid_yt = composites[intensity_bin][model][run_control].grid_yt.values
                        print(x_.shape, y_.shape)
                        
                        if x_.shape[1] < y_.shape[1]:
                            # if x is larger
                            x_ = np.hstack((np.full((x_.shape[0], 1), np.nan), x_[:, :y_.shape[1]], np.full((x_.shape[0], 1), np.nan)))
                            grid_xt = np.concatenate(([composites[intensity_bin][model][run_experiment].grid_xt.values[0]], composites[intensity_bin][model][run_control].grid_xt.values, [composites[intensity_bin][model][run_experiment].grid_xt.values[-1]]))
                        elif x_.shape[1] > y_.shape[1]:
                            y_ = np.vstack((np.full((y_.shape[0], 1), np.nan), y_[:, :x_.shape[1]], np.full((y_.shape[0], 1), np.nan)))
                            grid_xt = np.concatenate(([composites[intensity_bin][model][run_control].grid_xt.values[0]], composites[intensity_bin][model][run_experiment].grid_xt.values, [composites[intensity_bin][model][run_control].grid_xt.values[-1]]))
                        else:
                            grid_xt = composites[intensity_bin][model][run_control].grid_xt.values
                        print(x_.shape, y_.shape)
                            
                        delta = y_ - x_
                    
                    # Populate the dictionary with reconstructed xArray
                    composites[intensity_bin][model][run_difference] = xr.DataArray(data=delta, dims=['grid_yt', 'grid_xt'],
                                                                                    coords={'grid_yt': (['grid_yt'], grid_yt), 
                                                                                            'grid_xt': (['grid_xt'], grid_xt)})
                
    ''' Get normalizations and colormaps. '''
    # Initialize normalization and colormaps. Normalizations will be intensity_bin specific. Save into dictionary for future use in colorbars.
    norms = {intensity_bin: {run_control: None, run_experiment: None, run_difference: None} for intensity_bin in intensity_bins}
    bounds = 24 # number of levels to bin for the normalization
    # Iterate over intensity bins to populate norms. Perform this over all prescribed experiments per intensity bin.
    for intensity_bin in intensity_bins:
        ''' Normalization #1: raw data. '''
        # Initialize minimum and maximum values
        vmin, vmax = np.nan, np.nan
        # Iterate over each experiment
        for experiment_plot in [run_control, run_experiment]:
            # Find extrema for this experiment
            vmin_, vmax_ = [min([np.nanmin(sv) for k, v in composites[intensity_bin].items() for sk, sv in v.items() if sk in [run_control, run_experiment]]),
                            max([np.nanmax(sv) for k, v in composites[intensity_bin].items() for sk, sv in v.items() if sk in [run_control, run_experiment]])]
            # If the extrema are larger than their corresponding extrema, update
            vmin, vmax = vmin_ if (vmin_ <= np.nanmin([vmin, vmin_])) else vmin, vmax_ if (vmax_ > np.nanmax([vmax, 0])) else vmax
        # Subdivide normalization into sequential (all one sign) or diverging (two signs, including 0) bins
        if vmin < 0 and vmax > 0:
            # Get larger of the two bounds
            extremum = max([abs(vmin), abs(vmax)])
            bins = np.concatenate([np.linspace(-extremum, 0, int(bounds/2))[:-1], [0], np.linspace(0, extremum, int(bounds/2))[1:]])
        else:
            bins = np.linspace(vmin, vmax, bounds+1)
        # Assign to normalization dictionary
        norms[intensity_bin][run_control], norms[intensity_bin][run_experiment] = matplotlib.colors.BoundaryNorm(bins, 256), matplotlib.colors.BoundaryNorm(bins, 256)
            
        ''' Normalization #2: difference plot, if requested. '''
        if run_difference in experiments:
            # Find extrema for this experiment
            vmin, vmax = [min([np.nanmin(sv) for k, v in composites[intensity_bin].items() for sk, sv in v.items() if sk == run_difference]),
                        max([np.nanmax(sv) for k, v in composites[intensity_bin].items() for sk, sv in v.items() if sk == run_difference])]
            # Subdivide normalization into sequential (all one sign) or diverging (two signs, including 0) bins
            if vmin < 0 and vmax > 0:
                # Get larger of the two bounds
                extremum = max([abs(vmin), abs(vmax)])
                bins = np.concatenate([np.linspace(-extremum, 0, int(bounds/2))[:-1], [0], np.linspace(0, extremum, int(bounds/2))[1:]])
            else:
                bins = np.linspace(vmin, vmax, bounds+1)
            # Assign to normalization dictionary
            norms[intensity_bin][run_difference] = matplotlib.colors.BoundaryNorm(bins, 256)
    # Initialize colormap for use in colorbars
    cmaps = {intensity_bin: {run_control: None, run_experiment: None, run_difference: None} for intensity_bin in intensity_bins}
    
    ''' Begin plotting. '''
    # Note: number of rows is dictated by number of models (+1 for rows), 
    #       number of columns is dictated by number of intensity bins TIMES the number of experiments desired to be plotted.
    nrows, ncols = len(model_names) + 1, len(experiment_plots)*len(intensity_bins)
    # Define height ratios as a function of the number of models
    height_ratios = [1 if row != nrows-1 else 0.1 for row in range(0, nrows)]
    # Initialize figure and grid. Note that nrows-1 is used to ignore the colorbar height.
    fig, grid = [plt.figure(figsize=(2*ncols, 2*(nrows-1)), constrained_layout=True, dpi=dpi), 
                 matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, height_ratios=height_ratios)]
    # Letter list
    # Iterate over intensity bins (columns)
    for intensity_bin_index, intensity_bin in enumerate(composites.keys()):
        # Iterate over models (rows)
        for model_index, model_name in enumerate(composites[intensity_bin].keys()):
            # Go through each desired experiment plot
            for experiment_index, experiment_plot in enumerate(experiment_plots):
                # Define colormap. To be modified when field-specific colormaps are accessed.
                # Note: this is in here because the colormap is normalization-specific due to data being used to define the colormap
                cmaps[intensity_bin][experiment_plot] = get_cmap(field, norms[intensity_bin][experiment_plot])
                # Plot data
                ax = fig.add_subplot(grid[model_index, intensity_bin_index + experiment_index])
                im = ax.pcolormesh(composites[intensity_bin][model_name][experiment_plot].grid_xt, composites[intensity_bin][model_name][experiment_plot].grid_yt, 
                                   composites[intensity_bin][model_name][experiment_plot], norm=norms[intensity_bin][experiment_plot], cmap=cmaps[intensity_bin][experiment_plot])
                # Subplot labeling - only label first rows
                # Change for publication to hide the composite mean or standard deviation
                # if model_index == 0:
                if '-' in experiment_plot:
                    in_plot_stat = ''
                else:
                    in_plot_stats = {'min': np.nanmin(composites[intensity_bin][model_name][experiment_plot]), 
                                     'mean': np.nanmean(composites[intensity_bin][model_name][experiment_plot]), 
                                     'max': np.nanmax(composites[intensity_bin][model_name][experiment_plot]), 
                                     'sum': np.nansum(composites[intensity_bin][model_name][experiment_plot])}
                    # Get the order of magnitude (oom) of the extremum, use the reciprocal to determine number of decimal points
                    oom = np.ceil(np.log10(1/in_plot_stats['mean']))
                    print('Order of magnitude: {0}'.format(oom))
                    if oom < -2 or oom > 1: # for small values
                        in_plot_stat = '{0:.2e};\n({1:.3e}, {2:.3e}); {3:.3e}'.format(in_plot_stats['mean'], in_plot_stats['min'], in_plot_stats['max'], in_plot_stats['sum']) 
                    else: # other
                        in_plot_stat = '{0:.2f};\n({1:.2f}, {2:.2f}); {3:.2f}'.format(in_plot_stats['mean'], in_plot_stats['min'], in_plot_stats['max'], in_plot_stats['sum'])
                # Modify label color depending on darkness of the field for raw fields only, not difference plots
                label_color = 'white' if ((field in ['olr', 'h', 'h_anom', 'wind']) and (experiment_plot != run_difference)) else 'k'
                ax.annotate('({0})'.format(chr(ord('a') + intensity_bin_index + experiment_index)), 
                            xy=(0.05, 0.95), xycoords='axes fraction', fontsize=9, color=label_color, alpha=1, **{'va': 'top'})
                if inline_statistics:
                    ax.annotate('{0}'.format(in_plot_stat), xy=(0.05, 0.05), xycoords='axes fraction', fontsize=6, color=label_color, alpha=1, **{'va': 'bottom'})
                
                # Tick formatting
                # Strategy: place minor ticks at 1-degree intervals, but make sure upper and lower limits are integers
                dx, dy = 1, 1
                x_min, x_max = composites[intensity_bin][model_name][experiment_plot].grid_xt.min(), composites[intensity_bin][model_name][experiment_plot].grid_xt.max()
                y_min, y_max = composites[intensity_bin][model_name][experiment_plot].grid_yt.min(), composites[intensity_bin][model_name][experiment_plot].grid_yt.max()
                gridline_x_minor, gridline_y_minor = np.arange(np.floor(x_min), np.floor(x_max), 1), np.arange(np.floor(y_min), np.floor(y_max), 1)
                ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_x_minor))
                ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_y_minor))
                
                # Labeling: column operations for intensity bin
                if (intensity_bin_index + experiment_index) == 0:
                    # Add label for the model on the first column
                    ax.set_ylabel(model_name, fontsize=10, labelpad=10)
                # Labeling: row operations for model name
                if model_index == 0:
                    # Hide tick labels as necessary (if not on last column or first row, hide them)
                    ax.xaxis.tick_top()
                    # Add label for the intensity bin on the first row
                    title_y = 1.25
                    # Only label intensities if there is more than one intensity bin
                    if len(intensity_bins) > 1:
                        ax.annotate('{0}'.format(intensity_bin), (0.5, title_y), 
                                    va='baseline', ha='center', xycoords='axes fraction', fontsize=10)
                        
                # Tick positioning: move xticks to top and hide if not the first row
                ax.xaxis.tick_top()
                ax.yaxis.tick_right()
                if model_index != 0:
                    ax.set_xticklabels([])
                # Tick positioning: move yticks to right and hide if not the last column
                if (intensity_bin_index + experiment_index) != ncols-1:
                    ax.set_yticklabels([])
    
    ''' Create colorbars.'''
    # Note that because normalizations are intensity bin specific, three colorbars will be made.
    # Colorbar ticklabel formatter
    # Limit colorbar label print to 1 print to avoid repetition
    label_index = False
    # Iterate over intensity bins and create colorbars
    for intensity_bin_index, intensity_bin in enumerate(composites.keys()):
        for experiment_index, experiment in enumerate(experiment_plots):
            # Get minimum and maximum for the intensity bin normalization
            vmin, vmax = norms[intensity_bin][experiment].vmin, norms[intensity_bin][experiment].vmax
            extremum = max([abs(vmin), abs(vmax)])
            # Get the order of magnitude (oom) of the extremum, use the reciprocal to determine number of decimal points
            oom = np.ceil(np.log10(1/extremum))
            print('Extremum: {0}, order of magnitude: {1}'.format(extremum, oom))
            if oom < 0: # for large values
                fmt = '%.1e'
            elif oom > 2: # for small values
                fmt = '%.1e'.format(int(oom))
            else: # other
                fmt = '%.{0}f'.format(int(oom))
            # Build the colorbar
            cw = fig.add_subplot(grid[-1, experiment_index]) # 'cw' stands for 'cax_wrapper'
            cw.set_axis_off()
            cax = cw.inset_axes([0, 0, 1, 1])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norms[intensity_bin][experiment], cmaps[intensity_bin][experiment]), 
                                    orientation='horizontal', cax=cax, format=matplotlib.ticker.ScalarFormatter())
            # Colorbar tick scientific notation formatting (see the ScalarFormatterForceFormat() class)
            fmt = ScalarFormatterForceFormat()
            fmt.set_powerlimits((-2, 3))
            cax.xaxis.set_major_formatter(fmt)
            # Modify font size of ticklabels
            cax.tick_params(labelsize=8)
            # Hide every other ticklabel
            [l.set_visible(False) for (i, l) in enumerate(cax.xaxis.get_ticklabels()) if i % 2 != 0]
            # Create colorbar label
            long_name, units = field_properties(field)
            vertical_level = ' at {0} hPa '.format(pressure_level) if pressure_level else ' '
            colorbar.set_label('[{1}]'.format(long_name, units, vertical_level), labelpad=10)
            
    # Universal plot labels
    fig.supxlabel('degrees from TC center', y=0.975, fontsize=10)   
    fig.supylabel('degrees from TC center', x=1, rotation=270, fontsize=10)     

    ''' Output statistics to obtain sample sizes for each composite. '''
    # Iterate over the models provided in the dictionary
    for model in model_names:
        # Iterate over the experiments provided in the dictionary
        for experiment in data[model].keys():
            [print('Number of records used for the composite for {0}, {1}, {2}: {3}'.format(model, experiment, intensity_bin, 
                                                                                            data[model][experiment]['composite_records'][intensity_bin])) for intensity_bin in intensity_bins]
       
    # return data
  
def azimuthal_composite_1d(data, model_names, field, dpi=96):
    
    # Initialize counter dictionary to hold data statistics for each model and corresponding experiment
    counts = {model_name: {'control': {}, 'swishe': {}} for model_name in model_names}
    
    # Multiplication factor for precipitation and evaporation
    factor = 3600 if field in ['precip', 'evap'] else 1
    
    # Initialize dictionary to hold processed data for each model and corresponding experiment
    individual_storms = {model_name: {'control': {}, 'swishe': {}} for model_name in model_names}
    # Iterate over each model and experiment to get the azimuthal means for a given field
    for model in individual_storms.keys():
        for experiment in individual_storms[model].keys():
            print(model, experiment, field)
            individual_storms[model][experiment] = composite.azimuthal_compositor(model, data[model][experiment]['data'], field)

    # Initialize dictionary to hold composite data
    outputs = {}
    # Outer radial index size
    index_ro = 20
    # Iterate over each model and experiment to populate the array that will be averaged to get the azimuthal composite mean
    for model in individual_storms.keys():
        outputs[model] = {}
        for experiment in individual_storms[model].keys():
            # Initialize an array with pre-defined dimensions characteristic of TCs to be composited
            # Dimensions are [0]: radius and [1]: number of storms
            print('Processing: ', model, experiment, field)
            outputs[model][experiment] = np.full((index_ro, len(individual_storms[model]['control'])), np.nan)
            # For each individual storm, ensure that the DataArray populates some fraction of the initialized nan array
            for i, value in enumerate(individual_storms[model][experiment]):
                # Trim data radially, if needed
                if value.values.shape[0] > index_ro:
                    arr = value.values[:index_ro]
                else:
                    arr = value.values
                outputs[model][experiment][0:arr.shape[0], i] = arr
            # Get number of entries that are being averaged
            counts[model][experiment] = outputs[model][experiment].shape[1]
            # Get mean and standard deviations with respect to axis 1: number of storms
            outputs[model][experiment] = {'mean': xr.DataArray(dims=['radius'], coords={'radius': (['radius'], np.arange(0, index_ro*0.9375, 0.9375))},
                                                               data=np.nanmean(outputs[model][experiment], axis=1)),
                                          'std': xr.DataArray(dims=['radius'], coords={'radius': (['radius'], np.arange(0, index_ro*0.9375, 0.9375))},
                                                              data=np.nanstd(outputs[model][experiment], axis=1))}
                
        # Get differences
        outputs[model]['swishe-control'] = {}
        for subfield in outputs[model]['control'].keys():
            counts[model]['swishe-control'] = counts[model]['control']
            outputs[model]['swishe-control'][subfield] = outputs[model]['swishe'][subfield] - outputs[model]['control'][subfield]

    # Outer radial index
    index_ro = 15

    ''' Filled field normalization. '''
    # Get extrema for each model run 

    ''' Plot. '''
    nrows, ncols = len(outputs.keys()), 2

    fig, gs = [plt.figure(figsize=(6, 2*nrows), constrained_layout=True, dpi=dpi), 
               matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, wspace=0.3, hspace=0.25)]
    linecolors = ['b', 'r']
    
    # Initialize axes dictionary
    axes = {}
    # Initialize extrema to be assigned during plotting
    vmin, vmax = np.nan, np.nan
    for model_index, model_name in enumerate(outputs.keys()):
        ax = fig.add_subplot(gs[model_index, 0])
        for experiment_index, experiment_name in enumerate(['control', 'swishe']):
            
            # Smooth the data with a quadratic interpolation 
            # Substep 1: create a basis vector with some multiple number of spaces relative to the input
            interp_basis = np.linspace(outputs[model_name][experiment_name]['mean'][:index_ro].radius.min(),
                                       outputs[model_name][experiment_name]['mean'][:index_ro].radius.max(),
                                       8*len(outputs[model_name][experiment_name]['mean'][:index_ro].radius))
            # Substep 2a: perform the interpolation. Perform rolling mean.
            interp_composite = outputs[model_name][experiment_name]['mean'][::-1].rolling(radius=4).mean()[::-1]
            # interp_composite = outputs[model_name][experiment_name]['mean']
            # Substep 2b: perform a quadratic interpolation to smooth the data
            interp_composite = interp_composite[:index_ro].interp(radius=interp_basis, method='slinear')
            # interp_composite = interp_composite.where(interp_composite >= 0, 0)
            # Plot data. Ensure the label only gets printed for the first plot to avoid repetition.
            label = experiment_name
            im = ax.plot(interp_basis, interp_composite*factor, lw=3, color=linecolors[experiment_index], label=label)
            # Modify y-axis to fit data maxima
            ax.set_xlim([0, 10])
            # Find extrema after interpolation + smoothing
            if np.isnan(vmax):
                vmax = np.nanmax(interp_composite*factor)
            elif np.nanmax(interp_composite*factor) > vmax:
                vmax = np.nanmax(interp_composite*factor)
                
            if np.isnan(vmin):
                vmin = np.nanmin(interp_composite*factor)
            elif np.nanmin(interp_composite*factor) < vmin:
                vmin = np.nanmin(interp_composite*factor)
            
            # ax.set_title('{0}, N = {1}'.format(experiment_name, counts[model_name][experiment_name]), fontsize=10)
            
    
            ax.legend(loc='best', frameon=False)
            ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, outputs[model_name][experiment_name]['mean'][:index_ro].radius.max(), 0.5)))

            # Labeling: column operations for intensity bin
            if experiment_index == 0:
                # Add label for the model on the first column
                ax.set_ylabel(model_name, fontsize=10, labelpad=10)
                
            if model_index < nrows-1:
                ax.set_xticklabels([])
                
            axes[model_index] = {'0': ax}
    
    print(vmin, vmax)
    # Iterate over axes to set y-limits
    for model_index, model_name in enumerate(outputs.keys()):
        for experiment_index, experiment_name in enumerate(['control', 'swishe']):
            axes[model_index]['0'].set_ylim([vmin, vmax])
            
    ''' Difference plotting. '''
             
    vmin, vmax = np.nan, np.nan   
    for model_index, model_name in enumerate(outputs.keys()):
        for experiment_index, experiment_name in enumerate(['swishe-control']):
            ax = fig.add_subplot(gs[model_index, -1])
            ax.axhline(0, c='k', lw=1, alpha=0.5)
            
            # Smooth the data with a quadratic interpolation 
            # Substep 1: create a basis vector with some multiple number of spaces relative to the input
            interp_basis = np.linspace(outputs[model_name][experiment_name]['mean'][:index_ro].radius.min(),
                                       outputs[model_name][experiment_name]['mean'][:index_ro].radius.max(),
                                       8*len(outputs[model_name][experiment_name]['mean'][:index_ro].radius))
            # Substep 2a: perform the interpolation. Perform rolling mean.
            interp_composite = outputs[model_name][experiment_name]['mean'][::-1].rolling(radius=4).mean()[::-1]
            # interp_composite = outputs[model_name][experiment_name]['mean']
            # Substep 2b: perform a quadratic interpolation to smooth the data
            interp_composite = interp_composite[:index_ro].interp(radius=interp_basis, method='slinear')
            
            im = ax.plot(interp_composite.radius, interp_composite*factor, lw=3, color='k')
            
            # Find extrema after interpolation + smoothing
            if np.isnan(vmax):
                vmax = np.nanmax(interp_composite*factor)
            elif np.nanmax(interp_composite*factor) > vmax:
                vmax = np.nanmax(interp_composite*factor)
                
            if np.isnan(vmin):
                vmin = np.nanmin(interp_composite*factor)
            elif np.nanmin(interp_composite*factor) < vmin:
                vmin = np.nanmin(interp_composite*factor)
            
            # ax.set_title('{0}, N = {1}'.format(experiment_name, counts[model_name][experiment_name]), fontsize=10)
            ax.set_xlim([0, 10])
            ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, outputs[model_name][experiment_name]['mean'][:index_ro].radius.max(), 0.5)))

            if model_index < nrows-1:
                ax.set_xticklabels([])
            axes[model_index] = {'1': ax}
            
    # Iterate over axes to set y-limits
    for model_index, model_name in enumerate(outputs.keys()):
        for experiment_index, experiment_name in enumerate(['swishe-control']):
            axes[model_index]['1'].set_ylim([vmin, vmax])
            
    # Universal plot labels
    fig.supxlabel('degrees from TC center', y=0, fontsize=10)
    
    long_name, units = field_properties(field)
    fig.supylabel('{0} [{1}]'.format(long_name, units), x=-0.05, fontsize=10)
    
def azimuthal_composite_2d(data, model_names, fields, dpi=96):
    
    # Initialize counter dictionary to hold data statistics for each model and corresponding experiment
    counts = {model_name: {'control': {}, 'swishe': {}} for model_name in model_names}
    
    # Initialize dictionary to hold processed data for each model and corresponding experiment
    individual_storms = {model_name: {'control': {}, 'swishe': {}} for model_name in model_names}
    # Iterate over each model and experiment to get the azimuthal means for a given field
    for model in individual_storms.keys():
        for experiment in individual_storms[model].keys():
            for field in fields:
                print(model, experiment, field)
                individual_storms[model][experiment][field] = composite.azimuthal_compositor(model, data[model][experiment]['data'], field)

                
    # Initialize dictionary to hold composite data
    outputs = {}
    # Outer radial index size
    index_ro = 20
    # Iterate over each model and experiment to populate the array that will be averaged to get the azimuthal composite mean
    for model in individual_storms.keys():
        outputs[model] = {}
        for experiment in individual_storms[model].keys():
            # Initialize an array with pre-defined dimensions characteristic of TCs to be composited
            # Dimensions are [0]: pfull (pressure level) and [1]: radius and [2]: number of storms
            outputs[model][experiment] = {}
            for field in fields:
                print('Processing: ', model, experiment, field)
                outputs[model][experiment][field] = np.full((32, index_ro, len(individual_storms[model]['control'][field])), np.nan)
                # For each individual storm, ensure that the DataArray populates some fraction of the initialized nan array
                for i, value in enumerate(individual_storms[model][experiment][field]):
                    # Ensure dimensions are properly oriented
                    v = value.transpose('pfull', 'radius')
                    # Trim data radially, if needed
                    if v.values.shape[1] > index_ro:
                        arr = v.values[:, :index_ro]
                    else:
                        arr = v.values
                    outputs[model][experiment][field][:, 0:arr.shape[1], i] = arr
                # Get number of entries that are being averaged
                counts[model][experiment] = outputs[model][experiment][field].shape[2]
                # Get mean and standard deviations with respect to axis 2: number of storms
                outputs[model][experiment][field] = {'mean': xr.DataArray(dims=['pfull', 'radius'],
                                                                          coords={'pfull': (['pfull'], v.pfull.values),
                                                                                  'radius': (['radius'], np.arange(0, index_ro*0.9375, 0.9375))},
                                                                          data=np.nanmean(outputs[model][experiment][field], axis=2)),
                                                     'std': xr.DataArray(dims=['pfull', 'radius'],
                                                                          coords={'pfull': (['pfull'], v.pfull.values),
                                                                                  'radius': (['radius'], np.arange(0, index_ro*0.9375, 0.9375))},
                                                                          data=np.nanstd(outputs[model][experiment][field], axis=2))}
                
                
        # Get differences
        outputs[model]['swishe-control'] = {}
        for field in fields:
            outputs[model]['swishe-control'][field] = {}
            for subfield in outputs[model]['control'][field].keys():
                print(field, subfield)
                counts[model]['swishe-control'] = counts[model]['control']
                outputs[model]['swishe-control'][field][subfield] = outputs[model]['swishe'][field][subfield] - outputs[model]['control'][field][subfield]
            
    # Filled contour field
    field_fill, field_open = fields

    # Define number of contour levels
    levels = 17
    # Outer radial index
    index_ro = 15
    # Initialize normalization dictionaries
    norm_fill, norm_open = {}, {}

    ''' Filled field normalization. '''
    # Get extrema for each model run 
    vmin, vmax = [np.nanmin([np.nanmin(outputs[model][experiment][field_fill]['mean'][:, :index_ro]) for model in outputs.keys() for experiment in outputs[model].keys()]), 
                  np.nanmax([np.nanmax(outputs[model][experiment][field_fill]['mean'][:, :index_ro]) for model in outputs.keys() for experiment in outputs[model].keys()])]
    extremum = max([abs(vmin), abs(vmax)])
    # Define normalization for each model run
    norm_fill['control'] = matplotlib.colors.BoundaryNorm(np.linspace(-extremum, extremum, levels), 256)
    norm_fill['swishe'] = matplotlib.colors.BoundaryNorm(np.linspace(-extremum, extremum, levels), 256)

    # Get extrema for each model run 
    vmin, vmax = [np.nanmin([np.nanmin(outputs[model]['swishe-control'][field_fill]['mean'][:, :index_ro]) for experiment in outputs[model].keys()]), 
                np.nanmax([np.nanmax(outputs[model]['swishe-control'][field_fill]['mean'][:, :index_ro]) for experiment in outputs[model].keys()])]
    extremum = max([abs(vmin), abs(vmax)])
    # Define normalization for the difference plot
    norm_fill['swishe-control'] = matplotlib.colors.BoundaryNorm(np.linspace(-extremum, extremum, levels), 256)

    # Initialize colormap dictionaries
    cmaps = {}
    # Initialize list of processed experiments for use later
    experiment_names = []

    ''' Open field normalization. '''
    # Define custom bounds for different fields for reasonable contour levels
    if field_open == 'wind_radial':
        step_size = 1
        field_open_levels = np.arange(-6, 6+step_size, step_size)
    elif field_open == 'omega':
        step_size = 0.25
        field_open_levels = np.arange(-1.5, 1.5+step_size, step_size)
    elif field_open == 'temp_anom':
        step_size = 0.5
        field_open_levels = np.arange(-4, 4+step_size, step_size)
    elif field_open == 'rh':
        step_size = 1
        field_open_levels = np.arange(-5, 5+step_size, step_size)
    else:
        field_open_levels = np.linspace(-extremum, extremum, 13)
    field_open_levels = field_open_levels[np.where(field_open_levels != 0.)]

    ''' Plot. '''
    nrows, ncols = len(outputs.keys()), len(outputs[list(outputs.keys())[0]].keys())
    height_ratios = [0.05] + [1]*nrows

    fig, gs = [plt.figure(figsize=(8, 3*nrows), constrained_layout=True, dpi=dpi), 
               matplotlib.gridspec.GridSpec(nrows=nrows+1, ncols=ncols, height_ratios=height_ratios, wspace=0.1, hspace=0.2)]

    for model_index, model_name in enumerate(outputs.keys()):
        cmaps = {experiment: get_cmap(field_fill, norm_fill[experiment]) for experiment in outputs[model_name].keys()}
        
        for experiment_index, experiment_name in enumerate(outputs[model_name].keys()):
            ax = fig.add_subplot(gs[model_index+1, experiment_index])
            im_fill = ax.contourf(outputs[model_name][experiment_name][field_fill]['mean'][:, :index_ro].radius,
                                outputs[model_name][experiment_name][field_fill]['mean'][:, :index_ro].pfull,
                                outputs[model_name][experiment_name][field_fill]['mean'][:, :index_ro], norm=norm_fill[experiment_name], levels=levels, cmap=cmaps[experiment_name])
            im_open = ax.contour(outputs[model_name][experiment_name][field_open]['mean'][:, :index_ro].radius,
                                outputs[model_name][experiment_name][field_open]['mean'][:, :index_ro].pfull,
                                outputs[model_name][experiment_name][field_open]['mean'][:, :index_ro], levels=field_open_levels, colors='k', alpha=0.5, linewidths=1)
            ax.clabel(im_open, im_open.levels[::2])
            
            
            # Labeling: column operations for intensity bin
            if experiment_index == 0:
                # Add label for the model on the first column
                ax.set_ylabel(model_name, fontsize=10, labelpad=10)
                
            # Subplot labeling - only label first rows
            if model_index == 0:
                ax.annotate('({0})'.format(chr(ord('a') + experiment_index)), (0.98, 0.97), xycoords='axes fraction', 
                            fontsize=10, color='k', alpha=1, **{'ha': 'right', 'va': 'top'})
                
            # ax.set_title('{0}, N = {1}'.format(experiment_name, counts[model_name][experiment_name]), fontsize=9)
            
            ax.set_ylim([10, 1000])
            ax.set_ylim(ax.get_ylim()[::-1])
            ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(np.arange(0, outputs[model_name][experiment_name][field_fill]['mean'][:, :index_ro].radius.max(), 1)))

            if model_index < nrows-1:
                ax.set_xticklabels([])
            if experiment_index > 0:
                ax.set_yticklabels([])
                
            experiment_names.append(experiment_name)
     
    # Get unique experiment names       
    experiment_names = sorted(list(set(experiment_names)))

    # Define counter for raw data (limit to one shared colorbar for non-derived datasets)
    counter = 0
    # Iterate over all experiments, print shared colorbar for raw data and one for the difference colorbar
    for experiment_index, experiment_name in enumerate(experiment_names):
        
        # Get minimum and maximum for the intensity bin normalization
        vmin, vmax = norm_fill[experiment_name].vmin, norm_fill[experiment_name].vmax
        extremum = max([abs(vmin), abs(vmax)])
        # Get the order of magnitude (oom) of the extremum, use the reciprocal to determine number of decimal points
        oom = np.ceil(np.log10(1/extremum))
        
        if oom < 0: # for large values
            fmt = '%.0f'
        elif oom > 2: # for small values
            fmt = '%.1e'.format(int(oom)+1)
        else: # other
            fmt = '%.{0}f'.format(int(oom)+1)
       
        if '-' in experiment_name:
            cax_wrapper = fig.add_subplot(gs[0, -1])
            cax_wrapper.set_axis_off()
        else:
            cax_wrapper = fig.add_subplot(gs[0, 0:-1])
            cax_wrapper.set_axis_off()
            
        if ('-' not in experiment_name and counter == 1) or '-' in experiment_name:
            cax = cax_wrapper.inset_axes([0.125, 1.5, 0.75, 1]) if '-' not in experiment_name else cax_wrapper.inset_axes([0, 1.5, 1, 1])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm_fill[experiment_name], cmaps[experiment_name]), 
                                    cax=cax, orientation='horizontal', format=matplotlib.ticker.ScalarFormatter())
            cax.xaxis.set_ticks_position('top')
            
            # Colorbar tick scientific notation formatting (see the ScalarFormatterForceFormat() class)
            fmt = ScalarFormatterForceFormat()
            fmt.set_powerlimits((-2, 3))
            cax.xaxis.set_major_formatter(fmt)
            # Modify font size of ticklabels
            cax.tick_params(labelsize=9)
            # Hide every other ticklabel
            [l.set_visible(False) for (i, l) in enumerate(cax.xaxis.get_ticklabels()) if i % 2 != 0]
            # Create colorbar label
            long_name, units = field_properties(field_fill)
            colorbar.set_label('{0} [{1}]'.format(long_name, units), labelpad=10)
            
            cax.xaxis.set_label_position('top')
            
        if '-' not in experiment_name:
            counter += 1
            
    # Universal plot labels
    fig.supxlabel('degrees from TC center', y=0.05, fontsize=10)
    fig.supylabel('Pressure [hPa]', x=-0.025, fontsize=10)     
            
def density_grid(data, model_names=['AM2.5', 'FLOR', 'HIRAM'], bin_resolution=5, dpi=96):
    """
    Method to plot the daily spatial density of all TC occurrences given a track data dictionary. 
    See tc_analysis.py --> tc_track_data() for more information.

    Args:
        data (dictionary): 3-tiered dictionary.
        model_name (list): list of strings with names of model (usually 'AM2.5', 'HIRAM', 'FLOR'.)
        bin_resolution (int, optional): resolution at which to generate spatial density, in degrees. Defaults to 5.
    """
    
    ''' Data processing. '''
    # Initialize dictionary for spatial density
    density = {}
    # Iterate over each model provided
    for model in data.keys():
        # Initialize model-specific subdictionary
        density[model] = {}
        # Iterate over each experiment provided
        for experiment in data[model].keys():
            # Define storage list for TC data
            dataset = [] 
            # Iterate through each unique TC and resample by day to get number of TC days per grid point
            for storm_id in data[model][experiment]['raw']['storm_id'].unique():
                dataset.append(data[model][experiment]['raw'].loc[data[model][experiment]['raw']['storm_id'] == storm_id].resample('D', on='time').first().reset_index())
            dataset = pd.concat(dataset)
            # Define dictionary to hold data relevant to track density heatmap
            counts = {'lon': [], 'lat': [], 'count': [], 'num_storms': []}
            # Group DataFrame into defined bins and add storm counts to the dictionary
            for lon_g, lon_v in dataset.groupby(pd.cut(dataset['center_lon'], np.arange(0, 360+bin_resolution, bin_resolution))):
                for lat_g, lat_v in lon_v.groupby(pd.cut(lon_v['center_lat'], np.arange(-90, 90+bin_resolution, bin_resolution))):
                    counts['lon'].append(lon_g.left)
                    counts['lat'].append(lat_g.left)
                    counts['count'].append(len(lat_v))
                    counts['num_storms'].append(len(dataset))
            # Create DataFrame from the dictionary
            counts = pd.DataFrame(counts)
            # Add time metadata to the DataFrame for bin per year estimate
            counts['time_min'], counts['time_max'] = dataset['time'].min(), dataset['time'].max()
            # Concatenate to get comprehensive DataFrame
            density[model][experiment] = counts

    # Define control and experiment names (assume both are in the density dictionary)
    run_control, run_experiment = ['control', 'swishe']
    # Get minimum and maximum years based on the control experiment
    year_min, year_max = density[model_names[0]][run_control]['time_min'].unique().year.item(), density[model_names[0]][run_control]['time_max'].unique().year.item()
    
    ''' Begin plotting. '''
    # Note: number of rows is dictated by number of models + 1 (for colorbar) 
    # Note: number of columns is dictated by number of experiments + 1 (for difference plot)
    nrows, ncols = len(model_names) + 1, len(density[model_names[0]].keys()) + 1
    # Pre-define height ratios based on number of rows
    height_ratios = [1 if i < nrows-1 else 0.05 for i in range(0, nrows)]
    # Initialize figure and grid
    fig, grid = [plt.figure(figsize=(15, 2.25*(nrows-1)), dpi=dpi, constrained_layout=True), 
                 matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, height_ratios=height_ratios, wspace=0.1, hspace=0)]
    # Define longitudinal offset
    longitude_offset = 180
    # Define projections (working and reference projections)
    proj, proj_ref = ccrs.PlateCarree(central_longitude=longitude_offset), ccrs.PlateCarree()
    
    # Collect processed density maps. Note: the difference will be calculated for the last subplot.
    densities = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for model in model_names:
        # Initialize subdictionary for the iterand model
        densities[model] = {}
        # Iterate over each given experiment
        for experiment in [run_control, run_experiment]:
            # Get longitude and latitude bins
            x, y = density[model][experiment].lon.unique(), density[model][experiment].lat.unique()
            # Get density array
            v = np.reshape(density[model][experiment]['count'].values, (len(x), len(y)))
            # Assign to dictionary for the given model/experiment configuration. 
            # Normalize by number of years and by day (assume 6-hourly data, so 6 x 4 = 1 day)
            densities[model][experiment] = v.T/(year_max - year_min)
        
    # Initialize normalization and colormaps. Save into dictionary for future use in colorbars.
    norms, cmaps = {}, {}
    bounds = 12 # number of levels to bin for the normalization
    # Maximum value across 'run_control' and 'run_experiment' for all models. Scale down to improve visibility of bins with lower counts.
    vmax_runs = 0.75*max([np.nanmax(sv) for k, v in densities.items() for sk, sv in v.items()]) 
    # Calculate density differences for the models of interest (given by 'model_names' argument)
    run_difference = '{1} - {0}'.format(run_control, run_experiment)
    # Log extrema across all model difference plots for future normalization
    difference_extrema = {model_name: None for model_name in model_names}
    # Iterate over models to calculate differences for each
    for model_name in densities.keys():
        densities[model_name][run_difference] = densities[model_name][run_experiment] - densities[model_name][run_control]
        difference_extrema[model_name] = max([abs(np.nanmin(densities[model_name][run_difference])), abs(np.nanmax(densities[model_name][run_difference]))])
    # Get extreme values for difference runs and find largest magnitude for normalization
    vmax_difference = max([np.nanmax(abs(v)) for k, v in difference_extrema.items()]) # maximum value across 'run_control' and 'run_experiment' for all models
    # Assign normalization and colormap values
    for experiment in densities[model_name].keys():
        if '-' in experiment: 
            norm_bins = [0.5*round(value/0.5) for  value in np.linspace(-vmax_difference, vmax_difference, bounds)]
            norms[experiment] = matplotlib.colors.BoundaryNorm(norm_bins, 256)
            cmaps[experiment] = 'bwr'
        else:
            norm_bins = [0.25*round(value/0.25) for value in np.linspace(0, vmax_runs, bounds+1)]
            norms[experiment] = matplotlib.colors.BoundaryNorm(norm_bins, 256)
            cmaps[experiment] = cmap_white_adjust('Reds')
    # Pass 2: Iterate through density dictionary to plot each model and each experiment's density data.
    for model_index, model_name in enumerate(densities.keys()):
        for experiment_index, experiment in enumerate(densities[model_name].keys()):
            
            try:
                # Use the basemap() method to obtain the base map figure upon which data is plotted
                fig, gs, ax, proj = basemap(fig, gs, model_name, experiment, 
                                            (densities[model_name][experiment]['time'].min().dt.year, densities[model_name][experiment]['time'].min()), row_num=model_index)
            
                print('Basemap connection successful.')
            except:
                print('Basemap connection unsuccessful.')
                # Initialize subplots
                ax = fig.add_subplot(grid[model_index, experiment_index], projection=proj)
                # Plot the data
                im = ax.pcolormesh(x, y, densities[model_name][experiment], norm=norms[experiment], cmap=cmaps[experiment], transform=proj_ref)
                # Define extent
                ax.set_extent([0, 359, -60, 60])
                
                ax.add_feature(cartopy.feature.LAND, color=(0.5, 0.5, 0.5, 0.25), zorder=99)
                ax.coastlines()
                
                # Subplot labeling
                title_y = 1.075
                # Set left-hand side to be {model name}, {min year} to {max year}
                subplot_title_model = ax.annotate('{0}, {1} to {2}'.format(model_name, year_min, year_max), 
                                                (0, title_y), va='baseline', ha='left', xycoords='axes fraction', fontsize=10)
                # Set right-hand side to be {experiment name}
                subplot_title_experiment = ax.annotate('{0}'.format(experiment), (1, title_y), va='baseline', ha='right', xycoords='axes fraction', fontsize=10)

                ''' Gridlines and ticks. Consider making this its own function. '''
                gridline_x_step, gridline_y_step = 60, 20
                gridline_minor_step = 10
                # Define major ticks
                gridline_x_major, gridline_y_major = [np.arange(0 - longitude_offset, 360 - longitude_offset + gridline_x_step, gridline_x_step), 
                                                    np.arange(-60, 60 + gridline_y_step, gridline_y_step)]
                gridline_x_minor, gridline_y_minor = [np.arange(0 - longitude_offset, 360 - longitude_offset + gridline_minor_step, gridline_minor_step), 
                                                    np.arange(-60, 60 + gridline_minor_step, gridline_minor_step)]
                # Draw gridlines
                gl = ax.gridlines(draw_labels=True, xlocs=[], ylocs=[], linewidth=0.5, color='k', alpha=0.25, linestyle='-')
                gl.top_labels, gl.right_labels = False, False
                # Define ticks
                ax.set_xticks(gridline_x_major, crs=proj)
                ax.set_yticks(gridline_y_major, [], crs=proj)
                ax.set_xticklabels([str(x + 180) for x in gridline_x_major]) # override tick values that are projection-dependent
                # Set minor ticks
                ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_x_minor))
                ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_y_minor))
            
                ''' End replacement here. '''
            # Tick positioning: move xticks to top and hide if not the first row
            if model_index != nrows-2:
                ax.set_xticklabels([])
            # Tick positioning: move yticks to right and hide if not the last column
            if experiment_index != 0:
                ax.set_yticklabels([])
    
    ''' Create 2 colorbars: 1 for the experiments, another for the difference. '''
    # Colorbar ticklabel formatter
    fmt = lambda x, pos: '{:.2f}'.format(x)
    
    # Experiments
    cw_experiments = fig.add_subplot(grid[-1, 0:2]) # 'cw' stands for 'cax_wrapper'
    cw_experiments.set_axis_off()
    cax_experiments = cw_experiments.inset_axes([0.25, -3, 0.5, 1])
    colorbar_experiments = fig.colorbar(matplotlib.cm.ScalarMappable(norms[run_control], cmaps[run_control]), 
                                        orientation='horizontal', cax=cax_experiments, format=matplotlib.ticker.FuncFormatter(fmt))
    colorbar_experiments.set_label('TC days per year per {0}$\degree$ bin'.format(bin_resolution), labelpad=10)
    # Difference
    cw_difference = fig.add_subplot(grid[-1, -1]) # 'cw' stands for 'cax_wrapper'
    cw_difference.set_axis_off()
    cax_difference = cw_difference.inset_axes([0, -3, 1, 1])
    colorbar_difference = fig.colorbar(matplotlib.cm.ScalarMappable(norms[run_difference], cmaps[run_difference]), 
                                       orientation='horizontal', cax=cax_difference, format=matplotlib.ticker.FuncFormatter(fmt))
    colorbar_difference.set_label('$\\Delta$(TC days per year per {0}$\degree$ bin)'.format(bin_resolution), labelpad=10)

def pdf(data, models=['AM2.5', 'FLOR', 'HIRAM'], param='center_lat', num_bins=60, mode='pdf'):
    
    fig, ax = plt.subplots(figsize=(5, 3), dpi=300)
    
    # Define list of colors (one color per model) and line styles (one linestyle per experiment)
    colors = ['b', 'g', 'm']
    linestyles = ['-', '--']
    
    # Collect data for output counts
    counts = {model: {} for model in models}

    # Obtain maximum values from histogram distributions
    vmax = 0
    # Iterate over all models and experiments
    for model_index, model in enumerate(models):
        for experiment_index, experiment in enumerate(data[model].keys()):
            # Get unique values for the parameter passed
            out = data[model][experiment]['unique'][param].dropna()
            # Filter out nans and infs
            out = out.loc[np.isfinite(out.values)]
            # Get bin values and bin edges, don't plot
            n, bins, _ = ax.hist(out, bins=num_bins, histtype=u'step', density=True, lw=0) 
            # Incorporate mutiplicative factor - if mode == 'full_count', factor = length of data; else, factor = 1. The former gives you a full count.
            factor = len(out) if mode == 'full_count' else 1
            # Plot the distribution
            ax.plot(bins[:-1], n*factor, lw=2.5, color=colors[model_index], linestyle=linestyles[experiment_index], 
                    label='{0}, {1}'.format(model, experiment))
            # Check to see if number of entries for a given bin exceeds the maximum (vmax)
            vmax = np.nanmax(n*factor) if np.nanmax(n*factor) > vmax else vmax
            # Add data length to counts
            counts[model][experiment] = {'bins': bins, 'counts': n*len(out), 'median': np.median(out), 'std': np.nanstd(out)}
            
    # Adjust axis limits: handle limit differently for center_lat to allow for inline legend
    ax.set_ylim([0, 1.1*vmax])
    if param == 'center_lat':
        # Curb latitudinal extent of TC monitoring to midlatitude extrema
        ax.set_xlim([-60, 60])
    elif param == 'duration':
        # Storms lasting > 30 d are unrealistic
        ax.set_xlim([0, 30])
    
    # Case-by-case dictionary for improved plot labeling (log the long names)
    param_names = {'center_lat': 'Latitude', 'max_wind': 'Maximum 10 m wind speed [m s$^{-1}$]', 'min_slp': 'Minimum sea level pressure [hPa]'}
    # Plot labeling and formatting
    ax.set_xlabel(param_names[param], labelpad=10)
    ax.set_ylabel('Probability density', labelpad=10)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
    # Handle legend plotting differently for center_lat to allow for space inline
    if param == 'center_lat':
        ax.legend(frameon=False, loc='upper left', prop={'size': 8})
    else:
        ax.legend(frameon=False, loc='best', prop={'size': 8})
        
    # Print median and standard deviations for all model-experiment combinations
    for model in counts.keys():
        for experiment in counts[model].keys():
            print('{0} - {1}, {2}: median: {3:.2f}, std.dev: {4:.2f}'.format(param, model, experiment,
                                                                             counts[model][experiment]['median'], counts[model][experiment]['std']))
        
    return counts

def swishe_frequency(models=['HIRAM', 'AM2.5', 'FLOR'], dpi=144, set_visible=True):
    """Plot frequency of SWISHE filter application using the 'swfq' diagnostic on 'atmos_4xdaily' files for SWISHE runs.

    Args:
        models (list, optional): list of models to visualize data for. Defaults to ['FLOR', 'HIRAM'].
    """
    # Initialize figure and GridSpec object with the number of rows equal to the number of models plotted
    fig, gs = plt.figure(figsize=(6, 2*len(models)), dpi=dpi), matplotlib.gridspec.GridSpec(nrows=len(models), ncols=1, hspace=0.75)
    # Define colormap and normalization
    cmap, norm = cmap_white_adjust('Reds'), matplotlib.colors.BoundaryNorm(np.linspace(0, 2, 11), 256)
    # Filename dictionary
    filenames = {'HIRAM': 'model_HIRAM-exp_CTL1990s_swishe-type_atmos-var_swfq-mean_year-resample-101_150.nc',
                 'AM2.5': 'model_AM2.5-exp_CTL1990s_swishe-type_atmos-var_swfq-mean_year-resample-full-101_150.nc',
                 'FLOR': 'model_FLOR-exp_CTL1990s_swishe-type_atmos-var_swfq-mean_year-resample-full-2050_2100.nc'}
    # Iterate over each model, using 'row' as the indexer
    for row, model in enumerate(models):
        # Get directory name for data
        dirname = '/tigress/GEOCLIM/gr7610/analysis/model_out'
        # Get filenames for data
        filename = os.path.join(dirname, filenames[model])
        # Access the data, targeting 'swfq' to get SWISHE application frequency
        data = xr.open_dataset(filename)['swfq']
        # Use the basemap() method to obtain the base map figure upon which data is plotted
        fig, gs, ax, proj = basemap(fig, gs, model, 'SWISHE application frequency', (data.time.min().dt.year.item(), data.time.max().dt.year.item()), row_num=row, land=True)
        # Plot contoured data at 10 levels over a 2% range (intervals of 0.2% for legibility)
        im = ax.contourf(data.grid_xt, data.grid_yt, 100*data.sum(dim='time')/len(data.time.values), 
                         levels=10, cmap=cmap, norm=norm, vmin=0, vmax=2, transform=proj)
        # If row isn't the last one, hide the xlabel since axes are implicitly shared and equal
        if row != gs.nrows-1:
            ax.set_xlabel('')
        # Plot the colorbar
        cax = ax.inset_axes([1.03, 0, 0.03, 1])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
        
        long_name, units = field_properties('swfq')
        colorbar.set_label('{0}'.format(units), rotation=270, labelpad=20)
        
    if not set_visible:
        plt.clf()
        
    return data

def basemap(fig, gs, model_name, experiment, year_range=None, row_num=0, col_num=0, extent=[0, 359, -60, 60], land=False, xlabel=True, ylabel=True):
    """
    Function to provide template map upon which geospatial data can be plotted.

    Args:
        fig (matplotlib Figure object): Figure object.
        gs (matplotlib GridSpec object): GridSpec object.
        model_name (str): name of model, used for annotation.
        experiment (str): name of experiment, used for annotation.
        year_range (list or tuple): 2-item list with minimum and maximum years.
        row_num (int, optional): row number for use in a GridSpec object. Defaults to 0.
        col_num (int, optional): column number for use in a GridSpec object. Defaults to 0.
        extent (list, optional): 4-item geospatial extent list, with parameters of [minimum longitude, maximum longitude, minimum latitude, maximum latitude]. 
                                 Defaults to [0, 359, -60, 60].
        land (bool, optional): determine whether or not to add the land feature. Useful if land mask is used for data. Defaults to False.

    Returns:
        fig (matplotlib Figure object): Figure object.
        gs (matplotlib GridSpec object): GridSpec object.
        ax (matplotlib Subplot object): Subplot object.
        proj (Cartopy projection object): Cartopy projection object.
    """
    
    # Define longitudinal offset
    longitude_offset = 180
    # Define projections (working and reference projections)
    proj, proj_ref = ccrs.PlateCarree(central_longitude=longitude_offset), ccrs.PlateCarree()
    
    # Initialize subplot
    ax = fig.add_subplot(gs[row_num, col_num], projection=proj)
    # Define extent
    ax.set_extent(extent, crs=proj_ref)
    
    if land:
        ax.add_feature(cartopy.feature.LAND, color=(0.5, 0.5, 0.5, 0.25), zorder=99)
    ax.coastlines()
    
    # Subplot labeling
    title_y = 1.05
    # Define year range string plotting
    year_range_str = ', {0} to {1}'.format(min(year_range), max(year_range)) if year_range else ''
    # Set left-hand side to be {model name}, {min year} to {max year}
    subplot_title_model = ax.annotate('{0}{1}'.format(model_name, year_range_str),
                                      (0, title_y), va='baseline', ha='left', xycoords='axes fraction', fontsize=10)
    # Set right-hand side to be {experiment name}
    subplot_title_experiment = ax.annotate('{0}'.format(experiment), (1, title_y), va='baseline', ha='right', xycoords='axes fraction', fontsize=10)

    # Initialize gridline parameters
    gridline_x_step, gridline_y_step = 60, 20
    gridline_minor_step = 10
    # Define major ticks
    gridline_x_major, gridline_y_major = [np.arange(extent[0] - longitude_offset, extent[1] - longitude_offset + gridline_x_step, gridline_x_step), 
                                          np.arange(extent[2], extent[3] + gridline_y_step, gridline_y_step)]
    gridline_x_minor, gridline_y_minor = [np.arange(extent[0] - longitude_offset, extent[1] - longitude_offset + gridline_minor_step, gridline_minor_step), 
                                          np.arange(extent[2], extent[3] + gridline_minor_step, gridline_minor_step)]
    
    # Draw gridlines
    gl = ax.gridlines(draw_labels=True, xlocs=[], ylocs=[], linewidth=0.5, color='k', alpha=0.25, linestyle='-')
    gl.top_labels, gl.right_labels = False, False
    # Define ticks
    ax.set_xticks(gridline_x_major, crs=proj)
    ax.set_yticks(gridline_y_major, [], crs=proj)
    ax.set_xticklabels([str(x + 180) for x in gridline_x_major]) # override tick values that are projection-dependent
    # Set minor ticks
    ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_x_minor))
    ax.yaxis.set_minor_locator(matplotlib.ticker.FixedLocator(gridline_y_minor))
    # Axis labels
    if xlabel:
        ax.set_xlabel('Longitude')
    if ylabel:
        ax.set_ylabel('Latitude')
    ax.set_extent(extent, crs=proj_ref)

    return fig, gs, ax, proj_ref

def basins(visualize=False):
    
    basins = {'NI': [40, 100, 0, 30],
              'SI': [40, 120, -40, 0],
              'WP': [100, 180, 0, 40],
              'SP': [120, 260, -40, 0],
              'EP': [180, 260, 0, 40],
              'NA': [260, 360, 0, 40],
              'SA': [300, 360, -40, 0],
              'global': [0, 360, -60, 60]}

    if visualize:
        proj, proj_ref = ccrs.PlateCarree(central_longitude=180), ccrs.PlateCarree()
        fig, gs = plt.figure(figsize=(8, 3)), matplotlib.gridspec.GridSpec(nrows=1, ncols=1)

        fig, gs, ax, proj_ref = basemap(fig, gs, 'Basins', '', row_num=0, col_num=0, extent=[0, 359, -60, 60])
        
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        
        ax.coastlines()
        for i, basin_data in enumerate(basins.items()):
            basin, basin_extent = basin_data
            patch = matplotlib.patches.Rectangle(xy=(basin_extent[0], basin_extent[2]),
                                                 width=(basin_extent[1] - basin_extent[0]), height=(basin_extent[3] - basin_extent[2]), 
                                                 fc='none', lw=2, ec=colors[i], transform=proj_ref, zorder=i+10)
            ax.annotate(text=basin, xy=(basin_extent[0], basin_extent[2]), xytext=(5, 5), 
                        textcoords='offset points', transform=proj_ref)
            ax.add_patch(patch)
        ax.set_extent([0, 359, -60, 60])
        plt.show()

    return basins

def swishe_TC_overlay(model_name, years=None, bin_resolution=2.5, dpi=144):
    
    start_year, end_year = min(years), max(years)
    track_data = tc_analysis.tc_track_data([model_name], ['control', 'swishe'], year_range=(start_year, end_year), storm_type='TS')
    
    ''' Data processing. '''
    # Initialize dictionary for spatial density
    density = {}
    # Iterate over each model provided
    for model in track_data.keys():
        # Initialize model-specific subdictionary
        density[model] = {}
        # Iterate over each experiment provided
        for experiment in track_data[model].keys():
            # Define storage list for TC data
            dataset = [] 
            # Iterate through each unique TC and resample by day to get number of TC days per grid point
            for storm_id in track_data[model][experiment]['raw']['storm_id'].unique():
                dataset.append(track_data[model][experiment]['raw'].loc[track_data[model][experiment]['raw']['storm_id'] == storm_id].resample('D', on='time').first().reset_index())
            dataset = pd.concat(dataset)
            # Define dictionary to hold data relevant to track density heatmap
            counts = {'lon': [], 'lat': [], 'count': [], 'num_storms': []}
            # Group DataFrame into defined bins and add storm counts to the dictionary
            for lon_g, lon_v in dataset.groupby(pd.cut(dataset['center_lon'], np.arange(0, 360+bin_resolution, bin_resolution))):
                for lat_g, lat_v in lon_v.groupby(pd.cut(lon_v['center_lat'], np.arange(-90, 90+bin_resolution, bin_resolution))):
                    counts['lon'].append(lon_g.left)
                    counts['lat'].append(lat_g.left)
                    counts['count'].append(len(lat_v))
                    counts['num_storms'].append(len(dataset))
            # Create DataFrame from the dictionary
            counts = pd.DataFrame(counts)
            # Add time metadata to the DataFrame for bin per year estimate
            counts['time_min'], counts['time_max'] = dataset['time'].min(), dataset['time'].max()
            # Concatenate to get comprehensive DataFrame
            density[model][experiment] = counts
            
    # Collect processed density maps. Note: the difference will be calculated for the last subplot.
    densities = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for model in [model_name]:
        # Initialize subdictionary for the iterand model
        densities[model] = {}
        # Iterate over each given experiment
        for experiment in ['control', 'swishe']:
            # Get longitude and latitude bins
            x, y = density[model][experiment].lon.unique(), density[model][experiment].lat.unique()
            # Get density array
            v = np.reshape(density[model][experiment]['count'].values, (len(x), len(y)))
            # Assign to dictionary for the given model/experiment configuration. 
            # Normalize by number of years and by day (assume 6-hourly data, so 6 x 4 = 1 day)
            densities[model][experiment] = v.T/(end_year - start_year)
            
    fig, gs = plt.figure(figsize=(8, 3), dpi=dpi), matplotlib.gridspec.GridSpec(nrows=1, ncols=1)
    fig, gs, ax, proj_ref = basemap(fig, gs, model_name, 'Track density & SWISHE application', year_range=(start_year, end_year), land=True)

    contours = np.arange(0, 1.6 + 0.2, 0.2)
    norm = matplotlib.colors.BoundaryNorm(contours, 256)
    cmap = cmap_white_adjust('Blues', levels=len(contours)-1)

    im = ax.contourf(x, y, densities[model_name]['control'], transform=proj_ref, norm=norm, cmap=cmap, levels=len(contours)-1)

    swishe_contours = swishe_frequency(models=[model_name], dpi=144, set_visible=False)
    im_swishe = ax.contour(swishe_contours.grid_xt, swishe_contours.grid_yt, 100*swishe_contours.sum(dim='time')/len(swishe_contours.time.values), 
                           levels=contours[1:], cmap='Reds', norm=norm, vmin=0, vmax=max(contours), linewidths=1.5, transform=proj_ref)
    print(im_swishe.levels)

    cax = ax.inset_axes([1.05, 0, 0.03, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
    colorbar.set_label('% of time with TC', labelpad=15, rotation=270)