import cartopy, cartopy.crs as ccrs
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import matplotlib, matplotlib.pyplot as plt

import composite, utilities

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

def get_cmap(field, norm):
    """
    Method to pull colormap for a given field from a .csv with colormap info.

    Args:
        field (str): field of interest
        norm (matplotlib.colors.Normalize): normalization with data to extract values from

    Returns:
        colormap: matplotlib.colors.colormap
    """
    
    print('Colormap field: {0}'.format(field))
    # Get colormap DataFrame from .csv
    colormaps = pd.read_csv('/projects/GEOCLIM/gr7610/reference/param_colormaps.csv')
    # Get field-specific colormap. If field is not in colormaps, default to 'other'.
    if field in colormaps['param'].unique():
        colormap = colormaps.loc[colormaps['param'] == field]
    else:
        colormap = colormaps.loc[colormaps['param'] == 'other']
    
    return colormap['difference'].item() if (norm.vmin < 0 and norm.vmax > 0) else colormap['normal'].item()

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
                  'wind': {'long_name': '10m horizontal wind speed', 'units': 'm s$^{-1}$'},
                  'WVP': {'long_name': 'column-integrated water vapor', 'units': 'kg m$^{-2}$'},
                  'h': {'long_name': 'moist static energy', 'units': 'J kg$^{-1}$'},
                  'slp': {'long_name': 'sea level pressure', 'units': 'hPa'},}
    
    if field in properties.keys():
        return properties[field]['long_name'], properties[field]['units']
    else:
        return field, '-'

def planar_composite(data, model_names, intensity_bins, field, pressure_level, experiment_plots=['control']):
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
            # Create subdictionaries for prescribed intensity bins
            data[model][experiment]['composite'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Create dictionary to hold number of records found for each intensity bin per experiment per model
            data[model][experiment]['composite_records'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Iterate over intensity bins
            for intensity_bin in intensity_bins:
                # Pull composite and store
                key, N, composite_mean = composite.planar_compositor(model, data[model][experiment]['data'], intensity_bin, field, pressure_level=pressure_level)
                print('\t {0} records found for {1}, {2}, {3}'.format(len(N), model, experiment, intensity_bin))
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
                # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
                delta = composites[intensity_bin][model][run_experiment].values - composites[intensity_bin][model][run_control].values
                # Populate the dictionary with reconstructed xArray
                composites[intensity_bin][model][run_difference] = xr.DataArray(data=delta, dims=['grid_yt', 'grid_xt'],
                                                                            coords={'grid_yt': (['grid_yt'], composites[intensity_bin][model][run_control]['grid_yt'].values), 
                                                                                    'grid_xt': (['grid_xt'], composites[intensity_bin][model][run_control]['grid_xt'].values)})
            
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
    fig, grid = [plt.figure(figsize=(2*ncols, 2*(nrows-1)), constrained_layout=True), 
                 matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, height_ratios=height_ratios)]
    
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
            if oom < 0: # for large values
                fmt = '%.0f'
            elif oom > 2: # for small values
                fmt = '%.1e'.format(int(oom))
            else: # other
                fmt = '%.{0}f'.format(int(oom))
            # Build the colorbar
            cw = fig.add_subplot(grid[-1, experiment_index]) # 'cw' stands for 'cax_wrapper'
            cw.set_axis_off()
            cax = cw.inset_axes([0, 0, 1, 1])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norms[intensity_bin][experiment], cmaps[intensity_bin][experiment]), 
                                    orientation='horizontal', cax=cax, format=matplotlib.ticker.FormatStrFormatter(fmt))
            # Modify font size of ticklabels
            cax.tick_params(labelsize=8)
            # Hide every other ticklabel
            [l.set_visible(False) for (i, l) in enumerate(cax.xaxis.get_ticklabels()) if i % 2 != 0]
            # Create colorbar label
            long_name, units = field_properties(field)
            vertical_level = ' at {0} hPa '.format(pressure_level) if pressure_level else ' '
            colorbar.set_label('[{1}]'.format(long_name, units, vertical_level), labelpad=10)
            
    ''' Output statistics to obtain sample sizes for each composite. '''
    # Iterate over the models provided in the dictionary
    for model in model_names:
        # Iterate over the experiments provided in the dictionary
        for experiment in data[model].keys():
            [print('Number of records used for the composite for {0}, {1}, {2}: {3}'.format(model, experiment, intensity_bin, 
                                                                                            data[model][experiment]['composite_records'][intensity_bin])) for intensity_bin in intensity_bins]
            
def azimuthal_composite(data, model_names, intensity_bins, field, experiment_plots=['control']):
    """
    Method to obtain azimuthal composite for prescribed data for a given field.
    Note: input must be in the form of a 3-tiered dictionary that is the output of tc_analyis.tc_model_data().

    Args:
        data (dict): 3-tiered dictionary
        model_names (list of str): list of strs with model names
        intensity_bins (list of str): list of strs with intensity bin listings
        field (str): name of field to evaluate
        experiment_plots (list of str): experiments to plot (typically control, swishe, or swishe-control)
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
            # Create subdictionaries for prescribed intensity bins
            data[model][experiment]['composite'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Create dictionary to hold number of records found for each intensity bin per experiment per model
            data[model][experiment]['composite_records'] = {intensity_bin: None for intensity_bin in intensity_bins}
            # Iterate over intensity bins
            for intensity_bin in intensity_bins:
                # Pull composite and store
                key, N, composite_mean = composite.azimuthal_compositor(model, data[model][experiment]['data'], intensity_bin, field)
                # Populate the dictionary for the given intensity bin
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
    # Pass 2: Iterate through experiments to get the differences.
    if '-' in ''.join(experiment_plots):
        for intensity_bin in intensity_bins:
            for model in composites[intensity_bin].keys():
                try:
                    # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
                    delta = composites[intensity_bin][model][run_experiment].values - composites[intensity_bin][model][run_control].values
                    # Populate the dictionary with reconstructed xArray. Handle 1D and 2D data differently.
                    if 'pfull' in composites[intensity_bin][model][run_experiment].dims:
                        composites[intensity_bin][model][run_difference] = xr.DataArray(data=delta, dims=['pfull', 'radius'],
                                                                                        coords={'pfull': (['pfull'], composites[intensity_bin][model][run_control]['pfull'].values), 
                                                                                                'radius': (['radius'], composites[intensity_bin][model][run_control]['radius'].values)})
                    else:
                        composites[intensity_bin][model][run_difference] = xr.DataArray(data=delta, dims=['radius'], coords={'radius': (['radius'], composites[intensity_bin][model][run_control]['radius'].values)})
                except:
                    continue
            
    ''' Get normalizations and colormaps. '''
    # Initialize normalization and colormaps. Normalizations will be intensity_bin specific. Save into dictionary for future use in colorbars.
    norms = {intensity_bin: {run_control: None, run_experiment: None, run_difference: None} for intensity_bin in intensity_bins}
    bounds = 16 # number of levels to bin for the normalization
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
    fig, grid = [plt.figure(figsize=(2*ncols, 2*(nrows-1)), constrained_layout=True), 
                 matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, height_ratios=height_ratios)]
    
    # Iterate over intensity bins (columns)
    for intensity_bin_index, intensity_bin in enumerate(composites.keys()):
        # Iterate over models (rows)
        for model_index, model_name in enumerate(composites[intensity_bin].keys()):
            # Go through each desired experiment plot
            for experiment_index, experiment_plot in enumerate(experiment_plots):
                # Initialize subplot
                ax = fig.add_subplot(grid[model_index, intensity_bin_index + experiment_index])
                # A split occurs here: if data is 1D (azimuthally-composited planar values), process one ways. Else, process another.
                if 'pfull' in composites[intensity_bin][model_name][experiment_plot].dims:
                    # Define colormap. To be modified when field-specific colormaps are accessed.
                    # Note: this is in here because the colormap is normalization-specific due to data being used to define the colormap
                    cmaps[intensity_bin][experiment_plot] = get_cmap(field, norms[intensity_bin][experiment_plot])
                    # Plot data.
                    im = ax.contourf(composites[intensity_bin][model_name][experiment_plot].radius, composites[intensity_bin][model_name][experiment_plot].pfull, 
                                    composites[intensity_bin][model_name][experiment_plot], norm=norms[intensity_bin][experiment_plot], cmap=cmaps[intensity_bin][experiment_plot])
                    # Invert y-axis
                    ax.set_ylim(ax.get_ylim()[::-1])
                    # Tick positioning: move yticks to right and hide if not the last column
                    if intensity_bin_index != ncols-1:
                        ax.set_yticklabels([])
                    else:
                        ax.yaxis.tick_right()
                else:
                    # Create horizontal line to denote 0
                    ax.axhline([0], lw=0.5, c='k')
                    # Smooth the data with a quadratic interpolation 
                    # Substep 1: create a basis vector with some multiple number of spaces relative to the input
                    interp_basis = np.linspace(composites[intensity_bin][model_name][experiment_plot].radius.min(),
                                            composites[intensity_bin][model_name][experiment_plot].radius.max(),
                                            4*len(composites[intensity_bin][model_name][experiment_plot].radius))
                    # Substep 2: perform the interpolation
                    interp_composite = composites[intensity_bin][model_name][experiment_plot].interp(radius=interp_basis)
                    # Plot data.
                    im = ax.plot(interp_composite.radius, interp_composite, lw=2)
                    # Modify y-axis to fit data maxima
                    ax.set_ylim([norms[intensity_bin][experiment_plot].vmin, norms[intensity_bin][experiment_plot].vmax])
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
                if model_index == 0:
                    ax.xaxis.tick_top()
                else:
                    if model_index != 0:
                        ax.set_xticklabels([])
                    # Tick positioning: move yticks to right and hide if not the last column
                    if (intensity_bin_index + experiment_index) != ncols-1:
                        ax.set_yticklabels([])
    
    ''' Create colorbars.'''
    # Note that because normalizations are intensity bin specific, three colorbars will be made.
    # Colorbar ticklabel formatter
    # Iterate over intensity bins and create colorbars - only if pressure coordinates in data
    if 'pfull' in composites[intensity_bin][model_name][run_control].dims:
        for intensity_bin_index, intensity_bin in enumerate(composites.keys()):
            for experiment_index, experiment in enumerate(experiment_plots):
                # Get minimum and maximum for the intensity bin normalization
                vmin, vmax = norms[intensity_bin][experiment].vmin, norms[intensity_bin][experiment].vmax
                extremum = max([abs(vmin), abs(vmax)])
                # Get the order of magnitude (oom) of the extremum, use the reciprocal to determine number of decimal points
                oom = np.ceil(np.log10(1/extremum))
                if oom < 0: # for large values
                    fmt = '%.0f'
                elif oom > 2: # for small values
                    fmt = '%.1e'.format(int(oom))
                else: # other
                    fmt = '%.{0}f'.format(int(oom))
                # Build the colorbar
                cw = fig.add_subplot(grid[-1, experiment_index]) # 'cw' stands for 'cax_wrapper'
                cw.set_axis_off()
                cax = cw.inset_axes([0, 0, 1, 1])
                colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norms[intensity_bin][experiment], cmaps[intensity_bin][experiment]), 
                                        orientation='horizontal', cax=cax, format=matplotlib.ticker.FormatStrFormatter(fmt))
                # Modify font size of ticklabels
                cax.tick_params(labelsize=8)
                # Hide every other ticklabel
                [l.set_visible(False) for (i, l) in enumerate(cax.xaxis.get_ticklabels()) if i % 2 != 0]
                # Create colorbar label
                long_name, units = field_properties(field)
                colorbar.set_label('[{1}]'.format(long_name, units), labelpad=10)
            
    ''' Output statistics to obtain sample sizes for each composite. '''
    # Iterate over the models provided in the dictionary
    for model in model_names:
        # Iterate over the experiments provided in the dictionary
        for experiment in data[model].keys():
            [print('Number of records used for the composite for {0}, {1}, {2}: {3}'.format(model, experiment, intensity_bin, 
                                                                                            data[model][experiment]['composite_records'][intensity_bin])) for intensity_bin in intensity_bins]
            
def density_grid(data, model_names=['AM2.5', 'FLOR', 'HIRAM'], bin_resolution=5):
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
    fig, grid = [plt.figure(figsize=(15, 2.25*(nrows-1)), dpi=300, constrained_layout=True), 
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
            # Initialize subplot
            ax = fig.add_subplot(grid[model_index, experiment_index], projection=proj)
            # Plot the data
            im = ax.pcolormesh(x, y, densities[model_name][experiment], norm=norms[experiment], cmap=cmaps[experiment], transform=proj_ref)
            # Define extent
            ax.set_extent([0, 359, -60, 60])
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
    # Differnece
    cw_difference = fig.add_subplot(grid[-1, -1]) # 'cw' stands for 'cax_wrapper'
    cw_difference.set_axis_off()
    cax_difference = cw_difference.inset_axes([0, -3, 1, 1])
    colorbar_difference = fig.colorbar(matplotlib.cm.ScalarMappable(norms[run_difference], cmaps[run_difference]), 
                                       orientation='horizontal', cax=cax_difference, format=matplotlib.ticker.FuncFormatter(fmt))
    colorbar_difference.set_label('$\\Delta$(TC days per year per {0}$\degree$ bin)'.format(bin_resolution), labelpad=10)

def pdf(data, param='center_lat', num_bins=60):
    
    fig, ax = plt.subplots(figsize=(4, 3))
    
    # Define list of colors (one color per model) and line styles (one linestyle per experiment)
    colors = ['b', 'g', 'm']
    linestyles = ['-', '--']
    
    # Obtain maximum values from histogram distributions
    vmax = 0
    # Iterate over all models and experiments
    for model_index, model in enumerate(data.keys()):
        for experiment_index, experiment in enumerate(data[model].keys()):
            # Get unique values for the parameter passed
            out = data[model][experiment]['unique'][param].dropna()
            # Filter out nans and infs
            out = out.loc[np.isfinite(out.values)]
            # Get bin values and bin edges, don't plot
            n, bins, _ = ax.hist(out, bins=num_bins, histtype=u'step', density=True, lw=0) 
            # Plot the distribution
            ax.plot(bins[:-1], n, lw=2.5, color=colors[model_index], linestyle=linestyles[experiment_index], 
                    label='{0}, {1}'.format(model, experiment))
            # Check to see if number of entries for a given bin exceeds the maximum (vmax)
            vmax = np.nanmax(n) if np.nanmax(n) > vmax else vmax
            
    # Adjust axis limits: handle limit differently for center_lat to allow for inline legend
    ax.set_ylim([0, 1.1*vmax])
    if param == 'center_lat':
        # Curb latitudinal extent of TC monitoring to midlatitude extrema
        ax.set_xlim([-60, 60])
    elif param == 'duration':
        # Storms lasting > 30 d are unrealistic
        ax.set_xlim([0, 30])
    
    # Case-by-case dictionary for improved plot labeling (log the long names)
    param_names = {'center_lat': 'Latitude', 'max_wnd': 'Maximum 10 m wind speed [m s$^{-1}$]', 'min_slp': 'Minimum sea level pressure [hPa]'}
    # Plot labeling and formatting
    ax.set_xlabel(param_names[param], labelpad=10)
    ax.set_ylabel('Probability density', labelpad=10)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.2f'))
    # Handle legend plotting differently for center_lat to allow for space inline
    if param == 'center_lat':
        ax.legend(frameon=False, loc='upper left', prop={'size': 8})
    else:
        ax.legend(frameon=False, loc='best', prop={'size': 8})