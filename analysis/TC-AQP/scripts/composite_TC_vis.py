import matplotlib, matplotlib.pyplot as plt
import numpy as np
import sys
import xarray as xr

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import visualization

def multiplot_TC_GCM_composite(configuration_data: dict,
                               field_name: str,
                               contour_levels: int=12,
                               maximum_latitude: int=30):

    ''' 
    Method to plot composite means of TC and GCM fields, 
    as well as the composite mean anomalies of TCs from climatology.

    Only valid for a single model and multiple model-specific experiments.

    Note: `configuration_data` is provided in the following format:
    {MODEL_NAME-EXPERIMENT_NAME: {'TC': DATA_TC,
                                  'GCM': DATA_GCM}}
    '''

    ''' Define figure structure parameters. '''
    nrows = len(configuration_data.keys()) # each configuration has its own row
    ncols = 3 # composite TC mean, composite GCM mean, composite mean anomaly

    ax_width = 2 # subplot size (square)
    wspace, hspace = 0.125, 0.125 # horizontal and vertical padding for subplots
    width_ratios = [0.25, 1, 1, 1, 0.25] # one colorbar on each side, three central columns
    height_ratios = [1] * nrows # all rows of equal height
    ncols = ncols + 2 # add colorbars to column number

    fig_width = ax_width * np.sum(width_ratios) + (wspace * ncols * ax_width) # final figure width
    fig_height = ax_width * np.sum(height_ratios) + (hspace * nrows * ax_width) # final figure height

    fig, gs = [plt.figure(figsize=(fig_width, fig_height)),
               matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, 
                                            wspace=wspace, hspace=hspace, 
                                            width_ratios=width_ratios, height_ratios=height_ratios)]

    for config_index, config_name in enumerate(configuration_data.keys()):
        # Prep inputs for the modular TC-GCM compositing function
        entries_TC = configuration_data[config_name]['TC']
        entries_GCM = configuration_data[config_name]['GCM']
        model_name, experiment_name = config_name.split(':')
        # Plot the configuration composite sequence
        plot_TC_GCM_composite(entries_TC=entries_TC,
                              entries_GCM=entries_GCM,
                              model_name=model_name,
                              experiment_name=experiment_name,
                              field_name=field_name,
                              contour_levels=contour_levels,
                              maximum_latitude=maximum_latitude,
                              fig=fig,
                              gs=gs,
                              ax_index=config_index)

def plot_TC_GCM_composite(entries_TC: xr.Dataset,
                          entries_GCM: xr.Dataset,
                          model_name: str,
                          experiment_name: str,
                          field_name: str,
                          contour_levels: int=12,
                          maximum_latitude: int=30,
                          fig: matplotlib.figure.Figure|None=None,
                          gs: matplotlib.gridspec.GridSpec|None=None,
                          ax_index: int|None=None):

    # Filter entries by latitude
    TMP_TC, TMP_GCM = [entries_TC.where(abs(entries_TC['center_lat']) < maximum_latitude).dropna(dim='storm_id'),
                       entries_GCM.where(abs(entries_GCM['center_lat']) < maximum_latitude).dropna(dim='storm_id')]

    # Build configuration dictionary data structure
    configs = {'TC': TMP_TC[field_name],
               'GCM': TMP_GCM[field_name],
               'TC—GCM': TMP_TC[field_name] - TMP_GCM[field_name]}

    # Initialize the figure; if `fig` and `gs` are None, assume that a standalone figure is being requested
    if fig is None and gs is None:
        # Figure structure parameters
        ax_width = 2
        wspace = 0.125
        width_ratios = [0.2, 1, 1, 1, 0.2]
        ncols = len(configs) + 2
        fig_width = ax_width * np.sum(width_ratios) + (wspace * ncols * ax_width)
        # Initialize the figure
        fig, gs = [plt.figure(figsize=(fig_width, ax_width), dpi=144), 
                   matplotlib.gridspec.GridSpec(nrows=1, ncols=ncols, wspace=wspace, width_ratios=width_ratios)]

    # Initial pass: plot composite means for all configurations
    config_means = {} # store TC-averaged values for normalization processing
    config_extrema = {} # store extrema for experiment data
    for index, config_key in enumerate(configs.keys()):
        # Assign a temporary alias for the configuration value
        value = configs[config_key]
        # Get a TC-averaged value (i.e., TC composite)
        config_means[config_key] = value.mean(['storm_id', 'grid_yt'])
        # Store extrema
        config_extrema[config_key] = {'min': config_means[config_key].min().item(), 
                                      'max': config_means[config_key].max().item(), 
                                      'count': value['storm_id'].count().item()}
    
    # Define normalization parameters
    experiment_norm, experiment_cmap = visualization.norm_cmap([], 
                                                               field=field_name, 
                                                               extrema=(min([config_extrema['TC']['min'], config_extrema['GCM']['min']]),
                                                                        max([config_extrema['TC']['max'], config_extrema['GCM']['max']])),
                                                               num_bounds=contour_levels) # for TC and GCM
    difference_norm, difference_cmap = visualization.norm_cmap([], 
                                                               field=field_name, 
                                                               extrema=(config_extrema['TC—GCM']['min'], config_extrema['TC—GCM']['max']),
                                                               num_bounds=contour_levels) # for the TC-GCM difference
    
    # Second pass: plot data
    ims = {} # store plot data
    for col_index, config_key in enumerate(configs.keys()):
        row_index = 0 if ax_index is None else ax_index
        # Initialize the iterand configuration axis
        ax = fig.add_subplot(gs[row_index, col_index + 1])
    
        # Define normalization and colormap based on configuration type
        norm = difference_norm if config_key == 'TC—GCM' else experiment_norm
        cmap = difference_cmap if config_key == 'TC—GCM' else experiment_cmap
        
        # Plot the composite mean
        ims[config_key] = ax.pcolormesh(config_means[config_key].grid_xt_TC, config_means[config_key].grid_yt_TC, config_means[config_key], norm=norm, cmap=cmap)
    
        # Plot metadata
        if row_index == 0: ax.set_title(config_key, fontsize=11)
        if col_index > 0: ax.set_yticklabels([]) 
        if fig is not None and (row_index != (gs.nrows-1)): ax.set_xticklabels([]) 
        ax.set_aspect('equal')
    
        # Plot annotation with statistics
        statistics = f'N = {config_extrema[config_key]['count']}\n({config_extrema[config_key]['min']:.1f}, {config_extrema[config_key]['max']:.1f})'
        annotation = ax.annotate(statistics, xy=(0.03, 0.03), xycoords='axes fraction', 
                                 c='k', fontsize=8, ha='left', va='bottom')
        annotation.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground=(1, 1, 1, 0.75)),
                                     matplotlib.patheffects.Normal()])
    
    # Third pass: plot colorbars
    for index, config_key in enumerate(configs.keys()):
        if index == 1: continue # don't repeat the first colorbar
        # Define colorbar axis parameters based on configuration type
        cax_column = -1 if config_key == 'TC—GCM' else 0
        # Define normalization and colormap based on configuration type
        norm = difference_norm if config_key == 'TC—GCM' else experiment_norm
        cmap = difference_cmap if config_key == 'TC—GCM' else experiment_cmap

        # Initialize the colorbar axis
        colorbar_ax = fig.add_subplot(gs[row_index, cax_column]) 
        colorbar_ax.set_axis_off()
        cax = colorbar_ax.inset_axes([0, 0, 0.25, 1])
        # Plot the colorbar
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                cax=cax)
        # Format the tick formatting
        number_of_cax_ticks = (len(norm.boundaries) - 1) // 4 + 1
        cax.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(number_of_cax_ticks))
        # Apply labeling
        long_name, units = visualization.field_properties(field_name)
        long_name = f'$\delta$ (TC — GCM)' if config_key == 'TC—GCM' else long_name

        if config_key in ['TC', 'GCM']:
            cax.yaxis.set_ticks_position('left')
            colorbar.set_label(f'{long_name} [{units}]', labelpad=-60)
        else:
            cax.yaxis.set_ticks_position('right')
            colorbar.set_label(f'{long_name} [{units}]', rotation=270, labelpad=20)


    if fig is None:
        fig.tight_layout()
