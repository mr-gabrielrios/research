import matplotlib, matplotlib.pyplot as plt
import numpy as np
import os
import sys
import xarray as xr

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import visualization

def get_composite_latitude_distribution(configuration_data: dict,
                                        maximum_latitude: int=90,
                                        savefig: bool=False):

    ''' Method to get meridional distribution of TC snapshots used for compositing. '''

    # Initialize the figure
    dpi = 300 if savefig else 144
    dist_fig, dist_ax = plt.subplots(figsize=(6, 2.5), dpi=dpi)

    for config_index, config_name in enumerate(configuration_data.keys()):
        # Prep inputs for the modular TC-GCM compositing function
        entries_TC = configuration_data[config_name]['TC']
        entries_GCM = configuration_data[config_name]['GCM']
        # Get latitudes of the configuration snapshots
        snapshot_latitudes = entries_TC['center_lat'].where(abs(entries_TC['center_lat']) <= maximum_latitude, np.nan).values

        dist_ax.hist(snapshot_latitudes, bins=np.arange(-90, 90, 2.5), histtype='step', label=config_name)
    
    dist_ax.set_xlabel('Latitude')
    dist_ax.set_ylabel('Count')
    dist_ax.legend(frameon=False)

    # Save the figure, is the option is chosen
    if savefig:
        # Define storage parameters
        storage_model_name = [config_name.split(':')[0] for config_name in configuration_data.keys()][0]
        storage_experiment_names = ':'.join([config_name.split(':')[1].split('.')[1] for config_name in configuration_data.keys()])
        storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-AQP/figs'
        storage_filename = f'TC-GCM.meridional_distribution.{storage_model_name}-{storage_experiment_names}.pdf'
        # Save the figure
        plt.savefig(os.path.join(storage_dirname, storage_filename), dpi=dpi, bbox_inches='tight', format='pdf')

    
def multiplot_TC_GCM_composite(configuration_data: dict,
                               field_name: str,
                               contour_levels: int=12,
                               maximum_latitude: int=30,
                               plotting_method: str='pcolormesh',
                               savefig: bool=False):

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
    dpi = 300 if savefig else 144

    fig, gs = [plt.figure(figsize=(fig_width, fig_height), dpi=dpi),
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
                              plotting_method=plotting_method,
                              fig=fig,
                              gs=gs,
                              ax_index=config_index)

    # Add figure labels
    long_name, units = visualization.field_properties(field_name)
    fig.supxlabel('degrees from TC center', y=0.05)
    fig.supylabel(f'{long_name} [{units}]', fontsize=11) # left-hand y-label
    fig.text(0.95, 0.5, f'$\delta$ (TC — GCM) [{units}]', ha='center', va='center', rotation=270, fontsize=10) # right-hand y-label

    # Save the figure, is the option is chosen
    if savefig:
        # Define storage parameters
        storage_model_name = [config_name.split(':')[0] for config_name in configuration_data.keys()][0]
        storage_experiment_names = ':'.join([config_name.split(':')[1].split('.')[1] for config_name in configuration_data.keys()])
        storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-AQP/figs'
        storage_filename = f'TC-GCM.configuration.{storage_model_name}-{storage_experiment_names}.field_name.{field_name}.pdf'
        # Save the figure
        plt.savefig(os.path.join(storage_dirname, storage_filename), dpi=dpi, bbox_inches='tight', format='pdf')

    # Plot the meridional distribution of snapshots used for compositing
    get_composite_latitude_distribution(configuration_data=configuration_data,
                                        maximum_latitude=maximum_latitude,
                                        savefig=savefig)


def plot_TC_GCM_composite(entries_TC: xr.Dataset,
                          entries_GCM: xr.Dataset,
                          model_name: str,
                          experiment_name: str,
                          field_name: str,
                          contour_levels: int=12,
                          maximum_latitude: int=30,
                          plotting_method: str='pcolormesh',
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
    
    # Define normalization parameters for raw experiment values
    experiment_norm, experiment_cmap = visualization.norm_cmap([], 
                                                               field=field_name, 
                                                               extrema=(min([config_extrema['TC']['min'], config_extrema['GCM']['min']]),
                                                                        max([config_extrema['TC']['max'], config_extrema['GCM']['max']])),
                                                               num_bounds=contour_levels) # for TC and GCM
    
    # Define normalization parameters for difference experiment
    colormap_override = False # lazy flag to use to prevent confusing sequential colormaps from being used
    if field_name in ['netrad_toa']: 
        # Override difference extrema for net radiation to ensure centered extrema about zero
        extremum = max([abs(config_extrema['TC—GCM']['min']), abs(config_extrema['TC—GCM']['max'])]) 
        extrema = (-extremum, extremum)
        colormap_override = True # ensure a diverging colormap for a field with typically-divergent values
    else:
        extrema = (config_extrema['TC—GCM']['min'], config_extrema['TC—GCM']['max'])

    difference_norm, difference_cmap = visualization.norm_cmap([], 
                                                               field=field_name, 
                                                               extrema=extrema,
                                                               num_bounds=contour_levels) # for the TC-GCM difference
    
    # Second pass: plot data
    ims = {} # store plot data
    for col_index, config_key in enumerate(configs.keys()):
        row_index = 0 if ax_index is None else ax_index
        # Initialize the iterand configuration axis
        ax = fig.add_subplot(gs[row_index, col_index + 1])
    
        # Define normalization and colormap based on configuration type
        norm = difference_norm if (config_key == 'TC—GCM' or colormap_override) else experiment_norm
        cmap = difference_cmap if (config_key == 'TC—GCM' or colormap_override) else experiment_cmap
        
        # Plot the composite mean
        if plotting_method == 'contour':
            ims[config_key] = ax.contourf(config_means[config_key].grid_xt_TC, config_means[config_key].grid_yt_TC, config_means[config_key], norm=norm, cmap=cmap, levels=contour_levels)
        else:
            ims[config_key] = ax.pcolormesh(config_means[config_key].grid_xt_TC, config_means[config_key].grid_yt_TC, config_means[config_key], norm=norm, cmap=cmap, linewidth=0, rasterized=True)

    
        # Plot metadata
        if row_index == 0: ax.set_title(f'{model_name}, {config_key}', fontsize=10)
        if col_index > 0: ax.set_yticklabels([]) 
        if fig is not None and (row_index != (gs.nrows-1)): ax.set_xticklabels([]) 
        ax.set_aspect('equal')
    
        # Plot subfigure letter
        subplot_letter_factor = 1 if model_name == 'HIRAM' else 0
        subplot_figure_index = row_index * (gs.ncols - 2) + col_index + gs.nrows * (gs.ncols - 2) * subplot_letter_factor # calculate iterand subfigure index, subtract 2 for colorbars
        
        subplot_figure_letter = visualization.get_alphabet_letter(subplot_figure_index) # get subfigure letter
        letter_annotation = ax.annotate(f'({subplot_figure_letter}) {experiment_name.split('.')[1]}', xy=(0.03, 0.96), xycoords='axes fraction', 
                                        ha='left', va='top', fontsize=9)
        letter_annotation.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground=(1, 1, 1, 0.75)),
                                            matplotlib.patheffects.Normal()])

        # Plot annotation with statistics
        statistics = f'N = {config_extrema[config_key]['count']}\n({config_extrema[config_key]['min']:.1f}, {config_extrema[config_key]['max']:.1f})'
        statistics_annotation = ax.annotate(statistics, xy=(0.03, 0.03), xycoords='axes fraction',
                                            c='k', fontsize=8, ha='left', va='bottom')
        statistics_annotation.set_path_effects([matplotlib.patheffects.Stroke(linewidth=2, foreground=(1, 1, 1, 0.75)),
                                                matplotlib.patheffects.Normal()])
    
    # Third pass: plot colorbars
    for index, config_key in enumerate(configs.keys()):
        if index == 1: continue # don't repeat the first colorbar
        # Define colorbar axis parameters based on configuration type
        cax_column = -1 if config_key == 'TC—GCM' else 0
        # Define normalization and colormap based on configuration type
        norm = difference_norm if (config_key == 'TC—GCM' or colormap_override) else experiment_norm
        cmap = difference_cmap if (config_key == 'TC—GCM' or colormap_override) else experiment_cmap

        # Initialize the colorbar axis
        colorbar_ax = fig.add_subplot(gs[row_index, cax_column]) 
        colorbar_ax.set_axis_off()
        cax = colorbar_ax.inset_axes([0, 0, 0.25, 1])
        # Plot the colorbar
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                cax=cax, extend='both')
        # Format the tick formatting
        number_of_cax_ticks = (len(norm.boundaries) - 1) // 4 + 1
        cax.yaxis.set_major_locator(matplotlib.ticker.LinearLocator(number_of_cax_ticks))
        # Apply labeling
        long_name, units = visualization.field_properties(field_name)
        long_name = f'$\delta$ (TC — GCM)' if config_key == 'TC—GCM' else long_name

        if config_key in ['TC', 'GCM']:
            cax.yaxis.set_ticks_position('left')
            if ax_index is None: colorbar.set_label(f'{long_name} [{units}]', labelpad=-60)
        else:
            cax.yaxis.set_ticks_position('right')
            if ax_index is None: colorbar.set_label(f'{long_name} [{units}]', rotation=270, labelpad=20)


    if fig is None:
        fig.tight_layout()
