import matplotlib, matplotlib.pyplot as plt
import numpy as np
import sys
import xarray as xr

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import visualization

def plot_TC_GCM_composite(entries_TC: xr.Dataset,
                          entries_GCM: xr.Dataset,
                          model_name: str,
                          experiment_name: str,
                          field_name: str,
                          contour_levels: int=12,
                          maximum_latitude: int=30,
                          fig: matplotlib.figure.Figure|None=None,
                          gs: matplotlib.gridspec.GridSpec|None=None):

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
        ax_width = 2.5
        height_ratios = [1, 0.05]
        ncols = len(configs)
        wspace, hspace = 0.125, 0.25
        # Initialize the figure
        fig, gs = [plt.figure(figsize=(ax_width * ncols, sum(ax_width * np.array(height_ratios) * (1 + hspace))), dpi=144), 
                   matplotlib.gridspec.GridSpec(nrows=2, ncols=ncols, height_ratios=height_ratios)]

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
    for index, config_key in enumerate(configs.keys()):
        ax = fig.add_subplot(gs[0, index])
    
        # Define normalization and colormap based on configuration type
        norm = difference_norm if config_key == 'TC—GCM' else experiment_norm
        cmap = difference_cmap if config_key == 'TC—GCM' else experiment_cmap
        
        # Plot the composite mean
        ims[config_key] = ax.pcolormesh(config_means[config_key].grid_xt_TC, config_means[config_key].grid_yt_TC, config_means[config_key], norm=norm, cmap=cmap)
    
        # Plot metadata
        ax.set_title(config_key, fontsize=10)
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
        ax_column_span = -1 if config_key == 'TC—GCM' else slice(0, 2)
        cax_column_width = 1 if config_key == 'TC—GCM' else 0.5
        # Define normalization and colormap based on configuration type
        norm = difference_norm if config_key == 'TC—GCM' else experiment_norm
        cmap = difference_cmap if config_key == 'TC—GCM' else experiment_cmap

        # Initialize the colorbar axis
        ax = fig.add_subplot(gs[1, ax_column_span]) 
        ax.set_axis_off()
        cax = ax.inset_axes([(1 - cax_column_width)/2, 0, cax_column_width, 1])
        # Plot the colorbar
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), 
                                cax=cax,
                                orientation='horizontal')
        # Format the tick formatting
        number_of_cax_ticks = (len(norm.boundaries) - 1) // 4 + 1
        cax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(number_of_cax_ticks))
        # Apply labeling
        long_name, units = visualization.field_properties(field_name)
        long_name = f'$\delta$ (TC — GCM)' if config_key == 'TC—GCM' else long_name
        colorbar.set_label(f'{long_name} [{units}]')

    if fig is None:
        fig.tight_layout()
