''' Import packages. '''
# Time packages
import calendar, cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import functools, importlib, os, pickle, collections, sys

import pandas as pd, xarray as xr, nc_time_axis
xr.set_options(keep_attrs=True)
# Visualization tools
import cartopy, cartopy.util as cutil, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

import analysis

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import derived, tc_analysis, utilities, visualization, track_TCs
importlib.reload(utilities)
importlib.reload(analysis)
importlib.reload(visualization)
importlib.reload(tc_analysis)
importlib.reload(derived)
importlib.reload(track_TCs)

def get_experiment_alias(experiment_name: str) -> str:
    
    ''' Helper function to assign shorthand names for long-form experiment names. '''
    
    experiment_aliases = {'CTL1990_C192_CONST': 'CTL, const.',
                          'CTL1990_swishe_C192_CONST': 'SWISHE, const.',
                          'CTL1990_C192_FLAT': 'CTL, flat',
                          'CTL1990_swishe_C192_FLAT': 'SWISHE, flat',
                          'CTL1990_C192_FLAT75': 'CTL, flat + warm',
                          'CTL1990_swishe_C192_FLAT75': 'SWISHE, flat + warm',
                          'CTL1990_C192_AQP15': 'CTL, 15° N',
                          'CTL1990_swishe_C192_AQP15': 'SWISHE, 15° N',
                          'CTL1990s_tiger3': 'CTL, AGCM',
                          'CTL1990s_swishe_tiger3': 'SWISHE, AGCM',
                          'CTL1990_FA_tiger3': 'CTL, AOGCM',
                          'CTL1990_swishe_FA_tiger3': 'SWISHE, AOGCM',
                          'CTL1990.SLAB_40': 'HIRAM-CTL',
                          'CTL1990_SWISHE.SLAB_40': 'HIRAM-SWISHE'}
    
    # Catch the case in which the experiment name is not in the dictionary
    if experiment_name in experiment_aliases.keys():
        return experiment_aliases[experiment_name]
    else:
        return experiment_name

def get_field_name_alias(field_name: str) -> str:
    
    field_name_aliases = {'olr': 'OLR',
                          'netrad_toa': '$\mathrm{R_{TOA}}$',
                          'evap': 'E',
                          'precip': 'P',
                          'WVP': 'WVP',
                          'rh': 'RH',
                          'rh200': 'RH, 200 hPa',
                          'rh500': 'RH, 500 hPa',
                          'cld_amt': 'CC, 200 hPa',
                          't_surf': '$T_S$'}
    
    return field_name if field_name not in field_name_aliases.keys() else field_name_aliases[field_name]

def plot_global_mean(data: dict,
                     field_name: str):
    
    ''' Function to plot area-averaged global average for a quantity. '''

    # Define helper function for area-averaging
    weighted_mean = lambda x: x.weighted(np.cos(np.deg2rad(x['grid_yt']))).mean(['grid_xt', 'grid_yt'])
    # Get field name metadata
    long_name, units = visualization.field_properties(field_name)

    # Initialize figure
    fig, ax = plt.subplots(figsize=(5, 2), dpi=144)
    # Iterate over experiments
    for configuration_name in data.keys():
        assert field_name in data[configuration_name].data_vars, f'Data variable {field_name} not available. '
        # Perform area-averaging
        TMP = weighted_mean(data[configuration_name][field_name])
        # Plot the global mean
        ax.plot(TMP['time'], TMP, label=configuration_name.upper())
        
        # Add horizontal line
        if np.nanmin(TMP) < 0 and np.nanmax(TMP) > 0:
            ax.axhline(0, c='k', lw=0.5)

    ax.legend(frameon=False, fontsize=9)
    ax.set_title(long_name, ha='left', loc='left', fontsize=10)
    ax.set_xlabel('Model year')
    ax.set_ylabel(f'{units}')

def plot_experiment_means_1D(ax,
                             data: dict,
                             model_name: str,
                             experiment_names: list[str],
                             field_name: str,
                             dimension_x_axis: str='grid_xt',
                             dimension_y_axis: str='grid_yt',
                             experiment_comparison: str|None=None,
                             plot_TC_histogram: bool=False,
                             year_range: tuple[int, int] | None=None,
                             plot_std: bool=False,
                             savefig: bool=False):

    # Boolean to dictate whether plot is part of a larger figure, or a standalone
    individual_plot = True if ax is None else False

    ''' Method to plot results from two different experiments and obtain their difference. '''
    
    # Define rolling function
    rolling_degrees = 20 # in degrees
    model_resolution = 0.5 # in degrees
    apply_rolling_mean = lambda x: x.rolling({dimension_y_axis: int(np.round(rolling_degrees / model_resolution))}, 
                                              center=True).mean(dimension_y_axis)

    ''' Figure structure parameters. '''
    # Overall figure structure
    ncols = len(experiment_names) + 1
    subplot_width, subplot_height = 4, 1.5
    # Tick locations
    axis_tick_locations = {'grid_xt': np.arange(0, 360, 60),
                           'grid_yt': np.arange(-90, 90 + 30, 30),
                           'pfull': np.array([200, 500, 700, 850, 1000])}

    ''' Gather data and generate normalizations. '''
    axes, axes_data, axes_data_mean = {}, {}, {}
    # Define factors for certain fields
    factor = 100 if 'cld_amt' in field_name else 1
    # Load datasets into a temporary container
    TMP = {}
    # Gather data
    for experiment_name in experiment_names:
        TMP[experiment_name] = data[model_name][experiment_name]
    for col_number, experiment_name in enumerate(experiment_names):
        if dimension_x_axis == 'grid_xt' and dimension_y_axis == 'grid_yt':
            # Get experiment data
            axes_data[experiment_name] = TMP[experiment_name][field_name] * factor
            # Define weights
            weights = np.cos(np.deg2rad(axes_data[experiment_name][dimension_y_axis]))
            # Apply weights
            zonal_mean_weighted = axes_data[experiment_name].mean('time').weighted(weights)
            axes_data_mean[experiment_name] = zonal_mean_weighted.mean(dim=dimension_x_axis)
            # Apply rolling average
            axes_data_mean[experiment_name] = apply_rolling_mean(axes_data_mean[experiment_name])
        else:
            print('Please define axis dimension names correctly. Exiting.')
            sys.exit()
            

    # Obtain normalization and colormap values
    norm, cmap = visualization.norm_cmap(axes_data, field=field_name)

    ''' Experiment plots. '''
    dpi = 300 if savefig else 144 # Define resolution based on whether or not a figure is savedd
    # Initialize figure
    if individual_plot:
        fig, ax = plt.subplots(figsize=(subplot_width, subplot_height), dpi=dpi)
    
    # Define dashed line parameters
    linestyle_dash = (0, (5, 7))
    
    # Draw reference line
    ax.axvline(0, c='k', lw=0.5)
    # Iterate over subplots
    for col_number, experiment_name in enumerate(experiment_names):
        # Get experiment name shorthand
        experiment_name_alias = get_experiment_alias(experiment_name)
        im = ax.plot(axes_data_mean[experiment_name][dimension_y_axis], 
                     axes_data_mean[experiment_name],
                     label=experiment_name_alias,
                     lw=1,
                     ls='--',
                     zorder=1)

    # Plot metadata
    long_name, units = visualization.field_properties(field_name)
    long_name_alias = get_field_name_alias(field_name)
    
    if individual_plot:
        ax.set_xlabel('Latitude', labelpad=10)
        ax.set_ylabel(f'{long_name_alias} [{units}]', va='bottom', labelpad=10)
        legend_ncols = min([3, len(experiment_names)])
        ax.legend(bbox_to_anchor=(0.5, 1.15), ncols=legend_ncols, loc='lower center', frameon=False)
    
    ''' Difference plotting. '''
    # Automatically compare if only 2 experiments are provided
    experiment_comparison = True if len(experiment_names) == 2 else experiment_comparison
    if experiment_comparison:
        ax_difference = ax.twinx()
        # Obtain difference
        difference = axes_data[experiment_names[0]] - axes_data[experiment_names[1]]
        # Get mean difference
        difference_mean = difference.mean(['grid_xt', 'time'])
        difference_mean = np.cos(np.deg2rad(difference_mean[dimension_y_axis])) * difference_mean
        difference_mean = apply_rolling_mean(difference_mean)
        # Get standard deviation in differences
        difference_std = difference.mean('time').std('grid_xt')
        difference_std = np.cos(np.deg2rad(difference_std[dimension_y_axis])) * difference_std
        difference_std = apply_rolling_mean(difference_std)
        # Obtain normalization and colormap values
        norm, cmap = visualization.norm_cmap(difference, field=field_name)
        ax_difference.plot(difference_mean[dimension_y_axis],
                           difference_mean,
                           lw=2,
                           c='k',
                           zorder=2)
        # Show standard deviation
        if plot_std:
            ax_difference.fill_between(x=difference_std[dimension_y_axis],
                                    y1=difference_mean-difference_std,
                                    y2=difference_mean+difference_std,
                                    alpha=0.1,
                                    fc='k')
        
        ax_difference.axhline(0, c='k', lw=1)
        if individual_plot:
            ax_difference.set_ylabel(f'$\Delta$ ({long_name_alias}) [{units}]', va='bottom', rotation=270, labelpad=10)
        
    # Modify axis spine properties to differentiate plot types
    ax.spines['left'].set_visible(False)
    ax.yaxis.label.set_color('tab:blue')
    ax.tick_params(axis='y', colors='tab:blue')
    ax_difference.spines['left'].set_color('tab:blue')
    ax_difference.spines['left'].set_linestyle(linestyle_dash)
    ax.set_xticks(axis_tick_locations[dimension_y_axis])
    ax.set_xlim([-90, 90])
    
    # Get zonal mean distributions
    if plot_TC_histogram:
        visualization.TC_density_histogram(ax=ax,
                                           model_name=model_name,
                                           experiment_names=experiment_names,
                                           bin_size=2.5,
                                           axis_depth=0.15,
                                           year_range=year_range,
                                           orientation='vertical',
                                           basin_name=None)
        
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-CLIMO/figs'
    figure_identifier = 'zonal_time_mean'
    experiment_identifier = '-'.join([experiment_name for experiment_name in experiment_names])
    filename = f'{figure_identifier}.model_name.{model_name}.experiment_name.{experiment_identifier}.field_name.{field_name}.png'

    if savefig:
        storage_pathname = os.path.join(storage_dirname, filename)
        plt.savefig(storage_pathname, transparent=True, dpi=dpi, bbox_inches='tight')

def plot_experiment_means_2D(data: dict,
                             model_name: str,
                             experiment_names: list[str],
                             field_name: str,
                             dimension_x_axis: str='grid_xt',
                             dimension_y_axis: str='grid_yt',
                             plot_type: str='pcolormesh',
                             plot_TC_histogram: bool=False,
                             year_range: tuple[int, int] | None=None,
                             savefig: bool=False):


    ''' Method to plot results from two different experiments and obtain their difference. '''

    assert len(experiment_names) == 2, '2 experiments must be provided to generate a difference.'
    dpi = 300 if savefig else 144
    
    # Define rolling function
    rolling_degrees = 20 # in degrees
    model_resolution = 0.5 # in degrees
    apply_rolling_mean = lambda x: x.rolling({dimension_y_axis: int(np.round(rolling_degrees / model_resolution))}, 
                                              center=True).mean(dimension_y_axis)

    ''' Figure structure parameters. '''
    # Overall figure structure
    ncols = len(experiment_names) + 1
    subplot_width, subplot_height = 4, 3
    # Tick locations
    axis_tick_locations = {'grid_xt': np.arange(0, 360, 60),
                           'grid_yt': np.arange(-90, 90 + 30, 30),
                           'pfull': np.array([200, 500, 700, 850, 1000])}
    # Colorbar parameters
    cax_width = 0.5
    cax_diff_width = 1
    cax_y_position = -1

    ''' Gather data and generate normalizations. '''
    axes, axes_data = {}, {}
    # Load datasets into a temporary container
    TMP = {}
    # Gather data
    for experiment_name in experiment_names:
        TMP[experiment_name] = data[model_name][experiment_name]
    for col_number, experiment_name in enumerate(experiment_names):
        if dimension_x_axis == 'grid_xt' and dimension_y_axis == 'grid_yt':
            axes_data[experiment_name] = TMP[experiment_name][field_name].mean('time') 
        elif dimension_x_axis == 'grid_yt' and dimension_y_axis == 'pfull':
            axes_data[experiment_name] = TMP[experiment_name][field_name].mean(['grid_xt', 'time']).sel(pfull=slice(100, 1000)) 
        else:
            print('Please define axis dimension names correctly. Exiting.')
            sys.exit()
            
    # Obtain normalization and colormap values
    norm, cmap = visualization.norm_cmap(axes_data, field=field_name)

    ''' Experiment plots. '''
    # Initialize figure
    fig, gs = [plt.figure(figsize=(ncols * subplot_width, subplot_height), dpi=dpi), 
               matplotlib.gridspec.GridSpec(nrows=2, ncols=ncols, height_ratios=(1, 0.075))]
    # Iterate over subplots
    for col_number, experiment_name in enumerate(experiment_names):
        ax = fig.add_subplot(gs[0, col_number])
        axes[f'{field_name.upper()}-{experiment_name.upper()}'] = ax
        if plot_type == 'pcolormesh':
            im = ax.pcolormesh(axes_data[experiment_name][dimension_x_axis], 
                               axes_data[experiment_name][dimension_y_axis], 
                               axes_data[experiment_name],
                               norm=norm,
                               cmap=cmap)
        else:
            im = ax.contourf(axes_data[experiment_name][dimension_x_axis], 
                             axes_data[experiment_name][dimension_y_axis], 
                             axes_data[experiment_name],
                             norm=norm,
                             cmap=cmap,
                             levels=len(norm.boundaries))

    # Plot raw plot colorbar
    cax_wrapper = fig.add_subplot(gs[1, 0:2])
    cax_wrapper.set_axis_off()
    cax = cax_wrapper.inset_axes([(1 - cax_width)/2, cax_y_position, cax_width, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
    long_name, units = visualization.field_properties(field_name)
    colorbar.set_label(f'{long_name} [{units}]', labelpad=10)
    cax.ticklabel_format(scilimits=(-3, 3))
    
    ''' Difference plot. '''
    ax = fig.add_subplot(gs[0, -1])
    experiment_name = f'{experiment_names[0].upper()}-{experiment_names[1].upper()}'
    difference = axes_data[experiment_names[0]] - axes_data[experiment_names[1]]
    
    norm_difference, cmap_difference = visualization.norm_cmap(difference, field=field_name)
        
    if plot_type == 'pcolormesh':
        im = ax.pcolormesh(difference[dimension_x_axis], 
                           difference[dimension_y_axis], 
                           difference, 
                           norm=norm_difference, 
                           cmap=cmap_difference)
    else:
        im = ax.contourf(difference[dimension_x_axis], 
                         difference[dimension_y_axis], 
                         difference, 
                         norm=norm_difference, 
                         cmap=cmap_difference,
                         levels=len(norm.boundaries))
        
    axes[f'{field_name}-{experiment_name}'] = ax
    axes_data[f'{field_name}-{experiment_name}'] = difference
    
    # Zonal mean comparison
    ax_zonal_mean = ax.inset_axes([1, 0, 0.2, 1])
    zonal_mean_difference = difference.mean(dimension_x_axis)
    if dimension_y_axis == 'grid_yt':
        zonal_mean_difference_weighted = np.cos(np.deg2rad(difference[dimension_y_axis])) * zonal_mean_difference
        zonal_mean_difference_weighted_smoothed = apply_rolling_mean(zonal_mean_difference_weighted)
    else:
        zonal_mean_difference_weighted_smoothed = zonal_mean_difference
    
    ax_zonal_mean.plot(zonal_mean_difference_weighted_smoothed, zonal_mean_difference_weighted_smoothed[dimension_y_axis])
    
    if dimension_x_axis == 'grid_xt':
        ax_zonal_mean.axhline(0, c='k', lw=0.5)
    ax_zonal_mean.axvline(0, c='k', lw=0.5)
    
    ax_zonal_mean.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    ax_zonal_mean.ticklabel_format(scilimits=(-3, 3))
    [ax_zonal_mean.get_xticklabels()[i].set_visible(False) for i, l in enumerate(ax_zonal_mean.get_xticks()) if l == 0]
    ax_zonal_mean.set_yticks([])
    if dimension_y_axis == 'pfull':
        ax_zonal_mean.set_ylim([100, 1000])
        ax_zonal_mean.set_ylim(ax_zonal_mean.get_ylim()[::-1])
    
    # Get zonal mean distributions
    if plot_TC_histogram:
        visualization.TC_density_histogram(ax=ax_zonal_mean,
                                           model_name=model_name,
                                            experiment_names=experiment_names,
                                            bin_size=2.5,
                                            axis_depth=0.25,
                                            year_range=year_range,
                                            basin_name=None)
    
    # Colorbar for difference plot
    cax_diff_wrapper = fig.add_subplot(gs[1, -1])
    cax_diff_wrapper.set_axis_off()
    cax_diff = cax_diff_wrapper.inset_axes([(1 - cax_diff_width)/2, cax_y_position, cax_diff_width, 1])
    cax_diff.ticklabel_format(scilimits=(-3, 3))
    colorbar_diff = fig.colorbar(matplotlib.cm.ScalarMappable(norm_difference, cmap_difference), 
                                 cax=cax_diff, 
                                 orientation='horizontal')
    colorbar_diff.set_label(f'{long_name} [{units}]', labelpad=10)
    [l.set_visible(False) for i, l in enumerate(cax_diff.get_xticklabels()) if i % 2 == 1]
    
    for ax_index, ax_object in enumerate(axes.items()):
        ax_name, ax = ax_object
        ax_title_str = f"({visualization.get_alphabet_letter(ax_index)}) {'-'.join(ax_name.split('-')[1:])}"
        ax.set_title(ax_title_str, ha='left', loc='left', fontsize=10)

        # Draw equatorial reference line
        if dimension_x_axis == 'grid_xt':
            ax.axhline(0, c='k', lw=0.5)
        elif dimension_x_axis == 'grid_yt':
            ax.axvline(0, c='k', lw=0.5)
            
        ax.set_xticks(axis_tick_locations[dimension_x_axis])
        ax.set_yticks(axis_tick_locations[dimension_y_axis])
        
        if dimension_y_axis == 'pfull':
            ax.set_ylim(ax.get_ylim()[::-1])
    
    fig.tight_layout()
    
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-AQP/figs'
    figure_identifier = 'planar_time_mean'
    experiment_identifier = '-'.join([experiment_name for experiment_name in experiment_names])
    filename = f'{figure_identifier}.model_name.{model_name}.experiment_name.{experiment_identifier}.field_name.{field_name}.png'

    if savefig:
        storage_pathname = os.path.join(storage_dirname, filename)
        plt.savefig(storage_pathname, transparent=True, dpi=dpi, bbox_inches='tight')
    
def plot_multiexperiment_means_1D(experiment_configurations: dict,
                                  field_name: str,
                                  data_type: str='month',
                                  savefig: bool=False):
    
    ''' 
    Generalization of `plot_experiment_means_1D` to multiple experiments in a single figure. 
    
    Note: `experiment_configurations` must be a nested dictionary with format:
    {`EXPERIMENT_ALIAS`: {'name': `MODEL_NAME`:`EXPERIMENT_A-EXPERIMENT_B,
                          'year_range': (`YEAR_A`, `YEAR_B`)}}
    '''
    
    # Container list for configuration aliases
    configuration_alias_container = []
    
    ''' Figure definition. '''
    # Define figure structure parameters
    subplot_width, subplot_height = 4, 1.5 # subplot dimensions
    nrows = int(np.ceil(len(experiment_configurations.keys()) / 2)) # number of rows
    ncols = 2 # number of columns

    # Initialize figure
    dpi = 300 if savefig else 144
    fig, gs = [plt.figure(figsize=(ncols * subplot_width, nrows * subplot_height), dpi=dpi),
            matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, wspace=0.375, hspace=0.875)]
    
    # Iterate over each configuration
    for axis_index, experiment_configuration_item in enumerate(experiment_configurations.items()):
        
        # Get configuration name and values
        experiment_configuration_name, experiment_configuration = experiment_configuration_item

        # Get iterand subplot row and column number from index
        nrow = axis_index // ncols
        ncol = axis_index % ncols

        print(f'Axis index: {axis_index}, row number: {nrow}, column number: {ncol}')
        
        ax = fig.add_subplot(gs[nrow, ncol])
        # Scrape configuration name
        configuration_name = experiment_configuration['name']
        # Obtain model and experiment names by splitting repeatedly
        model_name, experiment_names = configuration_name.split(':')
        experiment_names = experiment_names.split('-')
        # Define year range
        year_range = experiment_configuration['year_range']

        data = analysis.load_GCM_data(model_name=model_name,
                                    experiment_names=experiment_names,
                                    field_names=[field_name],
                                    year_range=year_range,
                                    data_type=data_type)

        plot_experiment_means_1D(ax=ax,
                                 data=data,
                                 model_name=model_name,
                                 experiment_names=experiment_names,
                                 field_name=field_name,
                                 plot_TC_histogram=True,
                                 year_range=year_range)

        # Subplot annotation
        subplot_ax = ax.twiny() # create dummy axis for vertical stacking purposes
        subplot_ax.set_axis_off()
        subplot_identifier = visualization.get_alphabet_letter(axis_index)
        ann = subplot_ax.annotate(f'{subplot_identifier}) {experiment_configuration_name}', xy=(0, 1.15), xycoords='axes fraction', 
                                  fontsize=9, va='bottom', ha='left'),
                                #   bbox={'facecolor': 'white', 'alpha': 0.75, 'pad': 2, 'edgecolor': 'white'})
        
        # Append experiment name for future plot storage
        configuration_alias_container.append(experiment_configuration_name)

    # Plot metadata
    long_name, units = visualization.field_properties(field_name)
    long_name_alias = get_field_name_alias(field_name)

    fig.supxlabel('Latitude', y=0.025, va='top', fontsize=11)
    fig.supylabel(f'{long_name_alias} [{units}]', x=0.025, c='tab:blue', fontsize=11)
    ax.annotate(f'$\Delta$ ({long_name_alias}) [{units}]', xy=(0.96, 0.6), xycoords='figure fraction', va='center', rotation=270, fontsize=11)

    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-CLIMO/figs'
    figure_identifier = 'zonal_time_mean'
    experiment_identifier = '-'.join([configuration_alias for configuration_alias in configuration_alias_container])
    filename = f'{figure_identifier}.model_name.{model_name}.experiment_names.{experiment_identifier}.field_name.{field_name}.png'

    if savefig:
        storage_pathname = os.path.join(storage_dirname, filename)
        plt.savefig(storage_pathname, transparent=True, dpi=dpi, bbox_inches='tight')

def TC_parameter_distribution(track_data: dict,
                              experiment_names: str | list[str]=None,
                              distribution_parameter: str='min_slp',
                              sampling_method: str='all',
                              plot_median: bool=False,
                              normalize_by_year: bool=False,
                              savefig: bool=False):

    ''' Plot distribution histogram for a given parameter. Intensity parameters are shown at the provided sampling method. '''

    assert distribution_parameter in ['min_slp', 'max_wind', 'duration', 'center_lat'], 'Intensity parameter must be `min_slp`, `max_wind`, `center_lat`, or `duration`.'

    # Define TC sampling method and corresponding function
    # 'genesis' == TC genesis, which is assigned to the first timestamp of 'raw')
    # 'LMI'     == lifetime maximum intensity, which is assigned to 'unique')
    # 'lysis'   == TC lysis, which is assigned to the last timestamp of 'raw')
    # 'all'     == all TC timestamps, which is assigned to 'raw')
    sampling_types = {'genesis': lambda x, name, parameter: np.array([n.sort_values('time').iloc[0][parameter].item() 
                                                                      for _, n in x[name]['raw'].groupby('storm_id')]),
                      'lysis': lambda x, name, parameter: np.array([n.sort_values('time').iloc[-1][parameter].item() 
                                                                    for _, n in x[name]['raw'].groupby('storm_id')]),
                      'LMI': lambda x, name, parameter: x[name]['unique'][parameter],
                      'all': lambda x, name, parameter: x[name]['raw'][parameter]} 
    
    # Define metadata
    parameter_properties = {'min_slp': {'long_name': 'Minimum sea-level pressure [hPa]',
                                        'histogram_interval': 2.5},
                            'max_wind': {'long_name': 'Maximum 10 m wind speed [m s$^{-1}$]',
                                         'histogram_interval': 1},
                            'duration': {'long_name': 'Duration of TC lifetime [days]',
                                         'histogram_interval': 5},
                            'center_lat': {'long_name': 'Latitude',
                                           'histogram_interval': 2.5}}
    # Initialize extrema containers
    vmin, vmax = np.nan, np.nan
    # Container list for configuration aliases
    configuration_alias_container = []
    
    # Define colormap
    colormap = 'tab10'
    colors = matplotlib.colormaps[colormap](np.linspace(0.1, 0.9, len(track_data.keys())))

    # Initialize container for number of years per experiment, for temporal normalization purposes
    num_years = {}

    # Iterate over experiemnts to get extrema
    for configuration_name in track_data.keys():
        model_name, experiment_name = configuration_name.split('-')
        experiment_data = sampling_types[sampling_method](track_data, configuration_name, distribution_parameter)
        # Get extrema
        vmin = experiment_data.min() if (experiment_data.min() < vmin or np.isnan(vmin)) else vmin
        vmax = experiment_data.max() if (experiment_data.max() > vmax or np.isnan(vmax)) else vmax
        # Get number of years for each experiment
        num_years[experiment_name] = np.round((track_data[configuration_name]['unique']['time'].max() - track_data[configuration_name]['unique']['time'].min()).days / 365)
    
    # Get intensity bins
    distribution_bins = np.arange(vmin, 
                                  vmax + parameter_properties[distribution_parameter]['histogram_interval'], 
                                  parameter_properties[distribution_parameter]['histogram_interval'])

    # Initialize container to collect histogram counts for y-axis delimiting
    counts = np.array([])
    # Determine if a PDF or raw count will be output 
    density = True if normalize_by_year else False 
    # Plot histograms
    dpi = 300 if savefig else 144
    fig, ax = plt.subplots(figsize=(6, 2), dpi=dpi)
    
    for configuration_number, configuration_name in enumerate(track_data.keys()):
        model_name, experiment_name = configuration_name.split('-')
        # Get iterand experiment data
        experiment_data = sampling_types[sampling_method](track_data, configuration_name, distribution_parameter)
        # Define dataset label
        TC_counts = len(track_data[configuration_name]['raw']['storm_id'].unique())
        TC_counts = TC_counts / num_years[experiment_name] if normalize_by_year else TC_counts # Normalize by number of years, if chosen
        if distribution_parameter == 'center_lat':
            if normalize_by_year:
                label = f'{get_experiment_alias(experiment_name).upper()} (N={TC_counts:.1f} yr$^{{-1}}$)'
            else:
                label = f'{get_experiment_alias(experiment_name).upper()} (N={TC_counts:.1f})'
        else:
            label = f'{get_experiment_alias(experiment_name).upper()} ({np.median(experiment_data):.1f})'
        # Get data distribution
        experiment_counts, _, im = ax.hist(experiment_data, bins=distribution_bins, 
                                           histtype='step', label=label, density=density,
                                           color=colors[configuration_number], linewidth=1.5)
        # Plot distribution median
        if plot_median:
            ax.axvline(np.median(experiment_data), lw=0.5, ls='--', c=colors[configuration_number])
        # Extend the counts array to hold new counts
        counts = np.concatenate([counts, experiment_counts])
        print(f'Statistics for {experiment_name}: mean = {experiment_data.mean():.1f} +/- {experiment_data.std():.1f}; 10th percentile: {np.quantile(experiment_data, 0.1):.1f}; 90th percentile: {np.quantile(experiment_data, 0.9):.1f}')

    # Plot metadata
    ylabel_str = 'Probability density' if density else 'Count'
    ax.set_xlabel(parameter_properties[distribution_parameter]['long_name'], labelpad=10)
    ax.set_ylabel(ylabel_str, labelpad=10)
    ax.set_ylim([0, np.max(counts)*1.1])

    position_legend = (ax.get_position().x0, ax.get_position().y1)
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(0, 1.025), fontsize=7, labelspacing=0.25)
    
    storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-CLIMO/figs'
    figure_identifier = 'TC_distribution'
    configuration_identifier = '-'.join([configuration_alias for configuration_alias in configuration_alias_container])
    filename = f'{figure_identifier}.configuration_names.{configuration_identifier}.field_name.{distribution_parameter}.png'

    if savefig:
        storage_pathname = os.path.join(storage_dirname, filename)
        plt.savefig(storage_pathname, transparent=True, dpi=dpi, bbox_inches='tight')

def TC_tracks_raw(track_data: dict,
                  configuration_name: str,
                  track_format: str='plot',
                  projection=ccrs.PlateCarree(central_longitude=180)):

    ''' Plots TC tracks for a given model-experiment configuration. '''

    assert track_format in ['plot', 'scatter'], 'Plot format must be `plot` or `scatter`.'

    # configuration_name = list(track_data.keys())[0]
    model_name, experiment_name = configuration_name.split('-')
    
    track_data_reconfigured = {model_name: {experiment_name: track_data[configuration_name]}}
    
    # Normalize scatter plot colormap using minimum sea-level pressure as the intensity parameter
    norm = matplotlib.colors.BoundaryNorm(np.arange(900, 1000, 10), 256)
    cmap = 'viridis_r'
    
    projection, reference_projection = projection, ccrs.PlateCarree()
    figsize = (7, 3)
    histogram_axis_depth = 0.05
    
    
    linewidth = 0.5
    linecolor = 'k'
    
    fig, ax = plt.subplots(figsize=figsize, dpi=144, subplot_kw={'projection': projection})

    # Pull experiment-specific track data
    experiment_track_data = track_data[configuration_name]['raw']
    # Iterate over all storms and plot tracks
    for TC_ID, TC in experiment_track_data.groupby('storm_id'):
        TC = TC.sort_values('time')
    
        # If a TC crosses the "International Date Line" once, split the track and plot
        # Else if it does it more than twice, skip it
        # Else, plot normally

        if track_format == 'plot':
            if (abs(TC['center_lon'].diff()) > 180).sum() == 1:
                split_index = np.argmax(abs(TC['center_lon'].diff()) > 180)
                TC_A, TC_B = TC[:split_index], TC[split_index:]
                ax.plot(TC_A.center_lon, TC_A.center_lat, c=linecolor, lw=linewidth, transform=reference_projection)
                ax.plot(TC_B.center_lon, TC_B.center_lat, c=linecolor, lw=linewidth, transform=reference_projection)
            elif (abs(TC['center_lon'].diff()) > 180).sum() == 0:
                TC_im = ax.plot(TC.center_lon, TC.center_lat, lw=linewidth, c=linecolor, transform=ccrs.PlateCarree())
        else:
            ax.scatter(TC.center_lon, TC.center_lat, s=2.5, norm=norm, c=TC['min_slp'], cmap='viridis_r', transform=ccrs.PlateCarree())
    
    if projection == ccrs.PlateCarree(central_longitude=180):
        ax.set_extent([-180, 180, -90, 90])
    else:
        ax.set_extent([0, 360, -90, 90], crs=reference_projection)
    ax.gridlines(draw_labels=['bottom', 'left'], xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 15), ls=':', alpha=0.5)
    ax.set_title(f'{configuration_name.upper()}', ha='left', loc='left', fontsize=10)
    
    # Get zonal mean distribution
    if projection == ccrs.PlateCarree(central_longitude=180):
        visualization.TC_density_histogram(ax=ax,
                                        model_name=model_name,
                                        experiment_names=experiment_name,
                                        track_dataset=track_data_reconfigured,
                                        bin_size=5,
                                        axis_depth=histogram_axis_depth,
                                        basin_name=None)

    # Plot colorbar if the scatter method is chosen
    if track_format == 'scatter':
        cax = ax.inset_axes([1 + histogram_axis_depth + 0.05, 0, 0.025, 1])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
        colorbar.set_label('Minimum sea-level pressure [hPa]', labelpad=20, rotation=270)

def TC_track_density(track_data: dict,
                     experiment_comparison: str|None=None,
                     bin_resolution: int=5,
                     output_unit: str='raw',
                     add_land: bool=False,
                     diagnostic: bool=False):

    ''' Generate map of spatial frequency of TC occurrence. Provides a comparison option to compare two experiments. '''

    diagnostic_tag = f'[TC_track_density()]'
    zonal_mean_axis_depth = 0.05

    if experiment_comparison:
        # Ensure both experiments are in the provided track data
        experiment_names = experiment_comparison.split(':')
        if diagnostic:
            print(f'{diagnostic_tag} Experiment comparison input: {experiment_comparison}; experiment names: {experiment_names} of length {len(experiment_names)}')
        assert (len(experiment_names) == 2), 'Experiment names must be separated by a hyphen and not contain any hyphens.'
        assert (experiment_names[0] in track_data.keys()) & (experiment_names[1] in track_data.keys()), 'Experiment names must be in the track data.'

        axis_width, axis_height = 6, 3
        height_ratios = (1, 0.05)
        hspace = 0 if add_land else 0.125
        wspace = 0.125
        fig, gs = [plt.figure(figsize=(3 * axis_width, sum(axis_height * np.array(height_ratios) * (1 + hspace))), dpi=144, constrained_layout=True),
                   matplotlib.gridspec.GridSpec(nrows=2, ncols=3, hspace=hspace, wspace=wspace, height_ratios=height_ratios)]
        axes = {}

        ''' Generate density arrays for constituent experiments. '''
        # Re-format dictionary into the structure needed for `tc_analysis` scripts
        configuration_names = [s.split('-') for s in experiment_names]
        if diagnostic:
            print(f'{diagnostic_tag} Configuration names: {configuration_names}')
        
        # Initialize containers for track data
        # Assumes configurations being compared share a model name. This makes life easier for the rest of the script.
        model_name = experiment_names[0].split('-')[0]
        modified_track_data = {model_name: {}} 
        
        # Get configuration-specific TC activity
        TC_densities = {} 
        for configuration_name in experiment_names:
            if diagnostic:
                print(f'{diagnostic_tag} Experiment comparison, configuration {configuration_name}...')
        
            _, experiment_name = configuration_name.split('-')
            modified_track_data[model_name][experiment_name] = track_data[f'{model_name}-{experiment_name}']

            TC_densities[configuration_name] = tc_analysis.TC_density(model_names=[model_name],
                                                                      experiment_names=[experiment_name],
                                                                      track_dataset=modified_track_data,
                                                                      bin_resolution=bin_resolution)[model_name][experiment_name]
        
        ''' Experiment output plots. '''
        experiment_TC_densities = {}
        # Get data for both experiment to obtain a normalization and colormap
        # Normalize density values by TC count per day (divide by 365) as a percentage of the year (mult. by 100)
        experiment_TC_densities = {experiment_name: TC_densities[experiment_name] * (100 / 365) for experiment_name in experiment_names}
        # Get normalization values
        norm, _ = visualization.norm_cmap(experiment_TC_densities, field='density', num_bounds=8)
        cmap = visualization.cmap_white_adjust('Blues', levels=len(norm.boundaries))
        
        for experiment_index, experiment_name in enumerate(experiment_TC_densities.keys()):
            axes[experiment_name] = fig.add_subplot(gs[0, experiment_index], projection=ccrs.PlateCarree(central_longitude=180))
            
            # Plot the spatial distribution
            axes[experiment_name].pcolormesh(experiment_TC_densities[experiment_name].columns, 
                                              experiment_TC_densities[experiment_name].index, 
                                              experiment_TC_densities[experiment_name], 
                                              norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
            
        # Plot shared colorbar
        cax_wrapper = fig.add_subplot(gs[1, 0:2])
        cax_wrapper.set_axis_off()
        cax = cax_wrapper.inset_axes([0.25, -1, 0.5, 1])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
        colorbar.set_label('% of time with TC', labelpad=10)

        ''' Difference plot. '''
        
        # Get difference values
        experiment_TC_densities[experiment_comparison] = experiment_TC_densities[experiment_names[0]] - experiment_TC_densities[experiment_names[1]]
        # Get normalization values
        norm, cmap = visualization.norm_cmap(experiment_TC_densities[experiment_comparison], field='density', num_bounds=8)
        # Plot difference
        axes[experiment_comparison] = fig.add_subplot(gs[0, -1],  projection=ccrs.PlateCarree(central_longitude=180))
        # Plot the spatial distribution
        axes[experiment_comparison].pcolormesh(experiment_TC_densities[experiment_comparison].columns, 
                                              experiment_TC_densities[experiment_comparison].index, 
                                              experiment_TC_densities[experiment_comparison], 
                                              norm=norm, cmap=cmap, transform=ccrs.PlateCarree())

        # Plot difference colorbar
        cax_wrapper = fig.add_subplot(gs[1, -1])
        cax_wrapper.set_axis_off()
        cax = cax_wrapper.inset_axes([0, -1, 1, 1])
        cax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
        colorbar.set_label('$\Delta$ (% of time with TC) ', labelpad=10)
        
        # Perform common axis operations here
        histogram_colors = ['tab:blue', 'tab:green']
        for index, (ax_name, ax) in enumerate(axes.items()):

            draw_labels = ['bottom'] if index > 0 else ['bottom', 'left']

            # Draw gridlines
            ax.gridlines(draw_labels=draw_labels, 
                         xlocs=np.arange(-180, 180 + 90, 45),
                         ylocs=np.arange(-90, 90, 30), 
                         ls=':', alpha=0.5)
            ax.set_ylim([-60, 60]) if add_land else ax.set_ylim([-90, 90])
            ax.set_title(f'{ax_name.upper()}', ha='left', loc='left', fontsize=10)

            # Plot TC density per latitude bin
            if ':' not in ax_name:
                model_name, experiment_name = ax_name.split('-')
                visualization.TC_density_histogram(ax=ax,
                                                    model_name=model_name,
                                                    experiment_names=experiment_name,
                                                    track_dataset=modified_track_data,
                                                    bin_size=bin_resolution,
                                                    axis_depth=zonal_mean_axis_depth,
                                                    basin_name=None,
                                                    colors=histogram_colors[index])
            else:
                model_name = ax_name.split('-')[0] # assumes both model names are the same, which is likely the case
                experiment_names = [ax_name_entry.split('-')[1] for ax_name_entry in ax_name.split(':')]
                visualization.TC_density_histogram(ax=ax,
                                                   model_name=model_name,
                                                   experiment_names=experiment_names,
                                                   track_dataset=modified_track_data,
                                                   bin_size=bin_resolution,
                                                   axis_depth=zonal_mean_axis_depth,
                                                   basin_name=None,
                                                   colors=histogram_colors)


            # Perform additional cosmetic plot formatting
            if add_land: 
                ax.add_feature(cartopy.feature.LAND, fc=(0.875, 0.875, 0.875), zorder=8)
                ax.coastlines(zorder=9)
                # Ensure the axis spines stay above all other elements
                for k, spine in ax.spines.items():
                    spine.set_zorder(10)
            
        fig.tight_layout()
        
    else:
        
        for configuration_name in track_data.keys():
            model_name, experiment_name = configuration_name.split('-')
            
            # Re-format dictionary into the structure needed for `tc_analysis` scripts
            modified_track_data = {model_name: {experiment_name: track_data[configuration_name]}}
            
            TC_density = tc_analysis.TC_density(model_names=[model_name],
                                                experiment_names=[experiment_name],
                                                track_dataset=modified_track_data,
                                                bin_resolution=bin_resolution)[model_name]
        
            fig, ax = plt.subplots(figsize=(7, 2.5), dpi=144, subplot_kw={'projection': ccrs.PlateCarree()})
    
            # Normalize density values by TC count per day (divide by 365) as a percentage of the year (mult. by 100)
            experiment_TC_density = TC_density[experiment_name] * (100 / 365)
            
            # Get normalization values
            norm, cmap = visualization.norm_cmap(experiment_TC_density, field='density', num_bounds=8)
            cmap = visualization.cmap_white_adjust('Blues', levels=len(norm.boundaries))
    
            # Plot the spatial distribution
            ax.pcolormesh(experiment_TC_density.columns, experiment_TC_density.index, experiment_TC_density, 
                          norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
            # Get zonal mean distribution
            visualization.TC_density_histogram(ax=ax,
                                               model_name=model_name,
                                               experiment_names=experiment_name,
                                               track_dataset=modified_track_data,
                                               bin_size=bin_resolution,
                                               axis_depth=zonal_mean_axis_depth,
                                               basin_name=None)
    
            ax.gridlines(draw_labels=['bottom', 'left'], 
                         xlocs=np.arange(-180, 180, 45),
                         ylocs=np.arange(-90, 90, 30), 
                         ls=':', alpha=0.5)
            ax.set_ylim([-90, 90])
            
            ax.set_title(f'{model_name}-{experiment_name.upper()}', ha='left', loc='left', fontsize=10)
    
            cax = ax.inset_axes([1 + 2 * zonal_mean_axis_depth, 0, 0.025, 1])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
            
            colorbar.set_label('% of time with TC per bin', rotation=270, labelpad=20)