''' Import packages. '''
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import functools, importlib, os, pickle, collections, sys
import pandas as pd, xarray as xr, nc_time_axis
xr.set_options(keep_attrs=True)
# Visualization tools
import matplotlib, matplotlib.pyplot as plt
# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import derived, utilities, visualization, track_TCs, zonal_mean

def generate_field_plot_data(data,
                             field_name,
                             averaging_dimensions: str | list[str],
                             contour_levels: int=8,
                             extrema: tuple[int|float, int|float]|None=None):
    
    average = data[field_name].mean(averaging_dimensions)
    norm, cmap = visualization.norm_cmap(average, field=field_name, num_bounds=contour_levels, extrema=extrema)

    return average, norm, cmap

def plot_overturning_circulation(datasets: dict,
                                 configuration_name: str,
                                 contour_levels: int=8,
                                 single_level: int|float|None=None,
                                 single_level_color: str|None=None,
                                 extrema: tuple[int|float, int|float]|None=None,
                                 overlay: bool=False,
                                 ax: matplotlib.axes.Axes | None=None,
                                 **kwargs):

    # Flag to determine if this figure is recycled as part of a larger figure, or if it's a standalone
    standalone_fig = True if ax is None else False
    # Define a shorthand for the working dataset
    TMP = datasets[configuration_name]

    # Define the field to use for the main and secondary contour sets
    field_name, overlay_field_name = 'psi_m', 'ucomp'
    assert field_name in TMP.data_vars, 'Main field name (psi_m) not in the provided dataset.'

    ''' Load data. '''
    # Load data for the main contour set
    dataset, norm, cmap = generate_field_plot_data(TMP, 
                                                   field_name=field_name, 
                                                   averaging_dimensions=['time'],
                                                   contour_levels=contour_levels,
                                                   extrema=extrema)

    # Load data for the secondary (overlay) contour set
    overlay_dataset, overlay_norm, overlay_cmap = generate_field_plot_data(TMP, 
                                                                           field_name=overlay_field_name, 
                                                                           averaging_dimensions=['grid_xt', 'time'])

    ''' Figure parameters. '''
    # Determine if there should only be one level in the normalization depending on input
    if single_level is None:
        levels = len(norm.boundaries)
    else:
        levels = [single_level]
        norm, cmap = None, None

    # Initialize a new figure if one isn't provided
    if standalone_fig: fig, ax = plt.subplots(figsize=(5, 2.5), dpi=144)

    ''' Plotting. '''
    # Plot secondary contours (typically zonal wind)
    if overlay:
        overlay_contours = ax.contour(overlay_dataset.grid_yt, overlay_dataset.pfull, overlay_dataset, 
                                    norm=overlay_norm, colors=['grey'], alpha=0.5,
                                    linewidths=1, linestyles='--', levels=len(overlay_norm.boundaries))
        overlay_contour_labels = ax.clabel(overlay_contours)
    
    # Plot primary contour set
    if 'color' in kwargs: 
        kwargs['colors'] = kwargs['color']
    elif single_level is not None and single_level_color is not None:
        kwargs['colors'] = single_level_color
    im = ax.contour(dataset.grid_yt, dataset.pfull, dataset.T,
                    norm=norm, cmap=cmap, levels=levels,
                    **kwargs)

    # Plot reference equator line
    ax.axvline(0, c='k', lw=0.5)

    # Apply figure formatting if it's a standalone, else output the axis
    if standalone_fig:
        ax.set_xticks(np.arange(-60, 90, 30))
        ax.set_yticks([200, 500, 700, 850, 1000])
        
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel('Latitude [deg N]')
        ax.set_ylabel('Pressure level [hPa]')
        # ax.set_title(f'{format_experiment_name(configuration_name)}', loc='left', ha='left', fontsize=10)

        cax = ax.inset_axes([1.05, 0, 0.05, 1])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
        
        fig.tight_layout()

    return im

def multiplot_overturning_circulation(data: dict,
                                      configuration_names: dict[list[str]],
                                      psi_magnitudes: dict[str: int|float],
                                      colormap_name: str='viridis',
                                      diagnostic: bool=False):

    diagnostic_tag = f'[multiplot_overturning_circulation()]'

    ''' Method to generate multiple overlaid plots of overturning circulation. '''

    # Define colors to use for different configurations

    
    ''' Initialize plot. '''

    ncols = 2
    nrows = np.ceil(len(configuration_names) / ncols).astype(int)
    ax_width, ax_height = 4, 2
    dpi = 144
    fig, gs = [plt.figure(figsize=(ax_width * ncols, ax_height * nrows), dpi=dpi),
               matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols)]


    # Iterate over the selected configurations
    for index, config_name in enumerate(configuration_names):
                             
        nrow = index // ncols
        ncol = index - nrow * ncols
        subconfig_names = configuration_names[config_name]
        if diagnostic: print(f'Figure index: {index}; row: {nrow}, column: {ncol}, {subconfig_names}')

        # Initialize catalog to store plot data. This will be used for generating a legend.
        ims = {}

        ax = fig.add_subplot(gs[nrow, ncol])
        for subindex, subconfig_name in enumerate(subconfig_names):

            colormap_name = 'viridis' if ':CTL' in subconfig_name else 'bone'
            contour_colors = visualization.get_colormap_samples(number_of_samples=2,
                                                               colormap_name=colormap_name,
                                                               interval_minimum=0.25, 
                                                               interval_maximum=0.75)

            psi_magnitude = psi_magnitudes['CTL'] if ':CTL' in subconfig_name else psi_magnitudes['DIFF']

            if diagnostic:  print(f'\tPlotting subconfiguration: {subconfig_name}...')

            # Get configuration properties and make sure configuration names conform to structure and filter rules
            assert len(subconfig_name.split(':')) == 3, 'Revisit configuration name definitions such that the configuration is of form {model_name}:{experiment_name}:{experiment_type}.'
                        
            # Determine configuration-specific model and experiment name
            subconfig_model_name, subconfig_experiment_name = [subconfig_name.split(':')[0],
            subconfig_name.split(':')[1]]
        
            # Plot the negative streamfunction (southward transport)
            ims[subconfig_name] = plot_overturning_circulation(datasets=data,
                                                               configuration_name=subconfig_name,
                                                               single_level=-psi_magnitude,
                                                               single_level_color=contour_colors[subindex],
                                                               extrema=None,
                                                               ax=ax)  
            
            # Plot the positive streamfunction (southward transport)
            plot_overturning_circulation(datasets=data,
                                        configuration_name=subconfig_name,
                                        single_level=psi_magnitude,
                                        single_level_color=contour_colors[subindex],
                                        extrema=None,
                                        ax=ax)

            # Plot TC density histogram
            get_density_difference = True if 'DIFF' in subconfig_name else False
            model_name, experiment_names = [subconfig_name.split(':')[0],
                                            (f'CTL1990.{subconfig_name.split(':')[1]}', f'CTL1990_SWISHE.{subconfig_name.split(':')[1]}')]
            experiment_names = experiment_names if 'DIFF' in subconfig_name else experiment_names[0] # only use control experiment for non-difference plots
            zonal_mean.get_TC_track_histogram(ax=ax, 
                                              model_name=model_name, 
                                              experiment_names=experiment_names, 
                                              get_density_difference=get_density_difference, 
                                              colors=[contour_colors[subindex]],
                                              orientation='vertical',
                                              override_difference_colors=True)

            # Generate custom legend
            if nrow == 0: 
                handles = [matplotlib.lines.Line2D([0], [0], color=im.colors, lw=2) for im in ims.values()]
                labels = [config_name.split(':')[0] for config_name in ims.keys()]
                ax.legend(handles=handles, labels=labels, loc='upper right', frameon=False, fontsize=9)

        ''' Subplot metadata. '''
        # Plot ticks and limits
        ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        ax.set_yticks([200, 500, 700, 850, 1000])
        if ncol > 0: ax.set_yticklabels([])
        ax.set_xlim([-90, 90])
        ax.set_ylim(ax.get_ylim()[::-1])

        subplot_name = f'{subconfig_name.split(':')[1]}, CTL — SWISHE' if 'DIFF' in subconfig_name else f'{subconfig_name.split(':')[1]}, CTL'
        subplot_annotation = f'({visualization.get_alphabet_letter(index)}) {subplot_name}'
        ax.annotate(subplot_annotation, xy=(0.02, 0.95), xycoords='axes fraction', ha='left', va='top', fontsize=9)

    fig.tight_layout()

    # Plot labeling
    fig.supxlabel('Latitude [deg N]', y=-0.03, fontsize=11)
    fig.supylabel('Pressure level [hPa]', x=-0.03, fontsize=11)
