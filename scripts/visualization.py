import cartopy, cartopy.crs as ccrs
import numpy as np, pandas as pd, scipy as sp
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

def density_grid(data, model_name, bin_resolution=5):
    """
    Method to plot the spatial density of unique TC occurrences given a track data dictionary. 
    See tc_analysis.py --> tc_track_data() for more information.

    Args:
        data (dictionary): 3-tiered dictionary.
        model_name (str): name of model (usually 'AM2.5', 'HIRAM', 'FLOR'.)
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
            # Define shorthand for unique TC data
            dataset = data[model][experiment]['unique']
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
    year_min, year_max = density[model_name][run_control]['time_min'].unique().year.item(), density[model_name][run_control]['time_max'].unique().year.item()
    
    ''' Begin plotting. '''
    # Note: number of columns is dictated by number of experiments + 1 (for difference plot)
    ncols = len(density[model_name].keys()) + 1
    # Initialize figure and grid
    fig, grid = plt.figure(figsize=(14, 2.5), constrained_layout=True), matplotlib.gridspec.GridSpec(nrows=2, ncols=ncols, height_ratios=(0.95, 0.05))
    # Define longitudinal offset
    longitude_offset = 180
    # Define projections (working and reference projections)
    proj, proj_ref = ccrs.PlateCarree(central_longitude=longitude_offset), ccrs.PlateCarree()
    
    # Collect processed density maps. Note: the difference will be calculated for the last subplot.
    densities = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for model in density.keys():
        # Initialize subdictionary for the iterand model
        densities[model] = {}
        # Iterate over each given experiment
        for experiment in [run_control, run_experiment]:
            # Get longitude and latitude bins
            x, y = density[model][experiment].lon.unique(), density[model][experiment].lat.unique()
            # Get density array
            v = np.reshape(density[model][experiment]['count'].values, (len(x), len(y)))
            # Assign to dictionary for the given model/experiment configuration
            densities[model][experiment] = v.T/(year_max - year_min)
            
    # Initialize normalization and colormaps. Save into dictionary for future use in colorbars.
    norms, cmaps = {}, {}
    bounds = 8 # number of levels to bin for the normalization
    vmax = max([np.nanmax(sv) for k, v in densities.items() for sk, sv in v.items()]) # maximum value across 'run_control' and 'run_experiment' for all models
    # Calculate density differences for the model of interest (given by 'model_name' argument)
    run_difference = '{1} - {0}'.format(run_control, run_experiment)
    densities[model_name][run_difference] = densities[model_name][run_experiment] - densities[model_name][run_control]
    difference_extremum = max([abs(np.nanmin(densities[model_name][run_difference])), abs(np.nanmax(densities[model_name][run_difference]))])
    # Assign normalization and colormap values
    for experiment in densities[model_name].keys():
        if '-' in experiment:
            norms[experiment] = matplotlib.colors.BoundaryNorm(np.linspace(-difference_extremum, difference_extremum, bounds), 256)
            cmaps[experiment] = 'bwr'
        else:
            norms[experiment] = matplotlib.colors.BoundaryNorm(np.linspace(0, vmax, bounds+1), 256)
            cmaps[experiment] = cmap_white_adjust('Reds')
    # Pass 2: Iterate through density dictionary to plot each experiment's density data.
    for experiment_index, experiment in enumerate(densities[model_name].keys()):
        # Initialize subplot
        ax = fig.add_subplot(grid[0, experiment_index], projection=proj)
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
    
    ''' Create 2 colorbars: 1 for the experiments, another for the difference. '''
    # Colorbar ticklabel formatter
    fmt = lambda x, pos: '{:.1f}'.format(x)
    
    # Experiments
    cw_experiments = fig.add_subplot(grid[1, 0:2]) # 'cw' stands for 'cax_wrapper'
    cw_experiments.set_axis_off()
    cax_experiments = cw_experiments.inset_axes([0.25, 0, 0.5, 1])
    colorbar_experiments = fig.colorbar(matplotlib.cm.ScalarMappable(norms[run_control], cmaps[run_control]), 
                                        orientation='horizontal', cax=cax_experiments, format=matplotlib.ticker.FuncFormatter(fmt))
    # Differnece
    cw_difference = fig.add_subplot(grid[1, -1]) # 'cw' stands for 'cax_wrapper'
    cw_difference.set_axis_off()
    cax_difference = cw_difference.inset_axes([0, 0, 1, 1])
    colorbar_difference = fig.colorbar(matplotlib.cm.ScalarMappable(norms[run_difference], cmaps[run_difference]), 
                                        orientation='horizontal', cax=cax_difference, format=matplotlib.ticker.FuncFormatter(fmt))

def pdf(data, param='center_lat', num_bins=60):
    
    fig, ax = plt.subplots(figsize=(5, 2.5))
    
    # Define list of colors (one color per model) and line styles (one linestyle per experiment)
    colors = ['tab:blue', 'tab:green', 'tab:purple']
    linestyles = ['-', '--']
    
    # Hold list of PDF extrema for axis adjustment
    pdfs = []
    # Iterate over all models and experiments
    for model_index, model in enumerate(data.keys()):
        for experiment_index, experiment in enumerate(data[model].keys()):
            # Get unique values for the parameter passed
            out = data[model][experiment]['unique'][param].dropna()
            # Filter out nans and infs
            out = out.loc[np.isfinite(out.values)]
            # Get bin values and bin edges, don't plot
            n, bins, _ = ax.hist(out, bins=num_bins, histtype=u'step', density=True, lw=0) 
            # Build kernel density estimator to obtain smooth PDF curves
            kde = sp.stats.gaussian_kde(out)(bins)
            # Append to PDF list
            pdfs.append(np.nanmax(kde))
            # Plot the distribution
            ax.plot(bins, kde, lw=3, color=colors[model_index], linestyle=linestyles[experiment_index], 
                    label='{0}, {1}'.format(model, experiment))
            
    # Adjust axis limits
    ax.set_ylim([0, 1.1*max(pdfs)])
    if param == 'center_lat':
        # Curb latitudinal extent of TC monitoring to midlatitude extrema
        ax.set_xlim([-60, 60])
    elif param == 'duration':
        # Storms lasting > 30 d are unrealistic
        ax.set_xlim([0, 30])
    
    # Plot labeling
    ax.set_xlabel(param)
    ax.set_ylabel('probability density')
    ax.legend(frameon=False, bbox_to_anchor=(1, 1.025), loc='upper left')