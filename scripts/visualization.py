import cartopy
import cartopy.crs as ccrs
import numpy as np
import os
import pandas as pd
import scipy as sp
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

import utilities
import tc_analysis

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

def diverging_colormap_adjust(norm, cmap, additional=0):

    '''
    Adjust a diverging colormap to insert a white region in the middle for masking small values.
    The 'additional' argument represents the elements away from the center that are added for additional masking
    '''
    # Number of colors
    ncolors = 256
    # Get numerical array for colormap
    arr = matplotlib.cm.get_cmap(cmap)(np.linspace(0, 1, 256))
    # Get the center index for the numerical array (in other words, where the colormap diverges) and subtract by one 
    center = ncolors / 2
    # Offset
    additional = (additional + 1)*(center / (len(norm.boundaries) - 1))
    # If the number of colors is even, target the center entry and account for additional masking
    if center.is_integer():
        lower, upper = int(center - additional), int(center + additional)
        arr[lower:upper] = [1, 1, 1, 1]
    # Else, target the entries closest to center
    else:
        lower, upper = int(np.floor(center) - additional), int(np.ceil(center) + additional)
        arr[lower:(upper + 1)] = [1, 1, 1, 1], [1, 1, 1, 1]
    # Re-establish the colormap
    cmap = matplotlib.colors.ListedColormap(arr)

    return cmap

def norm_cmap(data, field, num_bounds=16, extrema=None, white_adjust=False, cmap_white_fraction=0.1, diagnostic=False):
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

    # Adjust vmin and vmax to nice numbers
    extremum = max([abs(vmin), abs(vmax)])
    order_of_magnitude = 10**(np.floor(np.log10(extremum)))
    
    num_bounds = int(4 * np.round(num_bounds/4))
    extremum = order_of_magnitude * np.round(extremum/order_of_magnitude)
        
    # Subdivide normalization into sequential (all one sign) or diverging (two signs, including 0) bins
    if vmin < 0 and vmax > 0:
        # Get larger of the two bounds
        bins = np.linspace(-extremum, extremum, num_bounds + 1)
    else:        
        # Adjust vmin and vmax to "nice" numbers
        delta = abs(vmin - vmax) # get difference to discern order of magnitude
        order_of_magnitude_delta = 10**(np.floor(np.log10(delta))) # get order of magnitude of difference
        # Redefine the extrema based on the new difference magnitude
        vmin = order_of_magnitude_delta * np.floor(vmin/order_of_magnitude_delta)
        vmax = order_of_magnitude_delta * np.ceil(vmax/order_of_magnitude_delta)
        bins = np.linspace(vmin, vmax, num_bounds+1)
        
        if diagnostic:
            print('Order of magnitude: extremum = {0}'.format(order_of_magnitude))
            print('Extrema: extremum = {0}, vmin = {1}, vmax = {2}'.format(extremum, vmin, vmax))
            print('Difference: {0} and order of magnitude: {1}'.format(delta, order_of_magnitude_delta))
            
    if diagnostic:
        print('Normalization bins: {0}'.format(bins))
        print('-----------------------------------------------------')
    
    # Assign to normalization dictionary
    norm = matplotlib.colors.BoundaryNorm(bins, 256)
    # Get colormap
    cmap = get_cmap(field, norm)
    # If the colormap is a "white-to-color" sequential colormap, apply an adjustment for a white basis
    sequential_cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds','YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 
                        'RdPu', 'BuPu', 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    # Perform white adjustment for colormaps (set lowest color to white for sequential, center to white for diverging)
    if white_adjust:
        cmap = cmap_white_adjust(cmap, levels=num_bounds) if cmap in sequential_cmaps else diverging_colormap_adjust(norm, cmap, additional=1)
        
    return norm, cmap

def field_properties(field):
    """
    Method to generate long name and units for a given field. Returns the field name and no units if field doesn't have a property entry.

    Args:
        field (str): GCM field name.

    Returns:
        long_name, units: str, str
    """

    long_name_addendum = ''
    if '_clr' in field:
        field = field.split('_clr')[0]
        long_name_addendum = ', clear-sky'
    
    properties = {'lhflx': {'long_name': 'latent heat flux', 'units': 'W m$^{-2}$'},
                  'shflx': {'long_name': 'sensible heat flux', 'units': 'W m$^{-2}$'},
                  'olr': {'long_name': 'outgoing longwave radiation', 'units': 'W m$^{-2}$'},
                  'netrad_toa': {'long_name': 'net radiation, TOA', 'units': 'W m$^{-2}$'},
                  'wind': {'long_name': '10m horizontal wind speed', 'units': 'm s$^{-1}$'},
                  'ucomp': {'long_name': 'zonal wind speed', 'units': 'm s$^{-1}$'},
                  'vcomp': {'long_name': 'meridional wind speed', 'units': 'm s$^{-1}$'},
                  'WVP': {'long_name': 'column-integrated water vapor\n', 'units': 'kg m$^{-2}$'},
                  'h': {'long_name': 'moist static energy', 'units': 'J kg$^{-1}$'},
                  'slp': {'long_name': 'sea level pressure', 'units': 'hPa'},
                  'omega': {'long_name': 'pressure velocity', 'units': 'Pa s$^{-1}$'},
                  'rh': {'long_name': 'relative humidity', 'units': '%'},
                  'vi_h': {'long_name': 'vertically-integrated MSE', 'units': 'J m$^{-2}$'},
                  'sphum': {'long_name': 'specific humidity', 'units': 'kg kg$^{-1}$'},
                  'precip': {'long_name': 'precipitation', 'units': 'mm d$^{-1}$'},
                  'evap': {'long_name': 'evaporation', 'units': 'mm d$^{-1}$'},
                  'h_anom': {'long_name': 'domainwise MSE anomaly\n', 'units': 'J kg$^{-1}$'},
                  'temp': {'long_name': 'temperature', 'units': 'K'},
                  'temp_anom': {'long_name': 'domainwise temperature anomaly\n', 'units': 'K'},
                  'sphum_anom': {'long_name': 'domainwise specific humidity anomaly\n', 'units': 'kg kg$^{-1}$'},
                  'wind_tangential': {'long_name': 'tangential velocity', 'units': 'm s$^{-1}$'},
                  'wind_radial': {'long_name': 'radial velocity', 'units': 'm s$^{-1}$'},
                  'tm': {'long_name': '300-500 hPa mean temperature\n', 'units': 'K'},
                  't_surf': {'long_name': 'surface temperature', 'units': 'K'},
                  'swfq': {'long_name': 'SWISHE frequency', 'units': '% of timestamps'},
                  'p-e': {'long_name': 'surface moisture flux', 'units': 'mm d$^{-1}$'},
                  'cld_amt': {'long_name': 'cloud amount', 'units': '%'},
                  'max_wind': {'long_name': 'maximum 10 meter wind speed', 'units': 'm s$^{-1}$'},
                  'min_slp': {'long_name': 'minimum sea level pressure', 'units': 'hPa'},
                  'theta_e': {'long_name': '$\\theta_e$', 'units': 'K'},
                  'olr_clr': {'long_name': '$\\theta_e$', 'units': 'K'},
                  'theta_e': {'long_name': '$\\theta_e$', 'units': 'K'},
                  'mld': {'long_name': 'mixed-layer depth', 'units': 'm'}}

    properties[field]['long_name'] = properties[field]['long_name'] + long_name_addendum
    
    if field in properties.keys():
        return properties[field]['long_name'], properties[field]['units']
    else:
        return field, '-'

def line_plot_interpolator(X, Y, interpolation_factor=10):

    # Remove nans and infs (note that the value array will have missing values)
    X, Y = X[np.isfinite(Y)], Y[np.isfinite(Y)]
    # Get new basis array with 'interpolation_factor' times as many points
    X_interp = np.linspace(X.min(), X.max(), len(X)*interpolation_factor)
    # Interpolate values using a BSpline
    spl = sp.interpolate.make_interp_spline(X, Y, k=3)
    Y_interp = spl(X_interp)

    return X_interp, Y_interp

def contour_tracing_levels(im, norm, diverging_mask_width=1, diagnostic=False):
    
    # Perform adjustment to remove white-masked values from the contours
    if np.sign(min(norm.boundaries)) != np.sign(max(norm.boundaries)):
        # Find index in levels where the value is 0
        center_index = np.nanmin(np.argwhere(im.levels == 0))
        # Redefine levels based on the zero-index and the width of the diverging white mask
        # This allows the new contour levels to identically trace the filled contour levels
        im_levels = np.concatenate([im.levels[0:(center_index - diverging_mask_width + 1)],
                                    im.levels[(center_index + diverging_mask_width):]])
        # Get normalization half-width
        norm_halfwidth = int(len(norm.boundaries) // 2)
        norm_mask_vmin = norm.boundaries[norm_halfwidth - diverging_mask_width]
        norm_mask_vmax = norm.boundaries[norm_halfwidth + diverging_mask_width]
        
        # Remove levels from the new contour levels if they fall within the normalization mask extrema
        im_levels = np.where((im_levels <= norm_mask_vmin) | (im_levels >= norm_mask_vmax), im_levels, np.nan)
    else:
        im_levels = norm.boundaries
        
    return im_levels

def subplot_statistics(arr, fig, ax, label_position_x=0.96, label_position_y=0.05, fontsize=7, text_color=None, text_offset=False):

    # Get relevant data statistics
    mean, std = np.nanmean(arr), np.nanstd(arr)
    vmin, vmax = np.nanmin(arr), np.nanmax(arr)
    # Get sample size from attributes, if available
    sample_size = arr.attrs['Sample size'] if 'Sample size' in arr.attrs.keys() else None
    
    # Get order of magnitude from the mean
    order_of_magnitude = np.floor(np.log10(mean))
    # Define string formatter
    formatter = '.2e' if abs(order_of_magnitude) > 3 else '.2f'
    
    # Define the string
    sample_size_str = 'N = {0}; '.format(sample_size) if sample_size else ''
    statistics_str = '{5}mean: {1:{0}} $\pm$ {2:{0}}\nmin: {3:{0}}, max: {4:{0}}'.format(formatter, mean, std, vmin, vmax, sample_size_str)
    # Adjust text alignment based on position such that it aligns towards a corner
    ha = 'left' if label_position_x < 0.5 else 'right'
    va = 'bottom' if label_position_y < 0.5 else 'top'
    # Define text color
    text_color = 'k' if not text_color else text_color

    # If a text offset is requested, do the following:
    # 1. Plot a dummy text box with transparent text.
    # 2. Get the height
    # 3. Offset the vertical position by the text height
    if text_offset:
        # 1. Plot a dummy text box with transparent text.
        ann = ax.annotate(statistics_str, xy=(label_position_x, label_position_y), xycoords='axes fraction',
                          ha=ha, va=va, fontsize=fontsize, color=text_color, alpha=0)
        # 2. Get the height relative to the axes fraction.
        transform = ax.transAxes.inverted()
        bounding_box = ann.get_window_extent(renderer=fig.canvas.get_renderer()).transformed(transform)
        text_height = bounding_box.y1 - bounding_box.y0
        direction_factor = -1 if va == 'top' else 1 # adjust direction of translation based on alignment
        spacing_factor = 1.25 # multiplying factor for vertical padding
        # 2.  Offset the vertical position by the text height
        ann = ax.annotate(statistics_str, xy=(label_position_x, label_position_y + spacing_factor*direction_factor*text_height), 
                          xycoords='axes fraction', ha=ha, va=va, fontsize=fontsize, color=text_color)
        
    # Else, just plot normally
    else:
        ann = ax.annotate(statistics_str, xy=(label_position_x, label_position_y), xycoords='axes fraction',
                          ha=ha, va=va, fontsize=fontsize, color=text_color)
    # White stroke on the text
    ann.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1, foreground='white'),
                          matplotlib.patheffects.Normal()])

    return ax

### Begin compositing support methods

def composite_tick_formatting(ax, plot_type, norm, nrows, ncols, row=None, col=None):

    # Define axis limits
    if plot_type == 'planar':
        min_x, max_x = [-10, 10]
        min_y, max_y = [-10, 10]
    elif plot_type == 'azimuthal_1D':
        min_x, max_x = [0, 7.5]
        min_y, max_y = min(norm.boundaries), max(norm.boundaries)
    else:
        min_x, max_x = [0, 10]
        min_y, max_y = ax.get_ylim()

    # Define tick intervals
    tick_interval_x = 1
    tick_interval_y = 1 if plot_type == 'planar' else 100

    # Define tick locations
    tick_locator_x = np.arange(min_x, max_x + tick_interval_x, tick_interval_x)
    tick_locator_y = np.arange(min_y, max_y + tick_interval_y, tick_interval_y)
    
    # Set x-axis ticks
    ax.set_xlim([min_x, max_x])
    ticks_x = ax.set_xticks(tick_locator_x)
    ax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    
    # Set y-axis ticks
    if plot_type in ['planar', 'azimuthal_2D']:
        ax.set_ylim([min_y, max_y])
        ticks_y = ax.set_yticks(tick_locator_y)
        ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    else:
        ax.set_ylim([min_y, max_y])
        ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    # Reverse y-axis direction for a 2D azimuthal plot
    if plot_type == 'azimuthal_2D':
        ax.set_ylim(ax.get_ylim()[::-1])

    ''' Control tick labeling locations. '''
    
    # Only plot tick labels on the top and right of a subplot
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()
        
    # Hide tick labels if the subplot is not an edge case
    # Only plot x-tick labels on the top row
    if row != nrows - 1:
        ax.set_xticklabels([])
    # Only plot y-tick labels on the top row
    if plot_type in ['planar', 'azimuthal_2D'] and col != 0:
        ax.set_yticklabels([])

    return ax

def composite_colorbar(fig, axes, field, experiments, norms, cmaps):
    
    field_name, field_units = field_properties(field)
    colorbar_offset_y = 1.1 # in units of axes fraction
    colorbar_height = 0.075 # in units of axes fraction
    
    for experiment_index, experiment in enumerate(experiments):

        ax = axes['{0}, {1}'.format(0, experiment_index)]
        norm, cmap = norms[experiment], cmaps[experiment]

        cax = ax.inset_axes([0, colorbar_offset_y + colorbar_height, 1, colorbar_height])
        colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')

        cax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(norm.boundaries, 5))
        cax.xaxis.set_label_position('top')
        cax.xaxis.set_ticks_position('top')

        # Change string formatting for tick notation (see the ScalarFormatterForceFormat() class)
        fmt = ScalarFormatterForceFormat()
        fmt.set_powerlimits((-2, 3))
        cax.xaxis.set_major_formatter(fmt)
        
        if experiment_index == int(np.floor(len(experiments)/2)):
            colorbar_label = '{0} [{1}]'.format(field_name, field_units)
            colorbar.set_label(colorbar_label, ha='center', labelpad=10)

        fig.canvas.draw()
        offset = cax.xaxis.get_major_formatter().get_offset()
        if len(offset) > 0:
            cax.xaxis.get_offset_text().set_visible(False)
            fontsize = cax.get_xticklabels()[0].get_fontsize()
            cax.annotate(offset, (1, -0.5), xycoords='axes fraction', ha='right', va='top', fontsize=fontsize - 1)

    return colorbar

def composite_normalization_colormap(dataset, experiments_in=['CTL1990s', 'CTL1990s_swishe'], field='WVP', 
                                     pressure_limits=None, contour_levels=16, diverging_mask_width=1, extremum_type='absolute', diagnostic=False):
    
    if pressure_limits:
        pressure_bottom, pressure_top = max(pressure_limits), min(pressure_limits)
    else:
        pressure_bottom, pressure_top = 1000, 100
        
    model_names = dataset.keys()
    run_control, run_experiment = experiments_in
    experiments = [run_control, run_experiment]
    difference_experiment = '{1}-{0}'.format(run_control, run_experiment)
    
    extrema = {'vmin': np.nan, 'vmax': np.nan}
    difference_extrema = {'vmin': np.nan, 'vmax': np.nan}
    
    for model in dataset.keys():
        for experiment in experiments:
            # Define iterand dataset for convenience
            arr = dataset[model][experiment]
            # Slice vertically if the data has a vertical component
            if 'pfull' in arr.dims:
                if diagnostic:
                    print('[composite_normalization_colormap()] Pressure selection being made.')
                arr = arr.sel(pfull=slice(pressure_top, pressure_bottom))
            # Get iterand dataset extrema
            if extremum_type == 'absolute':
                vmin, vmax = arr.min(), arr.max()
            elif extremum_type == 'std':
                num_deviations = 2
                mean, std = np.nanmean(arr), np.nanstd(arr)
                vmin, vmax = mean - num_deviations * std, mean + num_deviations * std
            # Update extrema if the iterand extrema superseded it/them
            extrema['vmin'] = vmin.item() if (np.isnan(extrema['vmin']) or vmin < extrema['vmin']) else extrema['vmin']
            extrema['vmax'] = vmax.item() if (np.isnan(extrema['vmax']) or vmax > extrema['vmax']) else extrema['vmax']
            if diagnostic:
                print('Model: {0:6s}; experiment: {1:15s}; field: {2} | minimum value: {3:.2f}; maximum value: {4:.2f}'.format(model, experiment, field, vmin, vmax))
    
        # Define iterand dataset for convenience
        arr = dataset[model][difference_experiment]
        # Slice vertically if the data has a vertical component
        if 'pfull' in arr:
            arr = arr.sel(pfull=slice(pressure_top, pressure_bottom), method='nearest')
        # Get iterand dataset extrema - use a standard deviation exceedance approach in case of anomalous values
        mean, std = np.nanmean(arr), np.nanstd(arr)
        num_deviations = 4
        if std > num_deviations * mean:
            vmin, vmax = mean - num_deviations * std, mean + num_deviations * std
        else:
            vmin, vmax = arr.min(), arr.max()
        # Update extrema if the iterand extrema superseded it/them
        difference_extrema['vmin'] = vmin.item() if (np.isnan(difference_extrema['vmin']) or vmin < difference_extrema['vmin']) else difference_extrema['vmin']
        difference_extrema['vmax'] = vmax.item() if (np.isnan(difference_extrema['vmax']) or vmax > difference_extrema['vmax']) else difference_extrema['vmax']
        
        if diagnostic:
            print('[difference] Model: {0:6s}; experiment: {1:15s}; field: {2} | minimum value: {3:.2f}; maximum value: {4:.2f}'.format(model, experiment, field, vmin, vmax))
    
    experiments.append(difference_experiment)
    norms, cmaps = {}, {}
    for model in dataset.keys():
        for experiment in experiments:
            extremum = extrema if experiment in experiments_in else difference_extrema
            norms[experiment], cmaps[experiment] = norm_cmap(0, field=field, extrema=(extremum['vmin'], extremum['vmax']), num_bounds=contour_levels)

            if diagnostic:
                print('Model: {0:6s}; experiment: {1:15s}; field: {2} | normalization levels: {3}'.format(model, experiment, field, norms[experiment].boundaries))
            
            if min(norms[experiment].boundaries) < 0 and max(norms[experiment].boundaries) > 0:
                cmaps[experiment] = diverging_colormap_adjust(norms[experiment], cmaps[experiment], additional=diverging_mask_width)

    return norms, cmaps

def composite_field_overlay(dataset, field, base_im):
    # Get normalization and colormap data for the iterand data
    norm, cmap = composite_normalization_colormap(dataset, field=field, contour_levels=base_im.levels)
    
def composite_plot(dataset, models, experiments, field, plot_type, contour_levels=16, format='contour', overlay=False, parent_fig=None, grid_overlay=False, dpi=96, diagnostic=False):

    # Get xarray-compatible field name for pressure-level slices
    field = field.split('hPa')[0] if 'hPa' in field else field
    # Define the width of the white masking for divergent data
    diverging_mask_width = 1
    # Get normalization and colormap data for the iterand data
    norms, cmaps = composite_normalization_colormap(dataset, field=field, contour_levels=contour_levels, diverging_mask_width=diverging_mask_width)
    
    # Set up basis vector names for plots
    coordinates = {'planar': {'x': 'grid_xt', 'y': 'grid_yt'},
                   'azimuthal_1D': {'x': 'radius'},
                   'azimuthal_2D': {'x': 'radius', 'y': 'pfull'}}

    # Define figure structure parameters
    width_factor, height_factor = 3, 3
    nrows, ncols = len(models), len(experiments)
    if plot_type == 'azimuthal_1D':
        ncols -= 1
        width_factor -= 0
        height_factor -= 1.5
    aspect_ratio = height_factor/width_factor
    
    if parent_fig and overlay:
        fig, gs, axes = parent_fig
    else:
        fig, gs = [plt.figure(figsize=(width_factor*ncols, height_factor*nrows), dpi=dpi),
                   matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols)]
        axes = {}
    
    for model_index, model in enumerate(models):
        for experiment_index, experiment in enumerate(experiments):

            # Get iteration-specific plot properties for line plotting
            cycle_properties = cycler(experiment_index)
            linecolor = cycle_properties['c']
            
            # Boolean to determine whether the row/column combination has been used yet.
            # This is needed to distinguish 1D from 2D plotting.
            preplotted, text_color = False, None
            
            # Determine column index for grid structure. 
            # Set to equal experiment_index if this is 2D, otherwise make it dependent on experiment type.
            if plot_type in ['planar', 'azimuthal_2D']:
                column_index = experiment_index 
            else:
                # This assumes that the difference experiment has a '-' in it
                column_index = 0 if '-' not in experiment else 1
                
            # Define identifying name for the specific subplot
            subplot_name = '{0}, {1}'.format(model_index, column_index)
            
            # If the subplot exists, use it. If not, define it
            if subplot_name not in axes.keys() and not overlay:
                ax = fig.add_subplot(gs[model_index, column_index])
            else:
                ax = axes[subplot_name]
                preplotted = True

            # Define a shorthand for the iterand dataset
            arr = dataset[model][experiment]
            # Get normalization and colormap
            norm, cmap = norms[experiment], cmaps[experiment]
            # If the data is sequential, white-out the near-zero values
            if (np.sign(min(norm.boundaries)) == np.sign(max(norm.boundaries))) or (min(norm.boundaries) == 0 or max(norm.boundaries) == 0):
                cmap = cmap_white_adjust(cmap)

            # Handle 1D plotting
            experiment_label = experiment if model_index == 0 else ''
            # Get the x-basis vector
            coordinate_x = coordinates[plot_type]['x'] 
            # Handle 2D plotting
            if plot_type in ['planar', 'azimuthal_2D']:
                coordinate_y = coordinates[plot_type]['y']
                if format == 'contour' and not overlay:
                    im = ax.contourf(arr[coordinate_x], arr[coordinate_y], arr, norm=norm, cmap=cmap, levels=contour_levels)
                    tracer_levels = contour_tracing_levels(im, norm, diverging_mask_width=diverging_mask_width, diagnostic=diagnostic)
                    # im_trace = ax.contour(arr[coordinate_x], arr[coordinate_y], arr, levels=tracer_levels, colors=['k'], alpha=0.25, linewidths=0.5)
                elif overlay:
                    im = ax.contourf(arr[coordinate_x], arr[coordinate_y], arr, norm=norm, cmap=cmap, alpha=0)
                    tracer_levels = contour_tracing_levels(im, norm, diverging_mask_width=diverging_mask_width, diagnostic=diagnostic)
                    im = ax.contour(arr[coordinate_x], arr[coordinate_y], arr, norm=norm, levels=tracer_levels, colors=['k'], alpha=0.25, linewidths=1)
                    print('[visualization.py(), composite_plot()] {0} {1} overlay contours: {2}'.format(model, experiment, im.levels))
                else:
                    im = ax.pcolormesh(arr[coordinate_x], arr[coordinate_y], arr, norm=norm, cmap=cmap)
                
            elif plot_type == 'azimuthal_1D':
                radius_interpolated, arr_interpolated = line_plot_interpolator(arr.radius, arr)
                im = ax.plot(radius_interpolated, arr_interpolated, linewidth=3, color=linecolor, label=experiment_label)
                # Plot reference 0 line if the difference has different signs
                if min(norm.boundaries) < 0 and max(norm.boundaries) > 0:
                    ax.axhline(0, c='k', lw=1, zorder=0)
                # Collect the text color for annotation
                text_color = im[0].get_color()
            
            # Tick formatting. Note the boolean use to avoid overwriting.
            if not preplotted:
                ax = composite_tick_formatting(ax, plot_type, norm, nrows, ncols, row=model_index, col=experiment_index)

            # Plot centerlines and set aspects equal if the requested plot is planar
            if plot_type == 'planar':
                ax.set_aspect('equal')
                ax.axhline(0, color='k', alpha=0.25, linestyle='--', linewidth=0.5)
                ax.axvline(0, color='k', alpha=0.25, linestyle='--', linewidth=0.5)            
            # Create grid
            if grid_overlay:
                grid = ax.grid(color='k', alpha=0.25, linestyle='--', linewidth=0.5)
            # Add inline statistics
            if not overlay:
                ax = subplot_statistics(arr, fig, ax, text_color=text_color, text_offset=preplotted)
            # Add subplot identifier
            if not preplotted:
                subplot_identifier = ax.annotate('({0})'.format(chr(ord('a') + len(axes))), 
                                                xy=(0.04*aspect_ratio, 0.96), xycoords='axes fraction', 
                                                fontsize=8, color='k', va='top', ha='left')
                subplot_identifier.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1, foreground='white'),
                                                    matplotlib.patheffects.Normal()])
            
            # Append subplot to collection
            axes[subplot_name] = ax

    ''' Construct colorbars. '''
    if plot_type in ['planar', 'azimuthal_2D'] and not overlay:
        colorbar = composite_colorbar(fig, axes, field, experiments, norms, cmaps)
    else:
        fig.tight_layout()
        fig.legend(ncols=len(experiments), bbox_to_anchor=[0.5, 1.05], loc='center', frameon=False)
        
    
    
    return fig, gs, axes

### End compositing support methods

def density_grid(data, model_names=['HIRAM', 'AM2.5', 'FLOR'], experiments=['control', 'swishe'], difference_order='rtl', bin_resolution=5, dpi=96):
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
    run_control, run_experiment = experiments
    # Get minimum and maximum years based on the control experiment
    year_min, year_max = density[model_names[0]][run_control]['time_min'].unique().year.item(), density[model_names[0]][run_control]['time_max'].unique().year.item()
    
    ''' Begin plotting. '''
    # Note: number of rows is dictated by number of models + 1 (for colorbar) 
    # Note: number of columns is dictated by number of experiments + 1 (for difference plot)
    nrows, ncols = len(model_names) + 1, 3
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
    run_difference = '{0} - {1}'.format(run_control, run_experiment)
    # Log extrema across all model difference plots for future normalization
    difference_extrema = {model_name: None for model_name in model_names}
    # Define direction for difference calculation
    difference_factor = -1 if difference_order == 'rtl' else 'ltr'
    # Iterate over models to calculate differences for each
    for model_name in densities.keys():
        densities[model_name][run_difference] = difference_factor*(densities[model_name][run_control] - densities[model_name][run_experiment])
        difference_extrema[model_name] = max([abs(np.nanmin(densities[model_name][run_difference])), abs(np.nanmax(densities[model_name][run_difference]))])
    # Get extreme values for difference runs and find largest magnitude for normalization
    vmax_difference = max([np.nanmax(abs(v)) for k, v in difference_extrema.items()]) # maximum value across 'run_control' and 'run_experiment' for all models
    # Assign normalization and colormap values
    for experiment in densities[model_name].keys():
        if '-' in experiment: 
            # This is old code to be deleted on next commit to git - GR (9/16)
            # norm_bins = [0.5*round(value/0.5) for  value in np.linspace(-vmax_difference, vmax_difference, bounds)]
            # norms[experiment] = matplotlib.colors.BoundaryNorm(norm_bins, 256)
            # cmaps[experiment] = 'bwr'
            norm, cmap = norm_cmap([0], 'density', extrema=(-vmax_difference, vmax_difference))
            cmap = diverging_colormap_adjust(norm, cmap, additional=1)
            norms[experiment] = norm
            cmaps[experiment] = cmap
        else:
            norm, cmap = norm_cmap([0], 'density', extrema=(0, vmax_runs))
            norms[experiment] = norm
            cmaps[experiment] = 'Reds'
            cmaps[experiment] = cmap_white_adjust(cmaps[experiment], len(norm.boundaries))
    # Pass 2: Iterate through density dictionary to plot each model and each experiment's density data.
    for model_index, model_name in enumerate(densities.keys()):
        for experiment_index, experiment in enumerate(densities[model_name].keys()):
            # Boolean to control whether or not model and years print. Turn off for the difference plot.
            subplot_title = False if "-" in experiment else True 
            # Set up the figure structure for the iterand subplot
            fig, grid, ax, proj = basemap(fig, grid, model_name, experiment, (year_min, year_max), 
                                          row_num=model_index, col_num=experiment_index, land=True, xlabel=False, ylabel=False, 
                                          subplot_title=subplot_title, label_fontsize=10)

            # Plot the data
            im = ax.pcolormesh(x, y, densities[model_name][experiment], norm=norms[experiment], cmap=cmaps[experiment], transform=proj_ref)
            
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

def pdf(data, models=['AM2.5', 'FLOR', 'HIRAM'], param='center_lat', num_bins=60, mode='pdf', dpi=144):
    
    fig, ax = plt.subplots(figsize=(4, 3), dpi=dpi)
    
    # Define list of colors (one color per model) and line styles (one linestyle per experiment)
    colors = ['b', 'g', 'm']
    linestyles = ['-', '--']
    
    # Collect data for output counts
    counts = {model: {} for model in models}

    # Obtain maximum values from histogram distributions
    vmax = 0
    # Save label names
    legend_labels = []
    # Iterate over all models and experiments
    for model_index, model in enumerate(models):
        for experiment_index, experiment in enumerate(data[model].keys()):
            # Get unique values for the parameter passed
            out = data[model][experiment]['unique'][param].dropna()
            # Filter out nans and infs
            out = out.loc[np.isfinite(out.values)]
            
            legend_labels.append('{0}, {1}'.format(model, experiment))
            bw_method = 0.1 if param == 'center_lat' else 0.25
            pd.DataFrame(out).plot.kde(bw_method=bw_method, lw=2.5, color=colors[model_index], linestyle=linestyles[experiment_index],
                                       legend=False, ax=ax, zorder=5)
            if param not in ['center_lat']:
                ax.axvline(np.nanmedian(out), lw=1, alpha=0.375, color=colors[model_index], linestyle=linestyles[experiment_index], zorder=1, label='_nolegend_')
            
            # Get bin values and bin edges, don't plot
            n, bins, _ = ax.hist(out, bins=num_bins, histtype=u'step', density=True, lw=0, label='_nolegend_') 
            # Incorporate mutiplicative factor - if mode == 'full_count', factor = length of data; else, factor = 1. The former gives you a full count.
            factor = len(out) if mode == 'full_count' else 1
            # Plot the distribution
            # ax.plot(bins[:-1], n*factor, lw=2.5, color=colors[model_index], linestyle=linestyles[experiment_index], 
            #         label='{0}, {1}'.format(model, experiment))
            # Check to see if number of entries for a given bin exceeds the maximum (vmax)
            vmax = np.nanmax(n*factor) if np.nanmax(n*factor) > vmax else vmax
            # Add data length to counts
            counts[model][experiment] = {'bins': bins, 'counts': n*len(out), 'median': np.median(out), 'std': np.nanstd(out)}
            
    # Adjust axis limits: handle limit differently for center_lat to allow for inline legend
    if param == 'max_wind':
        ax.set_xlim([10, 50])
    elif param == 'min_slp':
        ax.set_xlim([900, 1020])
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
        ax.legend(legend_labels, frameon=False, loc='upper left', prop={'size': 8})
    else:
        ax.legend(legend_labels, frameon=False, loc='best', prop={'size': 8})
        
    # Print median and standard deviations for all model-experiment combinations
    for model in counts.keys():
        for experiment in counts[model].keys():
            print('{0} - {1}, {2}: median: {3:.2f}, std.dev: {4:.2f}'.format(param, model, experiment,
                                                                             counts[model][experiment]['median'], counts[model][experiment]['std']))
        
    return counts

def intensity_distribution_plot(data, intensity_metric='max_wind', intensity_bin=None,
                                bin_width=4, rolling_mean_window=6, interpolation_factor=4, 
                                plot_histogram=False, dpi=96, print_statistics=False):
    
    def rolling_mean(a, window=3):
        ret = np.cumsum(a, dtype=float)
        ret[window:] = ret[window:] - ret[:-window]
        return ret[window - 1:] / window

    # Get model names from the input data
    models = data.keys()
    # Collect iterand extrema
    extrema = {'x': {'min': np.nan, 'max': np.nan},
               'y': {'min': np.nan, 'max': np.nan}}
    # Collect processed data
    output = {}
    
    # Define figure as a function of rows (different models)
    fig, axes = plt.subplots(figsize=(4, 2*len(models)), dpi=dpi, nrows=len(models))
    # Iterate over each model
    for model_index, model in enumerate(models):
        # Initialize dictionary at model level to collect processed data
        output[model] = {}
        # Get experiment names from the input data
        experiments = data[model].keys()
        # Define the model-specific axis
        ax = fig.axes[model_index]
        # Iterate over each experiment
        for experiment_index, experiment in enumerate(experiments):
            # Get iterand color scheme
            color, linestyle = cycler(experiment_index)['c'], cycler(experiment_index)['ls']
            # Get intensity data for a given model and experiment configuration for the chosen metric
            if intensity_bin:
                # Get all storm data
                entries = data[model][experiment]['raw']
                # Filter by intensity bin
                intensity_bin_entries = entries.loc[entries['intensity_bin'] == intensity_bin]
                # Initialize list to hold strongest storm entry for each storm, to be concatenated lated
                entry_list = []
                # For each unique storm, get the strongest entry corresponding to the given intensity bin
                for unique_storm_id in intensity_bin_entries['storm_id'].unique():
                    # Get unique storm entry
                    unique_storm = intensity_bin_entries.loc[intensity_bin_entries['storm_id'] == unique_storm_id]
                    # Get metric-dependent maximum value
                    
                    if intensity_metric == 'min_slp':
                        entry = unique_storm.loc[unique_storm[intensity_metric] == unique_storm[intensity_metric].min()][intensity_metric]
                    else:
                        entry = unique_storm.loc[unique_storm[intensity_metric] == unique_storm[intensity_metric].max()][intensity_metric]
                    entry_list.append(entry)
                # Concatenate entries
                x = pd.concat(entry_list)
            else:
                x = data[model][experiment]['unique'][intensity_metric]
            # Collect processed data for the iterand configuration
            output[model][experiment] = x
            # Define extrema for the intensity distribution. Make them integers for cleaner processing.
            experiment_vmin_x, experiment_vmax_x = np.round(np.nanmin(x)), np.round(np.nanmax(x))
            # Round to nearest bin width for cleaner divisions. 
            # Note the decrement and increment for floor and ceiling, respectively, to foresee domain reductions due to rolling mean.
            experiment_vmin_x = bin_width * np.floor(experiment_vmin_x/bin_width) - (rolling_mean_window // 2) * bin_width # set to floor to capture minima
            experiment_vmax_x = bin_width * np.ceil(experiment_vmax_x/bin_width) + (rolling_mean_window // 2) * bin_width # set to ceiling to capture maxima
            # Define bins for histogram definition
            bins = np.arange(experiment_vmin_x, experiment_vmax_x + bin_width, bin_width)
            # Plot histogram and extract outputs for future interpolation
            histogram_alpha = 0.25 if plot_histogram else 0
            histogram_output = ax.hist(x, density=True, bins=bins, histtype='stepfilled',
                                       facecolor=color, align='left', alpha=histogram_alpha, zorder=1)
            densities, bin_locations = histogram_output[0], histogram_output[1]
            
            ''' Derive an interpolated pseudo-probability density function. '''
            # Derive the interpolation number (length of bins multiplied by some factor)
            interpolation_number = len(bins) * interpolation_factor
            # Derive the interpolation bin locations
            interpolated_bins = np.arange(experiment_vmin_x, experiment_vmax_x + bin_width / interpolation_factor, bin_width / interpolation_factor)
            # Interpolate the histogram-derived densities to the interpolated bin locations
            interpolated_densities = np.interp(interpolated_bins, bin_locations[:-1], densities)
            # Redefine bins for the rolling mean
            interpolated_bins_rolling = interpolated_bins[(rolling_mean_window//2):-(rolling_mean_window//2 - 1)]
            # Apply a rolling mean (centered) to the densities
            interpolated_densities_rolling = rolling_mean(interpolated_densities, rolling_mean_window)
            # Plot the interpolated data. Shift the bin locations by half the bin width
            ax.plot(interpolated_bins_rolling, interpolated_densities_rolling,
                    color=color, linewidth=3, zorder=10, label=experiment)
            
            # Get median value
            median = np.median(x)
            median_density = interpolated_densities_rolling[np.argmin(np.abs(interpolated_bins_rolling - median))]
            # Plot the median value
            ax.plot([median, median], [0, median_density], lw=1, ls=':', color=color)
        
            # Check experiment extrema against overall extrema
            extrema['x']['min'] = experiment_vmin_x if (np.isnan(extrema['x']['min']) or experiment_vmin_x < extrema['x']['min']) else extrema['x']['min']
            extrema['x']['max'] = experiment_vmax_x if (np.isnan(extrema['x']['max']) or experiment_vmax_x > extrema['x']['max']) else extrema['x']['max']
            extrema['y']['min'] = min(densities) if (np.isnan(extrema['y']['min']) or min(densities) < extrema['y']['min']) else extrema['y']['min']
            extrema['y']['max'] = max(densities) if (np.isnan(extrema['y']['max']) or max(densities) > extrema['y']['max']) else extrema['y']['max']
    
    # Override extrema for winds.
    if intensity_metric == 'max_wind':
        extrema['x']['min'] = 0
        
    ''' Apply extrema to define axis limits. '''
    # Iterate over each model
    for model_index, model in enumerate(models):
        # Define the model-specific axis
        ax = fig.axes[model_index]
        # Add annotation to show model name
        ax.annotate(model, (0.03, 0.94), xycoords='axes fraction', fontsize=10, va='top')
        # Iterate over each experiment
        for experiment_index, experiment in enumerate(experiments):
            # Set minor tick locations automatically
            ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
            # Print statistics if chosen
            if print_statistics:
                sample_count = len(data[model][experiment]['unique'][intensity_metric])
                ax.annotate('N({0}) = {1}'.format(experiment, sample_count), (0.98, 0.95 - experiment_index*0.08), xycoords='axes fraction', 
                            fontsize=8, ha='right', va='top')
            # Ensure that subplot limits are defined by the data extrema
            ax.set_xlim([extrema['x']['min'], extrema['x']['max']])
            ax.set_ylim([extrema['y']['min'], extrema['y']['max']])
            # Only plot the legend for the top subplot
            if model_index == 0:
                ax.legend(frameon=False, ncols=len(experiments), loc='upper center', bbox_to_anchor=(0.5, 1.3))
            if model_index == len(models) - 1:
                long_name, units = field_properties(intensity_metric)
                ax.set_xlabel("{0} [{1}]".format(long_name, units))
                
    return output

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

def basemap(fig, gs, model_name, experiment, year_range=None, row_num=0, col_num=0, extent=[0, 359, -60, 60], land=False,
            xlabel=True, ylabel=True, subplot_title=True, label_fontsize=11):
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
    if subplot_title:
        subplot_title_model = ax.annotate('{0}{1}'.format(model_name, year_range_str),
                                          (0, title_y), va='baseline', ha='left', xycoords='axes fraction', fontsize=label_fontsize)
    # Set right-hand side to be {experiment name}
    subplot_title_experiment = ax.annotate('{0}'.format(experiment), (1, title_y), va='baseline', ha='right', 
                                           xycoords='axes fraction', fontsize=label_fontsize)

    # Initialize gridline parameters
    gridline_x_step, gridline_y_step = 60, 20
    gridline_minor_step = 10
    # Define major ticks
    gridline_x_major, gridline_y_major = [np.arange(extent[0] - longitude_offset, 
                                                    extent[1] - longitude_offset + gridline_x_step, 
                                                    gridline_x_step), 
                                          np.arange(extent[2], 
                                                    extent[3] + gridline_y_step, 
                                                    gridline_y_step)]
    gridline_x_minor, gridline_y_minor = [np.arange(extent[0] - longitude_offset, 
                                                    extent[1] - longitude_offset + gridline_minor_step, 
                                                    gridline_minor_step), 
                                          np.arange(extent[2], 
                                                    extent[3] + gridline_minor_step, 
                                                    gridline_minor_step)]
    
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
    
    basins = {'global': [1, 360, -60, 59],
             'NI': [40, 100, 0, 30],
             'SI': [40, 120, -40, 0],
             'WP': [100, 180, 0, 40],
             'SP': [120, 260, -40, 0],
             'EP': [180, 260, 0, 40],
             'NA': [260, 360, 0, 40],
             'SA': [300, 360, -40, 0],}
    
    # Generate masks to filter by basin
    basin_masks = {}
    grid = xr.open_dataset('/tigress/GEOCLIM/gr7610/tools/AM2.5_atmos_area.nc')
    for i, basin_data in enumerate(basins.items()):
        basin, basin_extent = basin_data
        mask = np.full(shape=(len(grid.grid_yt), len(grid.grid_xt)), fill_value=0)

        for j, y in enumerate(grid.grid_yt.values):
            for i, x in enumerate(grid.grid_xt.values):
                # Rectangular basins
                if (x >= basin_extent[0]) and (x < basin_extent[1]) and (y >= basin_extent[2]) and (y < basin_extent[3]):
                    mask[j, i] = 1
                # Special cases - handle the Central American isthmus
                if basin in ['EP', 'NA']:
                    x0, x1, y0, y1 = [260, 295, 0, 22]
                    m = (y0 - y1)/(x1 - x0)
                    f_y = m*(grid.grid_xt.values - x1) + y0
                    if (x >= x0) and (x < x1) and (y >= y0) and (y < f_y[i]):
                        mask[j, i] = 1 if basin == 'EP' else 0
            basin_masks[basin] = xr.DataArray(data=mask,
                                              dims=grid.dims,
                                              coords=grid.coords)

    if visualize:
        proj, proj_ref = ccrs.PlateCarree(central_longitude=180), ccrs.PlateCarree()
        fig, gs = plt.figure(figsize=(8, 3), dpi=300), matplotlib.gridspec.GridSpec(nrows=1, ncols=1)

        fig, gs, ax, proj_ref = basemap(fig, gs, 'Basins', '', row_num=0, col_num=0, extent=[0, 359, -60, 60], land=True)
        
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        
        ax.coastlines()
        for i, basin_data in enumerate(basins.items()):
            basin, basin_extent = basin_data
            text_offset_horizontal = 40 if basin == 'NA' else 0
            if basin in ['EP', 'NA']:
                if basin == 'EP':
                    xy = [(180, 40), (260, 40), (260, 22), (295, 0), (180, 0)]
                    patch = matplotlib.patches.Polygon(xy, closed=True, fc='none', lw=2, ec=colors[i], transform=proj_ref, zorder=i+10)
                else:
                    xy = [(260, 40), (360, 40), (360, 0), (295, 0), (260, 22)]
                    patch = matplotlib.patches.Polygon(xy, closed=True, fc='none', lw=2, ec=colors[i], transform=proj_ref, zorder=i+10)
            else:
                patch = matplotlib.patches.Rectangle(xy=(basin_extent[0], basin_extent[2]),
                                                    width=(basin_extent[1] - basin_extent[0]), height=(basin_extent[3] - basin_extent[2]), 
                                                    fc='none', lw=2, ec=colors[i], transform=proj_ref, zorder=i+10)
            ax.annotate(text=basin, xy=(basin_extent[0], basin_extent[2]), xytext=(5 + text_offset_horizontal, 6), 
                        textcoords='offset points', transform=proj_ref, bbox=dict(facecolor=(1, 1, 1, 0.9), ec='None'))
            ax.add_patch(patch)
        ax.set_extent([0, 359, -60, 60])
        plt.show()
        
    return basins, basin_masks

def tc_activity(track_data, model_name, experiments=['control', 'swishe'], storm_types=['TS', 'C15w'], 
                year_range=None, savefig=False):
    
    ''' Function to generate monthly climatologies of TC activity, subdivided per basin. See the script snippet below for recommended function use. '''
    
    ####### Begin recommended run snippet. ##########
    monthly_track_data = {}
    # for model in ['HIRAM', 'AM2.5', 'FLOR']:
    #     year_range = (101, 150) if model in ['HIRAM', 'AM2.5'] else (150, 200)
    #     track_data = {model: {}}
    #     for storm_type in ['TS', 'C15w']:
    #         temp = tc_analysis.tc_track_data([model], ['control', 'swishe'], storm_type=storm_type, year_range=year_range)
    #         track_data[model][storm_type] = temp[model]
    #         monthly_track_data[model] = track_data
    #     tc_activity(track_data, model, experiments=['control', 'swishe'], storm_types=['TS', 'C15w'], year_range=year_range, savefig=True)
    ####### End recommended run snippet. ##########
    
    ''' Data loading. '''
    # Import basin data
    basins, basin_masks = basins()

    ''' Data processing. '''
    # Set up dictionary to hold monthly-grouped track data for a given model
    monthly_track_data = {}
    
    # Iterate over each experiment
    for experiment in experiments:
        monthly_track_data[experiment] = {}
        # Iterate over each basin
        for basin in basins.keys():
            monthly_track_data[experiment][basin] = {}
            # Iterate over each storm type
            for storm_type in storm_types:
                monthly_track_data[experiment][basin][storm_type] = {}
                # Retrieve mask from the iterand basin to determine if TCs are being searched for here
                basin_mask = basin_masks[basin]
                # Get pseudonym for iterand basin data
                df = track_data[model_name][storm_type][experiment]['unique'].copy().sort_values('time')
                # Iterate over each month in the DataFrame
                for month, month_data in df.groupby(df['time'].dt.strftime('%m')):
                    # Define helper functions to retrieve the closest GCM data point to the iterand TC center
                    find_lon = lambda x: basin_mask.grid_xt.values[(np.abs(basin_mask.grid_xt.values - x)).argmin()]
                    find_lat = lambda y: basin_mask.grid_yt.values[(np.abs(basin_mask.grid_yt.values - y)).argmin()]
                    # Save the GCM data point coordinates
                    month_data['grid_xt'] = month_data['center_lon'].apply(find_lon)
                    month_data['grid_yt'] = month_data['center_lat'].apply(find_lat)
                    # Retrieve the mask data
                    month_data['mask'] = month_data.apply(lambda c: basin_mask.sel(grid_xt=c['grid_xt'], 
                                                                                   grid_yt=c['grid_yt'], method='nearest').item(), axis=1)
                    # If the mask data is 1, TC is in the iterand basin, and keep it. Else, drop it
                    monthly_track_data[experiment][basin][storm_type][month] = month_data.loc[month_data['mask'] == 1]

    ''' Visualization. '''
    fig, gs = plt.figure(figsize=(14, 4), dpi=144), matplotlib.gridspec.GridSpec(nrows=2, ncols=4, wspace=0.25, hspace=0.5)
    
    axes = {}
    ax_count = 0
    for i in range(0, gs.nrows):
        for j in range(0, gs.ncols):
            sharey = axes['(0, 0)'] if (i > 0 or j > 0) else None
            axes['({0}, {1})'.format(i, j)] = fig.add_subplot(gs[i, j], sharey=sharey)
            ax = axes['({0}, {1})'.format(i, j)]
            basin_name = list(basins.keys())[ax_count]

            statistics = {}
    
            for k, experiment in enumerate(experiments):
                props = cycler(k)
                statistics[experiment] = {storm_type: 0 for storm_type in storm_types}
                for m, storm_type in enumerate(storm_types):
                    months = [int(k) for k, v in monthly_track_data[experiment][basin_name][storm_type].items()]
                    counts = [v['storm_id'].nunique()/(max(year_range) - min(year_range)) 
                              for k, v in monthly_track_data[experiment][basin_name][storm_type].items()]
        
                    width = 0.4
                    hatch = 'xxxx' if storm_type == 'C15w' else ''
                    facecolor = 'white' if storm_type == 'C15w' else props['c']
                    alpha = 1 if storm_type == 'C15w' else 0.625
                    # Bar plot - first pass for hatching color
                    ax.bar(np.array(months) + width*k, counts, width=width, fc=facecolor, ec=props['c'], alpha=alpha, hatch=hatch)
                    # Redo for edges
                    ax.bar(np.array(months) + width*k, counts, width=width, fc='none', ec='k', alpha=1)

                    if m == 0:
                        ax.set_xticks(months, [utilities.month_letter(l) for l in months])

                    statistics[experiment][storm_type] = sum(counts)
            
            ax.set_title('{0}, {1}'.format(model_name, basin_name), loc='left', fontsize=10)
            ax_count += 1

            ax.annotate('{0}: CTL = {1:.1f}, EXP = {2:.1f}'.format('TS', statistics['control']['TS'], statistics['swishe']['TS']),
                        xy=(0.03, 0.96), xycoords='axes fraction', fontsize=7, ha='left', va='top')
            ax.annotate('{0}: CTL = {1:.1f}, EXP = {2:.1f}'.format('HU', statistics['control']['C15w'], statistics['swishe']['C15w']),
                        xy=(0.03, 0.84), xycoords='axes fraction', fontsize=7, ha='left', va='top')

    if savefig:
        filename = 'model_{0}-storm_count_monthly-{1}-year_s{2}_e{3}.png'.format(model_name, storm_type, min(year_range), max(year_range))
        plt.savefig(os.path.join('/tigress/GEOCLIM/gr7610/figs/climate', filename), dpi=300, bbox_inches='tight')

    return monthly_track_data

def swishe_TC_overlay(model_name, storm_type='TS', years=None, experiments=['CTL1990s', 'CTL1990s_swishe'], FLOR_year_adjustment=1900, bin_resolution=2.5, dpi=144):
    
    start_year, end_year = min(years), max(years)
    track_data = tc_analysis.tc_track_data([model_name], experiments, year_range=(start_year, end_year), storm_type=storm_type, FLOR_year_adjustment=FLOR_year_adjustment)
    
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
        for experiment in experiments:
            # Get longitude and latitude bins
            x, y = density[model][experiment].lon.unique(), density[model][experiment].lat.unique()
            # Get density array
            v = np.reshape(density[model][experiment]['count'].values, (len(x), len(y)))
            # Assign to dictionary for the given model/experiment configuration. 
            # Normalize by number of years and by day (assume 6-hourly data, so 6 x 4 = 1 day)
            densities[model][experiment] = v.T/(end_year - start_year)
            
    fig, gs = plt.figure(figsize=(8, 3), dpi=dpi), matplotlib.gridspec.GridSpec(nrows=1, ncols=1)
    fig, gs, ax, proj_ref = basemap(fig, gs, model_name, 'TC density (blue) & SWISHE application (red)', 
                                    year_range=(start_year + FLOR_year_adjustment, end_year + FLOR_year_adjustment), land=True)

    contours = np.arange(0, 1.6 + 0.2, 0.2)
    norm = matplotlib.colors.BoundaryNorm(contours, 256)
    cmap = cmap_white_adjust('Blues', levels=len(contours)-1)

    im = ax.contourf(x, y, densities[model_name][experiments[0]], transform=proj_ref, norm=norm, cmap=cmap, levels=len(contours)-1)

    swishe_contours = swishe_frequency(models=[model_name], dpi=144, set_visible=False)
    im_swishe = ax.contour(swishe_contours.grid_xt, swishe_contours.grid_yt, 100*swishe_contours.sum(dim='time')/len(swishe_contours.time.values), 
                           levels=contours[1:], cmap='Reds', norm=norm, vmin=0, vmax=max(contours), linewidths=1.5, transform=proj_ref)
    
    cax = ax.inset_axes([1.03, 0, 0.02, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
    colorbar.set_label('% of time with TC', labelpad=15, rotation=270)
    
def circular_subplot(fig, axes, data, norm, cmap, mask_degrees=10, plot_style='pcolormesh'):

    '''
    Constructs a circular subplot centered on TC centers, akin to Figure 2 in Fischer et al. (2018).
    Reference: 10.1175/MWR-D-17-0239.1
    '''

    # Construct a clipping patch to hide the pcolormesh
    mask_num_vertices = 100
    theta = np.linspace(0, 2*np.pi, mask_num_vertices)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    mask_patch = matplotlib.path.Path(verts * radius + center)

    # Plot the data and mask it
    if plot_style == 'pcolormesh':
        im = axes.pcolormesh(data.grid_xt, data.grid_yt, data, norm=norm, cmap=cmap, clip_path=(mask_patch, axes.transAxes))
    else:
        im = axes.contourf(data.grid_xt, data.grid_yt, data, norm=norm, cmap=cmap, clip_path=(mask_patch, axes.transAxes))
    axes.set_xlim([-mask_degrees, mask_degrees])
    axes.set_ylim(axes.get_xlim())
    axes.set_aspect('equal')
    axes.set_axis_off()

    # Build a polar grid
    ax_polar = axes.inset_axes([0, 0, 1, 1], polar=True, frameon=False)
    ax_polar.set_rmax(mask_degrees)
    ax_polar.set_yticks(np.arange(0, mask_degrees, 2))
    ax_polar.set_yticklabels([])
    ax_polar.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2], ['N', '', '', ''])
    ax_polar.set_theta_zero_location('N')
    ax_polar.set_theta_direction(-1)
    ax_polar.grid(True, c='k', linewidth=0.5, linestyle='--', alpha=0.25)
    
    # Construct a clipping patch to hide the pcolormesh
    mask_num_vertices = 100
    theta = np.linspace(0, 2*np.pi, mask_num_vertices)
    center, radius = [0, 0], mask_degrees
    verts = np.vstack([theta, [radius]*mask_num_vertices]).T
    ring_patch = matplotlib.path.Path(verts)
    patch = matplotlib.patches.PathPatch(ring_patch, facecolor='none', edgecolor='black', linewidth=1.5)
    # Add the patch to the axes
    ax_polar.add_patch(patch)
    # circ =  []
    # circ.append(patch)
    # coll = matplotlib.collections.PatchCollection(circ, zorder=10, linewidth=2, facecolor='none')
    # axes.add_collection(coll)
    
    return axes