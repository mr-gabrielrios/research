# %%

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt
import numba, numpy as np
import os
import pandas as pd
import xarray as xr
import time
# Import local scripts
import accessor, utilities, visualization

def month_letter(month):
    month_letters = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    return month_letters[month-1]

def month_selection(ds, start_month, end_month):
    return (ds >= start_month) & (ds <= end_month)

############################################################
# Begin numba-specific parallelization methods.

@numba.njit
def apply_along_axis_0(func1d, arr):
    """Like calling func1d(arr, axis=0)"""
    if arr.size == 0:
        raise RuntimeError("Must have arr.size > 0")
    ndim = arr.ndim
    if ndim == 0:
        raise RuntimeError("Must have ndim > 0")
    elif 1 == ndim:
        return func1d(arr)
    else:
        result_shape = arr.shape[1:]
        out = np.empty(result_shape, arr.dtype)
        _apply_along_axis_0(func1d, arr, out)
        return out

@numba.njit
def _apply_along_axis_0(func1d, arr, out):
    """Like calling func1d(arr, axis=0, out=out). Require arr to be 2d or bigger."""
    ndim = arr.ndim
    if ndim < 2:
        raise RuntimeError("_apply_along_axis_0 requires 2d array or bigger")
    elif ndim == 2:  # 2-dimensional case
        for i in range(len(out)):
            out[i] = func1d(arr[:, i])
    else:  # higher dimensional case
        for i, out_slice in enumerate(out):
            _apply_along_axis_0(func1d, arr[:, i], out_slice)

@numba.njit
def nb_mean_axis_0(arr):
    return apply_along_axis_0(np.nanmean, arr)
    
# End numba-specific parallelization methods.
############################################################

def xy_plot(model_name, data, median, significance_field, field, months=None, 
            extent=None, vertical_level=None, contour_levels=13, domain_average=True, extrema=None, dpi=96):
    
    factor = 86400 if field in ['precip', 'evap', 'p-e'] else 1
    year_min, year_max = 2001, 2050
    
    min_lon, max_lon, min_lat, max_lat = extent[0], extent[1], extent[2], extent[3]

    ncols = 3
    fig, grid = [plt.figure(dpi=dpi, figsize=(12, 2), constrained_layout=True), 
                 matplotlib.gridspec.GridSpec(nrows=1, ncols=ncols, width_ratios=(0.2, 0.2, 0.6), wspace=0.2, hspace=0)]
    # Define longitudinal offset
    longitude_offset = 180
    # Define projections (working and reference projections)
    proj, proj_ref = ccrs.PlateCarree(central_longitude=longitude_offset), ccrs.PlateCarree()

    data_control = {'mean': data['control'][field]*factor, 
                    'std': data['control'][field]*factor}
    data_swishe = {'mean': data['swishe'][field]*factor,
                    'std': data['swishe'][field]*factor}

    # Pre-define metadata
    long_name, units = data['control'][field].attrs['long_name'], data['control'][field].attrs['units'] if field not in ['precip', 'evap', 'p-e'] else 'mm d$^{-1}$'
    month_str = 'annual' if (months == None or (min(months) == 1 and max(months) == 12)) else ''.join([month_letter(month) for month in months])
    
    for i in range(0, ncols): 
        # Raw comparison
        if i == 0:
            ax = fig.add_subplot(grid[0, 0])
            ax.axhline(0, c='k')
            
            im = ax.plot(data_control['mean'].mean(dim='grid_xt').mean(dim='time'), data_control['mean'].grid_yt, lw=2, c='b', ls=':')
            # im = ax.fill_betweenx(data_control['mean'].grid_yt, data_control['mean']-data_control['std'], data_control['mean']+data_control['std'], fc='b', alpha=0.1)
            im = ax.plot(data_swishe['mean'].mean(dim='grid_xt').mean(dim='time'), data_swishe['mean'].grid_yt, lw=2, c='r', ls='--')
            # im = ax.fill_betweenx(data_swishe['mean'].grid_yt, data_swishe['mean']-data_swishe['std'], data_swishe['mean']+data_swishe['std'], fc='r', alpha=0.1)
            
            ax.set_ylim([min_lat, max_lat])
            ax.set_xlabel('{0} [{1}]'.format(field, units), labelpad=10)
        # Difference
        elif i == 1:
            ax = fig.add_subplot(grid[0, 1])
            ax.axhline(0, c='k')
            ax.axvline(0, c='k')
            # Get difference
            delta = data['swishe'][field].mean(dim='grid_xt').mean(dim='time')*factor - data['control'][field].mean(dim='grid_xt').mean(dim='time')*factor
            # 5-degree rolling mean for 0.5-degree resolution data
            # delta = delta.rolling(grid_yt=10).mean()
            ax.plot(delta, data['control'][field].grid_yt, lw=2, c='k')
            ax.set_ylim([min_lat, max_lat])
            ax.set_xlabel('{0} [{1}]'.format(field, units), labelpad=10)
        else:
            # Initialize subplot
            level_str = '' if not vertical_level else ', {0}'.format(vertical_level)
            fig, grid, ax, proj_ref = visualization.basemap(fig, grid, model_name, 'SWISHE - control{0}, {1}'.format(level_str, month_str), year_range=(year_min, year_max),
                                                            row_num=0, col_num=2, extent=extent, land=True)
            
            ax = fig.add_subplot(grid[0, -1], projection=proj)
            # norm, cmap = norm_cmap(median, field, difference=True, n_bounds=contour_levels, extrema=extrema)
            norm, cmap = visualization.norm_cmap(median, field, num_bounds=contour_levels, extrema=extrema, white_adjust=True)
    
            factor = 86400 if field in ['precip', 'evap', 'p-e'] else 1
            # Plot the data
            im = ax.contourf(data['control'][field].grid_xt, data['control'][field].grid_yt, median.T*factor, 
                             cmap=cmap, norm=norm, levels=contour_levels, transform=proj_ref)
            
            X, Y = np.meshgrid(data['control'][field].grid_xt.values.ravel(), data['control'][field].grid_yt.values.ravel())
            ax.scatter(X.ravel(), Y.ravel(), c=significance_field.T.ravel(), cmap='Greys_r', s=0.5, linewidth=0.5, alpha=0.1, transform=ccrs.PlateCarree())
            del X, Y
            
            cax = ax.inset_axes([1.03, 0, 0.03, 1])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
            
            colorbar.set_label('{0}\n[{1}]'.format(long_name, units), labelpad=30, rotation=270)
        
            # if np.log(abs(norm.vmax)) < -2:
            #     colorbar.ax.set_xticks(colorbar.get_ticks(), ['{0:.1e}'.format(label) for label in colorbar.get_ticks()])
            # else:
            #     colorbar.ax.set_xticks(colorbar.get_ticks(), ['{0:.2f}'.format(label) for label in colorbar.get_ticks()])
            
            ax.add_feature(cartopy.feature.LAND, color=(0.5, 0.5, 0.5, 0.25), zorder=99)
            ax.coastlines()
        
            # if domain_average:
            #     domain_average = np.nanmean(median)*factor
            #     ax.annotate('net $\Delta$ = {0:.2e}'.format(domain_average), (0.99, 0.03), xycoords='axes fraction', fontsize=8, ha='right')

def bootstrap_3d(a, b, plane=None, N=10000, level=0.95):
    '''
    Method to run bootstrap statistical testing on 3-dimensional model output. Dimensions include time, grid_xt (longitude), grid_yt (latitude).
    Current method diagnoses difference in means and establishes statistical significance for a given level from 0 to 1.
    Uses a 2-tailed approach.
    Parallel approach produces a speedup of ~4x.
    See Delsole and Tippett (2013), Chapter 3.6 for details.
    '''

    # Extract correct dimensions by checking to make sure data dimensions
    for dataset in [a, b]:
        # x-y planar data
        if ((plane == 'xy') or (plane is None)) and ('grid_xt' in dataset.dims) and ('grid_yt' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = a.dims.index('time'), a.dims.index('grid_xt'), a.dims.index('grid_yt')
        # y-z planar data 
        elif (plane == 'yz') and ('grid_yt' in dataset.dims) and ('pfull' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = a.dims.index('time'), a.dims.index('grid_yt'), a.dims.index('pfull')
        # x-z planar data
        elif (plane == 'xz') and ('grid_xt' in dataset.dims) and ('pfull' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = a.dims.index('time'), a.dims.index('grid_xt'), a.dims.index('pfull')
        else:
            print('Incorrect dimensions provided. Please check data for correct dimensions.')
            return None
            
    # Extract numeric values from the xArray Datasets
    x, y = a.values, b.values
    # Initialize empty arrays to contain output data
    out_x, out_y = [np.full((N, x.shape[y_axis], x.shape[x_axis]), np.nan), 
                    np.full((N, y.shape[y_axis], y.shape[x_axis]), np.nan)]

    @numba.njit(parallel=True)
    def bootstrap_resample(in_x, in_y, out_x, out_y, N):
        for repetition in numba.prange(0, N):
            out_x_temp, out_y_temp = np.full(in_x.shape, np.nan, dtype=np.float32), np.full(in_y.shape, np.nan, dtype=np.float32)
            for k in range(0, in_x.shape[time_axis]):
                out_x_temp[k] = in_x[np.random.randint(0, in_x.shape[time_axis]), :, :]
                out_y_temp[k] = in_y[np.random.randint(0, in_y.shape[time_axis]), :, :]
            out_x[repetition, :, :] = nb_mean_axis_0(out_x_temp)
            out_y[repetition, :, :] = nb_mean_axis_0(out_y_temp)
        return out_x, out_y

    out_x, out_y = bootstrap_resample(x, y, out_x, out_y, N)
    
    # Get difference between datasets
    delta = out_y - out_x
    # Get values at each respective tail 
    ci_min, ci_max = np.quantile(delta, (1 - level)/2, axis=0), np.quantile(delta, (1 + level)/2, axis=0)
    # Wherever the signs are equal, output the mean. 
    # This indicates that the confidence interval does not intersect 0, such that the null hypothesis is rejected.
    out_binary = np.where(np.sign(ci_min) == np.sign(ci_max), 1, np.nan).T
    out_full = np.where(np.sign(ci_min) == np.sign(ci_max), np.nanmedian(delta, axis=0), np.nan).T
    median = np.nanmedian(delta, axis=0).T

    return out_binary, out_full, median, delta, ci_min, ci_max

def data_access(model, domain, field, frequency='month', vertical_level=None, resample=True):
    dirname = '/projects/GEOCLIM/gr7610/analysis/model_out'
    experiments = ['control', 'swishe']

    data, derived_field_list = {}, derived_fields()
    
    derived = False
    if field in derived_field_list.keys():
        derived = True
        for experiment in experiments:
            data[experiment] = {}
            for subfield in derived_field_list[field]['fields']:
                data[experiment][subfield] = None
                if 'swishe' in experiment:
                    files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                             (model in file) & (domain in file) & (subfield in file) & (frequency in file) & 
                             ('swishe' in file) & ('resample' in file) & (file.endswith('nc'))]
                else:
                    files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                             (model in file) & (domain in file) & (subfield in file) & (frequency in file) & ('ens_' not in file) &
                             ('swishe' not in file) & ('resample' in file) & (file.endswith('nc'))]
                
                # Select vertical level
                if vertical_level:
                    data_levels = []
                    for file in files:
                        data_levels.append(int([substr.split('_')[-1][:-1] for substr in file.split('-') if 'level' in substr][0]))
                    file_index = (np.abs(np.array(data_levels) - vertical_level)).argmin()
                    filename = files[file_index]
                    print('File selection for vertical level {0}: {1}'.format(vertical_level, filename))
                else:
                    filename = files[0]
                    print('Reading data from: {0}'.format(files[0]))
                
                data[experiment][subfield] = xr.open_dataset(filename)
    else:
        for experiment in experiments:
            if 'swishe' in experiment:
                files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                         (model in file) & (domain in file) & (field in file) & (frequency in file) & ('ens_' not in file) & 
                         ('swishe' in file) & ('resample' in file) & (file.endswith('nc'))]
            else:
                files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                         (model in file) & (domain in file) & (field in file) & (frequency in file) & ('ens_' not in file) & 
                         ('swishe' not in file) & ('resample' in file) & (file.endswith('nc'))]
            
            # Select vertical level
            if vertical_level:
                data_levels = []
                for file in files:
                    data_levels.append(int([substr.split('_')[-1][:-1] for substr in file.split('-') if 'level' in substr][0]))
                file_index = (np.abs(np.array(data_levels) - vertical_level)).argmin()
                filename = files[file_index]
                print('File selection for vertical level {0}: {1}'.format(vertical_level, filename))
            else:
                filename = files[0]
                print('Reading data from: {0}'.format(files[0]))
            
            data[experiment] = xr.open_dataset(filename)
            if domain == 'ocean':
                data[experiment]['xt_ocean'] = np.where(data[experiment].xt_ocean < 0, 
                                                        data[experiment].xt_ocean + 360, 
                                                        data[experiment].xt_ocean)
                data[experiment] = data[experiment].sortby('xt_ocean')
                if 'st_ocean' in data[experiment].dims:
                    data[experiment] = data[experiment].rename({'xt_ocean': 'grid_xt', 'yt_ocean': 'grid_yt', 'st_ocean': 'pfull'})
                else:
                    data[experiment] = data[experiment].rename({'xt_ocean': 'grid_xt', 'yt_ocean': 'grid_yt'})
    
    return data, derived

def derived_fields():
    
    fields = {'p-e': {'fields': ['precip', 'evap'],
                      'attrs': {'long_name': 'precipitation - evaporation',
                                'units': 'mm d^-1'}},
              'theta_e': {'fields': ['sphum', 'temp', 'rh'],
                          'attrs': {'long_name': 'equivalent potential temperature',
                                    'units': 'K'}},
              'netsw_200hPa': {'fields': ['swdn_200hPa', 'swup_200hPa'],
                               'attrs': {'long_name': 'net SW flux at 200hPa',
                                         'units': 'watts/m2'}},
              'netrad_200hPa': {'fields': ['netlw_200hPa', 'swdn_200hPa', 'swup_200hPa'],
                                'attrs': {'long_name': 'net radiation at 200hPa',
                                          'units': 'watts/m2'}}}
    return fields 

def derived(data, field, plane='xy'):
    
    derived_fields_list = derived_fields()
    outputs = {}
    
    for experiment in data.keys():
        outputs[experiment] = {}
        if field == 'p-e':
            outputs[experiment][field] = data[experiment]['precip'] - data[experiment]['evap']
        if field == 'netsw_200hPa':
            outputs[experiment][field] = data[experiment]['swdn_200hPa'] - data[experiment]['swup_200hPa']
        if field == 'netrad_200hPa':
            outputs[experiment][field] = data[experiment]['swdn_200hPa'] - data[experiment]['swup_200hPa'] + data[experiment]['netlw_200hPa']
        if field == 'theta_e':
            R_d, R_v, L_v, c_pd, c_l = 287, 461.5, 2.501e6, 1005.7, 4190
            r = data[experiment]['sphum']/(1 - data[experiment]['sphum'])
            c_tot = (c_pd + c_l*r)
            a, b = R_d/c_tot, -r*R_v/c_tot
            outputs[experiment][field] = (data[experiment]['temp']*(1000/data[experiment]['temp'].pfull)**a) * (data[experiment]['rh']**b) * np.exp(L_v*r/(c_tot*data[experiment]['temp']))
            outputs[experiment][field] = outputs[experiment][field].where(outputs[experiment][field].pfull >= 85, np.nan)
            outputs[experiment][field] = outputs[experiment][field].dropna(dim='time', how='all').dropna(dim='pfull', how='all')
        try:
            outputs[experiment][field].attrs = derived_fields_list[field]['attrs']
        except:
            print(outputs)

    return outputs

def prune(data, field, vertical_level=None, plane='xy', extent=None, months=(1, 13), coarsen_factor=None):
    start_month, end_month = min(months), max(months)
    if not extent:
        extent = {'min_lon': 0, 'max_lon': 359, 'min_lat': -60, 'max_lat': 60}
    if not coarsen_factor:
        coarsen_factor = {'x': 2, 'y': 2, 'z': 1}
    if not vertical_level and plane == 'xy':
        vertical_level = 500
        print('Using default pressure level of 500 hPa...')

    outputs = {}
    derived_field_list = derived_fields()
    for experiment in data.keys():
        fields = [f if field in derived_field_list.keys() else field for f in data[experiment].keys()]
        print(fields, extent)
        outputs[experiment] = {}
        for subfield in fields:
            if 'pfull' in data['control'][subfield].dims:
                if plane == 'xy':
                    outputs[experiment][subfield] = data[experiment][subfield].sel(grid_xt=slice(extent['min_lon'], extent['max_lon'])).sel(grid_yt=slice(extent['min_lat'], extent['max_lat'])).sel(pfull=vertical_level, method='nearest').sel(time=month_selection(data[experiment][subfield]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor['x'], grid_yt=coarsen_factor['y']).mean()
                    outputs[experiment][subfield] = outputs[experiment][subfield].transpose('time', 'grid_yt', 'grid_xt').drop_duplicates('time')
                elif plane == 'yz':
                    outputs[experiment][subfield] = data[experiment][subfield].sel(grid_xt=slice(extent['min_lon'], extent['max_lon'])).mean(dim='grid_xt').sel(time=month_selection(data[experiment][subfield]['time.month'], start_month, end_month)).coarsen(grid_yt=coarsen_factor['y'], pfull=coarsen_factor['z']).mean()
                    outputs[experiment][subfield] = outputs[experiment][subfield].transpose('time', 'grid_yt', 'pfull').drop_duplicates('time')
                elif plane == 'xz':
                    outputs[experiment][subfield] = data[experiment][subfield].sel(grid_xt=slice(extent['min_lon'], extent['max_lon'])).mean(dim='grid_yt').sel(time=month_selection(data[experiment][subfield]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor['x'], pfull=coarsen_factor['z']).mean()
                    outputs[experiment][subfield] = outputs[experiment][subfield].transpose('time', 'grid_xt', 'pfull').drop_duplicates('time')
            else:
                outputs[experiment][subfield] = data[experiment][subfield].sel(grid_xt=slice(extent['min_lon'], extent['max_lon']), grid_yt=slice(extent['min_lat'], extent['max_lat'])).sel(time=month_selection(data[experiment][subfield]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor['x'], grid_yt=coarsen_factor['y'], boundary='trim').mean().drop_duplicates('time')
        outputs[experiment] = xr.merge(outputs[experiment].values())
    return outputs

def main(model_name, field, vertical_level=None, domain='atmos', 
         month_range=(1, 12), plane='xy', extent=[0, 359, -60, 60], coarsen_factor=[4, 4, 4],
         N=1000, ci=0.95):

    start_time = time.time()
    start_month, end_month = min(month_range), max(month_range)
    
    extent = {'min_lon': extent[0], 'max_lon': extent[1], 'min_lat': extent[2], 'max_lat': extent[3]}
    coarsen_factor = {'x': coarsen_factor[0], 'y': coarsen_factor[1], 'z': coarsen_factor[2]}

    raw, derived_boolean = data_access(model_name, domain, field, frequency='month', vertical_level=vertical_level, resample=True)
    pruned = prune(raw, field, vertical_level=vertical_level, plane=plane, 
                   extent=extent, months=(start_month, end_month), coarsen_factor=coarsen_factor)
    
    if derived_boolean:
        data = derived(pruned, field)
    else:
        data = pruned.copy()
        
    stipple, median_mask, median, delta, ci_min, ci_max = bootstrap_3d(data['control'][field], data['swishe'][field],
                                                                       plane=plane, N=N, level=ci)
    print('Runtime: {0:.2f} s'.format(time.time() - start_time))
    
    return data, stipple, median_mask, median, delta, ci_min, ci_max