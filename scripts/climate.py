# %%

import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt
import numba, numpy as np
import os
import pandas as pd
import xarray as xr
import time
import visualization

# x-axis scaling, inverse Mercator
def inverse(a):
    a = np.deg2rad(a)
    return np.rad2deg(np.log(np.abs(np.tan(a) + 1.0 / np.cos(a))))

def forward(a):
    a = np.deg2rad(a)
    return np.rad2deg(np.arctan(np.sinh(a)))

def month_letter(month):
    month_letters = ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D']
    return month_letters[month-1]

def month_selection(ds, start_month, end_month):
    return (ds >= start_month) & (ds <= end_month)

def bootstrap_3d(a, b, plane=None, N=10000, level=0.95):
    '''
    Method to run bootstrap statistical testing on 3-dimensional model output. Dimensions include time, grid_xt (longitude), grid_yt (latitude).
    Current method diagnoses difference in means and establishes statistical significance for a given level from 0 to 1.
    Uses a 2-tailed approach.
    See Delsole and Tippett (2013), Chapter 3.6 for details.
    '''

    
    # Check to make sure correct dimensions are in provided data
    for dataset in [a, b]:
        if ((plane == 'xy') or (plane is None)) and ('grid_xt' in dataset.dims) and ('grid_yt' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = x_.dims.index('time'), x_.dims.index('grid_xt'), x_.dims.index('grid_yt')
        elif (plane == 'yz') and ('grid_yt' in dataset.dims) and ('pfull' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = x_.dims.index('time'), x_.dims.index('grid_yt'), x_.dims.index('pfull')
        else:
            print('Incorrect dimensions provided. Please check data for correct dimensions.')
            return None
    
    # Obtain dimensional axes
    if plane == 'None' or plane == 'xy':
        X, Y = a.grid_xt, a.grid_yt
    elif plane == 'yz':
        X, Y = a.grid_yt, a.pfull
    # Extract numeric values from the xArray Datasets
    x, y = a.values, b.values
    # Initialize empty arrays to contain output data
    out_x, out_y = np.full((N, x.shape[x_axis], x.shape[y_axis]), np.nan), np.full((N, y.shape[x_axis], y.shape[y_axis]), np.nan)
    # Perform repetitive re-samplping
    repetition = 0
    while repetition < N:
        # Resample each dataset over its respective time axis
        out_x[repetition, :, :] = np.array([x[np.random.randint(0, x.shape[time_axis]), :, :] for k in range(0, x.shape[time_axis])]).mean(axis=0).T
        out_y[repetition, :, :] = np.array([y[np.random.randint(0, y.shape[time_axis]), :, :] for k in range(0, y.shape[time_axis])]).mean(axis=0).T
        # out_x[repetition, :, :] = np.array([x[k, :, :] for k in range(0, x.shape[time_axis])]).mean(axis=0).T
        # out_y[repetition, :, :] = np.array([y[k, :, :] for k in range(0, y.shape[time_axis])]).mean(axis=0).T
        repetition += 1
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

def bootstrap_3d_parallel(a, b, plane=None, N=10000, level=0.95):

    '''
    Method to run bootstrap statistical testing on 3-dimensional model output. Dimensions include time, grid_xt (longitude), grid_yt (latitude).
    Current method diagnoses difference in means and establishes statistical significance for a given level from 0 to 1.
    Uses a 2-tailed approach.
    Parallel approach produces a speedup of ~4x.
    See Delsole and Tippett (2013), Chapter 3.6 for details.
    '''
    
    # Check to make sure correct dimensions are in provided data
    for dataset in [a, b]:
        if ((plane == 'xy') or (plane is None)) and ('grid_xt' in dataset.dims) and ('grid_yt' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = x_.dims.index('time'), x_.dims.index('grid_xt'), x_.dims.index('grid_yt')
        elif (plane == 'yz') and ('grid_yt' in dataset.dims) and ('pfull' in dataset.dims) and ('time' in dataset.dims):
            time_axis, x_axis, y_axis = x_.dims.index('time'), x_.dims.index('grid_yt'), x_.dims.index('pfull')
        else:
            print('Incorrect dimensions provided. Please check data for correct dimensions.')
            return None
    
    # Obtain dimensional axes
    if plane == 'None' or plane == 'xy':
        X, Y = a.grid_xt, a.grid_yt
    elif plane == 'yz':
        X, Y = a.grid_yt, a.pfull
    # Extract numeric values from the xArray Datasets
    x, y = a.values, b.values
    # Initialize empty arrays to contain output data
    out_x, out_y = np.full((N, x.shape[x_axis], x.shape[y_axis]), np.nan), np.full((N, y.shape[x_axis], y.shape[y_axis]), np.nan)
    # Perform repetitive re-samplping
    # Perform repetitive re-samplping
    repetition = 0

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
    
    @numba.njit(parallel=True)
    def bootstrap_resample(in_x, in_y, out_x, out_y, N):
        # print(out_x.shape)
        for repetition in numba.prange(0, N):
            out_x_temp, out_y_temp = np.full(in_x.shape, np.nan, dtype=np.float32), np.full(in_y.shape, np.nan, dtype=np.float32)
            for k in range(0, in_x.shape[time_axis]):
                out_x_temp[k] = in_x[np.random.randint(0, in_x.shape[time_axis]), :, :]
                out_y_temp[k] = in_y[np.random.randint(0, in_y.shape[time_axis]), :, :]
                # out_x_temp[k] = in_x[k, :, :]
                # out_y_temp[k] = in_y[k, :, :]
            out_x[repetition, :, :] = nb_mean_axis_0(out_x_temp).transpose()
            out_y[repetition, :, :] = nb_mean_axis_0(out_y_temp).transpose()
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

def xy_plot(median, significance_field, field, contour_levels=13, domain_average=True):
    fig, ax = plt.subplots(figsize=(8, 3), subplot_kw={'projection': ccrs.Robinson(central_longitude=180)})

    norm, cmap = visualization.norm_cmap(median, field, difference=True, n_bounds=contour_levels)
    
    factor = 86400 if field in ['precip', 'evap'] else 1
    im = ax.contourf(x_.grid_xt, y_.grid_yt, median*factor, cmap=cmap, norm=norm, levels=contour_levels, transform=ccrs.PlateCarree())
    
    X, Y = np.meshgrid(x_.grid_xt.values.ravel(), y_.grid_yt.values.ravel())
    ax.scatter(X.ravel(), Y.ravel(), c=significance_field.ravel(), cmap='Greys_r', s=0.5, linewidth=0.25, alpha=0.25, transform=ccrs.PlateCarree())
    del X, Y
    
    gl = ax.gridlines(draw_labels=False, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator, gl.ylocator = matplotlib.ticker.FixedLocator([0]), matplotlib.ticker.FixedLocator([0])
    
    cax = ax.inset_axes([0, -0.15, 1, 0.05])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    long_name, units = x_.attrs['long_name'], x_.attrs['units'] if field not in ['precip', 'evap'] else 'mm d^-1'
    colorbar.set_label('{0} [{1}]'.format(long_name, units))

    if np.log(abs(norm.vmax)) < -2:
        colorbar.ax.set_xticks(colorbar.get_ticks(), ['{0:.1e}'.format(label) for label in colorbar.get_ticks()])
    else:
        colorbar.ax.set_xticks(colorbar.get_ticks(), ['{0:.2f}'.format(label) for label in colorbar.get_ticks()])
    
    month_str = [month_letter(m) for m in range(start_month, end_month+1)]
    month_title = ''.join(month_str) if len(month_str) < 12 else 'annually-averaged'
    if 'pfull' in data['control'].dims:
        ax.set_title('{0}, {1} mean\n{2} at {3} hPa, {4:.1f}% significance stippled'.format(model, month_title, field, pressure_level, level*100))
    else:
        ax.set_title('{0}, {1} mean\n{2}, {3:.1f}% significance stippled'.format(model, month_title, field, level*100))
    ax.coastlines()

    if domain_average:
        domain_average = np.nansum(median)
        ax.annotate('net $\Delta$ = {0:.2e}'.format(domain_average), (1, -0.05), xycoords='axes fraction', fontsize=8, ha='right')
    plt.show()

def yz_plot(median, median_masked, field, contour_levels=13, domain_average=True):
    fig, ax = plt.subplots(figsize=(5, 3))

    norm, cmap = norm_cmap(median, field, difference=True, n_bounds=contour_levels) 
    ax.contourf(x_.grid_yt, x_.pfull, median, norm=norm, cmap=cmap, levels=contour_levels)
    ax.axvline(0, color='gray', alpha=0.5, linestyle='--')
    
    X, Y = np.meshgrid(x_.grid_yt.values.ravel(), y_.pfull.values.ravel())
    ax.scatter(X.ravel(), Y.ravel(), c=median_masked.ravel(), cmap='Greys_r', s=1, linewidth=1, alpha=0.25)
    del X, Y
    
    ytick_step = 100
    yticks = np.arange(100, 1000+ytick_step, ytick_step)
    ax.set_xlim([extent['min_lat'], extent['max_lat']])
    ax.set_xscale('function', functions=(forward, inverse))
    ax.set_ylim([100, 1000])
    ax.set_yscale('log')
    ax.set_yticks(yticks, ['{0:0d}'.format(ytick) for ytick in yticks])
    ax.set_ylim(ax.get_ylim()[::-1])
    
    ax.set_xlabel('Latitude [deg N]')
    ax.set_ylabel('Pressure [hPa]')
    
    cax = ax.inset_axes([1.03, 0, 0.03, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    
    if np.log(abs(norm.vmax)) < -2:
        colorbar.ax.set_yticks(colorbar.get_ticks(), ['{0:.1e}'.format(label) for label in colorbar.get_ticks()])
    else:
        colorbar.ax.set_yticks(colorbar.get_ticks(), ['{0:.2f}'.format(label) for label in colorbar.get_ticks()])
        
    month_str = [month_letter(m) for m in range(start_month, end_month+1)]
    month_title = ''.join(month_str) if len(month_str) < 12 else 'annually-averaged'
    ax.set_title('{0}, {1} mean\n{2}, {3:.1f}% significance'.format(model, month_title, field, level*100))
    
    if domain_average:
        domain_average = np.nansum(median)
        ax.annotate('net $\Delta$ = {0:.2e}'.format(domain_average), (0.97, 0.92), xycoords='axes fraction', fontsize=8, ha='right')

def output_field(model, field, frequency='month', resample=True):
    dirname = '/projects/GEOCLIM/gr7610/analysis/model_out'
    experiments = ['control', 'swishe']
    
    data = {}
    for experiment in experiments:
        if 'swishe' in experiment:
            files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                     (model in file) & (field in file) & (frequency in file) & 
                     ('swishe' in file) & ('resample' in file) & (file.endswith('nc'))]
        else:
            files = [os.path.join(dirname, file) for file in os.listdir(dirname) if
                     (model in file) & (field in file) & (frequency in file) & 
                     ('swishe' not in file) & ('resample' in file) & (file.endswith('nc'))]
        print(files[0])
        data[experiment] = xr.open_dataset(files[0]).load()
    return data

if __name__ == '__main__':
    start_time = time.time()
    model, field = 'FLOR', 'temp'

    data = output_field(model, field, frequency='month', resample=True)
    # data = derived_field(model, field, frequency='month', resample=True)

    start_month, end_month = 1, 12
    plane = 'xy'
    extent = {'min_lon': 0, 'max_lon': 359, 'min_lat': -60, 'max_lat': 60}
    coarsen_factor_x, coarsen_factor_y, coarsen_factor_z = 4, 4, 1
    pressure_level = 500

    if 'pfull' in data['control'][field].dims:
        if plane == 'xy':
            x_ = data['control'][field].sel(pfull=pressure_level, method='nearest').sel(time=month_selection(data['control'][field]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor_x, grid_yt=coarsen_factor_y).mean()
            y_ = data['swishe'][field].sel(pfull=pressure_level, method='nearest').sel(time=month_selection(data['swishe'][field]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor_x, grid_yt=coarsen_factor_y).mean()
        elif plane == 'yz':
            x_ = data['control'][field].sel(grid_xt=slice(extent['min_lon'], extent['max_lon'])).mean(dim='grid_xt').sel(time=month_selection(data['control'][field]['time.month'], start_month, end_month)).coarsen(grid_yt=coarsen_factor_y, pfull=coarsen_factor_z).mean()
            y_ = data['swishe'][field].sel(grid_xt=slice(extent['min_lon'], extent['max_lon'])).mean(dim='grid_xt').sel(time=month_selection(data['swishe'][field]['time.month'], start_month, end_month)).coarsen(grid_yt=coarsen_factor_y, pfull=coarsen_factor_z).mean()
    else:
        x_ = data['control'][field].sel(time=month_selection(data['control'][field]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor_x, grid_yt=coarsen_factor_y).mean()
        y_ = data['swishe'][field].sel(time=month_selection(data['swishe'][field]['time.month'], start_month, end_month)).coarsen(grid_xt=coarsen_factor_x, grid_yt=coarsen_factor_y).mean()

    N, level = 100, 0.95
    stipple_field, median_masked, median, delta, ci_min, ci_max = bootstrap_3d_parallel(x_, y_, plane=plane, N=N, level=level)
    print('Runtime: {0:.2f} s'.format(time.time() - start_time))

    if plane == 'xy':
        xy_plot(median, stipple_field, field)
    else:
        yz_plot(median, median_masked, field, contour_levels=25, domain_average=False)