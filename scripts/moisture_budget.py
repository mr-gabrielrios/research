import numpy as np, xarray as xr
import importlib, os
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

import utilities, visualization
# Make any mathematical text render in LaTeX
matplotlib.rcParams['mathtext.fontset'] = 'cm'

def eddy(data, model, experiments, fields=['ucomp', 'vcomp', 'sphum']):
    """
    Derive eddy values from a monthly mean and output into xArray DataArray.
    
    Args:
    - data (dict): 3-tiered dictionary (model, experiment, field), with the field key containing an xArray Dataset
    - model (str): name of the model being processed
    - experiments (list): list of experiments being processed
    - fields (list): list of fields being processed
    Returns:
    - data (dict): 3-tiered dictionary (model, experiment, field), with the field key containing an xArray Dataset
    
    """
   # For each model, experiment, and field: compute the monthly mean and subtract it from the daily data
    for experiment in experiments:
        for field in fields:
            daily_anom, monthly_mean = [], []
            # Sample each month
            print(model, experiment, field)
            for month, monthly_data in data[model][experiment][field].resample(time='1M'): 
                # Get monthly mean
                m = monthly_data.mean(dim='time')
                # Get daily anomaly by subtracting the monthly mean that month's daily data 
                anom = data[model][experiment][field].sel(time=((data[model][experiment][field].time.dt.year == month.year) &
                                                                (data[model][experiment][field].time.dt.month == month.month))) - m
                monthly_mean.append(m)
                daily_anom.append(anom)
            data[model][experiment]['{0}_mean'.format(field)] = xr.concat(monthly_mean, dim='time').sortby('time')
            data[model][experiment]['{0}_eddy'.format(field)] = xr.concat(daily_anom, dim='time').sortby('time')
    
    return data

def mma(data, param='sphum'):
    ''' Derive mean moisture advection (MMA). '''
    dx_q = utilities.domain_differentiation(data['{0}_mean'.format(param)], dim='grid_xt')
    dy_q = utilities.domain_differentiation(data['{0}_mean'.format(param)], dim='grid_yt')
    data['mma'] = (data['ucomp_mean']*dx_q + data['vcomp_mean']*dy_q).rename({'time': 'month'})
    return data

def ema(data, param='sphum'):
    ''' Derive eddy moisture advection (EMA). '''
    dx_q = utilities.domain_differentiation(data['{0}_eddy'.format(param)], dim='grid_xt')
    dy_q = utilities.domain_differentiation(data['{0}_eddy'.format(param)], dim='grid_yt')
    data['ema'] = (data['ucomp_eddy']*dx_q + data['vcomp_eddy']*dy_q).groupby('time.month').mean()
    return data

def mmc(data, param='sphum'):
    ''' Derive mean moisture convergence (EMA). '''
    dx_u = utilities.domain_differentiation(data['ucomp_mean'], dim='grid_xt')
    dy_v = utilities.domain_differentiation(data['vcomp_mean'], dim='grid_yt')
    data['mmc'] = -(data['{0}_mean'.format(param)]*dx_u + data['{0}_mean'.format(param)]*dy_v).rename({'time': 'month'})
    return data

def emc(data, param='sphum'):
    ''' Derive eddy moisture convergence (EMC). '''
    dx_u = utilities.domain_differentiation(data['ucomp_eddy'], dim='grid_xt')
    dy_v = utilities.domain_differentiation(data['vcomp_eddy'], dim='grid_yt')
    data['emc'] = -(data['{0}_eddy'.format(param)]*dx_u + data['{0}_eddy'.format(param)]*dy_v).groupby('time.month').mean()
    return data

def vertical_integral(data, field):
    ''' Get a vertical integral for a given field. '''
    # Initialize empty list for future concatenation
    vi = []
    # Iterate over each level to get the proper integrand, starting from the surface moving up (the [::-1][:-1] reverses the list and trims the last index)
    # Use midpoint rule to integrate [f_{p+1/2} * dp], where f_{p+1/2} = (f_{p} + f_{p+1})/2
    for p, level in enumerate(data[field].pfull.values[::-1][:-1]):
        # print(p, level)
        # Get pressure level differential [dp], in Pa
        dp = 100*(data[field].isel(pfull=p+1).pfull.item() - data[field].isel(pfull=p).pfull.item())
        # Get the quantity at the integrand midpoint [f_{p+1/2}], in kg kg^-1 s^-1
        f_p12 = (data[field].isel(pfull=p+1) + data[field].isel(pfull=p))/2        
        # Get their product, [f_{p+1/2} * dp], in Pa kg kg^-1 s^-1 
        integrand = f_p12*dp
        # Build the xArray Dataset with the pressure value
        integrand = integrand.assign_coords(pfull=level).expand_dims('pfull')
        # Add to the list
        vi.append(integrand)
    # Concatenate the list and change units to mm d^-1
    data['vi_{0}'.format(field)] = xr.concat(integrand, dim='pfull').sum(dim='pfull') * 1000 * 86400 / (1000*9.81)
    data['vi_{0}'.format(field)].attrs = {'units': 'mm d^-1'}

    return data

def load(model, levels, years=(125, 150), extent=[0, 359, -40, 40], coarsen_factor=4):
    dirname = '/projects/GEOCLIM/gr7610/analysis/model_out'
    fields = ['ucomp', 'vcomp', 'sphum']
    experiments = {'control': 'CTL1990s_tigercpu_intelmpi_18_540PE',
                   'swishe': 'CTL1990s_swishe_tigercpu_intelmpi_18_540PE'}
    data = {}
    data[model] = {}
    for experiment in experiments.keys():
        data[model][experiment] = {}
        for field in fields:
            temp = []
            for level in levels:
                print(model, experiment, field, level)
                filenames = [os.path.join(dirname, file) for file in os.listdir(dirname)
                            if (model in file) and (str(level) in file) and (field in file) 
                            and (experiments[experiment] in file) and ('{0}_{1}'.format(min(years), max(years)) in file)] 
                if len(filenames) == 1:
                    temp.append(xr.open_dataset(filenames[0]).sel(grid_xt=slice(extent[0], extent[1]), grid_yt=slice(extent[2], extent[3]))[field])
                else:
                    print('Unable to process, make sure only one file exists in the directory that matches these criteria. Files found:')
                    print(filenames)
            data[model][experiment][field] = xr.concat(temp, dim='pfull').coarsen(grid_xt=coarsen_factor, grid_yt=coarsen_factor, boundary='trim').mean()
            
    return data
            
def plot_2D(data, model, experiment, start_year, end_year, months=[1, 12], basin='custom', extent=None, contour=False, dpi=96):
    
    terms = {'vi_mma': {'eq': "$\\frac{1}{\\rho_w g} \int \ \overline{\overline{\mathbf{u}} \cdot \\nabla \overline{q}} \ \mathrm{d}p $"},
             'vi_ema': {'eq': "$\\frac{1}{\\rho_w g} \int \ \overline{\mathbf{u'} \cdot \\nabla q'} \ \mathrm{d}p $"},
             'vi_mmc': {'eq': "$\\frac{1}{\\rho_w g} \int \ \overline{q \\nabla \cdot \mathbf{u}} \ \mathrm{d}p $"}, 
             'vi_emc': {'eq': "$\\frac{1}{\\rho_w g} \int \ \overline{q' \\nabla \cdot \mathbf{u}'} \ \mathrm{d}p $"}}
    fig, gs = [plt.figure(figsize=(12, 6), dpi=dpi), 
            matplotlib.gridspec.GridSpec(nrows=2, ncols=2, wspace=0.2, hspace=0.1)]

    month_slice = range(min(months), max(months))

    basins = visualization.basins()
    if basin in basins.keys():
        extent = basins[basin]
    else:
        extent = [0, 359, -40, 40] if not extent else extent

    levels = 16
    base_proj = ccrs.PlateCarree()
    counter = 0
    for i in range(0, gs.ncols):
        for j in range(0, gs.nrows):
            term = list(terms.keys())[counter]

            fig, gs, ax, proj_ref = visualization.basemap(fig, gs, model, '{0}, {1}'.format(term.split('_')[-1].upper(), experiment), year_range=(start_year, end_year), row_num=j, col_num=i, extent=extent)
            terms[term]['ax'] = ax

            if experiment == 'control':
                terms[term]['data'] = (data[model]['control'][term].sel(month=month_slice).mean(dim='month')).sel(grid_xt=slice(extent[0], extent[1]),
                                                                                                                  grid_yt=slice(extent[2], extent[3]))
            elif experiment == 'swishe':
                terms[term]['data'] = (data[model]['swishe'][term].sel(month=month_slice).mean(dim='month')).sel(grid_xt=slice(extent[0], extent[1]),
                                                                                                                 grid_yt=slice(extent[2], extent[3]))
            elif experiment == 'swishe-control':
                terms[term]['data'] = (data[model]['swishe'][term].sel(month=month_slice).mean(dim='month') - 
                                       data[model]['control'][term].sel(month=month_slice).mean(dim='month')).sel(grid_xt=slice(extent[0], extent[1]),
                                                                                                                  grid_yt=slice(extent[2], extent[3]))
            
            extrema = [-3, 3] if term in ['vi_mma', 'vi_mmc'] else None
            norm, cmap = visualization.norm_cmap(terms[term]['data'], None, num_bounds=levels, extrema=extrema)
            if contour:
                terms[term]['im'] = ax.contourf(terms[term]['data'].grid_xt, terms[term]['data'].grid_yt, terms[term]['data'], 
                                                levels=levels, norm=norm, cmap=cmap, transform=base_proj)
            else:
                terms[term]['im'] = ax.pcolormesh(terms[term]['data'].grid_xt, terms[term]['data'].grid_yt, 
                                                  terms[term]['data'], norm=norm, cmap=cmap, transform=base_proj)
                
            ax.axhline(0, c='k', lw=1, alpha=0.5)
            
            cax = ax.inset_axes([0, -0.625, 1, 0.075])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
            colorbar.set_label(data[model]['swishe'][term].attrs['units'])
            
            ax.coastlines()
            counter += 1
            
def plot_zonal_mean(data, model, experiments, start_year, end_year, months=[1, 12], term='vi_ema', basin='custom', extent=None):
    fig, axes = plt.subplots(figsize=(4, 3), ncols=2, sharey=True)

    output = {}

    term = 'vi_ema'
    basins = visualization.basins()
    if basin in basins.keys():
        extent = basins[basin]
    else:
        extent = [0, 359, -40, 40] if not extent else extent
    month_slice = range(min(months), max(months))
        
    colors = ['b', 'r']
    for experiment in experiments:
        output[experiment] = {'mean': data[model][experiment][term].sel(month=month_slice, 
                                                                        grid_xt=slice(extent[0], extent[1]),
                                                                        grid_yt=slice(extent[2], extent[3])).mean(dim='month').mean(dim='grid_xt'),
                            'std': data[model][experiment][term].sel(month=month_slice, 
                                                                        grid_xt=slice(extent[0], extent[1]),
                                                                        grid_yt=slice(extent[2], extent[3])).mean(dim='month').std(dim='grid_xt')}

    plot_std = True
    axes[0].axvline(0, c='k', alpha=0.5, lw=1)
    for i, experiment in enumerate(experiments):
        if plot_std:
            axes[0].fill_betweenx(output[experiment]['mean'].grid_yt, 
                                x1=output[experiment]['mean'] - output[experiment]['std'], 
                                x2=output[experiment]['mean'] + output[experiment]['std'], 
                            alpha=0.1, fc=colors[i])
        axes[0].plot(output[experiment]['mean'], output[experiment]['mean'].grid_yt, lw=2, c=colors[i])
        
    axes[0].set_ylim([-30, 30])
    axes[0].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    axes[0].grid(which='both', ls=':')
    axes[0].set_ylabel('Latitude')
    title_orig = axes[0].set_title('{0}'.format(term.split('_')[-1].upper()), ha='left', x=0, fontsize=10)

    axes[1].axvline(0, c='k', alpha=0.5, lw=1)
    axes[1].plot((output['swishe']['mean'] - output['control']['mean']).rolling(grid_yt=8, center=True).mean(), 
                output[experiment]['mean'].rolling(grid_yt=8).mean().grid_yt, lw=2, c='k')
    axes[1].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    axes[1].grid(which='both', ls=':')
    title_diff = axes[1].set_title('$\delta$ ({0})'.format(term.split('_')[-1].upper()), ha='left', x=0, fontsize=10)

    if basin == 'custom':
        basin = extent
    title_str = '{0}, {1} to {2}, {3}'.format(model, start_year, end_year, basin)
    fig.tight_layout()

    title_xpos, title_ypos = fig.transFigure.inverted().transform(title_orig._get_xy_display())
    fig.suptitle(title_str, x=title_xpos, y=title_ypos+0.125, fontsize=10, ha='left')
                
def main(model, levels, years=(125, 150), extent=[0, 359, -40, 40], coarsen_factor=4):
    experiments = ['control', 'swishe']
    
    data = load(model, levels, years, extent, coarsen_factor=coarsen_factor)
    data = eddy(data, model, experiments)
    
    for experiment in experiments:
        print(model, experiment)
        data[model][experiment] = mma(data[model][experiment], param='sphum')
        data[model][experiment] = ema(data[model][experiment], param='sphum')
        data[model][experiment] = mmc(data[model][experiment], param='sphum')
        data[model][experiment] = emc(data[model][experiment], param='sphum')
        
    land_mask = xr.open_dataset('/tigress/GEOCLIM/gr7610/tools/land_mask.nc')['land_mask'].sel(grid_xt=slice(extent[0], extent[1]), 
                                                                                               grid_yt=slice(extent[2], extent[3])).coarsen(grid_xt=coarsen_factor, grid_yt=coarsen_factor, boundary='trim').mean()
    
    integral_fields = ['mma', 'ema', 'mmc', 'emc']
    for experiment in experiments:
        for integral_field in integral_fields:
            print(model, experiment, integral_field)
            data[model][experiment] = vertical_integral(data[model][experiment], field=integral_field)
            data[model][experiment]['vi_{0}'.format(integral_field)] = data[model][experiment]['vi_{0}'.format(integral_field)].where(land_mask == 0)
            
    return data

# if __name__ == '__main__':
#     data = main(model='HIRAM', levels=[200, 400, 500, 700, 850], coarsen_factor=2)