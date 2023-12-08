import numpy as np
import xarray as xr
import utilities

def wind(data):
    """
    Horizontal wind.

    Args:
        data (xArray Dataset): Dataset from GCM output.
    Returns:
        data (xArray Dataset): Dataset from GCM output with wind data.
    """ 
    if 'u_ref' in data['tc_model_output'].data_vars and 'v_ref' in data['tc_model_output'].data_vars:
        data['tc_model_output']['wind'] = np.sqrt(data['tc_model_output']['u_ref']**2 + data['tc_model_output']['v_ref']**2)
        data['tc_model_output']['wind'].attrs = {'long_name': 'horizontal wind at 10 m', 'units': 'm/s'}
    
    return data

def rh(data):
    ''' 
    Calculate relative humidity (H) as a function of atmospheric pressure, specific humidity, and temperature.

    H = e(p, q)/e_s(T)
    e = p/(eps/q - e + 1) # see Emanuel (1994), Eq. 4.1.4
    e_s = 6.112*exp(17.67*t/(t+243.5)) # see Bolton (1980), Eq. 10. 
    
    Bolton (1980) doi: https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2
    ''' 

    # Define shorthand for variables
    p, q, t = data['tc_vertical_output'].pfull, data['tc_vertical_output']['sphum'], data['tc_vertical_output']['temp']

    eps = 0.622 # ratio of R_d to R_v
    # Get vapor pressure
    e = p/(eps/q - eps + 1)
    # Convert temperature
    tc = t - 273.16
    # Get saturation vapor pressure
    e_s = 6.112*np.exp(17.67*tc/(tc + 243.5))
    # Save to xArray DataArray
    data['tc_vertical_output']['rh'] = 100*e/e_s
    data['tc_vertical_output']['rh'].attrs = {'long_name': 'relative humidity', 'units': '%'}

    return data

def density(data):
    """Derive density of air using Bolton (1980) and Huang (2018).

    Args:
        data (xarray Dataset): Dataset

    Returns:
        data (xarray Dataset): Dataset
    """
    
    # Get relevant constants
    R_d = utilities.get_constants('R_d')
    R_v = utilities.get_constants('R_v')
    
    # Reference: Huang (2018), doi:10.1175/JAMC-D-17-0334.1
    e_s_huang_water = lambda t: np.exp(34.494 - 4924.99/((t-273.15)+237.1))/((t-273.15) + 105)**(1.57)
    e_s_huang_ice = lambda t: np.exp(43.494 - 6545.8/((t-273.15)+278))/((t-273.15) + 868)**2
    # Create shorthand name for vertical data
    iterdata = data['tc_vertical_output']
    # Modify saturation vapor pressure (e_s) approximation based on temperature
    e_s = xr.where(iterdata['temp'] >= 273.15, e_s_huang_water(iterdata['temp']), e_s_huang_ice(iterdata['temp']))
    # Get vapor pressure in Pa, then divide by 100 to get hPa
    iterdata['e'] = (iterdata['rh']*e_s/100)/100
    iterdata['e'].attrs = {'long_name': 'vapor pressure', 'units': 'hPa'}
    # Get density
    iterdata['rho'] = 100*(iterdata.pfull - iterdata['e'])/(R_d*iterdata['temp']) + 100*iterdata['e']/(R_v*iterdata['temp'])
    # Assign to Dataset
    data['tc_vertical_output']['rho'] = iterdata['rho']
    data['tc_vertical_output']['rho'].attrs = {'long_name': 'density of air', 'units': 'kg m^-3'}
    
    return data

def mse(data, benchmarking=False):
    
    """ Derive moist static energy (MSE).
    
    Methodology: use Eq. 4.5.23 from Emanuel (1994), assuming r_t = 0, and set gz = p/rho
    - Get height from density, and use vapor pressures to deduce density
    - Column-integrated MSE obtained 
    """    
    
    # Get relevant constants
    c_p = utilities.get_constants('c_p')
    L_v = utilities.get_constants('L_v')
    
    # Get density of air
    data = density(data)
    # Create shorthand name for vertical data
    iterdata = data['tc_vertical_output']
    # Get moist static energy
    iterdata['h'] = c_p*iterdata['temp'] + L_v*iterdata['sphum'] + (100*iterdata.pfull)/iterdata['rho']
    iterdata['h'].attrs = {'long_name': 'moist static energy', 'units': 'J kg^-1'}

    data['tc_vertical_output']['h_anom'] = iterdata['h'] - iterdata['h'].mean(dim='grid_xt').mean(dim='grid_yt')
    data['tc_vertical_output']['h_anom'].attrs = {'long_name': 'domainwise moist static energy anomaly', 'units': 'J kg^-1'}

    return data

def lhflx(data): 
    
    """Get latent heat flux from evaporation.
    
    Args:
        data (xArray Dataset): Dataset from GCM output.
    Returns:
        data (xArray Dataset): Dataset from GCM output with wind data.
    """ 
    
    if 'evap' in data['tc_model_output'].data_vars:
        # Multiply evaporation (kg m^-2 s^-1) by latent heat of vaporization (J kg^-1) to get latent heat flux (W m^-2)
        # Assume L_v = 2.5e6 J kg^-1, nominal value from Emanuel (1994)
        data['tc_model_output']['lhflx'] = 2.5e6*data['tc_model_output']['evap']
        data['tc_model_output']['lhflx'].attrs = {'long_name': 'latent heat flux', 'units': 'W m^-2'}
        
    return data

def thflx(data):
    
    """Get turbulent heat flux: thflx = shflx + lhflx.
    
    Args:
        data (xArray Dataset): Dataset from GCM output.
    Returns:
        data (xArray Dataset): Dataset from GCM output with wind data.
    """ 
    
    if 'lhflx' in data['tc_model_output'].data_vars and 'shflx' in data['tc_model_output'].data_vars:
        data['tc_model_output']['thflx'] = data['tc_model_output']['shflx'] + data['tc_model_output']['lhflx']
        data['tc_model_output']['thflx'].attrs = {'long_name': 'net upward surface heat flux', 'units': 'W m^-2'}
        
    return data

def net_lw(data): 
    """Net column longwave radiation. Convention is positive into column (upward at surface, downward at TOA).
    Assume that downward longwave radiation at TOA is equal to downward longwave radiation at surface

    Args:
        data (xArray Dataset): Dataset

    Returns:
        data (xArray Dataset): Dataset
    """
    
    data['tc_model_output']['net_lw'] =  data['tc_model_output']['lwup_sfc'] - data['tc_model_output']['lwdn_sfc'] - data['tc_model_output']['olr']
    data['tc_model_output']['net_lw'].attrs = {'long_name': 'column net longwave radiation', 'units': 'W m^-2'}
    
    return data
    
def net_sw(data):  
    """Net column shortwave radiation. Convention is positive up.

    Args:
        data (xArray Dataset): Dataset

    Returns:
        data (xArray Dataset): Dataset
    """
    
    net_sw_sfc = data['tc_model_output']['swup_sfc'] - data['tc_model_output']['swdn_sfc']
    net_sw_toa = data['tc_model_output']['swup_toa'] - data['tc_model_output']['swdn_toa']
    data['tc_model_output']['net_sw'] = net_sw_sfc - net_sw_toa
    data['tc_model_output']['net_sw'].attrs = {'long_name': 'column net shortwave radiation', 'units': 'W m^-2'}
    
    return data
  
def scalar_divergence(data, field, units=None):
    """Calculate the horizontal scalar divergence for a given field.
    
    For example, this translates to $\nabla \cdot (a b)$, where 'a' is a vector field (assumed to be velocity) and 'b' is a scalar field.
    
    Expanded, $\nabla \cdot (a b)$ = $b \nabla \cdot a + a \cdot \nabla b$
    
    Currently, 1st-order differencing is used for d/dx and d/dy.

    Args:
        data (xarray Dataset): Dataset
        field (str): scalar field
    """
    
    # Get grid distances
    distance = utilities.distance_grid(data['tc_vertical_output'])
    
    # Get advective term, $a \cdot \nabla b$ 
    db_dx = utilities.domain_differentiation(data['tc_vertical_output'], distance, field, 'grid_xt')
    db_dy = utilities.domain_differentiation(data['tc_vertical_output'], distance, field, 'grid_yt')
    adv = (data['tc_vertical_output']['ucomp']*db_dx + data['tc_vertical_output']['vcomp']*db_dy)
    # Append to input Dataset
    data['tc_vertical_output']['adv_{0}'.format(field)] = adv
    data['tc_vertical_output']['adv_{0}'.format(field)].attrs = {'long_name': 'horiz. advection of {0}'.format(data['tc_vertical_output'][field].attrs['long_name']),
                                                                 'units': units}
    
    # Get divergence term, $b \nabla \cdot a$ 
    du_dx = utilities.domain_differentiation(data['tc_vertical_output'], distance, 'ucomp', 'grid_xt')
    dv_dy = utilities.domain_differentiation(data['tc_vertical_output'], distance, 'vcomp', 'grid_yt')
    div = (data['tc_vertical_output'][field]*du_dx + data['tc_vertical_output'][field]*dv_dy)
    # Append to input Dataset
    data['tc_vertical_output']['div_{0}'.format(field)] = div
    data['tc_vertical_output']['div_{0}'.format(field)].attrs = {'long_name': 'horiz. divergence of {0}'.format(data['tc_vertical_output'][field].attrs['long_name']),
                                                                'units': units}
    
    # Append to input Dataset
    data['tc_vertical_output']['flux_{0}'.format(field)] = adv + div
    data['tc_vertical_output']['flux_{0}'.format(field)].attrs = {'long_name': 'horiz. flux of {0}'.format(data['tc_vertical_output'][field].attrs['long_name']),
                                                                  'units': units}
    
    return data    

def vertical_integral(data, field, bottom=950, top=100, benchmarking=False):
    """Vertically integrate a given field.

    Args:
        data (xArray Dataset): xArray Dataset
        field (str): field to integrate
        bottom (int, optional): lower integration bound in hPa. Defaults to 950.
        top (int, optional): upper integration bound in hPa. Defaults to 100.
        benchmarking (bool, optional): boolean to dictate performance benchmarking for debugging. Defaults to False.
    """
    
    # Beware of indexing here - indexing errors may arise from implictly-defined indices
    # for time, pfull, grid_yt, and grid_xt.
    
    # Obtain relevant constants
    g = utilities.get_constants('g')
     # Trim unused dimensions
    for drop_dim in ['bnds', 'phalf']:
        data['tc_vertical_output'] = data['tc_vertical_output'].drop_dims(drop_dim) if drop_dim in data['tc_vertical_output'].dims else data['tc_vertical_output']
    # Initialize container list to be concatenated over time later
    container = []
    # Iterate over each timestep
    for t, timestamp in enumerate(data['tc_model_output'].time.values):
        # Create shorthand name for vertical data
        iterdata = data['tc_vertical_output'].sel(time=timestamp, method='nearest').dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        # Align dimensions to ensure uniformity
        iterdata = iterdata.transpose('pfull', 'grid_yt', 'grid_xt')
        # Initialize container for column vertical integral
        cvi = np.full(iterdata[field].values.shape, np.nan)
        # Get pfull indices that are closest to specified levels
        index_top, index_bottom = [(np.abs(data['tc_vertical_output'].pfull.values - level)).argmin() for level in [top, bottom]]
        # Iterate over each vertical level - bottom pressure level is the second-nearest to surface
        for j in range(index_top, index_bottom+1):
            if benchmarking:
                print('\t Iterand level: {0} hPa'.format(iterdata.pfull.values[j]))
            # Mass-average (see Wing et al, 2019 - Eqn. A1) 
            cvi[j, :, :] = iterdata[field].isel(pfull=j).values * 100*(iterdata.pfull.values[j] - iterdata.pfull.values[j-1])
        # Sum in the vertical
        cvi = np.nansum(cvi, axis=0)/g
        # Remove zeroes to allow for nan dropping
        cvi = np.where(cvi != 0.0, cvi, np.nan)
        # Remove zeroes to allow for nan dropping
        cvi = np.expand_dims(cvi, axis=2)
        # Initiate an xArray data structure to store verrtically-integrated data
        output = xr.DataArray(dims=['grid_yt', 'grid_xt', 'time'], 
                              coords={'grid_yt': (['grid_yt'], iterdata[field].grid_yt.values),
                                      'grid_xt': (['grid_xt'], iterdata[field].grid_xt.values),
                                      'time': (['time'], [timestamp])},
                              data=cvi)
        
        # Plot testing to ensure output is reflected in appended quantity
        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots(figsize=(7, 2), ncols=2)
        # im = ax[0].pcolormesh(cvi[:, :, 0])
        # fig.colorbar(im)
        # fig.suptitle('{0}: {1}'.format(t, timestamp))
        # output.isel(time=0).plot(ax=ax[1])
        
        # Append
        container.append(output)
    # Place vertically-integrated field in the planar Dataset (tc_model_output)
    data['tc_model_output']['vi_{0}'.format(field)] = xr.concat(container, dim='time')
    
    return data
    
def radial_tangential_velocities(data, x, y, R): 
    '''
    Calculate radial and tangential velocity components from zonal and meridional velocities.
    '''

    u, v = data['ucomp'], data['vcomp']

    if np.nanmin(data.grid_yt) < 0:
        u, v = -u, -v
    
    theta = np.arctan(v/u)
    data['v_radial'] = (u*x + v*y)/R
    data['v_tangential'] = (v*x - u*y)/R
    
    data['v_radial'].attrs = {'long_name': 'radial velocity', 'units': 'm s$^{-1}$'}
    data['v_tangential'].attrs = {'v_tangential': 'tangential velocity', 'units': 'm s$^{-1}$'}
    
    return data