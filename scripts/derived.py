import numpy as np
import xarray as xr

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