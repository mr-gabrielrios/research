import xarray as xr
import os

path = '/projects/GEOCLIM/gr7610/analysis/model_out/model_AM2.5-exp_CTL1990s_tigercpu_intelmpi_18_540PE-var_temp-mean_month-resample-101_150.nc'
data = xr.open_dataset(path)

print(data.data_vars)