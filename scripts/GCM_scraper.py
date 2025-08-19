import argparse
import cftime
import datetime
import time
import numpy as np
import os
import pandas as pd
import xarray as xr

import utilities

import warnings
warnings.filterwarnings("ignore")


def access(filename='/projects/GEOCLIM/gr7610/reference/generals_data_processing.csv'):
    ''' Method to scrape CSV file for relevant data and format it in a DataFrame. '''

    # Read CSV into a Pandas DataFrame
    data = pd.read_csv(filename).drop(columns=['Status'])
    data.columns = map(str.lower, data.columns)

    data['start_year'] = data['year range'].str.split(' ').str[0].astype(int)
    data['end_year'] = data['year range'].str.split(' ').str[-1].astype(int)
    data = data.rename(columns={'data frequency': 'data type',
                                'resampling frequency': 'resample',
                                'pressure level': 'level'})
    data['level'] = data['level'].fillna('None')

    return data

def generator(data, experiment_names, diagnostic=False):

    # Initialize output list
    output_strs = []
    # Initialize output package for processing within a Jupyter notebook
    jupyter_pkg = []
    for i in range(0, len(data)):
        # Define specific row package
        row_pkg = {}
        row = data.iloc[i]
        # 1. Get model name
        model = row['model']
        if diagnostic:
            print('1. Accessing data for the {0} model...'.format(model))
        row_pkg['model'] = model
        # 2. Get experiment and corresponding experiment name
        experiment = row['experiment'].lower()
        experiment_name = experiment_names[model][experiment]
        if diagnostic:
            print('2. Data accessed for the {0} configuration...'.format(
                experiment_name))
        row_pkg['experiment'] = experiment_name
        # 3. Custom directory name
        custom_dir = utilities.directories(
            model, experiment_name, data_type='model_output')
        if diagnostic:
            print('3. Custom directory being used: {0}...'.format(custom_dir))
        # 4. Get data type
        data_type = row['domain'] + '_' + row['data type']
        if diagnostic:
            print('4. Data type extracted: {0}...'.format(data_type))
        row_pkg['data_type'] = data_type
        # 5. Get field name
        field = row['field']
        if diagnostic:
            print('5. Field extracted is: {0}...'.format(field))
        row_pkg['field'] = field
        # 6. Get level data
        level = row['level'].replace(' ', '')
        if 'm' in level:
            level = level.split('m')[0]
        elif 'hPa' in level:
            level = level.split('hPa')[0]
        if diagnostic:
            print('6. Data extracted at following level: {0}...'.format(level))
        row_pkg['level'] = level
        # 7. Get resample frequency
        resample_boolean = True
        time_mean = row['resample']
        if diagnostic:
            print('7. Resampling at {0} frequency...'.format(time_mean))
        row_pkg['resample'] = resample_boolean
        row_pkg['time_mean'] = time_mean
        # 8. Setting up script dates
        start_year, end_year = row['start_year'], row['end_year']
        row_pkg['start_year'] = start_year
        row_pkg['end_year'] = end_year

        # Output string to execute function
        if level == 'None':
            str_out = 'python -m extraction.extraction_submit --model={0} --data_type={1} --years={2}:{3} --experiment={4} --field={5} --resample={7} --time_mean={8} --custom_dir={9} &'.format(
                model, data_type, start_year, end_year, experiment_name, field, level, resample_boolean, time_mean, custom_dir)
        else:
            str_out = 'python -m extraction.extraction_submit --model={0} --data_type={1} --years={2}:{3} --experiment={4} --field={5} --level={6} --resample={7} --time_mean={8} --custom_dir={9} &'.format(
                model, data_type, start_year, end_year, experiment_name, field, level, resample_boolean, time_mean, custom_dir)
        # Append to output list
        output_strs.append(str_out)

        if diagnostic:
            print(str_out)
            print('\n')

        jupyter_pkg.append(row_pkg)

    return output_strs, jupyter_pkg

def model_out(model_name, experiment, data,
              data_var=None, level=None, months=None, resample=False,
              time_mean=None, spatial_mean=None,
              user='gr7610', savefile=False, check_performance=True,
              data_type=None, level_subset=None, overwrite=False):
    ''' Method to save select and/or averaged data to output netCDF file for faster data processing. 

    Assumptions:
    - data is output from a standard GFDL GCM (AM2.5, CM2.5, HIRAM, AM4, CM4)
    - data is for a given year

    Inputs:
    - model_name (str):            name of the model for which data is being processed
    - experiment (str):            name of the experiment for which data is being processed
    - data (xArray Dataset):       Dataset containing yearlong dataset of a given type (type inferred from attributes)
    - data_var (str):              name of data variables of interest
    - level (int or float):        pressure level at which data will be selected (in units of hPa)
    - months (tuple or list):      2-element iterable with a minimum and maximum month
    - resample (bool):             dictates whether averaged data will be resampled to coarser temporal frequency. If 
                                   If True, time mean will dictate resampling frequency. Else, mean will be applied over entire time range.
    - time_mean (str):             temporal frequency of averaging operation 
                                   - valid entries are: 'daily', 'month', 'season', 'year')
    - spatial_mean (str):           method of spatial averaging 
                                   - valid entries are: 'zonal', 'meridional', 'NA', 'EP', 'NI', 'SI', 'AU', 'EP', 'WP')
    - user (str):                  Princeton netID of the user 
                                   (defaults to 'gr7610' because he's likely the only one who'll ever use this)
    - savefile (bool)              boolean to specify whether data that is processed will be saved to Tiger
    - check_performance (bool)     boolean to specify whether performance time checks will be output to console
    '''

    if check_performance:
        counter = 0
        checkpoint = time.time()

    # Correct an unintentional data type change
    level = None if level == 'None' else level

    # Get benchmark Dataset size
    benchmark_size = data.nbytes

    # Get resample frequency. Note that frequency must be greater than or equal to data type temporal frequency.
    resample_frequencies = {'month': 'M',
                            'daily': 'D', 'year': 'Y', 'season': 'QS-DEC'}

    # Get start and end year in model years for output filename
    model_start_year, model_end_year = data.time.min().dt.year.item(
    ), data.isel(time=slice(0, -2)).time.max().dt.year.item()

    ''' Get number of time samples depending on averaging. '''
    # Pick start year to move dates into the Pandas range
    # https://pandas.pydata.org/docs/reference/api/pandas.Timestamp.min.html, https://pandas.pydata.org/docs/reference/api/pandas.Timestamp.max.html
    start_year = 1890 if (data.time.min().dt.year <
                          1678 or data.time.max().dt.year > 2361) else 0
    # Define minimum and maximum times
    min_time = datetime.datetime(year=data.time.min().dt.year.item() + start_year,
                                 month=data.time.min().dt.month.item(),
                                 day=data.time.min().dt.day.item(),
                                 minute=data.time.min().dt.minute.item(),
                                 second=data.time.min().dt.second.item())
    max_time = datetime.datetime(year=data.time.max().dt.year.item() + start_year,
                                 month=data.time.max().dt.month.item(),
                                 day=data.time.max().dt.day.item(),
                                 minute=data.time.max().dt.minute.item(),
                                 second=data.time.max().dt.second.item())

    print('Time bounds: {0} to {1}'.format(min_time, max_time))

    # Get number of time sample depending on averaging
    if time_mean:
        min_time, max_time = [cftime.datetime(min_time.year, min_time.month, min_time.day),
                              cftime.datetime(max_time.year, max_time.month, max_time.day)]
        # min_time, max_time = pd.to_datetime(min_time), pd.to_datetime(max_time)
        if time_mean == 'year':
            num_samples_time = (max_time.year - min_time.year) + \
                (max_time.month - min_time.month)
        elif time_mean == 'month':
            num_samples_time = 12 * \
                ((max_time.year - min_time.year) +
                 (max_time.month - min_time.month))
        else:
            num_samples_time = 365 * \
                ((max_time.year - min_time.year) +
                 (max_time.month - min_time.month))
    else:
        num_samples_time = len(data.time.values)

    ''' Filename naming convention definitions. '''
    # Initialize output filename components
    operation_substr, time_substr, level_substr = '', '', ''
    # Define substrings denoting the operations to be performed
    if time_mean or spatial_mean:
        operation_substr = '-mean_'
        if time_mean and spatial_mean:
            operation_substr = operation_substr + \
                '{0}_{1}'.format(time_mean, spatial_mean)
        elif time_mean:
            resample_suffix = '-resample' if resample else ''
            operation_substr = operation_substr + time_mean + resample_suffix
        elif spatial_mean:
            operation_substr = operation_substr + spatial_mean

    ''' Perform data paring and metadata-assignment operations. '''
    # Infer data type from input data attributes ("atmos_4xdaily", "atmos_daily", "atmos_month", etc.)
    if not data_type:
        data_type = data.attrs['filename'].split('.')[0]
    # Get keyword for climate domain (atmosphere or ocean)
    domain = data_type.split('_')[0]

    # Get dimensions dependent on data type (atmosphere or ocean)
    dim_zonal = 'grid_xt' if 'atmos' in data_type else 'xt_ocean'
    dim_meridional = 'grid_yt' if 'atmos' in data_type else 'yt_ocean'
    dim_vertical = 'pfull' if 'atmos' in data_type else 'st_ocean'

    # Define substrings denoting the vertical levels to be included
    if level:
        level_substr = '{0}hPa-'.format(
            level) if 'atmos' in data_type else '{0}m-'.format(level)
    else:
        level_substr = 'full-'
    # Define substrings denoting the times to be selected
    if months:
        time_substr = '-year_{0:04d}_{1:04d}-month_{2:02d}_{3:02d}'.format(data.time.min().dt.year.values,
                                                                           data.time.max().dt.year.values,
                                                                           min(months), max(months))
    else:
        time_substr = '-year_{0:04d}_{1:04d}'.format(data.time.min().dt.year.values,
                                                     data.time.max().dt.year.values)

    # Define resampling substring
    resample_str = 'resample-' if resample else ''
    # Define output directory
    out_dirname = '/projects/GEOCLIM/{0}/analysis/model_out'.format(user)
    # Save the output filename
    out_filename = 'model_{0}-exp_{1}-type_{2}-var_{3}-mean_{4}-{5}{6}'.format(
        model_name, experiment, domain, data_var, time_mean, resample_str, level_substr)

    # Define output path
    out_path = os.path.join(out_dirname, out_filename +
                            '{0}_{1}.nc'.format(model_start_year, model_end_year))

    print('Test output path: {0}'.format(out_path))

    if os.path.isfile(out_path) and not overwrite:
        print('File exists and will not be overwritten...')
        return

    if check_performance:
        print('Checkpoint {0}; elapsed time: {1:.4f} s'.format(
            counter, time.time() - checkpoint))
        checkpoint = time.time()
        counter += 1

    # Ensure that a variable is provided
    if data_var is None or data_var not in data.data_vars:
        print(data.data_vars)
        import sys
        print('Data variable specified: {0}'.format(data_var))
        print('Please specify a data variable when you re-run, or make sure that the variable you specified is in the input data. Exiting...')
        sys.exit()

    print(data[data_var].dims, dim_vertical, level == None, level == 'None')
    print('Condition 1:', (level != None or level != 'None')
          and (dim_vertical in data[data_var].dims))
    print('Condition 2:', (level == None or level == 'None')
          and (dim_vertical in data[data_var].dims))

    print(level, dim_vertical, data[data_var].dims, level_subset)
    # Pare down dataset by data variable and vertical level
    if (level != None) and (dim_vertical in data[data_var].dims):
        # For cases when a vertical level is explicitly requested
        data = data[data_var].sel({dim_vertical: level}, method='nearest')
    elif (level != None) and (dim_vertical in data[data_var].dims) and (level_subset is not None) and (level_subset != 'full'):
        # For cases when all vertical levels are needed but a subset is provided
        data = data[data_var].sel(
            {dim_vertical: level_subset}, method='nearest')
    elif (level == None) and (dim_vertical in data[data_var].dims) and (level_subset in 'full'):
        print('All levels being pulled!')
        # For cases when all vertical levels are needed
        data = data[data_var].isel({dim_vertical: range(
            0, len(data[data_var][dim_vertical].values))})
    else:
        # For cases where data is planar
        data = data[data_var]

    if dim_vertical in data.dims:
        print('Values: ', data[dim_vertical].values)
    # Perform month selection
    if months:
        min_month, max_month = min(months), max(months)
        data = data.sel(time=((data.time.dt.month >= min_month) &
                              (data.time.dt.month < max_month)))

    if check_performance:
        print('Checkpoint {0}; elapsed time: {1:.4f} s'.format(
            counter, time.time() - checkpoint))
        checkpoint = time.time()
        counter += 1

    ''' Perform user-defined operations. '''
    # Perform temporal averaging
    if time_mean:
        if resample:
            # Note: SWISHE frequency handled differently
            if data_var == 'swfq':
                print('Handling SWISHE data, with {0} time samples.'.format(
                    len(data.time.values)))
                temp = data.sum(dim='time')/len(data.time.values)
                data = temp.assign_coords(
                    {'time': data.time.values[0]}).expand_dims('time')
            else:
                # Note: this will not work with non-continuous years
                temp = data.resample(
                    time=resample_frequencies[time_mean]).mean()
                data = xr.concat(temp, dim='time').sortby('time')
                del temp
        else:
            data = data.groupby('time.{0}'.format(time_mean))
            data = data.mean(dim='time')
    # Perform temporal averaging
    num_samples_space = np.nan
    if spatial_mean:
        if spatial_mean == 'zonal':
            mean_dim_space = dim_zonal
        if spatial_mean == 'meridional':
            mean_dim_space = dim_meridional

        num_samples_space = np.prod(data[mean_dim_space].shape)
        data = data.mean(dim=mean_dim_space)
    else:
        # Ignoring phalf because none of the relevant variables are on the half-level vertical grid
        num_samples_space = np.array([len(data[data_dim]) for data_dim in data.dims if data_dim in [
                                     dim_zonal, dim_meridional, dim_vertical]])
        num_samples_space = np.prod(num_samples_space)

    print('Time sample size: {0}; Spatial sample size: {1}'.format(
        num_samples_time, num_samples_space))

    if check_performance:
        print('Checkpoint {0}; elapsed time: {1:.4f} s'.format(
            counter, time.time() - checkpoint))
        checkpoint = time.time()
        counter += 1

    print('Script statistics:')
    print('----------------------------------------------------------------------------')
    print('---> Name of output file: {0}\nSize of output file: {1:.3f} GB'.format(
        out_path, data.nbytes/1e9))
    print(
        '---> File reduction: {0:.2f} %'.format(100*(1 - data.nbytes/benchmark_size)))

    if savefile:
        if os.path.isfile(out_path):
            if overwrite:
                data.load().to_netcdf(out_path)
            else:
                print('File exists, revisit this combination...')
        else:
            data.load().to_netcdf(out_path)

    if check_performance:
        print('Checkpoint {0}; elapsed time: {1:.4f} s'.format(
            counter, time.time() - checkpoint))
        checkpoint = time.time()
        counter += 1

    print('\n')

def get_pressure_levels(data_types: str | list | None = None,
                        subset: list[int | float] | None = None):
    ''' Helper function that returns desired pressure levels to index at for data with a vertical component. '''

    if isinstance(data_types, str):
        data_types = [data_types]
    data_types = ['atmos_month'] if not data_types else data_types

    assert isinstance(data_types, list) and len(data_types) == 1
    data_types = data_types[0]
    domain = 'atmos' if 'atmos' in data_types else 'ocean'

    # Register the default pressure levels from GFDL atmospheric and ocean models
    pressure_levels = {'atmos': np.array([2.16404256,   5.84530754,  10.74508016,  17.10653726,
                                          25.11380513,  35.22119682,  48.13790369,  64.55837161,
                                          85.11248401, 110.41833291, 141.09225787, 177.73071875,
                                          220.89206496, 271.06313793, 328.51356472, 392.78924215,
                                          461.94997195, 532.46603469, 600.43331524, 663.10749157,
                                          719.30944054, 768.81882807, 811.8477175, 848.83778205,
                                          880.34978115, 906.99664502, 929.39135743, 948.12523145,
                                          963.73190073, 976.68620119, 987.38959146, 996.10994928]),
                       'ocean': np.array([5.00000000e+00, 1.50000000e+01, 2.50000000e+01, 3.50000000e+01,
                                          4.50000000e+01, 5.50000000e+01, 6.50000000e+01, 7.50000000e+01,
                                          8.50000000e+01, 9.50000000e+01, 1.05000000e+02, 1.15000000e+02,
                                          1.25000000e+02, 1.35000000e+02, 1.45000000e+02, 1.55000000e+02,
                                          1.65000000e+02, 1.75000000e+02, 1.85000000e+02, 1.95000000e+02,
                                          2.05000000e+02, 2.15000000e+02, 2.25000000e+02, 2.36122818e+02,
                                          2.50599976e+02, 2.70620819e+02, 2.98304932e+02, 3.35675629e+02,
                                          3.84634277e+02, 4.46936646e+02, 5.24170593e+02, 6.17736328e+02,
                                          7.28828491e+02, 8.58421509e+02, 1.00725708e+03, 1.17583484e+03,
                                          1.36440625e+03, 1.57297131e+03, 1.80127869e+03, 2.04882861e+03,
                                          2.31487915e+03, 2.59845630e+03, 2.89836523e+03, 3.21320581e+03,
                                          3.54138989e+03, 3.88116211e+03, 4.23062061e+03, 4.58774268e+03,
                                          4.95040869e+03, 5.31642871e+03])}

    if subset:
        subset_mask = np.where((pressure_levels[domain] >= min(
            subset)) and (pressure_levels[domain] <= max(subset)))
        return pressure_levels[domain][subset_mask]
    else:
        return pressure_levels[domain]

def process_entries(models: str | list[str],
                    experiments: str | list[str],
                    fields: str | list[str],
                    data_types: str,
                    levels: str,
                    time_mean: str,
                    year_range: tuple[int, int],
                    jupyter_pkg: list):

    diagnostic_tag = '[process_entries()]'

    # Ensure relevant inputs are in an iterable format
    if isinstance(models, str): models = [models]
    if isinstance(experiments, str): experiments = [experiments]
    if isinstance(fields, str): fields = [fields]

    start_year, end_year = min(year_range), max(year_range)
    pressure_level_subset = get_pressure_levels(data_types=data_types, subset=None) 
                                                if levels in [['None'], ['full']] 
                                                and set(fields).issubset(set(['temp', 'sphum', 'ucomp', 'vcomp', 'omega'])) else levels

    for entry in jupyter_pkg:
        for model in models:
            for field in fields:

                print('==============================================================================================')
                print(f'{diagnostic_tag} Processing entry: {entry}; model name: {model}; experiments: {experiments}; field name: {field}...')
                print(f'\tModel check: {model in entry['model']}; experiment check: {entry['experiment'] in experiments}; field check: {entry['field'] == field}')
                if model == entry['model'] and entry['experiment'] in experiments and entry['field'] == field:

                    print(f'{diagnostic_tag} Intaking {entry}...')
                    print(f'{diagnostic_tag} Data type: {entry['data_type'] in data_types}')
                    print(f'{diagnostic_tag} Field: {entry['field'] in fields}')
                    print(f'{diagnostic_tag} Level: {entry['level'] in levels}')
                    print(f'{diagnostic_tag} Time mean: {entry['time_mean'] not in time_mean}')
                    print(f'{diagnostic_tag} Start year: {int(entry['start_year']) != start_year}')
                    print(f'{diagnostic_tag} End year: {int(entry['end_year']) != end_year}')

                    if (entry['data_type'] not in data_types) or 
                       (entry['field'] not in fields) or 
                       (entry['level'] not in levels) or 
                       (entry['time_mean'] not in time_mean) or 
                       (int(entry['start_year']) != start_year) or 
                       (int(entry['end_year']) != end_year):
                        continue

                    print(f'{diagnostic_tag} Match found! Processing {entry}...')
                    dirname = utilities.directories(model, entry['experiment'])

                    filenames = [os.path.join(dirname, filename) for filename in os.listdir(dirname)
                                 if (('{0}.'.format(entry['data_type']) in filename) or
                                     ('{0}.'.format(entry['data_type'].split('_')[0]) in filename))
                                 and (int(filename[0:4]) >= int(entry['start_year']))
                                 and (int(filename[0:4]) <= int(entry['end_year']))]

                    print(
                        f'{diagnostic_tag} Processing data in directory {dirname} from files:')
                    [print(f) for f in sorted(filenames)]

                    vertical_dimension_name = 'pfull' if 'atmos' in entry['data_type'] else 'st_ocean'
                    temp = xr.open_mfdataset(filenames).sel({vertical_dimension_name: pressure_level_subset}, method='nearest') if pressure_level_subset is [
                        'None'] else xr.open_mfdataset(filenames)

                    # Define pressure level subset name based on how it aligns with the full subset for the given data type
                    pressure_level_subset = 'full' if len(pressure_level_subset) == len(
                        get_pressure_levels(data_types=data_types, subset=None)) else pressure_level_subset
                    model_out(entry['model'], entry['experiment'], temp, data_var=entry['field'], level=entry['level'], resample=entry['resample'], time_mean=entry['time_mean'],
                              data_type=entry['data_type'], level_subset=pressure_level_subset, savefile=True, overwrite=False)

                    del dirname, filenames, temp
                    print(
                        f'{diagnostic_tag} Entry {entry} processed! Onto the next one...')


def main(pathname: str,
         model_names: str | list[str],
         experiment_names: str | list[str],
         field_names: str | list[str],
         data_frequency: str,
         resampling_frequency: str,
         year_range: tuple[int, int]):

    pproc = access()

    reference_experiment_names = {'ERA5': {'reanalysis': 'reanalysis'},
                                  'AM2.5': {'control': 'CTL1990s',
                                            'swishe': 'CTL1990s_swishe',
                                            'control_tiger3': 'CTL1990s_tiger3',
                                            'swishe_tiger3': 'CTL1990s_swishe_tiger3',
                                            'swishe_allflux_tiger3': 'CTL1990s_swishe_allflux_tiger3',
                                            'control_dsst_tiger3': 'CTL1990s_dSST_tiger3',
                                            'ctl.const': 'CTL1990.CONST',
                                            'swishe.const': 'CTL1990_SWISHE.CONST',
                                            'ctl.0n': 'CTL1990.0N',
                                            'swishe.0n': 'CTL1990_SWISHE.0N',
                                            'ctl.15n': 'CTL1990.15N',
                                            'swishe.15n': 'CTL1990_SWISHE.15N',
                                            'ctl.timevar': 'CTL1990.TIMEVAR',
                                            'swishe.timevar': 'CTL1990_SWISHE.TIMEVAR',
                                            'ctl.amip': 'CTL1990.AMIP',
                                            'swishe.amip': 'CTL1990_SWISHE.AMIP'},
                                  'HIRAM': {'control': 'CTL1990s',
                                            'swishe': 'CTL1990s_swishe',
                                            'ctl.const': 'CTL1990.CONST',
                                            'swishe.const': 'CTL1990_SWISHE.CONST',
                                            'ctl.0n': 'CTL1990.0N',
                                            'swishe.0n': 'CTL1990_SWISHE.0N',
                                            'ctl.15n': 'CTL1990.15N',
                                            'swishe.15n': 'CTL1990_SWISHE.15N',
                                            'ctl.timevar': 'CTL1990.TIMEVAR',
                                            'swishe.timevar': 'CTL1990_SWISHE.TIMEVAR'},
                                  'FLOR': {'control': 'CTL1990s',
                                           'control_fa': 'CTL1990s_FA',
                                           'swishe': 'CTL1990s_swishe',
                                           'swishe_fa': 'CTL1990s_swishe_FA',
                                           'control_fa_tiger3': 'CTL1990_FA_tiger3',
                                           'swishe_fa_tiger3': 'CTL1990_swishe_FA_tiger3',
                                           'swishe_fa_ens01_tiger3': 'CTL1990_swishe_FA_ens01_tiger3',
                                           'swishe_fa_ens02_tiger3': 'CTL1990_swishe_FA_ens02_tiger3',
                                           '2xco2': 'CTL1990s_2xco2',
                                           'swishe_2xco2': 'CTL1990s_swishe_2xco2'}}

    data = access(filename=pathname)

    output_strs, jupyter_pkg = generator(
        data, reference_experiment_names, diagnostic=False)

    process_entries(models=model_names,
                    experiments=experiment_names,
                    fields=field_names,
                    data_types=data_frequency,
                    time_mean=resampling_frequency,
                    levels='None',
                    year_range=year_range,
                    jupyter_pkg=jupyter_pkg)


if __name__ == '__main__':

    ''' Collect and process arguments. '''
    parser = argparse.ArgumentParser()

    # Argument intake from command line
    parser.add_argument('--pathname', type=str,
                        help='Pathname to csv file containing data to be processed.')
    parser.add_argument('--model_names', type=str,
                        help='List of model names separated by a colon.')
    parser.add_argument('--experiment_names', type=str,
                        help='List of experiment names separated by a colon.')
    parser.add_argument('--field_names', type=str,
                        help='List of field names separated by a colon.')
    parser.add_argument('--data_frequency', type=str,
                        help='Frequency at which source data occurs.')
    parser.add_argument('--resampling_frequency', type=str,
                        help='Frequency at which resampled data is desired.')
    parser.add_argument('--year_range', type=str,
                        help='Years between which data are generated, separated by a colon.')

    args = parser.parse_args()

    # Convert colon-separated substrings into lists of strings
    args.model_names = args.model_names.split(':')
    args.experiment_names = args.experiment_names.split(':')
    args.field_names = args.field_names.split(':')
    args.year_range = tuple([int(year) for year in args.year_range.split(':')])

    # Show arguments
    print(args)

    ''' Run main function. '''

    main(pathname=args.pathname,
         model_names=args.model_names,
         experiment_names=args.experiment_names,
         field_names=args.field_names,
         data_frequency=args.data_frequency,
         resampling_frequency=args.resampling_frequency,
         year_range=args.year_range)
