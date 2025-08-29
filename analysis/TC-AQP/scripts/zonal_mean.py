''' Import packages. '''
# Time packages
import calendar, cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import functools, importlib, os, pickle, collections, sys
from multiprocess import Pool
import pandas as pd, xarray as xr, nc_time_axis
xr.set_options(keep_attrs=True)
# Visualization tools
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt

# Local imports
sys.path.insert(1, '/projects/GEOCLIM/gr7610/scripts')
import derived, utilities, visualization, track_TCs

def get_field_name_alias(field_name: str) -> str:

    field_name_aliases = {'WVP': 'PW',
                          'netrad_toa': 'R$_\mathrm{TOA}$',
                          'netrad_sfc': 'R$_\mathrm{sfc}$',
                          'olr': 'R$_\mathrm{TOA, LW, up}$',
                          'swabs_toa': 'R$_\mathrm{TOA, SW, abs}$',
                          'precip': 'P',
                          'evap': 'E',
                          'p-e': 'P-E',
                          'omega_daily': '$\\omega_{500}$',
                          'cld_amt': 'F$_\\text{cld, 200hPa}$',
                          'rh': 'RH$_\\text{500}$',
                          'trans_atm': '$\\langle Q \\rangle_\\mathrm{atm}$',
                          'trans_cpT': '$\\langle Q \\rangle_\\mathrm{sens}$',
                          'trans_lq': '$\\langle Q \\rangle_\\mathrm{lat}$',
                          'q_atm': 'Q$_\mathrm{atm}$',
                          'swfq': 'SWISHE frequency',
                          't_surf': 'SST'}
    
    return field_name_aliases[field_name] if field_name in field_name_aliases.keys() else field_name

def calculate_zonal_mean(dataset: dict,
                         config_names: list[str]=None,
                         field_name: str='WVP',
                         rolling_degrees: int=10,
                         diagnostic: bool=False):


    ''' Obtains GCM output and track data differences between two experiments. '''

    # Determine rolling mean parameters
    degrees_per_grid_point = 0.5 # grid cell size, in degrees
    rolling_points = np.round(rolling_degrees / degrees_per_grid_point).astype(int) # number of grid cells to average over

    # Metadata labeling
    long_name, units = visualization.field_properties(field_name)

    # Get configuration names to plot
    config_names = dataset.keys() if config_names is None else config_names

    # Get model and experiments names
    # Each unique model will be printed in its own row
    # Each unique experiment will be printed in its own column
    model_names = list(set([config_name.split(':')[0] for config_name in config_names]))
    experiment_names = list(sorted(set([config_name.split(':')[1] for config_name in config_names])))

    differences = {}
    
    ''' Pass 1: acquire data properties and plot. '''
    for index in range(len(config_names)):

        # Get configuration properties and make sure configuration names conform to structure and filter rules
        config_name = config_names[index]
        model_name, experiment_names = [config_name.split(':')[0], 
                                        (f'CTL1990.{config_name.split(':')[1]}', f'CTL1990_SWISHE.{config_name.split(':')[1]}')]
        assert len(config_name.split(':')) == 3, 'Revisit configuration name definitions such that the configuration is of form {model_name}:{experiment_name}:{experiment_type}.'
        if 'DIFF' not in config_name: continue
        
        differences[config_name] = {}

        # Determine configuration-specific model and experiment name
        config_model_name, config_experiment_name = [config_name.split(':')[0],
                                                     config_name.split(':')[1]]

        # Get zonal-time mean
        differences[config_name]['field_difference'] = dataset[config_name][field_name].mean('grid_xt').rolling(grid_yt=rolling_points, center=True).mean('grid_yt')
        # Get track data (setting `ax` to None to avoid plotting)
        differences[config_name]['track_frequency'] = get_TC_track_histogram(ax=None, model_name=model_name, experiment_names=experiment_names, 
                                                                get_density_difference=True)
    return differences

def inferred_heat_transport(dataset: xr.DataArray,
                            field_name: str='netrad_toa',
                            flux_adjustment: int|float=0,
                            normalize_units: str='PW'):
    '''
    Credit: Brian Rose, UAlbany; https://github.com/brian-rose/ClimateLaboratoryBook/blob/main/courseware/advanced-heat-transport.ipynb.
    Modifications made by G. Rios.
    
    Compute heat transport as integral of local energy imbalance.
    As either numpy array or xarray.DataArray, returns the heat transport in PW.
    Will attempt to return data in xarray.DataArray if possible.
    '''

    latitude_dimension_name = 'grid_yt'
    # Get averaging dimensions
    averaging_dimensions = [d for d in dataset.dims if d != latitude_dimension_name]
    # Get degrees latitude in units of radians
    lat_radians = dataset[latitude_dimension_name] * np.pi/180
    # Weigh the dataset by latitude
    
    dataset_weighted = (dataset[field_name] - flux_adjustment).mean(averaging_dimensions) * np.cos(lat_radians)
    # Get global integral
    integral = scipy.integrate.cumulative_trapezoid(dataset_weighted, 
                                                    x=lat_radians, 
                                                    initial=0., 
                                                    axis=dataset_weighted.get_axis_num(latitude_dimension_name))

    unit_prefix = {'PW': 1e-15, 'MW': 1e-6, 'GW': 1e-9, 'KW': 1e-3}
    
    result = (unit_prefix[normalize_units] * 2 * np.pi * utilities.get_constants('a')**2 * integral)
    
    if isinstance(dataset_weighted, xr.DataArray):
        result_xarray = dataset_weighted.copy()
        result_xarray.values = result
        return result_xarray
    else:
        return result

def get_derived_fields(dataset: xr.Dataset,
                       diagnostic: bool=False) -> xr.Dataset:

    diagnostic_tag = '[get_derived_fields()]'

    assert isinstance(dataset, xr.Dataset), f'Input data must be an xarray Dataset. It is currently of type {type(dataset)}.'

    # Define each derived field by its constituent fields and the derivation function.
    derived_fields = {'p-e': {'field_names': ['precip', 'evap'],
                              'func': lambda x: x['precip'] - x['evap'],
                              'long_name': 'surface moisture flux',
                              'units': 'mm d$^{-1}$'},
                      'swabs_toa': {'field_names': ['swdn_toa', 'swup_toa'],
                                    'func': lambda x: x['swdn_toa'] - x['swup_toa'],
                                    'long_name': 'absorbed shortwave radiation',
                                    'units': 'W m$^{-2}$'},
                      'netrad_sfc': {'field_names': ['swup_sfc', 'swdn_sfc', 'lwup_sfc', 'lwdn_sfc'],
                                      'func': lambda x: -x['swup_sfc'] + x['swdn_sfc'] - x['lwup_sfc'] + x['lwdn_sfc'],
                                      'long_name': 'net downward surface radiative flux (+ down)',
                                      'units': 'W m$^{-2}$'},
                      'q_atm': {'field_names': ['netrad_toa', 'netrad_sfc', 'shflx', 'precip'],
                                'func': lambda x: x['netrad_toa'] - x['netrad_sfc'] + x['shflx'] + (x['precip'] / 86400) * utilities.get_constants('L_v'),
                                'long_name': 'net atmospheric heating',
                                'units': 'W m$^{-2}$'},
                      'le': {'field_names': ['evap'],
                             'func': lambda x: (x['evap'] / 86400) * utilities.get_constants('L_v'),
                             'long_name': 'latent heat flux',
                             'units': 'W m$^{-2}$'},
                      'div_atm': {'field_names': ['netrad_toa', 'netrad_sfc', 'le', 'shflx'],
                                  'func': lambda x: x['netrad_toa'] - x['netrad_sfc'] +  x['le'] + x['shflx'],
                                  'long_name': 'horizontal divergence of atmospheric energy',
                                  'units': 'W m$^{-2}$'},
                      'trans_atm': {'field_names': ['div_atm'],
                                    'func': lambda x: inferred_heat_transport(dataset=x, field_name='div_atm'),
                                    'long_name': 'northward transport of atmospheric energy',
                                    'units': 'PW'},
                      'div_lq': {'field_names': ['p-e'],
                                 'func': lambda x: -(x['p-e'] / 86400) * utilities.get_constants('L_v'),
                                 'long_name': 'horizontal divergence of latent energy',
                                 'units': 'W m$^{-2}$'},
                      'trans_lq': {'field_names': ['div_lq'],
                                   'func': lambda x: inferred_heat_transport(dataset=x, field_name='div_lq'),
                                   'long_name': 'northward transport of latent energy',
                                   'units': 'PW'},
                      'trans_cpT': {'field_names': ['trans_lq', 'trans_atm'],
                                    'func': lambda x: x['trans_atm'] - x['trans_lq'],
                                    'long_name': 'northward transport of sensible heat',
                                    'units': 'PW'},
                      'omega_daily': {'field_names': ['omega'],
                                'func': lambda x: x['omega'] * 86400,
                                'long_name': 'pressure velocity',
                                'units': 'Pa d$^{-1}$'},
                      'psi_m': {'field_names': ['vcomp'],
                                'func': lambda x: utilities.meridional_overturning(x)['psi_m'],
                                'long_name': 'meridional streamfunction',
                                'units': 'kg s$^{-1}$'}}

    # Iterate over all fields and derived if the constituent fields are in the dataset.
    for derived_field, derived_dict in derived_fields.items():
        if diagnostic:
            print(f'{diagnostic_tag} Inspecting data to derive field {derived_field} using {derived_dict['field_names']}...')
            
        # Check that all fields for derived quantity are in dataset
        check_fields = [field in dataset.data_vars for field in derived_dict['field_names']]
        if not (sum(check_fields) == len(derived_dict['field_names'])):
            if diagnostic:
                print(f'{diagnostic_tag} not all fields found for data when deriving {derived_field}; only fields ...')
            continue
        # Derive the quantity
        dataset[derived_field] = derived_dict['func'](dataset)
        
    return dataset

def get_experiment(field_names: str|list[str],
                   data_type: str,
                   zonal_average: bool,
                   time_average: bool,
                   pressure_level: int|float|None,
                   configuration_name: str,
                   diagnostic: bool=False) -> xr.Dataset:

    ''' Function to load data and perform basic operations for zonal-mean analysis. '''
         
    diagnostic_tag = '[get_experiment_set]'

    if isinstance(field_names, str): field_names = [field_names]

    # Scrape configuration-specific parameters
    configuration_type = configuration_name.split('=')[0] # control (CTL) or experiment (EXP)
    model_name, experiment_name = configuration_name.split('=')[1].split(':')

    # Define directory name to pull data from
    dirname = f'/scratch/gpfs/GEOCLIM/gr7610/tiger3/{model_name}/work/{experiment_name}_tiger3_intelmpi_24_540PE/POSTP'

    # Obtain pathnames matching the iterand criteria
    min_year, max_year = 2, 125 # minimum and maximum potential years in the iterand datasets
    pathnames = [f'{dirname}/{year:04d}0101.{data_type}.nc' for year in range(min_year, max_year)]
    pathnames = [pathname for pathname in pathnames if pathname.split('/')[-1] in os.listdir(dirname)]
    # Hack: only get last N entries if the field names are used for overturning circulations
    pathnames = pathnames[-4:] if ('ucomp' in field_names) or ('vcomp' in field_names) else pathnames

    if diagnostic: print(f'{diagnostic_tag} evaluating configuration name {configuration_name} with pathname {pathnames}')
    # Pull data
    dataset = xr.open_mfdataset(pathnames)[field_names]
    # Find matching pressure level if requested, else load planar data
    # Obtain time-mean
    dataset = dataset.sel(pfull=pressure_level, method='nearest').load() if ('pfull' in dataset.dims) and (pressure_level is not None) else dataset.load()
    # Apply time-averaging, if chosen
    if time_average: dataset = dataset.mean('time')
    # Adjust values based on field name (for example, precipitation and evaporation should be 
    # multiplied by 86400 to convert kg/m^2/s to mm/d)
    dataset = utilities.field_correction(data=dataset)
    # Apply latitude weighting
    dataset = dataset.weighted(np.cos(np.deg2rad(dataset.grid_yt))).mean('grid_xt') if zonal_average else dataset * np.cos(np.deg2rad(dataset.grid_yt))
    # Acquire derived values
    dataset = get_derived_fields(dataset)

    experiment_name_abridged = experiment_name.split('.')[1]
    return {f'{model_name}:{experiment_name_abridged}:{configuration_type}': dataset}

def load(model_names: str|list[str],
         experiment_names: str|list[str],
         field_names: str|list[str],
         zonal_average: bool=True,
         time_average: bool=True,
         pressure_level: int|float|None=500,
         diagnostic: bool=False):

    if isinstance(model_names, str): model_names = [model_names]
    if isinstance(experiment_names, str): experiment_names = [experiment_names]

    data_type = 'atmos_month'

    # Define the configuration names
    configuration_names = {'CTL': [f'CTL={model_name}:CTL1990.{experiment_name}' for experiment_name in experiment_names for model_name in model_names],
                           'EXP': [f'EXP={model_name}:CTL1990_SWISHE.{experiment_name}' for experiment_name in experiment_names for model_name in model_names]}
    configuration_names = configuration_names['CTL'] + configuration_names['EXP']

    ''' Parallel processing. '''
    # Maximum number of processors for computation
    max_number_procs = 24
    # Specify number of processors to use
    number_procs = len(configuration_names) if len(configuration_names) < max_number_procs else max_number_procs
    # Define partial function to use in parallel processing method
    partial_get_experiment = functools.partial(get_experiment,
                                               field_names,
                                               data_type,
                                               zonal_average,
                                               time_average,
                                               pressure_level)
    # Load data in parallel
    with Pool(processes=number_procs) as pool:
        container = pool.map(partial_get_experiment, configuration_names)
        pool.close()
    # Assign distributed data to singular dictionary
    data = {}
    for entry in container: data.update(entry)

    # Iterate over all configurations to get differences for configuration pairs
    for model_name in model_names:
        for experiment_name in experiment_names:
            # Find matching pairs of control-experiment configuration names
            configuration_keys = sorted([k for k in data.keys() if model_name in k and experiment_name in k])
            assert len(configuration_keys) == 2, f'Only two entries are allowed for differencing. Number of entries found: {len(configuration_keys)}.'
            # Define configuration difference name and get simple difference
            config_difference = (data[configuration_keys[0]] - data[configuration_keys[1]])
            config_base_name = ':'.join(configuration_keys[0].split(':')[0:2])
            config_difference_name = f'{config_base_name}:DIFF.CTL-EXP' # 
            data[config_difference_name] = config_difference

    return data

def plot_zonal_mean(dataset: dict,
                    config_names: list[str]=None,
                    field_names: str | list[str]='WVP',
                    normalize: bool=False,
                    rolling_degrees: int=10,
                    colormap_name: str='viridis',
                    label_start_index: int=0,
                    savefig: bool=False,
                    diagnostic: bool=False):

    diagnostic_tag = '[plot_zonal_mean()]'

    # Ensure the `field_name` variables is in a list data type
    if isinstance(field_names, str): field_names = [field_names]
    elif isinstance(field_names, list): field_names = field_names
    principal_field_name = field_names[0]
    nonprincipal_field_names = [field for field in field_names if field != principal_field_name]

    # Determine rolling mean parameters
    degrees_per_grid_point = 0.5 # grid cell size, in degrees
    rolling_points = np.round(rolling_degrees / degrees_per_grid_point).astype(int) # number of grid cells to average over

    # Metadata labeling - using the principal (or only) field name
    long_name, units = visualization.field_properties(field_names[0])

    # Get configuration names to plot
    config_names = dataset.keys() if config_names is None else config_names

    # Get model and experiments names
    # Each unique model will be printed in its own row
    # Each unique experiment will be printed in its own column
    model_names = list(sorted(set([config_name.split(':')[0] for config_name in config_names]))) 
    experiment_names = list(set([config_name.split(':')[1] for config_name in config_names]))

    # MANUAL OVERRIDE FOR CUSTOM ORDER, DELETE LATER
    experiment_names = ['CONST', '0N', '15N', 'TIMEVAR', 'AMIP']

    # Factors to apply 
    
    ''' Initialize plot. '''
    
    nrows, ncols = len(model_names), len(experiment_names)
    dpi = 300 if savefig else 144
    fig, gs = [plt.figure(figsize=(1.5 * ncols, 2.5 * nrows), dpi=dpi),
               matplotlib.gridspec.GridSpec(nrows=nrows, ncols=ncols, hspace=0.375, wspace=0.75)]
    axes = {}

    # Define colors to use
    colors = visualization.get_colormap_samples(number_of_samples=len(field_names),
                                                colormap_name=colormap_name)

    # Initialize dictionaries to hold model-specific extrema
    vmin = {model_name: np.nan for model_name in model_names}
    vmax = {model_name: np.nan for model_name in model_names}
    
    # Initialize a dictionary to hold configuration-specific extrema
    extrema = {} 
    
    ''' Pass 1: acquire data properties and plot. '''
    TMP = {config_name: {} for config_name in config_names}
    for index in range(len(config_names)):

        # Get configuration properties and make sure configuration names conform to structure and filter rules
        config_name = config_names[index]
        assert len(config_name.split(':')) == 3, 'Revisit configuration name definitions such that the configuration is of form {model_name}:{experiment_name}:{experiment_type}.'
                       
        # Initialize configuration-specific extrema entry
        extrema[config_name] = {}
        # Determine configuration-specific model and experiment name
        config_model_name, config_experiment_name = [config_name.split(':')[0],
                                                     config_name.split(':')[1]]
        
        nrow, ncol = [list(model_names).index(config_model_name),
                      list(experiment_names).index(config_experiment_name)]
        
        ax = fig.add_subplot(gs[nrow, ncol])
        # Reference line for a zonal mean
        ax.axhline(0, c='k', lw=0.5)

        for field_index, field_name in enumerate(field_names):
            # Get zonal-time mean
            zonal_mean_TMP = dataset[config_name][field_name].mean('grid_xt') if 'grid_xt' in dataset[config_name][field_name].dims else dataset[config_name][field_name]
            TMP[config_name][field_name] = zonal_mean_TMP.rolling(grid_yt=rolling_points, center=True).mean('grid_yt')
            # Add extrema pre-normalization
            extrema[config_name][field_name] = {'min': zonal_mean_TMP.min(), 'mean': zonal_mean_TMP.mean(), 'max': zonal_mean_TMP.max()}

        for field_index, field_name in enumerate(field_names):
            # Derive parameters for plotting format
            linestyle, linewidth, linecolor, alpha = ('-', 2, 'k', 1) if field_name == principal_field_name else ('--', 1, colors[field_index], 0.75)
            # Apply normalization, if chosen
            field_normalizer = abs(TMP[config_name][principal_field_name]).max() if field_name == principal_field_name else max([abs(TMP[config_name][field]).max() for field in nonprincipal_field_names])
            # Normalize non-principal fields by the same denominator
            normalized_TMP = TMP[config_name][field_name] / field_normalizer if normalize else TMP[config_name][field_name]
            # Plot zonal-time mean
            ax.plot(normalized_TMP, TMP[config_name][field_name].grid_yt, 
                    linestyle=linestyle, linewidth=linewidth, c=linecolor, alpha=alpha)

            # Print the global mean if the field name is the principal and if the configuration is a difference
            # if 'DIFF' in config_name and field_name == principal_field_name:
            #     print(f'Global mean difference of {principal_field_name} for configuration {config_name}: {TMP[config_name][field_name].sum().item():.2f}')
    
            axes[config_name] = {'ax': ax, 'nrow': nrow, 'ncol': ncol}
            # Collect configuration/field-specific extrema
            vmin[config_model_name] = TMP[config_name][field_name].min() if TMP[config_name][field_name].min() < vmin[config_model_name] or np.isnan(vmin[config_model_name]) else vmin[config_model_name]
            vmax[config_model_name] = TMP[config_name][field_name].max() if TMP[config_name][field_name].max() > vmax[config_model_name] or np.isnan(vmax[config_model_name]) else vmax[config_model_name]
    del TMP

    ''' Pass 2: Set model-specific x-axis limits (with some padding). '''
    for model_name in model_names:    
        x_axis_range = vmax[model_name] - vmin[model_name]
        x_axis_padding = 0.1 * x_axis_range
        vmin[model_name], vmax[model_name] = vmin[model_name] - x_axis_padding, vmax[model_name] + x_axis_padding
    
    ''' Pass 3: now that data is acquired, apply axis formatting. '''
    for _, config_name in enumerate(config_names):
        ax = axes[config_name]['ax']
        nrow, ncol = axes[config_name]['nrow'], axes[config_name]['ncol']
        model_name, experiment_names = [config_name.split(':')[0], 
                                        (f'CTL1990.{config_name.split(':')[1]}', f'CTL1990_SWISHE.{config_name.split(':')[1]}')]
        
        index = ncol + nrow * ncols
        
        ax.set_xlim([-1.1, 1.1]) if normalize else ax.set_xlim([vmin[model_name], vmax[model_name]])
        # Plot a horizontal reference line if extrema are of opposite
        if vmin[model_name] < 0 and vmax[model_name] > 0:
            ax.axvline(0, c='k', lw=0.5, zorder=0)

        # Axis metadata formatting
        if axes[config_name]['ncol'] != 0: ax.set_yticklabels([])

        ax.yaxis.set_major_locator(matplotlib.ticker.AutoLocator())
        ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())

        ax.set_ylim([-90 + rolling_degrees/2, 90 - rolling_degrees/2])

        ''' Sample configuration TC activity histogram. '''
    
        if diagnostic:
            print(f'{diagnostic_tag} Configuration name: {config_name} --> model name: {model_name}; experiment_name: {experiment_name}')

        # Acquire TC frequency density histograms
        get_density_difference = True if 'DIFF' in ''.join(config_names) else False
        if diagnostic:
            print(f'{diagnostic_tag} density difference: {get_density_difference}; model name: {model_name}; experiment names: {experiment_names}')
        get_TC_track_histogram(ax=ax, model_name=model_name, experiment_names=experiment_names, get_density_difference=get_density_difference)

        # Ensure axis spines stay on top
        for k, spine in ax.spines.items(): spine.set_zorder(100)

        # Prevent overcrowding of x-axis tick labels
        if len(ax.get_xticklabels()) % 2 != 0:
            [l.set_visible(False) for (i, l) in enumerate(ax.xaxis.get_ticklabels()) if i % 2 == 0]

        # Add subplot title
        ax_title = f'{visualization.get_alphabet_letter(index + label_start_index)}) {'-'.join(config_name.split(':')[0:2])}'
        ax_title_spacing = 1.075 if 'DIFF' in ''.join(config_names) else 1
        ax.set_title(ax_title, fontsize=9, y=ax_title_spacing)
        # Add subplot annotation with extrema, if normalization applied
        if normalize:
            # Get statistics
            config_field_min, config_field_mean, config_field_max = [extrema[config_name][principal_field_name]['min'], # minimum
                                                                     extrema[config_name][principal_field_name]['mean'], # mean
                                                                     extrema[config_name][principal_field_name]['max']] # maximum
            statistics_annotation_text = f'{config_field_mean:.2f} ({config_field_min:.2f}, {config_field_max:.2f})'
            statistics_annotation = ax.annotate(statistics_annotation_text, (0.5, 1.025), xycoords='axes fraction', 
                                                 ha='center', va='bottom', fontsize=7)
            statistics_annotation.set_path_effects([matplotlib.patheffects.Stroke(linewidth=1.5, foreground='white'),
                                                    matplotlib.patheffects.Normal()])

    ''' Legend definition. '''

    fig.tight_layout()
    # Add a conditional prepend string in case the experiment is a difference
    label_prepend = '$\delta$' if 'diff' in ''.join(config_names).lower() else ''
    supxlabel = f'{label_prepend}{get_field_name_alias(principal_field_name)} [{units}]'
    if normalize: 
        supxlabel = f'$\\parallel${label_prepend}{get_field_name_alias(principal_field_name)}$\\parallel$'
    else:
        supxlabel = f'{label_prepend}{get_field_name_alias(principal_field_name)} [{units}]'
    fig.supxlabel(supxlabel, 
                  y=0, fontsize=10, wrap=True)
    
    fig.supylabel('Latitude', x=0.025, fontsize=10, wrap=True)

    if len(field_names) > 1:
        labels = [f'{label_prepend}{get_field_name_alias(field_name)}'for field_name in field_names]
        handles = [matplotlib.lines.Line2D([0], [0], color='k', lw=2) 
                   if field_name == principal_field_name 
                   else matplotlib.lines.Line2D([0], [0], color=colors[index], lw=2) 
                   for index, field_name in enumerate(field_names)]
        fig.legend(handles, labels, loc='upper left', frameon=False, bbox_to_anchor=(0.95, 0.9))

    if savefig:
        storage_dirname = '/projects/GEOCLIM/gr7610/analysis/TC-AQP/figs'
        figure_identifier = 'zonal_mean'
        figure_config_type = 'difference' if 'diff' in ''.join(config_names).lower() else 'climo'
        figure_config_names = sorted(list(set([name.split(':')[1] for name in config_names])))
        filename = f'{figure_identifier}.{figure_config_type}.configuration_names.{'-'.join(figure_config_names)}.field_name.{principal_field_name}.pdf'

        storage_pathname = os.path.join(storage_dirname, filename)
        plt.savefig(storage_pathname, transparent=True, dpi=dpi, format='pdf', bbox_inches='tight')

def get_TC_track_histogram(ax: matplotlib.axes.Axes,
                           model_name: str,
                           experiment_names: str|list[str],
                           colors: list[str],
                           bin_size: int|float=2.5,
                           get_density_difference: bool=True,
                           orientation: str='horizontal',
                           override_difference_colors: bool=False) -> dict:

    ''' Get histograms of TC track frequency and location. Allows individual densities or density differences to be generated. '''

    # Gets the difference in TC track density between two experiemnts
    if get_density_difference:
        
        assert len(experiment_names) == 2, 'Two experiment names must be provided to generate a difference. '''
        density_TCs = visualization.TC_density_histogram(ax=None,
                                                         model_name=model_name,
                                                         experiment_names=experiment_names,
                                                         bin_size=bin_size,
                                                         axis_depth=0.2,
                                                         orientation=orientation,
                                                         basin_name=None,
                                                         diagnostic=False)
            
        density_TCs = {'DIFF': {'counts': density_TCs[experiment_names[0]]['counts'] - density_TCs[experiment_names[1]]['counts'],
                                'bins': density_TCs[experiment_names[0]]['bins']}}

        # Determine difference colors
        colors = colors if override_difference_colors else 'tab:gray'
        
        _ = visualization.TC_density_histogram(ax=ax,
                                               model_name=model_name,
                                               experiment_names=experiment_names,
                                               density_data=density_TCs,
                                               bin_size=bin_size,
                                               axis_depth=0.2,
                                               orientation=orientation,
                                               basin_name=None,
                                               colors=colors,
                                               diagnostic=False)

    else:
        density_TCs = visualization.TC_density_histogram(ax=ax,
                                                         model_name=model_name,
                                                         experiment_names=experiment_names,
                                                         bin_size=bin_size,
                                                         axis_depth=0.2,
                                                         orientation=orientation,
                                                         colors=colors,
                                                         basin_name=None,
                                                         diagnostic=False)

    return density_TCs