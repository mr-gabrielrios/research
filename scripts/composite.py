import cftime, datetime
import multiprocessing
import matplotlib, matplotlib.pyplot as plt
import numpy as np, pandas as pd, scipy as sp, xarray as xr
import os, pickle, random

import concurrent.futures
import itertools

import importlib
import composite_weighting, rotate, utilities, visualization
importlib.reload(composite_weighting)

def azimuthal_composite(data, model_names, field, intensity_bin=None, pressure_level=None, experiments=['CTL1990s', 'CTL1990s_swishe'], 
                        track_data=None, weighting_intensity_metric=None, dpi=96, parallel=False, diagnostic=False):
    
    # Define experiment names as variables
    run_control, run_experiment = experiments
    run_difference = '{1}-{0}'.format(run_control, run_experiment)
    
    if not weighting_intensity_metric:
        print('[visualization.py, azimuthal_composite_2d()] Default weighting intensity metric being defined.')
        weighting_intensity_metric = 'min_slp'
    
    # Initialize counter dictionary to hold data statistics for each model and corresponding experiment
    counts = {model_name: {run_control: {}, run_experiment: {}} for model_name in model_names}
    
    # Initialize dictionary to hold processed data for each model and corresponding experiment
    individual_storms = {model_name: {run_control: {}, run_experiment: {}} for model_name in model_names}
    
    # Iterate over each model and experiment to get the azimuthal means for a given field
    for model in individual_storms.keys():
        for experiment in individual_storms[model].keys():
            print('[visualization.py, azimuthal_composite_2d()] Iterating over model {0} for experiment {1}'.format(model, experiment))
            individual_storms[model][experiment], counts[model][experiment] = compositor_preprocessor(model, data[model][experiment]['data'], intensity_bin, field, pressure_level=pressure_level,
                                                                                                      compositing_mode='azimuthal', track_data=track_data, weighting_intensity_metric=weighting_intensity_metric)
            # Populate the dictionaries for the given intensity bin, correct for factors or signs before loading
            individual_storms[model][experiment] = utilities.field_correction(data=individual_storms[model][experiment], field=field)   
            
    # Iterate over each model and experiment to get the azimuthal means for a given field
    for model in individual_storms.keys():
        for experiment in individual_storms[model].keys():        
            # Populate the dictionaries for the given intensity bin, correct for factors or signs before loading
            individual_storms[model][experiment] = utilities.field_correction(data=individual_storms[model][experiment], field=field)   
    
    
    for model in individual_storms.keys():
        individual_storms[model][run_difference] = individual_storms[model][run_experiment] - individual_storms[model][run_control]
    
    return individual_storms, counts

def planar_composite(data, model_names, intensity_bin, field, pressure_level, 
                     experiments=None, track_data=None, weighting_intensity_metric=None,
                     contour_levels=16, dpi=96, lon_nudge=0, lat_nudge=0, 
                     rotation=True, inline_statistics=False, subplot_format='circular', plot_style='contourf', diagnostic=False):
    """
    Method to obtain planar composite for prescribed data for a given field.
    Note: input must be in the form of a 3-tiered dictionary that is the output of tc_analyis.tc_model_data().

    Args:
        data (dict): 3-tiered dictionary
        model_names (list of str): list of strs with model names
        intensity_bins (list of str): list of strs with intensity bin listings
        field (str): name of field to evaluate
        pressure_level (numeric): level at which to evaluate field
        experiment_plot (str): experiment to plot (typically control, swishe, or swishe-control)
        dpi (int): dots per inch - 96 is default, 300 is recommended for publication-quality
        inline_statistics (bool): boolean to control wheether or not statistics are published in the figures
    """
 
    # Determine which experiments to evaluate; define a difference plot if 2 experiment names provided
    if len(experiments) == 2:
        # Define control and experiment names (assume both are in the density dictionary)
        run_control, run_experiment = experiments[0], experiments[1]
        run_difference = '{0}-{1}'.format(run_experiment, run_control)
        # Define experiment list
        experiments = [run_control, run_experiment, run_difference]
        
    ''' Generate composites. '''
    # Initialize counter dictionary to hold data statistics for each model and corresponding experiment
    counts = {model_name: {run_control: {}, run_experiment: {}} for model_name in model_names}
    # Iterate over the models provided in the dictionary
    for model in model_names:
        # Iterate over the experiments provided in the dictionary
        for experiment in data[model].keys():
            print('[visualization.py, planar_composite()] Processing model {0}, experiment: {1}'.format(model, experiment))
            # Create subdictionaries for prescribed intensity bins
            data[model][experiment]['composite'] = None
            # Determine experiments to undergo weighting
            composite_mean, counts[model][experiment] = compositor_preprocessor(model, data[model][experiment]['data'], intensity_bin=intensity_bin, field=field,
                                                                                pressure_level=pressure_level, compositing_mode='planar', track_data=track_data, weighting_intensity_metric='min_slp')
            # Populate the dictionaries for the given intensity bin, correct for factors or signs before loading
            data[model][experiment]['composite'] = utilities.field_correction(data=composite_mean, field=field)
            
    # Collect composites. Note: the differences will be calculated after this loop series.
    # Note: composite dictionary is structured top-down as: (1) intensity bin, (2) model, (3) experiment
    composites = {}
    # Pass 1: Iterate through experiments to get each experiment's density data.
    for model in model_names:
        # Initialize subdictionary for the iterand model
        composites[model] = {}
        # Iterate over each given experiment
        for experiment in data[model].keys():
            composites[model][experiment] = data[model][experiment]['composite']
                
    # Pass 2: Iterate through experiments to get the differences if a difference plot type (anything with a subtraction, or '-', is detected in the "experiment_plots" argument)
    for model in composites.keys():
        try:
            # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
            delta = composites[model][run_experiment].values - composites[model][run_control].values
            # Populate the dictionary with reconstructed xArray
            composites[model][run_difference] = xr.DataArray(data=delta, dims=['grid_yt', 'grid_xt'],
                                                                        coords={'grid_yt': (['grid_yt'], 
                                                                                            composites[model][run_control]['grid_yt'].values), 
                                                                                'grid_xt': (['grid_xt'], 
                                                                                            composites[model][run_control]['grid_xt'].values)})
        except:
            # Get raw difference (this is done the brute force way due to artificial xArray mismatches despite proper dimension alignment)
            # Methodology: trim to common centers of data. In other words, find the smaller composite and trim the bigger one to their dimensions.
            # Note: I wrote this stoned, so it can likely be optimized
            
            # Get shape
            shape_control, shape_exp = composites[model][run_control].shape, composites[model][run_experiment].shape
            # Round decimal places for the coordinates (floating point precision errors lead to weird indexing)
            composites[model][run_control]['grid_xt'] = composites[model][run_control]['grid_xt'].round(3)
            composites[model][run_control]['grid_yt'] = composites[model][run_control]['grid_yt'].round(3)
            composites[model][run_experiment]['grid_xt'] = composites[model][run_experiment]['grid_xt'].round(3)
            composites[model][run_experiment]['grid_yt'] = composites[model][run_experiment]['grid_yt'].round(3)
            
            # Compare latitudes:
            if shape_control[0] < shape_exp[0]:
                temp = composites[model][run_experiment].where((composites[model][run_experiment]['grid_yt'] >= composites[model][run_control]['grid_yt'].min())
                                                                            & (composites[model][run_experiment]['grid_yt'] <= composites[model][run_control]['grid_yt'].max()), drop=True)
                del composites[model][run_experiment]
                composites[model][run_experiment] = temp
            else:
                temp = composites[model][run_control].where((composites[model][run_control]['grid_yt'] >= composites[model][run_experiment]['grid_yt'].min())
                                                                            & (composites[model][run_control]['grid_yt'] <= composites[model][run_experiment]['grid_yt'].max()), drop=True)
                del composites[model][run_control]
                composites[model][run_control] = temp
            # now longitudes
            if shape_control[1] < shape_exp[1]:
                temp = composites[model][run_experiment].where((composites[model][run_experiment]['grid_xt'] >= composites[model][run_control]['grid_xt'].min())
                                                                            & (composites[model][run_experiment]['grid_xt'] <= composites[model][run_control]['grid_xt'].max()), drop=True)
                del composites[model][run_experiment]
                composites[model][run_experiment] = temp
            else:
                temp = composites[model][run_control].where((composites[model][run_control]['grid_xt'] >= composites[model][run_experiment]['grid_xt'].min())
                                                                            & (composites[model][run_control]['grid_xt'] <= composites[model][run_experiment]['grid_xt'].max()), drop=True)
                del composites[model][run_control]
                composites[model][run_control] = temp
            # Crude fix written when I was stoned, but it worked
            x_, y_ = composites[model][run_control].values, composites[model][run_experiment].values
            if x_.shape == y_.shape:
                delta = y_ - x_
                grid_xt = composites[model][run_control]['grid_xt'].values
                grid_yt = composites[model][run_control]['grid_yt'].values
            else:
                pass
                print('[visualization.py, planar_composite()] intake: sample A shape = {0}; sample B shape = {1}'.format(x_.shape, y_.shape))
                if x_.shape[0] > y_.shape[0]:
                    print('[visualization.py, planar_composite()] sample A y-axis > sample B y-axis')
                    # if first sample's y-axis is larger
                    x_ = np.vstack((np.full((1, x_.shape[1]), np.nan), x_[:y_.shape[0], :], np.full((1, x_.shape[1]), np.nan)))
                    grid_yt = np.concatenate(([composites[model][run_experiment].grid_yt.values[0]], composites[model][run_control].grid_yt.values, [composites[model][run_experiment].grid_yt.values[-1]]))
                elif y_.shape[0] > x_.shape[0]: 
                    print('[visualization.py, planar_composite()] sample B y-axis > sample A y-axis')
                    # if second sample's x-axis is larger
                    y_ = y_[:x_.shape[0], :]
                    grid_yt = composites[model][run_experiment].grid_yt.values[:x_.shape[0]]
                    # y_ = np.vstack((np.full((1, y_.shape[1]), np.nan), y_[:x_.shape[0], :], np.full((1, y_.shape[1]), np.nan)))
                    # grid_yt = np.concatenate(([composites[model][run_control].grid_yt.values[0]], composites[model][run_experiment].grid_yt.values, [composites[model][run_control].grid_yt.values[-1]]))
                else:
                    grid_yt = composites[model][run_control].grid_yt.values
                print('[visualization.py, planar_composite()] post y-axis filtering: sample A shape = {0}; sample B shape = {1}'.format(x_.shape, y_.shape))
                
                if x_.shape[1] < y_.shape[1]:
                    print('[visualization.py, planar_composite()] sample A y-axis > sample B y-axis')
                    # if first sample's y-axis is larger
                    x_ = np.hstack((np.full((x_.shape[0], 1), np.nan), x_[:, :y_.shape[1]], np.full((x_.shape[0], 1), np.nan)))
                    grid_xt = np.concatenate(([composites[model][run_experiment].grid_xt.values[0]], composites[model][run_control].grid_xt.values, [composites[model][run_experiment].grid_xt.values[-1]]))
                elif x_.shape[1] > y_.shape[1]:
                    print('[visualization.py, planar_composite()] sample A y-axis > sample B y-axis')
                    # if second sample's y-axis is larger
                    y_ = np.vstack((np.full((y_.shape[0], 1), np.nan), y_[:, :x_.shape[1]], np.full((y_.shape[0], 1), np.nan)))
                    grid_xt = np.concatenate(([composites[model][run_control].grid_xt.values[0]], composites[model][run_experiment].grid_xt.values, [composites[model][run_control].grid_xt.values[-1]]))
                else:
                    grid_xt = composites[model][run_control].grid_xt.values
                print('[visualization.py, planar_composite()] post x-axis filtering: sample A shape = {0}; sample B shape = {1}'.format(x_.shape, y_.shape))
                    
                if x_.shape != y_.shape:
                    continue
                else:
                    delta = y_ - x_
            
            # Populate the dictionary with reconstructed xArray
            composites[model][run_difference] = xr.DataArray(data=delta, dims=['grid_yt', 'grid_xt'],
                                                             coords={'grid_yt': (['grid_yt'], grid_yt),
                                                                     'grid_xt': (['grid_xt'], grid_xt)})
             
    return composites, counts
        
def azimuthal_processor(data, field='wind_tangential', pressure_level=None, visual_check=False, diagnostic=False):
    """
    Method to composite data azimuthally for a single storm.

    Args:
        data (dict): xArray DataArray for a single storm ID
        field (str, optional): _description_. Defaults to 'sphum'.
        visuals (bool, optional): _description_. Defaults to False.

    Returns:
        composite_azim (xArrray DataArray): _description_
    """

    ''' Methodology.
    1. Create radius field from storm-centered coordinates.
    2. Define grid resolution to define annular thickness.
    3. Define the starting radius to initialize annulus
    4. Create annular basis vector to determine where the annular radii will be defined.
    5. Initialize output array, with dimensions of (pfull, basis_vector).
    6. Create annulus, whose area will be averaged over.
    7. Average over annular area.
    8. Expand annulus by a grid cell and repeat 6-8 until domain edge is reached.
    '''
    
    if diagnostic:
        print('[composite.py, azimuthal_compositor()] Azimuthally processing {0}...'.format(data.attrs['storm_ID']))
    
    # Clean nans from the input data, if any exist
    dataset = data.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all') 
    visuals = {}
    
    # print(dataset['slp'].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all') .values)
    pressure_field = dataset['slp'].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')

    # Check to see if storm-centered coordinate system exists. If not, build it.
    if 'TC_xt' not in dataset.coords:
        # Get new midpoints for trimmed data 
        center_x, center_y = len(dataset.grid_xt) // 2, len(dataset.grid_yt) // 2
        print('From {1}, {0}...'.format(center_x, center_y))
        center_x, center_y = np.unravel_index(np.argmin(pressure_field.values), pressure_field.shape)
        print('... to {1}, {0}... to'.format(center_x, center_y))

        # Create coordinates for a TC-centric coordinate system (TC_xt, TC_yt)
        dataset = dataset.assign_coords({'TC_xt': dataset['grid_xt'] - dataset['grid_xt'].isel(grid_xt=center_x),
                                         'TC_yt': dataset['grid_yt'] - dataset['grid_yt'].isel(grid_yt=center_y)})

    # 1. Create radial field relative to storm center.
    X, Y = np.meshgrid(dataset['TC_xt'], dataset['TC_yt'])
    dataset['radius'] = xr.DataArray(data=np.sqrt(X**2 + Y**2), dims=('grid_yt', 'grid_xt'))
    # 2. Define grid resolution from coarser of the two domain dimensions ('TC_xt' or 'TC_yt')
    # Note: because the spatial extent is square, the dimension with fewer elements is coarser by definition
    limiting_dimension = 'TC_xt' if len(dataset['TC_xt']) < len(dataset['TC_yt']) else 'TC_yt'
    resolution = 1.5*np.diff(dataset[limiting_dimension])[0] # assume all elements are equal
    # 3. Create the first annulus inner radius and the radial limit
    r_i = 0
    r_limit = dataset[limiting_dimension].max().item()
    # 4. Create a basis vector over the radius to guide the annular expansion
    # Note: if the TC_xt/TC_yt storm-centered domain isn't centered on the storm for whatever reason, re-define the coordinate system
    #       This is taken from derived.py --> radial_tangential_velocities()
    if r_limit < 0:
        # Get new midpoints for trimmed data 
        center_x, center_y = len(dataset.grid_xt) // 2, len(dataset.grid_yt) // 2
        print('From {1}, {0}...'.format(center_x, center_y))
        center_x, center_y = np.unravel_index(np.argmin(pressure_field.values), pressure_field.shape)
        print('... to {1}, {0}... to'.format(center_x, center_y))
        # Create coordinates for a TC-centric coordinate system (TC_xt, TC_yt)
        dataset = dataset.assign_coords({'TC_xt': dataset['grid_xt'] - dataset['grid_xt'].isel(grid_xt=center_x),
                                         'TC_yt': dataset['grid_yt'] - dataset['grid_yt'].isel(grid_yt=center_y)})
        r_i = 0
        r_limit = dataset[limiting_dimension].max().item()
    rs = np.arange(r_i, r_limit + resolution, resolution)
    # 5. Create output array
    if not pressure_level and 'pfull' in dataset.dims:
        if diagnostic:
            print('[composite.py, azimuthal_compositor()] 2D data (radial x vertical) being processed...')
        # Handle data with a vertical component only (2D, radius and pressure level)
        out = np.full(shape=(len(dataset.pfull.values), len(rs)-1), fill_value=np.nan)
        # 6. Create the annulus using a loop based on the limiting dimension and resolution
        # Note: use a basis vector from the center (0) to the domain edge
        # Note: begin from index 1 to establish the outer annular radius
        for index, r_o in enumerate(rs[1:]):
            if diagnostic:
                print('[composite.py, azimuthal_compositor()] Annulus #{0}: inner radius = {1:.2f}, outer radius = {2:.2f}'.format(index, r_i, r_o))
        
            # 7. Average over annular area.
            average = dataset[field].where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o)).mean(dim='grid_xt').mean(dim='grid_yt')
            out[:, index] = average
            # Append this output to a storage list for easy debugging, if needed.
            visuals['Annulus #{0}:\ninner rad. = {1:.2f}, out. rad. = {2:.2f}'.format(index, r_i, r_o)] = \
            dataset.where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o)).sel(pfull=850, method='nearest')
            
            # 8. Expand annulus by a grid cell and repeat steps 6-8 until domain edge is reached.
            r_i += resolution

        # 9. Create xArray DataArray for the composite
        composite_azim = xr.DataArray(data=out, 
                                        dims=('pfull', 'radius'), 
                                        coords={'pfull': (['pfull'], dataset.pfull.values),
                                              'radius': (['radius'], rs[:-1])},
                                        attrs={'storm_ID': dataset.attrs['storm_ID'],
                                               'slp': dataset['slp'].min().values})
    else:
        if diagnostic:
            print('[composite.py, azimuthal_compositor()] 1D data (radial) being processed...')
        # Handle data with only a radial component
        out = np.full(shape=(len(rs)-1), fill_value=np.nan)
        # 6. Create the annulus using a loop based on the limiting dimension and resolution
        # Note: use a basis vector from the center (0) to the domain edge
        # Note: begin from index 1 to establish the outer annular radius
        for index, r_o in enumerate(rs[1:]):
            if diagnostic:
                print('[composite.py, azimuthal_compositor()] Annulus #{0}: inner radius = {1:.2f}, outer radius = {2:.2f}'.format(index, r_i, r_o))
        
            # 7. Average over annular area.
            average = dataset[field].where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o)).mean(dim='grid_xt').mean(dim='grid_yt')
            out[index] = average
            # Append this output to a storage list for easy debugging, if needed.
            visuals['Annulus #{0}:\ninner rad. = {1:.2f}, out. rad. = {2:.2f}'.format(index, r_i, r_o)] = dataset.where((dataset['radius'] >= r_i) & (dataset['radius'] < r_o))
            
            # 8. Expand annulus by a grid cell and repeat steps 6-8 until domain edge is reached.
            r_i += resolution

        # 9. Create xArray DataArray for the composite
        composite_azim = xr.DataArray(data=out, dims=('radius'), coords={'radius': (['radius'], rs[:-1])}, attrs=dataset.attrs)
        
    if visual_check:
        if pressure_level:
            fig, ax = plt.subplots(figsize=(4, 2))
            composite_azim.plot(ax=ax)
        else:
            fig, ax = plt.subplots(figsize=(3, 4))
            norm, cmap = visualization.norm_cmap(composite_azim, field=field)
            cmap = visualization.diverging_colormap_adjust(norm, cmap, additional=1)
            composite_azim.plot.contourf(levels=16, ax=ax, norm=norm, cmap=cmap)
            ax.set_ylim(ax.get_ylim()[::-1])

    composite_azim['time'] = dataset['time']

    return composite_azim

def compositor_preprocessor(model, datasets, intensity_bin, field, pressure_level=None, experiment_designator=False, rotation=False,
                            compositing_mode='azimuthal', track_data=None, weighting_intensity_metric=None, NH_TCs=False, diagnostic=False):
    
    # Define a sample decrement to subtract filtered data from the output data array for counts
    sample_count_decrement = 0
    # Define key with format {FIELD}{HEIGHT} to match GFDL diagnostics names.
    # For example, specific humidity at 700 hPa would be 'q700'.
    key = '{0}{1}'.format(field, pressure_level) if pressure_level else field
    # Grab minimum extents for the datasets for future trimming based on longitude (grid_xt) and latitude (grid_yt) extents
    min_x, min_y = {}, {}
    # Define buffer to prevent edge effects of domain trimming. Will be applied to each edge, such that the final domain will be 2*buffer less after trimming.
    buffer = 1
    
    # Initialize dictionary to hold metadata for future composite weighting
    weighting_metadata = {}
    
    ''' 
    Begin step 1. 
    '''
    
    ################################################################################################################
    # Iterate over each storm entry (iterates over each entry i in list data[MODEL][EXPERIMENT]['data'][i]). 
    # This step collects relevant data from each storm belonging to the input dataset.
    ################################################################################################################
    
    max_procs = 16
    num_procs = int(len(datasets)/4) if len(datasets) < max_procs else max_procs
    print('[composite.py, compositor_preprocessor()] Parallelizing on {0} processors...'.format(num_procs))
    # Prepare inputs for snapshot collection
    inputs = [[model, key, ds, field, pressure_level, intensity_bin, min_x, min_y, compositing_mode, track_data, weighting_intensity_metric, diagnostic]
              for ds in datasets]
    with multiprocessing.get_context("spawn").Pool(num_procs) as p: 
        storm_data = [result for result in p.starmap(snapshot_collection, inputs)]
    del inputs
    
    # Collect all model output associated with each storm ID and concatenate into a dictionary
    # Initialize dictionary to hold data for each storm corresponding to the intensity bin and field of choice
    
    # Auxiliary function to avoid repetition
    def item_input(entry, varname, vardict):
        if len(entry[varname]) == 1:
            key = list(entry[varname].keys())[0]
            vardict[key] = list(entry[varname].values())[0]
        elif len(entry[varname]) > 1:
            for k, v in entry[varname].items():
                vardict[k] = v
            if diagnostic:
                print('[composite.py, compositor_preprocessor()] More than 1 entry found after collecting snapshots.')
        else:
            if diagnostic:
                print('[composite.py, compositor_preprocessor()] Less than 1 entry found after collecting snapshots.')
            
        # Initialize output dictionary
        result = {}
        # Ensure unique values
        for k, v in vardict.items():
            if k not in result.keys():
                result[k] = v
            
        return result
    
    # Create collection data structure for all outputs from the loop
    collection_keys = ['data', 'weighting_metadata', 'min_x', 'min_y']
    collection =  {k: {} for k in collection_keys}

    # Iterate over the outputs
    for entry in storm_data:
        for varname in collection_keys:
            collection[varname] = item_input(entry, varname, collection[varname])
            
    data, weighting_metadata, min_x, min_y = collection['data'], collection['weighting_metadata'], collection['min_x'], collection['min_y']

    # Get minima for the longitude and latitude axes of the data
    min_x, min_y = np.nanmin(np.fromiter(min_x.values(), dtype=int)), np.nanmin(np.fromiter(min_y.values(), dtype=int))
    if diagnostic:
        print('[composite.py, compositor_preprocessor()] grid_xt minimum: {0}; grid_yt minimum: {1}'.format(min_x, min_y))
    
    ''' End step 1. '''
    
    ''' Begin step 2. '''

    ################################################################################################################
    # SUMMARY OF STEP
    ################################################################################################################
    
    if diagnostic:
        print('[composite.py, planar_compositor()] Weighting intensities: {0}'.format(weighting_metadata))
    # Get values for composite distribution weighting based on collected values in `weighting_metadata`
    if track_data and weighting_intensity_metric:
        if diagnostic:
            print('[composite.py, planar_compositor()] intensity weighting undergoing...')
        weighting_input_experiments = list(track_data[model].keys())
        weighting_inputs = {'CTL': {weighting_intensity_metric: track_data[model][weighting_input_experiments[0]]['unique'][weighting_intensity_metric]},
                            'EXP': {weighting_intensity_metric: np.array(list(weighting_metadata.values()))}}
        weighting_bins, weights, _ = composite_weighting.distribution_weighting(weighting_inputs, weighting_intensity_metric, diagnostic=False)
        # If the experiment is considered an experiment rather than a control, invert the weights for better normalization
        if experiment_designator:
            weights = {k: 1/v for k, v in weights.items()}
    else:
        weighting_bins, weights = None, None
            
    # Initialize dictionary to hold weights based on entries
    storm_weights = {}
    
    # Iterate over each storm and its sub-IDs and trim
    inputs = [[model, key, data, storm_id, field, pressure_level, intensity_bin, min_x, min_y, 
               weights, weighting_metadata, weighting_bins, rotation, 1, NH_TCs, diagnostic] 
              for storm_id in data.keys()]
    with multiprocessing.get_context("spawn").Pool(num_procs) as p: 
        data_weights = [result for result in p.starmap(snapshot_weight_trim, inputs)]
    del inputs
    
    # Create collection data structure for all outputs from the loop
    collection_keys = ['data', 'storm_weight']
    collection_initial = {k: {} for k in collection_keys}
    # Iterate over the outputs
    for entry in data_weights:
        for varname in collection_keys:
            collection_initial[varname] = item_input(entry, varname, collection_initial[varname])
            collection_initial[varname] = {k: v for k, v in collection_initial[varname].items()
                                   if v != np.isnan}
    
    collection = {k: {} for k in collection_keys}
    for varname in collection.keys():
        for k, v in collection_initial[varname].items():
            if not ('float' in str(type(v)) and np.isnan(v)):
                collection[varname][k] = v
           
    data, storm_weights = collection['data'], collection['storm_weight']
    
    # Get number of samples being processed for the composite mean (this is pre-weighting, where sample counts get modified)
    sample_count = len(data)
    print('Pre-weighting sample count: {0}'.format(sample_count))
    
    ''' End step 2. '''

    ''' Begin step 3. '''
    
    ################################################################################################################
    # SUMMARY OF STEP
    ################################################################################################################
    
    # Initialize list of snapshots, which will be stacked and averaged for a composite mean
    # The idea here is that if the experiment:control ratio is greater than 1, control entries are repeated some number of times. 
    # Otherwise, experiment entries are repeated by the inverse of the ratio.
    snapshot_stack = []
    # Iterate over each snapshot
    for k, v in data.items():
        # Use intensity-based weighting for each selected timestamp
        # Determine weight of the entry
        storm_weight = storm_weights[k]
        # Round to nearest defined interval to determine number of repetitions for the entry
        storm_weight_rounded = 0.25 * round(storm_weight/0.25)
        num_repetitions_rounded = int(storm_weight_rounded)
        print('[composite.py, compositor_preprocessor()] Number of repetitions for {0}: {1}'.format(k, num_repetitions_rounded))
        if diagnostic:
            print('[composite.py, planar_compositor()] Number of repetitions: {0}, interval-rounded: {1}, and fully rounded: {2}'.format(storm_weight, storm_weight_rounded, num_repetitions_rounded))
        # This loop controls the number of times a sample is added to the stack
        for repetition in range(num_repetitions_rounded):
            # v.attrs['snapshot_methodology'] = 
            snapshot_stack.append(v)
        
    
    # Get composite mean of all provided snapshots in 'snapshot_stack'
    if compositing_mode == 'planar':
        composite_mean = planar_composite_mean(data, snapshot_stack, field)
        return snapshot_stack, composite_mean, sample_count
    
    # Get azimuthal mean of all provided snapshots in 'snapshot_stack'
    elif compositing_mode == 'azimuthal':
        # Assemble the input lists for Pool.starmap
        inputs = [[snapshot, field, pressure_level, diagnostic, diagnostic] for snapshot in snapshot_stack]
        # Calculate azimuthal means in parallel
        print('[composite.py, compositor_preprocessor()] Processing in parallel over {0} processors'.format(num_procs))
        with multiprocessing.get_context("spawn").Pool(num_procs) as p:
            snapshots = [result for result in p.starmap(azimuthal_processor, inputs)]
            composite_mean = azimuthal_composite_mean(snapshots, outer_radius_index=20, diagnostic=diagnostic)
        
        return snapshots, composite_mean, sample_count
    
    else:
        print('[composite.py, planar_compositor()] Compositing mode not recognized - this script only handles planar and azimuthal composites.')
        return None, None
    
    # return key, data, len(data) + sample_count_decrement, composite_mean

def snapshot_collection(model, key, ds, field, pressure_level, intensity_bin, min_x, min_y, 
                        compositing_mode='azimuthal', track_data=None, weighting_intensity_metric=None, 
                        diagnostic=True):
        
    data = {} 
    weighting_metadata = {}
        
    # Get storm ID
    storm_id = ds['track_output']['storm_id'].unique().item()
    if diagnostic:
        print('[composite.py, snapshot_collection()] Processing {0}...'.format(storm_id))
    
    # Get timestamps relevant to the intensity bin at hand
    times = ds['track_output']['time'].loc[ds['track_output']['intensity_bin'] == intensity_bin]
    if diagnostic:
        if len(times) > 0:
            [print('[composite.py] Timestamps level 0:', t) for t in times]
        else:
            print('Intensity bin: {0}'.format(intensity_bin))
            print('Data sample: ', ds['track_output'])

    # Change from Pandas to cftime formats
    times = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in times if not ((t.month == 2) & (t.day == 29))]
    if diagnostic:
        [print('[composite.py] Timestamps level 1:', t) for t in times]
    
    # Determine which dictionary the field data is in
    subdict = None
    for sd in ['tc_model_output', 'tc_vertical_output']:
        # Remove duplicate entries along the time axis
        ds[sd] = ds[sd].drop_duplicates('time')
        
        # Wind-specific conditional: since 'wind' is a field in both types of datasets, handle manually based on pressure level input
        if field == 'wind':
            subdict = 'tc_vertical_output' if pressure_level is not None else 'tc_model_output'
        else:
            if field in ds[sd].data_vars:
                subdict = sd
    # Print error message
    if subdict is None:
        print('[composite.py, snapshot_collection()] Field not found in the data. Storm ID is: {0}. Continuing to the next storm...'.format(storm_id))
    
    # Check again for duplicate time values
    ds[subdict] = ds[subdict].drop_duplicates('time')
    # Get the matching timestamps between the track output and the corresponding dictionary with field data
    index_timestamps = list(sorted(set(ds['tc_model_output'].time.values) & set(ds['tc_vertical_output'].time.values)))

    if diagnostic:
        [print('[composite.py, snapshot_collection()] Timestamps level 2:', index_timestamp) for index_timestamp in index_timestamps]
    
    diagnostic = False
    # Initialize container dataset
    dataset = {}
    # Data merging
    if compositing_mode == 'planar':
        dataset = dataset_alignment(ds, subdict, field, times, pressure_level=pressure_level, 
                                    index_timestamps=index_timestamps, diagnostic=diagnostic)
    else:
        # Get data
        if pressure_level:
            dataset[field] = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') 
        else:
            dataset[field] = ds[subdict][field].sel(time=index_timestamps)
        # Get center data
        for planar_field in ['center_lon', 'center_lat', 'slp']:
            dataset[planar_field] = ds['tc_model_output'].sel(time=index_timestamps)[planar_field]
            
        # Merge dictionary into an xArray Dataset
        dataset = xr.merge(dataset.values())

    # If timestamps match, proceed. Else, continue.
    if len(dataset.time.values) > 0:
        if diagnostic:
            print('\t \t {0} has {1} matching timestamps! Further processing them now...'.format(storm_id, len(dataset.time.values)))
        
        # Methodology: extract all timestamps with matching intensity for the given storm and assign them a storm 'sub-ID', which preserves the storm ID but gives it a unique identifier.
        storm_subids = [] # create storage list for all storm sub-IDs to be generated
        
        ''' Populate dictionary with data. '''
        # Methodology:  (1) get matching timestamps between GCM output and TC tracker output. 
        #               (2) then, get the timestamp corresponding to LMI from the TC tracker output.
        #               (3) finally, get timestamps from +/- N days from the LMI timestamp, where N is defined in the script as 'n_days'
        n_days = 0
        # Get Pandas-friendly timestamps for DataFrame indexing
        timestamps = [utilities.time_adjust(model, t, method='cftime_to_pandas') 
                        for t in dataset.time.values if not ((t.month == 2) & (t.day == 29))] 
        # Get strongest timestamp
        timestamp_lmi = ds['track_output'].loc[ds['track_output']['time'].isin(timestamps)].sort_values('min_slp', ascending=True).iloc[0]['time'] 
        # Get the timestamps found from +/- n_days and revert to cftime for DataArray indexing
        timestamps = [utilities.time_adjust(model, t, method='pandas_to_cftime') for t in timestamps 
                        if abs(t - timestamp_lmi) <= pd.Timedelta(n_days, "d")]
        
        # Process data for each filtered timestamp
        for i, timestamp in enumerate(timestamps):   
            # Assign the storm sub-ID
            storm_subid = '{0}_{1:03d}'.format(storm_id, i)
            # Extract all nans
            data[storm_subid] = {key: dataset.sel(time=timestamp).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')}
            # Get time and intensity for the selected storm at the given timestamp
            if track_data and weighting_intensity_metric:
                timestamp_TC_tracker = utilities.time_adjust(model, timestamp, 'cftime_to_pandas') # convert GCM timestamp in cftime format to Pandas timestamp format
                timestamp_intensity = ds['track_output'].loc[ds['track_output']['time'] == timestamp_TC_tracker][weighting_intensity_metric] # get the TC tracker output corresponding to this timestamp
                if diagnostic:
                    print('[composite.py, planar_compositor()] weighting timestamp, GCM: {0}; weighting timestamp, tracker: {1}; weighting intensity: {2}'.format(timestamp, timestamp_TC_tracker, timestamp_intensity))
                weighting_metadata[storm_id] = timestamp_intensity.item() # append to the dictionary
            # Append to storage list
            storm_subids.append(storm_subid)
            if diagnostic:
                print('\t \t \t {0} created, continuing...'.format(storm_subid))
                
        # Iterate through all matching entries (also known as storm sub-IDs) for the iterand storm ID
        for storm_subid in storm_subids:
            # Get aspect ratio of clipped GCM data - used to see if an entry has invalid dimensions
            storm_aspect_ratio = len(data[storm_subid][key].grid_xt) / len(data[storm_subid][key].grid_yt)
            # Check to see if data proportions are too far from unity (typical values range from 0.75 to 1.35). If so, drop it.
            if (storm_aspect_ratio < (1/1.5) or storm_aspect_ratio > 1.5) and compositing_mode == 'planar':
                if diagnostic:
                    print('[composite.py, planar_compositor()] dimensions too skewed for this entry: {0} with shape {1}'.format(storm_subid, data[storm_subid][key][field].shape))
                del data[storm_subid]
                continue
            # Log the extents in longitude (grid_xt) and latitude (grid_yt)
            min_x[storm_subid] = len(data[storm_subid][key].grid_xt)
            min_y[storm_subid] = len(data[storm_subid][key].grid_yt)
            # Add storm ID to the dataset
            data[storm_subid][key].attrs['storm_ID'] = storm_subid
    else:
        print('\t \t {0} was not processed since no timestamps were found...'.format(storm_id))
     
    output_dict = {'data': data,
                   'weighting_metadata': weighting_metadata,
                   'min_x': min_x,
                   'min_y': min_y} 
        
    return output_dict

def snapshot_weight_trim(model, key, data, storm_id, field, pressure_level, intensity_bin, min_x, min_y,
                         weights=None, weighting_metadata=None, weighting_bins=None, rotation=False, 
                         buffer=1, NH_TCs=False, diagnostic=False):
    
    storm_weights = {}
    
    if diagnostic:
        print('[composite.py, snapshot_weight_trim()]: storm ID: {0}'.format(storm_id))

    # Make minimum domain lengths even by reducing frame size, if not even
    min_x = min_x - 1 if min_x % 2 == 1 else min_x
    min_y = min_y - 1 if min_y % 2 == 1 else min_y
    # Get midpoints
    center_x, center_y = int(len(data[storm_id][key].grid_xt) // 2), int(len(data[storm_id][key].grid_yt) // 2)
    # Define slices
    slice_x = slice(center_x - min_x//2 + buffer, center_x + min_x//2 - buffer)
    slice_y = slice(center_y - min_y//2 + buffer, center_y + min_y//2 - buffer)
    # Slice the domain
    data[storm_id][key] = data[storm_id][key].isel(grid_xt=slice_x, grid_yt=slice_y)
    
    # Redefine midpoints for meridional position filtering
    center_y = int(len(data[storm_id][key].grid_yt) // 2)
    # Latitude check - only use storms within 5 to 35 latitude
    if abs(data[storm_id][key]['center_lat']) > 35 or abs(data[storm_id][key]['center_lat']) < 5:
        print('[composite.py, snapshot_weight_trim()] Not using {0} due to TC center being too poleward at {1:.2f}.'.format(storm_id, data[storm_id][key]['center_lat'].values))
        data[storm_id][key] = np.nan
        storm_weights[storm_id] = 0
        # continue
    # Flip the data about the x-axis if the center latitude is < 0 (for Southern Hemisphere storms).
    elif data[storm_id][key].isel(grid_yt=center_y)['grid_yt'] < 0:
        if NH_TCs:
            print('[composite.py, snapshot_weight_trim()] Not using {0} due to TC being in the Southern Hemisphere at {1:.2f}.'.format(storm_id, data[storm_id][key].isel(grid_yt=center_y)['grid_yt']))
            data[storm_id][key] = np.nan
            storm_weights[storm_id] = 0
        else:
            data[storm_id][key][field] = np.flip(data[storm_id][key][field], axis=0)
    
    # Rotate the domain and load the data
    if rotation:
        print('[composite.py, planar_compositor()] Rotation enabled...')
        rotated = rotate.rotate_snapshot(data[storm_id][key], field, diagnostic=False)
        data[storm_id][key] = rotated['{0}_rotated'.format(field)]
    
    # Get intensity bin for the specific storm ID (find the bin that the storm ID intensity metric falls into)
    weighting_storm_id = storm_id.split('_')[0] if weighting_bins else 1
    # if diagnostic:
    #     print('[composite.py, planar_compositor()] Weighting bins: {0}; storm-specific intensity: {1}'.format(weighting_bins, weighting_metadata[weighting_storm_id]))
    weighting_index = np.nanmax(np.flatnonzero(np.array(weighting_bins) < weighting_metadata[weighting_storm_id])) if weighting_bins else 1
    if diagnostic:
        print('[composite.py, planar_compositor()] Weights: {0}; Weighting index: {1}'.format(weights, weighting_index))
    # Use intensity-based weighting for each selected timestamp
    if weights:
        storm_weights[storm_id] = list(weights.values())[weighting_index] if weighting_bins else 1
    else:
        storm_weights[storm_id] = 1
    
    # nan-out the storm weight if the data entry is invalid
    if 'float' in str(type(data[storm_id][key])):
        storm_weights[storm_id] = np.nan
    
    output_dict = {'data': {storm_id: data[storm_id][key]}, 'storm_weight': {storm_id: storm_weights[storm_id]}}
    
    return output_dict

def dataset_alignment(ds, subdict, field, times, pressure_level=None, index_timestamps=None, diagnostic=False):
    
    index_timestamps = [t for t in times if t in index_timestamps]
    if diagnostic:
        [print('[composite.py] Timestamps level 3:', index_timestamp) for index_timestamp in index_timestamps]
    
    if len(index_timestamps) == 0:
        print('No matching timestamps, continuing to the next TC...')
        # continue
    
    ''' Rotation support steps (aligns planar and vertical data to unite into single data structure). '''
    # Get planar data
    # and field + vertical level from the vertical data
    # This try/except is supposed to catch duplicate time indices in ds[subdict]
    dataset_planar = {}
    try:
        dataset_planar[field] = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[subdict][field].sel(time=index_timestamps)
    except:
        # Check for duplicate time values
        ds[subdict] = ds[subdict].drop_duplicates('time')
        # Now try reloading
        dataset_planar[field] = ds[subdict][field].sel(time=index_timestamps).sel(pfull=pressure_level, method='nearest') if pressure_level else ds[subdict][field].sel(time=index_timestamps)
    
        
    # Attempt to add center coordinates to support rotation function
    for planar_rotation_field in ['center_lon', 'center_lat', 'slp']:
        dataset_planar[planar_rotation_field] = ds['tc_model_output'].sel(time=index_timestamps)[planar_rotation_field]
    dataset_planar = xr.merge(dataset_planar.values())
    if diagnostic:
        print('Planar dataset merged.')
    
    # Get vertical data
    dataset_vertical = {}
    for vertical_rotation_field in ['ucomp', 'vcomp', 'sphum']:
        if vertical_rotation_field != field:
            dataset_vertical[vertical_rotation_field] = ds['tc_vertical_output'].sel(time=index_timestamps)[vertical_rotation_field]
    dataset_vertical = xr.merge(dataset_vertical.values())
    if diagnostic:
        print('Vertical dataset merged.')
    
    # Merge the two datasets
    dataset = xr.merge([dataset_planar, dataset_vertical])
    if diagnostic:
        print('Vertical and planar datasets merged.')
    for data_var in dataset_planar.data_vars:
        if 'grid_xt' in dataset_planar[data_var].dims and 'grid_yt' in dataset_planar[data_var].dims:
            dataset[data_var] = dataset[data_var].interpolate_na(dim='grid_xt').interpolate_na(dim='grid_yt')
    
    ''' End planar-specific block of code. '''

    return dataset

def planar_composite_mean(data, snapshot_stack, key):
    
    snapshot_stack_key = [s[key] for s in snapshot_stack]
    
    # This is the compositing ~magic~, where the arrays get flattened into a 2D array of axes 'grid_xt' and 'grid_yt'
    composite_mean = np.nanmean(np.stack(snapshot_stack_key), axis=0)
        
    # composite_mean = np.nanmean(np.stack([v[key] * storm_weights[k] for k, v in data.items()]), axis=0)
    # Get degree spacing for the data. Assumes all spacings are equal for each dimension
    dx, dy = [[entry[key] for entry in data.values()][0].grid_xt.diff(dim='grid_xt').values[0],
              [entry[key] for entry in data.values()][0].grid_yt.diff(dim='grid_yt').values[0]]
    # Define dictionary to hold distance data (x and y) from composite_mean shape
    # The domains are define as +/- 1/2 the width (centered on 0), stepping by 1 and scaling by dx and dy
    composite_mean_x = dx*np.arange(-composite_mean.shape[1]//2, composite_mean.shape[1]//2, 1)
    composite_mean_y = dy*np.arange(-composite_mean.shape[0]//2, composite_mean.shape[0]//2, 1)
    # Populate with the dimensional data
    composite_mean = {'grid_xt': composite_mean_x, 'grid_yt': composite_mean_y, 'data': composite_mean}    
    # Note: composite_mean can be made into an xArray DataArray directly from a dictionary.
    composite_mean = xr.DataArray(data=composite_mean['data'], dims=['grid_yt', 'grid_xt'], 
                                  coords={'grid_yt': (['grid_yt'], composite_mean['grid_yt']), 
                                          'grid_xt': (['grid_xt'], composite_mean['grid_xt'])})
    
    return composite_mean

def azimuthal_composite_mean(snapshots, outer_radius_index=20, diagnostic=False):
    
    print('[composite.py, azimuthal_composite_mean()] Processing azimuthal mean...')
    # Ensure all entries have like dimensions
    snapshot_dimensions = [snapshot.dims for snapshot in snapshots]
    if snapshot_dimensions.count(snapshot_dimensions[0]) != len(snapshot_dimensions):
        print('[composite.py, azimuthal_composite_mean()] Not all snapshots have the same dimensions. Check that out...')
    # For each individual storm, ensure that the DataArray populates some fraction of the initialized nan array
    
    # Handle 2D data (vertical and radial)
    if 'pfull' in snapshot_dimensions[0] and 'radius' in snapshot_dimensions[0]:
        # Initialize a container for all snapshots. This is made to grab a mean after iteration.
        container = np.full(shape=(32, outer_radius_index, len(snapshots)), fill_value=np.nan)
        for i, value in enumerate(snapshots):
            # Ensure dimensions are properly oriented
            v = value.transpose('pfull', 'radius')
            # Trim data radially, if needed
            if v.values.shape[1] > outer_radius_index:
                arr = v.values[:, :outer_radius_index]
            else:
                arr = v.values
            container[:, 0:arr.shape[1], i] = arr
            
        composite_mean = xr.DataArray(dims=['pfull', 'radius'],
                                      coords={'pfull': (['pfull'], v.pfull.values),
                                              'radius': (['radius'], np.arange(0, outer_radius_index*0.9375, 0.9375))},
                                      data=np.nanmean(container, axis=2))
        return composite_mean
    
    # Handle 1D data (radial)
    elif 'radius' in snapshot_dimensions[0]:
        # Initialize a container for all snapshots. This is made to grab a mean after iteration.
        container = np.full(shape=(outer_radius_index, len(snapshots)), fill_value=np.nan)
        for i, value in enumerate(snapshots):
            if value.values.shape[0] > outer_radius_index:
                arr = value.values[:outer_radius_index]
            else:
                arr = value.values
            container[0:arr.shape[0], i] = arr
            
        composite_mean = xr.DataArray(dims=['radius'], 
                                      coords={'radius': (['radius'], np.arange(0, outer_radius_index*0.9375, 0.9375))},
                                      data=np.nanmean(container, axis=1))
        
        return composite_mean
            
    else:
        print('[composite.py, azimuthal_composite_mean()] Dimensions not recognized and cannot be composited.')

