''' Import packages. '''
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy, numba
# Local data storage packages
import dill, os, pickle
# Data structure packages
import pandas as pd, xarray as xr
# Visualization tools
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt
# Local imports
import accessor, composite, composite_weighting, derived, utilities, socket, visualization, tc_analysis, tc_processing, rotate

import importlib
importlib.reload(composite)
importlib.reload(composite_weighting)
importlib.reload(utilities)
importlib.reload(tc_analysis)
importlib.reload(tc_processing)
importlib.reload(visualization)
importlib.reload(rotate)
importlib.reload(derived)


def data_access(models, experiments, ):
    models = ['HIRAM', 'AM2.5', 'FLOR']
    experiments = ['CTL1990s', 'CTL1990s_swishe']
    storm_type = 'TS'
    data = tc_analysis.tc_model_data(models, experiments, storm_type, num_storms=100, single_storm_id=None)
    track_data = tc_analysis.track_data_loading(models=models, experiments=experiments, storm_type=storm_type, year_range=(101, 150))

    return data, track_data

def snapshot_generator(data, track_data, models, experiments, intensity_bin='b1', fields=None, diagnostic=False):
    if not fields:
        fields = {'planar': {'precip': None, 'WVP': None, 'olr': None},
                  'azimuthal': {'lhflx': None, 'shflx': None, 'wind': None, 'omega': None, 'precip': None, 'WVP': None, 'olr': None, 'wind_tangential': None, 'wind_radial': None, 'h_anom': None}}

    for compositing_mode in ['planar', 'azimuthal']:
        for field, pressure_level in fields[compositing_mode].items():
            for model in ['HIRAM']:
                for experiment in experiments:
                    if diagnostic:
                        print('[snapshot_generator()] Processing {0} snapshots for {1} at {2} hPa for TCs in the {3} {4} configuration...'.format(compositing_mode, field, pressure_level, model, experiment))
                    
                    snapshots_experiment, _, _ = composite.compositor_preprocessor(model, data[model][experiment]['data'], intensity_bin, field, NH_TCs=True, 
                                                                                   compositing_mode=compositing_mode, track_data=track_data, weighting_intensity_metric='min_slp')
                    
                    snapshot_field_save(snapshots_experiment, compositing_mode=compositing_mode, model=model, experiment=experiment, storm_type='TS', field=field, intensity_bin=intensity_bin)

def snapshot_field_save(snapshots, compositing_mode, model, experiment, storm_type, field, intensity_bin, diagnostic=False):
    ''' Method to save a snapshot for a given configuration. This is supposed to accelerate compositing and bootstrapping. '''
    # Get times from each snapshot, assumes each snapshot has a single time entry
    timestamps = [item.time.item() for item in snapshots]
    # Get each unique snapshot and the number of times each unique snapshot is repeated
    snapshot_count = {timestamp: timestamps.count(timestamp) for timestamp in timestamps}
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/snapshots'
    # Iterate through each snapshot
    for snapshot in snapshots:
        # Pull the storm ID from the attributes
        storm_ID = snapshot.attrs['storm_ID']
        if diagnostic:
            print('[composite_snapshots.py, snapshot_field_save()] Iterand storm ID: {0}'.format(storm_ID))
        # Capture how many times the snapshot is repeated for the filename
        snapshot_repetitions = snapshot_count[snapshot.time.item()]
        # Construct a filename
        filename = 'TC_snapshot-{0}-{1}-{2}-{3}-{4}-{5}-{6}.nc'.format(compositing_mode, model, experiment, storm_type, storm_ID, field, intensity_bin)
        snapshot_path = os.path.join(dirname, filename)
        # List all files in directory
        filenames = os.listdir(dirname)
        # Save, if it does not already exist
        if not os.path.isfile(snapshot_path) and filename not in filenames:
            print('[composite_snapshots.py, snapshot_field_save()] Saving snapshot to {0}...'.format(snapshot_path))
            # Note: this is done for the azimuthal field because of theD
            if compositing_mode == 'azimuthal':
                snapshot.to_dataset(name=field, promote_attrs=True).to_netcdf(path=snapshot_path, mode='w', format='netcdf4')
            else:
                snapshot.to_netcdf(path=snapshot_path, mode='w', format='netcdf4')

def snapshot_load(compositing_mode, model, experiment, storm_type, field, intensity_bin, diagnostic=False):
    ''' Method to load a snapshot for a given configuration. This is supposed to accelerate compositing and bootstrapping. '''
    
    # List to collect snapshots that match specific criteria
    snapshot_stack = []
    # Define the directory name where snapshots are stored
    dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/snapshots'

    filenames = [filename for filename in os.listdir(dirname) if filename.endswith('.nc')]
    
    if diagnostic:
        print('[composite_snapshots, snapshot_load()] Looking for snapshots for model {0} in configuration {1} with field {2} for a {3} composite...'.format(model, experiment, field, compositing_mode))
    for filename in filenames:
        if diagnostic:
            print('[composite_snapshots, snapshot_load()] Processing filename {0}'.format(filename))
        # Decompose the filename
        filename_split = filename.split('.nc')[0].split('-')
        # Get file identifiers
        compositing_mode_substr = filename_split[1]
        model_substr = filename_split[2]
        experiment_substr = filename_split[3]
        storm_type_substr = filename_split[4]
        field_substr = filename_split[7]
        intensity_bin_substr = filename_split[8]
        # print(compositing_mode_substr, model_substr, experiment_substr, storm_type_substr, field_substr, intensity_bin_substr)
        # Get a boolean criterion to parse for matching files
        matching_criteria = (compositing_mode_substr == compositing_mode) & (model_substr == model) & (experiment_substr == experiment) & \
                            (storm_type_substr == storm_type) & (field_substr == field) & (intensity_bin_substr == intensity_bin)
        # diagnostic = True
        if diagnostic:
            print('[composite_snapshots, snapshot_load()] Matching criteria:\n', 
                  compositing_mode_substr == compositing_mode, 
                  model_substr == model,
                  experiment_substr == experiment,
                  storm_type_substr == storm_type,
                  field_substr == field,
                  intensity_bin_substr == intensity_bin)
        # Collect files that match the criteria
        if matching_criteria:
            if diagnostic:
                print('[snapshot_load()] Filename {0} matches the criteria!'.format(filename))

            # Load the data
            temp = xr.open_dataset(os.path.join(dirname, filename))

            # Crude modification due to failure of correcting it during composite generation
            # This corrects the values for precipitation (kg/m2/s to mm/d)
            if compositing_mode == 'planar' and field == 'precip':
                temp[field] = utilities.field_correction(temp[field], field)
            
            # Append to the output list of snapshots
            snapshot_stack.append(temp)
            del temp

    if diagnostic:
            print('[composite_snapshots.py, snapshot_load()] Number of loaded snapshots: {0}'.format(len(snapshot_stack)))
    
    return snapshot_stack

def snapshot_constructor(compositing_method, model, experiments, storm_type, field, intensity_bin, diagnostic=True):
    snapshots = {}
    for experiment in experiments:
        if diagnostic:
            print('[composite_snapshots, snapshot_constructor()] Loading snapshots for {0} in the {1} configuration.'.format(model, experiment))
        snapshots[experiment] = snapshot_load(compositing_method, model, experiment, storm_type, field, intensity_bin, diagnostic=diagnostic)
        if diagnostic:
            print('[composite_snapshots.py, snapshot_constructor()] Number of loaded snapshots: {0}'.format(len(snapshots[experiment])))
    
    return snapshots

def get_snapshot_weights(compositing_method, snapshots, experiments, intensity_bins=None, diagnostic=False):

    experiment_control, experiment_run = experiments

    intensity_bins = np.arange(900, 1020, 4)
    intensity_distributions, intensity_distributions_binned = {}, {}
    for experiment in experiments:
        if compositing_method == 'planar':
            [print(experiment, snapshot['slp'].min().item()) for snapshot in snapshots[experiment]]
            intensity_distributions[experiment] = np.array([snapshot['slp'].min().item() for snapshot in snapshots[experiment]])
        else:
            intensity_distributions[experiment] = np.array([snapshot.attrs['slp'] for snapshot in snapshots[experiment]])
        intensity_distributions_binned[experiment], _ = np.histogram(intensity_distributions[experiment], density=True, bins=intensity_bins)

    weights = {}
    for experiment_index, experiment in enumerate(experiments):
        if experiment == experiment_control:
            raw_weights = intensity_distributions_binned[experiment_run] / intensity_distributions_binned[experiment_control]
            # print(raw_weights)
            # print('-------------')
            filtered_weights = np.where(np.isfinite(raw_weights), raw_weights, 0)
            # print(filtered_weights)
            weights[experiment] = {intensity_bins[i]: filtered_weights[i] for i in range(len(filtered_weights))}
        else:
            filtered_weights = np.repeat(1, len(intensity_distributions_binned[experiment]))
            weights[experiment] = {intensity_bins[i]: filtered_weights[i] for i in range(len(filtered_weights))}
        
    diagnostic = True
    if diagnostic:
        fig, ax = plt.subplots(figsize=(4, 2))
        linestyles = ['-', '--']
        for experiment_index, experiment in enumerate(experiments):
            ax.step(intensity_bins[:-1], intensity_distributions_binned[experiment]*np.fromiter(weights[experiment].values(), 'float'), 
                    label=experiment, linestyle=linestyles[experiment_index])
        ax.legend(frameon=False)

    return intensity_bins, weights

def weighted_snapshot_constructor(compositing_method, snapshots, experiments, weights, intensity_bins, diagnostic=False):
    
    # Initialize size for minimum dimension
    min_x, min_y = np.nan, np.nan
    # Define dimensions for x and y
    dimension_x = 'radius' if 'azimuthal' in compositing_method else 'grid_xt'
    dimension_y = 'pfull' if 'azimuthal' in compositing_method else 'grid_yt'

    # Define interval to which the snapshot weight is rounded
    repetition_interval = 0.25
    # Weighted snapshot container
    weighted_snapshots = {}
    # Initialize container to hold total number of snapshots after repetitions
    repetition_sample_size = {}
    # Iterate over each experiment
    for experiment in experiments:
        weighted_snapshots[experiment] = {}
        repetition_sample_size[experiment] = 0
        if diagnostic:
            print('[composite_snapshots.py, weighted_snapshot_constructor()] Number of snapshots for experiment {0}: {1}'.format(experiment, len(snapshots[experiment])))
        # Iterate over each snapshot in the experiment
        for index, snapshot in enumerate(snapshots[experiment]):
            # Get storm ID for future dictionary assignment
            storm_ID = snapshot.attrs['storm_ID']
            # Get the minimum intensity from the snapshot
            intensity = snapshot['slp'].min().item() if compositing_method == 'planar' else snapshot.attrs['slp']
            # Find the closest matching intensity bin in the weights list
            closest_intensity = intensity_bins[(np.abs(intensity_bins - intensity)).argmin()]
            # Pull the weight matching the closest intensity bin
            snapshot_weight = weights[experiment][closest_intensity]
            # Get an integer value of repetitions for the given weight
            snapshot_repetitions = round(repetition_interval * np.round(snapshot_weight/repetition_interval))
            # Define new key:value pair where the key is the number of repetitions, and the value is the data
            weighted_snapshots[experiment][storm_ID] = {'repetitions': snapshot_repetitions, 'data': snapshot}
            # Add to the experiment-specific counter
            repetition_sample_size[experiment] += snapshot_repetitions
            if diagnostic:
                print('[] Snapshot intensity: {0:.2f}; closest bin intensity: {1:.2f}; corresponding weight: {2:.2f}; number of reps: {3}'.format(intensity, closest_intensity, snapshot_weight, snapshot_repetitions))
        if diagnostic:
            print('[composite_snapshots.py, weighted_snapshot_constructor()] Number of weighted snapshots for experiment {0}: {1}'.format(experiment, len(weighted_snapshots[experiment])))
            print('[] Number of unique snapshots: {0}; number of repetitions: {1}'.format(len(snapshots[experiment]), repetition_sample_size[experiment]))
        for storm_ID, _ in weighted_snapshots[experiment].items():
            # Trim the pressure levels for azimuthal dataset
            if compositing_method == 'azimuthal':
                weighted_snapshots[experiment][storm_ID]['data'] = weighted_snapshots[experiment][storm_ID]['data'].sel(pfull=slice(100, 1000))
            # New repetition value is equal to (number of snapshot repetitions / number of total experiment repetitions)
            weighted_snapshots[experiment][storm_ID]['weight'] = weighted_snapshots[experiment][storm_ID]['repetitions'] / repetition_sample_size[experiment]
            # Check to see if dimension length is shorter than the existing threshold
            length_x, length_y = len(weighted_snapshots[experiment][storm_ID]['data'][dimension_x]), len(weighted_snapshots[experiment][storm_ID]['data'][dimension_y])
            min_x = length_x if (np.isnan(min_x) or length_x < min_x) else min_x
            min_y = length_y if (np.isnan(min_y) or length_y < min_y) else min_y
            
    # Make another pass so that all snapshots are trimmed to the group minima
    for experiment in experiments:
        for storm_ID in weighted_snapshots[experiment].keys():
            weighted_snapshots[experiment][storm_ID]['data'] = weighted_snapshots[experiment][storm_ID]['data'].isel({dimension_x: slice(0, min_x)})
        if diagnostic:
            print('Snapshot stack length for experiment {0}: {1}'.format(experiment, len(weighted_snapshots[experiment])))

    return weighted_snapshots

def stack_constructor(weighted_snapshots, experiments, field, diagnostic=False):
    ''' Construct stacks for snapshots, bootstrapping snapshots (this is just repeated snapshots), and weights.'''
    # Get xarray-compatible field name for pressure-level slices
    xarray_field = field.split('hPa')[0] if 'hPa' in field else field
    # Initialize stacks
    snapshot_stack, bootstrap_stack, snapshot_values = {}, {}, {}
    # Iterate over each given experiment
    for experiment in experiments:
        # Initialize experiment-specific stack list
        snapshot_value_stack, snapshot_bootstrap_stack, snapshot_weight_stack = [], [], []
        print('[composite_snapshots.py, stack_constructor()] Snapshot stack length for {0}: {1}'.format(experiment, len(weighted_snapshots[experiment])))
        # Iterate through each storm in the constructed weighted snapshot stack
        for storm_ID in weighted_snapshots[experiment].keys():
            x_i = weighted_snapshots[experiment][storm_ID]['data'][xarray_field]
            w_i = weighted_snapshots[experiment][storm_ID]['weight']
            r_i = weighted_snapshots[experiment][storm_ID]['repetitions']
            snapshot_value_stack.append(w_i * x_i)
            snapshot_weight_stack.append(w_i)
            
            # Ensure minimum number of repetitions is 1
            r_i = 1 if r_i == 0 else r_i
            for i in range(r_i):
                snapshot_bootstrap_stack.append(x_i)
        
        snapshot_stack[experiment] = np.stack(snapshot_value_stack)
        bootstrap_stack[experiment] = np.stack(snapshot_bootstrap_stack)
        snapshot_values[experiment] = np.sum(np.stack(snapshot_value_stack), axis=0)
        if diagnostic:
            print('Bootstrap stack length: ', len(snapshot_bootstrap_stack))
        print('[composite_snapshots.py, stack_constructor()] Bootstrap stack length for {0}: {1}'.format(experiment, len(bootstrap_stack[experiment])))
        

    return snapshot_stack, bootstrap_stack, snapshot_values

def composite_bootstrap(compositing_method, bootstrap_stack, model, experiments, field, bootstrap_estimates=10000, confidence_interval=0.9, diagnostic=False):
    
    experiment_control, experiment_run = experiments
    a, b = bootstrap_stack[experiment_control], bootstrap_stack[experiment_run] 

    # Check for input array size mismatch along the x- and y-axes
    if a.shape[1:] != b.shape[1:]:
        print('[composite_snapshots.py, composite_bootstrap()] Input array shape mismatch for model {0} and field {1}'.format(model, field))
        print('[composite_snapshots.py, composite_bootstrap()] Control shape: {0}, experiment shape: {1}'.format(a.shape, b.shape))

    start_time = time.time()
    _, _, median, differences, lower_tail, upper_tail = utilities.bootstrap_3d(b, a, plane=None, N=bootstrap_estimates, level=confidence_interval)
    end_time = time.time()
    if diagnostic:
        print('[composite_snapshots.py, composite_bootstrap()] Elapsed bootstrap time: {0:.3f} s; elapsed time per bootstrap sample: {1:.3e}'.format((end_time - start_time),
                                                                                                                                                 (end_time - start_time)/bootstrap_estimates))

    # Get statistics
    bootstrap_sample_number = min(len(bootstrap_stack[experiment_control]), len(bootstrap_stack[experiment_run]))
    
    print('[composite_snapshots.py, composite_bootstrap()] Number of unique samples (minimum from both input snapshot stacks): {0}'.format(bootstrap_sample_number))

    significance_masking = np.where(np.sign(lower_tail) == np.sign(upper_tail), True, False)
    significant_difference = np.where(significance_masking, median.T, np.nan)
    outputs = [lower_tail, significant_difference, upper_tail]
    if diagnostic:
        print('[composite_snapshots.py, composite_bootstrap()] Domain-integrated difference: {0:.2f}'.format(np.nansum(significant_difference)))

    if diagnostic:
        fig, axes = plt.subplots(figsize=(8, 3), ncols=3)
        pfull = np.array([110.419627, 141.09261 , 177.729388, 220.892397, 271.066624, 328.516337, 
                          392.785273, 461.947262, 532.465907, 600.430867, 663.107383, 719.307118, 
                          768.814284, 811.846869, 848.836021, 880.346139, 906.995722, 929.394583, 
                          948.128523, 963.73257 , 976.687397, 987.392458, 996.109949])

        for index, ax in enumerate(fig.axes):
            vmin, vmax = np.nanmin(outputs[index]), np.nanmax(outputs[index])
            vmin = -1 if np.isnan(vmin) else vmin
            vmax = 1 if np.isnan(vmax) else vmax
            norm, cmap = visualization.norm_cmap(0, field=field, extrema=(vmin, vmax))
            if compositing_method == 'azimuthal':
                im = ax.contourf(range(outputs[index].shape[1]), pfull, outputs[index], norm=norm, cmap=cmap, levels=len(norm.boundaries))
                ax.set_ylim(ax.get_ylim()[::-1])
            else:
                im = ax.pcolormesh(outputs[index], norm=norm, cmap=cmap)

            cax = ax.inset_axes([0, 1.1, 1, 0.05])
            colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
            cax.xaxis.set_ticks_position('top')
            cax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        
        fig.tight_layout()
        
        title = 'Model: {0}; field: {1}'.format(model, field)
        fig.suptitle(title, y=1.1)

        fig, ax = plt.subplots(figsize=(3, 3))
        norm, cmap = visualization.norm_cmap(0, field=field, extrema=(np.nanmin(median), np.nanmax(median)))
        if compositing_method == 'azimuthal':
            X, Y = np.meshgrid(range(median.shape[0]), pfull)
            ax.contourf(X, Y, median.T, norm=norm, cmap=cmap, levels=len(norm.boundaries))
            ax.set_ylim(ax.get_ylim()[::-1])
        else:
            ax.pcolormesh(median, norm=norm, cmap=cmap)

        cax = ax.inset_axes([0, 1.1, 1, 0.05])
        fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
        cax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        cax.xaxis.set_ticks_position('top')

    return lower_tail, median, upper_tail, significance_masking

def main(mode, model, experiments, field, storm_type, intensity_bin, compositing_method=None, number_of_samples=None, diagnostic=False):

    if mode == 'snapshot_generation':
        model_output, track_output = data_access(model, experiments, storm_type)
        snapshot_generator(model_output, track_output, model, experiments, intensity_bin, diagnostic=True)
    
    else:
        # Build dictionary housing snapshot stacks for both experiments considered, and generate intensity distributions to inform sample weighting
        snapshots = snapshot_constructor(compositing_method, model, experiments, storm_type, field, intensity_bin, diagnostic=diagnostic)
        
        # If a sample number is given, randomly sample the snapshot stack and choose the number of sapmles
        if number_of_samples and number_of_samples < len(snapshots):
            sample_indices = np.random.randint(low=0, high=len(snapshots), shape=(number_of_samples))
            snapshots = snapshots[sample_indices]
        else:
            snapshots = snapshots

        # Get snapshot intensity bins and weights from the preloaded distributions
        intensity_bins, weights = get_snapshot_weights(compositing_method, snapshots, experiments, intensity_bins=None, diagnostic=diagnostic)

        # Construct a dictionary with weighted snapshots using the calculated intensity bins and weights
        weighted_snapshots = weighted_snapshot_constructor(compositing_method, snapshots, experiments, weights, intensity_bins, diagnostic=diagnostic)
        
        # Obtain snapshot stacks for bootstrapping
        _, bootstrapping_stack, _ = stack_constructor(weighted_snapshots, experiments, field, diagnostic=diagnostic)

        # Perform bootstrapping and return lower limit, mean, and upper limit of the distribution of mean differences between experiments. 
        # Note that the mean is masked for significance to the preset confidence interval.
        lower_limit, median, upper_limit, significance_masking = composite_bootstrap(compositing_method, bootstrapping_stack, model, experiments, field, 
                                                                   bootstrap_estimates=10000, confidence_interval=0.95, diagnostic=False)

        return bootstrapping_stack, lower_limit, median, upper_limit, significance_masking

if __name__ == '__main__':
    models = ['AM2.5']
    fields = ['WVP']
    compositing_method = 'planar'
    for field in fields:
        for model in models:
            stacks, lower_tail, median, upper_tail, significance_masking = main(mode='bootstrapping', model=model, experiments=['CTL1990s', 'CTL1990s_swishe'], field=field,
                                                                                storm_type='TS', intensity_bin='b1', compositing_method=compositing_method, diagnostic=False)