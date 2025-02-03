''' Import packages. '''
# Time packages
import cftime, datetime, time
# Numerical analysis packages
import numpy as np, random, scipy
# Local data storage packages
import dill, os, pickle
# Data structure packages

import pandas as pd, xarray as xr
# Visualization tools
import cartopy, cartopy.crs as ccrs, matplotlib, matplotlib.pyplot as plt
# Local imports
import accessor, composite, derived, utilities, socket, visualization, tc_analysis, tc_processing, rotate

import importlib
importlib.reload(composite)
importlib.reload(utilities)
importlib.reload(tc_analysis)
importlib.reload(tc_processing)
importlib.reload(visualization)
importlib.reload(rotate)
importlib.reload(derived)

def distribution_weighting(data, intensity_metric='max_wind', input_type='snapshots', diagnostic=True):

    '''
    This method takes two distributions of data and normalizes the first distribution to the second.
    This assumes that the first distribution is a subset of the second.
    The distribution is based on a given intensity metric for TCs, typically 'max_wind' or 'min_slp'.
    '''

    if len(data) != 2:
        print('[distribution_weighting()] Incorrect number of input experiments.')

    # Get names of the experiments
    CTL, EXP = list(data.keys())[0], list(data.keys())[1]

    # Define intensity bins common to both distributions
    bin_size = 2 if intensity_metric == 'max_wind' else 4
    bins = np.arange(10, 50, bin_size) if intensity_metric == 'max_wind' else np.arange(900, 1020, bin_size)
    
    # Get probability densities for each experiment dataset
    distributions = {}
    for experiment in [CTL, EXP]:

        distributions[experiment] = np.histogram(data[experiment][intensity_metric], density=True, bins=bins)
    
    # Get weights for the control distribution 
    # For each intensity bin, get the ratio of the experiment probability density to the control probability density
    weights = {}

    weighting_mode = 'density_ratio'
    if weighting_mode == 'density_ratio':
        for bin_index, bin in enumerate(bins[:-1]):
            density_CTL = distributions[CTL][0][bin_index]
            density_EXP = distributions[EXP][0][bin_index]
            density_ratio = density_EXP/density_CTL
            # Get density ratios
            if np.isinf(density_ratio):
                # If the ratio is infinite, that means no control entries exist in this intensity bin but experiment entries do
                if diagnostic:
                    print('Intensity bin: {0}; EXP to CTL density ratio: {1:.2f}'.format(bin, 0))
                weights[bin] = 0
            elif np.isnan(density_ratio):
                # If the ratio is nan, that means no entries exist in this intensity bin for either distribution
                if diagnostic:
                    print('Intensity bin: {0}; EXP to CTL density ratio: {1:.2f}'.format(bin, 0))
                weights[bin] = np.nan
            else:
                if diagnostic:
                    print('Intensity bin: {0}; EXP to CTL density ratio: {1:.2f}'.format(bin, density_ratio))
                weights[bin] = density_ratio
                
        # maximum_ratio = np.nanmax(list(weights.values()))
        # for bin_index, bin in enumerate(bins[:-1]):
        #     weights[bin] = weights[bin]/maximum_ratio if weights[bin] > 1 else weights[bin]       
    else:
        normal_distribution = np.random.normal(data[EXP][intensity_metric].mean(), data[EXP][intensity_metric].std(), size=10000)
        normal_weights, _ = np.histogram(normal_distribution, density=True, bins=bins)
        median_ratio = np.mean(data[CTL][intensity_metric])/np.mean(data[EXP][intensity_metric])
        print('[distribution_weighting()] ratio of CTL to EXP median: {0:.2f}'.format(median_ratio))
        for bin_index, bin in enumerate(bins[:-1]):
            weights[bin] = normal_weights[bin_index]/np.max(normal_weights)

    weighted_distributions = {}
    for index, experiment in enumerate([CTL, EXP]):
        arr, _ = np.histogram(data[experiment][intensity_metric], density=True, bins=bins)
        weighted_distributions[experiment] = arr * np.array(list(weights.values())) if experiment == CTL else arr

    # Plot the distribution densities before and after normalization
    if diagnostic:
        colors = ['tab:blue', 'tab:orange']
        linestyles = ['-', '--']
        fig, axes = plt.subplots(figsize=(8, 2), ncols=2)

        # # Pre-normalization
        # for index, experiment in enumerate(distributions.keys()):
        #     axes[0].step(bins[:-1], distributions[experiment][0], lw=1.5, c=colors[index], ls=linestyles[index])
        #     sample_count = len(data[experiment][intensity_metric])
        #     axes[0].annotate('N({0})={1}'.format(experiment, sample_count), (0.03, 0.96 - 0.1*index), xycoords='axes fraction', 
        #                      fontsize=8, va='top', ha='left')
        # axes[0].set_title('Before normalization')
            
        # # Post-normalization
        # for index, experiment in enumerate(weighted_distributions.keys()):
        #     axes[1].step(bins[:-1], weighted_distributions[experiment], lw=1.5, c=colors[index], ls=linestyles[index], label=experiment)
        # fig.tight_layout()
        # axes[1].set_title('After normalization')
        # fig.legend(ncols=2, frameon=False, loc='upper center', bbox_to_anchor=(0.55, 1.1))

        ######################################################################################################
        
        # Pre-normalization
        for index, experiment in enumerate(distributions.keys()):
            if experiment == CTL:
                axes[0].step(bins[:-1], weights.values(), lw=1.5, c=colors[index], ls=linestyles[index])
            else:
                axes[0].step(bins[:-1], 1/np.array(list(weights.values())), lw=1.5, c=colors[index], ls=linestyles[index])
                
            sample_count = len(data[experiment][intensity_metric])
            axes[0].annotate('N({0})={1}'.format(experiment, sample_count), (0.03, 0.96 - 0.1*index), xycoords='axes fraction', 
                             fontsize=8, va='top', ha='left')
        axes[0].set_title('Before normalization')
            
        # Post-normalization
        for index, experiment in enumerate(weighted_distributions.keys()):
            axes[1].step(bins[:-1], weighted_distributions[experiment], lw=1.5, c=colors[index], ls=linestyles[index], label=experiment)
        fig.tight_layout()
        axes[1].set_title('After normalization')
        fig.legend(ncols=2, frameon=False, loc='upper center', bbox_to_anchor=(0.55, 1.1))
    
    return sorted(bins), weights, weighted_distributions