''' Import packages. '''
# Analytical packages
import numpy as np, pandas as pd, xarray as xr
# Utilities
import importlib, datetime, os
# Visualization
import cartopy, cartopy.crs as ccrs
import matplotlib, matplotlib.pyplot as plt
# Internal methods
import derived, tc_processing, utilities, visualization

matplotlib.rcParams['mathtext.fontset'] = 'cm'

def storm_track(data, averaging_period=1, hourly_interval=3):
    track = data['track_output'].iloc[range(0, len(data['track_output']))]
    track['du_dt'] = track['max_wind'].diff()
    
    du_dt_period = {}
    for i in range(averaging_period*24 // hourly_interval, len(track)):
        dt = track.iloc[i]['time'] - track.iloc[i-8]['time']
        if dt == pd.Timedelta('1D'):
            # print(track.iloc[i]['time'], track.iloc[i]['max_wind'] - track.iloc[i-8]['max_wind'])
            du_dt_period[track.iloc[i]['time']] = track.iloc[i]['max_wind'] - track.iloc[i-8]['max_wind']
        else:
            du_dt_period[track.iloc[i]['time']] = np.nan
    du_dt_period = pd.DataFrame.from_dict(du_dt_period, columns=['du_dt_period'], orient='index').reset_index(names=['time'])
    
    track = track.merge(du_dt_period, on='time', how='left')

    data['track_output'] = track
    
    return data

def dh_dt(data):
    ddt = []
    for t, timestamp in enumerate(data['tc_model_output'].sortby('time', ascending=True).time.values):
        if (t >= 1) and (t < len(data['tc_model_output'].sortby('time', ascending=True).time.values)-1):
            dh_ = data['tc_model_output']['vi_h'].isel(time=t+1) - data['tc_model_output']['vi_h'].isel(time=t-1)
            dt_ = (data['tc_model_output']['time'].isel(time=t+1) - data['tc_model_output']['time'].isel(time=t-1)).data/np.timedelta64(1, 's')
            dh_dt_ = dh_/dt_
            dh_dt_['time'] = timestamp
            ddt.append(dh_dt_)
        else:
            ddt.append(data['tc_model_output']['vi_h'].isel(time=t)*np.nan)
    data['tc_model_output']['dh_dt'] = xr.concat(ddt, dim='time')

    return data

def storm_period_binning(data, diagnostic=False):

    # Define averaging parameters
    averaging_period = 1 # days
    hourly_interval = 3 # hours
    averaging_index_span = averaging_period*24 // hourly_interval
    # Get track output from the data
    data = storm_track(data, averaging_period=averaging_period)
    # Shorthand for track data
    track = data['track_output']
    track = track.reset_index()
    # Get storm ID
    storm_id = track['storm_id'].unique()[0]

    # Identify timestamp corresponding to LMI
    lmi = track['time'].loc[track['max_wind'] == track['max_wind'].max()]
    print('LMI occurs at {0}'.format(lmi))

    # Define absolute intensity bins (same for all storms)
    bins = [-np.inf, -2, 2, np.inf]
    # Define period names
    period_names = ['intensification', 'weakening', 'peak']

    # Initialize dictionary
    periods = {}
    # Sample count
    num_samples = 3
    # Define times for future visualization of samplin
    sampling_times = {period: {} for period in period_names}
    # DEfine sampling bins in decreasing order
    sampling_bins = [bin for bin in track.groupby(pd.cut(track['du_dt_period'], bins))][::-1]
    sampling_bins = [sampling_bins[0], sampling_bins[-1], sampling_bins[1]]
    # Iterate over each bin
    for i, bin in enumerate(sampling_bins):
        # Unpack data (bin and data)
        bin_name, bin_data = bin
        print(period_names[i], bin_name)
        
        # New method pseudocode
        # 1. identify LMI. This will tell me when I should look for peak/steady-state
        # 2. identify intensification period. This must be before LMI.
        # 3. identify peak/steady-state. This must be after intensification.
        # 4. identify weakening period. This must be after intensification and peak.
        
        # Iterate through data to find the most representative data point for each bin
        # If 'weakening', get most negative du_dt; if 'peak', get smallest du_dt; if 'strengthening', get maximum du_dt
        if period_names[i] == 'intensification':
            times = bin_data.sort_values('du_dt_period', ascending=False)
            times = times.loc[times['time'] < lmi.values[0]] 
        elif period_names[i] == 'weakening':
            times = bin_data.sort_values('du_dt_period', ascending=True)
            times = times.loc[times['time'] > lmi.values[0]]
        elif period_names[i] == 'peak':
            bin_data['du_dt_period'] = bin_data['du_dt_period'].abs()
            times = bin_data.sort_values('du_dt_period', ascending=True)
            times = times.loc[(times['time'] < min(periods['weakening']['time'])) & (times['time'] > min(periods['intensification']['time']))]['time']
        
        # If more sample are requested than available, use all available
        # Sampling method 1: get top 3 instances matching period-specific criteria
        times = times.iloc[0:-1] if len(times) < num_samples else times.iloc[0:num_samples]
        # Sampling method 2: get 3 timestamps surrounding the timestamp most matching the period-matching criteria
        timestamp_index = times.index.values[0]        
        times = track['time'].iloc[[timestamp_index-1, timestamp_index, timestamp_index+1]]
        
        # Populate sample times for visualization
        sampling_times[period_names[i]] = [t for t in times.values]
        
        # Adjust timestamps from Pandas datetime to cftime formats
        timestamps = [utilities.time_adjust(model='HIRAM', timestamp=timestamp) for timestamp in times]
        
        # Initialize and populate datasets
        planar, vertical = [], []
        for timestamp in timestamps:
            planar.append(data['tc_model_output'].sel(time=timestamp))
            vertical.append(data['tc_vertical_output'].sel(time=timestamp))
        
        # Concatenate data
        planar, vertical = xr.concat(planar, dim='time'), xr.concat(vertical, dim='time')
        
        # Build dictionary for the  given intensification period
        periods[period_names[i]] = {'bin': bin_name, 'time': times}
        periods[period_names[i]]['data'] = {'planar': planar, 'vertical': vertical}
        # Diagnostic print statement
        if diagnostic:
            print('{0}\ndu/dt bin: {1} m s^-1; timestamps: \n'.format(period_names[i], timestamps))
    
    return periods, sampling_times, storm_id

def lifetime_plots(data, sampling_times, storm_id):
    
    track = data['track_output']

    fig, axes = plt.subplots(figsize=(12, 2), ncols=3)
    
    axes[0].plot(track['time'], track['max_wind'])
    axes[0].xaxis.set_major_locator(matplotlib.dates.AutoDateLocator())
    axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=30, ha='right')
    axes[0].set_title('Maximum 10-meter winds', fontsize=10)
    axes[0].annotate('{0}'.format(storm_id), xy=(0.03, 0.95), xycoords='axes fraction', ha='left', va='top', fontsize=10)
    # Plot sampling points
    period_colors = ['r', 'b', 'g']
    for p, period in enumerate(sampling_times.keys()):
        for t in sampling_times[period]:
            axes[0].axvline(t, alpha=0.5, c=period_colors[p])
    
    axes[1].axhline(0, lw=0.5, c='k')
    axes[1].plot(track['time'], track['du_dt_period'])
    axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=30, ha='right')
    axes[1].set_title('{0}-h change in 10-meter winds'.format(24), fontsize=10)
    
    axes[2].axvline(0, lw=0.5, c='k')
    track['du_dt_period'].diff().hist(ax=axes[2], grid=False, bins=20)
    axes[2].set_title('Histogram of 24-h change in 10-meter winds', fontsize=10)
    
def term_selection(periods):
    
    for period in periods.keys():
        print('Processing data for the {0} period...'.format(period))
        periods[period]['data']['budget'] = {'dh_dt': [], 'thf': [], 'q_rad': [],}
        planar = periods[period]['data']['planar']
    
        for t in planar.time.values:
            p = planar.sel(time=t)
            p.load()
    
            # Get turbulent heat fluxes (THF)
            q_h = p['shflx'].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            q_l = p['lhflx'].dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
    
            if q_h.shape[0] > 30:
                q_h = q_h.isel(grid_yt=slice(0, 30))
            if q_h.shape[1] > 30:
                q_h = q_h.isel(grid_xt=slice(0, 30))
            
            if q_l.shape[0] > 30:
                q_l = q_l.isel(grid_yt=slice(0, 30))
            if q_l.shape[1] > 30:
                q_l = q_l.isel(grid_xt=slice(0, 30))
            
            periods[period]['data']['budget']['thf'].append(q_h + q_l)
            # Get radiative convergence in column (q_rad)
            q_sw = (-p['swup_toa'] + p['swdn_toa'] + p['swup_sfc'] - p['swdn_sfc']).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            q_lw = (p['lwup_sfc'] - p['olr']).dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
            
            if q_sw.shape[0] > 30:
                q_sw = q_sw.isel(grid_yt=slice(0, 30))
            if q_l.shape[1] > 30:
                q_sw = q_sw.isel(grid_xt=slice(0, 30))
            
            if q_lw.shape[0] > 30:
                q_lw = q_lw.isel(grid_yt=slice(0, 30))
            if q_lw.shape[1] > 30:
                q_lw = q_lw.isel(grid_xt=slice(0, 30))
                
            periods[period]['data']['budget']['q_rad'].append(q_sw + q_lw)
            # Get MSE time tendency
            ddt = p['dh_dt'].interp(grid_xt=periods[period]['data']['budget']['thf'][0].grid_xt, 
                                    grid_yt=periods[period]['data']['budget']['thf'][0].grid_yt)
            periods[period]['data']['budget']['dh_dt'].append(ddt)        
    
            # Turn on for troubleshooting
            # print(q_h.shape, q_l.shape, q_sw.shape, q_lw.shape, ddt.shape)
    
        # Iterate over each key to get time mean and build composite xArray
        for key in periods[period]['data']['budget'].keys():
            # Make shorthand for iterand term
            kd = periods[period]['data']['budget'][key]
            # Extract data values, stack arrays on each other, and get mean
            coord_x = kd[0].grid_xt - kd[0].isel(grid_xt=len(kd[0].grid_xt)//2).grid_xt
            coord_y = kd[0].grid_yt - kd[0].isel(grid_yt=len(kd[0].grid_xt)//2).grid_yt
            # Extract data values, stack arrays on each other, and get mean
            kdd = np.nanmean(np.stack([kd[i] for i in range(0, len(kd))]), axis=0)
            # Construct DataArray
            periods[period]['data']['budget'][key] = xr.DataArray(data=kdd,
                                                                  dims=['TC_xt', 'TC_yt'],
                                                                  coords={'TC_xt': (['TC_xt'], coord_x.data),
                                                                          'TC_yt': (['TC_yt'], coord_y.data)})
        # Get residual
        periods[period]['data']['budget']['-div_u_h'] = periods[period]['data']['budget']['dh_dt'] - periods[period]['data']['budget']['thf'] - periods[period]['data']['budget']['q_rad']

    return periods

def main(storms, troublshooting=False):
    output = []
    for storm in storms:
        print(storm)
        filename = storm
        # Read data int
        data = pd.read_pickle(filename)
        # Derive MSE time tendency
        data = dh_dt(data)
        # Sample the TC at points where intensification, steady-state (peak), and weakening occur
        if troublshooting:
            periods, sampling_times, storm_id = storm_period_binning(data)
        else:
            try:
                periods, sampling_times, storm_id = storm_period_binning(data)
            except:
                print('Not enough data found, proceeding to the next one...')
                continue
        # Generate plots to show where sampling occurs relative to the storm lifetime
        lifetime_plots(data, sampling_times, storm_id)
        # Build the MSE time tendency budget
        periods = term_selection(periods)
        output.append(periods)
    return output

N = 15
dirname = '/projects/GEOCLIM/gr7610/analysis/tc_storage/individual_TCs/processed'
outputs = {'control': [], 'swishe': []}
for experiment in outputs.keys():
    storms = [os.path.join(dirname, filename) for filename in os.listdir(dirname) if 'HIRAM-8xdaily' in filename and experiment in filename][:N]
    outputs[experiment] = main(storms)

outputs['swishe-control'] = []
for n in range(0, min(len(outputs['control']), len(outputs['swishe']))):
    temp = {}
    for period in outputs['control'][n].keys():
        temp[period] = {'data': {'budget': {}}}
        for term in outputs['control'][n][period]['data']['budget'].keys():
            # Derive raw difference from experiments
            delta = outputs['swishe'][n][period]['data']['budget'][term].values - outputs['control'][n][period]['data']['budget'][term].values
            # Construct DataArray for difference
            temp[period]['data']['budget'][term] = xr.DataArray(data=delta,
                                                                dims=['TC_xt', 'TC_yt'],
                                                                coords={'TC_xt': (['TC_xt'], outputs['control'][n][period]['data']['budget'][term].TC_xt.data),
                                                                        'TC_yt': (['TC_yt'], outputs['control'][n][period]['data']['budget'][term].TC_yt.data)})

    outputs['swishe-control'].append(temp)
    
importlib.reload(visualization)

experiment = 'control'

for experiment in outputs.keys():

    budget_terms = {'dh_dt': '$\partial_t h$',
                    'thf': '$\ddot{q}_{\mathrm{surf}} = \\ddot{q}_h$ + $\\ddot{q}_l$',
                    'q_rad': '$\ddot{q}_{\mathrm{rad}} = \\ddot{q}_{SW}$ + $\\ddot{q}_{LW}$',
                    '-div_u_h': '$-\\nabla_h \cdot \\left(\\widehat{ \\mathbf{u}h }\\right) $'}

    for period in ['intensification', 'peak', 'weakening']:

        sample = outputs[experiment][0][period]['data']['budget']
        
        ncols = len(sample.keys())
        fig, gs = plt.figure(figsize=(2*ncols, 4), dpi=96), matplotlib.gridspec.GridSpec(nrows=1, ncols=ncols)

        for i, term in enumerate(sample.keys()):
            print(period, term)
            ax = fig.add_subplot(gs[0, i])
            ax.set_title(budget_terms[term], y=1.35)

            out = np.nanmean(np.stack([outputs[experiment][i][period]['data']['budget'][term] for i in range(0, len(outputs[experiment]))]), axis=0)
        
            norm, cmap = visualization.norm_cmap(out, term,
                                                num_bounds=16, extrema=None, white_adjust=False) 
            
            im = ax.pcolormesh(sample[term].TC_xt, 
                            sample[term].TC_yt, 
                            out, norm=norm, cmap=cmap)
            
            cax = ax.inset_axes([0, 1.05, 1, 0.05])

            if sum(np.isnan(norm.boundaries)) == len(norm.boundaries):
                norm = matplotlib.colors.Normalize()
                colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
            else:
                colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax, orientation='horizontal')
            cax.xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
            cax.xaxis.tick_top()
            cax.xaxis.set_label_position('top')
        
            if i > 0:
                ax.set_yticklabels([])
            ax.set_aspect('equal')

        fig.tight_layout()
        fig.suptitle(period, y=1)