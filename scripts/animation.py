# %%

import multiprocessing, os
import numpy as np, pandas as pd
import xarray as xr

import cartopy, cartopy.crs as ccrs
from matplotlib import animation
import matplotlib, matplotlib.pyplot as plt

# Suppress warnings
import warnings
warnings.filterwarnings("ignore")

from colormap_normalization import get_cmap, norm_cmap

def retrieve_tracked_TCs(dirname, storm_type, year_range=(101, 125)):

    '''
    Function to collect tracked TC data and add derived data, such as duration and storm speed.
    
    Input(s):
    - dirname (str):              name of directory containing files of interest
    - storm_type (str):           type of storm to evaluate from TC tracks data ("TS" for all storms or "C15w" for hurricanes)
    - year_range (tuple of ints): 2-element tuple with a start and end year
    Output(s):
    - data (Pandas DataFrame):    Pandas DataFrame with tracked TC data
    '''
    
    ''' File collection. '''
    # Get filenames for all files within the specified directory 
    # Filenames will correspond to the determined storm type    
    fnames = [[os.path.join(dirname, file, 'Harris.TC', f) for f in os.listdir(os.path.join(dirname, file, 'Harris.TC')) 
               if '{0}.world'.format(storm_type) in f]
               for file in sorted(os.listdir(dirname))]
    # Compress 2D list to 1D list
    fnames = [item for sublist in fnames for item in sublist]

    # Select files with dates within 'year_range'
    # Note: the '+ 1900' is added because tracked TCs are on the 2000 year range, whereas model output is on the 100 year range
    # Note: conditional added as exception for provisional FLOR data ahead of 10th NE Tropical Workshop
    if dirname == '/tigress/GEOCLIM/grios/HIRAM/exp/CTL1990_v201905/analysis_lmh/cyclones_gav_ro110_1C_330k':
        fnames = [f for f in fnames]
    else:
        fnames = [f for f in fnames 
                  if min(year_range) + 1900 <= pd.to_datetime(f.split('.')[-2].split('-')[0]).year <= max(year_range) + 1900]
    
    # Concatenate all tracked TC data from the filename list
    data = pd.concat([lmh_parser(os.path.join(dirname, fname)) for fname in fnames])
    
    ''' Derived track-based data algorithm. Storm-specific derived properties will be generated in here. '''
    
    # Initialize empty duration column to populate iteratively
    data[['duration', 'speed', 'direction']] = np.nan
    # Initialize list to populate iteratively for each storm, then concatenate
    storms = []
    # Iterate through each unique storm (identify by 'storm_id') and get duration
    for storm_id in data['storm_id'].unique():
        # Define iterand storm
        storm = data.loc[data['storm_id'] == storm_id].copy().reset_index(drop=True)
        
        ''' Duration derivation. '''
        # Get difference between minimum and maximum timestamps
        dt = (storm['time'].max() - storm['time'].min())
        # Convert difference timedelta into hours
        dt = dt.days + dt.seconds/86400
        # Add duration to the outer DataFrame for the corresponding storm
        data.loc[data['storm_id'] == storm_id, 'duration'] = dt
        # Re-define iterand storm to incorporate duration
        storm = data.loc[data['storm_id'] == storm_id].copy().reset_index(drop=True)
        
        ''' Velocity (speed, direction) derivation. '''
        # Initialize dictionary for preliminary storage. Will be reassigned into the DataFrame by the join() method using time as the matching criterion.
        velocity = {'time': [storm.iloc[0]['time']], 'speed': [np.nan], 'direction': [np.nan]}
        # Iterate over all of the iterand storm timestamps
        for i in range(1, len(storm)):
            # Define coordinates for two points considered (i, i-1)
            lon_a, lat_a = [storm.iloc[i-1]['lon'], storm.iloc[i-1]['lat']]
            lon_b, lat_b = [storm.iloc[i]['lon'], storm.iloc[i]['lat']]
            # Determine timedelta between points (i, i-1)
            dt = storm.iloc[i]['time'] - storm.iloc[i-1]['time']
            # Derive speed (distance / time in m s^-1)
            speed = coords_to_dist((lon_b, lat_b), (lon_a, lat_a))/dt.seconds
            # Get changes in longtiude and latitude
            dlon, dlat = lon_b - lon_a, lat_b - lat_a
            # Derive direction relative to north (range of 0 to 360)
            direction = 180*np.arctan(dlon/dlat)/np.pi % 360
            # Append quantities to the 'velocity' dictionary
            velocity['time'].append(storm.iloc[i]['time'])    
            velocity['speed'].append(speed)    
            velocity['direction'].append(direction)
        # Build DataFrame
        velocity = pd.DataFrame(velocity)
        # Re-cast time column as a datetime object
        velocity['time'] = pd.to_datetime(velocity['time'])
        # Merge the storm and velocity DataFrames
        storm = storm.merge(velocity, how='left', on='time', suffixes=['_x', None]).drop(columns={'speed_x', 'direction_x'}).reset_index(drop=True)
        # Append to the list for future concatenation
        storms.append(storm)
        
    # Concatenate DataFrames
    data = pd.concat(storms)   
    # Rename columns for future addition into xArray Dataset, and reset index
    data = data.rename(columns={'lon': 'center_lon', 'lat': 'center_lat', 'flag': 'core_temp', 'slp': 'min_slp'}).reset_index(drop=True)
        
    return data

def lmh_parser(path):
    
    ''' 
    This method parses through text files from Lucas Harris' run outputs (held in directories titled 'analysis_lmh') 
    and produces an output DataFrame. 
    
    Input(s):
    - path (str):            path containing raw tracker data from Lucas Harris' runs.
    Output(s):
    - df (Pandas DataFrame): Pandas DataFrame containing tracked TC data
    '''
    
    # Create file object instance
    fobj = open(path, 'r').readlines()
    # Initialize dictionary to hold data
    data = {'storm_num': {}}
    # Initialize storm counter
    count = 1
    # Iterate through text file
    for line in fobj:
        # Extract information from the line
        content = line.strip()
        # Creates new storm-specific dict in the parent dict. The '+++' demarcates a new storm.
        if '+++' in line:
            storm_num = '{0:04d}'.format(count)
            data['storm_num'][storm_num] = {'storm_id': [], 'time': [], 'lon': [], 'lat': [], 'slp': [], 'max_wnd': [], 'flag': []}
            count += 1
        # Populates the storm-specific dict
        else:
            storm_num = '{0:04d}'.format(count-1) 
            tc_info = [x for x in content.split(' ') if x]
            year = tc_info[0][0:4] # get 4-digit year
            data['storm_num'][storm_num]['storm_id'].append('{0}-{1:04d}'.format(year, count-1))
            data['storm_num'][storm_num]['time'].append(tc_info[0])
            data['storm_num'][storm_num]['lon'].append(tc_info[1])
            data['storm_num'][storm_num]['lat'].append(tc_info[2])
            data['storm_num'][storm_num]['slp'].append(tc_info[3])
            data['storm_num'][storm_num]['max_wnd'].append(tc_info[4])
            data['storm_num'][storm_num]['flag'].append(tc_info[5])
    
    try:
        # Converts the dictionary into a DataFrame
        df = pd.concat({k: pd.DataFrame(v).T for k, v in data.items()}, axis=1)['storm_num']
        df = df.explode(df.columns.to_list()).reset_index().rename(columns={'index': 'storm_num'})
        # Re-cast column data types
        df = df.astype({'lon': 'float', 'lat': 'float', 'slp': 'float', 'max_wnd': 'float', 'flag': 'float'})
    except:
        df = pd.DataFrame(columns=['storm_id', 'time', 'lon', 'lat', 'slp', 'max_wnd', 'flag'])
    
    ''' DataFrame refinement. '''
    # Remove cold-core data points (flag == -1)
    # df = df.loc[df['flag'] != -1].reset_index(drop=True)
    # Convert timestamps to datetime objects
    df['time'] = pd.to_datetime(df['time'], format='%Y%m%d%H')
    
    return df

def coords_to_dist(a, b):
    ''' Convert coordinates to distance in meters. '''
    
    R = 6371e3
    
    lon_a, lat_a = np.array(a)*np.pi/180
    lon_b, lat_b = np.array(b)*np.pi/180
    
    dlon, dlat = lon_b - lon_a, lat_b - lat_a
    
    a = np.sin(dlat/2)**2 + np.cos(lat_a)*np.cos(lat_b)*np.sin(dlon/2)**2    
    c = 2*np.arctan2(np.sqrt(a), np.sqrt(1-a))
    
    distance = R*c
    
    return distance

def access(model_name='HIRAM', experiment='control', storm_category='C15w'):
    
    model_run_dirname = '/tiger/scratch/gpfs/GEOCLIM/gr7610/{0}/work'.format(model_name)
    
    file_type = 'atmos_4xdaily'
    dirname_modifier = ''
    if model_name in ['AM2.5', 'HIRAM']:
        npes = '540PE' 
    elif model_name in ['FLOR']:
        npes = '576PE'
        dirname_modifier = 'v201905_'
    else:
        npes = '1080PE'
    
    experiment_names = {'control': 'CTL1990s_{0}tigercpu_intelmpi_18_{1}'.format(dirname_modifier, npes),
                        'swishe': 'CTL1990s_swishe_{0}tigercpu_intelmpi_18_{1}'.format(dirname_modifier, npes)}
    
    storm_categories = ['TS', 'C15w']
    
    # Populate dictionary
    experiment_dirname = os.path.join(model_run_dirname, experiment_names[experiment], 'POSTP')
    track_dirname = [os.path.join(model_run_dirname, experiment_names[experiment], 'analysis_lmh', dirname) 
                     for dirname in os.listdir(os.path.join(model_run_dirname, experiment_names[experiment], 'analysis_lmh')) 
                     if 'cyclones_gav' in dirname][0]

    track_data = retrieve_tracked_TCs(track_dirname, storm_category)
    
    return track_data, experiment_dirname

def timestamp_conversion(input_time, mode='to_cftime', calendar='noleap'):
    import cftime
    cftimestamp = cftime.datetime(year=input_time.year-1900, month=input_time.month, day=input_time.day, hour=input_time.hour,
                                  calendar=calendar)
    return cftimestamp

def get_data(model_name='HIRAM', experiment='control', storm_category='C15w', storm_num=None):

    track_data, experiment_dirname = access(model_name, experiment, storm_category)

    valid_lon = False
    res = 30

    if not storm_num:
        while valid_lon is False:
            storm_num = np.random.choice(track_data['storm_id'].unique())
            print(storm_num)
            storm_data = track_data.loc[track_data['storm_id'] == storm_num]
            lon_filter = storm_data['center_lon'][(storm_data['center_lon'] < (180 + res/2)) & (storm_data['center_lon'] > (180 - res/2))]
            if len(lon_filter) == 0:
                valid_lon = True
    else:
        storm_data = track_data.loc[track_data['storm_id'] == storm_num]
        
    # Get the times corresponding to the storm
    timestamps = min(storm_data['time']), max(storm_data['time'])
    
    ''' Get 6-hourly data corresponding to storm. '''
    # Get filename matching corresponding model output.
    model_filename = [os.path.join(experiment_dirname, filename) for filename in os.listdir(experiment_dirname)
                      if ('atmos_4xdaily' in filename) and (str(min(timestamps).year - 1900) in filename.split('/')[-1][0:4])][0]    
    print([os.path.join(experiment_dirname, filename) for filename in os.listdir(experiment_dirname)
                      if ('atmos_4xdaily' in filename) and (str(min(timestamps).year - 1900) in filename)])
    # Convert timestamps to cftime format for output indexing.
    timestamps = [timestamp_conversion(t) for t in timestamps]
    # Get time indices to allow for storm loading before/after tracker start/end
    model_data = xr.open_dataset(model_filename).sel(time=slice(min(timestamps), max(timestamps))).load()

    ''' Filter data to storm extent. '''
    extent_lon_min, extent_lon_max = storm_data['center_lon'].min(), storm_data['center_lon'].max()
    extent_lat_min, extent_lat_max = storm_data['center_lat'].min(), storm_data['center_lat'].max()
    extent = [extent_lon_min - res, extent_lon_max + res, extent_lat_min - res, extent_lat_max + res]
    model_data = model_data.sel(grid_xt=slice(extent[0], extent[1]),
                                grid_yt=slice(extent[2], extent[3]))
    
    return storm_data, model_data

storm_data, model_data = get_data(model_name='HIRAM', experiment='control')

plt.rcParams["animation.html"] = "jshtml"
plt.ioff()
plt.cla()

param = 'olr'
proj, proj_pc = ccrs.PlateCarree(), ccrs.PlateCarree()

fig = plt.figure()
gs = matplotlib.gridspec.GridSpec(nrows=1, ncols=2, width_ratios=(1, 0.03), wspace=0)

geo_ax = fig.add_subplot(gs[0, 0], projection=proj)
cax_container = fig.add_subplot(gs[0, 1])

model_data['U'] = np.sqrt(model_data['u_ref']**2 + model_data['v_ref']**2)

def init_func():
    frame = 0
    geo_ax.clear()

    levels = 24
    norm, cmap = norm_cmap(model_data, param, bounds=True, n_bounds=levels+1)
    im = geo_ax.contourf(model_data.grid_xt, model_data.grid_yt, model_data[param].isel(time=frame),
                           norm=norm, cmap=cmap, levels=levels)
    geo_ax.contour(model_data.grid_xt, model_data.grid_yt, model_data['U'].isel(time=frame),
                   levels=[15, 29], cmap='Blues')

    res = 15
    center_lat, center_lon = storm_data.iloc[frame]['center_lat'], storm_data.iloc[frame]['center_lon']

    geo_ax.set_extent([center_lon - res, center_lon + res, center_lat - res, center_lat + res])
    
    gridlines = geo_ax.gridlines(ls='--', alpha=0.5)
    gridlines.bottom_labels = True
    
    geo_ax.coastlines()

    storm_id = storm_data.iloc[0]['storm_id']
    timestamp = model_data[param].isel(time=frame)['time'].values
    max_wind, slp = storm_data.iloc[frame]['max_wnd'], storm_data.iloc[frame]['min_slp']
    title = 'Storm ID: {0}\nTime:{1}; Coordinates: ({2:.1f}, {3:.1f})\nMaximum winds: {4:.1f} m s$^{{-1}}$; Minimum pressure: {5:.1f} hPa'.format(storm_id, timestamp,
                                                                                                                                          center_lat, center_lon, max_wind, slp)
    geo_ax.set_title(title, ha='left', x=0, fontsize=9)

    cax_container.clear()
    cax_container.set_axis_off()
    cax = cax_container.inset_axes([-2, 0, 1, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
    colorbar_label = '{0} [{1}]'.format(model_data[param].attrs['long_name'], model_data[param].attrs['units'])
    colorbar.set_label(colorbar_label, labelpad=20, rotation=270)

def animate(iter_frame):

    storm_frame, model_frame = iter_frame, iter_frame
    geo_ax.clear()

    levels = 24
    norm, cmap = norm_cmap(model_data, param, bounds=True, n_bounds=levels+1)
    im = geo_ax.contourf(model_data.grid_xt, model_data.grid_yt, model_data[param].isel(time=model_frame),
                           norm=norm, cmap=cmap, levels=levels)
    geo_ax.contour(model_data.grid_xt, model_data.grid_yt, model_data['U'].isel(time=model_frame),
                   levels=[15, 29], cmap='Blues')

    res = 15
    center_lat, center_lon = storm_data.iloc[storm_frame]['center_lat'], storm_data.iloc[storm_frame]['center_lon']

    geo_ax.set_extent([center_lon - res, center_lon + res, center_lat - res, center_lat + res])
    
    gridlines = geo_ax.gridlines(ls='--', alpha=0.5)
    gridlines.bottom_labels = True
    
    geo_ax.coastlines()

    storm_id = storm_data.iloc[0]['storm_id']
    timestamp = model_data[param].isel(time=model_frame)['time'].values
    max_wind, slp = storm_data.iloc[storm_frame]['max_wnd'], storm_data.iloc[storm_frame]['min_slp']
    title = 'Storm ID: {0}\nTime:{1}; Coordinates: ({2:.1f}, {3:.1f})\nMaximum winds: {4:.1f} m s$^{{-1}}$; Minimum pressure: {5:.1f} hPa'.format(storm_id, timestamp,
                                                                                                                                          center_lat, center_lon, max_wind, slp)
    geo_ax.set_title(title, ha='left', x=0, fontsize=9)

    cax_container.clear()
    cax_container.set_axis_off()
    cax = cax_container.inset_axes([-2, 0, 1, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
    colorbar_label = '{0} [{1}]'.format(model_data[param].attrs['long_name'], model_data[param].attrs['units'])
    colorbar.set_label(colorbar_label, labelpad=20, rotation=270)

def artist_animate(iter_frame):

    storm_frame, model_frame = iter_frame, iter_frame
    geo_ax.clear()

    levels = 24
    norm, cmap = norm_cmap(model_data, param, bounds=True, n_bounds=levels+1)
    im = geo_ax.contourf(model_data.grid_xt, model_data.grid_yt, model_data[param].isel(time=model_frame),
                           norm=norm, cmap=cmap, levels=levels)
    geo_ax.contour(model_data.grid_xt, model_data.grid_yt, model_data['U'].isel(time=model_frame),
                   levels=[15, 29], cmap='Blues')

    res = 15
    center_lat, center_lon = storm_data.iloc[storm_frame]['center_lat'], storm_data.iloc[storm_frame]['center_lon']

    geo_ax.set_extent([center_lon - res, center_lon + res, center_lat - res, center_lat + res])
    
    gridlines = geo_ax.gridlines(ls='--', alpha=0.5)
    gridlines.bottom_labels = True
    
    geo_ax.coastlines()

    storm_id = storm_data.iloc[0]['storm_id']
    timestamp = model_data[param].isel(time=model_frame)['time'].values
    max_wind, slp = storm_data.iloc[storm_frame]['max_wnd'], storm_data.iloc[storm_frame]['min_slp']
    title = 'Storm ID: {0}\nTime:{1}; Coordinates: ({2:.1f}, {3:.1f})\nMaximum winds: {4:.1f} m s$^{{-1}}$; Minimum pressure: {5:.1f} hPa'.format(storm_id, timestamp,
                                                                                                                                          center_lat, center_lon, max_wind, slp)
    geo_ax.set_title(title, ha='left', x=0, fontsize=9)

    cax_container.clear()
    cax_container.set_axis_off()
    cax = cax_container.inset_axes([-2, 0, 1, 1])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm, cmap), cax=cax)
    colorbar_label = '{0} [{1}]'.format(model_data[param].attrs['long_name'], model_data[param].attrs['units'])
    colorbar.set_label(colorbar_label, labelpad=20, rotation=270)
    
    return fig


# Animation generation and save settings
animation_type = 'artist'
frame_offset = 1
max_frame = 12
frame_limit = max_frame if (len(storm_data)-1 > max_frame) else len(storm_data)-1

ims = [artist_animate(i) for i in range(0, max_frame)]

if animation_type == 'func':
    matplotlib.animation.FuncAnimation(fig, animate, frames=range(frame_offset, frame_limit), init_func=init_func, interval=150, blit=False)
else:
    matplotlib.animation.ArtistAnimation(fig, ims, interval=150, blit=False)

