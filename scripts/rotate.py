import numpy as np, scipy as sp, xarray as xr
import os, pickle, random

import matplotlib.patheffects as pe
import matplotlib, matplotlib.pyplot as plt
import utilities
import visualization

def snapshot_visual(snapshot, fields=['ucomp', 'vcomp', 'vort850', 'slp'], vertical_level=850):
    '''
    Plot to show TC-relevant quantities used for filtering (wind speed at 850 hPa, vorticity at 850 hPa, sea-level pressure).
    '''

    # Obtain the filtered snapshot to superimpose on the raw fields to illustrate where filtering occurs.
    _, filtered_snapshot = mask_TC_core(snapshot)
    
    # Construct the plot to accommodate all fields
    fig, axes = plt.subplots(figsize=(4*len(fields), 3), ncols=len(fields))
    # Iterate over all fields
    for dataset_index, dataset in enumerate([snapshot, filtered_snapshot]):
        for field_index, field in enumerate(fields):
            # Define local subplot
            ax = axes[field_index]
            # Subsample at a given vertical level if the 
            if 'pfull' in dataset[field].dims:
                subplot_dataset = dataset[field].sel(pfull=vertical_level, method='nearest')
            else:
                subplot_dataset = dataset[field]
                
            # Plot the data - if dataset_index is 0, plot raw data. Else, get masked data and plot.
            if dataset_index == 0:
                subplot_dataset.plot(ax=ax)
            else:
                # Mask the data to binary values and plot a contour line
                xr.where(~subplot_dataset.isnull(), 1, 0).plot.contour(ax=ax, levels=[0], colors=['k'], linestyles=['--'])
            
            # Figure formatting
            ax.set_aspect('equal')
            ax.set_title('{0} at {1} hPa\n{2}'.format(subplot_dataset.attrs['long_name'], vertical_level, subplot_dataset.time.values))
    
    fig.tight_layout()
    
def mask_TC_core(snapshot):
    '''
    Mask the TC core at a given timestamp (snapshot in time). 
    Filters out grid cells with mean sea-level pressure < 1 standard deviation from the mean.
    
    Rationale:
    - This was chosen to filter out likely-TC areas and eliminate areas with multiple TCs in a given domain.
    '''

    mask = snapshot['slp'] > (snapshot['slp'].mean() - snapshot['slp'].std())
    snapshot = snapshot.where(mask)
    
    return mask, snapshot

def distance_from_center(snapshot, diagnostic=False):
    '''
    For a given TC snapshot, get distance in meters from center for all points in domain. 
    This can also work for an array with a time axis, but not recommended.
    '''

    # Build grid composed of longitude and latitude values
    X, Y = np.meshgrid(snapshot.grid_xt.values, snapshot.grid_yt.values)
    # Obtain center coordinate. Note that 'center_lon' and 'center_lat' are repeated values, but a uniqueness filter is applied to be sure.
    if 'Dataset' in str(type(snapshot)):
        center = (np.unique(snapshot['center_lon'])[0], np.unique(snapshot['center_lat'])[0])
    else:
        center = (0, 0)
    # Iterate through all domain points using pairings from the meshgrid to calculate distances
    distance_1d = [utilities.coords_to_dist((point[0], point[1]), center) for point in np.c_[X.ravel(), Y.ravel()]]
    # Reshape distances to the input domain shape
    distance = np.reshape(distance_1d, X.shape)
    # Reshape distance to a DataArray for a cohesive data structure
    snapshot['distance'] = xr.DataArray(data=distance, 
                                        coords={'grid_yt': (['grid_yt'], snapshot.grid_yt.values), 
                                                'grid_xt': (['grid_xt'], snapshot.grid_xt.values)})

    if diagnostic:
        fig, ax = plt.subplots(figsize=(3, 3))
        snapshot['distance'].plot.contourf(ax=ax)
        ax.set_aspect('equal')

    return snapshot

def vertical_wind_shear(snapshot, full_snapshot, upper_level=200, lower_level=850, diagnostic=False):
    '''
    Derive vertical wind shear, preferably from a filtered wind field to remove TC effects in some way.
    '''

    # Obtain zonal wind shear
    du_dp = snapshot['ucomp'].sel(pfull=upper_level, method='nearest') - snapshot['ucomp'].sel(pfull=lower_level, method='nearest')
    # Obtain meridional wind shear
    dv_dp = snapshot['vcomp'].sel(pfull=upper_level, method='nearest') - snapshot['vcomp'].sel(pfull=lower_level, method='nearest')
    # Get shear magnitude and domain mean
    shear_magnitude = np.sqrt(du_dp**2 + dv_dp**2) # keep this separate as a non-averaged field for plotting
    shear_magnitude_mean = shear_magnitude.mean().values
    # Get shear direction mean based on domain-averaged wind components
    shear_direction_mean = np.mod(180 + (180/np.pi)*np.arctan2(du_dp.mean().values, dv_dp.mean().values), 360) * np.pi/180

    if diagnostic:
        ''' Plot the raw wind data used for the calculation. '''
        # Define fields and levels to iterate over
        wind_fields, pressure_levels = ['ucomp', 'vcomp'], [upper_level, lower_level]
        # Obtain number of rows and columns as a function of fields and levels
        nrows, ncols = len(pressure_levels), len(wind_fields)
        # Plot the base wind data
        fig, axes = plt.subplots(figsize=(4*ncols, 3*nrows), nrows=nrows, ncols=ncols)
        for row_index, row in enumerate(range(nrows)):
            for col_index, col in enumerate(range(ncols)):
                field, level = wind_fields[col_index], pressure_levels[row_index]
                ax = axes[row, col]
                snapshot[field].sel(pfull=level, method='nearest').plot(ax=ax)
                ax.set_title("{0} at {1} hPa".format(field, level))
                ax.annotate('{0:.1f} m s$^{{-1}}$'.format(snapshot[field].sel(pfull=level, method='nearest').mean().values), 
                            xy=(0.04, 0.04), xycoords='axes fraction', color='k',
                            fontsize=8, path_effects=[pe.Stroke(linewidth=2, foreground='white'), pe.Normal()])
                ax.set_aspect('equal')

        ''' Plot the wind shear and shear vector overlaid on a tracer plot. '''
        fig, axes = plt.subplots(figsize=(8, 3), ncols=2)
        
        # Shear plot
        ax_wind_shear = axes[0]
        shear_magnitude.plot(ax=ax_wind_shear)
        
        # Moisture plot - the plot shows moisture distributions at 2 different levels as a visual check on whether the shear vector aligns with tracer displacement
        ax_shear_vector = axes[1]
        # Plot moisture at both levels
        full_snapshot['sphum'].sel(pfull=850, method='nearest').plot.contourf(cmap='Greens', ax=ax_shear_vector, levels=16)
        full_snapshot['sphum'].sel(pfull=200, method='nearest').plot.contour(cmap='Blues', ax=ax_shear_vector, levels=8)
        # Obtain center coordinate. Note that 'center_lon' and 'center_lat' are repeated values, but a uniqueness filter is applied to be sure.
        center_lon, center_lat = (np.unique(snapshot['center_lon'])[0], np.unique(snapshot['center_lat'])[0])
        # Plot the shear vector
        length_factor = -5
        ax_shear_vector.arrow(center_lon, center_lat, 
                              dx=length_factor*np.sin(shear_direction_mean), 
                              dy=length_factor*np.cos(shear_direction_mean),
                              width=0.25, head_width=1, color='k', zorder=99)
        ax_shear_vector.set_title('blue: $q_{{{0}}}$, green: $q_{{{1}}}$'.format(upper_level, lower_level))

        # Annotate with shear information
        ax_shear_vector.annotate('{0:.1f} m s$^{{-1}}$ at {1:.1f} deg'.format(shear_magnitude_mean, shear_direction_mean * 180/np.pi), 
                                 xy=(0.04, 0.04), xycoords='axes fraction', color='k',
                                 fontsize=8, path_effects=[pe.Stroke(linewidth=2, foreground='white'), pe.Normal()])
        # Even out the aspect raio
        for ax in axes:
            ax.set_aspect('equal')

    # Plot average values over the domain
    return shear_magnitude_mean, shear_direction_mean

def rotate_shear_vector(snapshot, field, shear_vector_direction, shear_alignment_direction=90, nudging=False, center_param='slp', diagnostic=False):
    '''
    Rotates a given field from an xArray Dataset snapshot by the difference between the shear alignment direction and the shear vector direction.
    The shear alignment direction is the angle along which the shear vector should be oriented.
    '''

    # Get rotation angle
    rotation_angle = shear_alignment_direction - shear_vector_direction*180/np.pi
    # Rotate the snapshot, replace empty areas with nans
    rotation_array = sp.ndimage.rotate(snapshot[field].dropna('grid_xt', how='all').dropna('grid_yt', how='all'), rotation_angle, reshape=False, cval=np.nan)
    # Rotated snapshot coordinates - defines which coordinates to build the rotated snapshot coordinates with.
    # Only accepts dimensions in 3 dimensions (ignores 'phalf')
    coord_names = [coord for coord in snapshot[field].coords if coord in ['grid_xt', 'grid_yt', 'pfull']]
    # Create new Dataset
    rotation_snapshot = snapshot.copy()
    # Add the rotated data to the Dataset. Assumes the dimensions are axes, so the only dimensions associated with the coordinates is the axis dimension.
    rotation_snapshot['{0}_rotated'.format(field)] = xr.DataArray(data=rotation_array,
                                                         coords={coord_name: ([coord_name], snapshot[coord_name].values) for coord_name in coord_names},
                                                         dims=coord_names)
    
    
    ''' Use nudging to realign minimum sea-level pressure location with domain centerpoint. '''
    if nudging:
        # Define name for rotated center-detecting field
        center_field_rotated = '{0}_rotated'.format(center_param)
        # Get domain centerpoint
        domain_center_x, domain_center_y = int(round(len(snapshot.grid_xt) / 2)), int(round(len(snapshot.grid_yt) / 2))
        # Get minimum value of field and ignore all others to get lon and lat
        center_field_value = xr.where(abs(rotation_snapshot[center_field_rotated]) == abs(rotation_snapshot[center_field_rotated]).min().values, 
                                      rotation_snapshot[center_field_rotated], np.nan)
        center_field_value = center_field_value.dropna(dim='grid_xt', how='all').dropna(dim='grid_yt', how='all')
        # Trim off nans for all values
        center_field_lon, center_field_lat = center_field_value['grid_xt'].values[0], center_field_value['grid_yt'].values[0] 
        # Get array indices for the location of the field value, with matching occurring within a given tolerance
        center_field_x, center_field_y = [np.where(np.isclose(snapshot['grid_xt'], center_field_lon, atol=1e-04))[0][0], 
                                          np.where(np.isclose(snapshot['grid_yt'], center_field_lat, atol=1e-04))[0][0]]
        # Nudge the data such that the minimum rotated value matches the domain center
        grid_yt_step, grid_xt_step = 0.5, 0.625
        index_nudge_x, index_nudge_y = domain_center_x - center_field_x, domain_center_y - center_field_y
        rotation_snapshot['grid_xt'] = rotation_snapshot['grid_xt'] + index_nudge_x*grid_xt_step
        rotation_snapshot['grid_yt'] = rotation_snapshot['grid_yt'] + index_nudge_y*grid_yt_step
        if diagnostic:
            print("Domain center x-position: {0}, Domain center y-position {1}\nMinimum value x-position: {2}, Minimum value y-position: {3}".format(domain_center_x, domain_center_y, center_field_x, center_field_y))
            print("Nudging, x: {0}; Nudging, y: {1}".format(index_nudge_x, index_nudge_y))
    
    if diagnostic:
        fig, axes = plt.subplots(figsize=(6, 3), ncols=2)
        # Plot the original field
        axes[0].pcolormesh(snapshot['grid_xt'], snapshot['grid_yt'], snapshot[field])
        axes[0].set_title('Original')
        # Plot the rotated field
        axes[1].pcolormesh(rotation_snapshot['grid_xt'], rotation_snapshot['grid_yt'], rotation_snapshot['{0}_rotated'.format(field)])
        axes[1].set_title('Rotated {0:.2f} degrees'.format(rotation_angle))
        # Plot centerlines for perspective
        for ax in fig.axes:
            ax.axhline(snapshot['grid_yt'].isel(grid_yt=int(np.round(len(snapshot['grid_yt']) / 2))), c='k', zorder=99)
            ax.axvline(snapshot['grid_xt'].isel(grid_xt=int(np.round(len(snapshot['grid_xt']) / 2))), c='k', zorder=99)
            ax.axhline(rotation_snapshot['grid_yt'].isel(grid_yt=int(np.round(rotation_array.shape[0] / 2))), c='r', ls='--', zorder=99)
            ax.axvline(rotation_snapshot['grid_xt'].isel(grid_xt=int(np.round(rotation_array.shape[1] / 2))), c='r', ls='--', zorder=99)
            ax.set_aspect('equal')
    
    return rotation_snapshot

def rotate_snapshot(snapshot, field, shear_alignment_direction=90, shear_level_upper=200, shear_lever_lower=850, diagnostic=False):
    # Get original planar dimensions
    snapshot_grid_xt_original, snapshot_grid_yt_original = snapshot.grid_xt, snapshot.grid_yt
    # Define lists to collect all dimension arrays for a union later
    snapshot_grid_xt_arr, snapshot_grid_yt_arr = [], []
    # Drop nans for all fields
    for data_var in snapshot.data_vars:
        if 'grid_xt' in snapshot[data_var].dims and 'grid_yt' in snapshot[data_var].dims:
            grid_xt_trim = snapshot[data_var].dropna('grid_xt', how='all').dropna('grid_yt', how='all').grid_xt.values
            grid_yt_trim = snapshot[data_var].dropna('grid_xt', how='all').dropna('grid_yt', how='all').grid_yt.values
            snapshot_grid_xt_arr.append(grid_xt_trim)
            snapshot_grid_yt_arr.append(grid_yt_trim)
    # Get sets of coordinates
    grid_xt_trimmed = sorted(list(set.intersection(*map(set, snapshot_grid_xt_arr))))
    grid_yt_trimmed = sorted(list(set.intersection(*map(set, snapshot_grid_yt_arr))))
    # Finalize the trim
    snapshot = snapshot.sel(grid_xt=grid_xt_trimmed).sel(grid_yt=grid_yt_trimmed)
    # Add distances from TC center as a field to the Dataset
    snapshot = distance_from_center(snapshot, diagnostic=False)
    # Perform masking, first by removing areas with sea-level pressure below 1-standard deviation (filters out any TCs in domain), 
    # then by taking storm-centered winds field at some radial distance out
    _, masked_snapshot = mask_TC_core(snapshot)
    masked_snapshot = snapshot.where((masked_snapshot['distance'] >= 600e3) & ((masked_snapshot['distance'] <= 1400e3)))
    # Get shear vector information
    shear_magnitude, shear_direction = vertical_wind_shear(masked_snapshot, snapshot, upper_level=shear_level_upper, lower_level=shear_lever_lower, diagnostic=diagnostic)
    # Perform rotation
    snapshot = rotate_shear_vector(snapshot, field, shear_direction, shear_alignment_direction=shear_alignment_direction, diagnostic=diagnostic)

    return snapshot