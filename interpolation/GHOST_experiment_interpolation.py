#WRITTEN BY DENE BOWDALO

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#GHOST_experiment_interpolation.py

#module which interpolates experiment output to surface observations 

###--------------------------------------------------------------------------------------------------###
###IMPORT MODULES
###--------------------------------------------------------------------------------------------------###

from calendar import monthrange
import cartopy.crs as ccrs
import glob
import os
from netCDF4 import Dataset
import numpy as np
import pyproj
import os
from scipy import spatial
from shapely.geometry import Polygon, Point
import sys
import time
import unit_converter

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

def process_monthly_interpolated_netCDF(experiment_to_process, exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth):

    '''function that processes monthly interpolated netCDF for an experiment with reference to observational points on the surface'''

    print('Interpolating in %s'%(yearmonth))

    #get year/month string
    year = yearmonth[:4]
    month = yearmonth[4:]

    #get relevant observational file
    obs_file = glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/%s/%s/%s/%s_%s*.nc'%(GHOST_network_to_interpolate_against, temporal_resolution_to_output, speci_to_process, speci_to_process, yearmonth))[0]

    #get relevant model files
    model_files = np.sort(glob.glob('%s/%s/%s/%s/%s_%s*.nc'%(exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, speci_to_process, yearmonth)))

    #get number of days in month processing
    days_in_month = monthrange(int(year),int(month))[1]

    #create descriptive temporal resolution variable
    if temporal_resolution_to_output == 'hourly':
        descriptive_temporal_resolution = 'hours'
    elif temporal_resolution_to_output == 'daily':
        descriptive_temporal_resolution = 'days'
    elif temporal_resolution_to_output == 'monthly':
        descriptive_temporal_resolution = 'months'

    #create time variable
    if temporal_resolution_to_output == 'hourly':
        yearmonth_time = np.arange(0,days_in_month*24.0)
    elif temporal_resolution_to_output == 'daily':
        yearmonth_time = np.arange(0,days_in_month)
    #monthly output resolution is initially in hours for processing reasons, this is modified later 
    elif temporal_resolution_to_output == 'monthly':
        yearmonth_time = np.arange(0,days_in_month*24.0)

    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    #read model data

    #take first valid model file in month and get grid dimension/coordinate information 
    #put initial object read in a  try/except to handle reading of corrupted files
    #iterate through files until have read a valid file
    #if do not read a valid file, skip month

    for model_file_ii, model_file in enumerate(model_files):
        try:
            #load instance of model file netCDF
            mod_nc_root = Dataset(model_file)
            
            #get all variable names in file
            mod_nc_varnames = list(mod_nc_root.variables.keys())

            #get instance of species variable
            mod_speci_obj = mod_nc_root[speci_to_process]
        
            #get species variable units
            mod_speci_units = mod_speci_obj.units
            
            #get model grid type
            mod_grid_type = mod_speci_obj.grid_mapping
            
            #get x/y grid dimension variable names
            grid_dimensions = mod_speci_obj.dimensions
    
        except:
            #if have got to last file of month and that is corrupted, return from function
            if model_file_ii == (len(model_files)-1):
                print('------ All model files corrupted in %s. Skipping month.'%(yearmonth))
                return
            #else, continue to next file in month
            else:
                continue 

        #get indivudual dimension variable names  
        #no vertical dimension
        if len(mod_speci_obj.shape) == 3:
            x_varname = mod_speci_obj.dimensions[2]
            y_varname = mod_speci_obj.dimensions[1]
        #with vertical dimension
        elif len(mod_speci_obj.shape) == 4:
            x_varname = mod_speci_obj.dimensions[3]
            y_varname = mod_speci_obj.dimensions[2]
            z_varname = mod_speci_obj.dimensions[1]
            #check if vertical dimension goes up or down to get correct index for surface
            mod_vert_obj = mod_nc_root[z_varname]
            direction = mod_vert_obj.positive
            #if direction == 'up', surface index is 0 
            if direction == 'up':
                z_index = 0
            #if direction == 'down', surface index is -1
            elif direction == 'down':
                z_index = -1
            #if cannot determine a surface index, terminate process
            else: 
                print('Cannot determine surface index in vertical dimension. Terminating process.')
                sys.exit()

        #check if species grid dimensions are named correctly, and in correct BSC standard order, if not terminate process
        #this is done by checking the variable names of the x, y (and z if required) dimensions
        #X dimension is valid if 'lon' is contained within name, or is == 'x'
        if ('lon' not in x_varname) & (x_varname != 'x'):
           print('X dimension incorrectly named. Terminating process.')
           sys.exit()
        #Y dimension is valid if 'lat' is contained within name, or is == 'y'
        if ('lat' not in y_varname) & (y_varname != 'y'):
           print('Y dimension incorrectly named. Terminating process.')
           sys.exit()
        #Z dimension is valid if == 'z' or 'lev' or 'alt'
        if len(mod_speci_obj.shape) == 4:
            if (z_varname != 'lev') & (z_varname != 'z') & (z_varname != 'alt'):
                print('Z dimension incorrectly named. Terminating process.')
                sys.exit()

        #get instances of x/y grid dimension variables
        mod_lon_obj = mod_nc_root[x_varname]
        mod_lat_obj = mod_nc_root[y_varname]
        #get size of x/y grid dimensions
        x_N = mod_lon_obj.size
        y_N = mod_lat_obj.size
 
        #get name of longitude/latitude grid centre coordinate variables
        grid_centre_coordinates = mod_speci_obj.coordinates.split(' ')
        lon_centre_varname = grid_centre_coordinates[1] 
        lat_centre_varname = grid_centre_coordinates[0]
     
        #check if species grid centre coordinate are named correctly, and in correct BSC standard order, if not terminate process
        #longitude coordinate is valid if 'lon' is contained within name
        if ('lon' not in lon_centre_varname):
            print('Longitude grid centre coordinate incorrectly named. Terminating process.')
            sys.exit()
        #latitude coordinate is valid if 'lat' is contained within name
        if ('lat' not in lat_centre_varname):
            print('Latitude grid centre coordinate incorrectly named. Terminating process.')
            sys.exit()

        #get longitude and latitude grid centre values
        mod_lons_centre = np.float32(mod_nc_root[lon_centre_varname][:])
        mod_lats_centre = np.float32(mod_nc_root[lat_centre_varname][:])

        #break out of for loop, now that have read a valid model file in the month
        break

    #correct model grid longitude/latitude centre variables to be between -180:180 and -90:90 respectively
    mod_lons_centre[np.where(mod_lons_centre > 180.0)] = mod_lons_centre[np.where(mod_lons_centre > 180.0)] - 360.0
    mod_lats_centre[np.where(mod_lats_centre > 90.0)] = mod_lats_centre[np.where(mod_lats_centre > 90.0)] - 180.0

    #create grid domain edge polygon
    #this is handled differtly for regular grids (i.e. following lines of longitude/latitude), and non-regular grids (e.g. lambert-conformal) 

    #if grid type == 'crs, then is regular grid
    if mod_grid_type == 'crs':

        #get longitude/latitude grid resolution (taken from average of increment between longitude/latitude grid centres)
        lon_res = np.mean(np.diff(mod_lons_centre))
        lat_res = np.mean(np.diff(mod_lats_centre))
        #if either of the longitude/latitude increment resolutions are negative, this is because they run unconventionally from east to west/south to north
        #set indices for indexing accordingly
        if lon_res < 0:
            lon_left_ind = -1
            lon_right_ind = 0
        else:
            lon_left_ind = 0
            lon_right_ind = -1
        if lat_res < 0:
            lat_bottom_ind = -1
            lat_top_ind = 0
        else:
            lat_bottom_ind = 0
            lat_top_ind = -1

        #force longitude/latitude resolutions increments to be positive
        lon_res = np.abs(lon_res)
        lat_res = np.abs(lat_res)

        #get geographic coordinates of grid edges
        mod_lons_edge = mod_lons_centre - (lon_res/2.)
        if lon_right_ind == -1:
            mod_lons_edge = np.append(mod_lons_edge, mod_lons_edge[lon_right_ind]+lon_res)
        else:
            mod_lons_edge = np.insert(mod_lons_edge, 0, mod_lons_edge[lon_right_ind]+lon_res)
        mod_lats_edge = mod_lats_centre - (lat_res/2.)
        if lat_top_ind == -1:
            mod_lats_edge = np.append(mod_lats_edge, mod_lats_edge[lat_top_ind]+lat_res)
        else:
            mod_lats_edge = np.insert(mod_lats_edge, 0, mod_lats_edge[lat_top_ind]+lat_res)

        #get longitude/latitudes coordinates around grid edges
        left_edge_lon = np.full(len(mod_lats_edge), mod_lons_edge[lon_left_ind])
        right_edge_lon = np.full(len(mod_lats_edge), mod_lons_edge[lon_right_ind])
        if lon_left_ind == 0:
            top_edge_lon = mod_lons_edge[1:-1]
            bottom_edge_lon = mod_lons_edge[-2::-1]
        else:
            top_edge_lon = mod_lons_edge[-2::-1]
            bottom_edge_lon = mod_lons_edge[1:-1]
        lon_grid_edge = np.concatenate((left_edge_lon,top_edge_lon,right_edge_lon,bottom_edge_lon))

        bottom_edge_lat = np.full(len(mod_lons_edge)-1, mod_lats_edge[lat_bottom_ind])
        top_edge_lat = np.full(len(mod_lons_edge)-2, mod_lats_edge[lat_top_ind])
        
        if lat_bottom_ind == 0:
            left_edge_lat = mod_lats_edge[:]
            right_edge_lat = mod_lats_edge[::-1]
        else:
            left_edge_lat = mod_lats_edge[::-1]
            right_edge_lat = mod_lats_edge[:]
        lat_grid_edge = np.concatenate((left_edge_lat,top_edge_lat,right_edge_lat,bottom_edge_lat))

        #convert longitude/latitude bounding edge coordinate pairs to a shapely polygon
        model_grid_outline = np.vstack((lon_grid_edge,lat_grid_edge)).T
        model_grid_outline_poly = Polygon(model_grid_outline)

        #get 2D mesh of model longitude/latitude gridcell centres
        mod_lons_centre, mod_lats_centre = np.meshgrid(mod_lons_centre, mod_lats_centre)

    #here handle non-regular grids
    #currently only rotated pole and lambert conformal are supported
    elif mod_grid_type in ['rotated_pole', 'Lambert_conformal']:
  
        #get instance of variable which defines key variables associated with the non-regular grid
        non_regular_grid_type_obj = mod_nc_root[mod_grid_type]

        #get non-regular grid-specific variables used for defining coordinate reference system
        if mod_grid_type == 'rotated_pole':
            pole_longitude = np.float32(non_regular_grid_type_obj.grid_north_pole_longitude)
            pole_latitude = np.float32(non_regular_grid_type_obj.grid_north_pole_latitude)
        elif mod_grid_type == 'Lambert_conformal':
            standard_parallels = np.float32(non_regular_grid_type_obj.standard_parallel.split(','))
            central_longitude = np.float32(non_regular_grid_type_obj.longitude_of_central_meridian)
            central_latitude  = np.float32(non_regular_grid_type_obj.latitude_of_projection_origin)

        #read in non-regular grid x/y grid centre coordinates 
        x_centre = np.float32(mod_nc_root[x_varname][:])
        y_centre = np.float32(mod_nc_root[y_varname][:])

        #get x/y grid resolution (taken from average of increment between x/y grid centres)
        x_res = np.mean(np.diff(x_centre))
        y_res = np.mean(np.diff(y_centre))
        #if either of the x/y latitude increment resolutions are negative, this is because they run unconventionally from east to west/south to north
        #set indices for indexing accordingly
        if x_res < 0:
            x_left_ind = -1
            x_right_ind = 0
        else:
            x_left_ind = 0
            x_right_ind = -1
        if y_res < 0:
            y_bottom_ind = -1
            y_top_ind = 0
        else:
            y_bottom_ind = 0
            y_top_ind = -1

        #force x/y resolutions increments to be positive
        x_res = np.abs(x_res)
        y_res = np.abs(y_res)

        #get x/y coordinates of grid edges
        x_edge = x_centre - (x_res/2.)
        if x_right_ind == -1:
            x_edge = np.append(x_edge, x_edge[x_right_ind]+x_res)
        else:
            x_edge = np.insert(x_edge, 0, x_edge[x_right_ind]+x_res)
        y_edge = y_centre - (y_res/2.)
        if y_top_ind == -1:
            y_edge = np.append(y_edge, y_edge[y_top_ind]+y_res)
        else:
            y_edge = np.insert(y_edge, 0, y_edge[y_top_ind]+y_res) 

        #get x/y coordinates around non-regular grid edges
        left_edge_x = np.full(len(y_edge), x_edge[x_left_ind])
        right_edge_x = np.full(len(y_edge), x_edge[x_right_ind])
        if x_left_ind == 0:
            top_edge_x = x_edge[1:-1]
            bottom_edge_x = x_edge[-2::-1]
        else:
            top_edge_x = x_edge[-2::-1]
            bottom_edge_x = x_edge[1:-1]
        x_grid_edge = np.concatenate((left_edge_x,top_edge_x,right_edge_x,bottom_edge_x))

        bottom_edge_y = np.full(len(x_edge)-1, y_edge[y_bottom_ind])
        top_edge_y = np.full(len(x_edge)-2, y_edge[y_top_ind])
        if y_bottom_ind == 0:
            left_edge_y = y_edge[:]
            right_edge_y = y_edge[::-1]
        else:
            left_edge_y = y_edge[::-1]
            right_edge_y = y_edge[:]
        y_grid_edge = np.concatenate((left_edge_y,top_edge_y,right_edge_y,bottom_edge_y))

        #create cartopy coordinate reference system for the specific type non-standard grid, on WGS84 ellipsoid
        if mod_grid_type == 'rotated_pole':
            non_regular_grid_crs = ccrs.RotatedPole(pole_longitude=pole_longitude, pole_latitude=pole_latitude)
        elif mod_grid_type == 'Lambert_conformal':
            non_regular_grid_crs = ccrs.LambertConformal(central_longitude=central_longitude, central_latitude=central_latitude, standard_parallels=standard_parallels)
        #define a regular gridded coordinate reference system (Plate Carree), on WGS84 ellipsoid
        plate_carree = ccrs.PlateCarree()

        #convert bounding coordinates of box defined in non-regular grid coordinates to regular grid coordinates
        #then convert longitude/latitude bounding coordinate pairs to a shapely polygon 
        model_grid_outline = plate_carree.transform_points(non_regular_grid_crs, x_grid_edge, y_grid_edge)[:,:2]
        model_grid_outline_poly = Polygon(model_grid_outline)

    #the grid type cannot be handled, therefore terminate process
    else:
        print('Cannot handle grid of type: %s. Terminating process'%(mod_grid_type))
        sys.exit()

    #close model netcdf - now have all neccessary grid information
    mod_nc_root.close()

    #create monthly array for concatenating monthly model data together
    monthly_model_data = np.full((len(yearmonth_time), y_N, x_N), np.NaN, dtype=np.float32)

    #iterate and read daily chunked model files
    for model_ii, model_file in enumerate(model_files):

        #put model file read in try/except to catch corrup model files
        try:

            #read in daily chunked netcdf
            mod_nc_root = Dataset(model_file)

            #check if have time dimension in daily file, if do not, do not process file
            if 'time' not in list(mod_nc_root.dimensions.keys()):
                print('------ File %s is corrupt. Skipping.'%(model_file))
                continue 

            #get day of month from filename
            day = int(model_file.split('_%s'%(yearmonth))[-1][:2])
    
            #get difference from start of current day to start of month (in hours)
            diff_hours_start_day = (day-1)*24
       
            #get how number of hours from start of month is the start of the next day
            diff_hours_start_next_day = day*24

            #cross compare monthly time with daily file provided time to get indices on where to place data in monthly netCDF
            #if any of hours are >= diff_hours_start_next_day, then these will not be processed onto monthly netcdf, to ensure no overlapping data between days is processed
            adjusted_daily_file_time = diff_hours_start_day + mod_nc_root['time'][:]
            valid_daily_file_time_indices = np.where(adjusted_daily_file_time < diff_hours_start_next_day)[0]
       
            #read valid daily chunked data for valid indices
            #has no vertical dimension?
            if len(mod_nc_root[speci_to_process].shape) == 3:
                daily_chunked_data = mod_nc_root[speci_to_process][valid_daily_file_time_indices,:,:]
            #has vertical dimension 
            elif len(mod_nc_root[speci_to_process].shape) == 4:
                daily_chunked_data = mod_nc_root[speci_to_process][valid_daily_file_time_indices,z_index,:,:] 

            #set any values <0.0 to be NaN
            daily_chunked_data[daily_chunked_data<0.0] = np.NaN 

            #put daily chunked data in monthly chunked array 
            #doing data averaging where necessary
            if (temporal_resolution_to_output == 'hourly') or (temporal_resolution_to_output == 'monthly'): 
                adjusted_daily_file_time = adjusted_daily_file_time[valid_daily_file_time_indices]
                monthly_time_indices = np.arange(len(yearmonth_time))[np.in1d(yearmonth_time,adjusted_daily_file_time)]
                monthly_model_data[monthly_time_indices,:,:] = daily_chunked_data
 
            elif temporal_resolution_to_output == 'daily':
                day_ind = day-1
                monthly_model_data[day_ind,:,:] = np.nanmean(daily_chunked_data, axis=0)
        
            #close daily model netcdf
            mod_nc_root.close()

        except:
            print('------ File %s is corrupt. Skipping.'%(model_file))

    #for monthly resolution output case, now take average and set output time array
    if temporal_resolution_to_output == 'monthly':
        yearmonth_time = np.arange(1)
        monthly_model_data = np.nanmean(monthly_model_data, axis=0).reshape((1, y_N, x_N))

    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#
    #read observational data

    #read monthly observational file
    obs_nc_root = Dataset(obs_file)
    obs_lon_obj = obs_nc_root['longitude']
    obs_lons = obs_lon_obj[:]
    obs_lat_obj = obs_nc_root['latitude']
    obs_lats = obs_lat_obj[:]
    obs_station_reference_obj = obs_nc_root['station_reference']
    obs_station_references = obs_station_reference_obj[:]
    obs_measured_var_obj = obs_nc_root[speci_to_process]
    obs_speci_units = obs_measured_var_obj.units

    #get units conversion factor between model and observations (go from model to observational units)
    conv_obj = unit_converter.convert_units(mod_speci_units, obs_speci_units, 1.0)
    conv_obj.do_conversion()  
    conversion_factor = conv_obj.conversion_factor
    
    #get interpolation weights of model grid to observational stations (using inverse distance weighting interpolation) 
    n_neighbours_to_find = 4
    nearest_neighbour_inds, inverse_dists = n_nearest_neighbour_inverse_distance_weights(obs_lons, obs_lats, mod_lons_centre, mod_lats_centre, model_grid_outline_poly, n_neighbours=n_neighbours_to_find)

    #-------------------------------------------------------------------#
    #-------------------------------------------------------------------#

    #create observational interpolated monthly model netcdf
    
    #create new netCDF frame
    output_dir = '/gpfs/projects/bsc32/AC_cache/recon/ghost_interp/%s/%s/%s/%s/%s'%(experiment_to_process,grid_type_to_process,temporal_resolution_to_output,speci_to_process,GHOST_network_to_interpolate_against)
    #if output directory does not exist yet, create it
    try:
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir) 
    except OSError as err:
        pass
    root_grp = Dataset('%s/%s_%s.nc'%(output_dir,speci_to_process,yearmonth), 'w', format="NETCDF4") 
    #auto mask arrays
    root_grp.set_auto_mask(True)	

    #file contents
    root_grp.title = 'Inverse distance weighting (%s neighbours) interpolated %s experiment data for the component %s with reference to the measurement stations in the %s network in %s-%s.'%(n_neighbours_to_find, experiment_to_process, speci_to_process, GHOST_network_to_interpolate_against, year, month)
    root_grp.institution = 'Barcelona Supercomputing Center'
    root_grp.source = 'Experiment %s'%(experiment_to_process)
    root_grp.creator_name = 'Dene R. Bowdalo'
    root_grp.creator_email = 'dene.bowdalo@bsc.es'
    root_grp.conventions = 'CF-1.7'
    root_grp.data_version= '1.0'

    #netcdf dimensions
    root_grp.createDimension('time', len(yearmonth_time))
    root_grp.createDimension('station', len(obs_lons))
    root_grp.createDimension('model_centre_longitude', x_N)
    root_grp.createDimension('model_centre_latitude', y_N)
    root_grp.createDimension('grid_edge', model_grid_outline.shape[0])

    #create and write time variable
    time_var = root_grp.createVariable('time', 'u4', ('time'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['time'].size])
    time_var.standard_name = 'time'
    time_var.long_name = 'time'
    time_var.units = '%s since %s-%s-01 00:00:00'%(descriptive_temporal_resolution,year,month)
    time_var.description = 'Time in %s since %s-%s-01 00:00 UTC. Time given refers to the start of the time window the measurement is representative of (temporal resolution).'%(descriptive_temporal_resolution,year,month)
    time_var.axis = 'T'
    time_var.calendar = 'standard'
    time_var.tz = 'UTC'
    time_var[:] = yearmonth_time

    #create and write observational equivalent station reference variable
    station_reference_var = root_grp.createVariable('station_reference', str, ('station'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['station'].size])
    station_reference_var.standard_name = obs_station_reference_obj.standard_name
    station_reference_var.long_name = obs_station_reference_obj.long_name
    station_reference_var.units = obs_station_reference_obj.units
    station_reference_var.description = obs_station_reference_obj.description
    station_reference_var[:] = obs_station_references

    #create and write observational equivalent longitude/latitude variables
    longitude_var = root_grp.createVariable('longitude', 'f8', ('station'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['station'].size])
    longitude_var.standard_name = obs_lon_obj.standard_name
    longitude_var.long_name = obs_lon_obj.long_name
    longitude_var.units = obs_lon_obj.units
    longitude_var.description = obs_lon_obj.description    
    longitude_var.axis = 'X'
    longitude_var[:] = obs_lons

    latitude_var = root_grp.createVariable('latitude', 'f8', ('station'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['station'].size])
    latitude_var.standard_name = obs_lat_obj.standard_name
    latitude_var.long_name = obs_lat_obj.long_name
    latitude_var.units = obs_lat_obj.units
    latitude_var.description = obs_lat_obj.description        
    latitude_var.axis = 'Y'
    latitude_var[:] = obs_lats

    #create and write 2D meshed longitude/latitude gridcell centre variables
    model_centre_longitude_var = root_grp.createVariable('model_centre_longitude', 'f8', ('model_centre_latitude','model_centre_longitude'), zlib=True, complevel=4, chunksizes=[y_N,x_N])
    model_centre_longitude_var.standard_name = 'model centre longitude'
    model_centre_longitude_var.long_name = 'model centre longitude'
    model_centre_longitude_var.units = obs_lon_obj.units
    model_centre_longitude_var.description = '2D meshed grid centre longitudes with %s longitudes in %s bands of latitude'%(x_N,y_N)
    model_centre_longitude_var.axis = 'X'
    model_centre_longitude_var[:] = mod_lons_centre
    
    model_centre_latitude_var = root_grp.createVariable('model_centre_latitude', 'f8', ('model_centre_latitude','model_centre_longitude'), zlib=True, complevel=4, chunksizes=[y_N,x_N])
    model_centre_latitude_var.standard_name = 'model centre latitude'
    model_centre_latitude_var.long_name = 'model centre latitude'
    model_centre_latitude_var.units = obs_lat_obj.units
    model_centre_latitude_var.description = '2D meshed grid centre latitudes with %s latitudes in %s bands of longitude'%(y_N,x_N)
    model_centre_latitude_var.axis = 'Y'
    model_centre_latitude_var[:] = mod_lats_centre

    #create and write grid domain edge longitude/latitude variables
    grid_edge_longitude_var = root_grp.createVariable('grid_edge_longitude', 'f8', ('grid_edge'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['grid_edge'].size])   
    grid_edge_longitude_var.standard_name = 'grid edge longitude'
    grid_edge_longitude_var.long_name = 'grid edge longitude'
    grid_edge_longitude_var.units = obs_lon_obj.units
    grid_edge_longitude_var.description = 'Longitude coordinate along edge of grid domain (going clockwise around grid boundary from bottom-left corner).'
    grid_edge_longitude_var.description = 'X'
    grid_edge_longitude_var[:] = model_grid_outline[:,0]

    grid_edge_latitude_var = root_grp.createVariable('grid_edge_latitude', 'f8', ('grid_edge'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['grid_edge'].size])
    grid_edge_latitude_var.standard_name = 'grid edge latitude'
    grid_edge_latitude_var.long_name = 'grid edge latitude'
    grid_edge_latitude_var.units = obs_lat_obj.units
    grid_edge_latitude_var.description = 'Latitude coordinate along edge of grid domain (going clockwise around grid boundary from bottom-left corner).'
    grid_edge_latitude_var.description = 'Y'
    grid_edge_latitude_var[:] = model_grid_outline[:,1]

    #create measured variable
    measured_var = root_grp.createVariable(speci_to_process, 'f4', ('time','station'), zlib=True, complevel=4, chunksizes=[root_grp.dimensions['time'].size, root_grp.dimensions['station'].size])
    measured_var.standard_name = obs_measured_var_obj.standard_name
    measured_var.long_name = obs_measured_var_obj.long_name
    measured_var.units = obs_measured_var_obj.units
    measured_var.descrption = 'Interpolated value of %s from the experiment %s with reference to the measurement stations in the %s network'%(obs_measured_var_obj.standard_name, experiment_to_process, GHOST_network_to_interpolate_against)

    #iterate through observational stations
    #use calculated interpolated weights per station to model grid to produce model reciprocal output
    for station_ii, (obs_lon, obs_lat) in enumerate(zip(obs_lons, obs_lats)):
        station_weights = inverse_dists[station_ii,:]
        #if all station weights are 0, station is outside model grid domain
        #set all values to be NaN
        if np.all(station_weights == 0):
            interp_vals = np.full(len(yearmonth_time), np.NaN, dtype=np.float32)
        else:
            #get reciprocal model data at N nearest neighbours to observational station 
            cut_model_data = monthly_model_data[:,nearest_neighbour_inds[station_ii,:n_neighbours_to_find],nearest_neighbour_inds[station_ii,n_neighbours_to_find:]]
            #create mask where data == NaN
            invalid_mask = np.isnan(cut_model_data)
            #create masked array
            cut_model_data = np.ma.MaskedArray(cut_model_data, mask=invalid_mask)
            #interpolate masked array across time dimension using interpolated weights per station
            interp_vals = np.ma.average(cut_model_data, weights = station_weights, axis=1)

        #write measured variable (multiplying with conversion factor in the process)
        measured_var[:,station_ii] = interp_vals * conversion_factor

    #close completed netCDF
    root_grp.close() 

    #close observational netcdf
    obs_nc_root.close()

    return

###--------------------------------------------------------------------------------------------------###
#define inverse distance weighting interpolation function
###--------------------------------------------------------------------------------------------------###

def n_nearest_neighbour_inverse_distance_weights(obs_lons,obs_lats,mod_lons_centre,mod_lats_centre,model_grid_outline_poly,n_neighbours=1):
   
    '''function that calculates N nearest neighbour inverse distance weights (and indices) of model gridcells centres from an array of geographic observational station coordinates.
       both observational and model geographic longitude/latitude coordinates are first converted to cartesian ECEF (Earth Centred, Earth Fixed) coordinates, before calculating distances.
       weights returned for obervational stations not contained within model grid extents are all zero.
    '''

    #for each pair of observational station geographic coordinates, test if station is inside model grid
    obs_inside_model_grid = np.array([model_grid_outline_poly.contains(Point(obs_lon,obs_lat)) for obs_lon,obs_lat in zip(obs_lons,obs_lats)]) 

    #flatten model centre lon/lat arrays
    #generate equal length altitude arrays
    obs_alts = np.zeros(len(obs_lons))
    flat_mod_lons_centre = mod_lons_centre.ravel()
    flat_mod_lats_centre = mod_lats_centre.ravel()
    flat_mod_alts = np.zeros(len(flat_mod_lats_centre))

    #convert observational/model geographic longitude/latitude coordinates to cartesian ECEF (Earth Centred, Earth Fixed) coordinates, assuming WGS84 datum and ellipsoid, and that all heights = 0
    #ECEF coordiantes represent positions (in meters) as X, Y, Z coordinates, approximating the earth surface as an ellipsoid of revolution.
    #This conversion is for the subsequent calculation of euclidean distances of the model gridcell centres from each observational station 
    #Defining the distance between two points on the earth's surface as simply the euclidean distance between the two lat/lon pairs could lead to inaccurate results depending on the distance between two points (i.e. 1 deg. of longitude varies with latitude)  
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    obs_x, obs_y, obs_z = pyproj.transform(lla, ecef, obs_lons, obs_lats, obs_alts, radians=False)
    mod_x, mod_y, mod_z = pyproj.transform(lla, ecef, flat_mod_lons_centre, flat_mod_lats_centre, flat_mod_alts, radians=False)

    #stack converted cartesian coordinates for preparation of calculation of nearest neighbour distances using Scipy cKDTree
    obs_lonlat = np.column_stack((obs_x,obs_y))
    mod_lonlat = np.column_stack((mod_x,mod_y))
    
    #generate KDtree
    tree = spatial.cKDTree(mod_lonlat)
    
    #get n-neighbour nearest distances/indices (ravel form) of model gridcell centres from each observational station  
    dists,idx = tree.query(obs_lonlat,k=n_neighbours)
    nearest_neighbour_inds = np.column_stack(np.unravel_index(idx,mod_lons_centre.shape))
      
    #take the reciprocals of the nearest neighbours distances from the observational points
    inverse_dists = np.reciprocal(dists)

    #set reciprocal distances for all observational points outside model grid extent to be 0
    inverse_dists[~obs_inside_model_grid,:] = 0.0

    return nearest_neighbour_inds, inverse_dists

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#get command line arguments variables passed upon execution
experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]

#read defined experiments dictionary
from defined_experiments import defined_experiments_dictionary

#get experiment specific directory (take P9 experiment directory preferentially over esarchive directory)
exp_dict = defined_experiments_dictionary[experiment_to_process]
if 'P9' in list(exp_dict.keys()):
   exp_dir = exp_dict['P9']
else:
   exp_dir = exp_dict['esarchive']

#pass variables to interpolation function and do interpolation 
process_monthly_interpolated_netCDF(experiment_to_process, exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth)

print('DONE')
