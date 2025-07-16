""" Module storing writing functions """

import copy
from datetime import datetime, timedelta
import os
import sys
import yaml

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray as xr

from providentia.auxiliar import CURRENT_PATH, join
from .configuration import write_conf

# get current path and providentia root path
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

# define possible temporal resolutions
possible_resolutions = ['hourly', 'hourly_instantaneous', '3hourly', '3hourly_instantaneous', 
                        '6hourly', '6hourly_instantaneous', 'daily', 'monthly', 'annual']

def export_data_npz(prv, fname, input_dialogue=False, set_in_memory=False):
    """ Function that writes out current data / ghost data / metadata 
        in memory to .npy file. 

        :prv: Instance of providentia 
        :type prv: instance of Dashboard / Library
        :fname: Name of the file to save
        :type fname: str
        :input_dialogue: boolean informing whether to open input prompt to ask if to filter data or not
        :type input_dialogue: boolean
        :set_in_memory: booolean informing if saved data is set in memory
        :type set_in_memory: boolean
    """
    # import function from dashboard_elements
    from .dashboard_elements import InputDialog

    # ensure fname has correct extension
    if fname[-4:] != '.npz':
        fname = '{}.npz'.format(fname)

    # open dialog to choose if data is filtered or not
    if input_dialogue:
        title = 'Export data'
        msg = 'Select option'
        options = ['Apply metadata filters, data filters, temporal colocation, calibration factor, and resampling to exported data', 
                   'Do not apply metadata filters, data filters, temporal colocation, calibration factor, and resampling to exported data']
        dialog = InputDialog(prv, title, msg, options)
        selected_option, okpressed = dialog.selected_option, dialog.okpressed
        if okpressed is False:
            return False
        if selected_option == options[0]:
            apply_filters = True
        elif selected_option == options[1]:
            apply_filters = False
    else:
        apply_filters = True

    # create dict to save data
    save_data_dict = {}

    # save data / ghost data / metadata
    for networkspeci in prv.networkspecies:

        # get species
        speci = networkspeci.split('|')[1]

        # get data array
        # if apply_filters is active get the filtered data array
        # also temporally resample data if required
        if apply_filters:
            valid_station_inds = prv.valid_station_inds_temporal_colocation[networkspeci][prv.observations_data_label]
            data_array = copy.deepcopy(prv.data_in_memory_filtered[networkspeci])
            
            # apply NaNs for temporal colocation
            data_array[:, prv.temporal_colocation_nans[networkspeci]] = np.NaN

            # cut data array for valid station inds
            data_array = np.take(data_array, valid_station_inds, axis=1)

        # otherwise, take unfiltered array
        else:
            valid_station_inds = np.arange(prv.data_in_memory[networkspeci].shape[1], dtype=np.int32)
            data_array = copy.deepcopy(prv.data_in_memory[networkspeci])

        # save data / metadata
        save_data_dict['{}_data'.format(networkspeci)] = data_array
        save_data_dict['{}_metadata'.format(networkspeci)] = np.take(prv.metadata_in_memory[networkspeci], 
                       valid_station_inds, axis=0)
        # save GHOST specific variables
        if prv.reading_ghost:
            save_data_dict['{}_ghost_data'.format(networkspeci)] = np.take(prv.ghost_data_in_memory[networkspeci], 
                           valid_station_inds, axis=1)
            save_data_dict['{}_qa'.format(networkspeci)] = np.unique(prv.qa_per_species[speci])
            save_data_dict['{}_flags'.format(networkspeci)] = np.unique(prv.flags)

    # save out key variables 
    save_data_dict['time'] = prv.time_array
    save_data_dict['resolution'] = prv.resolution
    save_data_dict['data_labels'] = prv.data_labels
    save_data_dict['start_date'] = prv.start_date
    save_data_dict['end_date'] = prv.end_date
    save_data_dict['temporal_colocation'] = prv.temporal_colocation
    save_data_dict['spatial_colocation'] = prv.spatial_colocation
    save_data_dict['filter_species'] = prv.filter_species
    save_data_dict['forecast'] = prv.forecast
    if prv.reading_ghost:
        save_data_dict['ghost_version'] = prv.ghost_version
        save_data_dict['ghost_data_variables'] = prv.ghost_data_vars_to_read

    # save resampled time and resolution
    if (prv.resampling_resolution in possible_resolutions) & (apply_filters):
        save_data_dict['time_resampled'] = prv.time_index
        save_data_dict['resampling_resolution'] = prv.resampling_resolution

    # save out statistical information
    save_data_dict['statistic_mode'] = prv.statistic_mode
    save_data_dict['statistic_aggregation'] = prv.statistic_aggregation
    save_data_dict['periodic_statistic_mode'] = prv.periodic_statistic_mode
    save_data_dict['periodic_statistic_aggregation'] = prv.periodic_statistic_aggregation
    save_data_dict['timeseries_statistic_aggregation'] = prv.timeseries_statistic_aggregation

    # save out miscellaneous variables
    calibration_factor = ''
    for factor_ii, (exp, factor) in enumerate(prv.calibration_factor.items()):
        if factor_ii == (len(prv.calibration_factor) - 1):
            calibration_factor += '{} ({})'.format(exp, factor)
        else:
            calibration_factor += '{} ({}), '.format(exp, factor)
    save_data_dict['calibration_factor'] = calibration_factor
    save_data_dict['remove_extreme_stations'] = prv.remove_extreme_stations
    save_data_dict['spatial_colocation_tolerance'] = prv.spatial_colocation_tolerance
    save_data_dict['spatial_colocation_station_reference'] = prv.spatial_colocation_station_reference
    save_data_dict['spatial_colocation_station_name'] = prv.spatial_colocation_station_name
    save_data_dict['spatial_colocation_longitude_latitude'] = prv.spatial_colocation_longitude_latitude
    save_data_dict['spatial_colocation_measurement_altitude'] = prv.spatial_colocation_measurement_altitude
    save_data_dict['spatial_colocation_validation'] = prv.spatial_colocation_validation
    save_data_dict['spatial_colocation_validation_tolerance'] = prv.spatial_colocation_validation_tolerance

    # save out dict to .npz file
    np.savez(fname, **save_data_dict)

    # if set_in_memory is active, load and return the variable in memory  
    if set_in_memory:   
        data = np.load(fname, allow_pickle=True)
        # delete temporary save file after load
        os.remove(fname)  
        return data        

    return True          

def export_netcdf(prv, fname, input_dialogue=False, set_in_memory=False, xarray=False):
    """ Write data and metadata to netcdf file. 
    
        :prv: Instance of providentia
        :type prv: instance of Dashboard / Library
        :fname: Name of the file to save
        :type fname: str
        :input_dialogue: boolean informing whether to open input prompt to ask if to filter data or not
        :type input_dialogue: boolean
        :set_in_memory: booolean informing if saved data is set in memory
        :type set_in_memory: boolean
        :xarray: booolean informing if data is read in xarray format or not
        :type xarray: boolean

    """
    # import function from dashboard_elements
    from .dashboard_elements import InputDialog

    # ensure fname has correct extension
    if fname[-3:] != '.nc':
        fname = '{}.nc'.format(fname)

    # open dialog to choose if data is filtered or not
    if input_dialogue:
        title = 'Export data'
        msg = 'Select option'
        options = ['Apply metadata filters, data filters, temporal colocation, calibration factor, and resampling to exported data', 
                   'Do not apply metadata filters, data filters, temporal colocation, calibration factor, and resampling to exported data']
        dialog = InputDialog(prv, title, msg, options)
        selected_option, okpressed = dialog.selected_option, dialog.okpressed
        if okpressed is False:
            return False
        if selected_option == options[0]:
            apply_filters = True
        elif selected_option == options[1]:
            apply_filters = False
    else:
        apply_filters = True

    # set up some structural variables
    from GHOST_standards import (standard_parameters, get_standard_data, get_standard_metadata,
                                 standard_QA_name_to_QA_code, standard_data_flag_name_to_data_flag_code)
    parameter_dictionary = {}
    for _, param_dict in standard_parameters.items():
        parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict

    # dictionary to map python types to netcdf types
    type_map = {np.uint8: 'u1', np.uint32: 'u4', object: str,
                np.float32: 'f4', np.float64: 'f8'}

    # start file
    fout = Dataset(fname, 'w', format="NETCDF4")

    # file contents
    fout.title = 'Saved data from Providentia'
    fout.institution = 'Barcelona Supercomputing Center'
    fout.source = 'Providentia'
    if prv.reading_ghost:
        fout.data_version = prv.ghost_version

    # netcdf dimensions
    fout.createDimension('data_label', len(prv.data_labels))
    fout.createDimension('time', len(prv.time_array))
    fout.createDimension('month', len(prv.yearmonths))
    # create resampling resolution if needed
    if (prv.resampling_resolution in possible_resolutions) & (apply_filters):
        fout.createDimension('time_resampled', len(prv.time_index))
    # create GHOST variables
    if prv.reading_ghost:
        fout.createDimension('qa', len(list(standard_QA_name_to_QA_code.keys())))
        fout.createDimension('flag', len(list(standard_data_flag_name_to_data_flag_code.keys())))
        fout.createDimension('ghost_data_variable', len(prv.ghost_data_vars_to_read))

    # iterate through networkspecies 
    for speci_ii, networkspeci in enumerate(prv.networkspecies):

        # get species
        speci = networkspeci.split('|')[1]

        # get prefix (name of networkspeci) to be added to variable names
        if prv.reading_ghost:
            var_prefix = networkspeci
        else:
            var_prefix = networkspeci.replace('/', '_')

        # get some key variables for speci
        parameter_details = parameter_dictionary[speci]
        metadata_format_dict = get_standard_metadata(parameter_details)
        data_format_dict = get_standard_data(parameter_details)
       
        # get data array
        # if apply_filters is active get the filtered data array
        if apply_filters:
            valid_station_inds = prv.valid_station_inds_temporal_colocation[networkspeci][prv.observations_data_label]
            data_array = copy.deepcopy(prv.data_in_memory_filtered[networkspeci])

            # apply NaNs for temporal colocation
            data_array[:, prv.temporal_colocation_nans[networkspeci]] = np.NaN

            # cut data array for valid station inds
            data_array = np.take(data_array, valid_station_inds, axis=1)

        # otherwise, take unfiltered array
        else:
            valid_station_inds = np.arange(prv.data_in_memory[networkspeci].shape[1], dtype=np.int32)
            data_array = copy.deepcopy(prv.data_in_memory[networkspeci])

        # set variables independent of network / speci on first pass
        if speci_ii == 0:

            # time
            current_data_type = type_map[data_format_dict['time']['data_type']]
            var = fout.createVariable('time', current_data_type, ('time',))
            # set attributes
            if 'hourly' in prv.resolution:
                res_str = 'hours'
            elif 'daily' in prv.resolution:
                res_str = 'days'
            elif 'monthly' in prv.resolution:
                res_str = 'months'
            var.standard_name = data_format_dict['time']['standard_name']
            var.long_name = data_format_dict['time']['long_name']
            var.units = '{} since {}-{}-01 00:00:00'.format(res_str, 
                                                            str(prv.start_date)[:4], 
                                                            str(prv.start_date)[4:6])
            msg = 'Time in {} since {}-{}-01 00:00 UTC. Time given refers '.format(res_str, 
                                                                                   str(prv.start_date)[:4], 
                                                                                   str(prv.start_date)[4:6])
            msg += 'to the start of the time window the measurement is representative of '
            msg += '(temporal resolution).'
            var.description = msg
            var.axis = 'T'
            var.calendar = 'standard'
            var.tz = 'UTC'
            if prv.resolution in ['3hourly', '6hourly']:
                # get indices of time_array in an array with all hours from start date to end date
                all_hours_array = pd.to_datetime(np.arange(datetime(int(str(prv.start_date)[0:4]), 
                                                                    int(str(prv.start_date)[4:6]), 
                                                                    int(str(prv.start_date)[6:8])), 
                                                           datetime(int(str(prv.end_date)[0:4]), 
                                                                    int(str(prv.end_date)[4:6]), 
                                                                    int(str(prv.end_date)[6:8])), 
                                                           timedelta(hours=1)))
                time_var = np.where(np.in1d(all_hours_array, prv.time_array))[0]
            else:
                time_var = np.arange(len(prv.time_array))
            var[:] = time_var

            # time resampled - create variable if have resampled data, and apply_filters is active
            if (prv.resampling_resolution in possible_resolutions) & (apply_filters):

                current_data_type = type_map[data_format_dict['time']['data_type']]
                var = fout.createVariable('time_resampled', current_data_type, ('time_resampled',))
                # set attributes
                if 'hourly' in prv.resampling_resolution:
                    res_str = 'hours'
                elif 'daily' in prv.resampling_resolution:
                    res_str = 'days'
                elif 'monthly' in prv.resampling_resolution:
                    res_str = 'months'
                var.standard_name = data_format_dict['time']['standard_name']
                var.long_name = data_format_dict['time']['long_name']
                var.units = '{} since {}-{}-01 00:00:00'.format(res_str, 
                                                                str(prv.start_date)[:4], 
                                                                str(prv.start_date)[4:6])
                msg = 'Time in {} since {}-{}-01 00:00 UTC. Time given refers '.format(res_str, 
                                                                                    str(prv.start_date)[:4], 
                                                                                    str(prv.start_date)[4:6])
                msg += 'to the start of the time window the measurement is representative of '
                msg += '(temporal resolution).'
                var.description = msg
                var.axis = 'T'
                var.calendar = 'standard'
                var.tz = 'UTC'
                if prv.resampling_resolution in ['3hourly', '6hourly']:
                    # get indices of time_index in an array with all hours from start date to end date
                    all_hours_array = pd.to_datetime(np.arange(datetime(int(str(prv.start_date)[0:4]), 
                                                                        int(str(prv.start_date)[4:6]), 
                                                                        int(str(prv.start_date)[6:8])), 
                                                            datetime(int(str(prv.end_date)[0:4]), 
                                                                     int(str(prv.end_date)[4:6]), 
                                                                     int(str(prv.end_date)[6:8])), 
                                                            timedelta(hours=1)))
                    time_var_resampled = np.where(np.in1d(all_hours_array, prv.time_index))[0]
                else:
                    time_var_resampled = np.arange(len(prv.time_index))
                var[:] = time_var_resampled
                
            # miscellaneous variables 
            var = fout.createVariable('data_labels', str, ('data_label',))
            var.standard_name = 'data_labels'
            var.long_name = 'data_labels'
            var.description = 'Labels associated with each data array, e.g. observations, experiment_1, etc.'
            var[:] = np.array(prv.data_labels)
            # GHOST variables
            if prv.reading_ghost:
                var = fout.createVariable('ghost_data_variables', str, ('ghost_data_variable',))
                var.standard_name = 'ghost_data_variables'
                var.long_name = 'ghost_data_variables'
                var.description = 'The names of the GHOST data variables used for additional filtering.'
                var[:] = np.array(prv.ghost_data_vars_to_read)
                 
        # set networkspeci station dimension
        station_dimension_var = 'station_{}'.format(var_prefix)
        fout.createDimension(station_dimension_var, len(valid_station_inds))

        # set data variable
        current_data_type = type_map[data_format_dict[speci]['data_type']]
        # set dimension to be time_resampled if are resampling and apply_filters is active
        if (prv.resampling_resolution in possible_resolutions) & (apply_filters):
            var = fout.createVariable('{}_data'.format(var_prefix), current_data_type, 
                                     ('data_label', station_dimension_var, 'time_resampled',))
        else:
            var = fout.createVariable('{}_data'.format(var_prefix), current_data_type, 
                                     ('data_label', station_dimension_var, 'time',))

        # set attributes
        var.standard_name = data_format_dict[speci]['standard_name']
        var.long_name = data_format_dict[speci]['long_name']
        var.units = data_format_dict[speci]['units']
        var.description = data_format_dict[speci]['description']
        var.resolution = str(prv.resolution)
        var.resampling_resolution = str(prv.resampling_resolution)
        var.start_date = str(prv.start_date)
        var.end_date = str(prv.end_date)
        var.temporal_colocation = str(prv.temporal_colocation)
        var.spatial_colocation = str(prv.spatial_colocation)
        var.filter_species = str(prv.filter_species)
        var.forecast = str(prv.forecast)
        var.calibration_factor = str(prv.calibration_factor)
        var.remove_extreme_stations = str(prv.remove_extreme_stations)
        var.statistic_mode = str(prv.statistic_mode)
        var.statistic_aggregation = str(prv.statistic_aggregation)
        var.periodic_statistic_mode = str(prv.periodic_statistic_mode)
        var.periodic_statistic_aggregation = str(prv.periodic_statistic_aggregation)
        var.timeseries_statistic_aggregation = str(prv.timeseries_statistic_aggregation)
        var.spatial_colocation_tolerance = str(prv.spatial_colocation_tolerance)
        var.spatial_colocation_station_reference = str(prv.spatial_colocation_station_reference)
        var.spatial_colocation_station_name = str(prv.spatial_colocation_station_name)
        var.spatial_colocation_longitude_latitude = str(prv.spatial_colocation_longitude_latitude)
        var.spatial_colocation_measurement_altitude = str(prv.spatial_colocation_measurement_altitude)
        var.spatial_colocation_validation = str(prv.spatial_colocation_validation)
        var.spatial_colocation_validation_tolerance = str(prv.spatial_colocation_validation_tolerance)
        if prv.reading_ghost:
            var.ghost_version = str(prv.ghost_version)
        var[:] = data_array

        # GHOST data
        if prv.reading_ghost:

            # set GHOST data variable (e.g. representativity)
            var = fout.createVariable('{}_ghost_data'.format(networkspeci), 'f4', 
                                     ('ghost_data_variable', station_dimension_var, 'time',))
            # set attributes and data
            var.standard_name = '{}_ghost_data'.format(networkspeci) 
            var.long_name = '{}_ghost_data'.format(networkspeci)
            var.description = 'GHOST data variables used for additional filtering.'
            var[:] = np.take(prv.ghost_data_in_memory[networkspeci], valid_station_inds, axis=1)

            # set qa variable
            var = fout.createVariable('{}_qa'.format(networkspeci), 'u1', ('qa',), fill_value=255)
            # set attributes and data
            var.standard_name = '{}_qa'.format(networkspeci)
            var.long_name = '{}_qa'.format(networkspeci)
            var.description = 'GHOST QA flag codes applied to filter data.'
            unique_qa = list(np.unique(prv.qa_per_species[speci]))
            padded_unique_qa = np.array([unique_qa + [255]*(fout.dimensions['qa'].size - len(unique_qa))], dtype=np.uint8)
            var[:] = padded_unique_qa

            # set flags variable
            var = fout.createVariable('{}_flags'.format(networkspeci), 'u1', ('flag',), fill_value=255)
            # set attributes and data
            var.standard_name = '{}_flags'.format(networkspeci)
            var.long_name = '{}_flags'.format(networkspeci)
            var.description = 'GHOST standardised data reporter flag codes applied to filter data.'
            unique_flag = list(np.unique(prv.flags))
            padded_unique_flag = np.array([unique_flag + [255]*(fout.dimensions['flag'].size - len(unique_flag))], dtype=np.uint8)
            var[:] = padded_unique_flag

        # save metadata (as individual variables)
        metadata_arr = np.take(prv.metadata_in_memory[networkspeci], valid_station_inds, axis=0)
        for metadata_var in metadata_arr.dtype.names:
            
            current_data_type = type_map[metadata_format_dict[metadata_var]['data_type']] 
            var = fout.createVariable('{}_{}'.format(var_prefix, metadata_var), 
                                      current_data_type, (station_dimension_var, 'month',))

            # set attributes
            var.standard_name = metadata_format_dict[metadata_var]['standard_name']
            var.long_name = metadata_format_dict[metadata_var]['long_name']
            var.units = metadata_format_dict[metadata_var]['units']
            var.description = metadata_format_dict[metadata_var]['description']
            if metadata_var == 'longitude':
                var.axis = 'X'
            elif metadata_var == 'latitude':
                var.axis = 'Y'
            if current_data_type == str:
                var[:] = metadata_arr[metadata_var].astype(str)
            else:
                var[:] = metadata_arr[metadata_var]

    # close writing to netCDF
    fout.close()

    # if set_in_memory is active, load and return the variable in memory  
    if set_in_memory:   
        if xarray:
            data = xr.load_dataset(fname)
        else:
            data = Dataset(fname)

        # delete temporary save file after load
        os.remove(fname)  

        return data
    
    return True

def export_configuration(prv, cname, separator="||"):
    """ Create all items to be written in configuration file
        and send them to write_conf.

        :prv: Instance of providentia
        :type prv: instance of Dashboard / Library
        :cname: Name for the configuration file
        :type cname: str
        :separator: delimiter for keep/remove fields
        :type separator: str
    """
    
    # if no data was loaded, there won't be any maximum nor minimum value
    if prv.le_minimum_value.text() == '' and  prv.le_minimum_value.text() == '':
        raise Exception("Error: No data available for writing. Please click on READ before trying to save any file.")
    
    # load initialisation defaults
    init_defaults = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'init_prov_dev.yaml')))
    # load variable defaults
    var_defaults = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'prov_defaults.yaml')))
    # load modifiable variable defaults
    modifiable_var_defaults = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'init_prov.yaml')))
    # merge defaults
    merged_defaults = init_defaults.copy()
    merged_defaults.update(var_defaults)
    merged_defaults.update(modifiable_var_defaults)

    # ensure cname has correct extension
    if cname[-5:] != '.conf':
        cname = '{}.conf'.format(cname)

    # set section and subsection names in config file
    if not hasattr(prv, 'section'):
        section = 'SECTION1'
        #subsection = '[SUBSECTION1]'
        subsection = None
    else:
        if '·' in prv.section:
            section = prv.section.split('·')[0]
            #subsection = '[' + prv.section.split('·')[1] + ']'
            subsection = None
        else:
            section = prv.section
            subsection = None

    options = {}
    options['section'] = {}
    options['subsection'] = {}

    # default
    options['section'] = {'network': prv.network[0],
                          'species': prv.species[0],
                          'resolution': prv.resolution,
                          'start_date': prv.start_date,
                          'end_date': prv.end_date}

    # experiments
    if prv.experiments_menu['checkboxes']['keep_selected']:
        options['section']['experiments'] = ",".join(str(i) for i in prv.experiments_menu['checkboxes']['keep_selected'])

    # colocation variables
    options['section'].update({'temporal_colocation': prv.temporal_colocation,
                               'spatial_colocation': prv.spatial_colocation})

    # filter species
    if len(prv.filter_species) > 0:
        filter_species = str(copy.deepcopy(prv.filter_species))
        filter_species = filter_species.replace("[", "(").replace("]", ")")
        filter_species = filter_species.replace("{", "").replace("}", "")
        filter_species = filter_species.replace("'", "")
        filter_species = filter_species.replace(":", "")
        filter_species = filter_species.replace("|", ":")
        filter_species = filter_species.replace("((", "(")
        filter_species = filter_species.replace("))", ")")
        options['section'].update({'filter_species': filter_species})

    # resampling_resolution
    if prv.resampling_resolution in possible_resolutions:
        options['section'].update({'resampling_resolution': prv.resampling_resolution})

    # statistic_mode
    if prv.statistic_mode != merged_defaults['statistic_mode']:
        options['section'].update({'statistic_mode': prv.statistic_mode})

    # statistic_aggregation 
    if prv.statistic_aggregation != merged_defaults['statistic_aggregation'][prv.statistic_mode]:
        options['section'].update({'statistic_aggregation': prv.statistic_aggregation})

    # calibration_factor
    if len(prv.calibration_factor) > 0:
        calibration_factor = ''
        for factor_ii, (exp, factor) in enumerate(prv.calibration_factor.items()):
            if factor_ii == (len(prv.calibration_factor) - 1):
                calibration_factor += '{} ({})'.format(exp, factor)
            else:
                calibration_factor += '{} ({}), '.format(exp, factor)
        options['section'].update({'calibration_factor': calibration_factor})

    # forecast
    if len(prv.forecast) > 0:
        options['section']['forecast'] = ",".join(str(i) for i in prv.forecast)

    # qa
    if set(prv.qa_menu['checkboxes']['remove_selected']) != set(prv.qa_menu['checkboxes']['remove_default']):
        options['section']['qa'] = ",".join(str(i) for i in prv.qa_menu['checkboxes']['remove_selected'])
        
    # flags
    if prv.flag_menu['checkboxes']['remove_selected']:
        options['section']['flags'] = ",".join(str(i) for i in prv.flag_menu['checkboxes']['remove_selected'])

    # representativity
    for i, label in enumerate(prv.representativity_menu['rangeboxes']['labels']):
        if 'max_gap' in label:
            if prv.representativity_menu['rangeboxes']['current_lower'][i] != '100':
                options['section'][label] = prv.representativity_menu['rangeboxes']['current_lower'][i]
        else:
            if prv.representativity_menu['rangeboxes']['current_lower'][i] != '0':
                options['section'][label] = prv.representativity_menu['rangeboxes']['current_lower'][i]

    # period
    if prv.period_menu['checkboxes']['keep_selected'] or prv.period_menu['checkboxes']['remove_selected']:
        period_k = "keep: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['keep_selected']) + separator
        period_r = " remove: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['remove_selected']) + separator
        options['section']['period'] = period_k + period_r
    
    # bounds
    if np.float32(prv.le_minimum_value.text()) != \
            np.float32(prv.parameter_dictionary[prv.species[0]]['extreme_lower_limit']):
        options['section']['lower_bound'] = prv.le_minimum_value.text()
    if np.float32(prv.le_maximum_value.text()) != \
            np.float32(prv.parameter_dictionary[prv.species[0]]['extreme_upper_limit']):
        options['section']['upper_bound'] = prv.le_maximum_value.text()

    # metadata
    for menu_type in prv.metadata_types:
        # treat ranges first
        for i, label in enumerate(prv.metadata_menu[menu_type]['rangeboxes']['labels']):
            lower_cur = prv.metadata_menu[menu_type]['rangeboxes']['current_lower'][i]
            lower_def = prv.metadata_menu[menu_type]['rangeboxes']['lower_default'][i]
            upper_cur = prv.metadata_menu[menu_type]['rangeboxes']['current_upper'][i]
            upper_def = prv.metadata_menu[menu_type]['rangeboxes']['upper_default'][i]
            # do not write nans
            if (pd.isnull(lower_cur)) or (pd.isnull(upper_cur)) or (lower_cur == 'nan') or (upper_cur == 'nan'):
                continue
            # write field if different from default values
            elif (lower_cur != lower_def) or (upper_cur != upper_def):
                options['section'][label] = lower_cur + ", " + upper_cur

        # and then treat the keep/remove
        for label in prv.metadata_menu[menu_type]['navigation_buttons']['labels']:
            keeps = prv.metadata_menu[menu_type][label]['checkboxes']['keep_selected']
            removes = prv.metadata_menu[menu_type][label]['checkboxes']['remove_selected']
            if keeps or removes:
                meta_keep = "keep: " + ",".join(str(i) for i in keeps) + separator
                meta_remove = " remove: " + ",".join(str(i) for i in removes) + separator
                options['section'][label] = meta_keep + meta_remove

    # map extent
    if prv.map_extent != [-180.0,180.0,-90.0,90.0]:
        options['section']['map_extent'] = ",".join(str(i) for i in prv.map_extent)

    # active dashboard plots 
    if prv.active_dashboard_plots != merged_defaults['active_dashboard_plots']:
        options['section']['active_dashboard_plots'] = ",".join(str(i) for i in prv.active_dashboard_plots)

    # miscellaneous fields that need string joining 
    #options['section'].update({'map_extent': ",".join(str(i) for i in prv.map_extent),
    #                           'active_dashboard_plots': ",".join(str(i) for i in prv.active_dashboard_plots)})

    # plot_characteristics_filename
    if ((prv.plot_characteristics_filename != join(PROVIDENTIA_ROOT, 'settings/plot_characteristics.yaml')) &
       (prv.plot_characteristics_filename != '')):
        options['section'].update({'plot_characteristics_filename': prv.plot_characteristics_filename})

    # add each variable to be writen out if: not a variable that should not be written, if value not default value, 
    # not already written out, not None or empty str, and not empty list/dict
    for var, default_val in merged_defaults.items():
        default_val = str(copy.deepcopy(default_val))
        current_val = str(copy.deepcopy(getattr(prv, var)))
        if ((var not in merged_defaults['non_writing_vars']) & (current_val != default_val) & 
           (var not in options['section']) & 
           ((current_val != "None") & (current_val != "") & (current_val != "[]") & (current_val != "{}"))):
            options['section'].update({var: current_val})
    
    # write .conf file
    write_conf(section, subsection, cname, options)

    return True