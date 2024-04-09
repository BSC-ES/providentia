""" Module storing writing functions """

import copy
from netCDF4 import Dataset
import numpy as np
import os
import pandas as pd
import sys
import xarray as xr

from .configuration import write_conf
from .dashboard_elements import InputDialog


def export_data_npz(prv, fname, input_dialogue=False, set_in_memory=False):
    """ Function that writes out current data / ghost data / metadata 
        in memory to .npy file. 

        :prv: Instance of providentia 
        :type prv: instance of ProvidentiaMainWindow / Interactive
        :fname: Name of the file to save
        :type fname: str
        :input_dialogue: boolean informing whether to open input prompt to ask if to filter data or not
        :type input_dialogue: boolean
        :set_in_memory: booolean informing if saved data is set in memory
        :type set_in_memory: boolean
    """

    # ensure fname has correct extension
    if fname[-4:] != '.npz':
        fname = '{}.npz'.format(fname)

    # open dialog to choose if data is filtered or not
    if input_dialogue:
        title = 'Export data'
        msg = 'Select option'
        options = ['Apply metadata filters and temporal colocation (if active) to exported data', 
                'Do not apply metadata filters and temporal colocation to exported data']
        dialog = InputDialog(prv, title, msg, options)
        selected_option, okpressed = dialog.selected_option, dialog.okpressed
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

        # get valid station indices (from observations because valid stations for the experiment is a 
        # subset of the observations)
        if prv.temporal_colocation:
            valid_station_inds = prv.valid_station_inds_temporal_colocation[networkspeci][prv.observations_data_label]
        else:
            valid_station_inds = prv.valid_station_inds[networkspeci][prv.observations_data_label]

        if apply_filters:
            if prv.reading_ghost:
                save_data_dict['{}_ghost_data'.format(networkspeci)] = np.take(prv.ghost_data_in_memory[networkspeci], 
                    valid_station_inds, axis=1)
            save_data_dict['{}_data'.format(networkspeci)] = np.take(prv.data_in_memory_filtered[networkspeci], 
                valid_station_inds, axis=1)
            save_data_dict['{}_metadata'.format(networkspeci)] = np.take(prv.metadata_in_memory[networkspeci], 
                valid_station_inds, axis=0)
        else:
            if prv.reading_ghost:
                save_data_dict['{}_ghost_data'.format(networkspeci)] = prv.ghost_data_in_memory[networkspeci]
            save_data_dict['{}_data'.format(networkspeci)] = prv.data_in_memory_filtered[networkspeci]
            save_data_dict['{}_metadata'.format(networkspeci)] = prv.metadata_in_memory[networkspeci]

    # save out miscellaneous variables 
    save_data_dict['time'] = prv.time_array
    save_data_dict['data_labels'] = prv.data_labels
    save_data_dict['resolution'] = prv.resolution
    save_data_dict['start_date'] = prv.start_date
    save_data_dict['end_date'] = prv.end_date
    save_data_dict['temporal_colocation'] = prv.temporal_colocation
    save_data_dict['spatial_colocation'] = prv.spatial_colocation
    save_data_dict['filter_species'] = prv.filter_species
    if prv.reading_ghost:
        save_data_dict['ghost_version'] = prv.ghost_version
        save_data_dict['ghost_data_variables'] = prv.ghost_data_vars_to_read

    # save out dict to .npz file
    np.savez(fname, **save_data_dict)

    # if set_in_memory is active, load and return the variable in memory  
    if set_in_memory:   
        data = np.load(fname, allow_pickle=True)
        # delete temporary save file after load
        os.remove(fname)  
        return data                  

def export_netcdf(prv, fname, input_dialogue=False, set_in_memory=False, xarray=False):
    """ Write data and metadata to netcdf file. 
    
        :prv: Instance of providentia
        :type prv: instance of ProvidentiaMainWindow / Interactive
        :fname: Name of the file to save
        :type fname: str
        :input_dialogue: boolean informing whether to open input prompt to ask if to filter data or not
        :type input_dialogue: boolean
        :set_in_memory: booolean informing if saved data is set in memory
        :type set_in_memory: boolean
        :xarray: booolean informing if data is read in xarray format or not
        :type xarray: boolean

    """

    # ensure fname has correct extension
    if fname[-3:] != '.nc':
        fname = '{}.nc'.format(fname)

    # open dialog to choose if data is filtered or not
    if input_dialogue:
        title = 'Export data'
        msg = 'Select option'
        options = ['Apply metadata filters and temporal colocation (if active) to exported data', 
                'Do not apply metadata filters and temporal colocation to exported data']
        dialog = InputDialog(prv, title, msg, options)
        selected_option, okpressed = dialog.selected_option, dialog.okpressed
        if selected_option == options[0]:
            apply_filters = True
        elif selected_option == options[1]:
            apply_filters = False
    else:
        apply_filters = True

    # set up some structural variables
    from GHOST_standards import standard_parameters, get_standard_data, get_standard_metadata
    parameter_dictionary = {}
    for _, param_dict in standard_parameters.items():
        parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
    
    # dictionary to map python types to netcdf types
    type_map = {np.uint8: 'u1', np.uint32: 'u4', object: str,
                np.float32: 'f4', np.float64: 'f8'}

    # start file
    fout = Dataset(fname, 'w', format="NETCDF4")

    # file contents
    if prv.interactive:
        fout.title = 'Saved data from Providentia interactive.'
    else:
        fout.title = 'Saved data from the Providentia dashboard.'
    fout.institution = 'Barcelona Supercomputing Center'
    fout.source = 'Providentia'
    if prv.reading_ghost:
        fout.data_version = prv.ghost_version

    # netcdf dimensions
    fout.createDimension('data_label', len(prv.data_labels))
    fout.createDimension('station', None)
    fout.createDimension('time', len(prv.time_array))
    fout.createDimension('month', len(prv.yearmonths))
    
    # create dimensions only for GHOST case
    if prv.reading_ghost:
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
            var[:] = np.arange(len(prv.time_array))
            
            # miscellaneous variables 
            var = fout.createVariable('data_labels', str, ('data_label',))
            var.standard_name = 'data_labels'
            var.long_name = 'data_labels'
            var.description = 'Labels associated with each data array, e.g. observations, experiment_1, etc.'
            var[:] = np.array(prv.data_labels)
            
            if prv.reading_ghost:
                var = fout.createVariable('ghost_data_variables', str, ('ghost_data_variable',))
                var.standard_name = 'ghost_data_variables'
                var.long_name = 'ghost_data_variables'
                var.description = 'The names of the GHOST data variables used for additional filtering.'
                var[:] = np.array(prv.ghost_data_vars_to_read)
         
        # data
        current_data_type = type_map[data_format_dict[speci]['data_type']]
        var = fout.createVariable('{}_data'.format(var_prefix), current_data_type, 
                                  ('data_label', 'station', 'time',))
        
        # get valid station indices (from observations because valid stations for the experiment is a 
        # subset of the observations)
        if prv.temporal_colocation:
            valid_station_inds = prv.valid_station_inds_temporal_colocation[networkspeci][prv.observations_data_label]
        else:
            valid_station_inds = prv.valid_station_inds[networkspeci][prv.observations_data_label]

        # set attributes
        var.standard_name = data_format_dict[speci]['standard_name']
        var.long_name = data_format_dict[speci]['long_name']
        var.units = data_format_dict[speci]['units']
        var.description = data_format_dict[speci]['description']
        var.resolution = str(prv.resolution)
        var.start_date = str(prv.start_date)
        var.end_date = str(prv.end_date)
        var.temporal_colocation = str(prv.temporal_colocation)
        var.spatial_colocation = str(prv.spatial_colocation)
        var.filter_species = str(prv.filter_species)
        if prv.reading_ghost:
            var.ghost_version = str(prv.ghost_version)
        if apply_filters:
            test = np.take(prv.data_in_memory_filtered[networkspeci], valid_station_inds, axis=1)
            var[:] = np.take(prv.data_in_memory_filtered[networkspeci], valid_station_inds, axis=1)
        else:
            var[:] = prv.data_in_memory[networkspeci]

        # GHOST data
        if prv.reading_ghost:
            var = fout.createVariable('{}_ghost_data'.format(networkspeci), 'f4', 
                                      ('ghost_data_variable', 'station', 'time',))
            # set attributes 
            var.standard_name = '{}_ghost_data'.format(networkspeci) 
            var.long_name = '{}_ghost_data'.format(networkspeci)
            var.description = 'GHOST data variables used for additional filtering.'
            if apply_filters:
                var[:] = np.take(prv.ghost_data_in_memory[networkspeci], valid_station_inds, axis=1)
            else:
                var[:] = prv.ghost_data_in_memory[networkspeci]

        # save metadata (as individual variables)
        if apply_filters:
            metadata_arr = np.take(prv.metadata_in_memory[networkspeci], valid_station_inds, axis=0)
        else:
            metadata_arr = prv.metadata_in_memory[networkspeci]
        
        for metadata_var in metadata_arr.dtype.names:
            
            current_data_type = type_map[metadata_format_dict[metadata_var]['data_type']] 
            var = fout.createVariable('{}_{}'.format(var_prefix, metadata_var), 
                                      current_data_type, ('station', 'month',))

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

def export_configuration(prv, cname, separator="||"):
    """ Create all items to be written in configuration file
        and send them to write_conf.

        :prv: Instance of providentia
        :type prv: instance of ProvidentiaMainWindow / Interactive
        :cname: Name for the configuration file
        :type cname: str
        :separator: delimiter for keep/remove fields
        :type separator: str
    """

    # ensure cname has correct extension
    if cname[-5:] != '.conf':
        cname = '{}.conf'.format(cname)

    # set section and subsection names in config file
    if not hasattr(prv, 'section'):
        section = 'SECTION1'
        subsection = '[SUBSECTION1]'
    else:
        if '·' in prv.section:
            section = prv.section.split('·')[0]
            subsection = '[' + prv.section.split('·')[1] + ']'
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
                          'end_date': prv.end_date,
                          'statistic_mode': prv.statistic_mode,
                          'periodic_statistic_mode': prv.periodic_statistic_mode,
                          'periodic_statistic_aggregation': prv.periodic_statistic_aggregation
                          }

    # statistic aggregation
    if len(prv.statistic_aggregation) > 0:
        options['section']['statistic_aggregation'] = prv.statistic_aggregation

    # experiments
    if prv.experiments_menu['checkboxes']['keep_selected']:
        options['section']['experiments'] = ",".join(str(i) for i in prv.experiments_menu['checkboxes']['keep_selected'])

    # add information about colocation
    options['section'].update({'temporal_colocation': prv.temporal_colocation,
                               'spatial_colocation': prv.spatial_colocation,
                              })

    # add information about filter species if any
    if len(prv.filter_species) > 0:
        filter_species = str(copy.deepcopy(prv.filter_species))
        filter_species = filter_species.replace("[", "(").replace("]", ")")
        filter_species = filter_species.replace("{", "").replace("}", "")
        filter_species = filter_species.replace("'", "")
        filter_species = filter_species.replace(":", "")
        filter_species = filter_species.replace("|", ":")
        filter_species = filter_species.replace("((", "(")
        filter_species = filter_species.replace("))", ")")
        options['section'].update({'filter_species': filter_species
                                  })

    # add information about report
    options['section'].update({'report_type': prv.report_type,
                               'report_summary': prv.report_summary,
                               'report_stations': prv.report_stations,
                               'report_title': prv.report_title,     
                               'report_filename': prv.report_filename 
                              })

    # add other miscellaneous fields
    options['section'].update({'map_extent': ",".join(str(i) for i in prv.map_extent),
                               'active_dashboard_plots': ",".join(str(i) for i in prv.active_dashboard_plots)
                              })
    
    if subsection != None:

        # QA
        if set(prv.qa_menu['checkboxes']['remove_selected']) != set(prv.qa_menu['checkboxes']['remove_default']):
            options['subsection']['QA'] = ",".join(str(i) for i in prv.qa_menu['checkboxes']['remove_selected'])
            
        # flags
        if prv.flag_menu['checkboxes']['remove_selected']:
            options['subsection']['flags'] = ",".join(str(i) for i in prv.flag_menu['checkboxes']['remove_selected'])

        # representativity
        for i, label in enumerate(prv.representativity_menu['rangeboxes']['labels']):
            if 'max_gap' in label:
                if prv.representativity_menu['rangeboxes']['current_lower'][i] != '100':
                    options['subsection'][label] = prv.representativity_menu['rangeboxes']['current_lower'][i]
            else:
                if prv.representativity_menu['rangeboxes']['current_lower'][i] != '0':
                    options['subsection'][label] = prv.representativity_menu['rangeboxes']['current_lower'][i]

        # period
        if prv.period_menu['checkboxes']['keep_selected'] or prv.period_menu['checkboxes']['remove_selected']:
            period_k = "keep: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['keep_selected']) + separator
            period_r = " remove: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['remove_selected']) + separator
            options['subsection']['period'] = period_k + period_r
     
        # bounds
        if np.float32(prv.le_minimum_value.text()) != \
                np.float32(prv.parameter_dictionary[prv.species[0]]['extreme_lower_limit']):
            options['subsection']['lower_bound'] = prv.le_minimum_value.text()
        if np.float32(prv.le_maximum_value.text()) != \
                np.float32(prv.parameter_dictionary[prv.species[0]]['extreme_upper_limit']):
            options['subsection']['upper_bound'] = prv.le_maximum_value.text()

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
                    options['subsection'][label] = lower_cur + ", " + upper_cur

            # and then treat the keep/remove
            for label in prv.metadata_menu[menu_type]['navigation_buttons']['labels']:
                keeps = prv.metadata_menu[menu_type][label]['checkboxes']['keep_selected']
                removes = prv.metadata_menu[menu_type][label]['checkboxes']['remove_selected']
                if keeps or removes:
                    meta_keep = "keep: " + ",".join(str(i) for i in keeps) + separator
                    meta_remove = " remove: " + ",".join(str(i) for i in removes) + separator
                    options['subsection'][label] = meta_keep + meta_remove
    
    write_conf(section, subsection, cname, options)
