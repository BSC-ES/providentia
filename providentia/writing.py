""" Module storing writing functions """

import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date
from .configuration import write_conf
from .dashboard_aux import InputDialog


def export_data_npz(canvas_instance, fname):
    """ Function that writes out current data / ghost data / metadata in memory to .npy file. """

    # open dialog to choose if data is filtered or not
    title = 'Export data'
    msg = 'Select option'
    options = ['Apply metadata filters to exported data', 
               'Do not apply metadata filters to exported data']
    dialog = InputDialog(canvas_instance.read_instance, title, msg, options)
    selected_option, okpressed = dialog.selected_option, dialog.okpressed
    if selected_option == options[0]:
        apply_filters = True
    elif selected_option == options[1]:
        apply_filters = False

    # create dict to save data
    save_data_dict = {}

    # save data / ghost data / metadata
    for networkspeci in canvas_instance.read_instance.networkspecies:
        if apply_filters:
            if canvas_instance.read_instance.reading_ghost:
                save_data_dict['{}_ghost_data'.format(networkspeci)] = canvas_instance.read_instance.ghost_data_in_memory[networkspeci][:, ~np.isnan(canvas_instance.read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))]
            save_data_dict['{}_data'.format(networkspeci)] = canvas_instance.read_instance.data_in_memory_filtered[networkspeci][:, ~np.isnan(canvas_instance.read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))]
            stations_after_filter_inds = np.array(np.where(~np.isnan(canvas_instance.read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))))
            stations_after_filter_inds = np.reshape(stations_after_filter_inds, stations_after_filter_inds.shape[1])
            save_data_dict['{}_metadata'.format(networkspeci)] = canvas_instance.read_instance.metadata_in_memory[networkspeci][stations_after_filter_inds, :]
        else:
            if canvas_instance.read_instance.reading_ghost:
                save_data_dict['{}_ghost_data'.format(networkspeci)] = canvas_instance.read_instance.ghost_data_in_memory[networkspeci]
            save_data_dict['{}_data'.format(networkspeci)] = canvas_instance.read_instance.data_in_memory_filtered[networkspeci]
            save_data_dict['{}_metadata'.format(networkspeci)] = canvas_instance.read_instance.metadata_in_memory[networkspeci]

    # save out miscellaneous variables 
    save_data_dict['time'] = canvas_instance.read_instance.time_array
    save_data_dict['data_labels'] = canvas_instance.read_instance.data_labels
    save_data_dict['resolution'] = canvas_instance.read_instance.resolution
    save_data_dict['start_date'] = canvas_instance.read_instance.start_date
    save_data_dict['end_date'] = canvas_instance.read_instance.end_date
    save_data_dict['temporal_colocation'] = canvas_instance.read_instance.temporal_colocation
    save_data_dict['spatial_colocation'] = canvas_instance.read_instance.spatial_colocation
    save_data_dict['filter_species'] = canvas_instance.read_instance.filter_species
    if canvas_instance.read_instance.reading_ghost:
        save_data_dict['ghost_version'] = canvas_instance.read_instance.ghost_version
        save_data_dict['ghost_data_variables'] = canvas_instance.read_instance.ghost_data_vars_to_read

    # save out dict to .npz file
    np.savez(fname, **save_data_dict)

def export_netcdf(canvas_instance, fname):
    """ Write data and metadata to netcdf file. """

    # open dialog to choose if data is filtered or not
    title = 'Export data'
    msg = 'Select option'
    options = ['Apply metadata filters to exported data', 
               'Do not apply metadata filters to exported data']
    dialog = InputDialog(canvas_instance.read_instance, title, msg, options)
    selected_option, okpressed = dialog.selected_option, dialog.okpressed
    if selected_option == options[0]:
        apply_filters = True
    elif selected_option == options[1]:
        apply_filters = False
        
    # set up some structural variables
    read_instance = canvas_instance.read_instance
    sys.path.append('/gpfs/projects/bsc32/AC_cache/obs/ghost/GHOST_standards/{}'
                    .format(read_instance.ghost_version))
    from GHOST_standards import standard_parameters, get_standard_data, get_standard_metadata
    parameter_dictionary = {}
    for _, param_dict in standard_parameters.items():
        parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
    
    # dictionary to map python types to netcdf types
    type_map = {np.uint8: 'u1', np.uint32: 'u4', np.object: str,
                np.float32: 'f4', np.float64: 'f8'}

    # start file
    fout = Dataset(fname+".nc", 'w', format="NETCDF4")

    # file contents
    fout.title = 'Saved data from the Providentia dashboard.'
    fout.institution = 'Barcelona Supercomputing Center'
    fout.source = 'Providentia'
    if read_instance.reading_ghost:
        fout.data_version = read_instance.ghost_version

    # netcdf dimensions
    fout.createDimension('data_label', len(read_instance.data_labels))
    fout.createDimension('station', None)
    fout.createDimension('time', len(read_instance.time_array))
    fout.createDimension('month', len(read_instance.yearmonths))
    
    # create dimensions only for GHOST case
    if read_instance.reading_ghost:
        fout.createDimension('ghost_data_variable', len(read_instance.ghost_data_vars_to_read))

    # iterate through networkspecies 
    for speci_ii, networkspeci in enumerate(read_instance.networkspecies):

        # get network / species
        network = networkspeci.split('|')[0]
        speci = networkspeci.split('|')[1]

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
            if 'hourly' in read_instance.resolution:
                res_str = 'hours'
            elif 'daily' in read_instance.resolution:
                res_str = 'days'
            elif 'monthly' in read_instance.resolution:
                res_str = 'months'
            var.standard_name = data_format_dict['time']['standard_name']
            var.long_name = data_format_dict['time']['long_name']
            var.units = '{} since {}-{}-01 00:00:00'.format(res_str, 
                                                            str(read_instance.start_date)[:4], 
                                                            str(read_instance.start_date)[4:6])
            msg = 'Time in {} since {}-{}-01 00:00 UTC. Time given refers '.format(res_str, 
                                                                                   str(read_instance.start_date)[:4], 
                                                                                   str(read_instance.start_date)[4:6])
            msg += 'to the start of the time window the measurement is representative of '
            msg += '(temporal resolution).'
            var.description = msg
            var.axis = 'T'
            var.calendar = 'standard'
            var.tz = 'UTC'
            var[:] = np.arange(len(read_instance.time_array))
            
            # miscellaneous variables 
            var = fout.createVariable('data_labels', str, ('data_label',))
            var.standard_name = 'data_labels'
            var.long_name = 'data_labels'
            var.description = 'Labels associated with each data array, e.g. observations, experiment_1, etc.'
            var[:] = np.array(read_instance.data_labels)
            
            var = fout.createVariable('ghost_data_variables', str, ('ghost_data_variable',))
            var.standard_name = 'ghost_data_variables'
            var.long_name = 'ghost_data_variables'
            var.description = 'The names of the GHOST data variables used for additional filtering.'
            var[:] = np.array(read_instance.ghost_data_vars_to_read)
         
        # data
        current_data_type = type_map[data_format_dict[speci]['data_type']]
        var = fout.createVariable('{}_data'.format(networkspeci), current_data_type, 
                                  ('data_label', 'station', 'time',))
        
        # set attributes
        var.standard_name = data_format_dict[speci]['standard_name']
        var.long_name = data_format_dict[speci]['long_name']
        var.units = data_format_dict[speci]['units']
        var.description = data_format_dict[speci]['description']
        var.resolution = str(read_instance.resolution)
        var.start_date = str(read_instance.start_date)
        var.end_date = str(read_instance.end_date)
        var.temporal_colocation = str(read_instance.temporal_colocation)
        var.spatial_colocation = str(read_instance.spatial_colocation)
        var.filter_species = str(read_instance.filter_species)
        if read_instance.reading_ghost:
            var.ghost_version = str(read_instance.ghost_version)
        if apply_filters:
            var[:] = read_instance.data_in_memory_filtered[networkspeci][:, ~np.isnan(read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))]
        else:
            var[:] = read_instance.data_in_memory[networkspeci]

        # GHOST data
        if read_instance.reading_ghost:
            var = fout.createVariable('{}_ghost_data'.format(networkspeci), 'f4', 
                                      ('ghost_data_variable', 'station', 'time',))
            # set attributes 
            var.standard_name = '{}_ghost_data'.format(networkspeci) 
            var.long_name = '{}_ghost_data'.format(networkspeci)
            var.description = 'GHOST data variables used for additional filtering.'
            if apply_filters:
                var[:] = read_instance.ghost_data_in_memory[networkspeci][:, ~np.isnan(read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))]
            else:
                var[:] = read_instance.ghost_data_in_memory[networkspeci]

        # save metadata (as individual variables)
        if apply_filters:
            stations_after_filter_inds = np.array(np.where(~np.isnan(read_instance.data_in_memory_filtered[networkspeci]).all(axis=(0, -1))))
            stations_after_filter_inds = np.reshape(stations_after_filter_inds, stations_after_filter_inds.shape[1])
            metadata_arr = read_instance.metadata_in_memory[networkspeci][stations_after_filter_inds, :]
        else:
            metadata_arr = read_instance.metadata_in_memory[networkspeci]
        for metadata_var in metadata_arr.dtype.names:
            current_data_type = type_map[metadata_format_dict[metadata_var]['data_type']]
            var = fout.createVariable('{}_{}'.format(networkspeci, metadata_var), 
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

def export_configuration(prv, cname, separator="||"):
    """ Create all items to be written in configuration file
        and send them to write_conf.

        :prv: Instance of providentia main window
        :type prv: instance of ProvidentiaMainWindow
        :cname: Name for the configuration file
        :type cname:
        :separator: delimiter for keep/remove fields
        :type separator: str
    """

    # set section and subsection names in config file
    if not hasattr(prv, 'section'):
        section = 'SECTION1'
        subsection = '[SUBSECTION1]'
    else:
        if '|' in prv.section:
            section = prv.section.split('|')[0]
            subsection = '[' + prv.section.split('|')[1] + ']'
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
                          }

    # experiments
    if prv.experiments_menu['checkboxes']['keep_selected']:
        options['section']['experiments'] = ",".join(str(i) for i in prv.experiments_menu['checkboxes']['keep_selected'])

    # add information about colocation
    options['section'].update({'temporal_colocation': prv.temporal_colocation,
                               'spatial_colocation': prv.spatial_colocation,
                              })

    # add information about filter species if any
    if len(prv.filter_species) > 0:
        options['section'].update({'filter_species': prv.filter_species
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
                               'active_dashboard_plots': ",".join(str(i) for i in prv.active_dashboard_plots),
                               'plot_characteristics_filename': prv.plot_characteristics_filename
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
    
    write_conf(section, subsection, cname + '.conf', options)