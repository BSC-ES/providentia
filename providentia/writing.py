""" Module storing writing functions """
import sys

import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date
from .configuration import write_conf


def export_data_npz(mpl_canvas, fname):
    """Function that writes out current data in memory to .npy file"""

    if mpl_canvas.read_instance.reading_nonghost:
        mdata = mpl_canvas.read_instance.datareader.nonghost_metadata
    else:
        mdata = mpl_canvas.read_instance.datareader.metadata_in_memory

    np.savez(fname, data=mpl_canvas.read_instance.data_in_memory_filtered,
             metadata=mdata,
             data_resolution=mpl_canvas.read_instance.active_resolution)


def export_netcdf(mpl_canvas, fname):
    """Write data and metadata to netcdf file"""

    instance = mpl_canvas.read_instance
    sys.path.append('/gpfs/projects/bsc32/AC_cache/obs/ghost/GHOST_standards/{}'
                    .format(instance.ghost_version))
    from GHOST_standards import standard_parameters, get_standard_data, get_standard_metadata
    parameter_dictionary = {}
    for _, param_dict in standard_parameters.items():
        parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict

    speci = instance.active_species
    network = instance.active_network
    start = instance.le_start_date.text()
    end = instance.le_end_date.text()
    relevant_yearmonths = instance.relevant_yearmonths

    # frequency for pandas
    fq = instance.datareader.active_frequency_code

    # create time array in selected resolution between start and end date
    pd_time = pd.date_range(start=start, end=end, freq=fq)[:-1]
    time = np.arange(len(pd_time))

    # dictionary to map python types to netcdf types
    type_map = {np.uint8: 'u1', np.uint32: 'u4', np.object: str,
                np.float32: 'f4', np.float64: 'f8'}

    parameter_details = parameter_dictionary[speci]
    metadata_format_dict = get_standard_metadata(parameter_details)
    data_format_dict = get_standard_data(parameter_details)

    metadata_keys = instance.metadata_vars_to_read
    data_arr = instance.data_in_memory_filtered['observations'][speci]
    metadata_arr = instance.datareader.metadata_in_memory
    expids = instance.experiments_menu['checkboxes']['keep_selected']
    exp_to_write = []
    # change some vars if we're treating nonghost
    if instance.reading_nonghost:
        network = instance.active_network.replace("*", "")
        # metadata_keys = ['station_name', 'latitude', 'longitude', 'altitude']
        metadata_arr = instance.datareader.nonghost_metadata
        metadata_keys = list(metadata_arr.dtype.names)

    # start file
    fout = Dataset(fname+".nc", 'w', format="NETCDF4")

    # file contents
    fout.title = 'Surface {} data in the {} network between {}-{}.'\
        .format(speci, network, start, end)
    fout.institution = 'Barcelona Supercomputing Center'
    fout.source = 'Surface observations'
    fout.conventions = 'CF-1.7'
    fout.data_version = instance.ghost_version

    # netcdf dimensions
    fout.createDimension('station', None)
    fout.createDimension('time', len(time))
    # create month dimension only for GHOST case
    if not instance.reading_nonghost:
        fout.createDimension('month', len(relevant_yearmonths))

    data_keys = ['time', speci]
    for data_key in data_keys:
        current_data_type = type_map[data_format_dict[data_key]['data_type']]
        if data_key == 'time':
            var = fout.createVariable('time', current_data_type, ('time',))
        else:
            var = fout.createVariable(data_key+"_"+network, current_data_type, ('station', 'time'))

        # set variable attributes
        var.standard_name = data_format_dict[data_key]['standard_name']
        var.long_name = data_format_dict[data_key]['long_name']
        var.units = data_format_dict[data_key]['units']
        var.description = data_format_dict[data_key]['description']
        # time variable specific attributes
        if data_key == 'time':
            var.units = 'hours since {}-{}-01 00:00:00'.format(start[:4], start[4:6])
            var.description = 'Time in hours since {}-{}-01 00:00 UTC. Time given refers ' \
                              'to the start of the time window the measurement is representative of ' \
                              '(temporal resolution).'.format(start[:4], start[4:6])
            var.axis = 'T'
            var.calendar = 'standard'
            var.tz = 'UTC'

    if mpl_canvas.temporal_colocation:
        for k in instance.data_in_memory_filtered.keys():
            if 'colocatedto' in k:
                expids.append(k)

    # create vars for exps
    if expids:
        for exp in expids:
            if mpl_canvas.temporal_colocation:
                key = speci + "_" + exp
            else:
                if 'colocatedto' not in exp:
                    key = speci + "_" + exp
                else:
                    continue
            exp_to_write.append(exp)
            var = fout.createVariable(key, current_data_type, ('station', 'time'))
            var.standard_name = data_format_dict[speci]['standard_name']
            var.long_name = data_format_dict[speci]['long_name']
            var.units = data_format_dict[speci]['units']
            var.description = data_format_dict[speci]['description']

    # write station data to netCDF
    for data_key in data_keys:
        if data_key == 'time':
            fout[data_key][:] = time
        else:
            fout[data_key+"_"+network][:, :] = data_arr

    for exp in exp_to_write:
        fout[speci+"_"+exp][:, :] = instance.data_in_memory_filtered[exp][speci]

    # metadata variables
    for metadata_key in metadata_keys:
        current_data_type = type_map[metadata_format_dict[metadata_key]['data_type']]
        if instance.reading_nonghost:
            var = fout.createVariable(metadata_key, current_data_type, ('station',))
        else:
            var = fout.createVariable(metadata_key, current_data_type, ('station', 'month'))

        # set variable attributes
        var.standard_name = metadata_format_dict[metadata_key]['standard_name']
        var.long_name = metadata_format_dict[metadata_key]['long_name']
        var.units = metadata_format_dict[metadata_key]['units']
        var.description = metadata_format_dict[metadata_key]['description']

        # variable specific attributes
        if metadata_key == 'longitude':
            var.axis = 'X'
        elif metadata_key == 'latitude':
            var.axis = 'Y'

    # write station metadata to netCDF
    for metadata_key in metadata_keys:
        if instance.reading_nonghost:
            if fout[metadata_key].dtype == str:
                    fout[metadata_key][:] = metadata_arr[metadata_key].astype(str)
            else:
                fout[metadata_key][:] = metadata_arr[metadata_key]
        else:
            if fout[metadata_key].dtype == str:
                fout[metadata_key][:, :] = metadata_arr[metadata_key].astype(str)
            else:
                fout[metadata_key][:, :] = metadata_arr[metadata_key]

    # close writing to netCDF
    fout.close()


def export_configuration(prv, cname, separator="||"):
    """
    Create all items to be written in configuration file
    and send them to write_conf

    :prv: Instance of providentia main window
    :type prv: instance of ProvidentiaMainWindow

    :cname: Name for the configuration file
    :type cname:

    :separator: delimiter for keep/remove fields
    :type separator: str
    """

    # default
    options = {'network': prv.network,
               'resolution': prv.resolution,
               'matrix': prv.matrix,
               'species': prv.species,
               'start_date': prv.start_date,
               'end_date': prv.end_date}

    # QA
    if set(prv.qa_menu['checkboxes']['remove_selected']) != set(prv.qa_menu['checkboxes']['remove_default']):
        options['qa'] = ",".join(str(i) for i in prv.qa_menu['checkboxes']['remove_selected'])
    # flags
    if prv.flag_menu['checkboxes']['remove_selected']:
        options['flags'] = ",".join(str(i) for i in prv.flag_menu['checkboxes']['remove_selected'])
    # experiments
    if prv.experiments_menu['checkboxes']['keep_selected']:
        options['experiments'] = ",".join(str(i) for i in prv.experiments_menu['checkboxes']['keep_selected'])

    # representativity
    for i, label in enumerate(prv.representativity_menu['rangeboxes']['labels']):
        if 'max_gap' in label:
            if prv.representativity_menu['rangeboxes']['current_lower'][i] != '100':
                options[label] = prv.representativity_menu['rangeboxes']['current_lower'][i]
        else:
            if prv.representativity_menu['rangeboxes']['current_lower'][i] != '0':
                options[label] = prv.representativity_menu['rangeboxes']['current_lower'][i]

    # period
    if prv.period_menu['checkboxes']['keep_selected'] or prv.period_menu['checkboxes']['remove_selected']:
        period_k = "keep: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['keep_selected']) + separator
        period_r = " remove: " + ",".join(str(i) for i in prv.period_menu['checkboxes']['remove_selected']) + separator
        options['period'] = period_k + period_r

    # bounds
    if np.float32(prv.le_minimum_value.text()) != \
            np.float32(prv.parameter_dictionary[prv.active_species]['extreme_lower_limit']):
        options['lower_bound'] = prv.le_minimum_value.text()
    if np.float32(prv.le_maximum_value.text()) != \
            np.float32(prv.parameter_dictionary[prv.active_species]['extreme_upper_limit']):
        options['upper_bound'] = prv.le_maximum_value.text()

    # metadata
    for menu_type in prv.metadata_types:
        # treat ranges first
        for i, label in enumerate(prv.metadata_menu[menu_type]['rangeboxes']['labels']):
            lower_cur = prv.metadata_menu[menu_type]['rangeboxes']['current_lower'][i]
            lower_def = prv.metadata_menu[menu_type]['rangeboxes']['lower_default'][i]
            upper_cur = prv.metadata_menu[menu_type]['rangeboxes']['current_upper'][i]
            upper_def = prv.metadata_menu[menu_type]['rangeboxes']['upper_default'][i]
            if (lower_cur != lower_def) or (upper_cur != upper_def):
                options[label] = lower_cur + ", " + upper_cur

        # and then treat the keep/remove
        for label in prv.metadata_menu[menu_type]['navigation_buttons']['labels']:
            keeps = prv.metadata_menu[menu_type][label]['checkboxes']['keep_selected']
            removes = prv.metadata_menu[menu_type][label]['checkboxes']['remove_selected']

            if keeps or removes:
                meta_keep = "keep: " + ",".join(str(i) for i in keeps) + separator
                meta_remove = " remove: " + ",".join(str(i) for i in removes) + separator
                options[label] = meta_keep + meta_remove

    # map z
    if prv.cb_z_stat.currentText() != prv.basic_z_stats[0]:
        options['map_z'] = prv.cb_z_stat.currentText()

    section_name = cname[cname.rfind("/")+1:]
    write_conf(section_name, cname+".conf", options)
