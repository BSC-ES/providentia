"""
Contains functions for the processing/calculation of statistics and colourbars
"""
from .calculate import Stats
from .calculate import ExpBias
from .read_aux import drop_nans
from .aux import exceedance_lim, get_relevant_temporal_resolutions

import copy
import json
import os
import sys
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))


def to_pandas_dataframe(read_instance, canvas_instance, networkspecies, 
                        station_index=False, data_range_min=False, data_range_max=False, stddev_max=False):
    """ Function that takes data in memory puts it in a pandas dataframe, per network / species, per data label.

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param networkspecies: List of networkspeci strings
        :type networkspecies: list
        :param station_index: Indices of stations to keep per network/species
        :type station_index: list
        :param data_range_min: current minimum of data range per networkspecies
        :type data_range_min: dict
        :param data_range_max: current maximum of data range per networkspecies
        :type data_range_max: dict
        :param stddev_max: current maximum of StdDev per networkspecies
        :type stddev_max: dict
    """

    if read_instance.resampling:

        # update relevant temporal resolutions 
        read_instance.relevant_temporal_resolutions = get_relevant_temporal_resolutions(read_instance.resampling_resolution)
                        
        # transform resolution to code for .resample function
        if read_instance.resampling_resolution in ['hourly', 'hourly_instantaneous']:
            temporal_resolution_to_output_code = 'H'
        elif read_instance.resampling_resolution in ['3hourly', '3hourly_instantaneous']:
            temporal_resolution_to_output_code = '3H'
        elif read_instance.resampling_resolution in ['6hourly', '6hourly_instantaneous']:
            temporal_resolution_to_output_code = '6H'
        elif read_instance.resampling_resolution == 'daily':
            temporal_resolution_to_output_code = 'D'
        elif read_instance.resampling_resolution == 'monthly':
            temporal_resolution_to_output_code = 'MS'
        elif read_instance.resampling_resolution == 'yearly':
            temporal_resolution_to_output_code = 'AS'
    
    else:
        # update relevant temporal resolutions 
        read_instance.relevant_temporal_resolutions = get_relevant_temporal_resolutions(read_instance.resolution)    

    # create new dictionaries to store selected station data by network / species, per data label
    canvas_instance.selected_station_data = {}
    canvas_instance.selected_station_data_min = {}
    canvas_instance.selected_station_data_max = {}
    canvas_instance.selected_station_stddev_max = {}

    # iterate through networks / species  
    for networkspeci_ii, networkspeci in enumerate(networkspecies):

        # add nested dictionary for networkspeci in selected station data dictionary
        canvas_instance.selected_station_data[networkspeci] = {}
        if data_range_min:
            canvas_instance.selected_station_data_min[networkspeci] = data_range_min[networkspeci]
        else:
            canvas_instance.selected_station_data_min[networkspeci] = np.inf
        if data_range_max:
            canvas_instance.selected_station_data_max[networkspeci] = data_range_max[networkspeci]
        else:
            canvas_instance.selected_station_data_max[networkspeci] = 0.0
        if stddev_max:
            canvas_instance.selected_station_stddev_max[networkspeci] = stddev_max[networkspeci]
        else:
            canvas_instance.selected_station_stddev_max[networkspeci] = 0.0      

        print('getting data array')
        start = time.time()

        # get data array for networkspeci
        data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][:,:,:])

        # temporally colocate data array
        if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
            data_array[:, read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
        
        # get selected station indices
        if station_index:
            read_instance.station_inds = np.array([station_index])
        else:
            if read_instance.offline:
                read_instance.station_inds = np.arange(data_array.shape[1])
            else:
                read_instance.station_inds = canvas_instance.relative_selected_station_inds
        
        # get data cut
        data_array = data_array[:,read_instance.station_inds,:]

        print(time.time()-start)
                    
        # if data array has no valid data for selected stations, do not create a pandas dataframe
        # data array has valid data and is not all nan?
        if data_array.size > 0 and not np.isnan(data_array).all():
                            
            print('creating selected station data')
            start = time.time()

            # apply calibration operation if required
            data_array = apply_calibration_factor(data_array, read_instance, canvas_instance, networkspeci_ii)

            # temporally resample data array if required
            if read_instance.resampling:
                # flatten networkspecies dimension for creation of pandas dataframe
                data_array_reduced = data_array.reshape(data_array.shape[0]*data_array.shape[1], data_array.shape[2])
                
                # create pandas dataframe of data array
                data_array_df = pd.DataFrame(data_array_reduced.transpose(), index=read_instance.time_array, 
                                             columns=np.arange(data_array_reduced.shape[0]), dtype=np.float32)
                # resample data array
                data_array_df_resampled = data_array_df.resample(temporal_resolution_to_output_code, axis=0).mean()

                # save back out as numpy array (reshaping to get back networkspecies dimension)
                data_array_resampled = data_array_df_resampled.to_numpy().transpose()
                data_array = data_array_resampled.reshape(data_array.shape[0],data_array.shape[1],data_array_resampled.shape[1])

            # save timeseries array
            if len(read_instance.station_inds) == 1:
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = data_array[:,0,:]
            else:
                aggregated_data = aggregation(data_array, read_instance.statistic_aggregation, axis=1)
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = aggregated_data

            print('ts', time.time()-start)

            # save data per station
            if read_instance.statistic_mode == 'Spatial|Temporal':
                canvas_instance.selected_station_data[networkspeci]['per_station'] = canvas_instance.selected_station_data[networkspeci]['timeseries'][:,np.newaxis,:]
            elif read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
                canvas_instance.selected_station_data[networkspeci]['per_station'] = data_array

            print('per station', time.time()-start)

            # flatten data across stations
            canvas_instance.selected_station_data[networkspeci]['flat'] = canvas_instance.selected_station_data[networkspeci]['per_station'].reshape(data_array.shape[0],1,data_array.shape[1]*data_array.shape[2])

            print('flat', time.time()-start)

            # set active data array for statistical mode
            if read_instance.statistic_mode in ['Spatial|Temporal', 'Temporal|Spatial']:
                canvas_instance.selected_station_data[networkspeci]['active_mode'] = canvas_instance.selected_station_data[networkspeci]['per_station']
            elif read_instance.statistic_mode == 'Flattened':
                canvas_instance.selected_station_data[networkspeci]['active_mode'] = canvas_instance.selected_station_data[networkspeci]['flat']

            print('active mode', time.time()-start)
                
            #if read_instance.statistic_mode == 'Spatial|Temporal':
            #current_min = np.nanmin(canvas_instance.selected_station_data[networkspeci]['flat'])
            #current_max = np.nanmax(canvas_instance.selected_station_data[networkspeci]['flat'])
            #elif read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
            #    current_min = np.nanpercentile(canvas_instance.selected_station_data[networkspeci][data_label]['flat'],q=5)
            #    current_max = np.nanpercentile(canvas_instance.selected_station_data[networkspeci][data_label]['flat'],q=95)
            canvas_instance.selected_station_data_min[networkspeci] = np.nanmin(canvas_instance.selected_station_data[networkspeci]['flat'])
            canvas_instance.selected_station_data_max[networkspeci] = np.nanmax(canvas_instance.selected_station_data[networkspeci]['flat'])
            canvas_instance.selected_station_stddev_max[networkspeci] = np.nanmax(np.nanstd(canvas_instance.selected_station_data[networkspeci]['flat'], axis=-1))

            print(time.time()-start)

            print('setting periodic chunks')
            start = time.time()

            # group data into periodic chunks
            group_periodic(read_instance, canvas_instance, networkspeci)
            print(time.time()-start)

def group_periodic(read_instance, canvas_instance, networkspeci):
    """ Function that groups data into periodic chunks

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param canvas_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type canvas_instance: object
        :param networkspeci: name of networkspeci str
        :type networkspeci: str
    """

    # iterate through all defined temporal aggregation resolutions
    for temporal_aggregation_resolution in read_instance.relevant_temporal_resolutions:

        # get all temporal periods for current resolution
        all_periods = getattr(read_instance.time_array, temporal_aggregation_resolution)
       
        # get all unique temporal periods for current resolution
        unique_periods = np.unique(all_periods)

        # create lists to store grouped data and valid ticks
        canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution] = {}
        canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['active_mode'] = []
        canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['valid_xticks'] = []

        # iterate through unique temporal periods and store associated data with each period, per data label
        for unique_period in unique_periods:

            # get mask for current period
            valid_period = all_periods == unique_period

            # get associated data with period
            period_data = canvas_instance.selected_station_data[networkspeci]['per_station'][:,:,valid_period]

            # if have valid data for period, append it and xtick for period
            if period_data.size > 0:
                # flatten group for flattened stat  mode
                if read_instance.statistic_mode == 'Flattened':
                    period_data = period_data.reshape(period_data.shape[0],1,period_data.shape[1]*period_data.shape[2])
                    
                canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['active_mode'].append(period_data)
                canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['valid_xticks'].append(unique_period)


def calculate_statistic(read_instance, canvas_instance, zstats, data_labels_a, data_labels_b, map=False, period=None):
    """Function that calculates a statistic for data labels, either absolute or bias, 
       for different aggregation modes.
    """

    start = time.time()

    #if data_labels_a, data_labels_b are strings then convert to lists
    if type(data_labels_a) != list:
        data_labels_a = [data_labels_a]
    if type(data_labels_b) != list:
        data_labels_b = [data_labels_b]

    #if have empty strings in lists then remove them
    data_labels_a = [label for label in data_labels_a if label != '']
    data_labels_b = [label for label in data_labels_b if label != '']

    #if zstats is str then make it a list
    if type(zstats) != list:
        zstats = [zstats]

    # iterate through zstats and calculate statistics
    stats_calc = {}
    for zstat in zstats:

        print('calculating statistic:', zstat)

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        # for map statistics, get active map valid station indices and then data_labels_a data 
        if map:
            # check if have valid station data first
            # if not update z statistic and active map valid station indices to be empty lists and return
            if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                n_valid_stations = len(read_instance.valid_station_inds_temporal_colocation[read_instance.networkspeci]['observations'])
            else:
                n_valid_stations = len(read_instance.valid_station_inds[read_instance.networkspeci]['observations'])
            if n_valid_stations == 0:             
                z_statistic = np.array([], dtype=np.float32)
                active_map_valid_station_inds = np.array([], dtype=np.int)
                return z_statistic, active_map_valid_station_inds 

            # get active map valid station indices (i.e. the indices of the stations data to plot on the map)
            # if only have data_labels_a, valid map indices are those simply for the data_labels_a array
            if len(data_labels_b) == 0:
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    active_map_valid_station_inds = read_instance.valid_station_inds_temporal_colocation[read_instance.networkspeci][data_labels_a[0]]
                else:
                    active_map_valid_station_inds = read_instance.valid_station_inds[read_instance.networkspeci][data_labels_a[0]]
            else:
                # if have data_labels_b, get intersection of data_labels_a and data_labels_b valid station indices
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    active_map_valid_station_inds = \
                        np.intersect1d(read_instance.valid_station_inds_temporal_colocation[read_instance.networkspeci][data_labels_a[0]],
                                    read_instance.valid_station_inds_temporal_colocation[read_instance.networkspeci][data_labels_b[0]])
                else:
                    active_map_valid_station_inds = \
                        np.intersect1d(read_instance.valid_station_inds[read_instance.networkspeci][data_labels_a[0]],
                                    read_instance.valid_station_inds[read_instance.networkspeci][data_labels_b[0]])

            # get data_label_a array data
            data_array_a = read_instance.data_in_memory_filtered[read_instance.networkspeci][read_instance.data_labels.index(data_labels_a[0]),:,:]

            # temporally colocate data (if active)
            if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                data_array_a[read_instance.temporal_colocation_nans[read_instance.networkspeci]] = np.NaN
            # cut for valid stations
            data_array_a = data_array_a[active_map_valid_station_inds,:]

        # for other cases, get cut of selected station data for data_labels_a
        else:
            indices_a = np.array([read_instance.data_labels.index(label) for label in data_labels_a])

            # periodic grouped data
            if period:
                data_array_a = [arr[indices_a] for arr in canvas_instance.selected_station_data[read_instance.networkspeci][period]['active_mode']]
            elif z_statistic_period:
                data_array_a = [arr[indices_a] for arr in canvas_instance.selected_station_data[read_instance.networkspeci][z_statistic_period]['active_mode']]
            # non-periodic grouped data
            else:
                data_array_a = canvas_instance.selected_station_data[read_instance.networkspeci]['active_mode'][indices_a]

        # get dictionary containing necessary information for calculation of selected statistic
        if z_statistic_type == 'basic':
            stats_dict = basic_stats[base_zstat]
        else:
            stats_dict = expbias_stats[base_zstat]

        # if have no data_labels_b, calculate 'absolute' basic statistic
        if len(data_labels_b) == 0:

            # load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            # if stat is exceedances then add threshold value (if available)  
            if base_zstat == 'Exceedances':
                function_arguments['threshold'] = exceedance_lim(read_instance.networkspeci)

            # calculate statistics
            
            # calculate statistics per periodic grouping per station
            if period:
                z_statistic = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments) 
                                        for group in data_array_a])

            # calculate periodic statistic per station
            elif z_statistic_period:
                # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                if read_instance.periodic_statistic_mode == 'Cycle':
                    # aggregation in each group, per station, by periodic statistic
                    z_statistic = np.array([aggregation(group, read_instance.periodic_statistic_aggregation, axis=-1)
                                        for group in data_array_a]).transpose()
                    # calculate statistic per station (removing period dimension)
                    z_statistic = np.array(getattr(Stats, stats_dict['function'])(z_statistic, **function_arguments))

                # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                # and then aggregate 
                elif read_instance.periodic_statistic_mode == 'Independent':
                    # calculate statistic per periodic grouping per station
                    z_statistic = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments)
                                            for group in data_array_a]).transpose()
                    # aggregate data per station (removing period dimension)
                    z_statistic = aggregation(z_statistic, read_instance.periodic_statistic_aggregation, axis=-1)

            # calculate statistics per station 
            else:
                z_statistic = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments))

        # else, get data_labels_b data then calculate 'difference' statistic
        else:

            # get data_labels_b data for map
            if map:
                data_array_b = \
                    copy.deepcopy(read_instance.data_in_memory_filtered[read_instance.networkspeci][read_instance.data_labels.index(data_labels_b[0]),:,:])
                # temporally colocate data (if active)
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    data_array_b[read_instance.temporal_colocation_nans[read_instance.networkspeci]] = np.NaN
                # cut for valid stations
                data_array_b = data_array_b[active_map_valid_station_inds,:]
            # for other cases, get cut of selected station data for data_labels_b
            else:
                indices_b = np.array([read_instance.data_labels.index(label) for label in data_labels_b])

                # periodic grouped data
                if period:
                    data_array_b = [arr[indices_b] for arr in canvas_instance.selected_station_data[read_instance.networkspeci][period]['active_mode']]
                elif z_statistic_period:
                    data_array_b = [arr[indices_b] for arr in canvas_instance.selected_station_data[read_instance.networkspeci][z_statistic_period]['active_mode']]
                # non-periodic grouped data
                else:
                    data_array_b = canvas_instance.selected_station_data[read_instance.networkspeci]['active_mode'][indices_b]
            
            # is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                # load default selected statistic arguments and make separate arguments
                # dictionaries for data_labels_a/data_labels_b calculations (as doing 2 separate calculations for data_labels_a/data_labels_b and subtracting)
                function_arguments_a = stats_dict['arguments']
                # if stat is exceedances then add threshold value (if available)  
                if base_zstat == 'Exceedances':
                    function_arguments_a['threshold'] = exceedance_lim(read_instance.networkspeci)
                function_arguments_b = copy.deepcopy(function_arguments_a)

                # calculate statistics for data_labels_a and data_labels_b, then subtract data_labels_b - data_labels_a
                
                # calculate statistics per periodic grouping per station
                if period:
                    statistic_a = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_a) 
                                            for group in data_array_a])
                    statistic_b = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_b) 
                                            for group in data_array_b])

                # calculate periodic statistic per station
                elif z_statistic_period:
                    # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                    if read_instance.periodic_statistic_mode == 'Cycle':
                        # aggregation in each group, per station, by periodic statistic
                        statistic_a = np.array([aggregation(group, read_instance.periodic_statistic_aggregation, axis=-1)
                                            for group in data_array_a]).transpose()
                        statistic_b = np.array([aggregation(group, read_instance.periodic_statistic_aggregation, axis=-1)
                                            for group in data_array_b]).transpose()
                        
                        # calculate statistic per station (removing period dimension)
                        statistic_a = np.array(getattr(Stats, stats_dict['function'])(statistic_a, **function_arguments_a))
                        statistic_b = np.array(getattr(Stats, stats_dict['function'])(statistic_b, **function_arguments_b))

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif read_instance.periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        statistic_a = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_a)
                                                for group in data_array_a]).transpose()
                        statistic_b = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_b)
                                                for group in data_array_b]).transpose()

                        # aggregate data per station (removing period dimension)
                        statistic_a = aggregation(statistic_a, read_instance.periodic_statistic_aggregation, axis=-1)
                        statistic_b = aggregation(statistic_b, read_instance.periodic_statistic_aggregation, axis=-1)

                # calculate statistics per station 
                else:
                    statistic_a = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments_a))
                    statistic_b = np.array(getattr(Stats, stats_dict['function'])(data_array_b, **function_arguments_b))

                # take difference: statistic_b - statistic_a
                z_statistic = statistic_b - statistic_a

            # else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'expbias':

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # calculate statistics
            
                # calculate statistics per periodic grouping per station
                if period:
                    z_statistic = np.array([getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':group_a,'exp':group_b}})
                                            for group_a, group_b in zip(data_array_a, data_array_b)])

                # calculate periodic statistic per station
                elif z_statistic_period:
                    # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                    if read_instance.periodic_statistic_mode == 'Cycle':
                        # aggregation in each group, per station, by periodic statistic
                        statistic_a = np.array([aggregation(group, read_instance.periodic_statistic_aggregation, axis=-1)
                                            for group in data_array_a]).transpose()
                        statistic_b = np.array([aggregation(group, read_instance.periodic_statistic_aggregation, axis=-1)
                                            for group in data_array_b]).transpose()
                        
                        # calculate statistic per station (removing period dimension)
                        z_statistic = np.array(getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':statistic_a,'exp':statistic_b}}))

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif read_instance.periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        z_statistic = np.array([getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':group_a,'exp':group_b}})
                                                for group_a, group_b in zip(data_array_a, data_array_b)]).transpose()

                        # aggregate data per station (removing period dimension)
                        z_statistic = aggregation(z_statistic, read_instance.periodic_statistic_aggregation, axis=-1)

                # calculate statistics per station 
                else:
                    z_statistic = np.array([getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':group_a,'exp':group_b}})
                                            for group_a, group_b in zip(data_array_a, data_array_b)])

        # if any calculated statistics are infinite, then set them to be NaNs 

        finite_boolean = np.isfinite(z_statistic)
        z_statistic[~finite_boolean] = np.NaN

        # return map statistics
        if map:
            # if any station z statistics come out as NaN/inf, cut z_statistic to remove invalid NaNs/infs, 
            # and also remove respective stations from active map valid station indices
            print(time.time()-start)
            return z_statistic[finite_boolean], active_map_valid_station_inds[finite_boolean] 

        # otherwise, save desired statistic for specific statistical calculation mode 
        else:
            if read_instance.statistic_mode == 'Temporal|Spatial':
                z_statistic = aggregation(z_statistic, read_instance.statistic_aggregation,axis=-1)
            elif read_instance.statistic_mode in ['Flattened', 'Spatial|Temporal']:
                z_statistic = np.squeeze(z_statistic, axis=-1)
            print(time.time()-start)
            stats_calc[zstat] = z_statistic

    #return statistics calculated (if just one statistic then remove dict)
    if len(zstats) == 1:
        stats_calc = stats_calc[zstats[0]]
    return stats_calc

def get_axes_vminmax(axs):
    """ Function that get minimum and maximum of plotted data across relevant axes.

        :param axs: list of relevant axes
        :type axs: list
        :return: minimum plotted value, maximum plotted value
        :rtype: np.float32, np.float32
    """

    # get axes plotted vmin/vmax
    ax_min = []
    ax_max = []
    for ax in axs:
        for collection in ax.collections:
            if ((isinstance(collection, matplotlib.collections.PathCollection)) or 
                (isinstance(collection, matplotlib.collections.QuadMesh))):
                col_array = collection.get_array()
                ax_min.append(np.nanmin(col_array))
                ax_max.append(np.nanmax(col_array))
    plotted_min = np.nanmin(ax_min)
    plotted_max = np.nanmax(ax_max)

    return plotted_min, plotted_max

def generate_colourbar_detail(read_instance, zstat, plotted_min, plotted_max, plot_characteristics, speci):
    """ Function that generates neccessary detail to crate colourbar.

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param zstat: name of statistic
        :type zstat: str
        :param plotted_min: minimum plotted value
        :type plotted_min: np.float32
        :param plotted_max: maximum plotted value
        :type plotted_max: np.float32
        :param plot_characteristics: dictionary of plot characteristics
        :type plot_characteristics: dict
        :param speci: speci to plot
        :type speci: str
        :return: cbar min, cbar max, cbar label, cbar cmap
        :rtype: np.float32, np.float32, str, str
    """

    # get zstat information
    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

    # get dictionary containing necessary information for calculation of selected statistic
    if z_statistic_type == 'basic':
        stats_dict = basic_stats[base_zstat]
    else:
        stats_dict = expbias_stats[base_zstat]
    label_units = stats_dict['units']
    if label_units == 'measurement_units':
        label_units = read_instance.measurement_units[speci]

    # generate z colourbar label
    # first check if have defined label (in this order: 1. configuration file 2. specific for z statistic)
    set_label = False

    #1. check configuration file
    if 'cb_label' in plot_characteristics:
        if plot_characteristics['cb_label']['label'] != '':
            z_label = plot_characteristics['cb_label']['label']
            set_label = True
    #2. get label specific for z statistic
    if not set_label:
        if z_statistic_sign == 'absolute':
            if label_units != '':
                z_label = '{} [{}]'.format(stats_dict['label'], label_units)
            else:
                z_label = copy.deepcopy(stats_dict['label'])
        else:
            if z_statistic_type == 'basic':
                if label_units != '':
                    z_label = '{} bias [{}]'.format(stats_dict['label'], label_units)
                else:
                    z_label = '{} bias'.format(stats_dict['label'])
            else:
                if label_units != '':
                    z_label = '{} [{}]'.format(stats_dict['label'], label_units)
                else:
                    z_label = copy.deepcopy(stats_dict['label'])        

    # set cmap for z statistic
    # first check if have defined cmap (in this order: 1. configuration file 2. specific for z statistic)
    set_cmap = False
    if z_statistic_sign == 'absolute':
        cmap_var_name = 'cmap_absolute'
    else:
        cmap_var_name = 'cmap_bias'
    #1. check configuration file
    if cmap_var_name in plot_characteristics['cb']:
        if plot_characteristics['cb'][cmap_var_name] != '':
            z_colourmap = plot_characteristics['cb'][cmap_var_name]
            set_cmap = True
    #2. get cmap specific for z statistic
    if not set_cmap:
        z_colourmap = stats_dict[cmap_var_name]

    # check if have defined vmin (in this order: 1. configuration file 2. specific for z statistic)
    # if have no defined vmin, then take vmin as minimum range value of calculated statistic
    set_vmin = False
    if z_statistic_sign == 'absolute':
        vmin_var_name = 'vmin_absolute'
    else:
        vmin_var_name = 'vmin_bias'
    #1. check configuration file
    if vmin_var_name in plot_characteristics['cb']:
        if plot_characteristics['cb'][vmin_var_name] != '':
            z_vmin = plot_characteristics['cb'][vmin_var_name]
            set_vmin = True
    #2. get vmin specific for z statistic
    if not set_vmin:
        if vmin_var_name in stats_dict:
            if stats_dict[vmin_var_name] != '':
                z_vmin = stats_dict[vmin_var_name]
                set_vmin = True
    # if have no defined vmin, take vmin as minimum range value of calculated statistic
    if not set_vmin:
        z_vmin = plotted_min

    # check if have defined vmax (in this order: 1. configuration file 2. specific for z statistic)
    # if have no defined vmax, then take vmax as maximum range value of calculated statistic
    set_vmax = False
    if z_statistic_sign == 'absolute':
        vmax_var_name = 'vmax_absolute'
    else:
        vmax_var_name = 'vmax_bias'
    #1. check configuration file
    if vmax_var_name in plot_characteristics['cb']:
        if plot_characteristics['cb'][vmax_var_name] != '':
            z_vmax = plot_characteristics['cb'][vmax_var_name]
            set_vmax = True
    #2. get vmax specific for z statistic
    if not set_vmax:
        if vmax_var_name in stats_dict:
            if stats_dict[vmax_var_name] != '':
                z_vmax = stats_dict[vmax_var_name]
                set_vmax = True
    # if have no defined vmax, take vmax as maximum range value of calculated statistic
    if not set_vmax:
        z_vmax = plotted_max

    # if z statistic is a bias stat, and one of vmin/vmax were not configured,
    # force vmin/vmax to be symmetrical across 0
    if (z_statistic_sign == 'bias') & (not set_vmin) & (not set_vmax):
        limit_stat = np.max(np.abs([z_vmin, z_vmax]))
        z_vmin = -limit_stat
        z_vmax = limit_stat

    return z_vmin, z_vmax, z_label, z_colourmap

def generate_colourbar(read_instance, axs, cb_axs, zstat, plot_characteristics, speci):
    """ Function that generates colourbar.

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param axs: list of relevant axes
        :type axs: list
        :param cb_axs: list of relevant colourbar axes
        :type cb_axs: list
        :param zstat: name of statistic
        :type zstat: str    
        :param plot_characteristics: dictionary of plot characteristics
        :type plot_characteristics: dict
        :param speci: speci to plot
        :type speci: str
    """

    # get plotted vmin and vmax over relevant axes
    plotted_min, plotted_max = get_axes_vminmax(axs)

    # get colourbar limits/label
    z_vmin, z_vmax, z_label, z_colourmap = generate_colourbar_detail(read_instance, zstat, plotted_min, plotted_max, 
                                                                     plot_characteristics, speci)

    # generate colourbar tick array
    tick_array = np.linspace(z_vmin, z_vmax, plot_characteristics['cb']['n_ticks'], endpoint=True)

    # normalise colourbar range (into the 0.0 - 1.0 interval)
    norm = matplotlib.colors.Normalize(vmin=z_vmin, vmax=z_vmax)

    # get cmap (handling discrete cases)
    if z_vmin != z_vmax:
        if plot_characteristics['cb']['discrete']:
            cmap = plt.get_cmap(z_colourmap, plot_characteristics['cb']['n_discrete'])
        else:
            cmap = plt.get_cmap(z_colourmap)
    else:
        cmap = plt.get_cmap(z_colourmap, 1)

    # create colourbar
    for cb_ax in cb_axs:

        # make colourbar on axis
        cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, 
                                              orientation=plot_characteristics['cb']['orientation'], 
                                              ticks=tick_array)

        # set colourbar label
        if 'cb_label' in plot_characteristics:
            cb_label_characteristics = copy.deepcopy(plot_characteristics['cb_label'])
            del cb_label_characteristics['label']
            if plot_characteristics['cb']['orientation'] == 'horizontal':
                cb_label_characteristics['xlabel'] = z_label
                cb_ax.set_xlabel(**cb_label_characteristics)
            else:
                cb_ax.yaxis.set_label_position("right")
                cb_label_characteristics['ylabel'] = z_label
                cb_ax.set_ylabel(**cb_label_characteristics)
           
        # set cb tick params
        if 'cb_tick_params' in plot_characteristics:
            # remove ticks for discrete colourbars
            # we do this because different screen resolutions slightly offset the tick position
            # https://earth.bsc.es/gitlab/ac/Providentia/-/issues/166
            if plot_characteristics['cb']['discrete']:
                plot_characteristics['cb_tick_params']['size'] = 0
            cb.ax.tick_params(**plot_characteristics['cb_tick_params'])

    # update plot axes (to take account of new colourbar vmin/vmax/cmap)
    for ax in axs:
        for collection in ax.collections:
            if ((isinstance(collection, matplotlib.collections.PathCollection)) or 
                (isinstance(collection, matplotlib.collections.QuadMesh))):
                collection.set_clim(vmin=z_vmin,vmax=z_vmax)
                collection.set_cmap(cmap=cmap)

def get_z_statistic_comboboxes(base_zstat, second_data_label=''):
    """ Function that gets appropriate zstat name for selected zstatistic comboboxes.

        :param base_zstat: name of statistic
        :type base_zstat: str   
        :param second_data_label: name if secondary data label (if exists)
        :type second_data_label: str
        :return: zstat name
        :rtype: str
    """
    
    # get zstat sign 
    # this is bias, if second data label has been provided
    if second_data_label == '':
        z_statistic_sign = 'absolute'
    else:
        z_statistic_sign = 'bias'

    # get zstat name
    zstat = copy.deepcopy(base_zstat)
    if z_statistic_sign == 'bias':
        if base_zstat in basic_stats:
            zstat = '{}_bias'.format(base_zstat)
        
    return zstat

def get_z_statistic_type(zstat):
    """ Function that checks if the z statistic is basic or expbias statistic.
    
        :param zstat: name of statistic
        :type zstat: str   
        :return: zstat type
        :rtype: str
    """

    # check if the chosen statistic is a basic statistic
    if zstat in basic_stats.keys():
        return 'basic'
    # if not a basic statistic, it must be an experiment bias statistic
    else:
        return 'expbias'

def get_z_statistic_sign(zstat, zstat_type):
    """ Function that checks if the z statistic is an absolute or bias statistic.

        :param zstat: name of statistic
        :type zstat: str   
        :param zstat_type: type of statistic
        :type zstat_type: str   
        :return: zstat sign
        :rtype: str
    """

    # statistic is bias?
    if ('_bias' in zstat) or (zstat_type == 'expbias'):
        return 'bias'
    # statistic is bias?
    else:
        return 'absolute'

def get_z_statistic_info(plot_type=None, zstat=None):
    """ Get z statistic name, type (basic or expbias), sign (absolute or bias), 
        base name (dropping '_bias' suffix) and period (if any)  
        from plot_type (or known zstat name).
    
        :param plot_type: plot type
        :type plot_type: str
        :param zstat: name of statistic
        :type plot_type: str
        :return zstat name, base zstat name, zstat type, zstat sign, zstat period
        :rtype: str, str, str, str, str
    """

    # have plot_type? Therefore need to extract zstat from plot_type name (if available)
    if plot_type:
        # have zstat in plot_type name?
        if ('-' in plot_type) & ('-violin' not in plot_type):
            # have other options in plot_type?
            if '_' in plot_type:
                # bias plot or not (if so, add bias to zstat)
                if '_bias' not in plot_type:
                    zstat = '-'.join(plot_type.split('_')[0].split('-')[1:])
                else:
                    zstat = '-'.join(plot_type.split('_')[0].split('-')[1:] + '_bias')
            # no other options in plot_type
            else:
                zstat = '-'.join(plot_type.split('-')[1:])
        # otherwise return None for all vars
        else:
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = None, None, None, None, None

    # zstat not 'None'? Then get information for it            
    if zstat:
        # get base name name of zstat, dropping 'bias' suffix, and dropping period
        base_zstat = zstat.split('_bias')[0].split('-')[0]
        # get zstat type (basic or expbias) 
        z_statistic_type = get_z_statistic_type(base_zstat)
        # get zstat sign (absolute or bias)
        z_statistic_sign = get_z_statistic_sign(zstat, z_statistic_type)
        # get zstat period (if any)
        if '-' in zstat:
            z_statistic_period = zstat.split('_bias')[0].split('-')[1]
        else:
            z_statistic_period = None

    return zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period

def apply_calibration_factor(data_array, read_instance, canvas_instance, networkspeci_ii):
    """ Apply calibration factor to add or subtract a number to the experiments, 
        multiply or divide the experiment data by a certain value.
    
        :param data_array: data array to perfrom calibration on
        :type data_array: ndarray
        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param canvas_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type canvas_instance: object
        :param networkspeci: position of networkspeci str in networkspecies
        :type networkspeci: int
    """

    if hasattr(read_instance, 'calibration_factor'):
                    
        # get calibration factor per experiment
        for data_label in read_instance.calibration_factor:
            calibration_factor = read_instance.calibration_factor[data_label]

            # get data label index
            data_label_index = read_instance.data_labels.index(data_label)

            # get calibration factor per networkspeci
            if (len(read_instance.networkspecies) > 1) and (',' in calibration_factor):
                calibration_factor = calibration_factor.split(',')[networkspeci_ii]
            
            print('{0} in {1}'.format(calibration_factor, data_label))
            
            # apply calibration factor
            if '*' in calibration_factor:
                data_array[data_label_index,:,:] *= \
                    float(calibration_factor.replace('*', ''))
            elif '/' in calibration_factor:
                data_array[data_label_index,:,:] /= \
                    float(calibration_factor.replace('/', ''))
            elif '-' in calibration_factor:
                data_array[data_label_index,:,:] -= \
                    float(calibration_factor.replace('-', ''))
            else:
                data_array[data_label_index,:,:] += \
                    float(calibration_factor)

    return data_array
    
def aggregation(data_array, statistic_aggregation, axis=0):
    """ Aggregate data across a the specific axis using a given statistic
    
        :param data_array: array of data
        :type data_array: numpy.ndarray
        :param statistic_aggregation: name of aggregation statistic
        :type statistic_aggregation: str
        :param axis: axis to aggregate across
        :type axis: int
    """

    if statistic_aggregation in ['Mean','']:
        aggregated_data = np.nanmean(data_array, axis=axis)
    elif statistic_aggregation == 'Median':
        aggregated_data = np.nanmedian(data_array, axis=axis)
    elif statistic_aggregation in ['p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']:
        aggregated_data = np.nanpercentile(data_array, 
                                           q=int(statistic_aggregation.split('p')[1]),
                                           axis=axis)
    else:
        error = 'Aggregation statistic {0} is not available. '.format(statistic_aggregation)
        error += 'The options are: Median, Mean, p1, p5, p10, p25, p75, p90, p95 and p99'
        sys.exit(error)

    return aggregated_data
