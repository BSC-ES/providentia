""" Functions for the processing/calculation of statistics and colourbars """

from calendar import monthrange
import copy
import datetime
import json
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

from .calculate import Stats, ExpBias
from .read_aux import (drop_nans, get_nonrelevant_temporal_resolutions, get_relevant_temporal_resolutions, 
                       get_frequency_code, get_lower_resolutions)

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
basic_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.json')))


def get_selected_station_data(read_instance, canvas_instance, networkspecies, 
                              station_index=False, data_range_min=False, data_range_max=False, stddev_max=False):
    """ Function that takes full data array and cuts it for selected stations, per network / species, per data label.

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

    possible_resolutions = ['hourly', 'hourly_instantaneous', '3hourly', '3hourly_instantaneous', 
                            '6hourly', '6hourly_instantaneous', 'daily', 'monthly', 'annual']
    
    if read_instance.resampling_resolution in possible_resolutions:
        
        # update relevant/nonrelevant temporal resolutions 
        read_instance.relevant_temporal_resolutions = get_relevant_temporal_resolutions(read_instance.resampling_resolution)
        read_instance.nonrelevant_temporal_resolutions = get_nonrelevant_temporal_resolutions(read_instance.resampling_resolution)

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
        elif read_instance.resampling_resolution == 'annual':
            temporal_resolution_to_output_code = 'AS'
    
    else:
        # update relevant/nonrelevant temporal resolutions 
        read_instance.relevant_temporal_resolutions = get_relevant_temporal_resolutions(read_instance.resolution)    
        read_instance.nonrelevant_temporal_resolutions = get_nonrelevant_temporal_resolutions(read_instance.resolution) 

    # create new dictionaries to store selected station data and metadata by network / species, per data label
    # and station inds per networkspeci
    canvas_instance.selected_station_data = {}
    canvas_instance.selected_station_metadata = {}
    canvas_instance.selected_station_data_labels = {}
    canvas_instance.selected_station_data_min = {}
    canvas_instance.selected_station_data_max = {}
    canvas_instance.selected_station_stddev_max = {}
    canvas_instance.station_inds = {}

    # iterate through networks / species  
    for networkspeci_ii, networkspeci in enumerate(networkspecies):
        
        # initialise data labels
        canvas_instance.selected_station_data_labels[networkspeci] = []

        # add nested dictionary for networkspeci in selected station data / metadata dictionaries
        canvas_instance.selected_station_data[networkspeci] = {}
        canvas_instance.selected_station_metadata[networkspeci] = {}
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

        # get data array for networkspeci
        read_instance.data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][:,:,:])

        # temporally colocate data array
        if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
            read_instance.data_array[:, read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
        
        # get selected station indices
        canvas_instance.station_inds[networkspeci] = get_station_inds(read_instance, canvas_instance, networkspeci, station_index)

        # get data cut for relevant stations
        read_instance.data_array = read_instance.data_array[:,canvas_instance.station_inds[networkspeci],:]
                    
        # get NaNs in data array
        nan_data_array = np.isnan(read_instance.data_array)

        # if data array has no valid data for selected stations, do not cut data array
        # data array has valid data and is not all nan?
        if read_instance.data_array.size > 0 and not np.all(nan_data_array):

            # set metadata cut for relevant stations
            canvas_instance.selected_station_metadata[networkspeci] = read_instance.metadata_in_memory[networkspeci][canvas_instance.station_inds[networkspeci],:]

            # get which data labels have some valid data
            valid_data_labels_mask = ~np.all(np.all(nan_data_array, axis=-1), axis=-1)
            canvas_instance.selected_station_data_labels[networkspeci] = list(np.array(read_instance.data_labels)[valid_data_labels_mask])

            # cut data array for valid data labels
            read_instance.data_array = read_instance.data_array[valid_data_labels_mask]

            # temporally resample data array if required
            if read_instance.resampling_resolution in possible_resolutions:
                # flatten networkspecies dimension for creation of pandas dataframe
                data_array_reduced = read_instance.data_array.reshape(read_instance.data_array.shape[0]*read_instance.data_array.shape[1], 
                                                                      read_instance.data_array.shape[2])
                
                # create pandas dataframe of data array
                data_array_df = pd.DataFrame(data_array_reduced.transpose(), index=read_instance.time_array, 
                                             columns=np.arange(data_array_reduced.shape[0]), dtype=np.float32)
                # resample data array
                data_array_df_resampled = data_array_df.resample(temporal_resolution_to_output_code, axis=0).mean()
                read_instance.time_index = data_array_df_resampled.index

                # save back out as numpy array (reshaping to get back networkspecies dimension)
                data_array_resampled = data_array_df_resampled.to_numpy().transpose()
                read_instance.data_array = data_array_resampled.reshape(read_instance.data_array.shape[0], read_instance.data_array.shape[1],
                                                                        data_array_resampled.shape[1])
            else:
                read_instance.time_index = read_instance.time_array

            # save timeseries array
            if len(canvas_instance.station_inds[networkspeci]) == 1:
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = read_instance.data_array[:,0,:]
            else:
                if (read_instance.offline) or (read_instance.interactive):
                    timeseries_stat = read_instance.statistic_aggregation
                else:
                    timeseries_stat = canvas_instance.timeseries_stat.currentText()
                aggregated_data = aggregation(read_instance.data_array, timeseries_stat, axis=1)
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = aggregated_data

            # save data per station
            if read_instance.statistic_mode == 'Spatial|Temporal':
                canvas_instance.selected_station_data[networkspeci]['per_station'] = canvas_instance.selected_station_data[networkspeci]['timeseries'][:,np.newaxis,:]
            elif read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
                canvas_instance.selected_station_data[networkspeci]['per_station'] = read_instance.data_array

            # transform timeseries to pandas dataframe
            canvas_instance.selected_station_data[networkspeci]['timeseries'] = pd.DataFrame(canvas_instance.selected_station_data[networkspeci]['timeseries'].T, 
                                                                                             columns=canvas_instance.selected_station_data_labels[networkspeci], 
                                                                                             index=read_instance.time_index)

            # flatten data across stations
            canvas_instance.selected_station_data[networkspeci]['flat'] = canvas_instance.selected_station_data[networkspeci]['per_station'].reshape(canvas_instance.selected_station_data[networkspeci]['per_station'].shape[0],
                                                                                                                                                     1,
                                                                                                                                                     canvas_instance.selected_station_data[networkspeci]['per_station'].shape[1]*canvas_instance.selected_station_data[networkspeci]['per_station'].shape[2])

            # set active data array for statistical mode
            if read_instance.statistic_mode in ['Spatial|Temporal', 'Temporal|Spatial']:
                canvas_instance.selected_station_data[networkspeci]['active_mode'] = canvas_instance.selected_station_data[networkspeci]['per_station']
            elif read_instance.statistic_mode == 'Flattened':
                canvas_instance.selected_station_data[networkspeci]['active_mode'] = canvas_instance.selected_station_data[networkspeci]['flat']

            # set lower/upper limits for specific plots
            # lower limit is always min of the data
            # The upper limit is set to be the inner Tukey fence, 
            # so that limits are not distorted by outlying extreme values
            current_min = np.nanmin(canvas_instance.selected_station_data[networkspeci]['flat'])
            if read_instance.statistic_mode == 'Spatial|Temporal':
                current_max = np.nanmax(canvas_instance.selected_station_data[networkspeci]['flat'])
            elif read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
                lower_inner_fence, upper_inner_fence = boxplot_inner_fences(canvas_instance.selected_station_data[networkspeci]['flat'])
                current_max = upper_inner_fence
            canvas_instance.selected_station_data_min[networkspeci] = current_min
            canvas_instance.selected_station_data_max[networkspeci] = current_max
            canvas_instance.selected_station_stddev_max[networkspeci] = np.nanmax(np.nanstd(canvas_instance.selected_station_data[networkspeci]['flat'], axis=-1))

            # group data into periodic chunks
            group_periodic(read_instance, canvas_instance, networkspeci)
            
            # group data into timeseries chunks
            group_timeseries(read_instance, canvas_instance, networkspeci)
            
def boxplot_inner_fences(data):

    ''' Using adjusted boxplot methodology, calaculate Tukey inner fences of data, which beyond these limits data are 
        considered 'possible outliers'. 

        check is only done when have >= 20 values to ensure have sufficient values to use methodology.
        otherwise the minimum nax maximum of the data are returned

        Tukey's boxplot is a very popular tool for detection of outliers. It reveals the location, spread and skewness of the data.
        The definition of the inner fences is such that the expected percentage values which exceed is close to 0.7% for a normal distribution.
        The method is only recommended to be used if a small number of outliers is presumed (at most 5%), and data is normally distributed.

        See References here:
        https://wis.kuleuven.be/stat/robust/papers/2008/adjboxplot-revision.pdf
        https://www.researchgate.net/publication/277943905_A_Modified_Approach_for_Detection_of_Outliers
        https://en.wikipedia.org/wiki/Medcouple
        https://en.wikipedia.org/wiki/Box_plot#Variations

    '''

    #if have < 20 points then simply return min/max of data 
    if data.size < 20:
        return np.nanmin(data), np.nanmax(data)

    #otherwise, calculated Tukey boxplot inner fences
    else:    

        #calculate the 25th percentile 
        p25 = np.nanpercentile(data, 25)
        #calculate the 75th percentile 
        p75 = np.nanpercentile(data, 75)

        #calculate the interquartile range
        iqr = p75-p25

        #calculate lower/upper inner fences and return values
        lower_inner_fence = p25 - 1.5 * iqr 
        upper_inner_fence = p75 + 1.5 * iqr

        return lower_inner_fence, upper_inner_fence

def get_station_inds(read_instance, canvas_instance, networkspeci, station_index):
    """ Get selected station indices """
        
    if station_index:
        station_inds = np.array([station_index])
    else:
        if (read_instance.offline) or (read_instance.interactive):
            if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                station_inds = read_instance.valid_station_inds_temporal_colocation[networkspeci][read_instance.observations_data_label]
            else:
                station_inds = read_instance.valid_station_inds[networkspeci][read_instance.observations_data_label] 
        else:
            station_inds = canvas_instance.relative_selected_station_inds

    return station_inds

def group_periodic(read_instance, canvas_instance, networkspeci):
    """ Function that groups data into periodic chunks

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param canvas_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type canvas_instance: object
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3)
        :type networkspeci: str
    """

    # iterate through all defined temporal aggregation resolutions
    for temporal_aggregation_resolution in read_instance.relevant_temporal_resolutions:

        # get all temporal periods for current resolution
        all_periods = getattr(read_instance.time_index, temporal_aggregation_resolution)
        
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
            
            # if have valid data for period, append it
            # otherwise, append empty list
            if period_data.size > 0:
                # flatten group for flattened stat mode
                if read_instance.statistic_mode == 'Flattened':
                    period_data = period_data.reshape(period_data.shape[0],1,period_data.shape[1]*period_data.shape[2])
            else:
                period_data = []
            
            canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['active_mode'].append(period_data)
            canvas_instance.selected_station_data[networkspeci][temporal_aggregation_resolution]['valid_xticks'].append(unique_period)


def group_timeseries(read_instance, canvas_instance, networkspeci):

    # update timeseries chunk resolution, to all higher resolutions
    if (read_instance.resampling_resolution == "None") or (read_instance.resampling_resolution is None):
        chunk_resolutions = list(get_lower_resolutions(read_instance.resolution))
    else:
        chunk_resolutions = list(get_lower_resolutions(read_instance.resampling_resolution))
    
    # init dictionary where all data will be saved in chunks
    canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'] = {}

    # iterate through available chunk resolutions  
    for chunk_resolution in chunk_resolutions:
        
        # get new frequency and apply to timeseries data
        timeseries_data = copy.deepcopy(canvas_instance.selected_station_data[networkspeci]['timeseries'])
        new_freq = get_frequency_code(chunk_resolution)
        new_index = timeseries_data.asfreq(new_freq).index

        # init arrays where chunk data will be saved
        canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution] = {}
        canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution]['valid_xticks'] = new_index
        canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution]['active_mode'] = []
        
        # get timeseries data for each chunk
        for i in range(len(new_index)):
            
            # is daily chunk resolution?
            if new_freq == "D":
                start_date = datetime.datetime(year=new_index[i].year, 
                                               month=new_index[i].month, 
                                               day=new_index[i].day, 
                                               hour=0)
                end_date = datetime.datetime(year=new_index[i].year, 
                                             month=new_index[i].month, 
                                             day=new_index[i].day, 
                                             hour=23)
            # is monthly chunk resolution?
            elif new_freq == "MS":
                start_date = datetime.datetime(year=new_index[i].year, 
                                               month=new_index[i].month, 
                                               day=1,
                                               hour=0)
                end_day = monthrange(new_index[i].year, new_index[i].month)[1]
                end_date = datetime.datetime(year=new_index[i].year, 
                                             month=new_index[i].month, 
                                             day=end_day, 
                                             hour=23)
                
            # is annual chunk resolution?
            elif new_freq == "AS":
                start_date = datetime.datetime(year=new_index[i].year, 
                                               month=1, 
                                               day=1,
                                               hour=0)
                end_date = datetime.datetime(year=new_index[i].year, 
                                             month=12, 
                                             day=31, 
                                             hour=23)
            
            time_indices = timeseries_data.index.get_indexer(timeseries_data[start_date:end_date].index)
            timeseries_per_station_data = copy.deepcopy(canvas_instance.selected_station_data[networkspeci]['per_station'])
            timeseries_per_station_data = np.take(timeseries_per_station_data, time_indices, axis=2)
            
            # if have valid data for period, append it
            # otherwise, append empty list
            if  timeseries_per_station_data.size > 0:
                # flatten group for flattened stat mode
                if read_instance.statistic_mode == 'Flattened':
                    timeseries_per_station_data =  timeseries_per_station_data.reshape(timeseries_per_station_data.shape[0], 
                                                                                       1, 
                                                                                       timeseries_per_station_data.shape[1] * timeseries_per_station_data.shape[2])
            else:
                 timeseries_per_station_data = []

            canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution]['active_mode'].append(
                timeseries_per_station_data
            )

def calculate_statistic(read_instance, canvas_instance, networkspeci, zstats, data_labels_a, 
                        data_labels_b, map=False, per_station=False, period=None, chunking=None, chunk_stat=None, chunk_resolution=None):
    """Function that calculates a statistic for data labels, either absolute or bias, 
       for different aggregation modes.
    """

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

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        # for map statistics, get active map valid station indices and then data_labels_a data 
        if (map) or (per_station):
            # check if have valid station data first
            # if not update z statistic and active map valid station indices to be empty lists and return
            if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                n_valid_stations = len(read_instance.valid_station_inds_temporal_colocation[networkspeci][read_instance.observations_data_label])
            else:
                n_valid_stations = len(read_instance.valid_station_inds[networkspeci][read_instance.observations_data_label])
            if n_valid_stations == 0:             
                z_statistic = np.array([], dtype=np.float32)
                if map:
                    active_map_valid_station_inds = np.array([], dtype=np.int64)
                    return z_statistic, active_map_valid_station_inds 
                elif per_station:
                    return z_statistic

            # get active map valid station indices (i.e. the indices of the stations data to plot on the map)
            # if only have data_labels_a, valid map indices are those simply for the data_labels_a array
            if len(data_labels_b) == 0:
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    active_map_valid_station_inds = read_instance.valid_station_inds_temporal_colocation[networkspeci][data_labels_a[0]]
                else:
                    active_map_valid_station_inds = read_instance.valid_station_inds[networkspeci][data_labels_a[0]]
            else:
                # if have data_labels_b, get intersection of data_labels_a and data_labels_b valid station indices
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    active_map_valid_station_inds = \
                        np.intersect1d(read_instance.valid_station_inds_temporal_colocation[networkspeci][data_labels_a[0]],
                                    read_instance.valid_station_inds_temporal_colocation[networkspeci][data_labels_b[0]])
                else:
                    active_map_valid_station_inds = \
                        np.intersect1d(read_instance.valid_station_inds[networkspeci][data_labels_a[0]],
                                    read_instance.valid_station_inds[networkspeci][data_labels_b[0]])

            # get data_label_a array data
            data_array_a = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(data_labels_a[0]),:,:])

            # temporally colocate data (if active)
            if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                data_array_a[read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
            # cut for valid stations
            data_array_a = data_array_a[active_map_valid_station_inds,:]

        # for other cases, get cut of selected station data for data_labels_a
        else:
            indices_a = np.array([canvas_instance.selected_station_data_labels[networkspeci].index(label) for label in data_labels_a])
            # periodic grouped data
            if period:
                data_array_a = [arr[indices_a] for arr in canvas_instance.selected_station_data[networkspeci][period]['active_mode']]
            elif z_statistic_period:
                data_array_a = [arr[indices_a] for arr in canvas_instance.selected_station_data[networkspeci][z_statistic_period]['active_mode']]
            # timeseries chunked data
            elif chunking:
                data_array_a = [arr[indices_a] for arr in canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution]['active_mode']]
            # non-periodic grouped data
            else:
                data_array_a = canvas_instance.selected_station_data[networkspeci]['active_mode'][indices_a]

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
                function_arguments['threshold'] = exceedance_lim(networkspeci)

            # need to do the aggregation inside the calculation of NStations, because if not we
            # get arrays with different time dimensions (e.g. different months have different number of hours)
            # and therefore np.array explodes
            if base_zstat == 'NStations':
                function_arguments['statistic_aggregation'] = read_instance.statistic_aggregation

            # calculate statistics
            
            # calculate statistics per periodic grouping per station
            if period or chunking:
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
                    z_statistic = np.array(getattr(Stats, stats_dict['function'])(z_statistic, **function_arguments)).transpose()

                # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                # and then aggregate 
                elif read_instance.periodic_statistic_mode == 'Independent':
                    # calculate statistic per periodic grouping per station
                    z_statistic = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments)
                                            for group in data_array_a]).transpose()
                    # aggregate data per station (removing period dimension)
                    z_statistic = aggregation(z_statistic, read_instance.periodic_statistic_aggregation, axis=-1).transpose()

            # calculate statistics per station 
            else:
                z_statistic = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments))

        # else, get data_labels_b data then calculate 'difference' statistic
        else:

            # get data_labels_b data for map
            if (map) or (per_station):
                data_array_b = \
                    copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(data_labels_b[0]),:,:])
                # temporally colocate data (if active)
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    data_array_b[read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                # cut for valid stations
                data_array_b = data_array_b[active_map_valid_station_inds,:]
            # for other cases, get cut of selected station data for data_labels_b
            else:
                indices_b = np.array([canvas_instance.selected_station_data_labels[networkspeci].index(label) for label in data_labels_b])

                # periodic grouped data
                if period:
                    data_array_b = [arr[indices_b] for arr in canvas_instance.selected_station_data[networkspeci][period]['active_mode']]
                elif z_statistic_period:
                    data_array_b = [arr[indices_b] for arr in canvas_instance.selected_station_data[networkspeci][z_statistic_period]['active_mode']]
                # timeseries chunked data
                elif chunking:
                    data_array_b = [arr[indices_b] for arr in canvas_instance.selected_station_data[networkspeci]['timeseries_chunks'][chunk_resolution]['active_mode']]
                # non-periodic grouped data
                else:
                    data_array_b = canvas_instance.selected_station_data[networkspeci]['active_mode'][indices_b]
            
            # is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                # load default selected statistic arguments and make separate arguments
                # dictionaries for data_labels_a/data_labels_b calculations (as doing 2 separate calculations for data_labels_a/data_labels_b and subtracting)
                function_arguments_a = stats_dict['arguments']
                # if stat is exceedances then add threshold value (if available)  
                if base_zstat == 'Exceedances':
                    function_arguments_a['threshold'] = exceedance_lim(networkspeci)
                function_arguments_b = copy.deepcopy(function_arguments_a)

                # calculate statistics for data_labels_a and data_labels_b, then subtract data_labels_b - data_labels_a
                
                # calculate statistics per periodic grouping per station
                if period or chunking:
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
                        statistic_a = np.array(getattr(Stats, stats_dict['function'])(statistic_a, **function_arguments_a)).transpose()
                        statistic_b = np.array(getattr(Stats, stats_dict['function'])(statistic_b, **function_arguments_b)).transpose()

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif read_instance.periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        statistic_a = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_a)
                                                for group in data_array_a]).transpose()
                        statistic_b = np.array([getattr(Stats, stats_dict['function'])(group, **function_arguments_b)
                                                for group in data_array_b]).transpose()

                        # aggregate data per station (removing period dimension)
                        statistic_a = aggregation(statistic_a, read_instance.periodic_statistic_aggregation, axis=-1).transpose()
                        statistic_b = aggregation(statistic_b, read_instance.periodic_statistic_aggregation, axis=-1).transpose()

                # calculate statistics per station 
                else:
                    statistic_a = np.array(getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments_a))
                    statistic_b = np.array(getattr(Stats, stats_dict['function'])(data_array_b, **function_arguments_b))

                # take difference: statistic_b - statistic_a
                z_statistic = statistic_b - statistic_a

            # else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'expbias':

                # temporal colocation must be turned on for calculation, so if not return NaNs
                if not read_instance.temporal_colocation:
                    if (map) or (per_station):
                        z_statistic = np.array([], dtype=np.float32)
                        if map:
                            active_map_valid_station_inds = np.array([], dtype=np.int64)
                            return z_statistic, active_map_valid_station_inds
                        elif per_station:
                            return z_statistic
                    else: 
                        if period or chunking: 
                            stats_calc[zstat] = np.full((len(data_array_b),len(data_labels_b)), np.NaN)
                        else:
                            stats_calc[zstat] = np.full((len(data_labels_b)), np.NaN)
                        continue

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # calculate statistics
            
                # calculate statistics per periodic grouping per station
                if period or chunking:
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
                        z_statistic = np.array(getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':statistic_a,'exp':statistic_b}})).transpose()

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif read_instance.periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        z_statistic = np.array([getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':group_a,'exp':group_b}})
                                                for group_a, group_b in zip(data_array_a, data_array_b)]).transpose()

                        # aggregate data per station (removing period dimension)
                        z_statistic = aggregation(z_statistic, read_instance.periodic_statistic_aggregation, axis=-1).transpose()

                # calculate statistics per station 
                else:
                    z_statistic = np.array(getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':data_array_a,'exp':data_array_b}}))

        # if any calculated statistics are infinite, then set them to be NaNs 
        finite_boolean = np.isfinite(z_statistic)
        z_statistic[~finite_boolean] = np.NaN

        # return map statistics
        if map:
            # if any station z statistics come out as NaN/inf, cut z_statistic to remove invalid NaNs/infs, 
            # and also remove respective stations from active map valid station indices
            return z_statistic[finite_boolean], active_map_valid_station_inds[finite_boolean] 

        # return per station statistics
        elif per_station:
            return z_statistic

        # otherwise, save desired statistic for specific statistical calculation mode 
        else:
            if (read_instance.statistic_mode == 'Temporal|Spatial') and (base_zstat != 'NStations'):
                z_statistic = aggregation(z_statistic, read_instance.statistic_aggregation,axis=-1)
            elif read_instance.statistic_mode in ['Flattened', 'Spatial|Temporal']:
                z_statistic = np.squeeze(z_statistic, axis=-1)
            stats_calc[zstat] = z_statistic

    # return statistics calculated (if just one statistic then remove dict)
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
                if (len(col_array) > 0):
                    ax_min.append(np.nanmin(col_array))
                    ax_max.append(np.nanmax(col_array))
        
    if (len(ax_min) > 0) and (len(ax_max) > 0):
        plotted_min = np.nanmin(ax_min)
        plotted_max = np.nanmax(ax_max)
        return plotted_min, plotted_max
    else:
        return np.NaN, np.NaN

def generate_colourbar_detail(read_instance, zstat, plotted_min, plotted_max, plot_characteristics, speci, 
                              only_label=False):
    """ Function that generates neccessary detail to crate colourbar.

        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param zstat: Statistic
        :type zstat: str
        :param plotted_min: minimum plotted value
        :type plotted_min: np.float32
        :param plotted_max: maximum plotted value
        :type plotted_max: np.float32
        :param plot_characteristics: dictionary of plot characteristics
        :type plot_characteristics: dict
        :param speci: speci to plot
        :type speci: str
        :param only_label: boolean if only to return label
        :type only_label: boolean
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
    if label_units == '[measurement_units]':
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
    # return label if only that is wanted
    if only_label:
        return z_label

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
        if isinstance(stats_dict[cmap_var_name], dict):
            if speci in stats_dict[cmap_var_name].keys():
                z_colourmap = stats_dict[cmap_var_name][speci]
            else:
                error = "Error: Colormap needs to be defined for all species, using a dictionary with the "
                error += f"cmap for each speci or a string for all. Colormap for {speci} has not been defined."
                sys.exit(error)
        else:
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
                # 3. get vmin specific for species
                set_vmin = True
                if isinstance(stats_dict[vmin_var_name], dict):
                    if speci in stats_dict[vmin_var_name].keys():
                        z_vmin = stats_dict[vmin_var_name][speci]
                    else:
                        set_vmin = False
                else:
                    z_vmin = stats_dict[vmin_var_name]
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
                set_vmax = True
                # 3. get vmax specific for species
                if isinstance(stats_dict[vmax_var_name], dict):
                    if speci in stats_dict[vmax_var_name].keys():
                        z_vmax = stats_dict[vmax_var_name][speci]
                    else:
                        set_vmax = False
                else:
                    z_vmax = stats_dict[vmax_var_name]
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
        :param zstat: Statistic
        :type zstat: str    
        :param plot_characteristics: dictionary of plot characteristics
        :type plot_characteristics: dict
        :param speci: speci to plot
        :type speci: str
    """

    # get plotted vmin and vmax over relevant axes
    plotted_min, plotted_max = get_axes_vminmax(axs)
    if np.isnan(plotted_min) or np.isnan(plotted_max):
        for cb_ax in cb_axs:
            cb_ax.axis('off')
            cb_ax.set_visible(False)
        return

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

def get_z_statistic_comboboxes(base_zstat, bias=False):
    """ Function that gets appropriate zstat name for selected zstatistic comboboxes.

        :param base_zstat: Statistic
        :type base_zstat: str   
        :param second_data_label: name if secondary data label (if exists)
        :type second_data_label: str
        :return: zstat name
        :rtype: str
    """
    
    # get zstat sign 
    if not bias:
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
    
        :param zstat: Statistic
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

def get_z_statistic_sign(zstat, zstat_type=None):
    """ Function that checks if the z statistic is an absolute or bias statistic.

        :param zstat: Statistic
        :type zstat: str   
        :param zstat_type: type of statistic
        :type zstat_type: str   
        :return: zstat sign
        :rtype: str
    """
    
    if zstat_type is None:
        zstat_type = get_z_statistic_type(zstat)

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
        :param zstat: Statistic
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
                    zstat = plot_type.split('_')[0].split('-')[1]
                else:
                    zstat = plot_type.split('_')[0].split('-')[1] + '_bias'
            # no other options in plot_type
            else:
                zstat = plot_type.split('-')[1]
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
            if z_statistic_period == 'diurnal':
                z_statistic_period = 'hour'
            elif z_statistic_period == 'weekly':
                z_statistic_period ='dayofweek'
            elif z_statistic_period == 'monthly':
                z_statistic_period = 'month'
        else:
            z_statistic_period = None

    return zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period
    
def aggregation(data_array, statistic_aggregation, axis=0):
    """ Aggregate data across a the specific axis using a given statistic
    
        :param data_array: array of data
        :type data_array: numpy.ndarray
        :param statistic_aggregation: name of aggregation statistic
        :type statistic_aggregation: str
        :param axis: axis to aggregate across
        :type axis: int
    """

    if statistic_aggregation in ['Median', '']:
        aggregated_data = np.nanmedian(data_array, axis=axis)
    elif statistic_aggregation == 'Mean':
        aggregated_data = np.nanmean(data_array, axis=axis)
    elif statistic_aggregation in ['p1', 'p5', 'p10', 'p25', 'p75', 'p90', 'p95', 'p99']:
        aggregated_data = np.nanpercentile(data_array, 
                                           q=int(statistic_aggregation.split('p')[1]),
                                           axis=axis)
    else:
        error = 'Aggregation statistic {0} is not available. '.format(statistic_aggregation)
        error += 'The options are: Mean, Median, p1, p5, p10, p25, p75, p90, p95 and p99'
        sys.exit(error)

    return aggregated_data


def exceedance_lim(networkspeci):
    """ Return the exceedance limit depending on the species input. 
        If species doesn't have a reported limit, returns np.NaN.

        Try to get limit for specific networkspeci first, and then species.

        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3)
        :type networkspeci: str
        :return: value of exceedance limit
        :rtype: int
    """

    # get speci
    speci = networkspeci.split('|')[1]

    exceedance_limits = json.load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/exceedances.json')))
    if networkspeci in exceedance_limits:
        return exceedance_limits[networkspeci]
    elif speci in exceedance_limits:
        return exceedance_limits[speci]
    else:
        return np.NaN
