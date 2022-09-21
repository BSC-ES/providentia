"""
Contains functions for the processing/calculation of statistics and colourbars
"""
from .calculate import Stats
from .calculate import ExpBias
from .read_aux import drop_nans
from .aux import exceedance_lim

import copy
import json
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
basic_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/basic_stats.json')))
expbias_stats = json.load(open(os.path.join(CURRENT_PATH, 'conf/experiment_bias_stats.json')))


def to_pandas_dataframe(read_instance, canvas_instance, networkspecies, station_index=False):
    """Function that takes data in memory puts it in a pandas dataframe, per network / species, per data label.
    For summary plots this involves take the median timeseries across the timeseries.
    For station plots it is just the station in question.
    Also temporally aggregate selected data dataframes (by hour, day of week, month),
    and if have some experiment data associated with selected stations, calculate
    temporally aggregated basic statistic differences and bias statistics between
    observations and experiment data arrays.

    :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type read_instance: object
    :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
    :type canvas_instance: object
    :param networkspecies: List of networkspeci strings
    :type networkspecies: list
    :param station_index: Indices of stations to keep per network/species
    :type station_index: list
    """

    # create new dictionary to store selection station data by network / species, per data label
    canvas_instance.selected_station_data = {}
    canvas_instance.selected_station_data_min = {}
    canvas_instance.selected_station_data_max = {}
    canvas_instance.selected_station_data_number_non_nan = {}

    # iterate through networks / species  
    for networkspeci_ii, networkspeci in enumerate(networkspecies):

        # add nested dictionary for networkspeci in selected station data dictionary
        canvas_instance.selected_station_data[networkspeci] = {}
        canvas_instance.selected_station_data_min[networkspeci] = np.inf
        canvas_instance.selected_station_data_max[networkspeci] = 0
        canvas_instance.selected_station_data_number_non_nan[networkspeci] = []

        # iterate through data labels
        for data_label in read_instance.data_labels:

            # get data for selected stations
            if station_index:
                data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(data_label),:,:])
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    data_array[read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                data_array = data_array[station_index,:]
            else:
                if read_instance.offline:
                    if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                        station_inds = read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label]
                    else:
                        station_inds = read_instance.valid_station_inds[networkspeci][data_label]
                else:
                    if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                        station_inds = np.intersect1d(canvas_instance.relative_selected_station_inds, read_instance.valid_station_inds_temporal_colocation[networkspeci][data_label])
                    else:
                        station_inds = np.intersect1d(canvas_instance.relative_selected_station_inds, read_instance.valid_station_inds[networkspeci][data_label])
                data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(data_label),:,:])
                if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
                    data_array[read_instance.temporal_colocation_nans[networkspeci]] = np.NaN
                data_array = data_array[station_inds,:]
                    
            # if data array has no valid data for selected stations, do not create a pandas dataframe
            # data array has valid data and is not all nan?
            if data_array.size and not np.isnan(data_array).all():

                # add nested dictionary for data label in selected station data dictionary
                canvas_instance.selected_station_data[networkspeci][data_label] = {}
                
                # take cross station median of selected data for data array, and place it in a pandas dataframe 
                if station_index:
                    canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df'] = pd.DataFrame(data_array, index=read_instance.time_array, columns=['data'])
                else:
                    canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df'] = pd.DataFrame(np.nanmedian(data_array, axis=0), index=read_instance.time_array, columns=['data'])

                # get min / max across all selected station data per network / species
                current_min = canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df']['data'].min()
                current_max = canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df']['data'].max()
                if current_min < canvas_instance.selected_station_data_min[networkspeci]:
                    canvas_instance.selected_station_data_min[networkspeci] = current_min
                if current_max > canvas_instance.selected_station_data_max[networkspeci]:
                    canvas_instance.selected_station_data_max[networkspeci] = current_max

                # get number of points per data label
                n_points = len(canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df'].dropna())
                canvas_instance.selected_station_data_number_non_nan[networkspeci].append(n_points)

        # temporally aggregate selected data dataframes (if have periodic plot to make) 
        # (by hour, day of week, month)
        pandas_temporal_aggregation(read_instance, canvas_instance, networkspeci)

        # if have some experiment data associated with selected stations, calculate
        # temporally aggregated basic statistic differences and bias statistics between
        # observations and experiment data arrays
        if len(read_instance.data_labels) > 1:
            calculate_temporally_aggregated_experiment_statistic(read_instance, canvas_instance, networkspeci)
            
def pandas_temporal_aggregation(read_instance, canvas_instance, networkspeci):
    """Function that aggregates pandas dataframe data, for all data arrays,
    into desired temporal groupings also calculates all defined basic
    statistics for each individual temporal grouping

    :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type read_instance: object
    :param canvas_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type canvas_instance: object
    :param networkspeci: name of networkspeci str
    :type networkspeci: str
    """

    # define statistics to calculate (all basic statistics)
    statistics_to_calculate = list(basic_stats.keys())

    # define all temporal aggregation resolutions that will be used to aggregate data
    relevant_temporal_resolutions = read_instance.relevant_temporal_resolutions + ['all']

    # iterate through all defined temporal aggregation resolutions
    for temporal_aggregation_resolution in relevant_temporal_resolutions:

        # define all possible xticks for temporal resolution 
        if temporal_aggregation_resolution != 'all':
            all_xticks = canvas_instance.periodic_xticks[temporal_aggregation_resolution]

        # iterate through data arrays names in selected station data dictionary
        for data_label in canvas_instance.selected_station_data[networkspeci]:

            # create nested dictionary inside selected station data dictionary for storing
            # aggregated data by data array label and temporal aggregation resolution
            canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution] = {}

            # if temporal group is all then simply take data in pandas df as is
            if temporal_aggregation_resolution == 'all':
                full_grouped_data = [canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df']['data'].dropna()]
            else:
                # else, aggregate data array into desired temporal groups (dropping NaNs)
                grouped_data = [g['data'].dropna() for n, g in
                                canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df'].groupby(
                                    getattr(canvas_instance.selected_station_data[networkspeci][data_label]['pandas_df'].index,
                                            temporal_aggregation_resolution))]
                # drop groups which have no data
                grouped_data = [group for group in grouped_data if len(group) > 0]
                # get xticks for groups which have valid data (i.e. which hours/days/months have valid data)
                valid_xticks = [getattr(group.index, temporal_aggregation_resolution)[0] for group in grouped_data]
                # create array of the size of the full range of each aggregation periods,
                # initialised with empty lists per element (i.e. 24 for hourly aggregation)
                full_grouped_data = [[] for _ in range(len(all_xticks))]
                # place valid grouped data in correct positions within full array
                full_group_indices_to_place = np.array([np.where(all_xticks == valid_xtick)[0][0]
                                                        for valid_xtick in valid_xticks], dtype=np.int)
                for grouped_data_ii, full_group_index_to_place in enumerate(full_group_indices_to_place):
                    full_grouped_data[full_group_index_to_place] = grouped_data[grouped_data_ii]
                # add valid xticks for group to selected data dictionary
                # (i.e. the group xtick indexes which have valid data)
                canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution]['valid_xticks'] = valid_xticks

            # add full grouped data to selected data dictionary
            canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution]['grouped_data'] = full_grouped_data

            # calculate basic statistics in each group and add them to selected station data dictionary
            for stat in statistics_to_calculate:
                # get specific statistic dictionary (containing necessary
                # information for calculation of selected statistic)
                stats_dict = basic_stats[stat]
                # load default statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']
                # if stat is exceedances then add threshold value (if available)  
                if stat == 'Exceedances':
                    function_arguments['threshold'] = exceedance_lim(networkspeci)
                # create empty array for storing calculated statistic by group
                stat_output_by_group = []
                # iterate through grouped data
                for group in full_grouped_data:
                    # only calculate statistic if have valid data in group
                    if len(group) > 0:
                        # calculate statistic (appending to all group statistic output array)
                        stat_output_by_group = \
                            np.append(stat_output_by_group,
                                        getattr(Stats, stats_dict['function'])(group, **function_arguments))
                    # if no valid data in group, append NaN
                    else:
                        stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                # save statistical output by group to selected station data dictionary
                canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution][stat] = stat_output_by_group

def calculate_temporally_aggregated_experiment_statistic(read_instance, canvas_instance, networkspeci):
    """Function that calculates temporally aggregated basic statistic
    differences and bias statistics between observations and experiment
    data arrays

    :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type read_instance: object
    :param canvas_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type canvas_instance: object
    :param networkspeci: name of networkspeci str
    :type networkspeci: str
    """

    # define all basic statistics that will be subtracted
    # (each experiment - observations) for each temporal aggregation
    basic_statistics = list(basic_stats.keys())
    # define all experiment bias statistics that will be calculated
    # between each experiment and observations for each temporal aggregation
    bias_statistics = list(expbias_stats.keys())

    # define all temporal aggregation resolutions that will be used to aggregate data
    relevant_temporal_resolutions = read_instance.relevant_temporal_resolutions + ['all']

    # iterate through all defined temporal aggregation resolutions
    for temporal_aggregation_resolution in relevant_temporal_resolutions:
        # iterate through data arrays names in selected station data dictionary
        for data_label in canvas_instance.selected_station_data[networkspeci]:
            # make sure the data array is an experimental one
            if data_label != 'observations':
                relevant_aggregated_observations_dict = \
                    canvas_instance.selected_station_data[networkspeci]['observations'][temporal_aggregation_resolution]

                # calculate temporally aggregated basic statistic differences between experiment and observations
                # iterate through basic statistics
                for basic_stat in basic_statistics:
                    # create empty array for storing calculated experiment-observations
                    # difference statistic by group
                    stat_diff_by_group = []
                    # iterate through aggregated index times
                    for group_ii in range(len(relevant_aggregated_observations_dict[basic_stat])):
                        # get observational and experiment group aggregated statistic
                        group_obs_stat = relevant_aggregated_observations_dict[basic_stat][group_ii]
                        group_exp_stat = canvas_instance.selected_station_data[networkspeci][data_label][
                            temporal_aggregation_resolution][basic_stat][group_ii]

                        # take difference between observations and experiment statistics, if both values finite
                        if (np.isfinite(group_obs_stat)) & (np.isfinite(group_exp_stat)):
                            # calculate difference statistic (experiment - observations)
                            stat_diff_by_group = np.append(stat_diff_by_group, group_exp_stat-group_obs_stat)
                        # else, if one (or both) of observations/experiment statistics are NaN, append NaN
                        else:
                            stat_diff_by_group = np.append(stat_diff_by_group, np.NaN)
                    # save statistical difference output by group to selected station data dictionary
                    canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution][
                        '{}_bias'.format(basic_stat)] = stat_diff_by_group

                # if colocation is active, calculate temporally aggregated experiment bias
                # statistical differences between experiment and observations
                if read_instance.temporal_colocation:
                    # iterate through bias statistics
                    for bias_stat in bias_statistics:
                        # get specific statistic dictionary (containing necessary information
                        # for calculation of selected statistic)
                        stats_dict = expbias_stats[bias_stat]
                        # load default statistic arguments for passing to statistical function
                        function_arguments = stats_dict['arguments']
                        # create empty array for storing calculated statistic by group
                        stat_output_by_group = []
                        # iterate through experimental grouped data
                        for group_ii, exp_group in enumerate(canvas_instance.selected_station_data[networkspeci][data_label]
                                                                [temporal_aggregation_resolution]['grouped_data']):
                            # add aggregated observational and experiment group data as
                            # arguments to pass to statistical function
                            function_arguments['obs'] = relevant_aggregated_observations_dict['grouped_data'][
                                group_ii]
                            function_arguments['exp'] = exp_group

                            # calculate experiment bias statistic between observations and
                            # experiment group data
                            try:
                                # calculate experiment bias statistic
                                stat_output_by_group = \
                                    np.append(stat_output_by_group,
                                                getattr(ExpBias, stats_dict['function'])(**function_arguments))
                            # if can not calculate statistic, append NaN
                            except Exception as e:
                                stat_output_by_group = np.append(stat_output_by_group, np.NaN)
                        # save experiment bias statistic by group to selected station data dictionary
                        canvas_instance.selected_station_data[networkspeci][data_label][temporal_aggregation_resolution][
                            '{}'.format(bias_stat)] = stat_output_by_group

def calculate_z_statistic(read_instance, z1, z2, zstat, networkspeci):
    """Function that calculates selected statistic across stations for map plots.
    
    :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type read_instance: object
    :param z1: name of first data array to plot
    :type z1: str
    :param z2: name of second data array to plot in case of bias plots (empty str if plot is absolute)
    :type z2: str
    :param zstat: name of statistic
    :type zstat: str
    :param networkspeci: name of networkspeci str
    :type networkspeci: str
    :return: calculated map statistic and active station indices on map
    :rtype: np.float32, np.int
    """

    # check if have valid station data first
    # if not update z statistic and active map valid station indices to be empty lists and return
    if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
        n_valid_stations = len(read_instance.valid_station_inds_temporal_colocation[networkspeci]['observations'])
    else:
        n_valid_stations = len(read_instance.valid_station_inds[networkspeci]['observations'])
    if n_valid_stations == 0:             
        z_statistic = np.array([], dtype=np.float32)
        active_map_valid_station_inds = np.array([], dtype=np.int)
        return z_statistic, active_map_valid_station_inds 

    # get zstat information 
    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)

    # get dictionary containing necessary information for calculation of selected statistic
    if z_statistic_type == 'basic':
        stats_dict = basic_stats[base_zstat]
    else:
        stats_dict = expbias_stats[base_zstat]

    # read selected data arrays (after subsetting arrays by intersection of z1/z2 valid station indices)
    # and calculate desired Z statistic (after removing NaNs from arrays)

    # get active map valid station indices (i.e. the indices of the stations data to plot on the map)
    # if only have z1, valid map indices are those simply for the z1 array
    if z2 == '':
        if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
            active_map_valid_station_inds = read_instance.valid_station_inds_temporal_colocation[networkspeci][z1]
        else:
            active_map_valid_station_inds = read_instance.valid_station_inds[networkspeci][z1]
    else:
        # if have z2 array, get intersection of z1 and z2 valid station indices
        if read_instance.temporal_colocation and len(read_instance.data_labels) > 1:
            active_map_valid_station_inds = \
                np.intersect1d(read_instance.valid_station_inds_temporal_colocation[networkspeci][z1],
                               read_instance.valid_station_inds_temporal_colocation[networkspeci][z2])
        else:
            active_map_valid_station_inds = \
                np.intersect1d(read_instance.valid_station_inds[networkspeci][z1],
                               read_instance.valid_station_inds[networkspeci][z2])

    # read z1 data
    z1_array_data = \
        read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(z1),active_map_valid_station_inds,:]
    # drop NaNs and reshape to object list of station data arrays (if not checking data %)
    if base_zstat != 'Data%':
        z1_array_data = drop_nans(z1_array_data)
    else:
        z1_array_data.tolist()

    # create empty array to store z statistic
    z_statistic = np.empty(len(z1_array_data))

    # if have no z2 data, calculate 'absolute' basic statistic
    if z2 == '':

        # load default selected z statistic arguments for passing to statistical function
        function_arguments = stats_dict['arguments']

        # iterate through stations calculating statistic
        for z_ii in range(len(z_statistic)):

            # and calculate its statistics
            z_statistic[z_ii] = \
                getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments)

    # else, read z2 data then calculate 'difference' statistic
    else:
        # read z2 data
        z2_array_data = \
            read_instance.data_in_memory_filtered[networkspeci][read_instance.data_labels.index(z2),active_map_valid_station_inds,:]
        # drop NaNs and reshape to object list of station data arrays (if not checking data %)
        if base_zstat != 'Data%':
            z2_array_data = drop_nans(z2_array_data)
        else:
            z2_array_data = z2_array_data.tolist()
        # is the difference statistic basic (i.e. mean)?
        if z_statistic_type == 'basic':

            # load default selected z statistic arguments and make separate arguments
            # dictionaries for z1/z2 calculations (as doing 2 separate calculations for z1/z2 and subtracting)
            function_arguments_z1 = stats_dict['arguments']
            function_arguments_z2 = copy.deepcopy(function_arguments_z1)

            # iterate through stations calculating statistic
            for z_ii in range(len(z_statistic)):

                # calculate statistics for z1/z2 arrays and subtract z2-z1
                z_statistic[z_ii] = \
                    getattr(Stats, stats_dict['function'])(z2_array_data[z_ii], **function_arguments_z2) - \
                    getattr(Stats, stats_dict['function'])(z1_array_data[z_ii], **function_arguments_z1)

        # else, is the difference statistic an experiment bias statistic (i.e. r)?
        elif z_statistic_type == 'expbias':

            # load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            # iterate through stations calculating statistic
            for z_ii in range(len(z_statistic)):

                # set station z1/z2 arrays as arguments in argument dictionary
                function_arguments['obs'] = z1_array_data[z_ii]
                function_arguments['exp'] = z2_array_data[z_ii]

                # calculate statistic
                z_statistic[z_ii] = getattr(ExpBias, stats_dict['function'])(**function_arguments)

    # if any station z statistics come out as NaN/inf, remove respective
    # stations from active map valid station indices
    # also cut z_statistic to remove invalid NaNs/infs
    valid_z_statistic_boolean = np.isfinite(z_statistic)
    active_map_valid_station_inds = active_map_valid_station_inds[valid_z_statistic_boolean]
    z_statistic = z_statistic[valid_z_statistic_boolean]

    return z_statistic, active_map_valid_station_inds

def get_axes_vminmax(axs):
    """Function that get minimum and maximum of plotted data across relevant axes

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
    """Function that generates neccessary detail to crate colourbar.

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

    #get zstat information
    zstat, base_zstat, z_statistic_type, z_statistic_sign = get_z_statistic_info(zstat=zstat)

    # get dictionary containing necessary information for calculation of selected statistic
    if z_statistic_type == 'basic':
        stats_dict = basic_stats[base_zstat]
        if base_zstat not in ['Data%','Exceedances']:
            label_units = ' ({})'.format(read_instance.measurement_units[speci])
        else:
            label_units = ''
    else:
        stats_dict = expbias_stats[base_zstat]
        label_units = ''

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
            z_label = '{} {}'.format(stats_dict['label'], label_units)
        else:
            if z_statistic_type == 'basic':
                z_label = '{} bias {}'.format(stats_dict['label'], label_units)
            else:
                z_label = '{} {}'.format(stats_dict['label'], label_units)
    
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
    """Function that generates colourbar.

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
    if plot_characteristics['cb']['discrete']:
        cmap = plt.get_cmap(z_colourmap, plot_characteristics['cb']['n_discrete'])
    else:
        cmap = plt.get_cmap(z_colourmap)

    # create colourbar
    for cb_ax in cb_axs:

        # make colourbar on axis
        cb = matplotlib.colorbar.ColorbarBase(cb_ax, cmap=cmap, norm=norm, 
                                              orientation=plot_characteristics['cb']['orientation'], 
                                              ticks=tick_array)

        # set colorbar label
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
            # remove ticks for discrete colorbars
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
    """Function that gets appropriate zstat name for selected zstatistic comboboxes 

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
    """Function that checks if the z statistic is basic or expbias statistic
    
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
    """Function that checks if the z statistic is an absolute or bias statistic

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
    """Get z statistic name, type (basic or expbias), sign (absolute or bias), base name (dropping '_bias' suffix) from plot_type (or known zstat name)
    
    :param plot_type: plot type
    :type plot_type: str
    :param zstat: name of statistic
    :type plot_type: str
    :return zstat name, base zstat name, zstat type, zstst sign
    :rtype: str, str, str, str
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
            zstat, base_zstat, z_statistic_type, z_statistic_sign = None, None, None, None

    # zstat not 'None'? Then get information for it            
    if zstat:
        # get base name name of zstat (dropping _bias suffix)
        base_zstat = zstat.split('_bias')[0]
        # get zstat type (basic or expbias) 
        z_statistic_type = get_z_statistic_type(base_zstat)
        # get zstat sign (absolute or bias)
        z_statistic_sign = get_z_statistic_sign(zstat, z_statistic_type)

    return zstat, base_zstat, z_statistic_type, z_statistic_sign