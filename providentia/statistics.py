""" Functions for the processing/calculation of statistics and colourbars """

from calendar import monthrange
import copy
import datetime
import json
import os
import sys
import time
import yaml

import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as st

from providentia.auxiliar import CURRENT_PATH, join, get_conversion_factor, get_standard_parameters_by_speci
from .calculate import Stats, ExpBias
from .read_aux import (get_frequency_code, get_chunk_size,
                       get_periodic_nonrelevant_temporal_resolutions, get_periodic_relevant_temporal_resolutions)
from .warnings_prv import show_message


PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
basic_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/basic_stats.yaml')))
expbias_stats = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/experiment_bias_stats.yaml')))
fairmode_settings = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/fairmode.yaml')))


def get_selected_station_data(read_instance, canvas_instance, networkspecies, 
                              station_index=None, data_range_min=None, data_range_max=None, stddev_max=None):
    """ Function that takes full data array and cuts it for selected stations, per network / species, per data label.

        :param read_instance: Instance of class Dashboard or Report
        :type read_instance: object
        :param canvas_instance: Instance of class Canvas or Report
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
        if data_range_min is not None:
            canvas_instance.selected_station_data_min[networkspeci] = data_range_min[networkspeci]
        else:
            canvas_instance.selected_station_data_min[networkspeci] = np.inf
        if data_range_max is not None:
            canvas_instance.selected_station_data_max[networkspeci] = data_range_max[networkspeci]
        else:
            canvas_instance.selected_station_data_max[networkspeci] = 0.0
        if stddev_max is not None:
            canvas_instance.selected_station_stddev_max[networkspeci] = stddev_max[networkspeci]
        else:
            canvas_instance.selected_station_stddev_max[networkspeci] = 0.0      

        # get data array for networkspeci
        data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci])

        # temporally colocate data array
        if read_instance.temporal_colocation:
            data_array[:, read_instance.temporal_colocation_nans[networkspeci]] = np.nan
        
        # get selected station indices
        canvas_instance.station_inds[networkspeci] = get_station_inds(read_instance, canvas_instance, networkspeci, station_index)

        # get data cut for relevant stations
        data_array = data_array[:,canvas_instance.station_inds[networkspeci],:]

        # get NaNs in data array
        nan_data_array = np.isnan(data_array)

        # if data array has no valid data for selected stations, do not cut data array
        # data array has valid data and is not all nan?
        if data_array.size > 0 and not np.all(nan_data_array):

            # set metadata cut for relevant stations
            canvas_instance.selected_station_metadata[networkspeci] = read_instance.metadata_in_memory[networkspeci][canvas_instance.station_inds[networkspeci],:]

            # get which data labels have some valid data
            valid_data_labels_mask = ~np.all(np.all(nan_data_array, axis=-1), axis=-1)
            canvas_instance.selected_station_data_labels[networkspeci] = list(np.array(read_instance.data_labels)[valid_data_labels_mask])

            # cut data array for valid data labels
            data_array = data_array[valid_data_labels_mask]

            # do resampling
            data_array = do_resampling(read_instance, data_array)

            # if have daily forecast active, reshape array to have forecast dimension
            if read_instance.daily_forecast:

                # group data per hour, shape after is (24, label, station, time_per_hour_tiled)
                data_array_forecast_grouped = group_periodic(read_instance, canvas_instance, networkspeci, 'hour', False, '', '', data_array)

                # reshape to add a forecast day dimension, shape after is (24, label, station, time_per_hour, fct)
                data_array_forecast_grouped_rs = np.reshape(data_array_forecast_grouped, 
                                                            (data_array_forecast_grouped.shape[0], data_array_forecast_grouped.shape[1], data_array_forecast_grouped.shape[2], -1, read_instance.max_forecast_days), 
                                                            order='F')

                # aggregate across the time dimension, after this shape is (24, label, station, fct)
                data_array_forecast_agg = aggregation(data_array_forecast_grouped_rs, read_instance.timeseries_statistic_aggregation, read_instance=read_instance, axis=-2)

                # move per hour dimension from first to second last
                data_array_forecast_agg = data_array_forecast_agg.transpose(1, 2, 0, 3)

                # reshape to remove the forecast day dimension, and just have a statistic per hour
                data_array_forecast_agg = np.reshape(data_array_forecast_agg,
                                                    (data_array_forecast_agg.shape[0], data_array_forecast_agg.shape[1], -1), 
                                                    order='F') 

                # set data array for timeseries
                data_array_ts = data_array_forecast_agg

                #set time_index
                # determine steps per day based on resolution
                if read_instance.active_resolution == 'hourly':
                    steps_per_day = 24
                elif read_instance.active_resolution == '3hourly':
                    steps_per_day = 8    # 24 / 3
                elif read_instance.active_resolution == '6hourly':
                    steps_per_day = 4    # 24 / 6
                elif read_instance.active_resolution == 'daily':
                    steps_per_day = 1

                # take the first 'steps_per_day' timestamps as the base block
                base_block = read_instance.time_index[:steps_per_day]

                # shift the block for each active forecast day
                time_index_blocks = [base_block + pd.Timedelta(days=(day-1)) for day in read_instance.active_forecast_days]

                # concatenate all blocks
                time_index = pd.DatetimeIndex(np.concatenate(time_index_blocks))

            # if have all combined forecast active, reshape array to have forecast dimension
            elif read_instance.combined_forecast:

                # reshape to add a forecast day dimension
                data_array_rs = np.reshape(data_array, 
                                           (data_array.shape[0], data_array.shape[1], -1, read_instance.max_forecast_days), 
                                           order='F')

                # aggregate across the forecast day dimension, after this shape is (label, station, time)
                data_array_forecast_agg = aggregation(data_array_rs, read_instance.timeseries_statistic_aggregation, read_instance=read_instance, axis=-1)

                # set data array for timeseries
                data_array_ts = data_array_forecast_agg

                # set time index, and set to self
                time_index = read_instance.time_index[:data_array_forecast_agg.shape[-1]]
                read_instance.time_index = time_index

            # otherwise, set data array for timeseries and time_index
            else:
                data_array_ts = data_array
                time_index = read_instance.time_index

            # save timeseries array
            if len(canvas_instance.station_inds[networkspeci]) == 1:
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = data_array_ts[:,0,:]
            else:
                canvas_instance.selected_station_data[networkspeci]['timeseries'] = aggregation(data_array_ts, read_instance.timeseries_statistic_aggregation, read_instance=read_instance, axis=1)

            # save data per station
            if read_instance.statistic_mode == 'Spatial|Temporal':
                # if have daily forecast active then do not want to use tiled timeseries array for per station data
                # rather we want the standard aggregated time series across the data period 
                if read_instance.daily_forecast:
                    # do aggregation across stations
                    if len(canvas_instance.station_inds[networkspeci]) == 1:
                        spatial_agg = data_array[:,0,:]
                    else:
                        spatial_agg = aggregation(data_array, read_instance.statistic_aggregation, read_instance=read_instance, axis=1)
                    canvas_instance.selected_station_data[networkspeci]['per_station'] = spatial_agg[:,np.newaxis,:]

                # non-daily forecast cases
                else:
                    # if statistic aggregation is the same as the timeseries statistic aggregation then can take the timeseries
                    if read_instance.statistic_aggregation == read_instance.timeseries_statistic_aggregation:
                        canvas_instance.selected_station_data[networkspeci]['per_station'] = canvas_instance.selected_station_data[networkspeci]['timeseries'][:,np.newaxis,:]
                    # otherwise do aggregation
                    else:
                        aggregated_data = aggregation(data_array_ts, read_instance.statistic_aggregation, read_instance=read_instance, axis=1)
                        canvas_instance.selected_station_data[networkspeci]['per_station'] = aggregated_data[:,np.newaxis,:]

            elif read_instance.statistic_mode in ['Temporal|Spatial', 'Flattened']:
                canvas_instance.selected_station_data[networkspeci]['per_station'] = data_array

            # transform timeseries to pandas dataframe
            canvas_instance.selected_station_data[networkspeci]['timeseries'] = pd.DataFrame(canvas_instance.selected_station_data[networkspeci]['timeseries'].T, 
                                                                                             columns=canvas_instance.selected_station_data_labels[networkspeci], 
                                                                                             index=time_index)

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

def do_resampling(read_instance, data_array, writing=False):

    """ Function which handles resampling of data """

    # if resampling resolution is None, then do not resample
    if read_instance.resampling_resolution == 'None':
        do_resampling = False
    else:
        do_resampling = True

    # update relevant/nonrelevant periodic temporal resolutions 
    if not writing:
        read_instance.periodic_relevant_temporal_resolutions = get_periodic_relevant_temporal_resolutions(read_instance.active_resolution)    
        read_instance.periodic_nonrelevant_temporal_resolutions = get_periodic_nonrelevant_temporal_resolutions(read_instance.active_resolution) 

    # need to do resampling
    if do_resampling:

        # transform resolution to code for .resample function
        temporal_resolution_to_output_code = get_frequency_code(read_instance.resampling_resolution)       

        # flatten data label dimension for creation of pandas dataframe
        data_array_reduced = data_array.reshape(data_array.shape[0]*data_array.shape[1], data_array.shape[2])

        # determine number of forecast days (if daily or combined forecast is active)
        if (read_instance.daily_forecast) or (read_instance.combined_forecast):
            forecast_days = read_instance.max_forecast_days
        # otherwise set it as 1
        else:
            forecast_days = 1

        # set n base times
        n_base_times = len(read_instance.time_array)
        resampled_dfs = []

        # iterate through each forecast day and resample data
        for day in range(forecast_days):
            # slice the portion of data for this forecast day
            start = day * n_base_times
            end = (day + 1) * n_base_times
            
            # create DataFrame for this forecast day
            df_day = pd.DataFrame(
                data_array_reduced[:, start:end].T,  # transpose so time is rows
                index=pd.DatetimeIndex(read_instance.time_array),
                columns=np.arange(data_array_reduced.shape[0]),
                dtype=np.float32
            )
            
            # resample for this forecast day
            df_day_resampled = df_day.resample(temporal_resolution_to_output_code).mean()
            
            # append df if have more than 1 forecast day
            if forecast_days > 1:
                resampled_dfs.append(df_day_resampled)

        # concatenate resampled dfs across each forecast day
        if forecast_days > 1:
            data_array_df_resampled = pd.concat(resampled_dfs)
        # otherwise just take only resampled df
        else:
            data_array_df_resampled = df_day_resampled

        # save back out as numpy array (reshaping to get back networkspecies dimension)
        data_array_resampled = data_array_df_resampled.to_numpy().transpose()
        data_array_resampled = data_array_resampled.reshape(data_array.shape[0], data_array.shape[1], data_array_resampled.shape[1])

        # update time index and return
        if writing:
            return data_array_resampled, data_array_df_resampled.index
        else:
            read_instance.time_index = data_array_df_resampled.index
            return data_array_resampled

    # resampling not neccessary?
    else:
        if writing:
            return data_array, read_instance.time_index_after_filter
        else:
            read_instance.time_index = read_instance.time_index_after_filter
            return data_array

def merge_forecast_days(read_instance, networkspeci, data_labels, unique_base_data_labels, data_array):
    """
    Function which joins different forecast days separated as different expertiments as 1 tiled experiment.
    Observations are repeatedly tiled to macth the tiled experiment shape.
    Non-forecst experiments also change to match the shape but only the first forecast day is filled.
    """

    # get n_labels and n_stations of data array
    n_labels, n_stations, _ = data_array.shape

    # set block size (len of original time array)
    block_size = len(read_instance.time_array)

    # Create a new array filled with NaNs to store updated forecast data
    new_data_in_memory = np.full(
        (len(unique_base_data_labels), n_stations, read_instance.max_forecast_days * block_size),
        np.nan, dtype=np.float32
    )

    # Iterate over each unique base label
    for base_label_ii, base_label in enumerate(unique_base_data_labels):
        # Find indices of all data_labels that match this base label
        relevant_inds = np.array([i for i, lbl in enumerate(data_labels) if base_label in lbl], dtype=np.int32)

        if len(relevant_inds) == 0:
            continue  # Skip if no matching labels

        for j, ind in enumerate(relevant_inds):
            data_block = data_array[ind, :, :]  # Extract data for this label
            if base_label == read_instance.observations_data_label:
                # Observations are repeated across all forecast days
                for day in range(read_instance.max_forecast_days):
                    new_data_in_memory[base_label_ii, :, day*block_size:(day+1)*block_size] = data_block
            else:
                # For forecast data, place each day in its corresponding block
                new_data_in_memory[base_label_ii, :, j*block_size:(j+1)*block_size] = data_block
                # If there are no forecast indices and this is the first block, stop
                if (len(read_instance.forecast_indices_per_data_label[networkspeci][base_label]) == 0) & (j == 0):
                    break

    return new_data_in_memory

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
        p25 = np.nanpercentile(data, 25, method='nearest')
        #calculate the 75th percentile 
        p75 = np.nanpercentile(data, 75, method='nearest')

        #calculate the interquartile range
        iqr = p75-p25

        #calculate lower/upper inner fences and return values
        lower_inner_fence = p25 - 1.5 * iqr 
        upper_inner_fence = p75 + 1.5 * iqr

        return lower_inner_fence, upper_inner_fence

def get_station_inds(read_instance, canvas_instance, networkspeci, station_index):
    """ Get selected station indices """
        
    if station_index is not None:
        station_inds = np.array([station_index])
    else:
        if read_instance.mode in ['report', 'library']:
            if read_instance.temporal_colocation:
                station_inds = read_instance.valid_station_inds_temporal_colocation[networkspeci][read_instance.observations_data_label]
            else:
                station_inds = read_instance.valid_station_inds[networkspeci][read_instance.observations_data_label] 
        else:
            station_inds = canvas_instance.relative_selected_station_inds

    return station_inds

def group_periodic(read_instance, canvas_instance, networkspeci, period_resolution, 
                   per_station, statistic_mode, base_zstat, data_array,
                   return_nan_padding_counts=False):
    """
    Function that groups data into periodic chunks.

    Input:  (label, station, time)
    Output: (chunk, label, station, chunk_time)

    Periodic grouping is done for hours of the day, days of the week, or months of the year.
    Applies NaN padding where chunk lengths differ, and tracks padding counts.
    """

    # ------------------------------------------------------------------
    # Get all periods for current resolution
    # ------------------------------------------------------------------
    # For example, for 'hour' resolution: all_periods = array of hours 0-23 repeated for all days
    all_periods = getattr(read_instance.time_index, period_resolution)

    # Extract unique periods (e.g., unique hours, weekdays, or months)
    canvas_instance.unique_periods = np.unique(all_periods)

    # ------------------------------------------------------------------
    # Iterate through each unique period and collect data
    # ------------------------------------------------------------------
    periodic_data = []   # list to store data arrays per period
    max_chunk_len = 0    # track maximum length of chunks (for NaN padding)

    for periodic_xtick in canvas_instance.unique_periods:
        # Create boolean mask for current period
        valid_period = all_periods == periodic_xtick

        # Extract data corresponding to current period
        period_data = data_array[:, :, valid_period]

        # Update max_chunk_len to know how much padding is needed
        chunk_len = period_data.shape[-1]
        if chunk_len > max_chunk_len:
            max_chunk_len = chunk_len

        periodic_data.append(period_data)

    # ------------------------------------------------------------------
    # Pad each period to max_chunk_len with NaNs
    # ------------------------------------------------------------------
    nan_padding_counts = []
    padded_groups = []
    for group in periodic_data:
        pad_len = max_chunk_len - group.shape[-1]  # number of NaNs to add
        if pad_len > 0:
            pad_shape = list(group.shape[:-1]) + [pad_len]  # shape for padding
            group = np.concatenate([group, np.full(pad_shape, np.nan)], axis=-1)
        padded_groups.append(group)
        nan_padding_counts.append(pad_len)

    # Stack into single numpy array: (chunk, label, station, chunk_time)
    periodic_data = np.stack(padded_groups, axis=0)
    nan_padding_counts = np.array(nan_padding_counts, dtype=np.int32)

    # ------------------------------------------------------------------
    # Flattened mode
    # ------------------------------------------------------------------
    if (statistic_mode == 'Flattened') and (not per_station) and (base_zstat not in ['NStations','MDA8']):

        # grouped_data shape before flattening: (n_chunks, n_labels, n_stations, chunk_time)
        n_chunks, n_labels, n_stations, chunk_time = periodic_data.shape

        # Flatten stations and chunk_time into single axis
        periodic_data = periodic_data.reshape(n_chunks, n_labels, 1, n_stations * chunk_time)

        # Scale nan_padding_counts to match the flattened last axis
        nan_padding_counts = nan_padding_counts * n_stations  # total padded NaNs per chunk

        # Reshape for broadcasting
        nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis]  # shape (n_chunks, 1, 1)

    else:
        # Non-flattened mode: reshape for broadcasting
        nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis]

    # ------------------------------------------------------------------
    # Return data
    # ------------------------------------------------------------------
    if return_nan_padding_counts:
        return periodic_data, nan_padding_counts
    return periodic_data

def group_temporal(read_instance, canvas_instance, networkspeci, chunk_resolution, 
                   per_station, statistic_mode, base_zstat, data_array,
                   prev_nan_padding_counts=None, return_nan_padding_counts=False):
    """
    Function that groups data into temporal chunks.

    Handles both non-forecast data (3D -> 4D) and forecast data (5D -> 5D).
    Applies NaN padding where chunk lengths differ, and tracks padding counts.
    """

    # get existing timeseries data 
    timeseries_data = canvas_instance.selected_station_data[networkspeci]['timeseries']

    # get frequency code for new resolution
    new_freq = get_frequency_code(chunk_resolution)
    
    # get new resampled timeseries indices
    if chunk_resolution == "monthly":
        canvas_instance.grouped_ts_index = timeseries_data.index.to_period("M").to_timestamp().unique().sort_values()
    elif chunk_resolution == "annual":
        canvas_instance.grouped_ts_index = timeseries_data.index.to_period("Y").to_timestamp().unique().sort_values()
    else:
        canvas_instance.grouped_ts_index = timeseries_data.index.floor(new_freq).unique().sort_values()

    # ------------------------------------------------------------------
    # 3D NON-FORECAST CASE
    # ------------------------------------------------------------------
    if data_array.ndim == 3:
        n_labels, n_stations, n_times = data_array.shape

        grouped_data = []
        max_chunk_len = 0  # track maximum length of chunks (for NaN padding)
        nan_padding_counts = []

        for index in canvas_instance.grouped_ts_index:
            # Determine start and end date for current chunk based on resolution
            if new_freq == "h":
                start_date = end_date = index
            elif new_freq == "3h":
                start_date = index
                end_date = index + datetime.timedelta(hours=2)
            elif new_freq == "6h":
                start_date = index
                end_date = index + datetime.timedelta(hours=5)
            elif new_freq == "D":
                start_date = datetime.datetime(index.year, index.month, index.day, 0)
                end_date   = datetime.datetime(index.year, index.month, index.day, 23)
            elif new_freq == "MS":
                start_date = datetime.datetime(index.year, index.month, 1, 0)
                end_day = monthrange(index.year, index.month)[1]
                end_date = datetime.datetime(index.year, index.month, end_day, 23)
            elif new_freq == "YS":
                start_date = datetime.datetime(index.year, 1, 1, 0)
                end_date   = datetime.datetime(index.year, 12, 31, 23)

            # Get the indices corresponding to this chunk
            time_indices = timeseries_data.index.get_indexer(timeseries_data[start_date:end_date].index)
            cut_data = np.take(data_array, time_indices, axis=-1)

            # Update max_chunk_len to know how much NaN padding is needed
            chunk_len = cut_data.shape[-1]
            if chunk_len > max_chunk_len:
                max_chunk_len = chunk_len

            grouped_data.append(cut_data)

        # Pad each chunk to max_chunk_len with NaNs if needed
        padded_groups = []
        for group in grouped_data:
            pad_len = max_chunk_len - group.shape[-1]
            if pad_len > 0:
                pad_shape = list(group.shape[:-1]) + [pad_len]
                group = np.concatenate([group, np.full(pad_shape, np.nan)], axis=-1)
            padded_groups.append(group)
            nan_padding_counts.append(pad_len)

        # Stack into single array with shape: (chunk, label, station, chunk_time)
        grouped_data = np.stack(padded_groups, axis=0)
        nan_padding_counts = np.array(nan_padding_counts, dtype=np.float32)

    # ------------------------------------------------------------------
    # 5D FORECAST CASE
    # ------------------------------------------------------------------
    elif data_array.ndim == 5:
        # Monthly/annual not allowed for forecast data
        if chunk_resolution in ["monthly", "annual"]:
            raise ValueError("Cannot group forecast data to monthly or annual resolution.")

        #get chunk size
        chunk_size = get_chunk_size(read_instance.active_resolution, chunk_resolution)

        n_chunks, n_labels, n_stations, n_forecast, n_times_in_chunk = data_array.shape
        n_out_chunks = int(np.ceil(n_chunks / chunk_size))

        # Pad the data array if the number of chunks is not divisible by chunk_size
        pad_chunks = n_out_chunks * chunk_size - n_chunks
        if pad_chunks > 0:
            pad_shape = (pad_chunks, n_labels, n_stations, n_forecast, n_times_in_chunk)
            pad_array = np.full(pad_shape, np.nan)
            data_array = np.concatenate([data_array, pad_array], axis=0)

        # Reshape and transpose to regroup chunks
        grouped_data = data_array.reshape(
            n_out_chunks, chunk_size, n_labels, n_stations, n_forecast, n_times_in_chunk
        )
        grouped_data = grouped_data.transpose(0, 2, 3, 4, 1, 5)
        grouped_data = grouped_data.reshape(
            n_out_chunks, n_labels, n_stations, n_forecast, chunk_size * n_times_in_chunk
        )

        # ----- NaN padding counts -----
        if prev_nan_padding_counts is not None:
            prev_nan_padding_counts = np.squeeze(prev_nan_padding_counts)
            if pad_chunks > 0:
                prev_nan_padding_counts = np.concatenate([
                    prev_nan_padding_counts,
                    np.zeros(pad_chunks, dtype=np.float32)
                ])
            prev_nan_padding_counts = prev_nan_padding_counts.reshape(n_out_chunks, chunk_size)
            nan_padding_counts = prev_nan_padding_counts.sum(axis=1)
        else:
            nan_padding_counts = np.zeros(n_out_chunks, dtype=np.float32)

        if pad_chunks > 0:
            # add NaNs for padded chunks
            nan_padding_counts[-1] += (
                pad_chunks * n_labels * n_stations * n_forecast * n_times_in_chunk
            )

    else:
        raise ValueError("Input array must be 3D or 5D.")

    # ------------------------------------------------------------------
    # Flattened mode
    # ------------------------------------------------------------------
    if (statistic_mode == 'Flattened') and (not per_station) and (base_zstat not in ['NStations','MDA8']):
        if grouped_data.ndim == 4:
            # grouped_data shape: (n_chunks, n_labels, n_stations, chunk_time)
            n_chunks, n_labels, n_stations, chunk_time = grouped_data.shape

            # flatten stations and chunk_time into one axis
            grouped_data = grouped_data.reshape(n_chunks, n_labels, 1, n_stations * chunk_time)

            # scale nan_padding_counts to match the new flattened last axis
            nan_padding_counts = nan_padding_counts * n_stations

            # reshape for broadcasting
            nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis]  # (n_chunks, 1, 1)

        elif grouped_data.ndim == 5:
            # grouped_data shape: (n_chunks, n_labels, n_stations, n_forecast, chunk_time)
            n_chunks, n_labels, n_stations, n_forecast, chunk_time = grouped_data.shape

            # reorder and flatten stations and chunk_time
            grouped_data = grouped_data.transpose(0, 1, 3, 2, 4)
            grouped_data = grouped_data.reshape(n_chunks, n_labels, 1, n_forecast, n_stations * chunk_time)

            # scale nan_padding_counts to match flattened last axis
            nan_padding_counts = nan_padding_counts * n_stations
            nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis, np.newaxis]

    else:
        # Non-flattened mode: reshape nan_padding_counts for broadcasting
        if grouped_data.ndim == 4:
            nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis]

        elif grouped_data.ndim == 5:
            nan_padding_counts = nan_padding_counts[:, np.newaxis, np.newaxis, np.newaxis]

    if return_nan_padding_counts:
        return grouped_data, nan_padding_counts

    return grouped_data



def calculate_statistic(read_instance, canvas_instance, networkspeci, zstats, data_labels_a, 
                        data_labels_b, map=False, per_station=False, period=None, chunk_resolution=None, 
                        reduction=True, mask=None, statistic_mode=None, statistic_aggregation=None, 
                        periodic_statistic_mode=None, periodic_statistic_aggregation=None,
                        forecast_type=None):
    """Function that calculates a statistic for data labels, either absolute or bias, 
       for different aggregation modes.
    """

    # if statistic mode is None, then take the global one
    if not statistic_mode:
        statistic_mode = read_instance.statistic_mode 

    # if statistic aggregation is None, then take the global one
    if not statistic_aggregation:
        statistic_aggregation = read_instance.statistic_aggregation 

    # if periodic statistic mode is None, then take the global one
    if not periodic_statistic_mode:
        periodic_statistic_mode = read_instance.periodic_statistic_mode

    # if periodic statistic aggregation is None, then take the global one
    if not periodic_statistic_aggregation:
        periodic_statistic_aggregation = read_instance.periodic_statistic_aggregation 

    # if data_labels_a, data_labels_b are strings then convert to lists
    if type(data_labels_a) != list:
        data_labels_a = [data_labels_a]
    if type(data_labels_b) != list:
        data_labels_b = [data_labels_b]

    # if have empty strings in lists then remove them
    data_labels_a = [label for label in data_labels_a if label != '']
    data_labels_b = [label for label in data_labels_b if label != '']

    # if zstats is str then make it a list
    if type(zstats) != list:
        zstats = [zstats]

    # iterate through zstats and calculate statistics
    stats_calc = {}
    for zstat in zstats:

        #print()

        # initialise nan_padding_counts as None (used to track padded NaNs in grouped arrays)
        nan_padding_counts_a = None
        nan_padding_counts_b = None

        # get zstat information 
        zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        if z_statistic_period is not None:
            pp = z_statistic_period
        elif period is not None:
            pp = period

        #if (z_statistic_period is not None) & (period is not None) & (chunk_resolution is not None):
        #    print('Stat:{}, Mode:{}, zstat period:{}, period:{}, resolution:{}'.format(zstat, statistic_mode, z_statistic_period, period, chunk_resolution))
        #elif (z_statistic_period is not None) & (period is not None):
        #    print('Stat:{}, Mode:{}, zstat period:{}, period:{}'.format(zstat, statistic_mode, z_statistic_period, period))
        #elif (z_statistic_period is not None) & (chunk_resolution is not None):
        #    print('Stat:{}, Mode:{}, zstat period:{}, resolution:{}'.format(zstat, statistic_mode, z_statistic_period, chunk_resolution))
        #elif (period is not None) & (chunk_resolution is not None):
        #    print('Stat:{}, Mode:{}, period:{}, resolution:{}'.format(zstat, statistic_mode, period, chunk_resolution))
        #elif (z_statistic_period is not None):
        #    print('Stat:{}, Mode:{}, zstat period:{}'.format(zstat, statistic_mode, z_statistic_period))
        #elif (period is not None):
        #    print('Stat:{}, Mode:{}, period:{}'.format(zstat, statistic_mode, period))
        #elif (chunk_resolution is not None):
        #    print('Stat:{}, Mode:{}, resolution:{}'.format(zstat, statistic_mode, chunk_resolution))
        #else:
        #    print('Stat:{}, Mode:{}'.format(zstat, statistic_mode))

        # for map statistics, get active map valid station indices and then data_labels_a data 
        if (map) or (per_station):

            # get relevant station indices
            if per_station:
                station_inds = canvas_instance.station_inds[networkspeci]

            elif map:
                if read_instance.temporal_colocation:
                    inds = read_instance.valid_station_inds_temporal_colocation[networkspeci]
                else:
                    inds = read_instance.valid_station_inds[networkspeci]
                # if have just data_labels_a, station indices are simply those relevant for data_labels_a
                if len(data_labels_b) == 0:
                    if data_labels_a[0] not in inds:
                        error = f'Data label {data_labels_a[0]} not in array. Options: {list(inds.keys())}'
                        read_instance.logger.error(error)
                        sys.exit(1)
                    station_inds = inds[data_labels_a[0]]
                # elif have data_labels_b, get intersection of data_labels_a and data_labels_b valid station indices
                else:
                    if data_labels_a[0] not in inds:
                        error = f'Data label {data_labels_a[0]} not in array. Options: {list(inds.keys())}'
                        read_instance.logger.error(error)
                        sys.exit(1)
                    if data_labels_b[0] not in inds:
                        error = f'Data label {data_labels_b[0]} not in array. Options: {list(inds.keys())}'
                        read_instance.logger.error(error)
                        sys.exit(1)
                    station_inds = np.intersect1d(inds[data_labels_a[0]], inds[data_labels_b[0]])

            # check if valid station data
            # if not return
            if len(station_inds) == 0:
                z_statistic = np.array([], dtype=np.float32)
                if map:
                    station_inds = np.array([], dtype=np.int32)
                    return z_statistic, station_inds 
                elif per_station:
                    return z_statistic
                
            # get data array_a
            data_label_a_indices = np.array([read_instance.data_labels.index(label) for label in data_labels_a], dtype=np.int32)
            data_array_a = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][data_label_a_indices,:,:])

            # temporally colocate data (if active)
            if read_instance.temporal_colocation:
                data_array_a[:, read_instance.temporal_colocation_nans[networkspeci]] = np.nan
                
            # get data cut for relevant stations
            data_array_a = data_array_a[:,station_inds,:]

        # for other cases, get cut of selected station data for data_labels_a
        else:
            data_label_a_indices = np.array([canvas_instance.selected_station_data_labels[networkspeci].index(label) for label in data_labels_a], dtype=np.int32)
                        
            # for grouping data take per_station array
            if (chunk_resolution is not None) or (period is not None) or (z_statistic_period is not None) or (base_zstat in ['NStations','MDA8']):
                data_array_a = canvas_instance.selected_station_data[networkspeci]['per_station'][data_label_a_indices]
            # otherwise take active mode
            else:
                data_array_a = canvas_instance.selected_station_data[networkspeci]['active_mode'][data_label_a_indices]

        #print('Data Array a: ', data_array_a.shape)

        # if need to mask data, then do so
        if mask is not None:
            data_array_a[mask[data_label_a_indices]] = np.nan

        # if need to reshape forecast data, then do so
        if forecast_type == 'daily':

            # group data per hour, shape after is (24, label, station, time_per_hour_tiled)
            data_array_a, nan_padding_counts_a = group_periodic(read_instance, canvas_instance, networkspeci, 'hour', False, '', '', data_array_a, return_nan_padding_counts=True)
                
            #print('Data Array a, forecast daily group periodic: ', data_array_a.shape)

            # reshape to add a forecast day dimension, shape after is (24, label, station, fct, time_per_hour)
            data_array_a = np.reshape(data_array_a, 
                                      (data_array_a.shape[0], data_array_a.shape[1], data_array_a.shape[2], read_instance.max_forecast_days, -1), 
                                      order='F')

            #print('Data Array a, reshape forecast period: ', data_array_a.shape)

        # if need to temporally chunk data, then do so
        if chunk_resolution is not None:
            if ((chunk_resolution == read_instance.active_resolution) & (statistic_mode in ['Temporal|Spatial', 'Spatial|Temporal']) & (forecast_type != 'daily')) or ((chunk_resolution == read_instance.active_resolution) & (statistic_mode == 'Flattened') & (base_zstat in ['NStations','MDA8']) & (forecast_type != 'daily')):
                data_array_a = np.expand_dims(np.transpose(data_array_a, (2,0,1)), -1)
                if len(data_labels_b) == 0:
                    chunk_resolution = None
            else:
                data_array_a, nan_padding_counts_a = group_temporal(read_instance, canvas_instance, networkspeci, chunk_resolution, per_station, statistic_mode, base_zstat, data_array_a, prev_nan_padding_counts=nan_padding_counts_a, return_nan_padding_counts=True)            
            #print('Data Array a, group temporal: ', data_array_a.shape)

        # if need to group data for a period, then do so
        # this can be for calculating statistics per period, or an integated periodic statistic
        elif period is not None:
            data_array_a, nan_padding_counts_a = group_periodic(read_instance, canvas_instance, networkspeci, period, per_station, statistic_mode, base_zstat, data_array_a, return_nan_padding_counts=True)
            #print('Data Array a, group periodic: ', data_array_a.shape)

        elif z_statistic_period is not None:
            data_array_a, nan_padding_counts_a = group_periodic(read_instance, canvas_instance, networkspeci, z_statistic_period, per_station, statistic_mode, base_zstat, data_array_a, return_nan_padding_counts=True)
            #print('Data Array a, group periodic: ', data_array_a.shape)

        # get dictionary containing necessary information for calculation of selected statistic
        if z_statistic_type == 'basic':
            stats_dict = copy.deepcopy(basic_stats[base_zstat])
        else:
            stats_dict = copy.deepcopy(expbias_stats[base_zstat])

        # if have no data_labels_b, calculate 'absolute' basic statistic
        if len(data_labels_b) == 0:

            # load default selected z statistic arguments for passing to statistical function
            function_arguments = stats_dict['arguments']

            # if stat is exceedances then add threshold value (if available)  
            if base_zstat == 'Exceedances':
                function_arguments['threshold'] = exceedance_lim(read_instance, networkspeci)

            # need to do the aggregation inside function for the calculation of NStations and MDA8
            # this is due to handling excepetions in how these are calculated across modes
            elif base_zstat in ['NStations','MDA8']:
                function_arguments['statistic_mode'] = statistic_mode
                function_arguments['statistic_aggregation'] = statistic_aggregation
                function_arguments['per_station'] = per_station
                if z_statistic_period is not None:
                    function_arguments['periodic_statistic_mode'] = periodic_statistic_mode
                    function_arguments['periodic_statistic_aggregation'] = periodic_statistic_aggregation

            # add argument to correct caculation of Data%, when using groups because of padded NaNs
            elif (base_zstat == 'Data%') & (nan_padding_counts_a is not None):
                function_arguments['nan_padding_counts'] = nan_padding_counts_a

            # calculate statistics
            
            # calculate periodic statistic per station
            if z_statistic_period is not None:
                # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                if periodic_statistic_mode == 'Cycle':
                    # aggregation in each group, per station, by periodic statistic
                    z_statistic = aggregation(data_array_a, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()

                    #print('Calculating Stat, Cycle Aggregation ', z_statistic.shape)

                    # need to reshape nan_padding_counts if set as argument
                    if 'nan_padding_counts' in function_arguments:
                        # sum over the chunk axis (first axis) to get total padded NaNs per label/station,
                        # then broadcast to match z_statistic non-time axes
                        function_arguments['nan_padding_counts'] = np.broadcast_to(function_arguments['nan_padding_counts'].sum(axis=0), z_statistic.shape[:-1])

                    # calculate statistic per station (removing period dimension)
                    z_statistic = getattr(Stats, stats_dict['function'])(z_statistic, **function_arguments).transpose()

                # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                # and then aggregate 
                elif periodic_statistic_mode == 'Independent':
                    # calculate statistic per periodic grouping per station
                    z_statistic = getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments).transpose()
                    
                    #print('Calculating Stat, Independent Aggregation ', z_statistic.shape)

                    # aggregate data per station (removing period dimension)
                    if base_zstat != 'NStations':
                        z_statistic = aggregation(z_statistic, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()

            # calculate statistics per station 
            else:
                z_statistic = getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments)
            
            #print('Calculated Stat: ', z_statistic.shape)

        # else, get data_labels_b data then calculate 'difference' statistic
        else:

            # get data_labels_b data for map
            if (map) or (per_station):

                # get data array_b
                data_label_b_indices = np.array([read_instance.data_labels.index(label) for label in data_labels_b], dtype=np.int32)
                data_array_b = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci][data_label_b_indices,:,:])

                # temporally colocate data (if active)
                if read_instance.temporal_colocation:
                    data_array_b[:, read_instance.temporal_colocation_nans[networkspeci]] = np.nan
                    
                # get data cut for relevant stations
                data_array_b = data_array_b[:,station_inds,:]

            # for other cases, get cut of selected station data for data_labels_b
            else:
                data_label_b_indices = np.array([canvas_instance.selected_station_data_labels[networkspeci].index(label) for label in data_labels_b], dtype=np.int32)

                # for grouping data take per_station array
                if (chunk_resolution is not None) or (period is not None) or (z_statistic_period is not None) or (base_zstat in ['NStations','MDA8']):
                    data_array_b = canvas_instance.selected_station_data[networkspeci]['per_station'][data_label_b_indices]
                # otherwise take active mode
                else:
                    data_array_b = canvas_instance.selected_station_data[networkspeci]['active_mode'][data_label_b_indices]
            # if need to mask data, then do so
            if mask is not None:
                data_array_b[mask[data_label_b_indices]] = np.nan

            # if need to reshape forecast data, then do so
            if forecast_type == 'daily':

                # group data per hour, shape after is (24, label, station, time_per_hour_tiled)
                data_array_b, nan_padding_counts_b = group_periodic(read_instance, canvas_instance, networkspeci, 'hour', False, '', '', data_array_b, return_nan_padding_counts=True)
                    
                # reshape to add a forecast day dimension, shape after is (24, label, station, fct, time_per_hour)
                data_array_b = np.reshape(data_array_b, 
                                          (data_array_b.shape[0], data_array_b.shape[1], data_array_b.shape[2], read_instance.max_forecast_days, -1), 
                                          order='F')

            # if need to temporally chunk data, then do so
            if chunk_resolution is not None:
                if ((chunk_resolution == read_instance.active_resolution) & (statistic_mode in ['Temporal|Spatial', 'Spatial|Temporal']) & (forecast_type != 'daily')) or ((chunk_resolution == read_instance.active_resolution) & (statistic_mode == 'Flattened') & (base_zstat in ['NStations','MDA8']) & (forecast_type != 'daily')):
                    data_array_b = np.expand_dims(np.transpose(data_array_b, (2,0,1)), -1)
                    chunk_resolution = None
                else:
                    data_array_b, nan_padding_counts_b = group_temporal(read_instance, canvas_instance, networkspeci, chunk_resolution, per_station, statistic_mode, base_zstat, data_array_b, prev_nan_padding_counts=nan_padding_counts_b, return_nan_padding_counts=True)

            # if need to group data for a period, then do so
            # this can be for calculating statistics per period, or an integated periodic statistic
            elif period is not None:
                data_array_b, nan_padding_counts_b = group_periodic(read_instance, canvas_instance, networkspeci, period, per_station, statistic_mode, base_zstat, data_array_b, return_nan_padding_counts=True)

            elif z_statistic_period is not None:
                data_array_b, nan_padding_counts_b = group_periodic(read_instance, canvas_instance, networkspeci, z_statistic_period, per_station, statistic_mode, base_zstat, data_array_b, return_nan_padding_counts=True)

            # is the difference statistic basic (i.e. mean)?
            if z_statistic_type == 'basic':

                # load default selected statistic arguments and make separate arguments
                # dictionaries for data_labels_a/data_labels_b calculations (as doing 2 separate calculations for data_labels_a/data_labels_b and subtracting)
                function_arguments_a = stats_dict['arguments']

                # if stat is exceedances then add threshold value (if available)  
                if base_zstat == 'Exceedances':
                    function_arguments_a['threshold'] = exceedance_lim(read_instance, networkspeci)

                # need to do the aggregation inside function for the calculation of NStations and MDA8
                # this is due to handling excepetions in how these are calculated across modes
                elif base_zstat in ['NStations','MDA8']:
                    function_arguments_a['statistic_mode'] = statistic_mode
                    function_arguments_a['statistic_aggregation'] = statistic_aggregation
                    function_arguments_a['per_station'] = per_station
                    if z_statistic_period is not None:
                        function_arguments_a['periodic_statistic_mode'] = periodic_statistic_mode
                        function_arguments_a['periodic_statistic_aggregation'] = periodic_statistic_aggregation

                function_arguments_b = copy.deepcopy(function_arguments_a)

                # add argument to correct caculation of Data%, when using groups because of padded NaNs
                if (base_zstat == 'Data%') & (nan_padding_counts_a is not None):
                    function_arguments_a['nan_padding_counts'] = nan_padding_counts_a
                if (base_zstat == 'Data%') & (nan_padding_counts_b is not None):
                    function_arguments_b['nan_padding_counts'] = nan_padding_counts_b

                # calculate statistics for data_labels_a and data_labels_b, then subtract data_labels_b - data_labels_a
                                
                # calculate periodic statistic per station
                if z_statistic_period is not None:
                    # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                    if periodic_statistic_mode == 'Cycle':
                        # aggregation in each group, per station, by periodic statistic
                        statistic_a = aggregation(data_array_a, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()
                        statistic_b = aggregation(data_array_b, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()
                        
                        # need to reshape nan_padding_counts if set as argument
                        if 'nan_padding_counts' in function_arguments_a:
                            # sum over the chunk axis (first axis) to get total padded NaNs per label/station,
                            # then broadcast to match z_statistic non-time axes
                            function_arguments_a['nan_padding_counts'] = np.broadcast_to(function_arguments_a['nan_padding_counts'].sum(axis=0), statistic_a.shape[:-1])
                            function_arguments_b['nan_padding_counts'] = np.broadcast_to(function_arguments_b['nan_padding_counts'].sum(axis=0), statistic_b.shape[:-1])

                        # calculate statistic per station (removing period dimension)
                        statistic_a = getattr(Stats, stats_dict['function'])(statistic_a, **function_arguments_a).transpose()
                        statistic_b = getattr(Stats, stats_dict['function'])(statistic_b, **function_arguments_b).transpose()

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        statistic_a = getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments_a).transpose()
                        statistic_b = getattr(Stats, stats_dict['function'])(data_array_b, **function_arguments_b).transpose()

                        # aggregate data per station (removing period dimension)
                        if base_zstat != 'NStations':
                            statistic_a = aggregation(statistic_a, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()
                            statistic_b = aggregation(statistic_b, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()

                # calculate statistics per station 
                else:
                    statistic_a = getattr(Stats, stats_dict['function'])(data_array_a, **function_arguments_a)
                    statistic_b = getattr(Stats, stats_dict['function'])(data_array_b, **function_arguments_b)

                # take difference: statistic_b - statistic_a
                z_statistic = statistic_b - statistic_a

            # else, is the difference statistic an experiment bias statistic (i.e. r)?
            elif z_statistic_type == 'expbias':

                # temporal colocation must be turned on for calculation, and have some experiments, if not return NaNs
                if (not read_instance.temporal_colocation) or (len(read_instance.data_labels) == 1):
                    if (map) or (per_station):
                        z_statistic = np.array([], dtype=np.float32)
                        if map:
                            station_inds = np.array([], dtype=np.int32)
                            return z_statistic, station_inds
                        elif per_station:
                            return z_statistic
                    else: 
                        if (period is not None) or (chunk_resolution is not None): 
                            stats_calc[zstat] = np.full((len(data_array_b),len(data_labels_b)), np.nan)
                        else:
                            stats_calc[zstat] = np.full((len(data_labels_b)), np.nan)
                        continue

                # load default selected z statistic arguments for passing to statistical function
                function_arguments = stats_dict['arguments']

                # calculate statistics
            
                # calculate periodic statistic per station
                if z_statistic_period is not None:
                    # if periodic statistic mode is cycle, then aggregate per periodic grouping, and then calculate stat
                    if periodic_statistic_mode == 'Cycle':
                        # aggregation in each group, per station, by periodic statistic
                        statistic_a = aggregation(data_array_a, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()
                        statistic_b = aggregation(data_array_b, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()
                        
                        # calculate statistic per station (removing period dimension)
                        z_statistic = getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':statistic_a,'exp':statistic_b}}).transpose()

                    # if periodic statistic mode is independent, then calculate stats independently per periodic grouping,
                    # and then aggregate 
                    elif periodic_statistic_mode == 'Independent':
                        # calculate statistic per periodic grouping per station
                        z_statistic = getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':data_array_a,'exp':data_array_b}}).transpose()

                        # aggregate data per station (removing period dimension)
                        z_statistic = aggregation(z_statistic, periodic_statistic_aggregation, read_instance=read_instance, axis=-1).transpose()

                # calculate statistics per station 
                else:
                    z_statistic = getattr(ExpBias, stats_dict['function'])(**{**function_arguments, **{'obs':data_array_a,'exp':data_array_b}})

        # if any calculated statistics are infinite, then set them to be NaNs 
        finite_boolean = np.isfinite(z_statistic)
        z_statistic[~finite_boolean] = np.nan

        # reshape forecast data
        if forecast_type == 'daily':
            if base_zstat in ['NStations','MDA8']:
                n_chunks, n_labels, n_forecast_days = z_statistic.shape
                z_statistic = z_statistic.transpose(0, 2, 1).reshape(-1, n_labels, order='F')
            else:
                n_chunks, n_labels, n_stations, n_forecast_days = z_statistic.shape
                z_statistic = z_statistic.transpose(0, 3, 1, 2).reshape(-1, n_labels, n_stations, order='F')
            #print('Calculated Stat Forecast reshape: ', z_statistic.shape)

        # return map statistics
        if map:
            # if any station z statistics come out as NaN/inf, cut z_statistic to remove invalid NaNs/infs, 
            # and also remove respective stations from station idncies
            finite_boolean = finite_boolean[0,:]
            z_statistic = z_statistic[0,:]
            return z_statistic[finite_boolean], station_inds[finite_boolean] 

        # return per station statistics
        elif per_station:
            return z_statistic

        # otherwise, save desired statistic for specific statistical calculation mode 
        else:
            if reduction:
                if (statistic_mode == 'Temporal|Spatial') & ((base_zstat not in ['NStations','MDA8']) or (z_statistic_period is not None)):
                    z_statistic = aggregation(z_statistic, statistic_aggregation, read_instance=read_instance, axis=-1)
                elif (statistic_mode in ['Flattened', 'Spatial|Temporal']) & ((base_zstat not in ['NStations','MDA8']) or (z_statistic_period is not None)):
                    if base_zstat in ['NStations','MDA8']:
                        z_statistic = aggregation(z_statistic, statistic_aggregation, read_instance=read_instance, axis=-1)
                    else:
                        z_statistic = np.squeeze(z_statistic, axis=-1)
            stats_calc[zstat] = z_statistic

            #print('Final: ', z_statistic.shape)

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
        return np.nan, np.nan

def generate_colourbar_detail(read_instance, zstat, plotted_min, plotted_max, plot_characteristics, speci, 
                              only_label=False):
    """ Function that generates neccessary detail to crate colourbar.

        :param read_instance: Instance of class Dashboard or Report
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
    # first check if have defined label (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    set_label = False
    #1. get label specific for z statistic
    if 'label' in stats_dict:
        if (stats_dict['label'] != '') and (stats_dict['label'] != {}):
            set_label = True
            # 2. get label specific for species
            if isinstance(stats_dict['label'], dict):
                if speci in stats_dict['label'].keys():
                    z_label = stats_dict['label'][speci]
                else:
                    set_label = False
            else:
                z_label = stats_dict['label']
            # adjust label to include units (if set)
            if set_label:
                if z_statistic_sign == 'absolute':
                    if label_units != '':
                        z_label = '{} [{}]'.format(z_label, label_units)
                    else:
                        z_label = copy.deepcopy(z_label)
                else:
                    if z_statistic_type == 'basic':
                        if label_units != '':
                            z_label = '{} bias [{}]'.format(z_label, label_units)
                        else:
                            z_label = '{} bias'.format(z_label)
                    else:
                        if label_units != '':
                            z_label = '{} [{}]'.format(z_label, label_units)
                        else:
                            z_label = copy.deepcopy(z_label)  
    # 3. check configuration file
    if not set_label:
        if 'cb_label' in plot_characteristics:
            if plot_characteristics['cb_label']['label'] != '':
                z_label = plot_characteristics['cb_label']['label']
                set_label = True
    # return label if only that is wanted
    if only_label:
        return z_label

    # set cmap for z statistic
    # first check if have defined cmap (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    set_cmap = False
    if z_statistic_sign == 'absolute':
        cmap_var_name = 'cmap_absolute'
    else:
        cmap_var_name = 'cmap_bias'
    #1. get cmap specific for z statistic
    if cmap_var_name in stats_dict:
        if (stats_dict[cmap_var_name] != '') and (stats_dict[cmap_var_name] != {}):
            set_cmap = True
            if isinstance(stats_dict[cmap_var_name], dict):
                if speci in stats_dict[cmap_var_name].keys():
                    z_colourmap = stats_dict[cmap_var_name][speci]
                    
                else:
                    error = f"Error: colourmap ({cmap_var_name}) is not defined for {speci}. "
                    error += f"{cmap_var_name} can be set as a string per statistic (for all species), or as a dict (per species)."
                    read_instance.logger.error(error)
                    sys.exit(1) 
            else:
                z_colourmap = stats_dict[cmap_var_name]
    #3. check configuration file
    if not set_cmap:
        if cmap_var_name in plot_characteristics['cb']:
            if (plot_characteristics['cb'][cmap_var_name] != '') and(plot_characteristics['cb'][cmap_var_name]):
                z_colourmap = plot_characteristics['cb'][cmap_var_name]
                set_cmap = True
    # if have no defined cmap, raise error
    if not set_cmap:
        error = f"Error: colourmap ({cmap_var_name}) for the colourbar needs to be defined, either in the "
        error += "configuration files for the map, or per statistic in 'basic_stats.yaml' or 'experiment_bias_stats.yaml'."
        read_instance.logger.error(error)
        sys.exit(1) 

    # check if have defined vmin (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    # if have no defined vmin, then take vmin as minimum range value of calculated statistic
    set_vmin = False
    if z_statistic_sign == 'absolute':
        vmin_var_name = 'vmin_absolute'
    else:
        vmin_var_name = 'vmin_bias'
    #1. get vmin specific for z statistic
    if vmin_var_name in stats_dict:
        if (stats_dict[vmin_var_name] != '') and (stats_dict[vmin_var_name] != {}):
            # 2. get vmin specific for species
            set_vmin = True
            if isinstance(stats_dict[vmin_var_name], dict):
                if speci in stats_dict[vmin_var_name].keys():
                    z_vmin = stats_dict[vmin_var_name][speci]
                else:
                    set_vmin = False
            else:
                z_vmin = stats_dict[vmin_var_name]
    #3. check configuration file
    if not set_vmin:
        if vmin_var_name in plot_characteristics['cb']:
            if (plot_characteristics['cb'][vmin_var_name] != '') and (plot_characteristics['cb'][vmin_var_name] != {}):
                z_vmin = plot_characteristics['cb'][vmin_var_name]
                set_vmin = True
    # if have no defined vmin, take vmin as minimum range value of calculated statistic
    if not set_vmin:
        z_vmin = plotted_min

    # check if have defined vmax (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    # if have no defined vmax, then take vmax as maximum range value of calculated statistic
    set_vmax = False
    if z_statistic_sign == 'absolute':
        vmax_var_name = 'vmax_absolute'
    else:
        vmax_var_name = 'vmax_bias'
    #1. get vmax specific for z statistic
    if vmax_var_name in stats_dict:   
        if (stats_dict[vmax_var_name] != '') and (stats_dict[vmax_var_name] != {}):
            set_vmax = True
            # 2. get vmax specific for species
            if isinstance(stats_dict[vmax_var_name], dict):
                if speci in stats_dict[vmax_var_name].keys():
                    z_vmax = stats_dict[vmax_var_name][speci]
                else:
                    set_vmax = False
            else:
                z_vmax = stats_dict[vmax_var_name]
    #3. check configuration file
    if not set_vmax:
        if vmax_var_name in plot_characteristics['cb']:
            if (plot_characteristics['cb'][vmax_var_name] != '') and (plot_characteristics['cb'][vmax_var_name] != {}):
                z_vmax = plot_characteristics['cb'][vmax_var_name]
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

    # check if have defined n_discrete (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    # if have no defined n_discrete, then take None
    set_n_discrete = False
    #1. get n_discrete specific for z statistic
    if 'n_discrete' in stats_dict:
        if (stats_dict['n_discrete'] != '') and (stats_dict['n_discrete'] != {}):
            set_n_discrete = True
            # 2. get n_discrete specific for species
            if isinstance(stats_dict['n_discrete'], dict):
                if speci in stats_dict['n_discrete'].keys():
                    n_discrete = stats_dict['n_discrete'][speci]
                else:
                    set_n_discrete = False
            else:
                n_discrete = stats_dict['n_discrete']    
    #3. check configuration file
    if not set_n_discrete:
        if 'n_discrete' in plot_characteristics['cb']:
            if (plot_characteristics['cb']['n_discrete'] != '') and (plot_characteristics['cb']['n_discrete'] != {}):
                n_discrete = plot_characteristics['cb']['n_discrete']
                set_n_discrete = True
    # if have no defined n_discrete, take None
    if not set_n_discrete:
        n_discrete = None
    
    # check if have defined n_ticks (in this order: 1. specific for z statistic 2. specific for species 3. configuration file)
    # if have no defined n_ticks, then raise error
    set_n_ticks = False
    #1. get n_ticks specific for z statistic
    if 'n_ticks' in stats_dict:
        if (stats_dict['n_ticks'] != '') and (stats_dict['n_ticks'] != {}):
            set_n_ticks = True
            # 2. get n_ticks specific for species
            if isinstance(stats_dict['n_ticks'], dict):
                if speci in stats_dict['n_ticks'].keys():
                    n_ticks = stats_dict['n_ticks'][speci]
                else:
                    set_n_ticks = False
            else:
                n_ticks = stats_dict['n_ticks']
    #3. check configuration file
    if not set_n_ticks:
        if 'n_ticks' in plot_characteristics['cb']:
            if (plot_characteristics['cb']['n_ticks'] != '') and (plot_characteristics['cb']['n_ticks'] != {}):
                n_ticks = plot_characteristics['cb']['n_ticks']
                set_n_ticks = True
    # if have no defined n_ticks, raise error
    if not set_n_ticks:
        error = 'Error: The number of ticks (n_ticks) in the colorbar need to be defined, either in the '
        error += 'configuration files for the map or per statistic.'
        read_instance.logger.error(error)
        sys.exit(1) 

    return z_vmin, z_vmax, z_label, z_colourmap, n_discrete, n_ticks

def generate_colourbar(read_instance, axs, cb_axs, zstat, plot_characteristics, speci):
    """ Function that generates colourbar.

        :param read_instance: Instance of class Dashboard or Report
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
    z_vmin, z_vmax, z_label, z_colourmap, n_discrete, n_ticks = generate_colourbar_detail(read_instance, zstat, plotted_min, plotted_max, 
                                                                                          plot_characteristics, speci)

    # generate colourbar tick array
    tick_array = np.linspace(z_vmin, z_vmax, n_ticks, endpoint=True)

    # normalise colourbar range (into the 0.0 - 1.0 interval)
    norm = matplotlib.colors.Normalize(vmin=z_vmin, vmax=z_vmax)

    # get color palettes
    color_palettes = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/color_palettes.yaml')))

    # get cmap (handling discrete cases)
    if z_vmin != z_vmax:
        if z_colourmap in color_palettes.keys():
            cmap = colors.LinearSegmentedColormap.from_list("providentia", color_palettes[z_colourmap], n_discrete)
        else:
            cmap = plt.get_cmap(z_colourmap, n_discrete)
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
            # https://github.com/BSC-ES/providentia/issues/159
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
    if plot_type is not None:
        # have zstat in plot_type name?
        if ('-' in plot_type) & ('-violin' not in plot_type) & ('-target' not in plot_type) & ('-statsummary' not in plot_type):
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
    if zstat is not None:
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
    
def aggregation(data_array, statistic_aggregation, read_instance=None, axis=0):
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
                                           axis=axis, method='nearest')
    else:
        if read_instance:
            error = 'Aggregation statistic {0} is not available. '.format(statistic_aggregation)
            error += 'The options are: Mean, Median, p1, p5, p10, p25, p75, p90, p95 and p99'
            read_instance.logger.error(error)
            sys.exit(1) 

    return aggregated_data

def exceedance_lim(read_instance, networkspeci):
    """ Return the exceedance limit depending on the species input. 
        If species doesn't have a reported limit, returns np.nan.

        Try to get limit for specific networkspeci first, and then species.

        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3)
        :type networkspeci: str
        :return: value of exceedance limit
        :rtype: int
    """

    # get speci
    speci = networkspeci.split('|')[1]

    exceedance_limits = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/exceedances.yaml')))
    if networkspeci in exceedance_limits:
        limit = exceedance_limits[networkspeci]['value']
        initial_units = exceedance_limits[networkspeci]['units']
    elif speci in exceedance_limits:
        limit = exceedance_limits[speci]['value']
        initial_units = exceedance_limits[speci]['units']
    else:
        limit = np.nan
    
    if limit is not np.nan:
        # get output units
        standard_parameter_speci = get_standard_parameters_by_speci(speci, read_instance.ghost_version)
        final_units = read_instance.measurement_units[speci]

        # convert units using conversion factor
        conversion_factor = get_conversion_factor(initial_units, final_units, standard_parameter_speci)
        if isinstance(conversion_factor, str):
            read_instance.logger.error(conversion_factor)
            sys.exit(1)
        limit *= conversion_factor

    return limit


def get_fairmode_data(read_instance, canvas_instance, networkspeci, data_labels):
    
    # get coverage
    speci = networkspeci.split('|')[1]
    coverage = fairmode_settings[speci]['coverage']

    # get data per station
    data_array = copy.deepcopy(read_instance.data_in_memory_filtered[networkspeci])

    # temporally colocate data (if active)
    if read_instance.temporal_colocation:
        data_array[:, read_instance.temporal_colocation_nans[networkspeci]] = np.nan

    # get data cut for relevant stations
    data_array = data_array[:,canvas_instance.station_inds[networkspeci],:]

    # if hourly data then make sure days with less than 75% coverage are nan
    if read_instance.active_resolution == 'hourly':

        # get number of days
        num_days = data_array.shape[-1] // 24

        # reshape data array to have days as time dimension
        daily_data = data_array.reshape(data_array.shape[0], data_array.shape[1], num_days, 24)  

        # get number of non-nan hours per day
        non_nan_count = np.count_nonzero(~np.isnan(daily_data), axis=-1)

        # get days that need to be nan because there are at least 25% of the hours per day that are nan
        threshold = (coverage / 100) * 24 # 18 hours out of 24 for 75% coverage
        days_to_nan = non_nan_count < threshold

        # convert days with less than 75% coverage to nan
        days_to_nan_expanded = np.repeat(days_to_nan, 24, axis=-1)
        data_array[days_to_nan_expanded] = np.nan

    # get indices of valid stations
    obs_representativity = Stats.calculate_data_avail_fraction(data_array[0, :, :])
    valid_station_idxs = obs_representativity >= coverage

    # do some extra processing for hourly resolution data
    if read_instance.active_resolution == 'hourly':

        # resample PM10/PM2.5 data to daily
        if speci in ['pm2p5', 'pm10']:

            # flatten networkspecies dimension for creation of pandas dataframe
            data_array_reduced = data_array.reshape(data_array.shape[0]*data_array.shape[1], data_array.shape[2])
            
            # create pandas dataframe of data array
            data_array_df = pd.DataFrame(data_array_reduced.transpose(), index=read_instance.time_index, 
                                         columns=np.arange(data_array_reduced.shape[0]), dtype=np.float32)
            # resample data array
            data_array_df_resampled = data_array_df.resample('D', axis=0).mean()

            # save back out as numpy array (reshaping to get back networkspecies dimension)
            data_array_resampled = data_array_df_resampled.to_numpy().transpose()
            data_array = data_array_resampled.reshape(data_array.shape[0], data_array.shape[1], 
                                                      data_array_resampled.shape[1])
        # calculate MDA8 for ozone
        elif speci in ['sconco3']:

            # get valid data labels for networkspeci
            valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

            # calculate MDA8 (applying mask for days with less than 75% coverage)
            stats_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, 'MDA8', 
                                             valid_data_labels, [], per_station=True, chunk_resolution='daily',
                                             mask=days_to_nan_expanded)

            # reshape data array
            data_array = np.transpose(stats_calc, (1, 2, 0))

    # filter by valid stations
    data_array = data_array[:, valid_station_idxs, :]

    return data_array, valid_station_idxs