""" Auxiliary reading functions """

import datetime
from glob import glob
import json
import os
import sys

import bisect
import cftime
from netCDF4 import Dataset, num2date, chartostring
import numpy as np
import pandas as pd

from providentia.auxiliar import CURRENT_PATH, join
from providentia.warnings_prv import show_message

_CFTIME_CLASS = {
    "standard": cftime.DatetimeGregorian,
    "gregorian": cftime.DatetimeGregorian,
    "proleptic_gregorian": cftime.DatetimeProlepticGregorian,
    "julian": cftime.DatetimeJulian,
    "noleap": cftime.DatetimeNoLeap,
    "365_day": cftime.DatetimeNoLeap,
    "all_leap": cftime.DatetimeAllLeap,
    "366_day": cftime.DatetimeAllLeap,
    "360_day": cftime.Datetime360Day,
}

# initialise dictionary for storing pointers to shared memory variables in read step
shared_memory_vars = {}

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])


def drop_nans(data):
    """
    Returns a 1D numpy array with all NaN values removed.

    Parameters
    ----------
    data : numpy.ndarray or pandas.Series
        The input data containing potential NaN values.

    Returns
    -------
    numpy.ndarray
        A filtered 1D array containing only valid numerical values.
    """
    return data[~pd.isnull(data)]


def init_shared_vars_read_netcdf_data(
    data_in_memory,
    data_in_memory_shape,
    ghost_data_in_memory,
    ghost_data_in_memory_shape,
    timestamp_array,
    qa,
    flags,
):
    """
    Initialises worker processes by granting access to shared memory variables before parallel NetCDF reading.

    Parameters
    ----------
    data_in_memory : multiprocessing.RawArray
        Shared memory buffer for storing the primary scientific data.
    data_in_memory_shape : tuple
        The dimensions of the primary data array.
    ghost_data_in_memory : multiprocessing.RawArray
        Shared memory buffer for storing metadata or GHOST-specific processed data.
    ghost_data_in_memory_shape : tuple
        The dimensions of the GHOST data array.
    timestamp_array : multiprocessing.RawArray
        Shared memory buffer containing the temporal coordinates for the data.
    qa : multiprocessing.RawArray
        Shared memory buffer for Quality Assurance scores.
    flags : multiprocessing.RawArray
        Shared memory buffer for data quality flags.
    """

    shared_memory_vars["data_in_memory"] = data_in_memory
    shared_memory_vars["data_in_memory_shape"] = data_in_memory_shape
    shared_memory_vars["ghost_data_in_memory"] = ghost_data_in_memory
    shared_memory_vars["ghost_data_in_memory_shape"] = ghost_data_in_memory_shape
    shared_memory_vars["timestamp_array"] = timestamp_array
    shared_memory_vars["qa"] = qa
    shared_memory_vars["flag"] = flags


def read_netcdf_data(tuple_arguments):
    """
    Handles the reading and filtering of observational or model netCDF data.

    Parameters
    ----------
    tuple_arguments : tuple
        A collection of arguments including file paths, station references, species names,
        shared memory handles, and quality control settings.

    Returns
    -------
    numpy.ndarray or list or None
        Returns a structured metadata array if reading observations; an empty list if no
        stations are found; or None for model data or missing files.
    """

    # assign arguments from tuple to variables
    (
        relevant_file,
        station_references,
        station_names,
        speci,
        observations_data_label,
        data_label,
        data_labels,
        reading_ghost,
        ghost_data_vars_to_read,
        metadata_dtype,
        metadata_vars_to_read,
        logger,
        default_qa,
        filter_read,
        network,
        forecast_indices,
    ) = tuple_arguments

    # wrap shared arrays as numpy arrays to more easily manipulate the data
    data_in_memory = np.frombuffer(
        shared_memory_vars["data_in_memory"], dtype=np.float32
    ).reshape(shared_memory_vars["data_in_memory_shape"][:])
    if (reading_ghost or network == "actris/actris") & (
        data_label == observations_data_label
    ):
        qa = np.frombuffer(shared_memory_vars["qa"], dtype=np.uint8)
        flags = np.frombuffer(shared_memory_vars["flag"], dtype=np.uint8)
        if reading_ghost:
            if not filter_read:
                ghost_data_in_memory = np.frombuffer(
                    shared_memory_vars["ghost_data_in_memory"], dtype=np.float32
                ).reshape(shared_memory_vars["ghost_data_in_memory_shape"][:])
    timestamp_array = np.frombuffer(
        shared_memory_vars["timestamp_array"], dtype=np.int64
    )

    # read netCDF frame
    ncdf_root = Dataset(relevant_file, mode="r")

    # get file time as integer timestamp
    time_var = ncdf_root.variables["time"]
    file_timestamp = time_var_to_asi8(time_var)

    # get valid file time indices (i.e. those times in active full time array)
    valid_file_time_indices = np.where(
        np.logical_and(
            file_timestamp >= timestamp_array[0], file_timestamp <= timestamp_array[-1]
        )
    )[0]

    # get indices relative to active full timestamp array
    full_array_time_indices = np.searchsorted(
        timestamp_array, file_timestamp[valid_file_time_indices]
    )

    # get all station references in file (do little extra work to get non-GHOST observational station references)
    if (not reading_ghost) & (data_label == observations_data_label):
        if "station_reference" in ncdf_root.variables:
            station_reference_var = "station_reference"
        elif "station_code" in ncdf_root.variables:
            station_reference_var = "station_code"
        elif "station_name" in ncdf_root.variables:
            station_reference_var = "station_name"
        else:
            error = "Error: {} cannot be read because it has no station_name.".format(
                relevant_file
            )
            logger.error(error)
            sys.exit(1)

        meta_shape = ncdf_root[station_reference_var].shape
        file_station_references = ncdf_root[station_reference_var][:]
        meta_val_dtype = np.array([file_station_references[0]]).dtype

        if len(meta_shape) == 2:
            if meta_val_dtype == np.dtype(object):
                file_station_references = np.array(
                    ["".join(val) for val in file_station_references]
                )
            else:
                file_station_references = chartostring(file_station_references)

    # GHOST and interpolated model data
    else:
        file_station_references = ncdf_root["station_reference"][:]

    # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
    refs = np.asarray(file_station_references).astype(str)
    original_indices = np.arange(refs.size, dtype=np.int32)
    mask_not_nan = np.char.lower(refs) != "nan"
    refs_non_nan = refs[mask_not_nan]
    non_nan_original_indices = original_indices[mask_not_nan]

    # get indices of file station station references that are contained in all unique station references array
    station_ref_set = set(station_references)
    mask_present = np.fromiter(
        (r in station_ref_set for r in refs_non_nan),
        dtype=bool,
        count=len(refs_non_nan),
    )
    matched_filtered_indices = np.where(mask_present)[0]
    current_file_station_indices = non_nan_original_indices[matched_filtered_indices]

    # for all unique station references that are contained within file station references array
    # get the index of the station reference in the unique station references array
    station_to_index = {ref: i for i, ref in enumerate(station_references)}
    valid_refs = refs_non_nan[matched_filtered_indices]
    full_array_station_indices = np.fromiter(
        (station_to_index[r] for r in valid_refs), dtype=np.int32, count=len(valid_refs)
    )

    # if have zero current_file_station_indices in all unique station references,
    # then check if it is because of old-style of Providentia-interpolation output,
    # where all station_references were for 'station_name'
    if (data_label != observations_data_label) & (
        len(current_file_station_indices) == 0
    ):
        # get indices of file station station references that are contained in all unique station references array
        current_file_station_indices = np.where(
            np.in1d(file_station_references, station_names)
        )[0]

        # for all unique station references that are contained within file station references array
        # get the index of the station reference in the unique station references array
        index = np.argsort(station_names)
        sorted_station_names = station_names[index]
        sorted_index = np.searchsorted(
            sorted_station_names, file_station_references[current_file_station_indices]
        )
        full_array_station_indices = np.take(index, sorted_index, mode="clip")

    # if still have zero current_file_station_indices in all unique station references (can happen due to station colocation)
    # then return from function without reading
    if len(current_file_station_indices) == 0:
        # return empty metadata list if reading observations
        if (data_label == observations_data_label) & (not filter_read):
            return []
        else:
            return

    # get first and last station index for current file to read only relevant subset of data,
    # and adjust current_file_station_indices to be relative to subset of data to read
    fsi = current_file_station_indices[0]
    lsi = current_file_station_indices[-1] + 1
    current_file_station_indices_adjusted = current_file_station_indices - fsi

    # read observations
    if data_label == observations_data_label:
        # read species variable
        # GHOST
        if reading_ghost:
            # if need to filter by qa load non-filtered array, otherwise load prefiltered array (if available)
            if (default_qa) & (
                "{}_prefiltered_defaultqa".format(speci)
                in list(ncdf_root.variables.keys())
            ):
                species_data = ncdf_root["{}_prefiltered_defaultqa".format(speci)][
                    fsi:lsi, valid_file_time_indices
                ]
                species_data = species_data[current_file_station_indices_adjusted]
                # set read_qa to False to not filter by them
                read_qa = False
            else:
                species_data = ncdf_root[speci][fsi:lsi, valid_file_time_indices]
                species_data = species_data[current_file_station_indices_adjusted]
                read_qa = True
        # non-GHOST
        else:
            # transpose array to swap station and time dimensions
            if ncdf_root[speci].dimensions == ("time", "station"):
                species_data = ncdf_root[speci][valid_file_time_indices, fsi:lsi].T
                species_data = species_data[current_file_station_indices_adjusted]
            # do not transpose
            else:
                species_data = ncdf_root[speci][fsi:lsi, valid_file_time_indices]
                species_data = species_data[current_file_station_indices_adjusted]

            # if network is actris then read qa
            if network == "actris/actris":
                read_qa = True

        # reading GHOST data?
        if reading_ghost:
            # read GHOST data variables
            if not filter_read:
                for ghost_data_var_ii, ghost_data_var in enumerate(
                    ghost_data_vars_to_read
                ):
                    ghost_data = ncdf_root[ghost_data_var][
                        fsi:lsi, valid_file_time_indices
                    ]
                    ghost_data = ghost_data[current_file_station_indices_adjusted]
                    ghost_data_in_memory[
                        ghost_data_var_ii,
                        full_array_station_indices[:, np.newaxis],
                        full_array_time_indices[np.newaxis, :],
                    ] = ghost_data

        if (reading_ghost) or (network == "actris/actris"):
            # if some qa flags selected then screen observations
            if read_qa:
                if len(qa) > 0:
                    # screen out observations which are associated with any of the selected qa flags
                    qa_data = ncdf_root["qa"][fsi:lsi, valid_file_time_indices, :]
                    qa_data = qa_data[current_file_station_indices_adjusted]
                    species_data[np.isin(qa_data, qa).any(axis=2)] = np.nan

            # if some data provider flags selected then screen observations
            if len(flags) > 0:
                # screen out observations which are associated with any of the selected data provider flags
                flag_data = ncdf_root["flag"][fsi:lsi, valid_file_time_indices, :]
                flag_data = flag_data[current_file_station_indices_adjusted]
                species_data[np.isin(flag_data, flags).any(axis=2)] = np.nan

        # write filtered species data to shared file data
        data_in_memory[
            data_labels.index(observations_data_label),
            full_array_station_indices[:, np.newaxis],
            full_array_time_indices[np.newaxis, :],
        ] = species_data

        # get file metadata
        if not filter_read:
            file_metadata = np.full(
                (len(station_references), 1), np.nan, dtype=metadata_dtype
            )
            for meta_var in metadata_vars_to_read:
                # do extra work for non-GHOST data
                if not reading_ghost:
                    # get correct variable name for .nc
                    if meta_var == "longitude":
                        if "longitude" in ncdf_root.variables:
                            meta_var_nc = "longitude"
                        else:
                            meta_var_nc = "lon"
                    elif meta_var == "latitude":
                        if "latitude" in ncdf_root.variables:
                            meta_var_nc = "latitude"
                        else:
                            meta_var_nc = "lat"
                    elif meta_var == "altitude":
                        if "altitude" in ncdf_root.variables:
                            meta_var_nc = "altitude"
                        else:
                            meta_var_nc = "alt"
                    elif meta_var == "station_reference":
                        if "station_reference" in ncdf_root.variables:
                            meta_var_nc = "station_reference"
                        elif "station_code" in ncdf_root.variables:
                            meta_var_nc = "station_code"
                        elif "station_name" in ncdf_root.variables:
                            meta_var_nc = "station_name"
                    elif meta_var == "area_classification":
                        if "area_classification" in ncdf_root.variables:
                            meta_var_nc = "area_classification"
                        else:
                            meta_var_nc = "station_area"
                    elif meta_var == "station_classification":
                        if "station_classification" in ncdf_root.variables:
                            meta_var_nc = "station_classification"
                        else:
                            meta_var_nc = "station_type"
                    else:
                        meta_var_nc = meta_var

                    # check meta variable is in netCDF
                    if meta_var_nc not in ncdf_root.variables:
                        continue

                    # get metadata
                    meta_shape = ncdf_root[meta_var_nc].shape
                    meta_val = ncdf_root[meta_var_nc][fsi:lsi]
                    meta_val = meta_val[current_file_station_indices_adjusted]
                    meta_val_dtype = np.array([meta_val[0]]).dtype

                    # do str formatting where neccessary
                    if meta_val_dtype not in [
                        np.int8,
                        np.int16,
                        np.int32,
                        np.int64,
                        np.uint8,
                        np.uint16,
                        np.uint32,
                        np.uint64,
                        np.float16,
                        np.float32,
                        np.float64,
                    ]:
                        if len(meta_shape) == 2:
                            if meta_val_dtype == np.dtype(object):
                                meta_val = np.array(["".join(val) for val in meta_val])
                            else:
                                meta_val = chartostring(meta_val)

                    # do str formatting (capitalisation) to the metadata
                    if isinstance(meta_val, str):
                        meta_val = np.char.capitalize(meta_val)

                # GHOST metadata
                else:
                    meta_var_nc = meta_var

                    # check meta variable is in netCDF
                    if meta_var_nc not in ncdf_root.variables:
                        continue

                    # get metadata
                    meta_val = ncdf_root[meta_var_nc][fsi:lsi]
                    meta_val = meta_val[current_file_station_indices_adjusted]

                # put metadata in array
                file_metadata[meta_var][full_array_station_indices, 0] = meta_val

    # model data
    else:
        # determine if data is structured as forecast data or not
        if "forecast_day" in ncdf_root[speci].dimensions:
            have_forecast = True
            # if so, check what type of forecast data
            # check if have daily forecast set
            daily_forecast = np.any(
                [True for data_label in data_labels if "-daily" in data_label]
            )
            # check if have combined forecast set
            combined_forecast = np.any(
                [True for data_label in data_labels if "-combined" in data_label]
            )
            # check if have day forecast set
            day_forecast = np.any(
                [True for data_label in data_labels if "-day" in data_label]
            )
        # data is not structured as forecast
        else:
            have_forecast = False

        # if have no passed forecast indices, then take first day preferentially (if data is structured as forecast data)
        if len(forecast_indices) == 0:
            forecast_indices = np.array([0], dtype=np.int32)

        # iterate through forecast indices
        for forecast_index in forecast_indices:
            # if want a specific forecast day, and it is available in the netCDF, then take it
            if have_forecast:
                # daily forecast
                if daily_forecast:
                    data_label_forecast = "{}-daily{}".format(
                        data_label, forecast_index + 1
                    )
                # combined forecast
                elif combined_forecast:
                    data_label_forecast = "{}-combined{}".format(
                        data_label, forecast_index + 1
                    )
                # N day forecast
                elif day_forecast:
                    data_label_forecast = "{}-day{}".format(
                        data_label, forecast_index + 1
                    )
                else:
                    data_label_forecast = data_label
                relevant_data = ncdf_root[speci][
                    fsi:lsi, valid_file_time_indices, forecast_index
                ]
                relevant_data = relevant_data[current_file_station_indices_adjusted]

            # else if forecast day not available in the netCDF, then just take the data as it is
            else:
                relevant_data = ncdf_root[speci][fsi:lsi, valid_file_time_indices]
                relevant_data = relevant_data[current_file_station_indices_adjusted]
                data_label_forecast = data_label

            # mask out fill values for parameter field
            relevant_data[relevant_data.mask] = np.nan

            # put data in array
            data_in_memory[
                data_labels.index(data_label_forecast),
                full_array_station_indices[:, np.newaxis],
                full_array_time_indices[np.newaxis, :],
            ] = relevant_data

    # close netCDF
    ncdf_root.close()

    # return metadata if reading observations
    if (data_label == observations_data_label) & (not filter_read):
        return file_metadata


def read_netcdf_metadata(tuple_arguments):
    """
    Extracts basic metadata from an observational NetCDF file.

    Parameters
    ----------
    tuple_arguments : tuple
        A collection containing the file path to be read, a boolean indicating if it is a GHOST
        format file, and a logger instance for error reporting.

    Returns
    -------
    list
        A list of NumPy arrays or lists containing the extracted metadata for
        'station_reference', 'longitude', 'latitude', 'station_name', and 'measurement_altitude'.
    """

    # assign arguments from tuple to variables
    relevant_file, reading_ghost, logger = tuple_arguments

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # set metadata variables to read
    metadata_vars_to_read = [
        "station_reference",
        "longitude",
        "latitude",
        "station_name",
        "measurement_altitude",
    ]
    metadata_read = []
    print(relevant_file)
    # iterate though metadata variables to read
    for meta_var in metadata_vars_to_read:
        # do extra work for non-GHOST data
        if not reading_ghost:
            # station reference
            if meta_var == "station_reference":
                if "station_reference" in ncdf_root.variables:
                    station_reference_var = "station_reference"
                elif "station_code" in ncdf_root.variables:
                    station_reference_var = "station_code"
                elif "station_name" in ncdf_root.variables:
                    station_reference_var = "station_name"
                else:
                    error = "Error: {} cannot be read because it has no station_name.".format(
                        relevant_file
                    )
                    logger.error(error)
                    sys.exit(1)

                meta_shape = ncdf_root[station_reference_var].shape
                meta_val = ncdf_root[station_reference_var][:]
                meta_val_dtype = np.array([meta_val[0]]).dtype
                if len(meta_shape) == 2:
                    if meta_val_dtype == np.dtype(object):
                        meta_val = np.array(["".join(val) for val in meta_val])
                    else:
                        meta_val = chartostring(meta_val)

                # get indices of all non-NaN stations (can be NaN for some non-GHOST files)
                non_nan_station_indices = np.array(
                    [
                        ref_ii
                        for ref_ii, ref in enumerate(meta_val)
                        if ref.lower() != "nan"
                    ]
                )

                # if have zero non-NaN station indices, then return from function without reading
                if len(non_nan_station_indices) == 0:
                    metadata_read = [np.array([]), np.array([]), np.array([]), [], []]
                    ncdf_root.close()
                    return metadata_read

                # get station references for non-NaN stations
                meta_val = meta_val[non_nan_station_indices]

                # get first and last non_nan_station_indices
                fi = non_nan_station_indices[0]
                li = non_nan_station_indices[-1] + 1
                non_nan_station_indices_adjusted = non_nan_station_indices - fi

            # longitude
            elif meta_var == "longitude":
                print(ncdf_root.variables.keys())
                if "longitude" in ncdf_root.variables:
                    meta_val = ncdf_root["longitude"][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]
                else:
                    meta_val = ncdf_root["lon"][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]

            # latitude
            elif meta_var == "latitude":
                if "latitude" in ncdf_root.variables:
                    meta_val = ncdf_root["latitude"][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]
                else:
                    meta_val = ncdf_root["lat"][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]

            # station name
            elif meta_var == "station_name":
                if "station_name" in ncdf_root.variables:
                    meta_shape = ncdf_root["station_name"].shape
                    meta_val = ncdf_root["station_name"][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]
                    meta_val_dtype = np.array([meta_val[0]]).dtype
                    if len(meta_shape) == 2:
                        if meta_val_dtype == np.dtype(object):
                            meta_val = np.array(["".join(val) for val in meta_val])
                        else:
                            meta_val = chartostring(meta_val)
                else:
                    meta_val = []

            # measurement altitude
            elif meta_var == "measurement_altitude":
                if meta_var in ncdf_root.variables:
                    meta_val = ncdf_root[meta_var][fi:li]
                    meta_val = meta_val[non_nan_station_indices_adjusted]
                else:
                    meta_val = []

        # GHOST metadata
        else:
            meta_val = ncdf_root[meta_var][:]

        # append read metadata
        metadata_read.append(meta_val)

    # close netCDF
    ncdf_root.close()

    # return read metadata
    return metadata_read


def check_forecast_dimension(tuple_arguments):
    """
    Determines if a NetCDF file contains a forecast dimension and returns the number of lead days.

    Parameters
    ----------
    tuple_arguments : str
        The file path to the NetCDF file being inspected.

    Returns
    -------
    int
        The size of the 'forecast_day' dimension, or 0 if the dimension is not present.
    """

    # assign arguments from tuple to variables
    relevant_file = tuple_arguments

    # read netCDF frame
    ncdf_root = Dataset(relevant_file)

    # check if forecast day dimension exists
    if "forecast_day" in ncdf_root.dimensions:
        n_forecast_days = ncdf_root.dimensions["forecast_day"].size
    else:
        n_forecast_days = 0

    # close netCDF
    ncdf_root.close()

    return n_forecast_days


def get_yearmonths_to_read(
    available_yearmonths, start_date_to_read, end_date_to_read, resolution
):
    """
    Determines which monthly data files fall within a specified date range based on temporal resolution.

    Parameters
    ----------
    available_yearmonths : list of str
        List of available file strings in 'YYYYMM' format.
    start_date_to_read : str or int
        The beginning of the requested period in 'YYYYMMDD' format.
    end_date_to_read : str or int
        The end of the requested period in 'YYYYMMDD' format.
    resolution : str
        The temporal resolution of the data (e.g., 'monthly', 'daily', or 'hourly').

    Returns
    -------
    list of str
        A subset of available_yearmonths representing the files required to cover the date range.
    """

    available_yearmonthdays = [
        int(yearmonth + "01") for yearmonth in available_yearmonths
    ]

    first_valid_file_ind = bisect.bisect_right(
        available_yearmonthdays, int(start_date_to_read)
    )
    if first_valid_file_ind != 0:
        first_valid_file_ind = first_valid_file_ind - 1
    last_valid_file_ind = bisect.bisect_left(
        available_yearmonthdays, int(end_date_to_read)
    )

    # read only complete months
    if resolution == "monthly":
        if str(end_date_to_read)[6:8] != "01":
            if str(end_date_to_read)[0:6] == str(available_yearmonths[-1]):
                last_valid_file_ind -= 1
        if str(start_date_to_read)[6:8] != "01":
            if str(start_date_to_read)[0:6] == str(available_yearmonths[0]):
                first_valid_file_ind += 1

    if first_valid_file_ind == last_valid_file_ind:
        return [available_yearmonths[first_valid_file_ind]]
    else:
        return available_yearmonths[first_valid_file_ind:last_valid_file_ind]


def get_default_qa(instance, speci):
    """
    Returns the standard GHOST quality assurance flags based on a species' physical limits.

    Parameters
    ----------
    instance : object
        An instance of the application class containing parameter metadata.
    speci : str
        The name of the chemical species or parameter to evaluate.

    Returns
    -------
    list
        A sorted list of integer codes representing the default QA flags.
    """

    if instance.parameter_dictionary[speci]["extreme_lower_limit"] < 0.0:
        return sorted(instance.default_qa_non_negative)
    else:
        return sorted(instance.default_qa_standard)


def get_frequency_code(resolution):
    """
    Maps human-readable temporal resolutions to Pandas offset aliases.

    Parameters
    ----------
    resolution : str
        The descriptive name of the temporal resolution (e.g., 'hourly', 'daily').

    Returns
    -------
    active_frequency_code : str
        The corresponding Pandas frequency string used for resampling and date ranges.
    """

    if resolution in ["hourly", "hourly_instantaneous"]:
        active_frequency_code = "h"
    elif resolution in ["3hourly", "3hourly_instantaneous"]:
        active_frequency_code = "3h"
    elif resolution in ["6hourly", "6hourly_instantaneous"]:
        active_frequency_code = "6h"
    elif resolution == "daily":
        active_frequency_code = "D"
    elif resolution == "monthly":
        active_frequency_code = "MS"
    elif resolution == "annual":
        active_frequency_code = "YS"

    return active_frequency_code


def get_chunk_size(active_resolution, chunk_resolution):
    """
    Calculates the number of data steps required to fill a specific temporal chunk.

    Parameters
    ----------
    active_resolution : str
        The base temporal resolution of the data (e.g., 'hourly', '3hourly').
    chunk_resolution : str
        The target resolution for data partitioning (e.g., 'daily').

    Returns
    -------
    int
        The ratio of the target duration to the base duration.
    """

    # convert resolutions to hours
    hours_dict = {
        "hourly": 1,
        "3hourly": 3,
        "6hourly": 6,
        "daily": 24,
    }

    base_hours = hours_dict[active_resolution]
    target_hours = hours_dict[chunk_resolution]

    return target_hours // base_hours


def check_for_ghost(network_name):
    """
    Determines if a given monitoring network comes from GHOST or not.

    Parameters
    ----------
    network_name : str
        The name of the atmospheric monitoring network to verify.

    Returns
    -------
    bool
        True if the network is a GHOST network, False otherwise.
    """

    if "/" in network_name:
        return False
    else:
        return True


def get_ghost_observational_tree(instance):
    """
    Creates a nested dictionary representing the GHOST observational data tree and exports it to JSON.

    Parameters
    ----------
    instance : object
        An instance of the application class containing GHOST configuration and versioning.

    Returns
    -------
    dict
        A nested dictionary structure: {network: {resolution: {matrix: {species: [YYYYMM, ...]}}}}.
    """

    # create dictionary for storing filetree
    ghost_observation_data = {}

    # iterate through available networks
    for network in instance.ghost_available_networks:
        # check if directory for network exists
        # if not, continue
        if not os.path.exists(
            "%s/%s/%s" % (instance.ghost_root, network, instance.ghost_version)
        ):
            continue

        # write empty dictionary for network
        ghost_observation_data[network] = {}

        # iterate through available resolutions
        for resolution in instance.ghost_available_resolutions:
            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists(
                "%s/%s/%s/%s"
                % (instance.ghost_root, network, instance.ghost_version, resolution)
            ):
                continue

            # write nested empty dictionary for resolution
            ghost_observation_data[network][resolution] = {}

            # get available species for network/resolution
            available_species = os.listdir(
                "%s/%s/%s/%s"
                % (instance.ghost_root, network, instance.ghost_version, resolution)
            )

            # iterate through available files per species
            for speci in available_species:
                # get all available netCDF files
                available_files = os.listdir(
                    "%s/%s/%s/%s/%s"
                    % (
                        instance.ghost_root,
                        network,
                        instance.ghost_version,
                        resolution,
                        speci,
                    )
                )

                # continue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted(
                    [f.split("_")[-1][:6] for f in available_files if f != "temporary"]
                )

                # get matrix for current species
                if speci in instance.parameter_dictionary:
                    matrix = instance.parameter_dictionary[speci]["matrix"]
                    if matrix not in ghost_observation_data[network][resolution]:
                        # write nested empty dictionary for matrix
                        ghost_observation_data[network][resolution][matrix] = {}

                    # write nested dictionary for species, with associated file yearmonths
                    ghost_observation_data[network][resolution][matrix][
                        speci
                    ] = file_yearmonths

    # save file tree out to yaml
    with open(
        join(
            PROVIDENTIA_ROOT,
            "settings/internal/ghost_filetree_{}.json".format(instance.ghost_version),
        ),
        "w",
    ) as json_file:
        json.dump(ghost_observation_data, json_file, indent=4)

    return ghost_observation_data


def get_nonghost_observational_tree(instance):
    """
    Creates a nested dictionary representing the non-GHOST observational data tree and exports it to JSON.

    Parameters
    ----------
    instance : object
        An instance of the application class containing non-GHOST configuration and directory roots.

    Returns
    -------
    dict
        A nested dictionary structure: {network: {resolution: {matrix: {species: [YYYYMM, ...]}}}}.
    """

    # create dictionary for storing filetree
    nonghost_observation_data = {}

    # iterate through networks
    for network in instance.nonghost_available_networks:
        # check if directory for network exists
        # if not, continue
        if not os.path.exists("%s/%s" % (instance.nonghost_root, network)):
            continue

        # write empty dictionary for network
        nonghost_observation_data[network] = {}

        # iterate through available resolutions
        for resolution in instance.nonghost_available_resolutions:
            # check if directory for resolution exists
            # if not, continue
            if not os.path.exists(
                "%s/%s/%s" % (instance.nonghost_root, network, resolution)
            ):
                continue

            # write nested empty dictionary for resolution
            nonghost_observation_data[network][resolution] = {}

            # get available species for network/resolution
            available_species = os.listdir(
                "%s/%s/%s" % (instance.nonghost_root, network, resolution)
            )

            # iterate through available files per species
            for speci in available_species:
                # get all available netCDF files
                available_files = glob(
                    "%s/%s/%s/%s/%s_??????.nc"
                    % (instance.nonghost_root, network, resolution, speci, speci)
                )

                # continue if have no files
                if len(available_files) == 0:
                    continue

                # get monthly start date (YYYYMM) of all files
                file_yearmonths = sorted(
                    [f.split("_")[-1][:6] for f in available_files]
                )

                # get matrix for current species
                if speci in instance.parameter_dictionary:
                    matrix = instance.parameter_dictionary[speci]["matrix"]
                    if matrix not in nonghost_observation_data[network][resolution]:
                        # write nested empty dictionary for matrix
                        nonghost_observation_data[network][resolution][matrix] = {}

                    # write nested dictionary for species, with associated file yearmonths
                    nonghost_observation_data[network][resolution][matrix][
                        speci
                    ] = file_yearmonths

    # save file tree out to yaml
    with open(
        join(PROVIDENTIA_ROOT, "settings/internal/nonghost_filetree.json"), "w"
    ) as json_file:
        json.dump(nonghost_observation_data, json_file, indent=4)

    return nonghost_observation_data


def get_valid_obs_files_in_date_range(instance, start_date, end_date):
    """
    Filters the full observational data tree to identify files available within a specific date range.

    Parameters
    ----------
    instance : object
        An instance of the application class containing the full observational data tree.
    start_date : str
        The start date in 'YYYYMMDD' format.
    end_date : str
        The end date in 'YYYYMMDD' format.
    """

    # create dictionary to store available observational data
    instance.available_observation_data = {}

    # check if start/end date are valid values, if not, return with no valid obs. files
    if (not valid_date(start_date)) or (not valid_date(end_date)):
        msg = (
            f"One of start date ({start_date}) or end date ({end_date}) are not valid."
        )
        show_message(instance, msg, print=True)
        return

    # check end date is > start date, if not, return with no valid obs. files
    if int(start_date) >= int(end_date):
        msg = f"Start date ({start_date}) exceeds end date ({end_date})."
        show_message(instance, msg, print=True)
        return

    # check start date and end date are both within if valid date range (19000101 - 20500101),
    # if not, return with no valid obs. files
    if (
        (int(start_date) < 19000101)
        or (int(end_date) < 19000101)
        or (int(start_date) >= 20500101)
        or (int(end_date) >= 20500101)
    ):
        msg = (
            f"One of start date ({start_date}) or end date ({end_date}) are not valid."
        )
        show_message(instance, msg, print=True)
        return

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + "01")

    # iterate through observational dictionary
    for network in instance.all_observation_data:
        for resolution in instance.all_observation_data[network]:
            for matrix in instance.all_observation_data[network][resolution]:
                for speci in instance.all_observation_data[network][resolution][matrix]:
                    # get available files
                    file_yearmonths = instance.all_observation_data[network][
                        resolution
                    ][matrix][speci]

                    # get file yearmonths within date range
                    valid_file_yearmonths = sorted(
                        [
                            ym
                            for ym in file_yearmonths
                            if (int("{}01".format(ym)) >= start_date_firstdayofmonth)
                            & (int("{}01".format(ym)) < int(end_date))
                        ]
                    )

                    # add yearmonths to available observation data dict
                    if len(valid_file_yearmonths) > 0:
                        if network not in instance.available_observation_data:
                            instance.available_observation_data[network] = {}
                        if (
                            resolution
                            not in instance.available_observation_data[network]
                        ):
                            instance.available_observation_data[network][
                                resolution
                            ] = {}
                        if (
                            matrix
                            not in instance.available_observation_data[network][
                                resolution
                            ]
                        ):
                            instance.available_observation_data[network][resolution][
                                matrix
                            ] = {}
                        instance.available_observation_data[network][resolution][
                            matrix
                        ][speci] = valid_file_yearmonths


def get_valid_models(instance, start_date, end_date, resolution, networks, species):
    """
    Identifies models within a given date range and chosen networks and species.

    Parameters
    ----------
    instance : object
        An instance of the application class containing directory roots and menu configurations.
    start_date : str
        The start date in 'YYYYMMDD' format.
    end_date : str
        The end date in 'YYYYMMDD' format.
    resolution : str
        The temporal resolution (e.g. 'hourly', 'daily').
    networks : list of str
        The monitoring networks to match against model data.
    species : list of str
        The chemical species or parameters to match against model data.
    """

    # get all different model names (from providentia-interpolation output dir)
    available_models = []
    if os.path.exists(join(instance.mod_root, instance.ghost_version)):
        available_models = os.listdir(
            "%s/%s" % (instance.mod_root, instance.ghost_version)
        )

    # create dictionary to store available model data
    instance.available_model_data = {}

    # list for saving models to add to models pop-up
    models_to_add = []

    # get start date on first of month
    start_date_firstdayofmonth = int(str(start_date)[:6] + "01")

    # iterate through networks and species
    for network, speci in zip(networks, species):
        # iterate through available models
        for model in available_models:
            # get folder where interpolated models are saved
            if "/" not in network:
                files_directory = "%s/%s/%s/%s/%s/%s" % (
                    instance.mod_root,
                    instance.ghost_version,
                    model,
                    resolution,
                    speci,
                    network,
                )
            else:
                files_directory = "%s/%s/%s/%s/%s/%s" % (
                    instance.mod_root,
                    instance.ghost_version,
                    model,
                    resolution,
                    speci,
                    network.replace("/", "-"),
                )

            # test if interpolated directory exists for model
            # if it does not exit, continue
            if not os.path.exists(files_directory):
                continue
            else:
                # get all available netCDF files (handling potential permissions issues)
                try:
                    available_files = os.listdir(files_directory)
                except PermissionError:
                    continue

            # get monthly start date (YYYYMM) of all files
            file_yearmonths = sorted([f.split("_")[-1][:6] for f in available_files])

            # write nested dictionary for model, with associated file yearmonths
            if len(file_yearmonths) > 0:
                # get file yearmonths within date range
                valid_file_yearmonths = sorted(
                    [
                        ym
                        for ym in file_yearmonths
                        if (int("{}01".format(ym)) >= start_date_firstdayofmonth)
                        & (int("{}01".format(ym)) < int(end_date))
                    ]
                )

                # if have valid files, then add model to pop-up menu,
                # and add yearmonths to available model data
                if len(valid_file_yearmonths) > 0:
                    models_to_add.append(model)

                    if network not in instance.available_model_data:
                        instance.available_model_data[network] = {}
                    if resolution not in instance.available_model_data[network]:
                        instance.available_model_data[network][resolution] = {}
                    if speci not in instance.available_model_data[network][resolution]:
                        instance.available_model_data[network][resolution][speci] = {}
                    if (
                        model
                        not in instance.available_model_data[network][resolution][speci]
                    ):
                        instance.available_model_data[network][resolution][speci][
                            model
                        ] = valid_file_yearmonths

    # set list of model names to add on models pop-up
    if instance.mode not in ["report", "library"]:
        models_to_add = np.array(sorted(models_to_add))
        instance.models_menu["models"]["labels"] = models_to_add
        instance.models_menu["models"]["map_vars"] = models_to_add


def get_possible_temporal_resolutions():
    """
    Returns a comprehensive list of all supported temporal resolutions for data processing.

    Returns
    -------
    list
        A list of strings representing temporal resolutions from hourly to annual.
    """

    # define possible temporal resolutions
    possible_resolutions = [
        "hourly",
        "hourly_instantaneous",
        "3hourly",
        "3hourly_instantaneous",
        "6hourly",
        "6hourly_instantaneous",
        "daily",
        "monthly",
        "annual",
    ]

    return possible_resolutions


def get_temporal_resolution_order():
    """
    Provides a dictionary mapping temporal resolution names to their preferred display order in menus.

    Returns
    -------
    dict
        A dictionary where keys are resolution strings and values are integer ranks.
    """

    resolution_order_dict = {
        "hourly": 1,
        "3hourly": 2,
        "6hourly": 3,
        "hourly_instantaneous": 4,
        "3hourly_instantaneous": 5,
        "6hourly_instantaneous": 6,
        "daily": 7,
        "monthly": 8,
    }

    return resolution_order_dict


def get_possible_resampling_resolutions(resolution, daily_forecast=False):
    """
    Identifies permissible lower temporal resolutions for resampling based on the base resolution.

    Parameters
    ----------
    resolution : str
        The current temporal resolution of the data (e.g. 'hourly', 'daily').
    daily_forecast : bool, optional
        Flag indicating if the data belongs to a daily forecast cycle (default is False).

    Returns
    -------
    list of str
        A list of target resolutions that are coarser than or equal to the input resolution.
    """

    if daily_forecast:
        if resolution in ["hourly", "hourly_instantaneous"]:
            resolutions = ["hourly", "3hourly", "6hourly", "daily"]
        elif resolution in ["3hourly", "3hourly_instantaneous"]:
            resolutions = ["3hourly", "6hourly", "daily"]
        elif resolution in ["6hourly", "6hourly_instantaneous"]:
            resolutions = ["6hourly", "daily"]
        elif resolution == "daily":
            resolutions = ["daily"]
        elif resolution == "monthly":
            resolutions = []
        elif resolution == "annual":
            resolutions = []
    else:
        if resolution in ["hourly", "hourly_instantaneous"]:
            resolutions = ["hourly", "3hourly", "6hourly", "daily", "monthly", "annual"]
        elif resolution in ["3hourly", "3hourly_instantaneous"]:
            resolutions = ["3hourly", "6hourly", "daily", "monthly", "annual"]
        elif resolution in ["6hourly", "6hourly_instantaneous"]:
            resolutions = ["6hourly", "daily", "monthly", "annual"]
        elif resolution == "daily":
            resolutions = ["daily", "monthly", "annual"]
        elif resolution == "monthly":
            resolutions = ["monthly", "annual"]
        elif resolution == "annual":
            resolutions = ["annual"]

    return resolutions


def get_periodic_relevant_temporal_resolutions(resolution):
    """
    Identifies the appropriate time cycles for periodic analysis based on a base temporal resolution.

    Parameters
    ----------
    resolution : str
        The name of the selected temporal resolution (e.g. 'hourly', 'daily').

    Returns
    -------
    relevant_temporal_resolutions : list
        A list of periodic groupings (e.g. 'hour', 'dayofweek') valid for the input resolution.
    """

    if "hourly" in resolution:
        relevant_temporal_resolutions = ["hour", "dayofweek", "month"]
    elif resolution == "daily":
        relevant_temporal_resolutions = ["dayofweek", "month"]
    elif resolution == "monthly":
        relevant_temporal_resolutions = ["month"]
    else:
        relevant_temporal_resolutions = []

    return relevant_temporal_resolutions


def get_periodic_nonrelevant_temporal_resolutions(resolution):
    """
    Identifies time cycles that cannot be analysed based on a base temporal resolution.

    Parameters
    ----------
    resolution : str
        The name of the selected temporal resolution (e.g. 'hourly', 'daily').

    Returns
    -------
    nonrelevant_temporal_resolutions : list
        A list of periodic groupings that are scientifically invalid for the input resolution.
    """

    if "hourly" in resolution:
        nonrelevant_temporal_resolutions = []
    elif resolution == "daily":
        nonrelevant_temporal_resolutions = ["hour"]
    elif resolution == "monthly":
        nonrelevant_temporal_resolutions = ["hour", "dayofweek"]
    else:
        nonrelevant_temporal_resolutions = ["hour", "dayofweek", "month"]

    return nonrelevant_temporal_resolutions


def valid_date(date_text):
    """
    Validates whether a given date string or integer conforms to the 'YYYYMMDD' format.

    Parameters
    ----------
    date_text : str or int
        The date value to be validated.

    Returns
    -------
    bool
        True if the date is valid and correctly formatted, False otherwise.
    """

    try:
        datetime.datetime.strptime(str(date_text), "%Y%m%d")
        return True
    except Exception:
        return False


def generate_file_trees(instance):
    """
    Handles the dynamic generation or loading of observational data catalogues.

    Parameters
    ----------
    instance : object
        An instance of the application class containing configuration flags and versioning.
    """

    # get dictionaries of observational GHOST and non-GHOST filetrees, either created dynamically or loaded
    # if have filetree flags, then these overwrite any defaults
    gft = False
    if instance.generate_file_tree:
        gft = True
    elif instance.disable_file_tree:
        gft = False
    # by default generate filetree on MN5
    elif instance.machine in ["mn5", "nord4"]:
        gft = True
    # by default generate filetree locally
    elif instance.filetree_type == "local":
        gft = True

    # generate file trees
    ghost_filetree_path = join(
        PROVIDENTIA_ROOT,
        "settings/internal/ghost_filetree_{}.json".format(instance.ghost_version),
    )
    nonghost_filetree_path = join(
        PROVIDENTIA_ROOT, "settings/internal/nonghost_filetree.json"
    )

    # generate file trees for ghost
    if gft or (not os.path.exists(ghost_filetree_path)):
        instance.logger.info("")
        if not os.path.exists(ghost_filetree_path):
            instance.logger.info(f"Generating file tree {ghost_filetree_path}...")
        else:
            instance.logger.info(f"Updating file tree {ghost_filetree_path}...")
        instance.all_observation_data = get_ghost_observational_tree(instance)
    # load file trees
    else:
        instance.logger.info(f"Loading file tree {ghost_filetree_path}...")
        try:
            instance.all_observation_data = json.load(
                open(
                    join(
                        PROVIDENTIA_ROOT,
                        "settings/internal/ghost_filetree_{}.json".format(
                            instance.ghost_version
                        ),
                    )
                )
            )
        except FileNotFoundError:
            error = "Error: Trying to load 'settings/internal/ghost_filetree_{}.json' but file does not exist. Run with the flag '--gft' to generate this file.".format(
                instance.ghost_version
            )
            instance.logger.error(error)
            sys.exit(1)

    if instance.nonghost_root is not None:
        # generate file trees for nonghost
        if gft or (not os.path.exists(nonghost_filetree_path)):
            if not os.path.exists(nonghost_filetree_path):
                instance.logger.info(
                    f"Generating file tree {nonghost_filetree_path}..."
                )
            else:
                instance.logger.info(f"Updating file tree {nonghost_filetree_path}...")
            nonghost_observation_data = get_nonghost_observational_tree(instance)
        # load file trees
        else:
            instance.logger.info(f"Loading file tree {nonghost_filetree_path}...")
            try:
                nonghost_observation_data = json.load(
                    open(
                        join(
                            PROVIDENTIA_ROOT, "settings/internal/nonghost_filetree.json"
                        )
                    )
                )
            except FileNotFoundError:
                error = "Error: Trying to load 'settings/internal/nonghost_filetree.json' but file does not exist. Run with the flag '--gft' to generate this file."
                instance.logger.error(error)
                sys.exit(1)

        # merge GHOST and non-GHOST filetrees
        instance.all_observation_data = {
            **instance.all_observation_data,
            **nonghost_observation_data,
        }


def get_valid_metadata(read_instance, variable, valid_station_idxs, networkspeci):
    """
    Extracts the first non-NaN metadata values for selected stations across all time steps.

    Parameters
    ----------
    read_instance : object
        An instance of the application class responsible for data reading operations.
    variable : str
        The specific metadata field to retrieve (e.g., 'station_name', 'latitude').
    valid_station_idxs : list of int
        Indices of the stations to be processed.
    networkspeci : str
        A combined string representing the network and species identifier.

    Returns
    -------
    valid_metadata : list
        A list of non-null metadata values, one per station.
    """

    stations_metadata = read_instance.canvas_instance.selected_station_metadata[
        networkspeci
    ][variable][valid_station_idxs, :]
    valid_metadata = []
    for station_metadata in stations_metadata:
        first_valid_station_metadata = next(
            (value for value in station_metadata if value == value), "nan"
        )
        # if metadata returns an array, get first
        if isinstance(first_valid_station_metadata, (np.ndarray)):
            first_valid_station_metadata = first_valid_station_metadata.item()
        valid_metadata.append(first_valid_station_metadata)

    return valid_metadata


def _is_leap_year(year):
    """
    Determine whether a given year is a leap year under the standard Gregorian rule.

    This helper is used when constructing monthly timestamps for calendars that
    require explicit month-length handling.

    Parameters
    ----------
    year : int
        Calendar year to test.

    Returns
    -------
    bool
        True if the year is a leap year under the standard Gregorian rule,
        otherwise False.
    """
    return (year % 4 == 0) and ((year % 100 != 0) or (year % 400 == 0))


def _days_in_month(year, month, calendar):
    """
    Return the number of days in a month for a specified calendar.

    The helper supports both standard/Gregorian-like calendars and non-standard
    netCDF calendars such as ``360_day``, ``365_day`` / ``noleap``, and
    ``366_day`` / ``all_leap``. It is primarily used when constructing dates
    from ``months since ...`` offsets.

    Parameters
    ----------
    year : int
        Calendar year of the target month.
    month : int
        Month number in the range 1-12.
    calendar : str
        NetCDF calendar name. Examples include ``'standard'``, ``'gregorian'``,
        ``'proleptic_gregorian'``, ``'360_day'``, ``'365_day'``, ``'noleap'``,
        ``'366_day'``, or ``'all_leap'``.

    Returns
    -------
    int
        Number of days in the requested month under the specified calendar.
    """
    if calendar == "360_day":
        return 30

    if month in (1, 3, 5, 7, 8, 10, 12):
        return 31
    if month in (4, 6, 9, 11):
        return 30

    # February
    if calendar in ("all_leap", "366_day"):
        return 29
    if calendar in ("noleap", "365_day"):
        return 28

    return 29 if _is_leap_year(year) else 28


def _parse_origin_parts(origin_str):
    """
    Parse a netCDF time-origin string into its date and time components.

    Parsing is delegated to ``pandas.Timestamp`` so that common ISO-like and
    netCDF-style origin strings are handled consistently.

    Parameters
    ----------
    origin_str : str
        Origin timestamp string extracted from a netCDF ``units`` attribute,
        for example ``'2018-12-01'`` or ``'2018-12-01 00:00:00'``.

    Returns
    -------
    year : int
        Parsed year component.
    month : int
        Parsed month component.
    day : int
        Parsed day component.
    hour : int
        Parsed hour component.
    minute : int
        Parsed minute component.
    second : int
        Parsed second component.
    microsecond : int
        Parsed microsecond component.
    """
    ts = pd.Timestamp(origin_str)
    return ts.year, ts.month, ts.day, ts.hour, ts.minute, ts.second, ts.microsecond


def _add_months_to_ym(year, month, delta_months):
    """
    Add an integer number of months to a ``(year, month)`` pair.

    This helper is used when converting ``months since ...`` time coordinates
    into concrete calendar dates.

    Parameters
    ----------
    year : int
        Starting year.
    month : int
        Starting month in the range 1-12.
    delta_months : int
        Integer number of months to add. May be positive, zero, or negative.

    Returns
    -------
    new_year : int
        Year after the month offset has been applied.
    new_month : int
        Month after the month offset has been applied, in the range 1-12.
    """
    total = year * 12 + (month - 1) + int(delta_months)
    new_year = total // 12
    new_month = (total % 12) + 1
    return new_year, new_month


def _num2date_months(values, origin_str, calendar):
    """
    Convert ``months since ...`` offsets into calendar-aware datetime objects.

    This helper explicitly handles monthly data, including non-standard netCDF
    calendars. The function assumes that the monthly offsets are integer-valued,
    which is the normal case for monthly data files. If the day-of-month in the
    origin timestamp is not valid in a target month, it is clamped to the last
    valid day of that month for the chosen calendar.

    Parameters
    ----------
    values : numpy.ndarray
        Numeric month offsets read from the netCDF ``time`` variable. These are
        expected to be integer month counts.
    origin_str : str
        Origin timestamp string extracted from the netCDF ``units`` attribute.
    calendar : str
        NetCDF calendar name associated with the time variable.

    Returns
    -------
    dates : list
        List of calendar-aware datetime objects, typically instances of the
        appropriate ``cftime`` datetime class for the requested calendar.
    """
    values = np.asarray(values, dtype=np.float64)

    if not np.allclose(values, np.round(values), equal_nan=False):
        raise ValueError(
            f"Monthly time coordinate contains non-integer month offsets: {values}"
        )

    offsets = np.round(values).astype(int)

    year, month, day, hour, minute, second, microsecond = _parse_origin_parts(
        origin_str
    )
    cls = _CFTIME_CLASS.get(calendar, cftime.DatetimeGregorian)

    out = []
    for delta_months in offsets:
        yy, mm = _add_months_to_ym(year, month, delta_months)
        dd = min(day, _days_in_month(yy, mm, calendar))
        out.append(cls(yy, mm, dd, hour, minute, second, microsecond))

    return out


def _seconds_to_asi8(seconds):
    """
    Convert integer Unix seconds to int64 nanosecond timestamps.

    This helper is used to preserve the historical Providentia behavior of
    working at whole-second resolution while still returning timestamps in the
    same integer-nanosecond format as ``DatetimeIndex.asi8``.

    Parameters
    ----------
    seconds : array-like
        Integer or float-like Unix seconds since 1970-01-01 00:00:00. Values
        are converted to int64 seconds before scaling to nanoseconds.

    Returns
    -------
    timestamps : numpy.ndarray
        One-dimensional array of int64 nanosecond timestamps.
    """
    seconds = np.asarray(seconds, dtype=np.int64)
    return seconds * 1_000_000_000


def _dates_to_asi8_seconds(dates, calendar):
    """
    Convert datetime-like objects to int64 nanosecond timestamps at second resolution.

    This helper intentionally discards sub-second precision. That behavior is
    designed to match the historical Providentia logic, which explicitly removed
    microseconds or cast decoded dates to ``datetime64[s]`` before obtaining
    integer timestamps. This avoids problems caused by tiny floating-point artefacts
    in decoded time coordinates.

    Parameters
    ----------
    dates : sequence
        Sequence of datetime-like objects. This may contain standard Python
        ``datetime`` objects, NumPy datetime-like objects, pandas-compatible
        datetime values, or ``cftime`` datetime instances.
    calendar : str
        NetCDF calendar name used to interpret ``cftime`` dates correctly.

    Returns
    -------
    timestamps : numpy.ndarray
        One-dimensional array of int64 nanosecond timestamps corresponding to
        whole-second times only. Any sub-second precision in the input is
        intentionally discarded.
    """
    dates = np.asarray(dates, dtype=object)

    if dates.size == 0:
        return np.array([], dtype=np.int64)

    first = dates.flat[0]

    # ------------------------------------------------------
    # cftime path (works for DatetimeGregorian too)
    # ------------------------------------------------------
    if isinstance(first, cftime.datetime):
        seconds = cftime.date2num(
            dates.tolist(), units="seconds since 1970-01-01 00:00:00", calendar=calendar
        )
        seconds = np.floor(np.asarray(seconds, dtype=np.float64)).astype(np.int64)
        return _seconds_to_asi8(seconds)

    # ------------------------------------------------------
    # Standard datetime-like path
    # ------------------------------------------------------
    dt64_s = pd.to_datetime(dates).to_numpy(dtype="datetime64[s]")
    seconds = dt64_s.astype(np.int64)
    return _seconds_to_asi8(seconds)


def time_var_to_asi8(time_var):
    """
    Convert a netCDF ``time`` variable to integer nanosecond timestamps.

    The function supports the common resolutions used in Providentia:
    hourly, 3-hourly, 6-hourly, daily, and monthly. It also respects the
    ``calendar`` attribute on the netCDF time variable.

    A fast vectorized conversion path is used for standard/Gregorian-like
    calendars when the units are expressed in seconds, minutes, hours, or
    days since a standard origin timestamp. Monthly data and non-standard
    calendars automatically fall back to a safer calendar-aware conversion path.

    IMPORTANT:
    ----------
    This function intentionally normalizes all timestamps to whole-second
    precision. That mirrors the historical Providentia behavior and avoids
    issues caused by tiny floating-point artefacts in decoded netCDF times.

    Parameters
    ----------
    time_var : netCDF4.Variable
        NetCDF time variable. The variable must contain numeric offsets and a
        ``units`` attribute of the form ``'<unit> since <origin>'``. If present,
        the ``calendar`` attribute is also used to ensure correct interpretation
        of the time axis.

    Returns
    -------
    file_timestamp : numpy.ndarray
        One-dimensional array of int64 nanosecond timestamps in Unix-epoch
        nanoseconds (``asi8`` style), suitable for comparison against
        pandas/NumPy ``datetime64[ns]``-based integer time axes.
    """
    units = time_var.units.strip()
    calendar = getattr(time_var, "calendar", "standard")
    values = np.asarray(time_var[:], dtype=np.float64)

    unit_part, origin_part = units.split("since", 1)
    unit = unit_part.strip().lower()
    origin_str = origin_part.strip()

    # ------------------------------------------------------
    # FAST PATH for standard/Gregorian-like calendars
    # and second/minute/hour/day units
    #
    # We intentionally normalize to whole seconds.
    # ------------------------------------------------------
    standard_cals = {"standard", "gregorian", "proleptic_gregorian"}

    sec_per_unit = {
        "second": 1,
        "seconds": 1,
        "sec": 1,
        "secs": 1,
        "minute": 60,
        "minutes": 60,
        "min": 60,
        "mins": 60,
        "hour": 3600,
        "hours": 3600,
        "hr": 3600,
        "hrs": 3600,
        "day": 86400,
        "days": 86400,
    }

    if calendar in standard_cals and unit in sec_per_unit:
        try:
            origin_s = np.datetime64(origin_str.replace(" ", "T"), "s").astype(np.int64)
            delta_s = np.floor(values * sec_per_unit[unit]).astype(np.int64)
            return _seconds_to_asi8(origin_s + delta_s)
        except Exception:
            # Fall through to safe path if the fast path fails to parse
            pass

    # ------------------------------------------------------
    # MONTHLY path
    # ------------------------------------------------------
    if unit in ("month", "months"):
        dates = _num2date_months(values, origin_str, calendar)
        return _dates_to_asi8_seconds(dates, calendar)

    # ------------------------------------------------------
    # SAFE GENERAL PATH for everything else
    # ------------------------------------------------------
    dates = num2date(values, units=units, calendar=calendar)
    return _dates_to_asi8_seconds(dates, calendar)
