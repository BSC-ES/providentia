from netCDF4 import Dataset
import numpy as np
import pandas as pd


def read_netcdf_station_information(tuple_to_read):

    """ Function that handles reading of observational
    desired station metadata in a netCDF file,
    returning a dictionary with read metadata.

    Parameters
    ----------
    tuple_to_read: tuple of the form (file_to_read, [metadata_vars_to_read])

    """

    # read netCDF frame
    ncdf_root = Dataset(tuple_to_read[0])

    # read all desired metadata, placing it within a dictionary by variable name
    read_metadata = {}
    for meta_var in tuple_to_read[1]:
        read_metadata[meta_var] = ncdf_root[meta_var][:]

    # close netCDF
    ncdf_root.close()

    return read_metadata


def drop_nans(data):

    """function that returns numpy object of lists of station data with NaNs removed"""

    # reshape numpy array to have lists of data per station
    data = data.tolist()
    # iterate through each list of station data and remove NaNs
    for station_ii, station_data in enumerate(data):
        data[station_ii] = np.array(station_data)[~np.isnan(station_data)]
    # return numpy object of lists of station data with NaNs removed
    return np.array(data)


# def read_netcdf_data(relevant_file):
#
#     """ Function that handles reading of observational/experiment
#     netCDF data also handles filtering of observational data based
#     on selected qa/flag/classification flags."""
#
#     # read netCDF frame
#     ncdf_root = netCDF4.Dataset(relevant_file)
#
#     # get time units
#     time_units = ncdf_root['time'].units
#
#     # get file time (handle monthly resolution data differently to hourly/daily
#     # as num2date does not support 'months since' units)
#     if 'months' in time_units:
#         monthly_start_date = time_units.split(' ')[2]
#         file_time = pd.date_range(start=monthly_start_date, periods=1, freq='MS')
#     else:
#         file_time = netCDF4.num2date(ncdf_root['time'][:], time_units)
#         # remove microseconds
#         file_time = pd.to_datetime([t.replace(microsecond=0) for t in file_time])
#
#     # get valid file time indices (i.e. those times in active full time array)
#     valid_file_time_indices = \
#         np.array([i for i, val in enumerate(file_time)
#                   if (val >= time_array[0]) & (val <= time_array[-1])], dtype=np.int)
#     # cut file time for valid indices
#     file_time = file_time[valid_file_time_indices]
#     # get indices relative to active full time array
#     full_array_time_indices = np.searchsorted(time_array, file_time)
#
#     # get all station references in file
#     file_station_references = ncdf_root['station_reference'][:]
#     # get indices of all unique station references that are contained within file station references array
#     full_array_station_indices = np.where(np.in1d(station_references, file_station_references))[0]
#     # get indices of file station station references that are contained in all unique station references array
#     current_file_station_indices = np.where(np.in1d(file_station_references, station_references))[0]
#
#     # read in species data
#     file_data = ncdf_root[active_species][:, valid_file_time_indices]
#     # get masked data
#     data_mask = file_data.mask
#     # set masked data as NaN
#     file_data[data_mask] = np.NaN
#
#     # for observations, set species data based on selected qa flags/standard data provider
#     # flags/classifications to retain or remove as NaN
#     if process_type == 'observations':
#         # if some qa flags selected then screen
#         if len(selected_qa) > 0:
#             # screen out observations which are associated with any of the selected qa flags
#             file_data[np.isin(ncdf_root['qa'][:, valid_file_time_indices, :], selected_qa).any(axis=2)] = np.NaN
#         # if some data provider flags selected then screen
#         if len(selected_flags) > 0:
#             # screen out observations which are associated with any of the selected data provider flags
#             file_data[np.isin(ncdf_root['flag'][:, valid_file_time_indices, :], selected_flags).any(axis=2)] = np.NaN
#         # if some classification flags (retain  or remove) selected then screen
#         if (len(selected_classifications_to_retain) > 0) or (len(selected_classifications_to_remove) > 0):
#             file_classifications = ncdf_root['classification'][:, valid_file_time_indices, :]
#             # screen out all observations that aren't associated with all of the selected classifications to retain
#             if len(selected_classifications_to_retain) > 0:
#                 file_data[np.isin(file_classifications, selected_classifications_to_retain,
#                                   invert=True).all(axis=2)] = np.NaN
#             # screen out all observations that are associated with any of the selected classifications to remove
#             if len(selected_classifications_to_remove) > 0:
#                 file_data[np.isin(file_classifications, selected_classifications_to_remove).any(axis=2)] = np.NaN
#
#     # close netCDF
#     ncdf_root.close()
#
#     # return valid species data, time indices relative to active full time array,
#     # file station indices relative to all unique station references array
#     return file_data[current_file_station_indices, :], full_array_time_indices, full_array_station_indices
