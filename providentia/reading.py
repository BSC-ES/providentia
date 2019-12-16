from netCDF4 import Dataset
import numpy as np


def read_netcdf_station_information(tuple_to_read):
    """Function that handles reading of observational
    desired station metadata in a netCDF file,
    returning a dictionary with read metadata.

    Args:
        tuple_to_read(tuple): file, [metadata_vars_to_read]
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
    """Function that returns numpy object of lists
    of station data with NaNs removed
    """

    # reshape numpy array to have lists of data per station
    data = data.tolist()
    # iterate through each list of station data and remove NaNs
    for station_ii, station_data in enumerate(data):
        data[station_ii] = np.array(station_data)[~np.isnan(station_data)]
    # return numpy object of lists of station data with NaNs removed
    return np.array(data)


# TODO: add def read_netcdf_data(relevant_file):
