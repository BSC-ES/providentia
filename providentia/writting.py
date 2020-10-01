""" Module storing writting functions """
import os

import numpy as np
from netCDF4 import Dataset


def export_data_npz(mpl_canvas, fname):
    """Function that writes out current data in memory to .npy file"""

    np.savez(fname, data=mpl_canvas.read_instance.data_in_memory_filtered,
             metadata=mpl_canvas.read_instance.metadata_in_memory)


def export_netcdf(mpl_canvas, fname):
    """Write data and metadata to netcdf file"""

    fillval = -9999.

    # nc_dims = {
    #     'time': None,
    #     'station': mpl_canvas.read_instance.station_references.size,
    #     'strlen': 75,
    # }

    nc_vars = {
        mpl_canvas.read_instance.active_species: {
            'dims': ('time', 'station'),
            'dtype': 'f',
            'attrs': {
                'missing_value': fillval,
                'units': '-',
                'long_name': "",
            },
            'val': None,
        }
    }

    # standard variables
    st_vars = {
        'lon': {
            'dims': ('station',),
            'dtype': 'f',
            'attrs': {
                'standard_name': "longitude",
                'units': "degrees_east",
            },
            'val': None,
        },
        'lat': {
            'dims': ('station',),
            'dtype': 'f',
            'attrs': {
                'standard_name': "latitude",
                'units': "degrees_north",
            },
            'val': None,
        },
        'altitude': {
            'dims': ('station',),
            'dtype': 'f',
            'attrs': {
                'standard_name': "altitude",
                'units': "meters",
            },
            'val': None,
        },
        'time': {
            'dims': ('time',),
            'dtype': 'f',
            'attrs': {
                'units': 'hours since %s',
            },
            'val': None,
        },
    }

    # let's say that we have the file prepared, more or less
    # start file
    fout = Dataset(fname, 'w', format="NETCDF4")

    fout.close()

    return
