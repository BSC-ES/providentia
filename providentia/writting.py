""" Module storing writting functions """
import numpy as np


def write_data_npz(mpl_canvas, fname):
    """Function that writes out current data in memory to .npy file"""

    np.savez(fname, data=mpl_canvas.read_instance.data_in_memory_filtered,
             metadata=mpl_canvas.read_instance.metadata_in_memory)


def write_netcdf():
    """Write data and metadata to netcdf file"""
    return
