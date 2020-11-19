""" Navigation toolbar and buttons/options functions"""
import os

from PyQt5 import QtWidgets
import matplotlib
from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from .writting import export_data_npz, export_netcdf


class NavigationToolbar(NavigationToolbar2QT):
    """Define class that updates available buttons on matplotlib toolbar"""

    # only display wanted buttons
    NavigationToolbar2QT.toolitems = (
        ('Home', 'Reset original view', 'home', 'home'),
        ('Back', 'Back to previous view', 'back', 'back'),
        ('Forward', 'Forward to next view', 'forward', 'forward'),
        (None, None, None, None),
        ('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'),
        ('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'),
        (None, None, None, None),
        ('Save', 'Save the figure', 'filesave', 'save_figure'),
        (None, None, None, None)
    )


def save_data(mpl_canvas):
    """Pops window for choosing directory, filename and type
    for saving data, metadata and configuration.
    Available filetypes: Numpy file: .npz, netCDF: .nc"""

    filetypes = {'Numpy file': 'npz', 'NetCDF': 'nc'}
    sorted_filetypes = sorted(filetypes.items())
    startpath = os.path.expanduser(matplotlib.rcParams['savefig.directory'])
    daterange = mpl_canvas.read_instance.le_start_date.text() + "_" \
                + mpl_canvas.read_instance.le_end_date.text()
    try:
        eg_name = "PRV_" + mpl_canvas.read_instance.active_species + "_" + daterange
    except:
        eg_name = "default_filename"
    start = os.path.join(startpath, eg_name)

    filter_ext = ['%s (%s)' % (name, '*.%s' % ext) for name, ext in sorted_filetypes]
    filter_ext = ';;'.join(filter_ext)
    # prompt with file name and extension
    fname, fext = qt_compat._getSaveFileName(None, "Choose a filename to save to", start, filter_ext)
    chose_npz = "npz" in fext
    if fname:
        # Save dir for next time, unless empty str (i.e., use cwd).
        if startpath != "":
            matplotlib.rcParams['savefig.directory'] = (os.path.dirname(fname))
            try:
                if chose_npz:
                    export_data_npz(mpl_canvas, fname)
                else:
                    export_netcdf(mpl_canvas, fname)
            except Exception as e:
                QtWidgets.QMessageBox.critical(mpl_canvas, "Error saving file", str(e),
                                               QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)
