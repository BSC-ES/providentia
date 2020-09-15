import os

from PyQt5 import QtWidgets
import matplotlib
from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from .prov_canvas import MPLCanvas


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

    def save_data_button(self):
        filetypes = {'Numpy file': ['npz']}
        sorted_filetypes = sorted(filetypes.items())
        startpath = os.path.expanduser(matplotlib.rcParams['savefig.directory'])
        start = os.path.join(startpath, 'default_filename')
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(['*.%s' % ext for ext in exts])
            filter = '%s (%s)' % (name, exts_list)
            filters.append(filter)
            filters = ';;'.join(filters)
            fname, filter = qt_compat._getSaveFileName(self.parent(), "Choose a filename to save to",
                                                       start, filters, selectedFilter)
            if fname:
                # Save dir for next time, unless empty str (i.e., use cwd).
                if startpath != "":
                    matplotlib.rcParams['savefig.directory'] = (os.path.dirname(fname))
                    try:
                        MPLCanvas.write_data_npz(self, fname)
                    except Exception as e:
                        QtWidgets.QMessageBox.critical(self, "Error saving file", str(e), QtWidgets.QMessageBox.Ok,
                                                       QtWidgets.QMessageBox.NoButton)
