""" Navigation toolbar and buttons/options functions"""
import os
import configparser

from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import QFileDialog
import matplotlib
from matplotlib.backends import qt_compat
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT

from .writing import export_data_npz, export_netcdf, export_configuration
from .aux import representativity_conf, period_conf, metadata_conf, load_conf

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

def save_data(canvas_instance):
    """Pops window for choosing directory, filename and type
    for saving data, metadata and configuration.
    Available filetypes: Numpy file: .npz, netCDF: .nc"""

    filetypes = {'NetCDF': 'nc', 'Numpy file': 'npz', 'Configuration': 'conf'}
    sorted_filetypes = sorted(filetypes.items())
    startpath = os.path.expanduser(matplotlib.rcParams['savefig.directory'])
    daterange = canvas_instance.read_instance.le_start_date.text() + "_" \
                + canvas_instance.read_instance.le_end_date.text()
    try:
        eg_name = "PRV_" + canvas_instance.read_instance.species + "_" + daterange
    except:
        eg_name = "default_filename"
    start = os.path.join(startpath, eg_name)

    filter_ext = ['%s (%s)' % (name, '*.%s' % ext) for name, ext in sorted_filetypes]
    filter_ext = ';;'.join(filter_ext)
    # prompt with file name and extension
    fname, fext = qt_compat._getSaveFileName(None, "Choose a filename to save to", start, filter_ext)
    chose_npz = "npz" in fext
    chose_conf = "conf" in fext
    if fname:
        # Save dir for next time, unless empty str (i.e., use cwd).
        if startpath != "":
            matplotlib.rcParams['savefig.directory'] = (os.path.dirname(fname))
            try:
                QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
                if chose_npz:
                    export_data_npz(canvas_instance, fname)
                elif chose_conf:
                    export_configuration(canvas_instance.read_instance, fname)
                else:
                    export_netcdf(canvas_instance, fname)
                QtWidgets.QApplication.restoreOverrideCursor()
            except Exception as e:
                QtWidgets.QMessageBox.critical(canvas_instance, "Error saving file", str(e),
                                               QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)

def conf_dialogs(instance):
    """Pops window for selecting configuration file. If file selcted, pops an
    input dialog for the user to select which section wants to load. Calls
    reload_conf to reset the fields"""

    conf_to_load = filename_dialog(instance)
    # is user pressed cancel
    if conf_to_load is None:
        return

    try:
        aux.load_conf(fpath=conf_to_load)
        all_sections = instance.sub_opts.keys()
        selected_section, okpressed = QtWidgets.QInputDialog.getItem(instance, 'Sections',
                                                                     'Select section to load',  
                                                                     all_sections, 0, False)
        if okpressed:
            reload_conf(instance, selected_section, conf_to_load)
    except Exception as e:
        QtWidgets.QMessageBox.critical(instance, "Error loading configuration file",
                                       str(e), QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)

def filename_dialog(instance):
    options = QFileDialog.Options()
    options |= QFileDialog.DontUseNativeDialog
    filename, _ = QFileDialog.getOpenFileName(instance, "Choose configuration file to load", "",
                                              "All Files (*);;Python Files (*.py)", options=options)
    if filename:
        return filename

def reload_conf(instance, section, fpath):
    """Resets previous selections, fills values according to new conf file,
    reads and filters."""

    # delete attributes from default config file
    CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))
    dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))
    dconf = configparser.RawConfigParser(empty_lines_in_values=False)
    dconf.read(dconf_path)
    if instance.from_conf:
        for k in dconf['default'].keys():
            delattr(instance, k)

    # update config and section attributes of instance
    instance.config = fpath
    instance.section = section
    instance.from_conf = True
    vars(instance).update({(k, instance.parse_parameter(k, val)) for k, val in instance.sub_opts[section].items()})

    # update species, experiments, qa & flags
    instance.config_bar_initialisation = True
    instance.update_configuration_bar_fields()
    instance.config_bar_initialisation = False

    # read
    instance.handle_data_selection_update()

    # reset the filter fields 
    instance.reset_options()

    # set fields available for filtering
    aux.representativity_conf(instance)
    aux.period_conf(instance)
    aux.metadata_conf(instance)
        
    # filter
    instance.mpl_canvas.handle_data_filter_update()