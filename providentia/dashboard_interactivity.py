""" Objects and functions to interact with the axes """

import copy
import datetime
import math
import time

from matplotlib.dates import num2date
from matplotlib.lines import Line2D
from matplotlib.widgets import _SelectorWidget
import numpy as np
from packaging.version import Version
import pandas as pd
from PyQt5.QtWidgets import QToolTip, QWidget

from .dashboard_elements import set_formatting
from .plot_aux import get_map_extent, get_hex_code
from .plot_formatting import harmonise_xy_lims_paradigm


class LassoSelector(_SelectorWidget):
    """
    Selection curve of an arbitrary shape.
    For the selector to remain responsive you must keep a reference to it.
    The selected path can be used in conjunction with `~.Path.contains_point`
    to select data points from an image.
    In contrast to `Lasso`, `LassoSelector` is written with an interface
    similar to `RectangleSelector` and `SpanSelector`, and will continue to
    interact with the Axes until disconnected.
    Example usage::
        ax = plt.subplot()
        ax.plot(x, y)
        def onselect(verts):
            print(verts)
        lasso = LassoSelector(ax, onselect)
    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The parent Axes for the widget.
    onselect : function
        Whenever the lasso is released, the *onselect* function is called and
        passed the vertices of the selected path.
    useblit : bool, default: True
        Whether to use blitting for faster drawing (if supported by the
        backend). See the tutorial :doc:`/tutorials/advanced/blitting`
        for details.
    props : dict, optional
        Properties with which the line is drawn, see `matplotlib.lines.Line2D`
        for valid properties. Default values are defined in ``mpl.rcParams``.
    button : `.MouseButton` or list of `.MouseButton`, optional
        The mouse buttons used for rectangle selection.  Default is ``None``,
        which corresponds to all buttons.
    """

    def __init__(self, ax, onselect, useblit=True, props=None, button=None):

        super().__init__(ax, onselect, useblit=useblit, button=button)
        self.verts = None
        props = {
            **(props if props is not None else {}),
            # Note that self.useblit may be != useblit, if the canvas doesn't
            # support blitting.
            'animated': self.useblit, 'visible': False,
        }
        line = Line2D([], [], **props)
        self.ax.add_line(line)
        self._selection_artist = line

        return None
    
    def _press(self, event):
            
        self.verts = [self._get_data(event)]
        self._selection_artist.set_visible(True)

        return None
    
    def _onmove(self, event):

        if self.verts is None:
            return
        self.verts.append(self._get_data(event))
        self._selection_artist.set_data(list(zip(*self.verts)))

        self.update()

        return None

    def _release(self, event):

        if self.verts is not None:
            self.verts.append(self._get_data(event))
            self.onselect(self.verts)
        self._selection_artist.set_data([[], []])
        self._selection_artist.set_visible(False)
        self.verts = None

        return None

def zoom_map_func(canvas_instance, event):
    """ Function to handle zoom on map using scroll wheel. """

    if event.inaxes == canvas_instance.plot_axes['map']:

        # lock canvas drawing if can, else return
        if canvas_instance.figure.canvas.widgetlock.locked():
            if not canvas_instance.figure.canvas.widgetlock.isowner(canvas_instance.zoom_map_event):
                return None
        else:        
            canvas_instance.figure.canvas.widgetlock(canvas_instance.zoom_map_event)

        # get the current x and y limits
        current_xlim = canvas_instance.plot_axes['map'].get_xlim()
        current_ylim = canvas_instance.plot_axes['map'].get_ylim()

        # get position of cursor
        xdata = event.xdata
        ydata = event.ydata
        base_scale = canvas_instance.plot_characteristics['map']['base_scale']

        if event.button == 'up':
            # deal with zoom in
            scale_factor = base_scale
        elif event.button == 'down':
            # deal with zoom out
            scale_factor = 1/base_scale
        else:
            # exceptions
            scale_factor = 1
        
        if event.button == 'up' or event.button == 'down':
            
            # set new limits
            canvas_instance.plot_axes['map'].set_xlim([xdata - (xdata - current_xlim[0]) / scale_factor, 
                                                        xdata + (current_xlim[1] - xdata) / scale_factor])
            canvas_instance.plot_axes['map'].set_ylim([ydata - (ydata - current_ylim[0]) / scale_factor, 
                                                        ydata + (current_ylim[1] - ydata) / scale_factor])
            
            # save map extent (in data coords)
            canvas_instance.read_instance.map_extent = get_map_extent(canvas_instance)
            
            # draw changes
            canvas_instance.figure.canvas.draw_idle()
        
            # update buttons (previous-forward) history
            canvas_instance.read_instance.navi_toolbar.push_current()
            canvas_instance.read_instance.navi_toolbar.set_history_buttons()

        # unlock canvas drawing
        if canvas_instance.figure.canvas.widgetlock.isowner(canvas_instance.zoom_map_event):
            canvas_instance.figure.canvas.widgetlock.release(canvas_instance.zoom_map_event)

    return None

def picker_block_func(canvas_instance, event):
    """ Block or unblock the station and legend pick functions
        to avoid interferences.
    """

    if event.inaxes == canvas_instance.plot_axes['legend']:
        # unblock legend picker in legend
        canvas_instance.lock_legend_pick = False
    
    else:
        # block legend picker
        canvas_instance.lock_legend_pick = True

    return None

def legend_picker_func(canvas_instance, event):
    """ Function to handle legend picker. """

    if canvas_instance.lock_legend_pick == False:
        if canvas_instance.plot_elements:
            # lock legend pick
            canvas_instance.lock_legend_pick = True
    
            # get legend label information
            legend_label = event.artist
            data_label = legend_label.get_text()

            if data_label not in canvas_instance.plot_elements['data_labels_active']:
                visible = True
                # put observations label always first in pop-ups on hover
                if data_label == canvas_instance.read_instance.observations_data_label:
                    canvas_instance.plot_elements['data_labels_active'].insert(0, data_label)
                # put experiment labels in the same order as in the legend
                else:
                    canvas_instance.plot_elements['data_labels_active'].insert(list(canvas_instance.read_instance.experiments.values()).index(data_label)+1, 
                                                                                data_label)
            else:
                visible = False
                canvas_instance.plot_elements['data_labels_active'].remove(data_label)

            # iterate through plot types stored in plot_elements (if have selected stations)
            if len(canvas_instance.relative_selected_station_inds) > 0:
                for plot_type in canvas_instance.plot_elements:  

                    if plot_type not in ['data_labels_active', 'metadata', 'map', 'heatmap', 
                                            'table', 'statsummary']:

                        # get currently selected options for plot
                        plot_options = canvas_instance.current_plot_options[plot_type]
                    
                        # get active (absolute / bias)
                        active = canvas_instance.plot_elements[plot_type]['active']

                        # change visibility of plot elements (if data label in plot elements dictionary)
                        if data_label in canvas_instance.plot_elements[plot_type][active]:
                            for element_type in canvas_instance.plot_elements[plot_type][active][data_label]:
                                for plot_element in canvas_instance.plot_elements[plot_type][active][data_label][element_type]:
                                    if visible:
                                        plot_element.set_visible(True)
                                    else:
                                        plot_element.set_visible(False)

                        # reset axes limits (harmonising across subplots for periodic plots) 
                        if plot_type not in ['taylor', "fairmode-statsummary"]:
                            if plot_type == 'scatter':
                                harmonise_xy_lims_paradigm(canvas_instance, canvas_instance.read_instance, 
                                                        canvas_instance.plot_axes[plot_type], plot_type, 
                                                        canvas_instance.plot_characteristics[plot_type], 
                                                        plot_options, relim=True)
                            else:
                                harmonise_xy_lims_paradigm(canvas_instance, canvas_instance.read_instance, 
                                                        canvas_instance.plot_axes[plot_type], plot_type, 
                                                        canvas_instance.plot_characteristics[plot_type], 
                                                        plot_options, relim=True, autoscale=True)

            # change font weight of label
            legend_label._fontproperties = canvas_instance.legend.get_texts()[0]._fontproperties.copy()
            if visible:
                legend_label.set_fontweight('bold')
            else:
                legend_label.set_fontweight('regular')

            # draw changes
            canvas_instance.figure.canvas.draw_idle()

            # unlock legend pick 
            canvas_instance.lock_legend_pick = False

    return None


class HoverAnnotation(object):

    def __init__(self, canvas_instance, add_vline=False):
        
        self.canvas_instance = canvas_instance

        # set up formatting for canvas annotations
        self.canvas_instance.figure.canvas = set_formatting(self.canvas_instance.figure.canvas, 
                                                            self.canvas_instance.read_instance.formatting_dict['canvas_annotation'])

        # set up formatting for canvas annotation vline
        self.canvas_instance.canvas_annotation_vline = set_formatting(QWidget(self.canvas_instance), self.canvas_instance.read_instance.formatting_dict['canvas_annotation_vline'])
        self.canvas_instance.canvas_annotation_vline.hide()

        return None

    def hover_annotation(self, event):
        """ Function that annotates on hover in timeseries, scatter, distribution, taylor and FAIRMODE target plots.
        """

        #hide annotation and vline
        self.canvas_instance.figure.canvas.setToolTip('')
        QToolTip.hideText()
        self.canvas_instance.canvas_annotation_vline.hide()

        # identify which axis is currently being hovered over
        plot_type = None
        for test_plot_type in self.canvas_instance.plot_axes:

            if not plot_type:
                if test_plot_type in ['periodic','periodic-violin']:
                    if hasattr(self.canvas_instance.read_instance, 'relevant_temporal_resolutions'):
                        for resolution in self.canvas_instance.read_instance.relevant_temporal_resolutions:
                            if event.inaxes == self.canvas_instance.plot_axes[test_plot_type][resolution]:
                                plot_type = copy.deepcopy(test_plot_type)
                                break
                elif test_plot_type == 'fairmode-statsummary':
                    for i in range(len(self.canvas_instance.plot_axes[test_plot_type])):
                        if event.inaxes == self.canvas_instance.plot_axes[test_plot_type][i]:
                            plot_type = copy.deepcopy(test_plot_type)
                            break
                else:
                    if event.inaxes == self.canvas_instance.plot_axes[test_plot_type]:
                        plot_type = copy.deepcopy(test_plot_type)
                        break

        # if an axis is being hovered over then now check if a point is being hovered over
        if plot_type:

            # add active axis to self
            if plot_type in ['periodic','periodic-violin']:
                self.ax = self.canvas_instance.plot_axes[plot_type][resolution]
            elif plot_type == 'fairmode-statsummary':
                self.ax = self.canvas_instance.plot_axes[plot_type][i]
            else:
                self.ax = self.canvas_instance.plot_axes[plot_type]

            # activate hover over plot
            if plot_type == 'map':
                search_plot = 'stations_scatter'
            elif plot_type == 'periodic':
                search_plot = 'periodic_plots'
            elif plot_type == 'periodic-violin':
                search_plot = 'violin_plot'
            else:
                search_plot = '{}_plot'.format(plot_type.replace('-','_'))

            # get plot element name
            if plot_type == 'periodic':
                plot_element_name = 'plot_{}'.format(resolution)
            elif plot_type == 'periodic-violin':
                plot_element_name = 'Median_plot_{}'.format(resolution)
            else:
                plot_element_name = 'plot'

            if ((hasattr(self.canvas_instance.plot, search_plot)) and (plot_type in self.canvas_instance.plot_elements)):

                is_contained = False

                for data_label in self.canvas_instance.plot_elements['data_labels_active']:

                    # skip observations for bias plot
                    if ((plot_type in ['timeseries', 'distribution','periodic','periodic-violin']) 
                        and (self.canvas_instance.plot_elements[plot_type]['active'] == 'bias')
                        and (data_label == self.canvas_instance.read_instance.observations_data_label)):
                        continue

                    # do not annotate if plot is cleared
                    if data_label not in self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']].keys():
                        continue
                    
                    if plot_type == 'map':
                        is_contained, annotation_index = self.canvas_instance.plot.stations_scatter.contains(event)
                    else:
                        # do no annotate if hidedata is active
                        if len(self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label][plot_element_name]) == 0:
                            continue
                        line = self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label][plot_element_name]
                        for n_point, point in enumerate(line):
                            is_contained, annotation_index = point.contains(event)
                            if is_contained:
                                if plot_type == 'fairmode-target':
                                    annotation_index = {'ind': np.array([n_point], dtype=np.int32)}
                                break

                    if is_contained:
                        break
                
                # point is being hovered over?
                if is_contained:
                    # add event coordinates to self (handling pixel scaling)
                    self.x = round(event.x / self.canvas_instance.read_instance.devicePixelRatio())
                    self.y = round(event.y / self.canvas_instance.read_instance.devicePixelRatio())

                    # update annotation and show vline
                    func = getattr(self, 'update_{}_annotation'.format(plot_type.replace('-','_')))
                    if plot_type in ['periodic','periodic-violin']:
                        func(annotation_index,resolution)
                    elif plot_type in ['fairmode-target','scatter','taylor']:
                        func(annotation_index,data_label)
                    else:
                        func(annotation_index)

        return None

    def update_map_annotation(self, annotation_index):
        """ Update annotation for each station that is hovered. """

        # retrieve stations references and coordinates
        station_name = self.canvas_instance.read_instance.station_names[self.canvas_instance.read_instance.networkspeci][self.canvas_instance.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_reference = self.canvas_instance.read_instance.station_references[self.canvas_instance.read_instance.networkspeci][self.canvas_instance.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_location = self.canvas_instance.plot.stations_scatter.get_offsets()[annotation_index['ind'][0]]
        station_value = self.canvas_instance.z_statistic[annotation_index['ind'][0]]

        # create annotation text
        text_label = ('Station: {0}\n').format(station_name)
        text_label += ('Reference: {0}\n').format(station_reference)
        text_label += ('Longitude: {0:.2f}\n').format(station_location[0])
        text_label += ('Latitude: {0:.2f}\n').format(station_location[1])
        text_label += ('{0}: {1:.{2}f}').format(self.canvas_instance.map_z_stat.currentText(), station_value, self.canvas_instance.plot_characteristics['map']['marker_annotate_rounding'])

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_timeseries_annotation(self, annotation_index):
        """ Update annotation for each timeseries point that is hovered. """

        # initialise annotation text
        text_label = ''

        #iterate through active data labels
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['timeseries']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue

            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['timeseries'][self.canvas_instance.plot_elements['timeseries']['active']].keys():
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements['timeseries'][self.canvas_instance.plot_elements['timeseries']['active']][data_label]['plot'][0]
            time = line.get_xdata()[annotation_index['ind'][0]]
            concentration = line.get_ydata()[annotation_index['ind'][0]]

            # first valid data label?
            if not text_label:

                # update vline position
                self.update_vline_position()

                # add time to annotation text
                text_label += ("<p style='white-space:pre'><i>Time: {0}</i>").format(time.astype('datetime64[us]').astype(datetime.datetime).strftime("%m/%d/%Y %H:%M:%S"))

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, data_label, concentration, self.canvas_instance.plot_characteristics['timeseries']['marker_annotate_rounding'])

        # end formatting of text label
        text_label += '</p>'

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_scatter_annotation(self, annotation_index, data_label):

        # initialise annotation text
        text_label = ''

        # do not annotate if plot is cleared
        if data_label not in self.canvas_instance.plot_elements['scatter'][self.canvas_instance.plot_elements['scatter']['active']].keys():
            return None

        # retrieve concentrations in x and y axis
        line = self.canvas_instance.plot_elements['scatter'][self.canvas_instance.plot_elements['scatter']['active']][data_label]['plot'][0]
        x = line.get_xdata()[annotation_index['ind'][0]]
        y = line.get_ydata()[annotation_index['ind'][0]]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(hex_colour, data_label)
        # observations label
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'x', x, self.canvas_instance.plot_characteristics['scatter']['marker_annotate_rounding'])
        # experiment label
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'y', y, self.canvas_instance.plot_characteristics['scatter']['marker_annotate_rounding'])

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_fairmode_target_annotation(self, annotation_index, data_label):

        # initialise annotation text
        text_label = ''

        # do not annotate if plot is cleared
        if data_label not in self.canvas_instance.plot_elements['fairmode-target'][self.canvas_instance.plot_elements['fairmode-target']['active']].keys():
            return None

        # retrieve CRMSE / β·RMSᵤ and Mean Bias / β·RMSᵤ
        line = self.canvas_instance.plot_elements['fairmode-target'][self.canvas_instance.plot_elements['fairmode-target']['active']][data_label]['plot'][annotation_index['ind'][0]]
        x = line.get_xdata()[0]
        y = line.get_ydata()[0]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(hex_colour, data_label)
        # CRMSE
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'CRMSE / β·RMSᵤ', x, self.canvas_instance.plot_characteristics['fairmode-target']['marker_annotate_rounding'])
        # MB
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'MB / β·RMSᵤ', y, self.canvas_instance.plot_characteristics['fairmode-target']['marker_annotate_rounding'])

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None
    
    def update_fairmode_statsummary_annotation(self, annotation_index, data_label):

        # initialise annotation text
        text_label = ''

        # do not annotate if plot is cleared
        if data_label not in self.canvas_instance.plot_elements['fairmode-target'][self.canvas_instance.plot_elements['fairmode-target']['active']].keys():
            return None

        # retrieve CRMSE / β·RMSᵤ and Mean Bias / β·RMSᵤ
        line = self.canvas_instance.plot_elements['fairmode-target'][self.canvas_instance.plot_elements['fairmode-target']['active']][data_label]['plot'][annotation_index['ind'][0]]
        x = line.get_xdata()[0]
        y = line.get_ydata()[0]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(hex_colour, data_label)
        # CRMSE
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'CRMSE / β·RMSᵤ', x, self.canvas_instance.plot_characteristics['fairmode-target']['marker_annotate_rounding'])
        # MB
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'MB / β·RMSᵤ', y, self.canvas_instance.plot_characteristics['fairmode-target']['marker_annotate_rounding'])

        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None
    
    def update_distribution_annotation(self, annotation_index):
        """ Update annotation for each distribution point that is hovered. """
        
        # initialise annotation text
        text_label = ''

        #iterate through active data labels
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # skip observations for bias plot
            if (self.canvas_instance.plot_elements['distribution']['active'] == 'bias') and (data_label == self.canvas_instance.read_instance.observations_data_label):
                continue
            
            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']].keys():
                continue

            # retrieve concentration and density
            line = self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']][data_label]['plot'][0]
            concentration = line.get_xdata()[annotation_index['ind'][0]]
            density = line.get_ydata()[annotation_index['ind'][0]]

             # first valid data label?
            if not text_label:

                # update vline position
                self.update_vline_position()

                # create annotation text
                text_label += ("<p style='white-space:pre'><i>{0}: {1:.{2}f}</i>").format(self.canvas_instance.read_instance.species[0], concentration, self.canvas_instance.plot_characteristics['distribution']['marker_annotate_rounding'])

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, data_label, density, self.canvas_instance.plot_characteristics['distribution']['marker_annotate_rounding'])

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None

    def update_taylor_annotation(self, annotation_index, data_label):
        
        # initialise annotation text
        text_label = ''
                
        # do not annotate if plot is cleared
        if data_label not in self.canvas_instance.plot_elements['taylor'][self.canvas_instance.plot_elements['taylor']['active']].keys():
            return None

        # retrieve time and concentration
        line = self.canvas_instance.plot_elements['taylor'][self.canvas_instance.plot_elements['taylor']['active']][data_label]['plot'][0]
        corr_stat = line.get_xdata()[annotation_index['ind'][0]]
        stddev = line.get_ydata()[annotation_index['ind'][0]]

        # get colour for data label
        colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

        # convert data label colour to hex code
        hex_colour = get_hex_code(colour)

        # add text label
        text_label += ('<font color="{0}">{1}</font>').format(hex_colour, data_label)
        # corr stat
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, self.canvas_instance.plot_characteristics['taylor']['corr_stat'], 
                                                                                np.cos(corr_stat), self.canvas_instance.plot_characteristics['taylor']['marker_annotate_rounding'])
        # stddev
        text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, 'StdDev', stddev, self.canvas_instance.plot_characteristics['taylor']['marker_annotate_rounding'])
    
        # update tooltip
        self.canvas_instance.figure.canvas.setToolTip(text_label)

        return None

    def update_periodic_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic point that is hovered. """

        # initialise annotation text
        text_label = ''

        #iterate through active data labels
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
                
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['periodic']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue
            
            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['periodic'][self.canvas_instance.plot_elements['periodic']['active']].keys():
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements['periodic'][self.canvas_instance.plot_elements['periodic']['active']][data_label]['plot_' + resolution][0]
            time = line.get_xdata()[annotation_index['ind'][0]]
            concentration = line.get_ydata()[annotation_index['ind'][0]]

            # first valid data label?
            if not text_label:

                # update vline position
                self.update_vline_position()

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.canvas_instance.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.canvas_instance.periodic_xticks[resolution]]
                    if resolution == 'dayofweek':
                        time_text = time_options[time]
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        time_text = time_options[time-1]
                        resolution_text = 'Month'
                text_label += ("<p style='white-space:pre'><i>{0}: {1}</i>").format(resolution_text, time_text)
        
            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, data_label, concentration, self.canvas_instance.plot_characteristics['periodic']['marker_annotate_rounding'])

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None
    
    def update_periodic_violin_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic violin point that is hovered. """

        # initialise annotation text
        text_label = ''

        #iterate through active data labels
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # skip observations for bias plot
            if self.canvas_instance.plot_elements['periodic-violin']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue
            
            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['periodic-violin'][self.canvas_instance.plot_elements['periodic-violin']['active']].keys():
                continue

            # retrieve time and concentration
            line = self.canvas_instance.plot_elements['periodic-violin'][self.canvas_instance.plot_elements['periodic-violin']['active']][data_label]['Median_plot_' + resolution][0]
            time = line.get_xdata()[annotation_index['ind'][0]]
            concentration = line.get_ydata()[annotation_index['ind'][0]]

            # first valid data label?
            if not text_label:

                # update vline position
                self.update_vline_position()

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.canvas_instance.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.canvas_instance.periodic_xticks[resolution]]
                    if resolution == 'dayofweek':
                        time_text = time_options[time]
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        time_text = time_options[time-1]
                        resolution_text = 'Month'
                text_label += ("<p style='white-space:pre'><i>{0}: {1}</i>").format(resolution_text, time_text)

            # get colour for data label
            colour = self.canvas_instance.read_instance.plotting_params[data_label]['colour']

            # convert data label colour to hex code
            hex_colour = get_hex_code(colour)

            # add text label
            text_label += ('<br><font color="{0}">{1}: {2:.{3}f}</font>').format(hex_colour, data_label, concentration, self.canvas_instance.plot_characteristics['periodic-violin']['marker_annotate_rounding'])

        # update tooltip and show vline
        self.canvas_instance.figure.canvas.setToolTip(text_label)
        self.canvas_instance.canvas_annotation_vline.show()

        return None
    
    def update_vline_position(self):
        """Function to update vline position 
        """

        # get current canvas width / height
        canvas_width, canvas_height = self.canvas_instance.figure.canvas.get_width_height()

        # get axis ylim (in data coordinates)
        ylim = self.ax.get_ylim()

        # transform matplotlib ylim data coordinates to display coordinates, handling pixel scaling
        ymin_display_mpl = round(self.ax.transData.transform([0,ylim[0]])[1] / self.canvas_instance.read_instance.devicePixelRatio())
        ymax_display_mpl = round(self.ax.transData.transform([0,ylim[1]])[1] / self.canvas_instance.read_instance.devicePixelRatio())
        
        # transform matplotlib display coordinates to Qt display coordinates                
        ymin_display_qt = int(canvas_height - ymin_display_mpl)
        ymax_display_qt = int(canvas_height - ymax_display_mpl)

        # calculate length of vline to plot
        vline_length = ymin_display_qt - ymax_display_qt 

        # update vline position 
        self.canvas_instance.canvas_annotation_vline.setGeometry(self.x, ymax_display_qt, 1, vline_length)

        return None