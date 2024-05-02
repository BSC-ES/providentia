""" Objects and functions to interact with the axes """

import copy
import datetime
import math
import time

import matplotlib
from matplotlib.lines import Line2D
from matplotlib.widgets import _SelectorWidget
import numpy as np
from packaging.version import Version
import pandas as pd

from .plot_aux import get_map_extent
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

                        #print(plot_type, active, visible)

                        # change visibility of plot elements (if data label in plot elements dictionary)
                        if data_label in canvas_instance.plot_elements[plot_type][active]:
                            for element_type in canvas_instance.plot_elements[plot_type][active][data_label]:
                                for plot_element in canvas_instance.plot_elements[plot_type][active][data_label][element_type]:
                                    #if plot_type == 'periodic':
                                        #print(plot_element)

                                    if visible:
                                        plot_element.set_visible(True)
                                    else:
                                        plot_element.set_visible(False)

                        # reset axes limits (harmonising across subplots for periodic plots) 
                        if plot_type not in ['taylor']:
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

    def __init__(self, canvas_instance, plot_type, ax, plot_characteristics, add_vline=False, resolution=None):
        
        self.canvas_instance = canvas_instance
        
        # get reference coordinates
        if plot_type == 'map':
            xycoords = canvas_instance.datacrs._as_mpl_transform(ax)
        else:
            xycoords = 'data'

        # create annotation
        # using matplotlib < 3.3.0, text corresponds to s
        annotation_dict = {"xy": [0, 0], 
                           "xycoords": xycoords, 
                           "bbox": {**plot_characteristics['marker_annotate_bbox']},
                           "arrowprops": {**plot_characteristics['marker_annotate_arrowprops']}}
        if Version(matplotlib.__version__) >= Version("3.3"):
            annotation_dict["text"] = ""
        else:
            annotation_dict["s"] = ""
        self.annotation = ax.annotate(**{**annotation_dict, **plot_characteristics['marker_annotate']})
        
        # hide when it is initialised
        self.annotation.set_visible(False)

        # add to elements for legend picking
        canvas_instance.annotation_elements.extend([self.annotation])

        # add vertical line
        if add_vline:
            self.vline = ax.axvline(0, 
                **plot_characteristics['marker_annotate_vline'])
            self.vline.set_visible(False)
            canvas_instance.annotation_elements.extend([self.vline])

        # initialise dict with values that are in the middle of x axis per plot type
        self.x_middle = {}
        self.x_middle[plot_type] = {}

        return None

    def hover_annotation(self, event, plot_type):
        """ Function that annotates on hover in timeseries, scatter, distribution and taylor plots.
        """

        # activate hover over plot
        if (plot_type in self.canvas_instance.read_instance.active_dashboard_plots):
            if event.inaxes == self.canvas_instance.plot_axes[plot_type]:

                if ((hasattr(self.canvas_instance.plot, plot_type + '_plot')) 
                    and (plot_type in self.canvas_instance.plot_elements)
                    and (self.canvas_instance.annotations_lock[plot_type] == False)):

                    # lock annotation
                    self.canvas_instance.annotations_lock[plot_type] = True
                    is_contained = False

                    for data_label in self.canvas_instance.plot_elements['data_labels_active']:

                        # skip observations for bias plot
                        if ((plot_type in ['timeseries', 'distribution']) 
                            and (self.canvas_instance.plot_elements[plot_type]['active'] == 'bias')
                            and (data_label == self.canvas_instance.read_instance.observations_data_label)):
                            continue

                        # do not annotate if plot is cleared
                        if data_label not in self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']].keys():
                            continue
                        
                        # do no annotate if hidedata is active
                        if len(self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label]['plot']) == 0:
                            continue
                        
                        line = self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label]['plot'][0]
                        is_contained, annotation_index = line.contains(event)
                        if is_contained:
                            self.annotate_data_label = data_label
                            break
                    
                    if is_contained:
                        # update annotation if hovered
                        func = getattr(self, 'update_' + plot_type + '_annotation')
                        func(annotation_index)
                        self.annotation.set_visible(True)
                        if hasattr(self, 'vline'):
                            self.vline.set_visible(True)
                    else:
                        # hide annotation if not hovered
                        if self.annotation.get_visible():
                            self.annotation.set_visible(False)
                            if hasattr(self, 'vline'):
                                self.vline.set_visible(False)

                    # draw changes
                    self.canvas_instance.figure.canvas.draw_idle()
                        
                    # unlock annotation 
                    self.canvas_instance.annotations_lock[plot_type] = False

        return None
    
    def update_x_middle(self, event, plot_type):
        """ Function to find middle value in x axis per plot type.
        """

        # do not annotate if plot has not been made yet
        if plot_type not in self.canvas_instance.plot_elements:
            return

        # get current limits on x axis
        xlim_range = event.get_xlim()

        # transform range into dates for timeseries
        if plot_type == 'timeseries':
            xdata_range = [matplotlib.dates.num2date(xlim) for xlim in xlim_range]
        else:
            xdata_range = xlim_range
            
        # get value/date in the middle of range
        x_middle = xdata_range[0] + (xdata_range[1] - xdata_range[0])/2

        # save into dictionary
        if 'periodic' in plot_type:
            for resolution in self.canvas_instance.read_instance.relevant_temporal_resolutions:
                if event == self.canvas_instance.plot_axes[plot_type][resolution]:
                    self.x_middle[plot_type][resolution] = x_middle
                    break
        else:
            self.x_middle[plot_type] = x_middle

    def update_timeseries_annotation(self, annotation_index):
        """ Update annotation for each timeseries point that is hovered. """
    
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
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

                # update location
                self.annotation.xy = (time, concentration)
                self.vline.set_xdata(time)

                # update bbox position
                if time > np.datetime64(self.x_middle['timeseries']):
                    self.annotation.set_x(-10)
                    self.annotation.set_ha('right')
                else:
                    self.annotation.set_x(10)
                    self.annotation.set_ha('left')

                # create annotation text
                text_label = ('Time: {0}').format(time.astype('datetime64[us]').astype(datetime.datetime).strftime("%m/%d/%Y %H:%M:%S"))
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['timeseries']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue

            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['timeseries'][self.canvas_instance.plot_elements['timeseries']['active']].keys():
                continue

            # retrieve concentration
            line = self.canvas_instance.plot_elements['timeseries'][self.canvas_instance.plot_elements['timeseries']['active']][data_label]['plot'][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]

            # for all labels if there is data
            if len(concentration) >= 1:
                text_label += ('\n{0}: {1:.2f}').format(data_label, concentration[0])

        self.annotation.set_text(text_label)

        return None

    def update_scatter_annotation(self, annotation_index):

        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
                # do not annotate if plot is cleared
                if data_label not in self.canvas_instance.plot_elements['scatter'][self.canvas_instance.plot_elements['scatter']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.canvas_instance.plot_elements['scatter'][self.canvas_instance.plot_elements['scatter']['active']][data_label]['plot'][0]
                concentration_x = line.get_xdata()[annotation_index['ind'][0]]
                concentration_y = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.annotation.xy = (concentration_x, concentration_y)

                # update bbox position
                if concentration_x > self.x_middle['scatter']:
                    self.annotation.set_x(-10)
                    self.annotation.set_ha('right')
                else:
                    self.annotation.set_x(10)
                    self.annotation.set_ha('left')

                # create annotation text
                # experiment label
                text_label = copy.deepcopy(data_label)
                # observations label
                text_label += ('\n{0}: {1:.2f}').format('x', concentration_x)
                # experiment label
                text_label += ('\n{0}: {1:.2f}').format('y', concentration_y)

        self.annotation.set_text(text_label)

        return None

    def update_distribution_annotation(self, annotation_index):
        """ Update annotation for each distribution point that is hovered. """
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
                # skip observations for bias plot
                if self.canvas_instance.plot_elements['distribution']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                    continue
                
                # do not annotate if plot is cleared
                if data_label not in self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']][data_label]['plot'][0]
                concentration = line.get_xdata()[annotation_index['ind'][0]]
                density = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.annotation.xy = (concentration, density)
                self.vline.set_xdata(concentration)

                # update bbox position
                if concentration > self.x_middle['distribution']:
                    self.annotation.set_x(-10)
                    self.annotation.set_ha('right')
                else:
                    self.annotation.set_x(10)
                    self.annotation.set_ha('left')

                # create annotation text
                text_label = ('{0}: {1:.3f}').format(self.canvas_instance.read_instance.species[0], concentration)
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['distribution']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue

            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']].keys():
                continue

            # retrieve density
            line = self.canvas_instance.plot_elements['distribution'][self.canvas_instance.plot_elements['distribution']['active']][data_label]['plot'][0]
            density = line.get_ydata()[np.where(line.get_xdata() == concentration)[0]]

            # for all labels if there is data
            if len(density) >= 1:
                text_label += ('\n{0}: {1:.3f}').format(data_label, density[0])
   
        self.annotation.set_text(text_label)

        return None

    def update_taylor_annotation(self, annotation_index):
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
                # do not annotate if plot is cleared
                if data_label not in self.canvas_instance.plot_elements['taylor'][self.canvas_instance.plot_elements['taylor']['active']].keys():
                    continue

                # retrieve time and concentration
                line = self.canvas_instance.plot_elements['taylor'][self.canvas_instance.plot_elements['taylor']['active']][data_label]['plot'][0]
                corr_stat = line.get_xdata()[annotation_index['ind'][0]]
                stddev = line.get_ydata()[annotation_index['ind'][0]]

                # update location
                self.annotation.xy = (corr_stat, stddev)

                # update bbox position
                corr_stat_middle = line.get_xdata()[math.floor((len(line.get_xdata()) - 1)/2)]
                if corr_stat > corr_stat_middle:
                    self.annotation.set_x(-10)
                    self.annotation.set_ha('right')
                else:
                    self.annotation.set_x(10)
                    self.annotation.set_ha('left')

                # create annotation text
                text_label = copy.deepcopy(data_label)
                text_label += ('\n{0}: {1:.2f}').format(self.canvas_instance.plot_characteristics['taylor']['corr_stat'], 
                                                        np.cos(corr_stat))
                text_label += ('\n{0}: {1:.2f}').format('StdDev', stddev)
        
        self.annotation.set_text(text_label)

        return None

    def hover_periodic_annotation(self, event, plot_type):
        """ Function that annotates on hover in periodic and periodic violin plots.
        """

        # activate hover over periodic plot
        if (plot_type in self.canvas_instance.read_instance.active_dashboard_plots):
            if hasattr(self.canvas_instance.read_instance, 'relevant_temporal_resolutions'):
                for resolution in self.canvas_instance.read_instance.relevant_temporal_resolutions:
                    if event.inaxes == self.canvas_instance.plot_axes[plot_type][resolution]:
                        search_plot = 'periodic_plots' if plot_type == 'periodic' else 'violin_plot'
                        if ((hasattr(self.canvas_instance.plot, search_plot)) and (plot_type in self.canvas_instance.plot_elements)
                            and (self.canvas_instance.annotations_lock[plot_type][resolution] == False)):

                            # lock annotation
                            self.canvas_instance.annotations_lock[plot_type][resolution] = True
                            is_contained = False

                            for data_label in self.canvas_instance.plot_elements['data_labels_active']:

                                # skip observations for bias plot
                                if self.canvas_instance.plot_elements[plot_type]['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                                    continue

                                # do not annotate if plot is cleared
                                if data_label not in self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']].keys():
                                    continue
                                
                                if plot_type == 'periodic':
                                    line = self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label]['plot_' + resolution][0]
                                else:
                                    line = self.canvas_instance.plot_elements[plot_type][self.canvas_instance.plot_elements[plot_type]['active']][data_label]['Median_plot_' + resolution][0]
                                is_contained, annotation_index = line.contains(event)
                                if is_contained:
                                    self.annotate_data_label = data_label
                                    break
                            
                            if is_contained:
                                # update annotation if hovered
                                if plot_type == 'periodic':
                                    self.update_periodic_annotation(annotation_index, resolution)
                                else:
                                    self.update_periodic_violin_annotation(annotation_index, resolution)
                                self.canvas_instance.annotations[plot_type][resolution].set_visible(True)
                                self.canvas_instance.annotations_vline[plot_type][resolution].set_visible(True)
                            else:
                                # hide annotation if not hovered
                                if self.canvas_instance.annotations[plot_type][resolution].get_visible():
                                    self.canvas_instance.annotations[plot_type][resolution].set_visible(False)
                                    self.canvas_instance.annotations_vline[plot_type][resolution].set_visible(False)
                                    
                            # draw changes
                            self.canvas_instance.figure.canvas.draw_idle()
                                
                            # unlock annotation 
                            self.canvas_instance.annotations_lock[plot_type][resolution] = False

        return None

    def update_periodic_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic point that is hovered. """

        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
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

                # update location
                self.canvas_instance.annotations['periodic'][resolution].xy = (time, concentration)
                self.canvas_instance.annotations_vline['periodic'][resolution].set_xdata(time)

                # update bbox position
                if time > self.x_middle['periodic'][resolution]:
                    self.canvas_instance.annotations['periodic'][resolution].set_x(-10)
                    self.canvas_instance.annotations['periodic'][resolution].set_ha('right')
                else:
                    self.canvas_instance.annotations['periodic'][resolution].set_x(10)
                    self.canvas_instance.annotations['periodic'][resolution].set_ha('left')

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.canvas_instance.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.canvas_instance.periodic_xticks[resolution]]
                    time_text = time_options[time-1]
                    if resolution == 'dayofweek':
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        resolution_text = 'Month'
                text_label = ('{0}: {1}').format(resolution_text, time_text)
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['periodic']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue

            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['periodic'][self.canvas_instance.plot_elements['periodic']['active']].keys():
                continue

            # retrieve concentration
            line = self.canvas_instance.plot_elements['periodic'][self.canvas_instance.plot_elements['periodic']['active']][data_label]['plot_' + resolution][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]
            
            # for all labels if there is data
            if len(concentration) >= 1:
                text_label += ('\n{0}: {1:.2f}').format(data_label, concentration[0])

        self.canvas_instance.annotations['periodic'][resolution].set_text(text_label)

        return None
    
    def update_periodic_violin_annotation(self, annotation_index, resolution):
        """ Update annotation for each periodic violin point that is hovered. """

        for data_label in self.canvas_instance.plot_elements['data_labels_active']:

            # for annotate data label
            if data_label == self.annotate_data_label:
                
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

                # update location
                self.canvas_instance.annotations['periodic-violin'][resolution].xy = (time, concentration)
                self.canvas_instance.annotations_vline['periodic-violin'][resolution].set_xdata(time)

                # update bbox position
                if time > self.x_middle['periodic-violin'][resolution]:
                    self.canvas_instance.annotations['periodic-violin'][resolution].set_x(-10)
                    self.canvas_instance.annotations['periodic-violin'][resolution].set_ha('right')
                else:
                    self.canvas_instance.annotations['periodic-violin'][resolution].set_x(10)
                    self.canvas_instance.annotations['periodic-violin'][resolution].set_ha('left')

                # create annotation text
                if resolution == 'hour':
                    resolution_text = 'Hour'
                    time_text = time
                else:
                    time_options = [self.canvas_instance.temporal_axis_mapping_dict['long'][resolution][xtick] 
                                    for xtick in self.canvas_instance.periodic_xticks[resolution]]
                    time_text = time_options[time-1]
                    if resolution == 'dayofweek':
                        resolution_text = 'Day'
                    elif resolution == 'month':
                        resolution_text = 'Month'
                text_label = ('{0}: {1}').format(resolution_text, time_text)
        
        for data_label in self.canvas_instance.plot_elements['data_labels_active']:
            
            # skip observations for bias plot
            if self.canvas_instance.plot_elements['periodic-violin']['active'] == 'bias' and data_label == self.canvas_instance.read_instance.observations_data_label:
                continue

            # do not annotate if plot is cleared
            if data_label not in self.canvas_instance.plot_elements['periodic-violin'][self.canvas_instance.plot_elements['periodic-violin']['active']].keys():
                continue

            # retrieve concentration
            line = self.canvas_instance.plot_elements['periodic-violin'][self.canvas_instance.plot_elements['periodic-violin']['active']][data_label]['Median_plot_' + resolution][0]
            concentration = line.get_ydata()[np.where(line.get_xdata() == time)[0]]
            
            # for all labels if there is data
            if len(concentration) >= 1:
                text_label += ('\n{0}: {1:.2f}').format(data_label, concentration[0])

        self.canvas_instance.annotations['periodic-violin'][resolution].set_text(text_label)

        return None

    def hover_map_annotation(self, event):
        """ Function that annotates on hover in map.
        """

        if event.inaxes == self.canvas_instance.plot_axes['map']:

            # activate hover over map
            if (hasattr(self.canvas_instance.plot, 'stations_scatter')):

                is_contained, annotation_index = self.canvas_instance.plot.stations_scatter.contains(event)
                
                if is_contained:
                    # update annotation if hovered
                    self.update_map_annotation(annotation_index)
                    self.annotation.set_visible(True)
                else:
                    # hide annotation if not hovered
                    if self.annotation.get_visible():
                        self.annotation.set_visible(False)

                # draw changes
                self.canvas_instance.figure.canvas.draw_idle()

        return None

    def update_map_annotation(self, annotation_index):
        """ Update annotation for each station that is hovered. """

        # retrieve stations references and coordinates
        station_name = self.canvas_instance.read_instance.station_names[self.canvas_instance.read_instance.networkspeci][self.canvas_instance.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_reference = self.canvas_instance.read_instance.station_references[self.canvas_instance.read_instance.networkspeci][self.canvas_instance.active_map_valid_station_inds[annotation_index['ind'][0]]]
        station_location = self.canvas_instance.plot.stations_scatter.get_offsets()[annotation_index['ind'][0]]
        station_value = self.canvas_instance.z_statistic[annotation_index['ind'][0]]

        # update location
        self.annotation.xy = station_location

        # update bbox position
        self.canvas_instance.read_instance.map_extent = get_map_extent(self.canvas_instance)
        lat_min = self.canvas_instance.read_instance.map_extent[2]
        lat_max = self.canvas_instance.read_instance.map_extent[3]
        if station_location[1] > ((lat_max + lat_min) / 2):
            self.annotation.set_y(-10)
            self.annotation.set_va('top')
        else:
            self.annotation.set_y(10)
            self.annotation.set_va('bottom')

        # create annotation text
        text_label = ('Station: {0}\n').format(station_name)
        text_label += ('Reference: {0}\n').format(station_reference)
        text_label += ('Longitude: {0:.2f}\n').format(station_location[0])
        text_label += ('Latitude: {0:.2f}\n').format(station_location[1])
        text_label += ('{0}: {1:.2f}').format(self.canvas_instance.map_z_stat.currentText(), station_value)
        self.annotation.set_text(text_label)

        return None
       