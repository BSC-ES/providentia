""" Objects and functions to interact with the axes """

from matplotlib.widgets import _SelectorWidget
from matplotlib.lines import Line2D
from PyQt5 import QtWidgets
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

    def __init__(self, canvas_instance, ax, onselect, useblit=True, props=None, button=None):
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
        self.canvas_instance = canvas_instance

    def _press(self, event):
        self.verts = [self._get_data(event)]
        self._selection_artist.set_visible(True)

    def _release(self, event):
        if self.verts is not None:
            self.verts.append(self._get_data(event))
            self.onselect(self.verts)
        self._selection_artist.set_data([[], []])
        self._selection_artist.set_visible(False)
        self.verts = None
        if self.canvas_instance.map_annotation_disconnect:
            self.canvas_instance.map_annotation_event = self.canvas_instance.figure.canvas.mpl_connect('motion_notify_event', self.canvas_instance.hover_map_annotation)
            self.canvas_instance.map_annotation_disconnect = False

    def _onmove(self, event):
        if not self.canvas_instance.map_annotation_disconnect:
            self.canvas_instance.figure.canvas.mpl_disconnect(self.canvas_instance.map_annotation_event)
            self.canvas_instance.map_annotation_disconnect = True
        if self.verts is None:
            return
        self.verts.append(self._get_data(event))
        self._selection_artist.set_data(list(zip(*self.verts)))

        self.update()


def zoom_map_func(canvas_instance, event):
    """ Function to handle zoom on map using scroll wheel. """

    if event.inaxes == canvas_instance.plot_axes['map']:
        
        if canvas_instance.lock_zoom == False:

            # lock zoom
            canvas_instance.lock_zoom = True

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

                # unlock zoom
                canvas_instance.lock_zoom = False

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

            # transform legend label into data labels
            for exp_label, exp_alias in canvas_instance.read_instance.experiments.items():
                if data_label == exp_alias:
                    data_label = exp_label
                    continue
            if data_label == canvas_instance.plot_characteristics['legend']['handles']['obs_label']:
                data_label = 'observations'

            if data_label not in canvas_instance.plot_elements['data_labels_active']:
                visible = True
                # put observations label always first in pop-ups on hover
                if data_label == 'observations':
                    canvas_instance.plot_elements['data_labels_active'].insert(0, data_label)
                # put experiment labels in the same order as in the legend
                else:
                    canvas_instance.plot_elements['data_labels_active'].insert(list(canvas_instance.read_instance.experiments.keys()).index(data_label)+1, 
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
                        if plot_type == 'scatter':
                            harmonise_xy_lims_paradigm(canvas_instance, canvas_instance.read_instance, 
                                                       canvas_instance.plot_axes[plot_type], plot_type, 
                                                       canvas_instance.plot_characteristics[plot_type], 
                                                       plot_options, relim=True)
                        elif plot_type != 'taylor':
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
