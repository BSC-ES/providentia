""" Functions to use plot options (annotate, log, etc.) """

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
import numpy as np

from .read_aux import drop_nans
from .statistics import calculate_statistic, get_z_statistic_info
from .warnings import show_message


def log_axes(relevant_axis, log_ax, plot_characteristics, undo=False):
    """ Log plot axes.

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param log_ax: which axis to log
        :type log_ax: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param undo: unlog plot axes
        :type undo: boolean
    """

    if not undo:
        if log_ax == 'logx':
            relevant_axis.set_xscale('log')
            
        if log_ax == 'logy':
            relevant_axis.set_yscale('log')

    else:
        if log_ax == 'logx':
            relevant_axis.set_xscale('linear')
        
        if log_ax == 'logy':
            relevant_axis.set_yscale('linear')
        

def linear_regression(canvas_instance, read_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
                      plot_characteristics, plot_options=[]):
    """ Add linear regression to plot.

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param networkspeci: str of currently active network and species 
        :type networkspeci: str
        :param data_labels: names of plotted data arrays  
        :type data_labels: list
        :param base_plot_type: plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
    """

    # get valid data labels for networkspeci
    valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

    # cut data_labels for those in valid data labels
    cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

    # get observations data (flattened and drop NaNs)
    observations_data = drop_nans(canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(read_instance.observations_data_label),0,:])

    # determine if number of points per data array exceeds max limit,
    # if so subset arrays
    subset = False
    data_array_size = observations_data.size
    if data_array_size > plot_characteristics['max_points']:
        subset = True
        inds_subset = np.random.choice(data_array_size, size=plot_characteristics['max_points'], replace=False)
        observations_data = observations_data[inds_subset]

    # iterate through experiment data, making regression line to observations
    for data_label in cut_data_labels:
        if data_label != read_instance.observations_data_label:
            # get experiement data (flattened and drop NaNs)
            experiment_data = drop_nans(canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])
            # subset data if neccessary
            if subset:
                experiment_data = experiment_data[inds_subset]
            m, b = np.polyfit(observations_data, experiment_data, deg=1)
            regression_line = relevant_axis.plot(observations_data, m*observations_data+b, 
                                                 color=read_instance.plotting_params[data_label]['colour'],
                                                 zorder=read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                                 **plot_characteristics['regression'])
            
            # track plot elements if using dashboard 
            if (not read_instance.offline) and (not read_instance.interactive) :
                canvas_instance.plot.track_plot_elements(data_label, base_plot_type, 'regression', regression_line, bias=False)


def smooth(canvas_instance, read_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
           plot_characteristics, plot_options=[]):
    """ Add smooth line to plot.

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param networkspeci: str of currently active network and species 
        :type networkspeci: str
        :param data_labels: names of plotted data arrays   
        :type data_labels: list
        :param base_plot_type: plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
    """

    # get valid data labels for networkspeci
    valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

    # cut data_labels for those in valid data labels
    cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

    # iterate through plotted data arrays making smooth line
    for data_label in cut_data_labels:

        # bias plot?
        if 'bias' in plot_options:
            # skip to next data label if making bias, and data label == observations
            if data_label == read_instance.observations_data_label:
                continue
            ts_obs = canvas_instance.selected_station_data[networkspeci]['timeseries'][read_instance.observations_data_label]
            ts_model = canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label] 
            ts = ts_model - ts_obs
            bias = True
        # normal plot?
        else:
            ts = canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label]
            bias = False

        # make smooth line
        smooth_line_data = ts.rolling(plot_characteristics['smooth']['window'], 
                                      min_periods=plot_characteristics['smooth']['min_points'], 
                                      center=True).mean()
        smooth_line = relevant_axis.plot(smooth_line_data,
                                         color=read_instance.plotting_params[data_label]['colour'],
                                         zorder=read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                         **plot_characteristics['smooth']['format'])

        # track plot elements if using dashboard 
        if (not read_instance.offline) and (not read_instance.interactive):
            canvas_instance.plot.track_plot_elements(data_label, base_plot_type, 'smooth', smooth_line, bias=bias)


def annotation(canvas_instance, read_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
               plot_characteristics, plot_options=[], plotting_paradigm=None):
    """ Add statistical annotations to plot.

        :param relevant_axis: axis to plot on 
        :type relevant_axis: object
        :param networkspeci: str of currently active network and species 
        :type networkspeci: str
        :param data_labels: names of plotted data arrays 
        :type data_labels: list
        :param base_plot_type: plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: list of options to configure plots
        :type plot_options: list
        :param plotting_paradigm: plotting paradigm (summary or station in offline reports)
        :type plotting_paradigm: str
    """

    # get stats wished to be annotated
    stats = plot_characteristics['annotate_stats']
    
    # if no stats defined, then return
    if len(stats) == 0:
        msg_dashboard = 'No annotation statistics are defined for {} in plot_characteristics_dashboard.json.'.format(base_plot_type)
        msg_offline = 'No annotation statistics are defined for {} in plot_characteristics_offline.json.'.format(base_plot_type)
        show_message(read_instance, msg=msg_dashboard, msg_offline=msg_offline)
        return

    # initialise list of strs to annotate, and colours of annotations
    str_to_annotate = []
    colours = []

    # get valid data labels for networkspeci
    valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

    # cut data_labels for those in valid data labels
    cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

    # bias plot?  
    if 'bias' in plot_options:
        bias = True
        if read_instance.observations_data_label in cut_data_labels:
            cut_data_labels.remove(read_instance.observations_data_label)
    else:
        bias = False

    # avoid plotting stats for observations data for scatter plots
    if base_plot_type == 'scatter':
        if read_instance.observations_data_label in cut_data_labels:
            cut_data_labels.remove(read_instance.observations_data_label)

    # generate annotation str to plot

    # show number of stations if defined
    if plot_characteristics['annotate_text']['n_stations']:
        colours.append('black')
        if (read_instance.offline) or (read_instance.interactive):
            if plotting_paradigm == 'station':
                str_to_annotate.append('Stations: 1')
            else:
                str_to_annotate.append('Stations: ' + str(len(canvas_instance.station_inds[networkspeci])))
        else:
            str_to_annotate.append('Stations: ' + str(len(canvas_instance.station_inds[networkspeci])))

    # generate annotation line by line (one line per data label, for all stats)
    for data_label_ii, data_label in enumerate(cut_data_labels):

        # get colour for data label
        colours.append(read_instance.plotting_params[data_label]['colour'])

        # iterate through stats to calculate
        stats_annotate = []
        for zstat in stats:

            # get zstat information
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

            # calculate stats
            if (bias) or (z_statistic_sign == 'bias'):
                if data_label != read_instance.observations_data_label:
                    stat_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, zstat, 
                                                    [read_instance.observations_data_label], [data_label])
                # skip bias stats for observations
                else:
                    continue
            else:
                stat_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, zstat, 
                                                [data_label], [])

            # format annotation line
            stats_annotate.append("{0}: {1:.{2}f}".format(zstat, stat_calc[0],
                                                          plot_characteristics['annotate_text']['round_decimal_places'])) 

        # append annotation line
        if (plot_characteristics['annotate_text']['exp_labels']):
            str_to_append = data_label + ' | ' + ', '.join(stats_annotate)
        else:
            str_to_append = ', '.join(stats_annotate)
        str_to_annotate.append(str_to_append)

    if len(str_to_annotate) != 0:
        if plot_characteristics['annotate_text']['color'] != "":
            colours = [plot_characteristics['annotate_text']['color']]*len(cut_data_labels)

        # add annotation to plot
        # see loc options at https://matplotlib.org/3.1.0/api/offsetbox_api.html
        lines = [TextArea(line, textprops=dict(color=colour, 
                                               size=plot_characteristics['annotate_text']['fontsize'])) 
                 for line, colour in zip(str_to_annotate, colours)]
        bbox = AnchoredOffsetbox(child=VPacker(children=lines, align='left', pad=0, sep=1),
                                 loc=plot_characteristics['annotate_text']['loc'],
                                 bbox_transform=relevant_axis.transAxes)
        bbox.zorder = plot_characteristics['annotate_bbox']['zorder']
        bbox.patch.set(**plot_characteristics['annotate_bbox'])
        relevant_axis.add_artist(bbox)

        # track plot elements if using dashboard 
        if (not read_instance.offline) and (not read_instance.interactive):
            canvas_instance.plot.track_plot_elements('ALL', base_plot_type, 'annotate', [bbox], bias=bias)
    else:
        msg = '{} could not be annotated'.format(base_plot_type)
        show_message(read_instance, msg)


def get_no_margin_lim(ax, lim):
    """ Get true limits of plot area. """

    # xlim
    if lim == 'xlim':
        xlim = ax.get_xlim()
        xwidth = xlim[1] - xlim[0]
        lower_lim = xlim[0] + (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth
        upper_lim = xlim[1] - (0.5 * ax.margins()[0]) / (0.5 + ax.margins()[0]) * xwidth

    # ylim
    if lim == 'ylim':
        ylim = ax.get_ylim()
        ywidth = ylim[1] - ylim[0]
        lower_lim = ylim[0] + (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth
        upper_lim = ylim[1] - (0.5 * ax.margins()[1]) / (0.5 + ax.margins()[1]) * ywidth

    return lower_lim, upper_lim


def log_validity(relevant_axis, log_ax):
    """ Determine if log operation for a given axes is valid (no values <= 0).
    
        :param relevant_axis: relevant axes
        :type relevant_axis: list
        :param log_ax: which axis to log
        :type log_ax: str
        :return: validity to log axis
        :rtype: boolean
    """

    if log_ax == 'logx':
        lower_lim, _ = get_no_margin_lim(relevant_axis, 'xlim')
        if round(lower_lim, 2) >= 0:
            validity = True
        else:
            validity = False
    
    if log_ax == 'logy':
        lower_lim, _ = get_no_margin_lim(relevant_axis, 'ylim')
        if round(lower_lim, 2) >= 0:
            validity = True
        else:
            validity = False

    return validity
