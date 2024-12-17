""" Functions to use plot options (annotate, log, etc.) """

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
import numpy as np

from .read_aux import drop_nans
from .statistics import calculate_statistic, get_z_statistic_info, exceedance_lim
from .warnings_prv import show_message
from .plot_aux import create_chunked_timeseries


def log_axes(relevant_axis, log_ax, plot_characteristics, undo=False):
    """ Log plot axes.

        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param log_ax: Axis to log
        :type log_ax: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param undo: Indicates if scale needs to be set to linear
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
                      plot_characteristics, plot_options):
    """ Add linear regression to plot.

        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
        :type networkspeci: str
        :param data_labels: Data arrays to plot  
        :type data_labels: list
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: Options to configure plots
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
    if 'max_points' in plot_characteristics:
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
           plot_characteristics, plot_options, chunk_stat=None, chunk_resolution=None):
    """ Add smooth line to plot.

        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
        :type networkspeci: str
        :param data_labels: Data arrays to plot   
        :type data_labels: list
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: Options to configure plots
        :type plot_options: list
        :param chunk_stat: Chunk statistic
        :type chunk_stat: str
        :param chunk_resolution: Chunk resolution
        :type chunk_resolution: str
    """

    # get valid data labels for networkspeci
    valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

    # cut data_labels for those in valid data labels
    cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

    # bias plot?
    if 'bias' in plot_options:
        bias = True
    else:
        bias = False

    # get chunking stat and resolution in dashboard
    if (not read_instance.offline) and (not read_instance.interactive):
        chunk_stat = canvas_instance.timeseries_chunk_stat.currentText()
        chunk_resolution = canvas_instance.timeseries_chunk_resolution.currentText()
        chunk_stat = None if chunk_stat == 'None' else chunk_stat
        chunk_resolution = None if chunk_resolution == 'None' else chunk_resolution
    
    # chunk timeseries
    if (chunk_stat is not None) and (chunk_resolution is not None):
        timeseries_data = create_chunked_timeseries(read_instance, canvas_instance, chunk_stat, 
                                                    chunk_resolution, networkspeci, cut_data_labels, 
                                                    bias)
    # normal timeseries
    else:
        timeseries_data = canvas_instance.selected_station_data[networkspeci]["timeseries"]
            
    # iterate through plotted data arrays making smooth line
    for data_label in cut_data_labels:

        # bias plot?
        if 'bias' in plot_options:
            
            # skip to next data label if making bias, and data label == observations
            if data_label == read_instance.observations_data_label:
                continue
            
            # chunk bias timeseries
            if (chunk_stat is not None) and (chunk_resolution is not None):
                ts = timeseries_data[data_label]
            # normal bias timeseries
            else:
                ts_obs = canvas_instance.selected_station_data[networkspeci]['timeseries'][read_instance.observations_data_label]
                ts_model = canvas_instance.selected_station_data[networkspeci]['timeseries'][data_label] 
                ts = ts_model - ts_obs

        # normal plot?
        else:
            ts = timeseries_data[data_label]

        # make smooth line
        min_points_percentage = plot_characteristics['smooth']['min_points_percentage'] / 100
        min_periods = int(round(plot_characteristics['smooth']['window'] * (min_points_percentage / 2)))
        smooth_line_data = ts.rolling(plot_characteristics['smooth']['window'], 
                                      min_periods=min_periods, 
                                      center=True, closed="both").mean()
        
        smooth_line = relevant_axis.plot(smooth_line_data,
                                         color=read_instance.plotting_params[data_label]['colour'],
                                         zorder=read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                         **plot_characteristics['smooth']['format'])

        # track plot elements if using dashboard 
        if (not read_instance.offline) and (not read_instance.interactive):
            canvas_instance.plot.track_plot_elements(data_label, base_plot_type, 'smooth', smooth_line, bias=bias)


def threshold(canvas_instance, read_instance, relevant_axis, networkspeci, base_plot_type, 
              plot_characteristics):
    """ Add threshold line/s to plot.

        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
        :type networkspeci: str
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
    """
    # get exceendance value
    threshold = exceedance_lim(networkspeci)

    # draw vertical line
    if base_plot_type in ['timeseries', 'scatter', 'periodic', 'periodic-violin', 'boxplot']:
        threshold_line = relevant_axis.axhline(y=threshold, 
                                               **plot_characteristics['threshold_line'])
    
    # draw horizontal line
    if base_plot_type in ['distribution', 'scatter']:
        threshold_line = relevant_axis.axvline(x=threshold, 
                                               **plot_characteristics['threshold_line'])

    # track plot elements if using dashboard 
    if (not read_instance.offline) and (not read_instance.interactive):
        canvas_instance.plot.track_plot_elements('ALL', base_plot_type, 'threshold', 
                                                 [threshold_line], bias=False)


def annotation(canvas_instance, read_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
               plot_characteristics, plot_options, plot_z_statistic_sign='absolute'):
    """ Add statistical annotations to plot.

        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
        :type read_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
        :type networkspeci: str
        :param data_labels: Data arrays to plot 
        :type data_labels: list
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: Options to configure plots
        :type plot_options: list
        :param plot_z_statistic_sign: sign of plotted z statistic (absolute or bias)
        :type plot_z_statistic_sign: str
    """

    # initialise list of strs to annotate
    str_to_annotate = []

    # add annotation text
    if base_plot_type == 'fairmode-target':
        bias = False
        if hasattr(canvas_instance.plot, 'faimode_target_annotate_text'):
            str_to_annotate = canvas_instance.plot.faimode_target_annotate_text
            colours = canvas_instance.plot.faimode_target_annotate_colour
    else:
        # get stats wished to be annotated
        stats = plot_characteristics['annotate_stats']
        
        # if no stats defined, then return
        if len(stats) == 0:
            msg = 'No annotation statistics are defined for {} in plot_characteristics.yaml.'.format(base_plot_type)
            show_message(read_instance, msg=msg)
            return

        # initialise colours of annotations
        colours = [] 

        # get valid data labels for networkspeci
        valid_data_labels = canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # bias plot? Then do not plot obs annotation label
        if 'bias' in plot_options:
            bias = True
            if read_instance.observations_data_label in cut_data_labels:
                cut_data_labels.remove(read_instance.observations_data_label)
        else:
            bias = False

        # making plot for a bias stat? Then do not plot obs annotation label
        if plot_z_statistic_sign == 'bias':
            if read_instance.observations_data_label in cut_data_labels:
                cut_data_labels.remove(read_instance.observations_data_label)

        # avoid plotting stats for observations data for scatter plots
        if base_plot_type == 'scatter':
            if read_instance.observations_data_label in cut_data_labels:
                cut_data_labels.remove(read_instance.observations_data_label)

        # generate annotation str to plot

        # show number of stations if defined
        if plot_characteristics['annotate_text']['n_stations']:
            colours.append('black')
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
                if (bias) or (plot_z_statistic_sign == 'bias') or (z_statistic_sign == 'bias'):
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

            if plot_characteristics['annotate_text']['color'] != "":
                colours = [plot_characteristics['annotate_text']['color']]*len(cut_data_labels)
    
    if len(str_to_annotate) != 0:
        # add annotation to plot
        lines = [TextArea(line, textprops=dict(color=colour, 
                                               size=plot_characteristics['annotate_text']['fontsize'])) 
                 for line, colour in zip(str_to_annotate, colours)]
        bbox = AnchoredOffsetbox(child=VPacker(children=lines, align='left', pad=0, 
                                               sep=plot_characteristics['annotate_text']['sep']),
                                 bbox_transform=relevant_axis.transAxes, 
                                 **plot_characteristics['annotate_offset'])
        # set zorder for plots that have the annotation box on top of the plot
        if base_plot_type != 'fairmode-target':
            bbox.zorder = plot_characteristics['annotate_bbox']['zorder']
        bbox.patch.set(**plot_characteristics['annotate_bbox'])
        relevant_axis.add_artist(bbox)

        # track plot elements if using dashboard 
        if (not read_instance.offline) and (not read_instance.interactive):
            canvas_instance.plot.track_plot_elements('ALL', base_plot_type, 'annotate', [bbox], bias=bias)
    else:
        msg = '{} could not be annotated'.format(base_plot_type)
        show_message(read_instance, msg)


def experiment_domain(canvas_instance, relevant_axis, data_labels, map_extent):
    """ Plot experiment domain extents on map

        :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
        :type canvas_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param data_labels: Data arrays to plot 
        :type data_labels: list
        :param map_extent: list of map extent bounds [lonmin, lonmax, latmin, latmax]
        :type map_extent: list
    """

    #get experiment domain polygons
    grid_edge_polygons = canvas_instance.plot.make_experiment_domain_polygons(data_labels=data_labels) 

    # plot grid edge polygons on map
    for grid_edge_polygon in grid_edge_polygons:
        relevant_axis.add_patch(grid_edge_polygon)

    # if map extent is not set then re-set automatic limits based now domain is plotted.
    if not map_extent:
        relevant_axis.relim(visible_only=True)
        relevant_axis.autoscale(tight=False)