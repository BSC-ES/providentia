""" Functions for applying plot options (annotate, log, etc.) """

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
import numpy as np

from .read_aux import drop_nans
from .statistics import calculate_statistic, get_z_statistic_info, exceedance_lim
from .warnings_prv import show_message
from .plot_aux import create_statistical_timeseries

def log_axes(relevant_axis, log_ax, plot_characteristics, undo=False):
    """
    Log or un-log a plot axis.

    Parameters
    ----------
    relevant_axis : object
        Axis to apply the log scale on.
    log_ax : str
        Which axis to log. Options are 'logx' or 'logy'.
    plot_characteristics : dict
        Plot characteristics.
    undo : bool, optional
        If True, sets the axis scale back to linear. Default is False.
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

def linear_regression(read_instance, canvas_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
                      plot_characteristics, plot_options):
    """
    Add linear regression to plot.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axis : object
        Axis to plot on.
    networkspeci : str
        Current networkspeci (e.g. EBAS|sconco3).
    data_labels : list
        Data arrays to plot.
    base_plot_type : str
        Plot type, without statistical information.
    plot_characteristics : dict
        Plot characteristics.
    plot_options : list
        Options to configure plots.
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

    # iterate through model data, making regression line to observations
    for data_label in cut_data_labels:
        if data_label != read_instance.observations_data_label:
            # get model data (flattened and drop NaNs)
            model_data = drop_nans(canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])
            # subset data if neccessary
            if subset:
                model_data = model_data[inds_subset]
            m, b = np.polyfit(observations_data, model_data, deg=1)
            regression_line = relevant_axis.plot(observations_data, m*observations_data+b, 
                                                 color=read_instance.plotting_params[data_label]['colour'],
                                                 zorder=read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels),
                                                 **plot_characteristics['regression'])
            
            # track plot elements if using dashboard 
            if read_instance.mode not in ['report', 'library']:
                canvas_instance.plotting.track_plot_elements(data_label, base_plot_type, 'regression', regression_line, bias=False)

def smooth(read_instance, canvas_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
           plot_characteristics, plot_options, chunk_stat=None, chunk_resolution=None):
    """
    Add smooth line to plot.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axis : object
        Axis to plot on.
    networkspeci : str
        Current networkspeci (e.g. EBAS|sconco3).
    data_labels : list
        Data arrays to plot.
    base_plot_type : str
        Plot type, without statistical information.
    plot_characteristics : dict
        Plot characteristics.
    plot_options : list
        Options to configure plots.
    chunk_stat : str, optional
        Chunk statistic for timeseries smoothing (default is None).
    chunk_resolution : str, optional
        Chunk resolution for timeseries smoothing (default is None).
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
    if read_instance.mode not in ['report', 'library']:
        chunk_stat = canvas_instance.timeseries_chunk_stat.currentText()
        chunk_resolution = canvas_instance.timeseries_chunk_resolution.currentText()
        chunk_stat = None if chunk_stat == 'None' else chunk_stat
        chunk_resolution = None if chunk_resolution == 'None' else chunk_resolution
    
    # chunk timeseries
    if (chunk_stat is not None) and (chunk_resolution is not None):
        timeseries_data = create_statistical_timeseries(read_instance, canvas_instance, chunk_stat, 
                                                        chunk_resolution, networkspeci, cut_data_labels, bias)
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
        if read_instance.mode not in ['report', 'library']:
            canvas_instance.plotting.track_plot_elements(data_label, base_plot_type, 'smooth', smooth_line, bias=bias)

def threshold(read_instance, canvas_instance, relevant_axis, networkspeci, base_plot_type, 
              plot_characteristics):
    """
    Add threshold line(s) to a plot.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axis : object
        Axis to plot on.
    networkspeci : str
        Current networkspeci (e.g. EBAS|sconco3).
    base_plot_type : str
        Plot type, without statistical information.
    plot_characteristics : dict
        Plot characteristics.
    """

    # get exceendance value
    threshold = exceedance_lim(read_instance, networkspeci)

    # draw vertical line
    if base_plot_type in ['timeseries', 'scatter', 'periodic', 'periodic-violin', 'boxplot']:
        threshold_line = relevant_axis.axhline(y=threshold, 
                                               **plot_characteristics['threshold_line'])
    
    # draw horizontal line
    if base_plot_type in ['distribution', 'scatter']:
        threshold_line = relevant_axis.axvline(x=threshold, 
                                               **plot_characteristics['threshold_line'])

    # track plot elements if using dashboard 
    if read_instance.mode not in ['report', 'library']:
        canvas_instance.plotting.track_plot_elements('ALL', base_plot_type, 'threshold', [threshold_line], bias=False)

def annotation(read_instance, canvas_instance, relevant_axis, networkspeci, data_labels, base_plot_type, 
               plot_characteristics, plot_options, plot_z_statistic_sign='absolute'):
    """
    Add statistical annotations to a plot.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axis : object
        Axis to plot on.
    networkspeci : str
        Current networkspeci (e.g. EBAS|sconco3).
    data_labels : list
        Data arrays to plot.
    base_plot_type : str
        Plot type, without statistical information.
    plot_characteristics : dict
        Plot characteristics.
    plot_options : list
        Options to configure plots.
    plot_z_statistic_sign : str, optional
        Sign of plotted z statistic, either 'absolute' or 'bias' (default is 'absolute').
    """

    # initialise list of strs to annotate
    str_to_annotate = []

    # add annotation text
    if base_plot_type == 'fairmode-target':
        bias = False
        if hasattr(canvas_instance.plotting, 'faimode_target_annotate_text'):
            str_to_annotate = canvas_instance.plotting.faimode_target_annotate_text
            colours = canvas_instance.plotting.faimode_target_annotate_colour
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
            if (plot_characteristics['annotate_text']['mod_labels']):
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
        if read_instance.mode not in ['report', 'library']:
            canvas_instance.plotting.track_plot_elements('ALL', base_plot_type, 'annotate', [bbox], bias=bias)
    else:
        msg = '{} could not be annotated'.format(base_plot_type)
        show_message(read_instance, msg)

def model_domain(canvas_instance, relevant_axis, data_labels, map_extent):
    """
    Plot model domain extents on a map.

    Parameters
    ----------
    canvas_instance : object
        Instance of class Canvas or Report.
    relevant_axis : object
        Axis to plot on.
    data_labels : list
        Data arrays to plot.
    map_extent : list
        List of map extent bounds [lonmin, lonmax, latmin, latmax].
    """

    #get model domain polygons
    grid_edge_polygons = canvas_instance.plotting.make_model_domain_polygons(data_labels=data_labels) 

    # plot grid edge polygons on map
    for grid_edge_polygon in grid_edge_polygons:
        relevant_axis.add_patch(grid_edge_polygon)

    # if map extent is not set then re-set automatic limits based now domain is plotted.
    if not map_extent:
        relevant_axis.relim(visible_only=True)
        relevant_axis.autoscale(tight=False)