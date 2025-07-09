""" Functions to format the axes """

import copy
import os
from PIL import Image

import cartopy
import cartopy.feature as cfeature
from datetime import datetime
import matplotlib
import matplotlib as mpl 
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.dates import num2date
from matplotlib import ticker
import numpy as np
from packaging.version import Version
import pandas as pd

from .plot_aux import get_land_polygon_resolution, set_map_extent
from .plot_options import annotation, experiment_domain, linear_regression, log_axes, smooth, threshold
from .statistics import get_z_statistic_info

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

Image.MAX_IMAGE_PIXELS = None

def set_equal_axes(ax, plot_options, plot_characteristics, base_plot_type):
    """ Set equal aspect and limits (useful for scatter plots). 

        :param ax: axis to set equal axes
        :type ax: object
        :param plot_options: active plot options 
        :type plot_options: list
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
    """

    # set equal aspect
    ax.set_aspect(aspect='equal', adjustable='box')

    if len(ax.lines) == 0 and base_plot_type == 'scatter':
        return None
    
    if ("xlim" not in plot_characteristics) and ("ylim" not in plot_characteristics):
        # get min and max values for axes from plotted data
        xmin, xmax = get_data_lims(ax, 'xlim', plot_options)
        ymin, ymax = get_data_lims(ax, 'ylim', plot_options)

        # compare min and max lims across axes
        if xmin < ymin:
            axmin = xmin
        else:
            axmin = ymin
        if xmax > ymax:
            axmax = xmax
        else:
            axmax = ymax
        
        # set equal lims
        ax.set_xlim(axmin, axmax)
        ax.set_ylim(axmin, axmax)
    elif ("xlim" not in plot_characteristics) and ("ylim" in plot_characteristics):
        # set xlim as ylim if ylim is passed
        if isinstance(plot_characteristics["ylim"], dict):
            ax.set_xlim(**plot_characteristics["ylim"])
        else:
            ax.set_xlim(plot_characteristics["ylim"])
    elif ("xlim" in plot_characteristics) and ("ylim" not in plot_characteristics):
        # set ylim as ylim if xlim is passed
        if isinstance(plot_characteristics["xlim"], dict):
            ax.set_ylim(**plot_characteristics["xlim"])
        else:
            ax.set_ylim(plot_characteristics["xlim"])
    elif ("xlim" in plot_characteristics) and ("ylim" in plot_characteristics):
        # set both limits
        if isinstance(plot_characteristics["xlim"], dict):
            ax.set_xlim(**plot_characteristics["xlim"])
        else:
            ax.set_xlim(plot_characteristics["xlim"])
        if isinstance(plot_characteristics["ylim"], dict):
            ax.set_ylim(**plot_characteristics["ylim"])
        else:
            ax.set_ylim(plot_characteristics["ylim"])

    return None


def harmonise_xy_lims_paradigm(canvas_instance, read_instance, relevant_axs, base_plot_type, plot_characteristics, 
                               plot_options, xlim=None, ylim=None, relim=False, autoscale=False, autoscale_x=False, 
                               autoscale_y=False, bias_centre=False, harmonise=True):
    """ Harmonise xy limits across paradigm of plot type, unless axis limits have been defined.
    
        :param canvas_instance: Instance of class Canvas or Report
        :type canvas_instance: object
        :param read_instance: Instance of class Dashboard or Report
        :type read_instance: object
        :param relevant_axs: relevant axes
        :type relevant_axs: list
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param plot_options: Options to configure plots
        :type plot_options: list
        :param xlim: xlimits to set
        :type xlim: dict
        :param ylim: ylimits to set
        :type ylim: dict
        :param relim: turn on relimiting of axes limits (when updating plotted data on axis)
        :type relim: boolean
        :param autoscale: Autoscale the axis view to the data (both x and y axes)
        :type autoscale: boolean
        :param autoscale_x: Autoscale the x axis view to the data
        :type autoscale_x: boolean
        :param autoscale_y: Autoscale the x axis view to the data
        :type autoscale_y: boolean
        :param bias_centre: centre bias plots at 0 on the y axis
        :type bias_centre: boolean
        :param harmonise: switch for harmonisation of axes 
        :type harmonise: boolean
    """
    
    # initialise arrays to save lower and upper limits in all axes
    all_xlim_lower = []
    all_xlim_upper = []
    all_ylim_lower = []
    all_ylim_upper = []

    # initialise variables for setting axis limits
    xlim_min = None
    xlim_max = None
    ylim_min = None
    ylim_max = None
    
    # transform axis dict or str to list
    if not isinstance(relevant_axs, list):
        # if changes only apply to one axis, put it in list
        if not isinstance(relevant_axs, dict):
            relevant_axs = [relevant_axs]
        # transform dictionaries into lists
        else:
            relevant_axs = [relevant_axs[relevant_temporal_resolution] for 
                            relevant_temporal_resolution in read_instance.relevant_temporal_resolutions]

    # get mapped resolution per axis for periodic plots
    if base_plot_type in ['periodic', 'periodic-violin']:
        mapped_resolutions = read_instance.relevant_temporal_resolutions*(int(len(relevant_axs)/len(
            read_instance.relevant_temporal_resolutions)))

    # remove any axes from relevant_axs which are not active (only for report and library)
    if (read_instance.report) or (read_instance.library):
        relevant_axs_active = []
        mapped_resolutions_active = []
        for ax_ii, ax in enumerate(relevant_axs):
            if ax.axison:
                relevant_axs_active.append(ax)
                if base_plot_type in ['periodic', 'periodic-violin']:
                    mapped_resolutions_active.append(mapped_resolutions[ax_ii])
    else:
        relevant_axs_active = relevant_axs
        if base_plot_type in ['periodic', 'periodic-violin']:
            mapped_resolutions_active = mapped_resolutions

    # get lower and upper limits across all relevant axes
    for ax in relevant_axs_active:
        if 'equal_aspect' in plot_characteristics:
            if plot_characteristics['equal_aspect']:
                set_equal_axes(ax, plot_options, plot_characteristics, base_plot_type)
                continue
        else:
            ax.set_aspect('auto')

        if relim:
            ax.relim(visible_only=True)
        if autoscale:
            ax.autoscale(tight=False)
        if autoscale_x:
            ax.autoscale(axis='x', tight=False)
        if autoscale_y:
            ax.autoscale(axis='y', tight=False)

        if (xlim is None) and ('xlim' not in plot_characteristics):
            if base_plot_type not in ['timeseries', 'boxplot', 'periodic', 'periodic-violin']:
                xlim_lower, xlim_upper = ax.get_xlim()
            elif base_plot_type == 'timeseries':
                xlim_lower, xlim_upper = get_no_margin_lim(ax, 'xlim')
                try:
                    xlim_lower = num2date(xlim_lower).replace(tzinfo=None)
                    xlim_upper = num2date(xlim_upper).replace(tzinfo=None)
                except ValueError:
                    continue

            if base_plot_type not in ['boxplot', 'periodic', 'periodic-violin']:
                all_xlim_lower.append(xlim_lower)
                all_xlim_upper.append(xlim_upper)

        if (ylim is None) and ('ylim' not in plot_characteristics):
            ylim_lower, ylim_upper = ax.get_ylim()
            all_ylim_lower.append(ylim_lower)
            all_ylim_upper.append(ylim_upper)

    # get minimum and maximum from all axes and set limits
    for ax_ii, ax in enumerate(relevant_axs_active):

        # get xlim
        if ax_ii == 0:
            set_xlim = False
            if ((xlim is None) and ('xlim' not in plot_characteristics) and (len(all_xlim_lower) > 0) and 
               (len(all_xlim_upper) > 0)):
                if base_plot_type not in ['boxplot', 'periodic', 'periodic-violin']:
                    xlim_min = np.min(all_xlim_lower)
                    xlim_max = np.max(all_xlim_upper)
                    xlim = {'left':xlim_min, 'right':xlim_max}
                    if harmonise:
                        set_xlim = True
            elif 'xlim' in plot_characteristics:
                xlim = plot_characteristics['xlim']
                set_xlim = True

        # set xlim
        if ((set_xlim) and ('equal_aspect' not in plot_characteristics) and
           (base_plot_type not in ['timeseries', 'boxplot', 'periodic', 'periodic-violin'])):
            if isinstance(xlim, dict):
                ax.set_xlim(**xlim)
            else:
                ax.set_xlim(xlim)

        # get ylim
        if ax_ii == 0:
            set_ylim = False
            if ((ylim is None) and ('ylim' not in plot_characteristics) and (len(all_ylim_lower) > 0) and 
               (len(all_ylim_upper) > 0)):
                ylim_min = np.min(all_ylim_lower) 
                ylim_max = np.max(all_ylim_upper)
                # if have bias_centre option, centre around zero
                if ('bias' in plot_options) & (bias_centre):                    
                    if np.abs(np.max(all_ylim_upper)) >= np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.max(all_ylim_upper))
                        ylim_max = np.abs(np.max(all_ylim_upper))
                    elif np.abs(np.max(all_ylim_upper)) < np.abs(np.min(all_ylim_lower)):
                        ylim_min = -np.abs(np.min(all_ylim_lower))
                        ylim_max = np.abs(np.min(all_ylim_lower))
                ylim = {'bottom':ylim_min, 'top':ylim_max}
                if harmonise:
                    set_ylim = True
            elif 'ylim' in plot_characteristics:
                ylim = plot_characteristics['ylim']
                set_ylim = True

        # set ylim
        if (set_ylim) and ('equal_aspect' not in plot_characteristics):
            if isinstance(ylim, dict):
                ax.set_ylim(**ylim)
            else:
                ax.set_ylim(ylim)
        
    # get minimum and maximum from all axes and set limits for periodic plots
    if base_plot_type in ['periodic', 'periodic-violin']:
        if (xlim is None) and ('xlim' not in plot_characteristics):
            for temporal_resolution, sub_ax in zip(mapped_resolutions_active, relevant_axs_active):
                # adjust plot x axis to have correct margin on edges
                xlim_lower, xlim_upper = sub_ax.get_xlim()
                first_valid_x = canvas_instance.periodic_xticks[temporal_resolution][(np.abs(
                    canvas_instance.periodic_xticks[temporal_resolution] - xlim_lower)).argmin()]
                last_valid_x = canvas_instance.periodic_xticks[temporal_resolution][(np.abs(
                    canvas_instance.periodic_xticks[temporal_resolution] - xlim_upper)).argmin()]
                if temporal_resolution == 'hour':
                    xlim_lower = first_valid_x - 0.65
                    xlim_upper = last_valid_x + 0.65
                elif temporal_resolution == 'dayofweek':
                    xlim_lower = first_valid_x - 0.55
                    xlim_upper = last_valid_x + 0.55
                elif temporal_resolution == 'month':
                    xlim_lower = first_valid_x - 0.55
                    xlim_upper = last_valid_x + 0.55
                xlim = {'left':xlim_lower, 'right':xlim_upper}
                sub_ax.set_xlim(**xlim)
        elif 'xlim' in plot_characteristics:
            xlim = plot_characteristics['xlim']
            for temporal_resolution, sub_ax in zip(mapped_resolutions_active, relevant_axs_active):
                if isinstance(xlim, dict):
                    sub_ax.set_xlim(**xlim)
                else:
                    sub_ax.set_xlim(xlim)

        # if harmonisation is off, and ylim not manually set, 
        # ensure harmonisation is at least done for a plot across resolutions
        if (not harmonise) and (not set_ylim):
            
            current_resolutions = []
            current_axs = []
            current_ylim_lower = []
            current_ylim_upper = []

            for ax_ii, (temporal_resolution, sub_ax) in enumerate(zip(mapped_resolutions_active, relevant_axs_active)):

                # get temporal resolution of next axis
                if ax_ii != (len(relevant_axs_active)-1):
                    next_temporal_resolution = mapped_resolutions_active[ax_ii+1]
                else:
                    next_temporal_resolution = None

                # if resolution not yet in current resolutions, then add information for it
                if temporal_resolution not in current_resolutions:
                     ylim_lower, ylim_upper = sub_ax.get_ylim()
                     current_resolutions.append(temporal_resolution)
                     current_axs.append(sub_ax)
                     current_ylim_lower.append(ylim_lower)
                     current_ylim_upper.append(ylim_upper)

                # if next resolution already in current resolutions, or on last axis, then set ylim for relevant axes
                if (next_temporal_resolution in current_resolutions) or (ax_ii == (len(relevant_axs_active)-1)):

                    ylim_min = np.min(current_ylim_lower) 
                    ylim_max = np.max(current_ylim_upper)
                    # if have bias_centre option, centre around zero
                    if ('bias' in plot_options) & (bias_centre):                    
                        if np.abs(np.max(current_ylim_upper)) >= np.abs(np.min(current_ylim_lower)):
                            ylim_min = -np.abs(np.max(current_ylim_upper))
                            ylim_max = np.abs(np.max(current_ylim_upper))
                        elif np.abs(np.max(current_ylim_upper)) < np.abs(np.min(current_ylim_lower)):
                            ylim_min = -np.abs(np.min(current_ylim_lower))
                            ylim_max = np.abs(np.min(current_ylim_lower))
                    ylim = {'bottom':ylim_min, 'top':ylim_max}
                    for current_ax in current_axs:
                        current_ax.set_ylim(**ylim)

                    # reset lists
                    current_resolutions = []
                    current_axs = []
                    current_ylim_lower = []
                    current_ylim_upper = []
      
    # get minimum and maximum from all axes and set limits for timeseries
    elif base_plot_type == 'timeseries':
        if (plot_characteristics['xtick_alteration']['define']) and (xlim):
            
            # get left and right
            if isinstance(xlim, dict):
                left = xlim['left']
                right = xlim['right']
            else:
                left = xlim[0]
                right = xlim[1]

            # get number of days
            n_days = (right - left).days

            first_step = plot_characteristics['xtick_alteration']['first_step']
            last_step = plot_characteristics['xtick_alteration']['last_step']
            n_slices = plot_characteristics['xtick_alteration']['n_slices']
            overlap = plot_characteristics['xtick_alteration']['overlap']

            # if there's more than 3 months, define time slices as the first day of the month
            if n_days >= 3 * 30:
                # remove hours, minutes and seconds from the right and end dates
                left, right = datetime(left.year, left.month, left.day), datetime(right.year, right.month, right.day)

                # get the first and last days of each month
                months_start = pd.date_range(left, right, freq='MS')
                if Version(matplotlib.__version__) < Version("3.9"):
                    months_end = pd.date_range(left, right, freq='M')
                else:
                    months_end = pd.date_range(left, right, freq='ME')

                # set steps as the start of the months
                steps = months_start

                # remove last day if there's less than a n days difference
                if last_step and 0 < (right - months_end[-1]).days <= overlap:
                    steps = steps[:-1]

                # remove first day if there's less than a n days difference
                if first_step and 0 < (months_start[0] - left).days <= overlap:
                    steps = steps[1:]

                # get xticks
                slices = int(np.ceil(len(steps) / int(n_slices+1)))
                xticks = steps[0::slices]

                # transform to numpy.datetime64
                if not isinstance(xticks[0], np.datetime64):
                    xticks = [np.datetime64(x, 'D') for x in xticks]
                if not isinstance(right, np.datetime64):
                    right = np.datetime64(right)

                # add last step to xticks
                if last_step and (xticks[-1] != right):
                    xticks = np.append(xticks, right)
                
                # add first step to xticks
                if first_step and (xticks[0] != left):
                    xticks = np.insert(xticks, 0, left)
        
            else:
                # round up the limit hours to the whole hour
                left = pd.to_datetime(left).ceil('H')   
                right = pd.to_datetime(right).floor('H')

                # set frequency to hourly when there's less than 7 days 
                freq = 'h' if n_days < 7 else 'D'

                # get all the dates in the frequency
                steps = pd.date_range(left, right, freq=freq)  

                # get n_periods dates from all_ticks
                periods = n_slices + int(first_step) + int(last_step) + 1

                # compute number of ticks to select, it can't exceed available steps
                n_ticks = min(periods, len(steps))
                xticks = steps[np.linspace(0, len(steps) - 1, n_ticks, dtype=int)]
               
            # show hours if number of days is less than 7
            if n_days < 7:
                ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d %Hh'))
            else:
                ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%Y-%m-%d'))
            
            # set modified xticks
            for ax in relevant_axs_active:
                ax.xaxis.set_ticks(xticks)

            # get the date format to create the margin
            clip_left = mdates.date2num(left)
            clip_right = mdates.date2num(right)

            # get the len of the original y axis
            ylen = ax.get_ylim()[1] - ax.get_ylim()[0]

            # create the rectangle that will define the margin
            clip_rect = mpatches.Rectangle(
                (clip_left, ax.get_ylim()[0] + ylen*0.05),
                clip_right - clip_left,
                ylen * 0.9,       
                transform=ax.transData
            )

            # set the margin
            for ts in canvas_instance.plotting.timeseries_plot:
                ts[0].set_clip_path(clip_rect)

def set_axis_title(read_instance, relevant_axis, title, plot_characteristics):
    """ Set title of plot axis.

        :param read_instance: Instance of class Dashboard or Report
        :type read_instance: object
        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param title: axis title
        :type title: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
    """    

    # return if title is empty str
    if title == '':
        return

    # get appropriate axis for plotting label for plots with multiple sub-axes
    axs_to_set_title = []
    if isinstance(relevant_axis, dict):
        # reorder dict to show axis title in monthly plot and not in DoW for daily plots
        relevant_dict = {key : relevant_axis[key] for key in ['hour', 'month', 'dayofweek']}
        for relevant_temporal_resolution, sub_ax in relevant_dict.items():
            if relevant_temporal_resolution in read_instance.relevant_temporal_resolutions:
                axs_to_set_title.append(sub_ax)
                break
    elif isinstance(relevant_axis, list):
        axs_to_set_title.append(relevant_axis[0])
    else:
        axs_to_set_title.append(relevant_axis)

    # set title for appropriate axes
    axis_title_characteristics = copy.deepcopy(plot_characteristics['axis_title'])
    axis_title_characteristics['label'] = title
    for relevant_axis in axs_to_set_title:
        relevant_axis.set_title(**axis_title_characteristics)


def set_axis_label(relevant_axis, label_ax, label, plot_characteristics, 
                   relevant_temporal_resolutions=None):
    """ Set label of plot axis.

        :param relevant_axis: Axis to plot on 
        :type relevant_axis: object
        :param label_ax: which axis to set label of
        :type label_ax: str
        :param label: axis label
        :type label: str
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param relevant_temporal_resolutions: list of relevant temporal resolutions  
        :type relevant_temporal_resolutions: list
    """

    # return if label is empty str
    if label == '':
        return

    # define default argument mutables
    if relevant_temporal_resolutions is None:
        relevant_temporal_resolutions = ['hour', 'month']

    # get appropriate axis for plotting label for plots with multiple sub-axes (hour and month axes)
    axs_to_set_label = []
    if isinstance(relevant_axis, dict):
        for relevant_temporal_resolution, sub_ax in relevant_axis.items():
            if relevant_temporal_resolution in relevant_temporal_resolutions:
                axs_to_set_label.append(sub_ax)
            # remove day of week axis label if setting ylabel
            if (relevant_temporal_resolution == 'dayofweek') & (label_ax == 'y'):                           
                sub_ax.yaxis.set_tick_params(which='both', labelleft=False)
                sub_ax.set_ylabel('')
    else:
        axs_to_set_label.append(relevant_axis)

    # set label for appropriate axes
    for relevant_axis in axs_to_set_label:
        if label_ax == 'x':
            axis_label_characteristics = copy.deepcopy(plot_characteristics['xlabel'])
            axis_label_characteristics['xlabel'] = label
            relevant_axis.set_xlabel(**axis_label_characteristics)
        elif label_ax == 'y':
            axis_label_characteristics = copy.deepcopy(plot_characteristics['ylabel'])
            axis_label_characteristics['ylabel'] = label
            relevant_axis.set_ylabel(**axis_label_characteristics)


def format_plot_options(canvas_instance, read_instance, relevant_axs, relevant_data_labels, networkspeci, 
                        base_plot_type, plot_type, plot_options, map_extent=False, chunk_stat=None, 
                        chunk_resolution=None):
    """ Function that handles formatting of a plot axis,
        based on given plot options.

        :param canvas_instance: Instance of class Canvas or Report
        :type canvas_instance: object
        :param read_instance: Instance of class Dashboard or Report
        :type read_instance: object
        :param relevant_axs: relevant axes
        :type relevant_axs: list
        :param relevant_data_labels: names of plotted data arrays per axis
        :type relevant_data_labels: list
        :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
        :type networkspeci: str
        :param base_plot_type: Plot type, without statistical information
        :type base_plot_type: str
        :param plot_type: plot type
        :type plot_type: str
        :param plot_options: Options to configure plots
        :type plot_options: list
        :param map_extent: list of map extent bounds [lonmin, lonmax, latmin, latmax]
        :type map_extent: list
        :param chunk_stat: Chunk statistic
        :type chunk_stat: str
        :param chunk_resolution: Chunk resolution
        :type chunk_resolution: str
    """

    # transform axis dict or str to list
    if not isinstance(relevant_axs, list):
        # if changes only apply to one axis, put it in list
        if not isinstance(relevant_axs, dict):
            relevant_axs = [relevant_axs]
        # transform dictionaries into lists
        else:
            relevant_axs = [relevant_axs[relevant_temporal_resolution] for 
                            relevant_temporal_resolution in read_instance.relevant_temporal_resolutions]
            relevant_data_labels = copy.deepcopy(relevant_data_labels) * len(read_instance.relevant_temporal_resolutions)

    # get zstat info (if any)
    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type=plot_type) 

    for relevant_ax_ii, relevant_ax in enumerate(relevant_axs):

        # log axes?
        if 'logx' in plot_options:            
            log_valid = log_validity(relevant_ax, 'logx')
            if log_valid:
                log_axes(relevant_ax, 'logx', canvas_instance.plot_characteristics[plot_type])
            else:
                msg = "It is not possible to log the x-axis "
                msg += "in {0} with negative values.".format(plot_type)
                show_message(read_instance, msg)

        if 'logy' in plot_options:
            log_valid = log_validity(relevant_ax, 'logy')
            if log_valid:
                log_axes(relevant_ax, 'logy', canvas_instance.plot_characteristics[plot_type])
            else:
                msg = "It is not possible to log the y-axis "
                msg += "in {0} with negative values.".format(plot_type)
                show_message(read_instance, msg)

        # domain
        if 'domain' in plot_options:
            if len(read_instance.data_labels) == 1:
                if read_instance.data_labels[0] == read_instance.observations_data_label:
                    msg = "'domain' plot option cannot be made as have no experiments."
                    show_message(read_instance, msg)
                    return
            experiment_domain(canvas_instance, relevant_ax, relevant_data_labels[relevant_ax_ii], map_extent)

        # annotation
        if 'annotate' in plot_options:
            if base_plot_type not in ['heatmap']:
                annotation(canvas_instance, read_instance, relevant_ax, networkspeci, 
                           relevant_data_labels[relevant_ax_ii], base_plot_type, 
                           canvas_instance.plot_characteristics[plot_type],
                           plot_options, plot_z_statistic_sign=z_statistic_sign)
                # annotate on first axis
                if base_plot_type in ['periodic', 'periodic-violin']:
                    break
        
        # regression line
        if 'regression' in plot_options:
            linear_regression(canvas_instance, read_instance, relevant_ax, networkspeci, 
                              relevant_data_labels[relevant_ax_ii], base_plot_type, 
                              canvas_instance.plot_characteristics[plot_type], plot_options)

        # smooth line
        if 'smooth' in plot_options:
            smooth(canvas_instance, read_instance, relevant_ax, networkspeci,
                   relevant_data_labels[relevant_ax_ii], base_plot_type, 
                   canvas_instance.plot_characteristics[plot_type], plot_options,
                   chunk_stat, chunk_resolution)

        # threshold line
        if 'threshold' in plot_options:
            threshold(canvas_instance, read_instance, relevant_ax, networkspeci, base_plot_type,
                      canvas_instance.plot_characteristics[plot_type])
            
            
def format_axis(canvas_instance, read_instance, ax, base_plot_type, plot_characteristics, col_ii=0, last_valid_row=True, 
                last_row_on_page=True, map_extent=False, relevant_temporal_resolutions=None):
    """ Format a plotting axis.
    
        :param canvas_instance: Instance of class Canvas or Report
        :type canvas_instance: object
        :param read_instance: Instance of class Dashboard or Report
        :type read_instance: object
        :param ax: axis object
        :type ax: object
        :param base_plot_type: plot to make, without statistical information
        :type base_plot_type: str  
        :param plot_characteristics: plot characteristics
        :type plot_characteristics: dict
        :param col_ii: column index (for reports)
        :type col_ii: int
        :param last_valid_row: boolean informing if last valid row to plot on (for reports)
        :type last_valid_row: boolean
        :param last_row_on_page: boolean informing if last valid row on page (for reports)
        :type last_row_on_page: boolean
        :param map_extent: list of map extent bounds [lonmin, lonmax, latmin, latmax]
        :type map_extent: list
        :param relevant_temporal_resolutions: list of relevant temporal resolutions
        :type relevant_temporal_resolutions: list
    """

    # define default argument mutables
    if relevant_temporal_resolutions is None:
        relevant_temporal_resolutions = ['hour', 'dayofweek', 'month']

    # get plot characteristics vars
    plot_characteristics_vars = list(plot_characteristics.keys())

    # get appropriate axes for nested axes
    axs_to_format = []
    temporal_resolutions_per_ax = []
    if isinstance(ax, dict):
        for relevant_temporal_resolution, sub_ax in ax.items():
            if relevant_temporal_resolution in relevant_temporal_resolutions:
                axs_to_format.append(sub_ax)
                temporal_resolutions_per_ax.append(relevant_temporal_resolution)
    elif isinstance(ax, list):
        axs_to_format = ax
        temporal_resolutions_per_ax = [""] * len(ax)
    else:
        axs_to_format.append(ax)
        temporal_resolutions_per_ax.append('')

    # iterate though relevant axes (and relevant temporal resolutions for periodic plots)
    for ax_to_format, relevant_temporal_resolution in zip(axs_to_format, temporal_resolutions_per_ax): 

        # set axis ticks and gridlines below all artists
        ax_to_format.set_axisbelow(True)

        # make axis ylabel (only on leftmost column of visible axes)?
        #if 'axis_title' in plot_characteristics_vars:
        #    ax_to_format.set_title(**plot_characteristics['axis_title'])

        # make axis xlabel?
        #if 'xlabel' in plot_characteristics_vars:
        #    ax_to_format.set_xlabel(**plot_characteristics['xlabel'])

        # make axis ylabel (only on leftmost column of visible axes)?
        #if 'ylabel' in plot_characteristics_vars:
        #    ax_to_format.set_ylabel(**plot_characteristics['ylabel'])

        # set xtick params ?
        if 'xtick_params' in plot_characteristics_vars:
            ax_to_format.xaxis.set_tick_params(**plot_characteristics['xtick_params'])

        # set ytick params ?
        if 'ytick_params' in plot_characteristics_vars:
            ax_to_format.yaxis.set_tick_params(**plot_characteristics['ytick_params'])

        # if are sharing xticks, and not on last row on page/last
        # valid row, then ensure current axis xticks are hidden
        if ('xtick_share' in plot_characteristics_vars) and (not last_valid_row) and (not last_row_on_page):
            plt.setp(ax_to_format.get_xticklabels(), visible=False)

        # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
        if ('ytick_share' in plot_characteristics_vars) and (col_ii != 0):
            plt.setp(ax_to_format.get_yticklabels(), visible=False)

        # set xlim?
        if 'xlim' in plot_characteristics_vars:
            if isinstance(plot_characteristics['xlim'], dict):
                ax_to_format.set_xlim(**plot_characteristics['xlim'])
            else:
                ax_to_format.set_xlim(plot_characteristics['xlim'])

        # set ylim? 
        if 'ylim' in plot_characteristics_vars:
            if isinstance(plot_characteristics['ylim'], dict):
                ax_to_format.set_ylim(**plot_characteristics['ylim'])
            else:
                ax_to_format.set_ylim(plot_characteristics['ylim'])

        # add gridlines (x and y)?
        if 'grid' in plot_characteristics_vars:
            ax_to_format.grid(**plot_characteristics['grid'])

        # add x gridlines?
        if 'xgrid' in plot_characteristics_vars:
            ax_to_format.xaxis.grid(**plot_characteristics['xgrid'])

        # add y gridlines?
        if 'ygrid' in plot_characteristics_vars:
            ax_to_format.yaxis.grid(**plot_characteristics['ygrid'])

        # set x axis decimal places?
        if 'round_decimal_places' in plot_characteristics_vars:
            if 'x' in plot_characteristics['round_decimal_places']:
                ax_to_format.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.{}f'.format(plot_characteristics['round_decimal_places']['x'])))

        # set y axis decimal places?
        if 'round_decimal_places' in plot_characteristics_vars:
            if 'y' in plot_characteristics['round_decimal_places']:
                ax_to_format.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.{}f'.format(plot_characteristics['round_decimal_places']['y'])))

        # remove spines?
        if 'remove_spines' in plot_characteristics_vars:
            for side in plot_characteristics['remove_spines']:
                ax_to_format.spines[side].set_visible(False)

            for side in list(set(['top', 'bottom', 'right', 'left']).symmetric_difference(plot_characteristics['remove_spines'])):
                ax_to_format.spines[side].set_visible(True)

        # handle formatting specific to plot types
        if base_plot_type in ['periodic','periodic-violin']:

            # add axis resolution label 
            ax_to_format.annotate(canvas_instance.periodic_labels[relevant_temporal_resolution], **plot_characteristics['label'])

            # set plotted x axis ticks/labels (if 'hour' aggregation --> a numeric tick every 3 hours)
            if relevant_temporal_resolution == 'hour':
                plot_characteristics['xticks'] = canvas_instance.periodic_xticks[relevant_temporal_resolution][::3]
                ax_to_format.set_xticks(plot_characteristics['xticks'])
            else:
                plot_characteristics['xticks'] = canvas_instance.periodic_xticks[relevant_temporal_resolution]
                ax_to_format.set_xticks(plot_characteristics['xticks'])
                ax_to_format.set_xticklabels([canvas_instance.temporal_axis_mapping_dict['short'][relevant_temporal_resolution][xtick] for xtick
                                              in canvas_instance.periodic_xticks[relevant_temporal_resolution]])
        
        elif base_plot_type == 'map':

            # set map background

            # providentia default background
            if plot_characteristics['background'] == 'providentia':
                feature = cfeature.NaturalEarthFeature(category='physical', name='land',
                                                       scale=get_land_polygon_resolution(canvas_instance.plot_characteristics_templates['map']['map_coastline_resolution']), 
                                                       **canvas_instance.plot_characteristics_templates['map']['land_polygon'])
                ax_to_format.add_feature(feature)

            # shaded relief (cartopy default)
            elif plot_characteristics['background'] == 'shaded_relief':
                ax_to_format.stock_img()

            # other type of map background
            else:
                # check file for background exists
                background_fname = join(CURRENT_PATH, "resources/{}.png".format(plot_characteristics['background']))
                if os.path.isfile(background_fname):
                    img = plt.imread(background_fname)
                    img_extent = (-180, 180, -90, 90)
                    ax_to_format.imshow(img, origin='upper', extent=img_extent, transform=canvas_instance.datacrs)
                else:
                    msg = "Specified map background file cannot be found."
                    show_message(read_instance, msg)

            # add gridlines ?
            if 'gridlines' in plot_characteristics_vars:
                gridlines_characteristics = plot_characteristics['gridlines']
                ax_to_format.gridlines(crs=canvas_instance.datacrs, 
                                       **gridlines_characteristics)

            # set map extent (if wanted)
            if map_extent:
                set_map_extent(canvas_instance, ax_to_format, map_extent)

        elif base_plot_type == 'fairmode-target':

            # update axis labels
            ax_to_format.set_xticks(**plot_characteristics['xticks'])
            ax_to_format.set_yticks(**plot_characteristics['yticks'])

def get_no_margin_lim(ax, lim):
    """ Get true limits of plot area (with no margins)

        :param ax: axis to get limits
        :type ax: object
        :param lim: xlim or ylim
        :type lim: str
        :return: lower_lim, upper_lim
        :rtype: float32, float32
    """

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


def get_data_lims(ax, lim, plot_options):
    """ Get x and y limits of a plot axis from plotted data

        :param ax: axis to get limits
        :type ax: object
        :param lim: xlim or ylim
        :type lim: str
        :param plot_options: Options to configure plots
        :type plot_options: list
        return: lower_lim, upper_lim
        :rtype: float32, float32
    """

    # get min and max values for axis
    lines = []
    for i, line in enumerate(ax.lines):
        if lim == 'xlim':
            line_data = line.get_xdata()
        elif lim == 'ylim':
            line_data = line.get_ydata()
        if ((list(line_data) == [0, 1]) 
            or (list(line_data) == [0, 0.5] 
            or (list(line_data) == [[1.0], [1.0]]))):
            continue
        lines.extend(line_data)

    # if log of an axis is active, remove values <= 0 to ensure limits are obtained correctly
    if (('logx' in plot_options) and (lim == 'xlim')) or (('logy' in plot_options) and (lim == 'ylim')):
        lines = np.array(lines)
        above_zero = lines > 0
        lines = lines[above_zero]

    # get min/max across line artists
    if len(lines) == 0:
        return np.NaN, np.NaN
    else:
        lower_lim = np.nanmin(lines)
        upper_lim = np.nanmax(lines)
        return lower_lim, upper_lim


def log_validity(ax, log_ax):
    """ Determine if log operation for a given axes is valid (no values < 0).
    
        :param ax: relevant axis
        :type ax: object
        :param log_ax: which axis to log
        :type log_ax: str
        :return: validity to log axis
        :rtype: boolean
    """

    if log_ax == 'logx':
        lower_lim, _ = get_data_lims(ax, 'xlim', ['logx'])
    elif log_ax == 'logy':
        lower_lim, _ = get_data_lims(ax, 'ylim', ['logy'])

    if lower_lim < 0:
        validity = False
    else:
        validity = True

    return validity
