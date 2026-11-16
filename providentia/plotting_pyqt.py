import copy

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np

from .calculate import ModBias
from .plot_aux import get_fairmode_RV_exceendance, get_taylor_diagram_ghelper
from .plot_formatting import (
    format_axis,
    format_plot_options,
    harmonise_xy_lims_paradigm,
    set_axis_label,
    set_axis_title,
)
from .plotting import Plotting, fairmode_settings
from .read_aux import get_valid_metadata
from .statistics import get_fairmode_data
from .warnings_prv import show_message


class PlottingPyQt(Plotting):
    """Class for plotting with matplotlib: PyQt dashboard (desktop app), report and library"""

    def create_figure(
        self, plot_type, base_plot_type, networkspeci, width=None, height=None, legend=True
    ):
        """
        Creates the matplotlib figure and axes of a library plot.

        Parameters
        ----------
        plot_type : str
            Plot type, with its options (e.g. timeseries_bias).
        base_plot_type : str
            Plot type without options (e.g. timeseries).
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        width : float, optional
            Width of the figure, by default None (the one in the plot characteristics).
        height : float, optional
            Height of the figure, by default None (the one in the plot characteristics).
        legend : bool, optional
            Not used, as the legend of matplotlib figures is made apart.

        Returns
        -------
        fig : matplotlib.figure.Figure
            Figure.
        relevant_ax : object
            Axis to plot on (a dictionary of axes per temporal resolution for periodic plots,
            a list of axes for fairmode-statsummary).
        main_gs : matplotlib.gridspec.GridSpec
            Grid of the figure: the plot is in its first row, the second is left for the legend.
        """

        plot_characteristics = self.canvas_instance.plot_characteristics

        if (width is not None) and (height is not None):
            fig = plt.figure(figsize=(width, height))
        else:
            msg = "Width and/or height have not been passed. The default values will be set."
            show_message(self.read_instance, msg)
            fig = plt.figure(figsize=plot_characteristics[plot_type]["figsize"])

        # create axes
        main_gs = gridspec.GridSpec(
            2, 1, **plot_characteristics[plot_type]["main_gs"]
        )
        if base_plot_type == "map":
            ax = fig.add_subplot(main_gs[0], projection=self.canvas_instance.plotcrs)
        elif base_plot_type == "taylor":
            reference_stddev = 7.5
            ghelper = get_taylor_diagram_ghelper(
                reference_stddev, plot_characteristics[plot_type]
            )
            ax = fig.add_subplot(
                main_gs[0], axes_class=fa.FloatingAxes, grid_helper=ghelper
            )
        elif base_plot_type == "fairmode-statsummary":
            # get current species
            speci = networkspeci.split("|")[1]

            # get number of rows and columns
            ncols = 4
            nrows = 8 if speci in ["sconco3", "sconcno2", "pm10"] else 7

            # create gridspec and add it to a list
            gs = gridspec.GridSpecFromSubplotSpec(
                nrows,
                ncols,
                subplot_spec=main_gs[0],
                **plot_characteristics["fairmode-statsummary"]["gridspec_kw"],
            )
            ax = [fig.add_subplot(gs[i, j]) for i in range(nrows) for j in range(ncols)]
        else:
            ax = fig.add_subplot(main_gs[0])

        if base_plot_type in ["periodic", "periodic-violin"]:
            gs = gridspec.GridSpecFromSubplotSpec(
                100, 100, subplot_spec=ax.get_subplotspec()
            )
            grid_dict = dict()
            grid_dict["hour"] = fig.add_subplot(gs[:46, :])
            grid_dict["dayofweek"] = fig.add_subplot(gs[54:, 64:])
            grid_dict["month"] = fig.add_subplot(gs[54:, :62])
            ax.spines["top"].set_color("none")
            ax.spines["bottom"].set_color("none")
            ax.spines["left"].set_color("none")
            ax.spines["right"].set_color("none")
            ax.tick_params(
                labelcolor="w", top=False, bottom=False, left=False, right=False
            )
            relevant_ax = grid_dict
        else:
            relevant_ax = ax

        # adjust margings and subplot spacing if defined
        if "subplots_adjust" in plot_characteristics[plot_type]:
            fig.subplots_adjust(
                **plot_characteristics[plot_type]["subplots_adjust"]
            )

        return fig, relevant_ax, main_gs

    def format_figure(
        self,
        fig,
        relevant_ax,
        main_gs,
        plot_type,
        base_plot_type,
        networkspeci,
        plot_options,
        title=None,
        xlabel=None,
        ylabel=None,
        relevant_data_labels=None,
        legend_data_labels=None,
        legend=True,
        set_obs_legend=True,
        cb=True,
        zstat=None,
        speci=None,
        labela="",
        map_extent=None,
        chunk_stat=None,
        chunk_resolution=None,
    ):
        """
        Sets the title and the axis labels of a library plot, formats its axes and adds its
        legend and colourbar. It is only for the library, whose legend and colourbar methods
        it calls.

        Parameters
        ----------
        fig : matplotlib.figure.Figure
            Figure of the plot.
        relevant_ax : object
            Axis of the plot.
        main_gs : matplotlib.gridspec.GridSpec
            Grid of the figure, whose second row is for the legend.
        plot_type : str
            Plot type, with its options (e.g. timeseries_bias).
        base_plot_type : str
            Plot type without options (e.g. timeseries).
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        plot_options : list
            Options to configure plot.
        title : str, optional
            Title of the plot.
        xlabel : str, optional
            Label of the x axis.
        ylabel : str, optional
            Label of the y axis.
        relevant_data_labels : list, optional
            Data labels plotted.
        legend_data_labels : list, optional
            Data labels to show in the legend.
        legend : bool, optional
            Indicates if legend appears on plot.
        set_obs_legend : bool, optional
            Indicates if observations appear on legend.
        cb : bool, optional
            Indicates if colourbar appears on plot.
        zstat : str, optional
            Statistic of the plot.
        speci : str, optional
            Current species.
        labela : str, optional
            Label of first dataset.
        map_extent : list, optional
            Extent of the map.
        chunk_stat : str, optional
            Chunk statistic.
        chunk_resolution : str, optional
            Chunk resolution.
        """

        plot_characteristics = self.canvas_instance.plot_characteristics[plot_type]

        # overwrite passed title, xlabels and ylabels
        if title:
            set_axis_title(self.read_instance, relevant_ax, title, plot_characteristics)
        if xlabel:
            set_axis_label(relevant_ax, "x", xlabel, plot_characteristics)
        if ylabel:
            set_axis_label(relevant_ax, "y", ylabel, plot_characteristics)

        # format plot axis/axes
        format_axis(
            self.read_instance,
            self.canvas_instance,
            relevant_ax,
            base_plot_type,
            plot_characteristics,
            map_extent=map_extent,
        )

        # format plot options
        format_plot_options(
            self.read_instance,
            self.canvas_instance,
            relevant_ax,
            [relevant_data_labels],
            networkspeci,
            base_plot_type,
            plot_type,
            plot_options,
            map_extent=map_extent,
            chunk_stat=chunk_stat,
            chunk_resolution=chunk_resolution,
        )

        # handle harmonisation of axes
        if base_plot_type == "scatter":
            harmonise_xy_lims_paradigm(
                self.read_instance,
                self.canvas_instance,
                relevant_ax,
                base_plot_type,
                plot_characteristics,
                plot_options,
                relim=True,
            )
        elif base_plot_type not in [
            "legend",
            "metadata",
            "map",
            "taylor",
            "fairmode-statsummary",
        ]:
            harmonise_xy_lims_paradigm(
                self.read_instance,
                self.canvas_instance,
                relevant_ax,
                base_plot_type,
                plot_characteristics,
                plot_options,
                relim=True,
                autoscale=True,
            )

        # make legend (embedded on plot axis)
        if (legend) & (base_plot_type != "legend"):
            if "legend" in plot_characteristics:
                # only make map legend in 'domain' plot option is a active
                # also remove observations from legend
                valid_legend = True
                if base_plot_type == "map":
                    if "domain" in plot_options:
                        set_obs_legend = False
                    else:
                        valid_legend = False

                if valid_legend:
                    legend_handles = self.canvas_instance.legend(
                        plot_type, data_labels=legend_data_labels, set_obs=set_obs_legend,
                        map_grid_data_label=labela if base_plot_type == "map" else None,
                    )
                    legend_ax = fig.add_subplot(main_gs[1])
                    legend_ax.axis("off")
                    legend_ax.legend(**legend_handles)

        # make colourbar (embedded on plot axis)
        if "cb" in plot_characteristics:
            cb_ax = self.canvas_instance.colourbar(fig, relevant_ax, zstat, speci, plot_type)
            # hide colourbar if requested, we still need to create it to get the correct cmap / bounds in the maps
            if not cb:
                cb_ax.set_visible(False)

    def make_map(
        self,
        relevant_axis,
        networkspeci,
        plot_characteristics,
        plot_options,
        zstat=None,
        labela="",
        labelb="",
        var=None,
        lat=None,
        lon=None,
        date_range=None,
        map_extent=None,
    ):
        """
        Renders a geospatial scatter plot of stations onto a map axis, coloured by a calculated statistical metric.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstat : str, optional
            Statistic to plot.
        labela : str, optional
            Label of first dataset.
        labelb : str, optional
            Label of second dataset (if defined then a bias plot is made).
        var : np.ndarray, optional
            Gridded model data, by default None (no grid).
        lat : np.ndarray, optional
            Latitudes of the gridded model data, by default None (no grid).
        lon : np.ndarray, optional
            Longitudes of the gridded model data, by default None (no grid).
        date_range : tuple, optional
            Start and end (inclusive) datetimes of period used for station statistics,
            by default None (all loaded period).
        map_extent : array-like, shape (4,), optional
            Extent the map will be shown at, used to size the markers against
            the stations that will actually be on show.
        """

        # get statistic per station (None if only model gridded data is loaded)
        z_statistic, active_map_valid_station_inds = self.get_map_data(
            networkspeci, zstat, labela=labela, labelb=labelb, date_range=date_range
        )
        plot_stations = z_statistic is not None

        # if not only model gridded data is loaded
        if plot_stations:
            # get marker size (for report and library). The size is restored
            # afterwards, as plot_characteristics is shared between every map a
            # report draws - left in place, the first map's size would be reused
            # for all the rest, whatever their own station count and extent
            original_markersize = plot_characteristics["plot"]["stations"]["s"]
            if self.read_instance.mode in ["report", "library"]:
                self.get_markersize(
                    relevant_axis,
                    "map",
                    networkspeci,
                    plot_characteristics,
                    active_map_valid_station_inds=active_map_valid_station_inds,
                    map_extent=map_extent,
                )
            # if using dashboard make z_statistic and active_map_valid_station_inds class variables
            else:
                self.canvas_instance.z_statistic = z_statistic
                self.canvas_instance.active_map_valid_station_inds = (
                    active_map_valid_station_inds
                )

        # plot model gridded data
        plot_grid = var is not None and lon is not None and lat is not None
        if plot_grid:
            grid_mesh = relevant_axis.pcolormesh(
                lon,
                lat,
                var,
                transform=self.canvas_instance.datacrs,
                **plot_characteristics["plot"]["grid"],
            )

            # save grid edges of gridded model (first gridded model is the one read),
            # from mesh cell corners (bottom, right, top, left edges), to plot model domain
            gridded_label = next(
                label
                for label, label_raw in zip(
                    self.read_instance.data_labels, self.read_instance.data_labels_raw
                )
                if label_raw.endswith("::gridded")
            )
            corners = grid_mesh.get_coordinates()
            grid_edges = np.concatenate(
                (
                    corners[0, :],
                    corners[1:, -1],
                    corners[-1, -2::-1],
                    corners[-2:0:-1, 0],
                )
            )
            self.read_instance.plotting_params[gridded_label]["grid_edge_longitude"] = grid_edges[:, 0]
            self.read_instance.plotting_params[gridded_label]["grid_edge_latitude"] = grid_edges[:, 1]

            # allow axis margins around the grid (pcolormesh sticks limits to the grid edges)
            grid_mesh.sticky_edges.x[:] = []
            grid_mesh.sticky_edges.y[:] = []
            relevant_axis.autoscale_view()

        # if not only model gridded data is loaded
        if plot_stations:
            # the stations_over_grid style (an edge, by default) only serves to separate
            # the stations from the gridded field underneath, so apply it only then.
            # copied, as plot_characteristics is shared between every map a report draws
            stations_kwargs = copy.deepcopy(plot_characteristics["plot"]["stations"])
            if plot_grid:
                stations_kwargs.update(
                    plot_characteristics["plot"].get("stations_over_grid", {})
                )

            # plot new station points on map - coloured by currently active z statisitic
            self.stations_scatter = relevant_axis.scatter(
                self.read_instance.station_longitudes[networkspeci][
                    active_map_valid_station_inds
                ],
                self.read_instance.station_latitudes[networkspeci][
                    active_map_valid_station_inds
                ],
                c=z_statistic,
                transform=self.canvas_instance.datacrs,
                **stations_kwargs,
            )
            plot_characteristics["plot"]["stations"]["s"] = original_markersize

            # track plot elements
            if self.read_instance.mode not in ["report"]:
                self.track_plot_elements(
                    self.read_instance.observations_data_label,
                    "map",
                    "plot",
                    [self.stations_scatter],
                    bias=False,
                )

    def make_timeseries(
        self,
        relevant_axis,
        networkspeci,
        data_labels,
        plot_characteristics,
        plot_options,
        chunk_stat=None,
        chunk_resolution=None,
    ):
        """
        Renders temporal data as a line plot, optionally applying statistical chunking or bias calculations.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        chunk_stat : str, optional
            Chunk statistic.
        chunk_resolution : str, optional
            Chunk resolution.
        """

        # create list for timeseries plot
        self.timeseries_plot = []

        # get timeseries per data label, and position of bias line (None if there is none)
        timeseries, bias_line_y, bias = self.get_timeseries_data(
            networkspeci,
            data_labels,
            plot_characteristics,
            plot_options,
            chunk_stat=chunk_stat,
            chunk_resolution=chunk_resolution,
        )

        # plot horizontal line across x axis if bias plot
        if bias_line_y is not None:
            bias_line_chars = {
                k: v for k, v in plot_characteristics["bias_line"].items() if k != "y"
            }
            bias_line = relevant_axis.axhline(y=bias_line_y, **bias_line_chars)
            # track plot elements
            if self.read_instance.mode not in ["report"]:
                self.track_plot_elements(
                    "ALL", "timeseries", "bias_line", [bias_line], bias=bias
                )

        # iterate through data labels
        for data_label, ts in timeseries.items():
            # get marker size (for report and library)
            if self.read_instance.mode in ["report", "library"]:
                self.get_markersize(
                    relevant_axis,
                    "timeseries",
                    networkspeci,
                    plot_characteristics,
                    data=ts,
                )

            # make timeseries plot
            self.timeseries_plot.append(
                relevant_axis.plot(
                    ts,
                    color=self.read_instance.plotting_params[data_label]["colour"],
                    **plot_characteristics["plot"],
                )
            )

            # update maximum smooth value
            if self.read_instance.mode not in ["report", "library"]:
                self.canvas_instance.timeseries_smooth_window_sl.setMaximum(len(ts))
                # the window is shown as the time it spans, which changes with
                # the resolution drawn at whether or not its step count does
                self.canvas_instance.timeseries_smooth_window_sl.refresh_readout()
                # To get straight line
                # if self.canvas_instance.timeseries_smooth_window_sl.value() != (len(ts)*2 - 1):
                #     self.canvas_instance.timeseries_smooth_window_sl.setMaximum(int(len(ts)*2 - 1))

            # track plot elements
            if self.read_instance.mode not in ["report"]:
                self.track_plot_elements(
                    data_label,
                    "timeseries",
                    "plot",
                    self.timeseries_plot[-1],
                    bias=bias,
                )

    def make_fairmode_statsummary(
        self,
        relevant_axis,
        networkspeci,
        data_labels,
        plot_characteristics,
        plot_options,
    ):
        """
        Renders a FAIRMODE statistical summary plot to provide a detailed breakdown of model performance indicators.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """

        # resample to daily for PM10 and PM2.5 if data is hourly
        # get MDA8 for ozone if data is hourly
        # finally filter by coverage
        data, valid_station_idxs = get_fairmode_data(
            self.read_instance, self.canvas_instance, networkspeci, data_labels
        )

        # skip making plot if there is no valid data
        # library and report modes are already handling this in advance
        if (self.read_instance.mode not in ["report", "library"]) and (
            not any(valid_station_idxs)
        ):
            msg = "No valid data to create FAIRMODE statistic summary plot after filtering by coverage."
            show_message(self.read_instance, msg)
            self.read_instance.handle_layout_update(
                "None",
                sender=self.canvas_instance.get_plot_type_position(
                    "fairmode-statsummary"
                ),
            )
            return

        observations_data = data[0, :, :]

        # get settings
        speci = networkspeci.split("|")[1]
        u_95r_RV = fairmode_settings[speci].get("u_95r_RV")
        RV = fairmode_settings[speci].get("RV")
        alpha = fairmode_settings[speci].get("alpha")
        beta = fairmode_settings[speci].get("beta")
        exc_threshold = fairmode_settings[speci].get("exc_threshold")
        percentile = fairmode_settings[speci].get("percentile")
        units = fairmode_settings[speci].get("units")

        # get RV and exceedance threshold per units
        RV, exc_threshold = get_fairmode_RV_exceendance(
            self.read_instance, speci, RV, exc_threshold, units
        )

        # get station references
        valid_station_references = get_valid_metadata(
            self, "station_reference", valid_station_idxs, networkspeci
        )
        valid_station_names = get_valid_metadata(
            self, "station_name", valid_station_idxs, networkspeci
        )
        # same area/station classifications the target plot's legend shows
        valid_area_classifications = None
        valid_station_type_classifications = None
        for other_type in ("area", "station"):
            try:
                fetched = get_valid_metadata(
                    self,
                    f"{other_type}_classification",
                    valid_station_idxs,
                    networkspeci,
                )
            except:
                fetched = np.full(
                    len(valid_station_references), np.nan, dtype=np.float32
                )
            if other_type == "area":
                valid_area_classifications = fetched
            else:
                valid_station_type_classifications = fetched

        # every row below is built from arrays in this same station order -
        # kept here so the hover annotation can name a hovered station
        self.canvas_instance.fairmode_statsummary_stations = (
            valid_station_names,
            valid_station_references,
            valid_area_classifications,
            valid_station_type_classifications,
        )

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[
            networkspeci
        ]

        # cut data_labels for those in valid data labels
        cut_data_labels = [
            data_label for data_label in data_labels if data_label in valid_data_labels
        ]

        self.fairmode_statsummary_row_titles = {}
        # per (data_label, row) station indices behind the middle-zone dots -
        # the collapsed left/right corner dots stand for more than one
        # station, so no single index is kept for those
        self.canvas_instance.fairmode_statsummary_row_station_indices = {}

        # iterate through data labels
        for data_label in cut_data_labels:
            self.fairmode_statsummary_row_titles[data_label] = []

            # get model data
            model_data = data[valid_data_labels.index(data_label), :, :]

            # calculate MQI for the current station
            exceedances = []
            means = []
            t_biases = []
            t_R_list = []
            t_sd_list = []
            h_perc_list = []

            for station_idx, station in enumerate(valid_station_references):
                st_observations_data = observations_data[station_idx, :]
                st_model_data = model_data[station_idx, :]

                mean, exc, t_bias, t_R, t_sd, h_perc = ModBias.calculate_fairmode_stats(
                    st_observations_data,
                    st_model_data,
                    u_95r_RV,
                    RV,
                    alpha,
                    beta,
                    exc_threshold,
                    percentile,
                    "summary",
                )

                means.append(mean)
                exceedances.append(exc)
                t_biases.append(t_bias)
                t_R_list.append(t_R)
                t_sd_list.append(t_sd)
                h_perc_list.append(h_perc)

            # plot data
            # join all statistics in one list
            statistics_list = [
                np.array(means),
                np.array(exceedances),
                np.array(t_biases),
                np.array(t_R_list),
                np.array(t_sd_list),
                np.array(h_perc_list),
                np.nanmean(t_sd_list),
                np.nanmean(t_R_list),
            ]

            # get subplots dictionary
            subplots = dict(plot_characteristics["auxiliar"]["subplots"])

            # get the variable that tells you if there are exceedances in this species
            has_exceedances = exc_threshold is not None

            # if there is no threshold don't create the exceedances row
            if not has_exceedances:
                subplots.pop("observed exceedances", None)
                statistics_list.pop(1)

            # create list for track_plot_elements
            self.fairmode_statsummary_plot = []

            # apply configuration to each row
            for i, (row, fairmode_data) in enumerate(zip(subplots, statistics_list)):
                # for the two firsts rows, skip the models
                if i < 2 and data_label != self.read_instance.observations_data_label:
                    continue

                # skip the observations all the rows but the two first ones
                if i >= 2 and data_label == self.read_instance.observations_data_label:
                    continue

                # get row dictionary
                plot_dict = subplots[row]

                # station indices behind fairmode_data, whittled down in
                # step with it below (None for the two scalar summary rows,
                # which are not per-station)
                if isinstance(fairmode_data, np.ndarray):
                    station_indices = np.arange(len(valid_station_references))
                else:
                    station_indices = None

                # remove axis from the dot on right side
                relevant_axis[i * 4 + 3].set_xticks([])

                # remove the axis of the dot
                for side in ["bottom", "top", "right"]:
                    relevant_axis[i * 4 + 3].spines[side].set_linestyle("none")

                # add dashed line on the left
                relevant_axis[i * 4 + 3].spines["left"].set_linestyle((10, (8, 5)))

                # configure color of the row and the dot on the right for rows 3 to 8
                if "range_style" in plot_dict:
                    # get the dictionary with the style of the row and its range
                    range_style_dict = plot_characteristics["auxiliar"]["range_style"][
                        plot_dict["range_style"]
                    ]
                    for span, color in zip(
                        range_style_dict["spans"], range_style_dict["colors"]
                    ):
                        relevant_axis[i * 4 + 1].axvspan(*span, color=color, lw=0)

                    # dot on the right configuration
                    # get the lowest and highest number on the range
                    min_span = range_style_dict["spans"][0][0]
                    max_span = range_style_dict["spans"][-1][-1]

                    # get the color of the dot on the right
                    arr = np.array(fairmode_data)[~np.isnan(fairmode_data)]
                    correct_arr = arr[(arr >= min_span) & (arr <= max_span)]

                    dot_color = plot_characteristics["auxiliar"]["right_dot_colors"][
                        "red"
                    ]
                    if len(arr) > 0 and len(correct_arr) / len(arr) >= 0.9:
                        dot_color = plot_characteristics["auxiliar"][
                            "right_dot_colors"
                        ]["green"]

                    # plot dot on the right
                    relevant_axis[i * 4 + 3].scatter(
                        **plot_characteristics["auxiliar"]["right_dot"],
                        color=dot_color,
                        edgecolor=dot_color,
                    )

                # y axis / grid
                # remove y axis ticks
                for j in range(4):
                    relevant_axis[i * 4 + j].set_yticks([])
                    relevant_axis[i * 4 + j].grid(False)

                # x axis
                # get the x axis limit for the current row
                x_limit = plot_dict["x_axis_limits"]

                # set x axis limits
                relevant_axis[i * 4 + 1].set_xlim(*x_limit)

                # right zone configuration
                # remove x ticks from the right dashed zone
                relevant_axis[i * 4 + 2].set_xticks([])

                # remove vertical lines separating middle and right dashed zone
                relevant_axis[i * 4 + 2].spines["left"].set_linestyle("none")
                relevant_axis[i * 4 + 1].spines["right"].set_linestyle("none")

                # change the linestyle to dashed
                for side in ["bottom", "top", "right"]:
                    relevant_axis[i * 4 + 2].spines[side].set_linestyle((10, (8, 5)))

                # right zone dot configuration
                # get the points that surpass the limit
                right_zone_mask = x_limit[1] < fairmode_data

                # set the right dashed zone range to (-1,1)
                relevant_axis[i * 4 + 2].set_xlim(-1, 1)

                # if there is a dot outside the limits
                if np.any(right_zone_mask):
                    # plot it in the middle of the right dashed zone
                    stations_dots = relevant_axis[i * 4 + 2].plot(
                        0,
                        0,
                        plot_characteristics["auxiliar"]["station_dots"]["marker"],
                        color=self.read_instance.plotting_params[data_label]["colour"],
                        markersize=plot_characteristics["auxiliar"]["station_dots"][
                            "markersize"
                        ],
                    )

                    # add the dots to the track plot elements list
                    self.fairmode_statsummary_plot.append(stations_dots[0])
                    self.fairmode_statsummary_row_titles[data_label].append(
                        f"{row}_right"
                    )

                    # remove it from the data plotted in the middle zone
                    fairmode_data = (
                        fairmode_data[~right_zone_mask]
                        if isinstance(fairmode_data, np.ndarray)
                        else np.array([])
                    )
                    if station_indices is not None:
                        station_indices = station_indices[~right_zone_mask]

                # left zone configuration
                # remove x ticks from the left dashed zone
                relevant_axis[i * 4 + 0].set_xticks([])

                # define left zone line style (can be dashed or no left zone)
                left_dashed_zone_linestyle = plot_dict["left_dashed_zone_linestyle"]

                # if left dashed zone exists in the current row
                if left_dashed_zone_linestyle != "none":
                    # convert to tuple [x,[x,x]] because yaml does not return python tuples
                    left_dashed_zone_linestyle = (
                        left_dashed_zone_linestyle[0],
                        tuple(left_dashed_zone_linestyle[1]),
                    )

                    # remove vertical lines separating middle and left dashed zone
                    relevant_axis[i * 4 + 1].spines["left"].set_linestyle("none")
                    relevant_axis[i * 4 + 0].spines["right"].set_linestyle("none")

                # change the linestyle to dashed or remove the dashed zone
                for side in ["bottom", "top", "left"]:
                    relevant_axis[i * 4 + 0].spines[side].set_linestyle(
                        left_dashed_zone_linestyle
                    )

                # left zone dot configuration
                # if left dashed zone exists in the current row
                if left_dashed_zone_linestyle != "none":
                    # get the points that surpass the limit
                    left_zone_mask = x_limit[0] > fairmode_data

                    # set the left dashed zone range to (-1,1)
                    relevant_axis[i * 4 + 0].set_xlim(-1, 1)

                    # if there is a dot outside the limits
                    if np.any(left_zone_mask):
                        # plot it in the middle of the left dashed zone
                        stations_dots = relevant_axis[i * 4 + 0].plot(
                            0,
                            0,
                            plot_characteristics["auxiliar"]["station_dots"]["marker"],
                            color=self.read_instance.plotting_params[data_label][
                                "colour"
                            ],
                            markersize=plot_characteristics["auxiliar"]["station_dots"][
                                "markersize"
                            ],
                        )

                        # add the dot to the track plot elements list
                        self.fairmode_statsummary_plot.append(stations_dots[0])
                        self.fairmode_statsummary_row_titles[data_label].append(
                            f"{row}_left"
                        )

                        # remove it from the data plotted in the middle zone
                        fairmode_data = (
                            fairmode_data[~left_zone_mask]
                            if isinstance(fairmode_data, np.ndarray)
                            else np.array([])
                        )
                        if station_indices is not None:
                            station_indices = station_indices[~left_zone_mask]

                # plot stations as dots
                stations_dots = relevant_axis[i * 4 + 1].plot(
                    fairmode_data,
                    np.zeros_like(fairmode_data),
                    plot_characteristics["auxiliar"]["station_dots"]["marker"],
                    color=self.read_instance.plotting_params[data_label]["colour"],
                    markersize=plot_characteristics["auxiliar"]["station_dots"][
                        "markersize"
                    ],
                )

                # change the size of the x-axis tick labels
                relevant_axis[i * 4 + 1].tick_params(
                    **plot_characteristics["auxiliar"]["station_dots"]["tick_params"]
                )

                # add the dots to the track plot elements list
                self.fairmode_statsummary_plot.append(stations_dots[0])
                self.fairmode_statsummary_row_titles[data_label].append(f"{row}")
                if station_indices is not None:
                    self.canvas_instance.fairmode_statsummary_row_station_indices[
                        (data_label, row)
                    ] = station_indices

            # track plot elements
            if self.read_instance.mode not in ["report"]:
                self.track_plot_elements(
                    data_label,
                    "fairmode-statsummary",
                    "plot",
                    self.fairmode_statsummary_plot,
                    bias=False,
                )

        # plot the row titles, descriptions, separators and units
        for i, (row, fairmode_data) in enumerate(zip(subplots, statistics_list)):
            # get row dictionary
            plot_dict = subplots[row]

            # get the row title
            row_title = plot_dict["title"]

            # write the threshold on the exceedances row title
            if row == "observed_exceedances":
                exc_threshold = round(
                    exc_threshold, plot_characteristics["round_decimal_places"]["title"]
                )
                row_title = row_title.format(exc_threshold)

            # add units to the first two rows
            if "units" in plot_dict:
                relevant_axis[i * 4 + 3].text(
                    *plot_characteristics["auxiliar"]["units"]["position"],
                    plot_dict["units"],
                    fontsize=plot_characteristics["auxiliar"]["units"]["fontsize"],
                )

            # add row title
            relevant_axis[i * 4 + 0].text(
                *plot_characteristics["auxiliar"]["row_title"]["position"],
                row_title,
                **plot_characteristics["auxiliar"]["row_title"],
                transform=relevant_axis[i * 4 + 0].transAxes,
            )

            # add left description
            relevant_axis[i * 4 + 0].text(
                *plot_characteristics["auxiliar"]["description"]["position"],
                plot_dict["description"],
                **plot_characteristics["auxiliar"]["description"],
                transform=relevant_axis[i * 4 + 0].transAxes,
            )

            # add separator
            if "separator" in plot_dict:
                relevant_axis[i * 4 + 0].text(
                    *plot_characteristics["auxiliar"]["separator"]["position"],
                    plot_characteristics["auxiliar"]["separator_text"],
                    **plot_characteristics["auxiliar"]["separator"],
                    transform=relevant_axis[i * 4 + 0].transAxes,
                )

        # add title if using dashboard
        if self.read_instance.mode not in ["report", "library"]:
            set_axis_title(
                self.read_instance,
                relevant_axis,
                fairmode_settings[speci]["title"],
                plot_characteristics,
            )
