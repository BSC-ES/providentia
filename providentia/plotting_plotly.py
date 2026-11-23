import cartopy.io.shapereader as shapereader
import matplotlib
from matplotlib.colors import to_hex
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from scipy.spatial import cKDTree
import yaml

from providentia.auxiliar import PROVIDENTIA_ROOT, join, get_map_colours
from .plot_aux import get_land_polygon_resolution
from .plot_options_plotly import format_plot_options
from .plotting import Plotting
from .statistics import generate_colourbar_detail



class PlottingPlotly(Plotting):
    """Class for plotting with plotly: Plotly Dash (web-based) and library"""

    def __init__(self, read_instance=None, canvas_instance=None):
        """
        Initialise the PlottingPlotly instance.

        Parameters
        ----------
        read_instance : object, optional
            The data management instance containing configuration and raw data.
        canvas_instance : object, optional
            The instance that holds the plot characteristics and the selected data.
        """

        super(PlottingPlotly, self).__init__(
            read_instance=read_instance, canvas_instance=canvas_instance
        )

        # land, coastlines and borders of the map in projected coordinates, kept for each
        # projection and resolution as projecting them takes about a second and they do not
        # change with the data
        self.projected_map_lines = {}

    def create_figure(
        self, plot_type, base_plot_type, networkspeci, width=None, height=None, legend=True
    ):
        """
        Creates the plotly figure of a library plot, which has no axes to pass around.

        Parameters
        ----------
        plot_type : str
            Plot type, with its options (not used, plotly figures are all created the same).
        base_plot_type : str
            Plot type without options (not used).
        networkspeci : str
            Current networkspeci (not used).
        width : float, optional
            Not used, the figure takes the size of where it is shown.
        height : float, optional
            Not used, the figure takes the size of where it is shown.
        legend : bool, optional
            Indicates if the legend is shown, by default True.

        Returns
        -------
        fig : plotly.graph_objects.Figure
            Figure.
        relevant_ax : None
            Plotly figures have no axes.
        main_gs : None
            Plotly figures have no grid for the legend, which plotly draws itself.
        """

        # plotly draws its own legend from the trace names, unless told not to
        fig = go.Figure(layout={"showlegend": legend})

        return fig, None, None

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
        Sets the title and the axis titles of a library plot, and applies its plot options
        (see plot_options_plotly.py). The other arguments are those that the matplotlib class
        needs to format its axes and add its legend and colourbar, which plotly figures do not
        have or draw themselves.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure of the plot.
        relevant_ax : object
            Axis of the plot (not used, None for plotly).
        main_gs : matplotlib.gridspec.GridSpec
            Grid of the figure, whose second row is for the legend (not used, None for plotly).
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

        # plotly figures have no matplotlib axes, so the title is set on the figure layout:
        # small and close to the plot, and centred over the plot as in matplotlib, not at the
        # left of the figure
        if title:
            fig.update_layout(
                title_text=title.replace("\n", "<br>"),
                title_font_size=plot_characteristics["axis_title"]["fontsize"],
                title_x=0.5,
                title_xanchor="center",
                title_xref="paper",
            )

        # axis titles, with the font sizes of the website block
        if xlabel:
            fig.update_xaxes(
                title_text=xlabel,
                title_font_size=plot_characteristics["xlabel"]["fontsize"],
            )
        if ylabel:
            fig.update_yaxes(
                title_text=ylabel,
                title_font_size=plot_characteristics["ylabel"]["fontsize"],
            )

        # apply the plot options (log axes, ...)
        format_plot_options(fig, plot_options)

    def make_map(
        self,
        fig,
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
        Renders the map on a plotly figure: cartopy projects everything into the map
        projection, and plotly draws it on plain x/y axes. The gridded model is a heatmap,
        resampled onto a regular raster in projected coordinates (any grid shape works, and
        its size does not depend on the model resolution), the stations are points on top,
        coloured by a calculated statistical metric, and the coastlines and borders are lines.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot (not used yet in plotly).
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
            Extent of the map (not used yet in plotly, which shows all the data).
        """

        # get statistic per station (None if only model gridded data is loaded)
        z_statistic, active_map_valid_station_inds = self.get_map_data(
            networkspeci, zstat, labela=labela, labelb=labelb, date_range=date_range
        )

        projection = self.canvas_instance.plotcrs
        datacrs = self.canvas_instance.datacrs
        speci = networkspeci.split("|")[1]
        map_template = self.canvas_instance.plot_characteristics_templates["map"]
        # decimals of the values on hover: only the dashboard has them in its characteristics
        # (2), so the same number is used when they are not given
        rounding = plot_characteristics.get("marker_annotate_rounding", 2)
        plot_grid = var is not None and lon is not None and lat is not None

        # points of the data in projected coordinates, to set the limits of the map
        extent_points = []
        # values plotted, to set the limits of the colourbar
        plotted_values = []

        # gridded model, under the stations
        cell_tree = None
        if plot_grid:
            if np.ndim(lon) == 1:
                lon, lat = np.meshgrid(lon, lat)
            cells = projection.transform_points(datacrs, lon, lat)[..., :2].reshape(-1, 2)
            cell_values = np.ma.filled(
                np.ma.masked_invalid(var).astype(np.float32), np.nan
            ).ravel()
            # cells that the projection does not show have no place in it
            shown = np.isfinite(cells).all(axis=1)
            cells, cell_values = cells[shown], cell_values[shown]

            # regular raster in projected coordinates, each pixel with the value of its
            # nearest cell, so that the values on hover are the ones of the model
            cell_tree = cKDTree(cells)
            raster_shape = plot_characteristics["raster_shape"]
            raster_x = np.linspace(cells[:, 0].min(), cells[:, 0].max(), raster_shape[1])
            raster_y = np.linspace(cells[:, 1].min(), cells[:, 1].max(), raster_shape[0])
            pixels_x, pixels_y = np.meshgrid(raster_x, raster_y)
            distance, nearest = cell_tree.query(
                np.column_stack([pixels_x.ravel(), pixels_y.ravel()]), workers=-1
            )
            # typical distance between cells, to leave empty the pixels outside the grid
            sample = cells[:: max(1, len(cells) // 2000)]
            cell_size = np.median(cell_tree.query(sample, k=2)[0][:, 1])
            raster = np.where(distance <= cell_size, cell_values[nearest], np.nan)
            raster = raster.reshape(raster_shape).astype(np.float32)
            # longitude and latitude of each pixel, for the hover box
            pixels_lonlat = datacrs.transform_points(projection, pixels_x, pixels_y)[..., :2]

            fig.add_trace(
                go.Heatmap(
                    x=raster_x,
                    y=raster_y,
                    z=raster,
                    customdata=pixels_lonlat.astype(np.float32),
                    coloraxis="coloraxis",
                    opacity=plot_characteristics["plot"]["grid"]["alpha"],
                    hovertemplate=(
                        "Grid: %{{z:.{}f}}<br>Longitude: %{{customdata[0]:.2f}}"
                        "<br>Latitude: %{{customdata[1]:.2f}}<extra></extra>"
                    ).format(rounding),
                )
            )
            extent_points.append(cells)
            plotted_values.append(cell_values)

        resolution = get_land_polygon_resolution(map_template["map_resolution"])

        # ocean and land in the colours of the map (those of its colour preset, if they are
        # not set), under the grid and the stations as in the dashboard. They are shapes of
        # the figure, which plotly draws under all the data
        land_colour, ocean_colour = get_map_colours(map_template)
        # edge of the world in the projection
        boundary_x, boundary_y = projection.boundary.xy
        shapes = []
        if map_template.get("ocean_polygon", {}).get("visible", True):
            shapes.append(
                {
                    "path": "M" + "L".join(
                        "{:.0f},{:.0f}".format(x, y) for x, y in zip(boundary_x, boundary_y)
                    ) + "Z",
                    "fillcolor": ocean_colour,
                }
            )
        # projected once for each projection, as the land does not change with the data
        key = (projection.proj4_init, resolution, "land")
        if key not in self.projected_map_lines:
            land_path = ""
            # the files cannot be downloaded on machines without internet
            try:
                path = shapereader.natural_earth(
                    resolution=resolution, category="physical", name="land"
                )
                geometries = list(shapereader.Reader(path).geometries())
            except Exception as error:
                self.read_instance.logger.warning(
                    "Map land ({}) not available, not drawn: {}".format(resolution, error)
                )
                geometries = []
            for geometry in geometries:
                projected = projection.project_geometry(geometry, datacrs).simplify(
                    plot_characteristics["land_tolerance"]
                )
                for polygon in getattr(projected, "geoms", [projected]):
                    if polygon.is_empty:
                        continue
                    # the outline of the land and then its holes (e.g. lakes)
                    for ring in [polygon.exterior] + list(polygon.interiors):
                        land_path += "M" + "L".join(
                            "{:.0f},{:.0f}".format(x, y) for x, y in ring.coords
                        ) + "Z"
            self.projected_map_lines[key] = land_path
        if self.projected_map_lines[key]:
            shapes.append({"path": self.projected_map_lines[key], "fillcolor": land_colour})
        for shape in shapes:
            shape.update(
                {"type": "path", "xref": "x", "yref": "y", "layer": "below", "line": {"width": 0}}
            )

        # coastlines, and borders if they are shown, at the resolution of the map
        map_lines = [("physical", "coastline", 0.8)]
        if map_template.get("borders", {}).get("visible", False):
            map_lines.append(("cultural", "admin_0_boundary_lines_land", 0.4))
        for category, name, width in map_lines:
            # projected once for each projection, as they do not change with the data
            key = (projection.proj4_init, resolution, name)
            if key not in self.projected_map_lines:
                line_x, line_y = [], []
                # the files cannot be downloaded on machines without internet
                try:
                    path = shapereader.natural_earth(
                        resolution=resolution, category=category, name=name
                    )
                    geometries = list(shapereader.Reader(path).geometries())
                except Exception as error:
                    self.read_instance.logger.warning(
                        "Map {} ({}) not available, not drawn: {}".format(name, resolution, error)
                    )
                    geometries = []
                for geometry in geometries:
                    projected = projection.project_geometry(geometry, datacrs)
                    for line in getattr(projected, "geoms", [projected]):
                        # a missing value separates lines in plotly
                        line_x.extend(list(line.xy[0]) + [np.nan])
                        line_y.extend(list(line.xy[1]) + [np.nan])
                # as arrays, which plotly writes in much less space than lists
                self.projected_map_lines[key] = (
                    np.array(line_x, dtype=np.float32),
                    np.array(line_y, dtype=np.float32),
                )
            line_x, line_y = self.projected_map_lines[key]
            fig.add_trace(
                go.Scattergl(
                    x=line_x,
                    y=line_y,
                    mode="lines",
                    hoverinfo="skip",
                    line={"color": "black", "width": width},
                )
            )

        # edge of the world, over the grid as the coastlines
        fig.add_trace(
            go.Scattergl(
                x=np.array(boundary_x, dtype=np.float32),
                y=np.array(boundary_y, dtype=np.float32),
                mode="lines",
                hoverinfo="skip",
                line={"color": "black", "width": 0.8},
            )
        )

        # stations, coloured by the statistic
        if z_statistic is not None:
            station_lon = self.read_instance.station_longitudes[networkspeci][
                active_map_valid_station_inds
            ]
            station_lat = self.read_instance.station_latitudes[networkspeci][
                active_map_valid_station_inds
            ]
            station_names = self.read_instance.station_names[networkspeci][
                active_map_valid_station_inds
            ]
            station_references = self.read_instance.station_references[networkspeci][
                active_map_valid_station_inds
            ]
            stations = projection.transform_points(datacrs, station_lon, station_lat)[:, :2]
            # stations that the projection does not show have no place in it
            shown = np.isfinite(stations).all(axis=1)

            # value of the grid under each station, added to its hover box
            station_grid = np.full(len(stations), np.nan)
            if cell_tree is not None:
                distance, nearest = cell_tree.query(stations[shown])
                station_grid[shown] = np.where(
                    distance <= cell_size, cell_values[nearest], np.nan
                )

            # hover box with the lines of the dashboard, and the grid value if there is one
            hover_texts = []
            for station_ii in np.flatnonzero(shown):
                text = "Station: {}<br>Reference: {}<br>Longitude: {:.2f}<br>Latitude: {:.2f}".format(
                    station_names[station_ii],
                    station_references[station_ii],
                    station_lon[station_ii],
                    station_lat[station_ii],
                )
                text += "<br>{}: {:.{}f}".format(zstat, z_statistic[station_ii], rounding)
                if np.isfinite(station_grid[station_ii]):
                    text += "<br>Grid: {:.{}f}".format(station_grid[station_ii], rounding)
                hover_texts.append(text)

            # an edge separates the stations from the gridded field underneath
            marker_line = {"width": 0}
            if plot_grid:
                stations_over_grid = plot_characteristics["plot"].get("stations_over_grid", {})
                marker_line = {
                    "color": stations_over_grid.get("edgecolor", "black"),
                    "width": 2 * stations_over_grid.get("linewidth", 0.4),
                }
            fig.add_trace(
                go.Scattergl(
                    x=stations[shown, 0],
                    y=stations[shown, 1],
                    mode="markers",
                    marker={
                        "color": np.asarray(z_statistic, dtype=np.float64)[shown],
                        "coloraxis": "coloraxis",
                        "size": plot_characteristics["marker_size"],
                        "line": marker_line,
                    },
                    hovertext=hover_texts,
                    hoverinfo="text",
                )
            )
            extent_points.append(stations[shown])
            plotted_values.append(np.asarray(z_statistic, dtype=np.float64)[shown])

        # axes without ticks and with the same scale, limited to the data with a margin
        layout = {
            "plot_bgcolor": "white",
            "shapes": shapes,
            # hover boxes with black text on white, as the tooltips of the dashboard, which do
            # not take the colour of the station or of the grid cell under the mouse
            "hoverlabel": {
                "bgcolor": "white",
                "bordercolor": "#757575",
                "font": {"color": "black"},
            },
            # room under the map, in pixels, for the colourbar and white space above it
            "margin": {**plot_characteristics["margin"], "b": plot_characteristics["colourbar_height"]},
            "xaxis": {"visible": False},
            "yaxis": {"visible": False, "scaleanchor": "x"},
        }
        if extent_points:
            extent_points = np.concatenate(extent_points)
            for axis, axis_ii in [("xaxis", 0), ("yaxis", 1)]:
                low, high = extent_points[:, axis_ii].min(), extent_points[:, axis_ii].max()
                margin = 0.05 * (high - low)
                layout[axis]["range"] = [low - margin, high + margin]

        # horizontal colourbar under the map, with the limits, colourmap and label that the
        # statistic has in the other backends
        plotted_values = np.concatenate(plotted_values) if plotted_values else np.array([np.nan])
        if np.isfinite(plotted_values).any():
            z_vmin, z_vmax, z_label, z_colourmap, n_discrete, n_ticks = generate_colourbar_detail(
                self.read_instance,
                zstat,
                np.nanmin(plotted_values),
                np.nanmax(plotted_values),
                plot_characteristics,
                speci,
            )
            color_palettes = yaml.safe_load(
                open(join(PROVIDENTIA_ROOT, "settings/color_palettes.yaml"))
            )
            if z_colourmap in color_palettes:
                cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                    "providentia", color_palettes[z_colourmap], n_discrete
                )
            else:
                cmap = plt.get_cmap(z_colourmap, n_discrete)
            # one step for each colour of a discrete colourmap
            colorscale = []
            for colour_ii in range(cmap.N):
                colour = to_hex(cmap(colour_ii))
                colorscale.append([colour_ii / cmap.N, colour])
                colorscale.append([(colour_ii + 1) / cmap.N, colour])
            layout["coloraxis"] = {
                "cmin": float(z_vmin),
                "cmax": float(z_vmax),
                "colorscale": colorscale,
                "colorbar": {
                    "orientation": "h",
                    "x": 0.5,
                    "xanchor": "center",
                    # at the bottom of the figure, apart from the map, nearly as wide as it
                    # and with the thickness it has in the dashboard
                    "yref": "container",
                    "y": 0,
                    "yanchor": "bottom",
                    "ypad": 2,
                    "len": 0.9,
                    "thickness": 14,
                    "tickvals": np.linspace(z_vmin, z_vmax, n_ticks, endpoint=True).tolist(),
                    "tickfont": {"size": 10},
                    "title": {"text": z_label, "side": "bottom", "font": {"size": 10}},
                },
            }
        fig.update_layout(**layout)

    def add_bias_line(self, fig, bias_line, y):
        """
        Draws the horizontal line of a bias plot across the x axis.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure to plot on.
        bias_line : dict
            Plot characteristics of the bias line, with matplotlib's line styles.
        y : float
            Position of the line in the y axis.
        """

        # the plot characteristics have matplotlib's line styles, which plotly names
        line_dashes = {"-": "solid", "--": "dash", ":": "dot", "-.": "dashdot"}
        fig.add_hline(
            y=y,
            line_dash=line_dashes[bias_line["linestyle"]],
            line_color=bias_line["color"],
            line_width=bias_line["linewidth"],
        )

    def make_timeseries(
        self,
        fig,
        networkspeci,
        data_labels,
        plot_characteristics,
        plot_options,
        chunk_stat=None,
        chunk_resolution=None,
    ):
        """
        Renders temporal data as a plotly scatter, optionally applying statistical chunking or bias calculations.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure to plot on.
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

        # get timeseries per data label, and position of bias line (None if there is none)
        timeseries, bias_line_y, bias = self.get_timeseries_data(
            networkspeci,
            data_labels,
            plot_characteristics,
            plot_options,
            chunk_stat=chunk_stat,
            chunk_resolution=chunk_resolution,
        )

        # plotly leaves 80px on the left, right and bottom and 100px on top by default: only the
        # room that the tick labels and the title need is kept, so that the plot takes more of
        # the figure
        fig.update_layout(margin=plot_characteristics["margin"])

        # plot horizontal line across x axis if bias plot
        if bias_line_y is not None:
            self.add_bias_line(fig, plot_characteristics["bias_line"], bias_line_y)

        # iterate through data labels
        for data_label, ts in timeseries.items():
            # drop missing values, otherwise plotly fits the time axis to the whole period
            # (matplotlib skips them when autoscaling)
            ts_valid = ts.dropna()
            fig.add_trace(
                go.Scatter(
                    x=ts_valid.index,
                    y=ts_valid.values,
                    mode="markers",
                    name=data_label,
                    # model colours are RGB tuples (0-1), which plotly does not accept
                    marker_color=to_hex(self.read_instance.plotting_params[data_label]["colour"]),
                )
            )

    def make_distribution(
        self,
        fig,
        networkspeci,
        data_labels,
        plot_characteristics,
        plot_options,
        data_range_min=None,
        data_range_max=None,
        zstat=None,
        station_statistic_range=None,
    ):
        """
        Renders the probability density functions of get_distribution_data as plotly lines.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        data_range_min : float, optional
            Minimum data range of distribution plot grid.
        data_range_max : float, optional
            Maximum data range of distribution plot grid.
        zstat : str, optional
            Statistic parsed from a "distribution-<stat>" plot type, see get_distribution_data.
        station_statistic_range : tuple of float, optional
            Minimum and maximum of the statistic's own per-station values, see get_distribution_data.
        """

        # get PDF per data label (None if the plot is not made)
        distribution_data = self.get_distribution_data(
            networkspeci,
            data_labels,
            plot_characteristics,
            plot_options,
            data_range_min=data_range_min,
            data_range_max=data_range_max,
            zstat=zstat,
            station_statistic_range=station_statistic_range,
        )
        if distribution_data is None:
            return
        x_grid, PDFs_sampled, bias = distribution_data

        # keep only the room that the tick labels and the titles need (see make_timeseries)
        fig.update_layout(margin=plot_characteristics["margin"])

        # plot horizontal line across x axis at 0 if bias plot
        if bias:
            self.add_bias_line(
                fig, plot_characteristics["bias_line"], plot_characteristics["bias_line"]["y"]
            )

        # iterate through data labels, with every point of the PDF as in matplotlib
        for data_label, PDF_sampled in PDFs_sampled.items():
            fig.add_trace(
                go.Scatter(
                    x=x_grid,
                    y=PDF_sampled,
                    mode="lines",
                    name=data_label,
                    # model colours are RGB tuples (0-1), which plotly does not accept
                    line_color=to_hex(self.read_instance.plotting_params[data_label]["colour"]),
                    line_width=plot_characteristics["plot"]["linewidth"],
                )
            )

    def make_table(
        self,
        fig,
        networkspeci,
        data_labels,
        plot_characteristics,
        plot_options,
        zstats=None,
        statsummary=False,
        subsection=None,
        plotting_paradigm=None,
        stats_df=None,
    ):
        """
        Renders a table of statistical metrics as a plotly table, with the
        cells or column labels coloured as the data labels.

        Parameters
        ----------
        fig : plotly.graph_objects.Figure
            Figure to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstats : list, optional
            Statistics.
        statsummary : boolean, optional
            To indicate if making alternative statistical summary table plot.
        subsection : str, optional
            Currently active subsection.
        plotting_paradigm : str, optional
            Plotting paradigm (summary or station report).
        stats_df : pandas dataframe, optional
            Dataframe of previously calculated statistics.
        """

        (
            stats_df,
            col_labels,
            networkspecies,
            subsections,
            cell_colours,
            col_colours,
            bias,
        ) = self.get_table_data(
            networkspeci,
            data_labels,
            plot_characteristics,
            plot_options,
            zstats=zstats,
            statsummary=statsummary,
            subsection=subsection,
            plotting_paradigm=plotting_paradigm,
            stats_df=stats_df,
        )

        # plotly tables are filled column by column, and show values as matplotlib does
        # when written as text
        columns = [[str(value) for value in stats_df[column]] for column in stats_df.columns]

        # plotly cannot merge cells: in the networkspecies and subsections columns (if any),
        # a value repeated in the rows below is left blank instead
        for column, rows in zip(
            ["networkspecies", "subsections"], (networkspecies, subsections)
        ):
            if column in stats_df.columns:
                column_ii = list(stats_df.columns).index(column)
                columns[column_ii] = [
                    "" if (row_ii > 0) and (row == rows[row_ii - 1]) else value
                    for row_ii, (row, value) in enumerate(zip(rows, columns[column_ii]))
                ]

        # model colours are RGB tuples (0-1), which plotly does not accept
        if cell_colours is not None:
            cells_fill = [
                [to_hex(colour) for colour in column] for column in cell_colours.T
            ]
        else:
            cells_fill = "white"
        if col_colours is not None:
            header_fill = [to_hex(colour) for colour in col_colours]
        else:
            header_fill = "white"

        # plotly leaves 80px on the left, right and bottom and 100px on top by default: only the
        # room that the title needs is kept, so that the table takes more of the figure
        fig.update_layout(margin=plot_characteristics["margin"])

        # plotly gives every column the same width: each one takes the share of the longest
        # text in it, header included, so that the labels of the first column fit
        column_widths = [
            max(len(str(text)) for text in [label, *column])
            for label, column in zip(col_labels, columns)
        ]

        fig.add_trace(
            go.Table(
                columnwidth=column_widths,
                header=dict(
                    values=col_labels,
                    fill_color=header_fill,
                    line_color="black",
                    font_size=plot_characteristics["fontsize"],
                    align="center",
                ),
                cells=dict(
                    values=columns,
                    fill_color=cells_fill,
                    line_color="black",
                    font_size=plot_characteristics["fontsize"],
                    align="center",
                ),
            )
        )
