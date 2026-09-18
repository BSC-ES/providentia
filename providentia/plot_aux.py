""" Auxiliary plotting functions """

import math
import os
import sys

import matplotlib
from matplotlib import transforms
from matplotlib.colors import cnames
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from scipy.signal import convolve
from scipy.signal.windows import gaussian
from scipy.sparse import coo_matrix
import seaborn as sns
import yaml

from providentia.auxiliar import (
    CURRENT_PATH,
    join,
    get_conversion_factor,
    get_standard_parameters_by_speci,
)
from .dashboard_elements import CheckDialog, MessageBox
from .statistics import (
    boxplot_inner_fences,
    calculate_statistic,
    get_z_statistic_info,
    get_z_statistic_sign,
    get_z_statistic_type,
)
from .warnings_prv import show_message

PROVIDENTIA_ROOT = "/".join(CURRENT_PATH.split("/")[:-1])


def get_AERONET_sizedist_bin_radius(species):
    """
    Get AERONET size distribution bin radius for a given species.

    Parameters
    ----------
    species : str
        Name of the species, e.g. 'vconcaerobin1'.

    Returns
    -------
    str
        The radius corresponding to the given species as a string.
    """

    radius_per_species = {
        "vconcaerobin1": "0.05",
        "vconcaerobin2": "0.066",
        "vconcaerobin3": "0.086",
        "vconcaerobin4": "0.113",
        "vconcaerobin5": "0.148",
        "vconcaerobin6": "0.194",
        "vconcaerobin7": "0.255",
        "vconcaerobin8": "0.335",
        "vconcaerobin9": "0.439",
        "vconcaerobin10": "0.576",
        "vconcaerobin11": "0.756",
        "vconcaerobin12": "0.992",
        "vconcaerobin13": "1.302",
        "vconcaerobin14": "1.708",
        "vconcaerobin15": "2.241",
        "vconcaerobin16": "2.940",
        "vconcaerobin17": "3.857",
        "vconcaerobin18": "5.061",
        "vconcaerobin19": "6.641",
        "vconcaerobin20": "8.713",
        "vconcaerobin21": "11.432",
        "vconcaerobin22": "15.00",
    }

    return radius_per_species[species]


def get_multispecies_aliases(networkspecies):
    """
    Map networkspecies to networkspecies aliases. Also get label for alias.

    Parameters
    ----------
    networkspecies : list of str
        List of networkspecies names to map.

    Returns
    -------
    networkspecies_aliases : list of str
        List of aliases corresponding to each networkspecies.
    unique_label : str
        Unique label corresponding to the mapped species.
    """

    multispecies_labels = {
        "vconcaerobin1": "Radius [µm]",
        "vconcaerobin2": "Radius [µm]",
        "vconcaerobin3": "Radius [µm]",
        "vconcaerobin4": "Radius [µm]",
        "vconcaerobin5": "Radius [µm]",
        "vconcaerobin6": "Radius [µm]",
        "vconcaerobin7": "Radius [µm]",
        "vconcaerobin8": "Radius [µm]",
        "vconcaerobin9": "Radius [µm]",
        "vconcaerobin10": "Radius [µm]",
        "vconcaerobin11": "Radius [µm]",
        "vconcaerobin12": "Radius [µm]",
        "vconcaerobin13": "Radius [µm]",
        "vconcaerobin14": "Radius [µm]",
        "vconcaerobin15": "Radius [µm]",
        "vconcaerobin16": "Radius [µm]",
        "vconcaerobin17": "Radius [µm]",
        "vconcaerobin18": "Radius [µm]",
        "vconcaerobin19": "Radius [µm]",
        "vconcaerobin20": "Radius [µm]",
        "vconcaerobin21": "Radius [µm]",
        "vconcaerobin22": "Radius [µm]",
    }

    networkspecies_aliases = [
        get_AERONET_sizedist_bin_radius(networkspeci)
        if networkspeci in multispecies_labels
        else networkspeci
        for networkspeci in networkspecies
    ]

    labels = np.unique(
        [
            multispecies_labels[networkspeci]
            for networkspeci in networkspecies
            if networkspeci in multispecies_labels
        ]
    )
    if len(labels) == 1:
        unique_label = labels[0]
    else:
        unique_label = ""

    return networkspecies_aliases, unique_label


def temp_axis_dict():
    """
    Return temporal mapping as a dictionary used for the plots.

    Returns
    -------
    map_dict : dict
        Numbering of months/days
    """

    map_dict = {
        "short": {
            "dayofweek": {0: "M", 1: "T", 2: "W", 3: "T", 4: "F", 5: "S", 6: "S"},
            "month": {
                1: "J",
                2: "F",
                3: "M",
                4: "A",
                5: "M",
                6: "J",
                7: "J",
                8: "A",
                9: "S",
                10: "O",
                11: "N",
                12: "D",
            },
        },
        "long": {
            "dayofweek": {
                0: "Monday",
                1: "Tuesday",
                2: "Wednesday",
                3: "Thursday",
                4: "Friday",
                5: "Saturday",
                6: "Sunday",
            },
            "month": {
                1: "January",
                2: "February",
                3: "March",
                4: "April",
                5: "May",
                6: "June",
                7: "July",
                8: "August",
                9: "September",
                10: "October",
                11: "November",
                12: "December",
            },
        },
    }

    return map_dict


def periodic_xticks():
    """
    Return xticks for periodic subplots based on temporal resolution.

    Returns
    -------
    xticks : dict
        Dictionary of xticks per temporal resolution
    """

    xticks = {
        "hour": np.arange(24, dtype=np.int8),
        "dayofweek": np.arange(7, dtype=np.int8),
        "month": np.arange(1, 13, dtype=np.int8),
    }

    return xticks


def periodic_labels():
    """
    Return axes labels for periodic subplots.

    Returns
    -------
    axes_labels : dict
        Axes labels per temporal resolution.
    """

    axes_labels = {"hour": "H", "dayofweek": "DoW", "month": "M"}

    return axes_labels


def get_land_polygon_resolution(selection):
    """
    Get resolution of land polygons to plot on map.

    Parameters
    ----------
    selection : str
        Selected temporal resolution.

    Returns
    -------
    resolution : str
        Selected land polygon resolution.
    """

    land_polygon_resolutions = {"low": "110m", "medium": "50m", "high": "10m"}
    resolution = land_polygon_resolutions[selection]

    return resolution


def update_plotting_parameters(
    instance, data_labels_to_remove=None, data_labels_to_add=None, daily_forecast=False
):
    """
    Update plotting parameters (colour, zorder, and grid edges)
    for data labels in a Report or Dashboard instance.

    Parameters
    ----------
    instance : object
        Instance of class Report or Dashboard.
    data_labels_to_remove : list, optional
        List of data labels to remove from plotting parameters.
    data_labels_to_add : list, optional
        List of data labels to add/update plotting parameters.
    daily_forecast : bool, optional
        Indicates whether daily forecast adjustments are needed.
    """

    # Reset plotting parameters if no labels to add or remove are specified
    if (data_labels_to_add is None) & (data_labels_to_remove is None):
        instance.plotting_params = {}
        data_labels_to_add = instance.data_labels  # Add all current data_labels

    # Add grid edge plotting parameters for data labels to add
    if data_labels_to_add is not None:
        # Extract unique base labels, stripping '-day', '-daily', '-combined' suffixes
        unique_base_data_labels = np.unique(
            [
                data_label.split("-day")[0].split("-daily")[0].split("-combined")[0]
                for data_label in data_labels_to_add
            ]
        )
        processed_data_labels = []  # Track labels that have already been processed

        for unique_base_data_label in unique_base_data_labels:
            if unique_base_data_label == instance.observations_data_label:
                # Observations do not require grid edges, but initialise plotting params dict
                instance.plotting_params[unique_base_data_label] = {}
            else:
                # Stop if all data labels have been processed
                if len(processed_data_labels) == len(data_labels_to_add):
                    break
                for valid_networkspeci in instance.networkspecies:
                    if len(processed_data_labels) == len(data_labels_to_add):
                        break

                    # read grid edge longitudes and latitudes if needed
                    if (not daily_forecast) & (
                        unique_base_data_label
                        in instance.files_to_read[valid_networkspeci]
                    ):
                        # Open netCDF file to extract grid edges
                        mod_nc_root = Dataset(
                            instance.files_to_read[valid_networkspeci][
                                unique_base_data_label
                            ][0]
                        )
                        grid_edge_longitude = mod_nc_root["grid_edge_longitude"][:]
                        grid_edge_latitude = mod_nc_root["grid_edge_latitude"][:]
                        # Close netCDF file
                        mod_nc_root.close()

                    # iterate through data labels to add
                    for data_label in data_labels_to_add:
                        base_data_label = (
                            data_label.split("-day")[0]
                            .split("-daily")[0]
                            .split("-combined")[0]
                        )

                        if (base_data_label == unique_base_data_label) & (
                            data_label not in processed_data_labels
                        ):
                            # initialise plotting params dict for data label to add
                            instance.plotting_params[data_label] = {}

                            # For daily forecast, copy grid edge info from removed labels
                            if daily_forecast:
                                relevant_inds = [
                                    ii
                                    for ii, data_label_to_remove in enumerate(
                                        data_labels_to_remove
                                    )
                                    if data_label in data_label_to_remove
                                ]
                                instance.plotting_params[data_label][
                                    "grid_edge_longitude"
                                ] = instance.plotting_params[
                                    data_labels_to_remove[relevant_inds[0]]
                                ][
                                    "grid_edge_longitude"
                                ]
                                instance.plotting_params[data_label][
                                    "grid_edge_latitude"
                                ] = instance.plotting_params[
                                    data_labels_to_remove[relevant_inds[0]]
                                ][
                                    "grid_edge_latitude"
                                ]
                                processed_data_labels.append(data_label)

                            # otherwise open netCDF file to extract grid edges
                            elif (
                                unique_base_data_label
                                in instance.files_to_read[valid_networkspeci]
                            ):
                                instance.plotting_params[data_label][
                                    "grid_edge_longitude"
                                ] = grid_edge_longitude
                                instance.plotting_params[data_label][
                                    "grid_edge_latitude"
                                ] = grid_edge_latitude
                                processed_data_labels.append(data_label)

    # Remove plotting parameters for labels marked for removal
    if data_labels_to_remove is not None:
        for data_label in data_labels_to_remove:
            del instance.plotting_params[data_label]

    # Add colour and zorder for observations
    instance.plotting_params[instance.observations_data_label][
        "colour"
    ] = instance.plot_characteristics_templates["general"]["obs_markerfacecolor"]
    instance.plotting_params[instance.observations_data_label][
        "zorder"
    ] = instance.plot_characteristics_templates["general"]["obs_zorder"]

    # Generate a list of RGB tuples for the number of models
    sns.reset_orig()  # Reset seaborn to default
    color_palette = instance.plot_characteristics_templates["general"][
        "legend_color_palette"
    ]
    color_palettes = yaml.safe_load(
        open(join(PROVIDENTIA_ROOT, "settings/color_palettes.yaml"))
    )

    if color_palette in color_palettes.keys():
        # Check that the number of colors matches the number of models
        if (len(instance.data_labels) - 1) > len(color_palettes[color_palette]):
            error = "Error: The number of models and palette colors should be equal. "
            error += f"Add more colors to your palette '{color_palette}' in settings/color_palettes.yaml "
            error += (
                "or change your legend_color_palette in the plot characteristics files."
            )
            instance.logger.error(error)
            sys.exit(1)
        else:
            clrs = sns.color_palette(color_palettes[color_palette])
    else:
        # If palette not in YAML, generate colors automatically
        clrs = sns.color_palette(color_palette, n_colors=len(instance.data_labels) - 1)

    # Add colours and zorder for each model (non-observations)
    model_ind = 1
    for data_label in instance.data_labels:
        if data_label != instance.observations_data_label:
            # Define colour for model
            instance.plotting_params[data_label]["colour"] = clrs[model_ind - 1]
            # Define zorder for model relative to observations
            instance.plotting_params[data_label]["zorder"] = (
                instance.plotting_params[instance.observations_data_label]["zorder"]
                + model_ind
            )
            # Update count of models
            model_ind += 1


def get_display_label(read_instance, data_label):
    """
    Get the text to actually show for a data label (legend, statsummary),
    applying any per-session rename set via double-clicking a legend entry
    (see legend_picker_func()/rename_legend_label() in
    dashboard_interactivity.py).

    data_label itself - the real identifier used for data selection,
    colour/style lookups, and everywhere else a data label is matched by
    value - is never touched by this; only what gets drawn as text. That's
    also why this doesn't affect the model pop-up menu, which shows
    data_label directly rather than going through this function.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard (or Report/Library, where this is
        always a no-op - legend_label_overrides is dashboard-only, reset on
        every data load, and never set outside the live dashboard session).
    data_label : str
        The data label's real identifier.

    Returns
    -------
    str
        data_label, or its override if one is set.
    """

    return getattr(read_instance, "legend_label_overrides", {}).get(
        data_label, data_label
    )


# fixed seed for get_deterministic_subsample() - not a secret or a tunable,
# just a constant so every call is reproducible
_DETERMINISTIC_SUBSAMPLE_SEED = 0


def get_deterministic_subsample(n, max_points):
    """
    Pick a fixed-size, reproducible random subset of indices into an array
    of length `n`, for plots that cap how many points they draw (a scatter
    cloud, a per-station Taylor diagram) once a selection is too large to
    render or read legibly in full.

    Seeded locally rather than drawn from the shared, global numpy random
    state: without a fixed seed, the exact same plot redrawn for a reason
    that has nothing to do with its own data (the window resizing, a
    sibling panel updating) silently swaps in a different random sample -
    a station a reader had spotted and was tracking moves or disappears for
    no reason connected to their own selection. A fixed seed instead
    reproduces the same subset every time the pool is the same size, and
    only ever changes because the underlying selection itself did. Using a
    local generator rather than reseeding the global one also means this
    can never be perturbed by unrelated code drawing random numbers
    elsewhere, nor perturb it in turn.

    Parameters
    ----------
    n : int
        Size of the full index range to sample from.
    max_points : int
        Number of indices to draw. If `n` is already this size or smaller,
        every index is returned and nothing is actually subsampled.

    Returns
    -------
    numpy.ndarray
        Sorted array of indices into an array of length `n`, `max_points`
        long (or `n` long, if there was nothing to cut down).
    """

    if n <= max_points:
        return np.arange(n)

    rng = np.random.default_rng(_DETERMINISTIC_SUBSAMPLE_SEED)
    return np.sort(rng.choice(n, size=max_points, replace=False))


def histogram_bin_target(data):
    """
    Work out what one set of data would want from a histogram's bins, as the
    widest bin it can carry and how far up the axis its values reach.

    Gathered across a report's subsections (the same way their data ranges
    are) so that every page can be drawn with one set of bins: bins following
    each subsection's own data would put every page on its own axis, and the
    counts on those pages are meant to be read against one another.

    Parameters
    ----------
    data : np.ndarray
        Values of one subsection, across all data labels

    Returns
    -------
    tuple of float
        Bin width the data would choose, and its upper inner Tukey fence.
        Both nan if the data cannot be binned
    """

    data = data[np.isfinite(data)]
    if (data.size == 0) or (np.nanmin(data) == np.nanmax(data)):
        return np.nan, np.nan

    _, upper_inner_fence = boxplot_inner_fences(data)

    # the bins this data would be given on its own (numpy's "auto" - the
    # larger of the Freedman-Diaconis and Sturges counts), as a width, which
    # unlike a count can be carried over to the range shared by every page
    edges = np.histogram_bin_edges(data, bins="auto")
    if len(edges) < 2:
        return np.nan, upper_inner_fence

    return float(edges[1] - edges[0]), float(upper_inner_fence)


def kde_fft(
    xin,
    gridsize=1024,
    extents=None,
    weights=None,
    adjust=1.0,
    bw="scott",
    xgrid=None,
    min_bandwidth=None,
):
    """
    A fft-based Gaussian kernel density estimate (KDE)
    for computing the KDE on a regular grid

    Note that this is a different use case than scipy's original
    scipy.stats.kde.gaussian_kde

    Parameters
    ----------
    xin : ndarray, shape (n,)
        1D array of data samples.
    gridsize : int, optional
        Size of the output grid (default is 1024).
    extents : tuple of float, optional
        Tuple of (xmin, xmax) specifying the range of the output grid.
        Defaults to the extent of the input data.
    weights : ndarray, shape (n,), optional
        Weights for each sample. If None, all samples are weighted equally.
    adjust : float, optional
        Adjustment factor for bandwidth. Bandwidth becomes `bw * adjust`.
    bw : {'scott', 'silverman'}, optional
        Method used to calculate bandwidth (default is 'scott').
    xgrid : ndarray, shape (n,), optional
        If provided, this grid will be used for KDE evaluation. Overrides `gridsize` and `extents`.
    min_bandwidth : float, optional
        Floor for the kernel's standard deviation, in the same units as xin
        (e.g. a species' instrument reporting resolution). Below it, the KDE
        reproduces the comb of rounded observations rather than smoothing
        over it, seen as ringing on the curve. None leaves the automatic
        bandwidth unmodified

    Returns
    -------
    grid_points : ndarray, shape (gridsize,)
        Grid points where the KDE is evaluated. Returned only if `xgrid` is None.
    kde : ndarray, shape (gridsize,)
        KDE evaluated at `grid_points` or `xgrid`.
    """
    # variable check
    x = np.squeeze(np.asarray(xin))

    # default extents are the extent of the data
    if xgrid is not None:
        xmin, xmax = xgrid.min(), xgrid.max()
        x = x[(x <= xmax) & (x >= xmin)]
    elif extents is None:
        xmin, xmax = x.min(), x.max()
    else:
        xmin, xmax = map(float, extents)
        x = x[(x <= xmax) & (x >= xmin)]

    n = x.size

    # apply weights
    if weights is None:
        # Default: Weight all points equally
        weights = np.ones(n)
    else:
        weights = np.squeeze(np.asarray(weights))
        if weights.size != x.size:
            error = "Input weights must be an array of the same size as xin. "
            return error

    # make grid and discretize the data and round it to the next power of 2
    # to optimize with the fft usage
    # ensure minimum gridsize of 1024 points
    if xgrid is None:
        if gridsize is None:
            gridsize = np.max((len(x), 1024.0))
        gridsize = 2 ** np.ceil(np.log2(gridsize))  # round to next power of 2
        nx = int(gridsize)
    else:
        nx = len(xgrid)

    # make the sparse histogram
    dx = (xmax - xmin) / (nx - 1)

    # basically, this is just doing what np.digitize does with one less copy
    xyi = x - xmin
    xyi /= dx
    xyi = np.floor(xyi, xyi)
    xyi = np.vstack((xyi, np.zeros(n, dtype=int)))

    # next, make a histogram of x
    # exploit a sparse coo_matrix avoiding np.histogram due to excessive
    # memory usage with many points
    try:
        grid = coo_matrix((weights, xyi), shape=(nx, 1)).toarray()
    except ValueError:
        error = "Too many zeros. "
        return error

    # Kernel Preliminary Calculations
    std_x = np.std(xyi[0])

    # Scaling factor for bandwidth
    if bw == "scott":
        bw_factor = (n ** (-1.0 / 5.0)) * adjust
    elif bw == "silverman":
        bw_factor = ((n * 3 / 4.0) ** (-1.0 / 5)) * adjust

    # kernel standard deviation, in grid cells
    kernel_std = bw_factor * std_x

    # floor the kernel width at min_bandwidth, converted from data units to
    # grid cells - a large enough sample count otherwise shrinks the automatic
    # bandwidth below the data's own reporting resolution
    if min_bandwidth is not None and min_bandwidth > 0:
        kernel_std = max(kernel_std, min_bandwidth / dx)

    # make the gaussian kernel
    # first, determine the bandwidth using defined bandwidth estimator rule
    kern_nx = int(np.round(kernel_std * 2 * np.pi))

    # If bandwidth is 0, skip plot for current data label
    if kern_nx == 0:
        error = "The kernel bandwidth is 0.  "
        error += "To change the bandwith, we recommend increasing the number of "
        error += "pdf_min_samples in the plot characteristics settings files. "
        return error

    # Then evaluate the gaussian function on the kernel grid
    kernel = np.reshape(gaussian(kern_nx, kernel_std), (kern_nx, 1))

    # convolve the histogram with the gaussian kernel
    # use symmetric padding to correct for data boundaries in the kde
    npad = np.min((nx, 2 * kern_nx))
    grid = np.vstack([grid[npad:0:-1], grid, grid[nx : nx - npad : -1]])
    grid = convolve(grid, kernel, mode="same")[npad : npad + nx]

    # normalization factor to divide result by so that units are in the same
    # units as scipy.stats.kde.gaussian_kde's output. Uses kernel_std (the
    # actual, possibly min_bandwidth-floored kernel width) rather than
    # bw_factor * std_x directly, so normalization stays consistent with
    # whichever kernel was actually used above.
    norm_factor = 2 * np.pi * kernel_std**2
    norm_factor = n * dx * np.sqrt(norm_factor)

    # normalize the result
    grid /= norm_factor

    # return grid points and estimated densities
    if xgrid is None:
        return np.linspace(xmin, xmax, nx), np.squeeze(grid)
    else:
        return np.squeeze(grid)


def round_decimal_places(x, decimal_places):
    """
    Round a number to a specified number of decimal places.

    Parameters
    ----------
    x : float
        Value to round.
    decimal_places : int
        Desired number of decimal places.

    Returns
    -------
    str
        Rounded value as a string. Uses scientific notation if necessary.
        Returns `'nan'` if input is NaN.
    """

    # if cell value is not nan
    if not np.isnan(x):
        # if number is zero, return int as str
        if x == 0:
            return "0"
        # if number of zeros is more than decimal places set by user, use scientific notation
        elif (-math.floor(math.log10(abs(x))) - 1) > decimal_places:
            return "{:0.{}e}".format(x, decimal_places)
        # if not, use float
        else:
            return "{:0.{}f}".format(x, decimal_places)
    # if cell value is nan, do nothing
    else:
        return "nan"


def merge_cells(table, cells, visibility=False):
    """
    Merge cells in a matplotlib Table.
    Reference: https://stackoverflow.com/a/53819765/12684122

    Parameters
    ----------
    table : matplotlib.Table
        The table object containing the cells.
    cells : list of tuple
        List of table coordinates to merge, e.g. [(0, 1), (0, 0), (0, 2)].
    visibility : bool
        Whether we want to show or hide the texts in the cells other than 0.
    """

    cells_array = [np.asarray(c) for c in cells]
    h = np.array(
        [cells_array[i + 1][0] - cells_array[i][0] for i in range(len(cells_array) - 1)]
    )
    v = np.array(
        [cells_array[i + 1][1] - cells_array[i][1] for i in range(len(cells_array) - 1)]
    )

    # if it's a horizontal merge, all values for `h` are 0
    if not np.any(h):
        # sort by horizontal coord
        cells = np.array(sorted(list(cells), key=lambda v: v[1]))
        edges = ["BTL"] + ["BT" for i in range(len(cells) - 2)] + ["BTR"]
    elif not np.any(v):
        cells = np.array(sorted(list(cells), key=lambda h: h[0]))
        edges = ["TRL"] + ["RL" for i in range(len(cells) - 2)] + ["BRL"]
    else:
        raise ValueError("Only horizontal and vertical merges allowed")

    for cell, e in zip(cells, edges):
        table[cell[0], cell[1]].visible_edges = e

    txts = [table[cell[0], cell[1]].get_text() for cell in cells]
    tpos = [np.array(t.get_position()) for t in txts]

    # transpose the text of the left cell
    trans = (tpos[-1] - tpos[0]) / 2

    # didn't had to check for ha because I only want ha='center'
    txts[0].set_transform(transforms.Affine2D().translate(*trans))
    if not visibility:
        for txt in txts[1:]:
            txt.set_visible(False)


def get_taylor_diagram_ghelper_info(
    reference_stddev, plot_characteristics, extend=False
):
    """
    Compute Taylor diagram axis extremes and labels.

    Parameters
    ----------
    reference_stddev : float
        Reference standard deviation used to set the limits.
    plot_characteristics : dict
        Dictionary containing plot characteristics.
    extend : bool, optional
        If True, the diagram shows negative correlations; otherwise, only positive correlations are shown (default is False).

    Returns
    -------
    tmin : float
        Minimum angle (radians) for the polar plot.
    tmax : float
        Maximum angle (radians) for the polar plot.
    smin : float
        Minimum standard deviation for the plot axis.
    smax : float
        Maximum standard deviation for the plot axis.
    gl1 : matplotlib.ticker.FixedLocator
        Locator for correlation labels.
    tf1 : matplotlib.ticker.DictFormatter
        Formatter for correlation labels.
    """

    tmin = 0

    # diagram extended to negative correlations
    if extend:
        tmax = np.pi
    # diagram limited to positive correlations
    else:
        tmax = np.pi / 2

    # get standard deviation axis extent
    srange = plot_characteristics["srange"]
    smin = srange[0] * reference_stddev
    smax = srange[1] * reference_stddev

    # correlation labels
    if extend:
        rlocs = np.array(plot_characteristics["rlocs"]["all"])
        rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
    else:
        rlocs = np.array(plot_characteristics["rlocs"]["positive"])

    # convert correlation values into polar angles
    tlocs = np.arccos(rlocs)
    gl1 = gf.FixedLocator(tlocs)
    tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    return tmin, tmax, smin, smax, gl1, tf1


def get_taylor_diagram_ghelper(reference_stddev, plot_characteristics, extend=False):
    """
    Create a Taylor diagram grid helper.

    Parameters
    ----------
    reference_stddev : float
        Reference standard deviation used to set the axis limits.
    plot_characteristics : dict
        Dictionary containing plot characteristics.
    extend : bool, optional
        If True, the diagram includes negative correlations,
        otherwise, only positive correlations are shown (default is False).

    Returns
    -------
    ghelper : mpl_toolkits.axisartist.grid_finder.GridHelperCurveLinear
        Grid helper object for the Taylor diagram.
    """

    # get axis extremes
    tmin, tmax, smin, smax, gl1, tf1 = get_taylor_diagram_ghelper_info(
        reference_stddev, plot_characteristics, extend
    )

    # get grid helper
    ghelper = fa.GridHelperCurveLinear(
        PolarAxes.PolarTransform(apply_theta_transforms=False),
        extremes=(tmin, tmax, smin, smax),
        grid_locator1=gl1,
        tick_formatter1=tf1,
    )

    return ghelper


# bounds on the automatic map marker size (see get_map_marker_size()).
#
# The floor is absolute: below it a point stops being visible at all, and that
# is just as true on a large figure as a small one. It is the one bound kept
# fixed across modes, so a dense map reads the same everywhere.
#
# The ceiling scales with the panel instead. Held fixed it bound far earlier on
# the library's large figure than on a report's small panel, so the same map
# came out relatively sparser in the library - the ceiling, not the sizing,
# was what made the modes disagree. It is only ever scaled up, never below
# MAP_MARKER_MAX_SIZE, so a small panel keeps a readable size, and never past
# MAP_MARKER_ABSOLUTE_MAX_SIZE, so a near-empty map cannot produce blobs
MAP_MARKER_MIN_SIZE = 6
MAP_MARKER_MAX_SIZE = 60
MAP_MARKER_ABSOLUTE_MAX_SIZE = 250
# panel area, in pixels squared, MAP_MARKER_MAX_SIZE is the ceiling for - a
# report map panel, the smallest of the three modes
MAP_MARKER_REFERENCE_AREA = 175000.0
# how much of the space each station has to itself its marker fills. The
# marker area is set to this share of the area per station, so points keep the
# same look however large the panel is and however many stations there are -
# raise it to make every map's points bigger
MAP_MARKER_AREA_FRACTION = 0.18


def get_map_marker_size(ax, datacrs, longitudes, latitudes, map_extent=None):
    """
    Get the marker size for a map, from how densely the stations being shown
    are packed into the axis. Shared by every mode, so the same map reads the
    same way in the dashboard, a report and the library.

    Counting only the stations on show is what makes this respond to zoom:
    zooming in leaves fewer of them on the same axis, so the markers grow,
    without needing a separate zoom term. Panel size is carried by the axis
    area, so the report's small 2x2 panels and the library's full-width figure
    each get a size suited to them.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Map axis the stations are plotted on
    datacrs : cartopy.crs.CRS
        CRS the station coordinates are given in
    longitudes : numpy array
        Station longitudes
    latitudes : numpy array
        Station latitudes
    map_extent : array-like, shape (4,), optional
        Extent the map will be shown at, as [lon_min, lon_max, lat_min,
        lat_max]. Given by report and library, which size the markers before
        the extent has been applied to the axis, so the axis cannot be asked.
        The dashboard leaves this as None and is counted against the axis's
        own current view, which is always live there.

    Returns
    -------
    float
        Marker size, in points squared
    """

    longitudes = np.asarray(longitudes)
    latitudes = np.asarray(latitudes)

    panel_area = ax.bbox.width * ax.bbox.height
    max_size = float(
        np.clip(
            MAP_MARKER_MAX_SIZE * panel_area / MAP_MARKER_REFERENCE_AREA,
            MAP_MARKER_MAX_SIZE,
            MAP_MARKER_ABSOLUTE_MAX_SIZE,
        )
    )

    if longitudes.size == 0 or panel_area <= 0:
        return max_size

    if map_extent is not None and len(map_extent) == 4:
        # count against the extent the map will be shown at, in lon/lat, which
        # is how map_extent is given
        lon_min, lon_max, lat_min, lat_max = map_extent
        coords_x, coords_y = longitudes, latitudes
        view_x = (min(lon_min, lon_max), max(lon_min, lon_max))
        view_y = (min(lat_min, lat_max), max(lat_min, lat_max))
    else:
        # count against the axis's current view, in the axis's own projected
        # coordinates rather than lon/lat, so this holds for every projection
        projected = ax.projection.transform_points(datacrs, longitudes, latitudes)
        coords_x, coords_y = projected[:, 0], projected[:, 1]
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        view_x = (min(xlim), max(xlim))
        view_y = (min(ylim), max(ylim))

    shown = (
        (coords_x >= view_x[0])
        & (coords_x <= view_x[1])
        & (coords_y >= view_y[0])
        & (coords_y <= view_y[1])
    )
    n_points = int(np.count_nonzero(shown))

    # a lone station should always be plainly visible
    if n_points == 1:
        return max_size

    if n_points == 0:
        # nothing on show either means the map genuinely holds no stations, or
        # the axis has not been given its limits yet - a fresh axis reports
        # (0, 1) on both, which no station falls inside. Fall back to counting
        # them all across the whole panel, rather than defaulting to the
        # sparsest possible size
        n_points = int(longitudes.size)
        area = panel_area
    else:
        # measure the density over the part of the panel the stations actually
        # occupy, not the whole of it - a network gathered in one region is
        # crowded there however much empty map surrounds it, and dividing by
        # the full panel made those points come out as large as a sparse
        # global network's. Spans are taken between the 2nd and 98th
        # percentile so a single distant station cannot stretch the area and
        # hide the crowding in the rest
        span_x = float(np.ptp(np.percentile(coords_x[shown], [2, 98])))
        span_y = float(np.ptp(np.percentile(coords_y[shown], [2, 98])))
        fraction_x = np.clip(span_x / (view_x[1] - view_x[0]), 0.0, 1.0)
        fraction_y = np.clip(span_y / (view_y[1] - view_y[0]), 0.0, 1.0)

        # keep a floor on the occupied area, so a very tight cluster still
        # divides by something rather than by (near) zero
        area = max(panel_area * fraction_x * fraction_y, 16.0)

    # give each marker a fixed share of the area its own station has to
    # itself, so the spacing between points reads the same whatever the panel
    # size or station count. Replaces the exponential of density this used to
    # apply (https://github.com/BSC-ES/providentia/issues/199), which fell
    # away far too steeply once a map held a few hundred stations - a Spanish
    # or European report map bottomed out at the floor while its points were
    # still visibly well separated
    area_per_station = area / n_points

    # matplotlib marker sizes are in points squared, the axis area in pixels
    dpi = ax.figure.dpi if ax.figure.dpi else 100.0
    area_per_station *= (72.0 / dpi) ** 2

    return float(
        np.clip(
            MAP_MARKER_AREA_FRACTION * area_per_station,
            MAP_MARKER_MIN_SIZE,
            max_size,
        )
    )


def set_map_extent(canvas_instance, ax, map_extent):
    """
    Set map extent, done set_xlim and set_ylim rather than set_extent
    to avoid axis cutting off slightly (https://github.com/SciTools/cartopy/issues/697).

    Parameters
    ----------
    canvas_instance : object
        Canvas instance containing the data CRS (`datacrs` attribute).
    ax : matplotlib.axes._subplots.AxesSubplot
        Axis on which to set the map extent.
    map_extent : array-like, shape (4,)
        Map extent in the form [lon_min, lon_max, lat_min, lat_max].
    """

    mlon = np.mean(map_extent[:2])
    mlat = np.mean(map_extent[2:])
    xtrm_data = np.array(
        [
            [map_extent[0], mlat],
            [mlon, map_extent[2]],
            [map_extent[1], mlat],
            [mlon, map_extent[3]],
        ]
    )
    proj_to_data = canvas_instance.datacrs._as_mpl_transform(ax) - ax.transData
    xtrm = proj_to_data.transform(xtrm_data)
    ax.set_xlim(xtrm[:, 0].min(), xtrm[:, 0].max())
    ax.set_ylim(xtrm[:, 1].min(), xtrm[:, 1].max())


def get_map_extent(canvas_instance):
    """
    Get the map extent from the xlim and ylim of a map axis.

    Parameters
    ----------
    canvas_instance : object
        Canvas instance containing the plot axes (`plot_axes`) and the data CRS (`datacrs`).

    Returns
    -------
    map_extent : list of float
        Map extent in the form [lon_min, lon_max, lat_min, lat_max].
    """

    # get plot extent
    ax = canvas_instance.plot_axes["map"]
    coords = np.array(ax.get_extent())
    current_xlim = coords[0:2]
    current_ylim = coords[2:4]

    # calculate means
    mlon = np.mean(current_xlim)
    mlat = np.mean(current_ylim)

    # get coordinates
    xcoords = np.array([current_xlim[0], mlon, current_xlim[1], mlon])
    ycoords = np.array([mlat, current_ylim[0], mlat, current_ylim[1]])

    # transform coordinates to projected data
    transformed_coords = canvas_instance.datacrs.transform_points(
        canvas_instance.plotcrs, xcoords, ycoords
    )[:, :2]

    # keep longitudes between -180 and 180
    lon_change = False
    if (np.isnan(transformed_coords[0, 0])) or (
        transformed_coords[0, 0] == -179.99999999999932
    ):
        transformed_coords[0, 0] = -180
        lon_change = True
    if (np.isnan(transformed_coords[2, 0])) or (
        transformed_coords[2, 0] == 179.99999999999932
    ):
        transformed_coords[2, 0] = 180
        lon_change = True

    # keep latitudes between -90 and 90
    lat_change = False
    if (np.isnan(transformed_coords[1, 1])) or (
        transformed_coords[1, 1] == -89.99999999999966
    ):
        transformed_coords[1, 1] = -90
        lat_change = True
    if (np.isnan(transformed_coords[3, 1])) or (
        transformed_coords[3, 1] == 89.99999999999966
    ):
        transformed_coords[3, 1] = 90
        lat_change = True

    # recalculate means
    if lon_change or lat_change:
        # recalculate longitude means
        mlon = np.mean(np.array([transformed_coords[0, 0], transformed_coords[2, 0]]))
        transformed_coords[1, 0] = mlon
        transformed_coords[3, 0] = mlon

        # recalculate latitude means
        mlat = np.mean(np.array([transformed_coords[1, 1], transformed_coords[3, 1]]))
        transformed_coords[0, 1] = mlat
        transformed_coords[2, 1] = mlat

    # get map extent
    map_extent = [
        transformed_coords[:, 0].min(),
        transformed_coords[:, 0].max(),
        transformed_coords[:, 1].min(),
        transformed_coords[:, 1].max(),
    ]

    return map_extent


def create_statistical_timeseries(
    read_instance,
    canvas_instance,
    chunk_stat,
    chunk_resolution,
    networkspeci,
    cut_data_labels,
    bias,
):
    """
    Create statistical timeseries data by chunk resolution.

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    canvas_instance : object
        Instance of class Canvas or Report.
    chunk_stat : str
        Statistic to calculate for each chunk.
    chunk_resolution : str
        Temporal resolution of the chunks.
    networkspeci : str
        Current network-species combination.
    cut_data_labels : list
        List of valid data labels to include in the calculation.
    bias : bool
        Indicates if bias statistics should be calculated.

    Returns
    -------
    timeseries_data : pandas.DataFrame
        DataFrame containing the calculated statistics for each chunk and data label.
    """

    z_statistic_sign = get_z_statistic_sign(chunk_stat)

    # if have forecast data in memory, set the forecast type
    forecast_type = None
    if read_instance.daily_forecast:
        forecast_type = "daily"
    elif read_instance.combined_forecast:
        forecast_type = "combined"

    if (z_statistic_sign == "bias") or (bias):
        if read_instance.observations_data_label in cut_data_labels:
            cut_data_labels.remove(read_instance.observations_data_label)
        stats_calc = calculate_statistic(
            read_instance,
            canvas_instance,
            networkspeci,
            chunk_stat,
            [read_instance.observations_data_label] * len(cut_data_labels),
            cut_data_labels,
            chunk_resolution=chunk_resolution,
            statistic_aggregation=read_instance.timeseries_statistic_aggregation,
            forecast_type=forecast_type,
        )
    else:
        stats_calc = calculate_statistic(
            read_instance,
            canvas_instance,
            networkspeci,
            chunk_stat,
            cut_data_labels,
            [],
            chunk_resolution=chunk_resolution,
            statistic_aggregation=read_instance.timeseries_statistic_aggregation,
            forecast_type=forecast_type,
        )

    if (chunk_resolution == read_instance.active_resolution) & (
        forecast_type == "combined"
    ):
        chunk_dates = read_instance.time_index_ts
    elif (chunk_resolution == read_instance.active_resolution) & (
        forecast_type != "daily"
    ):
        chunk_dates = read_instance.time_index
    else:
        chunk_dates = canvas_instance.grouped_ts_index

    timeseries_data = pd.DataFrame(
        index=chunk_dates, columns=cut_data_labels, dtype=np.float32
    )

    # if shape of stats_calc is not correct return error
    if (stats_calc.shape[0] != len(chunk_dates)) or (
        stats_calc.shape[1] != len(cut_data_labels)
    ):
        error = "Error: The shape of the calculated statistical timseseries does not match the expected shape."
        read_instance.logger.error(error)
        sys.exit(1)

    for chunk_date_idx, chunk_date in enumerate(chunk_dates):
        for label_idx, data_label in enumerate(cut_data_labels):
            timeseries_data.loc[chunk_date, data_label] = np.float32(
                stats_calc[chunk_date_idx, label_idx]
            )

    return timeseries_data


def get_hex_code(colour):
    """
    Convert a colour name or RGB decimal tuple to a hexadecimal colour code.

    Parameters
    ----------
    colour : str or tuple
        Colour input. Can be a string representing a named
        colour or a tuple of RGB decimals (each between 0 and 1).

    Returns
    -------
    hex_colour : str
        Hexadecimal colour code.
    """

    # convert from colour name to hex code
    if type(colour) is str:
        if colour[0] != "#":
            hex_colour = cnames[colour]
    # convert from rgb colour (as decimal) to hex code
    elif type(colour) is tuple:
        rgb_colour = tuple(round(255 * x) for x in colour)
        hex_colour = f"#{int(round(rgb_colour[0])):02x}{int(round(rgb_colour[1])):02x}{int(round(rgb_colour[2])):02x}"

    return hex_colour


def get_fairmode_RV_exceendance(read_instance, speci, RV, exc_threshold, units):
    """
    Convert standard GHOST units of RV and exceedances to actual ones

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    speci : str
        Speci to plot.
    RV : float
        Reference value.
    exc_threshold : float
        Exceedance threshold.
    units : str
        Units of the input RV and exceedance threshold.

    Returns
    -------
    RV : float
        Reference value in GHOST units.
    exc_threshold : float
        Exceedance threshold in GHOST units.
    """

    # get input and output units
    standard_parameter_speci = get_standard_parameters_by_speci(
        speci, read_instance.ghost_version
    )
    initial_units = units
    final_units = read_instance.measurement_units[speci]

    # convert units using conversion factor
    conversion_factor = get_conversion_factor(
        initial_units, final_units, standard_parameter_speci
    )
    if isinstance(conversion_factor, str):
        read_instance.logger.error(conversion_factor)
        sys.exit(1)
    RV *= conversion_factor

    # threshold can be None, as in the case of pm2p5
    if exc_threshold is not None:
        exc_threshold *= conversion_factor

    return RV, exc_threshold


def get_multispecies_conversion_factor(read_instance, speci):
    """
    Get conversion factor from original data units to multispecies units defined in configuration
    given a species

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    speci : str
        Speci to plot.

    Returns
    -------
    None, float
        Conversion factor
    """

    # get input and output units
    standard_parameter_speci = get_standard_parameters_by_speci(
        speci, read_instance.ghost_version
    )
    initial_units = read_instance.measurement_units[speci]
    final_units = read_instance.multispecies_units

    # do not convert if units do not change
    if initial_units == final_units:
        return None

    # convert units using conversion factor
    conversion_factor = get_conversion_factor(
        initial_units, final_units, standard_parameter_speci
    )
    if isinstance(conversion_factor, str):
        read_instance.logger.error(conversion_factor)
        sys.exit(1)

    # check if conversion factor is nan
    if np.isnan(conversion_factor):
        msg = f"Conversion factor is nan, units for {speci} could not be converted from {initial_units} to {final_units}, "
        msg += "multispecies plot will show data in different units."
        show_message(read_instance, msg)
        return None

    return conversion_factor


def convert_multispecies_df_units(read_instance, stats_df, zstats, base_plot_type):
    """
    Convert units in multispecies dataframe

    Parameters
    ----------
    read_instance : object
        Instance of class Dashboard or Report.
    stats_df : pandas.Dataframe
        Multispecies statistics dataframe

    Returns
    -------
    pandas.Dataframe
        Multispecies statistics dataframe with converted units
    """

    for networkspeci in read_instance.networkspecies:
        speci = networkspeci.split("|")[1]
        conversion_factor = get_multispecies_conversion_factor(read_instance, speci)
        if conversion_factor is not None:
            mask = stats_df.index.get_level_values("networkspecies") == networkspeci

            if base_plot_type in ["statsummary"]:
                # get stats that are in measurement units
                stat_cols = []
                stat_names = []
                for col in zstats:
                    # get units
                    base_zstat = col.split("_bias")[0].split("-")[0]
                    z_statistic_type = get_z_statistic_type(base_zstat)
                    if z_statistic_type == "basic":
                        units = read_instance.basic_stats[base_zstat]["units"]
                    else:
                        units = read_instance.modbias_stats[base_zstat]["units"]

                    # get stat position in dataframe and collect stats that are in measurement units
                    col_position = stats_df.columns.get_loc(col)
                    if units == "[measurement_units]":
                        stat_cols.append(col_position)
                        stat_names.append(col)

                # convert stats that are in measurement units
                if len(stat_cols) > 0:
                    msg = f"Converting units of {stat_names} for {networkspeci} to {read_instance.multispecies_units}."
                    show_message(read_instance, msg)
                    stats_df.iloc[mask, stat_cols] *= conversion_factor

            elif base_plot_type in ["heatmap", "table"]:
                # get units
                if isinstance(zstats, list):
                    zstat = zstats[0]
                else:
                    zstat = zstats
                base_zstat = zstat.split("_bias")[0].split("-")[0]
                z_statistic_type = get_z_statistic_type(zstat)
                if z_statistic_type == "basic":
                    units = read_instance.basic_stats[base_zstat]["units"]
                else:
                    units = read_instance.modbias_stats[base_zstat]["units"]

                # check if stat is in measurement units and convert stat if it is
                if units == "[measurement_units]":
                    msg = f"Converting units of {zstat} for {networkspeci} to {read_instance.multispecies_units}."
                    show_message(read_instance, msg)
                    stats_df.iloc[mask, :] *= conversion_factor

    return stats_df


def handle_test_or_save_df(
    read_instance, df, filename, path, tests_generate_output, msgs, decimal_places
):
    """
    Save dataframe or assert if dataframe generates the same outputs as the dataframes saved in tests folder

    Parameters
    ----------
    read_instance : object
        An instance of the application class responsible for data reading operations.
    df : pd.DataFrame
        Data to be downloaded
    filename : str
        Filename
    path : str
        Path to save file
    tests_generate_output : bool
        Indicates if we want to regenerate dataframes saved in tests folder
    msgs : list
        Text to show after downloading file
    decimal_places : int
        Decimal places to round the data to when saving dataframe
    """

    df = df.round(decimal_places)
    if read_instance.tests:
        generated_output = df
        generated_output = generated_output.replace("", np.nan)
        read_instance.logger.info("Generated_output")
        read_instance.logger.info(generated_output)
        if tests_generate_output:
            # during tests save empty spaces in dataframe as 'nan' because pd.read_csv creates nans on indices
            # and it is not possible to compare the dataframes
            # it is also possible to modify the behaviour of pd.read_csv to avoid reading nans defining keep_default_na=False
            # but in that way we hide the nans in the data and not only in the indices
            msgs = save_df(read_instance, df, filename, path, msgs, na_rep=np.nan)

        # read expected output
        parse_dates = []
        if "time" in generated_output.columns:
            parse_dates.append("time")
        expected_output = pd.read_csv(f"{path}/{filename}.csv", parse_dates=parse_dates)
        read_instance.logger.info(f'Expected_output ({f"{path}/{filename}.csv"})')
        read_instance.logger.info(expected_output)
        if "metadata" in filename:
            expected_output["value"] = expected_output["value"].astype(str)
        assert assert_frame_equal(generated_output, expected_output, atol=1e-5) is None

    else:
        msgs = save_df(read_instance, df, filename, path, msgs)

    return msgs


def save_df(read_instance, df, filename, path, msgs, na_rep=""):
    """
    Save dataframe to CSV file

    Parameters
    ----------
    read_instance : object
        An instance of the application class responsible for data reading operations.
    df : pd.DataFrame
        Data to be downloaded
    plot_type : str
        Plot type
    filename : str
        Filename
    data_label : str
        Data label
    element_type : str
        Element type.
    path : str
        Path to save file
    msgs : list
        Text to show after downloading file
    na_rep : str, np.nan, optional
        Representation of nans in saved dataframe
    """

    fname = join(PROVIDENTIA_ROOT, f"{path}/{filename}.csv")
    # in tests do not ask
    if os.path.exists(fname) and not read_instance.tests:
        overwrite_question = (
            f"File already exists in {fname}. Do you want to overwrite it?"
        )
        if read_instance.mode == "library":
            # ask the user whether they want to overwrite the file
            while True:
                overwrite_file = input(f"\n{overwrite_question} ([y]/n) ")
                if overwrite_file in ["", "y", "n"]:
                    if overwrite_file in ["n"]:
                        return msgs
                    break
        elif read_instance.mode == "dashboard":
            popup = MessageBox(overwrite_question, confirmation=True)
            if not popup.result:
                return msgs

    df.to_csv(fname, index=False, na_rep=na_rep)
    msg = f"\n- File {fname}"
    msgs.append(msg)

    return msgs


def download_plot_data_to_csv(
    read_instance,
    canvas_instance,
    base_plot_type,
    plot_type,
    plot_options,
    path,
    networkspeci,
    tests_generate_output=False,
    labela="",
    labelb="",
):
    """
    Extract data into dataframe and save to to CSV file

    Parameters
    ----------
    read_instance : object
        An instance of the application class responsible for data reading operations.
    canvas_instance : object
        An instance of the application class responsible for plotting operations.
    base_plot_type : str
        Base plot type without statistical information.
    plot_type : str
        Plot type
    labela : str
        Label of first dataset.
    labelb : str
        Label of second dataset (if defined then a bias plot is made).
    plot_options : list
        Plot options
    tests_generate_output : bool, optional
        Indicates if tests need to regenerate CSV files with plot data
    """

    if "bias" in plot_options:
        plot_element_varname = "bias"
    else:
        plot_element_varname = "absolute"
    if plot_element_varname not in canvas_instance.plot_elements[base_plot_type]:
        return

    # ask user which element types they want to download
    # keep only median for periodic violin plots
    element_types = list(
        {
            key
            for data_label in canvas_instance.plot_elements[base_plot_type][
                plot_element_varname
            ]
            for key in canvas_instance.plot_elements[base_plot_type][
                plot_element_varname
            ][data_label].keys()
            if not (
                (key.startswith("violin_plot"))
                and (base_plot_type == "periodic-violin")
            )
        }
    )

    element_types_to_save = []
    if read_instance.mode == "library":
        # in tests do not ask
        if read_instance.tests:
            element_types_to_save = element_types
        # ask the user whether they want to save specific element_types
        else:
            for element_type in element_types:
                while True:
                    create_file = input(
                        f"\nDo you want to save {element_type} data? ([y]/n) "
                    )
                    if create_file in ["", "y", "n"]:
                        if create_file in ["", "y"]:
                            element_types_to_save.append(element_type)
                        break
    elif read_instance.mode == "dashboard":
        # in dashboard open dialog to ask the user which element types they want to save
        # if there are more than 1 element types
        if len(element_types) > 1:
            dialog = CheckDialog(element_types)
            if dialog.exec_():
                element_types_to_save = dialog.get_checked_items()
        # do not open dialog if there is only one plot element
        elif len(element_types) == 1:
            element_types_to_save = element_types
        # do not continue if user selected no elements
        else:
            return

    # if no element types to save, do not continue
    if len(element_types_to_save) == 0:
        return

    # active statistic, from plot_characteristics (dashboard) or plot_type
    # (report/library, e.g. "distribution-r") - for a meaningful column name
    (parsed_zstat, parsed_base_zstat, _, _, _) = get_z_statistic_info(
        plot_type=plot_type
    )

    station_statistic = None
    if base_plot_type in ["distribution", "histogram"]:
        station_statistic = parsed_zstat or canvas_instance.plot_characteristics[
            plot_type
        ].get("station_statistic")
        if station_statistic in (None, "", "None"):
            station_statistic = None

    taylor_corr_stat = None
    if base_plot_type == "taylor":
        taylor_corr_stat = (
            canvas_instance.plot_characteristics[plot_type].get("corr_stat")
            or parsed_base_zstat
            or "r"
        )

    stats_dict = {**read_instance.basic_stats, **read_instance.modbias_stats}

    x_column = (
        "time"
        if base_plot_type == "timeseries"
        else read_instance.observations_data_label
        if base_plot_type == "scatter"
        else stats_dict[station_statistic]["label"]
        if (base_plot_type in ["distribution", "histogram"])
        and (station_statistic in stats_dict)
        else "concentration"
        if base_plot_type in ["distribution", "histogram"]
        else taylor_corr_stat
        if base_plot_type == "taylor"
        else canvas_instance.plot_characteristics[plot_type]["xlabel"]["xlabel"]
        if base_plot_type == "fairmode-target"
        else "x"
    )
    decimal_places = canvas_instance.plot_characteristics[plot_type][
        "round_decimal_places"
    ]["csv"]

    msgs = []
    combined_dfs = {}
    boxplot_accumulator = {}

    for data_label in canvas_instance.plot_elements[base_plot_type][
        plot_element_varname
    ]:
        for element_type in canvas_instance.plot_elements[base_plot_type][
            plot_element_varname
        ][data_label]:
            # ignore element types that user does not want to save
            if element_type not in element_types_to_save:
                continue

            plot_elements = canvas_instance.plot_elements[base_plot_type][
                plot_element_varname
            ][data_label][element_type]

            # for FAIRMODE target plot combine all dots (saved individually so that they can have different colors) in one Line2D
            if base_plot_type == "fairmode-target":
                x_array = []
                y_array = []

                for dot in plot_elements:
                    x_array.append(dot.get_xdata().tolist()[0])
                    y_array.append(dot.get_ydata().tolist()[0])

                plot_elements = [matplotlib.lines.Line2D(xdata=x_array, ydata=y_array)]

            for plot_element_i, plot_element in enumerate(plot_elements):
                if base_plot_type in [
                    "statsummary",
                    "heatmap",
                    "table",
                    "contingencytable",
                ]:
                    data = []

                    if base_plot_type in ["statsummary", "table", "contingencytable"]:
                        items = plot_element.get_celld().items()

                    elif base_plot_type == "heatmap":
                        items = np.ndenumerate(
                            plot_element.get_children()[0].get_array()
                        )

                    for (x, y), value in items:
                        data.append(
                            {
                                "x": x,
                                "y": y,
                                "z": value.get_text().get_text()
                                if base_plot_type != "heatmap"
                                else value,
                            }
                        )

                    df = (
                        pd.DataFrame(data, columns=["x", "y", "z"])
                        .pivot(index="x", columns="y", values="z")
                        .sort_index()
                        .rename_axis(index=None, columns=None)
                    )

                    df.columns = df.columns.astype(str)

                    # contingency plots have an index that starts at -1, reindex to start from 0
                    if base_plot_type == "contingencytable":
                        df.index = range(len(df))

                    filename = f"{plot_type}_{element_type}"

                    msgs = handle_test_or_save_df(
                        read_instance,
                        df,
                        filename,
                        path,
                        tests_generate_output,
                        msgs,
                        decimal_places,
                    )

                elif base_plot_type in [
                    "timeseries",
                    "distribution",
                    "histogram",
                    "scatter",
                    "fairmode-target",
                    "fairmode-statsummary",
                    "taylor",
                    "boxplot",
                    "periodic",
                    "periodic-violin",
                ]:
                    # skip periodic violin shapes
                    if isinstance(
                        plot_element, matplotlib.collections.FillBetweenPolyCollection
                    ):
                        continue

                    # extract annotations
                    if isinstance(plot_element, matplotlib.offsetbox.AnchoredOffsetbox):
                        data = []

                        for annotation in plot_element.get_child().get_children():
                            data.append(
                                {
                                    "dataset": annotation.get_text()
                                    .split("|")[0]
                                    .strip(),
                                    "annotation": annotation.get_text()
                                    .split("|")[1]
                                    .strip(),
                                }
                            )
                        df = pd.DataFrame(data)

                        filename = f"{plot_type}_{data_label}_{element_type}" + (
                            f"_{plot_element_i}" if len(plot_elements) > 1 else ""
                        )

                        msgs = handle_test_or_save_df(
                            read_instance,
                            df,
                            filename,
                            path,
                            tests_generate_output,
                            msgs,
                            decimal_places,
                        )

                    else:
                        if base_plot_type == "boxplot":
                            # skip patch
                            if isinstance(plot_element, matplotlib.patches.PathPatch):
                                continue

                            elif isinstance(plot_element, matplotlib.lines.Line2D):
                                y_value = plot_element.get_ydata()[0]

                            if data_label not in boxplot_accumulator:
                                boxplot_accumulator[data_label] = []
                            boxplot_accumulator[data_label].append(y_value)
                        else:
                            data = []
                            xy = plot_element.get_xydata()

                            # collapse the n+1 bin edges into centres, one
                            # per real bin, dropping the repeated closing edge
                            if base_plot_type == "histogram":
                                edges = xy[:, 0]
                                densities = xy[:-1, 1]
                                xy = np.column_stack(
                                    ((edges[:-1] + edges[1:]) / 2.0, densities)
                                )

                            for x, y in xy:
                                # no y axis on FAIRMODE statsummary
                                if base_plot_type == "fairmode-statsummary":
                                    data.append(
                                        {
                                            "value": x,
                                        }
                                    )
                                else:
                                    data.append(
                                        {
                                            # taylor's x is arccos(correlation) - undo it
                                            x_column: pd.to_datetime(
                                                x, unit="D", utc=True
                                            ).round("s")
                                            if base_plot_type == "timeseries"
                                            else math.cos(x)
                                            if base_plot_type == "taylor"
                                            else x,
                                            canvas_instance.plot_characteristics[
                                                plot_type
                                            ]["ylabel"]["ylabel"]
                                            if base_plot_type == "fairmode-target"
                                            else data_label: y,
                                        }
                                    )
                            df = pd.DataFrame(data)

                            # combine dataframes for some plots
                            if base_plot_type in [
                                "timeseries",
                                "scatter",
                                "distribution",
                                "histogram",
                                "periodic",
                                "periodic-violin",
                                "taylor",
                            ]:
                                # one dataframe per plot element
                                key = (element_type, plot_element_i)
                                df = df.set_index(x_column)
                                value_column = df.columns[0]

                                # column becomes the data label
                                df = df.rename(columns={value_column: data_label})

                                if key not in combined_dfs:
                                    combined_dfs[key] = df

                                else:
                                    combined_dfs[key] = pd.concat(
                                        [combined_dfs[key], df], axis=1
                                    )

                            # for other plot types save data per data label
                            else:
                                if base_plot_type == "fairmode-statsummary":
                                    plot_element_str = canvas_instance.plotting.fairmode_statsummary_row_titles[
                                        data_label
                                    ][
                                        plot_element_i
                                    ]
                                else:
                                    plot_element_str = (
                                        f"{plot_element_i}"
                                        if len(plot_elements) > 1
                                        else ""
                                    )
                                filename = f"{plot_type}_{data_label}_{element_type}_{plot_element_str}"
                                msgs = handle_test_or_save_df(
                                    read_instance,
                                    df,
                                    filename,
                                    path,
                                    tests_generate_output,
                                    msgs,
                                    decimal_places,
                                )

                elif base_plot_type == "metadata":
                    text = plot_element.get_text().split("\n")
                    column_names = []
                    values = []
                    for line in text:
                        if ":" in line:
                            title = line.split(":")[0]
                            value = " ".join(line.split(":")[1:]).strip()
                            column_names.append(title)
                            values.append(value)
                        else:
                            column_names.append("Stations")
                            values.append(line.split(" Stations")[0])

                    df = pd.DataFrame(columns=column_names)
                    df.loc[0] = values
                    df.columns = [str(col) for col in df.columns]
                    df = df.T.reset_index().rename(columns={"index": "key", 0: "value"})
                    filename = f"{plot_type}_{data_label}_{element_type}" + (
                        f"_{plot_element_i}" if len(plot_elements) > 1 else ""
                    )
                    msgs = handle_test_or_save_df(
                        read_instance,
                        df,
                        filename,
                        path,
                        tests_generate_output,
                        msgs,
                        decimal_places,
                    )

                elif base_plot_type == "map":
                    # extract annotations
                    if isinstance(plot_element, matplotlib.offsetbox.AnchoredOffsetbox):
                        data = []
                        for annotation in plot_element.get_child().get_children():
                            data.append(
                                {
                                    "dataset": annotation.get_text()
                                    .split("|")[0]
                                    .strip(),
                                    "annotation": annotation.get_text()
                                    .split("|")[1]
                                    .strip(),
                                }
                            )
                        filename = f"{plot_type}_{data_label}_{element_type}" + (
                            f"_{plot_element_i}" if len(plot_elements) > 1 else ""
                        )
                    # extract plot data
                    else:
                        coordinates = plot_element.get_offsets()
                        values = plot_element.get_array()

                        # extract data
                        data = []
                        for (lon, lat), val in zip(coordinates, values):
                            lon_ind = np.where(
                                read_instance.station_longitudes[networkspeci] == lon
                            )
                            lat_ind = np.where(
                                read_instance.station_latitudes[networkspeci] == lat
                            )
                            intersect_ind = np.intersect1d(lon_ind, lat_ind)

                            # if there is only one common lon and lat pair
                            # use intersection index to search for station name and reference
                            if intersect_ind.size == 1:
                                station_ind = intersect_ind[0]
                            else:
                                station_ind = None

                            if station_ind is not None:
                                station_name = read_instance.station_names[
                                    networkspeci
                                ][station_ind]
                                station_reference = read_instance.station_references[
                                    networkspeci
                                ][station_ind]
                            else:
                                station_name = np.nan
                                station_reference = np.nan

                            data.append(
                                {
                                    "Station name": station_name,
                                    "Station reference": station_reference,
                                    "Longitude": lat,
                                    "Latitude": lon,
                                    "Statistic": val,
                                }
                            )
                        label = ""
                        if (labela != "") and (labelb == ""):
                            label = labela
                        elif (labelb != "") and (labela == ""):
                            label = labelb
                        elif (labela != "") and (labelb != ""):
                            label = f"{labela}-{labelb}"
                        filename = f"{plot_type}_{element_type}_{label}"
                    df = pd.DataFrame(data)
                    msgs = handle_test_or_save_df(
                        read_instance,
                        df,
                        filename,
                        path,
                        tests_generate_output,
                        msgs,
                        decimal_places,
                    )

    # save combined dataframes into one file per plot element
    if base_plot_type in [
        "timeseries",
        "scatter",
        "distribution",
        "histogram",
        "periodic",
        "periodic-violin",
        "taylor",
        "boxplot",
    ]:
        if base_plot_type == "boxplot":
            stats = ["whisker_low", "q1", "median", "q3", "whisker_high"]
            data = {}

            for label, stats_list in boxplot_accumulator.items():
                stats_list_sorted = sorted(stats_list)
                data[label] = dict(zip(stats, stats_list_sorted))

            df = pd.DataFrame(data)
            df = df.reset_index()
            filename = "boxplot"

            msgs = handle_test_or_save_df(
                read_instance,
                df,
                filename,
                path,
                tests_generate_output,
                msgs,
                decimal_places,
            )

        else:
            for (element_type, plot_element_i), df in combined_dfs.items():
                df = df.reset_index()
                filename = f"{plot_type}_{element_type}" + (
                    f"_{plot_element_i}" if len(plot_elements) > 1 else ""
                )

                msgs = handle_test_or_save_df(
                    read_instance,
                    df,
                    filename,
                    path,
                    tests_generate_output,
                    msgs,
                    decimal_places,
                )

    if msgs:
        msg = f"Saving {plot_type} figure data to CSV:"
        msg += "".join(msgs)
        show_message(read_instance, msg)
