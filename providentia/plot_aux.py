""" Auxiliar functions to create plots """

import copy
import sys
import os
import math
import json
import scipy
import yaml

import matplotlib as mpl
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
import numpy as np
import pandas as pd
from pypdf import PdfReader, PdfWriter
from scipy.signal import convolve
from scipy.signal.windows import gaussian
from scipy.sparse import coo_matrix
import seaborn as sns
from .statistics import calculate_statistic, get_z_statistic_sign

from providentia.auxiliar import CURRENT_PATH

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

def get_multispecies_aliases(networkspecies):
    """ Map networkspecies to networkspecies aliases.
        Also get label for alias.   
    """

    multispecies_aliases = {'vconcaerobin1': '0.05',
                            'vconcaerobin2': '0.066',
                            'vconcaerobin3': '0.086',
                            'vconcaerobin4': '0.113',
                            'vconcaerobin5': '0.148',
                            'vconcaerobin6': '0.194',
                            'vconcaerobin7': '0.255',
                            'vconcaerobin8': '0.335',
                            'vconcaerobin9': '0.439',
                            'vconcaerobin10': '0.576',
                            'vconcaerobin11': '0.756',
                            'vconcaerobin12': '0.992',
                            'vconcaerobin13': '1.302',
                            'vconcaerobin14': '1.708',
                            'vconcaerobin15': '2.241',
                            'vconcaerobin16': '2.940',
                            'vconcaerobin17': '3.857',
                            'vconcaerobin18': '5.061',
                            'vconcaerobin19': '6.641',
                            'vconcaerobin20': '8.713',
                            'vconcaerobin21': '11.432',
                            'vconcaerobin22': '15.00'
                            }

    multispecies_labels =  {'vconcaerobin1': 'Radius [µm]',
                            'vconcaerobin2': 'Radius [µm]',
                            'vconcaerobin3': 'Radius [µm]',
                            'vconcaerobin4': 'Radius [µm]',
                            'vconcaerobin5': 'Radius [µm]',
                            'vconcaerobin6': 'Radius [µm]',
                            'vconcaerobin7': 'Radius [µm]',
                            'vconcaerobin8': 'Radius [µm]',
                            'vconcaerobin9': 'Radius [µm]',
                            'vconcaerobin10': 'Radius [µm]',
                            'vconcaerobin11': 'Radius [µm]',
                            'vconcaerobin12': 'Radius [µm]',
                            'vconcaerobin13': 'Radius [µm]',
                            'vconcaerobin14': 'Radius [µm]',
                            'vconcaerobin15': 'Radius [µm]',
                            'vconcaerobin16': 'Radius [µm]',
                            'vconcaerobin17': 'Radius [µm]',
                            'vconcaerobin18': 'Radius [µm]',
                            'vconcaerobin19': 'Radius [µm]',
                            'vconcaerobin20': 'Radius [µm]',
                            'vconcaerobin21': 'Radius [µm]',
                            'vconcaerobin22': 'Radius [µm]'
                            }
    
    networkspecies_aliases = [multispecies_aliases[networkspeci] 
                              if networkspeci in multispecies_aliases else networkspeci 
                              for networkspeci in networkspecies]

    labels = np.unique([multispecies_labels[networkspeci] 
                        for networkspeci in networkspecies if networkspeci in multispecies_labels])
    if len(labels) == 1:
        unique_label = labels[0]
    else:
        unique_label = ''

    return networkspecies_aliases, unique_label


def temp_axis_dict():
    """ Return temporal mapping as a dictionary used for the plots.

        :return: Numbering of months/days
        :rtype: dict
    """

    map_dict = {'short': {'dayofweek': {0: 'M', 1: 'T', 2: 'W', 3: 'T', 4: 'F', 5: 'S', 6: 'S'},
                          'month': {1: 'J', 2: 'F', 3: 'M', 4: 'A', 5: 'M', 6: 'J',
                                    7: 'J', 8: 'A', 9: 'S', 10: 'O', 11: 'N', 12: 'D'}},
                'long': {'dayofweek': {0: 'Monday', 1: 'Tuesday', 2: 'Wednesday', 3: 'Thursday', 
                                       4: 'Friday', 5: 'Saturday', 6: 'Sunday'},
                         'month': {1: 'January', 2: 'February', 3: 'March', 4: 'April', 5: 'May', 6: 'June',
                                   7: 'July', 8: 'August', 9: 'September', 10: 'October', 11: 'November', 
                                   12: 'December'}}
                }

    return map_dict


def periodic_xticks():
    """ Return xticks for periodic subplots.

        :return: Dictionary of xticks per temporal resolution
        :rtype: dict
    """

    return {'hour': np.arange(24, dtype=np.int64), 
            'dayofweek': np.arange(7, dtype=np.int64), 
            'month': np.arange(1, 13, dtype=np.int64)}


def periodic_labels():
    """ Return axes labels for periodic subplots.

        :return: Axes labels
        :rtype: dict
    """

    return {'hour':'H', 'dayofweek':'DoW', 'month':'M'}


def get_land_polygon_resolution(selection):
    """ Get resolution of land polygons to plot on map.

        :param selection: Selected temporal resolution
        :type resolution: str
        :return: Selected land polygon resolution
        :rtype: list
    """
    
    land_polygon_resolutions = {'low': '110m','medium': '50m','high': '10m'}
    
    return land_polygon_resolutions[selection]


def update_plotting_parameters(instance):
    """ Function that updates plotting parameters (colour and zorder) for data labels.

        :param instance: Instance of class ProvidentiaOffline or ProvidentiaMainWindow
        :type instance: object
    """

    # generate a list of RGB tuples for number of experiments there are
    sns.reset_orig()
    color_palette = instance.plot_characteristics_templates['general']['legend_color_palette']
    color_palettes = yaml.safe_load(open(os.path.join(PROVIDENTIA_ROOT, 'settings/color_palettes.yaml')))
    if color_palette in color_palettes.keys():
        if (len(instance.data_labels) - 1) > len(color_palettes[color_palette]):
            msg = "Error: The number of experiments and palette colors should be equal. "
            msg += f"Add more colors to your palette '{color_palette}' in settings/color_palettes.yaml "
            msg += "or change your legend_color_palette in the plot characteristics files."
            sys.exit()
        else:
            clrs = sns.color_palette(color_palettes[color_palette])
    else:
        clrs = sns.color_palette(color_palette, n_colors=len(instance.data_labels)-1)

    # add colour and zorder for observations
    instance.plotting_params[instance.observations_data_label]['colour'] = instance.plot_characteristics_templates['general']['obs_markerfacecolor']
    instance.plotting_params[instance.observations_data_label]['zorder'] = instance.plot_characteristics_templates['general']['obs_zorder']

    # add colours and zorder for each experiment
    experiment_ind = 1
    for data_label in instance.data_labels:
        if data_label != instance.observations_data_label:
            # define colour for experiment
            instance.plotting_params[data_label]['colour'] = clrs[experiment_ind-1]
            # define zorder for experiment (obs zorder + experiment_ind)
            instance.plotting_params[data_label]['zorder'] = \
                instance.plotting_params[instance.observations_data_label]['zorder'] + experiment_ind
            # update count of experiments
            experiment_ind += 1


def kde_fft(xin, gridsize=1024, extents=None, weights=None, adjust=1., bw='scott', xgrid=None):
    """
    A fft-based Gaussian kernel density estimate (KDE)
    for computing the KDE on a regular grid

    Note that this is a different use case than scipy's original
    scipy.stats.kde.gaussian_kde

    IMPLEMENTATION
    --------------

    Performs a gaussian kernel density estimate over a regular grid using a
    convolution of the gaussian kernel with a histogram of the data.

    It computes the sparse  histogram of the data samples where
    *x* is a 1-D sequence of the same length. If *weights* is None
    (default), this is a histogram of the number of occurences of the
    observations at (x[i]).
    histogram of the data is a faster implementation than numpy.histogram as it
    avoids intermediate copies and excessive memory usage!

    This function is typically *several orders of magnitude faster* than
    scipy.stats.kde.gaussian_kde.  For large (>1e7) numbers of points, it
    produces an essentially identical result.

    Boundary conditions on the data is corrected by using a symmetric /
    reflection condition. Hence the limits of the dataset does not affect the
    pdf estimate.

    INPUTS
    ------

        xin:  ndarray[ndim=1]
            The 1D data samples

        gridsize: int
            A nx integer of the size of the output grid (default: 1024)

        extents: (xmin, xmax) tuple
            tuple of the extents of output grid (default: extent of input data)

        weights: ndarray[ndim=1]
            An array of the same shape as x that weights each sample x_i
            by w_i.  Defaults to an array of ones the same size as x (default: None)

        adjust : float
            An adjustment factor for the bw. Bandwidth becomes bw * adjust.

        bw: str
            The method used to calculate the estimator bandwidth. 
            This can be 'scott' or 'silverman'(default: 'scott')

        xgrid: ndarray(ndim=1)
            The output grid (if this is provided gridsize and extents are ignored).

    OUTPUTS
    -------
        grid_points: ndarray[ndim=1]
            Grid points to sample density estimates

        kde: ndarray[ndim=1]
            A gridded 1D kernel density estimate of input data samples at grid points

    """
    # Variable check
    x = np.squeeze(np.asarray(xin))
    
    # Default extents are the extent of the data
    if xgrid is not None:
        xmin, xmax = xgrid.min(), xgrid.max()
        x = x[ (x <= xmax) & (x >= xmin) ]
    elif extents is None:
        xmin, xmax = x.min(), x.max()
    else:
        xmin, xmax = map(float, extents)
        x = x[ (x <= xmax) & (x >= xmin) ]

    n = x.size

    # apply weights
    if weights is None:
        # Default: Weight all points equally
        weights = np.ones(n)
    else:
        weights = np.squeeze(np.asarray(weights))
        if weights.size != x.size:
            raise ValueError('Input weights must be an array of the same size as xin!')

    # Optimize gridsize ------------------------------------------------------
    # Make grid and discretize the data and round it to the next power of 2
    # to optimize with the fft usage
    # ensure minimum gridsize of 1024 points
    if xgrid is None:
        if gridsize is None:
            gridsize = np.max((len(x), 1024.))
        gridsize = 2 ** np.ceil(np.log2(gridsize))  # round to next power of 2
        nx = int(gridsize)
    else:
        nx = len(xgrid)

    # Make the sparse histogram -------------------------------------------
    dx = (xmax - xmin) / (nx - 1)

    # Basically, this is just doing what np.digitize does with one less copy
    xyi = x - xmin
    xyi /= dx
    xyi = np.floor(xyi, xyi)
    xyi = np.vstack((xyi, np.zeros(n, dtype=int)))

    # Next, make a histogram of x
    # Exploit a sparse coo_matrix avoiding np.histogram due to excessive
    # memory usage with many points
    grid = coo_matrix((weights, xyi), shape=(nx, 1)).toarray()

    # Kernel Preliminary Calculations ---------------------------------------
    std_x = np.std(xyi[0])

    # Scaling factor for bandwidth
    if bw == 'scott': 
        bw_factor = (n ** (-1. / 5.)) * adjust
    elif bw == 'silverman':
        bw_factor =  ((n * 3 / 4.)**(-1. / 5)) * adjust 

    # Make the gaussian kernel ---------------------------------------------
    # First, determine the bandwidth using defined bandwidth estimator rule
    kern_nx = int(np.round(bw_factor * 2 * np.pi * std_x))

    # If bandwidth is 0, skip plot for current data label
    if kern_nx == 0:
        return None
    
    # Then evaluate the gaussian function on the kernel grid
    kernel = np.reshape(gaussian(kern_nx, bw_factor * std_x), (kern_nx, 1))

    #---- Produce the kernel density estimate --------------------------------

    # Convolve the histogram with the gaussian kernel
    # use symmetric padding to correct for data boundaries in the kde
    npad = np.min((nx, 2 * kern_nx))
    grid = np.vstack( [grid[npad: 0: -1], grid, grid[nx: nx - npad: -1]] )
    grid = convolve(grid, kernel, mode='same')[npad: npad + nx]

    # Normalization factor to divide result by so that units are in the same
    # units as scipy.stats.kde.gaussian_kde's output.
    norm_factor = 2 * np.pi * std_x * std_x * bw_factor ** 2
    norm_factor = n * dx * np.sqrt(norm_factor)

    # Normalize the result
    grid /= norm_factor

    # return grid points and estimated densities 
    if xgrid is None:
        return np.linspace(xmin,xmax,nx), np.squeeze(grid)
    else:
        return np.squeeze(grid)


def round_decimal_places(x, decimal_places):
    """ Round x to decimal places

    :param x: Value
    :type x: float
    :param decimal_places: Desired number of decimal places
    :type decimal_places: int
    :return: Rounded value
    :rtype: str
    """

    # if cell value is not nan
    if not np.isnan(x):
        # if number is zero, return int as str
        if x == 0:
            return '0'
        # if number of zeros is more than decimal places set by user, use scientific notation
        elif ((-math.floor(math.log10(abs(x))) - 1) > decimal_places):
            return '{:0.{}e}'.format(x, decimal_places) 
        # if not, use float
        else:
            return '{:0.{}f}'.format(x, decimal_places)
    # if cell value is nan, do nothing
    else:
        return 'nan'


def merge_cells(table, cells):
    """
    Merge cells in matplotlib.Table. Reference: https://stackoverflow.com/a/53819765/12684122.

    :param table: Table
    :type table: matplotlib.Table
    :param cells: Cells to merge in table coordinates, e.g. [(0,1), (0,0), (0,2)]
    :type cells: list
    """

    cells_array = [np.asarray(c) for c in cells]
    h = np.array([cells_array[i+1][0] - cells_array[i][0] for i in range(len(cells_array) - 1)])
    v = np.array([cells_array[i+1][1] - cells_array[i][1] for i in range(len(cells_array) - 1)])

    # if it's a horizontal merge, all values for `h` are 0
    if not np.any(h):
        # sort by horizontal coord
        cells = np.array(sorted(list(cells), key=lambda v: v[1]))
        edges = ['BTL'] + ['BT' for i in range(len(cells) - 2)] + ['BTR']
    elif not np.any(v):
        cells = np.array(sorted(list(cells), key=lambda h: h[0]))
        edges = ['TRL'] + ['RL' for i in range(len(cells) - 2)] + ['BRL']
    else:
        raise ValueError("Only horizontal and vertical merges allowed")

    for cell, e in zip(cells, edges):
        table[cell[0], cell[1]].visible_edges = e
        
    txts = [table[cell[0], cell[1]].get_text() for cell in cells]
    tpos = [np.array(t.get_position()) for t in txts]

    # transpose the text of the left cell
    trans = (tpos[-1] - tpos[0])/2

    # didn't had to check for ha because I only want ha='center'
    txts[0].set_transform(mpl.transforms.Affine2D().translate(*trans))
    for txt in txts[1:]:
        txt.set_visible(False)


def get_taylor_diagram_ghelper_info(reference_stddev, plot_characteristics, extend=False):
    """ Make Taylor diagram plot axis extremes and labels. 
        
        :param reference_stddev: Reference standard deviation to set the limits
        :type reference_stddev: float
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param extend: Indicates if plots needs to show negative correlations
        :type undo: boolean
    """

    tmin = 0

    # diagram extended to negative correlations
    if extend:
        tmax = np.pi
    # diagram limited to positive correlations
    else:
        tmax = np.pi/2

    # get standard deviation axis extent
    srange = plot_characteristics['srange']
    smin = srange[0] * reference_stddev
    smax = srange[1] * reference_stddev

    # correlation labels
    if extend:
        rlocs = np.array(plot_characteristics['rlocs']['all'])
        rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
    else:
        rlocs = np.array(plot_characteristics['rlocs']['positive'])

    # convert correlation values into polar angles
    tlocs = np.arccos(rlocs)
    gl1 = gf.FixedLocator(tlocs)
    tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    return tmin, tmax, smin, smax, gl1, tf1


def get_taylor_diagram_ghelper(reference_stddev, plot_characteristics, extend=False):
    """ Make Taylor diagram plot grid helper. 

        :param reference_stddev: Reference standard deviation to set the limits
        :type reference_stddev: float
        :param plot_characteristics: Plot characteristics  
        :type plot_characteristics: dict
        :param extend: Indicates if plots needs to show negative correlations
        :type undo: boolean
    """

    # get axis extremes
    tmin, tmax, smin, smax, gl1, tf1 = get_taylor_diagram_ghelper_info(reference_stddev, 
                                                                       plot_characteristics, 
                                                                       extend)

    # get grid helper
    ghelper = fa.GridHelperCurveLinear(PolarAxes.PolarTransform(apply_theta_transforms=False),
                                       extremes=(tmin, tmax, smin, smax),
                                       grid_locator1=gl1, tick_formatter1=tf1,)

    return ghelper


def set_map_extent(canvas_instance, ax, map_extent):
    """ Set map extent, done set_xlim and set_ylim rather than set_extent 
        to avoid axis cutting off slightly (https://github.com/SciTools/cartopy/issues/697).
    """

    mlon = np.mean(map_extent[:2])
    mlat = np.mean(map_extent[2:])
    xtrm_data = np.array([[map_extent[0], mlat], [mlon, map_extent[2]], 
                          [map_extent[1], mlat], [mlon, map_extent[3]]])
    proj_to_data = canvas_instance.datacrs._as_mpl_transform(ax) - ax.transData
    xtrm = proj_to_data.transform(xtrm_data)
    ax.set_xlim(xtrm[:,0].min(), xtrm[:,0].max())
    ax.set_ylim(xtrm[:,1].min(), xtrm[:,1].max())


def get_map_extent(canvas_instance):
    """ Get map extent from xlim and ylim. """ 

    # get plot extent
    ax = canvas_instance.plot_axes['map']
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
    transformed_coords = canvas_instance.datacrs.transform_points(canvas_instance.plotcrs, 
                                                                  xcoords, ycoords)[:, :2]

    # keep longitudes between -180 and 180
    lon_change = False
    if (np.isnan(transformed_coords[0, 0])) or (transformed_coords[0, 0] == -179.99999999999932):
        transformed_coords[0, 0] = -180
        lon_change = True
    if (np.isnan(transformed_coords[2, 0])) or (transformed_coords[2, 0] == 179.99999999999932):
        transformed_coords[2, 0] = 180  
        lon_change = True 

    # keep latitudes between -90 and 90
    lat_change = False
    if (np.isnan(transformed_coords[1, 1])) or (transformed_coords[1, 1] == -89.99999999999966):
        transformed_coords[1, 1] = -90
        lat_change = True  
    if (np.isnan(transformed_coords[3, 1])) or (transformed_coords[3, 1] == 89.99999999999966):
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
    map_extent = [transformed_coords[:,0].min(), transformed_coords[:,0].max(),
                  transformed_coords[:,1].min(), transformed_coords[:,1].max()]

    return map_extent


def create_chunked_timeseries(read_instance, canvas_instance, chunk_stat, chunk_resolution, 
                              networkspeci, cut_data_labels, bias):
    """ Create statistical timeseries data by chunk resolution

    :param read_instance: Instance of class ProvidentiaMainWindow or ProvidentiaOffline
    :type read_instance: object
    :param canvas_instance: Instance of class MPLCanvas or ProvidentiaOffline
    :type canvas_instance: object
    :param chunk_stat: Chunk statistic
    :type chunk_stat: str
    :param chunk_resolution: Chunk resolution
    :type chunk_resolution: str
    :param networkspeci: Current networkspeci (e.g. EBAS|sconco3) 
    :type networkspeci: str
    :param cut_data_labels: Valid data labels
    :type cut_data_labels: list
    :return: Dataframe with statistics per day / month / year
    :rtype: pandas dataframe
    """
    
    z_statistic_sign = get_z_statistic_sign(chunk_stat)
    if z_statistic_sign == 'bias' or bias:
        if read_instance.observations_data_label in cut_data_labels:
            cut_data_labels.remove(read_instance.observations_data_label)
        stats_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, chunk_stat, 
                                         [read_instance.observations_data_label]*len(cut_data_labels), 
                                         cut_data_labels, chunking=True, chunk_resolution=chunk_resolution)
        
    else:
        stats_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, 
                                         chunk_stat, cut_data_labels, [], chunking=True, 
                                         chunk_resolution=chunk_resolution)

    chunk_dates = canvas_instance.selected_station_data[networkspeci]["timeseries_chunks"][chunk_resolution]['valid_xticks']
    timeseries_data = pd.DataFrame(index=chunk_dates, columns=cut_data_labels, dtype=np.float64)

    for chunk_date_idx, chunk_date in enumerate(chunk_dates):
        for label_idx, data_label in enumerate(cut_data_labels):
            timeseries_data.loc[chunk_date, data_label] = stats_calc[chunk_date_idx][label_idx]
    
    return timeseries_data


def reorder_pdf_pages(input_pdf, output_pdf, summary_multispecies_pages, 
                      station_multispecies_pages, paradigm_break_page):
    """ Reorder PDF pages so that multispecies plots appear before other plots.

    :param input_pdf: Path to original PDF
    :type input_pdf: string
    :param output_pdf: Path to ordered PDF
    :type output_pdf: string
    :param summary_multispecies_pages: Pages that contain summary multispecies plots
    :type summary_multispecies_pages: list
    :param station_multispecies_pages: Pages that contain station multispecies plots
    :type station_multispecies_pages: list
    :param paradigm_break_page: Page where station plot start
    :type paradigm_break_page: int
    """

    # Get pages
    summary_multispecies_pages = np.array(summary_multispecies_pages)
    station_multispecies_pages = np.array(station_multispecies_pages)

    # Get original order
    input_pdf_file = PdfReader(open(input_pdf, "rb"))
    all_pages = np.arange(len(input_pdf_file.pages))

    # Initialise page order
    page_order = copy.deepcopy(all_pages)

    # Move summary pages after page 0 (cover)
    if len(summary_multispecies_pages) > 0:
        summary_multispecies_pages = np.concatenate((np.array([0]), summary_multispecies_pages))
        summary_other_plots_pages = all_pages[~np.isin(all_pages, summary_multispecies_pages)]
        page_order = np.concatenate((
            summary_multispecies_pages, 
            summary_other_plots_pages)).tolist()
    
    # Move station pages after paradigm break page (when we start to see station plots)
    if len(station_multispecies_pages) > 0:
        station_pages = all_pages[paradigm_break_page:]
        station_other_plots_pages = station_pages[~np.isin(station_pages, station_multispecies_pages)]
        page_order = np.concatenate((
            page_order[:paradigm_break_page], 
            station_multispecies_pages, 
            station_other_plots_pages)).tolist()

    # Reorder pages
    output_pdf_file = PdfWriter()
    for page_number in page_order:
        output_pdf_file.add_page(input_pdf_file.pages[page_number])

    # Write the rearranged pages to a new PDF file
    print(f'Writing {output_pdf}')
    with open(output_pdf, "wb") as outputStream:
        output_pdf_file.write(outputStream)
