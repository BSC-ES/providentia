""" Auxiliary plotting functions """

import copy
import math
import os
import sys

from matplotlib import transforms
from matplotlib.colors import cnames
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
from netCDF4 import Dataset
import numpy as np
import pandas as pd
from scipy.signal import convolve
from scipy.signal.windows import gaussian
from scipy.sparse import coo_matrix
import seaborn as sns
import yaml

from providentia.auxiliar import CURRENT_PATH, join, get_conversion_factor, get_standard_parameters_by_speci
from .statistics import calculate_statistic, get_z_statistic_sign

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])

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

    radius_per_species = {'vconcaerobin1': '0.05',
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
    
    networkspecies_aliases = [get_AERONET_sizedist_bin_radius(networkspeci) 
                              if networkspeci in multispecies_labels else networkspeci 
                              for networkspeci in networkspecies]

    labels = np.unique([multispecies_labels[networkspeci] 
                        for networkspeci in networkspecies if networkspeci in multispecies_labels])
    if len(labels) == 1:
        unique_label = labels[0]
    else:
        unique_label = ''

    return networkspecies_aliases, unique_label

def temp_axis_dict():
    """ 
    Return temporal mapping as a dictionary used for the plots.
        
    Returns
    -------
    map_dict : dict
        Numbering of months/days
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
    """
    Return xticks for periodic subplots based on temporal resolution.

    Returns
    -------
    xticks : dict
        Dictionary of xticks per temporal resolution
    """
    
    xticks = {'hour': np.arange(24, dtype=np.int8), 
            'dayofweek': np.arange(7, dtype=np.int8), 
            'month': np.arange(1, 13, dtype=np.int8)}

    return xticks


def periodic_labels():
    """
    Return axes labels for periodic subplots.

    Returns
    -------
    axes_labels : dict
        Axes labels per temporal resolution.
    """

    axes_labels = {'hour':'H', 'dayofweek':'DoW', 'month':'M'}

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
    
    land_polygon_resolutions = {'low': '110m','medium': '50m','high': '10m'}
    resolution = land_polygon_resolutions[selection]
    
    return resolution

def update_plotting_parameters(instance, data_labels_to_remove=None, data_labels_to_add=None, daily_forecast=False):
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
            [data_label.split('-day')[0].split('-daily')[0].split('-combined')[0] for data_label in data_labels_to_add]
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
                    if (not daily_forecast) & (unique_base_data_label in instance.files_to_read[valid_networkspeci]):
                        # Open netCDF file to extract grid edges
                        mod_nc_root = Dataset(instance.files_to_read[valid_networkspeci][unique_base_data_label][0])
                        grid_edge_longitude = mod_nc_root['grid_edge_longitude'][:]
                        grid_edge_latitude = mod_nc_root['grid_edge_latitude'][:]
                        # Close netCDF file
                        mod_nc_root.close() 

                    # iterate through data labels to add
                    for data_label in data_labels_to_add:

                        base_data_label = data_label.split('-day')[0].split('-daily')[0].split('-combined')[0]

                        if (base_data_label == unique_base_data_label) & (data_label not in processed_data_labels):

                            # initialise plotting params dict for data label to add
                            instance.plotting_params[data_label] = {}

                            # For daily forecast, copy grid edge info from removed labels
                            if daily_forecast:
                                
                                relevant_inds = [ii for ii, data_label_to_remove in enumerate(data_labels_to_remove) if data_label in data_label_to_remove]
                                instance.plotting_params[data_label]['grid_edge_longitude'] = instance.plotting_params[data_labels_to_remove[relevant_inds[0]]]['grid_edge_longitude'] 
                                instance.plotting_params[data_label]['grid_edge_latitude'] = instance.plotting_params[data_labels_to_remove[relevant_inds[0]]]['grid_edge_latitude'] 
                                processed_data_labels.append(data_label)
                            
                            # otherwise open netCDF file to extract grid edges
                            elif unique_base_data_label in instance.files_to_read[valid_networkspeci]:
                                
                                instance.plotting_params[data_label]['grid_edge_longitude'] = grid_edge_longitude
                                instance.plotting_params[data_label]['grid_edge_latitude'] = grid_edge_latitude
                                processed_data_labels.append(data_label)

    # Remove plotting parameters for labels marked for removal
    if data_labels_to_remove is not None:
        for data_label in data_labels_to_remove:
            del instance.plotting_params[data_label]

    # Add colour and zorder for observations
    instance.plotting_params[instance.observations_data_label]['colour'] = instance.plot_characteristics_templates['general']['obs_markerfacecolor']
    instance.plotting_params[instance.observations_data_label]['zorder'] = instance.plot_characteristics_templates['general']['obs_zorder']

    # Generate a list of RGB tuples for the number of models
    sns.reset_orig()  # Reset seaborn to default
    color_palette = instance.plot_characteristics_templates['general']['legend_color_palette']
    color_palettes = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/color_palettes.yaml')))

    if color_palette in color_palettes.keys():
        # Check that the number of colors matches the number of models
        if (len(instance.data_labels) - 1) > len(color_palettes[color_palette]):
            error = "Error: The number of models and palette colors should be equal. "
            error += f"Add more colors to your palette '{color_palette}' in settings/color_palettes.yaml "
            error += "or change your legend_color_palette in the plot characteristics files."
            instance.logger.error(error)
            sys.exit(1)
        else:
            clrs = sns.color_palette(color_palettes[color_palette])
    else:
        # If palette not in YAML, generate colors automatically
        clrs = sns.color_palette(color_palette, n_colors=len(instance.data_labels)-1)

    # Add colours and zorder for each model (non-observations)
    model_ind = 1
    for data_label in instance.data_labels:
        if data_label != instance.observations_data_label:
            # Define colour for model
            instance.plotting_params[data_label]['colour'] = clrs[model_ind-1]
            # Define zorder for model relative to observations
            instance.plotting_params[data_label]['zorder'] = \
                instance.plotting_params[instance.observations_data_label]['zorder'] + model_ind
            # Update count of models
            model_ind += 1

def kde_fft(xin, gridsize=1024, extents=None, weights=None, adjust=1., bw='scott', xgrid=None):
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
            error = 'Input weights must be an array of the same size as xin. '
            return error

    # make grid and discretize the data and round it to the next power of 2
    # to optimize with the fft usage
    # ensure minimum gridsize of 1024 points
    if xgrid is None:
        if gridsize is None:
            gridsize = np.max((len(x), 1024.))
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
        error = 'Too many zeros. '
        return error

    # Kernel Preliminary Calculations 
    std_x = np.std(xyi[0])

    # Scaling factor for bandwidth
    if bw == 'scott': 
        bw_factor = (n ** (-1. / 5.)) * adjust
    elif bw == 'silverman':
        bw_factor =  ((n * 3 / 4.)**(-1. / 5)) * adjust 

    # make the gaussian kernel
    # first, determine the bandwidth using defined bandwidth estimator rule
    kern_nx = int(np.round(bw_factor * 2 * np.pi * std_x))

    # If bandwidth is 0, skip plot for current data label
    if kern_nx == 0:
        error = 'The kernel bandwidth is 0.  '
        error += 'To change the bandwith, we recommend increasing the number of '
        error += 'pdf_min_samples in the plot characteristics settings files. '
        return error
    
    # Then evaluate the gaussian function on the kernel grid
    kernel = np.reshape(gaussian(kern_nx, bw_factor * std_x), (kern_nx, 1))

    # convolve the histogram with the gaussian kernel
    # use symmetric padding to correct for data boundaries in the kde
    npad = np.min((nx, 2 * kern_nx))
    grid = np.vstack( [grid[npad: 0: -1], grid, grid[nx: nx - npad: -1]] )
    grid = convolve(grid, kernel, mode='same')[npad: npad + nx]

    # normalization factor to divide result by so that units are in the same
    # units as scipy.stats.kde.gaussian_kde's output.
    norm_factor = 2 * np.pi * std_x * std_x * bw_factor ** 2
    norm_factor = n * dx * np.sqrt(norm_factor)

    # normalize the result
    grid /= norm_factor

    # return grid points and estimated densities 
    if xgrid is None:
        return np.linspace(xmin,xmax,nx), np.squeeze(grid)
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
    txts[0].set_transform(transforms.Affine2D().translate(*trans))
    if not visibility:
        for txt in txts[1:]:
            txt.set_visible(False)

def get_taylor_diagram_ghelper_info(reference_stddev, plot_characteristics, extend=False):
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
    tmin, tmax, smin, smax, gl1, tf1 = get_taylor_diagram_ghelper_info(reference_stddev, 
                                                                       plot_characteristics, 
                                                                       extend)

    # get grid helper
    ghelper = fa.GridHelperCurveLinear(PolarAxes.PolarTransform(apply_theta_transforms=False),
                                       extremes=(tmin, tmax, smin, smax),
                                       grid_locator1=gl1, tick_formatter1=tf1,)

    return ghelper

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
    xtrm_data = np.array([[map_extent[0], mlat], [mlon, map_extent[2]], 
                          [map_extent[1], mlat], [mlon, map_extent[3]]])
    proj_to_data = canvas_instance.datacrs._as_mpl_transform(ax) - ax.transData
    xtrm = proj_to_data.transform(xtrm_data)
    ax.set_xlim(xtrm[:,0].min(), xtrm[:,0].max())
    ax.set_ylim(xtrm[:,1].min(), xtrm[:,1].max())

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

def create_statistical_timeseries(read_instance, canvas_instance, chunk_stat, chunk_resolution, 
                                  networkspeci, cut_data_labels, bias):
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
        forecast_type = 'daily'
    elif read_instance.combined_forecast:
        forecast_type = 'combined'

    if (z_statistic_sign == 'bias') or (bias):
        if read_instance.observations_data_label in cut_data_labels:
            cut_data_labels.remove(read_instance.observations_data_label)
        stats_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, chunk_stat, 
                                         [read_instance.observations_data_label]*len(cut_data_labels), 
                                         cut_data_labels, chunk_resolution=chunk_resolution, 
                                         statistic_aggregation=read_instance.timeseries_statistic_aggregation,
                                         forecast_type=forecast_type)
    else:
        stats_calc = calculate_statistic(read_instance, canvas_instance, networkspeci, 
                                         chunk_stat, cut_data_labels, [], chunk_resolution=chunk_resolution,
                                         statistic_aggregation=read_instance.timeseries_statistic_aggregation,
                                         forecast_type=forecast_type) 

    if (chunk_resolution == read_instance.active_resolution) & (forecast_type == 'combined'):
        chunk_dates = read_instance.time_index_ts
    elif (chunk_resolution == read_instance.active_resolution) & (forecast_type != 'daily'):
        chunk_dates = read_instance.time_index
    else:
        chunk_dates = canvas_instance.grouped_ts_index

    timeseries_data = pd.DataFrame(index=chunk_dates, columns=cut_data_labels, dtype=np.float32)

    # if shape of stats_calc is not correct return error
    if (stats_calc.shape[0] != len(chunk_dates)) or (stats_calc.shape[1] != len(cut_data_labels)):
        error = "Error: The shape of the calculated statistical timseseries does not match the expected shape."
        read_instance.logger.error(error)
        sys.exit(1)

    for chunk_date_idx, chunk_date in enumerate(chunk_dates):
        for label_idx, data_label in enumerate(cut_data_labels):
            timeseries_data.loc[chunk_date, data_label] = np.float32(stats_calc[chunk_date_idx,label_idx])
    
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
    if type(colour) == str:
        if colour[0] != '#':
            hex_colour = cnames[colour]
    # convert from rgb colour (as decimal) to hex code
    elif type(colour) == tuple:
        rgb_colour = tuple(round(255*x) for x in colour)
        hex_colour = f'#{int(round(rgb_colour[0])):02x}{int(round(rgb_colour[1])):02x}{int(round(rgb_colour[2])):02x}'

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
    standard_parameter_speci = get_standard_parameters_by_speci(speci, read_instance.ghost_version)
    initial_units = units
    final_units = read_instance.measurement_units[speci]

    # convert units using conversion factor
    conversion_factor = get_conversion_factor(initial_units, final_units, standard_parameter_speci) 
    if isinstance(conversion_factor, str):
        read_instance.logger.error(conversion_factor)
        sys.exit(1)
    RV *= conversion_factor
    exc_threshold *= conversion_factor
        
    return RV, exc_threshold