import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.signal import convolve, gaussian
import seaborn as sns

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

        :return: numbering of months/days
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

        :return dictionary of xticks per temporal resolution
        :rtype dict
    """

    return {'hour': np.arange(24, dtype=np.int), 
            'dayofweek': np.arange(7, dtype=np.int), 
            'month': np.arange(1, 13, dtype=np.int)}


def periodic_labels():
    """ Return axes labels for periodic subplots.

        :return: axes labels
        :rtype: dict
    """

    return {'hour':'H', 'dayofweek':'DoW', 'month':'M'}


def get_land_polygon_resolution(selection):
    """ Get resolution of land polygons to plot on map.

        :param selection: name of selected temporal resolution
        :type resolution: str
        :return: selected land polygon resolution
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
    clrs = sns.color_palette(instance.plot_characteristics_templates['general']['legend_color_palette'], n_colors=len(instance.data_labels)-1)

    # add colour and zorder for observations
    instance.plotting_params['observations']['colour'] = instance.plot_characteristics_templates['general']['obs_markerfacecolor']
    instance.plotting_params['observations']['zorder'] = instance.plot_characteristics_templates['general']['obs_zorder']

    # add colours and zorder for each experiment
    experiment_ind = 1
    for data_label in instance.data_labels:
        if data_label != 'observations':
            # define colour for experiment
            instance.plotting_params[data_label]['colour'] = clrs[experiment_ind-1]
            # define zorder for experiment (obs zorder + experiment_ind)
            instance.plotting_params[data_label]['zorder'] = \
                instance.plotting_params['observations']['zorder'] + experiment_ind
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
    #Make grid and discretize the data and round it to the next power of 2
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

    Parameters
    ----------
    x : float
        Value
    decimal_places : int
        Desired number of decimal places
    """

    # if cell value is not nan
    if x == x:
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
    