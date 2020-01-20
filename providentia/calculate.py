"""
Provides functions for basic statistic
calculations and experiment bias evaluation

Mean, Percentile, Standard Deviation
Variance, Data Availability Fraction
"""

import numpy as np
import scipy

from providentia.configuration import sequential_colourmap_warm


def calculate_mean(data):
    """Calculate mean in a dataset"""
    return np.mean(data)


def calculate_percentile(data, percentile=50.0):
    """Calculate specific percentile in a dataset"""
    return np.percentile(data, percentile)


def calculate_standard_deviation(data):
    """Calculate standard deviation in a dataset"""
    return np.std(data)


def calculate_variance(data):
    """Calculate variance in a dataset"""
    return np.var(data)


def calculate_data_avail_fraction(data):
    """Calculate data availability fraction
    (i.e. fraction of total data array not equal to NaN)
    """
    return (100./data.shape[-1]) * (np.count_nonzero(~np.isnan(data), axis=-1))


def calculate_data_avail_number(data):
    """Calculate data availability absolute number
    (i.e. number of total data measurements not equal to NaN)
    """
    return np.count_nonzero(~np.isnan(data), axis=-1)

def calculate_apbe(obs, exp):
    """Calculate absolute percent bias error (APBE)
    between observations and experiment
    """
    return 100.0*np.sum(np.abs(exp-obs))/np.sum(obs)


def calculate_pbe(obs, exp):
    """Calculate percent bias error (PBE)
    between observations and experiment
    """
    return 100.0*np.sum(exp-obs)/np.sum(obs)


def calculate_coe(obs, exp):
    """Calculate coefficient of efficiency (COE) between observations and experiment,
       based on Legates and McCabe (1999, 2012). There have been many suggestions for
       measuring model performance over the years, but the COE is a simple formulation
       which is easy to interpret. A perfect model has a COE = 1. As noted by Legates
       although the COE has no lower bound, a value of COE = 0.0 has a fundamental meaning.
       It implies that the model is no more able to predict the observed values
       than does the observed mean. Therefore, since the model can explain no more of the
       variation in the observed values than can the observed mean, such a model can have
       no predictive advantage. For negative values of COE, the model is less effective than
       the observed mean in predicting the variation in the observations.
       References:
       Legates DR, McCabe GJ. (1999). Evaluating the use of goodness-of-fit measures in hydrologic
       and hydroclimatic model validation. Water Resources Research 35(1): 233-241.
       Legates DR, McCabe GJ. (2012). A refined index of model performance: a rejoinder,
       International Journal of Climatology.
    """
    return 1.0 - (np.mean(np.abs(exp-obs)) / np.mean(np.abs(obs-np.mean(obs))))


def calculate_ioa(obs, exp):
    """Calculate the Index of Agreement (IOA) between observations and experiment, based on Willmott et al. (2011)
       The metric spans between -1 and +1 with values approaching +1 representing better model performance.
       An IOA of 0.5, for example, indicates that the sum of the error-magnitudes is one half of the sum
       of the observed-deviation magnitudes.
       When IOA = 0.0, it signifies that the sum of the magnitudes of the errors
       and the sum of the observed-deviation magnitudes are equivalent.
       When IOA = -0.5, it indicates that the sum of the error-magnitudes is twice
       the sum of the perfect model-deviation and observed-deviation magnitudes.
       Values of IOA near -1.0 can mean that the model-estimated deviations about O
       are poor estimates of the observed deviations; but, they also can mean that there
       simply is little observed variability - so some caution is needed when the IOA approaches -1.
       References;
       Willmott, C.J., Robeson, S.M., Matsuura, K., 2011. A refined index of model performance. International
       Journal of Climatology.
    """
    return 1.0 - (np.sum((obs-exp)**2))/(np.sum((np.abs(exp-np.mean(obs)) +
                                                 np.abs(obs-np.mean(obs)))**2))


def calculate_mae(obs, exp, normalisation_type='none'):
    """Calculate mean absolute error (MAE)/ normalised mean absolute
     error (NMAE) between observations and experiment
     """
    mae = np.mean(np.abs(exp-obs))
    # handle normalisation if desired
    if normalisation_type == 'max_min':
        mae = mae / (np.max(obs) - np.min(obs))
    elif normalisation_type == 'mean':
        mae = mae / np.mean(obs)
    elif normalisation_type == 'iq':
        mae = mae / (np.percentile(obs, 75) - np.percentile(obs, 25))
    elif normalisation_type == 'stdev':
        mae = mae / np.std(obs)
    return mae


def calculate_mbe(obs, exp, normalisation_type='none'):
    """Calculate mean bias error (MBE)/ normalised mean bias
    error (NMBE) between observations and experiment
    """
    mbe = np.mean(exp-obs)
    # handle normalisation if desired
    if normalisation_type == 'max_min':
        mbe = mbe / (np.max(obs) - np.min(obs))
    elif normalisation_type == 'mean':
        mbe = mbe / np.mean(obs)
    elif normalisation_type == 'iq':
        mbe = mbe / (np.percentile(obs, 75) - np.percentile(obs, 25))
    elif normalisation_type == 'stdev':
        mbe = mbe / np.std(obs)
    return mbe


def calculate_rmse(obs, exp, normalisation_type='none'):
    """Calculate root mean squared error (RMSE) /
    normalised root mean squared error (NRMSE)
    between observations and experiment
    """
    rmse = np.sqrt(np.mean((exp-obs)**2))
    # handle normalisation if desired
    if normalisation_type == 'max_min':
        rmse = rmse / (np.max(obs) - np.min(obs))
    elif normalisation_type == 'mean':
        rmse = rmse / np.mean(obs)
    elif normalisation_type == 'iq':
        rmse = rmse / (np.percentile(obs, 75) - np.percentile(obs, 25))
    elif normalisation_type == 'stdev':
        rmse = rmse / np.std(obs)
    return rmse


def calculate_r(obs, exp):
    """Calculate the Pearson correlation coefficient (r) between observations and experiment
       The Pearson correlation coefficient measures the linear relationship between two datasets.
       Strictly speaking, Pearson’s correlation requires that each dataset be normally distributed.
       Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation.
       Correlations of -1 or +1 imply an exact linear relationship.
       Positive correlations imply that as x increases, so does y.
       Negative correlations imply that as x increases, y decreases.
    """
    return scipy.stats.pearsonr(obs, exp)[0]


def calculate_r_squared(obs, exp):
    """Calculate the coefficient of determination, r squared, between observations and experiment
       It is the proportion of the variance in the dependent variable
       that is predictable from the independent variable(s).
       In linear least squares multiple regression with an estimated intercept term,
       the r squared equals the square of the Pearson correlation coefficient
    """
    return calculate_r(obs, exp)**2


def calculate_fac2(obs, exp):
    """Calculate fraction of experiment values within
    a factor of two of observed values (FAC2)
    """
    frac = exp/obs
    return (100.0/len(frac)) * len(frac[(frac >= 0.5) & (frac <= 2.0)])


def calculate_upa(obs, exp):
    """Calculate unpaired peak accuracy (UPA)"""
    obs_max = np.max(obs)
    exp_max = np.max(exp)
    return (exp_max - obs_max) - obs_max
