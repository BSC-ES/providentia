"""
Provides functions for basic statistic
calculations and experiment bias evaluation

Mean, Percentile, Standard Deviation
Variance, Data Availability Fraction
"""

import numpy as np
import scipy


class Stats(object):

    @staticmethod
    def calculate_mean(data):
        """Calculate mean in a dataset"""
        return np.mean(data)

    @staticmethod
    def calculate_percentile(data, percentile=50.0):
        """Calculate specific percentile in a dataset"""
        return np.percentile(data, percentile)

    @staticmethod
    def calculate_standard_deviation(data):
        """Calculate standard deviation in a dataset"""
        return np.std(data)

    @staticmethod
    def calculate_variance(data):
        """Calculate variance in a dataset"""
        return np.var(data)

    @staticmethod
    def calculate_minimum(data):
        """Calculate minimum in a dataset"""
        return np.min(data)

    @staticmethod
    def calculate_maximum(data):
        """Calculate minimum in a dataset"""
        return np.max(data)

    @staticmethod
    def calculate_data_avail_fraction(data):
        """Calculate data availability fraction
        (i.e. fraction of total data array not equal to NaN)
        """
        return (100. / np.array(data).shape[-1]) * \
               (np.count_nonzero(~np.isnan(data), axis=-1))

    @staticmethod
    def calculate_data_avail_number(data):
        """Calculate data availability absolute number
        (i.e. number of total data measurements not equal to NaN)
        """
        return np.count_nonzero(~np.isnan(data), axis=-1)

    @staticmethod
    def max_repeated_nans_fraction(data):
        """Get % of total period of the maximum run of consecutive NaNs in array"""
        max_gap_pc = []

        for station_data in data:
            mask = np.concatenate(([False], np.isnan(station_data), [False]))
            if ~mask.any():
                max_gap_pc.append(0)
            else:
                idx = np.nonzero(mask[1:] != mask[:-1])[0]
                max_gap_pc.append((idx[1::2] - idx[::2]).max())

        return np.array(max_gap_pc) * (100. / data.shape[1])


class ExpBias(object):

    @staticmethod
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
        return 1.0 - (np.mean(np.abs(exp - obs)) /
                      np.mean(np.abs(obs - np.mean(obs))))

    @staticmethod
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
        return 1.0 - (np.sum((obs - exp) ** 2)) /\
                     (np.sum((np.abs(exp - np.mean(obs)) +
                              np.abs(obs - np.mean(obs))) ** 2))

    @staticmethod
    def calculate_mb(obs, exp, normalisation_type='none'):
        """Calculate mean bias (MB), or normalised derivation.
           The difference between a modelled and an observed value,
           𝑀𝑖 − 𝑂𝑖 , is referred to as the bias.
           The mean bias is simply the average bias between the modelled and observed values.
        """

        mbe = np.mean(exp - obs)
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

    @staticmethod
    def calculate_mae(obs, exp, normalisation_type='none'):
        """Calculate mean absolute error (MAE), or normalised derivation.
           It is calculated from the absolute of the difference between a modelled
           and an observed value,|𝑀𝑖 −𝑂𝑖|. Therefore the mean absolute error is always positive.
           This metric can highlight reveal somes biases not seen using the MB metric, where
           postive and negative biases can average out to be zero.
        """

        mae = np.mean(np.abs(exp - obs))
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

    @staticmethod
    def calculate_mnb(obs, exp):
        """Calculate mean normalised bias (MNB).
           The mean normalised bias (MNB) is calculated in a similar fashion to the mean bias.
           The mean normalised bias is calculated from the difference between the modelled and observed values
           (i.e., the bias, 𝑀𝑖 − 𝑂𝑖) is normalised (divided) by the observed value (𝑂𝑖).
           The mean normalised bias is reported as a percentage.
        """

        mnb = np.mean((exp - obs) / obs)
        return mnb

    @staticmethod
    def calculate_mnae(obs, exp):
        """Calculate mean normalised absolute error (MNAE).
           The mean normalised absolute error (MNAE) is calculated in a similar fashion to the mean absolute error.
           The mean normalised absolute error is calculated from the absolute of the bias, 𝑀𝑖 − 𝑂𝑖,
           normalised by the observed value, 𝑂𝑖. Therefore the mean normalised absolute error is always positive.
           The mean normalised absolute error is reported as a percentage.
        """

        mnae = np.mean((np.abs(exp - obs)) / obs)
        return mnae

    @staticmethod
    def calculate_mfb(obs, exp):
        """Calculate mean fractional bias (MFB).
           The mean fractional bias (MFB) is used as a substitute for the mean normalised bias (MNB),
           when the mean normalised bias becomes large.
           The mean normalised bias can become very large when a minimum threshold is not used for the observations.
           The fractional bias for cases with factors of 2 under-and over-prediction are -67 and +67%,
           respectively (as opposed to -50 and +100%, when using normalised bias).
           The mean fractional bias is a useful indicator because it has the advantage of equally weighting positive and
           negative bias estimates.
           It has also the advantage of not considering observations as the true value. The mean fractional bias can
           range in value from -200% to +200%.
        """

        mfb = np.mean((exp - obs) / ((exp + obs) / 2.))
        return mfb

    @staticmethod
    def calculate_mafb(obs, exp):
        """Calculate mean absolute fractional bias (MAFB).
        """

        mafb = np.mean(np.abs((exp - obs) / ((exp + obs) / 2.)))
        return mafb

    @staticmethod
    def calculate_rmse(obs, exp, normalisation_type='none'):
        """Calculate root mean squared error (RMSE) /
        normalised root mean squared error (NRMSE)
        between observations and experiment
        """
        rmse = np.sqrt(np.mean((exp - obs) ** 2))
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

    @staticmethod
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

    @staticmethod
    def calculate_r_squared(obs, exp):
        """Calculate the coefficient of determination, r squared, between observations and experiment
           It is the proportion of the variance in the dependent variable
           that is predictable from the independent variable(s).
           In linear least squares multiple regression with an estimated intercept term,
           the r squared equals the square of the Pearson correlation coefficient
        """
        return ExpBias.calculate_r(obs, exp) ** 2

    @staticmethod
    def calculate_fac2(obs, exp):
        """Calculate fraction of experiment values within
        a factor of two of observed values (FAC2)
        """
        frac = exp / obs
        return (100.0 / len(frac)) * len(frac[(frac >= 0.5) & (frac <= 2.0)])

    @staticmethod
    def calculate_upa(obs, exp):
        """Calculate unpaired peak accuracy (UPA)
           see here: https://gitlab.com/polyphemus/atmopy/-/blob/master/stat/measure.py
        """
        obs_max = np.max(obs)
        exp_max = np.max(exp)
        return (exp_max - obs_max) / obs_max
