"""
Provides functions for basic statistic
calculations and experiment bias evaluation.

Function defintions mainly stem from: 
https://www.tandfonline.com/doi/pdf/10.1080/10962247.2016.1265027 ,
https://www.cmascenter.org/conference/2003/session_poster/yu_abstract3.pdf ,
and https://github.com/davidcarslaw/openair/blob/HEAD/R/modStats.R 
"""

import numpy as np
import numpy.ma as ma
import scipy

class Stats(object):

    @staticmethod
    def calculate_mean(data):
        """ Calculate mean in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: mean value of data
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.nanmean(data, axis=-1)

    @staticmethod
    def calculate_median(data):
        """ Calculate median in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: median value of data
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.nanmedian(data, axis=-1)

    @staticmethod
    def calculate_percentile(data, percentile=50.0):
        """ Calculate specific percentile in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :param percentile: percentile, default = 50.0
            :type percentile: float
            :return: percentile of data
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.nanpercentile(data, percentile, axis=-1)

    @staticmethod
    def calculate_standard_deviation(data):
        """ Calculate standard deviation in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: standard deviation value of data
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.nanstd(data, axis=-1)

    @staticmethod
    def calculate_variance(data):
        """ Calculate variance in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: variance value of data
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.nanvar(data, axis=-1)

    @staticmethod
    def calculate_minimum(data):
        """ Calculate minimum in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: min value of data
            :rtype: numpy.float64
        """
        if len(data) == 0:
            return np.NaN
        else:
            return np.nanmin(data,axis=-1)

    @staticmethod
    def calculate_maximum(data):
        """ Calculate minimum in a dataset.

            :param data: array of data
            :type data: numpy.ndarray
            :return: max value of data
            :rtype: numpy.float64
        """
        if len(data) == 0:
            return np.NaN
        else:
            return np.nanmax(data,axis=-1)

    @staticmethod
    def calculate_data_avail_fraction(data):
        """ Calculate data availability fraction
            (i.e. fraction of total data array not equal to NaN).

            :param data: array of data
            :type data: numpy.ndarray
            :return: data availability percent
            :rtype: numpy.float64
        """
        if data.size == 0:
            return np.NaN
        else:
            return (100. / data.shape[-1]) * \
                   (np.count_nonzero(~np.isnan(data), axis=-1))

    @staticmethod
    def calculate_data_avail_number(data):
        """ Calculate data availability absolute number
            (i.e. number of total data measurements not equal to NaN).
        """
        if data.size == 0:
            return np.NaN
        else:
            return np.count_nonzero(~np.isnan(data), axis=-1).astype('float32')

    @staticmethod
    def max_repeated_nans_fraction(data):
        """ Get % of total period of the maximum run of consecutive NaNs in array. """
        if data.size == 0:
            return np.NaN
        else:
            max_gap_pc = []
            for station_ind in range(data.shape[-2]):
                mask = np.concatenate(([False],np.isnan(data[station_ind]),[False]))
                if ~mask.any():
                    max_gap_pc.append(0)
                else:
                    idx = np.nonzero(mask[1:] != mask[:-1])[0]
                    max_gap_pc.append((idx[1::2] - idx[::2]).max())

            return np.array(max_gap_pc) * (100. / data.shape[-1])

    @staticmethod
    def calculate_exceedances(data, threshold=0):
        """ Calculate number of data exceedances
            (i.e. number of measurements exceeding a set threshold).
        """
        if data.size == 0:
            return np.NaN
        else:
            n_exceed = np.nansum(data > threshold, axis=-1).astype('float32')
            return n_exceed

class ExpBias(object):

    @staticmethod
    def calculate_coe(obs, exp):
        """ Calculate coefficient of efficiency (COE) between observations and experiment,
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
        if obs.size == 0:
            return np.NaN
        else:
            return 1.0 - np.nansum(np.abs(exp - obs), axis=-1) / \
                   np.nansum(np.abs(obs - np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)), axis=-1) 

    @staticmethod
    def calculate_ioa(obs, exp):
        """ Calculate the Index of Agreement (IOA) between observations and experiment, based on Willmott et al. (2011)
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
        if obs.size == 0:
            return np.NaN
        else:
            lhs = np.nansum(np.abs(exp - obs), axis=-1)
            rhs = 2.0 * np.nansum(np.abs(obs - np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)), axis=-1)
            output = np.copy(lhs)
            lower_check = lhs <= rhs
            output[lower_check] = 1.0 - lhs[lower_check] / rhs[lower_check] 
            output[~lower_check] = rhs[~lower_check] / lhs[~lower_check] - 1.0
            return output

    @staticmethod
    def calculate_mb(obs, exp, normalisation_type='none'):
        """ Calculate mean bias (MB), or normalised derivation (NMB).
            The difference between a modelled and an observed value,
            𝑀𝑖 − 𝑂𝑖 , is referred to as the bias.
            The mean bias is simply the average bias between the modelled and observed values.
            This statistic is equivalent to the 'Mean_bias' when temporal_colocation is active.
        """
        if obs.size == 0:
            return np.NaN
        else:
            mb = np.nanmean(exp - obs, axis=-1)

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                mb = (mb / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0
            elif normalisation_type == 'mean':
                mb = (mb / np.nanmean(obs, axis=-1)) * 100.0
            elif normalisation_type == 'sum':
                mb = (mb / np.nansum(obs, axis=-1)) * 100.0
            elif normalisation_type == 'iq':
                mb = (mb / (np.nanpercentile(obs, 75, axis=-1) - np.nanpercentile(obs, 25, axis=-1))) * 100.0
            elif normalisation_type == 'stdev':
                mb = (mb / np.nanstd(obs, axis=-1)) * 100.0
            return mb

    @staticmethod
    def calculate_me(obs, exp, normalisation_type='none'):
        """ Calculate mean error (ME), or normalised derivation (NME).
            It is calculated from the absolute of the difference between a modelled
            and an observed value,|𝑀𝑖 −𝑂𝑖|. Therefore the mean error is always positive.
            This metric can highlight reveal somes biases not seen using the MB metric, where
            postive and negative biases can average out to be zero.
            Otherwise known as mean gross error (MGE), mean absolute error (MAE), 
            and mean absolute gross error (MAGE); 
            and normalised form as normalised mean gross error (NMGE) and normalised mean absolute error (NMAE).
        """
        if obs.size == 0:
            return np.NaN
        else:
            me = np.nanmean(np.abs(exp - obs), axis=-1)

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                me = (me / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0 
            elif normalisation_type == 'mean':
                me = (me / np.nanmean(obs, axis=-1)) * 100.0
            elif normalisation_type == 'sum':
                me = (me / np.nansum(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'iq':
                me = (me / (np.nanpercentile(obs, 75, axis=-1) - np.nanpercentile(obs, 25, axis=-1))) * 100.0 
            elif normalisation_type == 'stdev':
                me = (me / np.nanstd(obs, axis=-1)) * 100.0 
            return me

    @staticmethod
    def calculate_mnb(obs, exp):
        """ Calculate mean normalised bias (MNB).
            The mean normalised bias (MNB) is calculated in a similar fashion to the mean bias.
            The mean normalised bias is calculated from the difference between the modelled and observed values
            (i.e. the bias, 𝑀𝑖 − 𝑂𝑖) is normalised (divided) by the observed value (𝑂𝑖).
        """
        if obs.size == 0:
            return np.NaN
        else:
            mnb = np.nanmean(ma.masked_invalid((exp - obs) / obs), axis=-1) * 100.0
            return mnb

    @staticmethod
    def calculate_mne(obs, exp):
        """ Calculate mean normalised error (MNE).
            The mean normalised error (MNE) is calculated in a similar fashion to the mean error.
            The mean normalised error is calculated from the absolute of the bias, 𝑀𝑖 − 𝑂𝑖,
            normalised by the observed value, 𝑂𝑖. Therefore the mean normalised error is always positive.
            Otherwise known as mean normalised absolute error (MNAE).
        """
        if obs.size == 0:
            return np.NaN
        else:
            mne = np.nanmean(ma.masked_invalid((np.abs(exp - obs)) / obs), axis=-1) * 100.0
            return mne
    
    @staticmethod
    def calculate_mfb(obs, exp):
        """ Calculate mean fractional bias (MFB).
            The mean fractional bias (MFB) is used as a substitute for the mean normalised bias (MNB),
            when the MNB becomes large.
            The MNB can become very large when a minimum threshold is not used for the observations.
            The mean fractional bias for cases with factors of 2 under-and over-prediction are -67 and +67%,
            respectively (as opposed to -50 and +100%, when using normalised bias).
            The mean fractional bias is a useful indicator because it has the advantage of equally weighting positive and
            negative bias estimates.
            It has also the advantage of not considering observations as the true value. The mean fractional bias can
            range in value from -200% to +200%.
            Otherwise known as fractional bias (FB).
        """
        if obs.size == 0:
            return np.NaN
        else:
            mfb = np.nanmean((exp - obs) / ma.masked_invalid((exp + obs) / 2.0), axis=-1) * 100.0
            return mfb

    @staticmethod
    def calculate_mfe(obs, exp):
        """ Calculate mean fractional error (MFE).
            Otherwise known as fractional error (FE), fractional gross error (FGE), 
            or mean absolute fractional bias (MAFB).
        """
        if obs.size == 0:
            return np.NaN
        else:
            mfe = np.nanmean(np.abs(exp - obs) / ma.masked_invalid((exp + obs) / 2.0), axis=-1) * 100.0
            return mfe

    @staticmethod
    def calculate_rmse(obs, exp, normalisation_type='none'):
        """ Calculate root mean squared error (RMSE) /
            normalised root mean squared error (NRMSE)
            between observations and experiment.
        """

        if obs.size == 0:
            return np.NaN
        else:
            rmse = np.sqrt(np.nanmean((exp - obs) ** 2, axis=-1))

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                rmse = (rmse / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0 
            elif normalisation_type == 'mean':
                rmse = (rmse / np.nanmean(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'rmse':
                rmse = (rmse / np.nansum(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'iq':
                rmse = (rmse / (np.nanpercentile(obs, 75, axis=-1) - np.nanpercentile(obs, 25, axis=-1))) * 100.0 
            elif normalisation_type == 'stdev':
                rmse = (rmse / np.nanstd(obs, axis=-1)) * 100.0 
            return rmse

    @staticmethod
    def calculate_r(obs, exp):
        """ Calculate the Pearson correlation coefficient (r) between observations and experiment
            The Pearson correlation coefficient measures the linear relationship between two datasets.
            Strictly speaking, Pearson’s correlation requires that each dataset be normally distributed.
            Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation.
            Correlations of -1 or +1 imply an exact linear relationship.
            Positive correlations imply that as x increases, so does y.
            Negative correlations imply that as x increases, y decreases.
        """

        if obs.size == 0:
            return np.NaN
        else:
            mean_obs = np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)
            std_obs = np.expand_dims(np.nanstd(obs, axis=-1), axis=-1)
            mean_exp = np.expand_dims(np.nanmean(exp, axis=-1), axis=-1)
            std_exp = np.expand_dims(np.nanstd(exp, axis=-1), axis=-1)
            standard_score_obs = ma.masked_invalid((obs - mean_obs) / std_obs)
            standard_score_exp = ma.masked_invalid((exp - mean_exp) / std_exp)
            standard_score_mult = standard_score_obs*standard_score_exp
            # get number of non-masked values in the time dimension
            n = standard_score_mult.count(axis=-1)
            return np.nansum(standard_score_mult, axis=-1) / n

    @staticmethod
    def calculate_r_squared(obs, exp):
        """ Calculate the coefficient of determination, r squared, between observations and experiment
            It is the proportion of the variance in the dependent variable
            that is predictable from the independent variable(s).
            In linear least squares multiple regression with an estimated intercept term,
            the r squared equals the square of the Pearson correlation coefficient.
        """
        if obs.size == 0:
            return np.NaN
        else:
            return ExpBias.calculate_r(obs, exp) ** 2

    @staticmethod
    def calculate_fac2(obs, exp):
        """ Calculate fraction of experiment values within
            a factor of two of observed values (FAC2)
        """
        if obs.size == 0:
            return np.NaN
        else:
            frac = ma.masked_invalid(exp / obs)
            n = frac.count(axis=-1)
            return (100.0 / n) * np.nansum(((frac >= 0.5) & (frac <= 2.0)), axis=-1)

    @staticmethod
    def calculate_upa(obs, exp):
        """ Calculate unpaired peak accuracy (UPA).
            See here: https://gitlab.com/polyphemus/atmopy/-/blob/master/stat/measure.py.
        """
        if obs.size == 0:
            return np.NaN
        else:
            obs_max = np.nanmax(obs, axis=-1)
            exp_max = np.nanmax(exp, axis=-1)
            return (exp_max - obs_max) / obs_max
