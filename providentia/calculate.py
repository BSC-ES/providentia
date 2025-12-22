"""
Classes for the calculation of basic and model bias statistics.

Function defintions mainly stem from: 
https://www.tandfonline.com/doi/pdf/10.1080/10962247.2016.1265027 ,
https://www.cmascenter.org/conference/2003/session_poster/yu_abstract3.pdf ,
and https://github.com/davidcarslaw/openair/blob/HEAD/R/modStats.R 
"""

import copy
import numpy as np

def nansumwrapper(data, **kwargs):
    """Sum array elements over a given axis treating Not a Numbers as zero, returning NaN if all elements are NaN.

    Parameters
    ----------
    data : numpy.array
        Array containing numbers to sum.
    **kwargs : dict
        Arbitrary keyword arguments passed to the underlying numpy.nansum function.

    Returns
    -------
    numpy.array or scalar
        An array with the same shape as data, with the specified axis removed.
    """

    if np.isnan(data).all():
        return np.full(data.shape[:data.ndim-1], np.nan, dtype=np.float32)
    else:
        return np.nansum(data, **kwargs)

class Stats(object):
    """Class for the calculation of basic statistics."""

    @staticmethod
    def calculate_mean(data):
        """
        Calculate mean

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Mean
        """

        if data.size == 0:
            return np.nan
        else:
            return np.nanmean(data, axis=-1)

    @staticmethod
    def calculate_median(data):
        """
        Calculate median

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Median
        """

        if data.size == 0:
            return np.nan
        else:
            return np.nanmedian(data, axis=-1)

    @staticmethod
    def calculate_percentile(data, percentile=50.0):
        """
        Calculate specific percentile

        Parameters
        ----------
        data : numpy.array
            Data array
        percentile : float
            Percentile

        Returns
        -------
        float
            Percentile
        """

        if data.size == 0:
            return np.nan
        else:
            return np.nanpercentile(data, percentile, axis=-1, method='nearest')

    @staticmethod
    def calculate_standard_deviation(data):
        """
        Calculate standard deviation

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Standard deviation
        """

        if data.size == 0:
            return np.nan
        else:
            return np.nanstd(data, axis=-1)

    @staticmethod
    def calculate_variance(data):
        """
        Calculate variance

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Variance
        """

        if data.size == 0:
            return np.nan
        else:
            return np.nanvar(data, axis=-1)

    @staticmethod
    def calculate_minimum(data):
        """
        Calculate minimum

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Minimum
        """

        if len(data) == 0:
            return np.nan
        else:
            return np.nanmin(data,axis=-1)

    @staticmethod
    def calculate_maximum(data):
        """
        Calculate maximum

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Maximum
        """

        if len(data) == 0:
            return np.nan
        else:
            return np.nanmax(data,axis=-1)

    @staticmethod
    def calculate_data_avail_fraction(data, nan_padding_counts=None):
        """
        Calculate data availability fraction (% of non-NaN values).
        Corrects for padded NaNs if nan_padding_counts is provided.

        Parameters
        ----------
        data : numpy.array
            Input data array. Can be:
            - (label, station, time) for raw data
            - (chunk, label, station, chunk_time) for grouped data
        nan_padding_counts : np.ndarray or None, optional
            Number of padded NaNs per chunk (from grouping).
            Only relevant when data is grouped.

        Returns
        -------
        numpy.array
            Data availability percentage.
            Shape matches data[..., 0] (i.e. everything except time axis).
        """

        if data.size == 0:
            return np.nan

        # count valid values
        valid_counts = np.count_nonzero(~np.isnan(data), axis=-1)

        # denominator: correct for padded NaNs if provided
        if nan_padding_counts is not None:
            denom = data.shape[-1] - nan_padding_counts
        else:
            denom = data.shape[-1]

        # avoid division by zero
        denom = np.where(denom == 0, np.nan, denom)

        # return data %
        return (100.0 / denom) * valid_counts


    @staticmethod
    def calculate_stations_number(data, statistic_mode, statistic_aggregation, per_station,
                                  periodic_statistic_mode=None, periodic_statistic_aggregation=None):
        """
        Calculate number of stations

        Parameters
        ----------
        data : numpy.array
            Data array
        statistic_mode : str
            Statistic mode
        statistic_aggregation : str
            Statistic aggregation
        per_station : bool
            Per station
        periodic_statistic_mode : str
            Periodic statistic mode
        periodic_statistic_aggregation : str
            Periodic statistic aggregation
            
        Returns
        -------
        float
            Number of stations
        """

        if data.size == 0:
            return np.nan
        else:

            # if array is 5d are reading daily forecast data, so make neccesary adjustments to calculation
            if data.ndim == 5:
                relevant_dim = -3
            else:
                relevant_dim = -2

            # get number of stations that do not have missing values
            # if calculating periodic statistic for cycle mode, then aggregate across first dimension
            # this due to inconsistent way Nstations is calculated with respect to other statistics
            if periodic_statistic_mode == 'Cycle':
                stations_number = np.count_nonzero(~np.isnan(data), axis=0).astype('float32').transpose()
            else:
                stations_number = np.count_nonzero(~np.isnan(data), axis=relevant_dim).astype('float32')

            # if calculating periodic statistic for independent mode, then aggregate across last dimension
            # this due to inconsistent way Nstations is calculated with respect to other statistics
            if periodic_statistic_mode == 'Independent':
                from .statistics import aggregation
                stations_number = aggregation(stations_number, periodic_statistic_aggregation, axis=-1)

            # do aggregation (if not calculating periodic statistic, or per station)
            if (not periodic_statistic_mode) & (not per_station):
                from .statistics import aggregation
                if statistic_mode in ['Temporal|Spatial','Spatial|Temporal']:
                    stations_number = aggregation(stations_number, statistic_aggregation, axis=-1)
                elif statistic_mode == 'Flattened':
                    stations_number = aggregation(stations_number, 'Median', axis=-1)

            return stations_number


    @staticmethod
    def calculate_data_avail_number(data):
        """
        Calculate data availability absolute number
        (i.e. number of total data measurements not equal to NaN).

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Data availability absolute number
        """

        if data.size == 0:
            return np.nan
        else:
            return np.count_nonzero(~np.isnan(data), axis=-1).astype('float32')

    @staticmethod
    def calculate_exceedances(data, threshold=0):
        """
        Calculate number of data exceedances
        (i.e. number of measurements exceeding a set threshold).

        Parameters
        ----------
        data : numpy.array
            Data array

        Returns
        -------
        float
            Number of data exceedances
        """

        if data.size == 0:
            return np.nan
        else:
            return nansumwrapper(data > threshold, axis=-1).astype('float32')

    @staticmethod
    def calculate_mda8(data, statistic_mode, statistic_aggregation, per_station,
                       periodic_statistic_mode=None):
        """
        Calculate MDA8 (daily maximum 8 hour average) 

        Parameters
        ----------
        data : numpy.array
            Data array
        statistic_mode : str
            Statistic mode
        statistic_aggregation : str
            Statistic aggregation
        per_station : bool
            Per station
        periodic_statistic_mode : str
            Periodic statistic mode
            
        Returns
        -------
        float
            MDA8
        """

        from .statistics import aggregation

        if data.size == 0:
            return np.nan
        else:

            # if have periodic statistic, then cannot calculate MDA8, so return NaN
            if periodic_statistic_mode:
                return np.nan

            # if array is 5d are reading daily forecast data, so make neccesary adjustments to calculation
            if data.ndim == 5:
                agg_dim = -2
            else:
                agg_dim = -1

            # if last dimension is not 24, then reshape it to be so
            extra_agg = False
            if (data.shape[-1] != 24):
                data = data.reshape(*data.shape[:-1], -1, 24)
                extra_agg = True

            # calculate MDA8
            start_inds = np.arange(0,17)
            end_inds = np.arange(8,25)            
            mda8_arr = np.full((*data.shape[:-1], 17), np.nan, dtype=np.float32)

            for window_ind, (start_ind, end_ind) in enumerate(zip(start_inds, end_inds)):
                mda8_arr[..., window_ind] = np.nanmean(data[..., start_ind:end_ind], axis=-1)

            mda8 = np.nanmax(mda8_arr, axis=-1)

            # do aggregation (if not calculating periodic statistic, or per station)
            if (not periodic_statistic_mode) and (not per_station):
                if statistic_mode in ['Temporal|Spatial','Spatial|Temporal']:
                    mda8 = aggregation(mda8, statistic_aggregation, axis=-1)
                elif statistic_mode == 'Flattened':
                    mda8 = aggregation(mda8, 'Median', axis=-1)

            # do extra aggregation for daily forecast data (if not calculating periodic statistic, or per station)
            if (extra_agg) and (not periodic_statistic_mode) and (not per_station):
                if statistic_mode in ['Temporal|Spatial','Spatial|Temporal']:
                    mda8 = aggregation(mda8, statistic_aggregation, axis=agg_dim)
                elif statistic_mode == 'Flattened':
                    mda8 = aggregation(mda8, 'Median', axis=agg_dim)

            return mda8

    @staticmethod
    def calculate_rms_u(data, u_95r_RV, RV, alpha):
        """
        Calculate FAIRMODE RMSu
        See here: https://publications.jrc.ec.europa.eu/repository/handle/JRC129254

        Parameters
        ----------
        data : numpy.array
            Data array
        u_95r_RV : float
            Uncertainty around the reference value
        RV : float
            Reference value
        alpha : float
            Uncertainty parameter
            
        Returns
        -------
        float
            RMSu
        """

        if data.size == 0:
            return np.nan
        else:
            rms_u = u_95r_RV * np.sqrt(
                (1 - alpha ** 2) * (Stats.calculate_mean(data) ** 2 
                + Stats.calculate_standard_deviation(data) ** 2) 
                + (alpha * RV) ** 2)
            return rms_u


class ModBias(object):
    """Class for the calculation of model bias statistics."""

    @staticmethod
    def calculate_coe(obs, mod):
        """
        Calculate coefficient of efficiency (COE) between observations and model,
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

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            COE
        """

        if obs.size == 0:
            return np.nan
        else:
            return 1.0 - nansumwrapper(np.abs(mod - obs), axis=-1) / \
                   nansumwrapper(np.abs(obs - np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)), axis=-1) 

    @staticmethod
    def calculate_ioa(obs, mod):
        """
        Calculate the Index of Agreement (IOA) between observations and model, based on Willmott et al. (2011)
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

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            IOA
        """

        if obs.size == 0:
            return np.nan
        else:
            lhs = nansumwrapper(np.abs(mod - obs), axis=-1)
            rhs = 2.0 * nansumwrapper(np.abs(obs - np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)), axis=-1)
            output = np.copy(lhs)
            lower_check = lhs <= rhs
            output[lower_check] = 1.0 - lhs[lower_check] / rhs[lower_check] 
            output[~lower_check] = rhs[~lower_check] / lhs[~lower_check] - 1.0
            return output

    @staticmethod
    def calculate_mb(obs, mod, normalisation_type='none'):
        """
        Calculate mean bias (MB), or normalised derivation (NMB).
        The difference between a modelled and an observed value,
        𝑀𝑖 − 𝑂𝑖 , is referred to as the bias.
        The mean bias is simply the average bias between the modelled and observed values.
        This statistic is equivalent to the 'Mean_bias' when temporal_colocation is active.

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
        normalisation_type : str
            Normalisation type
            
        Returns
        -------
        float
            Mean bias
        """

        if obs.size == 0:
            return np.nan
        else:
            mb = np.nanmean(mod - obs, axis=-1)

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                mb = (mb / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0
            elif normalisation_type == 'mean':
                mb = (mb / np.nanmean(obs, axis=-1)) * 100.0
            elif normalisation_type == 'sum':
                mb = (mb / nansumwrapper(obs, axis=-1)) * 100.0
            elif normalisation_type == 'iq':
                mb = (mb / (np.nanpercentile(obs, 75, axis=-1, method='nearest') - np.nanpercentile(obs, 25, axis=-1, method='nearest'))) * 100.0
            elif normalisation_type == 'stdev':
                mb = (mb / np.nanstd(obs, axis=-1)) * 100.0
            return mb

    @staticmethod
    def calculate_me(obs, mod, normalisation_type='none'):
        """
        Calculate mean error (ME), or normalised derivation (NME).
        It is calculated from the absolute of the difference between a modelled
        and an observed value,|𝑀𝑖 −𝑂𝑖|. Therefore the mean error is always positive.
        This metric can highlight reveal somes biases not seen using the MB metric, where
        postive and negative biases can average out to be zero.
        Otherwise known as mean gross error (MGE), mean absolute error (MAE), 
        and mean absolute gross error (MAGE); 
        and normalised form as normalised mean gross error (NMGE) and normalised mean absolute error (NMAE).

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
        normalisation_type : str
            Normalisation type
            
        Returns
        -------
        float
            ME
        """

        if obs.size == 0:
            return np.nan
        else:
            me = np.nanmean(np.abs(mod - obs), axis=-1)

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                me = (me / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0 
            elif normalisation_type == 'mean':
                me = (me / np.nanmean(obs, axis=-1)) * 100.0
            elif normalisation_type == 'sum':
                me = (me / nansumwrapper(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'iq':
                me = (me / (np.nanpercentile(obs, 75, axis=-1, method='nearest') - np.nanpercentile(obs, 25, axis=-1, method='nearest'))) * 100.0 
            elif normalisation_type == 'stdev':
                me = (me / np.nanstd(obs, axis=-1)) * 100.0 
            return me

    @staticmethod
    def calculate_mnb(obs, mod):
        """
        Calculate mean normalised bias (MNB).
        The mean normalised bias (MNB) is calculated in a similar fashion to the mean bias.
        The mean normalised bias is calculated from the difference between the modelled and observed values
        (i.e. the bias, 𝑀𝑖 − 𝑂𝑖) is normalised (divided) by the observed value (𝑂𝑖).
        
        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array

        Returns
        -------
        float
            MNB
        """

        if obs.size == 0:
            return np.nan
        else:
            # to avoid ZeroDivisionError, replace all obs of 0 with NaN
            obs_nan = copy.deepcopy(obs)
            obs_nan[obs_nan == 0.0] = np.nan
            mnb = np.nanmean((mod - obs_nan) / obs_nan, axis=-1) * 100.0
            return mnb

    @staticmethod
    def calculate_mne(obs, mod):
        """
        Calculate mean normalised error (MNE).
        The mean normalised error (MNE) is calculated in a similar fashion to the mean error.
        The mean normalised error is calculated from the absolute of the bias, 𝑀𝑖 − 𝑂𝑖,
        normalised by the observed value, 𝑂𝑖. Therefore the mean normalised error is always positive.
        Otherwise known as mean normalised absolute error (MNAE).

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array

        Returns
        -------
        float
            MNE
        """

        if obs.size == 0:
            return np.nan
        else:
            # to avoid ZeroDivisionError, replace all obs of 0 with NaN
            obs_nan = copy.deepcopy(obs)
            obs_nan[obs_nan == 0.0] = np.nan
            mne = np.nanmean((np.abs(mod - obs_nan)) / obs_nan, axis=-1) * 100.0
            return mne
    
    @staticmethod
    def calculate_mfb(obs, mod):
        """
        Calculate mean fractional bias (MFB).
        The mean fractional bias (MFB) is used as a substitute for the mean normalised bias (MNB),
        when the MNB becomes large.
        The MNB can become very large when a minimum threshold is not used for the observations.
        The mean fractional bias for cases with factors of 2 under-and over-prediction are -67 and +67%,
        respectively (as opposed to -50 and +100%, when using normalised bias).
        The mean fractional bias is a useful indicator because it has the advantage of equally weighting positive and
        negative bias estimates.
        It has also the advantage of not considering observations as the true value. The mean fractional bias can
        range in value from -200% to +200%.
        Otherwise known as fractional bias (FB) or modified normalized mean bias (MNMB).

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array

        Returns
        -------
        float
            MFB
        """

        if obs.size == 0:
            return np.nan
        else:
            lhs = mod - obs
            rhs = (mod + obs) / 2.0
            # to avoid ZeroDivisionError, replace all rhs values of 0 with NaN
            rhs[rhs == 0.0] = np.nan
            mfb = np.nanmean(lhs / rhs, axis=-1) * 100.0
            return mfb

    @staticmethod
    def calculate_mfe(obs, mod):
        """
        Calculate mean fractional error (MFE).
        Otherwise known as fractional error (FE), fractional gross error (FGE), 
        or mean absolute fractional bias (MAFB).

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array

        Returns
        -------
        float
            MFE
        """

        if obs.size == 0:
            return np.nan
        else:
            lhs = np.abs(mod - obs)
            rhs = (mod + obs) / 2.0
            # to avoid ZeroDivisionError, replace all rhs values of 0 with NaN
            rhs[rhs == 0.0] = np.nan
            mfe = np.nanmean(lhs / rhs, axis=-1) * 100.0
            return mfe

    @staticmethod
    def calculate_rmse(obs, mod, normalisation_type='none'):
        """
        Calculate root mean squared error (RMSE) /
        normalised root mean squared error (NRMSE)
        between observations and model.

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
        normalisation_type : str
            Normalisation type
            
        Returns
        -------
        float
            RMSE
        """

        if obs.size == 0:
            return np.nan
        else:
            rmse = np.sqrt(np.nanmean((mod - obs) ** 2, axis=-1))

            # handle normalisation if desired
            if normalisation_type == 'max_min':
                rmse = (rmse / (np.nanmax(obs, axis=-1) - np.nanmin(obs, axis=-1))) * 100.0 
            elif normalisation_type == 'mean':
                rmse = (rmse / np.nanmean(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'rmse':
                rmse = (rmse / nansumwrapper(obs, axis=-1)) * 100.0 
            elif normalisation_type == 'iq':
                rmse = (rmse / (np.nanpercentile(obs, 75, axis=-1, method='nearest') - np.nanpercentile(obs, 25, axis=-1, method='nearest'))) * 100.0 
            elif normalisation_type == 'stdev':
                rmse = (rmse / np.nanstd(obs, axis=-1)) * 100.0 
            return rmse

    @staticmethod
    def calculate_crmse(obs, mod):
        """
        Calculate FAIRMODE cRMSE
        See here: https://publications.jrc.ec.europa.eu/repository/handle/JRC129254

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            cRMSE
        """

        if obs.size == 0:
            return np.nan
        else:
            crmse = np.sqrt(
                np.mean(((mod - np.mean(mod)) - (obs - np.mean(obs))) ** 2)
            )
            return crmse
    
    @staticmethod
    def calculate_r(obs, mod):
        """
        Calculate the Pearson correlation coefficient (r) between observations and model
        The Pearson correlation coefficient measures the linear relationship between two datasets.
        Strictly speaking, Pearson's correlation requires that each dataset be normally distributed.
        Like other correlation coefficients, this one varies between -1 and +1 with 0 implying no correlation.
        Correlations of -1 or +1 imply an exact linear relationship.
        Positive correlations imply that as x increases, so does y.
        Negative correlations imply that as x increases, y decreases.

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            r
        """

        if obs.size == 0:
            return np.nan
        else:
            mean_obs = np.expand_dims(np.nanmean(obs, axis=-1), axis=-1)
            std_obs = np.expand_dims(np.nanstd(obs, axis=-1), axis=-1)
            mean_mod = np.expand_dims(np.nanmean(mod, axis=-1), axis=-1)
            std_mod = np.expand_dims(np.nanstd(mod, axis=-1), axis=-1)
            # to avoid ZeroDivisionError, replace all stds of 0 with NaN
            std_obs[std_obs == 0.0] = np.nan
            std_mod[std_mod == 0.0] = np.nan
            standard_score_obs = (obs - mean_obs) / std_obs
            standard_score_mod = (mod - mean_mod) / std_mod
            standard_score_mult = standard_score_obs * standard_score_mod
            return nansumwrapper(standard_score_mult, axis=-1) / np.count_nonzero(~np.isnan(standard_score_mult), axis=-1)

    @staticmethod
    def calculate_r_squared(obs, mod):
        """
        Calculate the coefficient of determination, r squared, between observations and model
        It is the proportion of the variance in the dependent variable
        that is predictable from the independent variable(s).
        In linear least squares multiple regression with an estimated intercept term,
        the r squared equals the square of the Pearson correlation coefficient.

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            r2
        """

        if obs.size == 0:
            return np.nan
        else:
            return ModBias.calculate_r(obs, mod) ** 2

    @staticmethod
    def calculate_fac2(obs, mod):
        """
        Calculate fraction of model values within
        a factor of two of observed values (FAC2)

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            FAC2
        """

        if obs.size == 0:
            return np.nan
        else:
            # to avoid ZeroDivisionError, replace all obs of 0 with NaN
            obs_nan = copy.deepcopy(obs)
            obs_nan[obs_nan == 0.0] = np.nan
            frac = mod / obs_nan
            n = np.count_nonzero(~np.isnan(frac), axis=-1)
            return (100.0 / n) * nansumwrapper(((frac >= 0.5) & (frac <= 2.0)), axis=-1)

    @staticmethod
    def calculate_upa(obs, mod):
        """
        Calculate unpaired peak accuracy (UPA).
        See here: https://gitlab.com/polyphemus/atmopy/-/blob/master/stat/measure.py.

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
            
        Returns
        -------
        float
            UPA
        """

        if obs.size == 0:
            return np.nan
        else:
            obs_max = np.nanmax(obs, axis=-1)
            mod_max = np.nanmax(mod, axis=-1)
            return ((mod_max - obs_max) / obs_max) * 100.0

    @staticmethod
    def calculate_fairmode_stats(obs, mod, u_95r_RV, RV, alpha, beta, exc_threshold, percentile, plot):
        """ 
        Calculate FAIRMODE statistics
        See here: https://publications.jrc.ec.europa.eu/repository/handle/JRC129254

        Parameters
        ----------
        obs : numpy.array
            Observations data array
        mod : numpy.array
            Model data array
        u_95r_RV : float
            Uncertainty around the reference value
        RV : float
            Reference value
        alpha : float
            Uncertainty parameter
        beta : float
            Proportionality coefficient
        exc_threshold : int
            Threshold used to calculate the exceedances
        percentile : float
            Percentile
        plot : str
            Differenciates between fairmode target and statsummary
        """

        is_finite = np.isfinite(obs+mod)

        if np.any(is_finite):

            # remove missing data
            obs, mod = obs[is_finite], mod[is_finite]

            # calculate RMSu (root mean squared uncertainty)
            rms_u = Stats.calculate_rms_u(obs, u_95r_RV, RV, alpha)

            # calculate Mean Bias
            bias = ModBias.calculate_mb(obs, mod)
            
            # fairmode target plot
            if plot == 'target':
                        
                # calculate x-axis values (CRMSE/BETA*RMSu)
                crmse = ModBias.calculate_crmse(obs, mod)
                x = crmse / (beta * rms_u)
                
                # calculate y-axis values (Mean Bias/BETA*RMSu)
                y = bias / (beta * rms_u)
                
                # calculate ratio
                ratio = np.abs(
                    (Stats.calculate_standard_deviation(mod) - Stats.calculate_standard_deviation(obs))) / (
                        Stats.calculate_standard_deviation(obs) * np.sqrt(2 * (1 - ModBias.calculate_r(obs, mod))))
                
                # For ratios larger than one the σ error dominates and 
                # the station is represented on the right, whereas the reverse
                # applies for values smaller than one
                if ratio < 1:
                    x *= -1

                # calculate Modeling Quality Indicator (MQI)
                rmse = ModBias.calculate_rmse(obs, mod)
                mqi = rmse / (beta * rms_u)

                return x, y, mqi
            
            # fairmode summarystats plot
            elif plot == 'summary':

                # calculate exceedance
                exc = Stats.calculate_exceedances(obs,exc_threshold) if exc_threshold != None else None

                # calculate mean
                mean = Stats.calculate_mean(obs)
                
                # calculate Observation and Model Percentile
                obs_perc = Stats.calculate_percentile(obs, percentile=percentile)
                mod_perc = Stats.calculate_percentile(mod, percentile=percentile)

                # calculate Temporal Statistic for Bias
                t_bias = bias / (beta * rms_u)
                
                # calculate Temporal Statistic for Correlation
                t_R = (1 - ModBias.calculate_r(obs, mod)) / ((0.5 * (beta ** 2) * rms_u ** 2) / (Stats.calculate_standard_deviation(obs) * Stats.calculate_standard_deviation(mod)))
                
                # calculate Temporal Statistic for Standard Deviation
                t_sd = (np.abs(Stats.calculate_standard_deviation(mod) - Stats.calculate_standard_deviation(obs))) / (beta * rms_u)

                # calculate Uncertainty
                U = u_95r_RV * np.sqrt((1 - alpha ** 2) * obs_perc ** 2 + alpha ** 2 * RV ** 2)

                # calculate High Percentile
                h_perc = np.abs(mod_perc - obs_perc) / (beta * U)

                return mean, exc, t_bias, t_R, t_sd, h_perc

        else:
            if plot == 'target':
                return np.nan, np.nan, np.nan
            elif plot == 'summary':
                return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        