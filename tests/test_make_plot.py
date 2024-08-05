from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "ghost"),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "ghost"),
    (Interactive(conf='tests_ghost.conf'),
     "temporal_spatial", "ghost"),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost"),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost"),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost",
     )
]


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_timeseries(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'timeseries', [
              'annotate', 'smooth'])


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_distribution(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'distribution', [
              'annotate'])


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_statsummary(inst, statistic_mode, network_type):
    basic_stats = ["Mean", "StdDev", "Median", "Var", "Min", "Max",
                   "NData", "Data%", "Exceedances", "p1", "p5", "p10",
                   "p25", "p75", "p90", "p95", "p99"]
    make_plot(inst, statistic_mode, network_type,
              'statsummary', stats=basic_stats)

    expbias_stats = ["MB", "NMB", "ME", "NME", "MNB", "MNE", "MFB",
                     "MFE", "RMSE", "NRMSE", "COE", "FAC2", "IOA", "r", "r2", "UPA"]
    make_plot(inst, statistic_mode, network_type,
              'statsummary', ['bias'], stats=expbias_stats)


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_map(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'map-Median', [
              'domain'])


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_taylor(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type,
              'taylor-r', ['annotate'])


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_heatmap(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'heatmap-Median')
