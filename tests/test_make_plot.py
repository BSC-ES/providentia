from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation="",
                 tests=True),
     "flattened", "ghost"),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median",
                 tests=True),
     "spatial_temporal", "ghost"),
    (Interactive(conf='tests_ghost.conf',
                 tests=True),
     "temporal_spatial", "ghost"),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation="",
                 tests=True),
     "flattened", "nonghost"),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median",
                 tests=True),
     "spatial_temporal", "nonghost"),
    (Interactive(conf='tests_nonghost.conf',
                 tests=True),
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
    make_plot(inst, statistic_mode, network_type, 'statsummary')
    make_plot(inst, statistic_mode, network_type, 'statsummary', ['bias'])


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


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_table(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'table-RMSE')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_periodic(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'periodic-r')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_periodic_violin(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'periodic-violin')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_make_scatter(inst, statistic_mode, network_type):
    make_plot(inst, statistic_mode, network_type, 'scatter', ['regression'])
