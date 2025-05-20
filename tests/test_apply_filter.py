import warnings
from .aux_functions import check_filter_data
import providentia as prv
import pytest


possibilities = [
    (prv.load('tests_ghost.conf',
              statistic_mode="Flattened",
              statistic_aggregation="",
              tests=True),
     "flattened", "ghost"),
    (prv.load('tests_ghost.conf',
              statistic_mode="Spatial|Temporal",
              statistic_aggregation="Median",
              tests=True),
     "spatial_temporal", "ghost"),
    (prv.load('tests_ghost.conf',
              tests=True),
     "temporal_spatial", "ghost"),
    (prv.load('tests_nonghost.conf',
              statistic_mode="Flattened",
              statistic_aggregation="",
              tests=True),
     "flattened", "nonghost"),
    (prv.load('tests_nonghost.conf',
              statistic_mode="Spatial|Temporal",
              statistic_aggregation="Median",
              tests=True),
     "spatial_temporal", "nonghost"),
    (prv.load('tests_nonghost.conf',
              tests=True),
     "temporal_spatial", "nonghost")
]


@pytest.fixture(autouse=True)
def suppress_warnings():
    # Hide runtime errors coming from calculating statistics with nans
    # These nans appear because sometimes we don't have data for all stations
    # In the whole period of time
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        yield


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[0:3])
def test_apply_period(inst, statistic_mode, network_type):
    inst.reset(initialise=True)
    inst.filter('period', keep='Daytime')
    check_filter_data(inst, statistic_mode, network_type, filter='period')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_representativity(inst, statistic_mode, network_type):
    inst.reset(initialise=True)
    if network_type == 'ghost':
        inst.filter('all_representativity_percent', limit=50)
    else:
        inst.filter('all_representativity_percent', limit=20)
    check_filter_data(inst, statistic_mode, network_type,
                      filter='representativity')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_bounds(inst, statistic_mode, network_type):
    inst.reset(initialise=True)
    inst.filter('latitude', lower=50, upper=60)
    check_filter_data(inst, statistic_mode, network_type, filter='bounds')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_metadata(inst, statistic_mode, network_type):
    inst.reset(initialise=True)
    if network_type == 'ghost':
        value = ['AT0034G_UVP']
    else:
        value = ['Barcelona']
    inst.filter('station_reference', keep=value)
    check_filter_data(inst, statistic_mode, network_type, filter='keep')
    inst.filter('station_reference', remove=value)
    check_filter_data(inst, statistic_mode, network_type, filter='remove')
