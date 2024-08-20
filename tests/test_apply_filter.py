from .aux_functions import check_filter_data
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


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[0:3])
def test_apply_period(inst, statistic_mode, network_type):
    inst.reset_filter(initialise=True)
    inst.apply_filter('period', keep='Daytime')
    check_filter_data(inst, statistic_mode, network_type, filter='period')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_representativity(inst, statistic_mode, network_type):
    inst.reset_filter(initialise=True)
    if network_type == 'ghost':
        inst.apply_filter('all_representativity_percent', limit=50)
    else:
        inst.apply_filter('all_representativity_percent', limit=20)
    check_filter_data(inst, statistic_mode, network_type,
                      filter='representativity')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_bounds(inst, statistic_mode, network_type):
    inst.reset_filter(initialise=True)
    inst.apply_filter('latitude', lower=50, upper=60)
    check_filter_data(inst, statistic_mode, network_type, filter='bounds')


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_metadata(inst, statistic_mode, network_type):
    inst.reset_filter(initialise=True)
    if network_type == 'ghost':
        value = ['AT0034G_UVP']
    else:
        value = ['Barcelona']
    inst.apply_filter('station_reference', keep=value)
    check_filter_data(inst, statistic_mode, network_type, filter='keep')
    inst.apply_filter('station_reference', remove=value)
    check_filter_data(inst, statistic_mode, network_type, filter='remove')
