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


@ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_period(inst, statistic_mode, network_type):
    inst.apply_filter('period', keep='Daytime')
    check_filter_data(inst, statistic_mode, network_type, filter='period')


# @ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
# def test_apply_representativity(inst, statistic_mode, network_type):
#     inst.apply_filter('all_representativity_percent', limit=50)
#     check_filter_data(inst, statistic_mode, network_type,
#                       filter='representativity')


# @ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
# def test_apply_bounds(inst, statistic_mode, network_type):
#     inst.apply_filter('latitude', lower=50, upper=60)
#     check_filter_data(inst, statistic_mode, network_type, filter='bounds')


# @ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
# def test_apply_keep(inst, statistic_mode, network_type):
#     if network_type == 'ghost':
#         keep = ['AR0001R_UVP']
#     else:
#         keep = ['Barcelona']
#     inst.apply_filter('station_reference', keep=keep)
#     check_filter_data(inst, statistic_mode, network_type, filter='keep')


# @ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
# def test_apply_remove(inst, statistic_mode, network_type):
#     if network_type == 'ghost':
#         remove = ['Donon']
#     else:
#         remove = ['ATHENS-NOA', 'Zaragoza']
#     inst.apply_filter('station_name', remove=remove)
#     check_filter_data(inst, statistic_mode, network_type, filter='remove')
