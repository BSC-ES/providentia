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


@ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_read_data(inst, statistic_mode, network_type):
    read_data(inst, statistic_mode, network_type)
