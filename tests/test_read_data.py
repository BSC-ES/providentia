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
     "temporal_spatial", "nonghost")
]


@ pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_read_data(inst, statistic_mode, network_type):
    read_data(inst, statistic_mode, network_type)
