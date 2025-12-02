from .aux_functions import read_data
import providentia as prv
import pytest


possibilities = [
    (prv.Providentia('tests_ghost.conf',
                     statistic_mode="Flattened",
                     statistic_aggregation="",
                     tests=True),
     "flattened", "ghost"),
    (prv.Providentia('tests_ghost.conf',
                     statistic_mode="Spatial|Temporal",
                     statistic_aggregation="Median",
                     tests=True),
     "spatial_temporal", "ghost"),
    (prv.Providentia('tests_ghost.conf',
                     tests=True),
     "temporal_spatial", "ghost"),
    (prv.Providentia('tests_nonghost.conf',
                     statistic_mode="Flattened",
                     statistic_aggregation="",
                     tests=True),
     "flattened", "nonghost"),
    (prv.Providentia('tests_nonghost.conf',
                     statistic_mode="Spatial|Temporal",
                     statistic_aggregation="Median",
                     tests=True),
     "spatial_temporal", "nonghost"),
    (prv.Providentia('tests_nonghost.conf',
                     tests=True),
     "temporal_spatial", "nonghost"),
    (prv.Providentia('tests_nonghost_calibration.conf',
                     statistic_mode="Flattened",
                     statistic_aggregation="",
                     tests=True),
     "flattened", "nonghost"),
    (prv.Providentia('tests_nonghost_calibration.conf',
                     statistic_mode="Spatial|Temporal",
                     statistic_aggregation="Median",
                     tests=True),
     "spatial_temporal", "nonghost"),
    (prv.Providentia('tests_nonghost_calibration.conf',
                     tests=True),
     "temporal_spatial", "nonghost")
]


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[0:6])
def test_read_data(inst, statistic_mode, network_type):
    path = f'tests/reference/{network_type}/{statistic_mode}/data/data.npy'
    inst.load()
    read_data(inst, path, network_type)


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[6:9])
def test_calibration(inst, statistic_mode, network_type):
    path = f'tests/reference/{network_type}/{statistic_mode}/data/data_calibration.npy'
    inst.load()
    read_data(inst, path, network_type)
