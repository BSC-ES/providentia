from .aux_functions import read_data
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
     "temporal_spatial", "nonghost"),
    (prv.load('tests_nonghost_calibration.conf',
              statistic_mode="Flattened",
              statistic_aggregation="",
              tests=True),
     "flattened", "nonghost"),
    (prv.load('tests_nonghost_calibration.conf',
              statistic_mode="Spatial|Temporal",
              statistic_aggregation="Median",
              tests=True),
     "spatial_temporal", "nonghost"),
    (prv.load('tests_nonghost_calibration.conf',
              tests=True),
     "temporal_spatial", "nonghost")
]


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[0:6])
def test_read_data(inst, statistic_mode, network_type):
    path = f'tests/reference/{network_type}/{statistic_mode}/data/data.npy'
    read_data(inst, path)


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[6:9])
def test_calibration(inst, statistic_mode, network_type):
    path = f'tests/reference/{network_type}/{statistic_mode}/data/data_calibration.npy'
    read_data(inst, path)
