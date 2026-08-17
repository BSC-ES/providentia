from .aux_functions import save_data
import providentia as prv
import pytest


possibilities = [
    (
        prv.Providentia(
            "tests/tests_ghost.conf",
            statistic_mode="Flattened",
            statistic_aggregation="",
            tests=True,
        ),
        "flattened",
        "ghost",
    ),
    (
        prv.Providentia(
            "tests/tests_ghost.conf",
            statistic_mode="Spatial|Temporal",
            statistic_aggregation="Median",
            tests=True,
        ),
        "spatial_temporal",
        "ghost",
    ),
    (prv.Providentia("tests/tests_ghost.conf", tests=True), "temporal_spatial", "ghost"),
    (
        prv.Providentia(
            "tests/tests_nonghost.conf",
            statistic_mode="Flattened",
            statistic_aggregation="",
            tests=True,
        ),
        "flattened",
        "nonghost",
    ),
    (
        prv.Providentia(
            "tests/tests_nonghost.conf",
            statistic_mode="Spatial|Temporal",
            statistic_aggregation="Median",
            tests=True,
        ),
        "spatial_temporal",
        "nonghost",
    ),
    (
        prv.Providentia("tests/tests_nonghost.conf", tests=True),
        "temporal_spatial",
        "nonghost",
    ),
]


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_save_conf(inst, statistic_mode, network_type):
    inst.load()
    save_data(inst, "conf", "data", network_type, statistic_mode)


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_save_npz(inst, statistic_mode, network_type):
    inst.load()
    save_data(inst, "npz", "data", network_type, statistic_mode)


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_save_netcdf(inst, statistic_mode, network_type):
    inst.load()
    save_data(inst, "nc", "data", network_type, statistic_mode)
