import warnings
from .aux_functions import check_filter_data
import providentia as prv
import pytest


possibilities = [
    (
        prv.Providentia(
            "tests_ghost.conf",
            statistic_mode="Flattened",
            statistic_aggregation="",
            tests=True,
        ),
        "flattened",
        "ghost",
    ),
    (
        prv.Providentia(
            "tests_ghost.conf",
            statistic_mode="Spatial|Temporal",
            statistic_aggregation="Median",
            tests=True,
        ),
        "spatial_temporal",
        "ghost",
    ),
    (prv.Providentia("tests_ghost.conf", tests=True), "temporal_spatial", "ghost"),
    (
        prv.Providentia(
            "tests_nonghost.conf",
            statistic_mode="Flattened",
            statistic_aggregation="",
            tests=True,
        ),
        "flattened",
        "nonghost",
    ),
    (
        prv.Providentia(
            "tests_nonghost.conf",
            statistic_mode="Spatial|Temporal",
            statistic_aggregation="Median",
            tests=True,
        ),
        "spatial_temporal",
        "nonghost",
    ),
    (
        prv.Providentia("tests_nonghost.conf", tests=True),
        "temporal_spatial",
        "nonghost",
    ),
]


@pytest.fixture(autouse=True)
def suppress_warnings():
    # Hide runtime errors coming from calculating statistics with nans
    # These nans appear because sometimes we don't have data for all stations
    # In the whole period of time
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        yield


# Can only filter by period for GHOST network
@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities[0:3])
def test_apply_period(inst, statistic_mode, network_type):
    inst.load()

    # Filter by keeping Daytime
    inst.filter("period", keep="Daytime")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Daytime")

    # Filter by keeping Nighttime
    inst.reset(initialise=True)
    inst.filter("period", keep="Nighttime")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Nighttime")

    # Filter by removing Daytime
    inst.filter("period", remove="Daytime")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Daytime")

    # Filter by removing Nighttime
    inst.reset(initialise=True)
    inst.filter("period", remove="Nighttime")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Nighttime")

    # Filter by keeping Weekday
    inst.reset(initialise=True)
    inst.filter("period", keep="Weekday")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Weekday")

    # Filter by removing Weekday
    inst.reset(initialise=True)
    inst.filter("period", remove="Weekday")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Weekday")

    # Filter by keeping Weekend
    inst.reset(initialise=True)
    inst.filter("period", keep="Weekend")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Weekend")

    # Filter by removing Weekend
    inst.reset(initialise=True)
    inst.filter("period", remove="Weekend")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Weekend")

    # Filter by keeping Spring
    inst.reset(initialise=True)
    inst.filter("period", keep="Spring")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Spring")

    # Filter by removing Spring
    inst.reset(initialise=True)
    inst.filter("period", remove="Spring")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Spring")

    # Filter by keeping Winter
    inst.reset(initialise=True)
    inst.filter("period", keep="Winter")
    check_filter_data(inst, statistic_mode, network_type, filter="keep_Winter")

    # Filter by removing Winter
    inst.reset(initialise=True)
    inst.filter("period", remove="Winter")
    check_filter_data(inst, statistic_mode, network_type, filter="remove_Winter")


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_coverage(inst, statistic_mode, network_type):
    inst.load()
    if network_type == "ghost":
        inst.filter("total_coverage", limit=50)
    else:
        inst.filter("total_coverage", limit=20)
    check_filter_data(inst, statistic_mode, network_type, filter="coverage")


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_coords(inst, statistic_mode, network_type):
    inst.load()
    inst.filter("latitude", lower=50, upper=60)
    inst.filter("longitude", lower=10, upper=12)
    check_filter_data(inst, statistic_mode, network_type, filter="coords")


@pytest.mark.parametrize("inst, statistic_mode, network_type", possibilities)
def test_apply_station_reference(inst, statistic_mode, network_type):
    inst.load()
    if network_type == "ghost":
        value = ["AT0034G_UVP"]
    else:
        value = ["Barcelona"]

    # Filter by keeping station reference
    inst.filter("station_reference", keep=value)
    check_filter_data(
        inst, statistic_mode, network_type, filter="keep_station_reference"
    )

    # Filter by removing station reference
    inst.reset(initialise=True)
    inst.filter("station_reference", remove=value)
    check_filter_data(
        inst, statistic_mode, network_type, filter="remove_station_reference"
    )
