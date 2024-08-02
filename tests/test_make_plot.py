from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "ghost",
     {
        "timeseries": [
            "observations | Mean: 34.56",
            "MONARCH | Mean: 29.55, MB: -5.00, RMSE: 11.64, r: 0.69"
        ],
        "distribution": [
            "observations | Min: 0.00, Max: 119.35",
            "MONARCH | Min: 0.00, Max: 106.75"
        ],
    }
    ),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "ghost",
     {
        "timeseries": [
            "observations | Mean: 34.57",
            "MONARCH | Mean: 29.60, MB: -4.97, RMSE: 6.10, r: 0.89"
        ],
        "distribution": [
            "observations | Min: 17.00, Max: 57.78",
            "MONARCH | Min: 12.63, Max: 52.28"
        ],
    }
    ),
    (Interactive(conf='tests_ghost.conf'),
     "temporal_spatial", "ghost",
     {
        "timeseries": [
            "observations | Mean: 32.91",
            "MONARCH | Mean: 27.88, MB: -5.45, RMSE: 10.97, r: 0.68"
        ],
        "distribution": [
            "observations | Min: 1.78, Max: 80.38",
            "MONARCH | Min: 0.22, Max: 69.92"
        ],
    }),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost",
     {
        "timeseries": [
            "observations | Mean: 0.05",
            "osuite-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.05, r: 0.12",
            "cntrl-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.05, r: 0.03",
            "icap-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.05, r: 0.10"
        ],
        "distribution": [
            "observations | Min: 0.00, Max: 0.38",
            "osuite-global-000 | Min: 0.00, Max: 0.09",
            "cntrl-global-000 | Min: 0.00, Max: 0.10",
            "icap-global-000 | Min: 0.00, Max: 0.12"
        ],
    }
    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost", {
        "timeseries": [
            "observations | Mean: 0.04",
            "osuite-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.23",
            "cntrl-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.21",
            "icap-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.14"
        ],
        "distribution": [
            "observations | Min: 0.00, Max: 0.11",
            "osuite-global-000 | Min: 0.00, Max: 0.01",
            "cntrl-global-000 | Min: 0.00, Max: 0.01",
            "icap-global-000 | Min: 0.00, Max: 0.03"
        ],
    }),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost",
     {
        "timeseries": [
            "observations | Mean: 0.05",
            "osuite-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.05, r: 0.12",
            "cntrl-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.05, r: 0.14",
            "icap-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.05, r: 0.11"
        ],
        "distribution": [
            "observations | Min: 0.01, Max: 0.16",
            "osuite-global-000 | Min: 0.00, Max: 0.01",
            "cntrl-global-000 | Min: 0.00, Max: 0.01",
            "icap-global-000 | Min: 0.00, Max: 0.02"
        ],
    })
]


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_timeseries(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'timeseries', [
              'annotate'], expected_annotations)


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_distribution(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'distribution', [
              'annotate'], expected_annotations)


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_statsummary(inst, statistic_mode, network_type, expected_annotations):
    basic_stats = ["Mean", "StdDev", "Median", "Var", "Min", "Max",
                   "NData", "Data%", "Exceedances", "p1", "p5", "p10",
                   "p25", "p75", "p90", "p95", "p99"]
    make_plot(inst, statistic_mode, network_type,
              'statsummary', stats=basic_stats)

    expbias_stats = ["MB", "NMB", "ME", "NME", "MNB", "MNE", "MFB",
                     "MFE", "RMSE", "NRMSE", "COE", "FAC2", "IOA", "r", "r2", "UPA"]
    make_plot(inst, statistic_mode, network_type,
              'statsummary', ['bias'], stats=expbias_stats)


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_map(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'map-Median', [
              'annotate'], expected_annotations)