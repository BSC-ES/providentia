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
        "taylor": [
            "MONARCH | r: 0.69, StdDev: -2.11"
        ]
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
        "taylor": [
            "MONARCH | r: 0.89, StdDev: -0.09"
        ]
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
        "taylor": [
            "MONARCH | r: 0.68, StdDev: -1.43"
        ]
    }),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost",
     {
        "timeseries": [
            "observations | Mean: 0.05",
            "osuite-global-000 | Mean: 0.02, MB: -0.03, RMSE: 0.07, r: 0.64",
            "cntrl-global-000 | Mean: 0.02, MB: -0.03, RMSE: 0.08, r: 0.61",
            "icap-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.08, r: 0.49"
        ],
        "distribution": [
            "observations | Min: 0.00, Max: 1.94",
            "osuite-global-000 | Min: 0.00, Max: 1.08",
            "cntrl-global-000 | Min: 0.00, Max: 1.18",
            "icap-global-000 | Min: 0.00, Max: 1.14"
        ],
        "taylor": [
            "osuite-global-000 | r: 0.64, StdDev: -0.02",
            "cntrl-global-000 | r: 0.61, StdDev: -0.02",
            "icap-global-000 | r: 0.49, StdDev: -0.04"
        ]
    }
    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost", {
        "timeseries": [
            "observations | Mean: 0.03",
            "osuite-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.55",
            "cntrl-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.47",
            "icap-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.36"
        ],
        "distribution": [
            "observations | Min: 0.01, Max: 0.08",
            "osuite-global-000 | Min: 0.00, Max: 0.01",
            "cntrl-global-000 | Min: 0.00, Max: 0.01",
            "icap-global-000 | Min: 0.00, Max: 0.04"
        ],
        "taylor": [
            "osuite-global-000 | r: 0.55, StdDev: -0.01",
            "cntrl-global-000 | r: 0.47, StdDev: -0.01",
            "icap-global-000 | r: 0.36, StdDev: -0.01"
        ]
    }),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost",
     {
        "timeseries": [
            "observations | Mean: 0.03",
            "osuite-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.32",
            "cntrl-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.25",
            "icap-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.18"
        ],
        "distribution": [
            "observations | Min: 0.00, Max: 0.12",
            "osuite-global-000 | Min: 0.00, Max: 0.01",
            "cntrl-global-000 | Min: 0.00, Max: 0.01",
            "icap-global-000 | Min: 0.00, Max: 0.02"
        ],
        "taylor": [
            "osuite-global-000 | r: 0.32, StdDev: -0.02",
            "cntrl-global-000 | r: 0.25, StdDev: -0.02",
            "icap-global-000 | r: 0.18, StdDev: -0.01"
        ]
    })
]


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_timeseries(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'timeseries', [
              'annotate', 'smooth'], expected_annotations)


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
              'domain'], expected_annotations)


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_taylor(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type,
              'taylor-r', ['annotate'], expected_annotations)
