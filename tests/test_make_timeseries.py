from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "ghost", [
        "observations | Mean: 34.56",
        "MONARCH | Mean: 29.55, MB: -5.00, RMSE: 11.64, r: 0.69"
    ]
    ),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "ghost", [
        "observations | Mean: 34.57",
        "MONARCH | Mean: 29.60, MB: -4.97, RMSE: 6.10, r: 0.89"
    ]
    ),
    (Interactive(conf='tests_ghost.conf'),
     "temporal_spatial", "ghost", [
        "observations | Mean: 32.91",
        "MONARCH | Mean: 27.88, MB: -5.45, RMSE: 10.97, r: 0.68"
    ]),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost", [
        "observations | Mean: 0.05",
        "osuite-global-000 | Mean: 0.02, MB: -0.03, RMSE: 0.07, r: 0.64",
        "cntrl-global-000 | Mean: 0.02, MB: -0.03, RMSE: 0.08, r: 0.61",
        "icap-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.08, r: 0.49"
    ]
    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost", [
        "observations | Mean: 0.03",
        "osuite-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.55",
        "cntrl-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.47",
        "icap-global-000 | Mean: 0.00, MB: -0.03, RMSE: 0.03, r: 0.36"
    ]),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost",
     [
        "observations | Mean: 0.03",
        "osuite-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.32",
        "cntrl-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.25",
        "icap-global-000 | Mean: 0.00, MB: -0.02, RMSE: 0.03, r: 0.18"
    ])
]


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_timeseries(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'timeseries', [
              'annotate', 'smooth'], expected_annotations)
