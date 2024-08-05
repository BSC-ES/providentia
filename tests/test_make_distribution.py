from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "ghost", [
        "observations | Min: 0.00, Max: 119.35",
        "MONARCH | Min: 0.00, Max: 106.75"
    ]
    ),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "ghost", [
        "observations | Min: 17.00, Max: 57.78",
        "MONARCH | Min: 12.63, Max: 52.28"
    ]
    ),
    (Interactive(conf='tests_ghost.conf'),
     "temporal_spatial", "ghost", [
        "observations | Min: 1.78, Max: 80.38",
        "MONARCH | Min: 0.22, Max: 69.92"
    ]),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost", [
        "observations | Min: 0.00, Max: 0.38",
        "osuite-global-000 | Min: 0.00, Max: 0.09",
        "cntrl-global-000 | Min: 0.00, Max: 0.10",
        "icap-global-000 | Min: 0.00, Max: 0.12"
    ]
    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost", [
        "observations | Min: 0.00, Max: 0.11",
        "osuite-global-000 | Min: 0.00, Max: 0.01",
        "cntrl-global-000 | Min: 0.00, Max: 0.01",
        "icap-global-000 | Min: 0.00, Max: 0.03"
    ]),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost", [
        "observations | Min: 0.01, Max: 0.16",
        "osuite-global-000 | Min: 0.00, Max: 0.01",
        "cntrl-global-000 | Min: 0.00, Max: 0.01",
        "icap-global-000 | Min: 0.00, Max: 0.02"
    ])
]


@ pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_distribution(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type, 'distribution', [
              'annotate'], expected_annotations)
