
from .aux_functions import read_data, make_plot
from providentia import Interactive
import pytest


possibilities = [
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "ghost", [
        "MONARCH | r: 0.69, StdDev: -2.11"
    ]

    ),
    (Interactive(conf='tests_ghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "ghost", [
        "MONARCH | r: 0.89, StdDev: -0.09"
    ]

    ),
    (Interactive(conf='tests_ghost.conf'),
     "temporal_spatial", "ghost",
     [
        "MONARCH | r: 0.68, StdDev: -1.43"
    ]
    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Flattened",
                 statistic_aggregation=""),
     "flattened", "nonghost", [
        "osuite-global-000 | r: 0.12, StdDev: -0.02",
        "cntrl-global-000 | r: 0.03, StdDev: -0.02",
        "icap-global-000 | r: 0.10, StdDev: -0.02"
    ]

    ),
    (Interactive(conf='tests_nonghost.conf',
                 statistic_mode="Spatial|Temporal",
                 statistic_aggregation="Median"),
     "spatial_temporal", "nonghost", [
        "osuite-global-000 | r: 0.23, StdDev: -0.01",
        "cntrl-global-000 | r: 0.21, StdDev: -0.01",
        "icap-global-000 | r: 0.14, StdDev: -0.01"
    ]
    ),
    (Interactive(conf='tests_nonghost.conf'),
     "temporal_spatial", "nonghost",
     [
        "osuite-global-000 | r: 0.12, StdDev: -0.02",
        "cntrl-global-000 | r: 0.14, StdDev: -0.02",
        "icap-global-000 | r: 0.11, StdDev: -0.02"
    ]
    )
]


@pytest.mark.parametrize("inst, statistic_mode, network_type, expected_annotations", possibilities)
def test_make_taylor(inst, statistic_mode, network_type, expected_annotations):
    make_plot(inst, statistic_mode, network_type,
              'taylor-r', ['annotate'], expected_annotations)
