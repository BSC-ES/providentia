from .aux_functions import check_report
import providentia as prv
import pytest

possibilities = [
    (prv.Providentia("tests_ghost.conf", tests=True), "ghost"),
    (prv.Providentia("tests_nonghost.conf", tests=True), "nonghost"),
]


@pytest.mark.parametrize("inst, network_type", possibilities)
def test_statistics(inst, network_type):
    inst.report()
    check_report(network_type)
