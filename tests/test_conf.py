from .aux_functions import read_data
import providentia as prv
import pytest


valid_sections = [
    ("expID"),
    ("expID-domain"),
    ("expID-ensembleNum"),
    ("expID-forecast"),
    ("expID-domain-ensembleNum"),
    ("expID-domain-forecast"),
    ("expID-ensembleNum-forecast"),
    ("expID_domain"),
    ("expID_domain-ensembleNum"),
    ("expID_domain-forecast"),
    ("expID_domain-ensembleNum-forecast"),
    ("expID_ensembleNum"),
    ("expID_ensembleNum-domain"),
    ("expID_ensembleNum-forecast"),
    ("expID_ensembleNum-domain-forecast"),
    ("expID_domain_Num"),
    ("expID_domain_ensembleNum_forecast"),
    ("expID-expID"),
    ("expID-expID-domain"),
    ("expID-expID-ensembleNum"),
    ("expID-expID-forecast"),
    ("expID-expID-domain-ensembleNum"),
    ("expID-expID-domain-forecast"),
    ("expID-expID-domain-ensembleNum-forecast"),
    ("expID_domain-expID_domain"),
    ("expID_domain-expID_domain-ensembleNum"),
    ("expID_domain-expID_domain-forecast"),
    ("expID_domain-expID_domain-ensembleNum-forecast"),
    ("expID_domain_ensembleNum-expID_domain_ensembleNum"),
    ("expID_domain_ensembleNum_forecast-expID_domain_ensembleNum_forecast"),
    ("expID_domain_ensembleNum"),
    ("expID_ensembleNum_domain"),
]

@pytest.mark.parametrize("section", valid_sections)
def test_valid_config(section):
    path = f"tests/reference/conf/valid/{section}/data.npy"
    inst = prv.Providentia("tests/tests_conf_options_valid.conf", section=section)
    inst.load()
    read_data(inst, path)
