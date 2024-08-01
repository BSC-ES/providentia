from .aux_functions import read_data, make_plot
from providentia import Interactive

# initialise class (Spatial | Temporal)
network_type = "nonghost"
test_name = "local_spatial_temporal"
inst = Interactive(conf='tests_nonghost.conf', statistic_mode="Spatial|Temporal",
                   statistic_aggregation="Median")
expected_annotations = {
    "timeseries": [
        "observations | Mean: 0.04", 
        "osuite-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.23", 
        "cntrl-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.21", 
        "icap-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.04, r: 0.14"],
    "distribution": [
        "observations | Min: 0.00, Max: 0.11", 
        "osuite-global-000 | Min: 0.00, Max: 0.01", 
        "cntrl-global-000 | Min: 0.00, Max: 0.01", 
        "icap-global-000 | Min: 0.00, Max: 0.03"],
}

def test_read_data():
    read_data(inst, test_name, network_type)


def test_make_timeseries():
    make_plot(inst, test_name, network_type, 'timeseries', [
              'annotate'], expected_annotations)


def test_make_distribution():
    make_plot(inst, test_name, network_type, 'distribution', [
              'annotate'], expected_annotations)


def test_make_statsummary():
    stats = ["Mean", "StdDev", "Median", "Var", "Min", "Max",
             "NData", "Data%", "Exceedances", "p1", "p5", "p10",
             "p25", "p75", "p90", "p95", "p99"]
    make_plot(inst, test_name, network_type, 'statsummary', stats=stats)

    stats = ["MB", "NMB", "ME", "NME", "MNB", "MNE", "MFB",
             "MFE", "RMSE", "NRMSE", "COE", "FAC2", "IOA", "r", "r2", "UPA"]
    make_plot(inst, test_name, network_type,
              'statsummary', ['bias'], stats=stats)
