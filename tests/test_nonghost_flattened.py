from .aux_functions import read_data, make_plot
from providentia import Interactive

# initialise class (Flattened)
network_type = "nonghost"
test_name = "local_flattened"
inst = Interactive(conf='tests_nonghost.conf',
                   statistic_mode="Flattened", statistic_aggregation="")
expected_annotations = {
    "timeseries": [
        "observations | Mean: 0.05", 
        "osuite-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.05, r: 0.12", 
        "cntrl-global-000 | Mean: 0.00, MB: -0.04, RMSE: 0.05, r: 0.03", 
        "icap-global-000 | Mean: 0.01, MB: -0.04, RMSE: 0.05, r: 0.10"],
    "distribution": [
        "observations | Min: 0.00, Max: 0.38", 
        "osuite-global-000 | Min: 0.00, Max: 0.09", 
        "cntrl-global-000 | Min: 0.00, Max: 0.10", 
        "icap-global-000 | Min: 0.00, Max: 0.12"],
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
