from .aux_functions import read_data, make_plot
from providentia import Interactive

# initialise class (Spatial|Temporal)
test_name = "local_spatial_temporal"
inst = Interactive(conf='local.conf', statistic_mode="Spatial|Temporal",
                   statistic_aggregation="Median")
expected_annotations = {
    "timeseries": ["EBAS | Mean: 34.57", "MONARCH | Mean: 29.60, MB: -4.97, RMSE: 6.10, r: 0.89"],
    "distribution": ["EBAS | Min: 17.00, Max: 57.78", "MONARCH | Min: 12.63, Max: 52.28"],
}


def test_read_data():
    read_data(inst, test_name)


def test_make_timeseries():
    make_plot(inst, test_name, 'timeseries', [
              'annotate'], expected_annotations)


def test_make_distribution():
    make_plot(inst, test_name, 'distribution', [
              'annotate'], expected_annotations)


def test_make_statsummary():
    stats = ["Mean", "StdDev", "Median", "Var", "Min", "Max",
             "NData", "Data%", "Exceedances", "p1", "p5", "p10",
             "p25", "p75", "p90", "p95", "p99"]
    make_plot(inst, test_name, 'statsummary', stats=stats)

    stats = ["MB", "NMB", "ME", "NME", "MNB", "MNE", "MFB",
             "MFE", "RMSE", "NRMSE", "COE", "FAC2", "IOA", "r", "r2", "UPA"]
    make_plot(inst, test_name, 'statsummary', ['bias'], stats=stats)
