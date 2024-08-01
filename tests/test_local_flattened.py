from .aux_functions import read_data, make_plot
from providentia import Interactive

# initialise class (Flattened)
test_name = "local_flattened"
inst = Interactive(conf='local.conf', statistic_mode="Flattened", statistic_aggregation="")
expected_annotations = {
    "timeseries": ["EBAS | Mean: 34.56", "MONARCH | Mean: 29.55, MB: -5.00, RMSE: 11.64, r: 0.69"],
    "distribution": ["EBAS | Min: 0.00, Max: 119.35", "MONARCH | Min: 0.00, Max: 106.75"],
}

def test_read_data():
    read_data(inst, test_name)

def test_make_timeseries():
    make_plot(inst, test_name, 'timeseries', ['annotate'], expected_annotations)

def test_make_distribution():
    make_plot(inst, test_name, 'distribution', ['annotate'], expected_annotations)

def test_make_statsummary():
    make_plot(inst, test_name, 'statsummary')