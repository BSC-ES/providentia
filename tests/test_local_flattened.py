from .aux import read_data, make_plot
from providentia import Interactive

# initialise class (Flattened)
test_name = "local_flattened"
inst = Interactive(conf='local.conf', statistic_mode="Flattened", statistic_aggregation="")
expected_annotations = {
    "timeseries": ["EBAS | Mean: 34.56", "MONARCH | Mean: 29.55, MB: -5.00, RMSE: 11.64, r: 0.69"],
    "distribution": ["Min: 0.00, Max: 119.35", "MONARCH | Min: 0.00, Max: 106.75"],
    "scatter": [""]
}

def test_read_data():
    read_data(inst, test_name)

def test_make_timeseries():
    make_plot(inst, 'timeseries', expected_annotations, test_name)

def test_make_distribution():
    make_plot(inst, 'distribution', expected_annotations, test_name)

# def test_make_scatter():
#    make_plot('scatter')
