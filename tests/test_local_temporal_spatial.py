from .aux import read_data, make_plot
from providentia import Interactive

# initialise class (Temporal|Spatial)
inst = Interactive(conf='local.conf')
expected_annotations = {
    "timeseries": ["EBAS | Mean: 32.91", "MONARCH | Mean: 27.88, MB: -5.45, RMSE: 10.97, r: 0.68"],
    "distribution": ["EBAS | Min: 1.78, Max: 80.38", "MONARCH | Min: 0.22, Max: 69.92"],
    "scatter": ["MONARCH | r2: 0.46, RMSE: 10.97"]
}

def test_read_data():
    read_data(inst)

def test_make_timeseries():
    make_plot(inst, 'timeseries', expected_annotations)

def test_make_distribution():
    make_plot(inst, 'distribution', expected_annotations)

# def test_make_scatter():
#    make_plot('scatter')
