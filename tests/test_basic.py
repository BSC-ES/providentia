import matplotlib
import numpy as np

from providentia import Interactive

# initialise class
inst = Interactive(conf='local.conf')

expected_annotations = {
    "timeseries": ["EBAS | Mean: 32.91", "MONARCH | Mean: 27.88, MB: -5.45, RMSE: 10.97, r: 0.68"],
    "distribution": ["EBAS | Min: 1.78, Max: 80.38", "MONARCH | Min: 0.22, Max: 69.92"],
    "scatter": ["MONARCH | r2: 0.46, RMSE: 10.97"]
}

def test_data_read():

    # read expected output
    expected_output = np.load('tests/data/EBAS_sconco3_data.npy', allow_pickle=True)
    
    # get data in memory in xarray format
    data = inst.get_data(format='xr')
    generated_output = data['EBAS|sconco3_data'].values
    # np.save('tests/data/EBAS_sconco3_data', generated_output)

    # replace nans by -999
    expected_output = np.nan_to_num(expected_output, copy=True, nan=-999, posinf=None, neginf=None)
    generated_output = np.nan_to_num(generated_output, copy=True, nan=-999, posinf=None, neginf=None)

    assert (np.array_equal(generated_output, expected_output))

def make_plot(plot_type):

    # make plot
    fig = inst.make_plot(plot_type, annotate=True, return_plot=True)

    # check that a figure has been returned
    assert (type(fig) == matplotlib.figure.Figure)

    # check if annotations are correct
    annotations = [child for child in fig.axes[0].get_children() 
                   if type(child) == matplotlib.offsetbox.AnchoredOffsetbox][0]
    print(annotations.get_child().get_children()[0].get_text())
    for annotation, expected_annotation in zip(annotations.get_child().get_children(), 
                                               expected_annotations[plot_type]):
        assert annotation.get_text() == expected_annotation

    for line_i, line in enumerate(fig.axes[0].lines):
        
        #np.save(f'tests/data/{plot_type}_line_{line_i}', line.get_xydata())
            
        # read expected output
        expected_output = np.load(f'tests/data/{plot_type}_line_{line_i}.npy', allow_pickle=True)

        # get data
        generated_output = line.get_xydata()

        # replace nans by -999
        expected_output = np.nan_to_num(expected_output, copy=True, nan=-999, posinf=None, 
                                        neginf=None)
        generated_output = np.nan_to_num(generated_output, copy=True, nan=-999, posinf=None, 
                                         neginf=None)

        # check data for each line is correct
        assert (np.array_equal(generated_output, expected_output))

def test_make_timeseries():
    make_plot('timeseries')

def test_make_distribution():
    make_plot('distribution')

#def test_make_scatter():
#    make_plot('scatter')
