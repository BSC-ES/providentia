import matplotlib
import numpy as np


def read_data(inst, test_name):

    # read expected output
    expected_output = np.load(f'tests/data/{test_name}/EBAS_sconco3_data.npy', allow_pickle=True)
    
    # get data in memory in xarray format
    data = inst.get_data(format='xr')
    generated_output = data['EBAS|sconco3_data'].values
    np.save(f'tests/data/{test_name}/EBAS_sconco3_data', generated_output)

    # replace nans by -999
    expected_output = np.nan_to_num(expected_output, copy=True, nan=-999, posinf=None, neginf=None)
    generated_output = np.nan_to_num(generated_output, copy=True, nan=-999, posinf=None, neginf=None)

    assert (np.array_equal(generated_output, expected_output))


def make_plot(inst, plot_type, expected_annotations, test_name):

    # make plot
    fig = inst.make_plot(plot_type, annotate=True, return_plot=True)

    # check that a figure has been returned
    assert (type(fig) == matplotlib.figure.Figure)

    # check if annotations are correct
    annotations = [child for child in fig.axes[0].get_children() 
                   if type(child) == matplotlib.offsetbox.AnchoredOffsetbox][0]
    
    for annotation, expected_annotation in zip(annotations.get_child().get_children(), 
                                               expected_annotations[plot_type]):
        assert annotation.get_text() == expected_annotation

    for line_i, line in enumerate(fig.axes[0].lines):
        
        np.save(f'tests/data/{test_name}/{plot_type}_line_{line_i}', line.get_xydata())
            
        # read expected output
        expected_output = np.load(f'tests/data/{test_name}/{plot_type}_line_{line_i}.npy', allow_pickle=True)

        # get data
        generated_output = line.get_xydata()

        # replace nans by -999
        expected_output = np.nan_to_num(expected_output, copy=True, nan=-999, posinf=None, 
                                        neginf=None)
        generated_output = np.nan_to_num(generated_output, copy=True, nan=-999, posinf=None, 
                                         neginf=None)

        # check data for each line is correct
        assert (np.array_equal(generated_output, expected_output))
