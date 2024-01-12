import matplotlib
import numpy as np

from providentia import Interactive

# initialise class
inst = Interactive(conf='test_local.conf')


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

def test_make_timeseries():

    # make timeseries 
    ax = inst.make_plot('timeseries', return_plot=True)

    # check that an axis has been returned
    assert (type(ax) == matplotlib.axes._axes.Axes)

    for line_i, line in enumerate(ax.lines):

        # np.save(f'tests/data/timeseries_line_{line_i}', line.get_xydata())

        # read expected output
        expected_output = np.load(f'tests/data/timeseries_line_{line_i}.npy', allow_pickle=True)

        # get data in timeseries
        generated_output = line.get_xydata()

        # replace nans by -999
        expected_output = np.nan_to_num(expected_output, copy=True, nan=-999, posinf=None, neginf=None)
        generated_output = np.nan_to_num(generated_output, copy=True, nan=-999, posinf=None, neginf=None)

        # check data for each timeseries line is correct
        assert (np.array_equal(generated_output, expected_output))
