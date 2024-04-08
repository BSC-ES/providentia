import matplotlib
import numpy as np

from providentia import Interactive

# initialise class
inst = Interactive(conf='local.conf')


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

