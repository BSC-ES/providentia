import matplotlib
import numpy as np
import pandas as pd
import xarray as xr
from pandas.testing import assert_frame_equal

from providentia.statistics import get_z_statistic_info

GENERATE_OUTPUT = False

def read_data(inst, path):

    # get data in memory in xarray format
    inst.print_config()
    data = inst.data(format='xr')
    var = [v for v in list(data.data_vars) if v.endswith("_data")][0]
    generated_output = data[var].values

    # save data, uncomment if we want to update it
    if GENERATE_OUTPUT:
        np.save(path, generated_output)

    # read expected output
    expected_output = np.load(path, allow_pickle=True)

    assert (np.allclose(generated_output, expected_output, equal_nan=True))

    assert (generated_output.size != 0)


def plot(inst, statistic_mode, network_type, plot_type, plot_options=[]):

    # make plot
    inst.load()
    fig = inst.plot(plot_type, plot_options=plot_options, return_plot=True)

    # check that a figure has been returned
    assert (type(fig) == matplotlib.figure.Figure)

    # get zstat information from plot_type
    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(
        plot_type=plot_type)

    # get base plot type (without stat and options)
    if zstat:
        base_plot_type = plot_type.split('-')[0]
    else:
        base_plot_type = plot_type.split('_')[0]

    if base_plot_type in ['statsummary', 'heatmap', 'table']:

        # get table
        for child in fig.axes[0].get_children():
            if isinstance(child, (matplotlib.table.Table, matplotlib.collections.QuadMesh)):
                table = child
                break

        # extract data from the table/heatmap
        data = []
        if base_plot_type in ['statsummary', 'table']:
            for (row, col), cell in table.get_celld().items():
                data.append({
                    "row": row,
                    "col": col,
                    "value": cell.get_text().get_text()
                })
        else:
            for (x, y), value in np.ndenumerate(table.get_array()):
                data.append({
                    "x": x,
                    "y": y,
                    "value": value
                })
        generated_output = pd.DataFrame(data)

        # save data, uncomment if we want to update it
        if 'bias' in plot_options:
            path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{plot_type}_bias_values.csv'
        else:
            path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{plot_type}_values.csv'
        if GENERATE_OUTPUT:
            generated_output.to_csv(path, index=False)

        # read expected output
        expected_output = pd.read_csv(path, keep_default_na=False)

        assert assert_frame_equal(generated_output, expected_output) is None

    elif base_plot_type in ['timeseries', 'distribution', 'periodic', 'scatter', 
                            'periodic-violin', 'fairmode-target', 'fairmode-statsummary',
                            'taylor', 'boxplot']:

        for axis_i, axis in enumerate(fig.axes):
            if base_plot_type == 'taylor':
                lines = axis.parasites[0].lines
            else:
                lines = axis.lines
            for line_i, line in enumerate(lines):

                # extract data from each line
                data = []
                for x, y in line.get_xydata():
                    data.append({
                        "x": x,
                        "y": y,
                    })
                generated_output = pd.DataFrame(data)

                # save data, uncomment if we want to update it
                path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{plot_type}_{axis_i}_{line_i}.csv'
                if GENERATE_OUTPUT:
                    generated_output.to_csv(path, index=False)

                # read expected output
                expected_output = pd.read_csv(path)

                assert assert_frame_equal(
                generated_output, expected_output) is None

    elif base_plot_type in ['map']:

        # get coordinates and values
        for child in fig.axes[0].get_children():
            if isinstance(child, matplotlib.collections.PathCollection):
                coordinates = child.get_offsets()
                values = child.get_array()
                break

        # extract data from the table
        data = []
        for (lon, lat), val in zip(coordinates, values):
            data.append({
                "lon": lon,
                "lat": lat,
                "value": val
            })
        generated_output = pd.DataFrame(data)

        # save data, uncomment if we want to update it
        path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{plot_type}_values.csv'
        if GENERATE_OUTPUT:
            generated_output.to_csv(path, index=False)

        # read expected output
        expected_output = pd.read_csv(path, keep_default_na=False)

        assert assert_frame_equal(generated_output, expected_output) is None

    # check annotations
    if 'annotate' in plot_options:
        annotations = [child for child in fig.axes[0].get_children()
                       if type(child) == matplotlib.offsetbox.AnchoredOffsetbox][0]

        # extract annotations
        data = []
        for annotation in annotations.get_child().get_children():
            data.append({
                "dataset": annotation.get_text().split('|')[0].strip(),
                "annotation": annotation.get_text().split('|')[1].strip()
            })
        generated_output = pd.DataFrame(data)

        # save data, uncomment if we want to update it
        path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{plot_type}_annotations.csv'
        if GENERATE_OUTPUT:
            generated_output.to_csv(path, index=False)

        # read expected output
        expected_output = pd.read_csv(path, keep_default_na=False)

        assert assert_frame_equal(generated_output, expected_output) is None

def check_filter_data(inst, statistic_mode, network_type, filter):

    # Check filtered data
    filter_path = f'tests/reference/{network_type}/{statistic_mode}/data/data_{filter}.npy'
    read_data(inst, filter_path)

    # Reset filter and check original data
    inst.reset(initialise=True)
    orig_path = f'tests/reference/{network_type}/{statistic_mode}/data/data.npy'
    read_data(inst, orig_path)

    # Check filtered data is different from original data
    orig_output = np.load(orig_path, allow_pickle=True)
    filter_output = np.load(filter_path, allow_pickle=True)
    try:
        assert (not np.allclose(orig_output, filter_output, equal_nan=True))
    except ValueError as e:
        assert True

def save_data(inst, format, fname, network_type, statistic_mode):

    expected_path = f'tests/reference/{network_type}/{statistic_mode}/data/{fname}'
    if GENERATE_OUTPUT:
        inst.save(format=format, fname=expected_path)

    generated_path = f'saved_data/{fname}'
    inst.save(format=format, fname=generated_path)
    
    if format == 'npz':
        expected_data = np.load(f'{expected_path}.npz', allow_pickle=True)
        generated_data = np.load(f'{generated_path}.npz', allow_pickle=True)
        assert set(expected_data.files) == set(generated_data.files)
        for key in expected_data.files:
            a = expected_data[key]
            b = generated_data[key]
            # Test only numeric arrays
            if np.issubdtype(a.dtype, np.number) and np.issubdtype(b.dtype, np.number):
                assert np.allclose(a, b, equal_nan=True)
    
    elif format == 'nc':
        expected_data = xr.open_dataset(f'{expected_path}.nc')
        generated_data = xr.open_dataset(f'{generated_path}.nc')
        assert expected_data.equals(generated_data)

    elif format == 'conf':
        with open(f'{expected_path}.conf') as f:
            expected_conf = f.read()
        with open(f'{generated_path}.conf') as f:
            generated_conf = f.read()
        assert expected_conf == generated_conf

def check_unit_conversion(conv_obj, new_val, orig_val, orig_unit, new_unit, atol=1e-08):
    if np.isclose(conv_obj.converted_value, new_val, atol=atol) == False:
        print('{} {} to {} should be {}, is {}'.format(orig_val, orig_unit, new_unit, new_val, 
                                                       conv_obj.converted_value))
    else:
        print('{} {} to {} {} correctly converted'.format(orig_val, orig_unit, new_val, new_unit))
    assert np.isclose(conv_obj.converted_value, new_val, atol=atol)