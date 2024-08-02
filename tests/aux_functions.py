import matplotlib
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from providentia.statistics import get_z_statistic_info


def read_data(inst, statistic_mode, network_type):

    # get data in memory in xarray format
    data = inst.get_data(format='xr')
    if network_type == 'ghost':
        networkspeci = 'EBAS|sconco3_data'
    else:
        networkspeci = 'nasa-aeronet_oneill_v3-lev15|od500aerocoarse_data'
    generated_output = data[networkspeci].values

    # save data, uncomment if we want to update it
    path = f'tests/reference/{network_type}/{statistic_mode}/data/data.npy'
    # np.save(path, generated_output)

    # read expected output
    expected_output = np.load(path, allow_pickle=True)

    assert (np.allclose(generated_output, expected_output, equal_nan=True))


def make_plot(inst, statistic_mode, network_type, plot_type, plot_options=[], expected_annotations=[], stats=[]):

    # make plot
    fig = inst.make_plot(plot_type, plot_options=plot_options,
                         return_plot=True, stats=stats)

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

    if base_plot_type in ['statsummary']:

        # get table
        for artist in fig.axes[0].get_children():
            if isinstance(artist, matplotlib.table.Table):
                table = artist
                break

        # extract data from the table
        data = []
        for (row, col), cell in table.get_celld().items():
            data.append({
                "row": row,
                "col": col,
                "value": cell.get_text().get_text()
            })
        generated_output = pd.DataFrame(data)

        # save data, uncomment if we want to update it
        if 'bias' in plot_options:
            path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{base_plot_type}_bias_table_values.csv'
        else:
            path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{base_plot_type}_table_values.csv'
        # generated_output.to_csv(path, index=False)

        expected_output = pd.read_csv(path, keep_default_na=False)

        assert assert_frame_equal(generated_output, expected_output) is None

    elif base_plot_type in ['timeseries', 'distribution']:

        if 'annotate' in plot_options and expected_annotations:
            annotations = [child for child in fig.axes[0].get_children()
                           if type(child) == matplotlib.offsetbox.AnchoredOffsetbox][0]

            for annotation, expected_annotation in zip(annotations.get_child().get_children(),
                                                       expected_annotations[base_plot_type]):
                print(annotation.get_text())
                assert annotation.get_text() == expected_annotation

        # iterate through plotted lines
        for line_i, line in enumerate(fig.axes[0].lines):

            # extract data from each line
            data = []
            for x, y in line.get_xydata():
                data.append({
                    "x": x,
                    "y": y,
                })
            generated_output = pd.DataFrame(data)

            # save data, uncomment if we want to update it
            path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{base_plot_type}_line_{line_i}.csv'
            # generated_output.to_csv(path, index=False)

            # read expected output
            expected_output = pd.read_csv(path)

            assert assert_frame_equal(generated_output, expected_output) is None

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
        path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}/{base_plot_type}_values.csv'
        # generated_output.to_csv(path, index=False)

        # read expected output
        expected_output = pd.read_csv(path, keep_default_na=False)
        
        assert assert_frame_equal(generated_output, expected_output) is None
