import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def read_data(inst, test_name):

    # get data in memory in xarray format
    data = inst.get_data(format='xr')
    generated_output = data['EBAS|sconco3_data'].values

    # save data, uncomment if we want to update it
    # np.save(f'tests/reference/{test_name}/data/EBAS_sconco3_data', generated_output)

    # read expected output
    expected_output = np.load(f'tests/reference/{test_name}/data/EBAS_sconco3_data.npy', allow_pickle=True)
    
    assert (np.allclose(generated_output, expected_output, equal_nan=True))


def make_plot(inst, test_name, plot_type, plot_options=[], expected_annotations=[]):

    # make plot
    fig = inst.make_plot(plot_type, plot_options=plot_options, return_plot=True)

    # check that a figure has been returned
    assert (type(fig) == matplotlib.figure.Figure)

    if plot_type in ['statsummary']:

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
        # generated_output.to_csv(f'tests/reference/{test_name}/{plot_type}/{plot_type}_table_values.csv', index=False)

        expected_output = pd.read_csv(
            f'tests/reference/{test_name}/{plot_type}/{plot_type}_table_values.csv', keep_default_na=False)

        assert generated_output.equals(expected_output)

    elif plot_type in ['timeseries', 'distribution']:
        
        # check if annotations are correct
        if 'annotate' in plot_options and expected_annotations:
            annotations = [child for child in fig.axes[0].get_children() 
                        if type(child) == matplotlib.offsetbox.AnchoredOffsetbox][0]
            
            for annotation, expected_annotation in zip(annotations.get_child().get_children(), 
                                                    expected_annotations[plot_type]):
                assert annotation.get_text() == expected_annotation
        
        # iterate through plotted lines
        for line_i, line in enumerate(fig.axes[0].lines):
            
            # get data from line
            generated_output = line.get_xydata()

            # save data, uncomment if we want to update it
            # np.save(f'tests/reference/{test_name}/{plot_type}/{plot_type}_line_{line_i}', generated_output)
                
            # read expected output
            expected_output = np.load(f'tests/reference/{test_name}/{plot_type}/{plot_type}_line_{line_i}.npy', allow_pickle=True)

            # check data for each line is correct
            assert (np.allclose(generated_output, expected_output, equal_nan=True))
