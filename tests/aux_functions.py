import hashlib
import matplotlib
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from pypdf import PdfReader

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


def plot(inst, statistic_mode, network_type, plot_type, plot_options=[], load=True):

    # load data
    if load:
        inst.load()

    # get zstat information from plot_type
    zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(
        plot_type=plot_type)

    # get base plot type (without stat and options)
    if zstat:
        base_plot_type = plot_type.split('-')[0]
    else:
        base_plot_type = plot_type.split('_')[0]

    path = f'tests/reference/{network_type}/{statistic_mode}/{base_plot_type}'

    # make_plot
    fig = inst.plot(plot_type, plot_options=plot_options, return_plot=True, save_data=True, 
                    save_data_path=path, tests_generate_output=GENERATE_OUTPUT)

    # check that a figure has been returned
    assert (type(fig) == matplotlib.figure.Figure)


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


def extract_pdf_info(pdf_path):
    """Extract text and number of pages from a PDF file, and return a hash of the text for comparison.

    Parameters
    ----------
    pdf_path : str
        Path to file

    Returns
    -------
    dict
        Dictionary with information
    """

    reader = PdfReader(pdf_path)
    text = ""

    for page in reader.pages:
        extracted = page.extract_text()
        if extracted:
            text += extracted

    return {
        "pages": len(reader.pages),
        "text_hash": hashlib.md5(text.encode("utf-8")).hexdigest(),
        "text": text,
    }

def check_report(network_type):
    """Check if report is the same as reference

    Parameters
    ----------
    network_type : str
        Network type
    """

    # Read PDFs
    generated_report_info = extract_pdf_info(f'reports/tests_{network_type}.pdf')
    expected_report_info = extract_pdf_info(f'tests/reference/{network_type}/tests_{network_type}.pdf')

    if generated_report_info["text_hash"] != expected_report_info["text_hash"]:
       print("Different content detected") 

    if generated_report_info["pages"] != expected_report_info["pages"]:
        print("Different number of pages detected")

    print(f"Generated PDF pages: {generated_report_info['pages']}")
    print(f"Expected PDF pages: {expected_report_info['pages']}")

    # Compare
    assert generated_report_info["pages"] == expected_report_info["pages"]
    assert generated_report_info["text_hash"] == expected_report_info["text_hash"]