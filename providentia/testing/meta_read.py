import numpy as np
from netCDF4 import Dataset, chartostring

fnames = ['/esarchive/obs/nasa-aeronet/directsun_v3-lev15/3hourly/ae440-870aero/ae440-870aero_202303.nc',
          '/esarchive/obs/nasa-aeronet/directsun_v2-lev20/3hourly/ae440-870aero/ae440-870aero_201809.nc',
          '/esarchive/obs/nasa-aeronet/directsun_v3-lev20/3hourly_old/ae440-870aero/ae440-870aero_202012.nc',
          '/esarchive/obs/eea/eionet/hourly/sconco3/sconco3_202303.nc',
          '/esarchive/obs/cnrs-lisa/indaaf/hourly/pm10/pm10_201412.nc',
          '/esarchive/obs/port_barcelona/port-barcelona/hourly/sconcno2/sconcno2_202112.nc',
          '/esarchive/obs/csic/csic/monthly/sconcnh3/sconcnh3_201908.nc',
          '/esarchive/obs/MITECO_VOC/MITECO_VOC/hourly/sconcc6h6/sconcc6h6_202201.nc']

meta_vars = ['station_reference', 'station_name', 'station_classification', 'area_classification', 'longitude', 
           'latitude', 'altitude']

station_indices = np.array([0,2,7,9,11])

for fname in fnames:

    print(fname)

    root = Dataset(fname)

    for meta_var in meta_vars:

        if meta_var == 'station_reference':
            if "station_reference" not in root.variables:
                meta_var = 'station_code'
            else:
                meta_var = 'station_reference'

        # check meta variable is in netCDF
        if meta_var not in root.variables:
            continue

        meta_shape = root[meta_var].shape

        if meta_shape[0] < len(station_indices):
            station_indices = np.array([0])

        meta_val = root[meta_var][station_indices]
        meta_val_dtype = np.array([meta_val[0]]).dtype

        if len(meta_shape) == 2:
            if meta_val_dtype == np.dtype(object):
                meta_val =  np.array([''.join(val) for val in meta_val])
            else:
                meta_val = chartostring(meta_val)

print()
print('ALL FILES CAN BE READ CORRECTLY BY PROVIDENTIA')