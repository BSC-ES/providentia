#!/bin/sh
module purge

if [ "${BSC_MACHINE}" == "power" ]; then
    module load Python/3.7.0-foss-2018b \
     Cartopy/0.17.0-foss-2018b-Python-3.7.0 \
     cftime/1.0.3.4-foss-2018b-Python-3.7.0 \
     matplotlib/3.0.3-foss-2018b-Python-3.7.0 \
     netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0 \
     numpy/1.16.4-foss-2018b-Python-3.7.0 \
     pandas/0.24.2-foss-2018b-Python-3.7.0 \
     PyQt5/5.13.0-foss-2018b-Python-3.7.0 \
     scipy/1.3.0-foss-2018b-Python-3.7.0
#Marenustrum4
elif [ "${BSC_MACHINE}" == "mn4" ]; then
    module load bsc/1.0 gcc/5.4.0 intel/2017.4 impi/2017.4 mkl/2017.4 python/3.7.4_ES udunits/2.2.25 hdf5/1.10.5 netcdf/4.4.1.1_lf geos/3.6.1 proj/4.9.3 gdal/2.2.3
#Workstations/fatnodes
else
    module load Python/3.7.3-foss-2015a \
     ConfigArgParse/0.14.0-foss-2015a-Python-3.7.3 \
     Cartopy/0.17.0-foss-2015a-Python-3.7.3 \
     matplotlib/3.1.1-foss-2015a-Python-3.7.3 \
     netcdf4-python/1.5.1.2-foss-2015a-Python-3.7.3 \
     numpy/1.16.4-foss-2015a-Python-3.7.3 \
     pandas/0.24.2-foss-2015a-Python-3.7.3 \
     PyQt5/5.13.0-foss-2015a-Python-3.7.3 \
     seaborn/0.9.0-foss-2015a-Python-3.7.3 \
     scipy/1.3.0-foss-2015a-Python-3.7.3
fi
