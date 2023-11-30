#!/bin/sh
module purge

# CTE Power9
if [ "${BSC_MACHINE}" == "power" ]; then
    module load Python/3.7.0-foss-2018b \
        Cartopy/0.20.3-foss-2018b-Python-3.7.0 \
        cftime/1.0.3.4-foss-2018b-Python-3.7.0 \
        matplotlib/3.1.1-foss-2018b-Python-3.7.0 \
        netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0 \
        numpy/1.16.4-foss-2018b-Python-3.7.0 \
        pandas/0.24.2-foss-2018b-Python-3.7.0 \
        PyQt5/5.13.0-foss-2018b-Python-3.7.0 \
        scipy/1.3.0-foss-2018b-Python-3.7.0 \
        kdepy/1.1.1-fosscuda-2019b-Python-3.7.4 \
        Ghostscript/9.50-GCCcore-8.3.0

# Marenostrum4
elif [ "${BSC_MACHINE}" == "mn4" ]; then
    module load intel impi mkl \
        gcc/5.4.0 python/3.7.4_ES_test \
        sqlite/3.41.0 proj/8.2.1 hdf5 netcdf \
        geos/3.10.2 gdal/3.6.2 libtiff

# Nord3v2
elif [ "${BSC_MACHINE}" == "nord3v2" ]; then 
    module use /gpfs/projects/bsc32/software/suselinux/11/modules/all
    module load Python/3.7.4-GCCcore-8.3.0 \
        Ghostscript/9.50-GCCcore-8.3.0 \
        xarray/0.19.0-foss-2019b-Python-3.7.4 \
        jupyterlab/3.0.9-foss-2019b-Python-3.7.4 \
        matplotlib/3.1.1-foss-2019b-Python-3.7.4 \
        seaborn/0.9.0-foss-2019b-Python-3.7.4 \
        Cartopy/0.20.3-foss-2019b-Python-3.7.4 \
        netcdf4-python/1.5.3-foss-2019b-Python-3.7.4 \
        geopandas/0.7.0-foss-2019b-Python-3.7.4 \
        PyQt5/5.13.1-GCCcore-8.3.0-Python-3.7.4 \
        Qt5/5.14.1-GCCcore-8.3.0 \
        OpenMPI/4.0.5-GCC-8.3.0-nord3-v2 \
        pyproj/3.2.1-foss-2019b-Python-3.7.4 \
        kdepy/1.1.1-foss-2019b-Python-3.7.4 

# Hub
elif [ "${ip}" == "84.88.185.48" ]; then
    module load Python/3.9.6-GCCcore-11.2.0 \
        matplotlib/3.7.2-foss-2021b-Python-3.9.6 \
        netcdf4-python/1.6.1-foss-2021b-Python-3.9.6 \
        numpy/1.23.3-foss-2021b-Python-3.9.6 \
        pandas/1.5.3-foss-2021b-Python-3.9.6 \
        Qt5/5.15.5-GCCcore-11.2.0 \
        scipy/1.10.1-foss-2021b-Python-3.9.6 \
        pyproj/3.4.0-foss-2021b-Python-3.9.6 \
        ConfigArgParse/1.7-foss-2021b-Python-3.9.6 \
        PyQt5/5.15.5-GCCcore-11.2.0 \
        Cartopy/0.22.0-foss-2021b-Python-3.9.6 \
        Seaborn/0.12.2-foss-2021b-Python-3.9.6 \
        kdepy/1.1.8-foss-2021b-Python-3.9.6 \
        Ghostscript/10.01.2-GCCcore-11.2.0 \
        cftime/1.0.3.4-foss-2021b-Python-3.9.6
        
# Workstations
# TODO: Review modules
else
    module load Python/3.7.3-foss-2015a \
        ConfigArgParse/0.14.0-foss-2015a-Python-3.7.3 \
        Cartopy/0.20.3-foss-2015a-Python-3.7.3 \
        matplotlib/3.1.1-foss-2015a-Python-3.7.3 \
        netcdf4-python/1.5.1.2-foss-2015a-Python-3.7.3 \
        numpy/1.16.4-foss-2015a-Python-3.7.3 \
        pandas/0.24.2-foss-2015a-Python-3.7.3 \
        PyQt5/5.13.0-foss-2015a-Python-3.7.3 \
        Qt/5.13.0-foss-2015a \
        seaborn/0.9.0-foss-2015a-Python-3.7.3 \
        scipy/1.3.0-foss-2015a-Python-3.7.3 \
        pyproj/3.2.1-foss-2015a-Python-3.7.3
fi

export HDF5_USE_FILE_LOCKING=FALSE
