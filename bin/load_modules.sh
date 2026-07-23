#!/bin/sh
module purge

# Marenostrum5
if [ "${BSC_MACHINE}" == "mn5" ]; then
    module load intel/2024.1 impi/2021.12
    module load greasy/2.2.4.2
    module load hdf5/1.14.4.2 netcdf/c-4.9.2_hdf5-1.14.4.2 pnetcdf/1.12.3 libexpat udunits 

# Nord4
elif [ "${BSC_MACHINE}" = "nord4" ]; then
    module load intel/2021.4 impi/2017.4
    module load greasy/2.2.4
    module load mkl/2017.4 netcdf/4.4.1.1 udunits/2.2.25 gsl/2.7
    module load singularity
        
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
