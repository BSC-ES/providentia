#!/bin/sh
module purge

# CTE POWER
if [[ $BSC_MACHINE == "power" ]]; then
    module use /gpfs/projects/bsc32/software/rhel/7.4/ppc64le/POWER9/modules/all/
    module loaded gcc/6.4.0
    module loaded openmpi/3.0.0
    module load Python/3.7.0-foss-2018b
    module load Cartopy/0.17.0-foss-2018b-Python-3.7.0
    module load cftime/1.0.3.4-foss-2018b-Python-3.7.0
    module load netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0
    module load numpy/1.16.4-foss-2018b-Python-3.7.0
    module load pyproj/2.2.1-foss-2018b-Python-3.7.0
    module load scipy/1.3.0-foss-2018b-Python-3.7.0
    module load Shapely/1.6.4.post2-foss-2018b-Python-3.7.0
    module load xarray/0.15.1-foss-2018b-Python-3.7.0
    module load greasy/latest
    module load NCO/4.7.2-foss-2018b
# MareNostrum4
elif [[ $BSC_MACHINE == "mn4" ]]; then
    module load bsc/1.0 gcc/5.4.0 intel/2017.4 impi/2017.4 mkl/2017.4 python/3.7.4_ES udunits/2.2.25 hdf5/1.10.5_serial netcdf/4.4.1.1_serial geos/3.6.1 proj/4.9.3 gdal/2.2.3_py3 gsl/2.4 greasy/latest nco/4.6.7
# Nord3
elif [[ $BSC_MACHINE == "nord3" ]]; then
    module use /gpfs/projects/bsc32/software/suselinux/11/modules/all
    module load GCC/8.3.0 \
    Python/3.7.4-GCCcore-8.3.0 \
    geopandas/0.7.0-foss-2019b-Python-3.7.4 \
    Cartopy/0.18.0-foss-2019b-Python-3.7.4 \
    xarray/0.15.1-foss-2019b-Python-3.7.4 \
    netcdf4-python/1.5.3-foss-2019b-Python-3.7.4 \
    UDUNITS/2.2.26-GCCcore-8.3.0 \
    HDF5/1.10.5-gompi-2019b \
    netCDF/4.7.1-gompi-2019b \
    GEOS/3.7.2-foss-2019b-Python-3.7.4 \
    PROJ/6.2.1-GCCcore-8.3.0 \
    GSL/2.6-GCC-8.3.0 \
    NCO/4.9.2-foss-2019b
    module use /apps/modules/modulefiles/tools
    module load GREASY/latest
# Nord3v2
elif [ "${BSC_MACHINE}" == "nord3v2" ] || [ "${BSC_MACHINE}" = "amd" ]; then 
    module use /gpfs/projects/bsc32/software/suselinux/11/modules/all
    module load Python/3.7.4-GCCcore-8.3.0 \
    geopandas/0.7.0-foss-2019b-Python-3.7.4 \
    Cartopy/0.18.0-foss-2019b-Python-3.7.4 \
    xarray/0.15.1-foss-2019b-Python-3.7.4 \
    netcdf4-python/1.5.3-foss-2019b-Python-3.7.4 \
    UDUNITS/2.2.26-GCCcore-8.3.0 \
    HDF5/1.10.5-gompi-2019b \
    netCDF/4.7.1-gompi-2019b \
    GEOS/3.7.2-foss-2019b-Python-3.7.4 \
    PROJ/6.2.1-GCCcore-8.3.0 \
    GSL/2.6-GCC-8.3.0 \
    NCO/4.9.2-foss-2019b \
    greasy/2.2.3 \
    OpenMPI/4.0.5-GCC-8.3.0-nord3-v2
fi

export HDF5_USE_FILE_LOCKING=FALSE
