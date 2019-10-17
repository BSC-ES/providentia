#!/bin/sh
module purge

#CTE POWER
if [ $BSC_MACHINE == "power" ]; then
    module load Python/3.7.0-foss-2018b
    module load Cartopy/0.17.0-foss-2018b-Python-3.7.0
    module load cftime/1.0.3.4-foss-2018b-Python-3.7.0
    module load netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0
    module load numpy/1.16.4-foss-2018b-Python-3.7.0
    module load pyproj/2.2.1-foss-2018b-Python-3.7.0
    module load scipy/1.3.0-foss-2018b-Python-3.7.0
    module load Shapely/1.6.4.post2-foss-2018b-Python-3.7.0
#Marenustrum4
elif [ $BSC_MACHINE == "mn4" ]; then
    module load bsc/1.0 gcc/5.4.0 intel/2017.4 impi/2017.4 mkl/2017.4 python/3.7.4_ES udunits/2.2.25 hdf5/1.10.5 netcdf/4.4.1.1_lf geos/3.6.1 proj/4.9.3 gdal/2.2.3
fi
