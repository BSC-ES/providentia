#!/bin/sh
module purge

machine=$(hostname -a)
#CTE POWER
if [[ $machine == *"power.cte"* ]]; then
    module load Python/3.7.0-foss-2018b
    module load Cartopy/0.17.0-foss-2018b-Python-3.7.0
    module load cftime/1.0.3.4-foss-2018b-Python-3.7.0
    module load matplotlib/3.0.3-foss-2018b-Python-3.7.0
    module load netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0
    module load numpy/1.16.4-foss-2018b-Python-3.7.0
    module load pandas/0.24.2-foss-2018b-Python-3.7.0
    module load PyQt5/5.13.0-foss-2018b-Python-3.7.0
    module load scipy/1.3.0-foss-2018b-Python-3.7.0
    #get allocation of resources on CTE-POWER machine
    salloc -t 02:00:00 -n 1 -c 16 --mem=50Gb -J PRV -q debug --x11=first
#Marenustrum4
elif [[ $machine == *"bsc.mn"* ]]; then
    module load gcc/5.4.0 intel/2018.4 impi/2018.4 mkl/2018.4 python/3.7.4_ES udunits/2.2.25  hdf5/1.10.5  netcdf/4.4.1.1_lf geos/3.6.1 proj/4.9.3 gdal/2.2.3
    #get allocation of resources on MN4 machine
    salloc -t 02:00:00 -n 1 -c 16 -J PRV -q debug --x11=first
#Workstations
else
    module load Python/3.7.3-foss-2015a
    module load Cartopy/0.17.0-foss-2015a-Python-3.7.3
    module load matplotlib/3.1.1-foss-2015a-Python-3.7.3
    module load netcdf4-python/1.5.1.2-foss-2015a-Python-3.7.3
    module load numpy/1.16.4-foss-2015a-Python-3.7.3
    module load pandas/0.24.2-foss-2015a-Python-3.7.3
    module load PyQt5/5.13.0-foss-2015a-Python-3.7.3
    module load seaborn/0.9.0-foss-2015a-Python-3.7.3
    module load scipy/1.3.0-foss-2015a-Python-3.7.3
fi
