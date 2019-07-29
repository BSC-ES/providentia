module purge
module load Python/3.7.0-foss-2018b
module load Cartopy/0.17.0-foss-2018b-Python-3.7.0
module load cftime/1.0.3.4-foss-2018b-Python-3.7.0
module load matplotlib/3.0.3-foss-2018b-Python-3.7.0
module load netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0
module load numpy/1.16.4-foss-2018b-Python-3.7.0
module load pandas/0.24.2-foss-2018b-Python-3.7.0
module load PyQt5/5.13.0-foss-2018b-Python-3.7.0
module load scipy/1.3.0-foss-2018b-Python-3.7.0
#module load Qt5/5.13.0-foss-2018b
salloc -t 03:00:00 -n 1 -c 8 --mem 2000 -J GI --x11
python GHOST_interactive.py
