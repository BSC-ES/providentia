#!/bin/bash

#SBATCH --job-name=SubmitG                                         
#SBATCH --ntasks=1                                                
#SBATCH −−cpus−per−task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=01:00:00                                            

#module load Python/3.7.0-foss-2018b
#module load Cartopy/0.17.0-foss-2018b-Python-3.7.0
#module load netcdf4-python/1.5.1.2-foss-2018b-Python-3.7.0
#module load numpy/1.16.4-foss-2018b-Python-3.7.0
#module load pyproj/2.2.1-foss-2018b-Python-3.7.0
#module load scipy/1.3.0-foss-2018b-Python-3.7.0
#module load Shapely/1.6.4.post2-foss-2018b-Python-3.7.0
srun --cpu-bind=core python -u test.py 
