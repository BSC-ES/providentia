#WRITTEN BY DENE BOWDALO

###------------------------------------------------------------------------------------###
###IMPORT MODULES
###------------------------------------------------------------------------------------###

import glob
import numpy as np
import os
import subprocess
import sys
import time

###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

start = time.time()

#define current working directory and arguments/submit/interpolation log subdirectories
working_directory = os.getcwd()
arguments_dir = '%s/arguments'%(working_directory)
submit_dir = '%s/submit'%(working_directory)
interpolation_log_dir = '%s/interpolation_logs'%(working_directory)

#set SLURM job ID as unique ID for tracking tasks defined to process in the configuration file
unique_ID = sys.argv[1]

#Read configuration file
from configuration import start_date, end_date, experiments_to_process, species_to_process, grid_types_to_process, model_temporal_resolutions_to_process, GHOST_networks_to_interpolate_against, temporal_resolutions_to_output

#read defined experiments dictionary
from defined_experiments import defined_experiments_dictionary

#create argument text file
arguments_file= open("%s/%s.txt"%(arguments_dir, unique_ID),"w")

#create list to hold all .out files generated for each task processed
task_out_names = []

#iterate through desired experiment IDs
for experiment_to_process in experiments_to_process:

    #get experiment specific directory 
    exp_dir = defined_experiments_dictionary[experiment_to_process]

    #get all grid type subdirectories for current experiment
    available_grid_types = [name for name in os.listdir("%s/"%(exp_dir)) if os.path.isdir("%s/%s"%(exp_dir,name))]
    #get intersection between desired grid types to process and grid types available in directory
    available_grid_types = [x for x in available_grid_types if x in grid_types_to_process]

    #iterate through grid types to process
    for grid_type_to_process in available_grid_types:

        #get all temporal resolution subdirectories for current experiment/grid_type
        available_model_temporal_resolutions = [name for name in os.listdir("%s/%s/"%(exp_dir,grid_type_to_process)) if os.path.isdir("%s/%s/%s"%(exp_dir,grid_type_to_process,name))]
        #get intersection between desired temporal resolutions to process and temporal resolutions available in directory
        available_model_temporal_resolutions = [x for x in available_model_temporal_resolutions if x in model_temporal_resolutions_to_process]
        
        #iterate though temporal resolutions to process (should only be hourly)
        for model_temporal_resolution_to_process in available_model_temporal_resolutions:

            #get all species subdirectories for current experiment/grid_type/temporal_resolution
            available_species = [name for name in os.listdir("%s/%s/%s"%(exp_dir,grid_type_to_process,model_temporal_resolution_to_process)) if os.path.isdir("%s/%s/%s/%s"%(exp_dir,grid_type_to_process,model_temporal_resolution_to_process,name))]
            #get intersection between desired species to process and species available in directory
            available_species = [x for x in available_species if x in species_to_process]

            #iterate though species to process
            for speci_to_process in species_to_process:

                #iterate through GHOST observational networks to interpolate against
                for GHOST_network_to_interpolate_against in GHOST_networks_to_interpolate_against:
                    
                    #iterate through temporal_resolutions to output
                    for temporal_resolution_to_output in temporal_resolutions_to_output:

                        #get all GHOST network/species/temporal resolution observational files
                        obs_files = np.sort(glob.glob('/esarchive/obs/ghost/%s/%s/%s/%s*.nc'%(GHOST_network_to_interpolate_against, temporal_resolution_to_output, speci_to_process, speci_to_process)))

                        #get all relevant experiment files
                        exp_files = np.sort(glob.glob('%s/%s/%s/%s/%s*.nc'%(exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, speci_to_process)))
                        #remove new type of REDUCE output with member number (will be handled in future)
                        exp_files = np.array([exp_file for exp_file in exp_files if '-000_' not in exp_file]) 

                        #get all observational file start dates (year and month)
                        obs_files_dates = []
                        for obs_file in obs_files:
                            obs_file_date = obs_file.split('%s_'%(speci_to_process))[-1].split('.nc')[0]
                            obs_files_dates.append(obs_file_date[:6])
                        obs_files_dates = np.sort(obs_files_dates)

                        #get all experiment file start dates (year and month)
                        exp_files_dates = []   
                        for exp_file in exp_files:
                            exp_file_date = exp_file.split('%s_'%(speci_to_process))[-1].split('.nc')[0]
                            exp_files_dates.append(exp_file_date[:6])
                        exp_files_dates = np.sort(exp_files_dates) 

                        #remove observational/experiment files outside date ranges 
                        obs_files_ii = np.array([obs_files_ii for obs_files_ii, obs_files_date in enumerate(obs_files_dates) if (int(obs_files_date) >= int(start_date)) and (int(obs_files_date) < int(end_date))])       
                        if len(obs_files_ii) == 0:
                            continue
                        obs_files = obs_files[obs_files_ii]
                        obs_files_dates = obs_files_dates[obs_files_ii]
                        exp_files_ii = np.array([exp_files_ii for exp_files_ii, exp_files_date in enumerate(exp_files_dates) if (int(exp_files_date) >= int(start_date)) and (int(exp_files_date) < int(end_date))])      
                        if len(exp_files_ii) == 0:
                            continue
                        exp_files = exp_files[exp_files_ii]
                        exp_files_dates = exp_files_dates[exp_files_ii]

                        #get intersection of file yearmonths between observations and experiment
                        intersect_yearmonths = np.intersect1d(obs_files_dates, exp_files_dates)

                        #if have no intersecting months, continue
                        if len(intersect_yearmonths) == 0:
                            continue
    
                        #create directories to store slurm output/error logs for interpolation task of specific combination of iterated variables (if does not already exist)
                        if not os.path.exists('%s/%s/%s/%s/%s/%s/%s'%(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output)):
                            os.makedirs('%s/%s/%s/%s/%s/%s/%s'%(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output))

                        #iterate through intersecting yearmonths and write all current variable arguments to arguments file
                        for yearmonth in intersect_yearmonths:

                            #append current iterative arguments to arguments file
                            arguments_file.write("%s %s %s %s %s %s %s\n"%(experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth))                         

                            #append location of .out file that will be output for each processed task
                            task_out_names.append('%s/%s/%s/%s/%s/%s/%s/%s.out'%(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth))                           

                            #remove previous output logs
                            if os.path.exists('%s/%s/%s/%s/%s/%s/%s/%s.out'%(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth)):
                                os.remove('%s/%s/%s/%s/%s/%s/%s/%s.out'%(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth))

#close arguments file
arguments_file.close()

#get total number of tasks to submit
N_tasks = len(task_out_names)

#create job submit batch script
submit_file = open("%s/%s.sh"%(submit_dir, unique_ID),"w")
submit_file.write("#!/bin/bash\n")
submit_file.write("\n")
submit_file.write("#SBATCH --job-name=G%s\n"%(unique_ID))
submit_file.write("#SBATCH --ntasks=1\n")
submit_file.write("#SBATCH --cpus-per-task=1\n")
submit_file.write("#SBATCH --mem-per-cpu=2G\n") 
submit_file.write("#SBATCH --time=01:00:00\n")
submit_file.write("#SBATCH --array=1-%s%%1000\n"%(N_tasks))
submit_file.write("\n")
submit_file.write("arguments_store=%s/%s.txt\n"%(arguments_dir, unique_ID))
submit_file.write("argument_a=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n")
submit_file.write("argument_b=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')\n")  
submit_file.write("argument_c=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')\n")
submit_file.write("argument_d=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $4}')\n")
submit_file.write("argument_e=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $5}')\n")
submit_file.write("argument_f=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $6}')\n")
submit_file.write("argument_g=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $7}')\n")
submit_file.write("\n") 
submit_file.write("srun --cpu-bind=core --output=%s/$argument_a/$argument_b/$argument_c/$argument_d/$argument_e/$argument_f/$argument_g.out --error=%s/$argument_a/$argument_b/$argument_c/$argument_d/$argument_e/$argument_f/$argument_g.err python -u %s/GHOST_experiment_interpolation.py $argument_a $argument_b $argument_c $argument_d $argument_e $argument_f $argument_g"%(interpolation_log_dir, interpolation_log_dir, working_directory))

#close submit file
submit_file.close()

#submit job (make sure it is submitted successfully)
cmd = ['sbatch', '%s.sh'%(unique_ID)]

submit_complete = False
while submit_complete == False: 
    submit_process = subprocess.Popen(cmd, cwd=submit_dir, stdout=subprocess.PIPE)
    submit_status = submit_process.communicate()[0]
    submit_return_code = submit_process.returncode

    #if sbatch fails, time out for 20 seconds and then try again
    if submit_return_code != 0:
        time.sleep(20)
        continue
    else:
        submit_complete = True
        print('%s INTERPOLATION JOB/S SUBMITTED'%(N_tasks))

#now all interpoaltion tasks have been submitted
#monitor number of jobs in queue (every 30 seconds) until there are 0 left in the squeue
all_tasks_finished = False

while all_tasks_finished == False:
    cmd = ['squeue', '-h', '-n', 'G%s'%(unique_ID)]
    squeue_process =  subprocess.Popen(cmd, stdout=subprocess.PIPE, encoding='utf8')
    squeue_status = squeue_process.communicate()[0]
    n_jobs_in_queue = len(squeue_status.split('\n')[:-1])
    #if number of jobs in queue > 0, then sleep for 10 seconds and then check again how many jobs there are in queue
    if n_jobs_in_queue > 0:
        time.sleep(10)
        continue  

    #if no more jobs in the squeue, then now check the outcome of all the jobs  
    #if any jobs have failed, write them out to file
    failed_tasks=[]
    for task_out in task_out_names:
        line = subprocess.check_output(['tail', '-1', task_out], encoding='utf8').strip()
        if line != 'DONE':
            failed_tasks.append(task_out)

    #break out of while loop     
    all_tasks_finished = True

end = time.time()

#remove default .out files generated in submit directory (don't know why these are created --> probably to do with calling srun...)
for f in glob.glob("%s/*.out"%(submit_dir)):
    os.remove(f)

#remove default .out file generated in working directory (don't know why this is created --> probably to do with calling srun...)
for f in glob.glob("%s/*.out"%(working_directory)):
    os.remove(f)

if len(failed_tasks) == 0:
    print('ALL INTERPOLATIONS COMPLETED SUCCESSFULLY IN %s MINUTES'%((end-start)/60.))
else:
    print('%s INTERPOLATIONS FINISHED SUCCESSFULLY IN %s MINUTES'%(N_tasks-len(failed_tasks),(end-start)/60.))
    print('THE FOLLOWING INTERPOLATIONS FAILED:')
    print(failed_tasks)
