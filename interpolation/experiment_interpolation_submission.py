#WRITTEN BY DENE BOWDALO

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#experiment_interpolation_submission.py

#module which sets up parallel execution of interpolation jobs for defined configuration

###--------------------------------------------------------------------------------------------------###
###IMPORT MODULES
###--------------------------------------------------------------------------------------------------###

import copy
import glob
from netCDF4 import Dataset
import numpy as np
import os
import subprocess
import sys
import time

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

#start timer
start = time.time()

#define current working directory and arguments/submit/interpolation log subdirectories
working_directory = os.getcwd()
arguments_dir = '{}/arguments'.format(working_directory)
submit_dir = '{}/submit'.format(working_directory)
interpolation_log_dir = '{}/interpolation_logs'.format(working_directory)

#set SLURM job ID as unique ID for tracking tasks defined to process in the configuration file
unique_ID = sys.argv[1]

#Read configuration file
from configuration import start_date, end_date, experiments_to_process, species_to_process, grid_types_to_process, model_temporal_resolutions_to_process, GHOST_networks_to_interpolate_against, temporal_resolutions_to_output, n_neighbours_to_find, qos

#read defined experiments dictionary
from defined_experiments import defined_experiments_dictionary

#create arguments lost
arguments_list = []

#create list to hold all .out files generated for each task processed
task_out_names = []

#iterate through desired experiment IDs
for experiment_to_process in experiments_to_process:

    #get experiment specific directory (take gpfs directory preferentially over esarchive directory) 
    exp_dict = defined_experiments_dictionary[experiment_to_process]
    if 'gpfs' in list(exp_dict.keys()):
        exp_dir = exp_dict['gpfs']
    else:
        exp_dir = exp_dict['esarchive']

    #get all grid type subdirectories for current experiment
    available_grid_types = [name for name in os.listdir("{}/".format(exp_dir)) if os.path.isdir("{}/{}".format(exp_dir,name))]
    #get intersection between desired grid types to process and grid types available in directory
    available_grid_types = [x for x in available_grid_types if x in grid_types_to_process]

    #iterate through grid types to process
    for grid_type_to_process in available_grid_types:

        #get all temporal resolution subdirectories for current experiment/grid_type
        available_model_temporal_resolutions = [name for name in os.listdir("{}/{}/".format(exp_dir,grid_type_to_process)) if os.path.isdir("{}/{}/{}".format(exp_dir,grid_type_to_process,name))]
        #get intersection between desired temporal resolutions to process and temporal resolutions available in directory
        available_model_temporal_resolutions = [x for x in available_model_temporal_resolutions if x in model_temporal_resolutions_to_process]
        
        #iterate though temporal resolutions to process (should only be hourly)
        for model_temporal_resolution_to_process in available_model_temporal_resolutions:

            #get all species subdirectories for current experiment/grid_type/temporal_resolution
            available_species = [name for name in os.listdir("{}/{}/{}".format(exp_dir,grid_type_to_process,model_temporal_resolution_to_process)) if os.path.isdir("{}/{}/{}/{}".format(exp_dir,grid_type_to_process,model_temporal_resolution_to_process,name))]
            #get intersection between desired species to process and species available in directory
            available_species = [x for x in available_species if x in species_to_process]

            #iterate though species to process
            for speci_to_process in species_to_process:

                #iterate through GHOST observational networks to interpolate against
                for GHOST_network_to_interpolate_against in GHOST_networks_to_interpolate_against:
                    
                    #iterate through temporal_resolutions to output
                    for temporal_resolution_to_output in temporal_resolutions_to_output:

                        #get all GHOST network/species/temporal resolution observational files
                        obs_files = np.sort(glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/{}/{}/{}/{}*.nc'.format(GHOST_network_to_interpolate_against, temporal_resolution_to_output, speci_to_process, speci_to_process)))

                        #get all relevant experiment files
                        exp_files = np.sort(glob.glob('{}/{}/{}/{}/{}*.nc'.format(exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, speci_to_process)))
                        #remove new type of REDUCE output with member number (will be handled in future)
                        exp_files = np.array([exp_file for exp_file in exp_files if '-000_' not in exp_file]) 

                        #get all observational file start dates (year and month)
                        obs_files_dates = []
                        for obs_file in obs_files:
                            obs_file_date = obs_file.split('{}_'.format(speci_to_process))[-1].split('.nc')[0]
                            obs_files_dates.append(obs_file_date[:6])
                        obs_files_dates = np.sort(obs_files_dates)

                        #get all experiment file start dates (year and month)
                        exp_files_dates = []   
                        for exp_file in exp_files:
                            exp_file_date = exp_file.split('{}_'.format(speci_to_process))[-1].split('.nc')[0]
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
                        if not os.path.exists('{}/{}/{}/{}/{}/{}/{}'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output)):
                            os.makedirs('{}/{}/{}/{}/{}/{}/{}'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output))

                        #iterate through intersecting yearmonths and write all current variable arguments to arguments file
                        for yearmonth in intersect_yearmonths:

                            #append current iterative arguments to arguments list               
                            arguments_list.append("{} {} {} {} {} {} {} {}".format(experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth, n_neighbours_to_find))                         

                            #append location of .out file that will be output for each processed task
                            task_out_names.append('{}/{}/{}/{}/{}/{}/{}/{}.out'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth))                           

                            #remove previous output logs
                            if os.path.exists('{}/{}/{}/{}/{}/{}/{}/{}.out'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth)):
                                os.remove('{}/{}/{}/{}/{}/{}/{}/{}.out'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, yearmonth))

#initialise variable to store names of submit files created
submit_fnames = []

#split arguments list so that have no more than 999 arguments in each
arguments_split = np.array_split(arguments_list,int(np.ceil(len(arguments_list)/999)))

#loop through split arguments, creating a submit file per group
for file_ii, grouped_arguments in enumerate(arguments_split):

    #create arguments file/s
    arguments_file= open("{}/{}_{}.txt".format(arguments_dir, unique_ID, file_ii),"w")        

    for argument in grouped_arguments:
        arguments_file.write("{}\n".format(argument))                         

    #close arguments file
    arguments_file.close()
                        
    #create job submit batch script
    submit_fname = '{}_{}.sh'.format(unique_ID, file_ii)
    submit_file = open("{}/{}".format(submit_dir, submit_fname),"w")
    submit_file.write("#!/bin/bash\n")
    submit_file.write("\n")
    submit_file.write("#SBATCH --job-name=PRVI_{}\n".format(unique_ID))
    submit_file.write("#SBATCH --ntasks=1\n")
    submit_file.write("#SBATCH --time=01:00:00\n")
    submit_file.write("#SBATCH --array=1-{}\n".format(len(grouped_arguments)))
    submit_file.write("#SBATCH --qos={}\n".format(qos))
    submit_file.write("#SBATCH --output=/dev/null\n")
    submit_file.write("#SBATCH --error=/dev/null\n")
    submit_file.write("\n")
    submit_file.write("arguments_store={}/{}_{}.txt\n".format(arguments_dir, unique_ID, file_ii))
    submit_file.write("argument_a=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n")
    submit_file.write("argument_b=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')\n")  
    submit_file.write("argument_c=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $3}')\n")
    submit_file.write("argument_d=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $4}')\n")
    submit_file.write("argument_e=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $5}')\n")
    submit_file.write("argument_f=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $6}')\n")
    submit_file.write("argument_g=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $7}')\n")
    submit_file.write("argument_h=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $8}')\n")
    submit_file.write("\n")
    submit_file.write("source {}/load_modules.sh\n".format(working_directory)) 
    submit_file.write("\n")
    submit_file.write("srun --cpu-bind=core --output={}/$argument_a/$argument_b/$argument_c/$argument_d/$argument_e/$argument_f/$argument_g.out python -u {}/experiment_interpolation.py $argument_a $argument_b $argument_c $argument_d $argument_e $argument_f $argument_g $argument_h".format(interpolation_log_dir, working_directory))

    #close submit file
    submit_file.close()

    #add submit filename to list of submitted files
    submit_fnames.append(submit_fname)
    
#time start of interpolation jobs
interpolation_start = time.time()

#submit interpolation jobs (make sure they submitted successfully)
for submit_fname in submit_fnames:
    cmd = ['sbatch', '-A', 'bsc32', '{}'.format(submit_fname)]

    submit_complete = False
    while submit_complete == False: 
        submit_process = subprocess.Popen(cmd, cwd=submit_dir, stdout=subprocess.PIPE)
        submit_status = submit_process.communicate()[0]
        submit_return_code = submit_process.returncode

        #if sbatch fails, time out for 60 seconds and then try again
        if submit_return_code != 0:
            time.sleep(60)
            continue
        else:
            submit_complete = True

        #take a 1 second pause between submittals (to help slurm scheduler)
        time.sleep(1)

#now all interpolation jobs have been submitted
#monitor number of jobs in queue (every 10 seconds) until there are 0 left in the squeue
all_tasks_finished = False
 
while all_tasks_finished == False:
    cmd = ['squeue', '-h', '-n', 'PRVI_{}'.format(unique_ID)]
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
    process_times=[]
    for task_out in task_out_names:
        line = subprocess.check_output(['tail', '-1', task_out], encoding='utf8').strip()
        if 'DONE' not in line:
            failed_tasks.append(task_out)
        else:
            process_times.append(float(line.split('DONE:')[-1]))

    #break out of while loop     
    all_tasks_finished = True

#get interpolation time
interpolation_time = time.time() - interpolation_start

#remove default .out file generated in working directory (don't know why this is created --> probably to do with calling srun...)
for f in glob.glob("{}/*.out".format(working_directory)):
    os.remove(f)

#stop timer
total_time = time.time()-start

if len(failed_tasks) == 0:
    #get queue time 
    process_time = np.max(process_times)
    queue_time = interpolation_time - process_time 
    overhead_time = total_time - interpolation_time 
    print('ALL INTERPOLATIONS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES\n({:.2f} MINUTES PROCESING, {:.2f} MINUTES QUEUING, {:.2f} MINUTES ON OVERHEADS)'.format(total_time/60., process_time/60., queue_time/60., overhead_time/60.))
else:
    print('{}/{} INTERPOLATIONS FINISHED SUCCESSFULLY IN {:.2f} MINUTES'.format(len(task_out_names)-len(failed_tasks),len(task_out_names),total_time/60.))
    print('THE FOLLOWING INTERPOLATIONS FAILED: {}'.format(failed_tasks))