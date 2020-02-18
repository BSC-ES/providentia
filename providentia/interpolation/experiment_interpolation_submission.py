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
import numpy as np
import os
import random
import subprocess
import sys
import time

#Read configuration file
from configuration import *

###--------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------###

def gather_arguments(interpolation_log_dir):

    '''gather list of arguments for all unique tasks to process, as defined in the configuration file'''

    #read defined experiments dictionary
    from defined_experiments import defined_experiments_dictionary

    #create arguments list
    arguments_list = []

    #create list to hold root directory locations of .out log files generated for each task processed
    output_log_roots = []

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

                #iterate though available species to process
                for speci_to_process in available_species:

                    #iterate through GHOST observational networks to interpolate against
                    for GHOST_network_to_interpolate_against in GHOST_networks_to_interpolate_against:
                    
                        #iterate through temporal_resolutions to output
                        for temporal_resolution_to_output in temporal_resolutions_to_output:

                            #get all GHOST network/species/temporal resolution observational files
                            obs_files = np.sort(glob.glob('/gpfs/projects/bsc32/AC_cache/obs/ghost/{}/{}/{}/{}/{}*.nc'.format(GHOST_network_to_interpolate_against, GHOST_version, temporal_resolution_to_output, speci_to_process, speci_to_process)))

                            #get all relevant experiment files
                            exp_files_all = np.sort(glob.glob('{}/{}/{}/{}/{}*.nc'.format(exp_dir, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, speci_to_process)))
                            
                            #determine if simulation generates files with ensemble member numbers or not (test first file)
                            if exp_files_all[0].split('_')[0] != speci_to_process:
                                have_ensemble_members = True
                                #if have ensemble members in filename, get all unique numbers
                                unique_ensemble_members = np.unique([f.split('{}-'.format(speci_to_process))[0][:3] for f in exp_files_all])
                                #get intersection between desired ensemble members to process and those available in directory
                                #if no members defined explicitly to process, process them all 
                                if model_ensemble_members_to_process == []:
                                    available_ensemble_members = unique_ensemble_members
                                else:
                                    available_ensemble_members = [x for x in unique_ensemble_members if x in model_ensemble_members_to_process]
                            else:
                                have_ensemble_members = False
                                #if have defined ensemble members to process, then continue as no files in this directory have ensemble member number
                                if model_ensemble_members_to_process != []:
                                    continue
                                #otherwise, proceed (tag files as ensemble member '000' for sake of code)
                                else:
                                    unique_ensemble_members = ['000']
                                                                                     
                            #iterate through available ensemble members to process
                            for ensemble_member in model_ensemble_members_to_process:

                                #limit experiment files to be just those for current ensemble member
                                if have_ensemble_members == True:
                                    exp_files = [f for f in exp_files_all if '-{}_'.format(ensemble_member) in f]
                                else:
                                    exp_files = copy.deepcopy(exp_files_all)

                                #get all observational file start dates (year and month)
                                obs_files_dates = []
                                for obs_file in obs_files:
                                    obs_file_date = obs_file.split('{}_'.format(speci_to_process))[-1].split('.nc')[0]
                                    obs_files_dates.append(obs_file_date[:6])
                                obs_files_dates = np.sort(obs_files_dates)

                                #get all experiment file start dates (year and month)
                                exp_files_dates = []   
                                for exp_file in exp_files:
                                    exp_file_date = exp_file.split('.nc')[0][-10:]
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
                                if not os.path.exists('{}/{}/{}/{}/{}/{}/{}/{}'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, ensemble_member)):
                                    os.makedirs('{}/{}/{}/{}/{}/{}/{}/{}'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, ensemble_member))

                                #iterate through intersecting yearmonths and write all current variable arguments to arguments file
                                for yearmonth in intersect_yearmonths:

                                    #append current iterative arguments to arguments list               
                                    arguments_list.append("{} {} {} {} {} {} {} {}".format(experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, ensemble_member, yearmonth))                         

                                    #append root name of .out file that will be output for each processed task
                                    output_log_roots.append('{}/{}/{}/{}/{}/{}/{}/{}/{}'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, ensemble_member, yearmonth))                           

                                    #remove previous output logs
                                    previous_logs = glob.glob('{}/{}/{}/{}/{}/{}/{}/{}/{}*'.format(interpolation_log_dir, experiment_to_process, grid_type_to_process, model_temporal_resolution_to_process, speci_to_process, GHOST_network_to_interpolate_against, temporal_resolution_to_output, ensemble_member, yearmonth))
                                    for previous_log in previous_logs:
                                        os.remove(previous_log)

    return arguments_list, output_log_roots 

###--------------------------------------------------------------------------------------------------###

def create_greasy_arguments_file(unique_ID, arguments, working_directory, arguments_dir):

    '''function which creates greasy arguments text file storing all different tasks to run by greasy'''

    #define list to store chunked argument files (to be submitted using greasy)
    argument_files = []

    #get the CPU chunk size -- set initially as miniumum number of arguments per file 
    N_arguments_per_file_minimum = copy.deepcopy(chunk_size)
    N_arguments_per_file = copy.deepcopy(chunk_size)
    #divide the number of arguments by the CPU chunk size, to determine how many argument files will be needed to submit all jobs
    N_submit_files = int(np.ceil(len(arguments)/chunk_size))
    #set argument remainder as 0 initially
    argument_remainder = 0
    #if the number of argument files is greater than the job array limit (i.e the limit on the number of argument files that can be processed simultaneously)
    #then add adjust minimum N arguments per file appropriately (i.e. split extra arguments across the maximum number of argument files evenly)        
    if N_submit_files > job_array_limit:
        #update the minimum number of arguments per file
        N_arguments_per_file_minimum = int(np.floor(len(arguments)/job_array_limit))
        N_arguments_per_file = copy.deepcopy(N_arguments_per_file_minimum)
        #if the number of extra arguments does not divide equally between all files get the remainder
        argument_remainder = int(len(arguments)%job_array_limit)
        #if have argument remainder then update N_arguments_per_file variable to be 1 greater than minimum for first file written (and for all files thereafter until  remainder is accounted for)
        if argument_remainder > 0:
            N_arguments_per_file = N_arguments_per_file_minimum + 1
            #subtract 1 from the argument remainder
            argument_remainder -= 1
        #set N submit files as N of job array limit
        N_submit_files = copy.deepcopy(job_array_limit)

    #create file which will store a list of all chunked argument filenames
    greasy_file = open('{}/{}.grz'.format(arguments_dir, unique_ID), 'w')
    #create all chunked argument filenames
    for ii in range(N_submit_files):
        argument_fname = '{}/{}_{}.txt'.format(arguments_dir,unique_ID,ii)
        argument_files.append(argument_fname)
        greasy_file.write('{}\n'.format(argument_fname))
    greasy_file.close()

    #initialise variables to count N argument files written and N lines written
    file_ii = 0
    n_lines_written = 0
    #define special characters that will need to be escaped in string written        
    special_characters = ['(',')']

    #open first arguments file for writing
    arguments_file = open(argument_files[file_ii], 'w')

    #iterate through different arguments, writing line by line to arguments file
    for arguments_ii, str_to_write in enumerate(arguments):

        #escape certain special characters in str_to_write
        for ch in special_characters:
            str_to_write = str_to_write.replace(ch,'\{}'.format(ch))
        
        #write arguments str to current file
        arguments_file.write('python -u {}/experiment_interpolation.py {}\n'.format(working_directory,str_to_write))

        #iterate n lines written    
        n_lines_written += 1
    
        #if have written sufficient arguments to file (and not currently processing last argument)
        #then close current file and open another
        if (n_lines_written == N_arguments_per_file) & (arguments_ii < (len(arguments)-1)):
            #reset n lines written
            n_lines_written = 0
            #close current arguments file
            arguments_file.close()
            #iterate n arguments files written
            file_ii += 1
            #if have argument remainder, write an extra argument to the next file
            if argument_remainder > 0:
                N_arguments_per_file = N_arguments_per_file_minimum + 1
                #subtract 1 from the argument remainder
                argument_remainder -= 1
            #otherwise N arguments to write to next file should be minimum
            else:
                N_arguments_per_file = copy.deepcopy(N_arguments_per_file_minimum)
            #open new arguments file
            arguments_file = open(argument_files[file_ii], 'w')
                
    #close current arguments file
    arguments_file.close()

###--------------------------------------------------------------------------------------------------###

def create_slurm_submission_script(unique_ID, working_directory, arguments_dir, submit_dir):

    '''function that writes a slurm submission shell script that submits a greasy job'''

    #get current BSC machine
    machine = os.environ['BSC_MACHINE']

    #create slurm_fname (unique_ID + 'sh.')
    slurm_fname = unique_ID+'.sh'
    
    #get all argument files
    argument_files = sorted(glob.glob('{}/{}_*.txt'.format(arguments_dir,unique_ID)))

    #read how many lines are in first arguments file
    with open(argument_files[0]) as f: 
        for ii, line in enumerate(f):
            pass
        N_arguments = ii + 1

    #cap the number of simultaneously running tasks to be the defined CPU chunk size  
    max_tasks = copy.deepcopy(chunk_size)

    #if the number of arguments is > capped max tasks, then set N simultaneous tasks to be the max tasks permitted
    if N_arguments > max_tasks:
        n_simultaneous_tasks = copy.deepcopy(max_tasks) 
    #else, if the number of arguments is <= capped max tasks, then set N simultaneous tasks to be N arguments
    else:
        n_simultaneous_tasks = copy.deepcopy(N_arguments)

    #create slurm submission script
    submit_file = open(submit_dir+'/'+slurm_fname,"w")

    submit_file.write("#!/bin/bash\n")
    submit_file.write("\n")
    submit_file.write("#SBATCH --job-name=PRVI_{}\n".format(unique_ID))
    submit_file.write("#SBATCH --ntasks={}\n".format(n_simultaneous_tasks))
    #fix number of nodes to be 1 (for faster execution)
    submit_file.write("#SBATCH --nodes 1\n")
    submit_file.write("#SBATCH --time=02:00:00\n")
    submit_file.write("#SBATCH --array=1-{}\n".format(len(argument_files)))
    #if machine is power9, then there are 4 threads per physical CPU, which show as 4 separate CPUs. 
    #if multithreading is enabled, ensure 4 virtual cores are actually seen as 1 CPU.
    if (machine == 'power') & (multithreading == True):
        submit_file.write("#SBATCH --cpus-per-task=4\n")
    submit_file.write("#SBATCH --qos={}\n".format(qos))
    submit_file.write("#SBATCH --output=/dev/null\n")
    submit_file.write("#SBATCH --error=/dev/null\n")
    submit_file.write("\n")
    submit_file.write("source {}/load_modules.sh\n".format(working_directory))
    submit_file.write("export GREASY_NWORKERS=$SLURM_NPROCS\n") 
    submit_file.write("export GREASY_LOGFILE={}/{}_$SLURM_ARRAY_TASK_ID.log\n".format(submit_dir,unique_ID))
    submit_file.write("arguments_store={}/{}.grz\n".format(arguments_dir,unique_ID))
    submit_file.write("argument_file=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n")
    submit_file.write("\n")
    submit_file.write("greasy $argument_file")

    #close submit file
    submit_file.close()

    return slurm_fname

###--------------------------------------------------------------------------------------------------###

if __name__ == "__main__":

    #start timer
    start = time.time()

    #define current working directory and arguments/submit/interpolation log subdirectories
    working_directory = os.getcwd()
    arguments_dir = '{}/arguments'.format(working_directory)
    submit_dir = '{}/submit'.format(working_directory)
    interpolation_log_dir = '{}/interpolation_logs'.format(working_directory)

    #set SLURM job ID as unique ID for tracking tasks defined to process in the configuration file
    unique_ID = sys.argv[1]

    #get all unique arguments to process interpolation tasks 
    arguments, output_log_roots = gather_arguments(interpolation_log_dir)

    #randomise the order of the arguments list
    random.shuffle(arguments) 
    
    #create greasy arguments file
    create_greasy_arguments_file(unique_ID, arguments, working_directory, arguments_dir)

    #create slurm submission script
    slurm_fname = create_slurm_submission_script(unique_ID, working_directory, arguments_dir, submit_dir)

    #time start of interpolation jobs
    interpolation_start = time.time()

    #submit slurm script
    submit_complete = False
    while submit_complete == False: 
        
        submit_process = subprocess.Popen(['sbatch', '-A', 'bsc32', slurm_fname], cwd=submit_dir, stdout=subprocess.PIPE)
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
        #if any jobs have failed/not finished, write them out to file
        failed_tasks=[]
        not_finished_tasks=[]
        process_times=[]
        for output_log_root in output_log_roots:
            output_log_file = glob.glob('{}_*'.format(output_log_root))
            #have an output log file? (i.e. job has finished)
            if len(output_log_file) > 0:
                output_log_file = output_log_file[0]
                process_code = int(output_log_file.split('_')[-1].split('.out')[0])
                #if process code == 0, job completed successfully
                if process_code == 0:
                    #append interpolation job process time
                    process_times.append(float(subprocess.check_output(['tail', '-1', output_log_file], encoding='utf8').strip()))
                #otherwise, job failed --> append failed log file
                else:
                    failed_tasks.append(output_log_file)
            #no output log file, therefore append to not finished list
            else:
                not_finished_tasks.append(output_log_root)

        #break out of while loop     
        all_tasks_finished = True

    #get interpolation time
    interpolation_time = (time.time() - interpolation_start)/60.

    #stop timer
    total_time = (time.time()-start)/60.

    #have 0 failed/non-finished tasks?
    if (len(failed_tasks) == 0) & (len(not_finished_tasks) == 0):
        #get queue time 
        process_time = np.max(process_times)
        queue_time = interpolation_time - process_time 
        overhead_time = total_time - interpolation_time 
        print('ALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES\n({:.2f} MINUTES PROCESING, {:.2f} MINUTES QUEUING, {:.2f} MINUTES ON OVERHEADS)'.format(len(output_log_roots), total_time, process_time, queue_time, overhead_time))
    else:
        print('{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES'.format(len(output_log_roots)-(len(not_finished_tasks)+len(failed_tasks)),len(output_log_roots),total_time))
        if len(failed_tasks) > 0:
            print('THE FOLLOWING INTERPOLATION TASKS FAILED: {}'.format(failed_tasks))
        if len(not_finished_tasks) > 0:
            print('THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}'.format(not_finished_tasks))