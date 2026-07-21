""" Class for the submission of interpolation jobs """

import copy
import glob
import math
import multiprocessing
import os
import random
import subprocess
import sys
import time

from netCDF4 import Dataset
import numpy as np
import psutil
from pprint import pprint
import yaml

from providentia.auxiliar import CURRENT_PATH, join

# get current path and providentia root path
PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

sys.path.append(join(PROVIDENTIA_ROOT, "providentia", "interpolation"))
sys.path.append(join(PROVIDENTIA_ROOT, "providentia"))

from interpolation.aux_interp import (
    get_aeronet_bin_radius_from_bin_variable,
    get_model_bin_radii,
    check_for_ghost,
)
from configuration import ProvConfiguration, load_conf

# load the defined models and species yamls
interp_models = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "interp_models.yaml"))
)
mapping_species = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "mapping_species.yaml"))
)
temporal_resolution_map = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "temporal_resolution_map.yaml"))
)
interp_print_variables = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "interpolation_fields.yaml"))
)


class SubmitInterpolation(object):
    """Class that handles the interpolation submission."""

    def __init__(self, **kwargs):
        """
        Initialize Providentia interpolation mode.

        Parameters
        ----------
         **kwargs : dict
            Optional command-line arguments that override default configuration values.
        """

        # update self with command line arguments
        self.commandline_arguments = copy.deepcopy(kwargs)

        # start timer
        self.start = time.time()

        # define current working directory and
        # arguments/greasy/interpolation log subdirectories
        self.working_directory = CURRENT_PATH
        self.arguments_dir = join(PROVIDENTIA_ROOT, "logs/interpolation/arguments")
        self.submit_dir = join(PROVIDENTIA_ROOT, "logs/interpolation/greasy_logs")
        self.interpolation_log_dir = join(
            PROVIDENTIA_ROOT, "logs/interpolation/interpolation_logs"
        )

        # initialize commandline arguments, if given
        self.provconf = ProvConfiguration(self, **self.commandline_arguments)

        print("Starting Providentia interpolation...")

        # update variables from config file
        if self.config != "":
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            else:
                if os.path.exists(join(self.config_dir, self.config)):
                    self.config = join(self.config_dir, self.config)
                    read_conf = True
            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
            else:
                error = "Error: The path to the configuration file specified in the command line does not exist."
                sys.exit(error)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            sys.exit(error)

        # update variables from config file
        if self.config != "":
            read_conf = False
            if os.path.exists(self.config):
                read_conf = True
            elif os.path.exists(join(self.config_dir, self.config)):
                self.config = join(self.config_dir, self.config)
                read_conf = True

            if read_conf:
                load_conf(self, self.config)
                self.from_conf = True
                # get the section in case it was passed on the command line
                if "section" in self.commandline_arguments:
                    # config and section defined
                    if self.commandline_arguments["section"] in self.all_sections:
                        self.sections = [self.commandline_arguments["section"]]
                    else:
                        msg = "Error: The section specified in the command line ({0}) does not exist.".format(
                            self.commandline_arguments["section"]
                        )
                        msg += "\nTip: For subsections, add the name of the parent section followed by an interpunct (·) "
                        msg += "before the subsection name (e.g. SECTIONA·Spain). Available: {0}".format(
                            self.all_sections
                        )
                        self.logger.error(msg)
                        sys.exit(1)
                # if no section passed, then get all the parent sections
                else:
                    # if no parent section names are found throw an error
                    if len(self.parent_section_names) == 0:
                        error = "Error: No sections were found in the configuration file, make sure to name them using square brackets."
                        self.logger.error(error)
                        sys.exit(1)
                    self.sections = self.parent_section_names
            else:
                error = "Error: The path to the configuration file specified in the command line does not exist."
                self.logger.error(error)
                sys.exit(1)
        else:
            error = "Error: No configuration file found. The path to the config file must be added as an argument."
            self.logger.error(error)
            sys.exit(1)

    def run(self):
        """Execute the Providentia interpolation workflow for all configured sections."""

        for section_ind, section in enumerate(self.sections):
            print(f"\nStarting {section} section interpolation...\n")

            # update for new section parameters
            self.section = section
            self.section_opts = self.sub_opts[self.section]

            # dictionary that stores used interpolation variables
            self.interpolation_variables = {}

            # update self with section variables (if not passed via command line)
            for k, val in self.section_opts.items():
                # n_cpus is always passed via command line, either as default or explicit,
                # only read n_cpus from conf if it was not explicitly defined in command line
                if ((k not in self.commandline_arguments) 
                    or (k == "n_cpus" and self.commandline_arguments['n_cpus_explicit'] == 'false')):
                    setattr(self, k, self.provconf.parse_parameter(k, val))
            
            # if n_cpus are defined in configuration file, consider as explicit
            if "n_cpus" in self.section_opts.keys():
                setattr(self, "n_cpus_explicit", self.provconf.parse_parameter("n_cpus_explicit", True))
            
            # now all variables have been parsed, check validity of those, throwing errors where necessary
            self.provconf.check_validity()

            # print variables used, if all species are used print "All Species"
            print("Variables used for the interpolation:\n")
            for arg in interp_print_variables:
                if arg != "models":
                    print(f"{arg}: {getattr(self, arg)}")
                else:
                    print(f"{arg}:")
                    for mod, alias in getattr(self, "experiments").items():
                        if self.alias_flag:
                            print(f" - {mod} ({alias})")
                        else:
                            print(f" - {mod}")

            # define the QOS (Quality of Service) used to manage jobs on the SLURM system
            if self.machine == "mn5":
                self.qos = "gp_bsces"
            else:
                self.qos = "bsc_es"

            # initialise current line number for printing output
            self.current_line = -1

            # get all unique arguments to process interpolation tasks
            self.gather_arguments()

            # submit interpolation jobs
            if not self.interp_multiprocessing:

                # check if Greasy is installed
                is_greasy_installed = False
                try:
                    result = subprocess.run(
                        ["greasy", "-V"], stdout=subprocess.PIPE, text=True
                    )
                    if "greasy" in result.stdout:
                        is_greasy_installed = True
                except:
                    print("\nGreasy is not installed.")
                    pass
                
                if is_greasy_installed:
                    
                    print("\nUsing Greasy to manage the job submission.")

                    # create greasy arguments file
                    self.create_greasy_arguments_file()

                    # create submission script according to machine
                    self.create_slurm_submission_script()
                    self.submit_job_greasy()

                else:
                    print("Using multiprocessing to manage the job submission.")
                    self.submit_job_multiprocessing()

            else:
                print("\nUsing multiprocessing to manage the job submission.")
                self.submit_job_multiprocessing()

            # remove section variables from memory
            for k in self.section_opts:
                try:
                    vars(self).pop(k)
                except:
                    pass

            # reset domain and ensemble for new section
            self.domain = []
            self.ensemble = []

            # reinitialise default configuration variables
            # modified by commandline arguments, if given
            self.provconf = ProvConfiguration(self, **self.commandline_arguments)

            if section_ind != len(self.sections) - 1:
                print("\n" + "=" * 70)

    def gather_arguments(self):
        """Gather list of arguments for all unique tasks to process, as defined in the configuration file."""

        # create arguments list
        self.arguments = []

        # create list to hold root directory locations of .out log files generated for each task processed
        self.output_log_roots = []

        # create list of arguments to compare between the iterations
        last_arguments = None

        # get all networks in case the asterisk was used
        networks = self.get_all_networks() if self.network == ["*"] else self.network

        # get all models in case the asterisk was used
        models = (
            self.get_all_models()
            if self.experiments == {"*": "*"}
            else self.experiments
        )

        # iterate through desired model IDs and its types
        for mod_dom_ens, alias in models.items():
            model_to_process, grid_type, ensemble = mod_dom_ens.split("-")

            print("\nMODEL: {0}\n".format(alias))

            # initialise files lists
            obs_files = []
            mod_files = []
            obs_files_dates = []
            mod_files_dates = []

            # initialize model search variables
            self.mod_dir = None
            msg = ""

            # get model type and specific directories
            for model_type, model_dict in interp_models.items():
                if model_to_process in model_dict["models"]:
                    # take first functional directory
                    if self.machine != "local":
                        # for HPC machines, search in interp_models
                        for temp_mod_dir in model_dict["paths"]:
                            temp_mod_dir = join(temp_mod_dir, model_to_process)
                            if os.path.exists(temp_mod_dir):
                                self.mod_dir = temp_mod_dir
                                break
                        break
                    # for local machines, break to get model_type
                    else:
                        break

            # if local machine or if not mod_dir, get directory from data_paths
            if self.mod_dir is None:
                if self.machine != "local":
                    msg = f"The model '{model_to_process}' is in none of the model paths defined in settings/interp_models.yaml."
                    print(msg)

                mod_to_interp_path = join(self.mod_to_interp_root, model_to_process)
                if os.path.exists(mod_to_interp_path):
                    self.mod_dir = mod_to_interp_path
                else:
                    msg = f"The model '{model_to_process}' is not in {self.mod_to_interp_root}."
                    print(msg)
                    continue

            # get model bin edges
            r_edges, rho_bins = get_model_bin_radii(model_type)

            # get all resolutions in case the asterisk was used
            resolutions = (
                self.get_all_resolutions(model_to_process, grid_type)
                if self.resolution == ["*"]
                else self.resolution
            )

            # iterate through temporal_resolutions to output
            for temporal_resolution_to_output in resolutions:
                model_species_ensemblestat = []

                # get priority resolutions to look for in model files based on output resolution
                # look for same resolution first, then prioritise finer resolutions
                # if output resolution is non-instantaneous, prioritise non-instanstaneous resolutions first,
                # then instanstaneous
                # if output resolution is instantaneous, prioritise instanstaneous resolutions first,
                # then non-instanstaneous

                resolutions_to_keep = (
                    self.model_resolution
                    if self.model_resolution
                    else temporal_resolution_map[temporal_resolution_to_output]
                )

                # iterate through resolutions_to_keep until find one for which have speci_to_process (or mapped speci)
                searched_paths = {}
                for model_temporal_resolution in resolutions_to_keep:
                    searched_paths[model_temporal_resolution] = {}

                    # get all species in case the asterisk was used
                    species = (
                        self.get_all_species(
                            model_to_process, grid_type, model_temporal_resolution
                        )
                        if self.species == ["*"]
                        else self.species
                    )

                    # iterate through species to process
                    for speci_ii, speci_to_process in enumerate(species):
                        
                        have_valid_resolution = False
                        original_speci_to_process = copy.deepcopy(speci_to_process)
                        searched_paths[model_temporal_resolution][original_speci_to_process] = []

                        # before proceeding check if have already processed jobs for the model-grid_type/species/temporal_resolution_to_output pairing before
                        # if so, continue to next species
                        processed_species = False
                        for arg in self.arguments:
                            arg_list = arg.split(" ")
                            if (
                                (
                                    "-".join(arg_list[0].split("-")[:2])
                                    == "{}-{}".format(model_to_process, grid_type)
                                )
                                & (arg_list[2] == speci_to_process)
                                & (arg_list[4] == temporal_resolution_to_output)
                                & (arg_list[6] == original_speci_to_process)
                            ):
                                processed_species = True
                                break
                        if processed_species:
                            continue

                        # test if have directory for current speci_to_process
                        current_speci_path = "{}/{}/{}/{}".format(
                                self.mod_dir,
                                grid_type,
                                model_temporal_resolution,
                                speci_to_process,
                            )
                        searched_paths[model_temporal_resolution][original_speci_to_process].append(current_speci_path)

                        ensemble_stat_path = "{}/{}/{}/ensemble-stats".format(
                            self.mod_dir, grid_type, model_temporal_resolution
                        )
                        searched_paths[model_temporal_resolution][original_speci_to_process].append(ensemble_stat_path)

                        if os.path.isdir(current_speci_path):
                            have_valid_resolution = True

                        # test if have speci directory in ensemble-stats
                        elif os.path.isdir(ensemble_stat_path):
                            # get all ensemble-stats species
                            model_species_ensemblestat = list(
                                np.unique(
                                    [
                                        name.split("_")[0]
                                        for name in os.listdir(
                                            "{}/{}/{}/ensemble-stats".format(
                                                self.mod_dir,
                                                grid_type,
                                                model_temporal_resolution,
                                            )
                                        )
                                        if os.path.isdir(
                                            "{}/{}/{}/ensemble-stats/{}".format(
                                                self.mod_dir,
                                                grid_type,
                                                model_temporal_resolution,
                                                name,
                                            )
                                        )
                                    ]
                                )
                            )

                            # test if have speci_to_process in model_species_ensemblestat
                            if speci_to_process in model_species_ensemblestat:
                                have_valid_resolution = True
                            # test if have mapped speci_to_process in model_species_ensemblestat
                            elif speci_to_process in mapping_species:
                                for speci_to_map in mapping_species[speci_to_process]:
                                    if speci_to_map in model_species_ensemblestat:
                                        # if have a binned size distribution variable to map, first check if bin radius is within model's bin extents
                                        # if not, do not process species
                                        if ("vconcaerobin" in speci_to_process) or (
                                            "vconcaerobin" in speci_to_map
                                        ):
                                            # check if bin radius is within model's bin extents
                                            if "vconcaerobin" in speci_to_process:
                                                speci_to_check = copy.deepcopy(
                                                    speci_to_process
                                                )
                                            elif "vconcaerobin" in speci_to_map:
                                                speci_to_check = copy.deepcopy(
                                                    speci_to_map
                                                )
                                            bin_radius = get_aeronet_bin_radius_from_bin_variable(
                                                speci_to_check
                                            )
                                            if (bin_radius >= r_edges[0]) & (
                                                bin_radius <= r_edges[-1]
                                            ):
                                                speci_to_process = copy.deepcopy(
                                                    speci_to_map
                                                )
                                                print(
                                                    "Found speci to process on mapping",
                                                    speci_to_process,
                                                )
                                                have_valid_resolution = True
                                        else:
                                            speci_to_process = copy.deepcopy(
                                                speci_to_map
                                            )
                                            print(
                                                "Found speci to process on mapping",
                                                speci_to_process,
                                            )
                                            have_valid_resolution = True

                        # for some variables it is possible to extract the variable information by mapping to a different variable name, with a higher dimensionality
                        # this currently is implemented for 2 cases:
                        # -- 4D binned size distribution
                        # -- 4D gas variables
                        # first check if speci_to_process can be mapped
                        # then check if the variable to map to exists for the model/grid_type/resolution (these can be multiple, list order sets the priority)
                        else:
                            # get species that can do mapping for
                            available_species_to_map_from = list(mapping_species.keys())
                            # check if speci_to_process can be mapped
                            if speci_to_process in available_species_to_map_from:
                                # if it can be then check then if the variable to map to exists for the model/grid_type/resolution
                                # (these can be multiple, list order sets the priority)
                                for speci_to_map in mapping_species[speci_to_process]:
                                    speci_to_map_path = "{}/{}/{}/{}".format(
                                            self.mod_dir,
                                            grid_type,
                                            model_temporal_resolution,
                                            speci_to_map,
                                        )
                                    searched_paths[model_temporal_resolution][original_speci_to_process].append(speci_to_map_path)

                                    if os.path.isdir(speci_to_map_path):
                                        # if have a binned size distribution variable to map, first check if bin radius is within model's bin extents
                                        # if not, do not process species
                                        if ("vconcaerobin" in speci_to_process) or (
                                            "vconcaerobin" in speci_to_map
                                        ):
                                            # check if bin radius is within model's bin extents
                                            if "vconcaerobin" in speci_to_process:
                                                speci_to_check = copy.deepcopy(
                                                    speci_to_process
                                                )
                                            elif "vconcaerobin" in speci_to_map:
                                                speci_to_check = copy.deepcopy(
                                                    speci_to_map
                                                )
                                            bin_radius = get_aeronet_bin_radius_from_bin_variable(
                                                speci_to_check
                                            )
                                            if (bin_radius >= r_edges[0]) & (
                                                bin_radius <= r_edges[-1]
                                            ):
                                                speci_to_process = copy.deepcopy(
                                                    speci_to_map
                                                )
                                                print(
                                                    "Found speci to process on mapping",
                                                    speci_to_process,
                                                )
                                                have_valid_resolution = True
                                        else:
                                            speci_to_process = copy.deepcopy(
                                                speci_to_map
                                            )
                                            print(
                                                "Found speci to process on mapping",
                                                speci_to_process,
                                            )
                                            have_valid_resolution = True

                        # only proceed if have valid resolution
                        if have_valid_resolution:
                            # iterate through observational networks to interpolate against
                            for network_to_interpolate_against in networks:
                                # define if network is in GHOST format
                                self.reading_ghost = check_for_ghost(
                                    network_to_interpolate_against
                                )

                                # get all relevant observational files
                                # GHOST
                                if self.reading_ghost:
                                    obs_path = (
                                        self.ghost_root
                                        + "/{}/{}/{}/{}/{}*.nc".format(
                                            network_to_interpolate_against,
                                            self.ghost_version,
                                            temporal_resolution_to_output,
                                            original_speci_to_process,
                                            original_speci_to_process,
                                        )
                                    )
                                    obs_files = np.sort(glob.glob(obs_path))
                                # non-GHOST
                                else:
                                    obs_path = (
                                        self.nonghost_root
                                        + "/{}/{}/{}/{}*.nc".format(
                                            network_to_interpolate_against,
                                            temporal_resolution_to_output,
                                            original_speci_to_process,
                                            original_speci_to_process,
                                        )
                                    )
                                    obs_files = np.sort(glob.glob(obs_path))
                                
                                # if have no observational files then continue
                                if len(obs_files) == 0:
                                    print(
                                        f"Observation files for {network_to_interpolate_against} cannot be found for {temporal_resolution_to_output} resolution in {os.path.dirname(obs_path)}."
                                    )
                                    continue
                                else:
                                    print(
                                        f"{len(obs_files)} observation file(s) for {network_to_interpolate_against} and {temporal_resolution_to_output} resolution were found in {os.path.dirname(obs_path)}."
                                    )

                                # determine if ensemble is member or emsemble stat
                                ensemble_member = ensemble.isdigit()

                                # check if ensemble is ensemble stat and get all relevant model files
                                if not ensemble_member:
                                    ensemble_stat = True
                                    mod_path = (
                                        "{}/{}/{}/ensemble-stats/{}_{}/{}*{}.nc".format(
                                            self.mod_dir,
                                            grid_type,
                                            model_temporal_resolution,
                                            speci_to_process,
                                            ensemble,
                                            speci_to_process,
                                            ensemble,
                                        )
                                    )
                                    mod_files_all = np.sort(glob.glob(mod_path))
                                else:
                                    ensemble_stat = False
                                    mod_path = "{}/{}/{}/{}/{}*.nc".format(
                                        self.mod_dir,
                                        grid_type,
                                        model_temporal_resolution,
                                        speci_to_process,
                                        speci_to_process,
                                    )
                                    mod_files_all = np.sort(glob.glob(mod_path))

                                    # drop all analysis files ending with '_an.nc' which are not in ensemble-stats
                                    mod_files_all = np.array(
                                        [f for f in mod_files_all if "_an.nc" not in f]
                                    )

                                # if have no relevant model files then continue
                                if len(mod_files_all) == 0:
                                    print(
                                        f"Model files cannot be found for {temporal_resolution_to_output} resolution in {os.path.dirname(mod_path)}."
                                    )
                                    continue
                                else:
                                    print(
                                        f"{len(mod_files_all)} model file(s) for {temporal_resolution_to_output} resolution were found in {os.path.dirname(mod_path)}."
                                    )

                                # ensemble stat?
                                if ensemble_stat:
                                    available_ensemble = [ensemble]

                                # not ensemble stat?
                                else:
                                    # determine if simulation generates files with ensemble member numbers or not
                                    # (test first file)
                                    # if the ensemble number is in the nc file name next to species
                                    if (
                                        mod_files_all[0]
                                        .split("/")[-1]
                                        .rsplit("_", 1)[0]
                                        != speci_to_process
                                    ):
                                        have_ensemble_members = True
                                        # if have ensemble members in filename, get all unique numbers
                                        unique_ensemble_members = np.unique(
                                            [
                                                f.split(
                                                    "/{}-".format(speci_to_process)
                                                )[-1][:3]
                                                for f in mod_files_all
                                            ]
                                        )
                                        # get intersection between desired ensemble members to process and those
                                        # available in directory
                                        # if no members defined explicitly to process, process them all
                                        if ensemble == "allmembers":
                                            available_ensemble = unique_ensemble_members
                                        else:
                                            if ensemble in unique_ensemble_members:
                                                available_ensemble = [ensemble]
                                            else:
                                                continue
                                    # if there's no ensemble number in the file name
                                    else:
                                        have_ensemble_members = False
                                        # if have defined ensemble members to process, then continue as no files in this
                                        # directory have ensemble member number
                                        if ensemble not in ["allmembers", "000"]:
                                            continue
                                        # otherwise, proceed (tag files as ensemble member '000' for sake of operation)
                                        else:
                                            available_ensemble = ["000"]

                                # iterate through available ensemble to process
                                for available_ens in available_ensemble:
                                    # limit model files to be just those for specific ensemble member
                                    # (where neccessary)
                                    mod_file_speci = copy.deepcopy(speci_to_process)
                                    mod_files = copy.deepcopy(mod_files_all)
                                    if not ensemble_stat:
                                        if have_ensemble_members:
                                            mod_file_speci = "{}-{}".format(
                                                speci_to_process, available_ens
                                            )
                                            mod_files = np.sort(
                                                [
                                                    f
                                                    for f in mod_files_all
                                                    if "{}_".format(mod_file_speci) in f
                                                ]
                                            )

                                    # get all observational file start dates (year and month)
                                    obs_files_dates = []
                                    for obs_file in obs_files:
                                        obs_file_date = (
                                            obs_file.split(
                                                "{}_".format(original_speci_to_process)
                                            )[-1]
                                            .split("_")[0]
                                            .split(".nc")[0]
                                        )
                                        obs_files_dates = np.append(
                                            obs_files_dates, obs_file_date[:6]
                                        )

                                    # get all model file start dates (year and month)
                                    mod_files_dates = []
                                    for mod_file in mod_files:
                                        mod_file_date = (
                                            mod_file.split(
                                                "{}_".format(mod_file_speci)
                                            )[-1]
                                            .split("_")[0]
                                            .split(".nc")[0]
                                        )
                                        mod_files_dates = np.append(
                                            mod_files_dates, mod_file_date[:6]
                                        )

                                    # remove observational files outside date ranges
                                    obs_files_ii = np.array(
                                        [
                                            obs_files_ii
                                            for obs_files_ii, obs_files_date in enumerate(
                                                obs_files_dates
                                            )
                                            if (
                                                (
                                                    self.start_date == "*"
                                                    or int(obs_files_date)
                                                    >= int(self.start_date)
                                                )
                                                and (
                                                    self.end_date == "*"
                                                    or int(obs_files_date)
                                                    < int(self.end_date)
                                                )
                                            )
                                        ]
                                    )
                                    if len(obs_files_ii) == 0:
                                        obs_files = []
                                        obs_files_dates = []
                                        continue
                                    obs_files = obs_files[obs_files_ii]
                                    obs_files_dates = obs_files_dates[obs_files_ii]

                                    # remove model files outside date ranges
                                    mod_files_ii = np.array(
                                        [
                                            mod_files_ii
                                            for mod_files_ii, mod_files_date in enumerate(
                                                mod_files_dates
                                            )
                                            if (
                                                (
                                                    self.start_date == "*"
                                                    or int(mod_files_date)
                                                    >= int(self.start_date)
                                                )
                                                and (
                                                    self.end_date == "*"
                                                    or int(mod_files_date)
                                                    < int(self.end_date)
                                                )
                                            )
                                        ]
                                    )
                                    if len(mod_files_ii) == 0:
                                        mod_files = []
                                        mod_files_dates = []
                                        continue
                                    mod_files = mod_files[mod_files_ii]
                                    mod_files_dates = mod_files_dates[mod_files_ii]

                                    # get intersection of file yearmonths between observations and model
                                    intersect_yearmonths = np.intersect1d(
                                        obs_files_dates, mod_files_dates
                                    )

                                    # if have no intersecting months, continue
                                    if len(intersect_yearmonths) == 0:
                                        continue

                                    # create Providentia model code (modid-region-ensembleoption)
                                    prov_mod_code = "{}-{}-{}".format(
                                        model_to_process, grid_type, available_ens
                                    )

                                    # create directories to store slurm output/error logs for interpolation task of
                                    # specific combination of iterated variables (if does not already exist)
                                    if not os.path.exists(
                                        "{}/{}/{}/{}/{}".format(
                                            self.interpolation_log_dir,
                                            prov_mod_code,
                                            original_speci_to_process,
                                            network_to_interpolate_against,
                                            temporal_resolution_to_output,
                                        )
                                    ):
                                        os.makedirs(
                                            "{}/{}/{}/{}/{}".format(
                                                self.interpolation_log_dir,
                                                prov_mod_code,
                                                original_speci_to_process,
                                                network_to_interpolate_against,
                                                temporal_resolution_to_output,
                                            )
                                        )

                                    # iterate through intersecting yearmonths and write all current variable arguments
                                    # to arguments file
                                    for yearmonth in intersect_yearmonths:
                                        # append current iterative arguments to arguments list
                                        self.arguments.append(
                                            "{} {} {} {} {} {} {} {}".format(
                                                prov_mod_code,
                                                model_temporal_resolution,
                                                speci_to_process,
                                                network_to_interpolate_against,
                                                temporal_resolution_to_output,
                                                yearmonth,
                                                original_speci_to_process,
                                                self.slurm_job_id,
                                            )
                                        )

                                        # append root name of .out file that will be output for each processed task
                                        self.output_log_roots.append(
                                            "{}/{}/{}/{}/{}/{}".format(
                                                self.interpolation_log_dir,
                                                prov_mod_code,
                                                original_speci_to_process,
                                                network_to_interpolate_against,
                                                temporal_resolution_to_output,
                                                yearmonth,
                                            )
                                        )

                                        # remove previous output logs
                                        previous_logs = glob.glob(
                                            "{}/{}/{}/{}/{}/{}*".format(
                                                self.interpolation_log_dir,
                                                prov_mod_code,
                                                original_speci_to_process,
                                                network_to_interpolate_against,
                                                temporal_resolution_to_output,
                                                yearmonth,
                                            )
                                        )
                                        for previous_log in previous_logs:
                                            os.remove(previous_log)

            # boolean that says if there was any new valid data in the iteration
            new_arguments = last_arguments != self.arguments

            # if list is empty or have no arguments after iteration, return message stating that
            if len(self.arguments) == 0 or not new_arguments:
                msg = "\nNO INTERSECTING OBSERVATIONAL AND EXPERIMENT DATA FOR INTERPOLATION. \n"
                if not have_valid_resolution:
                    # show paths where we have searched for data
                    msg += f"\nSearched paths for model data: {pprint(searched_paths, width=120)}"
            else:
                msg = "\n***INTERSECTING OBSERVATIONAL AND EXPERIMENTAL DATA IS AVAILABLE FOR INTERPOLATION.***"

                # add a warning if nco is not installed
                _, output = subprocess.getstatusoutput("ncks")
                if "not found" in output:
                    msg += "\n\nNCO could not be found, please install it in your system with sudo apt install nco (Debian/Ubuntu) or brew install nco (macOS)."

            print(msg)

            # get the arguments from the last iteration
            last_arguments = copy.deepcopy(self.arguments)

        # if have no arguments for all models, return message stating that
        if len(self.arguments) == 0:
            error = "\nINTERPOLATION CANNOT BE DONE FOR ANY EXPERIMENT"
            sys.exit(error)

        # randomise the order of the arguments list
        random.shuffle(self.arguments)

    def create_greasy_arguments_file(self):
        """Create greasy arguments text file storing all different tasks to run by greasy."""

        # define list to store chunked argument files (to be submitted using greasy)
        argument_files = []

        # get the CPU chunk size -- set initially as miniumum number of arguments per file
        n_arguments_per_file_minimum = copy.deepcopy(int(self.interp_chunk_size))
        n_arguments_per_file = copy.deepcopy(int(self.interp_chunk_size))

        # divide the number of arguments by the CPU chunk size, to determine how many argument files will be needed to
        # submit all jobs
        n_submit_files = int(np.ceil(len(self.arguments) / int(self.interp_chunk_size)))

        # set argument remainder as 0 initially
        argument_remainder = 0

        # if the number of argument files is greater than the job array limit (i.e the limit on the number of argument
        # files that can be processed simultaneously)
        # then add adjust minimum N arguments per file appropriately (i.e. split extra arguments across the maximum
        # number of argument files evenly)
        if n_submit_files > int(self.interp_job_array_limit):
            # update the minimum number of arguments per file
            n_arguments_per_file_minimum = int(
                np.floor(len(self.arguments) / int(self.interp_job_array_limit))
            )
            n_arguments_per_file = copy.deepcopy(n_arguments_per_file_minimum)

            # if the number of extra arguments does not divide equally between all files get the remainder
            argument_remainder = int(
                len(self.arguments) % int(self.interp_job_array_limit)
            )

            # if have argument remainder then update n_arguments_per_file variable to be 1 greater than minimum for
            # first file written (and for all files thereafter until  remainder is accounted for)
            if argument_remainder > 0:
                n_arguments_per_file = n_arguments_per_file_minimum + 1
                # subtract 1 from the argument remainder
                argument_remainder -= 1

            # set N submit files as N of job array limit
            n_submit_files = copy.deepcopy(int(self.interp_job_array_limit))

        # create file which will store a list of all chunked argument filenames
        greasy_file = open(
            "{}/{}.grz".format(self.arguments_dir, self.slurm_job_id), "w"
        )

        # create all chunked argument filenames
        for ii in range(n_submit_files):
            argument_fname = "{}/{}_{}.txt".format(
                self.arguments_dir, self.slurm_job_id, ii
            )
            argument_files.append(argument_fname)
            greasy_file.write("{}\n".format(argument_fname))
        greasy_file.close()

        # initialise variables to count N argument files written and N lines written
        file_ii = 0
        n_lines_written = 0

        # define special characters that will need to be escaped in string written
        special_characters = ["(", ")"]

        # open first arguments file for writing
        arguments_file = open(argument_files[file_ii], "w")

        # iterate through different arguments, writing line by line to arguments file
        for arguments_ii, str_to_write in enumerate(self.arguments):
            # escape certain special characters in str_to_write
            for ch in special_characters:
                str_to_write = str_to_write.replace(ch, "\{}".format(ch))

            # write arguments str to current file
            command = (
                "python -u {}/interpolation/experiment_interpolation.py {}\n".format(
                    self.working_directory, str_to_write
                )
            )
            if self.machine == "nord4":
                command = "nord3_singu_es {}".format(command)
            arguments_file.write(command)

            # iterate n lines written
            n_lines_written += 1

            # if have written sufficient arguments to file (and not currently processing last argument)
            # then close current file and open another
            if (n_lines_written == n_arguments_per_file) & (
                arguments_ii < (len(self.arguments) - 1)
            ):
                # reset n lines written
                n_lines_written = 0

                # close current arguments file
                arguments_file.close()

                # iterate n arguments files written
                file_ii += 1

                # if have argument remainder, write an extra argument to the next file
                if argument_remainder > 0:
                    n_arguments_per_file = n_arguments_per_file_minimum + 1
                    # subtract 1 from the argument remainder
                    argument_remainder -= 1
                # otherwise N arguments to write to next file should be minimum
                else:
                    n_arguments_per_file = copy.deepcopy(n_arguments_per_file_minimum)

                # open new arguments file
                arguments_file = open(argument_files[file_ii], "w")

        # close current arguments file
        arguments_file.close()

    def create_slurm_submission_script(self):
        """Write a slurm submission shell script that submits a greasy job."""

        # create job_fname (unique_ID + 'sh.')
        self.job_fname = self.slurm_job_id + ".sh"

        # get all argument files
        argument_files = sorted(
            glob.glob("{}/{}_*.txt".format(self.arguments_dir, self.slurm_job_id))
        )

        # read how many lines are in first arguments file
        with open(argument_files[0]) as f:
            for ii, line in enumerate(f):
                pass
            N_arguments = ii + 1

        # cap the number of simultaneously running tasks to be the defined CPU chunk size
        max_tasks = copy.deepcopy(int(self.interp_chunk_size))

        # if the number of arguments is > capped max tasks, then set N simultaneous tasks to be the max tasks permitted
        if N_arguments > max_tasks:
            n_simultaneous_tasks = copy.deepcopy(max_tasks)
        # else, if the number of arguments is <= capped max tasks, then set N simultaneous tasks to be N arguments
        else:
            n_simultaneous_tasks = copy.deepcopy(N_arguments)

        # create slurm submission script
        submit_file = open(self.submit_dir + "/" + self.job_fname, "w")

        submit_file.write("#!/bin/bash\n")
        submit_file.write("\n")
        submit_file.write("#SBATCH --job-name=PRVI_{}\n".format(self.slurm_job_id))
        submit_file.write("#SBATCH --ntasks={}\n".format(n_simultaneous_tasks))
        # fix number of nodes to be 1 (for faster execution)
        submit_file.write("#SBATCH --nodes=1\n")
        submit_file.write("#SBATCH --time=48:00:00\n")
        submit_file.write("#SBATCH --array=1-{}\n".format(len(argument_files)))
        submit_file.write("#SBATCH --qos={}\n".format(self.qos))
        # submit_file.write("#SBATCH --output=/dev/null\n") # decomment when debugging
        # submit_file.write("#SBATCH --error=/dev/null\n")
        if self.machine in [
            "mn5",
            "nord4",
        ]:  # TODO when checking if debug works check this
            submit_file.write("#SBATCH --account=bsc32\n")
            submit_file.write(
                "#SBATCH --ntasks-per-node={}\n".format(n_simultaneous_tasks)
            )
            submit_file.write("#SBATCH --cpus-per-task=1\n")
            submit_file.write("\n")
        else:
            submit_file.write("\n")
            submit_file.write(
                "source {}/bin/load_modules.sh\n".format(PROVIDENTIA_ROOT)
            )
        submit_file.write("export GREASY_NWORKERS=$SLURM_NPROCS\n")
        submit_file.write(
            "export GREASY_LOGFILE={}/{}_$SLURM_ARRAY_TASK_ID.log\n".format(
                self.submit_dir, self.slurm_job_id
            )
        )
        submit_file.write("export SLURM_CPU_BIND=none\n")
        submit_file.write(
            "arguments_store={}/{}.grz\n".format(self.arguments_dir, self.slurm_job_id)
        )
        submit_file.write(
            "argument_file=$(cat $arguments_store | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')\n"
        )
        submit_file.write("\n")
        submit_file.write("greasy $argument_file")

        # close submit file
        submit_file.close()

    def submit_job_greasy(self):
        """Submit and monitor interpolation jobs using Greasy/SLURM."""

        # time start of interpolation jobs
        interpolation_start = time.time()

        # submit slurm script
        submit_complete = False
        while not submit_complete:
            submit_process = subprocess.Popen(
                ["sbatch", self.job_fname],
                cwd=self.submit_dir,
                stdout=subprocess.PIPE,
            )
            submit_return_code = submit_process.returncode

            # if sbatch fails, time out for 60 seconds and then try again
            if submit_return_code != 0:
                time.sleep(60)
                continue
            else:
                submit_complete = True

            # take a 1 second pause between submittals (to help slurm scheduler)
            time.sleep(1)

        # monitor number of jobs in queue (every 10 seconds) until there are 0 left in the squeue
        all_tasks_finished = False
        # flag for monitoring jobs that were submitted in n3
        job_entered = False

        while not all_tasks_finished:
            cmd = ["squeue", "-h", "-n", "PRVI_{}".format(self.slurm_job_id)]
            squeue_process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, encoding="utf8"
            )
            squeue_status = squeue_process.communicate()[0]
            n_jobs_in_queue = len(squeue_status.split("\n")[:-1])
            # if number of jobs in queue > 0, then sleep for 10
            # seconds and then check again how many jobs there are in queue
            if (self.machine in ("mn5", "nord4")) and (n_jobs_in_queue > 0):
                time.sleep(10)
                continue

            # if no more jobs in the squeue, then now check the outcome of all the jobs
            # if any jobs have failed/not finished, write them out to file
            failed_tasks = []
            not_finished_tasks = []
            process_times = []
            for output_log_root in self.output_log_roots:
                output_log_file = glob.glob("{}_*".format(output_log_root))
                # have an output log file? (i.e. job has finished)
                if len(output_log_file) > 0:
                    output_log_file = output_log_file[0]
                    process_code = int(output_log_file.split("_")[-1].split(".out")[0])
                    # if process code == 0, job completed successfully
                    if process_code == 0:
                        # append interpolation job process time
                        process_times.append(
                            float(
                                subprocess.check_output(
                                    ["tail", "-1", output_log_file], encoding="utf8"
                                ).strip()
                            )
                        )
                    # otherwise, job failed --> append failed log file
                    else:
                        failed_tasks.append(output_log_file)
                # no output log file, therefore append to not finished list
                else:
                    not_finished_tasks.append(output_log_root)

            # break out of while loop
            all_tasks_finished = True

        # get interpolation time
        interpolation_time = (time.time() - interpolation_start) / 60.0

        # stop timer
        total_time = (time.time() - self.start) / 60.0

        # have 0 failed/non-finished tasks?
        if (len(failed_tasks) == 0) & (len(not_finished_tasks) == 0):
            # get queue time
            process_time = np.max(process_times)
            queue_time = interpolation_time - process_time
            overhead_time = total_time - interpolation_time
            print(
                "\nALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES"
                "({:.2f} MINUTES PROCESSING, {:.2f} MINUTES QUEUING, {:.2f} MINUTES ON OVERHEADS)".format(
                    len(self.output_log_roots),
                    total_time,
                    process_time,
                    queue_time,
                    overhead_time,
                )
            )
        else:
            print(
                "\n{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES".format(
                    len(self.output_log_roots)
                    - (len(not_finished_tasks) + len(failed_tasks)),
                    len(self.output_log_roots),
                    total_time,
                )
            )
            if len(failed_tasks) > 0:
                print(
                    "THE FOLLOWING INTERPOLATION TASKS FAILED: {}".format(failed_tasks)
                )
            if len(not_finished_tasks) > 0:
                print(
                    "THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}".format(
                        not_finished_tasks
                    )
                )

    def submit_job_multiprocessing(self):
        """Submit interpolation jobs using multiprocessing pool."""

        # set resource usage parameters for estimating safe number of pool workers
        self.per_worker_mem_gb = self.guess_peak_RSS_per_worker()
        self.cpu_fraction_limit = 0.75
        self.mem_fraction_limit = 0.75

        # set run commands
        self.commands = [
            "python -u {}/interpolation/experiment_interpolation.py {}".format(
                self.working_directory, argument
            )
            for argument in self.arguments
        ]

        # if n_cpus hasn't been defined, use 1 or half of the available CPUS to
        # avoid having to kill other processes locally
        if self.machine == "local":
            # use value passed through --n_cpus in terminal or configuration file
            if self.n_cpus_explicit:
                n_cpus = self.n_cpus
                msg = f"Using {n_cpus} CPUs."
            # otherwise estimate the number of safe pool workers
            else:
                n_cpus = self.guess_pool_workers()
                msg = f"Using {n_cpus} out of {self.n_cpus} available CPUs to"
                msg += " ensure that other processes keep running. \nIf you encounter any problems"
                msg += (
                    " consider reducing the number of CPUS by running Providentia using"
                )
                msg += " --n_cpus (./bin/providentia --n_cpus=2)."
        else:
            n_cpus = self.n_cpus
            msg = f"Using {n_cpus} CPUs."

        print(msg)

        # cap number of cpus to not be larger than number of tasks
        if n_cpus > len(self.commands):
            old_n_cpus = copy.deepcopy(n_cpus)
            n_cpus = len(self.commands)
            print(
                f"Capping {old_n_cpus} CPUs to {n_cpus} since there are only {n_cpus} tasks to process."
            )

        # add print space
        print()

        # get total number of jobs
        self.total_jobs = len(self.commands)

        # launch interpolation
        # if have no swap memory or have just 1 cpu, run in serial without multiprocessing
        if n_cpus == 1:
            for i, cmd in enumerate(self.commands, start=1):
                self.run_command(cmd)
                print(f"{i}/{self.total_jobs} jobs completed", flush=True)
        # otherwise, use multiprocessing pool
        else:
            with multiprocessing.Pool(n_cpus) as pool:
                for i, r in enumerate(
                    pool.imap_unordered(self.resource_safe_job, self.commands), start=1
                ):
                    print(f"{i}/{self.total_jobs} jobs completed", flush=True)

        # stop timer
        total_time = (time.time() - self.start) / 60.0

        # if no more jobs in the squeue, then now check the outcome of all the jobs
        # if any jobs have failed/not finished, write them out to file
        failed_tasks = []
        not_finished_tasks = []
        for output_log_root in self.output_log_roots:
            output_log_file = glob.glob("{}_*".format(output_log_root))
            # have an output log file? (i.e. job has finished)
            if len(output_log_file) > 0:
                output_log_file = output_log_file[0]
                process_code = int(output_log_file.split("_")[-1].split(".out")[0])
                # if process code is not 0, job failed
                if process_code != 0:
                    failed_tasks.append(output_log_file)
            # no output log file, therefore append to not finished list
            else:
                not_finished_tasks.append(output_log_root)

        # have 0 failed/non-finished tasks?
        if (len(failed_tasks) == 0) & (len(not_finished_tasks) == 0):
            print(
                "\nALL {} INTERPOLATION TASKS COMPLETED SUCCESSFULLY IN {:.2f} MINUTES".format(
                    len(self.output_log_roots), total_time
                )
            )
        else:
            print(
                "\n{}/{} INTERPOLATION TASKS FINISHED SUCCESSFULLY IN {:.2f} MINUTES".format(
                    len(self.output_log_roots)
                    - (len(not_finished_tasks) + len(failed_tasks)),
                    len(self.output_log_roots),
                    total_time,
                )
            )
            if len(failed_tasks) > 0:
                print(
                    "THE FOLLOWING INTERPOLATION TASKS FAILED: {}".format(failed_tasks)
                )
            if len(not_finished_tasks) > 0:
                print(
                    "THE FOLLOWING INTERPOLATION TASKS DID NOT FINISH: {}".format(
                        not_finished_tasks
                    )
                )

    def guess_peak_RSS_per_worker(self):
        """Get conservative estimate of peak RSS per worker in GB."""

        # define conservative factors for estimation for peak RSS in GB
        array_overhead_factor = 1.5
        temp_array_factor = 3.0
        if self.operating_system == "Mac":
            process_overhead_gb = 0.15
        else:
            process_overhead_gb = 0.08
        hdf5_cache_gb = 0.1

        # initialise peak RSS as 0
        max_worker_rss_gb = 0.0

        # iterate through all jobs to submit and get peak RSS for a job
        for argument in self.arguments:
            worker_payload_gb = 0.0

            argument = argument.split(" ")
            prov_mod_code = argument[0]
            model_temporal_resolution = argument[1]
            speci_to_process = argument[2]
            yearmonth = argument[5]
            model_to_process, grid_type, ensemble = prov_mod_code.split("-")
            ensemble_member = ensemble.isdigit()

            if ensemble_member:
                all_model_files = np.sort(
                    glob.glob(
                        f"{self.mod_dir}/{grid_type}/{model_temporal_resolution}/"
                        f"{speci_to_process}/{speci_to_process}*{yearmonth}*.nc"
                    )
                )
                all_model_files = [f for f in all_model_files if "_an.nc" not in f]
                model_files = [
                    f for f in all_model_files if f"{speci_to_process}-{ensemble}_" in f
                ]
                if not model_files:
                    model_files = all_model_files
            else:
                model_files = np.sort(
                    glob.glob(
                        f"{self.mod_dir}/{grid_type}/{model_temporal_resolution}/"
                        f"ensemble-stats/{speci_to_process}_{ensemble}/"
                        f"{speci_to_process}*{yearmonth}*{ensemble}.nc"
                    )
                )

            for model_file in model_files:
                with Dataset(model_file, "r") as nc:
                    var = nc.variables[speci_to_process]
                    dims = [nc.dimensions[d].size for d in var.dimensions]
                    n_elements = math.prod(dims)
                    dtype_size = np.dtype(var.dtype).itemsize
                    worker_payload_gb += (n_elements * dtype_size) / (1024**3)

            # Inflate to peak RSS
            peak_worker_gb = (
                worker_payload_gb * array_overhead_factor * temp_array_factor
                + process_overhead_gb
                + hdf5_cache_gb
            )

            # get peak worker RSS in GB
            max_worker_rss_gb = max(max_worker_rss_gb, peak_worker_gb)

        return max_worker_rss_gb

    def guess_pool_workers(self):
        """
        Estimate a safe number of multiprocessing workers for the system.

        Returns
        -------
        n_workers : int
            Estimated number of worker processes that can safely run in parallel
            Always at least 1.
        """

        # static system capacity
        physical_cores = psutil.cpu_count(logical=False) or 1
        vm = psutil.virtual_memory()

        # per-worker estimates
        per_worker_mem_gb = self.per_worker_mem_gb

        # CPU-based limit
        cpu_capacity_workers = int(physical_cores * self.cpu_fraction_limit)
        cpu_capacity_workers = min(cpu_capacity_workers, physical_cores)
        cpu_capacity_workers = max(1, cpu_capacity_workers)

        # memory-based limit (peak RSS)
        available_mem_gb = vm.available / (1024**3)
        usable_mem_gb = max(0.0, available_mem_gb)

        if per_worker_mem_gb > 0:
            mem_capacity_workers = int(usable_mem_gb / per_worker_mem_gb)
        else:
            mem_capacity_workers = physical_cores

        mem_capacity_workers = max(1, mem_capacity_workers)

        # final worker count
        n_workers = min(cpu_capacity_workers, mem_capacity_workers, physical_cores)

        return max(1, n_workers)

    def resource_safe_job(self, cmd, check_interval=1.0, stagger_range=(0.0, 5.0)):
        """
        Manage multiprocessing worker submission with hysteresis.
        If system is overloaded it will wait until the system is stable before submitting more jobs.

        Parameters
        ----------
        cmd : str
            Command to execute as the job.
        check_interval : float, optional
            Time in seconds between resource checks (default is 1.0).
        stagger_range : tuple of float, optional
            Range (min, max) in seconds to randomly stagger job start
            to avoid simultaneous spikes (default is (0.0, 5.0)).
        """

        # set cpu, memory and swap limits
        cpu_enter = self.cpu_fraction_limit
        cpu_exit = self.cpu_fraction_limit * 0.9
        mem_enter = self.mem_fraction_limit
        mem_exit = self.mem_fraction_limit * 0.85
        swap_panic_limit = 0.95

        # initialise hysteresis counters
        overload_enter_count = 3
        recovery_exit_count = 3
        overload_count = 0
        recovery_count = 0
        paused = False
        max_pause_seconds = 120
        pause_start = None

        # while loop to trap jobs while there is an overload
        while True:
            # get current swap, cpu and memory percent
            swap_frac = psutil.swap_memory().percent / 100.0
            cpu = psutil.cpu_percent(interval=None) / 100.0
            mem = 1.0 - (
                psutil.virtual_memory().available / psutil.virtual_memory().total
            )

            # If have extremely high swap (and are not on Mac) then pause for 30 seconds
            if (self.operating_system != "Mac") and (swap_frac > swap_panic_limit):
                print(
                    f"Warning: Swap memory is extremely high: {swap_frac*100:.1f}% — sleeping 30s",
                    flush=True,
                )
                time.sleep(30)
                continue

            # check if system is overloaded
            # 1. system exceeds cpu limits
            # 2. system exeeds memory limits
            overloaded = (cpu > cpu_enter) or (mem > mem_enter)

            # hysteresis counters
            if overloaded:
                overload_count += 1
                recovery_count = 0
            else:
                overload_count = 0
                recovery_count += 1

            # if system is overloaded then pause submissions
            if (not paused) and (overload_count >= overload_enter_count):
                paused = True
                pause_start = time.time()
                print(
                    "Warning: Sustained overload. Pausing job submissions.", flush=True
                )

            # check if can resume submissions
            if paused:
                # Force progress after max pause
                if (pause_start) and (time.time() - pause_start > max_pause_seconds):
                    print(
                        "Warning: Forcing job submission to avoid deadlock.", flush=True
                    )
                    paused = False
                    break
                healthy = (cpu < cpu_exit and mem < mem_exit) or (cpu < 0.25)
                if healthy:
                    recovery_count += 1
                else:
                    recovery_count = 0
                if recovery_count >= recovery_exit_count:
                    paused = False
                    print(
                        "Warning: System has recovered. Resuming job submissions.",
                        flush=True,
                    )

            # break out of loop if not paused
            if not paused:
                break

            # if paused then sleep a little
            time.sleep(check_interval)

        # Stagger before starting job
        time.sleep(random.uniform(*stagger_range))
        self.run_command(cmd)

    def run_command(self, commands):
        """
        Function to submit interpolation job.

        Parameters
        ----------
        commands : str
            Command string to be executed for job submission.
        """

        arguments_list = commands.strip().split()
        print(
            "Submitting job with arguments: [{} {} {} {} {}]".format(
                arguments_list[3],
                arguments_list[5],
                arguments_list[6],
                arguments_list[7],
                arguments_list[8],
            ),
            flush=True,
        )
        if self.machine == "nord4":
            arguments_list.insert(0, "nord3_singu_es")
        result = subprocess.run(arguments_list, capture_output=True, text=True)

        if result.returncode != 0:
            error = result.stderr
            if error == "":
                error = "Unknown error"
            print(
                f"Error in submission using the arguments: {result.args[3:-1]}: "
                f"\nError: {error}"
                f"\nReturn code: {result.returncode}",
                flush=True,
            )

    def get_all_models(self):
        """
        Retrieve all available models.

        Returns
        -------
        models_dict : dict
            A dictionary of models.
        """

        models = []

        # from interp_models (remothe machine)
        if self.machine != "local":
            for model_dict in interp_models.values():
                for mod_id in model_dict["models"]:
                    for temp_mod_dir in model_dict["paths"]:
                        mod_to_interp_path = join(temp_mod_dir, mod_id)
                        if os.path.exists(mod_to_interp_path):
                            for dom in os.listdir(mod_to_interp_path):
                                models.append(
                                    f"{mod_id}-{dom}-{self.default_values['ensemble'][0]}"
                                )

        # from data_paths (local and remote machine)
        for mod_id in os.listdir(self.mod_to_interp_root):
            if not mod_id.startswith("."):
                mod_to_interp_path = join(self.mod_to_interp_root, mod_id)
                for dom in os.listdir(mod_to_interp_path):
                    models.append(
                        f"{mod_id}-{dom}-{self.default_values['ensemble'][0]}"
                    )

        # create dictionary from list
        models_dict = dict(zip(models, models))

        return models_dict

    def get_all_resolutions(self, mod_id, domain):
        """
        Retrieve all available resolutions.

        Parameters
        ----------
        mod_id : str
            Identifier of the model.
        domain : str
            Domain of the model to retrieve resolutions for.

        Returns
        -------
        resolutions : list
            A list of resolutions.
        """

        resolutions = []

        # from interp_models (remothe machine)
        if self.machine != "local":
            for model_dict in interp_models.values():
                if mod_id in model_dict["models"]:
                    for temp_mod_dir in model_dict["paths"]:
                        domain_dir = join(temp_mod_dir, mod_id, domain)
                        if os.path.exists(domain_dir):
                            resolutions += os.listdir(domain_dir)

        # from data_paths (local and remote machine)
        domain_dir = join(self.mod_to_interp_root, mod_id, domain)
        if os.path.exists(domain_dir):
            resolutions += os.listdir(domain_dir)

        return resolutions

    def get_all_species(self, mod_id, domain, resolution):
        """
        Retrieve all available species.

        Parameters
        ----------
        mod_id : str
            Identifier of the model.
        domain : str
            Domain of the model to retrieve species for.
        resolution : str
            Resolution of the model to retrieve species for.

        Returns
        -------
        species : list
            A list of species.
        """

        species = []

        # from interp_models (remothe machine)
        if self.machine != "local":
            for model_dict in interp_models.values():
                if mod_id in model_dict["models"]:
                    for temp_mod_dir in model_dict["paths"]:
                        resolution_dir = join(temp_mod_dir, mod_id, domain, resolution)
                        if os.path.exists(resolution_dir):
                            species += os.listdir(resolution_dir)

        # from data_paths (local and remote machine)
        resolution_dir = join(self.mod_to_interp_root, mod_id, domain, resolution)
        if os.path.exists(resolution_dir):
            species += os.listdir(resolution_dir)

        return species

    def get_all_networks(self):
        """
        Retrieve all available networks.

        Returns
        -------
        networks : list
            A list of networks.
        """

        if self.network_type:
            # exit if the value is neither GHOST or non-GHOST
            network_type = str(self.network_type).lower()
            if network_type not in ["ghost", "non-ghost"]:
                error = f"Error: Invalid 'network_type': '{self.network_type}'. Expected 'ghost' or 'non-ghost'."
                self.logger.error(error)
                sys.exit(1)

            self.reading_ghost = network_type == "ghost"
        else:
            # get user input to know which kind of network wants
            while True:
                download_source = input(
                    "\nDo you want to interpolate all the GHOST networks? (Otherwise all the non-GHOST networks will be interpolated) ([y]/n): "
                ).lower()
                if download_source in ["", "y", "n"]:
                    break

            # get the boolean value from the answer of the user
            self.reading_ghost = download_source in ["", "y"]

        networks = []

        if self.reading_ghost:
            # from data_paths (local and remote machine)
            if os.path.exists(self.ghost_root):
                networks = os.listdir(self.ghost_root)
        else:
            # from data_paths (local and remote machine)
            if os.path.exists(self.nonghost_root):
                for dir1 in os.listdir(self.nonghost_root):
                    for dir2 in os.listdir(join(self.nonghost_root, dir1)):
                        networks.append(f"{dir1}/{dir2}")

        return networks


def main(**kwargs):
    """
    Initialize the interpolation submission environment and launch the interpolation workflow.

    Parameters
    ----------
    **kwargs : dict
        Optional command-line arguments that override default configuration values.
    """

    interpolation = SubmitInterpolation(**kwargs)
    interpolation.run()
