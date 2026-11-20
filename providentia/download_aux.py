""" Auxiliary downloading functions """

import os
import time

from tqdm import tqdm
import yaml

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

# load the defined models paths, agrupations yaml and mapping species
data_paths = yaml.safe_load(open(join(PROVIDENTIA_ROOT, "settings/data_paths.yaml")))
interp_models = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "interp_models.yaml"))
)
mapping_species = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "mapping_species.yaml"))
)
dl_hpc = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "dl_hpc.yaml"))
)
temporal_resolution_map = yaml.safe_load(
    open(join(PROVIDENTIA_ROOT, "settings", "internal", "temporal_resolution_map.yaml"))
)


def find_model(download_instance, mod_id, domain, initial_check):
    # initialise warning message and model exists boolean
    msg = ""
    model_exists = False
    model_paths = []
    remote_dir = None

    # check model is in any of the interp_models.yaml lists
    for model_type, model_dict in interp_models.items():
        if mod_id in model_dict["models"]:
            model_paths = model_dict["paths"]
            break

    # if in the list, check if the paths work
    if model_paths:
        # stores all working paths
        mod_dir_functional_list = []

        for mod_dir in model_paths:
            # esarchive in transfer5 is located inside gpfs
            if "/esarchive/" == mod_dir[
                :11
            ] and download_instance.remote_hostname.startswith("transfer"):
                mod_dir = join("/gpfs/archive/bsc32/", mod_dir[1:])

            # add directory if it exists in the remote machine
            try:
                download_instance.sftp.stat(mod_dir)
                mod_dir_functional_list.append(mod_dir)
            except FileNotFoundError:
                pass

        # if none of the paths are in this current machine, break
        if not mod_dir_functional_list:
            msg += f"None of the paths specified in {join('settings', 'interp_models.yaml')} are available on the remote machine ({download_instance.remote_machine}). "

        else:
            # get first functional directory that has the model
            for mod_dir in mod_dir_functional_list:
                try:
                    remote_dir = join(mod_dir, mod_id, domain)
                    download_instance.sftp.stat(remote_dir)
                    model_exists = True
                    break
                except FileNotFoundError:
                    pass

            # if the model-domain combination is not possible, show the warning
            if model_exists is False:
                msg += f"There is no data available for the {mod_id} model with the {domain} domain in none of the paths specified in {join('settings', 'interp_models.yaml')} in the remote machine ({download_instance.remote_machine}). "

    # if no valid path model, search in the data_paths.yaml directory
    if model_exists is False:
        # get all possible models
        try:
            remote_dir = join(
                download_instance.mod_to_interp_remote_path, mod_id, domain
            )
            download_instance.sftp.stat(remote_dir)
            model_exists = True
        except FileNotFoundError:
            pass

    # if the model-domain combination is not possible, break
    if model_exists is False:
        # add to the message if model was not found in the gpfs remote directory
        msg += f"Cannot find the {mod_id} model with the {domain} domain in '{download_instance.mod_to_interp_remote_path}'."
        show_message(download_instance, msg, deactivate=initial_check)

    return model_exists, remote_dir


def find_available_resolutions(
    download_instance, remote_dir, mod_id, domain, initial_check
):
    resolution_not_found = False
    valid_resolutions = []

    # iterate through the resolutions
    for resolution in download_instance.model_resolution:
        # check existence of resolution directory
        try:
            download_instance.sftp.stat(join(remote_dir, resolution))
        except FileNotFoundError:

            if not resolution_not_found:
                # obtain the possible resolutions
                possible_resolutions = download_instance.sftp.listdir(remote_dir)

                # get the order of the possible resolutions from temporal_resolution_map
                resolution_order = temporal_resolution_map[resolution]

                # sort the possible resolutions in the temporal_resolution_map order
                sorted_resolutions = []
                for res in resolution_order:
                    if res in possible_resolutions:
                        sorted_resolutions.append(res)

                # get resolution input from user
                if sorted_resolutions:
                    while True:
                        user_resolution = input(
                        f"\nNo non-interpolated model data is available for {resolution} resolution. "
                        f"Available temporal resolutions for this model: {', '.join(sorted_resolutions)}. "
                        f"Please specify the desired resolution or press Enter to automatically select "
                        f"the first available option ({sorted_resolutions[0]}). "
                        f"Otherwise enter 'n' if you do not want to download this model at another resolution. "
                    ).lower()

                        if (user_resolution in possible_resolutions or user_resolution in ["", "n"]):
                            if user_resolution == "":
                                valid_resolutions.append(sorted_resolutions[0])
                            elif user_resolution != "n":
                                valid_resolutions.append(user_resolution)

                            # block input from the user happening again
                            resolution_not_found = True

                            break

        else:
            # add resolution to the return list
            valid_resolutions.append(resolution)

    return valid_resolutions


def find_available_species(
    download_instance,
    valid_resolutions,
    remote_dir,
    ensemble,
    mod_id,
    domain,
    initial_check,
):
    valid_resolutions_species_dir = []

    for resolution in valid_resolutions:
        species_exists = False

        # iterate through the species
        for original_species in download_instance.species:
            species_to_process = [original_species]

            # iterate through the original species and if not found through the mapped species
            if original_species in mapping_species:
                species_to_process = (
                    species_to_process + mapping_species[original_species]
                )

            for species in species_to_process:
                # check existance of species directory
                try:
                    # ensemble member
                    if ensemble.isdigit() or ensemble == "allmembers":
                        resolution_species_dir = join(remote_dir, resolution, species)
                    # ensemble statistic
                    else:
                        resolution_species_dir = join(
                            remote_dir,
                            resolution,
                            "ensemble-stats",
                            species + "_" + ensemble,
                        )

                    download_instance.sftp.stat(resolution_species_dir)

                    # store the resolution species pair
                    valid_resolutions_species_dir.append(
                        [
                            resolution_species_dir,
                            mod_id,
                            domain,
                            resolution,
                            species,
                            ensemble,
                        ]
                    )

                    # do not show warning
                    species_exists = True

                    # get the first species found
                    break

                # check existance of mapped species directory
                except FileNotFoundError:
                    pass

        # if no species were found, then show the message
        if species_exists is False:
            msg = (
                f"There is no data available in {download_instance.remote_machine} for the {mod_id} model with the "
                f"{domain} domain for {species} species at {resolution} resolution."
            )
            show_message(download_instance, msg, deactivate=initial_check)
            continue

    return valid_resolutions_species_dir


def find_available_nc_files(
    download_instance, valid_resolutions_species_dir, initial_check
):
    # initialise list with all the nc files to be downloaded
    path_files_dict = {}

    # get all the nc files in the date range
    for (
        resolution_species_dir,
        mod_id,
        domain,
        resolution,
        species,
        ensemble,
    ) in valid_resolutions_species_dir:
        # get nc files
        nc_files = download_instance.sftp.listdir(resolution_species_dir)

        if nc_files:
            # ensemble member
            if ensemble.isdigit() or ensemble == "allmembers":
                # identify format of the directory
                # the format is a tuple of how many '-' and how many '_' are there, e.g.: (0,1)
                # the directory format is choosen by popularity
                formats_list = [(file.count("-"), file.count("_")) for file in nc_files]
                number_of_formats_dict = {
                    format: formats_list.count(format) for format in set(formats_list)
                }
                format = max(number_of_formats_dict, key=number_of_formats_dict.get)

                # filter and get only the files that follow the format choosen
                nc_files = list(
                    filter(
                        lambda x: (x.count("-"), x.count("_")) == format
                        and x.endswith(".nc"),
                        nc_files,
                    )
                )

                # example: od550du_2019040212.nc (0,1)
                if format == (0, 1):
                    # if no ensemble in the name only allmembers and 000 are valid
                    if ensemble == "000" or ensemble == "allmembers":
                        nc_files = list(
                            filter(lambda x: x.split("_")[0] == species, nc_files)
                        )

                # example: od550du-000_2021020812.nc (1,1)
                elif format == (1, 1):
                    # filter by ensemble in case that ensemble is not allmembers
                    if ensemble != "allmembers":
                        nc_files = list(
                            filter(
                                lambda x: x.split("_")[0] == species + "-" + ensemble,
                                nc_files,
                            )
                        )

                # unknown format
                else:
                    msg = f"It is not possible to download this nc file type yet. Please, contact the developers. Files to download: {nc_files}"
                    show_message(download_instance, msg, deactivate=initial_check)
                    continue

                local_dir = join(
                    download_instance.mod_to_interp_root,
                    mod_id,
                    domain,
                    resolution,
                    species,
                )

            # ensemble statistic
            else:
                # example: sconco3_2025120300_av_an.nc
                nc_files = list(
                    filter(
                        lambda x: x.split("_")[0] == species
                        and "_".join(x[:-3].split("_")[2:]) == ensemble,
                        nc_files,
                    )
                )

                local_dir = join(
                    download_instance.mod_to_interp_root,
                    mod_id,
                    domain,
                    resolution,
                    "ensemble-stats",
                    species + "_" + ensemble,
                )

            # store the available files for each path
            if nc_files:
                path_files_dict[local_dir] = {
                    "nc_files": nc_files,
                    "remote_dir": resolution_species_dir,
                    "mod_id": mod_id,
                    "species": species,
                    "resolution": resolution,
                }

            else:
                msg = f"There is no data available in {download_instance.remote_machine} for the {mod_id} model with the {domain} domain with the {ensemble} ensemble."
                show_message(download_instance, msg, deactivate=initial_check)

    return path_files_dict


def build_nc_file_paths_in_range(download_instance, path_files_dict, initial_check):
    initial_check_nc_files = {}

    for dir, model_dict in path_files_dict.items():
        nc_files = model_dict["nc_files"]

        valid_nc_files = download_instance.get_valid_nc_files_in_date_range(nc_files)
        valid_nc_files.sort()

        # warning if model + species + resolution + network + date range combination gets no matching results
        if not valid_nc_files:
            msg = (
                f"There is no data available in {download_instance.remote_machine} from {download_instance.start_date} to {download_instance.end_date} "
                f"for {model_dict['mod_id']} model {model_dict['species']} species at {model_dict['resolution']} resolution."
            )
            show_message(download_instance, msg, deactivate=initial_check)
            continue

        for nc_file in valid_nc_files:
            if dir not in initial_check_nc_files:
                initial_check_nc_files[dir] = {
                    "remote_dir": model_dict["remote_dir"],
                    "nc_files": [nc_file],
                }
            else:
                initial_check_nc_files[dir]["nc_files"].append(nc_file)

    return initial_check_nc_files


# TODO change name
def download_non_interpolated_sftp(download_instance, files_to_download_info):
    for local_dir, files_to_download_dict in files_to_download_info.items():
        remote_dir = files_to_download_dict["remote_dir"]

        # print source and destination
        download_instance.logger.info(
            f"\n  - {remote_dir}, source: {local_dir} ({download_instance.remote_machine})"
        )

        # create directories if they don't exist
        os.makedirs(local_dir, exist_ok=True)

        # download each individual nc file using sftp protocol
        for nc_file in tqdm(
            files_to_download_dict["nc_files"],
            bar_format="{l_bar}{bar}|{n_fmt}/{total_fmt}",
            desc=f"    Downloading files ({len(files_to_download_dict['nc_files'])})",
        ):
            local_path = join(local_dir, nc_file)
            remote_path = join(remote_dir, nc_file)

            # get last downloaded file in case there was a keyboard interrupt
            download_instance.latest_nc_file_path = local_path

            # initialize the timeout and get the file
            download_instance.ncfile_dl_start_time = time.time()
            download_instance.sftp.get(
                remote_path, local_path, callback=download_instance.check_time
            )

            # change the last downloaded file
            download_instance.latest_nc_file_path = "/path/to/file"
