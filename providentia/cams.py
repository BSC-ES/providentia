import sys
import os

import cdsapi
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from netCDF4 import Dataset
import numpy as np
import re
import requests
import shutil
import yaml
import zipfile 

from providentia.auxiliar import CURRENT_PATH, join
from .warnings_prv import show_message

PROVIDENTIA_ROOT = os.path.dirname(CURRENT_PATH)

cams_options = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_dataset.yaml')))
cams_variables_level = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_variables_level.yaml')))
cams_formatting = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_formatting.yaml')))
ghost_cams_variables = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'ghost_cams_variables.yaml')))
cams_stream = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings', 'internal', 'cams', 'cams_stream.yaml')))

class Cams:

    def __init__(self, download_instance):
        self.download_instance = download_instance

    def print_request(self, dataset, request):
        
        self.download_instance.logger.info(f"dataset = '{dataset}'")
        self.download_instance.logger.info('request = {')
        for k,v in request.items():
            if type(v) == str:
                v = f"'{v}'"
            self.download_instance.logger.info(f"'{k}' : {v},")
        self.download_instance.logger.info('}\n')

    def fetch_cams_dates(self, url, cams_dict):
        
        # send HTTP GET request and get the text
        response = requests.get(url)
        r = response.text
        
        # do the webscrapping depending if there is whole dates or only month
        if cams_dict['month_names'] is False:
            # get the minstart and maxend dictionary
            minstart_dict = re.findall(r'"minStart":".*?"', r, re.DOTALL)
            maxend_dict = re.findall(r'"maxEnd":".*?"', r, re.DOTALL)
            
            # get the value from the dictionary
            minstart = datetime.strptime(minstart_dict[0].split('"')[-2], '%Y-%m-%d')
            maxend = datetime.strptime(maxend_dict[0].split('"')[-2], '%Y-%m-%d')
        else:
            # get the interval dictionary
            match = re.search(r'"interval":\[\["(.*?)","(.*?)"\]\]', r)

            # get the date value
            minstart = datetime.strptime(match.group(1), "%Y-%m-%dT%H:%M:%SZ")
            maxend = datetime.strptime(match.group(2), "%Y-%m-%dT%H:%M:%SZ")

        return minstart, maxend
        
    def control_dates(self, url, cams_dict):
        # check whether end date is bigger than start date, if not, return
        if int(self.download_instance.start_date) >= int(self.download_instance.end_date):
            msg = f'Start date ({self.download_instance.start_date}) exceeds end date ({self.download_instance.end_date}).'
            show_message(self.download_instance, msg, print=True)
            return None, None

        # get minimum and maximum possible dates
        min_start_date, max_end_date = self.fetch_cams_dates(url, cams_dict)
        
        # convert the selected dates to datetetime
        cams_start_date = datetime.strptime(self.download_instance.start_date, "%Y%m%d")
        cams_end_date = datetime.strptime(self.download_instance.end_date, "%Y%m%d") - timedelta(days=1)

        # if the minimum date is over the end date
        if min_start_date > cams_end_date or max_end_date < cams_start_date:
            msg = f"The selected dates are unavailable. Please choose dates between {min_start_date.strftime('%Y-%m-%d')} and {max_end_date.strftime('%Y-%m-%d')}."
            show_message(self.download_instance, msg)
            return None, None

        # check if the start date is within limits
        if min_start_date > cams_start_date:
            cams_start_date = min_start_date
            
        # check if the end date is within limits
        if max_end_date < cams_end_date:
            cams_end_date = max_end_date

        return cams_start_date, cams_end_date

    # TODO check that message in stream does not get repeated in the loop
    # TODO check if it is plaussible to keep this code as a separate function since there is a continue
    def create_request(self, cams_species, cams_dict, current_cams_date, next_cams_date, level, stream, url, exp_id):
        # create the request
        request = {"variable" : cams_species}

        # add leadtime hour to the request if the dataset has it
        if 'leadtime_hour' in cams_dict:
            request["leadtime_hour"] = cams_dict['leadtime_hour'][level]

        # add type to the request if the dataset has it
        if 'type' in cams_dict:
            request["type"] = cams_dict['type']

        # if it's forecast one file per day, analysis one file per month
        if cams_dict['month_names'] is False:
            request["date"] = f"{current_cams_date.strftime('%Y-%m-%d')}/{next_cams_date.strftime('%Y-%m-%d')}"
        else:
            request["year"] = str(current_cams_date.year)
            request["month"] = str(current_cams_date.strftime('%m'))

        # add interim and or validated stream to the request
        if cams_dict['stream'] is True:
            # check whether the stream is available for the year
            if stream not in cams_stream[request["year"]]:
                msg = f"The current stream '{stream}' is not available for the current date: {request['month']}-{request['year']}. Continuing...."          
                show_message(self.download_instance, msg)
                # continue
                
            request["type"] = stream

        # add type to the request if the dataset has it
        if 'type' in cams_dict:
            request["type"] = cams_dict['type']

        # add the experiment if models are available in the dataset
        if 'models' in cams_dict:
            request["model"] = exp_id

        # get the level and apply it if the species is multi level
        level_variable = 'level' if 'level' in cams_dict else 'model_level'
        if cams_species in cams_variables_level[url]['multi']:
            # choose the maximum level if there was a level increase at some point
            if 'level_boundary' in cams_dict:
                boundary_date = datetime.combine(cams_dict['level_boundary'], datetime.min.time())
                if current_cams_date <= boundary_date:
                    request[level_variable] = cams_dict[level_variable]['before_increase']
                else: 
                    request[level_variable] = cams_dict[level_variable]['after_increase']
            # without a level increase, just get the level
            else:
                request[level_variable] = cams_dict[level_variable]
        
        # add time to the request
        if 'time' in cams_dict:
            request['time'] = cams_dict['time']

        # add data_format to the request
        if 'data_format' in cams_dict:
            request['data_format'] = cams_dict['data_format']

        return request

    def create_cdsapirc(self, cdsapirc_path):      
        # ask the user whether they want to create the file in the home directory
        create_file = input(f"'.cdsapirc' file not found. Creating it at {cdsapirc_path}. Do you agree? ([y]/n) ").lower()
        while create_file not in ['','y','n']:
            create_file = input(f"'.cdsapirc' file not found. Creating it at {cdsapirc_path}. Do you agree? ([y]/n) ").lower()

        # create file if user agreed with it
        if create_file in ['', 'y']: 
            # ask the user for the personal access token
            personal_access_token = input("Enter your personal access token, which you can find at https://cds.climate.copernicus.eu/how-to-api: ")
            # create the .cdsapirc file with the user's acces token
            with open(cdsapirc_path, "w") as f:
                f.write("url: https://ads.atmosphere.copernicus.eu/api\n")
                f.write(f"key: {personal_access_token}\n")
        else:
            self.download_instance.logger.error("Error: Cannot proceed without '.cdsapirc'. CAMS experiment data download requires this file.")
            sys.exit(1)

    def get_experiment(self, cams_dict, u_count, config_expid, dataset, ensemble_options):
        # determine experiment ID or stream
        exp_id, stream = None, None
        
        if u_count == 1: # e.g. cams_forecast
            if 'models' in cams_dict:
                msg = f"The experiment '{config_expid}' is missing the model. Please add one (e.g., '{config_expid}_ensemble')."
                show_message(self.download_instance, msg)
                return None, None, True

        elif u_count == 2:
            # extract last element
            last_element = config_expid.rsplit('_', 1)[1]

            if cams_dict['stream'] is True and 'models' in cams_dict: # e.g. cams_reanalysis_ensemble-regional
                msg = f"The '{dataset}' dataset needs a model and a stream. E.g. 'cams_reanalysis_ensemble_interim')."    
                show_message(self.download_instance, msg)
                return None, None, True
            elif 'models' in cams_dict: # e.g. cams_analysis_ensemble
                exp_id = last_element
                # make sure the experiment is available in the dataset
                if exp_id not in cams_dict["models"]:
                    msg = f"Cannot find the {exp_id} model in the '{dataset}' dataset."    
                    show_message(self.download_instance, msg)
                    return None, None, True
            else:
                # if there are three elements and they
                msg = f"The '{dataset}' dataset does not admit models or streams."    
                show_message(self.download_instance, msg)
                return None, None, True
                
        elif u_count == 3:

            if not (cams_dict['stream'] is True and 'models' in cams_dict): 
                # if there are three elements and they
                msg = f"The '{dataset}' dataset does not admit models and streams, change the experiment in the configuration file."    
                show_message(self.download_instance, msg)
                return None, None, True
            
            # extract the last two elements
            _, exp_id, stream = config_expid.rsplit('_', 2)

            # make sure the experiment is available in the dataset
            if exp_id not in cams_dict["models"]:
                msg = f"Cannot find the {exp_id} model in the '{dataset}' dataset."    
                show_message(self.download_instance, msg)
                return None, None, True
            
            # make sure the stream is valid
            if stream not in ['validated','interim']:
                msg = f"'{stream}' is not a valid stream. Availabe streams: validated, interim."    
                show_message(self.download_instance, msg)
                return None, None, True
            
            # add reanalysis sufix to the stream
            stream += '_reanalysis'

        # get error if it does not get any of the conditions
        else:
            msg = f"The '{config_expid}' format is not valid."    
            show_message(self.download_instance, msg)
            return None, None, True

        # only ensemble options allmembers and 000 are valid
        if ensemble_options not in ['000', 'allmembers']:
            msg = (
            f"The current ensemble option '{ensemble_options}' is not valid for the CAMS '{dataset}' dataset."
            f"It must be '000' or 'allmembers'.")            
            show_message(self.download_instance, msg)
            return None, None, True
        
        return exp_id, stream, False

    def extract_date(self, input_file, prefix, domain):
        
        if prefix in ['cams_analysis','cams_forecast'] and domain == 'regional':
            time = input_file['time'].long_name.split()[-1]
            time = datetime.strptime(time, '%Y%m%d')
        elif prefix == 'cams_reanalysis' and domain == 'regional':
            time = input_file['time'].units.split()[2]
            time = datetime.strptime(time, '%Y-%m-%d')
        elif prefix == 'cams_forecast' and domain == 'global':
            time = input_file['forecast_reference_time'][0]
            time = datetime.fromtimestamp(int(time))
        elif prefix == 'cams_reanalysis' and domain == 'global':
            time = input_file['valid_time'][0]
            time = datetime.fromtimestamp(int(time))

        return time.strftime('%Y-%m-%d')

    def format_data(self, input_filepath, output_filepath, species, prefix, domain, resolution, final_path, cams_species, url): 

        self.download_instance.logger.info(f"Formatting {final_path}\n") 

        # get file formatting 
        cams_providentia_map = cams_formatting[prefix][domain]

        # open original netcdf file      
        input_file = Dataset(input_filepath, 'r', format="NETCDF4") 

        # extract date 
        date_str = self.extract_date(input_file, prefix, domain)

        # create new netcdf file
        output_file = Dataset(output_filepath, 'w', format="NETCDF4") 
        output_file.set_auto_mask(True)	

        for input_dim_name, output_dim_name in cams_providentia_map.items():
            # skip single species
            if output_dim_name == "level" and 'single' in cams_variables_level[url] and cams_species in cams_variables_level[url]['single']:
                pass
            # get dimension
            dim = input_file.dimensions[input_dim_name]
            # create the dimension with the new name 
            output_file.createDimension(output_dim_name, len(dim))

        # get species name in the input file
        input_species = list(set(input_file.variables) - set(input_file.dimensions))[0] # TODO take into account that cams_forecast_global has two elements

        # copy variables
        for input_var_name in input_file.variables:
            # get the output var name
            if input_var_name in cams_providentia_map:
                output_var_name = cams_providentia_map[input_var_name]
            elif input_var_name == input_species:
                output_var_name = species
            else:
                continue
            
            # get the variable
            input_var = input_file[input_var_name]
            
            # change the name of the dimensions into the providentia name
            output_var_dims = [cams_providentia_map.get(name, name) for name in input_var.dimensions if name in cams_providentia_map]
            
            # create the variable
            output_var = output_file.createVariable(output_var_name, input_var.datatype, output_var_dims)

            # add calendar and units attributes to the time variable
            if output_var_name == "time":
                output_var.setncattr('calendar', cams_options[prefix][domain]['calendar'])
                time_units = f"hours since {date_str}"
                output_var.setncattr('units', time_units)

            # add to level the priority of the level
            elif output_var_name == "level":
                output_var.positive = 'up'

            # add coordinates, grid_mapping and units to the species variable
            elif input_var_name == input_species:               
                output_var.setncattr('coordinates', 'latitude longitude')
                output_var.setncattr('grid_mapping', 'crs')
                output_var.setncattr('units', input_var.units)
            
            # get the data from 
            if output_var_name == "time":  
                data = np.arange(len(input_var[:]))
                if resolution == "3hourly":
                    data *= 3
            else:
                data = input_var[:]

            # add the data to the variable
            output_var[:] = data
        
        # add grid_mapping
        output_file[species].setncattr('grid_mapping', 'crs')
        
        # add coordinates
        output_file[species].setncattr('coordinates', 'latitude longitude')

        # add crs
        crs_var = output_file.createVariable('crs', 'u1')  
        crs_var.setncatts({
            'grid_mapping_name': 'latitude_longitude',
            'semi_major_axis': 6371000.0,
            'inverse_flattening': 0.0
        })
        
        # close the original and new netcdf files
        output_file.close()
        input_file.close()           

    def download_cams_experiment(self, experiment): 
        # print current_experiment
        self.download_instance.logger.info('\n'+'-'*40)
        self.download_instance.logger.info(f"\nDownloading {experiment} experiment data from the Atmosphere Data Store...")

        # create cdsapirc file in case it was not created
        cdsapirc_path = join(os.getenv("HOME"),'.cdsapirc')

        # create csapirc file necessary for the download
        if not os.path.isfile(cdsapirc_path):
            self.create_cdsapirc(cdsapirc_path)

        # get experiment id and the domain
        config_expid, domain, ensemble_options = experiment.split("-")
        
        # count underscores to determine format
        u_count = config_expid.count('_')

        # get the prefix
        prefix = config_expid.rsplit('_', u_count-1)[0]

        # check if the domain is the correct one for the dataset
        if domain not in cams_options[prefix]:
            possible_domains = "', '".join(cams_options[prefix])
            msg = (
            f"The current domain '{domain}' is not valid for the CAMS dataset."
            f"It must be '{possible_domains}'.")            
            show_message(self.download_instance, msg)
            return

        # get the dictionary with the dataset characteristics
        cams_dict = cams_options[prefix][domain]

        # get CAMS dataset and url
        dataset = cams_dict["dataset"]
        url = cams_dict['url']
        
        # make the necessary checks to the experiment
        exp_id, stream, invalid_experiment = self.get_experiment(cams_dict, u_count, config_expid, dataset, ensemble_options)
    
        # stop download if the experiment format is not correct
        if invalid_experiment:
            return
    
        # make the necessary checks to the dates
        cams_start_date, cams_end_date = self.control_dates(url, cams_dict)

        # stop download if the dates are not correct
        if cams_start_date is None and cams_end_date is None:
            return

        # iterate through the species
        for species in self.download_instance.species: 
            # check if species is in the ghost_cams_variables file
            if species not in ghost_cams_variables:
                msg = f"The species '{species}' is not available in CAMS."
                show_message(self.download_instance, msg)
                continue
        
            # get the species in the cams vocabulary
            cams_species = ghost_cams_variables[species]

            # check if the mapped species are available in the dataset
            if cams_species not in cams_dict['variable']:
                msg = f"Mapped species '{cams_species}' for input species '{species}' is not available in the CAMS '{dataset}' dataset."          
                show_message(self.download_instance, msg)
                continue
            
            # get the species' level
            level = 'multi' if cams_species in cams_variables_level[url]['multi'] else 'single'

            # iterate throught the resolutions
            for resolution in self.download_instance.resolution:
                # get the resolution for the cams dataset
                correct_resolution = cams_dict["resolution"] if type(cams_dict["resolution"]) == str else cams_dict["resolution"][level]
                
                # check if the resolution is the correct one for the dataset
                if resolution != correct_resolution:
                    msg = (
                    f"The current resolution '{resolution}' is not valid. It must be '{correct_resolution}'.")            
                    show_message(self.download_instance, msg)
                    continue

                # get directory structure
                dir_tail = join(config_expid, domain, resolution, species)

                # get temporal and final dir
                temp_dir = join(self.download_instance.exp_to_interp_root,'.temp', dir_tail)
                final_dir = join(self.download_instance.exp_to_interp_root, dir_tail)

                # create dir
                os.makedirs(final_dir, exist_ok=True)

                # create client
                self.download_instance.logger.info('')
                client = cdsapi.Client(retry_max=1, quiet=True)

                # iterate through the dates
                current_cams_date = cams_start_date
                while current_cams_date <= cams_end_date:

                    # add one day or one month depending if it is forecast or analysis
                    if cams_dict['forecast'] is False:
                        next_cams_date = current_cams_date.replace(day=1) + relativedelta(months=1) - timedelta(days=1)
                        next_cams_date = cams_end_date if next_cams_date > cams_end_date else next_cams_date
                    else:
                        next_cams_date =  current_cams_date

                    # create dictionary to do the request
                    request = self.create_request(cams_species, cams_dict, current_cams_date, next_cams_date, level, stream, url, exp_id)

                    # get filename depending whether it is a download for the whole month or just a day
                    date_format = '%Y%m' if cams_dict['forecast'] is False else '%Y%m%d'
                    
                    # get final path
                    file_name = f"{species}_{current_cams_date.strftime(date_format)}.nc"
                    final_path = join(final_dir, file_name)
                    
                    # create temporal dir to store the middle zip file with its directories
                    os.makedirs(temp_dir, exist_ok=True)

                    # get temporal path
                    temp_path = join(temp_dir, 'zip_file')

                    # print the request
                    self.print_request(cams_dict['dataset'], request)

                    # get last downloaded file in case there was a keyboard interrupt
                    self.download_instance.latest_nc_file_path = final_path

                    # make the request
                    try:
                        self.download_instance.logger.info(f"Downloading {final_path}") # TODO change message
                        client.retrieve(dataset, request, target=temp_path)
                    except requests.exceptions.HTTPError as err:
                        # invalid credential on .cdsapirc
                        if err.response.status_code == 401: 
                            self.download_instance.logger.info(
                                "\nBad request (401): Client Error. Invalid credentials in the .cdsapirc file. "
                                "Removing authentication file...\n"
                                "Please run the program again so Providentia can recreate the file automatically, "
                                "or manually create a new .cdsapirc file by following the instructions at: "
                                "https://cds.climate.copernicus.eu/how-to-api"
                            )
                            os.remove(cdsapirc_path)
                            return
                        # bad request
                        if err.response.status_code == 400: 
                            self.download_instance.logger.info("\nBad request (400): The server could not understand the request.")
                            self.download_instance.logger.info(f"Details: {err}")
                        # connection error
                        elif err.response.status_code == 500: 
                            self.download_instance.logger.info("\nServer error (500): The server encountered an error while processing the request.")
                            self.download_instance.logger.info(f"Details: {err}")
                            self.download_instance.logger.info("Please try again later.")
                            return
                        else:
                            self.download_instance.logger.info(f"\nUnexpected error ({err.response.status_code}):")
                            self.download_instance.logger.info(f"Details: {err}")
                        # make the next download
                        continue

                    # extract file 
                    with zipfile.ZipFile(temp_path, 'r') as zip_ref:
                        zip_file_name = zip_ref.namelist()[0]
                        zip_ref.extractall(temp_dir)

                    # format the cams files and move them to the corresponding folder
                    self.format_data(join(temp_dir,zip_file_name), final_path, species, prefix, 
                                     domain, resolution, final_path, cams_species, url)

                    # add one day to the date
                    current_cams_date = next_cams_date + timedelta(days=1)    

                    # remove the temp directory tail
                    shutil.rmtree(temp_dir)
