# BSC HPC machines

Providentia's download mode supports downloading **GHOST** and **non-GHOST** observational data, as well as model outputs, directly from BSC HPC systems. 

In order to do this type of download, a BSC HPC account is required.

## Supported HPC login nodes

Providentia can download GHOST and non-GHOST networks from the BSC HPC environment using available login nodes. By default, it connects to `storage5` (or `mn5` if `storage5` is unavailable).  

The login nodes used are defined in `settings/dl_hpc.yaml`:

- `transfer1.bsc.es` (storage5)
- `transfer2.bsc.es` (storage5)
- `transfer3.bsc.es` (storage5)  
- `transfer4.bsc.es` (storage5)  
- `glogin4.bsc.es` (mn5)  

The download mode will attempt to access the nodes in the order listed. If the node is unavailable, the system will automatically try the next node in the list, and so on.

## Download of networks from HPC machines

The network is essential to generate a report, as it provides access to the real observational data. At BSC, a large number of observational datasets are already stored on `gpfs`.

**How to enable:**  

- You must include at least one network in your configuration.  
- Set `dl_mode` to `both` or `obs` (or answer `both`/`obs` to the prompt _"Which type of data do you want to download? Observational, modelled or both?"_).  
- For GHOST network downloads, answer `y` to the prompt:  
  _"Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo)"_  
  or set `dl_ghost_source = bsc` in your configuration.
- **If the configuration file contains only network data and no models, network data will be downloaded automatically.**

**Download source paths:**

- Networks are **saved** under the `ghost_root` and `nonghost_root` paths specified in the `local` key of `settings/data_paths.yaml`.  
- Networks are **retrieved** from the `ghost_root` and `nonghost_root` paths specified in the `storage5` or `mn5` key of `settings/data_paths.yaml`.

If your network data to be **retrieved** is stored in a different directory, you can update the corresponding path in the `storage5` or `mn5` key to point to the correct location.

**Data format requirements:**

To be detected by Providentia, network data must follow these folder structures:

- GHOST networks: `{network}/{ghost_version}/{resolution}/{species}/{species_YYYYMM.nc}`

- Non-GHOST networks: `{network}/{resolution}/{species}/{species_YYYYMM.nc}`

## Download of interpolated models

An interpolated model consists of model output that have already been spatially interpolated to the observations. These models are ready to be used in modes such as report and dashboard. Most interpolated model data is stored on `gpfs`.

**How to enable:**  

- You must include at least one network and one model in your configuration.  
- Set `dl_mode` to `both` or `mod` (or answer `both`/`mod` to the prompt).  
- Answer `y` to the prompt:  
  _"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_  
  or set `dl_interpolated = True` in your configuration.

**Download source paths:**

- Models are **saved** under the `mod_root` path specified in the `local` key of `settings/data_paths.yaml`.  
- Models are **retrieved** from the `mod_root` path specified in the `storage5` or `mn5` key of `settings/data_paths.yaml`.

your model data to be **retrieved** in a different directory, you can update the corresponding path in the `storage5` or `mn5` key to point to the correct location.

**Data format requirements:**

To be detected by Providentia, interpolated model data must follow this folder structure:

`{ghost_version}/{model_id}-{domain}-{ensemble}/{resolution}/{species}/{network}/{species_YYYYMM.nc}`

## Download of non-interpolated models

### Local non-interpolated downloads

Non-interpolated model data refers to model outputs that are ready to be interpolated against a network using interpolation mode. Most non-interpolated datasets are stored in `esarchive`, although some may already exist in `gpfs`.

**How to enable:**  

- You must include at least one model in your configuration.  
- Set `dl_mode` to `both` or `mod` (or answer `both`/`mod` to the prompt).  
- Answer `n` to the prompt:  
  _"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_  
  or set `dl_interpolated = False` in your configuration.
- **If the configuration file contains only model data and no networks, non-interpolated model data will be downloaded automatically.**

**Download source paths:**

- Models are **saved** under the `mod_to_interp_root` path specified in the `local` key in `settings/data_paths.yaml`.  
- Firstly, models are **retrieved** from the paths specified in `settings/interp_models.yaml`. 
- If the model is not in the path, the system **retrieves** from the `mod_to_interp_root` path under the `storage5` or `mn5` key in `settings/data_paths.yaml`.  

If your data is stored in a different directory you can update paths directly in `settings/interp_models.yaml`. To learn how to define models, please see the [Defining models in interp_models.yaml](define-models) section in Interpolation.

You can also update the corresponding `storage5` or `mn5` path in `mod_to_interp_root` to point to the correct location.

**Data format requirements:**

To be detected by Providentia, interpolated model data must follow this folder structure:

- Standard non-interpolated models: `{model_id}/{domain}/{resolution}/{species}/{species_YYYYMM[DD].nc}`

- Ensemble-stats species non-interpolated models: `{model_id}/{domain}/{resolution}/ensemble-stats/{species}_{ensemble_stat}_an/{species_YYYYMM[DD]_{ensemble_stat}_an.nc}`

### HPC non-interpolated downloads
In some cases, model data exists in `esarchive`, but not all HPC machines have direct access to it. When interpolation needs to be performed on a machine without `esarchive` access, the model data must first be copied from `esarchive` to `gpfs`.

**How to enable:**  

- The download **must be performed from the `storage5` machine**.

**Download source paths:**

- Models are copied from `esarchive` to the `gpfs` `mod_to_interp_root` folder defined under the `storage5` key in `settings/data_paths.yaml`.  
- Only copies from paths specified in `settings/interp_models.yaml`. To learn how to define models, please see the [Defining models in interp_models.yaml](define-models) section in Interpolation.

## .env file

An `.env` file will appear in the Providentia root directory when using the download mode. It is designed to store specific user preferences.

   - **PRV_USER:** This setting specifies the username used to connect to the remote machines. It can be any valid username, e.g.: `bsc000000`.
   - **PRV_PWD:** This setting allows you to save the password needed for connecting to remote machines.  
  Note that the password is not required if you have configured a passwordless connection to the different servers.  
  Tutorial: [SSH Key Autologon](https://earth.bsc.es/wiki/doku.php?id=computing:sshkeyautologon&s%5B%5D=id_rsa.pub#ssh_key_autologon) _Only accessible for users with a BSC CAS account._

These values can be changed directly on the `.env` file and also be updated by Providentia during the next run.