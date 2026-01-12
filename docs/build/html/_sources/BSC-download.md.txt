# BSC HPC Machines

Providentia's download mode supports downloading **GHOST** and **non-GHOST** observational data, as well as model outputs, directly from BSC HPC systems. 

In order to do this type of download, a BSC HPC account is required.

## Supported HPC Login Nodes

Providentia can download GHOST and non-GHOST networks from the BSC HPC environment using available login nodes. By default, it connects to `storage5` (or `mn5` if `storage5` is unavailable).  

The login nodes used are defined in `settings/dl_hpc.yaml`:

- `transfer1.bsc.es` (primary)  
- `transfer2.bsc.es`  
- `transfer3.bsc.es`  
- `transfer4.bsc.es`  
- `glogin4.bsc.es`  

Downloads will attempt the nodes in the order listed. If the primary node is unavailable, the system automatically tries the next one in the list.
## Download of Networks from HPC Machines

- Saves networks in the `ghost_root` and `nonghost_root` paths specified in the `local` key of `settings/data_paths.yaml`.  
- Retrieves networks from the `ghost_root` and `nonghost_root` paths specified in the `storage5` or `mn5` key of `settings/data_paths.yaml`.

**How to enable:**  
- You must include at least one network in your configuration.  
- Set `dl_mode` to `both` or `obs` (or answer `both`/`obs` to the prompt _"Which type of data do you want to download? Observational, modelled or both?"_).  
- For GHOST network downloads, answer `y` to the prompt:  
  _"Do you want to download observational data from the BSC remote machine? (Otherwise, GHOST observational data will be retrieved from Zenodo)"_  
  or set `dl_ghost_source = bsc` in your configuration.

## Download of Interpolated Models

- Saves models in the `mod_root` path specified in the `local` key of `settings/data_paths.yaml`.  
- Retrieves models from the `mod_root` path specified in the `storage5` or `mn5` key of `settings/data_paths.yaml`.

**How to enable:**  
- You must include at least one network **and** one model in your configuration.  
- Set `dl_mode` to `both` or `mod` (or answer `both`/`mod` to the prompt).  
- Answer `y` to the prompt:  
  _"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_  
  or set `dl_interpolated = True` in your configuration.

## Download of Non-Interpolated Models

### Local Non-Interpolated Downloads

- Saves models in `mod_to_interp_root` paths defined in `settings/data_paths.yaml`.  
- Retrieves models from paths specified in `settings/interp_models.yaml`. To learn how to define models, please see the [Defining models in interp_models.yaml](define-models) section in Interpolation.
- If the model is not in the path, the system looks in the `mod_to_interp_root` paths under the `storage5` or `mn5` key in `settings/data_paths.yaml`.  

**How to enable:**  
- You must include at least one model in your configuration.  
- Set `dl_mode` to `both` or `mod` (or answer `both`/`mod` to the prompt).  
- Answer `n` to the prompt:  
  _"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_  
  or set `dl_interpolated = False` in your configuration.

### HPC Non-Interpolated Downloads

- This is used when experiments are available in `esarchive` and cannot be used directly for interpolation on machines without `esarchive` access.  
- Copies non-interpolated models from `esarchive` to the `gpfs` folder defined in the `storage5` key `mod_to_interp_root` in `settings/data_paths.yaml`.  
- Only copies from paths specified in `settings/interp_models.yaml`. To learn how to define models, please see the [Defining models in interp_models.yaml](define-models) section in Interpolation.

**How to enable:**  
- The download **must be performed from the `storage5` machine**.

## .env file

An `.env` file will appear in the Providentia root directory when using the download mode. It is designed to store specific user preferences.

   - **PRV_USER:** This setting specifies the username used to connect to the remote machines. It can be any valid username, e.g.: `bsc000000`.
   - **PRV_PWD:** This setting allows you to save the password needed for connecting to remote machines.  
  Note that the password is not required if you have configured a passwordless connection to the different servers.  
  Tutorial: [SSH Key Autologon](https://earth.bsc.es/wiki/doku.php?id=computing:sshkeyautologon&s%5B%5D=id_rsa.pub#ssh_key_autologon) _Only accessible for users with a BSC CAS account._

These values can be changed directly on the `.env` file and also be updated by Providentia during the next run.