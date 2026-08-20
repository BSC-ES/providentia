# ICAP

Providentia's download mode supports downloading model output data provided by the **International Cooperative for Aerosol Prediction** (ICAP) from the **U.S. Navy / Naval Research Laboratory** (NRL) ICAP-MME data archive:

[https://nrlgodae1.nrlmry.navy.mil/cgi-bin/datalist.pl?dset=nrl_icap_mme&summary=Go](https://nrlgodae1.nrlmry.navy.mil/cgi-bin/datalist.pl?dset=nrl_icap_mme&summary=Go)

No account is required to download the files.

## Available CAMS datasets

Currently, only one ICAP dataset is available for download:

1. [**ICAP MME - The International Cooperative for Aerosol Prediction (ICAP) Multi-Model Ensemble (MME)**](#icap-mme---the-international-cooperative-for-aerosol-prediction-icap-multi-model-ensemble-mme)

(how-to-enable-icap-download)=
## How to enable ICAP download 

To download ICAP data, include one of the following model names along with the domain:

* `icap_ensemble-global` - [**ICAP MME - The International Cooperative for Aerosol Prediction (ICAP) Multi-Model Ensemble (MME)**](#icap-mme---the-international-cooperative-for-aerosol-prediction-icap-multi-model-ensemble-mme)

* `icap_ensemble_C4-global` - [**ICAP MME - The International Cooperative for Aerosol Prediction (ICAP) Multi-Model Ensemble (MME)**](#icap-mme---the-international-cooperative-for-aerosol-prediction-icap-multi-model-ensemble-mme)

Once you have selected the dataset and specified the model and domain, make sure to set:

```ini
dl_interpolated = False
```

or

Answer `n` to the prompt:

_"Model data was detected in the configuration file. Do you want to download the interpolated version? (Otherwise, the non-interpolated model data will be downloaded)"_

## ICAP MME - The International Cooperative for Aerosol Prediction (ICAP) Multi-Model Ensemble (MME) 

[Dataset Link](https://nrlgodae1.nrlmry.navy.mil/ftp/outgoing/nrl/ICAP-MME/)

Using this dataset you can download **Global Multi-Model Ensemble data**.

Only 6hourly global data is available in this dataset. You must set `resolution = 6hourly` and `domain = global` in your configuration file; otherwise, the download will not work.

### Mandatory fields for Multi-Model Ensemble/Consensus (MME/MMC)

```ini
model = icap_ensemble
domain = global
resolution = 6hourly
```

### Mandatory fields for 4-model consensus (C4)

```ini
model = icap_ensemble_C4
domain = global
resolution = 6hourly
```

### Available variables (Species)

The only available species is:

* `od550du`

You may also specify related species in the [mapping_species.yaml](interpolation-mapping-species) file (e.g. `od550aero`, `od500aerocoarse` or other mapped aliases). Regardless of the requested alias, the downloaded data will always correspond to `od550du`.

### Dataset temporal coverage

The dataset covers the period from November 2014 to the present.

### Product availability

The ICAP MME archive changes its naming convention over time.

- **From November 2014 up to and including September 2022**, only the **MME** total dust product is available, meaning that all the files have the following format:
  - `icap_yyyymmddhh_MME_totaldustaod550.nc`

- **From September 2022 onwards**, the archive introduces separate **C4** and **MMC** dust products, meaning that all the files have the following format:
  - `icap_yyyymmddhh_C4_dustaod550.nc`
  - `icap_yyyymmddhh_MMC_dustaod550.nc`

When `model = icap_ensemble` is specified, Providentia selects the product based on the requested date:

* Dates up to and including **August 2022**: download the **MME** total dust product.

* Dates from **September 2022** onwards: download the **MMC** dust product.

When `model = icap_ensemble_C4` is specified, Providentia downloads the **C4** dust product for dates from **September 2022** onwards.

## Example configuration files

#### Example 1: ICAP Multi-Model Ensemble/Consensus (MME/MMC)
```ini
[ICAP-MME_MMC] 
start_date = 20220701
end_date = 20230101
species = od550aero
resolution = 6hourly
model = icap_ensemble-global
dl_interpolated = False
```

#### Example 2: ICAP 4-model consensus (C4)
```ini
[ICAP-C4] 
start_date = 20220901
end_date = 20230101
species = od550aero
resolution = 6hourly
model = icap_ensemble_C4-global
dl_interpolated = False
```
