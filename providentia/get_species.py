import yaml

ghost_cams_variables = yaml.safe_load(
    open(
        "/home/pserrano/format_bug/providentia/settings/internal/cams/ghost_cams_variables.yaml"
    )
)
cams_dataset = yaml.safe_load(
    open(
        "/home/pserrano/format_bug/providentia/settings/internal/cams/cams_dataset.yaml"
    )
)

ghost_list = []

for ecmwf_var in cams_dataset["era5_reanalysis"]["global"]["variable"]:
    for ghost_var, ecmwf_var2 in ghost_cams_variables.items():
        if ecmwf_var == ecmwf_var2:
            ghost_list.append(ghost_var)
            break

for ghost_var, ecmwf_var2 in ghost_cams_variables.items():
    if type(ecmwf_var2) is list:
        for var in ecmwf_var2:
            if var in ghost_list:
                ghost_list.append(ghost_var)
                break

for i in sorted(ghost_list):
    print(f"{i} : {ghost_cams_variables[i]}")
