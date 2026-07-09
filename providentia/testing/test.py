import yaml
import os

files = [
    "/home/pserrano/push_pr/providentia/settings/basic_stats.yaml",
    "/home/pserrano/push_pr/providentia/settings/color_palettes.yaml",
    "/home/pserrano/push_pr/providentia/settings/exceedances.yaml",
    "/home/pserrano/push_pr/providentia/settings/fairmode.yaml",
    "/home/pserrano/push_pr/providentia/settings/interp_models.yaml",
    "/home/pserrano/push_pr/providentia/settings/model_bias_stats.yaml",
    "/home/pserrano/push_pr/providentia/settings/plot_characteristics.yaml",
    "/home/pserrano/push_pr/providentia/settings/remove_extreme_stations.yaml",
    "/home/pserrano/push_pr/providentia/settings/report_plots.yaml",
]

path = "/home/pserrano/push_pr/providentia/settings/temp"

for file in files:
    d = yaml.safe_load(open(file))
    final_path = os.path.join(path, os.path.basename(file))
    yaml.safe_dump(d, open(final_path, "w"))
