import yaml
from variable_mapping import variable_mapping

result = {
    value['preferred_term'].replace('"', ''): {'var': key[2], 'units': key[0]}
    for key, value in variable_mapping.items()
}

with open('variable_mapping.yaml', 'w') as file:
    yaml.dump(result, file, default_flow_style=False)