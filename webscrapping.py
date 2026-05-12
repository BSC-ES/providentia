import requests
from bs4 import BeautifulSoup
import json

url = "https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=overview"
page = requests.get(url)

soup = BeautifulSoup(page.text, "html.parser")

# գտ the Next.js JSON
script = soup.find("script", {"id": "__NEXT_DATA__"})

data = json.loads(script.string)

def find_all_matching(obj, results=None, path=""):
    if results is None:
        results = []

    if isinstance(obj, dict):
        # condición: tiene las 3 claves
        if all(k in obj for k in ["name", "description", "units"]):
            results.append((path, obj))

        # seguir buscando (sin return!)
        for k, v in obj.items():
            new_path = f"{path}.{k}" if path else k
            find_all_matching(v, results, new_path)

    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            find_all_matching(item, results, f"{path}[{i}]")

    return results

vars = ["10m_u-component_of_wind",
"10m_v-component_of_wind",
"2m_dewpoint_temperature",
"2m_temperature",
"boundary_layer_height",
"cloud_base_height",
"high_cloud_cover",
"instantaneous_10m_wind_gust",
"mean_sea_level_pressure",
"sea_surface_temperature",
"snow_depth",
"snowfall",
"surface_pressure",
"total_cloud_cover",
"total_precipitation"]

for i, var in enumerate(vars):
    if '_' in var:
        vars[i] = var.replace("_", " ")

possible_units = {}
for i, d in find_all_matching(data):
    symbols = ["<sup>", "</sup>"]
    unit = d['units']
    for s in symbols:
        if s in unit:
            unit = unit.replace(s, "")
    
    if d['name'].lower() in vars:
        if unit in possible_units:
            possible_units[unit].append(d['name'].lower())
        else:
            possible_units[unit] = [d['name'].lower()]

for u in possible_units:
    print(u)
    for var in possible_units[u]:
        print("   -",var)