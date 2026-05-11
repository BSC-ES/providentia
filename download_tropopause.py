# import requests
# import re
# from urllib.parse import urljoin

# base_url = "https://datapub.fz-juelich.de/slcs/tropopause/data/v2/era5/2000/"

# # Get directory HTML
# response = requests.get(base_url)
# response.raise_for_status()


# print(f"Found {len(files)} files")

# # Download
# for f in files:
#     file_url = urljoin(base_url, f)
#     print("Downloading:", file_url)

#     r = requests.get(file_url)
#     with open(f, "wb") as out:
#         out.write(r.content)