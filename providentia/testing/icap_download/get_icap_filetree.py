import requests
from bs4 import BeautifulSoup

dict = {}

og_url = "https://nrlgodae1.nrlmry.navy.mil/ftp/outgoing/nrl/ICAP-MME/"

for i in range(2014, 2027):
    for j in range(1, 13):
        if i == 2014 and j not in [11, 12]:
            continue

        last_url = f"{i}{j:02d}"

        url = f"{og_url}{i}/{last_url}"

        print(last_url)

        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, "html.parser")

        files = [
            a["href"]
            for a in soup.find_all("a", href=True)
            if a["href"].startswith("icap_")
        ]

        dict[last_url] = files

        print(dict)

        # with open("urls.yaml", "w") as f:
        #     yaml.dump(dict, f, sort_keys=False)

        

