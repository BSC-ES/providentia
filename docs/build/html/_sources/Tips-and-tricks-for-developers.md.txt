# Tips and tricks for developers

This page can be used by developers to better understand certain parts of the code and perform certain actions:

- [Migrate repository from Gitlab to Github](#migrate-repository-from-gitlab-to-github)
- [Run tests](#run-tests)
- [Upload to PyPI](#upload-to-pypi)
- [Check that all reports are created](#check-that-all-reports-are-created)
- [Run Providentia inside Docker container](#run-providentia-inside-docker-container)
- [Create the docs for the first time and upload to readthedocs](#create-the-docs-for-the-first-time-and-upload-to-readthedocs)
- [Generate the docs](#generate-the-docs)
- [Create conda environments](#create-conda-environments)
- [Memory Profiling Code](#memory-profiling-code)
- [Obtain Zenodo Filetrees](#obtain-zenodo-filetrees)

## Migrate repository from Gitlab to Github

Things to do in advance:

- Create a file that maps all users from Gitlab to Github. For example: avilanov in Gitlab is albavilanova in Github. You can use this file as reference https://docs.google.com/spreadsheets/d/1O7pKzNNfRjM8O-iP_fB2HM5AvJiBXq31FO40BLUuu_M/edit?gid=0#gid=0.

- Request permission to Francesco or Albert to become an owner in https://github.com/BSC-ES.

- Invite users who created issues in Gitlab as team members to BSC-ES group in Github.

- Create empty repository in Github.

Once you have completed these steps:

- Prepare to use node-gitlab-2-github repository, which will be used to transfer all issues and milestones:

```bash
git clone --mirror https://earth.bsc.es/gitlab/ac/providentia.git
cd providentia.git
git push --no-verify --mirror git@github.com:BSC-ES/providentia.git
git remote set-url --push origin git@github.com:BSC-ES/providentia.git
git fetch -p origin
git push --no-verify --mirror
```

- Install node-gitlab-2-github: You can try to install the original repository but I don't recommend it since so many things are failing. I fixed them in a forked version at https://github.com/albavilanova/node-gitlab-2-github.

```bash
git clone https://github.com/albavilanova/node-gitlab-2-github
cd node-gitlab-2-github
```

- Edit the settings.ts within the repository that you've just cloned: Leave projectId as 0 the first time you run `npm run start`, that way it will throw an error showing you all the available project ids and their corresponding project names. Then update `settings.ts` with your id and run again. Our `settings.ts` was:

```
import Settings from './src/settings';
export default {
  gitlab: {
    url: 'https://earth.bsc.es/gitlab',
    token: 'YOUR-GITLAB-ACCOUNT-TOKEN',
    projectId: 574,
    listArchivedProjects: true,
    sessionCookie: "***",
  },
  github: {
    owner: 'BSC-ES',
    ownerIsOrg: true,
    token: 'YOUR-GITHUB-ACCOUNT-TOKEN',
    token_owner: 'albavilanova',
    repo: 'providentia',
    recreateRepo: false,
  },
  s3: {
    accessKeyId: null,
    secretAccessKey: null,
    bucket: null,
    region: null,
  },
  usermap: {
    'dbowdalo': 'denebowdalo',
    'avilanov': 'albavilanova',
    'pserrano': 'PaulaSerranoSierra',
    'avilamir': 'albertvilabsc',
    'avradi': 'amaliavr',
    'cmeikle': 'cmeikle',
    'ctena': 'charlio86',
    'cgile': 'carlottagile',
    'cferruz': 'cferruz',
    'eemili': 'emanueleemili',
    'fmacchia': 'f-macchia',
    'fbeninca': 'fbeninca',
    'hpetetin': 'hervepetetin',
    'jescriba': 'jeronimoescribano',
    'jmassagu': 'jomassa',
    'jyun': 'JayoungYun',
    'kdeoliv': 'kdeoliv',
    'kserrade': 'kserradell',
    'lilic': 'atmosphericdust',
    'mguevara': 'mguevarabsc',
    'molid': 'molidg',
    'ojorba': 'orioljorba',
    'pcamps': 'paulacamps',
    'rcruzalv': 'Rafaelaalves15',
    'rgrodofz': 'raphaelgrzg',
    'rgaratac': 'rgaratac1',
    'rsousse': 'rsousse',
    'tvintimi': 'tito-vintimilla',
    'yyousef': 'yarayo',
    'hnavarro': 'hectornav',
  },
  projectmap: {},
  conversion: {
    useLowerCaseLabels: true,
    addIssueInformation: true,
  },
  transfer: {
    description: true,
    milestones: true,
    labels: true,
    issues: true,
    mergeRequests: true,
    releases: true,
  },
  dryRun: false,
  exportUsers: true,
  useIssueImportAPI: true,
  usePlaceholderMilestonesForMissingMilestones: true,
  usePlaceholderIssuesForMissingIssues: true,
  useReplacementIssuesForCreationFails: true,
  useIssuesForAllMergeRequests: false,
  filterByLabel: undefined,
  trimOversizedLabelDescriptions: false,
  skipMergeRequestStates: [],
  skipMatchingComments: [],
  mergeRequests: {
    logFile: './merge-requests.json',
    log: false,
  },
  commitMap: {
  }
} as Settings;
```

Session cookie is obtained from Developer tools -> Application -> Cookies -> _gitlab_session while being in the Gitlab page. Tokens for Gitlab and Github need to be obtained from your account.

- Run code

```bash
npm run start
Transferring Description
Transferring Milestones
Transferring Labels
Transferring Releases
Transferring Issues
...
```

- Once it is finished, go to Github and check if everything is correct. You might see that the commits are not linked to your accounts, that's because they are associated to your BSC emails and not your Gitlab emails. To solve this, create a .mailmap file in your repository. Our .mailmap-gitlab contained:

```
# FORMAT: <replace-with--name>  <replace-with-email>  <commit-name>  <commit-email>
# Omit commit-name or commit-email if same as replace-with.
# git log --pretty="%aN <%aE>%n%cN <%cE>" | sort | uniq

Alba Vilanova               <alba.vilanova@outlook.com>          Alba Vilanova             <alba.vilanova@bsc.es>
Dene Bowdalo                <denebowdalo@googlemail.com>         Dene Bowdalo              <dene.bowdalo@bsc.es>
Paula Serrano Sierra        <paulaserranosierra@gmail.com>       Paula Serrano Sierra      <paula.serrano@bsc.es>
Amalia Vradi                <amalia.vradi@bsc.es>                Amalia Vradi              <amalia.vradi@bsc.es>
Francesco Benincasa         <fbeninca@gmail.com>                 Francesco Benincasa       <francesco.benincasa@bsc.es>
```

Then:

```bash
git clone https://github.com/BSC-ES/providentia.git
cd providentia
git filter-repo --mailmap .mailmap --force
git remote add origin https://github.com/BSC-ES/providentia.git
git push --set-upstream origin master  --force
git push --mirror --force origin
```

It is also possible that the default branch is not set to be the master, you can change that from Settings -> General.

## Run tests

To run all the pipeline tests in your local machine and read the data from the folder `tests/data`, first, you will need to change your `get_machine()` function in `auxiliar.py` to this:

```python
def get_machine():
    """
    Get machine where code is running

    Returns
    -------
    str
        Machine
    """

    return "github"
```

Next, you will need to run:

```bash
conda activate providentia-env_v[version]
pytest tests
```

If you want to run a specific set of tests (from: test_apply_filter, test_make_plot, test_read_data, test_save, test_unit_converter), you can do so by specifying the set:

```bash
pytest tests/test_read_data.py
```

To run a specific test, you will need to edit these files and comment the functions you are not interested in testing.

If you want to recreate the expected data in the tests, set the variable `GENERATE_OUTPUT` in `tests/aux_functions.py` to True temporarily and run the tests once before reverting it back to False.

If you want to see the coverage report use:

```bash
coverage report -i -m
```

## Upload to PyPI

Every time a version is released, we should update Providentia in PyPI. To do so, we will first create a source distribution in the folder `dist``:

```bash
python setup.py sdist
```

We can check if it can be installed doing:

```bash
pip install dist/providentia-X.X.X.tar.gz
```

If everything is correct, we can proceed to upload our distribution to the website.

To upload the package we need to install twine:

```bash
pip install twine
```

And get an API key in PyPI. I recommend creating a folder in your home directory called .pypirc, with the following content:

```
[pypi]
  username = __token__
  password = pypi-[key]
```

That way you won't be asked for the your credentials anymore. Make sure the version name has been updated in __init__.py, and no more changes are needed since **it is not allowed to update the repository for that version once it is uploaded even if deleted**. Be very sure that everything is working fine at this point, as the uploaded version cannot be overwritten, only deleted.

To upload it, you can do:

```
twine upload dist/providentia-X.X.X.tar.gz
```

## Check that all reports are created

To generate reports for all configuration files under `configurations` folder:

```bash
#!/bin/bash
folder_path="configurations"
error_log="error_log2.txt"
> "$error_log"
for file in "$folder_path"/*; do
        command="./bin/providentia --report --config="$file
  echo
  echo
  echo
  echo "ejecutando: $command"
  output=$(eval $command 2>&1)
  if [ $? -ne 0 ]; then
  	error_type=$(echo "$output" | grep -o -E '[a-zA-Z]*Error')
  	echo "Command failed with exit code $?."
  	echo "File: $file Error: $error_type" >> "$error_log"
  else
    echo "Command executed successfully."
  fi
done
```

## Run Providentia inside Docker container

We need to first build our image:

```bash
docker build -t "providentia-image" .
```

To be able to display the dashboard we will also run these commands:
```bash
export DISPLAY=:1.0
xhost +local:*
```

We can then start the service using compose.yml:
```bash
docker compose run providentia bash
```

This will open a "session" inside our container, from which we can run Providentia.

Prior to that you should update the paths to the volumes (/data, /home/avilanov/software/Providentia) in `compose.yml`, used to access the data and repositories that are found in your local machine (outside the Docker container). The changes you do to your local files will be automatically reflected in the container.

## Create the docs for the first time and upload to readthedocs

Install `myst-parser`, `Sphinx` and `Sphinx-rtd-theme` in your environment and add the modules to your `requirements.txt` file.

```bash
conda activate providentia-env_v[version]
pip install myst-parser Sphinx Sphinx-rtd-theme
```

Create a folder called docs inside your repository and do:

```bash
cd docs
sphinx-quickstart
```

This will ask you a few questions (we answered with yes to the question of whether we want to separate source and build folders) and a `conf.py` file, together with other files, will be generated in the `docs/source` folder. You can edit the configuration file. In our case, we added some extensions and plugins:

```
# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["myst_parser", "sphinx.ext.todo", "sphinx.ext.viewcode", "sphinx.ext.autodoc", "sphinx_rtd_theme"]
templates_path = ['_templates']
exclude_patterns = []
myst_heading_anchors = 2

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
```

Add your markdown files to `docs/source` and add the .md filenames to `index.rst` in order to show their contents. Then run:

```bash
make clean
make html
```

Here you might get some errors related to the format of your .md files, fix them and run both commands again. This will generate the HTML files under `docs/build`, you can open them from the browser to make sure they look correct.

Now create an account and project in https://app.readthedocs.org/, and associate it with your project in Github, which must be open source. To show the documentation, you will need to create the file .readthedocs.yaml in your repository. Ours looks like this:

```
# Read the Docs configuration file
# See https://docs.readthedocs.io/en/stable/config-file/v2.html for details

# Required
version: 2

# Set the OS, Python version, and other tools you might need
build:
  os: ubuntu-24.04
  tools:
    python: "mambaforge-22.9"

# Build documentation in the "docs/" directory with Sphinx
sphinx:
   configuration: docs/source/conf.py

# Create environment
conda:
  environment: environment.yaml
```

Trigger a build from your latest commit containing this file and you should be ready to go.

## Generate the docs

Edit the files under `docs/source`, then navigate to docs and simply run:

```bash
conda activate providentia-env_v3.1.0
cd docs
make clean
make html
```

Do not edit anything under `docs/build` as it gets deleted everytime you run `make clean`.

## Create conda environments

### Create local environment

If for some reason you want to create the environment from scratch, you can use:

```
conda env create -f environment.yaml
```

You might get a warning like:

```
WARNING conda.models.version:get_matcher(556): Using .* with relational operator is superfluous and deprecated and will be removed in a future version of conda. Your spec was 1.6.0.*, but conda is ignoring the .* and treating it as 1.6.0
```

This can be removed by updating conda:

```
conda update conda
conda install -n base conda=24.4.0 conda-build=24.3.0
```

Check what the latest versions of [conda](https://github.com/conda/conda/releases) and [conda-build](https://github.com/conda/conda-build/releases) are.

What the first command does is creating an environment called `providentia-env_v3.1.0` with the Python version 3.11.5, and installing `cartopy`, `jupyterlab`, `ghostscript`, `dask`, `nco`, `pyproj` and `proj` with conda, and the Python packages from `requirements.txt` using pip. The equivalent would be:

```bash
conda create -n providentia-env_v3.1.0 python=3.11.5 -c conda-forge --override-channels
conda activate providentia-env_v3.1.0
conda install -c conda-forge cartopy=0.25.0 jupyterlab==4.5.3 ghostscript==10.06.0 dask==2026.1.2 nco==5.3.6 pyproj=3.7.2 proj=9.7.1 --override-channels
pip install -r requirements.txt
```

### Create providentia-env_v3.1.0-nord4 in Nord4

```bash
module unload intel
module load GCC/10.2.0
conda create -p /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.1.0-nord4 -y python=3.11.5 -c conda-forge --override-channels
conda activate /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.1.0-nord4
conda install -c conda-forge cartopy=0.25.0 jupyterlab==4.5.3 ghostscript==10.06.0 dask==2026.1.2 nco==5.3.6 pyproj=3.7.2 proj=9.7.1 --override-channels
pip install -r requirements.txt
```

### Create providentia-env_v3.1.0 in MN5

```bash
conda create -p /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.1.0 -y python=3.11.5 -c conda-forge --override-channels
conda activate /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.1.0
conda install -c conda-forge cartopy=0.25.0 jupyterlab==4.5.3 ghostscript==10.06.0 dask==2026.1.2 nco==5.3.6 pyproj=3.7.2 proj=9.7.1 --override-channels
pip install -r requirements.txt
```

## Memory Profiling Code

To memory-profile a function, you can use the `Tracker` method from the `memray` python module.

```python
import memray

def function_to_profile1():
    pass

with memray.Tracker("output_file"):
    function_to_profile1()

```
Use the `memray` command in the terminal to create a human-readable output.

### Setup

- If working **locally**, install Memray with:
```bash
python3 -m pip install memray
```
- If working on **MN5**, load the module with:
```bash
module load memray/1.17.1-foss-2023b-Python-3.11.5
```

### Generate Outputs

To create an interactive flamegraph, use:

```bash
memray flamegraph output_file 
```

To create a text summary, use:

```bash
memray table output_file 
```

To get the detailed stats, use:

```bash
memray stats output_file 
```

## Obtain Zenodo Filetrees

Here the code used to obtain the Zenodo filetrees for the release of the next GHOST versions in Zenodo. 

In order to create individual filetrees, run a modified version of the Providentia Zenodo download. 

First of all, the `download_ghost_network_zenodo` function only has to be runned once with `initial_check = False`. To do so comment these lines on `configuration.py`:

```python
# download GHOST network
#initial_check_nc_files = download_fun(network, initial_check=True)
#files_to_download = self.select_files_to_download(initial_check_nc_files)
#if not initial_check_nc_files or files_to_download:
    download_fun(network, initial_check=False, files_to_download=files_to_download)
```

This makes `download_ghost_network_zenodo` run the next code:

```python
# get the GHOST artifact value for the corresponding network
artifact_network = self.artifact_mapping[network]    

# create temporal dir to store the zip file and its tar components
self.temp_dir = os.path.join(self.download_instance.ghost_root, ".temp")
os.makedirs(self.temp_dir, exist_ok=True)

# download zip on the temporal directory
zip_path = self.download_zip(network, artifact_network)

# extract zip on the temporal directory
valid_files_info = self.extract_zip(files_to_download, zip_path, initial_check)

# extract tar on the temporal directory
self.extract_tar(valid_files_info)          
```

Before that modify the `extract_zip` and `extract_tar` functions from the `zenodo.py` module.

```python
def extract_zip(self, network, zip_path):        
  # open the ZIP file
  with ZipFile(zip_path, "r") as zipf:            
      zipf.extractall(self.temp_dir)

      # list all files inside the ZIP
      self.file_paths = [join(self.temp_dir,f) for f in zipf.namelist() if not f.endswith("/")]

def extract_tar(self, network):
  import tarfile
  import json

  json_dict = {}
  for tar_path in tqdm(self.file_paths,desc=f"    Checking {network} tars",):
      with tarfile.open(tar_path, "r:*") as tar:
          _, _, _, _, _, _, _, _, network, resolution, species_tar = tar_path.split('/')
          species = species_tar[:-7]

          if network not in json_dict:
              json_dict[network] = {}
          
          if resolution not in json_dict[network]:
              json_dict[network][resolution] = {}

          if species not in json_dict[network][resolution]:
              json_dict[network][resolution][species] = []

          for member in tar.getmembers():
              file_name = member.name.split('/')[-1]
              if file_name.endswith('.nc'):
                  json_dict[network][resolution][species].append(file_name)

          json_dict[network][resolution][species] = list(sorted(json_dict[network][resolution][species]))
      
  final_path = f"/home/pserrano/providentia/settings/internal/zenodo/zenodo_networks/{self.download_instance.ghost_version}/{network}.json"
  with open(final_path, "w", encoding="utf-8") as f:
      json.dump(json_dict, f, indent=2)

  with open(final_path, "r", encoding="utf-8") as f:
      data = json.load(f)
```

After all the networks have generated individual filetrees, join all of them by running this script:

```python
import os
import json
import yaml
import sys

ghost_version = '1.5.1'

path = f"/home/pserrano/providentia/settings/internal/zenodo/zenodo_networks/{ghost_version}"
ld = os.listdir(path)

yaml_path = f"/home/pserrano/providentia/settings/internal/zenodo/zenodo_{ghost_version}.yaml"
yaml_obj = yaml.safe_load(open(yaml_path))

valid_networks = yaml_obj.keys()
actual_networks = [i[:-5] for i in ld]

if set(valid_networks)-set(actual_networks) or set(actual_networks)-set(valid_networks):
    print(set(valid_networks)-set(actual_networks))
    print()
    print(set(actual_networks)-set(valid_networks))
    sys.exit()

data = {}

for json_file in sorted(ld):
    full_path = os.path.join(path, json_file)
    with open(full_path) as f:
        d = json.load(f)
    n = list(d.keys())[0]
    data[n] = d[n]

final_path = os.path.join(os.path.dirname(os.path.dirname(path)), f'zenodo_ghost_filetree_{ghost_version}.json')
with open(final_path, 'w') as f:
    json.dump(data, f,indent=4)
```

Too see more useful scripts regarding the filetrees [this](https://github.com/BSC-ES/providentia/issues/812) issue.