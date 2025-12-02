# Tips and tricks for developers

This page can be used by developers to better understand certain parts of the code:

- [Migrate from Gitlab to Github](#migrate-from-gitlab-to-github)
- [Run tests](#run-tests)
- [Upload to PyPI](#upload-to-pypi)
- [Check that all reports are created](#check-that-all-reports-are-created)
- [Run Providentia inside Docker container](#run-providentia-inside-docker-container)
- [Create the docs for the first time and upload to readthedocs](#create-the-docs-for-the-first-time-and-upload-to-readthedocs)
- [Generate the docs](#generate-the-docs)
- [Create providentia-env_v3.0.0 in MN5](#create-providentia-env_v300-in-mn5)
- [Create local environment](#create-local-environment)
- [Memory Profiling Code](#memory-profiling-code)

## Migrate from Gitlab to Github

Things to do in advance:

- Create a file that maps all users from Gitlab to Github. For example: avilanov is albavilanova. You can use this file as reference https://docs.google.com/spreadsheets/d/1O7pKzNNfRjM8O-iP_fB2HM5AvJiBXq31FO40BLUuu_M/edit?gid=0#gid=0.

- Request permission to Francesco or Albert to become an owner in https://github.com/BSC-ES.

- Invite users who created issues in Gitlab as team members to BSC-ES group in Github.

- Create empty repository in Github

Once you have done completed these steps:

- Prepare to use node-gitlab-2-github repository

```
git clone --mirror https://earth.bsc.es/gitlab/ac/providentia.git
cd providentia.git
git push --no-verify --mirror git@github.com:BSC-ES/providentia.git
git remote set-url --push origin git@github.com:BSC-ES/providentia.git
git fetch -p origin
git push --no-verify --mirror
```

- Install node-gitlab-2-github: You can try to install the original repository but I don't recommend it since so many things are failing. I fixed them in a forked version.

```
git clone https://github.com/albavilanova/node-gitlab-2-github
cd node-gitlab-2-github
```

- Create settings.ts: Leave projectId as 0 the first time you do `npm run start`so that it throws an error showing you all the project ids and their corresponding project names. Then update `settings.ts`. Our file was:

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

Session cookie is obtained from Developer tools -> Application -> Cookies -> _gitlab_session while being in the Gitlab page.

- Run code

```
npm run start
Transferring Description
Transferring Milestones
Transferring Labels
Transferring Releases
Transferring Issues
```

- Once it is finished, go to Github and check if everything is correct. You might see that the commits are not linked to your accounts, that's because they are associated to your BSC emails and not your Gitlab emails. To solve this, create a .mailmap and do the following:

```
git clone https://github.com/BSC-ES/providentia.git
cd providentia
git filter-repo --mailmap .mailmap --force
git remote add origin https://github.com/BSC-ES/providentia.git
git push --set-upstream origin master  --force
git push --mirror --force origin
```

Our .mailmap-gitlab contains:
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

It is also possible that the default branch is not set to be the master, you can change that from Settings -> General.

## Run tests

To run all the pipeline tests in your local machine, you will need to add `return "github"` in the function `get_machine()` in `auxiliar.py` and then run:

```
conda activate providentia-env_v3.0.0
pytest tests
```

If you want to run a specific set of tests (there are three sets: test_apply_filter, test_make_plot, test_read_data), you can do so by specifying the set:

```
pytest tests/test_read_data.py
```

To run a specific test, you will need to edit these files and comment the functions you are not interested in testing.

If you want to see the coverage report use:

```
coverage report -i -m
```

## Upload to PyPI

Every time a version is released, we should update Providentia in PyPI. To do so, we will first create a source distribution:

```
python setup.py sdist
```

We can check if it can be installed doing:
```
pip install dist/providentia-X.X.X.tar.gz
```

If everything is correct, we can proceed to upload our distribution to the website. For this we need to install twine:

```
pip install twine
```

And get an API key in PyPI. I recommend creating a folder in your home directory called .pypirc, with the following content:

```
[pypi]
  username = __token__
  password = pypi-[key]
```

That way you won't be asked for the your credentials anymore. Make sure the version name has been updated in __init__.py, and no more changes are needed since **it is not allowed to update the repository for that version once it is uploaded even if deleted**.

To upload it, you can do:

```
twine upload dist/providentia-X.X.X.tar.gz
```

## Check that all reports are created

```
#!/bin/bash
folder_path="configurations"
error_log="error_log2.txt"
> "$error_log"
for file in "$folder_path"/*; do
        command="./bin/providentia --offline --config="$file
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

```
docker build -t "providentia-image" .
```

To be able to display the dashboard we will also run these commands:
```
export DISPLAY=:1.0
xhost +local:*
```

We can then start the service without using compose.yml:
```
docker run -v /data:/data -v /home/avilanov/software/Providentia:/tmp/Providentia -v /home/avilanov/software/providentia-interpolation:/tmp/providentia-interpolation -it -p 8888:8888 providentia-image /bin/bash
```

Or simply doing
```
docker compose run providentia bash
```

This will open a "session" inside our container, from which we can run Providentia.

Prior to that you should update the paths to the volumes (/data, /home/avilanov/software/Providentia, /home/avilanov/software/providentia-interpolation), used to access the data and repositories that are found in your local machine (outside the Docker container). The changes you do to your local files will be automatically reflected in the container.

## Create the docs for the first time and upload to readthedocs

Install `myst-parser==3.0.1`, `Sphinx==7.2.6` and `Sphinx-rtd-theme==2.0.0` in your environment and add the modules to your `requirements.txt` file.

```
conda activate providentia-env_v3.0.0
pip install myst-parser==3.0.1 Sphinx==7.2.6 Sphinx-rtd-theme==2.0.0
```

Create a folder called docs inside your repository and do:
```
cd docs
sphinx-quickstart
```

This will ask you a few questions (we answered with yes to the question of whether we want to separate source and build folders) and a `conf.py` file together with other files will be generated in the `docs/source` folder. You can edit the configuration file. In our case, we added some extensions and plugins:

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

Add your markdown files to `docs/source` and add the filenames to `index.rst` in order to show their contents. Then run:

```
make clean
make html
```

Here you might get some errors related to the format of your .md files, fix them and run the commands again. This will generate the HTML files, you can open them from the browser to make sure they look correct.

Now create an account and project in https://app.readthedocs.org/, and associate it with your project, which must be open source. To show the documentation, you will need to create the file .readthedocs.yaml in your repository. Ours looks like this:

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

Edit the files under `docs/source``, then navigate to docs and simply run:
```
conda activate providentia-env_v3.0.0
cd docs
make clean
make html
```

## Create providentia-env_v3.0.0-nord4 in Nord4

```
conda create -p /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.0.0-nord4 -y python=3.11.5 -c conda-forge --override-channels
conda activate /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.0.0-nord4
conda install -c conda-forge cartopy jupyterlab ghostscript dask --override-channels
pip install -r requirements.txt
```

## Create providentia-env_v3.0.0 in MN5

```
conda create -p /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.0.0 -y python=3.11.5 -c conda-forge --override-channels
conda activate /gpfs/projects/bsc32/repository/apps/conda_envs/providentia-env_v3.0.0
conda install -c conda-forge cartopy jupyterlab ghostscript dask --override-channels
pip install -r requirements.txt
```

## Create local environment

In order to test modules with pip, you need to create an environment. Once activated, you can start installing modules using either `pip` or `conda`, as in this example:

```
conda create -n providentia-env_v3.0.0 python=3.11.5 -c conda-forge --override-channels
conda activate providentia-env_v3.0.0
conda install -c conda-forge cartopy jupyterlab ghostscript dask --override-channels
pip install -r requirements.txt
```

Here we create an environment called `providentia-env_v3.00` with the Python version 3.11.5, and we install the latest version of Cartopy (with conda, with pip it gives problems), and the Python packages from `requirements.txt` using pip.


pip install memray

memray run -o result.bin python your_script.py
memray flamegraph result.bin  # Creates a flamegraph
memray table result.bin       # Text summary
memray stats result.bin       # Detailed stats


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