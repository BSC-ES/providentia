# Git Bash

## Steps to install Git Bash

Download the .exe file from https://git-scm.com/install/windows (Git-2.53.0-64-bit.exe) and execute. Click Yes if prompted and follow the installation steps. This will install Git CMD, Git Bash and and Git GUI, from these we will be using Git Bash.

## Install Providentia on Git Bash

Install Conda for Windows following the official guide: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend Miniconda because it is lightweight. If you get an error and your VPN is active, make sure to deactivate it before downloading it.

After downloading the .exe file, execute and follow the steps. In the step `Advanced Installation Options` check `Add installation to my PATH environment variable`.

Open the Git Bash application as an administrator and run:
```
conda init
```

Close terminal and reopen.

Now install Providentia:

```bash
cd Desktop (your preferred location)
git clone https://github.com/BSC-ES/providentia.git
```

Open the dasboard doing:
```bash
cd providentia
./bin/providentia
```

It will be blank as you don't have any data on the WSL. If you have data in your Windows machine, change the local paths in `providentia/settings/data_paths.yaml` to point to your local directories:

```
local:
  ghost_root: /Users/avilanov/Desktop/data/providentia/obs/ghost
  mod_root: /Users/avilanov/Desktop/data/providentia/mod
  mod_to_interp_root: /Users/avilanov/Desktop/data/providentia/mod_to_interp
  nonghost_root: /Users/avilanov/Desktop/data/providentia/obs/nonghost
```

This will download to the disk C.

Remember that you can use the [download mode](Download) of Providentia to get data from Zenodo, ACTRIS and CAMS, and also from the MN5 if you have access to the machine.
