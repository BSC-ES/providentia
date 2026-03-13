# Git Bash

## Install Git Bash

Download the .exe file from https://git-scm.com/install/windows (Git-2.53.0-64-bit.exe) and execute. Click Yes if prompted and follow the installation steps. This will install Git CMD, Git Bash and and Git GUI, from these we will be using Git Bash.

## Install conda on Git Bash

Install Conda for Windows following the official guide: [https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend Miniconda because it is lightweight. If you get an error and your VPN is active, make sure to deactivate it before downloading it.

After downloading the .exe file, execute and follow the steps. In the step `Advanced Installation Options` check `Add installation to my PATH environment variable`.

Open the Git Bash application as an administrator and run:
```bash
conda init
```

Close terminal and reopen.

## Install git on Git Bash

Install git using:

```bash
sudo apt install git
```

## Install Providentia on Git Bash

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

## Download data onto Windows

If you want to store the data used by Providentia on your Windows machine, change the local paths in `providentia/settings/data_paths.yaml` to point to your local directories. For example:

```yaml
local:
  ghost_root: /Users/avilanov/Desktop/data/providentia/obs/ghost
  mod_root: /Users/avilanov/Desktop/data/providentia/mod
  mod_to_interp_root: /Users/avilanov/Desktop/data/providentia/mod_to_interp
  nonghost_root: /Users/avilanov/Desktop/data/providentia/obs/nonghost
```

This setup will make Providentia download to and read the data from the folder `data` inside Desktop, which will be created the first time you run the tool. By default, if you don't change the paths, it will get the data from `/home/{user}/data`, a folder inside the WSL. 

Remember that you can use the [download mode](Download) of Providentia to get data from Zenodo, ACTRIS and CAMS, and also from the MN5 if you have access to the machine.

If you have access to MN5 and want to download data from there without using a password, you will need to create a public SSH key. In order to do this, run this command from the from WSL:

```bash
ssh-keygen -t ed25519
```

You can click enter to save file to default location and keep the paraphrase empty. Now copy the content within `~/.ssh/id_ed25519.pub` on WSL inside `~/.ssh/autorized_keys` on MN5 (you can connect to it as you would usually do). Once that's done, you can download data easily.
