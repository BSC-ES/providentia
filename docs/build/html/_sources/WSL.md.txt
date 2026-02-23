# Windows Subsystem for Linux

## Prerequisites

Go to System Information and make sure your Windows OS is in the Pro version, and not Home, otherwise you won't have access to the virtualisation feature and will need to use a Virtual Machine instead of using the WSL.

Go to Control panel -> Programs and Features -> Turn Windows features on and off and check if the options `Windows subsystem for Linux` and `Hyper-V` are checked, if they are not check them, apply the changes and restart your computer.

## Steps to install WSL

### 1. Open the Windows Powershell

As administrator, click Yes if prompted.

### 2. Run `wsl`

If you have checked the options Windows subsystem for linux and Hyper-V following the prerequisites instructions, there should be no problem here. If the command is not recognised, check your Windows version using the command `winver` and make sure your Windows installation is new. Otherwise, go to Update Windows through the control panel and proceed to update the system.

### 3. Run `wsl --status`

This should show:

```bash
Default Version: 2
```

It might not return anything, that's OK.

If you have updated your Windows version in step 2, it might suggest that you upgrade wsl by doing:

```bash
wsl --update
```

If when running `wsl --status` this error is shown:

```
This application requires the Windows Subsystem for Linux Optional Component.
Install it by running: wsl.exe --install --no-distribution
The system may need to be restarted so the changes can take effect.
Error code: Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED
```

Try installing the distribution doing `wsl.exe --install Ubuntu-22.04`. If it fails again with same error, open the control panel, navigate to Programs and Features and click on Turn Windows features on and off and uncheck the options Windows subsystem for linux and Hyper-V, click OK and restart. Then, open Turn Windows features on and off again and check the options and restart your computer. Now, if you open the Powershell and type `wsl --status` you should finally get:

```bash
Default Version: 2
```

It might not return anything, that's OK.

### 4. Install your Linux distribution

If the installation failed before or you haven't done this yet, run:

```bash
wsl.exe --install Ubuntu-22.04
```

If it returns the Arguments for running Linux binaries, try running this instead:

```bash
wsl --install -d Ubuntu-22.04
```

This will ask you to create a username and password, remember them because you will need them to install the packages required by Providentia.

If you get the error:

```bash
wsl: Using legacy distribution registration. Consider using a tar based distribution instead.
Downloading: Ubuntu 22.04 LTS
Ubuntu 22.04 LTS has been downloaded.
Distribution successfully installed. It can be launched via 'wsl.exe -d Ubuntu 22.04 LTS'
Launching Ubuntu 22.04 LTS...
Installing, this may take a few minutes...
WslRegisterDistribution failed with error: 0x80370102
Please enable the Virtual Machine Platform Windows feature and ensure virtualization is enabled in the BIOS.
For information please visit https://aka.ms/enablevirtualization
Press any key to continue...
The installation process for distribution 'Ubuntu-22.04' failed with exit code: 1.
```

It means that your virtualisation is not activated. This is probably due to your Windows not being the Pro version and therefore not having access to the Hyper-V feature that had to be checked in the prerequisites instructions.

Once it is installed, you can open the terminal running `wsl` from PowerShell or by opening the application named WSL. From the terminal you can access your Windows directories doing `cd /mnt/`. From outside, in Windows, you can see the folders of the subsystem in a section called Linux in your file explorer.

## 5. Install conda on WSL

Once WSL is installed, we need to install conda. We recommend Miniconda because it is lightweight. You can download the .sh file **for Linux** from [the official website](https://www.anaconda.com/docs/getting-started/miniconda/install#linux-2) after creating an account. If you get an error and your VPN is active, make sure to deactivate it before downloading it.

From the terminal on the WSL run:

```bash
cd /mnt/c/Users/{user}/Downloads
bash Miniconda3-latest-Linux-x86_64.sh
```

If conda is not detected after the installation finishes, run:

```bash
echo 'export PATH="/home/{user}/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
conda init
```

Close your WSL and reopen it.

## 6. Install Providentia on WSL

We can now install Providentia:

```bash
cd /mnt/c/Users/{user}/Desktop (your preferred location)
git clone https://github.com/BSC-ES/providentia.git
```

Navigate to your Providentia directory and launch the dashboard:

```bash
cd providentia
./bin/providentia
```

If you get the error `/usr/bin/env: 'bash\r': No such file or directory`, run:

```bash
sudo apt update
sudo apt install dos2unix
dos2unix ./bin/providentia
```

If you get the error `Could not load the Qt platform plugin "xcb" on local machine`, run:

```bash
sudo apt update
sudo apt install --reinstall libxcb-xinerama0 libxkbcommon-x11-0 libx11-xcb1 libxrender1 libxi6 libxcb1 libxext6
sudo apt install qtbase5-dev qtbase5-dev-tools qt5-qmake libqt5gui5 libqt5widgets5 libqt5core5a
```

Now you should be able to open the dashboard:

```bash
./bin/providentia
```

It will be blank as you don't have any data on the WSL. 

## 7. Download data onto your WSL

If you want to store the data used by Providentia on your Windows machine, change the local paths in `providentia/settings/data_paths.yaml` to point to your local directories under `/mnt`. For example:

```yaml
local:
  ghost_root: /mnt/c/Users/{user}/Desktop/data/providentia/obs/ghost
  mod_root: ~/mnt/c/Users/{user}/Desktop/data/providentia/mod
  mod_to_interp_root: /mnt/c/Users/{user}/Desktop/data/providentia/mod_to_interp
  nonghost_root: /mnt/c/Users/{user}/Desktop/data/providentia/obs/nonghost
```

This setup will make Providentia download to and read the data from the folder `data` inside Desktop, which will be created the first time you run the tool. By default, if you don't change the paths, it will get the data from `/home/{user}/data`, a folder inside the WSL. 

Remember that you can use the [download mode](Download) of Providentia to get data from Zenodo, ACTRIS and CAMS, and also from the MN5 if you have access to the machine.

If you have access to MN5 and want to download data from there without using a password, you will need to create a public SSH key. In order to do this from WSL run:

```bash
ssh-keygen -t ed25519
```

You can keep the paraphrase empty. Now copy the content within `~/.ssh/{yourkeyname}.pub` on WSL inside `~/.ssh/autorized_keys` on MN5 (you can connect to it as you would usually do). Once that's done, you can download data easily.