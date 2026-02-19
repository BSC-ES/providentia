# Windows Subsystem for Linux

## Steps to install WSL

### 1. Open the Windows Powershell

As administrator, click Yes if prompted.

### 2. Run `wsl`

If command is not recognised open the control panel, navigate to Programs and Features and click on Turn Windows features on and off, search for the options Windows subsystem for linux and Hyper-V, check them, click OK and restart your laptop. Run Windows Powershell as administrator and try to run `wsl` again. If it still says that it is not recognised, check your Windows version using the command `winver` and make sure your Windows installation is new. Otherwise, go to Update Windows through the control panel and proceed to update the system.

### 3. Run `wsl --status`

This should show:
```
Default Version: 2
```

If you have updated your Windows version in step 2, it might suggest that you upgrade wsl by doing `wsl --update`.

If when running `wsl --status` this error is shown:

```
This application requires the Windows Subsystem for Linux Optional Component.
Install it by running: wsl.exe --install --no-distribution
The system may need to be restarted so the changes can take effect.
Error code: Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED
```

Try installing the distribution doing `wsl.exe --install Ubuntu-22.04`. If it fails again with same error, open the control panel, navigate to Programs and Features and click on Turn Windows features on and off and uncheck the options Windows subsystem for linux and Hyper-V, click OK and restart. Then, open Turn Windows features on and off again and check the options and restart your computer. Now, if you open the Powershell and type `wsl --status` you should finally get:

```
Default Version: 2
```

### 4. Install your Linux distribution:
If the installation failed before, rerun:

```
wsl.exe --install Ubuntu-22.04
```

This will ask you to create a username and password, remember them because you will need them to install the packages required by Providentia.

If you get the error:

```
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
It means that your virtualisation is not activated.

6. Once it is installed, you can open the terminal running `wsl` from PowerShell or by opening the application named WSL. From inside you can see your Windows directories doing `cd /mnt/`. From outside, in Windows, you can see the folders of the subsystem in a section called Linux in your file explorer.

## Install Providentia on WSL

We can now install Providentia:

```bash
cd /mnt/c/Users/YOURUSER/Desktop (your preferred location)
git clone https://github.com/BSC-ES/providentia.git
```

We also need to install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend Miniconda because it is lightweight. After downloading the .sh file for Linux, run the installer from your terminal:

```bash
cd /mnt/c/Users/YOURUSER/Downloads
bash Miniconda3-latest-Linux-x86_64.sh
```

If conda is not detected after the installation finishes run:

```bash
echo 'export PATH="/home/YOURUSERINWSL/miniconda3/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
conda init
```
Close your WSL and reopen it.

Navigate to your Providentia directory and launch the dashboard:

```bash
cd /mnt/c/Users/YOURUSER/Desktop/providentia (your preferred location)
./bin/providentia
```

If you get the error `Could not load the Qt platform plugin "xcb" on local machine`, run:

```
sudo apt update
sudo apt install --reinstall libxcb-xinerama0 libxkbcommon-x11-0 libx11-xcb1 libxrender1 libxi6 libxcb1 libxext6
sudo apt install qtbase5-dev qtbase5-dev-tools qt5-qmake libqt5gui5 libqt5widgets5 libqt5core5a
```

Edit your .bashrc file and add the line `export QT_QPA_PLATFORM_PLUGIN_PATH=/usr/lib/x86_64-linux-gnu/qt5/plugins/platform`:

```bash
vi ~/.bashrc
source ~/.bashrc
```

Now you should be able to open the dashboard:

```
./bin/providentia
```

It will be blank as you don't have any data on the WSL. If you have data in your Windows machine, change the paths in `providentia/settings/data_paths.yaml` to point to your local directories under `/mnt`. Remember that you can use the [download mode](Download) of Providentia to get data from Zenodo, ACTRIS and CAMS, and also from the MN5 if you have access to the machine.