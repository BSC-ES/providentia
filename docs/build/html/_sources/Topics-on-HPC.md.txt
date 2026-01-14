# Topics on HPC

(Connection-setup)=
## Connection setup

### Remote HPC connection setup

For users who do not wish to run Providentia locally, there is the option to run it on a remote HPC.

For BSC users, Providentia has been designed to be run optimally on the **MareNostrum5** and **Nord4** HPC machines. It is recommended that users run Providentia on **MareNostrum5** as it is guaranteed that there will be no lack of computational resources. 

As users connect to these remote machines by ssh, the graphical refresh speed of Providentia is limited by x11 forwarding, thus in order to exploit the full potential of the Providentia, the connections to these machines should be setup in an optimal fashion.

Note that not all remote BSC HPC machines have access to the same data locations. In case you want to load observations that are located on the `/esarchive`, then you should run Providentia on **Nord4**, as **MareNostrum5** does not have access to the `/esarchive`.

#### .ssh/config

To optimise the ssh connection to the remote HPC machine, you will modify the `.ssh/config` file in the `$HOME` directory of the local machine used to connect to the remote machines following the examples for the connection to the **MareNostrum5** and **Nord4** BSC machines: 

```
Host mn5
    HostName glogin4.bsc.es
    User bscXXXXXX
    IdentityFile ~/.ssh/id_rsa
    ForwardX11 yes
    ForwardX11Trusted yes
    Compression yes
    Ciphers aes128-gcm@openssh.com
    ForwardX11Timeout 7d
```

```
Host nord4
    HostName n4login0.bsc.es
    User bscXXXXXX
    IdentityFile ~/.ssh/id_rsa
    ForwardX11 yes
    ForwardX11Trusted yes
    Compression yes
    ForwardX11Timeout 7d
    Ciphers aes128-gcm@openssh.com
```

`ForwardX11 yes` ensures x11 forwarding is active.  
`ForwardX11Trusted yes` ensures X11 forwarding is trusted.            
`Compression yes` compresses the connection (i.e. speeds up).          
`Ciphers aes128-gcm@openssh.com` sets a speed-optimised cipher for the connection to ensure less time is spent on encryption/decryption.                 
`ForwardX11Timeout 7d` ensures the x11 forwarding does not time out (it times out by default after 15 minutes on some machines). 

#### .bashrc for BSC HPC machines

For BSC users, you should also have these lines in your `.bashrc` file in the `$HOME` directory of the remote HPC machine:

```
if [ "$BSC_MACHINE" == "mn5" ]; then
    module use /gpfs/projects/bsc32/software/rhel/9.2/modules/all
elif [ $BSC_MACHINE == "nord4" ]; then
   module load nord3-singu
   module load bsc/current
   module use /gpfs/projects/bsc32/software/suselinux/11/modules/all
   unset PYTHONSTARTUP
fi
```

#### External connections (VPN)

If you are connecting to a HPC machine externally (i.e. via VPN), then the ssh connection configuration on the local machine to the gateway machine you are connecting to should be setup in the exact same way.

(Running-the-tool-on-debug)=
## Running the tool on debug

If you are on Nord4 or Marenostrum 5, this will request an interactive session:

```
./bin/providentia
``` 

After some seconds you will have entered onto an allocated node and the dashboard will be launched. If you want to launch it multiple times and avoid waiting in the queue you can use the debug mode. To do this, you will need to run the tool using the `--debug` flag and then rerun it:

```
./bin/providentia --debug
./bin/providentia
```

(Mounting-filesystems)=
## Mounting HPC filesystems

It is possible to mount filesystems, such as esarchive from MN5, on your personal computer. This can then be used to load observations/models that are stored remotely on these filesystems while using Providentia on your local machine. The Providentia data paths should be updated to point to the local mounted directory in `settings/data_paths.yaml`.

How this is done varies by operating system. In this section you will find information on how to mount esarchive, you can use this as inspiration for other filesystems.

### Mounting esarchive

#### Linux

On the Linux the esarchive can be setup permanently. 

Firstly, create a local directory where you want the esarchive to be mounted, e.g. `/localdir/tomount/esarchive`.

Next, add the following line to the `/etc/fstab` file (changing your BSC username and local directory):

`sshfs#bscXXXXXX@dt01.bsc.es:/gpfs/archive/bsc32/esarchive /localdir/tomount/esarchive fuse defaults,allow_other,reconnect,ro 0 0`

Install sshfs if it is not in your computer yet:

`sudo apt-get install sshfs`

The esarchive can finally be mounted by:

`sudo mount -av `

To unmount the esarchive, the following command can be used:

`sudo umount /localdir/tomount/esarchive`

#### Mac

On Mac, the mounting cannot be setup permanently, but will stay open while the terminal session is active. Before trying to setup the connection you must install FUSE and SSHFS, please follow the instructions here: https://earth.bsc.es/wiki/doku.php?id=computing:addesarchive&s[]=sshfs.

After that is done, create a local directory where you want the esarchive to be mounted, e.g. `/localdir/tomount/esarchive`.

Then run the following command (changing your BSC username and local directory): 

`sshfs bscXXXXXX@dt01.bsc.es:/gpfs/archive/bsc32/esarchive /localdir/tomount/esarchive -o volname=esarchive,defer_permissions,allow_other,IdentityFile=$HOME/.ssh/id_rsa,reconnect,ro`

To unmount the esarchive, the following command can be used:

`diskutil umount force /localdir/tomount/esarchive`

#### Windows

As present, the esarchive cannot be mounted on Windows machines.
