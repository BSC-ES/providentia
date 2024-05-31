# Mounting esarchive

It is possible to mount the esarchive on your personal computer. This can then be used to load observations/experiments that are stored remotely on the esarchive while using Providentia on your local machine. How this is done varies by operating system, as detailed below. The Providentia data paths should be updated to point to the local mounted directory in `settings/data_paths.yaml` 

## Linux

On the Linux the esarchive can be setup permanently. 

Firstly, create a local directory where you want the esarchive to be mounted, e.g. `/localdir/tomount/esarchive`

Next, add the following line to the `/etc/fstab` file (changing your BSC username and local directory):

`sshfs#bscXXXXXX@dt01.bsc.es:/gpfs/archive/bsc32/esarchive /localdir/tomount/esarchive fuse defaults,allow_other,reconnect,ro 0 0`

The esarchive can finally be mounted by:

`sudo mount -av `

To unmount the esarchive, the following command can be used:

`sudo umount /localdir/tomount/esarchive`

## Mac

On Mac, the mounting cannot be setup permanently, but will stay open while the terminal session is active. Before trying to setup the connection you must install FUSE and SSHFS, please follow the instructions here: https://earth.bsc.es/wiki/doku.php?id=computing:addesarchive&s[]=sshfs

After that is done, create a local directory where you want the esarchive to be mounted, e.g. `/localdir/tomount/esarchive`

Then run the following command (changing your BSC username and local directory): 

`sshfs bscXXXXXX@dt01.bsc.es:/gpfs/archive/bsc32/esarchive /localdir/tomount/esarchive -o volname=esarchive,defer_permissions,allow_other,IdentityFile=$HOME/.ssh/id_rsa,reconnect,ro`

To unmount the esarchive, the following command can be used:

`diskutil umount force /localdir/tomount/esarchive`

## Windows

As present, the esarchive cannot be mounted on Windows machines.