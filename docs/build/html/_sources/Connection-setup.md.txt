# Remote HPC connection setup

For users who do not wish to run Providentia locally, there is the option to run it on a remote HPC.

For BSC users, Providentia has been designed to be run optimally on the **MareNostrum5** and **Nord4** HPC machines. It is recommended that users run Providentia on **MareNostrum5** as it is guaranteed that there will be no lack of computational resources. 

As users connect to these remote machines by ssh, the graphical refresh speed of Providentia is limited by x11 forwarding, thus in order to exploit the full potential of the Providentia, the connections to these machines should be setup in an optimal fashion.

Note that not all remote BSC HPC machines have access to the same data locations. In case you want to load observations that are located on the `/esarchive`, then you should run Providentia on **Nord4**, as **MareNostrum5** does not have access to the `/esarchive`.

## .ssh/config

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

## .bashrc for BSC HPC machines

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

## External connections (VPN)

If you are connecting to a HPC machine externally (i.e. via VPN), then the ssh connection configuration on the local machine to the gateway machine you are connecting to should be setup in the exact same way.