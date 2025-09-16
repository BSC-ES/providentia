# Connection setup

Providentia has been designed to be run optimally on MareNostrum 5 and Nord4 machines. Therefore, in order to exploit the full potential of the Providentia tool, the connections to these machines should be setup in an optimal fashion.

As the users connect to these machines by ssh, the graphical refreshing speed of the tool is limited by the x11 forwarding.

It is recommended that users run Providentia on MareNostrum 5 as it is guaranteed that there will be no lack of computational resources. 

Note that not all the machines have access to the same data locations. In case you want to load observations that are located in `/esarchive`, then you should run Providentia in Nord4, as MareNostrum 5 does not have access to the `/esarchive`.

## .ssh/config Optimal Setup

To optimise the ssh connection to MareNostrum 5 and Nord4, you will modify the `.ssh/config` file in the `$HOME` directory of the gateway machine used to connect to the machines. For all gateway machines the `.ssh/config` connection setup should follow the following templates:

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

## .bashrc

You should have this in your .bashrc to be able to run Providentia in MN5 and Nord4 without problems:

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

If you are connecting to a BSC machine externally (i.e. via VPN), then the ssh connection configuration on the local machine to the gateway machine you are connecting to (e.g. to a workstation) should be setup in the exact same way.