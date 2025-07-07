# Notebook

As the library mode allows the importing of Providentia as a module, it follows that Providentia can therefore also be used in a Jupyter notebook. This type of interactive working is growing ever more popular, and thus it is important that Providentia can adapt to this type of operation.

## Launching a Jupyter notebook
A Jupyter notebook can be launched with the following command:
```
./bin/providentia --notebook
```         
This will either start a job if you are using a machine with SLURM, or directly open Jupyter notebooks. In the first case, a file named `notebook.out` will be created. Information in the file will then allow for a connection to a Jupyter notebook to be set up. Firstly, an SSH tunnel from the local machine needs to be set up, by pasting a given command into a Linux/MacOS terminal or equivalent (i.e. PuTTY on Windows), e.g.:

```
Create an SSH tunnel via terminal on your local machine:
ssh -N -L 8825:s04r1b14:8825 bsc32025@nord4.bsc.es
```
Secondly, the notebook can than be accessed via a web browser following information from the bottom of the file, e.g.:     
```
To access the server, open this file in a browser:
        file:///.statelite/tmpfs/gpfs/home/bsc32/bsc32025/.local/share/jupyter/runtime/jpserver-800644-open.html
    Or copy and paste one of these URLs:
        http://s04r1b14:8825/lab?token=57a994dfaecf92e5aefd1d5dbc5a9831fb8060f873fd8181
     or http://127.0.0.1:8825/lab?token=57a994dfaecf92e5aefd1d5dbc5a9831fb8060f873fd8181
```
One additional step is also needed inside the notebook if wishing to embed plots:     
```
%matplotlib inline
```              
This is a magic function that allows for the rendering of figures directly in the notebook.        
NOTE: This must set after the importing of the module.

A template Jupyter notebook, demonstrating some of the Providentia interactive features can be found in `notebooks/interactive_template.ipynb`.

## Importing Providentia
 
If the notebook you are running is not inside the Providentia home directory, in order to import from Providentia, then it is necessary to add these lines to your code first:            
```
import sys
sys.path.append(provdir)
```      
where `provdir` is the path where your Providentia code exists.         

You will also need to load the necessary modules for it to function, for this you can run:

```
source provdir/bin/load_modules.sh
```

After the import you will have full access to the backend of Providentia and all of its functions as described in the following sections.

NOTE: In the future we will allow for Providentia to be imported as a module on the BSC machines, making this step redundant. 
