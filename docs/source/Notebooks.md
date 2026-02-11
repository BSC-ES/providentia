# Notebooks

As the library mode allows the importing of Providentia as a module, it follows that Providentia can therefore also be used in a Jupyter notebook.

## Launching Jupyter notebook 
A Jupyter notebook can be launched with the following command:
```
./bin/providentia --notebook
```         
If you are on a local machine, this will directly open Jupyter in your web browser or if you are on a HPC it will start a job on the machine with SLURM. 

For the HPC case, a file named `notebook.out` will be created. Information in the file will then allow for a connection to a Jupyter notebook to be set up. Firstly, an SSH tunnel from the local machine needs to be set up, by pasting a given command into a Linux/MacOS terminal or equivalent (i.e. PuTTY on Windows), e.g.:

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

## Importing Providentia
 
If the notebook you are running is not inside the Providentia home directory, in order to import from Providentia, then it is necessary to add these lines to your code at the start:            
```
import sys
sys.path.append(provdir)
```      
where `provdir` is the path where your Providentia code exists.         

You will also need to load the necessary modules for it to function, for this you can activate the local environment that is created the first time you run the tool:

```
conda activate providentia-env_v[version]
```

where `version` is the relevant version of Providentia, e.g. `3.0.2`.

Then the Providentia library can be safely imported as a module:

```
import providentia as prv
```

Please see the [library](Library) page for a full description of the features of the providentia module.

## Embedding plots

If you wish to embed plots within the notebook, you must set the following code, **AFTER** importing Providentia: 
     
```
%matplotlib inline
```              
This is a magic function that allows for the rendering of figures directly in the notebook.        

## Tutorials

We have a number of tutorial notebooks that can be used to learn about Providentia's features, these can all be found in the [tutorials](https://github.com/BSC-ES/providentia/tree/master/tutorials) subdirectory. 