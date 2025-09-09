# Running the tool on debug

If you are on Nord4 or Marenostrum 5, this will request an interactive session:

```
./bin/providentia
``` 

After some seconds you will have entered onto an allocated node and the dashboard will be launched. If you want to launch it multiple times and avoid waiting in the queue you can use the debug mode. To do this, you will need to run the tool using the `--debug` flag and then rerun it:

```
./bin/providentia --debug
./bin/providentia
```