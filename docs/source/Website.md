# Website

The website mode runs Providentia's website on your machine, together with the backend that reads the data and makes the plots.

## Getting started

To start the website, run the following command in the command line:

```
./bin/providentia --website
```

It starts the backend, which loads Providentia and can take a few minutes, and then opens the dashboard of the website in your browser (`http://127.0.0.1:7782/...`). The configuration file and section are not used in this mode.

To stop the website and the backend, press `Ctrl+C` in the command line. If the backend stops on its own, the page tells you to stop the command and start it again. The `LOGS` button of the page shows what the backend wrote, to see why it stopped.

## Requirements

- The website mode is only available on local machines, not on MN5, Nord4 or the operational machine.
- Port 7782 must be free on your machine.
- The website packages (`ines-core` and `ines-providentia`) are installed with the `providentia-env_v4.0.0` environment from BSC's GitLab (https://gitlab.earth.bsc.es), so creating the environment needs access to it.
