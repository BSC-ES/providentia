# Website

The website mode runs Providentia's website on your machine, together with the backend that reads the data and makes the plots.

## Getting started

The website mode is not available to users yet; without `WEB=1`, `--website` only says so. To start the website, run the following command in the command line:

```
WEB=1 ./bin/providentia --website
```

It runs the checkout of ines-providentia at `WEBSITE_PATH`, at the top of `bin/providentia` (`${HOME}/software/ines-providentia` by default; change it if yours is elsewhere), from its folder: its `src` is put on `PYTHONPATH`, also for the backend it starts, and the ines-providentia installed in the environment is not used. Its dependencies must be in the environment, as they are when it was created with `requirements.txt`.

It starts the backend, which loads Providentia and can take a few minutes, and then opens the dashboard of the website in your browser (`http://127.0.0.1:7782/...`). The configuration file and section are not used in this mode.

To stop the website and the backend, press `Ctrl+C` in the command line. If the backend stops on its own, the page tells you to stop the command and start it again. The `LOGS` button of the page shows what the backend wrote, to see why it stopped.

## Requirements

- The website mode is only available on local machines, not on MN5, Nord4 or the operational machine.
- Port 7782 must be free on your machine.
- The website package (`ines-providentia`) is installed with the `providentia-env_v4.0.0` environment from BSC's GitLab (public, https://gitlab.earth.bsc.es).
