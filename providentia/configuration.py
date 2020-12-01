""" Providentia Configuration Module """

import os
import re
import subprocess


MACHINE = os.environ.get('BSC_MACHINE', '')

def parse_path(dir, f):
    #print "Opening data file", f
    if os.path.isabs(f):
        return f
    else:
        #log.info("Input: %s", f)
        return os.path.join(dir, f)


class ProvConfiguration(object):
    """ Configuration parameters definitions """

    def __init__(self, **kwargs):
        self.config_dir = kwargs.get('config_dir', \
             os.path.join(os.environ['HOME'], '.providentia'))
        self.ghost_version = kwargs.get('ghost_version', '1.3.1')
        self.cartopy_data_dir = kwargs.get('cartopy_data_dir', '')
        self.available_cpus = kwargs.get('available_cpus', '')
        self.n_cpus = kwargs.get('n_cpus', '')
        self.obs_root = kwargs.get('obs_root', '')
        self.nonghost_root = kwargs.get('nonghost_root', '')
        self.exp_root = kwargs.get('exp_root', '')
        self.sequential_colourmap = kwargs.get('sequential_colourmap',
                                               'viridis')
        self.sequential_colourmap_warm = \
                kwargs.get('sequential_colourmap_warm', 'Reds')

        self.diverging_colourmap = kwargs.get('diverging_colourmap', 'bwr')
        self.unsel_station_markersize = \
                                   kwargs.get('unsel_station_markersize', 3)
        self.sel_station_markersize = \
                                   kwargs.get('sel_station_markersize', 8)
        self.legend_markersize = kwargs.get('legend_markersize', 11)
        self.time_series_markersize = \
                                   kwargs.get('time_series_markersize', 1.1)
        self.temp_agg_markersize = \
                                   kwargs.get('temp_agg_markersize', 3)
        self.temp_agg_expbias_markersize = \
                kwargs.get('temp_agg_expbias_markersize', 3)
        self.map_coastline_resolution = \
                kwargs.get('map_coastline_resolution', 'low')
        self.available_networks = \
                kwargs.get('available_networks', "['AERONET_v3','EBAS','EEA_AQ_eReporting','NCDC_ISD','NCDC_ISD_EU','NCDC_ISD_IP','NCDC_ISD_NA']")

    def __setattr__(self, key, value):
        super(ProvConfiguration, self).__setattr__(key, self.parse_parameter(key, value))

    def parse_parameter(self, key, value):
        """ parse parameters """

        # get available N CPUs
        if key == 'available_cpus':
            if (MACHINE == 'power') or (MACHINE == 'mn4'):
                bash_command = 'squeue -h -o "%C"'
                process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
                output, _ = process.communicate()
                return int(re.findall(r'\d+', str(output))[0])

            value = int(os.cpu_count())

        elif key == 'cartopy_data_dir':
            # set cartopy data directory (needed on CTE-POWER/MN4 as has no external
            # internet connection)

            if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3'):
                value = '/gpfs/projects/bsc32/software/rhel/7.5/ppc64le/POWER9/software/Cartopy/0.17.0-foss-2018b-Python-3.7.0/lib/python3.7/site-packages/Cartopy-0.17.0-py3.7-linux-ppc64le.egg/cartopy/data'
            # on all machines except CTE-POWER/MN4, pull from internet

        elif key == 'n_cpus':
            # Define number of CPUs to process on (leave empty to automatically
            # utilise all available CPUs) NOTE: if this value is set higher than the
            # actual number of CPUs available, then the max number of CPUs is used.

            if (value == '') or (int(value) > self.available_cpus):
                value = self.available_cpus

        elif key == 'obs_root':
            # Define observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set observational root data directory if left undefined
            if value == '':
                # running on CTE-POWER/MN4?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3'):
                    value = '/gpfs/projects/bsc32/AC_cache/obs/ghost'
                else:
                    # running on workstation?
                    value = '/esarchive/obs/ghost'

        elif key == 'nonghost_root':
            # Define observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set observational root data directory if left undefined
            if value == '':
                # running on CTE-POWER/MN4?
                if MACHINE == 'nord3':
                    value = '/esarchive/obs'
                elif (MACHINE == 'power') or (MACHINE == 'mn4'):
                    value = None
                else:
                    # running on workstation?
                    value = '/esarchive/obs'

        elif key == 'exp_root':
            # Define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                # running on CTE-POWER?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3'):
                    return '/gpfs/projects/bsc32/AC_cache/recon/exp_interp'

                # running on workstation?
                value = '/esarchive/recon/ghost_interp_new'

        return value
