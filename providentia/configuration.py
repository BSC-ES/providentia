""" Providentia Configuration Module """

import configparser
import os
import re
import subprocess

MACHINE = os.environ.get('BSC_MACHINE', '')


def parse_path(dir, f):
    # print "Opening data file", f
    if os.path.isabs(f):
        return f
    else:
        return os.path.join(dir, f)


class ProvConfiguration(object):
    """ Configuration parameters definitions """

    def __init__(self, **kwargs):
        self.config_dir = kwargs.get('config_dir',
                                     os.path.join(os.environ['HOME'], '.providentia'))
        self.ghost_version = kwargs.get('ghost_version', '1.3.3')
        self.cartopy_data_dir = kwargs.get('cartopy_data_dir', '')
        self.available_cpus = kwargs.get('available_cpus', '')
        self.n_cpus = kwargs.get('n_cpus', '')
        self.obs_root = kwargs.get('obs_root', '')
        self.nonghost_root = kwargs.get('nonghost_root', '')
        self.exp_root = kwargs.get('exp_root', '')
        self.offline = kwargs.get('offline', '')
        self.available_networks =\
            kwargs.get('available_networks',
                       "['AERONET_v3_lev1.5','AERONET_v3_lev2.0','CANADA_NAPS','CAPMoN','CHILE_SINCA',"
                       "'EANET','EBAS','EEA_AIRBASE','EEA_AQ_eReporting','JAPAN_NIES','MEXICO_CDMX',"
                       "'MITECO','NOAA_ISD','NOAA_ISD_EU','NOAA_ISD_IP','NOAA_ISD_NA'," 
                       "'SEARCH','UK_AIR','US_EPA_AQS','US_EPA_CASTNET','US_NADP_AMNet','US_NADP_AMoN','WMO_WDCGG']")
        self.selected_species = kwargs.get('species', '')
        self.selected_network = kwargs.get('network', '')
        self.selected_matrix = kwargs.get('matrix', '')
        self.selected_resolution = kwargs.get('resolution', '')
        self.start_date = kwargs.get('start_date', '')
        self.end_date = kwargs.get('end_date', '')

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
            # set cartopy data directory (needed on CTE-POWER/MN4/N3 as has no external
            # internet connection)

            if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3') or (MACHINE == 'nord3v2'):
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
                # running on CTE-POWER/MN4/N3?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3') or (MACHINE == 'nord3v2'):
                    value = '/gpfs/projects/bsc32/AC_cache/obs/ghost'
                else:
                    # running on workstation?
                    value = '/esarchive/obs/ghost'

        elif key == 'nonghost_root':
            # Define observational root data directory (if undefined it is
            # automatically taken from the BSC machine the tool is ran on)

            # set observational root data directory if left undefined
            if value == '':
                if (MACHINE == 'power') or (MACHINE == 'nord3') or (MACHINE == 'nord3v2'):
                    value = '/esarchive/obs'
                elif MACHINE == 'mn4':
                    value = None
                else:
                    # running on workstation?
                    value = '/esarchive/obs'

        elif key == 'exp_root':
            # Define experiment root data directory
            # set experiment root data directory if left undefined
            if value == '':
                # not running on workstation?
                if (MACHINE == 'power') or (MACHINE == 'mn4') or (MACHINE == 'nord3') or (MACHINE == 'nord3v2'):
                    return '/gpfs/projects/bsc32/AC_cache/recon/exp_interp'

                # running on workstation?
                value = '/esarchive/recon/prov_interp'

        return value


def read_conf(section=None, fpath=None):
    """Read configuration"""

    config = configparser.RawConfigParser()
    config.read(fpath)
    #if no section defined, but just 1 section in file then set that as section   
    if (section is None) and (len(config.sections()) == 1):   
        section = config.sections()[0]
    #if section is undefined then cannot read
    if section is None:
        print('*** WARNING!!! CANNOT LOAD CONFIGURATION FILE AS NO SECTION DEFINED.')
        return None

    #convert numeric information appropriate types
    res = {}
    for k, val in config.items(section):
        try:
            res[k] = eval(val)
        except:
            res[k] = val
    return res


def read_offline_conf(fpath=None):
    """Read configuration files when running Providentia
    offline. When running offline, having a 'DEFAULT' section
    is mandatory. The 'DEFAULT' section contains options that
    are applied across all remaining sections."""

    config = configparser.RawConfigParser()
    config.read(fpath)

    # check if Default section has options
    if not config.defaults():
        return None

    #convert numeric information to appropriate types
    defaults = {}
    for k, val in config.items(config.default_section):
        try:
            defaults[k] = eval(val)
        except:
            defaults[k] = val

    # store remaining sections into dict
    res = {}
    for section in config.sections():
        res[section] = read_conf(section, fpath)
    return defaults, res


def write_conf(section, fpath, opts):
    """Write configurations on file. """

    config = configparser.RawConfigParser()

    # check if file exists
    if os.path.exists(fpath):
        config.read(fpath)

    # check if section exists
    if not config.has_section(section):
        config.add_section(section)

    # update configuration
    for item in opts:
        val = opts[item]
        config.set(section, item, val)

    # write configuration
    with open(fpath, 'w') as configfile:
        config.write(configfile)


def split_options(conf_string, separator="||"):
    """For the options in the configuration that define the keep and remove
    options. Returns the values in two lists, the keeps and removes"""
    keeps, removes = [], []
    if "keep:" in conf_string:
        keep_start, keep_end = conf_string.find("keep:"), conf_string.find(separator)
        keeps = conf_string[keep_start+5:keep_end]
        keeps = keeps.split(",")
        keeps = [k.strip() for k in keeps]

    if "remove:" in conf_string:
        remove_start = conf_string.find("remove:")
        removes = conf_string[remove_start+7:-1]
        removes = removes.split(",")
        removes = [r.strip() for r in removes]

    return keeps, removes
