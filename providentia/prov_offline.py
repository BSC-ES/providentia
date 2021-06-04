import os
import sys
import json
import copy

import numpy as np
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from providentia import aux

from .prov_read import DataReader
from .filter import DataFilter
from .configuration import parse_path
from .configuration import ProvConfiguration

CURRENT_PATH = os.path.abspath(os.path.dirname(__file__))


class ProvidentiaOffline(ProvConfiguration):
    """Run Providentia offline reports"""

    def __init__(self, read_type='parallel', **kwargs):
        print("starting Providentia offline")
        # portrait/landscape page figsize
        self.portrait_figsize = (8.27, 11.69)
        self.landscape_figsize = (11.69, 8.27)
        self.dpi = 200
        # super(ProvidentiaMainWindow, self).__init__()  # what was that about?
        ProvConfiguration.__init__(self, **kwargs)

        # put read_type into self
        self.read_type = read_type

        dconf_path = (os.path.join(CURRENT_PATH, 'conf/default.conf'))
        # update from config file
        if ('config' in kwargs) and ('section' in kwargs):
            self.load_conf(kwargs['section'], kwargs['config'])

        # update from command line
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in kwargs.items()})

        self.parameter_dictionary = {}
        self.standard_metadata = {}
        self.metadata_vars_to_read = {}
        self.metadata_dtype = {}
        self.standard_data_flag_name_to_data_flag_code = {}
        self.standard_QA_name_to_QA_code = {}
        self.init_standards()

        # load necessary dictionaries
        self.basic_stats_dict = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/basic_stats_dict.json')))
        self.expbias_dict = json.load(open(os.path.join(
            CURRENT_PATH, 'conf/experiment_bias_stats_dict.json')))

        self.station_subset_names = list('Barcelona')

        self.bounding_box = {'longitude': {'min': -12, 'max': 34}, 'latitude': {'min': 30, 'max': 46}}
        self.active_qa = aux.which_qa(self)
        self.active_flags = aux.which_flags(self)
        self.minimum_value, self.maximum_value = aux.which_bounds(self, self.selected_species)

        # get all netCDF monthly files per species
        species_files = os.listdir(
            '%s/%s/%s/%s/%s' % (self.obs_root, self.selected_network, self.ghost_version,
                                self.selected_resolution, self.selected_species))

        # get monthly start date (YYYYMM) of all species files
        species_files_yearmonths = \
            [int(f.split('_')[-1][:6] + '01') for f in species_files if f != 'temporary']

        # initialize structure to store all obs
        self.all_observation_data = {}
        self.all_observation_data[self.selected_network] = {}
        self.all_observation_data[self.selected_network][self.selected_resolution] = {}
        self.all_observation_data[self.selected_network][self.selected_resolution][self.selected_matrix] = {}
        self.all_observation_data[self.selected_network][
            self.selected_resolution][self.selected_matrix][self.selected_species] = species_files_yearmonths

        self.representativity_menu = init_representativity(self.selected_resolution)
        self.representativity_conf()

        # initialize DataReader
        self.datareader = DataReader(self)
        # read
        self.datareader.get_valid_obs_files_in_date_range(self.start_date, self.end_date)
        self.datareader.get_valid_experiment_files_in_date_range()

        self.datareader.read_setup(self.selected_resolution, self.start_date, self.end_date,
                                   self.selected_network, self.selected_species, self.selected_matrix)

        # create data dictionaries to fill
        self.datareader.reset_data_in_memory()
        self.metadata_inds_to_fill = np.arange(len(self.relevant_yearmonths))

        # read observations
        self.datareader.read_data('observations', self.start_date, self.end_date,
                                  self.selected_network, self.selected_resolution,
                                  self.selected_species, self.selected_matrix)
        # read selected experiments (iterate through)
        for exp in [exp.strip() for exp in self.experiments.split(",")]:
            self.datareader.read_data(exp, self.selected_start_date, self.selected_end_date, self.selected_network,
                                      self.selected_resolution, self.selected_species, self.selected_matrix)

        # update dictionary of plotting parameters (colour and zorder etc.) for each data array
        self.datareader.update_plotting_parameters()

        print(list(self.datareader.plotting_params.keys()))

        self.start_pdf()

    def init_standards(self):
        """Read and initialise GHOST standards."""
        sys.path.insert(1, '{}/GHOST_standards/{}'.
                        format(self.obs_root, self.ghost_version))
        from GHOST_standards import standard_parameters, \
            get_standard_metadata, standard_data_flag_name_to_data_flag_code, \
            standard_QA_name_to_QA_code
        # modify standard parameter dictionary to have BSC standard parameter
        # names as keys (rather than GHOST)
        for _, param_dict in standard_parameters.items():
            self.parameter_dictionary[param_dict['bsc_parameter_name']] = param_dict
        # get standard metadata dictionary
        self.standard_metadata = get_standard_metadata({'standard_units':''})
        # create list of metadata variables to read (make global)
        self.metadata_vars_to_read = [key for key in self.standard_metadata.keys()
                                      if pd.isnull(self.standard_metadata[key]['metadata_type'])
                                      == False]
        self.metadata_dtype = [(key, self.standard_metadata[key]['data_type'])
                               for key in self.metadata_vars_to_read]
        self.standard_data_flag_name_to_data_flag_code = \
            standard_data_flag_name_to_data_flag_code
        self.standard_QA_name_to_QA_code = standard_QA_name_to_QA_code
        self.qa_exceptions = ['dir10', 'spd10', 'rho2', 'acprec', 'acsnow', 'si',
                              'cldbot', 'vdist', 'ccovmean', 'cfracmean']
        self.specific_qa, self.general_qa, self.qa_diff = aux.get_qa_codes(self)

    def load_conf(self, section=None, fpath=None):
        """ Load existing configurations from file. """

        from .config import read_conf

        if fpath is None:
            fpath = parse_path(self.config_dir, self.config_file)

        if not os.path.isfile(fpath):
            print(("Error %s" % fpath))
            return

        opts = read_conf(section, fpath)
        if section is None:
            return opts

        self.opts = opts
        vars(self).update({(k, self.parse_parameter(k, val)) for k, val in opts.items()})

    def representativity_conf(self):
        """Comes here if there is a configuration loaded. Checks if there is a
        representative field loaded in the object from the conf and if there is
        assigns the value in the representativity menu"""
        for i, label in enumerate(self.representativity_menu['rangeboxes']['labels']):
            if hasattr(self, label):
                self.representativity_menu['rangeboxes']['current_lower'][i] = str(getattr(self, label))

    def start_pdf(self):
        filename = "test_pdf.pdf"

        # open new PDF file
        with PdfPages(filename) as pdf:
            self.pdf = pdf
            self.make_header()

            # define dictionary to store plot figures
            self.plot_dictionary = {}

            # make summary plots?
            if self.summary_pages:
                # setup plotting geometry of summary plots
                self.setup_summary_plots_geometry()

            # iterate through station_subset_names
            for station_subset_ind, station_subset in enumerate(self.station_subset_names):

                print('Filtering Data for {} Subset'.format(station_subset))

                # filter dataset for current station_subset
                # filter_instance = DataFilter(self, station_subset=station_subset)
                filter_instance = DataFilter(self)
                # TODO: initialize all the menus that are filtered separately in filter
                # in simple lists (?) and refactor the calls to the filtering module
                # maybe instead of calling
                # self.datareader.plotting_params = filter_instance.plotting_params ????
                # self.data_in_memory_filtered = filter_instance.data_in_memory_filtered

                print('Placing Data Arrays in Pandas Dataframes')

                self.station_subset = station_subset
                # convert filtered dataset to pandas dataframe
                self.to_pandas_dataframe()

                print('Doing Temporal Aggregation on Dataframes')

                # temporally aggregate selected data dataframes (by hour, day of week, month)
                self.pandas_temporal_aggregation()

                # if have some experiment data associated with selected stations, calculate
                # temporally aggregated basic statistic differences and bias statistics between
                # observations and experiment data arrays
                if len(self.offline_configuration.experiments) > 0:
                    print('Calculating Temporally Aggregated Bias Statisitics')
                    self.calculate_temporally_aggregated_experiment_bias_statistics()

                print('Making {} Subset Plots'.format(station_subset))

                # make summary plots?
                if self.offline_configuration.summary_pages:

                    # iterate through each summary plot to make
                    for plot_type in self.summary_plots_to_make:

                        # heatmap plot?
                        if 'heatm-' in plot_type:
                            heatmap_types = plot_type.split('-')[1:]
                            if station_subset_ind == 0:
                                heatmap_dict = {}
                                for heatmap_type in heatmap_types:
                                    heatmap_dict[heatmap_type] = {}
                                    for original_data_label in self.data_in_memory.keys():
                                        if original_data_label != 'observations':
                                            heatmap_dict[heatmap_type][
                                                self.offline_configuration.experiments[original_data_label][
                                                    'name']] = []
                            for heatmap_type in heatmap_types:
                                for original_data_label in self.data_in_memory.keys():
                                    if original_data_label != 'observations':
                                        if self.offline_configuration.temporal_colocation:
                                            data_label = '{}_colocatedto_observations'.format(original_data_label)
                                        else:
                                            data_label = original_data_label
                                        if data_label in list(self.selected_station_data.keys()):
                                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                                if heatmap_type in list(self.basic_stats_dict.keys()):
                                                    heatmap_dict[heatmap_type][
                                                        self.offline_configuration.experiments[original_data_label][
                                                            'name']].append(
                                                        self.selected_station_data[data_label]['all'][
                                                            '{}_bias'.format(heatmap_type)][0])
                                                else:
                                                    heatmap_dict[heatmap_type][
                                                        self.offline_configuration.experiments[original_data_label][
                                                            'name']].append(
                                                        self.selected_station_data[data_label]['all'][heatmap_type][
                                                            0])
                                            else:
                                                heatmap_dict[heatmap_type][
                                                    self.offline_configuration.experiments[original_data_label][
                                                        'name']].append(np.NaN)
                                        else:
                                            heatmap_dict[heatmap_type][
                                                self.offline_configuration.experiments[original_data_label][
                                                    'name']].append(np.NaN)

                            if station_subset_ind == (len(self.station_subset_names) - 1):
                                for heatmap_type_ii, heatmap_type in enumerate(heatmap_types):
                                    relevant_axis = self.get_relevant_axis('heatm', heatmap_type_ii)
                                    self.characteristics_per_plot_type['heatm']['axis_title'][
                                        'label'] = heatmap_type
                                    relevant_axis.set_title(
                                        **self.characteristics_per_plot_type['heatm']['axis_title'])
                                    if heatmap_dict:
                                        heatmap_df = pd.DataFrame(data=heatmap_dict[heatmap_type],
                                                                  index=self.station_subset_names)
                                        self.make_heatmap(relevant_axis, heatmap_type, heatmap_df)
                        # other plot type?
                        else:
                            # count how many plots are made per plot type
                            current_plot_ind = 0

                            # iterate through all data arrays
                            original_data_array_labels = list(self.data_in_memory.keys())
                            for original_data_label in original_data_array_labels:

                                # get zstat for plot type (if exists)
                                if '-' in plot_type:
                                    zstat = plot_type.split('-')[1]
                                    # get base name name of zstat (dropping _bias suffix)
                                    base_zstat = zstat.split('_bias')[0]
                                    # get zstat type (basic or expbias)
                                    z_statistic_type = self.get_z_statistic_type(base_zstat)
                                    # get zstat sign (absolute or bias)
                                    z_statistic_sign = self.get_z_statistic_sign(zstat, z_statistic_type)

                                # map plots (1 plot per data array/s (1 array if absolute plot,
                                # 2 arrays if making bias plot), per subset)
                                if 'map-' in plot_type:

                                    # get necessary data arrays
                                    if '-obs' in plot_type:
                                        if original_data_label != 'observations':
                                            continue
                                        if self.offline_configuration.temporal_colocation:
                                            z1 = 'observations_colocatedto_experiments'
                                        else:
                                            z1 = 'observations'
                                        z2 = ''
                                    elif z_statistic_sign == 'bias':
                                        if original_data_label == 'observations':
                                            continue
                                        if self.offline_configuration.temporal_colocation:
                                            z1 = 'observations_colocatedto_{}'.format(original_data_label)
                                            z2 = '{}_colocatedto_observations'.format(original_data_label)
                                        else:
                                            z1 = 'observations'
                                            z2 = original_data_label
                                    else:
                                        if original_data_label == 'observations':
                                            if self.offline_configuration.temporal_colocation:
                                                z1 = 'observations_colocatedto_experiments'
                                            else:
                                                z1 = 'observations'
                                        else:
                                            if self.offline_configuration.temporal_colocation:
                                                z1 = '{}_colocatedto_observations'.format(original_data_label)
                                            else:
                                                z1 = original_data_label
                                        z2 = ''

                                    # get relevant axis to plot on
                                    relevant_axis = self.get_relevant_axis(plot_type, (current_plot_ind * len(
                                        self.station_subset_names)) + station_subset_ind)

                                    # make map plot
                                    n_stations = self.make_map(relevant_axis, z1, z2, zstat)

                                    # set axis title
                                    if plot_type == 'map-Data %-obs':
                                        self.characteristics_per_plot_type[plot_type]['axis_title'][
                                            'label'] = '{} {}\n({} Stations)'.format(original_data_label,
                                                                                     station_subset, n_stations)
                                    else:
                                        self.characteristics_per_plot_type[plot_type]['axis_title'][
                                            'label'] = '{}\n{}'.format(original_data_label, station_subset)
                                    relevant_axis.set_title(
                                        **self.characteristics_per_plot_type[plot_type]['axis_title'])

                                # other plots (1 plot per subset, with multiple data arrays)
                                else:
                                    # get relevant axis to plot on
                                    relevant_axis = self.get_relevant_axis(plot_type, station_subset_ind)

                                    # get necessary data array to plot
                                    if original_data_label == 'observations':
                                        if self.offline_configuration.temporal_colocation:
                                            data_label = 'observations_colocatedto_experiments'
                                        else:
                                            data_label = 'observations'
                                    else:
                                        if self.offline_configuration.temporal_colocation:
                                            data_label = '{}_colocatedto_observations'.format(original_data_label)
                                        else:
                                            data_label = original_data_label

                                    # periodic plots
                                    if 'periodic-' in plot_type:
                                        # skip observational array if plotting bias stat
                                        if (z_statistic_sign == 'bias') & (original_data_label == 'observations'):
                                            continue
                                        func = getattr(self, 'make_periodic')
                                        if data_label in list(self.selected_station_data.keys()):
                                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                                func(relevant_axis, data_label, zstat)
                                        # set axis title
                                        self.characteristics_per_plot_type[plot_type]['axis_title'][
                                            'label'] = station_subset
                                        relevant_axis['hour'].set_title(
                                            **self.characteristics_per_plot_type[plot_type]['axis_title'])
                                    # other plot types (distribution, timeseries etc.)
                                    else:
                                        # determine if plotting bias stat
                                        bias_stat = False
                                        plot_type_split = plot_type.split('_')
                                        if len(plot_type_split) > 1:
                                            bias_stat = True
                                        # skip observational array if plotting bias stat
                                        if (bias_stat) & (original_data_label == 'observations'):
                                            continue
                                        # make plot
                                        func = getattr(self, 'make_{}'.format(plot_type_split[0]))

                                        if data_label in list(self.selected_station_data.keys()):
                                            if len(self.selected_station_data[data_label]['pandas_df']['data']) > 0:
                                                # bias plot or standard?
                                                if bias_stat:
                                                    func(relevant_axis, data_label, bias=True)
                                                else:
                                                    func(relevant_axis, data_label)
                                        # set axis title
                                        self.characteristics_per_plot_type[plot_type_split[0]]['axis_title'][
                                            'label'] = station_subset
                                        relevant_axis.set_title(
                                            **self.characteristics_per_plot_type[plot_type_split[0]]['axis_title'])

                                # iterate number of plots have made for current type of plot
                                current_plot_ind += 1

                # make individual station pages?

            # generate colourbars
            if self.offline_configuration.summary_pages:
                relevant_axs = []
                cb_axs = []
                previous_plot_type = ''
                for figure_n in self.plot_dictionary.keys():
                    plot_type = self.plot_dictionary[figure_n]['plot_type']
                    if (previous_plot_type != plot_type) & (len(relevant_axs) > 0):
                        self.generate_colourbar(cb_axs, relevant_axs, previous_plot_type)
                        relevant_axs = []
                        cb_axs = []
                    if 'cb' in list(self.characteristics_per_plot_type[plot_type].keys()):
                        relevant_axs.extend(self.plot_dictionary[figure_n]['axs'])
                        cb_axs.append(self.plot_dictionary[figure_n]['cb_ax'])
                    previous_plot_type = copy.deepcopy(plot_type)
                if len(relevant_axs) > 0:
                    self.generate_colourbar(cb_axs, relevant_axs, previous_plot_type)

                #
                for plot_type in self.summary_plots_to_make:
                    if ('timeseries' not in plot_type) & ('map-' not in plot_type) & ('heatm-' not in plot_type):
                        # if ('map-' not in plot_type) & ('heatm-' not in plot_type):
                        # if 'heatm-' not in plot_type:
                        self.harmonise_xy_lims(plot_type)

            # save page figures
            print('WRITING PDF')
            for figure_n in self.plot_dictionary.keys():
                fig = self.plot_dictionary[figure_n]['fig']
                self.pdf.savefig(fig, dpi=self.dpi)
                plt.close(fig)

    def make_header(self):
        # set tile
        page = plt.figure(figsize=self.portrait_figsize)
        if hasattr(self, 'report_title'):
            txt = self.report_title
        else:
            txt = 'An example report'
        # TODO: define default title in case no title in conf
        page.text(0.5, 0.9, txt, transform=page.transFigure,
                  size=20, ha="center", va='top', wrap=True)

        experiment_labels = [exp.strip() for exp in self.experiments.split(",")]
        txt = 'Network = {}\nTemporal Resolution = {}\n' \
              'Species = {}\nDate Range = {} - {}\nExperiments = {}\n'\
            .format(self.selected_network,
                    self.selected_resolution,
                    self.selected_species,
                    self.start_date,
                    self.end_date, experiment_labels)
        # txt += 'Station Subset/s = {}\n'.format(self.station_subset_names)
        if hasattr(self, 'bounding_box'):
            txt += 'Bounding Box = [{}E:{}E, {}N:{}N]\n'.\
                format(self.bounding_box['longitude']['min'],
                       self.bounding_box['longitude']['max'],
                       self.bounding_box['latitude']['min'],
                       self.bounding_box['latitude']['max'])

        """This statement cannot be applied. We do not
        have a data_availability_filter field nor group"""
        # if hasattr(self, 'data_availability_filter'):
        #     txt += 'Data Availability = {}'.format(self.data_availability_filter)
        # page.text(0.5, 0.82, txt, transform=page.transFigure,
        #           weight='light', size=15, ha="center", va='top', wrap=True)

        self.pdf.savefig(page, dpi=self.dpi)
        plt.close(page)

    def setup_summary_plots_geometry(self):
        """setup summary plots geometry"""

        # define projections for map plot and actual geographic coordinates
        self.datacrs = ccrs.PlateCarree()
        self.plotcrs = ccrs.Robinson()

        land_polygon_resolutions = {'low': '110m',
                                    'medium': '50m',
                                    'high': '10m'}

        feature = cfeature.NaturalEarthFeature('physical', 'land',
                                               land_polygon_resolutions[
                                                   self.map_coastline_resolution],
                                               facecolor='0.85')

        # define plot types per type of report
        # TODO: to what depend on these plot? what is the name-convention? should they be configurable?
        # TODO: I think it should be moved to a configuration file. Not everyone will want to make cams50 reports
        plot_types_per_report_type = {
            'cams50': ['map-Data %-obs', 'map-p5', 'map-p50',
                       'map-p95', 'map-p5_bias', 'map-p50_bias',
                       'map-p95_bias', 'map-MB', 'map-RMSE', 'map-MFB',
                       'map-MAFB', 'map-r', 'timeseries', 'distribution',
                       'distribution_bias', 'periodic-p5', 'periodic-p50',
                       'periodic-p95', 'periodic-Min', 'periodic-Max',
                       'periodic-p5_bias', 'periodic-p50_bias',
                       'periodic-p95_bias', 'periodic-Min_bias',
                       'periodic-Max_bias', 'periodic-MB', 'periodic-RMSE',
                       'periodic-MFB', 'periodic-MAFB', 'periodic-r',
                       'heatm-p5-p50-p95-Min-Max-MB-RMSE-MFB-MAFB-r']}

        # define characteristics per plot type
        self.characteristics_per_plot_type = {
            'timeseries': {'pages': [], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 3},
                           'xtick_share': True, 'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                           'page_title': {'t': 'Time Series', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                           'axis_title': {'label': '', 'fontsize': 8}, 'axis_xlabel': {'xlabel': 'Time', 'fontsize': 8},
                           'axis_ylabel': {'ylabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                           'xticks': {'labelsize': 7, 'rotation': 0}, 'yticks': {'labelsize': 7},
                           'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                           'subplots_adjust': {'top': 0.90, 'bottom': 0.08}},
            'distribution': {'pages': [], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                             'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                             'page_title': {'t': 'Distribution', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                             'axis_title': {'label': '', 'fontsize': 10},
                             'axis_xlabel': {'xlabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                             'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                             'yticks': {'labelsize': 7}, 'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0},
                             'tightlayout': True, 'subplots_adjust': {'top': 0.90}},
            'distribution_bias': {'pages': [], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                                  'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                                  'page_title': {'t': 'Distributional bias', 'fontsize': 18, 'ha': 'left', 'x': 0.05,
                                                 'y': 0.98}, 'axis_title': {'label': '', 'fontsize': 10},
                                  'axis_xlabel': {'xlabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                                  'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                                  'yticks': {'labelsize': 7},
                                  'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                                  'subplots_adjust': {'top': 0.90}},

            'heatm': {'pages': [], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 1},
                      'page_title': {'t': 'Statistical Heatmap', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                      'axis_title': {'label': '', 'fontsize': 8}, 'xticks': {'labelsize': 7, 'rotation': -270},
                      'yticks': {'labelsize': 7, 'rotation': -315}, 'tightlayout': True,
                      'subplots_adjust': {'top': 0.90}, 'cb_xlabel': {'xlabel': '', 'fontsize': 8},
                      'cb_xticks': {'labelsize': 8}, 'annot': True}
        }
        # add all types of map plot types
        map_basicstat_obs_template = {'pages': [], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                                                              'subplot_kw': {'projection': self.plotcrs}},
                                      'page_title': {'t': '', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                                      'axis_title': {'label': '', 'fontsize': 8},
                                      'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                                      'cb': {'position': [0.5, 0.95, 0.4, 0.04]},
                                      'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}}
        map_basicstat_template = {'pages': [], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                                                          'subplot_kw': {'projection': self.plotcrs}},
                                  'page_title': {'t': '', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                                  'axis_title': {'label': '', 'fontsize': 8},
                                  'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                                  'cb': {'position': [0.5, 0.95, 0.4, 0.04]},
                                  'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}}
        map_biasstat_template = {'pages': [], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                                                         'subplot_kw': {'projection': self.plotcrs}},
                                 'page_title': {'t': '', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                                 'axis_title': {'label': '', 'fontsize': 8},
                                 'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                                 'cb': {'position': [0.5, 0.95, 0.4, 0.04]}, 'cb_xlabel': {'xlabel': '', 'fontsize': 8},
                                 'cb_xticks': {'labelsize': 8}}
        periodic_basicstat_template = {'pages': [],
                                       'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                                       'page_title': {'t': '', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                                       'axis_title': {'label': '', 'fontsize': 8},
                                       'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                                       'yticks': {'labelsize': 7},
                                       'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                                       'subplots_adjust': {'top': 0.85}}
        periodic_biasstat_template = {'pages': [],
                                      'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                                      'page_title': {'t': '', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                                      'axis_title': {'label': '', 'fontsize': 8},
                                      'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                                      'yticks': {'labelsize': 7},
                                      'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                                      'subplots_adjust': {'top': 0.85}}

        for bstat in self.basic_stats_dict.keys():
            self.characteristics_per_plot_type['map-{}-obs'.format(bstat)] = copy.deepcopy(map_basicstat_obs_template)
            self.characteristics_per_plot_type['map-{}'.format(bstat)] = copy.deepcopy(map_basicstat_template)
            self.characteristics_per_plot_type['map-{}_bias'.format(bstat)] = copy.deepcopy(map_biasstat_template)
            self.characteristics_per_plot_type['periodic-{}'.format(bstat)] = copy.deepcopy(periodic_basicstat_template)
            self.characteristics_per_plot_type['periodic-{}_bias'.format(bstat)] = copy.deepcopy(
                periodic_biasstat_template)
            self.characteristics_per_plot_type['map-{}-obs'.format(bstat)]['page_title']['t'] = \
            self.basic_stats_dict[bstat]['label']
            self.characteristics_per_plot_type['map-{}'.format(bstat)]['page_title']['t'] = \
            self.basic_stats_dict[bstat]['label']
            self.characteristics_per_plot_type['map-{}_bias'.format(bstat)]['page_title']['t'] = '{} bias'.format(
                self.basic_stats_dict[bstat]['label'])
            self.characteristics_per_plot_type['periodic-{}'.format(bstat)]['page_title']['t'] = \
            self.basic_stats_dict[bstat]['label']
            self.characteristics_per_plot_type['periodic-{}_bias'.format(bstat)]['page_title']['t'] = '{} bias'.format(
                self.basic_stats_dict[bstat]['label'])

        for estat in self.expbias_dict.keys():
            self.characteristics_per_plot_type['map-{}'.format(estat)] = copy.deepcopy(map_biasstat_template)
            self.characteristics_per_plot_type['periodic-{}'.format(estat)] = copy.deepcopy(periodic_biasstat_template)
            self.characteristics_per_plot_type['map-{}'.format(estat)]['page_title']['t'] = self.expbias_dict[estat][
                'label']
            self.characteristics_per_plot_type['periodic-{}'.format(estat)]['page_title']['t'] = \
            self.expbias_dict[estat]['label']

        # iterate through needed plot types, creating
        self.summary_plots_to_make = plot_types_per_report_type[self.report_type]
        current_page_n = 1

        # iterate through plot types to make
        for plot_type in self.summary_plots_to_make:
            if 'heatm-' in plot_type:
                n_plots_per_plot_type = len(plot_type.split('-')) - 1
                plot_type = 'heatm'
            plot_characteristics = self.characteristics_per_plot_type[plot_type]
            plot_characteristics_vars = list(plot_characteristics.keys())

            if '-' in plot_type:
                zstat = plot_type.split('-')[1]
                # get base name name of zstat (dropping _bias suffix)
                base_zstat = zstat.split('_bias')[0]
                # get zstat type (basic or expbias)
                z_statistic_type = get_z_statistic_type(self.basic_stats_dict, base_zstat)
                # get zstat sign (absolute or bias)
                z_statistic_sign = get_z_statistic_sign(zstat, z_statistic_type)

            if 'map-' in plot_type:
                if '-obs' in plot_type:
                    n_plots_per_plot_type = len(self.station_subset_names)
                elif z_statistic_sign == 'bias':
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            (len(list(self.datareader.data_in_memory.keys())) - 1)
                else:
                    n_plots_per_plot_type = len(self.station_subset_names) * \
                                            len(list(self.datareader.data_in_memory.keys()))
            elif 'periodic-' in plot_type:
                if z_statistic_type == 'basic':
                    if base_zstat not in ['Data %', 'Exceedances']:
                        plot_characteristics['axis_ylabel']['ylabel'] = self.datareader.measurement_units  # 'µg m⁻³'
                    else:
                        plot_characteristics['axis_ylabel']['ylabel'] = self.basic_stats_dict[base_zstat]['label']
                else:
                    plot_characteristics['axis_ylabel']['ylabel'] = self.expbias_dict[base_zstat]['label']
            elif plot_type != 'heatm':
                n_plots_per_plot_type = len(self.station_subset_names)

            # get n pages per plot type
            n_pages_per_plot_type = int(np.ceil(n_plots_per_plot_type / (
                        plot_characteristics['figure']['ncols'] * plot_characteristics['figure']['nrows'])))
            plot_ii_per_type = 0
            # iterate through n pages per plot
            for page_n in range(n_pages_per_plot_type):
                fig, axs = plt.subplots(**plot_characteristics['figure'])
                self.plot_dictionary[current_page_n] = {'fig': fig, 'axs': [], 'plot_type': plot_type}

                # make page title?
                if 'page_title' in plot_characteristics_vars:
                    st = fig.suptitle(**plot_characteristics['page_title'])

                # iterate through axes (by row, then column)
                row_ii = -1
                col_ii = copy.deepcopy(plot_characteristics['figure']['ncols'])
                for ax in axs.flatten():
                    if col_ii == plot_characteristics['figure']['ncols']:
                        row_ii += 1
                        col_ii = 0
                    if row_ii == (plot_characteristics['figure']['nrows'] - 1):
                        last_row_on_page = True
                    else:
                        last_row_on_page = False

                    ax.set_rasterized(True)
                    if plot_ii_per_type < n_plots_per_plot_type:

                        # determine if are on last valid row to plot
                        if (n_plots_per_plot_type - plot_ii_per_type) <= plot_characteristics['figure']['ncols']:
                            last_valid_row = True
                        else:
                            last_valid_row = False

                            # setup periodic plot type gridspec
                        if 'periodic-' in plot_type:
                            gs = gridspec.GridSpecFromSubplotSpec(20, 20, subplot_spec=ax)
                            grid_dict = dict()
                            grid_dict['hour'] = fig.add_subplot(gs[:9, :])
                            grid_dict['month'] = fig.add_subplot(gs[11:, :11])
                            grid_dict['dayofweek'] = fig.add_subplot(gs[11:, 13:])
                            self.plot_dictionary[current_page_n]['axs'].append(grid_dict)
                            ax.spines['top'].set_color('none')
                            ax.spines['bottom'].set_color('none')
                            ax.spines['left'].set_color('none')
                            ax.spines['right'].set_color('none')
                            ax.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
                            ax.set_ylabel(**plot_characteristics['axis_ylabel'])
                        else:
                            self.plot_dictionary[current_page_n]['axs'].append(ax)

                            # make axis title?
                        if 'axis_title' in plot_characteristics_vars:
                            ax.set_title(**plot_characteristics['axis_title'])

                        # make axis xlabel (only on last row on page/last valid row of visible axes)?
                        if 'axis_xlabel' in plot_characteristics_vars:
                            if last_valid_row or last_row_on_page:
                                ax.set_xlabel(**plot_characteristics['axis_xlabel'])

                        # make axis ylabel (only on leftmost column of visible axes)?
                        if ('axis_ylabel' in plot_characteristics_vars) & (col_ii == 0):
                            ax.set_ylabel(**plot_characteristics['axis_ylabel'])

                        # if are sharing xticks, and not on last row on page/last
                        # valid row, then ensure current axis xticks are hidden
                        if ('xtick_share' in plot_characteristics_vars) & (not last_valid_row) & (not last_row_on_page):
                            plt.setp(ax.get_xticklabels(), visible=False)

                        # if are sharing yticks, and not on left column, then ensure current axis yticks are hidden
                        if ('ytick_share' in plot_characteristics_vars) & (col_ii != 0):
                            plt.setp(ax.get_yticklabels(), visible=False)

                        # add gridlines?
                        if 'grid' in plot_characteristics_vars:
                            ax.grid(**plot_characteristics['grid'])

                        # add coastlines/gridlines to map?
                        if 'map-' in plot_type:
                            ax.add_feature(feature)
                            ax.gridlines(linestyle='-', alpha=0.4)
                            if hasattr(self, 'bounding_box'):
                                if 'longitude' in self.bounding_box.keys():
                                    if 'min' in self.bounding_box['longitude']:
                                        min_lon, null = self.plotcrs.transform_point(
                                            self.bounding_box['longitude']['min'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(left=min_lon)
                                    if 'max' in self.bounding_box['longitude']:
                                        max_lon, null = self.plotcrs.transform_point(
                                            self.bounding_box['longitude']['max'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(right=max_lon)
                                if 'latitude' in self.bounding_box.keys():
                                    if 'min' in self.bounding_box['latitude']:
                                        null, min_lat = self.plotcrs.transform_point(0, self.bounding_box[
                                            'latitude']['min'], src_crs=self.datacrs)
                                        ax.set_ylim(bottom=min_lat)
                                    if 'max' in self.bounding_box['latitude']:
                                        null, max_lat = self.plotcrs.transform_point(0, self.bounding_box[
                                            'latitude']['max'], src_crs=self.datacrs)
                                        ax.set_ylim(top=max_lat)

                    else:
                        ax.set_visible(False)
                    plot_ii_per_type += 1
                    col_ii += 1

                    # tight layout?
                if 'tightlayout' in plot_characteristics_vars:
                    fig.tight_layout()
                # adjust subplots?
                if 'subplots_adjust' in plot_characteristics_vars:
                    fig.subplots_adjust(**plot_characteristics['subplots_adjust'])

                    # make legend?
                if 'legend' in plot_characteristics_vars:
                    leg_dict = plot_characteristics['legend']
                    leg_dict['handles'] = self.make_legend_elements()
                    fig.legend(**leg_dict)

                # add colourbar axis to plot dictionary?
                if 'cb' in plot_characteristics_vars:
                    self.plot_dictionary[current_page_n]['cb_ax'] = fig.add_axes(plot_characteristics['cb']['position'])
                    self.plot_dictionary[current_page_n]['cb_ax'].set_rasterized(True)

                # add current page number
                self.characteristics_per_plot_type[plot_type]['pages'].append(current_page_n)
                current_page_n += 1

    def make_legend_elements(self):
        """Function that makes legend elements"""

        # create legend elements
        # add observations element
        legend_elements = [Line2D([0], [0], marker='o', color='white',
                                  markerfacecolor=self.datareader.plotting_params['observations']['colour'],
                                  markersize=self.legend_markersize, label='observations')]
        # add element for each experiment
        for experiment_ind, experiment in enumerate(sorted(list(self.datareader.data_in_memory.keys()))):
            if experiment != 'observations':
                # add experiment element
                legend_elements.append(Line2D([0], [0], marker='o', color='white',
                                              markerfacecolor=self.datareader.plotting_params[experiment]['colour'],
                                              markersize=self.legend_markersize,
                                              label=experiment))
        return legend_elements


def get_z_statistic_type(stats_dict, zstat):
    """Function that checks if the z statistic is basic or expbias statistic"""

    # check if the chosen statistic is a basic statistic
    if zstat in stats_dict.keys():
        return 'basic'
    # if not a basic statistic, it must be an experiment bias statistic
    else:
        return 'expbias'


def get_z_statistic_sign(zstat, zstat_type):
    """Function that checks if the z statistic is an absolute or bias statistic"""

    # statistic is bias?
    if ('_bias' in zstat) or (zstat_type == 'expbias'):
        return 'bias'
    # statistic is bias?
    else:
        return 'absolute'


def init_representativity(resolution):
    """Initialize representativity structure, similar to the one used
    in Providentia dashboard, keeping only the necessary fields"""

    # initialize representativity menu, inly with necessary substructures
    representativity = {'rangeboxes': {'labels': [], 'current_lower': []}}

    if (resolution == 'hourly') or (resolution == 'hourly_instantaneous'):
        representativity['rangeboxes']['labels'] = ['hourly_native_representativity_percent',
                                                    'hourly_native_max_gap_percent',
                                                    'daily_native_representativity_percent',
                                                    'daily_representativity_percent',
                                                    'daily_native_max_gap_percent',
                                                    'daily_max_gap_percent',
                                                    'monthly_native_representativity_percent',
                                                    'monthly_representativity_percent',
                                                    'monthly_native_max_gap_percent',
                                                    'monthly_max_gap_percent',
                                                    'all_representativity_percent', 'all_max_gap_percent']
    # daily temporal resolution?
    elif (resolution == 'daily') or (resolution == '3hourly') or \
            (resolution == '6hourly') or (resolution == '3hourly_instantaneous') or \
            (resolution == '6hourly_instantaneous'):
        representativity['rangeboxes']['labels'] = ['daily_native_representativity_percent',
                                                    'daily_native_max_gap_percent',
                                                    'monthly_native_representativity_percent',
                                                    'monthly_representativity_percent',
                                                    'monthly_native_max_gap_percent',
                                                    'monthly_max_gap_percent',
                                                    'all_representativity_percent', 'all_max_gap_percent']
    # monthly temporal resolution?
    elif resolution == 'monthly':
        representativity['rangeboxes']['labels'] = ['monthly_native_representativity_percent',
                                                    'monthly_native_max_gap_percent',
                                                    'all_representativity_percent', 'all_max_gap_percent']

    # initialise rangebox values --> for data representativity fields
    # the default is 0%, for max gap fields % the default is 100%
    representativity['rangeboxes']['current_lower'] = []
    for label_ii, label in enumerate(representativity['rangeboxes']['labels']):
        if 'max_gap' in label:
            representativity['rangeboxes']['current_lower'].append('100')
        else:
            representativity['rangeboxes']['current_lower'].append('0')

    return representativity


def main_offline(**kwargs):
    """Main function when running offine reports"""
    ProvidentiaOffline(**kwargs)
