import os
import sys
import json
import copy

import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages

from .prov_read import DataReader
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
        # initialize DataReader
        self.datareader = DataReader(self)

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
        self.get_qa_codes()

    def get_qa_codes(self):
        """Retrieve QA codes from GHOST_standards using the qa flags' names.

        Specific flags are defined for the following species:
        ['WND_DIR_10M','WND_SPD_10M','RH_2M','PREC_ACCUM','SNOW_ACCUM',
        'SNOW_DEPTH','CEILING_HEIGHT','VIS_DIST','CLOUD_CVG','CLOUD_CVG_FRAC']"""

        # get names from json files
        specific_qa_names = json.load(open(
            "providentia/conf/default_flags.json"))['specific_qa']
        general_qa_names = json.load(open(
            "providentia/conf/default_flags.json"))['general_qa']
        # get codes
        self.specific_qa = [self.standard_QA_name_to_QA_code[qa_name]
                            for qa_name in specific_qa_names]
        self.general_qa = [self.standard_QA_name_to_QA_code[qa_name]
                           for qa_name in general_qa_names]
        # get difference of flags, needed later for updating default selection
        self.qa_diff = list(set(self.general_qa) - set(self.specific_qa))

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
                           'axis_ylabel': {'ylabel': '{}'.format(self.measurement_units), 'fontsize': 8},
                           'xticks': {'labelsize': 7, 'rotation': 0}, 'yticks': {'labelsize': 7},
                           'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                           'subplots_adjust': {'top': 0.90, 'bottom': 0.08}},
            # 'timeseries':  {'pages':[], 'figure':{'figsize':self.landscape_figsize,  'ncols':2, 'nrows':3}, 'xtick_share':True, 'grid':{'axis':'both', 'color':'lightgrey', 'alpha':0.8}, 'page_title':{'t':'', 'fontsize':18, 'ha':'left', 'x':0.05, 'y':0.98},         'axis_title':{'label':'', 'fontsize':8},  'axis_xlabel':{'xlabel':'Time','fontsize':8},                              'axis_ylabel':{'ylabel':'µg m⁻³','fontsize':8}, 'xticks':{'labelsize':7, 'rotation':0}, 'yticks':{'labelsize':7},                  'legend':{'loc':'upper right', 'ncol':3, 'fontsize':8.0}, 'tightlayout':True, 'subplots_adjust':{'top':0.90,'bottom':0.08},                                                                        'trend':{'n_points':24, 'min_points':6}},
            'distribution': {'pages': [], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                             'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                             'page_title': {'t': 'Distribution', 'fontsize': 18, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                             'axis_title': {'label': '', 'fontsize': 10},
                             'axis_xlabel': {'xlabel': '{}'.format(self.measurement_units), 'fontsize': 8},
                             'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                             'yticks': {'labelsize': 7}, 'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0},
                             'tightlayout': True, 'subplots_adjust': {'top': 0.90}},
            'distribution_bias': {'pages': [], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                                  'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                                  'page_title': {'t': 'Distributional bias', 'fontsize': 18, 'ha': 'left', 'x': 0.05,
                                                 'y': 0.98}, 'axis_title': {'label': '', 'fontsize': 10},
                                  'axis_xlabel': {'xlabel': '{}'.format(self.measurement_units), 'fontsize': 8},
                                  'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                                  'yticks': {'labelsize': 7},
                                  'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                                  'subplots_adjust': {'top': 0.90}},
            # 'distribution':{'pages':[], 'figure':{'figsize':self.portrait_figsize,  'ncols':2, 'nrows':4},                     'grid':{'axis':'both', 'color':'lightgrey', 'alpha':0.8}, 'page_title':{'t':'Distribution', 'fontsize':18, 'ha':'left', 'x':0.05, 'y':0.98},        'axis_title':{'label':'', 'fontsize':10}, 'axis_xlabel':{'xlabel':'µg m⁻³','fontsize':8}, 'axis_ylabel':{'ylabel':'Density','fontsize':8},                           'xticks':{'labelsize':7},                  'yticks':{'labelsize':7},                  'legend':{'loc':'upper right', 'ncol':3, 'fontsize':8.0}, 'tightlayout':True, 'subplots_adjust':{'top':0.90}},
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
        self.summary_plots_to_make = plot_types_per_report_type[self.offline_configuration.report_type]
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
                z_statistic_type = self.get_z_statistic_type(base_zstat)
                # get zstat sign (absolute or bias)
                z_statistic_sign = self.get_z_statistic_sign(zstat, z_statistic_type)

            if 'map-' in plot_type:
                if '-obs' in plot_type:
                    n_plots_per_plot_type = len(self.station_subset_names)
                elif z_statistic_sign == 'bias':
                    n_plots_per_plot_type = len(self.station_subset_names) * (len(list(self.data_in_memory.keys())) - 1)
                else:
                    n_plots_per_plot_type = len(self.station_subset_names) * len(list(self.data_in_memory.keys()))
            elif 'periodic-' in plot_type:
                if z_statistic_type == 'basic':
                    if base_zstat not in ['Data %', 'Exceedances']:
                        plot_characteristics['axis_ylabel']['ylabel'] = self.measurement_units  # 'µg m⁻³'
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
                            grid_dict = {}
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
                            if (last_valid_row) or (last_row_on_page):
                                ax.set_xlabel(**plot_characteristics['axis_xlabel'])

                        # make axis ylabel (only on leftmost column of visible axes)?
                        if ('axis_ylabel' in plot_characteristics_vars) & (col_ii == 0):
                            ax.set_ylabel(**plot_characteristics['axis_ylabel'])

                        # if are sharing xticks, and not on last row on page/last valid row, then ensure current axis xticks are hidden
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
                            if hasattr(self.offline_configuration, 'bounding_box'):
                                if 'longitude' in self.offline_configuration.bounding_box.keys():
                                    if 'min' in self.offline_configuration.bounding_box['longitude']:
                                        min_lon, null = self.plotcrs.transform_point(
                                            self.offline_configuration.bounding_box['longitude']['min'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(left=min_lon)
                                    if 'max' in self.offline_configuration.bounding_box['longitude']:
                                        max_lon, null = self.plotcrs.transform_point(
                                            self.offline_configuration.bounding_box['longitude']['max'], 0,
                                            src_crs=self.datacrs)
                                        ax.set_xlim(right=max_lon)
                                if 'latitude' in self.offline_configuration.bounding_box.keys():
                                    if 'min' in self.offline_configuration.bounding_box['latitude']:
                                        null, min_lat = self.plotcrs.transform_point(0,
                                                                                     self.offline_configuration.bounding_box[
                                                                                         'latitude']['min'],
                                                                                     src_crs=self.datacrs)
                                        ax.set_ylim(bottom=min_lat)
                                    if 'max' in self.offline_configuration.bounding_box['latitude']:
                                        null, max_lat = self.plotcrs.transform_point(0,
                                                                                     self.offline_configuration.bounding_box[
                                                                                         'latitude']['max'],
                                                                                     src_crs=self.datacrs)
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


def main_offline(**kwargs):
    """Main function when running offine reports"""
    ProvidentiaOffline(**kwargs)
