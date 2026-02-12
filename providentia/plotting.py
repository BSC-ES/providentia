""" Class for plotting """

import copy
from datetime import datetime
from itertools import groupby
import math
import sys

import cartopy
import cartopy.crs as ccrs
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.projections import PolarAxes
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.floating_axes as fa
import numpy as np
from packaging.version import Version
import pandas as pd
from PIL import Image
import pyproj
import seaborn as sns
import yaml
import xarray as xr
import xskillscore as xs

from providentia.auxiliar import CURRENT_PATH, join, get_conversion_factor, get_standard_parameters_by_speci
from .calculate import ModBias
from .statistics import (boxplot_inner_fences, calculate_statistic, group_periodic,
                         get_fairmode_data, get_z_statistic_info, get_z_statistic_type)
from .read_aux import drop_nans, get_valid_metadata
from .plot_aux import (create_statistical_timeseries, get_multispecies_aliases, get_AERONET_sizedist_bin_radius,
                       get_taylor_diagram_ghelper_info, kde_fft, merge_cells, periodic_labels, 
                       periodic_xticks, round_decimal_places, temp_axis_dict, get_fairmode_RV_exceendance)
from .plot_formatting import set_axis_title
from .warnings_prv import show_message

# speed up transformations in cartopy
pyproj.set_use_global_context()

PROVIDENTIA_ROOT = '/'.join(CURRENT_PATH.split('/')[:-1])
fairmode_settings = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/fairmode.yaml')))
contingency_settings = yaml.safe_load(open(join(PROVIDENTIA_ROOT, 'settings/contingency.yaml')))


class Plotting:
    """ Class for plotting """

    def __init__(self, read_instance=None, canvas_instance=None):
        """
        Initialise the Plotting instance.

        Parameters
        ----------
        read_instance : object, optional
            The data management instance containing configuration and raw data.
        canvas_instance : object, optional
            The Matplotlib/Cartopy canvas where visualisations are rendered.
        """
        
        self.read_instance = read_instance
        self.canvas_instance = canvas_instance

        # set cartopy data directory
        cartopy.config['pre_existing_data_dir'] = self.read_instance.cartopy_data_dir

        # define projections for map plot and actual geographic coordinates
        self.canvas_instance.datacrs = ccrs.PlateCarree()
        proj_class = getattr(ccrs, self.canvas_instance.plot_characteristics_templates['map']['projection']) 
        self.canvas_instance.plotcrs = proj_class()
        
        # set miscellaneous vars
        self.canvas_instance.temporal_axis_mapping_dict = temp_axis_dict()
        self.canvas_instance.periodic_xticks = periodic_xticks()
        self.canvas_instance.periodic_labels = periodic_labels()

    def set_plot_characteristics(self, plot_types, zstat=False, data_labels=None, format={}):
        """
        Validates requested plot types against available data and configures their visual attributes.

        Parameters
        ----------
        plot_types : list
            Plot types to create.
        zstat : str, optional
            Statistic to be plotted.
        data_labels : list, optional
            Data arrays to plot.
        format : dict, optional
            Dictionary to overwrite default formatting.

        Returns
        -------
        valid_plot_type : bool
            Returns a boolean indicating plot validity when operating in library mode.
        """

        # if data_labels are not defined, take all in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

        # add all valid defined plots to plot_characteristics
        for plot_type in plot_types:
            
            # initialize condition
            valid_plot_type = True

            # do not create empty plots
            if plot_type == 'None':
                continue
            
            # get options defined to configure plot (e.g. bias, individual, annotate, etc.)
            plot_options = plot_type.split('_')[1:]
            
            # get zstat information from plot_type
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(plot_type)
            
            # check if plot type is correct for report and library modes
            if self.read_instance.mode in ['report', 'library']:

                # remove plots where setting 'obs' and 'bias' options together
                if ('obs' in plot_options) & ('bias' in plot_options): 
                    msg = f"{plot_type} cannot not be created as 'obs' and 'bias' options set together."
                    show_message(self.read_instance, msg)
                    valid_plot_type = False

                # remove plots when calculating bias stat but temporal_colocation is not active
                elif (z_statistic_type == 'modbias') & (not self.read_instance.temporal_colocation):
                        msg = f'To calculate the model bias stat {zstat}, temporal_colocation must be set to True, so {plot_type} plot cannot be created.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                # if no models are defined, remove all bias plots, or plots with bias statistics 
                elif ('bias' in plot_options) or (z_statistic_sign == 'bias'):
                    if len(data_labels) == 1:
                        msg = f'No models defined, so {plot_type} plot cannot be created.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                # break loop if the plot type is not valid and remove plot type from lists
                if not valid_plot_type:
                    if self.read_instance.mode == 'report':
                        if plot_type in self.read_instance.summary_plots_to_make:
                            self.read_instance.summary_plots_to_make.remove(plot_type)
                        if plot_type in self.read_instance.station_plots_to_make:
                            self.read_instance.station_plots_to_make.remove(plot_type)
                    elif self.read_instance.mode == 'library':
                        return valid_plot_type
                    continue

            # add new keys to make plots with stats (map, periodic, heatmap, table)
            if zstat:

                # get base plot type (without stat and options)
                base_plot_type = plot_type.split('-')[0] 
                # combine basic and modbias stats dicts together
                stats_dict = {**self.read_instance.basic_stats, **self.read_instance.modbias_stats}
                
                # check if plot type is correct for report and library modes
                if self.read_instance.mode in ['report', 'library']:

                    # check all defined plot options are allowed for current plot type
                    invalid_plot_options = [plot_option for plot_option in plot_options if plot_option not in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options']]
                    if len(invalid_plot_options) > 0:
                        msg = f'{plot_type} cannot be created as {invalid_plot_options} plot options are not valid.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False
                    
                    # check desired statistic is defined in stats dict
                    elif base_zstat not in stats_dict:
                        msg = f"{plot_type} cannot be created as {base_zstat} not defined in Providentia's statistical library."
                        show_message(self.read_instance, msg)
                        valid_plot_type = False
                    
                    # remove plots where setting 'obs', but z_statistic_sign is 'bias'
                    elif ('obs' in plot_options) & (z_statistic_sign == 'bias'):
                        msg = f"{plot_type} cannot be created as are plotting a bias statistic but 'obs' option is set."
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for taylor plot if have no models
                    elif (base_plot_type in ['taylor']) & (len(data_labels) == 1):
                        msg = f'No models defined, so {plot_type} cannot be created.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for timeseries bias plot if the temporal colocation is not active
                    elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (not self.read_instance.temporal_colocation):
                        msg = f'{plot_type} cannot be created as temporal colocation is not active.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for periodic plot when active resolution is annual
                    elif ('periodic' == base_plot_type) & (self.read_instance.active_resolution == 'annual'): 
                        msg = f'{plot_type} cannot be created when active temporal resolution is annual'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # break loop if the plot type is not valid and remove plot type from lists
                    if not valid_plot_type:
                        if self.read_instance.mode == 'report':
                            if plot_type in self.read_instance.summary_plots_to_make:
                                self.read_instance.summary_plots_to_make.remove(plot_type)
                            if plot_type in self.read_instance.station_plots_to_make:
                                self.read_instance.station_plots_to_make.remove(plot_type)
                        elif self.read_instance.mode == 'library':
                            return valid_plot_type
                        continue

                # add information for plot type 
                # first try get it from custom plot charcteristics, and then from base plot type template 
                if plot_type in self.canvas_instance.plot_characteristics_templates:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[plot_type])
                else:
                    try:
                        self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])
                    except KeyError:
                        error = f'Error: Plot type {plot_type} is not available. Remove from settings/report_plots.yaml'
                        self.read_instance.logger.error(error)
                        sys.exit(1)

                # overwrite default plot characteristics with custom formatting
                for format_var in format:
                    self.canvas_instance.plot_characteristics[plot_type][format_var] = format[format_var]

                # set page title 
                if 'page_title' in self.canvas_instance.plot_characteristics[plot_type]:
                    if 't' in self.canvas_instance.plot_characteristics[plot_type]['page_title'].keys():
                        if 'bias' in plot_options:
                            self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'],
                                                                                                                          stats_dict[base_zstat]['label'])
                        else:
                            self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} {}'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'], 
                                                                                                                     stats_dict[base_zstat]['label'])
            # add new keys for plots without stats
            else:

                # get base plot type (without options)
                base_plot_type = plot_type.split('_')[0] 

                # check if plot type is correct for report and library modes
                if self.read_instance.mode in ['report', 'library']:

                    # check all defined plot options are allowed for current plot type
                    invalid_plot_options = [plot_option for plot_option in plot_options if plot_option not in self.canvas_instance.plot_characteristics_templates[base_plot_type]['plot_options']]
                    if len(invalid_plot_options) > 0:
                        msg = f'{plot_type} cannot be created as {invalid_plot_options} plot options are not valid.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for scatter and fairmode plots if have no models
                    elif (base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) & (len(data_labels) == 1):
                        msg = f'No models defined, so {plot_type} cannot be created.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for scatter and fairmode plots if the temporal colocation is not active 
                    elif (base_plot_type in ['scatter', 'fairmode-target', 'fairmode-statsummary']) & (not self.read_instance.temporal_colocation):
                        msg = f'{plot_type} cannot be created as temporal colocation is not active.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for timeseries bias plot if the temporal colocation is not active
                    elif ('timeseries' == base_plot_type) & ('bias' in plot_options) & (not self.read_instance.temporal_colocation):
                        msg = f'{plot_type} cannot be created as temporal colocation is not active.'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # warning for periodic-violin plot when active resolution is annual
                    elif ('periodic-violin' == base_plot_type) & (self.read_instance.active_resolution == 'annual'): 
                        msg = f'{plot_type} cannot be created when active temporal resolution is annual'
                        show_message(self.read_instance, msg)
                        valid_plot_type = False

                    # break loop if the plot type is not valid and remove plot type from lists
                    if not valid_plot_type:
                        if self.read_instance.mode == 'report':
                            if plot_type in self.read_instance.summary_plots_to_make:
                                self.read_instance.summary_plots_to_make.remove(plot_type)
                            if plot_type in self.read_instance.station_plots_to_make:
                                self.read_instance.station_plots_to_make.remove(plot_type)
                        elif self.read_instance.mode == 'library':
                            return valid_plot_type
                        continue

                # add information for plot type for base plot type 
                if plot_type in self.canvas_instance.plot_characteristics_templates:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[plot_type])
                else:
                    self.canvas_instance.plot_characteristics[plot_type] = copy.deepcopy(self.canvas_instance.plot_characteristics_templates[base_plot_type])

                # overwrite default plot characteristics with custom formatting
                for format_var in format:
                    self.canvas_instance.plot_characteristics[plot_type][format_var] = format[format_var]

                # change page title if have 'bias' option
                if 'page_title' in self.canvas_instance.plot_characteristics[plot_type]:
                    if 'bias' in plot_options:
                        self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'] = '{} bias'.format(self.canvas_instance.plot_characteristics[plot_type]['page_title']['t'])

            # add figsize for plot orientation
            if 'orientation' in self.canvas_instance.plot_characteristics[plot_type]:
                if self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'landscape':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.landscape_figsize
                elif self.canvas_instance.plot_characteristics[plot_type]['orientation'] == 'portrait':
                    self.canvas_instance.plot_characteristics[plot_type]['figure']['figsize'] = self.canvas_instance.portrait_figsize

        # return valid plot type if library mode
        if self.read_instance.mode == 'library':
            return valid_plot_type

    def make_legend_handles(self, plot_characteristics_legend, data_labels=None, set_obs=True):
        """
        Constructs Matplotlib legend handle objects for observational and model data labels.

        Parameters
        ----------
        plot_characteristics_legend : dict
            Plot characteristics for relevant legend.
        data_labels : list, optional
            Data arrays to plot.
        set_obs : bool, optional
            Indicates if to set observations in legend or not.

        Returns
        -------
        plot_characteristics_legend : dict
            Updated plot characteristics dictionary containing the generated handle objects.
        """

        # if data_labels are not defined, take all in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)

        # create legend elements
        legend_elements = []

        # add observations element, if available, and set_obs == True
        if (self.read_instance.observations_data_label in data_labels) and (set_obs):
            legend_elements.append(Line2D([0], [0], 
                                marker=plot_characteristics_legend['handles']['marker'], 
                                color=plot_characteristics_legend['handles']['color'],
                                markerfacecolor=self.read_instance.plotting_params[self.read_instance.observations_data_label]['colour'],
                                markersize=plot_characteristics_legend['handles']['markersize'], 
                                label=self.read_instance.observations_data_label))
                                  
        # add element for each model
        for model in data_labels:
            if model != self.read_instance.observations_data_label:
                # add model element
                legend_elements.append(Line2D([0], [0], 
                                              marker=plot_characteristics_legend['handles']['marker'],  
                                              color=plot_characteristics_legend['handles']['color'],
                                              markerfacecolor=self.read_instance.plotting_params[model]['colour'],
                                              markersize=plot_characteristics_legend['handles']['markersize'],
                                              label=model))
        
        plot_characteristics_legend['plot']['handles'] = legend_elements
        
        return plot_characteristics_legend

    def make_model_domain_polygons(self, data_labels=None):
        """
        Creates polygon objects representing the geographical boundaries of model domains.

        Parameters
        ----------
        data_labels : list, optional
            Data arrays to plot.

        Returns
        -------
        grid_edge_polygons : list
            List of Matplotlib Polygon objects representing the grid boundaries.
        """

        # if data_labels are not defined, or are just observations, then take all labels in memory
        if data_labels is None:
            data_labels = copy.deepcopy(self.read_instance.data_labels)
        elif len(data_labels) == 1:
            if data_labels[0] == self.read_instance.observations_data_label:
                data_labels = copy.deepcopy(self.read_instance.data_labels)

        grid_edge_polygons = []

        # iterate through read models and plot grid domain edges on map
        for model in data_labels:
            if model != self.read_instance.observations_data_label:
                # create matplotlib polygon object from model grid edge map projection coordinates
                grid_edge_outline_poly = \
                    Polygon(np.vstack((self.read_instance.plotting_params[model]['grid_edge_longitude'],
                                       self.read_instance.plotting_params[model]['grid_edge_latitude'])).T,
                                       edgecolor=self.read_instance.plotting_params[model]['colour'],
                                       transform=self.canvas_instance.datacrs,
                                       **self.canvas_instance.model_domain_polygon)
                # append polygon
                grid_edge_polygons.append(grid_edge_outline_poly)
            
        return grid_edge_polygons

    def draw_line(self, fig, key, value, plot_characteristics, y, spacing, max_len):
        """
        Renders a key-value pair as text on a PDF cover page with automatic line wrapping for long values.

        Parameters
        ----------
        fig : plt.figure
            Cover page.
        key : str
            Key of line.
        value : str
            Value of line.
        plot_characteristics : dict
            Formatting specifications for keys and values.
        y : float
            Y position of last line.
        spacing : float
            Spacing between lines.
        max_len : int
            Maximum characters in value.

        Returns
        -------
        n_lines : int
            The number of lines used to render the value.
        """

        # calculate number of lines for value if it is more than max length characters
        value = str(value)
        lines = [value[i:i + max_len] for i in range(0, len(value), max_len)]
        n_lines = len(lines)

        # add key
        fig.text(y=y + (n_lines - 1) * spacing, s=f"{key}:", **plot_characteristics['keys'])

        # add value, wrapped when too long        
        for j, line in enumerate(reversed(lines)):
            fig.text(y=y + j * spacing, s=line, **plot_characteristics['values'])

        return n_lines
    
    def format_value(self, value):
        """
        Serialises lists or NumPy arrays into comma-separated strings for textual display.

        Parameters
        ----------
        value : list, np.ndarray, str
            The value to be formatted.

        Returns
        -------
        formatted_value : str
            The input value converted to a string.
        """

        if isinstance(value, (list, np.ndarray)):
            return ', '.join(map(str, value))
        return value

    def get_key_by_bsc_name(self, standard_parameters, bsc_name):
        """
        Retrieves the GHOST standard parameter key that corresponds to a specific BSC species name.

        Parameters
        ----------
        standard_parameters : dict
            GHOST standard parameters dictionary.
        bsc_name : str
            BSC species name.

        Returns
        -------
        key : str
            The matching GHOST parameter key, or None if no match is found.
        """

        for key, value in standard_parameters.items():
            if value.get('bsc_parameter_name') == bsc_name:
                return key
        return None

    def make_header(self, pdf, plot_characteristics):
        """
        Creates and renders a formatted cover page for the PDF report, including background, logos, and metadata.

        Parameters
        ----------
        pdf : matplotlib.backends.backend_pdf.PdfPages
            The PDF object used to save the generated page.
        plot_characteristics : dict
            Configuration for visual styling, including dark mode settings and logo placement.
        """

        # load external image as background
        if plot_characteristics['dark_mode']:
            bg_image_path = 'assets/report_dark.png'
            for key in ['keys', 'values']:
                plot_characteristics[key]['color'] = plot_characteristics[key]['color']['dark_mode']
            logo_image_path = plot_characteristics['logo']['path']['dark_mode']
        else:
            bg_image_path = 'assets/report.png'
            for key in ['keys', 'values']:
                plot_characteristics[key]['color'] = plot_characteristics[key]['color']['light_mode']
            logo_image_path = plot_characteristics['logo']['path']['light_mode']
        
        # add background image
        bg_img = Image.open(join(PROVIDENTIA_ROOT, bg_image_path))
        dpi = self.canvas_instance.dpi
        fig_width, fig_height = bg_img.size[0] / dpi, bg_img.size[1] / dpi
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')
        ax.imshow(bg_img)

        # add logo image if path is defined
        if logo_image_path:
            logo_img = Image.open(join(PROVIDENTIA_ROOT, logo_image_path))

            # resize logo to scale of the background width (1/4 by default)
            bg_width, bg_height = bg_img.size
            scale = plot_characteristics['logo']['scale']
            logo_new_width = int(bg_width * scale)
            aspect_ratio = logo_img.height / logo_img.width
            logo_new_height = int(logo_new_width * aspect_ratio)
            logo_img_resized = logo_img.resize((logo_new_width, logo_new_height))
            logo_array = np.asarray(logo_img_resized)
            logo_ax = fig.add_axes([
                plot_characteristics['logo']['x'],
                plot_characteristics['logo']['y'],
                scale,
                logo_new_height / bg_height
            ])
            logo_ax.imshow(logo_array)
            logo_ax.axis('off')

        # if len of network or species uniques is 1, set that instead of long list of duplicates
        _, idx = np.unique(self.read_instance.network, return_index=True)
        network_to_write = np.array(self.read_instance.network)[np.sort(idx)]
        _, idx = np.unique(self.read_instance.species, return_index=True)
        species_to_write = np.array(self.read_instance.species)[np.sort(idx)]

        # get key from standard parameters as name for species
        from GHOST_standards import standard_parameters
        species = []
        for speci in species_to_write:
            species.append(self.get_key_by_bsc_name(standard_parameters, speci))
        
        # format date
        start_date = datetime.strptime(self.read_instance.start_date, "%Y%m%d").strftime("%d/%m/%Y")
        end_date = datetime.strptime(self.read_instance.end_date, "%Y%m%d").strftime("%d/%m/%Y")
        dates = f"{start_date} - {end_date}"

        variable_values = {
            'network': self.format_value(network_to_write),
            'species': self.format_value(species),
            'resolution': self.read_instance.resolution.capitalize().replace('_', ' '),
            'dates': dates,
            'models': self.format_value(list(self.read_instance.experiments.values())),
            'temporal_colocation': self.read_instance.temporal_colocation,
            'spatial_colocation': self.read_instance.spatial_colocation,
            'filter_species': self.format_value(list(self.read_instance.filter_species.keys())),
            'calibration': getattr(self.read_instance, 'calibration_factor', ''),
            'subsections': self.format_value(self.read_instance.subsections)
        }

        visible_items = []
        for key, value in variable_values.items():
            if not value:
                continue
            if plot_characteristics['variables'][key]['active']:
                if plot_characteristics['variables'][key]['value']:
                    value = plot_characteristics['variables'][key]['value']
                visible_items.append((key, value))

        # add key: value as line
        spacing = plot_characteristics['spacing']
        max_len = plot_characteristics['max_len']
        y = plot_characteristics['y_end']
        for key, value in reversed(visible_items):
            n_lines = self.draw_line(
                fig,
                plot_characteristics['variables'][key]['key'],
                value,
                plot_characteristics,
                y,
                spacing,
                max_len
            )
            y += n_lines * spacing 

        # save and close
        pdf.savefig(fig)
        plt.close(fig)

    def make_doi_pdf(self, plot_characteristics, reports_doi_path_temp):
        """
        Create temporary PDF with DOIs per station reference to be added after reordering pages.

        Parameters
        ----------
        plot_characteristics : dict
            Plot characteristics
        reports_doi_path_temp : str
            PDF temporary path

        Returns
        -------
        PdfPages
            PDF
        """

        with PdfPages(reports_doi_path_temp) as pdf:

            for key in ['keys', 'values']:
                plot_characteristics[key]['color'] = plot_characteristics[key]['color']['light_mode']
                
            # Get plot characteristics
            title_spacing = plot_characteristics['title_spacing']
            spacing = plot_characteristics['spacing']
            chunk_size = plot_characteristics['chunk_size']
                
            # Set figure dimensions
            dpi = self.canvas_instance.dpi
            bg_img = [1414, 2000]
            fig_width, fig_height = bg_img[0] / dpi, bg_img[1] / dpi
        
            # Get DOI per station
            for subsection in self.canvas_instance.station_reference_data:
                for networkspeci in self.read_instance.networkspecies:
                    # DOI pages can only be created for ACTRIS
                    if networkspeci.split('|')[0] != 'actris/actris':
                        continue

                    reference_data = self.canvas_instance.station_reference_data[subsection][networkspeci]
                    doi_data = self.canvas_instance.station_doi_data[subsection][networkspeci]

                    visible_items = []
                    extra_rows = 0
                    
                    for station_reference, station_doi in zip(reference_data, doi_data):
                        # Show DOIs only if station reference is consistent across time
                        unique_station_references = list(
                            np.unique(station_reference)
                        )
                        if len(unique_station_references) > 1:
                            error = (
                                f'Error: Station references change with time '
                                f'({unique_station_references}), cannot create DOI page.'
                            )
                            self.read_instance.logger.error(error)
                            return None

                        # Get DOIs per station reference
                        unique_station_dois = list(
                            np.unique(station_doi)
                        )
                        if len(unique_station_dois) > 1:
                            extra_rows += len(unique_station_dois) - 1

                        visible_items.append(
                            (unique_station_references, unique_station_dois)
                        )
                    
                    # If there are extra rows, chunk size needs to be lower to see all DOIs in a page
                    if extra_rows > 2:
                        msg = (
                            f'Warning: Some stations have multiple DOIs, '
                            f'lower chunk size in plot_characteristics.yaml '
                            f'(current: {chunk_size}).'
                        )
                        self.read_instance.logger.info(msg)

                    # Create a page for every 30 stations 
                    num_pages = math.ceil(len(visible_items) / chunk_size)

                    for page_idx in range(num_pages):
                        y = plot_characteristics['y_end']
                        fig = plt.figure(figsize=(fig_width, fig_height), dpi=dpi)
                        ax = fig.add_axes([0, 0, 1, 1])
                        ax.axis('off')

                        # Page title
                        fig.text(
                            y=y,
                            s=(
                                f"List of DOIs - {subsection} - {networkspeci} "
                                f"(Page {page_idx+1}/{num_pages})"
                            ),
                            fontweight='bold',
                            **plot_characteristics['keys']
                        )
                        y -= title_spacing

                        # Slice visible_items for this page
                        start = page_idx * chunk_size
                        end = start + chunk_size
                        page_items = visible_items[start:end]

                        # Show station (key) and DOI data (value)
                        for stations, dois in page_items:
                            
                            # Save starting height of this block
                            block_y = y

                            # Write stations
                            y = block_y
                            for i, station in enumerate(stations):
                                fig.text(y=y, s=station, **plot_characteristics['keys'])
                                if i < len(stations) - 1:
                                    y -= spacing

                            # Write DOIs
                            y = block_y
                            for j, doi in enumerate(dois):
                                fig.text(y=y, s=doi, **plot_characteristics['values'])
                                if j < len(dois) - 1:
                                    y -= spacing

                            # After finishing block, lower height once
                            y = min(y, block_y) - spacing

                        # Save page to PDF
                        pdf.savefig(fig)
                        plt.close(fig)

        return pdf

    def make_metadata(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """
        Generates and renders a text-based summary of station metadata, including statistical aggregations of attributes.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : str
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """

        # initialise string to plot on axis
        str_to_plot = ''
        
        # get number of selected stations
        n_stations = self.canvas_instance.selected_station_metadata[networkspeci].shape[0]

        # set first line of metadata print to be either selected stations or number of stations 
        if n_stations == 1:
            var_str_name = 'name_one' 
            station_references = self.canvas_instance.selected_station_metadata[networkspeci][
                                 'station_reference'].flatten()
            station_reference = station_references[~pd.isnull(station_references)].astype(str)[0]
            str_to_plot += 'Station: {}'.format(station_reference)
        else:
            var_str_name = 'name_multiple' 
            str_to_plot += '{} Stations'.format(n_stations)
        
        # iterate n vars per line and add spacing
        current_n_vars_per_line = 1

        # set up read for GHOST and non-GHOST data types 
        if self.read_instance.reading_ghost:
            char_var = 'ghost_vars'
        else:
            char_var = 'non-ghost_vars'
    
        # iterate through defined variables and add them
        for ghost_var, ghost_var_dict in plot_characteristics[char_var].items():

            # check var str name exists in ghost var dict, if not move to next var
            if var_str_name not in ghost_var_dict:
                continue
            
            # check if var is in metadata
            if ghost_var not in self.read_instance.metadata_vars_to_read:
                continue

            # if are on limit of vars allowed per line then break to new line
            if current_n_vars_per_line == plot_characteristics['max_vars_per_row']:
                str_to_plot += '\n'
                current_n_vars_per_line = 0
            # otherwise, add spacing between variables on line
            else:
                str_to_plot += (' '*plot_characteristics['var_spacing'])

            # round decimal places if float
            if 'dp' in ghost_var_dict:
                str_to_plot += '{}: {:.{}f}'.format(ghost_var_dict[var_str_name],
                                                    np.nanmedian(self.canvas_instance.selected_station_metadata[
                                                                 networkspeci][ghost_var].astype(np.float32)),
                                                    ghost_var_dict['dp'])

            # if str then get unique elements or percentage dependent on n uniques
            else:
                # gather all selected station metadata for current meta variable
                all_current_meta = self.canvas_instance.selected_station_metadata[networkspeci][
                                   ghost_var].flatten()
                # remove NaNs
                all_current_meta = all_current_meta[~pd.isnull(all_current_meta)].astype(str)

                # get counts of all unique metadata elements across selected stations
                unique_meta, meta_counts = np.unique(all_current_meta, return_counts=True)
                # get number of unique metadata elements across selected stations
                n_unique_meta = len(unique_meta)

                # 0 unique metadata elements, then must be all NaN, so set as NaN
                if n_unique_meta == 0:
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], 'nan')
                # 1 unique metadata element? then return it
                elif n_unique_meta == 1:
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], unique_meta[0])
                # if have > N unique metadata elements, just return count of the elements across the selected stations
                elif n_unique_meta > plot_characteristics['max_uniques']:
                    str_to_plot += '{}: {} uniques'.format(ghost_var_dict[var_str_name], n_unique_meta)
                # otherwise, get percentage of unique metadata elements across selected stations
                else:
                    meta_pc = (100. / len(all_current_meta)) * meta_counts
                    meta_pc = ['{:.1f}%'.format(meta) for meta in meta_pc]
                    # create string for variable to plot
                    str_to_plot += '{}: {}'.format(ghost_var_dict[var_str_name], ', '.join(
                        [':'.join([str(var), pc]) for var, pc in zip(unique_meta, meta_pc)]))

            # add units
            if 'units' in ghost_var_dict:
                if ghost_var_dict['units'] != '':
                    str_to_plot += ' {}'.format(ghost_var_dict['units'])

            # iterate current_n_vars_per_line
            current_n_vars_per_line += 1

        # plot string to axis
        plot_txt = relevant_axis.text(0.0, 1.0, str_to_plot, transform=relevant_axis.transAxes, **plot_characteristics['plot'])

        # modify limit to wrap text as axis width in pixels
        if self.read_instance.mode == 'report': 
            
            # get axis dimensions in pixels
            ax_width_px = relevant_axis.bbox.width * plot_characteristics['figure']['nrows']
            
        elif self.read_instance.mode == 'library':

            # get axis dimensions in pixels
            ax_width_px = relevant_axis.bbox.width

        else:
            # get axis bounding box
            ax_bbox = relevant_axis.get_window_extent().transformed(self.canvas_instance.figure.dpi_scale_trans.inverted())
            
            # get axis dimensions in inches
            ax_width_inches = ax_bbox.width

            # get axis dimensions in pixels
            ax_width_px = ax_width_inches * self.canvas_instance.figure.dpi

        # automatically sets limit as figure width
        plot_txt._get_wrap_line_width = lambda: ax_width_px

        # track plot elements if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            self.track_plot_elements(self.read_instance.observations_data_label, 'metadata', 'plot', [plot_txt], bias=False)

    def make_map(self, relevant_axis, networkspeci, plot_characteristics, plot_options, zstat=None, labela='', 
                 labelb=''):
        """
        Renders a geospatial scatter plot of stations onto a map axis, coloured by a calculated statistical metric.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstat : str, optional
            Statistic to plot.
        labela : str, optional
            Label of first dataset.
        labelb : str, optional
            Label of second dataset (if defined then a bias plot is made).
        """

        # calculate statistic
        z_statistic, active_map_valid_station_inds = calculate_statistic(self.read_instance, self.canvas_instance, 
                                                                         networkspeci, zstat, [labela], [labelb], 
                                                                         map=True)

        # get marker size (for report and library)
        if self.read_instance.mode in ['report', 'library']:
            self.get_markersize(relevant_axis, 'map', networkspeci, plot_characteristics, 
                                active_map_valid_station_inds=active_map_valid_station_inds)
        # if using dashboard make z_statistic and active_map_valid_station_inds class variables
        else:
            self.canvas_instance.z_statistic = z_statistic
            self.canvas_instance.active_map_valid_station_inds = active_map_valid_station_inds

        # plot new station points on map - coloured by currently active z statisitic
        self.stations_scatter = relevant_axis.scatter(self.read_instance.station_longitudes[networkspeci][active_map_valid_station_inds], 
                                                      self.read_instance.station_latitudes[networkspeci][active_map_valid_station_inds], 
                                                      c=z_statistic, transform=self.canvas_instance.datacrs,
                                                      **plot_characteristics['plot'])
        
        # track plot elements if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            self.track_plot_elements(self.read_instance.observations_data_label, 'map', 'plot', [self.stations_scatter], bias=False)

    def make_timeseries(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, 
                        chunk_stat=None, chunk_resolution=None):
        """
        Renders temporal data as a line plot, optionally applying statistical chunking or bias calculations.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        chunk_stat : str, optional
            Chunk statistic.
        chunk_resolution : str, optional
            Chunk resolution.
        """

        # create list for timeseries plot
        self.timeseries_plot = []

        # skip making timeseries (points) for report and library mode
        # we do not apply this in the dashboard to avoid being unable to see the points on certain changes
        if (self.read_instance.mode in ['report', 'library']) and ('hidedata' in plot_options):
            return

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # get bias and set if bias line will be added
        if 'bias' in plot_options:
            bias =  True
            add_bias_line = True
        else:
            bias = False
            add_bias_line = False

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get chunking stat and resolution in dashboard
        if self.read_instance.mode not in ['report', 'library']:
            chunk_stat = self.canvas_instance.timeseries_chunk_stat.currentText()
            chunk_resolution = self.canvas_instance.timeseries_chunk_resolution.currentText()
            chunk_stat = None if chunk_stat == 'None' else chunk_stat
            chunk_resolution = None if chunk_resolution == 'None' else chunk_resolution
        
        # chunk timeseries
        if (chunk_stat is not None) and (chunk_resolution is not None):
            timeseries_data = create_statistical_timeseries(self.read_instance, self.canvas_instance, chunk_stat, 
                                                            chunk_resolution, networkspeci, cut_data_labels, bias)

            # if it is a bias chunk statistic, add bias line
            z_statistic_type = get_z_statistic_type(chunk_stat)
            if z_statistic_type == 'modbias':
                add_bias_line = True
        # normal timeseries
        else:
            timeseries_data = self.canvas_instance.selected_station_data[networkspeci]["timeseries"]
        
        # plot horizontal line across x axis at 0 if bias plot
        if add_bias_line:
            if 'bias' not in plot_options:
               plot_characteristics['bias_line']['y'] = self.read_instance.modbias_stats[chunk_stat]['minimum_bias']
            bias_line = relevant_axis.axhline(**plot_characteristics['bias_line'])
            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements('ALL', 'timeseries', 'bias_line', [bias_line], bias=bias)

        # iterate through data labels
        for data_label in cut_data_labels:

            # bias plot?
            if bias:

                # skip if data label is for observations
                if data_label == self.read_instance.observations_data_label:
                    continue
                
                # chunk bias timeseries
                if (chunk_stat is not None) and (chunk_resolution is not None):
                    ts = timeseries_data[data_label] 
                # normal bias timeseries
                else:
                    ts_obs = timeseries_data[self.read_instance.observations_data_label]
                    ts_model = timeseries_data[data_label] 
                    ts = ts_model - ts_obs

            else:
                ts = timeseries_data[data_label]
            
            # get marker size (for report and library)
            if self.read_instance.mode in ['report', 'library']:
                self.get_markersize(relevant_axis, 'timeseries', networkspeci, plot_characteristics, data=ts)

            # make timeseries plot
            self.timeseries_plot.append(relevant_axis.plot(ts, 
                                                      color=self.read_instance.plotting_params[data_label]['colour'], 
                                                      **plot_characteristics['plot']))

            # update maximum smooth value
            if self.read_instance.mode not in ['report', 'library']:
                self.canvas_instance.timeseries_smooth_window_sl.setMaximum(len(ts))
                # To get straight line
                # if self.canvas_instance.timeseries_smooth_window_sl.value() != (len(ts)*2 - 1):
                #     self.canvas_instance.timeseries_smooth_window_sl.setMaximum(int(len(ts)*2 - 1))

            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements(data_label, 'timeseries', 'plot', self.timeseries_plot[-1], bias=bias)

    def make_periodic(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, zstat=None):
        """
        Renders periodic cycles as either line plots or split-violin distributions across temporal resolutions.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstat : str, optional
            Statistic to plot.
        """

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # get zstat information
        if zstat is not None:
            zstat, base_zstat, z_statistic_type, z_statistic_sign, z_statistic_period = get_z_statistic_info(zstat=zstat)

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]
        
        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get number of models in data  labels
        if self.read_instance.observations_data_label in cut_data_labels:
            n_mods = len(cut_data_labels) - 1
        else:
            n_mods = len(cut_data_labels) 

        # hide non-relevant resolution axes
        for nonrelevant_temporal_resolution in self.read_instance.periodic_nonrelevant_temporal_resolutions:
            # get subplot axis
            relevant_sub_ax = relevant_axis[nonrelevant_temporal_resolution]
            relevant_sub_ax.axis('off')
            relevant_sub_ax.set_visible(False)

        # iterate through all relevant temporal aggregation 
        for relevant_temporal_resolution in self.read_instance.periodic_relevant_temporal_resolutions:

            # get subplot axis
            relevant_sub_ax = relevant_axis[relevant_temporal_resolution]

            # un-hide relevant resolution axis
            relevant_sub_ax.axis('on')
            relevant_sub_ax.set_visible(True)

            # violin plot type?
            if zstat is None:

                # calculate medians
                medians = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'Median', 
                                              cut_data_labels, [], period=relevant_temporal_resolution)

                # calculate PDF for data label
                period_x_grid, PDFs_sampled = self.make_distribution(relevant_axis, networkspeci, data_labels, 
                                                                     plot_characteristics, plot_options,
                                                                     violin_resolution=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # list to save all violins per data label
                    violins = []

                    # get median zorder
                    median_zorder = self.read_instance.plotting_params[data_label]['zorder']+len(cut_data_labels)
                    
                    # get alpha and violin fill information 
                    if data_label == self.read_instance.observations_data_label:
                        alpha = plot_characteristics['violin_alphas']['alpha_obs']
                        violin_fill = plot_characteristics['violin_fill_obs']
                    else:
                        alpha = plot_characteristics['violin_alphas']['alpha_mod']
                        if (n_mods == 1) or ('individual' in plot_options):
                            violin_fill = plot_characteristics['violin_fill_1model']
                        else:
                            violin_fill = plot_characteristics['violin_fill_2+models']

                    # make plot of median
                    median_plots = relevant_sub_ax.plot(self.canvas_instance.unique_periods, medians[:, data_label_ii], 
                                                        color=self.read_instance.plotting_params[data_label]['colour'], 
                                                        zorder=median_zorder, 
                                                        **plot_characteristics['plot']['median'])

                    # make violin plot
                    for period_ii in range(len(self.canvas_instance.unique_periods)):

                        # get x_grid for period
                        x_grid = period_x_grid[period_ii]

                        # set xaxis position of each violin
                        if relevant_temporal_resolution == 'hour':
                            if self.read_instance.active_resolution == '3hourly':
                                period_pos = period_ii * 3
                            elif self.read_instance.active_resolution  == '6hourly':
                                period_pos = period_ii * 6
                            else:
                                period_pos = period_ii
                        elif relevant_temporal_resolution == 'month':
                            period_pos = period_ii + 1
                        else:
                            period_pos = period_ii

                        PDF_sampled = PDFs_sampled[data_label_ii, period_ii, :]
                        if not np.all(np.isnan(PDF_sampled)):
                            
                            PDF_sampled = 0.5 * plot_characteristics['violin_widths'] * PDF_sampled / PDF_sampled.max()

                            # make violin plot (filled or unfilled)
                            if violin_fill:
                                self.violin_plot = relevant_sub_ax.fill_betweenx(x_grid, -PDF_sampled + period_pos, PDF_sampled + period_pos,
                                                                                 facecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                                                 alpha=alpha,
                                                                                 **plot_characteristics['plot']['violin'])
                            else:
                                self.violin_plot = relevant_sub_ax.fill_betweenx(x_grid, -PDF_sampled + period_pos, PDF_sampled + period_pos,
                                                                                 facecolor='None', edgecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                                                                 **plot_characteristics['plot']['violin'])

                            # if have more than 1 valid data array (both obs and model), 
                            # split the violin plot across the horizontal
                            # (observations on left, model violin_plots on right)
                            if ((n_mods > 0) and (self.read_instance.observations_data_label in cut_data_labels) and
                               ('individual' not in plot_options)):
                                m = np.mean(self.violin_plot.get_paths()[0].vertices[:, 0])
                                # observations on left
                                if data_label == self.read_instance.observations_data_label:
                                    self.violin_plot.get_paths()[0].vertices[:, 0] = np.clip(self.violin_plot.get_paths()[0].vertices[:, 0], -np.inf, m)
                                # models on right
                                else:
                                    self.violin_plot.get_paths()[0].vertices[:, 0] = np.clip(self.violin_plot.get_paths()[0].vertices[:, 0], m, np.inf)

                            # save violin
                            violins.append(self.violin_plot)
                    
                    # add hidden poins for data limits per data label to allow for limit harmonisation
                    if len(violins) > 0:
                        limit_plot = relevant_sub_ax.plot([period_pos,period_pos], 
                                                          [np.min(period_x_grid),np.max(period_x_grid)],
                                                           alpha=0.0)
                        violins.extend(limit_plot)

                    # track plot elements if using dashboard 
                    if self.read_instance.mode not in ['report', 'library']:
                        self.track_plot_elements(data_label, 'periodic-violin', 'violin_plot_{}'.format(relevant_temporal_resolution), violins, bias=False)
                        self.track_plot_elements(data_label, 'periodic-violin', 'Median_plot_{}'.format(relevant_temporal_resolution), median_plots, bias=False)

            # periodic plot type
            else:

                # plot horizontal line/s across x axis at value/s of minimum model bias (if bias statistic)
                if z_statistic_sign == 'bias':
                    # get value/s of minimum bias for statistic
                    if z_statistic_type == 'basic':
                        minimum_bias = self.read_instance.basic_stats[base_zstat]['minimum_bias']
                    else:
                        minimum_bias = self.read_instance.modbias_stats[base_zstat]['minimum_bias']
                    bias_lines = []
                    for mb in minimum_bias:
                        bias_lines += [relevant_sub_ax.axhline(y=mb, **plot_characteristics['bias_line'])]
                    # track plot elements if using dashboard 
                    if self.read_instance.mode not in ['report', 'library']:
                        self.track_plot_elements('ALL', 'periodic', 'bias_line_{}'.format(relevant_temporal_resolution), 
                                                 bias_lines, bias=bias)

                # calculate statistic in each of periodic groups per data label
                if z_statistic_sign == 'bias':
                    if self.read_instance.observations_data_label in cut_data_labels:
                        cut_data_labels.remove(self.read_instance.observations_data_label)
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                                    [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                    cut_data_labels, period=relevant_temporal_resolution)
                else:
                    stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 
                                                     zstat, cut_data_labels, [], period=relevant_temporal_resolution)

                # iterate through data labels and plot violins
                for data_label_ii, data_label in enumerate(cut_data_labels): 

                    # skip observational array if bias stat
                    if (z_statistic_sign == 'bias') & (data_label == self.read_instance.observations_data_label):
                        continue

                    # make plot
                    self.periodic_plots = relevant_sub_ax.plot(self.canvas_instance.unique_periods, stats_calc[:, data_label_ii], 
                                                               color=self.read_instance.plotting_params[data_label]['colour'], 
                                                               zorder=self.read_instance.plotting_params[data_label]['zorder'], 
                                                               **plot_characteristics['plot'])

                    # track plot elements if using dashboard 
                    if self.read_instance.mode not in ['report', 'library']:
                        self.track_plot_elements(data_label, 'periodic', 'plot_{}'.format(relevant_temporal_resolution), self.periodic_plots, bias=bias)

    def make_distribution(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                          data_range_min=None, data_range_max=None, violin_resolution=None):
        """
        Computes and renders probability density functions using Fast Fourier Transform-based Kernel Density Estimation.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        data_range_min : float, optional
            Minimum data range of distribution plot grid.
        data_range_max : float, optional
            Maximum data range of distribution plot grid.
        violin_resolution : int, optional
            If are calculating distribution for violin plot, this is set to temporal resolution of groupings.

        Returns
        -------
        period_x_grid : list
            The grid of x-values used for periodic distribution calculations (only if violin_resolution is set).
        PDFs_sampled : np.ndarray
            The sampled PDF values for each data label and period (only if violin_resolution is set).
        """

        # determine if 'bias' in plot_options
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if 'obs' in plot_options, set data labels to just observations label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # set data ranges for distribution plot grid if not set explicitly
        if data_range_min is None:
            data_range_min = self.canvas_instance.selected_station_data_min[networkspeci]
        
        if data_range_max is None:
            data_range_max = self.canvas_instance.selected_station_data_max[networkspeci]

        # set xgrid for calculating distribution
        # if calculating for period n_samples is set to pdf_min_samples
        # otherwise it is inferred from data (if above minimum value)
        if violin_resolution is not None:
            n_samples = plot_characteristics['pdf_min_samples']
        else:
            minmax_diff = data_range_max - data_range_min
            if pd.isnull(self.read_instance.parameter_dictionary[networkspeci.split('|')[1]]['minimum_resolution']):
                n_samples = plot_characteristics['pdf_min_samples']
            else:
                n_samples = int(np.around(minmax_diff/(self.read_instance.parameter_dictionary[networkspeci.split('|')[1]]['minimum_resolution']/100.0),0))
                if n_samples < plot_characteristics['pdf_min_samples']:
                    n_samples = plot_characteristics['pdf_min_samples']

        # round n_samples to next next power of 2 (for fft optimisation)
        n_samples = 2 ** np.ceil(np.log2(n_samples)) 
        # set x_grid
        x_grid = np.linspace(data_range_min,data_range_max,int(n_samples))

        # plot horizontal line across x axis at 0 if bias plot
        # also remove observations from cut_data_labels
        if bias:
            bias_line = [relevant_axis.axhline(**plot_characteristics['bias_line'])]
            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements('ALL', 'distribution', 'bias_line', bias_line, bias=bias)
            if self.read_instance.observations_data_label in cut_data_labels:
                cut_data_labels.remove(self.read_instance.observations_data_label)

        # if violin plot setup arrays for saving data
        if violin_resolution is not None:
            PDFs_sampled = np.full((len(cut_data_labels), len(self.canvas_instance.periodic_xticks[violin_resolution]), int(n_samples)), np.nan, dtype=np.float32)

        # iterate through data labels
        for data_label_ii, data_label in enumerate(cut_data_labels):
        
            PDF_sampled_calculated = False

            # setup bias plot
            if bias:

                # calculate obs PDF on first pass
                if data_label_ii == 0:
                    kde_data_obs = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(self.read_instance.observations_data_label),0,:])
                    
                    #filter out data outside data range bounds
                    #kde_data_obs = kde_data_obs[(kde_data_obs > data_range_min) & (kde_data_obs < data_range_max)]
                    
                    # check if all values are equal in the dataframe
                    if kde_data_obs.size == 0:
                        msg = 'The kernel density cannot be calculated because there are no valid observational values.'
                        show_message(self.read_instance, msg)
                        return
                    elif np.all(kde_data_obs == kde_data_obs[0]):
                        msg = 'The kernel density cannot be calculated because all observational values are equal.'
                        show_message(self.read_instance, msg)
                        return
                    else:
                        PDF_obs_sampled = kde_fft(kde_data_obs, xgrid=x_grid)

                        if isinstance(PDF_obs_sampled, str):
                            msg = PDF_obs_sampled
                            msg += f'The distribution plot will be created and not include data for this label ({data_label}). '
                            show_message(self.read_instance, msg)
                            continue

                # calculate model PDF
                kde_data_model = drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])
                
                #filter out data outside data range bounds
                #kde_data_model = kde_data_model[(kde_data_model > data_range_min) & (kde_data_model < data_range_max)]
                
                # check if all values are equal in the dataframe
                if kde_data_model.size == 0:
                    msg = 'The kernel density cannot be calculated because there are no valid values for {} model.'.format(data_label)
                    show_message(self.read_instance, msg)
                    continue
                elif np.all(kde_data_model == kde_data_model[0]):
                    msg = 'The kernel density cannot be calculated because all values for {} model are equal.'.format(data_label)
                    show_message(self.read_instance, msg)
                    continue
                # calculate PDF
                PDF_model_sampled = kde_fft(kde_data_model, xgrid=x_grid)
                if isinstance(PDF_model_sampled, str):
                    msg = PDF_model_sampled
                    msg += f'The distribution plot will be created and not include data for this label ({data_label}). '
                    show_message(self.read_instance, msg)
                    continue

                # calculate difference
                PDF_sampled = PDF_model_sampled - PDF_obs_sampled
                PDF_sampled_calculated = True

            # setup standard plot
            else:
                
                # if first data label and calculating distributions for violin plot,
                # calculate the x_grid / data ranges per period  
                # use min for min data range and upper inner Tukey fence for max data range
                if (violin_resolution is not None) & (data_label_ii == 0):
                    period_data_range_min = []
                    period_data_range_max = []
                    period_x_grid = []

                    # get grouped data per period
                    grouped_data = group_periodic(self.read_instance, self.canvas_instance, networkspeci,
                                                  violin_resolution, False, self.read_instance.statistic_mode, '', 
                                                  self.canvas_instance.selected_station_data[networkspeci]['per_station'])

                    # iterate through periods
                    for group in grouped_data:
                        lower_inner_fence, upper_inner_fence = boxplot_inner_fences(group)
                        min_data = np.nanmin(group)
                        period_data_range_min.append(min_data)
                        period_data_range_max.append(upper_inner_fence)
                        period_x_grid.append(np.linspace(min_data,upper_inner_fence,int(n_samples)))
                
                # get data (flattened and drop NaNs)
                if violin_resolution is not None:
                    kde_data_grouped = [drop_nans(group[valid_data_labels.index(data_label)].flatten()) for group in grouped_data]
                else:    
                    kde_data_grouped = [drop_nans(self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:])]

                # iterate through kde data groups
                for period_ii, kde_data in enumerate(kde_data_grouped):

                    # get relevant data ranges / x_grid, for violin period distribution calculation
                    if violin_resolution is not None:
                        data_range_min = period_data_range_min[period_ii]
                        data_range_max = period_data_range_max[period_ii]
                        x_grid = period_x_grid[period_ii]
                    
                    #filter out data outside data range bounds
                    #kde_data = kde_data[(kde_data > data_range_min) & (kde_data < data_range_max)]

                    # check if all values are equal in the dataframe
                    if kde_data.size == 0:
                        if violin_resolution is None:
                            msg = 'The kernel density cannot be calculated because there are no valid values for {}.'.format(data_label)
                            show_message(self.read_instance, msg)
                        continue
                    elif np.all(kde_data == kde_data[0]):
                        if violin_resolution is None:
                            msg = 'The kernel density cannot be calculated because all {} values are equal.'.format(data_label)
                            show_message(self.read_instance, msg)
                        continue
                    else:
                        PDF_sampled = kde_fft(kde_data, xgrid=x_grid)
                        if isinstance(PDF_sampled, str):
                            msg = PDF_sampled
                            msg += f'The distribution plot will be created and not include data for this label ({data_label}). '
                            show_message(self.read_instance, msg)
                            continue
                        
                        # save PDF for violin plot
                        if violin_resolution is not None:
                            PDFs_sampled[data_label_ii, period_ii, :] = PDF_sampled
                        else:
                            PDF_sampled_calculated = True

            # make plot (if not making PDF for violin plot)
            if PDF_sampled_calculated:
                # make plot
                self.distribution_plot = relevant_axis.plot(x_grid, PDF_sampled, 
                                                            color=self.read_instance.plotting_params[data_label]['colour'], 
                                                            **plot_characteristics['plot'])

                # track plot elements if using dashboard 
                if self.read_instance.mode not in ['report', 'library']:
                    self.track_plot_elements(data_label, 'distribution', 'plot', self.distribution_plot, bias=bias)

        # if have made PDFs for violin plot then return it
        if violin_resolution is not None:
            return period_x_grid, PDFs_sampled

    def make_scatter(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """
        Renders a scatter plot comparing model predictions against observations, including reference ratio lines.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """
        
        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # add 1:1 line (if in plot_characteristics)
        if '1:1_line' in plot_characteristics:
            relevant_axis.plot([0, 1], [0, 1], transform=relevant_axis.transAxes, 
                               **plot_characteristics['1:1_line'])
        # add 1:2 line (if in plot_characteristics)
        if '1:2_line' in plot_characteristics:
            relevant_axis.plot([0, 1], [0, 0.5], transform=relevant_axis.transAxes, 
                                **plot_characteristics['1:2_line'])     
        # add 2:1 line (if in plot_characteristics)
        if '2:1_line' in plot_characteristics:
            relevant_axis.plot([0, 0.5], [0, 1], transform=relevant_axis.transAxes, 
                               **plot_characteristics['2:1_line'])

        # skip making scatter for report and library mode
        # we do not apply this in the dashboard to avoid being unable to see the points on certain changes
        if (self.read_instance.mode in ['report', 'library']) and ('hidedata' in plot_options):
            return
        
        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get observations data (flattened)
        observations_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(self.read_instance.observations_data_label),0,:]

        # determine if number of points per data array exceeds max limit,
        # if so subset arrays
        subset = False
        data_array_size = observations_data.size
        if 'max_points' in plot_characteristics:
            if data_array_size > plot_characteristics['max_points']:
                subset = True
                inds_subset = np.random.choice(data_array_size, size=plot_characteristics['max_points'], replace=False)
                observations_data = observations_data[inds_subset]

        # iterate through data labels
        for data_label in cut_data_labels:

            # continue for observations data label
            if data_label == self.read_instance.observations_data_label:
                continue

            # get model data (flattened)
            model_data = self.canvas_instance.selected_station_data[networkspeci]['flat'][valid_data_labels.index(data_label),0,:]

            # subset data if neccessary
            if subset:
                model_data = model_data[inds_subset]

            # get marker size (for report and library)
            if self.read_instance.mode in ['report', 'library']:
                self.get_markersize(relevant_axis, 'scatter', networkspeci, plot_characteristics, data=observations_data)

            # create scatter plot
            self.scatter_plot = relevant_axis.plot(observations_data, model_data, 
                                                   color=self.read_instance.plotting_params[data_label]['colour'],
                                                   **plot_characteristics['plot'])

            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements(data_label, 'scatter', 'plot', self.scatter_plot, bias=False)

    def make_boxplot(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """
        Renders box-and-whisker plots to visualise data distributions across observations and models.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """

        # if 'obs' in plot_options, set data labels to just observations data label
        if 'obs' in plot_options:
            data_labels = [self.read_instance.observations_data_label]

        # if multispecies in plot options then make plot for all networkspecies
        if 'multispecies' in plot_options:
            networkspecies = self.read_instance.networkspecies
            species = self.read_instance.species
        else:
            networkspecies = [networkspeci]
            species = [networkspeci.split('|')[1]]

        # determine if have all species are size distribution species
        if np.all([True if 'vconcaerobin' in speci else False for speci in species]):
            sizedist = True
            radii_bins = np.array([get_AERONET_sizedist_bin_radius(speci) for speci in species]).astype(np.float32)
            bin_widths = np.diff(np.log(radii_bins))
            bin_widths = np.append(bin_widths, bin_widths[-1])
        else:
            sizedist = False  

        # if normalise in plot options, then get factor for normalisation (per observations and model)
        if ('normalise' in plot_options): 

            # initialise dict to store normalisation factors
            norm_factor = {}

            # iterate through data labels
            for data_label in data_labels:
                
                # initialse list to store bin means
                y_bins = [] 
                
                # iterate through networkspecies
                for ns in networkspecies:
            
                    # get valid data labels for networkspeci
                    valid_data_labels = self.canvas_instance.selected_station_data_labels[ns]
                
                    # if data label not in valid data labels, continue
                    if data_label not in valid_data_labels:
                        continue

                    # get data (flattened and drop NaNs)
                    data_array = drop_nans(self.canvas_instance.selected_station_data[ns]['flat'][valid_data_labels.index(data_label),0,:])

                    # For area calculation, take mean over time for each bin
                    y_bins.append(np.nanmean(data_array))

                # create numpy array of bin means
                y_bins = np.array(y_bins)

                # calculate area under curve for normalisation
                if sizedist:
                    widths_used = bin_widths[:len(y_bins)]
                    area = np.sum(y_bins * widths_used)
                else:               
                    area = np.nansum(y_bins)
                norm_factor[data_label] = area
            
        # iterate through networkspecies
        ns_current = 0
        for ns in networkspecies:

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[ns]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # only proceed if have some data labels to plot
            if len(cut_data_labels) > 0:

                # get data label width and spacing
                if (('individual' in plot_options) or ('obs' in plot_options) 
                    or (len(self.read_instance.networkspecies) == 1) or (len(cut_data_labels) == 1)):
                    widths = plot_characteristics['group_widths']['singlespecies']
                else:
                    available_width = plot_characteristics['group_widths']['multispecies']
                    remainder_width = 1.0 - available_width
                    start_point = -0.5 + (remainder_width / 2.0)
                    widths = available_width / (len(cut_data_labels) + 0.15)
                    spacing = (available_width - (widths * len(cut_data_labels))) / (len(cut_data_labels) - 1)

                # get plot positions
                if ('individual' in plot_options) or ('obs' in plot_options) or (len(cut_data_labels) == 1):
                    positions = [ns_current]
                elif (len(self.read_instance.networkspecies) == 1):
                    positions = np.arange(len(cut_data_labels))
                else:
                    positions = [((start_point + (widths/2.0)) + (spacing * data_label_ii) + (widths * data_label_ii)) + ns_current 
                                for data_label_ii in range(len(cut_data_labels))]  
                
                # iterate ns_current
                ns_current += 1

                # iterate through cut data labels making plot
                for data_label_ii, data_label in enumerate(cut_data_labels):

                    # get data (flattened and drop NaNs)
                    data_array = drop_nans(self.canvas_instance.selected_station_data[ns]['flat'][valid_data_labels.index(data_label),0,:])

                    # normalise data array
                    if 'normalise' in plot_options:
                        data_array = data_array / norm_factor[data_label]

                    # make boxplot
                    boxplot = relevant_axis.boxplot(x=data_array, positions=[positions[data_label_ii]], widths=widths, 
                                                    **plot_characteristics['plot'])

                    # set box colour
                    for element in ['boxes', 'whiskers', 'fliers', 'medians', 'caps']:
                        plt.setp(boxplot[element], color=self.read_instance.plotting_params[data_label]['colour'])
                    # set fill colour to be white
                    for patch in boxplot['boxes']:
                        patch.set(facecolor='white')

                    # track plot elements if using dashboard 
                    if self.read_instance.mode not in ['report', 'library']:
                        self.track_plot_elements(data_label, 'boxplot', 'plot', boxplot, bias=False)

        # set xticklabels 
        # labels for multispecies plot
        xtick_params = copy.deepcopy(plot_characteristics['xtick_params'])
        xticklabel_params = copy.deepcopy(plot_characteristics['xticklabels'])
        if ('multispecies' in plot_options) & (len(self.read_instance.networkspecies) > 1):
            xticks = np.arange(len(self.read_instance.networkspecies))
            #if all networks or species are same, drop them from xtick label
            if len(np.unique(self.read_instance.network)) == 1:
                xtick_labels = copy.deepcopy(self.read_instance.species)
            elif len(np.unique(self.read_instance.species)) == 1:
                xtick_labels = copy.deepcopy(self.read_instance.network)
            else:
                xtick_labels = copy.deepcopy(self.read_instance.networkspecies)
            # get aliases for multispecies (if have any)
            xtick_labels, xlabel = get_multispecies_aliases(xtick_labels)

            # set xlabel if xlabels have changed due to alias, and have one to set
            if xlabel != '':
                plot_characteristics['xlabel']['xlabel'] = xlabel
                relevant_axis.set_xlabel(**plot_characteristics['xlabel'])
        # labels for standard plot
        else:
            xticks = positions
            xtick_labels = copy.deepcopy(cut_data_labels)

        # modify xticks to be horizontal as just have 1 label
        if len(xtick_labels) == 1:
            xtick_params['rotation'] = 0
            xticklabel_params = {}

        # set xticks / xticklabels
        relevant_axis.set_xticks(xticks)
        relevant_axis.xaxis.set_tick_params(**xtick_params)
        relevant_axis.set_xticklabels(xtick_labels, **xticklabel_params)

    def make_heatmap(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                     zstat=None, subsection=None, plotting_paradigm=None, stats_df=None):
        """
        Renders a statistical heatmap using Seaborn to visualise performance metrics across observations and models.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstat : str, optional
            Statistic.
        subsection : str, optional
            Currently active subsection.
        plotting_paradigm : str, optional
            Plotting paradigm (summary or station report).
        stats_df : pandas dataframe, optional
            Dataframe of previously calculated statistics.
        """

        # bias plot?
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if statistical dataframe is not provided then create it
        if not isinstance(stats_df, pd.DataFrame):

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # calculate statistics
            if bias:
                if self.read_instance.observations_data_label in cut_data_labels:
                    cut_data_labels.remove(self.read_instance.observations_data_label)
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, [zstat], 
                                                [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, [zstat], 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float64)

        # get subsections
        subsections = list(np.unique(stats_df.index.get_level_values(1)))

        # get relevant data
        if plotting_paradigm == 'station':
            stats_df = stats_df.iloc[stats_df.index.get_level_values('subsections') == subsection]
        if 'multispecies' not in plot_options:
            stats_df = stats_df.iloc[stats_df.index.get_level_values('networkspecies') == networkspeci]
        else:
            # replace subsection name by networkspecies if there is only one
            if (len(subsections) == 1) or (plotting_paradigm == 'station'):
                stats_df = stats_df.droplevel(level='subsections')

        # determine if want to add annotations or not from plot_options
        if 'annotate' in plot_options:
            # get rounded labels
            decimal_places = plot_characteristics['round_decimal_places']['table']
            if Version(pd.__version__) >= Version("2.1.0"):
                annotate = stats_df.map(lambda x: round_decimal_places(x, decimal_places))
            else:
                annotate = stats_df.applymap(lambda x: round_decimal_places(x, decimal_places))
        else:
            annotate = False

        # plot heatmap
        heatmap = sns.heatmap(stats_df, 
                              ax=relevant_axis, 
                              annot=annotate,
                              fmt='',
                              **plot_characteristics['plot'])
        
        # remove networkspecies-subsections label from y-axis
        relevant_axis.set_ylabel("")

        # if there is only one subsection or station data
        if (plotting_paradigm == 'station') or (len(subsections) == 1):
            # for multispecies, remove network names from labels
            if ('multispecies' in plot_options) and (not plot_characteristics['multispecies']['network_names']):
                if not plot_characteristics['multispecies']['network_names']:
                    yticklabels = [networkspeci.split('|')[1] if '|' in networkspeci else networkspeci 
                                   for networkspeci in stats_df.index]
            # for non multispecies, remove subsection names from labels
            elif ('multispecies' not in plot_options) and (not plot_characteristics['parent_section_names']):
                yticklabels = []
                for subsection_label in stats_df.index.get_level_values(1):
                    if "·" in subsection_label:
                        subsection_label = subsection_label.split('·')[1]
                    yticklabels.append(subsection_label)
            # keep original labels
            else:
                yticklabels = stats_df.index
        # if there is summary data for more than one subsection
        elif (plotting_paradigm == 'summary') and (len(subsections) > 1):
            # remove parent names from subsections
            if not plot_characteristics['parent_section_names']:
                yticklabels = []
                for subsection_label in stats_df.index.get_level_values(1):
                    if "·" in subsection_label:
                        subsection_label = subsection_label.split('·')[1]
                    yticklabels.append(subsection_label)
            # keep original labels
            else:
                yticklabels = stats_df.index.get_level_values(1)
        relevant_axis.set_yticklabels(yticklabels, **plot_characteristics['yticklabels'])

        # set xticklabels
        relevant_axis.set_xticklabels(stats_df.columns, **plot_characteristics['xticklabels'])

        # axis cuts off due to bug in matplotlib 3.1.1 - hack fix
        if Version(matplotlib.__version__) <= Version("3.1.1"):
            if len(stats_df.index) > 1:
                bottom, top = relevant_axis.get_ylim()
                relevant_axis.set_ylim(bottom + 0.5, top - 0.5)

        networkspecies = list(stats_df.index.get_level_values(0)[::(len(subsections))])
        n_rows = len(subsections)*len(networkspecies)
        n_cols = len(data_labels)

        # format for multispecies
        if 'multispecies' in plot_options:
            
            if plot_characteristics['multispecies']['xmin']:
                xmin = plot_characteristics['multispecies']['xmin']
            else:
                if n_rows == n_cols:
                    xmin = -0.35 * n_rows - 0.35 * n_cols
                else:
                    xmin = -0.35 * n_rows + 0.1 * n_cols

            # if we have more than one subsection and we are plotting summaries
            if (len(subsections) > 1) and (plotting_paradigm == 'summary'):
                # add horizontal lines to separate networkspecies
                y_separators = []
                for networkspeci_ii in range(len(networkspecies)+1):
                    y_separators.append(len(subsections)*networkspeci_ii)
                relevant_axis.hlines(y=y_separators, xmin=xmin, 
                                     xmax=0, clip_on=False, **plot_characteristics['multispecies']['hlines'])

                # annotate networkspecies names on the left
                for networkspeci, y_separator in zip(networkspecies, y_separators[:-1]):
                    y_position = (len(subsections)/2) + y_separator
                    # remove network names from networkspecies
                    if not plot_characteristics['multispecies']['network_names']:
                        networkspeci_label = networkspeci.split('|')[1]
                    else:
                        networkspeci_label = networkspeci
                    
                    if Version(matplotlib.__version__) >= Version("3.3"):
                        relevant_axis.annotate(text=networkspeci_label, annotation_clip=False,
                                               xy=(xmin, y_position), 
                                               **plot_characteristics['multispecies']['yticklabels'])
                    else:
                        relevant_axis.annotate(s=networkspeci_label, annotation_clip=False,
                                               xy=(xmin, y_position), 
                                               **plot_characteristics['multispecies']['yticklabels'])        
        
        # format for non multispecies
        else:

            # vertically align yticklabels due to bug again in matplotlib - hack fix. Remove in Future!
            for tick in relevant_axis.get_yticklabels():
                tick.set_verticalalignment("center")

        # track plot elements if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            self.track_plot_elements(self.read_instance.observations_data_label, 'heatmap', 'plot', heatmap, bias=bias)

    def make_table(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options,
                   zstats=None, statsummary=False, subsection=None, plotting_paradigm=None, stats_df=None):
        """
        Constructs a formatted table of statistical metrics, featuring merged cells and dynamic colour-coding.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstats : list, optional
            Statistics.
        statsummary : boolean, optional
            To indicate if making alternative statistical summary table plot.
        subsection : str, optional
            Currently active subsection.
        plotting_paradigm : str, optional
            Plotting paradigm (summary or station report).
        stats_df : pandas dataframe, optional
            Dataframe of previously calculated statistics.
        """

        # turn off axis to make table
        relevant_axis.axis('off')

        # bias plot?
        if 'bias' in plot_options:
            bias = True
        else:
            bias = False

        # if statistical dataframe is not provided then create it
        if not isinstance(stats_df, pd.DataFrame):

            # get valid data labels for networkspeci
            valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

            # cut data_labels for those in valid data labels
            cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

            # calculate statistics
            if bias:
                if self.read_instance.observations_data_label in cut_data_labels:
                    cut_data_labels.remove(self.read_instance.observations_data_label)
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                [self.read_instance.observations_data_label]*len(cut_data_labels), 
                                                cut_data_labels)
            else:
                stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstats, 
                                                 cut_data_labels, [])

            # create stats dataframe
            stats_df = pd.DataFrame(data=stats_calc, 
                                    index=cut_data_labels,
                                    dtype=np.float64)
        
        # when we have 1 stat in the statsummary, the column name is 0
        # we need to rename it to the stat name
        if (len(stats_df.columns) == 1) and (stats_df.columns[0] == 0):
            stats_df = stats_df.rename(columns={stats_df.columns[0]: zstats[0]})
        
        # rename columns to save space
        columns = {}
        for column in stats_df.columns:
            new_colum_name = copy.deepcopy(column)
            # rename columns to replace Diurnal, Weekly and Monthly by D, W, M
            if 'diurnal' in column:
                new_colum_name = column.replace('diurnal', 'D')
            elif 'weekly' in column:
                new_colum_name = column.replace('weekly', 'W')
            elif 'monthly' in column:
                new_colum_name = column.replace('monthly', 'M')
            # remove _bias from columns
            if '_bias' in new_colum_name:
                new_colum_name = new_colum_name.replace('_bias', '')   
            columns[column] = new_colum_name
        stats_df = stats_df.rename(columns=columns)
  
        # get column and row labels
        col_labels = stats_df.columns.tolist()
        row_labels = stats_df.index.tolist()

        # round dataframe
        decimal_places = plot_characteristics['round_decimal_places']['table']
        if Version(pd.__version__) >= Version("2.1.0"):
            stats_df = stats_df.map(lambda x: round_decimal_places(x, decimal_places))
        else: 
            stats_df = stats_df.applymap(lambda x: round_decimal_places(x, decimal_places))

        # reports
        if self.read_instance.mode in ['report', 'library']:
            
            # get relevant data
            if 'multispecies' not in plot_options:
                stats_df = stats_df.iloc[stats_df.index.get_level_values('networkspecies') == networkspeci]
            if plotting_paradigm == 'station':
                stats_df = stats_df.iloc[stats_df.index.get_level_values('subsections') == subsection]

            # get labels
            networkspecies = list(stats_df.index.get_level_values('networkspecies'))
            subsections = list(stats_df.index.get_level_values('subsections'))
            if statsummary:
                data_labels = list(stats_df.index.get_level_values('labels'))
                stats = list(stats_df.columns)
            else:
                data_labels = list(stats_df.columns)

            # reset index after filtering
            stats_df = stats_df.reset_index()

            # hide subsections from station plots or if there is only 1 section
            if self.read_instance.mode in ['report', 'library']:
                if plotting_paradigm == 'station' or len(np.unique(subsections)) == 1:
                    stats_df = stats_df.drop(columns='subsections')
        
            # hide networkspecies from plots that are not multispecies
            if 'multispecies' not in plot_options:
                stats_df = stats_df.drop(columns='networkspecies')
        
            # remove parent names from subsections
            if ('subsections' in stats_df.columns) and (not plot_characteristics['parent_section_names']):
                stats_df['subsections'] = [subsection_label.split('·')[1] if '·' in subsection_label 
                                           else subsection_label for subsection_label in subsections]

            # remove network names from networkspecies
            if (('multispecies' in plot_options) and ('networkspecies' in stats_df.columns) and 
                (not plot_characteristics['multispecies']['network_names'])):
                stats_df['networkspecies'] = [networkspeci_label.split('|')[1] for networkspeci_label in networkspecies]
            
            # get number of "empty" cells (without stats) and 
            # column labels (hide networkspecies, subsections and data labels)
            if statsummary:
                empty_cells = len(stats_df.columns) - len(stats)
                col_labels = ['']*empty_cells + stats
            else:
                empty_cells = len(stats_df.columns) - len(data_labels)
                col_labels = ['']*empty_cells + data_labels

        # dashboard
        else:
            # there is only statsummary
            if statsummary:

                # get labels
                data_labels = list(stats_df.index)
                stats = list(stats_df.columns)
                
                # reset index
                stats_df = stats_df.reset_index()

                # get number of "empty" cells (without stats)
                empty_cells = 1
                col_labels = ['']*empty_cells + stats

        # set cell colors
        if statsummary:
            if 'cell_colours' in plot_characteristics:
                if plot_characteristics['cell_colours']:
                    cell_colours = [[]] * (stats_df.shape[1])
                    for col in range(stats_df.shape[1]):
                        # custom colors for data labels cells
                        if col == (empty_cells-1):
                            for data_label in data_labels:
                                # observations in white
                                if data_label == self.read_instance.observations_data_label:
                                    color = 'white'
                                # models in legend colors
                                else:
                                    color = self.read_instance.plotting_params[data_label]['colour']
                                cell_colours[col].append(color)
                        # white for other cells
                        else:
                            cell_colours[col] = ['white'] * stats_df.shape[0]
                    if stats_df.shape == (1, 1):
                        plot_characteristics['plot']['cellColours'] = np.array(cell_colours, dtype=object)
                    else:
                        plot_characteristics['plot']['cellColours'] = np.array(cell_colours, dtype=object).T
        else:
            if 'col_colours' in plot_characteristics:
                if plot_characteristics['col_colours']:
                    col_colours = []
                    for data_label in data_labels:
                        # observations in white
                        if data_label == self.read_instance.observations_data_label:
                            color = 'white'
                        # models in legend colors
                        else:
                            color = self.read_instance.plotting_params[data_label]['colour']
                        col_colours.extend([color])
                    plot_characteristics['plot']['colColours'] = ['white']*empty_cells + col_colours

        # make table
        table = relevant_axis.table(cellText=stats_df.values, 
                                    colLabels=col_labels,
                                    **plot_characteristics['plot'])

        # merge cells in networkspecies and subsections columns (if any)
        if self.read_instance.mode in ['report', 'library']:
            column_ii = 0
            for column, rows in zip(['networkspecies', 'subsections'], (networkspecies, subsections)):
                if column in stats_df.columns:
                    # count consecutive duplicates
                    count_dups = [sum(1 for _ in group) for _, group in groupby(rows)]
                    
                    # merge cells that have consecutive duplicates
                    current_row = 0
                    for count_ii, count in enumerate(count_dups):
                        cells_to_merge = [(current_row + i, column_ii) for i in range(1, count+1)]
                        merge_cells(table, cells_to_merge)
                        current_row += count
                    column_ii += 1

        # adjust cell height
        if 'cell_height' in plot_characteristics:
            table.scale(1, plot_characteristics['cell_height'])

        # adjust fontsize
        if 'fontsize' in plot_characteristics:
            table.auto_set_font_size(False)
            table.set_fontsize(plot_characteristics['fontsize'])
            table.auto_set_column_width(np.arange(-1, len(col_labels)+1))

        # track plot elements if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            if statsummary:
                self.track_plot_elements(self.read_instance.observations_data_label, 'statsummary', 'plot', [table], bias=bias)
            else:
                self.track_plot_elements(self.read_instance.observations_data_label, 'table', 'plot', [table], bias=bias)
    
    def make_taylor(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options, zstat,
                    stddev_max=None):
        """
        Renders a Taylor diagram to evaluate model performance based on correlation, root-mean-square error, and standard deviation.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        zstat : str
            Statistic.
        stddev_max : float, optional
            Maximum standard deviation.

        Returns
        -------
        bool
            Returns True upon successful completion of the plot.
        """

        if self.read_instance.mode in ['report', 'library']:
            self.taylor_polar_relevant_axis = relevant_axis.get_aux_axes(
                PolarAxes.PolarTransform(apply_theta_transforms=False))

        # add observations to labels to get all standard deviation when plotting per label
        if 'individual' in plot_options:
            data_labels.insert(0, self.read_instance.observations_data_label)

        # calculate statistics
        stats_dict = {}

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # get data labels without observations
        obs_index = valid_data_labels.index(self.read_instance.observations_data_label)
        data_labels_sans_obs = copy.deepcopy(cut_data_labels)
        if self.read_instance.observations_data_label in data_labels_sans_obs:
            data_labels_sans_obs.remove(self.read_instance.observations_data_label)

        # standard deviation - absolute calculations of observations and models
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, 'StdDev', 
                                         data_labels, [])
        stats_dict['StdDev'] = stats_calc   

        # correlation - between observations and model
        stats_calc = calculate_statistic(self.read_instance, self.canvas_instance, networkspeci, zstat, 
                                         [self.read_instance.observations_data_label]*len(data_labels_sans_obs), 
                                         data_labels_sans_obs)
        stats_calc = np.insert(stats_calc, obs_index, np.nan)
        stats_dict[zstat] = stats_calc
        
        # get maximum stddev in dataframe (if not defined)
        if stddev_max is None:
            stddev_max = np.nanmax(stats_dict["StdDev"])

        # create stats dataframe
        stats_df = pd.DataFrame(data=stats_dict, 
                                index=data_labels,
                                dtype=np.float64)

        # get labels
        rlabel = stats_df.columns[1]
        xylabel = stats_df.columns[0]

        # check if we need to extend the axis to negative correlations
        extend = False
        if np.nanmin(stats_df[rlabel]) < 0:
            extend = True

        # update axis extremes and labels
        tmin, tmax, smin, smax, gl1, tf1 = get_taylor_diagram_ghelper_info(stddev_max, 
                                                                           plot_characteristics,
                                                                           extend)
        relevant_axis.get_grid_helper().update_grid_finder(
            extreme_finder=fa.ExtremeFinderFixed((tmin, tmax, smin, smax)),
            grid_locator1=gl1, tick_formatter1=tf1)

        # update axis position and size in dashboard
        if self.read_instance.mode not in ['report', 'library']:

            # find Taylor plot position in layout
            for plot_position in range(2, 6):
                plot_type = getattr(self.read_instance, 'position_{}'.format(plot_position))
                if plot_type == 'taylor':
                    break
            
            # changing the extend reduces the size of the plot and changes its start position
            if extend:
                old_position = relevant_axis.get_position().bounds
                if plot_position == 2:
                    new_position = (0.60, 0.42, 0.288, 0.594)
                elif plot_position == 3:
                    new_position = (0.03, 0, 0.256, 0.56)
                elif plot_position == 4:
                    new_position = (0.37, 0, 0.256, 0.56)
                elif plot_position == 5:
                    new_position = (0.70, 0, 0.256, 0.56)
            else:
                if plot_position == 2:
                    new_position = (0.64, 0.57, 0.16, 0.33)
                elif plot_position == 3:
                    new_position = (0.08, 0.08, 0.16, 0.35)
                elif plot_position == 4:
                    new_position = (0.41, 0.08, 0.16, 0.35)
                elif plot_position == 5:
                    new_position = (0.69, 0.08, 0.16, 0.35)
            relevant_axis.set_position(new_position)

        # clear axis, add grid and adjust limits 
        # as suggested by the Matpotlib devs in https://github.com/matplotlib/matplotlib/issues/25426
        relevant_axis.clear()
        relevant_axis.grid(**plot_characteristics['grid'])
        relevant_axis.adjust_axes_lim()

        # adjust top axis (curve)
        relevant_axis.axis['top'].set_axis_direction('bottom')
        relevant_axis.axis['top'].toggle(ticklabels=True, label=True)
        relevant_axis.axis['top'].major_ticklabels.set_axis_direction('top')
        relevant_axis.axis['top'].major_ticklabels.set(**plot_characteristics['rtick_params'])
        relevant_axis.axis['top'].label.set_text(rlabel)
        relevant_axis.axis['top'].label.set_fontsize(plot_characteristics['rlabel']['fontsize'])
        relevant_axis.axis['top'].label.set_axis_direction('top')

        # adjust right axis (y axis)
        relevant_axis.axis['right'].set_axis_direction('top')
        relevant_axis.axis['right'].toggle(ticklabels=True)
        relevant_axis.axis['right'].major_ticklabels.set_axis_direction('bottom' if extend else 'left')
        relevant_axis.axis['right'].major_ticklabels.set(**plot_characteristics['xytick_params'])

        # adjust left axis (x axis)
        relevant_axis.axis['left'].set_axis_direction('bottom')
        relevant_axis.axis['left'].major_ticklabels.set(**plot_characteristics['xytick_params'])

        # hide bottom axis ticks and tick labels
        relevant_axis.axis['bottom'].major_ticklabels.set_visible(False)
        relevant_axis.axis['bottom'].major_ticks.set_visible(False)

        # show label
        if extend:
            relevant_axis.axis['bottom'].label.set_text(xylabel)
            relevant_axis.axis['bottom'].label.set_pad(10)
            relevant_axis.axis['bottom'].label.set_fontsize(plot_characteristics['xylabel']['fontsize'])
        else:
            relevant_axis.axis['left'].label.set_text(xylabel)
            relevant_axis.axis['left'].label.set_fontsize(plot_characteristics['xylabel']['fontsize'])  

        # add contours around observations standard deviation
        reference_stddev = stats_df[xylabel].iloc[obs_index]
        num_levels = plot_characteristics['contours']['levels']['number']
        rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0, tmax))
        rms = np.sqrt(reference_stddev**2 + rs**2 - 2*reference_stddev*rs*np.cos(ts))
        contours = self.taylor_polar_relevant_axis.contour(ts, rs, rms, num_levels,
            **plot_characteristics['contours']['style']['general'])

        # add contour labels
        self.taylor_polar_relevant_axis.clabel(contours, contours.levels, inline=True, fmt = "%.2f", fontsize=6)

        # add reference contour of observational standard deviation
        ref_x = np.linspace(0, tmax)
        ref_y = np.zeros_like(ref_x) + reference_stddev
        self.taylor_polar_relevant_axis.plot(ref_x, ref_y, **plot_characteristics['contours']['style']['obs'])

        # add models
        for data_label, stddev, corr_stat in zip(stats_df.index, 
                                                 stats_df[xylabel], 
                                                 stats_df[rlabel]):
            if data_label == self.read_instance.observations_data_label:
                continue
            self.taylor_plot = self.taylor_polar_relevant_axis.plot(np.arccos(corr_stat), stddev,
                                                                    **plot_characteristics['plot'],
                                                                    mfc=self.read_instance.plotting_params[data_label]['colour'], 
                                                                    mec=self.read_instance.plotting_params[data_label]['colour'],
                                                                    label=data_label) 

            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements(data_label, 'taylor', 'plot', self.taylor_plot, bias=False)

        return True
    
    def make_fairmode_target(self, relevant_axis, networkspeci, data_labels, plot_characteristics, plot_options):
        """
        Renders a FAIRMODE Target plot to assess model quality objectives (MQO) and uncertainty.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """

        # resample to daily for PM10 and PM2.5 if data is hourly
        # get MDA8 for ozone if data is hourly
        # finally filter by coverage
        data, valid_station_idxs = get_fairmode_data(self.read_instance, self.canvas_instance, networkspeci, data_labels)
        
        # skip making plot if there is no valid data
        # library and report modes are already handling this in advance
        if (self.read_instance.mode not in ['report', 'library']) and (not any(valid_station_idxs)):
            msg = 'No valid data to create FAIRMODE target plot after filtering by coverage.'
            show_message(self.read_instance, msg)
            self.read_instance.handle_layout_update('None', sender=self.canvas_instance.get_plot_type_position('fairmode-target'))
            return
        
        observations_data = data[0, :, :]

        # get settings
        speci = networkspeci.split('|')[1]
        u_95r_RV = fairmode_settings[speci].get('u_95r_RV')
        RV = fairmode_settings[speci].get('RV')
        alpha = fairmode_settings[speci].get('alpha')
        beta = fairmode_settings[speci].get('beta')
        coverage = fairmode_settings[speci].get('coverage')
        exc_threshold = fairmode_settings[speci].get('exc_threshold')
        percentile = fairmode_settings[speci].get('percentile')
        units = fairmode_settings[speci].get('units')

        # get RV and exceedance threshold per units
        RV, exc_threshold = get_fairmode_RV_exceendance(self.read_instance, speci, RV, exc_threshold, units)
            
        # add target
        main_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['main'])
        relevant_axis.add_patch(main_circle)

        # add a black circle with radius 1 (continuous line)
        big_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['big'])
        relevant_axis.add_patch(big_circle)

        # add a black circle with radius 0.5 (dotted line)
        small_circle = plt.Circle(**plot_characteristics['auxiliar']['circle']['small'])
        relevant_axis.add_patch(small_circle)

        # add text in the sides
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['top'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['bottom'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['left'], transform=relevant_axis.transAxes)
        relevant_axis.text(**plot_characteristics['auxiliar']['sides']['right'], transform=relevant_axis.transAxes)

        # add diagonal lines (y = x and y = -x)
        xmin = np.min(plot_characteristics['xticks']['ticks'])
        xmax = np.max(plot_characteristics['xticks']['ticks'])
        ymin = np.min(plot_characteristics['yticks']['ticks'])
        ymax = np.max(plot_characteristics['yticks']['ticks'])
        relevant_axis.plot([xmin, xmax], [ymin, ymax], **plot_characteristics['auxiliar']['crosses']['increasing'])
        relevant_axis.plot([xmin, xmax], [ymax, ymin], **plot_characteristics['auxiliar']['crosses']['decreasing'])

        # get metadata without nans
        classification_type = plot_characteristics['markers']['type'].lower()
        valid_station_references = get_valid_metadata(self, 'station_reference', 
                                                      valid_station_idxs, networkspeci)
        try:
            valid_station_classifications = get_valid_metadata(self, f'{classification_type}_classification', 
                                                               valid_station_idxs, networkspeci)
        except:
            valid_station_classifications = np.full(len(valid_station_references), np.nan, dtype=np.float32)
            self.read_instance.logger.info(f'Data for {classification_type}_classification is not available and will not be shown in the legend')

        # get number of stations
        n_stations = len(valid_station_references)

        # initialise annotation text
        self.faimode_target_annotate_text = []
        self.faimode_target_annotate_colour = []
        if "parameters" in plot_characteristics['annotate_options']:
            self.faimode_target_annotate_text.append(f"α={alpha}")
            self.faimode_target_annotate_text.append(f"\nβ={beta}")
            self.faimode_target_annotate_text.append(f"\nRV={RV}")
            self.faimode_target_annotate_text.append(f"\nU₉₅,ᵣᴿⱽ={u_95r_RV}")
            self.faimode_target_annotate_colour.extend(['black']*4)
        if "covered_stations" in plot_characteristics['annotate_options']:
            self.faimode_target_annotate_text.append(f"\n\n{n_stations} stations with\ncoverage above {coverage}%")
            self.faimode_target_annotate_colour.append('black')

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # iterate through data labels
        for data_label in cut_data_labels:

            # continue for observations data label
            if data_label == self.read_instance.observations_data_label:
                continue

            # get model data
            model_data = data[valid_data_labels.index(data_label), :, :]
            
            # calculate MQI for the current station
            x_points = []
            y_points = []
            stations = []
            bad_stations = []
            classifications = []

            # get FAIRMODE statistics per station
            mqi_array = np.full(n_stations, np.nan)
         
            for station_idx, (station, classification) in enumerate(
                zip(valid_station_references, valid_station_classifications)):

                st_observations_data = observations_data[station_idx, :]
                st_model_data = model_data[station_idx, :]
      
                x, y, mqi = ModBias.calculate_fairmode_stats(
                    st_observations_data, st_model_data, 
                    u_95r_RV, RV, alpha, beta, exc_threshold, percentile, 'target')

                x_points.append(x)
                y_points.append(y)
                mqi_array[station_idx] = mqi
                stations.append(station)
                classifications.append(classification)
                if mqi > 1:
                    bad_stations.append(station)

            # create list for track_plot_elements
            self.fairmode_target_plot = []

            # plot data
            # we need to create the plot point by point to be able to set the marker
            # depending on the area classification since Matplotlib doesn't have a way to 
            # set different markers at the same time
            for x, y, classification in (zip(x_points, y_points, classifications)):
                if classification not in plot_characteristics['markers'][f'{classification_type}_classification'].keys():
                    marker = 'h'
                else:
                    marker = plot_characteristics['markers'][f'{classification_type}_classification'][classification]
                stations_dots = relevant_axis.plot(x, y, markeredgecolor=self.read_instance.plotting_params[data_label]['colour'], 
                                   marker=marker, **plot_characteristics['plot'])
                self.fairmode_target_plot.append(stations_dots[0])

            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements(data_label, 'fairmode-target', 'plot', self.fairmode_target_plot, bias=False)

            # add MQI90
            self.faimode_target_annotate_text.append(f"\n\n{data_label}")
            self.faimode_target_annotate_colour.append('black')
            if "MQI90" in plot_characteristics['annotate_options']:
                # calculate MQI90
                mqi_sorted = sorted(mqi_array[~np.isnan(mqi_array)])
                i_90 = int(0.9 * len(mqi_sorted)) - 1
                MQI90 = mqi_sorted[i_90]
                MQI90_formatted = '{0:.{1}f}'.format(
                    MQI90, plot_characteristics['annotate_text']['round_decimal_places'])
                self.faimode_target_annotate_text.append(f"\nMQI₉₀ = {MQI90_formatted}")
                if MQI90 > 1:
                    mqi_color = 'red'
                else:
                    mqi_color = 'green'
                self.faimode_target_annotate_colour.append(mqi_color)

            # add bad stations
            if "N_bad_stations" in plot_characteristics['annotate_options']:
                if len(bad_stations) == 1:
                    stations_name = 'station'
                else:
                    stations_name = 'stations'
                self.faimode_target_annotate_text.append(f'\n{len(bad_stations)} {stations_name} with MQI > 1')
                self.faimode_target_annotate_colour.append('black')
            if 'bad_stations' in plot_characteristics['annotate_options']:
                if (len(bad_stations) > 0):
                    formatted_stations = "\n".join(station.replace('_', '') for station in bad_stations)
                    self.faimode_target_annotate_text.append(f':\n{formatted_stations}')
                    self.faimode_target_annotate_colour.append('black')
                if data_label != cut_data_labels[-1]:
                    self.faimode_target_annotate_text += '\n\n'
                    self.faimode_target_annotate_colour.append('black')

        # strip empty characters from start
        self.faimode_target_annotate_text[0] = self.faimode_target_annotate_text[0].strip() 

        # create legend
        legend_elements = []
        for classification in np.unique(classifications):
            if classification not in plot_characteristics['markers'][f'{classification_type}_classification']:
                marker = 'h'
            else:
                marker = plot_characteristics['markers'][f'{classification_type}_classification'][classification]
            legend_element = mlines.Line2D([], [], marker=marker, 
                                           label=classification, 
                                           **plot_characteristics['markers']['plot'])
            legend_elements.append(legend_element)
        relevant_axis.legend(handles=legend_elements, 
                             **plot_characteristics['markers']['legend'])

        # add title if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            set_axis_title(self.read_instance, relevant_axis, fairmode_settings[speci]['title'], 
                           plot_characteristics)

    def make_fairmode_statsummary(self, relevant_axis, networkspeci, data_labels, plot_characteristics, 
                                  plot_options):
        """
        Renders a FAIRMODE statistical summary plot to provide a detailed breakdown of model performance indicators.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """

        # resample to daily for PM10 and PM2.5 if data is hourly
        # get MDA8 for ozone if data is hourly
        # finally filter by coverage
        data, valid_station_idxs = get_fairmode_data(self.read_instance, self.canvas_instance, networkspeci, data_labels)

        # skip making plot if there is no valid data
        # library and report modes are already handling this in advance
        if (self.read_instance.mode not in ['report', 'library']) and (not any(valid_station_idxs)):
            msg = 'No valid data to create FAIRMODE statistic summary plot after filtering by coverage.'
            show_message(self.read_instance, msg)
            self.read_instance.handle_layout_update('None', sender=self.canvas_instance.get_plot_type_position('fairmode-statsummary'))
            return

        observations_data = data[0, :, :]

        # get settings
        speci = networkspeci.split('|')[1]
        u_95r_RV = fairmode_settings[speci].get('u_95r_RV')
        RV = fairmode_settings[speci].get('RV')
        alpha = fairmode_settings[speci].get('alpha')
        beta = fairmode_settings[speci].get('beta')
        exc_threshold = fairmode_settings[speci].get('exc_threshold')
        percentile = fairmode_settings[speci].get('percentile')
        units = fairmode_settings[speci].get('units')

        # get RV and exceedance threshold per units
        RV, exc_threshold = get_fairmode_RV_exceendance(self.read_instance, speci, RV, exc_threshold, units)
        
        # get station references
        valid_station_references = get_valid_metadata(self, 'station_reference', 
                                                      valid_station_idxs, networkspeci)
        
        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # iterate through data labels
        for data_label in cut_data_labels:

            # get model data
            model_data = data[valid_data_labels.index(data_label), :, :]
            
            # calculate MQI for the current station
            exceedances = []
            means = []
            t_biases = []
            t_R_list = []
            t_sd_list = []
            h_perc_list = []
        
            for station_idx, station in enumerate(valid_station_references):

                st_observations_data = observations_data[station_idx, :]
                st_model_data = model_data[station_idx, :]

                mean, exc, t_bias, t_R, t_sd, h_perc = ModBias.calculate_fairmode_stats(
                    st_observations_data, st_model_data, 
                    u_95r_RV, RV, alpha, beta, exc_threshold, percentile, 'summary')

                means.append(mean)
                exceedances.append(exc)
                t_biases.append(t_bias)
                t_R_list.append(t_R)
                t_sd_list.append(t_sd)
                h_perc_list.append(h_perc)

            # plot data
            # join all statistics in one list
            statistics_list = [np.array(means), np.array(exceedances), np.array(t_biases), 
                               np.array(t_R_list), np.array(t_sd_list), np.array(h_perc_list), 
                               np.nanmean(t_sd_list), np.nanmean(t_R_list)]

            # get subplots dictionary
            subplots = dict(plot_characteristics["auxiliar"]["subplots"])

            # get the variable that tells you if there are exceedances in this species
            has_exceedances = exc_threshold != None

            # if there is no threshold don't create the exceedances row
            if not has_exceedances:
                subplots.pop("observed exceedances", None)
                statistics_list.pop(1)

            # create list for track_plot_elements
            self.fairmode_statsummary_plot = []
            
            # apply configuration to each row
            for i, (row, fairmode_data) in enumerate(zip(subplots,statistics_list)):

                # for the two firsts rows, skip the models
                if i < 2 and data_label != self.read_instance.observations_data_label:
                    continue

                # skip the observations all the rows but the two first ones
                if i >= 2 and data_label == self.read_instance.observations_data_label:
                    continue

                # get row dictionary
                plot_dict = subplots[row]

                # remove axis from the dot on right side
                relevant_axis[i*4 + 3].set_xticks([])

                # remove the axis of the dot
                for side in ['bottom', 'top', 'right']:
                    relevant_axis[i*4 + 3].spines[side].set_linestyle('none')                  
                
                # add dashed line on the left
                relevant_axis[i*4 + 3].spines['left'].set_linestyle((10, (8, 5)))                  
                
                # configure color of the row and the dot on the right for rows 3 to 8
                if 'range_style' in plot_dict:
                    # get the dictionary with the style of the row and its range
                    range_style_dict = plot_characteristics["auxiliar"]["range_style"][plot_dict['range_style']]
                    for span, color in zip(range_style_dict["spans"],range_style_dict["colors"]):
                        relevant_axis[i*4 + 1].axvspan(*span, color=color, lw=0)

                    # dot on the right configuration
                    # get the lowest and highest number on the range
                    min_span = range_style_dict["spans"][0][0]
                    max_span = range_style_dict["spans"][-1][-1]

                    # get the color of the dot on the right
                    arr = np.array(fairmode_data)[~np.isnan(fairmode_data)]
                    correct_arr = arr[(arr >= min_span) & (arr <= max_span)]
                    
                    dot_color = plot_characteristics["auxiliar"]["right_dot_colors"]['red']
                    if len(arr) > 0 and len(correct_arr)/len(arr) >= .9:
                        dot_color = plot_characteristics["auxiliar"]["right_dot_colors"]["green"]
                    
                    # plot dot on the right
                    relevant_axis[i*4 + 3].scatter(**plot_characteristics["auxiliar"]["right_dot"], 
                                                   color=dot_color, edgecolor=dot_color)

                # y axis / grid
                # remove y axis ticks
                for j in range(4):
                    relevant_axis[i*4 + j].set_yticks([])
                    relevant_axis[i*4 + j].grid(False)

                # x axis
                # get the x axis limit for the current row
                x_limit = plot_dict['x_axis_limits']
                
                # set x axis limits
                relevant_axis[i*4 + 1].set_xlim(*x_limit)

                # right zone configuration
                # remove x ticks from the right dashed zone
                relevant_axis[i*4 + 2].set_xticks([])
                
                # remove vertical lines separating middle and right dashed zone
                relevant_axis[i*4 + 2].spines['left'].set_linestyle('none')
                relevant_axis[i*4 + 1].spines['right'].set_linestyle('none')

                # change the linestyle to dashed
                for side in ['bottom', 'top', 'right']:
                    relevant_axis[i*4 + 2].spines[side].set_linestyle((10, (8, 5)))

                # right zone dot configuration
                # get the points that surpass the limit
                right_zone_mask = x_limit[1] < fairmode_data
                
                # set the right dashed zone range to (-1,1)
                relevant_axis[i*4 + 2].set_xlim(-1,1)
                
                # if there is a dot outside the limits
                if np.any(right_zone_mask):
                    # plot it in the middle of the right dashed zone
                    stations_dots = relevant_axis[i*4 + 2].plot(
                        0, 0, 
                        plot_characteristics["auxiliar"]["station_dots"]["marker"], 
                        color=self.read_instance.plotting_params[data_label]['colour'], 
                        markersize=plot_characteristics["auxiliar"]["station_dots"]["markersize"])
                    
                    # add the dots to the track plot elements list
                    self.fairmode_statsummary_plot.append(stations_dots[0])
                    
                    # remove it from the data plotted in the middle zone
                    fairmode_data = fairmode_data[~right_zone_mask] if isinstance(fairmode_data,np.ndarray) else np.array([])
         
                # left zone configuration
                # remove x ticks from the left dashed zone
                relevant_axis[i*4 + 0].set_xticks([])
                
                # define left zone line style (can be dashed or no left zone)
                left_dashed_zone_linestyle = plot_dict['left_dashed_zone_linestyle']
                
                # if left dashed zone exists in the current row
                if left_dashed_zone_linestyle != 'none':
                    # convert to tuple [x,[x,x]] because yaml does not return python tuples
                    left_dashed_zone_linestyle = (left_dashed_zone_linestyle[0],tuple(left_dashed_zone_linestyle[1]))
                    
                    # remove vertical lines separating middle and left dashed zone
                    relevant_axis[i*4 + 1].spines['left'].set_linestyle('none')
                    relevant_axis[i*4 + 0].spines['right'].set_linestyle('none')
                
                # change the linestyle to dashed or remove the dashed zone
                for side in ['bottom', 'top', 'left']:
                    relevant_axis[i*4 + 0].spines[side].set_linestyle(left_dashed_zone_linestyle)

                # left zone dot configuration               
                # if left dashed zone exists in the current row
                if left_dashed_zone_linestyle != 'none':
                    
                    # get the points that surpass the limit
                    left_zone_mask = x_limit[0] > fairmode_data
                    
                    # set the left dashed zone range to (-1,1)
                    relevant_axis[i*4 + 0].set_xlim(-1,1)
                   
                    # if there is a dot outside the limits
                    if np.any(left_zone_mask):
                        # plot it in the middle of the left dashed zone
                        stations_dots = relevant_axis[i*4 + 0].plot(
                            0, 0, 
                            plot_characteristics["auxiliar"]["station_dots"]["marker"], 
                            color=self.read_instance.plotting_params[data_label]['colour'], 
                            markersize=plot_characteristics["auxiliar"]["station_dots"]["markersize"])
                        
                        # add the dot to the track plot elements list
                        self.fairmode_statsummary_plot.append(stations_dots[0])
                        
                        # remove it from the data plotted in the middle zone
                        fairmode_data = fairmode_data[~left_zone_mask] if isinstance(fairmode_data,np.ndarray) else np.array([])

                # plot stations as dots
                stations_dots = relevant_axis[i*4 + 1].plot(
                    fairmode_data, np.zeros_like(fairmode_data), 
                    plot_characteristics["auxiliar"]["station_dots"]["marker"], 
                    color=self.read_instance.plotting_params[data_label]['colour'], 
                    markersize=plot_characteristics["auxiliar"]["station_dots"]["markersize"])
                
                # change the size of the x-axis tick labels
                relevant_axis[i*4 + 1].tick_params(**plot_characteristics["auxiliar"]["station_dots"]["tick_params"])
                
                # add the dots to the track plot elements list
                self.fairmode_statsummary_plot.append(stations_dots[0])
                
            # track plot elements if using dashboard 
            if self.read_instance.mode not in ['report', 'library']:
                self.track_plot_elements(data_label, 'fairmode-statsummary', 'plot', 
                                         self.fairmode_statsummary_plot, bias=False)
                        
        # plot the row titles, descriptions, separators and units
        for i, (row, fairmode_data) in enumerate(zip(subplots,statistics_list)):
            # get row dictionary
            plot_dict = subplots[row]

            # get the row title
            row_title = plot_dict["title"]

            # write the threshold on the exceedances row title
            if row == "observed_exceedances":
                row_title = row_title.format(exc_threshold)

            # add units to the first two rows
            if 'units' in plot_dict:
                relevant_axis[i*4 + 3].text(
                    *plot_characteristics["auxiliar"]["units"]["position"], 
                    plot_dict['units'], 
                    fontsize=plot_characteristics["auxiliar"]["units"]["fontsize"])
                
            # add row title
            relevant_axis[i*4 + 0].text(
                *plot_characteristics["auxiliar"]["row_title"]["position"], 
                row_title, 
                **plot_characteristics["auxiliar"]["row_title"], 
                transform=relevant_axis[i*4 + 0].transAxes)
        
            # add left description
            relevant_axis[i*4 + 0].text(
                *plot_characteristics["auxiliar"]["description"]["position"], 
                plot_dict["description"], 
                **plot_characteristics["auxiliar"]["description"], 
                transform=relevant_axis[i*4 + 0].transAxes)

            # add separator
            if "separator" in plot_dict:
                relevant_axis[i*4 + 0].text(
                    *plot_characteristics["auxiliar"]["separator"]["position"], 
                    plot_characteristics["auxiliar"]["separator_text"], 
                    **plot_characteristics["auxiliar"]["separator"], 
                    transform=relevant_axis[i*4 + 0].transAxes)

        # add title if using dashboard 
        if self.read_instance.mode not in ['report', 'library']:
            set_axis_title(self.read_instance, relevant_axis, fairmode_settings[speci]['title'], plot_characteristics)

    def make_contingencytable(self, relevant_axis, networkspeci, data_labels, plot_characteristics, 
                              plot_options):
        """
        Renders a contingency table or a gerrity score table per station

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        data_labels : list
            Data arrays to plot.
        plot_characteristics : dict
            Plot characteristics.
        plot_options : list
            Options to configure plot.
        """
        
        # get limits
        speci = networkspeci.split('|')[1]
        if speci in list(contingency_settings.keys()):
            limits = contingency_settings[speci]['limits']
            limits_units = contingency_settings[speci]['units']

        # get input and output units
        standard_parameter_speci = get_standard_parameters_by_speci(speci, self.read_instance.ghost_version)
        initial_units = limits_units
        final_units = self.read_instance.measurement_units[speci]

        # convert units using conversion factor
        conversion_factor = get_conversion_factor(initial_units, final_units, standard_parameter_speci) 
        if isinstance(conversion_factor, str):
            self.read_instance.logger.error(conversion_factor)
            sys.exit(1)
        limits = [limit*conversion_factor for limit in limits]
        
        # define categories
        index_levels = plot_characteristics["index_levels"]
        levels = plot_characteristics["levels"]
        edges = plot_characteristics["edges"] 
        gerrity_row_limit = plot_characteristics["gerrity_row_limit"]

        # get valid data labels for networkspeci
        valid_data_labels = self.canvas_instance.selected_station_data_labels[networkspeci]

        # cut data_labels for those in valid data labels
        cut_data_labels = [data_label for data_label in data_labels if data_label in valid_data_labels]

        # check if there is data
        if 'per_station' not in self.canvas_instance.selected_station_data[networkspeci]:
            msg = 'Contingency table cannot be calculated without data.'
            show_message(self.read_instance, msg)
            return
                        
        # get observations data (flattened)
        observations_data = self.canvas_instance.selected_station_data[networkspeci]['per_station'][0,:,:]
        n_stations = observations_data.shape[0]

        # if gerrity score is in plot options or we are in dashboard and we select more than one station, then make gerrity plot
        make_gerrity = True if 'gerrity' in plot_options or (self.read_instance.mode not in ['report', 'library'] and n_stations > 1) else False

        # iterate through data labels
        results = []
        for data_label in cut_data_labels:

            # continue for observations data label
            if data_label == self.read_instance.observations_data_label:
                continue

            # get model data (flattened)# If gerrity score is in plot options or we are in dashboard and we select more than one station
            model_data = self.canvas_instance.selected_station_data[networkspeci]['per_station'][valid_data_labels.index(data_label),:,:]
            
            for station_ind in np.arange(n_stations):

                st_observations_data = observations_data[station_ind, :]
                st_model_data = model_data[station_ind, :]

                # calculate rolling average for PM10 and PM2.5
                if speci in ['pm2p5', 'pm10']:
                    st_observations_data = pd.Series(st_observations_data).rolling(window=24, min_periods=18).mean().to_numpy()
                    st_model_data = pd.Series(st_model_data).rolling(window=24, min_periods=18).mean().to_numpy()

                # calculate contingency table
                contingency_table = ModBias.calculate_contingency_table(st_observations_data, st_model_data, 
                                                                        limits, index_levels, 
                                                                        self.read_instance.time_index, edges)
                
                # calculate gerrity score per station
                if make_gerrity:
                    gerrity_score = round(ModBias.calculate_gerrity_score(contingency_table), 2)
                    current_lon = round(self.canvas_instance.selected_station_metadata[networkspeci]['longitude'][station_ind][0], 2)
                    current_lat = round(self.canvas_instance.selected_station_metadata[networkspeci]['latitude'][station_ind][0], 2)
                    current_station_reference = self.canvas_instance.selected_station_metadata[networkspeci]['station_reference'][station_ind][0]
                    
                    # do not make table / append data if station reference is nan
                    if pd.isna(current_station_reference):
                        continue
                    current_station_name = self.canvas_instance.selected_station_metadata[networkspeci]['station_name'][station_ind][0]
                    
                    results.append({
                        'Station reference': current_station_reference,
                        'Station name': current_station_name,
                        'Latitude': current_lat,
                        'Longitude': current_lon,
                        'Gerrity Score': gerrity_score
                    })

                    # hide last rows if there is a limit
                    if gerrity_row_limit is not None:
                        if station_ind > gerrity_row_limit-1 and n_stations != gerrity_row_limit:
                            results.append({
                                'Station reference': '...',
                                'Station name': '...',
                                'Latitude': '...',
                                'Longitude': '...',
                                'Gerrity Score': '...'
                            })
                            break

                # show contingency table (only available per station)
                else:
                    results_df = pd.DataFrame(contingency_table.table.values, index=levels, columns=levels)

            # make table
            if make_gerrity:
                results_df = pd.DataFrame(results)

            if not results_df.empty:
                table = relevant_axis.table(cellText=results_df.values, 
                                            colLabels=results_df.columns,
                                            rowLabels=None if make_gerrity else results_df.index,
                                            **plot_characteristics['plot'])

                if not make_gerrity:
                    # Show data label on a row above value columns
                    for col in range(len(results_df.columns)):
                        # Center horizontally
                        if col == 3:
                            text = data_label
                        else:
                            text = ""
                        table.add_cell(
                            row=-1,
                            col=col,
                            width=1,
                            height=table[0, 0].get_height(),
                            text=text,
                            loc="center"
                        )
                    cells_to_merge = [(-1, col) for col in range(len(results_df.columns))]
                    merge_cells(table, cells_to_merge, visibility=True)

                    # Show observations label above index column
                    table.add_cell(
                        row=0,
                        col=-1,
                        width=table[0, 0].get_width(),
                        height=table[0, 0].get_height(),
                        text=self.read_instance.observations_data_label,
                        loc="center"
                    )

                # adjust cell height
                if 'cell_height' in plot_characteristics:
                    table.scale(1, plot_characteristics['cell_height'])

                # adjust fontsize
                if 'fontsize' in plot_characteristics:
                    table.auto_set_font_size(False)
                    table.set_fontsize(plot_characteristics['fontsize'])
                    table.auto_set_column_width(np.arange(-1, len(results_df.columns)+1))
                    
                # track plot elements if using dashboard 
                if self.read_instance.mode not in ['report', 'library']:
                    self.track_plot_elements(data_label, 'contingencytable', 'plot', [table], bias=False)
            
    def track_plot_elements(self, data_label, base_plot_type, element_type, plot_object, bias=False):
        """
        Registers plotted artists and collections to manage their visibility dynamically on the dashboard.

        Parameters
        ----------
        data_label : str
            Data array to plot.
        base_plot_type : str
            Plot type, without statistical information.
        element_type : str
            Element type.
        plot_object : object
            Plotted element object.
        bias : boolean, optional
            Indicates if plot is a bias plot.
        """

        # set variable name to access plot elements (absolute or bias versions)
        if not bias:
            plot_element_varname = 'absolute'
        else:
            plot_element_varname = 'bias'

        # add dictionary for plot_type elements if does not yet exist
        if base_plot_type not in self.canvas_instance.plot_elements:
            self.canvas_instance.plot_elements[base_plot_type] =  {'active': 'absolute', 'absolute': {}} 

        # add plot_element_varname if does not yet exist
        if plot_element_varname not in self.canvas_instance.plot_elements[base_plot_type]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname] = {}

        # add dictionary for data label if does not yet exist
        if data_label not in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label] = {}

        # add list for element type if does not yet exist
        if element_type not in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label]:
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] = []
        # if does exist already then remove plot element return, as element has already been plotted
        else:
            return 

        # track plot elements
        # periodic plot specific elements 
        if (base_plot_type in ['periodic', 'periodic-violin']) & (data_label != 'ALL'):
            # add list of collections
            if isinstance(plot_object, dict):
                self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['bodies']
            # add list of lines
            elif isinstance(plot_object, list):
                self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object
        # boxplot plot specific elements
        elif (base_plot_type == 'boxplot') & (data_label != 'ALL'):
            # add lines for all boxplot elements
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['boxes']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['medians']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['whiskers']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['caps']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['fliers']
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object['means']
        # do not save elements for plot objecplot type, without statistical information
        # all other plot elements
        else:
            # add list of lines
            self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type] += plot_object
            
        # set element visibility
        if (data_label not in self.canvas_instance.plot_elements['data_labels_active']) & (data_label != 'ALL'):
            for element in self.canvas_instance.plot_elements[base_plot_type][plot_element_varname][data_label][element_type]:
                element.set_visible(False)

    def get_markersize(self, relevant_axis, base_plot_type, networkspeci, plot_characteristics, 
                       data=None, active_map_valid_station_inds=[]):
        """
        Determines and updates the marker size within plot_characteristics based on data density and plot type.

        Parameters
        ----------
        relevant_axis : object
            Axis to plot on.
        base_plot_type : str
            Plot type, without statistical information.
        networkspeci : str
            Current networkspeci (e.g. EBAS|sconco3).
        plot_characteristics : dict
            Plot characteristics.
        data : numpy array, optional
            Data array to be plotted.
        active_map_valid_station_inds : numpy array, optional
            Valid map indices to plot.
        """

        if base_plot_type in ['timeseries', 'scatter']:
            
            if plot_characteristics['plot']['markersize'] == '':

                # get minimum number of non-NaN data points for plot type across data labels
                min_points = np.count_nonzero(~np.isnan(data))

                # configure size of plots if have very few points
                if min_points < plot_characteristics['markersize_npoints_threshold']:
                    markersize = plot_characteristics['markersize']['few_points'] 
                else:
                    markersize = plot_characteristics['markersize']['standard'] 

                # add to plot_characteristics yaml
                plot_characteristics['plot']['markersize'] = markersize
        
        elif base_plot_type == 'map':

            if plot_characteristics['plot']['s'] == '':

                # calculate marker size considering points density
                n_points = len(self.read_instance.station_longitudes[networkspeci][active_map_valid_station_inds])
                
                # calculate figure area and density
                # divide area by 1000 so the function below makes sense
                area = (relevant_axis.bbox.width * relevant_axis.bbox.height) / 1000
                density = n_points / area

                # marker size is calculated using an exponential equation
                # the maximum size is 40 (very low densities)
                # see https://github.com/BSC-ES/providentia/issues/199
                plot_characteristics['plot']['s'] = 1.2**(-density)*40