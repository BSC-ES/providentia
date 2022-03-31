def get_characteristics_per_plot_type(self):

    characteristics_per_plot_type = {
        'timeseries': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 3},
                       'xtick_share': True, 'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                       'page_title': {'t': 'Time Series', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                       'axis_title': {'label': '', 'fontsize': 8}, 'axis_xlabel': {'xlabel': 'Time', 'fontsize': 8},
                       'axis_ylabel': {'ylabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                       'xticks': {'labelsize': 7, 'rotation': 0}, 'yticks': {'labelsize': 7},
                       'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                       'subplots_adjust': {'top': 0.90, 'bottom': 0.08}},
            
        'distribution': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                         'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                         'page_title': {'t': 'Distribution', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                         'axis_title': {'label': '', 'fontsize': 8},
                         'axis_xlabel': {'xlabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                         'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                         'yticks': {'labelsize': 7}, 'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0},
                         'tightlayout': True, 'subplots_adjust': {'top': 0.90}},

        'distribution_bias': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                              'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                              'page_title': {'t': 'Distributional bias', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98}, 
                              'axis_title': {'label': '', 'fontsize': 8}, 
                              'axis_xlabel': {'xlabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                              'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                              'yticks': {'labelsize': 7},
                              'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 'tightlayout': True,
                              'subplots_adjust': {'top': 0.90}},

        'heatmap': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 1},
                    'page_title': {'t': 'Statistical Heatmap', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                    'axis_title': {'label': '', 'fontsize': 8}, 'xticks': {'labelsize': 7, 'rotation': -270},
                    'yticks': {'labelsize': 7, 'rotation': -315}, 'tightlayout': True,
                    'subplots_adjust': {'top': 0.90}, 'cb_xlabel': {'xlabel': '', 'fontsize': 8},
                    'cb_xticks': {'labelsize': 8}, 'annot': True},

        'periodic-violin': {'summary_pages':[], 'station_pages':[],
                            'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                            'page_title': {'t': 'Violin plots', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                            'axis_title': {'label': '', 'fontsize': 8},
                            'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                            'yticks': {'labelsize': 7},
                            'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                            'subplots_adjust': {'top': 0.85}},    

        'scatter': {'summary_pages':[], 'station_pages':[],'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                    'page_title': {'t': 'Scatter plots', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                    'axis_title': {'label': '', 'fontsize': 8}, 'xticks': {'labelsize': 7}, 'yticks': {'labelsize': 7},
                    'axis_xlabel': {'xlabel': 'Observations', 'fontsize': 8}, 'axis_ylabel': {'ylabel': 'Experiment', 'fontsize': 8},
                    'markers': {'character': 'o', 'size': 2},
                    'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0}}
    }

    return characteristics_per_plot_type

def get_characteristics_per_plot_type_templates(self):

    characteristics_per_plot_type_templates = {
        'map-basicstat-obs': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                              'subplot_kw': {'projection': self.plotcrs}},
                              'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                              'axis_title': {'label': '', 'fontsize': 8},
                              'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                              'cb': {'position': [0.5, 0.95, 0.4, 0.04]},
                              'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}},

        'map-basicstat': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                          'subplot_kw': {'projection': self.plotcrs}},
                          'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                          'axis_title': {'label': '', 'fontsize': 8},
                          'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                          'cb': {'position': [0.5, 0.95, 0.4, 0.04]},
                          'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}},

        'map-biasstat': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                         'subplot_kw': {'projection': self.plotcrs}},
                         'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                         'axis_title': {'label': '', 'fontsize': 8},
                         'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                         'cb': {'position': [0.5, 0.95, 0.4, 0.04]}, 'cb_xlabel': {'xlabel': '', 'fontsize': 8},
                         'cb_xticks': {'labelsize': 8}},

        'periodic-basicstat': {'summary_pages':[], 'station_pages':[],
                               'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                               'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                               'axis_title': {'label': '', 'fontsize': 8},
                               'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                               'yticks': {'labelsize': 7},
                               'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                               'subplots_adjust': {'top': 0.85}},

        'periodic-biasstat': {'summary_pages':[], 'station_pages':[],
                              'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                              'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                              'axis_title': {'label': '', 'fontsize': 8},
                              'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                              'yticks': {'labelsize': 7},
                              'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                              'subplots_adjust': {'top': 0.85}}
    }
    return characteristics_per_plot_type_templates
