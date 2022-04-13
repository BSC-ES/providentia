def get_characteristics_per_plot_type(self):

    characteristics_per_plot_type = {
        'timeseries': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 3},
                       'xtick_share': True, 'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                       'page_title': {'t': 'Timeseries', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                       'axis_title': {'label': '', 'fontsize': 8}, 'axis_xlabel': {'xlabel': 'Time', 'fontsize': 8},
                       'axis_ylabel': {'ylabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                       'xticks': {'labelsize': 7, 'rotation': 0}, 'yticks': {'labelsize': 7},
                       'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0}, 
                       'tightlayout': True, 'subplots_adjust': {'top': 0.90, 'bottom': 0.08}, 
                       'bias': {'title': 'Timeseries bias'},
                       'annotate_stats': ['Mean', 'Min'],
                       'annotate_text': {'s': '', 'loc': 'upper right', 'fontsize': 7},
                       'annotate_bbox': {'facecolor': 'white', 'edgecolor': 'gainsboro', 'alpha': 1}},
        
        'distribution': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                         'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                         'page_title': {'t': 'Distribution', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                         'axis_title': {'label': '', 'fontsize': 8},
                         'axis_xlabel': {'xlabel': '{}'.format(self.datareader.measurement_units), 'fontsize': 8},
                         'axis_ylabel': {'ylabel': 'Density', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                         'yticks': {'labelsize': 7}, 'legend': {'loc': 'upper right', 'ncol': 3, 'fontsize': 8.0},
                         'tightlayout': True, 'subplots_adjust': {'top': 0.90},
                         'bias': {'title': 'Distributional bias'},
                         'annotate_stats': ['Max', 'Min'],
                         'annotate_text': {'s': '', 'loc': 'upper right', 'fontsize': 7},
                         'annotate_bbox': {'facecolor': 'white', 'edgecolor': 'gainsboro', 'alpha': 1}},

        'heatmap': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 1},
                    'page_title': {'t': 'Statistical Heatmap', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                    'axis_title': {'label': '', 'fontsize': 8}, 
                    'xticks': {'labelsize': 7, 'rotation': -270}, 'yticks': {'labelsize': 7, 'rotation': -315}, 
                    'tightlayout': True, 'subplots_adjust': {'top': 0.90}, 
                    'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}, 'annot': True},

        'periodic-violin': {'summary_pages':[], 'station_pages':[],
                            'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                            'page_title': {'t': 'Violin plots', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                            'axis_title': {'label': '', 'fontsize': 8},
                            'xticks': {'labelsize': 7}, 'yticks': {'labelsize': 7}, 
                            'axis_ylabel': {'ylabel': '', 'fontsize': 8},
                            'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                            'subplots_adjust': {'top': 0.80, 'hspace': 0.30},
                            'bias': {'title': 'Periodic violin bias'},
                            'annotate_stats': ['Max', 'Min'],
                            'annotate_text': {'s': '', 'loc': 'upper right', 'fontsize': 7},
                            'annotate_bbox': {'facecolor': 'white', 'edgecolor': 'gainsboro', 'alpha': 1}},    

        'scatter': {'summary_pages':[], 'station_pages':[],'figure': {'figsize': self.portrait_figsize, 'ncols': 2, 'nrows': 4},
                    'grid': {'axis': 'both', 'color': 'lightgrey', 'alpha': 0.8},
                    'page_title': {'t': 'Scatter plots', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                    'axis_title': {'label': '', 'fontsize': 8}, 
                    'xticks': {'labelsize': 7}, 'yticks': {'labelsize': 7},
                    'axis_xlabel': {'xlabel': 'Observations', 'fontsize': 8}, 'axis_ylabel': {'ylabel': 'Experiment', 'fontsize': 8},
                    'markers': {'character': 'o', 'size': 2},
                    'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0}, 
                    'subplots_adjust': {'top': 0.90, 'hspace': 0.40},
                    'bias': {'title': 'Scatter plot bias'},
                    'annotate_stats': ['r2', 'RMSE'],
                    'annotate_text': {'s': '', 'loc': 'upper left', 'fontsize': 7},
                    'annotate_bbox': {'facecolor': 'white', 'edgecolor': 'gainsboro', 'alpha': 1}}
    }

    return characteristics_per_plot_type

def get_characteristics_per_plot_type_templates(self):

    characteristics_per_plot_type_templates = {
        
        'periodic': {'summary_pages':[], 'station_pages':[],
                     'figure': {'figsize': self.landscape_figsize, 'ncols': 2, 'nrows': 2},
                     'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                     'axis_title': {'label': '', 'fontsize': 8},
                     'axis_ylabel': {'ylabel': '', 'fontsize': 8}, 'xticks': {'labelsize': 7},
                     'yticks': {'labelsize': 7},
                     'legend': {'loc': 'upper right', 'ncol': 4, 'fontsize': 8.0},
                     'subplots_adjust': {'top': 0.85},
                     'annotate_stats': ['Max', 'Min'],
                     'annotate_text': {'s': '', 'loc': 'upper right', 'fontsize': 7},
                     'annotate_bbox': {'facecolor': 'white', 'edgecolor': 'gainsboro', 'alpha': 1}},

        'map': {'summary_pages':[], 'station_pages':[], 'figure': {'figsize': self.landscape_figsize, 'ncols': 4, 'nrows': 4,
                'subplot_kw': {'projection': self.plotcrs}},
                'page_title': {'t': '', 'fontsize': 15, 'ha': 'left', 'x': 0.05, 'y': 0.98},
                'axis_title': {'label': '', 'fontsize': 8},
                'subplots_adjust': {'top': 0.85, 'hspace': 0.28, 'wspace': 0.28},
                'cb': {'position': [0.5, 0.95, 0.4, 0.04]},
                'cb_xlabel': {'xlabel': '', 'fontsize': 8}, 'cb_xticks': {'labelsize': 8}},
    }
    
    return characteristics_per_plot_type_templates
