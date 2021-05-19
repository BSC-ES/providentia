import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# portrait/landscape page figsize
portrait_figsize = (8.27, 11.69)
landscape_figsize = (11.69, 8.27)
dpi = 200


class ProvidentiaOffline:
    """Run Providentia offline reports"""

    def __init__(self, initialised, read_type='parallel'):
        self.initialised = initialised
        self.read_type = read_type
        print("you are here!")

        filename = "test_pdf.pdf"

        # open new PDF file
        with PdfPages(filename) as pdf:
            self.pdf = pdf
            self.make_header()

        exit()

    def make_header(self):
        # set tile
        page = plt.figure(figsize=portrait_figsize)
        if hasattr(self.initialised, 'report_title'):
            txt = self.initialised.report_tile
        else:
            txt = 'An example report'
        # TODO: define default title in case no title in conf
        page.text(0.5, 0.9, txt, transform=page.transFigure,
                  size=20, ha="center", va='top', wrap=True)

        experiment_labels = [exp.strip() for exp in self.initialised.experiments.split(",")]
        txt = 'Network = {}\nTemporal Resolution = {}\n' \
              'Species = {}\nDate Range = {} - {}\nExperiments = {}\n'\
            .format(self.initialised.selected_network, self.initialised.selected_resolution,
                    self.initialised.selected_species, self.initialised.start_date,
                    self.initialised.end_date, experiment_labels)
        # txt += 'Station Subset/s = {}\n'.format(self.station_subset_names)
        if hasattr(self.initialised, 'bounding_box'):
            txt += 'Bounding Box = [{}E:{}E, {}N:{}N]\n'.format(self.initialised.bounding_box['longitude']['min'],
                                                                self.initialised.bounding_box['longitude']['max'],
                                                                self.initialised.bounding_box['latitude']['min'],
                                                                self.initialised.bounding_box['latitude']['max'])
        if hasattr(self.initialised, 'data_availability_filter'):
            txt += 'Data Availability = {}'.format(self.initialised.data_availability_filter)
        page.text(0.5, 0.82, txt, transform=page.transFigure, weight='light', size=15, ha="center", va='top', wrap=True)
        self.pdf.savefig(page, dpi=dpi)
        plt.close(page)
