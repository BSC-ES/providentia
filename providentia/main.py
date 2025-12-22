""" Class for launching Providentia """

from __future__ import print_function
import logging
import sys

from .argument_parser import ProvArgumentParser

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)

class Providentia(object):
    """ Interface class for Providentia. """

    def __init__(self, parser):
        self.parser = parser

    def getargs(self, args):
        """ Return arguments to be passed to the different Providentia modes """

        req = vars(args)
        # print help if no args
        if req.values() == [None for _ in range(len(req.values()))]:
            self.parser.parser.print_help()
            return False
        # pass only valid values and cast boolean strings to boolean
        res = {k: bool(v) if v in ('True', 'False') else v for k, v in req.items() if v}
        return res

    def main(self):
        """ Main functionality of the tool. """
        
        try:
            args = self.parser.parse_args()
            LOG.info(args)
            res = self.getargs(args)
            if res is False:
                return res

            LOG.info(res)

            if args.report:
                from .report import main
            elif args.interpolation:
                from .interpolation.experiment_interpolation_submission import main
            elif args.download:
                from .download import main
            else:
                from .dashboard import main
            main(**res)

        except Exception as err:
            LOG.error('Unhandled exception on Providentia: %s' % err, exc_info=True)
            return False

def main():
    """ Main function. """

    if Providentia(ProvArgumentParser()).main() is False:
        sys.exit(1)
    sys.exit(0)

