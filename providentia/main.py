""" Class for launching Providentia """

from __future__ import print_function
import logging
import sys

from .argument_parser import ProvArgumentParser

logging.basicConfig(level=logging.WARNING)
LOG = logging.getLogger(__name__)

class Providentia(object):
    """
    Class that handles parsing command-line arguments, 
    selecting the appropriate mode.
    """

    def __init__(self, parser):
        """
        Initialize the Providentia launcher with a parser.

        Parameters
        ----------
        parser : ProvArgumentParser
            Argument parser instance used to parse command-line inputs.
        """

        self.parser = parser

    def getargs(self, args):
        """
        Process parsed command-line arguments and return a dictionary of valid values.

        Parameters
        ----------
        args : argparse.Namespace
            Parsed arguments from the command line.

        Returns
        -------
        dict or bool
            Dictionary of valid argument values, or False if no arguments were provided.
        """

        req = vars(args)
        # print help if no args
        if req.values() == [None for _ in range(len(req.values()))]:
            self.parser.parser.print_help()
            return False
        # pass only valid values and cast boolean strings to boolean
        res = {k: bool(v) if v in ('True', 'False') else v for k, v in req.items() if v}
        return res

    def main(self):
        """
        Execute Providentia based on the mode parsed command-line arguments.

        Returns
        -------
        bool
            The return value of the executed mode's main function.
        """
        
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
    """
    Create a Providentia instance with the argument parser 
    and execute the main workflow.
    
    This function is called by the Providentia binary.    
    """

    if Providentia(ProvArgumentParser()).main() is False:
        sys.exit(1)
    sys.exit(0)

