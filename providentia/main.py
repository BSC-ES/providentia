""" Main Providentia module """
# Providentia is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# Providentia is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with Providentia. If not, see <http://www.gnu.org/licenses/>.

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
        """ Return arguments to be passed to the dashboard. """

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

            if args.offline:
                from . import offline as offline
                offline.main(**res)
            elif args.interpolation:
                print(res)
            else:
                from . import dashboard as dashboard
                dashboard.main(**res)

        except Exception as err:
            LOG.error('Unhandled exception on Providentia: %s' % err, exc_info=True)
            return False

def main():
    """ Main function. """

    if Providentia(ProvArgumentParser()).main() is False:
        sys.exit(1)
    sys.exit(0)

