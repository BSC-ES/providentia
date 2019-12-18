#!/usr/bin/env python

"""
Main module for Providentia. Only contains an interface class to all
functionality implemented on Providentia.
"""

# Copyright 2016 Earth Sciences Department, BSC-CNS

# This file is part of Providentia.

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
from .config import ArgumentParser
from . import dashboard

import sys
import logging
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger(__name__)


class Providentia(object):
    """
    Interface class for Providentia
    """
    def __init__(self, parser):
        self.parser = parser

    def main(self):
        """
        Main functionality of the tool
        """
        try:
            args = self.parser.parse_args()
            log.info(args)
            req = vars(args)
            # print help if no args
            if req.values() == [None for i in range(len(req.values()))]:
                self.parser.parser.print_help()
                return False
            # pass only valid values and cast boolean strings to boolean
            res = {k: eval(v) if v in ('True', 'False') else v for k, v in req.items() if v}
            log.info(res)

            # call dashboard
            dashboard.main(**res)

        except Exception as err:
            log.error('Unhandled exception on Providentia: %s' % err, exc_info=True)
            return False


def main():
    """ Main function """
    if Providentia(ArgumentParser()).main() is False:
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    if Providentia(ArgumentParser()).main() is False:
        sys.exit(1)
    sys.exit(0)
