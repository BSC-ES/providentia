#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2018 Barcelona Supercomputing Center - Centro Nacional de
# Supercomputación (BSC-CNS)

# This file is part of Providentia 

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

from os import path
from setuptools import setup

from providentia import __version__

REQUIREMENTS = {
    'test': [
        'pytest>=3.9',
        'pytest-cov',
        'pytest-env',
        'pytest-flake8',
        'pytest-html',
        'pytest-metadata>=1.5.1',
    ],
    'setup': [
#        'pytest-runner',
        'setuptools_scm',
    ],
}

setup(
    # Application name:
    name="providentia",
    license='GNU GPL v3',
    # Version number (initial):
    version=__version__,

    # Application author details:
    author="Dene Bowdalo, Amalia Vradi, Alba Vilanova Cortezón, Francesco Benincasa",
    author_email="dene.bowdalo@bsc.es, amalia.vradi@bsc.es, alba.vilanova@bsc.es, francesco.benincasa@bsc.es",

    # Packages
    packages=['providentia', 'providentia.interpolation'],

    # Include additional files into the package
    include_package_data=True,
    scripts=['bin/providentia'],

    # Details
    url="https://earth.bsc.es/gitlab/ac/Providentia",

    keywords=['earth sciences', 'evaluation', 'verification', 'observations', 'NWP models',
              'air quality'],
    description="Providentia is designed to allow on-the-fly and offline analysis of experiment outputs, with respect to processed observational data.",
    #    long_description=open("README.rst").read(),
    #    long_description_content_type='text/x-rst',

    # Dependent packages (distributions)
    setup_requires=REQUIREMENTS['setup'],
    install_requires=[
        "matplotlib",
        "pandas",
        "Cartopy",
        "numpy",
        "netCDF4",
        "seaborn",
        "PyQt5",
        "scipy",
    ],
    #tests_require=REQUIREMENTS['test'],

)
