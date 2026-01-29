#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

REQUIREMENTS = {
    'install': [
        'jupyterlab',
        'cartopy',
        'ghostscript',
        'matplotlib==3.10.7',
        'netCDF4==1.7.2',
        'numpy==2.3.5',
        'pandas==2.3.3',
        'pyproj==3.7.2',
        'PyQt5==5.15.11',
        'seaborn==0.13.2',
        'scipy==1.16.3',
        'ConfigArgParse==1.7.1',
        'cftime==1.6.5',
        'xarray==2025.11.0',
        'Sphinx==8.2.3',
        'Sphinx-rtd-theme==3.0.2',
        'coverage==7.12.0',
        'coverage-badge==1.1.2',
        'pytest==9.0.1',
        'pypdf==6.3.0',
        'bottleneck==1.6.0',
        'tqdm==4.67.1',
        'remotezip==0.12.3',
        'myst-parser==4.0.1',
        'PyYAML==6.0.1',
        'python-dotenv==1.2.1',
        'paramiko==4.0.0',
        'memray==1.19.1',
        'cdsapi==0.7.7',
        'pycountry==24.6.1'
    ],
    'setup': [
        'setuptools_scm',
    ],
}

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="providentia",
    license="GNU GPL v3",
    version="3.0.1",
    author="Dene Bowdalo, Alba Vilanova Cortezón, Paula Serrano Sierra, Francesco Benincasa, Amalia Vradi",
    author_email="dene.bowdalo@bsc.es, alba.vilanova@bsc.es, paula.serrano@bsc.es, francesco.benincasa@bsc.es, amalia.vradi@bsc.es",
    packages=find_packages(),
    url="https://github.com/BSC-ES/providentia",
    keywords=["earth sciences", "atmospheric composition",
              "evaluation", "verification", "observations", "air quality"],
    description="Providentia is designed to allow on-the-fly, offline and interactive analysis of model outputs, with respect to processed observational data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "Programming Language :: Python :: 3.11",
        "Operating System :: POSIX :: Linux",
        "Operating System :: MacOS",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Intended Audience :: Science/Research",
        "Natural Language :: English"
    ],
    include_package_data=True,
    package_data={"": [
        "README.md",
        "LICENSE"
    ]},
    setup_requires=REQUIREMENTS["setup"],
    install_requires=REQUIREMENTS["install"],
    python_requires=">=3.11",
)