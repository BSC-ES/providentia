#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

REQUIREMENTS = {
    'install': [
        'jupyterlab',
        'cartopy',
        'matplotlib==3.9.1',
        'netCDF4==1.6.5',
        'numpy==1.26.4',
        'pandas==2.2.2',
        'pyproj==3.6.1',
        'PyQt5==5.15.11',
        'seaborn==0.13.2',
        'scipy==1.12.0',
        'ConfigArgParse==1.7',
        'cftime==1.6.3',
        'KDEpy==1.1.9',
        'xarray==2024.11.0',
        'Sphinx==7.2.6',
        'Sphinx-rtd-theme==2.0.0',
        'coverage==7.4.3',
        'pytest==8.1.1',
        'pypdf==3.4.1',
        'bottleneck==1.3.8',
        'tqdm==4.66.4',
        'remotezip==0.12.3',
        'myst-parser==3.0.1',
        'PyYAML==6.0.1',
        'python-dotenv==1.0.1',
        'paramiko==3.4.0'
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
    version="2.4.0",
    author="Dene Bowdalo, Alba Vilanova Cortezón, Paula Serrano Sierra, Amalia Vradi, Francesco Benincasa",
    author_email="dene.bowdalo@bsc.es, alba.vilanova@bsc.es, paula.serrano@bsc.es, amalia.vradi@bsc.es, francesco.benincasa@bsc.es",
    packages=find_packages(),
    url="https://earth.bsc.es/gitlab/ac/Providentia",
    keywords=["earth sciences", "atmospheric composition",
              "evaluation", "verification", "observations", "air quality"],
    description="Providentia is designed to allow on-the-fly, offline and interactive analysis of experiment outputs, with respect to processed observational data.",
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
    package_data={"": [
        "README.md",
        "LICENSE",
    ]},
    setup_requires=REQUIREMENTS["setup"],
    install_requires=REQUIREMENTS["install"],
    python_requires=">=3.7",
)
