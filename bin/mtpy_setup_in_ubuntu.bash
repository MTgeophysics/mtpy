#!/usr/bin/env bash

# for Travis-ci, the following script sets installs all dependencies
set -ex

# Here is a user guide for Ubuntu 16.04 desktop users.
# Assume the system default python 2.7 is used.
# We plan to migrate mtpy to python3 later. But Not at this point of time.
# check the version by command
# python -V

# Note that it's easy to be confused if your system has multiple version of python.
# we recommend to use Ananconda python which is independent of the operating system:
# See wiki page:
# https://github.com/MTgeophysics/mtpy/wiki/MTPy-installation-guide-for-Linux-system

# Keep this script for reference.

sudo apt -y install python-pip


sudo add-apt-repository -y ppa:ubuntugis/ubuntugis-unstabl
sudo apt -y install gdal-bin python-gdal
sudo apt -y install libgdal-dev

# gdal-config --datadir

export GDAL_DATA=$(gdal-config --datadir) && echo GDAL_DATA=$GDAL_DATA


# do not use sudo for pip install !
 pip install --upgrade pandas
 pip install --upgrade geopandas
 pip install --upgrade pyyaml

 pip install pytest-xdist  # add xdist for distributing tests
#pip install pytest-xvfb  # run xvfb automatically
 pip install pytest-cov  # code coverage
 pip install coveralls

 pip install -q -r requirements.txt

echo "************************** show installed python package versions *********************** "
pip freeze


# OK python -c "import geopandas"

# apt install python-pytest

#pytest  tests/core/test_edi.py
#pytest  tests/core/test_ediCollection.py

# pytest tests/core/
# py.test tests/core/
# py.test tests/analysis/
# py.test tests/modeling/
# py.test tests/imaging/
