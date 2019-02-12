#!/usr/bin/env bash

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

sudo apt-get install python-pip

sudo apt-get install --upgrade spyder
# spyder  windows


sudo add-apt-repository -y ppa:ubuntugis/ubuntugis-unstabl
sudo apt install gdal-bin python-gdal
sudo apt-get install libgdal-dev

#  gdal-config --datadir

export GDAL_DATA=$(gdal-config --datadir)
#export GDAL_DATA=/usr/share/gdal/2.2/

sudo pip install --upgrade pandas
sudo pip install --upgrade geopandas
sudo pip install --upgrade pyyaml

pip install pytest-xdist  # add xdist for distributing tests
#pip install pytest-xvfb  # run xvfb automatically
pip install pytest-cov  # code coverage
pip install coveralls

pip install -q -r requirements.txt

#cd $HOME
#
#git clone https://github.com/MTgeophysics/mtpy.git
#
#cd mtpy
#
#sudo pip install -r requirements.txt


# OK python -c "import geopandas"

#sudo apt install python-pytest

#pytest  tests/core/test_edi.py
#pytest  tests/core/test_ediCollection.py

# pytest tests/core/
# py.test tests/core/
# py.test tests/analysis/
# py.test tests/modeling/
# py.test tests/imaging/
