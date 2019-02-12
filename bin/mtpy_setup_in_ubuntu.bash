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

sudo -H apt-get install python-pip

# spyder  windows


sudo -H add-apt-repository -y ppa:ubuntugis/ubuntugis-unstabl
sudo -H apt install gdal-bin python-gdal
sudo -H apt-get install libgdal-dev

#  gdal-config --datadir

export GDAL_DATA=$(gdal-config --datadir)
#export GDAL_DATA=/usr/share/gdal/2.2/

sudo -H pip install --upgrade pandas
sudo -H pip install --upgrade geopandas
sudo -H pip install --upgrade pyyaml

sudo -H pip install pytest-xdist  # add xdist for distributing tests
#pip install pytest-xvfb  # run xvfb automatically
sudo -H pip install pytest-cov  # code coverage
sudo -H pip install coveralls

sudo -H pip install -q -r requirements.txt


echo "************************** show installed python package versions *********************** "
pip freeze


export CACHED=TRUE


cd $HOME

git clone https://github.com/MTgeophysics/mtpy.git

cd mtpy

sudo -H pip install -r requirements.txt


# OK python -c "import geopandas"

#sudo -H apt install python-pytest

#pytest  tests/core/test_edi.py
#pytest  tests/core/test_ediCollection.py

# pytest tests/core/
# py.test tests/core/
# py.test tests/analysis/
# py.test tests/modeling/
# py.test tests/imaging/
