#!/usr/bin/env bash

# in Ubuntu 16.04 desktop
# IF the system default python 2.7 is used.

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

cd $HOME

git clone https://github.com/MTgeophysics/mtpy.git

cd mtpy

sudo pip install -r requirements.txt


# OK python -c "import geopandas"

sudo apt install python-pytest

pytest  tests/core/test_edi.py
pytest  tests/core/test_ediCollection.py

# pytest tests/core/
# py.test tests/core/
# py.test tests/analysis/
# py.test tests/modeling/
# py.test tests/imaging/
