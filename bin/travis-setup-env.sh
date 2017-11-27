#!/usr/bin/env bash

# for Travis-ci, the following script sets installs all dependencies
set -ex

sudo rm /etc/apt/sources.list.d/mongodb-3.2.list  # temp fix for issue travis-ci/travis-ci#8554
sudo add-apt-repository --yes ppa:ubuntu-sdk-team/ppa  # qt5 source
sudo apt update -qq

# command to install dependencies
sudo apt install libproj-dev  # undocumented gdal dependency
  # install qt4 and dependencies
if [ $QT_VERSION == 4 ]; then sudo apt -y install -qq "libqt4-dev"; fi  # install pyqt dependencies
if [ $QT_VERSION == 5 ]; then sudo apt -y install -qq qtdeclarative5-dev libqt5svg5-dev qtmultimedia5-dev && export QMAKE=/usr/lib/x86_64-linux-gnu/qt5/bin/qmake; fi  # install qt5 and switch to qt5 env
./bin/install_pyqt.sh
  # build gdal from source as the one on ubuntu trusty is out-of-date
#sudo apt -y install -qq libgdal-dev=2.2.2  # dependencies of gdal
./bin/install_libgdal.sh
export CPLUS_INCLUDE_PATH=/usr/include/gdal
export C_INCLUDE_PATH=/usr/include/gdal
export LD_LIBRARY_PATH="/usr/local/gdal/lib:$LD_LIBRARY_PATH"
export PATH="/usr/local/gdal/bin:$PATH"
export GDAL_DATA=$(gdal-config --datadir) && echo GDAL_DATA=$GDAL_DATA
pip install pytest-xdist  # add xdist for distributing tests
#pip install pytest-xvfb  # run xvfb automatically
pip install pytest-cov  # code coverage
pip install coveralls

# manually build and install gdal
#pip install GDAL==$(gdal-config --version | awk -F'[.]' '{print $1"."$2}') --global-option=build_ext --global-option="-I/usr/include/gdal/"
pip install -q matplotlib==$MATPLOTLIB_VERSION
pip install -q -r requirements.txt

# show installed package versions for debugging
pip freeze
#  # install mtpy2 in dev mode
#pip install .

export CACHED=TRUE
