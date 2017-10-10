#!/usr/bin/env bash
# building gdal from source

BUILDING_DIR=cache
GDAL_VERSION=2.2.2

set -ex

if [ ! -d "$BUILDING_DIR" ]; then
    mkdir "$BUILDING_DIR"
fi
pushd "$BUILDING_DIR"
if [ ! -d "gdal-$GDAL_VERSION" ]; then
    wget "http://download.osgeo.org/gdal/CURRENT/gdal-$GDAL_VERSION.tar.gz" -O "gdal-$GDAL_VERSION.tar.gz"
    tar xfz "gdal-$GDAL_VERSION.tar.gz"
fi
pushd "gdal-$GDAL_VERSION"
./configure --with-python --prefix=/usr/local/gdal
make
sudo make install
popd
popd
ldconfig
