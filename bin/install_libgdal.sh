#!/usr/bin/env bash
# building gdal from source

set -ex

BUILDING_DIR=builds
GDAL_VERSION=2.2.2

if [ ! -d "$BUILDING_DIR" ]; then
    mkdir "$BUILDING_DIR"
fi

pushd "$BUILDING_DIR"
if [ ! -d "gdal-$GDAL_VERSION" ]; then
    wget "http://download.osgeo.org/gdal/$GDAL_VERSION/gdal-$GDAL_VERSION.tar.gz" -O "gdal-$GDAL_VERSION.tar.gz"
    tar xfz "gdal-$GDAL_VERSION.tar.gz"
fi
pushd "gdal-$GDAL_VERSION"
./configure --with-python --prefix=/usr/local/gdal
make
sudo make install
sudo ldconfig
popd
popd
