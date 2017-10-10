#!/usr/bin/env bash
# building gdal from source

BUILDING_DIR=cache
GDAL_VERSION=2.2.2

set -ex

if [ ! -d "$BUILDING_DIR" ]; then
    mkdir "$BUILDING_DIR"
fi
if [ ! -f "$BUILDING_DIR/gdal-$GDAL_VERSION.tar.gz" ]; then
    wget "http://download.osgeo.org/gdal/CURRENT/gdal-$GDAL_VERSION.tar.gz" -O "$BUILDING_DIR/gdal-$GDAL_VERSION.tar.gz"
    tar xvfz "gdal-$GDAL_VERSION.tar.gz"
fi
pushd "$BUILDING_DIR/gdal-$GDAL_VERSION"
./configure --with-python
make
sudo make install
popd
