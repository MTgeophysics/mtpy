#!/usr/bin/env bash

set -ex

BUILDING_DIR=builds

if [ -z "$QT_VERSION" ]; then
    QT_VERSION=4  # install pyqt 4 by default
fi

if [ ! -d "$BUILDING_DIR" ]; then
    mkdir "$BUILDING_DIR"
fi

pushd "$BUILDING_DIR"

# SIP
if [ ! -d "sip-4.19.3" ];
    then curl -L -o sip-4.19.3.tar.gz https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.3/sip-4.19.3.tar.gz;
    tar xzf sip-4.19.3.tar.gz
fi
pushd "sip-4.19.3"
python configure.py
make
sudo make install
#make install
popd

# PYQT
if [ $QT_VERSION == 4 ]; then
    if [ ! -d PyQt4_gpl_x11-4.12.1 ];
        then curl -L -o PyQt4_gpl_x11-4.12.1.tar.gz https://sourceforge.net/projects/pyqt/files/PyQt4/PyQt-4.12.1/PyQt4_gpl_x11-4.12.1.tar.gz;
        tar xzf PyQt4_gpl_x11-4.12.1.tar.gz
    fi
    pushd PyQt4_gpl_x11-4.12.1
    python configure.py -c --confirm-license --no-designer-plugin -e QtCore -e QtGui -e QtTest
    make
    sudo make install
#    make install
    popd
    else
        echo "QT version $QT_VERSION is not yet supported"
        false
fi

popd
