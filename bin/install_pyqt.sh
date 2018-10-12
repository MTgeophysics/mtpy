#!/usr/bin/env bash

set -ex

BUILDING_DIR=builds

if [ -z "$QT_VERSION" ]; then
    QT_VERSION=5  # install pyqt 5 by default
fi

if [ ! -d "$BUILDING_DIR" ]; then
    mkdir "$BUILDING_DIR"
fi

pushd "$BUILDING_DIR"

# SIP
if [ ! -d "sip-4.19.3" ]; then
    curl -L -o sip-4.19.3.tar.gz https://sourceforge.net/projects/pyqt/files/sip/sip-4.19.3/sip-4.19.3.tar.gz
    echo "4708187f74a4188cb4e294060707106f sip-4.19.3.tar.gz" | md5sum -c -
    tar xzf sip-4.19.3.tar.gz
fi
pushd "sip-4.19.3"
if [ ! -f "BUILT" ] || [ $1 == rebuilt ]; then
    python configure.py
    make -j 4
    touch BUILT
fi
sudo make install
#make install
popd

# PYQT
if [ $QT_VERSION == 4 ]; then
    if [ ! -d PyQt4_gpl_x11-4.12.1 ]; then 
        curl -L -o PyQt4_gpl_x11-4.12.1.tar.gz https://sourceforge.net/projects/pyqt/files/PyQt4/PyQt-4.12.1/PyQt4_gpl_x11-4.12.1.tar.gz
        echo "0112e15858cd7d318a09e7366922f874 PyQt4_gpl_x11-4.12.1.tar.gz" | md5sum -c -
        tar xzf PyQt4_gpl_x11-4.12.1.tar.gz
    fi
    pushd PyQt4_gpl_x11-4.12.1
    if [ ! -f "BUILT" ] || [ $1 == rebuilt ]; then
        python configure.py -c --confirm-license --no-designer-plugin -e QtCore -e QtGui -e QtTest
        make -j 4
        touch BUILT
    fi
    sudo make install
#    make install
    popd
elif [ $QT_VERSION == 5 ]; then
    if [ ! -d "PyQt5_gpl-5.9" ]; then
	    curl -L -o PyQt5_gpl-5.9.tar.gz https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.9/PyQt5_gpl-5.9.tar.gz
	    echo "a409ac0d65ead9178b90c2822759a84b PyQt5_gpl-5.9.tar.gz" | md5sum -c -
	    tar xzf PyQt5_gpl-5.9.tar.gz
    fi
    pushd PyQt5_gpl-5.9
    if [ ! -f "BUILT" ] || [ $1 == rebuilt ]; then
        python configure.py -c --confirm-license --no-designer-plugin -e QtCore -e QtGui -e QtWidgets -e QtTest --qmake=/usr/lib/x86_64-linux-gnu/qt5/bin/qmake
        make -j 4
        touch BUILT
    fi
    sudo make install
    popd
else
    echo "QT version $QT_VERSION is not yet supported"
    false
fi

popd
