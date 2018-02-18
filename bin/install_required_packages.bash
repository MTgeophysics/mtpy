#! /usr/bin/bash

# run this script after click installed Anaconda2-5.1.0-Windows-x86_64.exe into c:/anaconda2
# and having configured Git/etc/bash.bashrc to use the anaconda paths
conda config --add channels conda-forge
conda install --file mtpy/requirements.txt
conda install -y pyqt=5
