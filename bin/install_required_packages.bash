#! /usr/bin/bash

# run this script after click installed Anaconda2-5.1.0-Windows-x86_64.exe into c:/anaconda2
# and having configured Git/etc/bash.bashrc to use the anaconda paths
conda config --add channels conda-forge
conda install -y --file requirements.txt
# If have problem connecting Internet, you can download the bz2 files from  https://repo.continuum.io/pkgs/
# drop into the folder anaconda2/pkgs, and use conda --offline install option. This function may NOT always work.
# conda install -y --offline --file requirements.txt
conda install -y pyqt=5
conda install -y basemap

