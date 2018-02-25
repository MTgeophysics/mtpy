#! /usr/bin/bash

# Purpose: Run this script on the command line to see your software environment
# This script can run in both Linux and Windows git-bash terminal.

# get into the mtpy repository dir
# cd mtpy

echo "check path to python, conda, and spyder"

which python
which conda
which pip
which spyder


echo "check git version  ...."
git --version

echo "Check python version ... "
python -V

echo

echo "Check conda configuration..."

conda info

echo "Testing some plotting functions. Close each figure to continue"

echo "Running python mtpy/imaging/plot_response.py "
python mtpy/imaging/plot_response.py

echo "Running python examples/scripts/plot_edis.py" 
python examples/scripts/plot_edis.py

echo " Strating the MTPY GUi...."
python -OO mtpy/gui/SmartMT/start.py &

#echo "Staring Spyder......................"
#spyder &



