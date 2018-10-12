#! /usr/bin/env bash

# Purpose: Run this script on the command line to see your software environment
# This script can run in Windows git-bash terminal.
#
# How to run:
#
# cd mtpy   # into the repository root dir mtpy
# bin/check_mtpy_env.bash

echo "check git version  ...."
git --version

echo "Check python version ...(should be python 2.7.??) "
python -V

echo "check path to python, conda, and spyder"

which python
which conda
which pip
which spyder

echo

echo "Check conda configuration..."
conda info

echo "Testing some plotting functions. Close the popup figure to continue ....."

echo "Try to run a simple script: examples/scripts/plot_edis.py"
python examples/scripts/plot_edis.py


#echo " Check strating the MTPY GUi...."
#python -OO mtpy/gui/SmartMT/start.py &

#echo "Check staring spyder......................"
#spyder &



