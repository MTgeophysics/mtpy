#! /usr/bin/bash

# get into the mtpy dir
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

echo

echo "Running python mtpy/imaging/plot_response.py "
python mtpy/imaging/plot_response.py


echo "Running python examples/scripts/plot_edis.py" 
python examples/scripts/plot_edis.py

echo " Staring Spyder"

spyder &



