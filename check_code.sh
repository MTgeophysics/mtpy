#!/usr/bin/env bash
# script for running code checks.

#set -eu
#set -x

#pep8 tests mtpy tests examples --max-line-length 120
pep8 tests  exp  --max-line-length 120

echo "pylint .... "
pylint -j 2 --reports=no tests 

echo "pylint ....end "

# Run tests, taking coverage.
# Users can specify extra folders as arguments.
#? py.test -r sx --cov datacube --durations=5 datacube tests datacube_apps $@

