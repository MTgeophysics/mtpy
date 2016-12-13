#!/usr/bin/env bash
# script for running code checks.

#set -eu
#set -x

#pep8 tests mtpy tests examples --max-line-length 120
pep8 tests  exp  --max-line-length 120

echo "pylint .... "
#pylint -j 2 --reports=no tests  #j2 cpu run no report
pylint -E tests
echo "pylint ....end "

py.test tests
# Run tests, taking coverage.
# Users can specify extra folders as arguments.
#? py.test -r sx --cov datacube --durations=5 datacube tests datacube_apps $@

