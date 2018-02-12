#!/usr/bin/env bash
# script for running code checks.

#set -eu
#set -x

#pep8 tests mtpy  --max-line-length 120
pep8 mtpy --max-line-length 120

echo "pylint .... "
#pylint -j 2 --reports=no tests  #j2 cpu run no report
pylint -E mtpy
echo "pylint ....end "

# do pytest and generate coverage report
py.test tests
# Run tests, taking coverage.
# Users can specify extra folders as arguments.
#? py.test -r sx --cov datacube --durations=5 datacube tests datacube_apps $@

# send/publish coverage stats to coveralls.io
# and get an URL like https://coveralls.io/jobs/30383099
# coveralls
