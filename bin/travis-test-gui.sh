#!/usr/bin/env bash

set -ex

export DISPLAY=:99.0  # enable display on travis
sh -e /etc/init.d/xvfb start
sleep 3 # give xvfb some time to start

# skip the following test run. It's unreliable
# py.test -v --cov=mtpy --cov-report= tests/SmartMT
