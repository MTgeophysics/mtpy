#!/usr/bin/env bash

export DISPLAY=:99.0  # enable display on travis
sh -e /etc/init.d/xvfb start
sleep 3 # give xvfb some time to start
py.test -v --cov=mtpy --cov-report= tests/SmartMT
