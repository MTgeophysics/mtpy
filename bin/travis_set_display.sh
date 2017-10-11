#!/usr/bin/env bash

set -ex
# setting up display on travis
export DISPLAY=:99.0
sh -e /etc/init.d/xvfb start
sleep 3 # give xvfb some time to start
