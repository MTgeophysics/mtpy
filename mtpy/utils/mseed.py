#!/usr/bin/env python

"""
/mtpy/utils/mseed.py

This modules contains functions for the conversion of raw time series: ASCII to/from miniSeed. 

The functionality is based on obspy.mseed/pyrocko


@UofA, 2013
(LK)

"""

#=================================================================

import obspy.mseed as omseed
#import pyrocko as pmseed

from mtpy.utils.exceptions import *

#=================================================================