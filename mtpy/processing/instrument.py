
"""
mtpy/processing/instrument.py

Functions for the removal of instrument response from time series data.

Either working on ASCII data or on miniSeed



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import os.path as op

import copy

import mtpy.utils.mseed as MTmseed
import obspy.mseed as Omseed
import mtpy.utils.exceptions as MTexceptions

#=================================================================
