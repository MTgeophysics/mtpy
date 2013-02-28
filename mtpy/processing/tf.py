
"""
mtpy/processing/tf.py

Functions for the time-frequency analysis of time series data.

Output can be visualised with the help of mtpy/imaging/spectrogram.py



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import os.path as op

import copy


import  mtpy.utils.exceptions as MTexceptions

#=================================================================
