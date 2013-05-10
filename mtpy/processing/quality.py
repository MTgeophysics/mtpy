
"""
mtpy/processing/quality.py

Functions for the analysis of time series data quality.

Output can be visualised with the help of mtpy/imaging/plotquality.py



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import os.path as op

import copy


import  mtpy.utils.exceptions as MTex
#=================================================================
