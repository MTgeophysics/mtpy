"""
mtpy/utils/metadata.py

Helper functions for the handling of meta data. 

E.g. converting field book xls-sheet information into a config-style configuration file, or check entries for consistency.



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import ConfigParser

import math, cmath

import mtpy.utils.exceptions as MTex


#=================================================================