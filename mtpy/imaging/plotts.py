#!/usr/bin/env python

"""
mtpy/utils/plotts.py


Class and functions for the imaging of a (set of) time series


Class:
    TimeSeries()

Functions:
    Batch processing 
    plot setup
    

Output via matplotlib



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np

import mtpy.core.z as MTz
import mtpy.core.edi as MTedi
import mtpy.analysis.pt as MTpt

from mtpy.utils.exceptions import *


#=================================================================