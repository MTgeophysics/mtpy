#!/usr/bin/env python

"""
mtpy/utils/plotcoherence.py


Class and functions for the imaging of coherences of one or more Z/PT/Edi element(s), as well as time series data coherence.


Class:
    Coherence()

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

import mtpy.processing.coherence as MTcoh
from mtpy.utils.exceptions import *


#=================================================================