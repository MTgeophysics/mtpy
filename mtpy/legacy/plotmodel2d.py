#!/usr/bin/env python

"""
mtpy/utils/plotmodel2d.py


Class and functions for the imaging of a 2D model from WingLink, ModEM, Occam


Class:
    Model2d()

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