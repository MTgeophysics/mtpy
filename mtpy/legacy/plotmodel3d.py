#!/usr/bin/env python

"""
mtpy/utils/plotmodel3d.py


Class and functions for the imaging of a 3D model from WingLink, ModEM, Ws3Dinv


Class:
    Model3d()

Functions:
    Batch processing 
    plot setup
    

Output to be sent to mayavi or paraview



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