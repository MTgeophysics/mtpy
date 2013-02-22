#!/usr/bin/env python

"""
mtpy/mtpy/analysis/geometry.py

Contains classes and functions for handling geometry analysis of impedance tensors:

dimensionality, strike directions, alphas/betas/...
 


    Class:

        Methods:



    Functions:


@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np
import os
import sys
import os.path as op
import math, cmath
import time, calendar 

import mtpy.core.edi as MTedi 
import mtpy.core.z as MTz 
import mtpy.utils.format as MTformat
import mtpy.utils.exceptions as MTexceptions
import mtpy.utils.calculator as MTc

reload(MTexceptions)
reload(MTedi)
reload(MTz)
reload(MTformat)
reload(MTc)


#=================================================================