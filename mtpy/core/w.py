#!/usr/bin/env python

"""
mtpy/mtpy/core/w.py

Contains classes and functions for handling GDS transfer function tensors (W). 
 
    Class:
    "W" contains information about a GDS W-tensor. 

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
reload (MTedi)

import mtpy.utils.format as MTformat
reload(MTformat)

import mtpy.utils.exceptions as MTexceptions
reload(MTexceptions)


#=================================================================


#------------------------
class W(object):
    """
        W class - generates a GDS W-tensor object.


    """

    def __init__(self, *a, **b):

        pass
    
