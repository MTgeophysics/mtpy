#!/usr/bin/env python

"""
mtpy/utils/roseplot.py


Class and functions for the general imaging of roseplots.
To be represented via matplotlib.

Class:
    Roseplot - not necessarily restricted to strike angles


Functions:

    Generate roseplot for Strike angle analysis from Z/PT/Edi
    Batch processing of different sites
    Returning plottable set of roseplot information, to be piped into GMT or other mapping tools


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