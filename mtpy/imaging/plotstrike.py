#!/usr/bin/env python

"""
mtpy/utils/plotstrike.py


Functions for the imaging of strike analysis of one or more Z/PT/Edi element.

Use roseplot.py module

Output via matplotlib or perhaps embedded in a map



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np

import mtpy.core.z as MTz
import mtpy.core.edi as MTedi
import mtpy.analysis.pt as MTpt
import mtpy.imaging.roseplot as MTrose


from mtpy.utils.exceptions import *


#=================================================================