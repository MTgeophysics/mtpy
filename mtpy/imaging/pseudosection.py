#!/usr/bin/env python

"""
mtpy/utils/pseudosections.py


Class and functions for the imaging of MT pseudosections.
To be represented via matplotlib.

Class:
    Pseudosection - not necessarily restricted to PT ellipses


Functions:

    Generate pseudosection for PhaseTensorEllipses
    Batch processing of different sites
    Returning plottable set of pseudosection information, to be piped into GMT or other mapping tools for


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