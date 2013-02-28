#!/usr/bin/env python

"""
mtpy/utils/inductionarrows.py


Class and functions for the imaging of induction arrows.
To be embedded in local area rectangular plots, GMT maps, or Python basemap

Class:
    Inductionarrows (from Z object/array or PT object/array or Edi object/file)


Functions:

    Batch processing of different sites
    Returning plottable vector information, to be piped into GMT or other mapping tools


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