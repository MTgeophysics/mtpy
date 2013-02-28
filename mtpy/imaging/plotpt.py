
"""
mtpy/utils/plotpt.py


Class and functions for the imaging of a (set of) Phase Tensor(s)


Class:
    PlotPT() -- generated with a MTpt.PhaseTensor instance

    calculates ellipses (incl. error-ellipses), to be used in pseudosections or in maps


Functions:
    Batch processing 
    plot setup
    convert to vector information (e.g. GMTs psxy style)
    

Output as ellipses directly via matplotlib, or ellipses vector information to be piped to GMT or basemap



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