#!/usr/bin/env python

"""
mtpy/mtpy/analysis/distortion.py

Contains functions for the determination of (galvanic) distortion of impedance tensors.
The methods used follow Bibby et al 2005.
As it has been pointed out in that paper, there are various possibilities for constraining the solution, esp. in the 2D case.

Here we just implement the 'most basic' variety for the calculation of the distortion tensor.
Other mehtods can be implemented, but since the optimal assumtions and constraints depend on the application, the actual place for further functions is in an independent, personalised module.



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
import mtpy.analysis.pt as MTpt 
import mtpy.utils.format as MTformat
import mtpy.utils.exceptions as MTexceptions
import mtpy.utils.calculator as MTc

reload(MTexceptions)
reload(MTedi)
reload(MTz)
reload(MTformat)
reload(MTc)
reload(MTpt)


#=================================================================
#Finding the distortion of a Z array. Using the phase tensor (so, Z arrays are transformed into PTs first), following Bibby et al. 2005. #

# First, try to find periods that indicate 1D. From them determine D incl. the g-factor by calculatiing a weighted mean. The g is assumed in order to cater for the missing unknown in the system, it is here set to det(X)^0.5. 
# After that is found, the function no_distortion from the Z module can be called to obtain the unperturbated regional impedance tensor.

#Second, if there are no 1D sections: Find the strike angle, then rotate the Z to the principal axis. In order to do that, use the rotate(-strike) method of the Z module. Then take the real part of the rotated Z. As in the 1D case, we need an assumption to get rid of the (2) unknowns:
# set det(D) = P and det(D) = T, where P,T can be chosen. Common choice is to set  one of P,T to an arbitrary value (e.g. 1). Then check, for which values of the other parameter  S^2 = T^2+4*P*X_12*X_21/det(X) > 0 holds.

def find_distortion():



    dist = np.zeros((2,2))


    return dist

