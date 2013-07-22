#!/usr/bin/env python

"""
mtpy/mtpy/analysis/niblettbostick.py

Contains functions for the calculation of the Niblett-Bostick transformation of
impedance tensors.

The methods follow
- Niblett 
- Bostick  
- Jones
- J. RODRIGUEZ, F.J. ESPARZA, E. GOMEZ-TREVINO

Niblett-Bostick transformations are possible in 1D and 2D.


    Functions:


@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np

import mtpy.core.z as MTz 
import mtpy.analysis.geometry as MTge 
import mtpy.utils.exceptions as MTex
import mtpy.utils.calculator as MTcc




def rhophi2rhodepth(rho, phase, period):

	"""
	Convert a period-dependent pair of rho/phase (Ohm meters/rad) 
	into rho/depth (Ohm meters/meters)

	The conversion uses the simplified transformation without derivatives.

	Input:
	- apparent resistivity (Ohm meters
	- phase angle (rad)
	- period (seconds)

	Output:
	- resistivity estimate (Ohm meters)
	- depth (meters) 

	"""

	depth = np.sqrt(rho*period/2/np.pi/MTcc.mu0)

	rho_nb = rho * (np.pi/2/phase - 1)




	return rho_nb, depth



