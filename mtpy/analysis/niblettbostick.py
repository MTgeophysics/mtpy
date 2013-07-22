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
	- phase angle (degrees)
	- period (seconds)

	Output:
	- resistivity estimate (Ohm meters)
	- depth (meters) 

	"""

	depth = np.sqrt(rho*period/2/np.pi/MTcc.mu0)

	# phase angle needed in rad
	rho_nb = rho * (np.pi/2/np.deg2rad(phase) - 1)

	return rho_nb, depth


def calculate_znb(z_object = None, z_array = None, periods = None):
	"""
	Determine an array of Z_nb (depth dependent Niblett-Bostick transformed Z) 
	from the 1D and 2D parts of an impedance tensor array Z.

	input:
	- Z

	output:
	- Z_nb

	The calculation of the Z_nb needs 6 steps:

	1) Determine the dimensionality of the Z(T), discard all 3D parts
	2) Rotate all Z(T) to TE/TM setup (T_parallel/T_ortho)
	3) Transform every component individually by Niblett-Bostick
	4) collect the respective 2 components each for equal/similar depths
	5) interprete them as TE_nb/TM_nb
	6) set up Z_nb(depth)

	If 1D layers occur inbetween 2D layers, the strike angle is undefined therein.
	We take an - arbitrarily chosen - linear interpolation of strike angle for 
	these layers, with the values varying between the angles of the bounding 
	upper and lower 2D layers (linearly w.r.t. the periods).   

	Use the output for instance for the determination of 
	NB-transformed phase tensors.

	Note:
	No propagation of errors implemented yet!

	"""


	#deal with inputs
	#if zobject:
	z = z_object.z
	periods = 1./z_object.freq
	#else:
	z = z 
	periods = periods


	dimensions = MTge.dimensionality()
	angles = MTge.strike_angle(z)

	#reduce actual Z by the 3D layers:
	z2 = z[np.where(dimensions != 3)[0]]
	angles2 = angles[np.where(dimensions != 3)[0]]
	periods2 = periods[np.where(dimensions != 3)[0]]



	return Z_nb
