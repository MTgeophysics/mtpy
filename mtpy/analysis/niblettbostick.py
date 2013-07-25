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
import copy

reload(MTz)


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
	rho_nb = rho * (np.pi/2/np.deg2rad(phase%90) - 1)
	#print rho,period,depth,rho_nb
	#print 'rho: {0:.1f} \t-\t rhoNB: {3:.1f}\t-\t period: {1:.1f} \t-\t depth: {2:.1f}'.format(
	#												rho,period,depth,rho_nb)

	return rho_nb, depth


def calculate_znb(z_object = None, z_array = None, periods = None):
	"""
	Determine an array of Z_nb (depth dependent Niblett-Bostick transformed Z) 
	from the 1D and 2D parts of an impedance tensor array Z.

	input:
	- Z (object or array)
	- periods (mandatory, if Z is just array)

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
	#z = z_object.z
	#periods = 1./z_object.freq
	#else:
	z = z_array 
	periods = periods


	dimensions = MTge.dimensionality(z)
	angles = MTge.strike_angle(z)

	#reduce actual Z by the 3D layers:
	z2 = z[np.where(dimensions != 3)[0]]
	angles2 = angles[np.where(dimensions != 3)[0]]
	periods2 = periods[np.where(dimensions != 3)[0]]

	angles_incl1D = interpolate_strike_angles(angles2[:,0],periods2)
	
	z3 = MTz.rotate_z(z2,-angles_incl1D)[0]

	#at this point we assume that the two modes are the off-diagonal elements!!
	#TE is element (1,2), TM at (2,1)
	lo_nb_max = []
	lo_nb_min = []

	app_res = MTz.z2resphi(z3,periods2)[0]
	phase = MTz.z2resphi(z3,periods2)[1]


	for i,per in enumerate(periods):
		
		te_rho, te_depth = rhophi2rhodepth(app_res[i][0,1], phase[i][0,1], per)
		tm_rho, tm_depth = rhophi2rhodepth(app_res[i][1,0], phase[i][1,0], per)
		
		if te_rho > tm_rho:
			lo_nb_max.append([te_depth, te_rho])
			lo_nb_min.append([tm_depth, tm_rho])
		else:
			lo_nb_min.append([te_depth, te_rho])
			lo_nb_max.append([tm_depth, tm_rho])




	return np.array(lo_nb_max), np.array(lo_nb_min)


def calculate_rho_minmax(z_object = None, z_array = None, periods = None):
	"""
	Determine 2 arrays of Niblett-Bostick transformed aparent resistivities: 
	minumum and maximum values for respective periods. 

	Values are calculated from the 1D and 2D parts of an impedance tensor array Z.

	input:
	- Z (object or array)
	- periods (mandatory, if Z is just array)

	output:
	- n x 3 array, depth/rho_nb/angle for rho_nb max
	- n x 3 array, depth/rho_nb/angle for rho_nb min

	The calculation is carried out by :
	
	1) Determine the dimensionality of the Z(T), discard all 3D parts
	2) loop over periods 
       * rotate Z and calculate app_res_NB for off-diagonal elements
       * find maximum and minimum values
       * write out respective depths and rho values  


	Note:
	No propagation of errors implemented yet!

	"""

	#deal with inputs
	#if zobject:
	#z = z_object.z
	#periods = 1./z_object.freq
	#else:
	z = z_array 
	periods = periods


	dimensions = MTge.dimensionality(z)
	angles = MTge.strike_angle(z)

	#reduce actual Z by the 3D layers:
	z2 = z[np.where(dimensions != 3)[0]]
	angles2 = angles[np.where(dimensions != 3)[0]]
	periods2 = periods[np.where(dimensions != 3)[0]]

	lo_nb_max = []
	lo_nb_min = []

	rotsteps = 360
	rotangles = np.arange(rotsteps)*180./rotsteps


	for i,per in enumerate(periods2):
		z_curr = z2[i]
		temp_vals = np.zeros((rotsteps,4))

		for j,d in enumerate(rotangles):
			new_z = MTcc.rotatematrix_incl_errors(z_curr, d)[0]
			#print i,per,j,d

			res = MTz.z2resphi(new_z,per)[0]
			phs = MTz.z2resphi(new_z,per)[1]

			te_rho, te_depth = rhophi2rhodepth(res[0,1], phs[0,1], per)
			tm_rho, tm_depth = rhophi2rhodepth(res[1,0], phs[1,0], per)

			temp_vals[j,0] = te_depth
			temp_vals[j,1] = te_rho
			temp_vals[j,2] = tm_depth
			temp_vals[j,3] = tm_rho
		
		column = (np.argmax([ np.max(temp_vals[:,1]),
								np.max(temp_vals[:,3])]))*2 + 1

		maxidx = np.argmax(temp_vals[:,column])
		max_rho = temp_vals[maxidx,column]
		max_depth = temp_vals[maxidx,column-1]
		max_ang = rotangles[maxidx]
		
		#alternative 1
		min_column = (np.argmin([ np.max(temp_vals[:,1]),
								np.max(temp_vals[:,3])]))*2 + 1
		if max_ang <= 90:
			min_ang = max_ang + 90
		else:
			min_ang = max_ang - 90
		minidx = np.argmin(np.abs(rotangles-min_ang))

		min_rho = temp_vals[minidx,min_column]
		min_depth = temp_vals[minidx,min_column-1]


		lo_nb_max.append([max_depth, max_rho, max_ang])
		lo_nb_min.append([min_depth, min_rho])
		

	return np.array(lo_nb_max), np.array(lo_nb_min)


def interpolate_strike_angles(angles,in_periods):

	"""
	expect 2 arrays

	1. sort ascending by periods
	2. loop over angles to find 'nan' values (i.e. 1D layers)
	3. determine linear interpolation between bounding 2D strike angles
	4. if 1D on top or bottom, set to 0 degrees
	"""

	new_angles = copy.copy(angles)

	#sort in ascending order: 
	orig_sorting = np.argsort(in_periods)
	back_sorting = np.argsort(orig_sorting)
	angles = angles[orig_sorting]
	periods = in_periods[orig_sorting]

	in_line = 0 
	while in_line < len(angles):
		curr_per = periods[in_line]
		curr_ang = angles[in_line]
		if np.isnan(curr_ang):

			if in_line in [0,len(angles)-1]:
				new_angles[in_line] = 0.
				in_line += 1
				continue
		
			#otherwise:

			#use value before current one:
			ang1 = new_angles[in_line - 1]
			per1 = periods[in_line - 1]
			#check for next non-nan:
			ang2 = None

			j = in_line
			while j < len(angles):
				j += 1 
				if not np.isnan(angles[j]):
					per2 = periods[j]
					ang2 = angles[j]
					break
			#catch case of all nan:
			if ang2 is None:
				ang2 = 0.  

			delta_per = per2-per1
			delta_ang = ang2-ang1
			per_step = periods[in_line] - per1
			new_angles[in_line] = ang1 + delta_ang/delta_per * per_step
												
		in_line += 1

	#asserting correct order (same as input) of the angles:
	return new_angles[back_sorting]