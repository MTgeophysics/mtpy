#!/usr/bin/env python

"""
This module contains helper functions for standard calculations. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time
import ConfigParser
import math, cmath

from mtpy.utils.exceptions import *
import mtpy.utils.format as MTformat
#=================================================================

#define uncertainty for differences between time steps
epsilon = 1e-9


#=================================================================


def invertmatrix_incl_errors(inmatrix, inmatrix_err = None):

    if inmatrix is None:
        raise MTexceptions.MTpyError_inputarguments('Matrix must be defined')

    if (inmatrix_err is not None) and (inmatrix.shape != inmatrix_err.shape):
        raise MTexceptions.MTpyError_inputarguments('Matrix and err-matrix shapes do not match: %s - %s'%(str(inmatrix.shape), str(inmatrix_err.shape)))

    if (inmatrix.shape[-2] != inmatrix.shape[-1]):
        raise MTexceptions.MTpyError_inputarguments('Matrices must be square!')

    if (inmatrix_err is not None) and (inmatrix_err.shape[-2] != inmatrix_err.shape[-1]) :
        raise MTexceptions.MTpyError_inputarguments('Matrices must be square!')

    dim = inmatrix.shape[-1]


    det = np.linalg.det(inmatrix)

    if det == 0:
        raise MTexceptions.MTpyError_inputarguments('Matrix is singular - I cannot invert that!')

    inv_matrix = np.zeros_like(inmatrix)



    if dim != 2:
        raise MTexceptions.MTpyError_inputarguments('Only 2D matrices supported yet')

    inv_matrix = np.linalg.inv(inmatrix)

    inv_matrix_err = None

    if (inmatrix_err is not None):
        inmatrix_err = np.real(inmatrix_err)
        inv_matrix_err = np.zeros_like(inmatrix_err)

        for i in range(2):
            for j in range(2):
                #looping over the entries of the error matrix
                err = 0.
                for k in range(2):
                    for l in range(2):
                        #each entry has 4 summands

                        err += np.abs (- inv_matrix[i,k]  * inv_matrix[l,j]  * inmatrix_err[k,l])


                inv_matrix_err[i,j] = err

 
    return inv_matrix, inv_matrix_err


def propagate_error_rect2polar(x,x_error,y, y_error):
    
    # x_error, y_error define a  rectangular uncertainty box  
    
    # rho error is the difference between the closest and furthest point of the box (w.r.t. the origin)
    # approximation: just take corners and midpoint of edges 
    lo_points = [ (x + x_error, y), (x - x_error, y), (x, y - y_error ), (x, y - y_error ),\
                  (x - x_error, y - y_error) ,(x + x_error, y - y_error) ,(x + x_error, y + y_error) ,(x - x_error, y + y_error) ]


    #check, if origin is within the box:
    origin_in_box = False
    if x_error >= np.abs(x) and y_error >= np.abs(y):
        origin_in_box = True

    lo_polar_points = [ cmath.polar(np.complex(*i)) for i in lo_points ]
    lo_rho = [i[0] for i in lo_polar_points ]
    lo_phi = [i[1] for i in lo_polar_points ]

    rho_err = lo_rho.max() 
    phi_err = math.degrees(  lo_phi.max() - lo_phi.min() )

    if origin_in_box is True:
        rho_err -= lo_rho.min()
        phi_err = 360.

    return rho_error, phi_error



#rotation:
#1. rotation positive in clockwise direction
#2. orientation of new X axis X' given by rotation angle
#3. express contents of Z/tipper (points P) in this new system (points P')
#4. rotation for points calculated as P' = ([cos , sin ],[-sin, cos]) * P <=> P' = R * P
#5. => B' = R * B and E' = R * E
# (Rt is the inverse rotation matrix)
#6. E = Z * B => Rt * E' = Z * Rt * B' => E' = (R*Z*Rt) * B' => Z' = (R*Z*Rt)  

#7. Bz = T * B => Bz = T * Rt * B' => T' = (T * Rt)
#(Since Tipper is a row vector)
#7.a for general column vectors v:  v' = R * v

# Rotation of the uncertainties:
# a) rotate Z into Z'
# b) use propagation of errors on Z' to obtain the rotated Z'err
# That is NOT the same as the rotated error matrix Zerr (although the result is similar)


def rotatematrix_incl_errors(inmatrix, angle, inmatrix_err = None) :
   
    if inmatrix is None :
        raise MTexceptions.MTpyError_inputarguments('Matrix AND eror matrix must be defined')

    if (inmatrix_err is not None) and (inmatrix.shape != inmatrix_err.shape):
        raise MTexceptions.MTpyError_inputarguments('Matrix and err-matrix shapes do not match: %s - %s'%(str(inmatrix.shape), str(inmatrix_err.shape)))


    try:
        degreeangle = angle%360
    except:
        raise MTexceptions.MTpyError_inputarguments('"Angle" must be a valid number (in degrees)')

    phi = math.radians(degreeangle)
    
    cphi = np.cos(phi)
    sphi = np.sin(phi)

    rotmat = np.matrix([[ cphi,-sphi],[sphi,cphi] ])

    rotated_matrix = np.dot(np.dot( rotmat, inmatrix ),rotmat.I  ) 


    errmat  = None
    if (inmatrix_err is not None) :
        err_orig = np.real(inmatrix_err) 
        errmat = np.zeros_like(inmatrix_err)

        # squared propagation of errors
        # zerr_rot[idx_freq,0,0] = np.sqrt( (cphi**2 * zerr_orig[0,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[0,1])**2 + (zerr_orig[1,0])**2) + (sphi**2 * zerr_orig[1,1])**2)

        # zerr_rot[idx_freq,0,1] = np.sqrt( (cphi**2 * zerr_orig[0,1])**2 + cphi**2 * sphi**2 * ( (zerr_orig[1,1])**2 + (zerr_orig[0,0])**2) + (sphi**2 * zerr_orig[1,0])**2) 

        # zerr_rot[idx_freq,1,0] = np.sqrt( (cphi**2 * zerr_orig[1,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[1,1])**2 + (zerr_orig[0,0])**2) + (sphi**2 * zerr_orig[0,1])**2) 

        # zerr_rot[idx_freq,1,1] = np.sqrt( (sphi**2 * zerr_orig[0,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[0,1])**2 + (zerr_orig[1,0])**2) + (cphi**2 * zerr_orig[1,1])**2) 


        errmat[0,0] = cphi**2 * err_orig[0,0] + np.abs(cphi * sphi) * (err_orig[0,1] + err_orig[1,0]) + sphi**2 * err_orig[1,1]
        errmat[0,1] = cphi**2 * err_orig[0,1] + np.abs(cphi * sphi) * (err_orig[1,1] + err_orig[0,0]) + sphi**2 * err_orig[1,0]
        errmat[1,0] = cphi**2 * err_orig[1,0] + np.abs(cphi * sphi) * (err_orig[1,1] + err_orig[0,0]) + sphi**2 * err_orig[0,1]
        errmat[1,1] = cphi**2 * err_orig[1,1] + np.abs(cphi * sphi) * (err_orig[0,1] + err_orig[1,0]) + sphi**2 * err_orig[0,0] 

    return rotated_matrix, errmat


def rotatevector_incl_errors(invector, angle, invector_err = None):
    #check for row or column vector 
    
    if invector is None :
        raise MTexceptions.MTpyError_inputarguments('Vector AND error-vector must be defined')

    if (invector_err is not None) and (invector.shape != invector_err.shape):
        raise MTexceptions.MTpyError_inputarguments('Vector and errror-vector shapes do not match: %s - %s'%(str(invector.shape), str(invector_err.shape)))

    try:
        degreeangle = angle%360
    except:
        raise MTexceptions.MTpyError_inputarguments('"Angle" must be a valid number (in degrees)')

    phi = math.radians(degreeangle)
    
    cphi = np.cos(phi)
    sphi = np.sin(phi)

    rotmat = np.matrix([[ cphi,-sphi],[sphi,cphi] ])

    rotated_vector = np.dot( rotmat, invector )
    
    errvec = None
    if (invector_err is not None) :   
        err_orig = np.real(invector_err)
        errvec = np.zeros_like(invector_err)

        errvec[0,0] = cphi**2 * err_orig[0,0] + np.abs(cphi * sphi) * (err_orig[0,1] + err_orig[1,0]) + sphi**2 * err_orig[1,1]
        errvec[0,1] = cphi**2 * err_orig[0,1] + np.abs(cphi * sphi) * (err_orig[1,1] + err_orig[0,0]) + sphi**2 * err_orig[1,0]


    return rotated_matrix, errvec
