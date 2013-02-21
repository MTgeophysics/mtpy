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


def invertmatrix_incl_errors(inmatrix, inmatrix_err):

    if inmatrix is None or inmatrix_err is None:
        raise MTexceptions.MTpyError_inputarguments('Matrix AND eror matrix must be defined')

    if inmatrix.shape != inmatrix_err.shape:
        raise MTexceptions.MTpyError_inputarguments('Matrix and err-matrix shapes do not match: %s - %s'%(str(inmatrix.shape), str(inmatrix_err.shape)))

    if (inmatrix.shape[-2] != inmatrix.shape[-1]) or (inmatrix_err.shape[-2] != inmatrix_err.shape[-1]) :
        raise MTexceptions.MTpyError_inputarguments('Matrices must be square!')

    dim = inmatrix.shape[-1]

    inmatrix_err = np.real(inmatrix_err)

    det = np.linalg.det(inmatrix)

    if det == 0:
        raise MTexceptions.MTpyError_inputarguments('Matrix is singular - I cannot invert that!')

    inv_matrix = np.zeros_like(inmatrix)
    inv_matrix_err = np.zeros_like(inmatrix_err)

    if dim != 2:
        raise MTexceptions.MTpyError_inputarguments('Only 2D matrices supported yet')

    inv_matrix = np.linalg.inv(inmatrix)

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

    if origin_in_box is False:
        rho_err -= lo_rho.min()
        phi_err = 360.

    return rho_error, phi_error
