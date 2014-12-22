#!/usr/bin/env python

"""
mtpy/mtpy/analysis/geometry.py

Contains classes and functions for handling geometry analysis of impedance tensors:

dimensionality, strike directions, alphas/betas/...
 


    Class:

        Methods:



    Functions:


@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np

import mtpy.core.z as MTz 
import mtpy.analysis.pt as MTpt 
import mtpy.utils.exceptions as MTex

# reload(MTex)
# reload(MTz)
reload(MTpt)


#=================================================================


# 1d - 2d : excentricity of ellipses
# 2d - 3d : beta < threshold (to be given as argument)
# strike: frequency - depending angle (incl. 90degree ambiguity)
# input: PT object, Z object (,edi object) 



def dimensionality(z_array = None, z_object = None, pt_array= None, 
                    pt_object = None, beta_threshold = 5, 
                    eccentricity_threshold = 0.1):
    """
    beta_threshold: angle in degrees - if beta is smaller than this, it's 2d
    
    eccentricity_threshold: fraction of eccentricity (0: circle - 1: line) -
    if eccentricity (ellipticity) is small than this, it's a 1D geometry.

    """

    lo_dimensionality = []
    
    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array = z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z('Input argument is not an instance of the Z class')        
        pt_obj = MTpt.PhaseTensor(z_object = z_object)
    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array= pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT('Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object
    

    #use criteria from Bibby et al. 2005 for determining the dimensionality for each frequency of the pt/z array:
    for idx_f in range(len(pt_obj.pt)):
        #1. determine beta value...
        beta = pt_obj.beta[0][idx_f]
            #compare with threshold for 3D
        if beta > beta_threshold:
            lo_dimensionality.append(3)
        else:
            #2.check for eccentricity:
            ecc = pt_obj._pi1()[0][idx_f] / pt_obj._pi2()[0][idx_f]
            if ecc > eccentricity_threshold:
                lo_dimensionality.append(2)
            else:
                lo_dimensionality.append(1)

    return np.array(lo_dimensionality)



def strike_angle(z_array = None, z_object = None, pt_array= None, 
                    pt_object = None, beta_threshold = 5, 
                    eccentricity_threshold = 0.1):

    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array = z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z('Input argument is not an instance of the Z class')
        pt_obj = MTpt.PhaseTensor(z_object = z_object)

    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array= pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT('Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object

    lo_dims =  dimensionality(pt_object = pt_obj, beta_threshold =beta_threshold , eccentricity_threshold = eccentricity_threshold )

    lo_strikes = []

    for idx, dim in enumerate(lo_dims):
        if dim == 1:
            lo_strikes.append((np.nan, np.nan))
            continue

        a = pt_obj.alpha[0][idx]
        b = pt_obj.beta[0][idx]

        strike1 = (a - b)%90
        if 0 < strike1 < 45 :
            strike2 = strike1 + 90
        else:
            strike2 = strike1 - 90
        
        s1 = min(strike1,strike2)
        s2 = max(strike1,strike2)

        lo_strikes.append(( s1,s2) )
        
        

    return np.array(lo_strikes)



def eccentricity(z_array = None, z_object = None, pt_array= None, pt_object = None):


    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array = z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z('Input argument is not an instance of the Z class')
        pt_obj = MTpt.PhaseTensor(z_object = z_object)
    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array= pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT('Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object

    lo_ecc = []
    lo_eccerr = []

    if not isinstance(pt_obj, MTpt.PhaseTensor):
        raise MTex.MTpyError_PT('Input argument is not an instance of the PhaseTensor class')
   

    for idx_f in range(len(pt_obj.pt)):
        lo_ecc.append( pt_obj._pi1()[0][idx_f] / pt_obj._pi2()[0][idx_f] )

        ecc_err = None
        if (pt_obj._pi1()[1] is not None) and (pt_obj._pi2()[1] is not None):
            ecc_err = np.sqrt( (pt_obj._pi1()[1][idx_f] /pt_obj._pi1()[0][idx_f] )**2 + (pt_obj._pi2()[1][idx_f] /pt_obj._pi2()[0][idx_f])**2)  


        lo_eccerr.append(ecc_err)

    return np.array(lo_ecc), np.array(lo_eccerr) 