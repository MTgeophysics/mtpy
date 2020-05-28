#!/usr/bin/env python

"""
mtpy/mtpy/analysis/geometry.py

Contains classes and functions for handling geometry analysis of impedance tensors:

dimensionality, strike directions, alphas/skews/...

    * 1d - 2d : excentricity of ellipses
    * 2d - 3d : skew < threshold (to be given as argument)
    * strike: frequency - depending angle (incl. 90degree ambiguity)

@UofA, 2013(LK)

Edited by JP, 2016

"""

# =================================================================
import numpy as np

import mtpy.analysis.pt as MTpt
import mtpy.core.z as MTz
import mtpy.utils.exceptions as MTex


# =================================================================

def dimensionality(z_array=None, z_object=None, pt_array=None,
                   pt_object=None, skew_threshold=5,
                   eccentricity_threshold=0.1):
    """
    Esitmate dimensionality of an impedance tensor, frequency by frequency.

    Dimensionality is estimated from the phase tensor given the threshold
    criteria on the skew angle and eccentricity following Bibby et al., 2005
    and Booker, 2014.

    Arguments
    ------------

        **z_array** : np.ndarray(nf, 2, 2)
                      numpy array of impedance elements
                      *default* is None

        **z_object** : mtpy.core.z.Z
                       z_object
                       *default* is None

        **pt_array** : np.ndarray(nf, 2, 2)
                       numpy array of phase tensor elements
                       *default* is None

        **pt_object** : mtpy.analysis.pt.PT
                        phase tensor object
                        *default* is None

        **skew_threshold** : float
                             threshold on the skew angle in degrees, anything
                             above this value is 3-D or azimuthally anisotropic
                             *default* is 5 degrees

        **eccentricity_threshold** : float
                                     threshold on eccentricty in dimensionaless
                                     units, anything below this value is 1-D
                                     *default* is 0.1

    Returns
    ----------

        **dimensions** : np.ndarray(nf, dtype=int)
                         an array of dimesions for each frequency
                         the values are [ 1 | 2 | 3 ]


    Examples
    ----------
        :Estimate Dimesions: ::

            >>> import mtpy.analysis.geometry as geometry
            >>> dim = geometry.dimensionality(z_object=z_obj,
            >>>                               skew_threshold=3)


    """

    lo_dimensionality = []

    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array=z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z(
                'Input argument is not an instance of the Z class')
        pt_obj = MTpt.PhaseTensor(z_object=z_object)
    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array=pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT(
                'Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object

    # use criteria from Bibby et al. 2005 for determining the dimensionality
    # for each frequency of the pt/z array:
    for idx_f in range(len(pt_obj.pt)):
        #1. determine skew value...
        skew = pt_obj.beta[idx_f]
            #compare with threshold for 3D
        if np.abs(skew) > skew_threshold:
            lo_dimensionality.append(3)
        else:
            # 2.check for eccentricity:
            ecc = pt_obj._pi1()[0][idx_f] / pt_obj._pi2()[0][idx_f]
            if ecc > eccentricity_threshold:
                lo_dimensionality.append(2)
            else:
                lo_dimensionality.append(1)

    return np.array(lo_dimensionality)


def strike_angle(z_array=None, z_object=None, pt_array=None,
                 pt_object=None, skew_threshold=5,
                 eccentricity_threshold=0.1):
    """
    Estimate strike angle from 2D parts of the phase tensor given the
    skew and eccentricity thresholds

        Arguments
    ------------

        **z_array** : np.ndarray(nf, 2, 2)
                      numpy array of impedance elements
                      *default* is None

        **z_object** : mtpy.core.z.Z
                       z_object
                       *default* is None

        **pt_array** : np.ndarray(nf, 2, 2)
                       numpy array of phase tensor elements
                       *default* is None

        **pt_object** : mtpy.analysis.pt.PT
                        phase tensor object
                        *default* is None

        **skew_threshold** : float
                             threshold on the skew angle in degrees, anything
                             above this value is 3-D or azimuthally anisotropic
                             *default* is 5 degrees

        **eccentricity_threshold** : float
                                     threshold on eccentricty in dimensionaless
                                     units, anything below this value is 1-D
                                     *default* is 0.1

    Returns
    ----------

        **strike** : np.ndarray(nf)
                         an array of strike angles in degrees for each frequency
                         assuming 0 is north, and e is 90.  There is a 90
                         degree ambiguity in the angle.


    Examples
    ----------
        :Estimate Dimesions: ::

            >>> import mtpy.analysis.geometry as geometry
            >>> strike = geometry.strike_angle(z_object=z_obj,
            >>>                                skew_threshold=3)

    """

    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array=z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z(
                'Input argument is not an instance of the Z class')
        pt_obj = MTpt.PhaseTensor(z_object=z_object)

    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array=pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT(
                'Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object

    lo_dims = dimensionality(pt_object=pt_obj,
                             skew_threshold=skew_threshold,
                             eccentricity_threshold=eccentricity_threshold)

    lo_strikes = []

    for idx, dim in enumerate(lo_dims):
        if dim == 1:
            lo_strikes.append((np.nan, np.nan))
#            continue
        
        elif dim == 3:
            lo_strikes.append((np.nan, np.nan))

        else:
            a = pt_obj.alpha[idx]
            b = pt_obj.beta[idx]
    
            strike1 = (a - b) % 180
            
            # change so that values range from -90 to +90
            # add alternative strikes to account for ambiguity
            if strike1 > 90:
                strike1 -= 180
                strike2 = strike1 + 90
            else:
                strike2 = strike1 - 90
        
            lo_strikes.append((strike1, strike2))

    return np.array(lo_strikes)


def eccentricity(z_array=None, z_object=None, pt_array=None, pt_object=None):
    """
    Estimate eccentricy of a given impedance or phase tensor object


    Arguments
    ------------

        **z_array** : np.ndarray(nf, 2, 2)
                      numpy array of impedance elements
                      *default* is None

        **z_object** : mtpy.core.z.Z
                       z_object
                       *default* is None

        **pt_array** : np.ndarray(nf, 2, 2)
                       numpy array of phase tensor elements
                       *default* is None

        **pt_object** : mtpy.analysis.pt.PT
                        phase tensor object
                        *default* is None


    Returns
    ----------

        **eccentricity** : np.ndarray(nf)


        **eccentricity_err** : np.ndarray(nf)



    Examples
    ----------
        :Estimate Dimesions: ::

            >>> import mtpy.analysis.geometry as geometry
            >>> ec, ec_err= geometry.eccentricity(z_object=z_obj)
    """

    if z_array is not None:
        pt_obj = MTpt.PhaseTensor(z_array=z_array)
    elif z_object is not None:
        if not isinstance(z_object, MTz.Z):
            raise MTex.MTpyError_Z(
                'Input argument is not an instance of the Z class')
        pt_obj = MTpt.PhaseTensor(z_object=z_object)
    elif pt_array is not None:
        pt_obj = MTpt.PhaseTensor(pt_array=pt_array)
    elif pt_object is not None:
        if not isinstance(pt_object, MTpt.PhaseTensor):
            raise MTex.MTpyError_PT(
                'Input argument is not an instance of the PhaseTensor class')
        pt_obj = pt_object

    lo_ecc = []
    lo_eccerr = []

    if not isinstance(pt_obj, MTpt.PhaseTensor):
        raise MTex.MTpyError_PT(
            'Input argument is not an instance of the PhaseTensor class')

    for idx_f in range(len(pt_obj.pt)):
        lo_ecc.append(pt_obj._pi1()[0][idx_f] / pt_obj._pi2()[0][idx_f])

        ecc_err = None
        if (pt_obj._pi1()[1] is not None) and (pt_obj._pi2()[1] is not None):
            ecc_err = np.sqrt((pt_obj._pi1()[1][idx_f] / pt_obj._pi1()[0][idx_f]) ** 2 +\
                              (pt_obj._pi2()[1][idx_f] / pt_obj._pi2()[0][idx_f]) ** 2)

        lo_eccerr.append(ecc_err)

    return np.array(lo_ecc), np.array(lo_eccerr)*np.array(lo_ecc)
