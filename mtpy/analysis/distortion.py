#!/usr/bin/env python

"""
mtpy/analysis/distortion.py

Contains functions for the determination of (galvanic) distortion of impedance tensors.
The methods used follow Bibby et al 2005.
As it has been pointed out in that paper, there are various possibilities for
constraining the solution, esp. in the 2D case.

Here we just implement the 'most basic' variety for the calculation of the 
distortion tensor. Other methods can be implemented, but since the optimal assumptions and
constraints depend on the application, the actual place for further functions
is in an independent, personalised module.

Algorithm Details:
Finding the distortion of a Z array. Using the phase tensor
so, Z arrays are transformed into PTs first), following Bibby et al. 2005.

First, try to find periods that indicate 1D. From them determine D incl.
the g-factor by calculatiing a weighted mean. The g is assumed in order to
cater for the missing unknown in the system, it is here set to det(X)^0.5.
After that is found, the function no_distortion from the Z module can be
called to obtain the unperturbated regional impedance tensor.

Second, if there are no 1D sections: Find the strike angle, then rotate the
Z to the principal axis. In order to do that, use the rotate(-strike) method
of the Z module. Then take the real part of the rotated Z. As in the 1D case,
we need an assumption to get rid of the (2) unknowns:
set det(D) = P and det(D) = T, where P,T can be chosen. Common choice is to
set  one of P,T to an arbitrary value (e.g. 1). Then check, for which values
of the other parameter  S^2 = T^2+4*P*X_12*X_21/det(X) > 0 holds.

@UofA, 2013  (LK)

Edited by JP, 2016

"""

import copy

import numpy as np

import mtpy.analysis.geometry as MTge
import mtpy.core.z as MTz
import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex


def find_distortion(z_object, g='det', num_freq=None, lo_dims=None):
    """
    find optimal distortion tensor from z object

    automatically determine the dimensionality over all frequencies, then find
    the appropriate distortion tensor D
    
    Parameters
    ----------
    
        **z_object** : mtpy.core.z object
                       
        **g** : [ 'det' | '01' | '10 ]
                type of distortion correction
                *default* is 'det'
                
        **num_freq** : int
                       number of frequencies to look for distortion from 
                       the index 0
                       *default* is None, meaning all frequencies are used
                       
        **lo_dims** : list
                      list of dimensions for each frequency
                      *default* is None, meaning calculated from data
                      
    Returns
    -------
    
        **distortion** : np.ndarray(2, 2)
                         distortion array all real values
        
        **distortion_err** : np.ndarray(2, 2)
                             distortion error array


    Examples
    --------

        :Estimate Distortion: ::
        
            >>> import mtpy.analysis.distortion as distortion
            >>> dis, dis_err = distortion.find_distortion(z_obj, num_freq=12)
            
    """

    if num_freq is not None:
        if num_freq > z_object.freq.size:
            num_freq = z_object.freq.size
            print('Number of frequencies to sweep over is too high for z')
            print('setting num_freq to {0}'.format(num_freq))
    else:
        num_freq = z_object.freq.size


    z_obj = MTz.Z(z_object.z[0:num_freq],
                  z_object.z_err[0:num_freq],
                  z_object.freq[0:num_freq])

    g = 'det'

    dim_arr = MTge.dimensionality(z_object=z_obj)
    st_arr = -1 * MTge.strike_angle(z_object=z_obj)[:, 0]

    dis = np.zeros_like(z_obj.z, dtype=np.float)
    dis_err = np.ones_like(z_obj.z, dtype=np.float)

    #dictionary of values that should be no distortion in case distortion
    #cannot be calculated for that component
    rot_mat = np.matrix([[0, -1], [1, 0]])
    for idx, dim in enumerate(dim_arr):
        if np.any(z_obj.z[idx] == 0.0 + 0.0j) == True:
            dis[idx] = np.identity(2)
            print('Found a zero in z at {0}, skipping'.format(idx))
            continue

        if dim == 1:

            if g in ['01', '10']:
                gr = np.abs(z_obj.z.real[idx, int(g[0]), int(g[1])])
                gi = np.abs(z_obj.z.imag[idx, int(g[0]), int(g[1])])
            else:
                gr = np.sqrt(np.linalg.det(z_obj.z.real[idx]))
                gi = np.sqrt(np.linalg.det(z_obj.z.imag[idx]))

            dis[idx] = np.mean(np.array([(1. / gr * np.dot(z_obj.z.real[idx],
                                                           rot_mat)),
                                         (1. / gi * np.dot(z_obj.z.imag[idx],
                                                           rot_mat))]),
                               axis=0)

            if z_obj.z_err is not None:
                # find errors of entries for calculating weights

                gr_err = 1. / gr * np.abs(z_obj.z_err[idx])
                gr_err[np.where(gr_err == 0.0)] = 1.0

                gi_err = 1. / gi * np.abs(z_obj.z_err[idx])
                gi_err[np.where(gi_err == 0.0)] = 1.0

                dis_err[idx] = np.mean(np.array([gi_err, gr_err]),
                                       axis=0)

        elif dim == 2:
            P = 1
            strike_ang = st_arr[idx]
            if np.isnan(strike_ang):
                strike_ang = 0.0

            if z_obj.z_err is not None:
                err_arr = z_obj.z_err[idx]
                err_arr[np.where(err_arr == 0.0)] = 1.0
            else:
                err_arr = None

            tetm_arr, tetm_err = MTcc.rotatematrix_incl_errors(z_obj.z[idx],
                                                               strike_ang,
                                                               inmatrix_err=err_arr)

            tetm_r = tetm_arr.real
            tetm_i = tetm_arr.imag
            t_arr_r = -4 * P * tetm_r[0, 1] * tetm_r[1, 0] / np.linalg.det(tetm_r)
            t_arr_i = -4 * P * tetm_i[0, 1] * tetm_i[1, 0] / np.linalg.det(tetm_i)

            try:
                T = np.sqrt(max([t_arr_r, t_arr_i])) + .001
            except ValueError:
                T = 2

            sr = np.sqrt(T ** 2 + 4 * P * tetm_r[0, 1] * tetm_r[1, 0] / np.linalg.det(tetm_r))
            si = np.sqrt(T ** 2 + 4 * P * tetm_i[0, 1] * tetm_i[1, 0] / np.linalg.det(tetm_i))

            par_r = 2 * tetm_r[0, 1] / (T - sr)
            orth_r = 2 * tetm_r[1, 0] / (T + sr)
            par_i = 2 * tetm_i[0, 1] / (T - si)
            orth_i = 2 * tetm_i[1, 0] / (T + si)

            mat2_r = np.matrix([[0, 1. / orth_r], [1. / par_r, 0]])
            mat2_i = np.matrix([[0, 1. / orth_i], [1. / par_i, 0]])

            avg_mat = np.mean(np.array([np.dot(tetm_r, mat2_r),
                                        np.dot(tetm_i, mat2_i)]),
                              axis=0)

            dis[idx] = avg_mat

            if err_arr is not None:
                # find errors of entries for calculating weights
                sigma_sr = np.sqrt((-(2 * P * tetm_r[0, 1] * tetm_r[1, 0] * \
                                      tetm_r[1, 1] * err_arr[0, 0]) / \
                                    (np.linalg.det(tetm_r) ** 2 * sr)) ** 2 + \
                                   ((2 * P * tetm_r[0, 0] * tetm_r[1, 0] *
                                     tetm_r[1, 1] * err_arr[0, 1]) /
                                    (np.linalg.det(tetm_r) ** 2 * sr)) ** 2 + \
                                   ((2 * P * tetm_r[0, 0] * tetm_r[0, 1] *
                                     tetm_r[1, 1] * err_arr[1, 0]) / \
                                    (np.linalg.det(tetm_r) ** 2 * sr)) ** 2 + \
                                   (-(2 * P * tetm_r[0, 1] * tetm_r[1, 0] * \
                                      tetm_r[0, 0] * err_arr[1, 1]) / \
                                    (np.linalg.det(tetm_r) ** 2 * sr)) ** 2)

                sigma_dr_11 = 0.5 * sigma_sr
                sigma_dr_22 = 0.5 * sigma_sr

                sigma_dr_12 = np.sqrt((mat2_r[0, 1] / tetm_r[0, 0] * err_arr[0, 0]) ** 2 + \
                                      (mat2_r[0, 1] / tetm_r[1, 0] * err_arr[1, 0]) ** 2 + \
                                      (0.5 * tetm_r[0, 0] / tetm_r[1, 0] * sigma_sr) ** 2)
                sigma_dr_21 = np.sqrt((mat2_r[1, 0] / tetm_r[1, 1] * err_arr[1, 1]) ** 2 + \
                                      (mat2_r[1, 0] / tetm_r[0, 1] * err_arr[0, 1]) ** 2 + \
                                      (0.5 * tetm_r[1, 1] / tetm_r[0, 1] * sigma_sr) ** 2)

                dis_err_r = np.array([[sigma_dr_11, sigma_dr_12],
                                      [sigma_dr_21, sigma_dr_22]])

                sigma_si = np.sqrt((-(2 * P * tetm_i[0, 1] * tetm_i[1, 0] * \
                                      tetm_i[1, 1] * err_arr[0, 0]) / \
                                    (np.linalg.det(tetm_i) ** 2 * sr)) ** 2 + \
                                   ((2 * P * tetm_i[0, 0] * tetm_i[1, 0] * \
                                     tetm_i[1, 1] * err_arr[0, 1]) / \
                                    (np.linalg.det(tetm_i) ** 2 * sr)) ** 2 + \
                                   ((2 * P * tetm_i[0, 0] * tetm_i[0, 1] * \
                                     tetm_i[1, 1] * err_arr[1, 0]) / \
                                    (np.linalg.det(tetm_i) ** 2 * sr)) ** 2 + \
                                   (-(2 * P * tetm_i[0, 1] * tetm_i[1, 0] * \
                                      tetm_i[0, 0] * err_arr[1, 1]) / \
                                    (np.linalg.det(tetm_i) ** 2 * sr)) ** 2)

                sigma_di_11 = 0.5 * sigma_si
                sigma_di_22 = 0.5 * sigma_si
                sigma_di_12 = np.sqrt((mat2_i[0, 1] / tetm_i[0, 0] * err_arr[0, 0]) ** 2 + \
                                      (mat2_i[0, 1] / tetm_i[1, 0] * err_arr[1, 0]) ** 2 + \
                                      (0.5 * tetm_i[0, 0] / tetm_i[1, 0] * sigma_si) ** 2)
                sigma_di_21 = np.sqrt((mat2_i[1, 0] / tetm_i[1, 1] * err_arr[1, 1]) ** 2 + \
                                      (mat2_i[1, 0] / tetm_i[0, 1] * err_arr[0, 1]) ** 2 + \
                                      (0.5 * tetm_i[1, 1] / tetm_i[0, 1] * sigma_si) ** 2)

                dis_err_i = np.array([[sigma_di_11, sigma_di_12],
                                      [sigma_di_21, sigma_di_22]])

                dis_err[idx] = np.mean(np.array([dis_err_r, dis_err_i]))
        else:
            dis[idx] = np.identity(2)

    nonzero_idx = np.array(list(set(np.nonzero(dis)[0])))

    dis_avg, weights_sum = np.average(dis[nonzero_idx],
                                      axis=0,
                                      weights=(1. / dis_err[nonzero_idx]) ** 2,
                                      returned=True)

    dis_avg_err = np.sqrt(1. / weights_sum)

    return dis_avg, dis_avg_err


def find_1d_distortion(z_object, include_non1d=False):
    """
    find 1D distortion tensor from z object

    ONly use the 1D part of the Z to determine D. 
    Treat all frequencies as 1D, if  "include_non1d = True".

    
    """

    if not isinstance(z_object, MTz.Z):
        raise MTex.MTpyError_inputarguments('first argument must be an '
                                            'instance of the Z class')

    z_obj = z_object

    lo_dims = MTge.dimensionality(z_object=z_obj)

    if include_non1d is True:
        lo_dims = [1 for i in lo_dims]

    if len(list(np.where(np.array(lo_dims) == 1))) == 0:
        raise MTex.MTpyError_inputarguments('Z object does not have '
                                            'frequencies with spatial 1D characteristic')

    print(lo_dims)

    return find_distortion(z_obj, lo_dims=lo_dims)


def find_2d_distortion(z_object, include_non2d=False):
    """
    find 2D distortion tensor from z object

    ONly use the 2D part of the Z to determine D. 
    Treat all frequencies as 2D, if  "include_non2d = True".
    
    """

    if not isinstance(z_object, MTz.Z):
        raise MTex.MTpyError_inputarguments('first argument must be an '
                                            'instance of the Z class')

    z_obj = z_object

    lo_dims = MTge.dimensionality(z_object=z_obj)

    # avoid the (standard) 1D distortion call -> remove all 1
    lo_dims = [4 if i == 1 else i for i in lo_dims]

    if include_non2d is True:
        lo_dims = [2 for i in lo_dims]

    if len(list(np.where(np.array(lo_dims) == 2))) == 0:
        raise MTex.MTpyError_inputarguments('Z object does not have'
                                            ' frequencies with spatial 2D characteristic')

    return find_distortion(z_obj, lo_dims=lo_dims)


def remove_distortion(z_array=None, z_object=None, num_freq=None, g='det'):
    """
    remove distortion from an impedance tensor using the method outlined by 
    Bibby et al., [2005].
    
    Parameters
    -----------
    
        **z_array** : np.ndarray((nf, 2, 2))
                      numpy array of impedance tensor
                      *default* is None
                      
        **z_object** : mtpy.core.z object
                       *default* is None

        **num_freq** : int
                       number of frequecies to look for distortion
                       *default* is None, meaning look over all frequencies

        **g** : [ 'det' | '01' | '10 ]
                type of distortion to look for
                *default* is 'det'

    Returns
    --------

        **distortion** : np.ndarray (2, 2)
                         distortion array

        **new_z_obj** : mtpy.core.z
                        z object with distortion removed and error calculated

    Examples
    -------------

        :Remove Distortion: ::

            >>> import mtpy.analysis.distortion as distortion
            >>> d, new_z = distortion.remove_distortion(z_object=z_obj)                         
    """

    if z_array is not None:
        z_obj = MTz.Z(z_array=z_array)

    elif z_object is not None:
        z_obj = z_object

    zero_idx = np.where(z_obj.z == 0 + 0j)

    # 0. generate a Z object
    # 1. find distortion via function above,
    # 2. remove distortion via method of z object

    dis, dis_err = find_distortion(z_obj, num_freq=num_freq, g=g)

    try:
        # distortion_tensor, zd, zd_err = z_obj.no_distortion(dis, distortion_err_tensor=dis_err)
        distortion_tensor, zd, zd_err = z_obj.remove_distortion(dis, distortion_err_tensor=dis_err)

        zd_err = np.nan_to_num(zd_err)
        zd_err[np.where(zd_err == 0.0)] = 1.0
        distortion_z_obj = z_obj
        distortion_z_obj.z = zd
        distortion_z_obj.z[zero_idx] = 0.0 + 0.0j
        distortion_z_obj.z_err = zd_err

        return distortion_tensor, distortion_z_obj

    except MTex.MTpyError_Z:
        print('Could not compute distortion tensor')

        return np.identity(2), z_obj
