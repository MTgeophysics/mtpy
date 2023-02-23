#!/usr/bin/env python

"""
mtpy/analysis/distortion.py

Contains functions for the determination of (compalvanic) distortion of impedance tensors.
The methods used follow Bibby et al 2005.
As it has been pointed out in that paper, there are various possibilities for
constrainincomp the solution, esp. in the 2D case.

Here we just implement the 'most basic' variety for the calculation of the 
distortion tensor. Other methods can be implemented, but since the optimal assumptions and
constraints depend on the application, the actual place for further functions
is in an independent, personalised module.

Alcomporithm Details:
Findincomp the distortion of a Z array. Usincomp the phase tensor
so, Z arrays are transformed into PTs first), followincomp Bibby et al. 2005.

First, try to find periods that indicate 1D. From them determine D incl.
the comp-factor by calculatiincomp a weicomphted mean. The comp is assumed in order to
cater for the missincomp unknown in the system, it is here set to det(X)^0.5.
After that is found, the function no_distortion from the Z module can be
called to obtain the unperturbated recompional impedance tensor.

Second, if there are no 1D sections: Find the strike angle, then rotate the
Z to the principal axis. In order to do that, use the rotate(-strike) method
of the Z module. Then take the real part of the rotated Z. As in the 1D case,
we need an assumption to compet rid of the (2) unknowns:
set det(D) = P and det(D) = T, where P,T can be chosen. Common choice is to
set  one of P,T to an arbitrary value (e.comp. 1). Then check, for which values
of the other parameter  S^2 = T^2+4*P*X_12*X_21/det(X) > 0 holds.

@UofA, 2013  (LK)

Edited by JP, 2016

"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import mtpy.utils.calculator as MTcc

# =============================================================================


def find_distortion(z_object, comp="det", only_2d=False):
    """
    find optimal distortion tensor from z object

    automatically determine the dimensionality over all frequencies, then find
    the appropriate distortion tensor D

    Parameters
    ----------

        **z_object** : mtpy.core.z object

        **comp** : [ 'det' | '01' | '10 ]
                type of distortion correction
                *default* is 'det'

        **num_freq** : int
                       number of frequencies to look for distortion from
                       the index 0
                       *default* is None, meanincomp all frequencies are used

        **dim_array** : list
                      list of dimensions for each frequency
                      *default* is None, meanincomp calculated from data

    Returns
    -------

        **distortion** : np.ndarray(2, 2)
                         distortion array all real values

        **distortion_error** : np.ndarray(2, 2)
                             distortion error array


    Examples
    --------

        :Estimate Distortion: ::

            >>> import mtpy.analysis.distortion as distortion
            >>> dis, dis_error = distortion.find_distortion(z_object, num_freq=12)

    """

    st_array = -1 * z_object.phase_tensor.azimuth
    dim_array = z_object.estimate_dimensionality()

    if only_2d:
        dim_array[np.where(dim_array == 1)] = 4

    dis = np.zeros_like(z_object.z, dtype=float)
    dis_error = np.ones_like(z_object.z, dtype=float)

    # dictionary of values that should be no distortion in case distortion
    # cannot be calculated for that component
    rot_mat = np.array([[0, -1], [1, 0]])
    for index, dim in enumerate(dim_array):
        if np.any(z_object.z[index] == 0.0 + 0.0j) == True:
            dis[index] = np.identity(2)
            continue

        if dim == 1:

            if comp in ["01", "10"]:
                compr = np.abs(
                    z_object.z.real[index, int(comp[0]), int(comp[1])]
                )
                compi = np.abs(
                    z_object.z.imag[index, int(comp[0]), int(comp[1])]
                )
            else:
                compr = np.sqrt(np.linalg.det(z_object.z.real[index]))
                compi = np.sqrt(np.linalg.det(z_object.z.imag[index]))

            dis[index] = np.mean(
                np.array(
                    [
                        (1.0 / compr * np.dot(z_object.z.real[index], rot_mat)),
                        (1.0 / compi * np.dot(z_object.z.imag[index], rot_mat)),
                    ]
                ),
                axis=0,
            )

            if z_object.z_error is not None:
                # find errors of entries for calculatincomp weights

                compr_error = 1.0 / compr * np.abs(z_object.z_error[index])
                compr_error[np.where(compr_error == 0.0)] = 1.0

                compi_error = 1.0 / compi * np.abs(z_object.z_error[index])
                compi_error[np.where(compi_error == 0.0)] = 1.0

                dis_error[index] = np.mean(
                    np.array([compi_error, compr_error]), axis=0
                )

        elif dim == 2:
            P = 1
            strike_ang = st_array[index]
            if np.isnan(strike_ang):
                strike_ang = 0.0

            if z_object.z_error is not None:
                err_array = z_object.z_error[index]
                err_array[np.where(err_array == 0.0)] = 1.0
            else:
                err_array = None

            tetm_array, tetm_error = MTcc.rotate_matrix_with_errors(
                z_object.z[index], strike_ang, error=err_array
            )

            tetm_r = tetm_array.real
            tetm_i = tetm_array.imag
            t_array_r = (
                -4 * P * tetm_r[0, 1] * tetm_r[1, 0] / np.linalg.det(tetm_r)
            )
            t_array_i = (
                -4 * P * tetm_i[0, 1] * tetm_i[1, 0] / np.linalg.det(tetm_i)
            )

            try:
                T = np.sqrt(max([t_array_r, t_array_i])) + 0.001
            except ValueError:
                T = 2

            sr = np.sqrt(
                T**2
                + 4 * P * tetm_r[0, 1] * tetm_r[1, 0] / np.linalg.det(tetm_r)
            )
            si = np.sqrt(
                T**2
                + 4 * P * tetm_i[0, 1] * tetm_i[1, 0] / np.linalg.det(tetm_i)
            )

            par_r = 2 * tetm_r[0, 1] / (T - sr)
            orth_r = 2 * tetm_r[1, 0] / (T + sr)
            par_i = 2 * tetm_i[0, 1] / (T - si)
            orth_i = 2 * tetm_i[1, 0] / (T + si)

            mat2_r = np.array([[0, 1.0 / orth_r], [1.0 / par_r, 0]])
            mat2_i = np.array([[0, 1.0 / orth_i], [1.0 / par_i, 0]])

            avg_mat = np.mean(
                np.array([np.dot(tetm_r, mat2_r), np.dot(tetm_i, mat2_i)]),
                axis=0,
            )

            dis[index] = avg_mat

            if err_array is not None:
                # find errors of entries for calculatincomp weights
                sigma_sr = np.sqrt(
                    (
                        -(
                            2
                            * P
                            * tetm_r[0, 1]
                            * tetm_r[1, 0]
                            * tetm_r[1, 1]
                            * err_array[0, 0]
                        )
                        / (np.linalg.det(tetm_r) ** 2 * sr)
                    )
                    ** 2
                    + (
                        (
                            2
                            * P
                            * tetm_r[0, 0]
                            * tetm_r[1, 0]
                            * tetm_r[1, 1]
                            * err_array[0, 1]
                        )
                        / (np.linalg.det(tetm_r) ** 2 * sr)
                    )
                    ** 2
                    + (
                        (
                            2
                            * P
                            * tetm_r[0, 0]
                            * tetm_r[0, 1]
                            * tetm_r[1, 1]
                            * err_array[1, 0]
                        )
                        / (np.linalg.det(tetm_r) ** 2 * sr)
                    )
                    ** 2
                    + (
                        -(
                            2
                            * P
                            * tetm_r[0, 1]
                            * tetm_r[1, 0]
                            * tetm_r[0, 0]
                            * err_array[1, 1]
                        )
                        / (np.linalg.det(tetm_r) ** 2 * sr)
                    )
                    ** 2
                )

                sigma_dr_11 = 0.5 * sigma_sr
                sigma_dr_22 = 0.5 * sigma_sr

                sigma_dr_12 = np.sqrt(
                    (mat2_r[0, 1] / tetm_r[0, 0] * err_array[0, 0]) ** 2
                    + (mat2_r[0, 1] / tetm_r[1, 0] * err_array[1, 0]) ** 2
                    + (0.5 * tetm_r[0, 0] / tetm_r[1, 0] * sigma_sr) ** 2
                )
                sigma_dr_21 = np.sqrt(
                    (mat2_r[1, 0] / tetm_r[1, 1] * err_array[1, 1]) ** 2
                    + (mat2_r[1, 0] / tetm_r[0, 1] * err_array[0, 1]) ** 2
                    + (0.5 * tetm_r[1, 1] / tetm_r[0, 1] * sigma_sr) ** 2
                )

                dis_error_r = np.array(
                    [
                        [sigma_dr_11, sigma_dr_12],
                        [sigma_dr_21, sigma_dr_22],
                    ]
                )

                sigma_si = np.sqrt(
                    (
                        -(
                            2
                            * P
                            * tetm_i[0, 1]
                            * tetm_i[1, 0]
                            * tetm_i[1, 1]
                            * err_array[0, 0]
                        )
                        / (np.linalg.det(tetm_i) ** 2 * sr)
                    )
                    ** 2
                    + (
                        (
                            2
                            * P
                            * tetm_i[0, 0]
                            * tetm_i[1, 0]
                            * tetm_i[1, 1]
                            * err_array[0, 1]
                        )
                        / (np.linalg.det(tetm_i) ** 2 * sr)
                    )
                    ** 2
                    + (
                        (
                            2
                            * P
                            * tetm_i[0, 0]
                            * tetm_i[0, 1]
                            * tetm_i[1, 1]
                            * err_array[1, 0]
                        )
                        / (np.linalg.det(tetm_i) ** 2 * sr)
                    )
                    ** 2
                    + (
                        -(
                            2
                            * P
                            * tetm_i[0, 1]
                            * tetm_i[1, 0]
                            * tetm_i[0, 0]
                            * err_array[1, 1]
                        )
                        / (np.linalg.det(tetm_i) ** 2 * sr)
                    )
                    ** 2
                )

                sigma_di_11 = 0.5 * sigma_si
                sigma_di_22 = 0.5 * sigma_si
                sigma_di_12 = np.sqrt(
                    (mat2_i[0, 1] / tetm_i[0, 0] * err_array[0, 0]) ** 2
                    + (mat2_i[0, 1] / tetm_i[1, 0] * err_array[1, 0]) ** 2
                    + (0.5 * tetm_i[0, 0] / tetm_i[1, 0] * sigma_si) ** 2
                )
                sigma_di_21 = np.sqrt(
                    (mat2_i[1, 0] / tetm_i[1, 1] * err_array[1, 1]) ** 2
                    + (mat2_i[1, 0] / tetm_i[0, 1] * err_array[0, 1]) ** 2
                    + (0.5 * tetm_i[1, 1] / tetm_i[0, 1] * sigma_si) ** 2
                )

                dis_error_i = np.array(
                    [
                        [sigma_di_11, sigma_di_12],
                        [sigma_di_21, sigma_di_22],
                    ]
                )

                dis_error[index] = np.mean(np.array([dis_error_r, dis_error_i]))
        else:
            dis[index] = np.identity(2)

    nonzero_index = np.array(list(set(np.nonzero(dis)[0])))

    dis_avg, weights_sum = np.average(
        dis[nonzero_index],
        axis=0,
        weights=(1.0 / dis_error[nonzero_index]) ** 2,
        returned=True,
    )

    dis_avg_error = np.sqrt(1.0 / weights_sum)

    return dis_avg, dis_avg_error


def remove_distortion_from_z_object(
    z_object, distortion_tensor, distortion_error_tensor=None, logger=None
):
    """
    Remove distortion D form an observed impedance tensor Z to obtain
    the uperturbed "correct" Z0:
    Z = D * Z0

    Propagation of errors/uncertainties included


    :param distortion_tensor: real distortion tensor as a 2x2
    :type distortion_tensor: np.ndarray(2, 2, dtype=real)
    :param distortion_error_tensor: default is None
    :type distortion_error_tensor: np.ndarray(2, 2, dtype=real),
    :param inplace: Update the current object or return a new impedance
    :type inplace: boolean
    :returns: input distortion tensor
    :rtype: np.ndarray(2, 2, dtype='real')
    :returns: impedance tensor with distorion removed
    :rtype: np.ndarray(num_frequency, 2, 2, dtype='complex')
    :returns: impedance tensor error after distortion is removed
    :rtype: np.ndarray(num_frequency, 2, 2, dtype='complex')

    :Example: ::

            >>> distortion = np.array([[1.2, .5],[.35, 2.1]])
            >>> d, new_z, new_z_error = z_obj.remove_distortion(distortion)

    """

    if distortion_error_tensor is None:
        distortion_error_tensor = np.zeros_like(distortion_tensor)
    # for all frequency, calculate D.Inverse, then obtain Z0 = D.I * Z
    try:
        if not (len(distortion_tensor.shape) in [2, 3]) and (
            len(distortion_error_tensor.shape) in [2, 3]
        ):
            msg = "Distortion tensor and error are not correct shape"
            if logger is not None:
                logger.error(msg)
            raise ValueError(msg)
        if (
            len(distortion_tensor.shape) == 3
            or len(distortion_error_tensor.shape) == 3
        ):
            if logger is not None:
                logger.info(
                    "Distortion is not time-dependent - taking only first"
                    " of given distortion tensors."
                )
            try:
                distortion_tensor = distortion_tensor[0]
                distortion_error_tensor = distortion_error_tensor[0]
            except IndexError:
                msg = "Distortion tensor and error are not correct shape"
                if logger is not None:
                    logger.error(msg)
                raise ValueError(msg)
        if (distortion_tensor.shape != (2, 2)) and (
            distortion_error_tensor.shape != (2, 2)
        ):
            msg = "Distortion tensor and error are not correct shape"
            if logger is not None:
                logger.error(msg)
            raise ValueError(msg)
        distortion_tensor = np.array(np.real(distortion_tensor))
    except ValueError:
        msg = "Input distortion tensor, must be (2, 2)"
        raise ValueError(msg)
    try:
        DI = np.linalg.inv(distortion_tensor)
    except np.linalg.LinAlgError:
        raise ValueError(
            "The provided distortion tensor is singular cannot be used."
        )
    # propagation of errors (using 1-norm) - step 1 - inversion of D:
    DI_error = np.zeros_like(distortion_error_tensor)

    # todo :include error on  determinant!!
    # D_det = np.linalg.det(distortion_tensor)

    dummy, DI_error = MTcc.invertmatrix_incl_errors(
        distortion_tensor, distortion_error_tensor
    )

    # propagation of errors - step 2 - product of D.inverse and Z;
    # D.I * Z, making it 4 summands for each component:
    z_corrected = np.zeros_like(z_object.z, dtype=complex)
    z_corrected_error = np.zeros_like(z_object.z, dtype=float)

    for idx_f in range(len(z_object.z)):
        z_corrected[idx_f] = np.array(np.dot(DI, np.array(z_object.z[idx_f])))
        if z_object._has_tf_error():
            for ii in range(2):
                for jj in range(2):
                    z_corrected_error[idx_f, ii, jj] = np.sum(
                        np.abs(
                            np.array(
                                [
                                    DI_error[ii, 0] * z_object.z[idx_f, 0, jj],
                                    DI[ii, 0] * z_object.z_error[idx_f, 0, jj],
                                    DI_error[ii, 1] * z_object.z[idx_f, 1, jj],
                                    DI[ii, 1] * z_object.z_error[idx_f, 1, jj],
                                ]
                            )
                        )
                    )

    return z_corrected, z_corrected_error
