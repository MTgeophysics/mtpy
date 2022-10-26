#!/usr/bin/env python

"""
.. module:: Z
   :synopsis: Deal with MT responses Z and Tipper

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov> 
.. moduleauthor:: Lars Krieger

Updated 11/2020 for logging and formating (J. Peacock).
    - ToDo: add functionality for covariance matrix
"""

# =============================================================================
# Imports
# =============================================================================
import copy

import numpy as np

import mtpy.utils.calculator as MTcc
from mtpy.utils.exceptions import (
    MTpyError_Z,
)

from .res_phase import ResPhase
from .base import TFBase

# ==============================================================================
# Impedance Tensor Class
# ==============================================================================
class Z(TFBase, ResPhase):
    """
    Z class - generates an impedance tensor (Z) object.

    Z is a complex array of the form (n_frequency, 2, 2),
    with indices in the following order:

        - Zxx: (0,0)
        - Zxy: (0,1)
        - Zyx: (1,0)
        - Zyy: (1,1)

    All errors are given as standard deviations (sqrt(VAR))

    :param z: array containing complex impedance values
    :type z: numpy.ndarray(n_frequency, 2, 2)


    :param z_error: array containing error values (standard deviation)
                        of impedance tensor elements
    :type z_error: numpy.ndarray(n_frequency, 2, 2)

    :param frequency: array of frequencyuency values corresponding to impedance tensor
                 elements.
    :type frequency: np.ndarray(n_frequency)

    :Example: ::

        >>> import mtpy.core.z as mtz
        >>> import numpy as np
        >>> z_test = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
        >>> z_object = mtz.Z(z=z_test, frequency=[1])
        >>> z_object.rotate(45)
        >>> z_object.resistivity


    """

    def __init__(
        self,
        z=None,
        z_error=None,
        frequency=None,
        z_model_error=None,
    ):
        """
        Initialise an instance of the Z class.

        :param z: array containing complex impedance values
        :type z: numpy.ndarray(n_frequency, 2, 2)

        :param z_error: array containing error values (standard deviation)
                            of impedance tensor elements
        :type z_error: numpy.ndarray(n_frequency, 2, 2)

        :param frequency: array of frequencyuency values corresponding to impedance
                     tensor elements.
        :type frequency: np.ndarray(n_frequency)

        Initialises the attributes with None

        """

        TFBase.__init__(
            self,
            tf=z,
            tf_error=z_error,
            tf_model_error=z_model_error,
            frequency=frequency,
            _name="impedance",
        )

    @property
    def z(self):
        """
        Impedance tensor

        np.ndarray(nfrequency, 2, 2)
        """
        if self._has_tf():
            return self._dataset.transfer_function.values

    @z.setter
    def z(self, z):
        """
        Set the attribute 'z'.


        :param z: complex impedance tensor array
        :type z: np.ndarray(nfrequency, 2, 2)

        Test for shape, but no test for consistency!

        Nulling the rotation_angle
        """

        old_shape = None
        if self._has_tf():
            old_shape = self._dataset.transfer_function.shape
        z = self._validate_array_input(z, "complex", old_shape)

        comps = dict(input=self.inputs, output=self.outputs)
        if self._is_empty():
            self._dataset = self._initialize(tf=z)
        else:
            self._dataset["transfer_function"].loc[comps] = z

    # ----impedance error-----------------------------------------------------
    @property
    def z_error(self):
        if self._has_tf_error():
            return self._dataset.transfer_function_error.values

    @z_error.setter
    def z_error(self, z_error):
        """
        Set the attribute z_error

        :param z_error: error of impedance tensor array as standard
                            deviation
        :type z_error: np.ndarray(nfrequency, 2, 2)
        """
        old_shape = None
        if not self._has_tf_error():
            old_shape = self._dataset.transfer_function_error.shape
        self._dataset["transfer_function_error"] = self._validate_array_input(
            z_error, "float", old_shape
        )

    # ----impedance model error-----------------------------------------------------
    @property
    def z_model_error(self):
        if self._has_tf_model_error():
            return self._dataset.transfer_function_model_error.values

    @z_model_error.setter
    def z_model_error(self, z_model_error):
        """
        Set the attribute z_model_error

        :param z_model_error: error of impedance tensor array as standard
         deviation
        :type z_model_error: np.ndarray(nfrequency, 2, 2)
        """

        old_shape = None
        if not self._has_tf_error():
            old_shape = self._dataset.transfer_function_error.shape
        self._dataset[
            "transfer_function_model_error"
        ] = self._validate_array_input(z_model_error, "float", old_shape)

    def remove_ss(
        self, reduce_res_factor_x=1.0, reduce_res_factor_y=1.0, inplace=False
    ):
        """
        Remove the static shift by providing the respective correction factors
        for the resistivity in the x and y components.
        (Factors can be determined by using the "Analysis" module for the
        impedance tensor)

        Assume the original observed tensor Z is built by a static shift S
        and an unperturbated "correct" Z0 :

             * Z = S * Z0

        therefore the correct Z will be :
            * Z0 = S^(-1) * Z


        :param reduce_res_factor_x: static shift factor to be applied to x
                                    components (ie z[:, 0, :]).  This is
                                    assumed to be in resistivity scale
        :type reduce_res_factor_x: float or iterable list or array

        :param reduce_res_factor_y: static shift factor to be applied to y
                                    components (ie z[:, 1, :]).  This is
                                    assumed to be in resistivity scale
        :type reduce_res_factor_y: float or iterable list or array

        :returns: static shift matrix,
        :rtype: np.ndarray ((2, 2))

        :returns: corrected Z
        :rtype: mtpy.core.z.Z

        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots!

        """

        # check for iterable list/set of reduce_res_factor_x - if so, it must
        # have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_x) == 0:
            try:
                x_factor = float(reduce_res_factor_x)
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
            lo_x_factors = np.repeat(x_factor, len(self.z))
        elif len(reduce_res_factor_x) == 1:
            try:
                x_factor = float(reduce_res_factor_x)
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
            lo_x_factors = np.repeat(x_factor, len(self.z))
        else:
            try:
                lo_x_factors = np.repeat(x_factor, len(reduce_res_factor_x))
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
        if len(lo_x_factors) != len(self.z):
            msg = (
                f"Length of reduce_res_factor_x needs to be {len(self.z)}"
                + f" not {len(lo_x_factors)}"
            )
            self.logger.error(msg)
            raise ValueError(msg)
        # check for iterable list/set of reduce_res_factor_y - if so,
        # it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_y) == 0:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                msg = "reduce_res_factor_y must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
            lo_y_factors = np.repeat(y_factor, len(self.z))
        elif len(reduce_res_factor_y) == 1:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                msg = "reduce_res_factor_y must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
            lo_y_factors = np.repeat(y_factor, len(self.z))
        else:
            try:
                lo_y_factors = np.repeat(y_factor, len(reduce_res_factor_y))
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self.logger.error(msg)
                raise ValueError(msg)
        if len(lo_y_factors) != len(self.z):
            msg = (
                f"Length of reduce_res_factor_x needs to be {len(self.z)}"
                + f" not {len(lo_x_factors)}"
            )
            self.logger.error(msg)
            raise ValueError(msg)
        z_corrected = copy.copy(self.z)
        static_shift = np.zeros((len(self.z), 2, 2))

        for idx_f in range(len(self.z)):
            # correct for x-direction
            z_corrected[idx_f, 0, :] = self.z[idx_f, 0, :] / np.sqrt(
                lo_x_factors[idx_f]
            )
            # correct for y-direction
            z_corrected[idx_f, 1, :] = self.z[idx_f, 1, :] / np.sqrt(
                lo_y_factors[idx_f]
            )
            # make static shift array
            static_shift[idx_f, 0, 0] = np.sqrt(lo_x_factors[idx_f])
            static_shift[idx_f, 1, 1] = np.sqrt(lo_y_factors[idx_f])
        return static_shift, z_corrected

    def remove_distortion(self, distortion_tensor, distortion_err_tensor=None):
        """
        Remove distortion D form an observed impedance tensor Z to obtain
        the uperturbed "correct" Z0:
        Z = D * Z0

        Propagation of errors/uncertainties included


        :param distortion_tensor: real distortion tensor as a 2x2
        :type distortion_tensor: np.ndarray(2, 2, dtype=real)


        :param distortion_err_tensor: default is None
        :type distortion_err_tensor: np.ndarray(2, 2, dtype=real),

                :returns: input distortion tensor
        :rtype: np.ndarray(2, 2, dtype='real')

                :returns: impedance tensor with distorion removed
        :rtype: np.ndarray(num_frequency, 2, 2, dtype='complex')


                :returns: impedance tensor error after distortion is removed
        :rtype: np.ndarray(num_frequency, 2, 2, dtype='complex')


                :Example: ::

                        >>> import mtpy.core.z as mtz
                        >>> distortion = np.array([[1.2, .5],[.35, 2.1]])
                        >>> d, new_z, new_z_error = z_obj.remove_distortion(distortion)

        """

        if distortion_err_tensor is None:
            distortion_err_tensor = np.zeros_like(distortion_tensor)
        # for all frequency, calculate D.Inverse, then obtain Z0 = D.I * Z
        try:
            if not (len(distortion_tensor.shape) in [2, 3]) and (
                len(distortion_err_tensor.shape) in [2, 3]
            ):
                msg = "Distortion tensor and error are not correct shape"
                self.logger.error(msg)
                raise ValueError(msg)
            if (
                len(distortion_tensor.shape) == 3
                or len(distortion_err_tensor.shape) == 3
            ):
                self.logger.info(
                    "Distortion is not time-dependent - taking only first"
                    + "of given distortion tensors"
                )
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_err_tensor = distortion_err_tensor[0]
                except IndexError:
                    msg = "Distortion tensor and error are not correct shape"
                    self.logger.error(msg)
                    raise ValueError(msg)
            if not (distortion_tensor.shape == (2, 2)) and (
                distortion_err_tensor.shape == (2, 2)
            ):
                msg = "Distortion tensor and error are not correct shape"
                self.logger.error(msg)
                raise ValueError(msg)
            distortion_tensor = np.matrix(np.real(distortion_tensor))
        except ValueError:
            msg = "Input distortion tensor, must be (2, 2)"
            raise MTpyError_Z(msg)
        try:
            DI = distortion_tensor.I
        except np.linalg.LinAlgError:
            raise MTpyError_Z(
                "The provided distortion tensor is singular cannot be used."
            )
        # propagation of errors (using 1-norm) - step 1 - inversion of D:
        DI_err = np.zeros_like(distortion_err_tensor)

        # todo :include error on  determinant!!
        # D_det = np.linalg.det(distortion_tensor)

        dummy, DI_err = MTcc.invertmatrix_incl_errors(
            distortion_tensor, distortion_err_tensor
        )

        # propagation of errors - step 2 - product of D.inverse and Z;
        # D.I * Z, making it 4 summands for each component:
        z_corrected = np.zeros_like(self.z)
        z_corrected_err = np.zeros_like(self.z_error)

        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot(DI, np.matrix(self.z[idx_f])))
            for ii in range(2):
                for jj in range(2):
                    z_corrected_err[idx_f, ii, jj] = np.sum(
                        np.abs(
                            np.array(
                                [
                                    DI_err[ii, 0] * self.z[idx_f, 0, jj],
                                    DI[ii, 0] * self.z_error[idx_f, 0, jj],
                                    DI_err[ii, 1] * self.z[idx_f, 1, jj],
                                    DI[ii, 1] * self.z_error[idx_f, 1, jj],
                                ]
                            )
                        )
                    )
        return distortion_tensor, z_corrected, z_corrected_err

    @property
    def only_1d(self):
        """
        Return Z in 1D form.

        If Z is not 1D per se, the diagonal elements are set to zero,
        the off-diagonal elements keep their signs, but their absolute
        is set to the mean of the original Z off-diagonal absolutes.
        """

        z1d = copy.copy(self.z)

        for ii in range(len(z1d)):
            z1d[ii, 0, 0] = 0
            z1d[ii, 1, 1] = 0
            sign01 = np.sign(z1d[ii, 0, 1])
            sign10 = np.sign(z1d[ii, 1, 0])
            mean1d = 0.5 * (z1d[ii, 1, 0] + z1d[ii, 0, 1])
            z1d[ii, 0, 1] = sign01 * mean1d
            z1d[ii, 1, 0] = sign10 * mean1d
        return z1d

    @property
    def only_2d(self):
        """
        Return Z in 2D form.

        If Z is not 2D per se, the diagonal elements are set to zero.
        """
        if self.z is not None:
            z2d = np.zeros_like(self.z, dtype=complex)
            z2d[:, 0, 1] = self.z[:, 0, 1]
            z2d[:, 1, 0] = self.z[:, 1, 0]
            return z2d

    @property
    def trace(self):
        """
        Return the trace of Z

        :returns: Trace(z)
        :rtype: np.ndarray(nfrequency, 2, 2)

        """

        if self.z is not None:
            tr = np.array([np.trace(ii) for ii in self.z])

            return tr

    @property
    def trace_err(self):
        """
        Return the trace of Z

        :returns: Trace(z)
        :rtype: np.ndarray(nfrequency, 2, 2)

        """

        tr_err = None
        if self.z_error is not None:
            tr_err = np.zeros_like(self.trace, dtype=np.float)
            tr_err[:] = self.z_error[:, 0, 0] + self.z_error[:, 1, 1]
        return tr_err

    @property
    def skew(self):
        """
        Returns the skew of Z as defined by Z[0, 1] + Z[1, 0]

        .. note:: This is not the MT skew, but simply the linear algebra skew


        :returns: skew
        :rtype: np.ndarray(nfrequency, 2, 2)
        """

        if self.z is not None:
            skew = np.array([ii[0, 1] - ii[1, 0] for ii in self.z])

            return skew

    @property
    def skew_err(self):
        """
        Returns the skew error of Z as defined by z_error[0, 1] + z_error[1, 0]

        .. note:: This is not the MT skew, but simply the linear algebra skew

        :returns: skew_err
        :rtype: np.ndarray(nfrequency, 2, 2)
        """

        skew_err = None
        if self.z_error is not None:
            skew_err = np.zeros_like(self.skew, dtype=np.float)
            skew_err[:] = self.z_error[:, 0, 1] + self.z_error[:, 1, 0]
        return skew_err

    @property
    def det(self):
        """
        Return the determinant of Z

        :returns: det_Z
        :rtype: np.ndarray(nfrequency)
        """
        if self.z is not None:
            det_Z = np.array([np.linalg.det(ii) for ii in self.z])

            return det_Z

    @property
    def det_err(self):
        """
        Return the determinant of Z error

        :returns: det_z_error
        :rtype: np.ndarray(nfrequency)
        """
        det_z_error = None
        if self.z_error is not None:
            det_z_error = np.zeros_like(self.det, dtype=np.float)
            # components of the impedance tensor are not independent variables
            # so can't use standard error propagation
            # calculate manually:
            # difference of determinant of z + z_error and z - z_error then divide by 2
            det_z_error[:] = (
                np.abs(
                    np.linalg.det(self.z + self.z_error)
                    - np.linalg.det(self.z - self.z_error)
                )
                / 2.0
            )
        return det_z_error

    @property
    def norm(self):
        """
        Return the 2-/Frobenius-norm of Z

        :returns: norm
        :rtype: np.ndarray(nfrequency)
        """

        if self.z is not None:
            norm = np.array([np.linalg.norm(ii) for ii in self.z])

            return norm

    @property
    def norm_err(self):
        """
        Return the 2-/Frobenius-norm of Z  error

        :returns: norm_err
        :rtype: np.ndarray(nfrequency)
        """
        norm_err = None

        if self.z_error is not None:
            norm_err = np.zeros_like(self.norm, dtype=np.float)
            for idx, z_tmp in enumerate(self.z):
                value = self.norm[idx]
                error_matrix = self.z_error[idx]
                radicand = 0.0
                for ii in range(2):
                    for jj in range(2):
                        radicand += (
                            error_matrix[ii, jj] * np.real(z_tmp[ii, jj])
                        ) ** 2
                        radicand += (
                            error_matrix[ii, jj] * np.imag(z_tmp[ii, jj])
                        ) ** 2
                norm_err[idx] = 1.0 / value * np.sqrt(radicand)
        return norm_err

    @property
    def invariants(self):
        """
        Return a dictionary of Z-invariants.

        Contains
                -----------
                        * z1
                        * det
                        * det_real
                        * det_imag
                        * trace
                        * skew
                        * norm
                        * lambda_plus/minus,
                        * sigma_plus/minus
        """

        invariants_dict = {}

        if self.z is not None:
            z1 = (self.z[:, 0, 1] - self.z[:, 1, 0]) / 2.0
            invariants_dict["z1"] = z1

            invariants_dict["det"] = self.det[0]

            det_real = np.array([np.linalg.det(ii) for ii in np.real(self.z)])
            invariants_dict["det_real"] = det_real

            det_imag = np.array([np.linalg.det(ii) for ii in np.imag(self.z)])
            invariants_dict["det_imag"] = det_imag

            invariants_dict["trace"] = self.trace

            invariants_dict["skew"] = self.skew

            invariants_dict["norm"] = self.norm

            invariants_dict["lambda_plus"] = z1 + np.sqrt(z1 * z1 / self.det)

            invariants_dict["lambda_minus"] = z1 - np.sqrt(z1 * z1 / self.det)

            invariants_dict["sigma_plus"] = (
                0.5 * self.norm**2
                + np.sqrt(0.25 * self.norm**4)
                + np.abs(self.det**2)
            )

            invariants_dict["sigma_minus"] = (
                0.5 * self.norm**2
                - np.sqrt(0.25 * self.norm**4)
                + np.abs(self.det**2)
            )

        return invariants_dict
