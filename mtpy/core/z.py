#!/usr/bin/env python

"""
.. module:: Z
   :synopsis: Deal with MT responses Z and Tipper

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov> 
.. moduleauthor:: Lars Krieger

Updated 11/2020 for logging and formating (J. Peacock).
    - ToDo: add functionality for covariance matrix
"""

# =================================================================
import cmath
import copy
import math

import numpy as np

import mtpy.utils.calculator as MTcc
from mtpy.utils.exceptions import (
    MTpyError_Z,
    MTpyError_Tipper,
    MTpyError_input_arguments,
)
from mtpy.utils.mtpy_logger import get_mtpy_logger
from .res_phase import ResPhase


# ==============================================================================
# Impedance Tensor Class
# ==============================================================================
class Z(ResPhase):
    """
    Z class - generates an impedance tensor (Z) object.

    Z is a complex array of the form (n_frequency, 2, 2),
    with indices in the following order:

        - Zxx: (0,0)
        - Zxy: (0,1)
        - Zyx: (1,0)
        - Zyy: (1,1)

    All errors are given as standard deviations (sqrt(VAR))

    :param z_array: array containing complex impedance values
    :type z_array: numpy.ndarray(n_frequency, 2, 2)


    :param z_err_array: array containing error values (standard deviation)
                        of impedance tensor elements
    :type z_err_array: numpy.ndarray(n_frequency, 2, 2)

    :param frequency: array of frequencyuency values corresponding to impedance tensor
                 elements.
    :type frequency: np.ndarray(n_frequency)

    :Example: ::

        >>> import mtpy.core.z as mtz
        >>> import numpy as np
        >>> z_test = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
        >>> z_object = mtz.Z(z_array=z_test, frequency=[1])
        >>> z_object.rotate(45)
        >>> z_object.resistivity


    """

    def __init__(
        self,
        z_array=None,
        z_err_array=None,
        frequency=None,
        z_model_err_array=None,
    ):
        """
        Initialise an instance of the Z class.

        :param z_array: array containing complex impedance values
        :type z_array: numpy.ndarray(n_frequency, 2, 2)

        :param z_err_array: array containing error values (standard deviation)
                            of impedance tensor elements
        :type z_err_array: numpy.ndarray(n_frequency, 2, 2)

        :param frequency: array of frequencyuency values corresponding to impedance
                     tensor elements.
        :type frequency: np.ndarray(n_frequency)

        Initialises the attributes with None

        """
        self.rotation_angle = 0.0

        super().__init__()
        self.z = z_array
        self.z_err = z_err_array
        self.z_model_err = z_model_err_array
        self.frequency = frequency

        if self.z is not None:
            self.rotation_angle = np.zeros((len(self.z)))
        if self.z is not None and self.frequency is not None:
            self.compute_resistivity_phase()

    def __str__(self):
        lines = ["Impedance Tensor", "-" * 30]
        if self.frequency is not None:
            lines.append(
                f"\tNumber of frequencyuencies:  {self.frequency.size}"
            )
            lines.append(
                f"\tfrequencyuency range:        {self.frequency.min():.5E} -- {self.frequency.max():.5E} Hz"
            )
            lines.append(
                f"\tPeriod range:           {1/self.frequency.max():.5E} -- {1/self.frequency.min():.5E} s"
            )
            lines.append("")
            lines.append("\tElements:")
            for zz, ff in zip(self.z, self.frequency):
                lines.append(f"\tfrequencyuency: {ff:5E} Hz -- Period {1/ff} s")
                lines.append(
                    "\t\t"
                    + np.array2string(
                        zz,
                        formatter={
                            "complex_kind": lambda x: f"{x.real:.4e}{x.imag:+.4E}"
                        },
                    ).replace("\n", "\n\t\t")
                )
        else:
            if self.z is not None:
                lines.append("Elements:")
                for ff, zz in enumerate(self.z):
                    lines.append(f"\tIndex {ff}")
                    lines.append(
                        "\t\t"
                        + np.array2string(
                            zz,
                            formatter={
                                "complex_kind": lambda x: f"{x.real:.4e}{x.imag:+.4E}"
                            },
                        ).replace("\n", "\n\t\t")
                    )
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, Z):
            msg = f"Cannot compare {type(other)} with Z"
            self._logger.error(msg)
            raise MTpyError_Z(msg)
        if (self.z != other.z).all():
            return False
        if (self.frequency != other.frequency).all():
            return False
        if (self.z_err != other.z_err).all():
            return False
        return True

    def copy(self):
        return copy.deepcopy(self)

    # ---frequencyuency-------------------------------------------------------------
    @property
    def frequency(self):
        """
        frequencyuencies for each impedance tensor element

        Units are Hz.
        """
        return self._frequency

    @frequency.setter
    def frequency(self, frequency_arr):
        """
        Set the array of frequency.

        :param frequency_arr: array of frequencyunecies (Hz)
        :type frequency_arr: np.ndarray
        """

        if frequency_arr is None:
            self._frequency = None
            return

        self._frequency = np.array(frequency_arr, dtype="float")

        if self.z is not None:
            if self.z.shape[0] != len(self._frequency):
                msg = (
                    "New frequency array is not correct shape for existing z. "
                    + f"new: {self._frequency.size} != old: {self.z.shape[0]}"
                )
                self._logger.error(msg)
                raise MTpyError_Z(msg)
            self.compute_resistivity_phase()

    @property
    def period(self):
        """
        periods in seconds
        """

        return 1.0 / self.frequency

    @period.setter
    def period(self, value):
        """
        setting periods will set the frequencyuencies
        """

        self.frequency = 1.0 / value

    # ----impedance tensor -----------------------------------------------------
    def _validate_impedance_input(self, z_array, dtype, old_shape=None):
        """
        Validate an input impedance array

        :param array: DESCRIPTION
        :type array: TYPE
        :param dtype: DESCRIPTION
        :type dtype: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if z_array is None:
            return
        if not isinstance(z_array, np.ndarray):
            z_array = np.array(z_array, dtype=dtype)
        if z_array.dtype not in [dtype]:
            z_array = z_array.astype(dtype)
        # check to see if the new z array is the same shape as the old
        if old_shape is not None and old_shape != z_array.shape:
            msg = (
                "Shape of new array does not match old.  "
                + f"new shape {z_array.shape} != "
                + f"old shape {self._z.shape}. "
                + "Make a new Z instance to be safe."
            )
            self._logger.error(msg)
            raise MTpyError_Z(msg)

        if len(z_array.shape) == 3:
            if z_array.shape[1:3] == (2, 2):
                return z_array
            else:
                msg = f"Input array must be shape (n, 2, 2) not {z_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Z(msg)
        elif len(z_array.shape) == 2:
            if z_array.shape == (2, 2):
                return z_array.reshape((1, 2, 2))
                self._logger.debug(
                    "setting input z with shape (2, 2) to (1, 2, 2)"
                )
            else:
                msg = f"Input array must be shape (n, 2, 2) not {z_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Z(msg)
        else:
            msg = f"{z_array.shape} are not the correct dimensions, must be (n, 2, 2)"
            self._logger.error(msg)
            raise MTpyError_Z(msg)

    @property
    def z(self):
        """
        Impedance tensor

        np.ndarray(nfrequency, 2, 2)
        """
        return self._z

    @z.setter
    def z(self, z_array):
        """
        Set the attribute 'z'.


        :param z_array: complex impedance tensor array
        :type z_array: np.ndarray(nfrequency, 2, 2)

        Test for shape, but no test for consistency!

        Nulling the rotation_angle
        """
        old_shape = None
        if self._z is not None:
            old_shape = self._z.shape
        self._z = self._validate_impedance_input(z_array, "complex", old_shape)

        if self._z is not None:
            if isinstance(self.rotation_angle, float):
                self.rotation_angle = np.repeat(
                    self.rotation_angle, len(self._z)
                )

        # for consistency recalculate resistivity and phase
        try:
            self.compute_resistivity_phase()
        except ValueError as error:
            self.logger.debug(error)

    # ----impedance error-----------------------------------------------------
    @property
    def z_err(self):
        return self._z_err

    @z_err.setter
    def z_err(self, z_err_array):
        """
        Set the attribute z_err

        :param z_err_array: error of impedance tensor array as standard
                            deviation
        :type z_err_array: np.ndarray(nfrequency, 2, 2)
        """
        old_shape = None
        if self._z is not None:
            old_shape = self._z.shape

        self._z_err = self._validate_impedance_input(
            z_err_array, "float", old_shape
        )

        # for consistency recalculate resistivity and phase
        try:
            self.compute_resistivity_phase()
        except ValueError as error:
            self.logger.debug(error)

    # ----impedance model error-----------------------------------------------------
    @property
    def z_model_err(self):
        return self._z_model_err

    @z_model_err.setter
    def z_model_err(self, z_model_err_array):
        """
        Set the attribute z_model_err

        :param z_model_err_array: error of impedance tensor array as standard
                            deviation
        :type z_model_err_array: np.ndarray(nfrequency, 2, 2)
        """
        old_shape = None
        if self._z is not None:
            old_shape = self._z.shape

        self._z_model_err = self._validate_impedance_input(
            z_model_err_array, "float", old_shape
        )

        # for consistency recalculate resistivity and phase
        try:
            self.compute_resistivity_phase()
        except ValueError as error:
            self.logger.debug(error)

    @property
    def inverse(self):
        """
        Return the inverse of Z.

        (no error propagtaion included yet)

        """

        if self.z is None:
            self._logger.warn('z array is "None" - I cannot invert that')
            return
        inverse = copy.copy(self.z)
        for idx_f in range(len(inverse)):
            try:
                inverse[idx_f, :, :] = np.array(
                    (np.matrix(self.z[idx_f, :, :])).I
                )
            except:
                msg = f"The {idx_f + 1}ith impedance tensor cannot be inverted"
                raise MTpyError_Z(msg)
        return inverse

    def rotate(self, alpha):
        """
        Rotate Z array by angle alpha.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North, positive in clockwise direction.
        (Mathematically negative!)

        In non-rotated state, X refs to North and Y to East direction.

        Updates the attributes
            - *z*
            - *z_err*
            - *zrot*
            - *resistivity*
            - *phase*
            - *resistivity_err*
            - *phase_err*

        """

        if self.z is None:
            self._logger.warning('Z array is "None" and cannot be rotated')
            return
        # check for iterable list/set of angles - if so, it must have length
        # 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                msg = f"Angle must be a valid number (in degrees) not {alpha}"
                self._logger.error(msg)
                raise MTpyError_Z(msg)
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.z]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha % 360)
                except ValueError:
                    msg = (
                        f"Angle must be a valid number (in degrees) not {alpha}"
                    )
                    self._logger.error(msg)
                    raise MTpyError_Z(msg)
                # make an n long list of identical angles
                lo_angles = [degreeangle for ii in self.z]
            else:
                try:
                    lo_angles = [float(ii % 360) for ii in alpha]
                except ValueError:
                    msg = (
                        f"Angle must be a valid number (in degrees) not {alpha}"
                    )
                    self._logger.error(msg)
                    raise MTpyError_Z(msg)
        self.rotation_angle = np.array(
            [
                (oldangle + lo_angles[ii]) % 360
                for ii, oldangle in enumerate(self.rotation_angle)
            ]
        )

        if len(lo_angles) != len(self.z):
            msg = f"Wrong number of angles, need {len(self.z)}"
            self._logger.error(msg)
            raise MTpyError_Z(msg)
        z_rot = copy.copy(self.z)
        z_err_rot = copy.copy(self.z_err)

        for idx_frequency in range(len(self.z)):

            angle = lo_angles[idx_frequency]
            if np.isnan(angle):
                angle = 0.0
            if self.z_err is not None:
                (
                    z_rot[idx_frequency],
                    z_err_rot[idx_frequency],
                ) = MTcc.rotate_matrix_with_errors(
                    self.z[idx_frequency, :, :],
                    angle,
                    self.z_err[idx_frequency, :, :],
                )
            else:
                (
                    z_rot[idx_frequency],
                    z_err_rot,
                ) = MTcc.rotate_matrix_with_errors(
                    self.z[idx_frequency, :, :], angle
                )
        self.z = z_rot
        if self.z_err is not None:
            self.z_err = z_err_rot
        # for consistency recalculate resistivity and phase
        self.compute_resistivity_phase()

    def remove_ss(self, reduce_res_factor_x=1.0, reduce_res_factor_y=1.0):
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
                self._logger.error(msg)
                raise ValueError(msg)
            lo_x_factors = np.repeat(x_factor, len(self.z))
        elif len(reduce_res_factor_x) == 1:
            try:
                x_factor = float(reduce_res_factor_x)
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self._logger.error(msg)
                raise ValueError(msg)
            lo_x_factors = np.repeat(x_factor, len(self.z))
        else:
            try:
                lo_x_factors = np.repeat(x_factor, len(reduce_res_factor_x))
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self._logger.error(msg)
                raise ValueError(msg)
        if len(lo_x_factors) != len(self.z):
            msg = (
                f"Length of reduce_res_factor_x needs to be {len(self.z)}"
                + f" not {len(lo_x_factors)}"
            )
            self._logger.error(msg)
            raise ValueError(msg)
        # check for iterable list/set of reduce_res_factor_y - if so,
        # it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_y) == 0:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                msg = "reduce_res_factor_y must be a valid number"
                self._logger.error(msg)
                raise ValueError(msg)
            lo_y_factors = np.repeat(y_factor, len(self.z))
        elif len(reduce_res_factor_y) == 1:
            try:
                y_factor = float(reduce_res_factor_y)
            except ValueError:
                msg = "reduce_res_factor_y must be a valid number"
                self._logger.error(msg)
                raise ValueError(msg)
            lo_y_factors = np.repeat(y_factor, len(self.z))
        else:
            try:
                lo_y_factors = np.repeat(y_factor, len(reduce_res_factor_y))
            except ValueError:
                msg = "reduce_res_factor_x must be a valid number"
                self._logger.error(msg)
                raise ValueError(msg)
        if len(lo_y_factors) != len(self.z):
            msg = (
                f"Length of reduce_res_factor_x needs to be {len(self.z)}"
                + f" not {len(lo_x_factors)}"
            )
            self._logger.error(msg)
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
                        >>> d, new_z, new_z_err = z_obj.remove_distortion(distortion)

        """

        if distortion_err_tensor is None:
            distortion_err_tensor = np.zeros_like(distortion_tensor)
        # for all frequency, calculate D.Inverse, then obtain Z0 = D.I * Z
        try:
            if not (len(distortion_tensor.shape) in [2, 3]) and (
                len(distortion_err_tensor.shape) in [2, 3]
            ):
                msg = "Distortion tensor and error are not correct shape"
                self._logger.error(msg)
                raise ValueError(msg)
            if (
                len(distortion_tensor.shape) == 3
                or len(distortion_err_tensor.shape) == 3
            ):
                self._logger.info(
                    "Distortion is not time-dependent - taking only first"
                    + "of given distortion tensors"
                )
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_err_tensor = distortion_err_tensor[0]
                except IndexError:
                    msg = "Distortion tensor and error are not correct shape"
                    self._logger.error(msg)
                    raise ValueError(msg)
            if not (distortion_tensor.shape == (2, 2)) and (
                distortion_err_tensor.shape == (2, 2)
            ):
                msg = "Distortion tensor and error are not correct shape"
                self._logger.error(msg)
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
        z_corrected_err = np.zeros_like(self.z_err)

        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot(DI, np.matrix(self.z[idx_f])))
            for ii in range(2):
                for jj in range(2):
                    z_corrected_err[idx_f, ii, jj] = np.sum(
                        np.abs(
                            np.array(
                                [
                                    DI_err[ii, 0] * self.z[idx_f, 0, jj],
                                    DI[ii, 0] * self.z_err[idx_f, 0, jj],
                                    DI_err[ii, 1] * self.z[idx_f, 1, jj],
                                    DI[ii, 1] * self.z_err[idx_f, 1, jj],
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
        if self.z_err is not None:
            tr_err = np.zeros_like(self.trace, dtype=np.float)
            tr_err[:] = self.z_err[:, 0, 0] + self.z_err[:, 1, 1]
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
        Returns the skew error of Z as defined by Z_err[0, 1] + Z_err[1, 0]

        .. note:: This is not the MT skew, but simply the linear algebra skew

        :returns: skew_err
        :rtype: np.ndarray(nfrequency, 2, 2)
        """

        skew_err = None
        if self.z_err is not None:
            skew_err = np.zeros_like(self.skew, dtype=np.float)
            skew_err[:] = self.z_err[:, 0, 1] + self.z_err[:, 1, 0]
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

        :returns: det_Z_err
        :rtype: np.ndarray(nfrequency)
        """
        det_Z_err = None
        if self.z_err is not None:
            det_Z_err = np.zeros_like(self.det, dtype=np.float)
            # components of the impedance tensor are not independent variables
            # so can't use standard error propagation
            # calculate manually:
            # difference of determinant of z + z_err and z - z_err then divide by 2
            det_Z_err[:] = (
                np.abs(
                    np.linalg.det(self.z + self.z_err)
                    - np.linalg.det(self.z - self.z_err)
                )
                / 2.0
            )
        return det_Z_err

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

        if self.z_err is not None:
            norm_err = np.zeros_like(self.norm, dtype=np.float)
            for idx, z_tmp in enumerate(self.z):
                value = self.norm[idx]
                error_matrix = self.z_err[idx]
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


# ======================================================================
#                               TIPPER
# ======================================================================
class Tipper(object):
    """
    Tipper class --> generates a Tipper-object.

    Errors are given as standard deviations (sqrt(VAR))

    :param tipper_array: tipper array in the shape of [Tx, Ty]
                         *default* is None
    :type tipper_array: np.ndarray((nf, 1, 2), dtype='complex')


    :param tipper_err_array: array of estimated tipper errors
                             in the shape of [Tx, Ty].
                             Must be the same shape as tipper_array.
                             *default* is None
    :type tipper_err_array: np.ndarray((nf, 1, 2))


    :param frequency: array of frequencyuencies corresponding to the tipper elements.
                 Must be same length as tipper_array.
                 *default* is None
    :type frequency: np.ndarray(nf)


    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    frequency            array of frequencyuencies corresponding to elements of z
    rotation_angle  angle of which data is rotated by

    tipper          tipper array
    tipper_err       tipper error array
    =============== ===========================================================

    =============== ===========================================================
    Methods         Description
    =============== ===========================================================
    mag_direction   computes magnitude and direction of real and imaginary
                    induction arrows.
    amp_phase       computes amplitude and phase of Tx and Ty.
    rotate          rotates the data by the given angle
    =============== ===========================================================
    """

    def __init__(
        self, tipper_array=None, tipper_err_array=None, frequency=None
    ):
        """
        initialize
        """
        self._logger = get_mtpy_logger(self.__class__.__name__)
        self._tipper = tipper_array
        self._tipper_err = tipper_err_array
        self._frequency = frequency

        self.rotation_angle = 0.0
        if self.tipper is not None:
            self.rotation_angle = np.zeros((len(self.tipper)))
        self._amplitude = None
        self._amplitude_err = None
        self._phase = None
        self._phase_err = None

        self._mag_real = None
        self._mag_imag = None
        self._angle_real = None
        self._angle_imag = None
        self._mag_err = None
        self._angle_err = None

        if self._tipper is not None and self._frequency is not None:
            self.compute_amp_phase()
            self.compute_mag_direction()

    def __str__(self):
        lines = ["Induction Vector (Tippers)", "-" * 40]
        if self.frequency is not None:
            lines.append(
                f"\tNumber of frequencyuencies:  {self.frequency.size}"
            )
            lines.append(
                f"\tfrequencyuency range:        {self.frequency.min():.5E} -- {self.frequency.max():.5E} Hz"
            )
            lines.append(
                f"\tPeriod range:           {1/self.frequency.max():.5E} -- {1/self.frequency.min():.5E} s"
            )
            lines.append("")
            lines.append("\tElements:")
            for zz, ff in zip(self.tipper, self.frequency):
                lines.append(f"\tfrequencyuency: {ff:5E} Hz -- Period {1/ff} s")
                lines.append(
                    "\t\t"
                    + np.array2string(
                        zz,
                        formatter={
                            "complex_kind": lambda x: f"{x.real:.4e}{x.imag:+.4E}"
                        },
                    ).replace("\n", "\n\t\t")
                )
        else:
            if self.tipper is not None:
                lines.append("Elements:")
                for ff, zz in enumerate(self.tipper):
                    lines.append(f"\tIndex {ff}")
                    lines.append(
                        "\t\t"
                        + np.array2string(
                            zz,
                            formatter={
                                "complex_kind": lambda x: f"{x.real:.4e}{x.imag:+.4E}"
                            },
                        ).replace("\n", "\n\t\t")
                    )
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, Tipper):
            msg = f"Cannot compare {type(other)} with Tipper"
            self._logger.error(msg)
            raise MTpyError_Tipper(msg)
        if (self.tipper != other.tipper).all():
            return False
        if (self.frequency != other.frequency).all():
            return False
        if (self.tipper_err != other.tipper_err).all():
            return False
        return True

    def copy(self):
        return copy.deepcopy(self)

    # ==========================================================================
    # Define get/set and properties
    # ==========================================================================
    # ----frequency----------------------------------------------------------
    @property
    def frequency(self):
        return self._frequency

    @frequency.setter
    def frequency(self, frequency_arr):
        """
        Set the array of frequency.

        :param frequency_arr: array of frequencyunecies (Hz)
        :type frequency_arr: np.ndarray(num_frequencyuencies)
        """
        if frequency_arr is None:
            return
        self._frequency = np.array(frequency_arr, dtype="float")

        if self.tipper is not None:
            if self.tipper.shape[0] != len(self._frequency):
                msg = (
                    "New frequency array is not correct shape for existing z"
                    + "new: {self._frequency.size} != old: {self.tipper.shape[0]}"
                )
                self._logger.error(msg)
                raise MTpyError_Tipper
        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    @property
    def period(self):
        """
        periods in seconds
        """

        return 1.0 / self.frequency

    @period.setter
    def period(self, value):
        """
        setting periods will set the frequencyuencies
        """

        self.frequency = 1.0 / value

    # ---tipper--------------------------------------------------------------
    @property
    def tipper(self):
        return self._tipper

    @tipper.setter
    def tipper(self, tipper_array):
        """
        Set the attribute *tipper*

        :param tipper_array: tipper array in the shape of [Tx, Ty]
                             *default* is None
        :type tipper_array: np.ndarray((nf, 1, 2), dtype='complex')
        """
        if tipper_array is None:
            return
        if not isinstance(tipper_array, np.ndarray):
            tipper_array = np.array(tipper_array, dtype="complex")
        if not tipper_array.dtype in ["complex"]:
            tipper_array = tipper_array.astype("complex")
        # check to see if the new tipper array is the same shape as the old
        if (
            self._tipper is not None
            and self._tipper.shape != tipper_array.shape
        ):
            msg = (
                "Shape of new array does not match old.  "
                + f"new shape {tipper_array.shape} != "
                + f"old shape {self._tipper.shape}. "
                + "Make a new Tipper instance to be save."
            )
            self._logger.error(msg)
            raise MTpyError_Tipper(msg)
        if len(tipper_array.shape) == 3:
            if tipper_array.shape[1:3] == (1, 2):
                self._tipper = tipper_array
            else:
                msg = f"Input array must be shape (n, 1, 2) not {tipper_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Tipper(msg)
        elif len(tipper_array.shape) == 2:
            if tipper_array.shape == (1, 2):
                self._tipper = tipper_array.reshape((1, 1, 2))
                self._logger.debug(
                    "setting input tipper with shape (1, 2) to (1, 1, 2)"
                )
            else:
                msg = f"Input array must be shape (n, 1, 2) not {tipper_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Tipper(msg)
        else:
            msg = f"{tipper_array.shape} are not the correct dimensions, must be (n, 1, 2)"
            self._logger.error(msg)
            raise MTpyError_Tipper(msg)
        # neeed to set the rotation angle such that it is an array
        if self.rotation_angle is float:
            self.rotation_angle = np.repeat(
                self.rotation_angle, len(self._tipper)
            )
        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----tipper error---------------------------------------------------------
    @property
    def tipper_err(self):
        return self._tipper_err

    @tipper_err.setter
    def tipper_err(self, tipper_err_array):
        """
        Set the attribute *tipper_err*.

        :param tipper_err_array: array of estimated tipper errors
                                 in the shape of [Tx, Ty].
                                 Must be the same shape as tipper_array.
                                 *default* is None
        :type tipper_err_array: np.ndarray((nf, 1, 2))
        """
        if tipper_err_array is None:
            return
        if not isinstance(tipper_err_array, np.ndarray):
            tipper_err_array = np.array(tipper_err_array, dtype="float")
        if not tipper_err_array.dtype in ["float"]:
            tipper_err_array = tipper_err_array.astype("float")
        if len(tipper_err_array.shape) == 3:
            if not tipper_err_array.shape[1:3] == (1, 2):
                msg = f"Input array must be shape (n, 1, 2) not {tipper_err_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Tipper(msg)
        elif len(tipper_err_array.shape) == 2:
            if tipper_err_array.shape == (1, 2):
                tipper_err_array = tipper_err_array.reshape((1, 1, 2))
                self._logger.debug(
                    "setting input tipper with shape (1, 2) to (1, 1, 2)"
                )
            else:
                msg = f"Input array must be shape (n, 1, 2) not {tipper_err_array.shape}"
                self._logger.error(msg)
                raise MTpyError_Tipper(msg)
        else:
            msg = f"{tipper_err_array.shape} are not the correct dimensions, must be (n, 1, 2)"
            self._logger.error(msg)
            raise MTpyError_Tipper(msg)
        # check to see if the new tipper array is the same shape as the old
        if (
            self._tipper is not None
            and self._tipper.shape != tipper_err_array.shape
        ):
            raise MTpyError_Tipper(
                "Shape of new error array does not match old"
                + f"new shape {tipper_err_array.shape} != old shape {self._tipper.shape}"
            )
        self._tipper_err = tipper_err_array

        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----amplitude and phase
    def compute_amp_phase(self):
        """
        Sets attributes:
                        * *amplitude*
                        * *phase*
                        * *amplitude_err*
                        * *phase_err*

        values for resistivity are in in Ohm m and phase in degrees.
        """

        if self.tipper is None:
            # logging.error( 'tipper array is None - cannot calculate rho/phi')
            # print 'tipper array is None - cannot calculate rho/phi'
            return None
        self._amplitude_err = None
        self._phase_err = None
        if self.tipper_err is not None:
            self._amplitude_err = np.zeros(self.tipper_err.shape)
            self._phase_err = np.zeros(self.tipper_err.shape)
        self._amplitude = np.abs(self.tipper)
        self._phase = np.rad2deg(np.angle(self.tipper))

        if self.tipper_err is not None:
            for idx_f in range(len(self.tipper)):
                for jj in range(2):
                    if self.tipper_err is not None:
                        if type(self.tipper) == np.ma.core.MaskedArray:
                            if self.tipper.mask[idx_f, 0, jj]:
                                continue
                        r_err, phi_err = MTcc.propagate_error_rect2polar(
                            np.real(self.tipper[idx_f, 0, jj]),
                            self.tipper_err[idx_f, 0, jj],
                            np.imag(self.tipper[idx_f, 0, jj]),
                            self.tipper_err[idx_f, 0, jj],
                        )

                        self.amplitude_err[idx_f, 0, jj] = r_err
                        self._phase_err[idx_f, 0, jj] = phi_err

    def set_amp_phase(self, r_array, phi_array):
        """
        Set values for amplitude(r) and argument (phi - in degrees).

        Updates the attributes:
                        * tipper
                        * tipper_err

        """

        if self.tipper is not None:

            tipper_new = copy.copy(self.tipper)

            if self.tipper.shape != r_array.shape:
                self._logger.error(
                    'Error - shape of "r" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(r_array.shape), str(self.tipper.shape))
                )
                return
            if self.tipper.shape != phi_array.shape:
                self._logger.error(
                    'Error - shape of "phi" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(phi_array.shape), str(self.tipper.shape))
                )
                return
        else:

            tipper_new = np.zeros(r_array.shape, "complex")

            if r_array.shape != phi_array.shape:
                self._logger.error(
                    'Error - shape of "phi" array does not match shape '
                    + 'of "r" array: %s ; %s'
                    % (str(phi_array.shape), str(r_array.shape))
                )
                return
        # assert real array:
        if np.linalg.norm(np.imag(r_array)) != 0:
            self._logger.error('Error - array "r" is not real valued !')
            return
        if np.linalg.norm(np.imag(phi_array)) != 0:
            self._logger.error('Error - array "phi" is not real valued !')
            return
        for idx_f in range(len(r_array)):
            for jj in range(2):
                tipper_new[idx_f, 0, jj] = cmath.rect(
                    r_array[idx_f, 0, jj],
                    math.radians(phi_array[idx_f, 0, jj]),
                )
        self.tipper = tipper_new

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()
        self.compute_mag_direction()

    # ---------------------------------
    # properties
    @property
    def amplitude(self):
        return self._amplitude

    @property
    def phase(self):
        return self._phase

    @property
    def amplitude_err(self):
        return self._amplitude_err

    @property
    def phase_err(self):
        return self._phase_err

    # ----magnitude and direction----------------------------------------------
    def compute_mag_direction(self):
        """
        computes the magnitude and direction of the real and imaginary
        induction vectors.
        """

        if self.tipper is None:
            return None
        self._mag_real = np.sqrt(
            self.tipper[:, 0, 0].real ** 2 + self.tipper[:, 0, 1].real ** 2
        )
        self._mag_imag = np.sqrt(
            self.tipper[:, 0, 0].imag ** 2 + self.tipper[:, 0, 1].imag ** 2
        )

        self._mag_err = None
        self._angle_err = None
        # get the angle, need to make both parts negative to get it into the
        # parkinson convention where the arrows point towards the conductor

        self._angle_real = np.rad2deg(
            np.arctan2(-self.tipper[:, 0, 1].real, -self.tipper[:, 0, 0].real)
        )

        self._angle_imag = np.rad2deg(
            np.arctan2(-self.tipper[:, 0, 1].imag, -self.tipper[:, 0, 0].imag)
        )

        ## estimate error: THIS MAYBE A HACK
        if self.tipper_err is not None:
            self._mag_err = np.sqrt(
                self.tipper_err[:, 0, 0] ** 2 + self.tipper_err[:, 0, 1] ** 2
            )
            self._angle_err = (
                np.rad2deg(
                    np.arctan2(
                        self.tipper_err[:, 0, 0], self.tipper_err[:, 0, 1]
                    )
                )
                % 45
            )

    def set_mag_direction(self, mag_real, ang_real, mag_imag, ang_imag):
        """
        computes the tipper from the magnitude and direction of the real
        and imaginary components.

        Updates tipper

        No error propagation yet
        """

        self.tipper[:, 0, 0].real = np.sqrt(
            (mag_real**2 * np.arctan(ang_real) ** 2)
            / (1 - np.arctan(ang_real) ** 2)
        )

        self.tipper[:, 0, 1].real = np.sqrt(
            mag_real**2 / (1 - np.arctan(ang_real) ** 2)
        )

        self.tipper[:, 0, 0].imag = np.sqrt(
            (mag_imag**2 * np.arctan(ang_imag) ** 2)
            / (1 - np.arctan(ang_imag) ** 2)
        )

        self.tipper[:, 0, 1].imag = np.sqrt(
            mag_imag**2 / (1 - np.arctan(ang_imag) ** 2)
        )
        # for consistency recalculate mag and angle
        self.compute_mag_direction()
        self.compute_amp_phase()

    @property
    def mag_real(self):
        return self._mag_real

    @property
    def mag_imag(self):
        return self._mag_imag

    @property
    def angle_real(self):
        return self._angle_real

    @property
    def angle_imag(self):
        return self._angle_imag

    @property
    def mag_err(self):
        return self._mag_err

    @property
    def angle_err(self):
        return self._angle_err

    # ----rotate---------------------------------------------------------------
    def rotate(self, alpha):
        """
        Rotate  Tipper array.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North=0, positive in clockwise direction.
        (Mathematically negative!)

        In non-rotated state, 'X' refs to North and 'Y' to East direction.

        Updates the attributes:
                        * *tipper*
                        * *tipper_err*
                        * *rotation_angle*

        """

        if self.tipper is None:
            self._logger.error('tipper array is "None" - I cannot rotate that')
            return
        # check for iterable list/set of angles - if so, it must have length 1
        # or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self._logger.error(
                    '"Angle" must be a valid number (in degrees)'
                )
                return
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        elif len(alpha) == 1:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self._logger.error(
                    '"Angle" must be a valid number (in degrees)'
                )
                return
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        else:
            try:
                lo_angles = [float(ii % 360) for ii in alpha]
            except ValueError:
                self._logger.error(
                    '"Angles" must be valid numbers (in degrees)'
                )
                return
        self.rotation_angle = np.array(
            [
                (oldangle + lo_angles[ii]) % 360
                for ii, oldangle in enumerate(self.rotation_angle)
            ]
        )

        if len(lo_angles) != len(self.tipper):
            self._logger.error(
                'Wrong number Number of "angles" - need %ii '
                % (len(self.tipper))
            )
            self.rotation_angle = 0.0
            return
        tipper_rot = copy.copy(self.tipper)
        tipper_err_rot = copy.copy(self.tipper_err)

        for idx_frequency in range(len(tipper_rot)):
            angle = lo_angles[idx_frequency]

            if self.tipper_err is not None:
                (
                    tipper_rot[idx_frequency],
                    tipper_err_rot[idx_frequency],
                ) = MTcc.rotate_vector_with_errors(
                    self.tipper[idx_frequency, :, :],
                    angle,
                    self.tipper_err[idx_frequency, :, :],
                )
            else:
                (
                    tipper_rot[idx_frequency],
                    tipper_err_rot,
                ) = MTcc.rotate_vector_with_errors(
                    self.tipper[idx_frequency, :, :], angle
                )
        self.tipper = tipper_rot
        self.tipper_err = tipper_err_rot

        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()


# ------------------------
def correct4sensor_orientation(
    Z_prime, Bx=0, By=90, Ex=0, Ey=90, Z_prime_error=None
):
    """
    Correct a Z-array for wrong orientation of the sensors.

    Assume, E' is measured by sensors orientated with the angles
        E'x: a
        E'y: b

    Assume, B' is measured by sensors orientated with the angles
        B'x: c
        B'y: d

    With those data, one obtained the impedance tensor Z':
        E' = Z' * B'

    Now we define change-of-basis matrices T,U so that
        E = T * E'
        B = U * B'

    =>   T contains the expression of the E'-basis in terms of E
    (the standard basis)
    and  U contains the expression of the B'-basis in terms of B
    (the standard basis)
    The respective expressions for E'x-basis vector and E'y-basis
    vector are the columns of T.
    The respective expressions for B'x-basis vector and B'y-basis
    vector are the columns of U.

    We obtain the impedance tensor in default coordinates as:

    E' = Z' * B' => T^(-1) * E = Z' * U^(-1) * B
                 => E = T * Z' * U^(-1) * B
                 => Z = T * Z' * U^(-1)

    :param Z_prime: impedance tensor to be adjusted
    :dtype Z_prime: np.ndarray(num_frequency, 2, 2, dtype='complex')


    :param Bx: orientation of Bx relative to geographic north (0)
                                   *default* is 0
    :type Bx: float (angle in degrees)

    :param By:
    :type By: float (angle in degrees)
                         orientation of By relative to geographic north (0)
                                 *default* is 90

    :param Ex: orientation of Ex relative to geographic north (0)
                                   *default* is 0
    :type Ex: float (angle in degrees)

    :param Ey: orientation of Ey relative to geographic north (0)
                                  *default* is 90
    :type Ey: float (angle in degrees)

    :param Z_prime_error: impedance tensor error (std)
                                                 *default* is None
    :type Z_prime_error: np.ndarray(Z_prime.shape)

    :returns: adjusted impedance tensor
    :rtype: np.ndarray(Z_prime.shape, dtype='complex')

    :returns: impedance tensor standard deviation in
                                        default orientation
    :rtype: np.ndarray(Z_prime.shape, dtype='real')
    """
    try:
        if len(Z_prime.shape) != 2:
            raise
        if Z_prime.shape != (2, 2):
            raise
        if Z_prime.dtype not in ["complex", "float", "int"]:
            raise
        Z_prime = np.matrix(Z_prime)
    except:
        raise MTpyError_input_arguments(
            "ERROR - Z array not valid!" + "Must be 2x2 complex array"
        )
    if Z_prime_error is not None:
        try:
            if len(Z_prime_error.shape) != 2:
                raise
            if Z_prime_error.shape != (2, 2):
                raise
            if Z_prime_error.dtype not in ["float", "int"]:
                raise
        except:
            raise MTpyError_input_arguments(
                "ERROR - Z-error array not" + "valid! Must be 2x2 real array"
            )
    T = np.matrix(np.zeros((2, 2)))
    U = np.matrix(np.zeros((2, 2)))

    dummy1 = cmath.rect(1, math.radians(Ex))

    T[0, 0] = np.real(dummy1)
    T[1, 0] = np.imag(dummy1)
    dummy2 = cmath.rect(1, math.radians(Ey))
    T[0, 1] = np.real(dummy2)
    T[1, 1] = np.imag(dummy2)

    dummy3 = cmath.rect(1, math.radians(Bx))
    U[0, 0] = np.real(dummy3)
    U[1, 0] = np.imag(dummy3)
    dummy4 = cmath.rect(1, math.radians(By))
    U[0, 1] = np.real(dummy4)
    U[1, 1] = np.imag(dummy4)

    try:
        z_arr = np.array(np.dot(T, np.dot(Z_prime, U.I)))
    except:
        raise MTpyError_input_arguments(
            "ERROR - Given angles do not"
            + "define basis for 2 dimensions - cannot convert Z'"
        )
    z_err_arr = copy.copy(Z_prime_error)

    # TODO: calculate error propagation

    return z_arr, z_err_arr
