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
            self.logger.error(msg)
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
                self.logger.error(msg)
                raise MTpyError_Z(msg)

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
            self.logger.error(msg)
            raise MTpyError_Z(msg)

        if len(z_array.shape) == 3:
            if z_array.shape[1:3] == (2, 2):
                return z_array
            else:
                msg = f"Input array must be shape (n, 2, 2) not {z_array.shape}"
                self.logger.error(msg)
                raise MTpyError_Z(msg)
        elif len(z_array.shape) == 2:
            if z_array.shape == (2, 2):
                return z_array.reshape((1, 2, 2))
                self.logger.debug(
                    "setting input z with shape (2, 2) to (1, 2, 2)"
                )
            else:
                msg = f"Input array must be shape (n, 2, 2) not {z_array.shape}"
                self.logger.error(msg)
                raise MTpyError_Z(msg)
        else:
            msg = f"{z_array.shape} are not the correct dimensions, must be (n, 2, 2)"
            self.logger.error(msg)
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

    @property
    def inverse(self):
        """
        Return the inverse of Z.

        (no error propagtaion included yet)

        """

        if self.z is None:
            self.logger.warn('z array is "None" - I cannot invert that')
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
            self.logger.warning('Z array is "None" and cannot be rotated')
            return
        # check for iterable list/set of angles - if so, it must have length
        # 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                msg = f"Angle must be a valid number (in degrees) not {alpha}"
                self.logger.error(msg)
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
                    self.logger.error(msg)
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
                    self.logger.error(msg)
                    raise MTpyError_Z(msg)
        self.rotation_angle = np.array(
            [
                (oldangle + lo_angles[ii]) % 360
                for ii, oldangle in enumerate(self.rotation_angle)
            ]
        )

        if len(lo_angles) != len(self.z):
            msg = f"Wrong number of angles, need {len(self.z)}"
            self.logger.error(msg)
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
