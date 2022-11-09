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
from .base import TFBase
from .pt import PhaseTensor

# ==============================================================================
# Impedance Tensor Class
# ==============================================================================
class Z(TFBase):
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

        super().__init__(
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
        if z is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf=z)
        else:
            self._dataset["transfer_function"].loc[self.comps] = z

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

        z_error = self._validate_array_input(z_error, "float", old_shape)
        if z_error is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf_error=z_error)
        else:
            self._dataset["transfer_function_error"].loc[self.comps] = z_error

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

        z_model_error = self._validate_array_input(
            z_model_error, "float", old_shape
        )

        if z_model_error is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf_error=z_model_error)
        else:
            self._dataset["transfer_function_model_error"].loc[
                self.comps
            ] = z_model_error

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
         components (ie z[:, 0, :]).  This is assumed to be in resistivity scale
        :type reduce_res_factor_x: float or iterable list or array
        :param reduce_res_factor_y: static shift factor to be applied to y
         components (ie z[:, 1, :]).  This is assumed to be in resistivity scale
        :type reduce_res_factor_y: float or iterable list or array
        :param inplace: Update the current object or return a new impedance
        :type inplace: boolean
        :returns: static shift matrix,
        :rtype: np.ndarray ((2, 2))
        :returns: corrected Z if inplace is False
        :rtype: mtpy.core.z.Z

        .. note:: The factors are in resistivity scale, so the entries of the
         matrix "S" need to be given by their square-roots!

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
                f"Length of reduce_res_factor_x needs to be {len(self.z)} "
                f" not {len(lo_x_factors)}"
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
                f"Length of reduce_res_factor_x needs to be {len(self.z)} "
                f" not {len(lo_x_factors)}"
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

        if inplace:
            self.z = z_corrected
            return static_shift
        else:

            return static_shift, z_corrected

    def remove_distortion(
        self, distortion_tensor, distortion_error_tensor=None, inplace=False
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
                self.logger.error(msg)
                raise ValueError(msg)
            if (
                len(distortion_tensor.shape) == 3
                or len(distortion_error_tensor.shape) == 3
            ):
                self.logger.info(
                    "Distortion is not time-dependent - taking only first"
                    + "of given distortion tensors"
                )
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_error_tensor = distortion_error_tensor[0]
                except IndexError:
                    msg = "Distortion tensor and error are not correct shape"
                    self.logger.error(msg)
                    raise ValueError(msg)
            if not (distortion_tensor.shape == (2, 2)) and (
                distortion_error_tensor.shape == (2, 2)
            ):
                msg = "Distortion tensor and error are not correct shape"
                self.logger.error(msg)
                raise ValueError(msg)
            distortion_tensor = np.matrix(np.real(distortion_tensor))
        except ValueError:
            msg = "Input distortion tensor, must be (2, 2)"
            raise ValueError(msg)
        try:
            DI = distortion_tensor.I
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
        z_corrected = np.zeros_like(self.z)
        z_corrected_error = np.zeros_like(self.z_error)

        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot(DI, np.matrix(self.z[idx_f])))
            for ii in range(2):
                for jj in range(2):
                    z_corrected_error[idx_f, ii, jj] = np.sum(
                        np.abs(
                            np.array(
                                [
                                    DI_error[ii, 0] * self.z[idx_f, 0, jj],
                                    DI[ii, 0] * self.z_error[idx_f, 0, jj],
                                    DI_error[ii, 1] * self.z[idx_f, 1, jj],
                                    DI[ii, 1] * self.z_error[idx_f, 1, jj],
                                ]
                            )
                        )
                    )
        return distortion_tensor, z_corrected, z_corrected_error

    @property
    def resistivity(self):
        if self.z is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2, 0, self.z
            )

    @property
    def phase(self):
        if self.z is not None:
            return np.rad2deg(np.angle(self.z))

    @property
    def resistivity_error(self):
        if self.z is not None and self.z_error is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2,
                0,
                self.z_error,
            )

    @property
    def phase_error(self):
        if self.z is not None and self.z_error is not None:
            return np.degrees(
                np.arctan(self.resistivity_error / self.resistivity)
            )

    @property
    def resistivity_model_error(self):
        if self.z is not None and self.z_model_error is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2,
                0,
                self.z_model_error,
            )

    @property
    def phase_model_error(self):
        if self.z is not None and self.z_model_error is not None:
            return np.degrees(
                np.arctan(self.resistivity_model_error / self.resistivity)
            )

    def _compute_z_error(self, res, res_error, phase, phase_error):
        if res_error is None:
            return None
        return abs(
            np.sqrt(5.0 * self.frequency * (res_error.T)).T
            * np.exp(1j * np.radians(phase_error))
        )

    def set_resistivity_phase(
        self,
        resistivity,
        phase,
        frequency,
        res_error=None,
        phase_error=None,
        res_model_error=None,
        phase_model_error=None,
    ):
        """
        Set values for resistivity (res - in Ohm m) and phase
        (phase - in degrees), including error propagation.


        :param resistivity: resistivity array in Ohm-m
        :type resistivity: np.ndarray(num_frequency, 2, 2)

        :param phase: phase array in degrees
        :type phase: np.ndarray(num_frequency, 2, 2)

        :param frequency: frequency array in Hz
        :type frequency: np.ndarray(num_frequency)

        :param res_error: resistivity error array in Ohm-m
        :type res_error: np.ndarray(num_frequency, 2, 2)

        :param phase_error: phase error array in degrees
        :type phase_error: np.ndarray(num_frequency, 2, 2)


        """

        if resistivity is None or phase is None or frequency is None:
            self.logger.debug(
                "Cannot estimate resitivity and phase if resistivity, "
                "phase, or frequency is None."
            )
            return

        self.frequency = self._validate_frequency(frequency)
        resistivity = self._validate_array_input(resistivity, float)
        phase = self._validate_array_input(phase, float)

        res_error = self._validate_array_input(res_error, float)
        phase_error = self._validate_array_input(phase_error, float)
        res_model_error = self._validate_array_input(res_model_error, float)
        phase_model_error = self._validate_array_input(
            phase_model_error, float
        )

        abs_z = np.sqrt(5.0 * self.frequency * (resistivity.T)).T
        self.z = abs_z * np.exp(1j * np.radians(phase))

        self.z_error = self._compute_z_error(
            resistivity, res_error, phase, phase_error
        )
        self.z_model_error = self._compute_z_error(
            resistivity, res_model_error, phase, phase_model_error
        )

    @property
    def det(self):
        """
        Return the determinant of Z

        :returns: det_Z
        :rtype: np.ndarray(nfrequency)
        """
        if self.z is not None:
            det_z = np.array([np.linalg.det(ii) ** 0.5 for ii in self.z])

            return det_z

    @property
    def det_error(self):
        """
        Return the determinant of Z error

        :returns: det_z_error
        :rtype: np.ndarray(nfrequency)
        """
        det_z_error = None
        if self.z_error is not None:
            det_z_error = np.zeros_like(self.det, dtype=np.float)
            with np.errstate(invalid="ignore"):
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
                ) ** 0.5
        return det_z_error

    @property
    def det_model_error(self):
        """
        Return the determinant of Z error

        :returns: det_z_error
        :rtype: np.ndarray(nfrequency)
        """
        det_z_error = None
        if self.z_model_error is not None:
            det_z_error = np.zeros_like(self.det, dtype=np.float)
            with np.errstate(invalid="ignore"):
                # components of the impedance tensor are not independent variables
                # so can't use standard error propagation
                # calculate manually:
                # difference of determinant of z + z_error and z - z_error then divide by 2
                det_z_error[:] = (
                    np.abs(
                        np.linalg.det(self.z + self.z_model_error)
                        - np.linalg.det(self.z - self.z_model_error)
                    )
                    / 2.0
                ) ** 0.5
        return det_z_error

    @property
    def phase_det(self):
        if self.det is not None:
            return np.rad2deg(np.arctan2(self.det.imag, self.det.real))

    @property
    def phase_error_det(self):
        if self.det is not None:
            return np.rad2deg(np.arcsin(self.det_error / abs(self.det)))

    @property
    def phase_model_error_det(self):
        if self.det is not None:
            return np.rad2deg(np.arcsin(self.det_model_error / abs(self.det)))

    @property
    def res_det(self):
        if self.det is not None:
            return 0.2 * (1.0 / self.frequency) * abs(self.det) ** 2

    @property
    def res_error_det(self):
        if self.det_error is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self.det + self.det_error) ** 2
                - self.res_det
            )

    @property
    def res_model_error_det(self):
        if self.det_model_error is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self.det + self.det_model_error) ** 2
                - self.res_det
            )

    def _get_component(self, comp, array):
        if array is not None:
            index_dict = {"x": 0, "y": 1}
            ii = index_dict[comp[-2]]
            jj = index_dict[comp[-1]]

            return array[:, ii, jj]

    @property
    def res_xx(self):
        return self._get_component("xx", self.resistivity)

    @property
    def res_xy(self):
        return self._get_component("xy", self.resistivity)

    @property
    def res_yx(self):
        return self._get_component("yx", self.resistivity)

    @property
    def res_yy(self):
        return self._get_component("yy", self.resistivity)

    @property
    def res_error_xx(self):
        return self._get_component("xx", self.resistivity_error)

    @property
    def res_error_xy(self):
        return self._get_component("xy", self.resistivity_error)

    @property
    def res_error_yx(self):
        return self._get_component("yx", self.resistivity_error)

    @property
    def res_error_yy(self):
        return self._get_component("yy", self.resistivity_error)

    @property
    def res_model_error_xx(self):
        return self._get_component("xx", self.resistivity_model_error)

    @property
    def res_model_error_xy(self):
        return self._get_component("xy", self.resistivity_model_error)

    @property
    def res_model_error_yx(self):
        return self._get_component("yx", self.resistivity_model_error)

    @property
    def res_model_error_yy(self):
        return self._get_component("yy", self.resistivity_model_error)

    @property
    def phase_xx(self):
        return self._get_component("xx", self.phase)

    @property
    def phase_xy(self):
        return self._get_component("xy", self.phase)

    @property
    def phase_yx(self):
        return self._get_component("yx", self.phase)

    @property
    def phase_yy(self):
        return self._get_component("yy", self.phase)

    @property
    def phase_error_xx(self):
        return self._get_component("xx", self.phase_error)

    @property
    def phase_error_xy(self):
        return self._get_component("xy", self.phase_error)

    @property
    def phase_error_yx(self):
        return self._get_component("yx", self.phase_error)

    @property
    def phase_error_yy(self):
        return self._get_component("yy", self.phase_error)

    @property
    def phase_model_error_xx(self):
        return self._get_component("xx", self.phase_model_error)

    @property
    def phase_model_error_xy(self):
        return self._get_component("xy", self.phase_model_error)

    @property
    def phase_model_error_yx(self):
        return self._get_component("yx", self.phase_model_error)

    @property
    def phase_model_error_yy(self):
        return self._get_component("yy", self.phase_model_error)

    @property
    def phase_tensor(self):
        return PhaseTensor(
            z=self.z,
            z_error=self.z_error,
            z_model_error=self.z_model_error,
            frequency=self.frequency,
        )
