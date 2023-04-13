#!/usr/bin/env python

"""
Z
===

Container for the Impedance Tensor 

Originally written by Jared Peacock Lars Krieger
Updated 2022 by J. Peacock to work with new framework

"""

# =============================================================================
# Imports
# =============================================================================
import copy
import numpy as np

import mtpy.utils.calculator as MTcc
from .base import TFBase
from .pt import PhaseTensor
from .z_analysis import (
    ZInvariants,
    find_distortion,
    remove_distortion_from_z_object,
    calculate_depth_of_investigation,
)

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
    :param frequency: array of frequencyuency values corresponding to impedance
     tensor elements.
    :type frequency: np.ndarray(n_frequency)

    Create Impedance from scracth
    ------------------------------

    >>> import mtpy.core import Z
    >>> import numpy as np
    >>> z_test = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
    >>> z_object = Z(z=z_test, frequency=[1])
    >>> z_object.rotate(45)
    >>> z_object.resistivity

    Create from resistivity and phase
    -----------------------------------

    >>> z_object = Z()
    >>> z_object.set_resistivity_phase(
        np.array([[5, 100], [100, 5]]),
        np.array([[90, 45], [-135, -90]]),
        np.array([1])
        )
    >>> z_object.z
    array([[[ 3.06161700e-16 +5.j, 1.58113883e+01+15.8113883j],
            [-1.58113883e+01-15.8113883j, 3.06161700e-16 -5.j ]]])

    """

    def __init__(
        self,
        z=None,
        z_error=None,
        frequency=None,
        z_model_error=None,
    ):
        """
        Initialize an instance of the Z class.

        :param z: array containing complex impedance values
        :type z: numpy.ndarray(n_frequency, 2, 2)
        :param z_error: array containing error values (standard deviation)
         of impedance tensor elements
        :type z_error: numpy.ndarray(n_frequency, 2, 2)
        :param frequency: array of frequencyuency values corresponding to impedance
         tensor elements.
        :type frequency: np.ndarray(n_frequency)

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
        elif self._has_frequency():
            old_shape = (
                self.frequency.size,
                self._expected_shape[0],
                self._expected_shape[1],
            )
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
        """error of impedance tensor array as standard deviation"""
        if self._has_tf_error():
            return self._dataset.transfer_function_error.values

    @z_error.setter
    def z_error(self, z_error):
        """
        Set the attribute z_error

        :param z_error: error of impedance tensor array as standard deviation
        :type z_error: np.ndarray(nfrequency, 2, 2)
        """
        old_shape = None
        if not self._has_tf_error():
            old_shape = self._dataset.transfer_function_error.shape
        elif self._has_frequency():
            old_shape = (
                self.frequency.size,
                self._expected_shape[0],
                self._expected_shape[1],
            )

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
        """model error of impedance tensor array as standard deviation"""
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

        elif self._has_frequency():
            old_shape = (
                self.frequency.size,
                self._expected_shape[0],
                self._expected_shape[1],
            )

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

        def _validate_factor_single(factor):
            try:
                x_factor = float(factor)
            except ValueError:
                msg = f"factor must be a valid number not {factor}"
                self.logger.error(msg)
                raise ValueError(msg)
            return np.repeat(x_factor, len(self.z))

        def _validate_ss_input(factor):
            if not np.iterable(factor):
                x_factor = _validate_factor_single(factor)

            elif len(reduce_res_factor_x) == 1:
                x_factor = _validate_factor_single(factor)
            else:
                x_factor = np.array(factor, dtype=float)

            if len(x_factor) != len(self.z):
                msg = (
                    f"Length of reduce_res_factor_x needs to be {len(self.z)} "
                    f" not {len(x_factor)}"
                )
                self.logger.error(msg)
                raise ValueError(msg)
            return x_factor

        x_factors = np.sqrt(_validate_ss_input(reduce_res_factor_x))
        y_factors = np.sqrt(_validate_ss_input(reduce_res_factor_y))

        z_corrected = copy.copy(self.z)

        z_corrected[:, 0, 0] = self.z[:, 0, 0] / x_factors
        z_corrected[:, 0, 1] = self.z[:, 0, 1] / x_factors
        z_corrected[:, 1, 0] = self.z[:, 1, 0] / y_factors
        z_corrected[:, 1, 1] = self.z[:, 1, 1] / y_factors

        if inplace:
            self.z = z_corrected
        else:
            return Z(
                z=z_corrected,
                z_error=self.z_error,
                frequency=self.frequency,
                z_model_error=self.z_model_error,
            )

    def remove_distortion(
        self,
        distortion_tensor=None,
        distortion_error_tensor=None,
        n_frequencies=None,
        comp="det",
        only_2d=False,
        inplace=False,
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

        if distortion_tensor is None:
            (
                distortion_tensor,
                distortion_error_tensor,
            ) = self.estimate_distortion(
                n_frequencies=n_frequencies, comp=comp, only_2d=only_2d
            )

        z_corrected, z_corrected_error = remove_distortion_from_z_object(
            self, distortion_tensor, distortion_error_tensor, self.logger
        )

        if inplace:
            self.z = z_corrected
            self.z_error = z_corrected_error
        else:
            return Z(
                z=z_corrected,
                z_error=z_corrected_error,
                frequency=self.frequency,
                z_model_error=self.z_model_error,
            )

    @property
    def resistivity(self):
        """resistivity of impedance"""
        if self.z is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2, 0, self.z
            )

    @property
    def phase(self):
        """phase of impedance"""
        if self.z is not None:
            return np.rad2deg(np.angle(self.z))

    @property
    def resistivity_error(self):
        """resistivity error of impedance"""
        if self.z is not None and self.z_error is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2,
                0,
                self.z_error,
            )

    @property
    def phase_error(self):
        """phase error of impedance"""
        if self.z is not None and self.z_error is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.degrees(
                    np.arctan(self.resistivity_error / self.resistivity)
                )

    @property
    def resistivity_model_error(self):
        """resistivity model error of impedance"""
        if self.z is not None and self.z_model_error is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2,
                0,
                self.z_model_error,
            )

    @property
    def phase_model_error(self):
        """phase model error of impedance"""
        if self.z is not None and self.z_model_error is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.degrees(
                    np.arctan(self.resistivity_model_error / self.resistivity)
                )

    def _compute_z_error(self, res, res_error, phase, phase_error):
        """
        Compute z error from apparent resistivity and phase.

        :param res: resistivity array
        :type res: np.ndarray
        :param res_error: resistivity error array
        :type res_error: np.ndarray
        :param phase: phase array in degrees
        :type phase: np.ndarray
        :param phase_error: phase error array in degrees
        :type phase_error: np.ndarray
        :return: impedance error as a float
        :rtype: np.ndarray

        """
        if res_error is None:
            return None
        return np.abs(
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

        .. note:: The error propogation is causal, meaning the apparent
         resistivity error and phase error are linked through a Taylor exampsion
         approximation where the phase error is estimated from the apparent
         resistivity error.  Therefore if you set the phase error
         you will likely not get back the same phase error.

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
        phase_model_error = self._validate_array_input(phase_model_error, float)

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
        """determinant of impedance"""
        if self.z is not None:
            det_z = np.array([np.linalg.det(ii) ** 0.5 for ii in self.z])

            return det_z

    @property
    def det_error(self):
        """
        Return the determinant of impedance error
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
        Return the determinant of impedance model error
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
        """phase determinant"""
        if self.det is not None:
            return np.rad2deg(np.arctan2(self.det.imag, self.det.real))

    @property
    def phase_error_det(self):
        """phase error determinant"""
        if self.det is not None:
            return np.rad2deg(np.arcsin(self.det_error / abs(self.det)))

    @property
    def phase_model_error_det(self):
        """phase model error determinant"""
        if self.det is not None:
            return np.rad2deg(np.arcsin(self.det_model_error / abs(self.det)))

    @property
    def res_det(self):
        """resistivity determinant"""
        if self.det is not None:
            return 0.2 * (1.0 / self.frequency) * abs(self.det) ** 2

    @property
    def res_error_det(self):
        """resistivity error determinant"""
        if self.det_error is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self.det + self.det_error) ** 2
                - self.res_det
            )

    @property
    def res_model_error_det(self):
        """resistivity model error determinant"""
        if self.det_model_error is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self.det + self.det_model_error) ** 2
                - self.res_det
            )

    def _get_component(self, comp, array):
        """
        Get the correct component from an array

        :param comp: [ xx | xy | yx | yy ]
        :type comp: string
        :param array: impedance array
        :type array: np.ndarray
        :return: array component
        :rtype: np.ndarray

        """
        if array is not None:
            index_dict = {"x": 0, "y": 1}
            ii = index_dict[comp[-2]]
            jj = index_dict[comp[-1]]

            return array[:, ii, jj]

    @property
    def res_xx(self):
        """resistivity of xx component"""
        return self._get_component("xx", self.resistivity)

    @property
    def res_xy(self):
        """resistivity of xy component"""
        return self._get_component("xy", self.resistivity)

    @property
    def res_yx(self):
        """resistivity of yx component"""
        return self._get_component("yx", self.resistivity)

    @property
    def res_yy(self):
        """resistivity of yy component"""
        return self._get_component("yy", self.resistivity)

    @property
    def res_error_xx(self):
        """resistivity error of xx component"""
        return self._get_component("xx", self.resistivity_error)

    @property
    def res_error_xy(self):
        """resistivity error of xy component"""
        return self._get_component("xy", self.resistivity_error)

    @property
    def res_error_yx(self):
        """resistivity error of yx component"""
        return self._get_component("yx", self.resistivity_error)

    @property
    def res_error_yy(self):
        """resistivity error of yy component"""
        return self._get_component("yy", self.resistivity_error)

    @property
    def res_model_error_xx(self):
        """resistivity model error of xx component"""
        return self._get_component("xx", self.resistivity_model_error)

    @property
    def res_model_error_xy(self):
        """resistivity model error of xy component"""
        return self._get_component("xy", self.resistivity_model_error)

    @property
    def res_model_error_yx(self):
        """resistivity model error of yx component"""
        return self._get_component("yx", self.resistivity_model_error)

    @property
    def res_model_error_yy(self):
        """resistivity model error of yy component"""
        return self._get_component("yy", self.resistivity_model_error)

    @property
    def phase_xx(self):
        """phase of xx component"""
        return self._get_component("xx", self.phase)

    @property
    def phase_xy(self):
        """phase of xy component"""
        return self._get_component("xy", self.phase)

    @property
    def phase_yx(self):
        """phase of yx component"""
        return self._get_component("yx", self.phase)

    @property
    def phase_yy(self):
        """phase of yy component"""
        return self._get_component("yy", self.phase)

    @property
    def phase_error_xx(self):
        """phase error of xx component"""
        return self._get_component("xx", self.phase_error)

    @property
    def phase_error_xy(self):
        """phase error of xy component"""
        return self._get_component("xy", self.phase_error)

    @property
    def phase_error_yx(self):
        """phase error of yx component"""
        return self._get_component("yx", self.phase_error)

    @property
    def phase_error_yy(self):
        """phase error of yy component"""
        return self._get_component("yy", self.phase_error)

    @property
    def phase_model_error_xx(self):
        """phase model error of xx component"""
        return self._get_component("xx", self.phase_model_error)

    @property
    def phase_model_error_xy(self):
        """phase model error of xy component"""
        return self._get_component("xy", self.phase_model_error)

    @property
    def phase_model_error_yx(self):
        """phase model error of yx component"""
        return self._get_component("yx", self.phase_model_error)

    @property
    def phase_model_error_yy(self):
        """phase model error of yy component"""
        return self._get_component("yy", self.phase_model_error)

    @property
    def phase_tensor(self):
        """Phase tensor object based on impedance"""
        return PhaseTensor(
            z=self.z,
            z_error=self.z_error,
            z_model_error=self.z_model_error,
            frequency=self.frequency,
        )

    @property
    def invariants(self):
        """Weaver Invariants"""
        return ZInvariants(z=self.z)

    def estimate_dimensionality(
        self, skew_threshold=5, eccentricity_threshold=0.1
    ):
        """
        Estimate dimensionality of the impedance tensor from parameters such
        as strike and phase tensor eccentricity

        :return: DESCRIPTION
        :rtype: TYPE

        """

        dimensionality = np.ones(self.period.size, dtype=int)

        # need to get 2D first then 3D
        dimensionality[
            np.where(self.phase_tensor.eccentricity > eccentricity_threshold)
        ] = 2
        dimensionality[
            np.where(np.abs(self.phase_tensor.skew) > skew_threshold)
        ] = 3

        return dimensionality

    def estimate_distortion(
        self,
        n_frequencies=None,
        comp="det",
        only_2d=False,
    ):
        """

        :param n_frequencies: DESCRIPTION, defaults to 20
        :type n_frequencies: TYPE, optional
        :param comp: DESCRIPTION, defaults to "det"
        :type comp: TYPE, optional
        :param : DESCRIPTION
        :type : TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if n_frequencies is None:
            nf = self.frequency.size
        else:
            nf = n_frequencies

        if self._has_tf():
            new_z_object = Z(
                z=self.z[0:nf, :, :],
                frequency=self.frequency[0:nf],
            )
            if self._has_tf_error():
                new_z_object.z_error = self.z_error[0:nf]

        return find_distortion(new_z_object, comp=comp, only_2d=only_2d)

    def estimate_depth_of_investigation(self):
        """
        estimate depth of investigation

        :return: DESCRIPTION
        :rtype: TYPE

        """

        return calculate_depth_of_investigation(self)
