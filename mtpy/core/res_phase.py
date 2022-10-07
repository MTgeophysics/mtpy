# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 10:25:59 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import mtpy.utils.calculator as MTcc
from mtpy.utils.exceptions import (
    MTpyError_Z,
    MTpyError_input_arguments,
)
from mtpy.utils.mtpy_logger import get_mtpy_logger


# ==============================================================================
# Resistivity and phase object
# ==============================================================================
class ResPhase:
    """
    resistivity and phase container with convenience property attributes to
    access the different components.

    """

    def __init__(
        self,
        z_array=None,
        z_err_array=None,
        frequency=None,
        z_model_err=None,
        **kwargs,
    ):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self.resistivity = None
        self.phase = None

        self.resistivity_err = None
        self.phase_err = None

        self.resistivity_model_err = None
        self.phase_model_err = None

        self.frequency = frequency

        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        lines = ["Resistivity and Phase", "-" * 30]
        if self.frequency is not None:
            lines.append(
                f"\tNumber of frequencyuencies:  {self.frequency.size}"
            )
            lines.append(
                f"\tfrequency range:        {self.frequency.min():.5E} -- {self.frequency.max():.5E} Hz"
            )
            lines.append(
                f"\tPeriod range:           {1/self.frequency.max():.5E} -- {1/self.frequency.min():.5E} s"
            )
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __getattr__(self, name):
        """
        overload get attribute for components
        """

        index_dict = {"x": 0, "y": 1}
        if name in ["res_xx", "res_xy", "res_yx", "res_yy"]:
            if self.resistivity is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.resistivity[:, ii, jj]

        elif name in ["res_err_xx", "res_err_xy", "res_err_yx", "res_err_yy"]:
            if self.resistivity_err is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.resistivity_err[:, ii, jj]

        elif name in [
            "res_model_err_xx",
            "res_model_err_xy",
            "res_model_err_yx",
            "res_model_err_yy",
        ]:
            if self.resistivity_model_err is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.resistivity_model_err[:, ii, jj]

        elif name in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
            if self.phase is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.phase[:, ii, jj]

        elif name in [
            "phase_err_xx",
            "phase_err_xy",
            "phase_err_yx",
            "phase_err_yy",
        ]:
            if self.phase_err is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.phase_err[:, ii, jj]

        elif name in [
            "phase_model_err_xx",
            "phase_model_err_xy",
            "phase_model_err_yx",
            "phase_model_err_yy",
        ]:
            if self.phase_model_err is not None:
                ii = index_dict[name[-2]]
                jj = index_dict[name[-1]]
                return self.phase_model_err[:, ii, jj]

        else:
            return super().__getattribute__(name)

    def __setattr__(self, name, value):
        """
        overload get attribute for components
        """

        index_dict = {"x": 0, "y": 1}
        if name in ["res_xx", "res_xy", "res_yx", "res_yy"]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.resistivity[:, ii, jj] = self._validate_input_component(value)

        elif name in ["res_err_xx", "res_err_xy", "res_err_yx", "res_err_yy"]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.resistivity_err[:, ii, jj] = self._validate_input_component(
                value
            )

        elif name in [
            "res_model_err_xx",
            "res_model_err_xy",
            "res_model_err_yx",
            "res_model_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.resistivity_model_err[
                :, ii, jj
            ] = self._validate_input_component(value)

        elif name in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.phase[:, ii, jj] = self._validate_input_component(value)

        elif name in [
            "phase_err_xx",
            "phase_err_xy",
            "phase_err_yx",
            "phase_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.phase_err[:, ii, jj] = self._validate_input_component(value)

        elif name in [
            "phase_model_err_xx",
            "phase_model_err_xy",
            "phase_model_err_yx",
            "phase_model_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            self.phase_model_err[:, ii, jj] = self._validate_input_component(
                value
            )

        elif name in [
            "resistivity",
            "resistivity_err",
            "resistivity_model_err",
            "phase",
            "phase_err",
            "phase_model_err",
        ]:
            value = self._validate_real_valued(
                self._validate_input_array(value)
            )
            super().__setattr__(name, value)

        else:
            return super().__setattr__(name, value)

    def _validate_frequency(self, value):
        if value is None:
            return None

        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) > 1:
            value = value.flatten()

    def _validate_input_array(self, value):
        """
        validate an input array
        """
        if value is None:
            return None
        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) != 3 and len(value.shape) == 2:
            value = value.reshape((1, 2, 2))

        if value.shape[1] != 2 or value.shape[2] != 2:
            raise ValueError(
                f"Input array must have shape (nf, 2, 2) not {value.shape}"
            )
        if self.frequency is not None:
            if value.shape[0] != self.frequency.size:
                raise ValueError(
                    "Input array must have same number of frequencyuencies "
                    f"({self.frequency.size}) as current array "
                    f"not {value.shape[0]}"
                )
        return value

    def _validate_input_component(self, value):
        """
        validate a single component input

        needs to be same shape as frequency if present and 1D

        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if value is None:
            return None

        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) > 1:
            value = value.flatten()

        if self.frequency is not None:
            if value.shape != self.frequency.shape:
                raise ValueError(
                    "Input array must have same number of frequencyuencies "
                    f"({self.frequency.size}) as current array "
                    f"not {value.shape[0]}"
                )

        return value

    def _validate_real_valued(self, value):
        """
        make sure resistivity is real valued

        :param res: DESCRIPTION
        :type res: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        # assert real array:
        if np.linalg.norm(np.imag(value)) != 0:
            msg = "Array is not real valued"
            self.logger.error(msg)
            raise ValueError(msg)
        return value

    @property
    def frequency(self):
        if self._frequency is not None:
            return self._frequency

    @frequency.setter
    def frequency(self, value):
        if value is None:
            self._frequency = value
            return

        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) > 1:
            value = value.flatten()

        if self._frequency is not None:
            if value.shape != self._frequency.shape:
                raise ValueError(
                    "Input array must have same number of frequencyuencies "
                    f"({self.frequency.size}) as current array "
                    f"not {value.shape[0]}"
                )
        self._frequency = value
        self.set_res_phase()

    @property
    def period(self):
        if self.frequency is not None:
            return 1.0 / self.frequency

    @period.setter
    def period(self, value):
        self.frequency = 1.0 / value

    def from_impedance(
        self, z_array, z_err_array, frequency, z_model_err_array
    ):
        """

        :param z_array: DESCRIPTION
        :type z_array: TYPE
        :param z_err_array: DESCRIPTION
        :type z_err_array: TYPE
        :param frequency: DESCRIPTION
        :type frequency: TYPE
        :param z_model_err_array: DESCRIPTION
        :type z_model_err_array: TYPE
        :raises ValueError: DESCRIPTION
        :return: DESCRIPTION
        :rtype: TYPE

        """

        self.frequency = self._validate_frequency(frequency)

        # The _z_err can be None!!!
        if z_array is None or frequency is None:
            if z_array is None:
                msg = "z values are None, cannot compute parameters"
            elif frequency is None:
                msg = "frequency values are None, cannot compute parameters"
            self.logger.debug(msg)
            raise ValueError(msg)

        # if we set these to self., then we will be constantly updating, so
        # compute first and then set.
        self.resistivity = np.apply_along_axis(
            lambda x: np.abs(x) ** 2 / frequency * 0.2, 0, z_array
        )
        self.phase = np.rad2deg(np.angle(z_array))

        # Compute Errors
        self.resistivity_err = np.zeros_like(self.resistivity, dtype=np.float)
        self.phase_err = np.zeros_like(self.phase, dtype=np.float)
        self.resistivity_model_err = np.zeros_like(
            self.resistivity, dtype=np.float
        )
        self.phase_model_err = np.zeros_like(phase, dtype=np.float)

        # calculate resistivity and phase
        if z_err_array is not None:
            for idx_f in range(frequency.size):
                for ii in range(2):
                    for jj in range(2):
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            z_array[idx_f, ii, jj].real,
                            z_array[idx_f, ii, jj].imag,
                            z_err_array[idx_f, ii, jj],
                        )
                        self.resistivity_err[idx_f, ii, jj] = (
                            resistivity[idx_f, ii, jj] * r_err
                        )

                        self.phase_err[idx_f, ii, jj] = phi_err

        if z_model_err_array is not None:
            for idx_f in range(frequency.size):
                for ii in range(2):
                    for jj in range(2):
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            z_array[idx_f, ii, jj].real,
                            z_array[idx_f, ii, jj].imag,
                            z_model_err_array[idx_f, ii, jj],
                        )
                        self.resistivity_model_err[idx_f, ii, jj] = (
                            resistivity[idx_f, ii, jj] * r_err
                        )

                        self.phase_model_err[idx_f, ii, jj] = phi_err

    def compute_impedance(
        self,
        res_array=None,
        phase_array=None,
        frequency=None,
        res_err_array=None,
        phase_err_array=None,
        res_model_err=None,
        phase_model_err=None,
    ):
        """
        Set values for resistivity (res - in Ohm m) and phase
        (phase - in degrees), including error propagation.


        :param res_array: resistivity array in Ohm-m
        :type res_array: np.ndarray(num_frequency, 2, 2)

        :param phase_array: phase array in degrees
        :type phase_array: np.ndarray(num_frequency, 2, 2)

        :param frequency: frequency array in Hz
        :type frequency: np.ndarray(num_frequency)

        :param res_err_array: resistivity error array in Ohm-m
        :type res_err_array: np.ndarray(num_frequency, 2, 2)

        :param phase_err_array: phase error array in degrees
        :type phase_err_array: np.ndarray(num_frequency, 2, 2)


        """

        self.logger.debug("Resetting z, z_err, z_model_err")

        if res_array is not None:
            self.resistivity = res_array
        if phase_array is not None:
            self.phase = phase_array
        if frequency is not None:
            self.frequency = frequency
        if res_err_array is not None:
            self.resistivity_err = res_err_array
        if phase_err_array is not None:
            self.phase_err = phase_err_array
        if res_model_err is not None:
            self.resistivity_model_err = res_model_err
        if phase_model_err is not None:
            self.phase_model_err = phase_model_err

        if (
            self.resistivity is None
            or self.phase is None
            or self.frequency is None
        ):
            self.logger.debug(
                "Cannot estimate resitivity and phase if resistivity, "
                "phase, or frequency is None."
            )
            return

        abs_z = np.sqrt(5.0 * self.frequency * (self.resistivity.T)).T
        self._z = abs_z * np.exp(1j * np.radians(self.phase))

        self._z_err = np.zeros_like(self._z, dtype=np.float)
        self._z_model_err = np.zeros_like(self._z, dtype=np.float)
        # ---------------------------
        # error propagation:
        if self.resistivity_err is not None or self.phase_err is not None:
            for idx_f in range(self.frequency.shape[0]):
                for ii in range(2):
                    for jj in range(2):
                        abs_z = np.sqrt(
                            5
                            * self.frequency[idx_f]
                            * self.resistivity[idx_f, ii, jj]
                        )
                        rel_error_res = (
                            self.resistivity_err[idx_f, ii, jj]
                            / self.resistivity[idx_f, ii, jj]
                        )
                        # relative error varies by a factor of 0.5, which is the
                        # exponent in the relation between them:
                        abs_z_error = 0.5 * abs_z * rel_error_res

                        self._z_err[idx_f, ii, jj] = max(
                            MTcc.propagate_error_polar2rect(
                                abs_z,
                                abs_z_error,
                                self.phase[idx_f, ii, jj],
                                self.phase_err[idx_f, ii, jj],
                            )
                        )
        if (
            self.resistivity_model_err is not None
            or self.phase_model_err is not None
        ):
            for idx_f in range(self.frequency.shape[0]):
                for ii in range(2):
                    for jj in range(2):
                        abs_z = np.sqrt(
                            5
                            * self.frequency[idx_f]
                            * self.resistivity[idx_f, ii, jj]
                        )
                        rel_error_res = (
                            self.resistivity_model_err[idx_f, ii, jj]
                            / self.resistivity[idx_f, ii, jj]
                        )
                        # relative error varies by a factor of 0.5, which is the
                        # exponent in the relation between them:
                        abs_z_error = 0.5 * abs_z * rel_error_res

                        self._z_model_err[idx_f, ii, jj] = max(
                            MTcc.propagate_error_polar2rect(
                                abs_z,
                                abs_z_error,
                                self.phase[idx_f, ii, jj],
                                self.phase_model_err[idx_f, ii, jj],
                            )
                        )

    # calculate determinant values
    @property
    def _zdet(self):
        if self._z is not None:
            return np.array([np.linalg.det(zz) ** 0.5 for zz in self._z])

    @property
    def _zdet_var(self):
        if self._z_err is not None:
            return np.array(
                [abs(np.linalg.det(zzv)) ** 0.5 for zzv in self._z_err]
            )
        else:
            return np.ones_like(self._zdet, dtype=np.float)

    @property
    def _zdet_model_var(self):
        if self._z_model_err is not None:
            return np.array(
                [abs(np.linalg.det(zzv)) ** 0.5 for zzv in self._z_model_err]
            )
        else:
            return np.ones_like(self._zdet, dtype=np.float)

    @property
    def phase_det(self):
        if self._zdet is not None:
            return np.rad2deg(np.arctan2(self._zdet.imag, self._zdet.real))

    @property
    def phase_err_det(self):
        if self._zdet is not None:
            return np.rad2deg(np.arcsin(self._zdet_var / abs(self._zdet)))

    @property
    def phase_model_err_det(self):
        if self._zdet is not None:
            return np.rad2deg(np.arcsin(self._zdet_model_var / abs(self._zdet)))

    @property
    def res_det(self):
        if self._zdet is not None:
            return 0.2 * (1.0 / self.frequency) * abs(self._zdet) ** 2

    @property
    def res_err_det(self):
        if self._zdet is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self._zdet + self._zdet_var) ** 2
                - self.res_det
            )

    @property
    def res_model_err_det(self):
        if self._zdet is not None:
            return (
                0.2
                * (1.0 / self.frequency)
                * np.abs(self._zdet + self._zdet_model_var) ** 2
                - self.res_det
            )
