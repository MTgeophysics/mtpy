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
from mtpy.utils.mtpy_logger import get_mtpy_logger


# ==============================================================================
# Resistivity and phase object
# ==============================================================================
class ResPhase:
    """
    resistivity and phase container with convenience property attributes to
    access the different components.

    All data are stored as private variables _z, _z_err, _z_model_err and
    resistivity and phase attributes are built from them.

    To set resistivity and phase from arrays of resistivity and phase use

    `ResPhase.set_resistivity_phase()`


    """

    def __init__(
        self,
        z=None,
        z_err=None,
        frequency=None,
        z_model_err=None,
        **kwargs,
    ):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self._z = z
        self._z_err = z_err
        self._z_model_err = z_model_err

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

    def _validate_frequency(self, value):
        if value is None:
            return None

        if not isinstance(value, np.ndarray):
            value = np.array(value)

        if len(value.shape) > 1:
            value = value.flatten()

    def _validate_input(self, value):
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
        return self._validate_real_valued(value)

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

    def _get_component(self, comp, array):
        if array is not None:
            index_dict = {"x": 0, "y": 1}
            ii = index_dict[comp[-2]]
            jj = index_dict[comp[-1]]

            return array[:, ii, jj]

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

    @property
    def resistivity(self):
        if self._z is not None:
            return np.apply_along_axis(
                lambda x: np.abs(x) ** 2 / self.frequency * 0.2, 0, self._z
            )

    @property
    def phase(self):
        if self._z is not None:
            return np.rad2deg(np.angle(self._z))

    @property
    def resistivity_err(self):
        if self._z is not None and self._z_err is not None:
            res_err, phi_err = self._compute_rp_error(self._z, self._z_err)
            return res_err

    @property
    def phase_err(self):
        if self._z is not None and self._z_err is not None:
            res_err, phi_err = self._compute_rp_error(self._z, self._z_err)
            return phi_err

    @property
    def resistivity_model_err(self):
        if self._z is not None and self._z_model_err is not None:
            res_err, phi_err = self._compute_rp_error(
                self._z, self._z_model_err
            )
            return res_err

    @property
    def phase_model_err(self):
        if self._z is not None and self._z_err is not None:
            res_err, phi_err = self._compute_rp_error(
                self._z, self._z_model_err
            )
            return phi_err

    def _compute_rp_error(self, z, z_err):
        # calculate resistivity and phase
        res_err = None
        phase_err = None
        if z is not None and z_err is not None:
            res_err = np.zeros_like(self._z, dtype=np.float)
            phase_err = np.zeros_like(self._z, dtype=np.float)

            for idx_f in range(self.frequency.size):
                for ii in range(2):
                    for jj in range(2):
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            z[idx_f, ii, jj].real,
                            z[idx_f, ii, jj].imag,
                            z_err[idx_f, ii, jj],
                        )
                        res_err[idx_f, ii, jj] = r_err
                        phase_err[idx_f, ii, jj] = phi_err
        return res_err, phase_err

    def _compute_z_error(self, res, res_err, phase, phase_err):
        z_err = None
        if res_err is not None or phase_err is not None:
            z_err = np.zeros_like(self._z, dtype=np.float)
            for idx_f in range(self.frequency.shape[0]):
                for ii in range(2):
                    for jj in range(2):
                        abs_z = np.sqrt(
                            5 * self.frequency[idx_f] * res[idx_f, ii, jj]
                        )
                        rel_error_res = (
                            res_err[idx_f, ii, jj] / res[idx_f, ii, jj]
                        )
                        # relative error varies by a factor of 0.5, which is the
                        # exponent in the relation between them:
                        abs_z_error = 0.5 * abs_z * rel_error_res

                        z_err[idx_f, ii, jj] = max(
                            MTcc.propagate_error_polar2rect(
                                abs_z,
                                abs_z_error,
                                phase[idx_f, ii, jj],
                                phase_err[idx_f, ii, jj],
                            )
                        )
        return z_err

    def set_resistivity_phase(
        self,
        resistivity,
        phase,
        frequency,
        res_err=None,
        phase_err=None,
        res_model_err=None,
        phase_model_err=None,
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

        :param res_err: resistivity error array in Ohm-m
        :type res_err: np.ndarray(num_frequency, 2, 2)

        :param phase_err: phase error array in degrees
        :type phase_err: np.ndarray(num_frequency, 2, 2)


        """

        if resistivity is None or phase is None or frequency is None:
            self.logger.debug(
                "Cannot estimate resitivity and phase if resistivity, "
                "phase, or frequency is None."
            )
            return

        self.frequency = self._validate_frequency(frequency)
        resistivity = self._validate_input(resistivity)
        phase = self._validate_input(phase)

        res_err = self._validate_input(res_err)
        phase_err = self._validate_input(phase_err)
        res_model_err = self._validate_input(res_model_err)
        phase_model_err = self._validate_input(phase_model_err)

        abs_z = np.sqrt(5.0 * self.frequency * (resistivity.T)).T
        self._z = abs_z * np.exp(1j * np.radians(phase))

        self._z_err = self._compute_z_error(
            resistivity, res_err, phase, phase_err
        )
        self._z_model_err = self._compute_z_error(
            resistivity, res_model_err, phase, phase_model_err
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
    def res_err_xx(self):
        return self._get_component("xx", self.resistivity_err)

    @property
    def res_err_xy(self):
        return self._get_component("xy", self.resistivity_err)

    @property
    def res_err_yx(self):
        return self._get_component("yx", self.resistivity_err)

    @property
    def res_err_yy(self):
        return self._get_component("yy", self.resistivity_err)

    @property
    def res_model_err_xx(self):
        return self._get_component("xx", self.resistivity_model_err)

    @property
    def res_model_err_xy(self):
        return self._get_component("xy", self.resistivity_model_err)

    @property
    def res_model_err_yx(self):
        return self._get_component("yx", self.resistivity_model_err)

    @property
    def res_model_err_yy(self):
        return self._get_component("yy", self.resistivity_model_err)

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
    def phase_err_xx(self):
        return self._get_component("xx", self.phase_err)

    @property
    def phase_err_xy(self):
        return self._get_component("xy", self.phase_err)

    @property
    def phase_err_yx(self):
        return self._get_component("yx", self.phase_err)

    @property
    def phase_err_yy(self):
        return self._get_component("yy", self.phase_err)

    @property
    def phase_model_err_xx(self):
        return self._get_component("xx", self.phase_model_err)

    @property
    def phase_model_err_xy(self):
        return self._get_component("xy", self.phase_model_err)

    @property
    def phase_model_err_yx(self):
        return self._get_component("yx", self.phase_model_err)

    @property
    def phase_model_err_yy(self):
        return self._get_component("yy", self.phase_model_err)
