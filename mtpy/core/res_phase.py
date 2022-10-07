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
        freq=None,
        z_model_err=None,
        **kwargs,
    ):
        self._logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self._z = z_array
        self._z_err = z_err_array
        self._z_model_err = z_model_err

        self._resistivity = None
        self._phase = None

        self._resistivity_err = None
        self._phase_err = None

        self._freq = freq

        for key in kwargs:
            setattr(self, key, kwargs[key])

    def __str__(self):
        lines = ["Resistivity and Phase", "-" * 30]
        if self.freq is not None:
            lines.append(f"\tNumber of frequencies:  {self.freq.size}")
            lines.append(
                f"\tFrequency range:        {self.freq.min():.5E} -- {self.freq.max():.5E} Hz"
            )
            lines.append(
                f"\tPeriod range:           {1/self.freq.max():.5E} -- {1/self.freq.min():.5E} s"
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
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.resistivity[:, ii, jj]

        elif name in ["res_err_xx", "res_err_xy", "res_err_yx", "res_err_yy"]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.resistivity_err[:, ii, jj]

        elif name in [
            "res_model_err_xx",
            "res_model_err_xy",
            "res_model_err_yx",
            "res_model_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.resistivity_model_err[:, ii, jj]

        elif name in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.phase[:, ii, jj]

        elif name in [
            "phase_err_xx",
            "phase_err_xy",
            "phase_err_yx",
            "phase_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.phase_err[:, ii, jj]

        elif name in [
            "phase_model_err_xx",
            "phase_model_err_xy",
            "phase_model_err_yx",
            "phase_model_err_yy",
        ]:
            ii = index_dict[name[-2]]
            jj = index_dict[name[-1]]
            return self.phase_model_err[:, ii, jj]

        else:
            return super().__getattr__(name)

    @property
    def freq(self):
        return self._freq

    @freq.setter
    def freq(self, value):
        self._freq = value

    @property
    def resistivity(self):
        return self._resistivity

    @resistivity.setter
    def resistivity(self, res_array):
        self._resistivity = res_array

    @property
    def resistivity_err(self):
        return self._resistivity_err

    @resistivity_err.setter
    def resistivity_err(self, res_err_array):
        self._resistivity_err = res_err_array

    @property
    def phase(self):
        return self._phase

    @phase.setter
    def phase(self, phase_array):
        self._phase = phase_array

    @property
    def phase_err(self):
        return self._phase_err

    @phase_err.setter
    def phase_err(self, phase_err_array):
        self._phase_err = phase_err_array

    def compute_resistivity_phase(
        self, z_array=None, z_err_array=None, freq=None
    ):
        """
        compute resistivity and phase from z and z_err
        """

        if z_array is not None:
            self._z = z_array
        if z_err_array is not None:
            self._z_err = z_err_array
        if freq is not None:
            self.freq = freq
        # The _z_err can be None!!!
        if self._z is None or self.freq is None:
            if self._z is None:
                msg = "z values are None, cannot compute parameters"
            elif self._freq is None:
                msg = "freq values are None, cannot compute parameters"
            self._logger.error(msg)
            raise MTpyError_Z(msg)
        self._resistivity = np.apply_along_axis(
            lambda x: np.abs(x) ** 2 / self.freq * 0.2, 0, self._z
        )
        self._phase = np.rad2deg(np.angle(self._z))

        self._resistivity_err = np.zeros_like(self._resistivity, dtype=np.float)
        self._phase_err = np.zeros_like(self._phase, dtype=np.float)

        # calculate resistivity and phase
        if self._z_err is not None:
            for idx_f in range(self.freq.size):
                for ii in range(2):
                    for jj in range(2):
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            self._z[idx_f, ii, jj].real,
                            self._z[idx_f, ii, jj].imag,
                            self._z_err[idx_f, ii, jj],
                        )
                        self._resistivity_err[idx_f, ii, jj] = (
                            self._resistivity[idx_f, ii, jj] * r_err
                        )

                        self._phase_err[idx_f, ii, jj] = phi_err

        if self._z_model_err is not None:
            for idx_f in range(self.freq.size):
                for ii in range(2):
                    for jj in range(2):
                        r_err, phi_err = MTcc.z_error2r_phi_error(
                            self._z[idx_f, ii, jj].real,
                            self._z[idx_f, ii, jj].imag,
                            self._z_model_err[idx_f, ii, jj],
                        )
                        self._resistivity_model_err[idx_f, ii, jj] = (
                            self._resistivity[idx_f, ii, jj] * r_err
                        )

                        self._phase_model_err[idx_f, ii, jj] = phi_err

    def set_res_phase(
        self,
        res_array,
        phase_array,
        freq,
        res_err_array=None,
        phase_err_array=None,
        res_model_err=None,
        phase_model_err=None,
    ):
        """
        Set values for resistivity (res - in Ohm m) and phase
        (phase - in degrees), including error propagation.


        :param res_array: resistivity array in Ohm-m
        :type res_array: np.ndarray(num_freq, 2, 2)

        :param phase_array: phase array in degrees
        :type phase_array: np.ndarray(num_freq, 2, 2)

        :param freq: frequency array in Hz
        :type freq: np.ndarray(num_freq)

        :param res_err_array: resistivity error array in Ohm-m
        :type res_err_array: np.ndarray(num_freq, 2, 2)

        :param phase_err_array: phase error array in degrees
        :type phase_err_array: np.ndarray(num_freq, 2, 2)


        """

        self._logger.debug("Resetting z and z_err")

        self._resistivity = res_array
        self._phase = phase_array
        self.freq = freq
        self._resistivity_err = res_err_array
        self._phase_err = phase_err_array
        self._resistivity_model_err = res_model_err
        self._phase_model_err = phase_model_err

        # assert real array:
        if np.linalg.norm(np.imag(res_array)) != 0:
            msg = "Resistivity is not real valued"
            self._logger.error(msg)
            raise MTpyError_input_arguments(msg)
        if np.linalg.norm(np.imag(phase_array)) != 0:
            msg = "Phase is not real valued"
            self._logger.error(msg)
            raise MTpyError_input_arguments(msg)
        abs_z = np.sqrt(5.0 * self.freq * (self.resistivity.T)).T
        self._z = abs_z * np.exp(1j * np.radians(self.phase))

        self._z_err = np.zeros_like(self._z, dtype=np.float)
        # ---------------------------
        # error propagation:
        if self._resistivity_err is not None or self._phase_err is not None:
            for idx_f in range(self.freq.shape[0]):
                for ii in range(2):
                    for jj in range(2):
                        abs_z = np.sqrt(
                            5
                            * self.freq[idx_f]
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
            self._resistivity_model_err is not None
            or self._phase_model_err is not None
        ):
            for idx_f in range(self.freq.shape[0]):
                for ii in range(2):
                    for jj in range(2):
                        abs_z = np.sqrt(
                            5
                            * self.freq[idx_f]
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

    # @property
    # def res_xx(self):
    #     return self._resistivity[:, 0, 0]

    # @property
    # def res_xy(self):
    #     return self._resistivity[:, 0, 1]

    # @property
    # def res_yx(self):
    #     return self._resistivity[:, 1, 0]

    # @property
    # def res_yy(self):
    #     return self._resistivity[:, 1, 1]

    # @property
    # def phase_xx(self):
    #     return self._phase[:, 0, 0]

    # @property
    # def phase_xy(self):
    #     return self._phase[:, 0, 1]

    # @property
    # def phase_yx(self):
    #     return self._phase[:, 1, 0]

    # @property
    # def phase_yy(self):
    #     return self._phase[:, 1, 1]

    # @property
    # def res_err_xx(self):
    #     return self._resistivity_err[:, 0, 0]

    # @property
    # def res_err_xy(self):
    #     return self._resistivity_err[:, 0, 1]

    # @property
    # def res_err_yx(self):
    #     return self._resistivity_err[:, 1, 0]

    # @property
    # def res_err_yy(self):
    #     return self._resistivity_err[:, 1, 1]

    # @property
    # def phase_err_xx(self):
    #     return self._phase_err[:, 0, 0]

    # @property
    # def phase_err_xy(self):
    #     return self._phase_err[:, 0, 1]

    # @property
    # def phase_err_yx(self):
    #     return self._phase_err[:, 1, 0]

    # @property
    # def phase_err_yy(self):
    #     return self._phase_err[:, 1, 1]

    # calculate determinant values
    @property
    def _zdet(self):
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
    def phase_det(self):
        return np.arctan2(self._zdet.imag, self._zdet.real) * (180 / np.pi)

    @property
    def phase_err_det(self):
        return np.arcsin(self._zdet_var / abs(self._zdet)) * (180 / np.pi)

    @property
    def res_det(self):
        return 0.2 * (1.0 / self.freq) * abs(self._zdet) ** 2

    @property
    def res_err_det(self):
        return (
            0.2 * (1.0 / self.freq) * np.abs(self._zdet + self._zdet_var) ** 2
            - self.res_det
        )
