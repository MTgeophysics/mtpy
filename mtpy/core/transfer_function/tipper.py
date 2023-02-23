# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 23:25:57 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import cmath
import copy
import math

import numpy as np

import mtpy.utils.calculator as MTcc
from .base import TFBase

# =============================================================================


class Tipper(TFBase):
    """
    Tipper class --> generates a Tipper-object.

    Errors are given as standard deviations (sqrt(VAR))

    :param tipper: tipper array in the shape of [Tx, Ty]
                         *default* is None
    :type tipper: np.ndarray((nf, 1, 2), dtype='complex')


    :param tipper_error: array of estimated tipper errors
                             in the shape of [Tx, Ty].
                             Must be the same shape as tipper.
                             *default* is None
    :type tipper_error: np.ndarray((nf, 1, 2))


    :param frequency: array of frequencyuencies corresponding to the tipper elements.
                 Must be same length as tipper.
                 *default* is None
    :type frequency: np.ndarray(nf)


    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    frequency            array of frequencyuencies corresponding to elements of z
    rotation_angle  angle of which data is rotated by

    tipper          tipper array
    tipper_error       tipper error array
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
        self,
        tipper=None,
        tipper_error=None,
        frequency=None,
        tipper_model_error=None,
    ):
        """
        initialize
        """
        super().__init__(
            tf=tipper,
            tf_error=tipper_error,
            tf_model_error=tipper_model_error,
            frequency=frequency,
            _name="tipper",
            _expected_shape=(1, 2),
            inputs=["x", "y"],
            outputs=["z"],
        )

    # --- tipper ----
    @property
    def tipper(self):
        if self._has_tf():
            return self._dataset.transfer_function.values

    @tipper.setter
    def tipper(self, tipper):
        """
        Set the attribute *tipper*

        :param tipper: tipper array in the shape of [Tx, Ty]
                             *default* is None
        :type tipper: np.ndarray((nf, 1, 2), dtype='complex')
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

        tipper = self._validate_array_input(tipper, "complex", old_shape)
        if tipper is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf=tipper)
        else:
            self._dataset["transfer_function"].loc[self.comps] = tipper

    # ----tipper error---------------
    @property
    def tipper_error(self):
        if self._has_tf_error():
            return self._dataset.transfer_function_error.values

    @tipper_error.setter
    def tipper_error(self, tipper_error):
        """
        Set the attribute *tipper_error*.

        :param tipper_error: array of estimated tipper errors
                                 in the shape of [Tx, Ty].
                                 Must be the same shape as tipper.
                                 *default* is None
        :type tipper_error: np.ndarray((nf, 1, 2))
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

        tipper_error = self._validate_array_input(
            tipper_error, "float", old_shape
        )
        if tipper_error is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf_error=tipper_error)
        else:
            self._dataset["transfer_function_error"].loc[
                self.comps
            ] = tipper_error

    # ----tipper model error---------------------------------------------------------
    @property
    def tipper_model_error(self):
        if self._has_tf_model_error():
            return self._dataset.transfer_function_model_error.values

    @tipper_model_error.setter
    def tipper_model_error(self, tipper_model_error):
        """
        Set the attribute *tipper_model_error*.

        :param tipper_model_error: array of estimated tipper errors
                                 in the shape of [Tx, Ty].
                                 Must be the same shape as tipper.
                                 *default* is None
        :type tipper_model_error: np.ndarray((nf, 1, 2))
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
        tipper_model_error = self._validate_array_input(
            tipper_model_error, "float", old_shape
        )
        if tipper_model_error is None:
            return

        if self._is_empty():
            self._dataset = self._initialize(tf_error=tipper_model_error)
        else:
            self._dataset["transfer_function_model_error"].loc[
                self.comps
            ] = tipper_model_error

    # ----amplitude and phase
    def _compute_amp_phase_error(self, error):
        """
        Sets attributes:
                        * *amplitude*
                        * *phase*
                        * *amplitude_error*
                        * *phase_error*

        values for resistivity are in in Ohm m and phase in degrees.
        """
        amplitude_error = None
        phase_error = None
        if error is not None:
            amplitude_error = np.zeros(self.tipper_error.shape)
            phase_error = np.zeros(self.tipper_error.shape)
            for idx_f in range(len(self.tipper)):
                for jj in range(2):
                    if type(self.tipper) == np.ma.core.MaskedArray:
                        if self.tipper.mask[idx_f, 0, jj]:
                            continue
                    r_error, phi_error = MTcc.propagate_error_rect2polar(
                        np.real(self.tipper[idx_f, 0, jj]),
                        error[idx_f, 0, jj],
                        np.imag(self.tipper[idx_f, 0, jj]),
                        error[idx_f, 0, jj],
                    )

                    amplitude_error[idx_f, 0, jj] = r_error
                    phase_error[idx_f, 0, jj] = phi_error

        return amplitude_error, phase_error

    def set_amp_phase(self, r, phi):
        """
        Set values for amplitude(r) and argument (phi - in degrees).

        Updates the attributes:
                        * tipper
                        * tipper_error

        """

        if self.tipper is not None:

            tipper_new = copy.copy(self.tipper)

            if self.tipper.shape != r.shape:
                self.logger.error(
                    'Error - shape of "r" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(r.shape), str(self.tipper.shape))
                )
                return
            if self.tipper.shape != phi.shape:
                self.logger.error(
                    'Error - shape of "phi" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(phi.shape), str(self.tipper.shape))
                )
                return
        else:

            tipper_new = np.zeros(r.shape, "complex")

            if r.shape != phi.shape:
                self.logger.error(
                    'Error - shape of "phi" array does not match shape '
                    + 'of "r" array: %s ; %s' % (str(phi.shape), str(r.shape))
                )
                return
        # assert real array:
        if np.linalg.norm(np.imag(r)) != 0:
            self.logger.error('Error - array "r" is not real valued !')
            return
        if np.linalg.norm(np.imag(phi)) != 0:
            self.logger.error('Error - array "phi" is not real valued !')
            return
        for idx_f in range(len(r)):
            for jj in range(2):
                tipper_new[idx_f, 0, jj] = cmath.rect(
                    r[idx_f, 0, jj],
                    math.radians(phi[idx_f, 0, jj]),
                )
        self.tipper = tipper_new

    # ---------------------------------
    # properties
    @property
    def amplitude(self):
        if self._has_tf():
            return np.abs(self.tipper)

    @property
    def phase(self):
        if self._has_tf():
            return np.rad2deg(np.angle(self.tipper))

    @property
    def amplitude_error(self):
        if self._has_tf_error():
            return self._compute_amp_phase_error(self.tipper_error)[0]

    @property
    def phase_error(self):
        if self._has_tf_error():
            return self._compute_amp_phase_error(self.tipper_error)[1]

    @property
    def amplitude_model_error(self):
        if self._has_tf_model_error():
            return self._compute_amp_phase_error(self.tipper_model_error)[0]

    @property
    def phase_model_error(self):
        if self._has_tf_model_error():
            return self._compute_amp_phase_error(self.tipper_model_error)[1]

    # ----magnitude and direction----------------------------------------------

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
        if self._has_tf():
            return np.sqrt(
                self.tipper[:, 0, 0].real ** 2 + self.tipper[:, 0, 1].real ** 2
            )

    @property
    def mag_imag(self):
        if self._has_tf():
            return np.sqrt(
                self.tipper[:, 0, 0].imag ** 2 + self.tipper[:, 0, 1].imag ** 2
            )

    @property
    def angle_real(self):
        if self._has_tf():
            return np.rad2deg(
                np.arctan2(
                    self.tipper[:, 0, 1].real, self.tipper[:, 0, 0].real
                )
            )

    @property
    def angle_imag(self):
        if self._has_tf():
            return np.rad2deg(
                np.arctan2(
                    self.tipper[:, 0, 1].imag, self.tipper[:, 0, 0].imag
                )
            )

    @property
    def mag_error(self):
        if self._has_tf_error():
            return np.sqrt(
                self.tipper_error[:, 0, 0] ** 2
                + self.tipper_error[:, 0, 1] ** 2
            )

    @property
    def angle_error(self):
        if self._has_tf_error():
            return np.abs(
                np.rad2deg(
                    np.arctan(
                        self.tipper_error[:, 0, 0] / self.tipper_error[:, 0, 1]
                    )
                )
                - 45
            )

    @property
    def mag_model_error(self):
        if self._has_tf_model_error():
            return np.sqrt(
                self.tipper_model_error[:, 0, 0] ** 2
                + self.tipper_model_error[:, 0, 1] ** 2
            )

    @property
    def angle_model_error(self):
        if self._has_tf_model_error():
            return np.abs(
                np.rad2deg(
                    np.arctan(
                        self.tipper_model_error[:, 0, 0]
                        / self.tipper_model_error[:, 0, 1]
                    )
                )
                - 45
            )
