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
from mtpy.utils.exceptions import (
    MTpyError_Tipper,
    MTpyError_input_arguments,
)
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================


class Tipper:
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
        self,
        tipper_array=None,
        tipper_err_array=None,
        frequency=None,
        tipper_model_err_array=None,
    ):
        """
        initialize
        """
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")
        self._tipper = tipper_array
        self._tipper_err = tipper_err_array
        self._frequency = frequency
        self._tipper_model_err = tipper_model_err_array

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
            self.logger.error(msg)
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
                self.logger.error(msg)
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

    def _validate_input_array(self, t_array, dtype, old_shape=None):
        """

        :param z_array: DESCRIPTION
        :type z_array: TYPE
        :param dtype: DESCRIPTION
        :type dtype: TYPE
        :param old_shape: DESCRIPTION, defaults to None
        :type old_shape: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if t_array is None:
            return
        if not isinstance(t_array, np.ndarray):
            t_array = np.array(t_array, dtype=dtype)
        if not t_array.dtype in [dtype]:
            t_array = t_array.astype(dtype)
        # check to see if the new tipper array is the same shape as the old
        if self._tipper is not None and self._tipper.shape != t_array.shape:
            msg = (
                "Shape of new array does not match old.  "
                + f"new shape {t_array.shape} != "
                + f"old shape {self._tipper.shape}. "
                + "Make a new Tipper instance to be save."
            )
            self.logger.error(msg)
            raise MTpyError_Tipper(msg)
        if len(t_array.shape) == 3:
            if t_array.shape[1:3] == (1, 2):
                return t_array
            else:
                msg = f"Input array must be shape (n, 1, 2) not {t_array.shape}"
                self.logger.error(msg)
                raise MTpyError_Tipper(msg)
        elif len(t_array.shape) == 2:
            if t_array.shape == (1, 2):
                self.logger.debug(
                    "setting input tipper with shape (1, 2) to (1, 1, 2)"
                )
                return t_array.reshape((1, 1, 2))

            else:
                msg = f"Input array must be shape (n, 1, 2) not {t_array.shape}"
                self.logger.error(msg)
                raise MTpyError_Tipper(msg)
        else:
            msg = f"{t_array.shape} are not the correct dimensions, must be (n, 1, 2)"
            self.logger.error(msg)
            raise MTpyError_Tipper(msg)

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

    # --- tipper ----
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
        old_shape = None
        if self._tipper is not None:
            old_shape = self._tipper.shape
        self._tipper = self._validate_input_array(
            tipper_array, "complex", old_shape
        )

        # neeed to set the rotation angle such that it is an array
        if self.rotation_angle is float:
            self.rotation_angle = np.repeat(
                self.rotation_angle, len(self._tipper)
            )
        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----tipper error---------------
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

        old_shape = None
        if self._tipper is not None:
            old_shape = self._tipper.shape
        self._tipper_err = self._validate_input_array(
            tipper_err_array, "float", old_shape
        )

        # for consistency recalculate mag and angle
        self.compute_mag_direction()

        # for consistency recalculate amplitude and phase
        self.compute_amp_phase()

    # ----tipper model error---------------------------------------------------------
    @property
    def tipper_model_err(self):
        return self._tipper_model_err

    @tipper_model_err.setter
    def tipper_model_err(self, tipper_model_err_array):
        """
        Set the attribute *tipper_err*.

        :param tipper_model_err_array: array of estimated tipper errors
                                 in the shape of [Tx, Ty].
                                 Must be the same shape as tipper_array.
                                 *default* is None
        :type tipper_model_err_array: np.ndarray((nf, 1, 2))
        """

        old_shape = None
        if self._tipper is not None:
            old_shape = self._tipper.shape
        self._tipper_model_err = self._validate_input_array(
            tipper_model_err_array, "float", old_shape
        )

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

        if self.tipper_model_err is not None:
            for idx_f in range(len(self.tipper)):
                for jj in range(2):
                    if self.tipper_model_err is not None:
                        if type(self.tipper) == np.ma.core.MaskedArray:
                            if self.tipper.mask[idx_f, 0, jj]:
                                continue
                        r_err, phi_err = MTcc.propagate_error_rect2polar(
                            np.real(self.tipper[idx_f, 0, jj]),
                            self.tipper_model_err[idx_f, 0, jj],
                            np.imag(self.tipper[idx_f, 0, jj]),
                            self.tipper_model_err[idx_f, 0, jj],
                        )

                        self.amplitude_model_err[idx_f, 0, jj] = r_err
                        self.phase_model_err[idx_f, 0, jj] = phi_err

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
                self.logger.error(
                    'Error - shape of "r" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(r_array.shape), str(self.tipper.shape))
                )
                return
            if self.tipper.shape != phi_array.shape:
                self.logger.error(
                    'Error - shape of "phi" array does not match shape of '
                    + "tipper array: %s ; %s"
                    % (str(phi_array.shape), str(self.tipper.shape))
                )
                return
        else:

            tipper_new = np.zeros(r_array.shape, "complex")

            if r_array.shape != phi_array.shape:
                self.logger.error(
                    'Error - shape of "phi" array does not match shape '
                    + 'of "r" array: %s ; %s'
                    % (str(phi_array.shape), str(r_array.shape))
                )
                return
        # assert real array:
        if np.linalg.norm(np.imag(r_array)) != 0:
            self.logger.error('Error - array "r" is not real valued !')
            return
        if np.linalg.norm(np.imag(phi_array)) != 0:
            self.logger.error('Error - array "phi" is not real valued !')
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

    @property
    def amplitude_model_err(self):
        return self._amplitude_model_err

    @property
    def phase_model_err(self):
        return self._phase_model_err

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

        ## estimate error: THIS MAYBE A HACK
        if self.tipper_model_err is not None:
            self._mag_model_err = np.sqrt(
                self.tipper_model_err[:, 0, 0] ** 2
                + self.tipper_model_err[:, 0, 1] ** 2
            )
            self._angle_model_err = (
                np.rad2deg(
                    np.arctan2(
                        self.tipper_model_err[:, 0, 0],
                        self.tipper_model_err[:, 0, 1],
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

    @property
    def mag_model_err(self):
        return self._mag_model_err

    @property
    def angle_model_err(self):
        return self._angle_model_err

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
            self.logger.error('tipper array is "None" - I cannot rotate that')
            return
        # check for iterable list/set of angles - if so, it must have length 1
        # or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self.logger.error('"Angle" must be a valid number (in degrees)')
                return
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        elif len(alpha) == 1:
            try:
                degreeangle = float(alpha % 360)
            except ValueError:
                self.logger.error('"Angle" must be a valid number (in degrees)')
                return
            # make an n long list of identical angles
            lo_angles = [degreeangle for ii in self.tipper]
        else:
            try:
                lo_angles = [float(ii % 360) for ii in alpha]
            except ValueError:
                self.logger.error('"Angles" must be valid numbers (in degrees)')
                return
        self.rotation_angle = np.array(
            [
                (oldangle + lo_angles[ii]) % 360
                for ii, oldangle in enumerate(self.rotation_angle)
            ]
        )

        if len(lo_angles) != len(self.tipper):
            self.logger.error(
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
