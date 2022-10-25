#!/usr/bin/env python

"""
.. module:: TFBase
   :synopsis: Generic Transfer Function object

.. moduleauthor:: Jared Peacock <jpeacock@usgs.gov> 

Updated 11/2020 for logging and formating (J. Peacock).
    - ToDo: add functionality for covariance matrix
"""

# =============================================================================
# Imports
# =============================================================================
import copy

import numpy as np
import xarray as xr

import mtpy.utils.calculator as MTcc
from mtpy.utils.mtpy_logger import get_mtpy_logger

# ==============================================================================
# Impedance Tensor Class
# ==============================================================================
class TFBase:
    """ """

    def __init__(
        self,
        tf=None,
        tf_error=None,
        frequency=None,
        tf_model_error=None,
        tf_dataset=None,
    ):
        """ """
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")
        self.rotation_angle = 0.0
        self.inputs = ["x"]
        self.outputs = ["y"]
        self._expected_shape = (1, 1)
        self._name = "base transfer function"

        frequency = self._validate_frequency(frequency)

        self._dataset = self._initialize(
            periods=1.0 / frequency,
            tf=tf,
            tf_error=tf_error,
            tf_model_error=tf_model_error,
        )

    def __str__(self):
        lines = [f"Transfer Function {self._name}", "-" * 30]
        if self.frequency is not None:
            lines.append(f"\tNumber of periods:  {self.frequency.size}")
            lines.append(
                f"\tFrequency range:        {self.frequency.min():.5E} -- "
                f"{self.frequency.max():.5E} Hz"
            )
            lines.append(
                f"\tPeriod range:           {1/self.frequency.max():.5E} -- "
                f"{1/self.frequency.min():.5E} s"
            )
            lines.append("")
            lines.append(f"\tHas {self._name}:              {self._has_tf()}")
            lines.append(
                f"\tHas {self._name}_error:        {self._has_tf_error()}"
            )
            lines.append(
                f"\tHas {self._name}_model_error:  {self._has_tf_model_error()}"
            )
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, TFBase):
            msg = f"Cannot compare {type(other)} with TFBase"
            self.logger.error(msg)
            raise ValueError(msg)
        if self._tf_dataset != other._tf_dataset:
            return False
        return True

    def copy(self):
        return copy.deepcopy(self)

    def _initialize(
        self, periods=[1], tf=[0 + 0j], tf_error=[0], tf_model_error=[0]
    ):
        """
        initialized based on input channels, output channels and period

        """

        if tf is not None:
            tf = self._validate_array_input(tf, complex)
            periods = self._validate_frequency(periods, tf.shape[0])
            if tf_error is not None:
                self._validate_array_shape(tf_error, tf.shape)
            else:
                tf_error = np.zeros_like(tf, dtype=float)

            if tf_model_error is not None:
                self._validate_array_shape(tf_model_error, tf.shape)
            else:
                tf_model_error = np.zeros_like(tf, dtype=float)

        elif tf_error is not None:
            tf_error = self._validate_array_input(tf_error, float)
            periods = self._validate_frequency(periods, tf_error.shape[0])
            tf = np.zeros_like(tf_error, dtype=complex)

            if tf_model_error is not None:
                self._validate_array_shape(tf_model_error, tf_error.shape)
            else:
                tf_model_error = np.zeros_like(tf_error, dtype=float)

        elif tf_model_error is not None:
            tf_model_error = self._validate_array_input(tf_model_error, float)
            tf = np.zeros_like(tf_model_error, dtype=complex)
            tf_error = np.zeros_like(tf_model_error, dtype=float)
            periods = self._validate_frequency(periods, tf_model_error.shape[0])

        else:
            periods = self._validate_frequency(periods)

        tf = xr.DataArray(
            data=tf,
            dims=["period", "output", "input"],
            coords={
                "period": periods,
                "output": self.outputs,
                "input": self.inputs,
            },
            name="transfer_function",
        )
        tf_err = xr.DataArray(
            data=tf_error,
            dims=["period", "output", "input"],
            coords={
                "period": periods,
                "output": self.outputs,
                "input": self.inputs,
            },
            name="transfer_function_error",
        )
        tf_model_err = xr.DataArray(
            data=tf_model_error,
            dims=["period", "output", "input"],
            coords={
                "period": periods,
                "output": self.outputs,
                "input": self.inputs,
            },
            name="transfer_function_model_error",
        )

        return xr.Dataset(
            {
                tf.name: tf,
                tf_err.name: tf_err,
                tf_model_err.name: tf_model_err,
            }
        )

    def _is_empty(self):
        """
        Check to see if the data set is empty, default settings

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if (
            (self._dataset.transfer_function.values == 0).all()
            and (self._dataset.transfer_function_error.values == 0).all()
            and (self._dataset.transfer_function_model_error.values == 0).all()
        ):
            return True
        return False

    def _has_tf(self):
        if not self._is_empty():
            return not (self._dataset.transfer_function.values == 0).all()
        return False

    def _has_tf_error(self):
        if not self._is_empty():
            return not (self._dataset.transfer_function_error.values == 0).all()
        return False

    def _has_tf_model_error(self):
        if not self._is_empty():
            return not (
                self._dataset.transfer_function_model_error.values == 0
            ).all()
        return False

    # ---frequencyuency-------------------------------------------------------------
    @property
    def frequency(self):
        """
        frequencyuencies for each impedance tensor element

        Units are Hz.
        """
        return 1.0 / self._dataset.period.values

    @frequency.setter
    def frequency(self, frequency):
        """
        Set the array of frequency.

        :param frequency: array of frequencyunecies (Hz)
        :type frequency: np.ndarray
        """

        if frequency is None:
            return

        if self._is_empty():
            frequency = self._validate_frequency(frequency)
        else:
            frequency = self._validate_frequency(
                frequency, n_frequencies=self._dataset.period.shape[0]
            )

        self._dataset = self._dataset.assign_coords({"period": 1.0 / frequency})

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

    @property
    def n_periods(self):
        return self.period.size

    def _validate_frequency(self, frequency, n_frequencies=None):
        """
        validate frequency
        """

        if frequency is None:
            return np.array([1])

        frequency = np.array(frequency, dtype=float)
        if len(frequency) > 1:
            frequency = frequency.flatten()

        if n_frequencies is not None:
            if frequency.size == 1:
                if (frequency == np.array([1])).all():
                    return np.arange(1, n_frequencies + 1, 1)
            if frequency.size != n_frequencies:
                raise ValueError(
                    f"input frequencies must have shape {n_frequencies} not "
                    f"{frequency.size}. "
                    "Use tf._dataset = TFBase._initialize(1./new_frequencies) "
                    "or make a new transfer function object"
                )

        return frequency

    def _validate_array_input(self, tf_array, expected_dtype, old_shape=None):
        """
        Validate an input impedance array

        :param array: DESCRIPTION
        :type array: TYPE
        :param dtype: DESCRIPTION
        :type dtype: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if tf_array is None:
            return
        if not isinstance(tf_array, np.ndarray):
            if isinstance(tf_array, (float, int, complex)):
                tf_array = [tf_array]
            tf_array = np.array(tf_array, dtype=expected_dtype)
        if tf_array.dtype not in [expected_dtype]:
            tf_array = tf_array.astype(expected_dtype)

        if len(tf_array.shape) == 3:
            if tf_array.shape[1:3] == self._expected_shape:
                return tf_array
            else:
                msg = (
                    f"Input array must be shape (n, "
                    f"{self.expected_shape[0]}, {self.expected_shape[1]}) "
                    f"not {tf_array.shape}"
                )
                self.logger.error(msg)
                raise ValueError(msg)
        elif len(tf_array.shape) == 2:
            if tf_array.shape == self._expected_shape:
                return tf_array.reshape(
                    (1, self._expected_shape[0], self._expected_shape[1])
                )
                self.logger.debug(
                    f"setting input tf with shape {self._expected_shape} "
                    f"to (1, self._expected_shape[0], self._expected_shape[1])"
                )
            else:
                msg = (
                    f"Input array must be shape (n, "
                    f"{self._expected_shape[0]}, {self._expected_shape[1]}) "
                    f"not {tf_array.shape}"
                )
                self.logger.error(msg)
                raise ValueError(msg)
        else:
            msg = (
                f"{tf_array.shape} are not the correct dimensions, "
                f"must be (n, {self._expected_shape[0]}, {self._expected_shape[1]})"
            )
            self.logger.error(msg)
            raise ValueError(msg)

        if old_shape is not None:
            self._validate_array_shape(tf_array, old_shape)

    def _validate_array_shape(self, array, expected_shape):
        """
        Check array for expected shape

        :param array: DESCRIPTION
        :type array: TYPE
        :param expected_shape: DESCRIPTION
        :type expected_shape: TYPE
        :raises ValueError: DESCRIPTION
        :return: DESCRIPTION
        :rtype: TYPE

        """
        # check to see if the new z array is the same shape as the old
        if array.shape != expected_shape:
            msg = (
                "Shape of array does not match expected.  "
                f"array shape {array.shape} != "
                f"expected shape {expected_shape}."
            )
            self.logger.error(msg)
            raise ValueError(msg)

    @property
    def inverse(self):
        """
        Return the inverse of transfer function.

        (no error propagtaion included yet)

        """

        if self.has_tf():
            inverse = self._dataset.copy()

            try:
                inverse.transfer_function = np.linalg.inv(
                    inverse.transfer_function
                )

            except np.linalg.LinAlgError:
                raise ValueError(
                    "Transfer Function is a singular matrix cannot invert"
                )

            return inverse

    def rotate(self, alpha, inplace=False):
        """
        Rotate transfer function array by angle alpha.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North, positive in clockwise direction.
        (Mathematically negative!)

        In non-rotated state, X refs to North and Y to East direction.

        """

        if not self._has_tf():
            self.logger.warning(
                "transfer function array is empty and cannot be rotated"
            )
            return

        def _validate_angle(self, angle):
            """validate angle to be a valid float"""
            try:
                return float(angle % 360)
            except ValueError:
                msg = f"Angle must be a valid number (in degrees) not {alpha}"
                self.logger.error(msg)
                raise ValueError(msg)

        if isinstance(alpha, (float, int, str)):
            degree_angle = _validate_angle(self, alpha)

        elif isinstance(alpha, (list, tuple, np.ndarray)):
            if len(alpha) == 1:
                degree_angle = _validate_angle(self, alpha[0])
            else:
                degree_angle = np.array(alpha, dtype=float) % 360
                if degree_angle.size != self.n_periods:
                    raise ValueError(
                        "angles must be the same size as periods "
                        f"{self.n_periods} not {degree_angle.size}"
                    )

        def _rotation_matrix(angle):
            phi = np.deg2rad(angle)
            cphi = np.cos(phi)
            sphi = np.sin(phi)

            return np.array([[cphi, -sphi], [sphi, cphi]])

        self.rotation_angle = self.rotation_angle + degree_angle

        ds = self._dataset.copy()

        if np.iterate(degree_angle):
            for index, angle in enumerate(degree_angle):
                r = _rotation_matrix(angle)
                if self._has_tf():
                    ds.transfer_function[index].values.dot(r)
                if self._has_tf_error():
                    ds.transfer_function_error[index].values.dot(r)
                if self._has_tf():
                    ds.transfer_function_model_error[index].values.dot(r)

        else:
            r = _rotation_matrix(degree_angle)
            if self._has_tf():
                ds.transfer_function.values.dot(r)
            if self._has_tf_error():
                ds.transfer_function_error.values.dot(r)
            if self._has_tf():
                ds.transfer_function_model_error.values.dot(r)

        if inplace:
            self._dataset = ds
        else:
            return ds

    def interpolate(self, new_periods):
        """
        interpolate onto a new period range

        :param new_periods: DESCRIPTION
        :type new_periods: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        pass