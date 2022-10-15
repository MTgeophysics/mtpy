# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 16:01:37 2022

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np

# =============================================================================


class ModelErrors:
    def __init__(self, array=None, **kwargs):

        self._functions = {
            "egbert": self.compute_geometric_mean_error,
            "geometric_mean": self.compute_geometric_mean_error,
            "arithmetic_mean": self.compute_off_diagonal_mean_error,
            "off-diagonals": self.compute_off_diagonals_error,
            "mean_od": self.compute_off_diagonal_mean_error,
            "median": self.compute_median_error,
            "eigen": self.compute_eigen_value_error,
            "percent": self.compute_percent_error,
            "absolute": self.compute_absolute_error,
            "abs": self.compute_absolute_error,
        }

        self._array_shapes = {
            "impedance": (2, 2),
            "z": (2, 2),
            "transfer_function": (3, 2),
            "tipper": (1, 2),
            "t": (1, 2),
        }

        self.error_value = 5
        self.error_type = "percent"
        self.floor = True
        self.mode = "impedance"
        self.array = array

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        lines = ["Model Errors:", "-" * 20]
        lines += [f"\terror_type:    {self.error_type}"]
        lines += [f"\terror_value:   {self.error_value}"]
        lines += [f"\tfloor:         {self.floor}"]
        lines += [f"\tmode:          {self.mode}"]

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def validate_percent(self, value):
        """
        Make sure the percent is a decimal

        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if value > 1:
            value /= 100.0

        return value

    @property
    def error_parameters(self):
        return {
            "error_value": self.error_value,
            "error_type": self.error_type,
            "floor": self.floor,
        }

    @property
    def error_type(self):
        return self._error_type

    @error_type.setter
    def error_type(self, value):
        if value not in self._functions.keys():
            raise NotImplementedError(f"Error Type {value} not supported.")
        self._error_type = value

    @property
    def floor(self):
        return self._floor

    @floor.setter
    def floor(self, value):
        if value not in [False, True]:
            raise ValueError("Floor must be True or False")
        self._floor = value

    @property
    def error_value(self):
        return self._error_value

    @error_value.setter
    def error_value(self, value):
        self._error_value = self.validate_percent(value)

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        if value not in self._array_shapes.keys():
            raise NotImplementedError(f"Mode {value} not supported.")
        self._mode = value

    def _get_shape(self):
        try:
            return self._array_shapes[self.mode]

        except KeyError:
            raise NotImplementedError(f"Mode {self.mode} not supported.")

    def validate_array_shape(self, array):
        """

        :param array: DESCRIPTION
        :type array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not isinstance(array, np.ndarray):
            array = np.array(array)

        expected_shape = self._get_shape()
        if array.shape == expected_shape:
            array = array.reshape((1, expected_shape[0], expected_shape[1]))

        return array

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, value):
        if value is not None:
            self._array = self.validate_array_shape(value)
        else:
            self._array = None

    def mask_zeros(self, array):
        """
        mask zeros

        :param array: DESCRIPTION
        :type array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return np.ma.masked_equal(array, 0)

    def compute_error(
        self, array=None, error_type=None, error_value=None, floor=None
    ):
        """

        :param array: DESCRIPTION, defaults to None
        :type array: TYPE, optional
        :param error_type: DESCRIPTION, defaults to None
        :type error_type: TYPE, optional
        :param error_value: DESCRIPTION, defaults to None
        :type error_value: TYPE, optional
        :param floor: DESCRIPTION, defaults to None
        :type floor: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if array is not None:
            self.array = array
        if error_type is not None:
            self.error_type = error_type
        if error_value is not None:
            self.error_value = error_value
        if floor is not None:
            self.floor = floor

        return self._functions[self.error_type]()

    def compute_percent_error(self):
        """
        Percent error

        :param array: DESCRIPTION
        :type array: TYPE
        :param percent: DESCRIPTION
        :type percent: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return self.error_value * np.abs(self.array)

    def set_floor(self, error_array):
        """
        Set error floor

        :param array: DESCRIPTION
        :type array: TYPE
        :param floor: DESCRIPTION
        :type floor: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        error_array[np.where(error_array < self.error_value)] = self.error_value
        return error_array

    def compute_off_diagonal_mean_error(self):
        """
        error_value * (Zxy + Zyx) / 2


        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        od = self.mask_zeros(
            np.array([self.array[:, 0, 1], self.array[:, 1, 0]])
        )
        err = self.error_value * np.ma.abs(np.ma.mean(od, axis=0))

        if self.floor:
            err = self.set_floor(err)

        if isinstance(err, np.ma.core.MaskedArray):
            return err.data

        return err

    def compute_median_error(self):
        """
        median(array) * error_value

        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        array = self.mask_zeros(self.array)
        err = np.abs(np.ma.median(array, axis=(1, 2))) * self.error_value

        if self.floor:
            err = self.set_floor(err)

        if isinstance(err, np.ma.core.MaskedArray):
            return err.data

        return err

    def compute_eigen_value_error(self):
        """
        error_value * eigen(array).mean()

        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        array = self.mask_zeros(self.array)

        try:
            err = self.error_value * np.abs(np.linalg.eigvals(array)).mean(
                axis=1
            )
        except np.Exception:
            err = self.error_value * np.abs(np.linalg.eigvals(array)).mean()

        if np.atleast_1d(err).sum(axis=0) == 0:
            err = self.error_value * array[np.nonzero(array)].mean()

        if self.floor:
            err = self.set_floor(err)
        return err

    def compute_geometric_mean_error(self):
        """
        error_value * sqrt(Zxy * Zyx)

        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        array = self.array.copy()

        zero_xy = np.where(array[:, 0, 1] == 0)
        array[zero_xy, 0, 1] = array[zero_xy, 1, 0]

        zero_yx = np.where(array[:, 1, 0] == 0)
        array[zero_yx, 1, 0] = array[zero_yx, 0, 1]

        array = self.mask_zeros(array)

        err = self.error_value * np.ma.sqrt(
            np.ma.abs(array[:, 0, 1] * array[:, 1, 0])
        )

        if self.floor:
            err = self.set_floor(err)

        if isinstance(err, np.ma.core.MaskedArray):
            return err.data

        return err

    def compute_off_diagonals_error(self):
        """
        set zxx and zxy the same error and zyy and zyx the same error

        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :param floor: DESCRIPTION, defaults to True
        :type floor: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        err_xy = self.compute_percent_error(self.array[:, 0, 1])
        err_yx = self.compute_percent_error(self.array[:, 1, 0])

        err = np.zeros_like(self.array, dtype=float)
        err[:, 0, :] = err_xy
        err[:, 1, :] = err_yx

        return err

    def compute_absolute_error(self):
        """

        :param array: DESCRIPTION
        :type array: TYPE
        :param error_value: DESCRIPTION
        :type error_value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        err = np.zeros_like(self.array, dtype=float)
        err[:] = self.error_value
        return err
