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


def validate_percent(value):
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


def validate_z_array_shape(array):
    """

    :param array: DESCRIPTION
    :type array: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    if array.shape == (2, 2):
        array = array.reshape((1, 2, 2))

    return array


def mask_zeros(array):
    """
    mask zeros

    :param array: DESCRIPTION
    :type array: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    return np.ma.masked_equal(array, 0)


def compute_percent_error(array, percent):
    """
    Percent error

    :param array: DESCRIPTION
    :type array: TYPE
    :param percent: DESCRIPTION
    :type percent: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    percent = validate_percent(percent)

    return percent * np.abs(array)


def set_floor(error_array, floor):
    """
    Set error floor

    :param array: DESCRIPTION
    :type array: TYPE
    :param floor: DESCRIPTION
    :type floor: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    floor = validate_percent(floor)

    error_array[np.where(error_array < floor)] = floor
    return error_array


def compute_off_diagonal_mean_error(z_array, error_value, floor=True):
    """
    error_value * (Zxy + Zyx) / 2


    :param z_array: DESCRIPTION
    :type z_array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    error_value = validate_percent(error_value)
    z_array = validate_z_array_shape(z_array)

    od = mask_zeros(np.array([z_array[:, 0, 1], z_array[:, 1, 0]]))
    err = error_value * np.ma.abs(np.ma.mean(od, axis=0))

    if floor:
        err = set_floor(err, error_value)

    if isinstance(err, np.ma.core.MaskedArray):
        return err.data

    return err


def compute_median_error(array, error_value, floor=True):
    """
    median(array) * error_value

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    error_value = validate_percent(error_value)
    array = mask_zeros(validate_z_array_shape(array))
    err = np.abs(np.ma.median(array, axis=(1, 2))) * error_value

    if floor:
        err = set_floor(err, error_value)

    if isinstance(err, np.ma.core.MaskedArray):
        return err.data

    return err


def compute_eigen_value_error(array, error_value, floor=True):
    """
    error_value * eigen(array).mean()

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    error_value = validate_percent(error_value)
    array = mask_zeros(validate_z_array_shape(array))

    try:
        err = error_value * np.abs(np.linalg.eigvals(array)).mean(axis=1)
    except np.Exception:
        err = error_value * np.abs(np.linalg.eigvals(array)).mean()

    if np.atleast_1d(err).sum(axis=0) == 0:
        err = error_value * array[np.nonzero(array)].mean()

    if floor:
        err = set_floor(err, error_value)
    return err


def compute_geometric_mean_error(array, error_value, floor=True):
    """
    error_value * sqrt(Zxy * Zyx)

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    error_value = validate_percent(error_value)
    array = validate_z_array_shape(array)

    zero_xy = np.where(array[:, 0, 1] == 0)
    array[zero_xy, 0, 1] = array[zero_xy, 1, 0]

    zero_yx = np.where(array[:, 1, 0] == 0)
    array[zero_yx, 1, 0] = array[zero_yx, 0, 1]

    array = mask_zeros(validate_z_array_shape(array))

    err = error_value * np.ma.sqrt(np.ma.abs(array[:, 0, 1] * array[:, 1, 0]))

    if floor:
        err = set_floor(err, error_value)

    if isinstance(err, np.ma.core.MaskedArray):
        return err.data

    return err


def compute_off_diagonals_error(array, error_value, floor=True):
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

    error_value = validate_percent(error_value)
    z_array = validate_z_array_shape(z_array)

    err_xy = compute_percent_error(array[:, 0, 1], error_value, floor=floor)
    err_yx = compute_percent_error(array[:, 1, 0], error_value, floor=floor)

    err = np.zeros_like(array, dtype=float)
    err[:, 0, :] = err_xy
    err[:, 1, :] = err_yx

    return err


def compute_absolute_error(array, error_value):
    """

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    error_value = validate_percent(error_value)
    err = np.zeros_like(array, dtype=float)
    err[:] = error_value
    return err


# =============================================================================
# Error dictionary
# =============================================================================
ERROR_DICT = {
    "egbert": compute_geometric_mean_error,
    "geometric_mean": compute_geometric_mean_error,
    "arithmetic_mean": compute_off_diagonal_mean_error,
    "off-diagonals": compute_off_diagonals_error,
    "mean_od": compute_off_diagonal_mean_error,
    "median": compute_median_error,
    "eigen": compute_eigen_value_error,
    "percent": compute_percent_error,
    "absolute": compute_absolute_error,
}
