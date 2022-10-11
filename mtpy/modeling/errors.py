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

    error_array[np.where(error_array) < floor] = floor
    return error_array


def compute_off_diagonal_mean_error(z_array, error_value):
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
    err = error_value * np.abs(np.ma.mean(od, axis=0))

    return err


def compute_median_error(array, error_value):
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
    array = mask_zeros(array)
    err = np.ma.abs(np.ma.median(array)) * error_value

    return err


def compute_eigen_value_error(array, error_value):
    """
    error_value * eigen(array).mean()

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    err = err_value * np.abs(np.linalg.eigvals(array)).mean()
    if np.atleast_1d(err).sum() == 0:
        err = error_value * array[np.nonzero(array)].mean()
    return err


def compute_geometric_mean_error(array, error_value):
    """
    error_value * sqrt(Zxy * Zyx)

    :param array: DESCRIPTION
    :type array: TYPE
    :param error_value: DESCRIPTION
    :type error_value: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    array = validate_z_array_shape(array)
    # if both components masked, then take error floor from
    # max of z_xx or z_yy
    zero_xy = np.where(array[:, 0, 1] == 0)
    array[zero_xy] = array[zero_xy[0], 1, 0]

    zero_yx = np.where(array[:, 1, 0] == 0)
    array[zero_yx] = array[zero_yx[0], 0, 1]

    return err_value * np.sqrt(np.abs(array[:, 0, 1] * array[:, 1, 0]))
