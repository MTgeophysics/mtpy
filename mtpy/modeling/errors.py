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

    od = np.array([z_array[:, 0, 1], z_array[:, 1, 0]])
    nz = np.nonzeros(od)

    return error_value * np.mean(od[nz])

def compute_median_error(array, error_value)