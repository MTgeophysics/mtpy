# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 13:20:28 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import pandas as pd
import numpy as np

from .z import Z, Tipper

# =============================================================================


df_data_dtypes = dict(
    [
        ("station", "U25"),
        ("latitude", float),
        ("longitude", float),
        ("elevation", float),
        ("east", float),
        ("north", float),
        ("utm_epsg", "U6"),
        ("model_east", float),
        ("model_north", float),
        ("model_elevation", float),
        ("period", float),
        ("zxx", complex),
        ("zxx_error", float),
        ("zxx_model_error", float),
        ("zxy", complex),
        ("zxy_error", float),
        ("zxy_model_error", float),
        ("zyx", complex),
        ("zyx_error", float),
        ("zyx_model_error", float),
        ("zyy", complex),
        ("zyy_error", float),
        ("zyy_model_error", float),
        ("tzx", complex),
        ("tzx_error", float),
        ("tzx_model_error", float),
        ("tzy", complex),
        ("tzy_error", float),
        ("tzy_model_error", float),
    ]
)


def make_empty_entry(n_entries):
    return dict(
        [
            (col, np.zeros(n_entries, dtype))
            for col, dtype in df_data_dtypes.items()
        ]
    )


def fill_impedance(impedance, impedance_error, entry):
    """
    Fill impedance
    :param impedance: DESCRIPTION
    :type impedance: TYPE
    :param index: DESCRIPTION
    :type index: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    z_dict = {
        "zxx": {"ii": 0, "jj": 0},
        "zxy": {"ii": 0, "jj": 1},
        "zyx": {"ii": 1, "jj": 0},
        "zyy": {"ii": 1, "jj": 1},
    }

    for z_key, z_index in z_dict.items():
        entry[z_key][:] = impedance[:, z_index["ii"], z_index["jj"]].to_numpy()
        entry[f"{z_key}_error"][:] = impedance_error[
            :, z_index["ii"], z_index["jj"]
        ].to_numpy()

    return entry


def fill_tipper(tipper, tipper_error, entry):
    """
    Fill tipper
    :param tipper: DESCRIPTION
    :type tipper: TYPE
    :param tipper_error: DESCRIPTION
    :type tipper_error: TYPE
    :param index: DESCRIPTION
    :type index: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    t_dict = {
        "tzx": {"ii": 0, "jj": 0},
        "tzy": {"ii": 0, "jj": 1},
    }

    for t_key, t_index in t_dict.items():
        entry[t_key][:] = tipper[:, t_index["ii"], t_index["jj"]].to_numpy()
        entry[f"{t_key}_error"] = tipper_error[
            :, t_index["ii"], t_index["jj"]
        ].to_numpy()

    return entry


def get_frequency(df):
    """
    Get frequencies

    :param df: DESCRIPTION
    :type df: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """

    return 1.0 / df.period


def to_z_object(df):
    """
    fill z_object from dataframe

    Need to have the components this way for transposing the elements so
    that the shape is (nf, 2, 2)
    """

    z = np.array(
        [
            [df.zxx, df.zyx],
            [df.zxy, df.zyy],
        ],
        dtype=complex,
    ).T
    z_err = np.array(
        [
            [df.zxx_error, df.zyx_error],
            [df.zxy_error, df.zyy_error],
        ],
        dtype=float,
    ).T

    return Z(z, z_err, get_frequency(df))


def to_z_model_object(df):
    """
    fill model z object

    :param df: DESCRIPTION
    :type df: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    z = np.array(
        [
            [df.zxx, df.zyx],
            [df.zxy, df.zyy],
        ],
        dtype=complex,
    ).T

    z_model_err = np.array(
        [
            [
                df.zxx_model_error,
                df.zyx_model_error,
            ],
            [
                df.zxy_model_error,
                df.zyy_model_error,
            ],
        ],
        dtype=float,
    ).T

    return Z(z, z_model_err, get_frequency(df))


def to_t_object(df):
    """
    To a tipper object

    :return: DESCRIPTION
    :rtype: TYPE

    """
    t = np.array(
        [
            [df.tzx],
            [df.tzy],
        ],
        dtype=complex,
    ).T
    t_err = np.array(
        [
            [df.tzx_error],
            [df.tzy_error],
        ],
        dtype=float,
    ).T
    return Tipper(t, t_err, get_frequency(df))


def to_t_model_object(df):
    """
    to model tipper object

    :param df: DESCRIPTION
    :type df: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    t = np.array(
        [
            [df.tzx],
            [df.tzy],
        ],
        dtype=complex,
    ).T
    t_model_err = np.array(
        [
            [df.tzx_model_error],
            [df.tzy_model_error],
        ],
        dtype=float,
    ).T
    return Tipper(t, t_model_err, get_frequency(df))
