# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 13:20:28 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

from . import Z, Tipper

# =============================================================================


class MTDataFrame:
    def __init__(self, **kwargs):
        self.base_df_dtypes = dict(
            [
                ("station", "U25"),
                ("latitude", float),
                ("longitude", float),
                ("elevation", float),
                ("datum_epsg", "U6"),
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
                ("res_xx", complex),
                ("res_xx_error", float),
                ("res_xx_model_error", float),
                ("res_xy", complex),
                ("res_xy_error", float),
                ("res_xy_model_error", float),
                ("res_yx", complex),
                ("res_yx_error", float),
                ("res_yx_model_error", float),
                ("res_yy", complex),
                ("res_yy_error", float),
                ("res_yy_model_error", float),
                ("phase_xx", complex),
                ("phase_xx_error", float),
                ("phase_xx_model_error", float),
                ("phase_xy", complex),
                ("phase_xy_error", float),
                ("phase_xy_model_error", float),
                ("phase_yx", complex),
                ("phase_yx_error", float),
                ("phase_yx_model_error", float),
                ("phase_yy", complex),
                ("phase_yy_error", float),
                ("phase_yy_model_error", float),
                ("ptxx", complex),
                ("ptxx_error", float),
                ("ptxx_model_error", float),
                ("ptxy", complex),
                ("ptxy_error", float),
                ("ptxy_model_error", float),
                ("ptyx", complex),
                ("ptyx_error", float),
                ("ptyx_model_error", float),
                ("ptyy", complex),
                ("ptyy_error", float),
                ("ptyy_model_error", float),
            ]
        )

        self.df_dtypes = self._get_dtypes(
            [
                "station",
                "latitude",
                "longitude",
                "elevation",
                "datum_epsg",
                "east",
                "north",
                "utm_epsg",
                "model_east",
                "model_north",
                "model_elevation",
                "period",
                "zxx",
                "zxx_error",
                "zxx_model_error",
                "zxy",
                "zxy_error",
                "zxy_model_error",
                "zyx",
                "zyx_error",
                "zyx_model_error",
                "zyy",
                "zyy_error",
                "zyy_model_error",
                "tzx",
                "tzx_error",
                "tzx_model_error",
                "tzy",
                "tzy_error",
                "tzy_model_error",
            ]
        )

        self._index_dict = {
            "xx": {"ii": 0, "jj": 0},
            "xy": {"ii": 0, "jj": 1},
            "yx": {"ii": 1, "jj": 0},
            "yy": {"ii": 1, "jj": 1},
            "zx": {"ii": 0, "jj": 0},
            "zy": {"ii": 0, "jj": 1},
        }

    def _get_dtypes(self, cols):
        """
        Get dtypes for the dataframe

        :param cols: columns to use in data types
        :type cols: list of strings
        :return: dictionary of data types
        :rtype: dictionary

        """

        dtypes = {}
        for key in cols:
            try:
                dtypes[key] = self.base_df_dtypes[key]
            except KeyError:
                raise KeyError(f"Could not find {key} in df_dtypes")

        return dtypes

    def _get_index(self, key):
        """ """

        if key.startswith("z") or key.startswith("t"):
            return self._index_dict[key[1:3]]

        elif key.startswith("res"):
            return self._index_dict[key[4:6]]
        elif key.startswith("phase"):
            return self._index_dict[key[6:8]]
        elif key.startswith("pt"):
            return self._index_dict[key[2:4]]
        else:
            return None

    def make_empty_entry(self, n_entries):
        """

        :param n_entries: number of expected rows
        :type n_entries: int
        :param dtypes: DESCRIPTION
        :type dtypes: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return dict(
            [
                (col, np.zeros(n_entries, dtype))
                for col, dtype in self.df_dtypes.items()
            ]
        )

    def from_z_object(self, z_object, entry):
        """
        Fill impedance
        :param impedance: DESCRIPTION
        :type impedance: TYPE
        :param index: DESCRIPTION
        :type index: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for key in self.df_dtypes.keys():
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                entry[key][:] = z_object.z[:, index["ii"], index["jj"]]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                entry[key][:] = z_object.z_err[:, index["ii"], index["jj"]]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                entry[key][:] = z_object.z_model_err[
                    :, index["ii"], index["jj"]
                ]
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                entry[key][:] = z_object.resistivity[
                    :, index["ii"], index["jj"]
                ]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                entry[key][:] = z_object.resistivity_err[
                    :, index["ii"], index["jj"]
                ]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                entry[key][:] = z_object.resistivity_model_err[
                    :, index["ii"], index["jj"]
                ]

            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                entry[key][:] = z_object.phase[:, index["ii"], index["jj"]]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                entry[key][:] = z_object.phase_err[:, index["ii"], index["jj"]]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                entry[key][:] = z_object.phase_model_err[
                    :, index["ii"], index["jj"]
                ]

        return entry

    def from_t_object(self, t_object, entry):
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
        for key in self.df_dtypes.keys():
            index = self._get_index(key)
            if index is None:
                continue
            if key in ["tzx", "tzy"]:
                entry[key][:] = t_object.tipper[:, index["ii"], index["jj"]]
            elif key in ["tzx_error", "tzy_error"]:
                entry[key][:] = t_object.tipper_err[
                    :, index["ii"], index["jj"]
                ]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                entry[key][:] = t_object.tipper_model_err[
                    :, index["ii"], index["jj"]
                ]

        return entry

    def get_frequency(self, df):
        """
        Get frequencies

        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        try:
            return 1.0 / df.period
        except KeyError:
            return None

    def get_period(self, df):
        """
        Get frequencies

        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        try:
            return df.period
        except KeyError:
            return None

    def to_z_object(self, df):
        """
        fill z_object from dataframe

        Need to have the components this way for transposing the elements so
        that the shape is (nf, 2, 2)
        """

        nf = self.get_frequency(df).size
        z = np.zeros((nf, 2, 2), dtype=complex)
        z_err = np.zeros((nf, 2, 2), dtype=float)
        z_model_err = np.zeros((nf, 2, 2), dtype=float)

        res = np.zeros((nf, 2, 2), dtype=float)
        res_err = np.zeros((nf, 2, 2), dtype=float)
        res_model_err = np.zeros((nf, 2, 2), dtype=float)

        phase = np.zeros((nf, 2, 2), dtype=float)
        phase_err = np.zeros((nf, 2, 2), dtype=float)
        phase_model_err = np.zeros((nf, 2, 2), dtype=float)

        for key in df.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                z[:, index["ii"], index["jj"]] = df[key][:]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                z_err[:, index["ii"], index["jj"]] = df[key][:]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                z_model_err[:, index["ii"], index["jj"]] = df[key][:]

            ### resistivity
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                res[:, index["ii"], index["jj"]] = df[key][:]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                res_err[:, index["ii"], index["jj"]] = df[key][:]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                res_model_err[:, index["ii"], index["jj"]] = df[key][:]

            ### Phase
            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                phase[:, index["ii"], index["jj"]] = df[key][:]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                phase_err[:, index["ii"], index["jj"]] = df[key][:]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                phase_model_err[:, index["ii"], index["jj"]] = df[key][:]

        z_object = Z(z, z_err, self.get_frequency(df), z_model_err)

        if (z == 0).all():
            # only load in resistivity and phase if impedance is 0, otherwise
            # its recreated from z.
            if (res != 0).all():
                if (phase != 0).all():
                    z_object.set_resistivity_phase(
                        res,
                        phase,
                        self.get_frequency(df),
                        res_err=res_err,
                        phase_err=phase_err,
                        res_model_err=res_model_err,
                        phase_model_err=phase_model_err,
                    )
                else:
                    raise ValueError(
                        "cannot estimate Z without phase information"
                    )

        return z_object

    def to_t_object(self, df):
        """
        To a tipper object

        :return: DESCRIPTION
        :rtype: TYPE

        """

        nf = self.get_frequency(df).size
        t = np.zeros((nf, 1, 2), dtype=complex)
        t_err = np.zeros((nf, 1, 2), dtype=float)
        t_model_err = np.zeros((nf, 1, 2), dtype=float)

        for key in df.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["tzx", "tzy"]:
                t[:, index["ii"], index["jj"]] = df[key][:]
            elif key in ["tzx_err", "tzy_err"]:
                t_err[:, index["ii"], index["jj"]] = df[key][:]
            elif key in ["tzx_model_err", "tzy_model_err"]:
                t_model_err[:, index["ii"], index["jj"]] = df[key][:]

        return Tipper(t, t_err, self.get_frequency(df), t_model_err)
