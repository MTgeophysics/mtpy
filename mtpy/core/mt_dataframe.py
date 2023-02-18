# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 13:20:28 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd

from . import Z, Tipper

# =============================================================================


class MTDataFrame(pd.DataFrame):
    def __init__(
        self,
        n_entries=0,
        impedance=True,
        tipper=True,
        res_phase=False,
        pt=False,
        *args,
        **kwargs,
    ):

        super().__init__(*args, **kwargs)

    def _get_initial_array(
        self,
        n_entries=0,
        impedance=True,
        tipper=True,
        res_phase=False,
        pt=False,
    ):

        base_columns = self._basic_columns
        if impedance:
            base_columns += self._impedance_columns
        if tipper:
            base_columns += self._tipper_columns
        if res_phase:
            base_columns += self._res_phase_columns
        if pt:
            base_columns += self._pt_columns

        return np.empty(n_entries, dtype=self._get_dtypes(base_columns))

    @property
    def _constructor(self):
        """
        Creates a self object that is basically a pandas.Dataframe.
        self is a dataframe-like object inherited from pandas.DataFrame
        self behaves like a dataframe + new custom attributes and methods.
        """
        return MTDataFrame

    @property
    def _constructor_sliced(self):
        return MTDataFrame

    @property
    def _basic_columns(self):
        return [
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
        ]

    @property
    def _impedance_columns(self):
        return [
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
        ]

    @property
    def _tipper_columns(self):
        return [
            "tzx",
            "tzx_error",
            "tzx_model_error",
            "tzy",
            "tzy_error",
            "tzy_model_error",
        ]

    @property
    def _res_phase_columns(self):
        return [
            "res_xx",
            "res_xx_error",
            "res_xx_model_error",
            "res_xy",
            "res_xy_error",
            "res_xy_model_error",
            "res_yx",
            "res_yx_error",
            "res_yx_model_error",
            "res_yy",
            "res_yy_error",
            "res_yy_model_error",
            "phase_xx",
            "phase_xx_error",
            "phase_xx_model_error",
            "phase_xy",
            "phase_xy_error",
            "phase_xy_model_error",
            "phase_yx",
            "phase_yx_error",
            "phase_yx_model_error",
            "phase_yy",
            "phase_yy_error",
            "phase_yy_model_error",
        ]

    @property
    def _pt_columns(self):
        return [
            "ptxx",
            "ptxx_error",
            "ptxx_model_error",
            "ptxy",
            "ptxy_error",
            "ptxy_model_error",
            "ptyx",
            "ptyx_error",
            "ptyx_model_error",
            "ptyy",
            "ptyy_error",
            "ptyy_model_error",
        ]

    @property
    def _base_df_dtypes(self):

        return dict(
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

    @property
    def _index_dict(self):
        return {
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

        dtypes = []
        for key in cols:
            try:
                dtypes.append((key, self._base_df_dtypes[key]))
            except KeyError:
                raise KeyError(f"Could not find {key} in df_dtypes")

        return np.dtype(dtypes)

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
                for col, dtype in self.dtypes.items()
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

        for key in self.dtypes.keys():
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                if z_object._has_tf():
                    entry[key][:] = z_object.z[:, index["ii"], index["jj"]]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                if z_object._has_tf_error():
                    entry[key][:] = z_object.z_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    entry[key][:] = z_object.z_model_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                if z_object._has_tf():
                    entry[key][:] = z_object.resistivity[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                if z_object._has_tf_error():
                    entry[key][:] = z_object.resistivity_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    entry[key][:] = z_object.resistivity_model_error[
                        :, index["ii"], index["jj"]
                    ]

            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                if z_object._has_tf():
                    entry[key][:] = z_object.phase[:, index["ii"], index["jj"]]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                if z_object._has_tf_error():
                    entry[key][:] = z_object.phase_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    entry[key][:] = z_object.phase_model_error[
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
        for key in self.dtypes.keys():
            index = self._get_index(key)
            if index is None:
                continue
            if key in ["tzx", "tzy"]:
                if t_object._has_tf():
                    entry[key][:] = t_object.tipper[:, index["ii"], index["jj"]]
            elif key in ["tzx_error", "tzy_error"]:
                if t_object._has_tf_error():
                    entry[key][:] = t_object.tipper_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                if t_object._has_tf_model_error():
                    entry[key][:] = t_object.tipper_model_error[
                        :, index["ii"], index["jj"]
                    ]

        return entry

    @property
    def frequency(self):
        """
        Get frequencies

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if hasattr(self, "period"):
            return 1.0 / self.period

    def get_station_df(self, station):
        """
        Get station data frame
        :param station: DESCRIPTION
        :type station: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if hasattr(self, "station"):
            if station in self.station:
                return self.loc[self.station == station]

    def to_z_object(self, station):
        """
        fill z_object from dataframe

        Need to have the components this way for transposing the elements so
        that the shape is (nf, 2, 2)
        """

        sdf = self.get_station_df(station)

        nf = self.period.size
        z = np.zeros((nf, 2, 2), dtype=complex)
        z_err = np.zeros((nf, 2, 2), dtype=float)
        z_model_err = np.zeros((nf, 2, 2), dtype=float)

        res = np.zeros((nf, 2, 2), dtype=float)
        res_err = np.zeros((nf, 2, 2), dtype=float)
        res_model_err = np.zeros((nf, 2, 2), dtype=float)

        phase = np.zeros((nf, 2, 2), dtype=float)
        phase_err = np.zeros((nf, 2, 2), dtype=float)
        phase_model_err = np.zeros((nf, 2, 2), dtype=float)

        for key in sdf.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                z[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                z_err[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                z_model_err[:, index["ii"], index["jj"]] = sdf[key][:]

            ### resistivity
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                res[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                res_err[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                res_model_err[:, index["ii"], index["jj"]] = sdf[key][:]

            ### Phase
            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                phase[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                phase_err[:, index["ii"], index["jj"]] = sdf[key][:]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                phase_model_err[:, index["ii"], index["jj"]] = sdf[key][:]

        z_object = Z(z, z_err, self.get_frequency(sdf), z_model_err)

        if (z == 0).all():
            # only load in resistivity and phase if impedance is 0, otherwise
            # its recreated from z.
            if (res != 0).all():
                if (phase != 0).all():
                    z_object.set_resistivity_phase(
                        res,
                        phase,
                        self.frequency,
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
            elif key in ["tzx_error", "tzy_error"]:
                t_err[:, index["ii"], index["jj"]] = df[key][:]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                t_model_err[:, index["ii"], index["jj"]] = df[key][:]

        return Tipper(t, t_err, self.get_frequency(df), t_model_err)
