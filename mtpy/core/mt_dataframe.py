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


class MTStationDataFrame:
    """
    Dataframe for a single station

    Tried subclassing pandas.DataFrame, but that turned out to not be straight
    forward, so when with compilation instead.

    Think about having period as an index?
    """

    def __init__(self, data=None, n_entries=0, **kwargs):
        self._dtype_list = [
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

        if data is not None:
            self.dataframe = self._validate_data(data)

        else:
            self.dataframe = self._get_initial_df(n_entries)

        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        if self._has_data():
            return self.dataframe.__str__()

        else:
            return "Empty MTStationDataFrame"

    def __repr__(self):
        if self._has_data():
            return self.dataframe.__repr__()
        else:
            return "MTStationDataFrame()"

    def __eq__(self, other):
        other = self._validata_data(other)
        return self.dataframe == other

    def _validate_data(self, data):
        """

        :param data: DESCRIPTION
        :type data: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if isinstance(data, (dict, np.ndarray)):
            df = pd.DataFrame(data)

        elif isinstance(data, pd.DataFrame):
            df = data

        else:
            raise TypeError("Input data must be a pandas.DataFrame")

        for col in ["station", "period", "latitude", "longitude", "elevation"]:
            if col not in df.columns:
                raise ValueError(f"Input missing column {col}.")

        return df

    def _get_initial_df(self, n_entries=0):

        return pd.DataFrame(
            np.empty(n_entries, dtype=np.dtype(self._dtype_list))
        )

    def _has_data(self):
        if self.dataframe is None:
            return False
        elif self.dataframe.shape[0] > 0:
            return True
        return False

    @property
    def size(self):
        if self._has_data():
            return self.period.size

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

    @property
    def period(self):
        """
        Get frequencies

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self._has_data():
            return self.dataframe.period

    @property
    def frequency(self):
        """
        Get frequencies

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self._has_data():
            return 1.0 / self.dataframe.period

    @property
    def station(self):
        """station name"""
        if self._has_data():
            return self.dataframe.station.unique()[0]

    @station.setter
    def station(self, value):
        """station name"""
        if self._has_data():
            self.dataframe.loc[:, "station"] = value

    @property
    def latitude(self):
        """latitude"""
        if self._has_data():
            return self.dataframe.latitude.unique()[0]

    @latitude.setter
    def latitude(self, value):
        """latitude"""
        if self._has_data():
            self.dataframe.loc[:, "latitude"] = value

    @property
    def longitude(self):
        """longitude"""
        if self._has_data():
            return self.dataframe.longitude.unique()[0]

    @longitude.setter
    def longitude(self, value):
        """longitude"""
        if self._has_data():
            self.dataframe.loc[:, "longitude"] = value

    @property
    def elevation(self):
        """elevation"""
        if self._has_data():
            return self.dataframe.elevation.unique()[0]

    @elevation.setter
    def elevation(self, value):
        """elevation"""
        if self._has_data():
            self.dataframe.loc[:, "elevation"] = value

    @property
    def datum_epsg(self):
        """datum_epsg"""
        if self._has_data():
            return self.dataframe.datum_epsg.unique()[0]

    @datum_epsg.setter
    def datum_epsg(self, value):
        """datum_epsg"""
        if self._has_data():
            self.dataframe.loc[:, "datum_epsg"] = value

    @property
    def east(self):
        """station"""
        if self._has_data():
            return self.dataframe.east.unique()[0]

    @east.setter
    def east(self, value):
        """east"""
        if self._has_data():
            self.dataframe.loc[:, "east"] = value

    @property
    def north(self):
        """north"""
        if self._has_data():
            return self.dataframe.north.unique()[0]

    @north.setter
    def north(self, value):
        """north"""
        if self._has_data():
            self.dataframe.loc[:, "north"] = value

    @property
    def utm_epsg(self):
        """utm_epsg"""
        if self._has_data():
            return self.dataframe.utm_epsg.unique()[0]

    @utm_epsg.setter
    def utm_epsg(self, value):
        """utm_epsg"""
        if self._has_data():
            self.dataframe.loc[:, "utm_epsg"] = value

    @property
    def model_east(self):
        """model_east"""
        if self._has_data():
            return self.dataframe.model_east.unique()[0]

    @model_east.setter
    def model_east(self, value):
        """model_east"""
        if self._has_data():
            self.dataframe.loc[:, "model_east"] = value

    @property
    def model_north(self):
        """model_north"""
        if self._has_data():
            return self.dataframe.model_north.unique()[0]

    @model_north.setter
    def model_north(self, value):
        """model_north"""
        if self._has_data():
            self.dataframe.loc[:, "model_north"] = value

    @property
    def model_elevation(self):
        """model_elevation"""
        if self._has_data():
            return self.dataframe.model_elevation.unique()[0]

    @model_elevation.setter
    def model_elevation(self, value):
        """model_elevation"""
        if self._has_data():
            self.dataframe.loc[:, "model_elevation"] = value

    def from_z_object(self, z_object):
        """
        Fill impedance
        :param impedance: DESCRIPTION
        :type impedance: TYPE
        :param index: DESCRIPTION
        :type index: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for key in self.dataframe.dtypes.keys():
            if key in ["period"]:
                self.dataframe.loc[:, "period"] = z_object.period

            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                if z_object._has_tf():
                    self.dataframe.loc[:, key] = z_object.z[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                if z_object._has_tf_error():
                    self.dataframe.loc[:, key] = z_object.z_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    self.dataframe.loc[:, key] = z_object.z_model_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                if z_object._has_tf():
                    self.dataframe.loc[:, key] = z_object.resistivity[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                if z_object._has_tf_error():
                    self.dataframe.loc[:, key] = z_object.resistivity_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    self.dataframe.loc[
                        :, key
                    ] = z_object.resistivity_model_error[
                        :, index["ii"], index["jj"]
                    ]

            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                if z_object._has_tf():
                    self.dataframe.loc[:, key] = z_object.phase[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                if z_object._has_tf_error():
                    self.dataframe.loc[:, key] = z_object.phase_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    self.dataframe.loc[:, key] = z_object.phase_model_error[
                        :, index["ii"], index["jj"]
                    ]

    def from_t_object(self, t_object):
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
        for key in self.dataframe.dtypes.keys():
            if key in ["period"]:
                self.dataframe.loc[:, "period"] = t_object.period

            index = self._get_index(key)
            if index is None:
                continue
            if key in ["tzx", "tzy"]:
                if t_object._has_tf():
                    self.dataframe.loc[:, key] = t_object.tipper[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["tzx_error", "tzy_error"]:
                if t_object._has_tf_error():
                    self.dataframe.loc[:, key] = t_object.tipper_error[
                        :, index["ii"], index["jj"]
                    ]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                if t_object._has_tf_model_error():
                    self.dataframe.loc[:, key] = t_object.tipper_model_error[
                        :, index["ii"], index["jj"]
                    ]

    def to_z_object(self):
        """
        fill z_object from dataframe

        Need to have the components this way for transposing the elements so
        that the shape is (nf, 2, 2)
        """

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

        for key in self.dataframe.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                z[:, index["ii"], index["jj"]] = self.dataframe.loc[:, key]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                z_err[:, index["ii"], index["jj"]] = self.dataframe.loc[:, key]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                z_model_err[:, index["ii"], index["jj"]] = self.dataframe[key][
                    :
                ]

            ### resistivity
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                res[:, index["ii"], index["jj"]] = self.dataframe.loc[:, key]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                res_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    :, key
                ]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                res_model_err[:, index["ii"], index["jj"]] = self.dataframe[
                    key
                ][:]

            ### Phase
            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                phase[:, index["ii"], index["jj"]] = self.dataframe.loc[:, key]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                phase_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    :, key
                ]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                phase_model_err[:, index["ii"], index["jj"]] = self.dataframe[
                    key
                ][:]

        z_object = Z(z, z_err, self.frequency, z_model_err)

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

    def to_t_object(self):
        """
        To a tipper object

        :return: DESCRIPTION
        :rtype: TYPE

        """

        nf = self.dataframe.period.size
        t = np.zeros((nf, 1, 2), dtype=complex)
        t_err = np.zeros((nf, 1, 2), dtype=float)
        t_model_err = np.zeros((nf, 1, 2), dtype=float)

        for key in self.dataframe.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["tzx", "tzy"]:
                t[:, index["ii"], index["jj"]] = self.dataframe[key][:]
            elif key in ["tzx_error", "tzy_error"]:
                t_err[:, index["ii"], index["jj"]] = self.dataframe[key][:]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                t_model_err[:, index["ii"], index["jj"]] = self.dataframe[key][
                    :
                ]

        return Tipper(t, t_err, self.frequency, t_model_err)
