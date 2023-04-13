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


class MTDataFrame:
    """
    Dataframe for a single station

    Tried subclassing pandas.DataFrame, but that turned out to not be straight
    forward, so when with compilation instead.

    Think about having period as an index?
    """

    def __init__(self, data=None, n_entries=0, **kwargs):
        self._dtype_list = [
            ("survey", "U25"),
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
            ("profile_offset", float),
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
            ("res_xx", float),
            ("res_xx_error", float),
            ("res_xx_model_error", float),
            ("res_xy", float),
            ("res_xy_error", float),
            ("res_xy_model_error", float),
            ("res_yx", float),
            ("res_yx_error", float),
            ("res_yx_model_error", float),
            ("res_yy", float),
            ("res_yy_error", float),
            ("res_yy_model_error", float),
            ("phase_xx", float),
            ("phase_xx_error", float),
            ("phase_xx_model_error", float),
            ("phase_xy", float),
            ("phase_xy_error", float),
            ("phase_xy_model_error", float),
            ("phase_yx", float),
            ("phase_yx_error", float),
            ("phase_yx_model_error", float),
            ("phase_yy", float),
            ("phase_yy_error", float),
            ("phase_yy_model_error", float),
            ("ptxx", float),
            ("ptxx_error", float),
            ("ptxx_model_error", float),
            ("ptxy", float),
            ("ptxy_error", float),
            ("ptxy_model_error", float),
            ("ptyx", float),
            ("ptyx_error", float),
            ("ptyx_model_error", float),
            ("ptyy", float),
            ("ptyy_error", float),
            ("ptyy_model_error", float),
            ("rms_zxx", float),
            ("rms_zxy", float),
            ("rms_zyx", float),
            ("rms_zyy", float),
            ("rms_tzx", float),
            ("rms_tzy", float),
        ]

        if data is not None:
            self.dataframe = self._validate_data(data)

        else:
            self.dataframe = self._get_initial_df(n_entries)

        self.working_survey = None
        self.working_station = None

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

    @property
    def _column_names(self):
        return [col[0] for col in self._dtype_list]

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

        if data is None:
            return

        if isinstance(data, (dict, np.ndarray, pd.DataFrame)):
            df = pd.DataFrame(data)

        elif isinstance(data, (MTDataFrame)):
            df = data.dataframe

        else:
            raise TypeError(
                f"Input data must be a pandas.DataFrame not {type(data)}"
            )

        for col in self._dtype_list:
            if col[0] not in df.columns:

                df[col[0]] = np.zeros(df.shape[0], dtype=col[1])

        # resort to the desired column order
        if df.columns.to_list() != self._column_names:
            df = df[self._column_names]

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

    def get_station_df(self, station=None):
        """
        get a single station df

        :return: DESCRIPTION
        :rtype: TYPE

        """
        if station is not None:
            self.working_station = station
        if self._has_data():
            if self.working_station is None:
                self.working_station = self.dataframe.station.unique()[0]

            if self.working_station not in self.dataframe.station.values:
                raise ValueError(
                    f"Could not find station {self.working_station} in dataframe."
                )

            return self.dataframe[
                self.dataframe.station == self.working_station
            ]

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
            return np.sort(self.dataframe.period.unique())

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
    def survey(self):
        """survey name"""
        if self._has_data():
            if self.working_survey is None:
                self.working_survey = self.dataframe.survey.unique()[0]
            return self.working_survey

    @survey.setter
    def survey(self, value):
        """survey name"""
        if self._has_data():
            if self.working_survey in [None, ""]:
                self.dataframe.loc[
                    self.dataframe.survey == "", "survey"
                ] = value
                self.working_survey = value

    @property
    def station(self):
        """station name"""
        if self._has_data():
            if self.working_station is None:
                self.working_station = self.dataframe.station.unique()[0]
            return self.working_station

    @station.setter
    def station(self, value):
        """station name"""
        if self._has_data():
            if self.working_station in [None, ""]:
                self.dataframe.loc[
                    self.dataframe.station == "", "station"
                ] = value
                self.working_station = value

    @property
    def latitude(self):
        """latitude"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "latitude"
            ].unique()[0]

    @latitude.setter
    def latitude(self, value):
        """latitude"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "latitude"
            ] = value

    @property
    def longitude(self):
        """longitude"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "longitude"
            ].unique()[0]

    @longitude.setter
    def longitude(self, value):
        """longitude"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "longitude"
            ] = value

    @property
    def elevation(self):
        """elevation"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "elevation"
            ].unique()[0]

    @elevation.setter
    def elevation(self, value):
        """elevation"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "elevation"
            ] = value

    @property
    def datum_epsg(self):
        """datum_epsg"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "datum_epsg"
            ].unique()[0]

    @datum_epsg.setter
    def datum_epsg(self, value):
        """datum_epsg"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "datum_epsg"
            ] = value

    @property
    def east(self):
        """station"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "east"
            ].unique()[0]

    @east.setter
    def east(self, value):
        """east"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "east"
            ] = value

    @property
    def north(self):
        """north"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "north"
            ].unique()[0]

    @north.setter
    def north(self, value):
        """north"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "north"
            ] = value

    @property
    def utm_epsg(self):
        """utm_epsg"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "utm_epsg"
            ].unique()[0]

    @utm_epsg.setter
    def utm_epsg(self, value):
        """utm_epsg"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "utm_epsg"
            ] = value

    @property
    def model_east(self):
        """model_east"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "model_east"
            ].unique()[0]

    @model_east.setter
    def model_east(self, value):
        """model_east"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "model_east"
            ] = value

    @property
    def model_north(self):
        """model_north"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "model_north"
            ].unique()[0]

    @model_north.setter
    def model_north(self, value):
        """model_north"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station, "model_north"
            ] = value

    @property
    def model_elevation(self):
        """model_elevation"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "model_elevation"
            ].unique()[0]

    @model_elevation.setter
    def model_elevation(self, value):
        """model_elevation"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station,
                "model_elevation",
            ] = value

    @property
    def profile_offset(self):
        """profile_offset"""
        if self._has_data():
            return self.dataframe.loc[
                self.dataframe.station == self.station, "profile_offset"
            ].unique()[0]

    @profile_offset.setter
    def profile_offset(self, value):
        """profile_offset"""
        if self._has_data():
            self.dataframe.loc[
                self.dataframe.station == self.station,
                "profile_offset",
            ] = value

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
                self.dataframe.loc[
                    self.dataframe.station == self.station, "period"
                ] = z_object.period

            index = self._get_index(key)
            if index is None:
                continue

            if key in ["zxx", "zxy", "zyx", "zyy"]:
                if z_object._has_tf():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.z[:, index["ii"], index["jj"]]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                if z_object._has_tf_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.z_error[:, index["ii"], index["jj"]]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.z_model_error[:, index["ii"], index["jj"]]
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                if z_object._has_tf():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.resistivity[:, index["ii"], index["jj"]]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                if z_object._has_tf_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.resistivity_error[:, index["ii"], index["jj"]]
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
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.phase[:, index["ii"], index["jj"]]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                if z_object._has_tf_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.phase_error[:, index["ii"], index["jj"]]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                if z_object._has_tf_model_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = z_object.phase_model_error[:, index["ii"], index["jj"]]

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
                self.dataframe.loc[
                    self.dataframe.station == self.station, "period"
                ] = t_object.period

            index = self._get_index(key)
            if index is None:
                continue
            if key in ["tzx", "tzy"]:
                if t_object._has_tf():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = t_object.tipper[:, index["ii"], index["jj"]]
            elif key in ["tzx_error", "tzy_error"]:
                if t_object._has_tf_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = t_object.tipper_error[:, index["ii"], index["jj"]]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                if t_object._has_tf_model_error():
                    self.dataframe.loc[
                        self.dataframe.station == self.station, key
                    ] = t_object.tipper_model_error[:, index["ii"], index["jj"]]

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
                z[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in ["zxx_error", "zxy_error", "zyx_error", "zyy_error"]:
                z_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in [
                "zxx_model_error",
                "zxy_model_error",
                "zyx_model_error",
                "zyy_model_error",
            ]:
                z_model_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]

            ### resistivity
            elif key in ["res_xx", "res_xy", "res_yx", "res_yy"]:
                res[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in [
                "res_xx_error",
                "res_xy_error",
                "res_yx_error",
                "res_yy_error",
            ]:
                res_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in [
                "res_xx_model_error",
                "res_xy_model_error",
                "res_yx_model_error",
                "res_yy_model_error",
            ]:
                res_model_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]

            ### Phase
            elif key in ["phase_xx", "phase_xy", "phase_yx", "phase_yy"]:
                phase[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in [
                "phase_xx_error",
                "phase_xy_error",
                "phase_yx_error",
                "phase_yy_error",
            ]:
                phase_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in [
                "phase_xx_model_error",
                "phase_xy_model_error",
                "phase_yx_model_error",
                "phase_yy_model_error",
            ]:
                phase_model_err[
                    :, index["ii"], index["jj"]
                ] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]

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

        nf = self.period.size
        t = np.zeros((nf, 1, 2), dtype=complex)
        t_err = np.zeros((nf, 1, 2), dtype=float)
        t_model_err = np.zeros((nf, 1, 2), dtype=float)

        for key in self.dataframe.columns:
            index = self._get_index(key)
            if index is None:
                continue

            if key in ["tzx", "tzy"]:
                t[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in ["tzx_error", "tzy_error"]:
                t_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]
            elif key in ["tzx_model_error", "tzy_model_error"]:
                t_model_err[:, index["ii"], index["jj"]] = self.dataframe.loc[
                    self.dataframe.station == self.station, key
                ]

        return Tipper(t, t_err, self.frequency, t_model_err)
