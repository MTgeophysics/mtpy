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
from mtpy.utils.mtpy_logger import get_mtpy_logger
from mt_metadata.transfer_functions.tf import Location

# =============================================================================


class MTDataFrame:
    """
    dataframe for MT data with some convenience properties

    """

    def __init__(self, **kwargs):
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self._data_dtypes = dict(
            [
                ("station", "U25"),
                ("latitude", float),
                ("longitude", float),
                ("elevation", float),
                ("utm_east", float),
                ("utm_north", float),
                ("utm_zone", "U4"),
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

        self._z_object = Z()
        self._z_model_object = Z()
        self._t_object = Tipper()
        self._t_model_object = Tipper()

        self.station = None
        self.latitude = 0
        self.longitude = 0
        self.elevation = 0
        self.utm_east = 0
        self.utm_north = 0
        self.utm_zone = None
        self.model_east = 0
        self.model_north = 0
        self.model_elevation = 0

    def __getattr__(self, name):
        """
        Over loat getattr to make the code more compact

        :param name: DESCRIPTION
        :type name: TYPE
        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if name in [
            "station",
            "latitude",
            "longitude",
            "elevation",
            "utm_east",
            "utm_north",
            "utm_zone",
            "model_east",
            "model_north",
            "model_elevation",
        ]:
            if self.has_data():
                return self.mt_dataframe[name].unique()[0]

        elif name in ["period"]:
            if self.has_data():
                return self.mt_dataframe[name]

        elif name in ["frequency"]:
            if self.has_data():
                return 1.0 / self.mt_dataframe["period"]

        elif name in [
            "zxx",
            "zxy",
            "zyx",
            "zyy",
            "tzx",
            "tzy",
            "zxx_error",
            "zxy_error",
            "zyx_error",
            "zyy_error",
            "tzx_error",
            "tzy_error",
            "zxx_model_error",
            "zxy_model_error",
            "zyx_model_error",
            "zyy_model_error",
            "tzx_model_error",
            "tzy_model_error",
        ]:
            if self.has_data():
                return self.mt_dataframe.loc[:, name]

        elif name in [
            "res_xx",
            "res_xy",
            "res_yx",
            "res_yy",
            "res_xx_error",
            "res_xy_error",
            "res_yx_error",
            "res_yy_error",
            "phase_xx",
            "phase_xy",
            "phase_yx",
            "phase_yy",
            "phase_xx_error",
            "phase_xy_error",
            "phase_yx_error",
            "phase_yy_error",
        ]:
            if self.has_data():
                return getattr(self._z_object, name.replace("error", "err"))

        elif name in [
            "res_xx_model_error",
            "res_xy_model_error",
            "res_yx_model_error",
            "res_yy_model_error",
            "phase_xx_model_error",
            "phase_xy_model_error",
            "phase_yx_model_error",
            "phase_yy_model_error",
        ]:
            if self.has_data():
                return getattr(
                    self._z_model_object, name.replace("model_error", "err")
                )

        else:
            return super().__getattr__(name)

    def __setattr__(self, name, value):
        """
        Over loat setattr to make the code more compact

        :param name: DESCRIPTION
        :type name: TYPE
        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if name in ["station", "utm_zone"]:
            if self.has_data():
                self.mt_dataframe.loc[:, name] = str(value)
            else:
                self.logger.warning(
                    f"Cannot set {name} value of empty dataframe."
                )

        elif name in [
            "latitude",
            "longitude",
            "elevation",
        ]:
            if self.has_data():
                location = Location()
                location.set_attr_from_name(name, value)
                self.mt_dataframe.loc[:, name] = location.get_attr_from_name(
                    name
                )
            else:
                self.logger.warning(
                    f"Cannot set {name} value of empty dataframe."
                )

        elif name in [
            "utm_east",
            "utm_north",
            "utm_zone",
            "model_east",
            "model_north",
            "model_elevation",
        ]:
            if self.has_data():
                self.mt_dataframe.loc[:, name] = float(value)
            else:
                self.logger.warning(
                    f"Cannot set {name} value of empty dataframe."
                )

        elif name in [
            "period",
            "zxx",
            "zxy",
            "zyx",
            "zyy",
            "tzx",
            "tzy",
            "zxx_error",
            "zxy_error",
            "zyx_error",
            "zyy_error",
            "tzx_error",
            "tzy_error",
            "zxx_model_error",
            "zxy_model_error",
            "zyx_model_error",
            "zyy_model_error",
            "tzx_model_error",
            "tzy_model_error",
        ]:
            if self.has_data():
                if isinstance(value, (pd.Series, np.ndarray)):
                    if value.size != self.mt_dataframe.shape[0]:
                        raise ValueError(
                            f"Input must have same size as original {self.mt_dataframe.shape[0]} not {value.size}"
                        )
                    self.mt_dataframe.loc[:, name] = value.astype(
                        self._data_dtypes[name]
                    )

                elif isinstance(value, (float, int, complex)):
                    self.logger.warning(
                        f"input is a number, setting all values in {name} to "
                        f"{value} as type {self._data_dtypes[name]}"
                    )
                    self.mt_dataframe.loc[:, name] = self._data_dtypes[name](
                        value
                    )
                else:
                    raise TypeError(
                        f"Cannot set {name} with type {type(value)}"
                    )
            else:
                self.logger.warning(
                    f"Cannot set {name} value of empty dataframe."
                )

        elif name in [
            "res_xx",
            "res_xy",
            "res_yx",
            "res_yy",
            "res_xx_error",
            "res_xy_error",
            "res_yx_error",
            "res_yy_error",
            "phase_xx",
            "phase_xy",
            "phase_yx",
            "phase_yy",
            "phase_xx_error",
            "phase_xy_error",
            "phase_yx_error",
            "phase_yy_error",
        ]:
            if self.has_data():
                raise AttributeError(f"Cannot currently set {name}")

        elif name in [
            "res_xx_model_error",
            "res_xy_model_error",
            "res_yx_model_error",
            "res_yy_model_error",
            "phase_xx_model_error",
            "phase_xy_model_error",
            "phase_yx_model_error",
            "phase_yy_model_error",
        ]:
            if self.has_data():
                raise AttributeError(f"Cannot currently set {name}")

        else:
            super().__setattr__(name, value)

    def _make_empty_entry(self, n_entries):
        return dict(
            [
                (col, np.zeros(n_entries, dtype))
                for col, dtype in self._data_dtypes.items()
            ]
        )

    def _fill_impedance(self, impedance, impedance_error, entry):
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
            entry[z_key][:] = impedance[
                :, z_index["ii"], z_index["jj"]
            ].to_numpy()
            entry[f"{z_key}_error"][:] = impedance_error[
                :, z_index["ii"], z_index["jj"]
            ].to_numpy()

        return entry

    def _fill_tipper(self, tipper, tipper_error, entry):
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
            entry[t_key][:] = tipper[
                :, t_index["ii"], t_index["jj"]
            ].to_numpy()
            entry[f"{t_key}_error"] = tipper_error[
                :, t_index["ii"], t_index["jj"]
            ].to_numpy()

        return entry

    def from_tf(self, tf):
        """
        fill dataframe from a TF object

        """
        n_entries = tf.period.size
        entry = self._make_empty_entry(n_entries)

        entry["station"][:] = tf.station
        entry["latitude"][:] = tf.latitude
        entry["longitude"][:] = tf.longitude
        entry["elevation"][:] = tf.elevation

        tf.project_to_utm(epsg=self.data_epsg, utm_zone=self.data_utm_zone)
        entry["utm_east"][:] = tf.east
        entry["utm_north"][:] = tf.north
        entry["utm_zone"][:] = tf.utm_zone
        entry["period"][:] = tf.period
        if tf.has_impedance():
            entry = self._fill_impedance(
                tf.impedance, tf.impedance_error, entry
            )
        if tf.has_tipper():
            entry = self._fill_tipper(tf.tipper, tf.tipper_error, entry)

        self._mt_dataframe = pd.DataFrame(entry)

    def has_data(self):
        if self.mt_dataframe.size > 0:
            return True
        return False

    @property
    def mt_dataframe(self):
        """dataframe of data"""
        return self._mt_dataframe

    @mt_dataframe.setter
    def mt_dataframe(self, value):
        """
        Check to make sure the input dataframe is of proper types

        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not isinstance(value, pd.DataFrame):
            msg = f"Input data must be a valid pandas.DataFrame not {type(value)}"
            self.logger.exception(msg)

            raise TypeError(msg)

        empty_df = pd.DataFrame(self._make_empty_entry(0))
        test_dtypes = value.dtypes == empty_df.dtypes

        if not test_dtypes.all():
            for row in test_dtypes[test_dtypes == False].index:
                try:
                    value[row] = value[row].astype(self._data_dtypes[row])
                except ValueError:
                    self.logger.exception(
                        f"{row} cannot be set to {self._data_dtypes[row]}"
                    )

                    raise ValueError(
                        f"{row} cannot be set to {self._data_dtypes[row]}"
                    )

        for key in [
            "station",
            "latitude",
            "longitude",
            "elevation",
            "utm_east",
            "utm_north",
            "utm_zone",
            "model_east",
            "model_north",
            "model_elevation",
        ]:
            if len(value[key].unique()) > 1:
                raise ValueError(
                    f"Input must contain only one {key} value not {len(value.station.unique())}"
                )

        self._mt_dataframe = value

    @property
    def _z_object(self):
        """
        fill z_object from dataframe

        Need to have the components this way for transposing the elements so
        that the shape is (nf, 2, 2)
        """

        if self.has_data():
            z = np.array(
                [
                    [self.mt_dataframe.zxx, self.mt_dataframe.zyx],
                    [self.mt_dataframe.zxy, self.mt_dataframe.zyy],
                ],
                dtype=complex,
            ).T
            z_err = np.array(
                [
                    [self.mt_dataframe.zxx_error, self.mt_dataframe.zyx_error],
                    [self.mt_dataframe.zxy_error, self.mt_dataframe.zyy_error],
                ],
                dtype=float,
            ).T

            return Z(z, z_err, self.frequency)

    @property
    def _z_model_object(self):

        if self.has_data():
            z = np.array(
                [
                    [self.mt_dataframe.zxx, self.mt_dataframe.zyx],
                    [self.mt_dataframe.zxy, self.mt_dataframe.zyy],
                ],
                dtype=complex,
            ).T

            z_model_err = np.array(
                [
                    [
                        self.mt_dataframe.zxx_model_error,
                        self.mt_dataframe.zyx_model_error,
                    ],
                    [
                        self.mt_dataframe.zxy_model_error,
                        self.mt_dataframe.zyy_model_error,
                    ],
                ],
                dtype=float,
            ).T

            return Z(z, z_model_err, self.frequency)

    @property
    def _t_object(self):
        """
        To a tipper object

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.has_data():
            t = np.array(
                [
                    [self.mt_dataframe.tzx],
                    [self.mt_dataframe.tzy],
                ],
                dtype=complex,
            ).T
            t_err = np.array(
                [
                    [self.mt_dataframe.tzx_error],
                    [self.mt_dataframe.tzy_error],
                ],
                dtype=float,
            ).T
            return Tipper(t, t_err, self.frequency)

    @property
    def _t_model_object(self):
        if self.has_data():
            t = np.array(
                [
                    [self.mt_dataframe.tzx],
                    [self.mt_dataframe.tzy],
                ],
                dtype=complex,
            ).T
            t_model_err = np.array(
                [
                    [self.mt_dataframe.tzx_model_error],
                    [self.mt_dataframe.tzy_model_error],
                ],
                dtype=float,
            ).T
            return Tipper(t, t_model_err, self.frequency)

    @property
    def impedance(self):
        if self.has_data():
            return self._z_object.z

    @property
    def impedance_error(self):
        if self.has_data():
            return self._z_object.z_err

    @property
    def impedance_model_error(self):
        if self.has_data():
            return self._z_model_object.z_err
