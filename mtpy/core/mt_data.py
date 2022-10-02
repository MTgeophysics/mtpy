# -*- coding: utf-8 -*-
"""
MT Data

Created on Sat Oct  1 17:47:19 2022

:author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd


from mtpy import MT
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================


class MTData:
    """
    MTData will hold transfer function information for input into models

    The underlying object will be a :class:`pandas.DataFrame` with columns
    for the parameters

     - station
     - latitude
     - longitude
     - elevation
     - utm_east
     - utm_north
     - utm_zone
     - model_east
     - model_north
     - model_elevation
     - period
     - zxx
     - zxx_error
     - zxx_model_error
     - zxy
     - zxy_error
     - zxy_model_error
     - zyx
     - zyx_error
     - zyx_model_error
     - zyy
     - zyy_error
     - zyy_model_error
     - tzx
     - tzx_error
     - zxx_model_error
     - tzy
     - tzy_error
     - tzy_model_error


    properties that can be built from the impedance tensor should be

     - res_xx
     - res_xy
     - res_yy
     - res_yy
     - res_det
     - phase_xx
     - phase_xy
     - phase_yy
     - phase_yy
     - phase_det
     - phase_tensor

    """

    def __init__(self, tf_list, **kwargs):
        self.logger = get_mtpy_logger(
            f"{self.__class__}.{self.__class__.__name__}"
        )

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

        self._tf_dataframe = None

        self.data_epsg = None
        self.data_utm_zone = None

        self.tf_list = tf_list

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def tf_list(self):
        """list of mtpy.core.MT transfer function objects"""
        return self._tf_list

    @tf_list.setter
    def tf_list(self, value):
        """set tf list making sure each element is of proper type"""

        if isinstance(value, MT):
            self._tf_list = [value]
        elif isinstance(value, (list, tuple)):
            tf_list = []
            msg = []
            for ii, tf in enumerate(value):
                if not isinstance(tf, MT):
                    self.logger.error(
                        f"Entry {ii} must be a mtpy.core.MT object not {type(tf)}"
                    )
                else:
                    tf_list.append(tf)

            self._tf_list = tf_list
            if msg != []:
                raise TypeError(
                    "One or more entries in the TF list is not of type mtpy.core.MT"
                )

        else:
            raise TypeError(
                "Input must be a list or tuple of MT objects, or an MT object"
            )

        self._tf_dataframe = self._fill_dataframe(self._tf_list)

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
            entry[t_key][:] = tipper[:, t_index["ii"], t_index["jj"]].to_numpy()
            entry[f"{t_key}_error"] = tipper_error[
                :, t_index["ii"], t_index["jj"]
            ].to_numpy()

        return entry

    def _fill_entry(self, tf):
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

        return pd.DataFrame(entry)

    def _fill_dataframe(self, tf_list):
        """
        Fill the data frame

        :param tf_list: DESCRIPTION
        :type tf_list: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        df_list = []

        for tf in tf_list:
            df_list.append(self._fill_entry(tf))

        return pd.concat(df_list)

    @property
    def data(self):
        """dataframe of data"""
        return self._tf_dataframe

    @data.setter
    def data(self, value):
        """
        Check to make sure the input dataframe is of proper types

        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not isinstance(value, pd.DataFrame):
            msg = (
                f"Input data must be a valid pandas.DataFrame not {type(value)}"
            )
            self.logger.exception(msg)

            raise TypeError(msg)

        empty_df = pd.DataFrame(self._make_empty_entry(0))
        test_dtypes = value.dtypes == empty_df.dtypes

        if not test_dtypes.all():
            for row in test_dtypes[test_dtypes == False].index:
                try:
                    value[row] = value[row].astype(self._data_dtypes[row])
                except ValueError as error:
                    self.logger.exception(
                        f"{row} cannot be set to {self._data_dtypes[row]}"
                    )

                    raise ValueError(
                        f"{row} cannot be set to {self._data_dtypes[row]}"
                    )

        self._tf_dataframe = value
