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

# =============================================================================


class MTDataFrame:
    """
    dataframe for MT data with some convenience properties

    """

    def __init__(self, **kwargs):

        self._data_dtypes = [
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

        self.data_epsg = None
        self.data_utm_zone = None

        self._mt_dataframe = pd.DataFrame(self._make_empty_entry(0))
        self.station = None

    def _make_empty_entry(self, n_entries):
        return dict(
            [
                (col, np.zeros(n_entries, dtype))
                for col, dtype in self._data_dtypes
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

    def _set_component(self, component, value):
        """
        Set a component of the dataframe

        :param component: DESCRIPTION
        :type component: TYPE
        :param value: DESCRIPTION
        :type value: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
