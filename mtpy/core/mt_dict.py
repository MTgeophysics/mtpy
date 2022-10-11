# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:58:56 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from collections import OrderedDict

import pandas as pd
import numpy as np

from .mt import MT
from .mt_stations import MTStations

# =============================================================================


class MTDict(OrderedDict, MTStations):
    def __init__(self, tf_list, **kwargs):

        for tf in tf_list:
            tf = self._validate_item(tf)
            self.update(OrderedDict([(tf.station, tf)]))

        MTStations.__init__(self, tf_list, None, **kwargs)

    def _validate_item(self, tf):
        if not isinstance(tf, MT):
            raise TypeError(
                f"entry must be a mtpy.core.MT object not type({type(tf)})"
            )
        return tf

    def to_dataframe(self, utm_crs=None, cols=None):
        """

        :param utm_crs: DESCRIPTION, defaults to None
        :type utm_crs: TYPE, optional
        :param cols: DESCRIPTION, defaults to None
        :type cols: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        df_list = [
            tf.to_dataframe(utm_crs=utm_crs, cols=cols) for tf in self.values()
        ]

        return pd.concat(df_list)

    def from_dataframe(self, df):
        """
        Create an dictionary of MT objects from a dataframe

        :param df: dataframe of mt data
        :type df: `pandas.DataFrame`
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for station in df.station.unique():
            sdf = df.loc[df.station == station]
            mt_object = MT()
            mt_object.from_dataframe(sdf)
            self.update(OrderedDict([(mt_object.station, mt_object)]))

    def interpolate(self, new_periods):
        """
        Interpolate onto common period range

        :param new_periods: DESCRIPTION
        :type new_periods: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for mt_obj in self.values():
            interp_periods = new_periods[
                np.where(
                    (new_periods >= mt_obj.period.max())
                    & (new_periods <= mt_obj.period.min())
                )
            ]

            interp_z, interp_t = mt_obj.interpolate(1.0 / interp_periods)

            mt_obj.Z = interp_z
            mt_obj.Tipper = interp_t
