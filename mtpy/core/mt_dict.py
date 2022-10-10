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

from .mt import MT

# =============================================================================


class MTDict(OrderedDict):
    def __init__(self, tf_list, **kwargs):

        for tf in tf_list:
            tf = self._validate_item(tf)
            self.update(OrderedDict([(tf.station, tf)]))

    def _validate_item(self, tf):
        if not isinstance(tf, MT):
            raise TypeError(
                f"entry must be a mtpy.core.MT object not type({type(tf)})"
            )
        return tf

    def to_dataframe(self):
        """
        create a dataframe from the tf_objects for modeling or plotting

        :return: DESCRIPTION
        :rtype: TYPE

        """

        df_list = [tf.to_dataframe() for tf in self.values()]

        return pd.concat(df_list)
