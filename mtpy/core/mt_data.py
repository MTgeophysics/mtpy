# -*- coding: utf-8 -*-
"""
MT Data

Created on Sat Oct  1 17:47:19 2022

:author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
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
        self.tf_list = tf_list
        
        self.data_epsg = None
        self.data_utm_zone = None
        
        self._data_columns = [     
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
             "zxx_model_error",
             "tzy",
             "tzy_error",
             "tzy_model_error",]

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
            
    def _make_empty_entry(self):
        return dict([(col, None) for col in self._data_columns])
    
    def _fill_impedance(self, impedance, impedance_error, index):
        """
        Fill impedance 
        :param impedance: DESCRIPTION
        :type impedance: TYPE
        :param index: DESCRIPTION
        :type index: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        
        z_dict = {"zxx": {"ii": 0, "jj": 0},
                  "zxy": {"ii": 0, "jj": 1},
                  "zyx": {"ii": 1, "jj": 0},
                  "zyy": {"ii": 1, "jj": 1}}
        
        entry = {}
        for z_key, z_index in z_dict.items():
            entry[z_key] = impedance[index, z_index["ii"], z_index[jj]]
            entry[f"{z_key}_error"] = impedance_error[index, z_index["ii"], z_index[jj]]

        return entry
    
    def _fill_tipper(self, tipper, tipper_error, index):
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

        
        t_dict = {"tzx": {"ii": 0, "jj": 0},
                  "tzy": {"ii": 0, "jj": 1},}
        
        entry = {}
        for t_key, t_index in t_dict.items():
            entry[t_key] = tipper[index, t_index["ii"], t_index[jj]]
            entry[f"{t_key}_error"] = tipper_error[index, t_index["ii"], t_index[jj]]

        return entry
        
        
    
    def _fill_entry(self, tf):
        entry = self._make_empty_entry()
        entry["station"] = tf.station
        entry["latitude"] = tf.latitude
        entry["longitude"]
        for index, period in enumerate(tf.period):
            if tf.has_impedance():
                entry["zxx"] = tf.impedance[index, 0, 0]
                
                
        return entry
            
    def _fill_dataframe(self, tf_list):
        """
        Fill the data frame
        
        :param tf_list: DESCRIPTION
        :type tf_list: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        
        entries = []
        
        for tf in tf_list:
            entry = self._make_empty_entry()
            entry["station"] = tf.station
            entry["latitude"] = tf.latitude
            entry["longitude"]
            for ff in tf.period:
                if tf.has_impedance():
                    entry["zxx"] = 
                    
                    
                
                
            

    def data(self):
        pass
        