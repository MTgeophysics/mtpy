# -*- coding: utf-8 -*-
"""
TransferFunction
====================

    * Deals with transfer functions and their covariance if given
    * Z and Tipper will be based on TransferFunction
    
Created on Tue Feb  2 10:49:18 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np

from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================
class TransferFunction:
    """
    Holds transfer function information and their covariances if given
    
    :param tf: Transfer function of shape (num_periods, num_inputs, num_outputs)
    :param sigma_e: Residual covariance matrix with shape 
    (n_periods, n_predicted_channels, n_predicted_channels) namely ex, ey, hz
    :param sigma_s: inverse coherent signal power with shape
    (n_periods, n_predictor_channels, n_predictor_channels) namely hx, hy. 
    For MT should always be a (nf, 2, 2)
    
    """
    
    def __init__(self, tf=None, sigma_s=None, sigma_e=None, periods=None,
                 z_err=None, t_err=None):
        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")
        self._tf = None
        self._sigma_e = None
        self._sigma_s = None
        self._z_err = None
        self._t_err = None
        self._periods = None
        
    @property
    def periods(self):
        return self._periods
    
    @periods.setter
    def periods(self, value):
        value = np.array(value, dtype=np.float)
        if len(value.shape) > 1:
            msg = f"Input must be 1-D, not shape {value.shape}"
            self.logger.error(msg)
            raise ValueError(msg)
        else:
            self._periods = value
            
        
    