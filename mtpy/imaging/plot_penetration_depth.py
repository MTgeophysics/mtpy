# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:58:58 2022

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata

import mtpy.modeling.occam2d_rewrite as occam2d
from mtpy.imaging.mtplot_tools import PlotBase
from mtpy.utils import MU0

# =============================================================================


class PlotPenetrationDepth(PlotBase):
    """
    Plot the depth of penetration based on the Niblett-Bostick approximation.


    """

    def __init__(self, tf_list, **kwargs):

        self.tf_list = tf_list

        super().__init__(self, **kwargs)

        self.z_mode = "det"
        self.scale_parameter = 2.0 * np.pi * MU0

    def calculate_niblett_bostick_depth(self, resistivity, period):
        """
        Use the Niblett-Bostick approximation for depth of penetration

        :param resistivity: DESCRIPTION
        :type resistivity: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return np.sqrt(resistivity * period / self.scale_parameter)
