# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase


class PenetrationDepth3D(VisualizationBase):
    @staticmethod
    def get_plot_name():
        return "Penetration Depth (3D)"

    def __init__(self):
        VisualizationBase.__init__(self)
        print "hi there"
