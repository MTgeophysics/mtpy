# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase


class PhaseTensorMap(VisualizationBase):
    @staticmethod
    def plot_description():
        return """
<p>Plot phase tensor map in Lat-Lon Coordinate System.</p>
        """

    def update_ui(self):
        pass

    def parameter_ui(self):
        pass

    @staticmethod
    def plot_name():
        return "Phase Tensor Map"

    def plot(self):
        # get data
        # NOTE: this is a hack because the existing bug(s) in the PlotPhaseTensorMaps class, that is the
        # constructor of the class only populate all the necessary information correctly when reads from file
        # this is the only way before this bug(s) is fixed
        file_list = []
        for mt_obj in self._mt_objs:
            file_list.append(mt_obj.fn)

        raise NotImplemented()

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # set up ui
