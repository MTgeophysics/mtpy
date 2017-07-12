# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.gui.plot_parameter import FrequencySingle, Ellipse, FrequencyTolerance
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase


class PhaseTensorMap(VisualizationBase):
    @staticmethod
    def plot_description():
        return """
<p>Plot phase tensor map in Lat-Lon Coordinate System.</p>
        """

    def update_ui(self):
        self._frequency_ui.set_data(self._mt_objs)

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
        self._frequency_ui = FrequencySingle(self._parameter_ui)
        self._frequency_ui.setTitle("Frequency (Hz)")
        self._parameter_ui.add_parameter_groubox(self._frequency_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._tolerance_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        # resize
        self._parameter_ui.resize(self._parameter_ui.width(),
                                  self._parameter_ui.sizeHint().height())

        self.update_ui()
