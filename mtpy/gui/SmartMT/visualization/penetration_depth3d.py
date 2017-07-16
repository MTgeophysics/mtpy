# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtCore

from mtpy.gui.SmartMT.gui.plot_parameter import PlotParameter, ZComponentSingle, FrequencySingle, FrequencyTolerance
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.penetration import Depth3D


class PenetrationDepth3D(VisualizationBase):
    def get_parameter_str(self):
        return "z-component=%s, period=%.5f, tolerance=%.2f%%" % (self._zcomponent, self._period, self._tolerance*100)

    def plot(self):
        # get parameters
        try:
            self._zcomponent = self._z_component_ui.get_selection()
            self._period = self._frequency_period_ui.get_frequency()
            self._tolerance = self._tolerance_ui.get_tolerance_in_float()
            self._plotting_object = Depth3D(self._mt_objs, self._period, self._zcomponent, self._tolerance)
            self._plotting_object.plot()
            self._fig = self._plotting_object.get_figure()

        except Exception as e:
            self._logger.warning(e.message)

    def update_ui(self):
        self._frequency_period_ui.set_data(self._mt_objs)

    @staticmethod
    def plot_description():
        return """
<p>For a batch of MT_stations (input as a list of MT objects):</p>

<ul>
<li>plot the Penetration Depth vs the station_location, for a given period value or index (1/freq)</li>
</ul>

<p><strong>Note:</strong> that the values of periods within10% tolerance (ptol=0.1) are considered as equal. Setting a smaller value for ptol(=0.05)may result less MT sites data included.</p>
"""

    @staticmethod
    def plot_name():
        return "Penetration Depth (3D)"

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # add parameter sub component
        self._z_component_ui = ZComponentSingle(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._z_component_ui)

        self._frequency_period_ui = FrequencySingle(self._parameter_ui, unit="seconds", distribution="Period", inverse=True)
        self._frequency_period_ui.setTitle("Frequency Period (seconds)")
        self._parameter_ui.add_parameter_groubox(self._frequency_period_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._tolerance_ui.setTitle("Period Tolerance")
        self._parameter_ui.add_parameter_groubox(self._tolerance_ui)

        # resize
        # self._parameter_ui.resize(self._parameter_ui.width(),
        #                           self._parameter_ui.sizeHint().height())

        self.update_ui()

        # plot parameters
        self._zcomponent = None
        self._period = None
        self._tolerance = None

# register subclass
VisualizationBase.register(PenetrationDepth3D)
