# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtCore

from mtpy.gui.SmartMT.gui.plot_parameter import PlotParameter, ZComponentSingle, FrequencySingle
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.penetration import Depth3D


class PenetrationDepth3D(VisualizationBase):
    def plot(self):
        # get parameters
        try:
            zcomponent = self._z_component_ui.get_selection()
            period = self._frequency_period_ui.get_frequency()
            self._fig = Depth3D(self._mt_objs, period, zcomponent)
            self._fig.plot()
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
        self._frequency_period_ui.setTitle("Frequency Period")
        self._parameter_ui.add_parameter_groubox(self._frequency_period_ui)

        self.update_ui()

# register subclass
VisualizationBase.register(PenetrationDepth3D)
