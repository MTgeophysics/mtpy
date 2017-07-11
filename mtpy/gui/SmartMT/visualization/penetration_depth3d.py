# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtCore

from mtpy.gui.SmartMT.visualization.plot_parameter import PlotParameter
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.penetration import Depth3D


class PenetrationDepth3D(VisualizationBase):
    def plot(self):
        # get parameters
        zcomponent = None
        if self.parameter_ui.ui.radioButton_det.isChecked():
            zcomponent = 'det'
        elif self.parameter_ui.ui.radioButton_zxy.isChecked():
            zcomponent = 'zxy'
        elif self.parameter_ui.ui.radioButton_zyx.isChecked():
            zcomponent = 'zyx'
        try:
            period = float(self.parameter_ui.ui.comboBoxPeriod.currentText())
            self._fig = Depth3D(self._mt_objs, period, zcomponent)
            self._fig.plot()
        except Exception as e:
            self._logger.warning(e.message)

    def update_ui(self):
        self._parameter_ui.set_data(self._mt_objs)
        self._parameter_ui.update_ui()

    @property
    def parameter_ui(self):
        return self._parameter_ui

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
        self._parameter_ui = PlotParameter(self._parent)
        # set to use single z-component group
        self._parameter_ui.enable_single_z_component_selection()
        # connect plot button
        QtCore.QObject.connect(self._parameter_ui.ui.pushButtonPlot, QtCore.SIGNAL("clicked()"), self.show_figure)
        # setup periods
        self.update_ui()

# register subclass
VisualizationBase.register(PenetrationDepth3D)
