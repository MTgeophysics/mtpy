# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui
import numpy as np


from mtpy.gui.SmartMT.visualization.plot_parameter import PlotParameter
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase


class PenetrationDepth3D(VisualizationBase):
    def update_ui(self):
        self.update_periods()

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
        # setup periods
        self.update_ui()

    def update_periods(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))
            all_periods = 1.0/np.array(all_freqs)
            # sort all frequencies in ascending order
            all_unique_periods = sorted(list(set(all_periods)))
            self._parameter_ui.slider_min = np.min(all_unique_periods)
            self._parameter_ui.slider_max = np.max(all_unique_periods)
            self._parameter_ui.slider_tick_size = (self._parameter_ui.slider_max -
                                                   self._parameter_ui.slider_min) / 100.0  # 100 is the default ticks of the slidere
            # self._parameter_ui.period_histogram.set_data(all_periods)
            for period in all_unique_periods:
                self._parameter_ui.ui.comboBoxPeriod.addItem("%.5f" % period)

# register subclass
VisualizationBase.register(PenetrationDepth3D)
