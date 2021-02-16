# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""

import mtpy.imaging.penetration
from mtpy.gui.SmartMT.Components.PlotParameter import (
    ZComponentMultiple,
    ZComponentSingle,
    FrequencyTolerance,
    FrequencySelection,
    FrequencyIndex,
    StationSelection,
    ZUnit,
)
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase


class Depth1D(VisualizationBase):
    def update_ui(self):
        self._station_ui.set_data(self._mt_objs)

    @staticmethod
    def plot_name():
        return "Penetration Depth (1D)"

    def get_plot_tooltip(self):
        return "station=%s, rho=%s" % (self._station.station, str(self._rhos))

    @staticmethod
    def plot_description():
        return """
<p>plot the penetration depth vs all the periods (1/frequency) of 1 station.</p>
<p><strong>Note:</strong> This only plot one station per plot.</p>
        """

    def plot(self):
        # get parameters
        self._rhos = self._z_component_ui.get_selection()
        self._station = self._station_ui.get_station()
        self._plotting_object = mtpy.imaging.penetration.Depth1D(
            self._station, self._rhos
        )
        self._plotting_object.plot()
        self._fig = self._plotting_object.get_figure()

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        self._rhos = None
        self._station = None
        # add parameter gui sub components
        self._station_ui = StationSelection(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._station_ui)
        self._z_component_ui = ZComponentMultiple(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._z_component_ui)

        self._parameter_ui.end_of_parameter_components()


class Depth2D(VisualizationBase):
    def plot(self):
        # get parameters
        self._rho = self._z_component_ui.get_selection()
        self._period_index = self._frequency_period_ui.get_index_list()
        self._plotting_object = mtpy.imaging.penetration.Depth2D(
            self._mt_objs, self._period_index, self._rho
        )
        self._plotting_object.plot()
        self._fig = self._plotting_object.get_figure()

    @staticmethod
    def plot_name():
        return "Penetration Depth (2D)"

    @staticmethod
    def plot_description():
        return """
<p>plot the Penetration Depth profile at the given periods vs the stations locations.</p>
<p><strong>Note:</strong> This plot requires the identical list of frequencies from all stations.</p>
        """

    def update_ui(self):
        self._frequency_period_ui.set_data(self._mt_objs)

    def get_plot_tooltip(self):
        return "rho=%s, period_index=%s" % (self._rho, self._period_index)

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        self._rho = None
        self._period_index = None
        # add parameter sub component
        self._z_component_ui = ZComponentSingle(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._z_component_ui)

        self._frequency_period_ui = FrequencyIndex(self._parameter_ui, use_period=True)
        self._parameter_ui.add_parameter_groupbox(self._frequency_period_ui)

        self._parameter_ui.end_of_parameter_components()


class Depth3D(VisualizationBase):
    def get_plot_tooltip(self):
        return "z-component=%s, period=%.5f, tolerance=%.2f%%, z_unit=%s" % (
            self._zcomponent,
            self._period,
            self._tolerance * 100,
            self._z_unit,
        )

    def plot(self):
        # get parameters
        self._zcomponent = self._z_component_ui.get_selection()
        self._period = self._frequency_period_ui.get_frequencies()
        self._tolerance = self._tolerance_ui.get_tolerance_in_float()
        self._plotting_object = mtpy.imaging.penetration.Depth3D(
            self._mt_objs, self._period, self._zcomponent, self._tolerance
        )
        self._z_unit = self._z_unit_ui.get_unit()
        self._plotting_object.plot(z_unit=self._z_unit)
        self._fig = self._plotting_object.get_figure()

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
        self._parameter_ui.add_parameter_groupbox(self._z_component_ui)

        # self._frequency_period_ui = FrequencySingle(self._parameter_ui, use_period=True)
        self._frequency_period_ui = FrequencySelection(
            self._parameter_ui,
            show_frequency=False,
            allow_range_select=False,
            select_multiple=False,
        )
        self._parameter_ui.add_parameter_groupbox(self._frequency_period_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._tolerance_ui.setTitle("Period Tolerance")
        self._parameter_ui.add_parameter_groupbox(self._tolerance_ui)

        self._z_unit_ui = ZUnit(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._z_unit_ui)

        self._parameter_ui.end_of_parameter_components()

        # resize
        # self._parameter_ui.resize(self._parameter_ui.width(),
        #                           self._parameter_ui.sizeHint().height())

        self.update_ui()

        # plot parameters
        self._zcomponent = None
        self._period = None
        self._tolerance = None


# register subclass
# VisualizationBase.register(Depth3D)
