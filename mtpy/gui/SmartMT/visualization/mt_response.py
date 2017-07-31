# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtCore

from mtpy.gui.SmartMT.gui.plot_parameter import StationSelection, Rotation, PlotControlMTResponse, Arrow, Ellipse, \
    CommonSettings
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.plot_mt_response import PlotMTResponse
from mtpy.imaging.plotnresponses import PlotMultipleResponses


class MTResponse(VisualizationBase):
    def plot(self):
        # get parameters
        self._station = self._station_ui.get_station()

        self._params = {
            'fn': self._station.fn,
            'rot_z': self._rotation_ui.get_rotation_in_degree(),
            'plot_num': self._plot_control_ui.get_plot_num(),
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            'plot_strike': self._plot_control_ui.get_strike(),
            'plot_skew': self._plot_control_ui.get_skew(),
            'plot_pt': self._plot_control_ui.get_ellipses(),
            # 'plot_title': self._common_ui.get_title()
        }

        if self._arrow_ui.ui.groupBox_advanced_options.isChecked():
            self._params['arrow_dict'] = self._arrow_ui.get_arrow_dict()

        if not self._ellipse_ui.isHidden():
            self._params['ellipse_dict'] = self._ellipse_ui.get_ellipse_dict()

        # plot
        self._plotting_object = PlotMTResponse(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # setup attributes here
        self._station = None
        self._params = None

        # set up parameter GUIs here

        self._station_ui = StationSelection(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._station_ui)
        self._station_ui.station_changed.connect(self._station_changed)

        # self._title_ui = CommonSettings(self._parameter_ui)
        # self._parameter_ui.add_parameter_groubox(self._title_ui)

        self._plot_control_ui = PlotControlMTResponse(self._parameter_ui)
        self._plot_control_ui.hide_plot_style()
        self._parameter_ui.add_parameter_groubox(self._plot_control_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._ellipse_ui.setHidden(True)
        # make the radio button toggle hidden of the ellipses groupbox
        self._plot_control_ui.ui.radioButton_ellipses_y.toggled.connect(self._ellipse_radio_button_toggled)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._arrow_ui)

        self._rotation_ui = Rotation(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._rotation_ui)

        self._parameter_ui.end_of_parameter_components()

    def _station_changed(self):
        # change title
        self._common_ui.set_title(self._station_ui.get_station().station)

    def _ellipse_radio_button_toggled(self, b):
        self._ellipse_ui.setHidden(not self._ellipse_ui.isHidden())

    def update_ui(self):
        self._station_ui.set_data(self._mt_objs)

    @staticmethod
    def plot_name():
        return "MT Response"

    def get_parameter_str(self):
        return "station=%s, rotation=%.2f°, plot_num=%d" % (
            self._station.station,
            self._params['rot_z'],
            self._params['plot_num']
        )

    @staticmethod
    def plot_description():
        return """
<p>Plots Resistivity and phase for the different modes of the MT response.</p>
<p>The plot places the apparent resistivity in log scale in the top panel(s),
depending on the plot_num. The phase is below this, note that 180 degrees 
has been added to the yx phase so the xy and yx phases plot in the same quadrant.
Both the resistivity and phase share the same x-axis which is in log period, 
short periods on the left to long periods on the right. So if you zoom in on the 
plot both plots will zoom in to the same x-coordinates.  If there is tipper 
information, you can plot the tipper as a third panel at the bottom, and also 
shares the x-axis.  The arrows are in the convention of pointing towards a 
conductor. The xx and yy components can be plotted as well, this adds two panels 
on the right. Here the phase is left unwrapped. Other parameters can be added as 
subplots such as strike, skew and phase tensor ellipses.</p>
<p><strong>Note:</strong> This only plot one station per plot.</p>
        """


class MultipleMTResponses(VisualizationBase):
    def get_parameter_str(self):
        return "rotation=%.2f°, plot_num=%d" % (
            self._params['rot_z'],
            self._params['plot_num'],
        )

    @staticmethod
    def plot_description():
        return """
<p>Plots multiple MT responses simultaneously in one plot of sub-figures 
or in a single plot with subfigures for comparing each component</p>
<p>The plot places the apparent resistivity in log scale in the top panel(s),
depending on the plot_num. The phase is below this, note that 180 degrees 
has been added to the yx phase so the xy and yx phases plot in the same quadrant.
Both the resistivity and phase share the same x-axis which is in log period, 
short periods on the left to long periods on the right. So if you zoom in on the 
plot both plots will zoom in to the same x-coordinates.  If there is tipper 
information, you can plot the tipper as a third panel at the bottom, and also 
shares the x-axis.  The arrows are in the convention of pointing towards a 
conductor. The xx and yy components can be plotted as well, this adds two panels 
on the right. Here the phase is left unwrapped. Other parameters can be added as 
subplots such as strike, skew and phase tensor ellipses.</p>
<p><strong>Warning</strong> This plot could take very long time.</p>
        """

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # setup attributes here
        self._params = None

        # set up parameter GUIs here
        # self._title_ui = CommonSettings(self._parameter_ui)
        # self._parameter_ui.add_parameter_groubox(self._title_ui)

        self._plot_control_ui = PlotControlMTResponse(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._plot_control_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._ellipse_ui.setHidden(True)
        # make the radio button toggle hidden of the ellipses groupbox
        self._plot_control_ui.ui.radioButton_ellipses_y.toggled.connect(self._ellipse_radio_button_toggled)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._arrow_ui)

        self._rotation_ui = Rotation(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._rotation_ui)

        self._parameter_ui.end_of_parameter_components()

    def _ellipse_radio_button_toggled(self, b):
        self._ellipse_ui.setHidden(not self._ellipse_ui.isHidden())

    def update_ui(self):
        pass

    @staticmethod
    def plot_name():
        return "Multiple MT responses"

    def plot(self):
        # get parameters
        self._params = {
            'fn_list': [mt_obj.fn for mt_obj in self._mt_objs],
            'rot_z': self._rotation_ui.get_rotation_in_degree(),
            'plot_num': self._plot_control_ui.get_plot_num(),
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            'plot_strike': self._plot_control_ui.get_strike(),
            'plot_skew': self._plot_control_ui.get_skew(),
            'plot_pt': self._plot_control_ui.get_ellipses(),
            # 'plot_title': self._common_ui.get_title(),
            'plot_style': self._plot_control_ui.get_style(),
            'plot_yn': 'n'
        }

        if self._arrow_ui.ui.groupBox_advanced_options.isChecked():
            self._params['arrow_dict'] = self._arrow_ui.get_arrow_dict()

        if not self._ellipse_ui.isHidden():
            self._params['ellipse_dict'] = self._ellipse_ui.get_ellipse_dict()

        # plot
        self._plotting_object = PlotMultipleResponses(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig
