# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 21/08/2017
"""
from mtpy.gui.SmartMT.Components.FigureSetting import Font
from mtpy.gui.SmartMT.Components.PlotParameter import FrequencyTolerance, Rotation
from mtpy.gui.SmartMT.gui.plot_control_guis import PlotControlStrike
from mtpy.gui.SmartMT.utils.matplotlib_utils import get_next_fig_num
from mtpy.gui.SmartMT.visualization import VisualizationBase
from mtpy.imaging.plotstrike import PlotStrike


class Strike(VisualizationBase):
    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # setup ui
        self._plot_control_ui = PlotControlStrike(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._plot_control_ui)

        self._rotation_ui = Rotation(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._rotation_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._tolerance_ui.ui.doubleSpinBox.setValue(0.05)  # set default value
        self._parameter_ui.add_parameter_groupbox(self._tolerance_ui)

        self._font_ui = Font(self._parameter_ui)
        self._font_ui.hide_weight()
        self._font_ui.hide_color()
        self._parameter_ui.add_figure_groupbox(self._font_ui)

        self._parameter_ui.end_of_parameter_components()

        self.update_ui()
        self._params = None

    def plot(self):
        # set up params
        self._params = {
            'fn_list': [mt_obj.fn for mt_obj in self._mt_objs],
            'rot_z': self._rotation_ui.get_rotation_in_degree(),
            'period_tolerance': self._tolerance_ui.get_tolerance_in_float(),
            'plot_range': self._plot_control_ui.get_plot_range(),
            'plot_type': self._plot_control_ui.get_plot_type(),
            'plot_tipper': self._plot_control_ui.get_plot_tipper(),
            'pt_error_floor': self._plot_control_ui.get_error_floor(),
            'fold': self._plot_control_ui.get_fold(),
            'fig_dpi': 100,
            "plot_yn": 'n',
            "fig_num": get_next_fig_num()
        }
        param = self._font_ui.get_size()
        if param is not None:
            self._params['font_size'] = param
        self._plotting_object = PlotStrike(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    def update_ui(self):
        pass

    @staticmethod
    def plot_name():
        return "Strike"

    @staticmethod
    def plot_description():
        return """
        <p>This plots the estrike estimated from the invariants, 
        phase tensor and the tipper in either a rose diagram of 
        xy plot</p>
        <p>plots the strike angle as determined by phase tensor 
        azimuth (Caldwell et al. [2004]) and invariants of the 
        impedance tensor (Weaver et al. [2003]).</p>
        <p>The data is split into decades where the histogram 
        for each is plotted in the form of a rose diagram with a
        range of 0 to 180 degrees. Where 0 is North and 90 is 
        East.   The median angle of the period band is set in 
        polar diagram.  The top row is the strike estimated from
        the invariants of the impedance tensor.  The bottom row 
        is the azimuth estimated from the phase tensor.  If 
        tipper is plotted then the 3rd row is the strike determined
        from the tipper, which is orthogonal to the induction 
        arrow direction.</p>
        """

    def get_plot_tooltip(self):
        pass
