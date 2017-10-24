# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
import numpy as np

from mtpy.gui.SmartMT.Components.FigureSetting import ColorBar, Font, AspectRatio, TextBox
from mtpy.gui.SmartMT.gui.plot_control_guis import PlotControlResistivityPhasePseudoSection
from mtpy.gui.SmartMT.gui.plot_parameter_guis import Ellipse, FrequencyTolerance, Arrow, Padding, \
    Scale, Stretch, LineDir, MeshGrid, UniqueFrequencies, FrequencySelect
from mtpy.gui.SmartMT.utils.matplotlib_utils import get_next_fig_num
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection
from mtpy.imaging.plotpseudosection import PlotResPhasePseudoSection


class PhaseTensorMap(VisualizationBase):
    def get_plot_tooltip(self):
        ellipse = self._params['ellipse_dict']
        tipper = self._params['plot_tipper']
        return "freq=%.5f, tolerance=%.2f%%, ellipse_size=%.2f, real_induction=%s, imaginary_induction=%s" % (
            self._params['plot_freq'],
            self._params['ftol'] * 100,
            ellipse['size'],
            'on' if tipper.find('r') >= 0 else 'off',
            'on' if tipper.find('i') >= 0 else 'off'
        )

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

        self._params = {
            'fn_list': [mt_obj.fn for mt_obj in self._mt_objs],
            'plot_freq': self._frequency_ui.get_frequencies(),
            'ftol': self._tolerance_ui.get_tolerance_in_float(),
            'ellipse_dict': self._ellipse_ui.get_ellipse_dict(),
            'mapscale': self._scale_ui.get_mapscale(),
            'tscale': self._scale_ui.get_tscale(),
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            # 'rot_z': self._rotation_ui.get_rotation_in_degree(),  # not implemented in PlotPhaseTensorMaps
            'station_id': (0, 10),
            'fig_num': get_next_fig_num()
        }

        if self._colorbar_ui.isChecked():
            self._params['cb_dict'] = self._colorbar_ui.get_colorbar_dict()

        # arrow_dict = {
        #         'size': 0.5,
        #         'lw': 0.2,
        #         'head_width': 0.04,
        #         'head_length': 0.04,
        #         'threshold': 0.8,
        #         'direction': 0}
        if self._arrow_ui.isChecked():
            self._params['arrow_dict'] = self._arrow_ui.get_arrow_dict()

        if self._padding_ui.isChecked():
            self._params['xpad'] = self._padding_ui.get_x_pad()
            self._params['ypad'] = self._padding_ui.get_y_pad()

        if self._label_font_ui.ui.checkBox_size.isChecked():
            self._params['font_size'] = self._label_font_ui.get_size()

        station_dict = {}
        if self._station_font_ui.ui.checkBox_size.isChecked():
            station_dict['size'] = self._station_font_ui.get_size()
        if self._station_font_ui.ui.checkBox_weight.isChecked():
            station_dict['weight'] = self._station_font_ui.get_weight()
        if self._station_font_ui.ui.checkBox_color.isChecked():
            station_dict['color'] = self._station_font_ui.get_color()
        if station_dict:
            self._params['station_dict'] = station_dict

        self._plotting_object = PlotPhaseTensorMaps(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)
        # set up ui
        self._scale_ui = Scale(self._parameter_ui)
        self._frequency_ui = FrequencySelect(self._parameter_ui,
                                             show_frequency=True,
                                             show_period=False,
                                             allow_range_select=False,
                                             select_multiple=False)
        # self._scale_ui.ui.comboBox_time.currentIndexChanged.connect(
        #     lambda index: self._frequency_ui.show_period() if index == 0
        #     else self._frequency_ui.show_frequency()
        # )
        self._parameter_ui.add_parameter_groupbox(self._scale_ui)
        self._parameter_ui.add_parameter_groupbox(self._frequency_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._tolerance_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._ellipse_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._arrow_ui)

        # not implemented in PlotPhaseTensorMaps
        # self._rotation_ui = Rotation(self._parameter_ui)
        # self._parameter_ui.add_parameter_groubox(self._rotation_ui)

        self._padding_ui = Padding(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._padding_ui)

        self._colorbar_ui = ColorBar(self._parameter_ui)
        self._parameter_ui.add_figure_groupbox(self._colorbar_ui)

        self._label_font_ui = Font(self._parameter_ui)
        self._label_font_ui.setToolTip("Font of the plot labels")
        self._label_font_ui.hide_color()
        self._label_font_ui.hide_weight()
        self._parameter_ui.add_figure_groupbox(self._label_font_ui)

        self._station_font_ui = Font(self._parameter_ui, simple_color=False)
        self._station_font_ui.setTitle('Station Label Font')
        self._parameter_ui.add_figure_groupbox(self._station_font_ui)

        self._parameter_ui.end_of_parameter_components()

        # resize
        # self._parameter_ui.resize(self._parameter_ui.width(),
        #                           self._parameter_ui.sizeHint().height())

        self.update_ui()
        self._params = None


class PhaseTensorPseudoSection(VisualizationBase):
    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)

        # setup gui
        self._plot_control = PlotControlResistivityPhasePseudoSection(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._plot_control)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._ellipse_ui)

        self._linedir_ui = LineDir(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._linedir_ui)

        self._stretch_ui = Stretch(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._stretch_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._arrow_ui)

        self._scale_ui = Scale(self._parameter_ui)
        self._scale_ui.hide_mapscale()
        self._parameter_ui.add_parameter_groupbox(self._scale_ui)

        # this is not implemented in PlotPhaseTensorPseudoSection
        # self._rotation_ui = Rotation(self._parameter_ui)
        # self._parameter_ui.add_parameter_groubox(self._rotation_ui)

        self._colorbar_ui = ColorBar(self._parameter_ui)
        self._parameter_ui.add_figure_groupbox(self._colorbar_ui)

        self._label_font_ui = Font(self._parameter_ui)
        self._label_font_ui.setToolTip("Font of the plot labels")
        self._label_font_ui.hide_color()
        self._label_font_ui.hide_weight()
        self._parameter_ui.add_figure_groupbox(self._label_font_ui)

        self._parameter_ui.end_of_parameter_components()

        self.update_ui()
        self._params = None

    def get_plot_tooltip(self):
        ellipse = self._params['ellipse_dict']
        tipper = self._params['plot_tipper']
        stretch = self._params['stretch']
        return "ellipse_size=%.2f, stretch=(%.2f,%.2f), linedir=%s, real_induction=%s, imaginary_induction=%s" % (
            ellipse['size'],
            stretch[0],
            stretch[1],
            self._params['linedir'],
            'on' if tipper.find('r') >= 0 else 'off',
            'on' if tipper.find('i') >= 0 else 'off'
        )

    def update_ui(self):
        # nothing to be done here
        pass

    def plot(self):
        # get parameters
        self._params = {
            'fn_list': [mt_obj.fn for mt_obj in self._mt_objs],
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            'tscale': self._scale_ui.get_tscale(),
            'ellipse_dict': self._ellipse_ui.get_ellipse_dict(),
            'stretch': self._stretch_ui.get_stretch(),
            'linedir': self._linedir_ui.get_linedir(),
            # 'rotz': self._rotation_ui.get_rotation_in_degree(), # this is not implemented in PlotPhaseTensorPseudoSection
            # default for testing
            'station_id': (0, 10),  # indices for showing station names,
            'fig_num': get_next_fig_num()
        }
        if self._arrow_ui.isChecked():
            self._params['arrow_dict'] = self._arrow_ui.get_arrow_dict()

        if self._colorbar_ui.isChecked():
            self._params['cb_dict'] = self._colorbar_ui.get_colorbar_dict()

        if self._label_font_ui.ui.checkBox_size.isChecked():
            self._params['font_size'] = self._label_font_ui.get_size()

        if self._stretch_ui.ui.checkBox_x_range.isChecked():
            self._params['xlimits'] = self._stretch_ui.get_x_limits()

        if self._stretch_ui.ui.checkBox_y_range.isChecked():
            self._params['ylimits'] = self._stretch_ui.get_y_limits()

        self._plotting_object = PlotPhaseTensorPseudoSection(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    @staticmethod
    def plot_description():
        return """
<p>plot the phase tensor ellipses in a pseudo section format</p>
        """

    @staticmethod
    def plot_name():
        return "Phase Tensor Pseudo Section"


class ResistivityPhasePseudoSection(VisualizationBase):
    @staticmethod
    def plot_description():
        return """
        <p> Plot resistivity and phase pseudo section for different components</p>
        """

    def plot(self):
        self._params = {
            'fn_list': [mt_obj.fn for mt_obj in self._mt_objs],
            'plot_style': self._mesh_grid_ui.get_grid_type(),
            'imshow_interp': self._mesh_grid_ui.get_interpolation_method(),
            'ftol': self._tolerance_ui.get_tolerance_in_float(),
            'linedir': self._linedir_ui.get_linedir(),
            'aspect': self._aspect_ui.get_aspect(),
            'plot_xx': self._plot_control.get_plot_xx(),
            'plot_xy': self._plot_control.get_plot_xy(),
            'plot_yx': self._plot_control.get_plot_yx(),
            'plot_yy': self._plot_control.get_plot_yy(),
            'res_cmap': self._plot_control.get_res_cmap(),
            'phase_cmap': self._plot_control.get_phase_cmap(),
            'xtickspace': self._plot_control.get_tickspace(),
            'stationid': (0, 20),
            'plot_yn': 'n',  # do not plot on class creation
            'fig_num': get_next_fig_num()
        }

        param = self._font_ui.get_size()
        if param is not None:
            self._params['font_size'] = param
        param = self._text_box.get_size()
        if param is not None:
            self._params['text_size'] = param
        param = self._text_box.get_weight()
        if param is not None:
            self._params['text_weight'] = param
        param = self._text_box.get_location()
        if param is not None:
            self._params['text_location'] = param
        if self._text_box.ui.groupBox_padding.isChecked():
            self._params['text_xpad'] = self._text_box.get_xpad()
            self._params['text_ypad'] = self._text_box.get_ypad()
        param = self._plot_control.get_period_limit()
        if param is not None:
            self._params['period_limits'] = param
        param = self._plot_control.get_phase_limit()
        if param is not None:
            self._params['phase_limits'] = param
        param = self._plot_control.get_resistivity_limits()
        if param is not None:
            self._params['res_limits'] = param
        if self._period_ui.isChecked():
            self._params['plot_period'] = np.array(self._period_ui.get_frequency_list())

        self._plotting_object = PlotResPhasePseudoSection(**self._params)
        self._plotting_object.plot(show=False)
        self._fig = self._plotting_object.fig

    def update_ui(self):
        self._period_ui.set_data(self._mt_objs)

    @staticmethod
    def plot_name():
        return "Resistivity and Phase Pseudo Section"

    def get_plot_tooltip(self):
        pass

    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)

        # setup gui
        self._plot_control = PlotControlResistivityPhasePseudoSection(self._parameter_ui)
        self.parameter_ui.add_parameter_groupbox(self._plot_control)

        self._mesh_grid_ui = MeshGrid(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._mesh_grid_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._parameter_ui.add_parameter_groupbox(self._tolerance_ui)

        self._linedir_ui = LineDir(self._parameter_ui)
        self._linedir_ui.ui.radioButton_ew.setChecked(True)
        self._parameter_ui.add_parameter_groupbox(self._linedir_ui)

        self._period_ui = UniqueFrequencies(self._parameter_ui, use_period=True)
        self._period_ui.setCheckable(True)
        self._period_ui.setChecked(False)
        self._parameter_ui.add_parameter_groupbox(self._period_ui)

        # figure settings
        self._aspect_ui = AspectRatio(self._parameter_ui)
        self._parameter_ui.add_figure_groupbox(self._aspect_ui)

        self._font_ui = Font(self._parameter_ui, point_size=True)
        self._font_ui.hide_color()
        self._font_ui.hide_weight()
        self._parameter_ui.add_figure_groupbox(self._font_ui)

        self._text_box = TextBox(self._parameter_ui, point_size=True, key_size=True)
        self._parameter_ui.add_figure_groupbox(self._text_box)

        self._parameter_ui.end_of_parameter_components()

        self.update_ui()
        self._params = None
