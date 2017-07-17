# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.gui.plot_parameter import FrequencySingle, Ellipse, FrequencyTolerance, ColorBar, Arrow, Padding, \
    Scale, Font, Stretch, LineDir
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.imaging.phase_tensor_maps import PlotPhaseTensorMaps
from mtpy.imaging.phase_tensor_pseudosection import PlotPhaseTensorPseudoSection


class PhaseTensorMap(VisualizationBase):
    def get_parameter_str(self):
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
        file_list = []
        for mt_obj in self._mt_objs:
            file_list.append(mt_obj.fn)

        self._params = {
            'fn_list': file_list,
            'plot_freq': self._frequency_ui.get_frequency(),
            'ftol': self._tolerance_ui.get_tolerance_in_float(),
            'ellipse_dict': self._ellipse_ui.get_ellipse_dict(),
            'mapscale': self._scale_ui.get_mapscale(),
            'tscale': self._scale_ui.get_tscale(),
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            'station_id': (0, 10)
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
        self._frequency_ui = FrequencySingle(self._parameter_ui)
        self._frequency_ui.setTitle("Frequency (Hz)")
        self._parameter_ui.add_parameter_groubox(self._frequency_ui)

        self._tolerance_ui = FrequencyTolerance(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._tolerance_ui)

        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._arrow_ui)

        self._scale_ui = Scale(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._scale_ui)

        self._padding_ui = Padding(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._padding_ui)

        self._colorbar_ui = ColorBar(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._colorbar_ui)

        self._label_font_ui = Font(self._parameter_ui)
        self._label_font_ui.setToolTip("Font of the plot labels")
        self._label_font_ui.hide_color()
        self._label_font_ui.hide_weight()
        self._parameter_ui.add_parameter_groubox(self._label_font_ui)

        self._station_font_ui = Font(self._parameter_ui, simple_color=False)
        self._station_font_ui.setTitle('Station Label Font')
        self._parameter_ui.add_parameter_groubox(self._station_font_ui)

        # resize
        self._parameter_ui.resize(self._parameter_ui.width(),
                                  self._parameter_ui.sizeHint().height())

        self.update_ui()
        self._params = None


class PhaseTensorPseudoSection(VisualizationBase):
    def __init__(self, parent):
        VisualizationBase.__init__(self, parent)

        # setup gui
        self._ellipse_ui = Ellipse(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._ellipse_ui)

        self._linedir_ui = LineDir(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._linedir_ui)

        self._stretch_ui = Stretch(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._stretch_ui)

        self._arrow_ui = Arrow(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._arrow_ui)

        self._scale_ui = Scale(self._parameter_ui)
        self._scale_ui.hide_mapscale()
        self._parameter_ui.add_parameter_groubox(self._scale_ui)

        self._colorbar_ui = ColorBar(self._parameter_ui)
        self._parameter_ui.add_parameter_groubox(self._colorbar_ui)

        self._label_font_ui = Font(self._parameter_ui)
        self._label_font_ui.setToolTip("Font of the plot labels")
        self._label_font_ui.hide_color()
        self._label_font_ui.hide_weight()
        self._parameter_ui.add_parameter_groubox(self._label_font_ui)

        self.update_ui()
        self._params = None

    def get_parameter_str(self):
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
        file_list = []
        for mt_obj in self._mt_objs:
            file_list.append(mt_obj.fn)
        self._params = {
            'fn_list': file_list,
            'plot_tipper': self._arrow_ui.get_plot_tipper(),
            'tscale': self._scale_ui.get_tscale(),
            'ellipse_dict': self._ellipse_ui.get_ellipse_dict(),
            'stretch': self._stretch_ui.get_stretch(),
            'linedir': self._linedir_ui.get_linedir(),
            # default for testing
            'station_id': (0, 10),  # indices for showing station names,
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
