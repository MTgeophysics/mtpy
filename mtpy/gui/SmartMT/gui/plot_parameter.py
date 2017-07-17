# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
import six
from matplotlib import colors as mcolors
import numpy as np
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.ui_asset.groupbox_arrow import Ui_GroupBox_Arrow
from mtpy.gui.SmartMT.ui_asset.groupbox_color_bar import Ui_GroupBox_ColorBar
from mtpy.gui.SmartMT.ui_asset.groupbox_ellipse import Ui_GroupBoxEllipse
from mtpy.gui.SmartMT.ui_asset.groupbox_font import Ui_GroupBox_Font
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_single import Ui_groupBoxFrequency_pereiod_single
from mtpy.gui.SmartMT.ui_asset.groupbox_linedir import Ui_GroupBox_Linedir
from mtpy.gui.SmartMT.ui_asset.groupbox_padding import Ui_GroupBox_Padding
from mtpy.gui.SmartMT.ui_asset.groupbox_scale import Ui_GroupBox_Scale
from mtpy.gui.SmartMT.ui_asset.groupbox_stretch import Ui_GroupBox_Stretch
from mtpy.gui.SmartMT.ui_asset.groupbox_tolerance import Ui_GroupBoxTolerance
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_multiple import Ui_groupBoxZ_Component_Multiple
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_single import Ui_groupBoxZ_Component_Single
from mtpy.gui.SmartMT.ui_asset.plot_parameters import Ui_GroupBoxParameters


COLORS = list(six.iteritems(mcolors.cnames))
# # add the single letter colors
# for name, rgb in six.iteritems(mcolors.ColorConverter.colors):
#     hex_ = mcolors.rgb2hex(rgb)
#     COLORS.append((name, hex_))
# sort by name
COLORS.sort(key=lambda c: c[0])

SIMPLE_COLORS = ['b',  # blue
                 'g',  # green
                 'r',  # red
                 'c',  # cyan
                 'm',  # magenta
                 'y',  # yellow
                 'k',  # black
                 'w'  # white
              ]


class PlotParameter(QtGui.QGroupBox):
    _slider_tick_size = 10.0
    _slider_min = 0.0
    _slider_max = 100.0

    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxParameters()
        self.ui.setupUi(self)

    def add_parameter_groubox(self, groupbox):
        self.ui.verticalLayout_2.addWidget(groupbox, QtCore.Qt.AlignLeft)
        self.resize(self.sizeHint())


class ZComponentMultiple(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Multiple()
        self.ui.setupUi(self)
        # z-component checkbox logic
        self.ui.checkBox_zyx.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_zxy.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_det.stateChanged.connect(self._multiple_zcomponent_logic)

    def _multiple_zcomponent_logic(self, int):
        """
        set up at least one component selected
        :return:
        """
        counter = 0
        if self.ui.checkBox_det.isChecked() + self.ui.checkBox_zxy.isChecked() + self.ui.checkBox_zyx.isChecked() == 1:
            # only one checkbox is checked, lock the checked box
            if self.ui.checkBox_det.isChecked():
                self.ui.checkBox_det.setEnabled(False)
            elif self.ui.checkBox_zxy.isChecked():
                self.ui.checkBox_zxy.setEnabled(False)
            else:
                self.ui.checkBox_zyx.setEnabled(False)
        else:
            self.ui.checkBox_det.setEnabled(True)
            self.ui.checkBox_zxy.setEnabled(True)
            self.ui.checkBox_zyx.setEnabled(True)

    def get_selection(self):
        zcomponent = []
        if self.parameter_ui.ui.radioButton_det.isChecked():
            zcomponent.append('det')
        if self.parameter_ui.ui.radioButton_zxy.isChecked():
            zcomponent.append('zxy')
        if self.parameter_ui.ui.radioButton_zyx.isChecked():
            zcomponent.append('zyx')
        return zcomponent


class ZComponentSingle(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Single()
        self.ui.setupUi(self)

    def get_selection(self):
        if self.ui.radioButton_det.isChecked():
            return 'det'
        elif self.ui.radioButton_zxy.isChecked():
            return 'zxy'
        elif self.ui.radioButton_zyx.isChecked():
            return 'zyx'


class FrequencySingle(QtGui.QGroupBox):
    """
    Frequency selection (single frequency)
    """

    _unit_period = 'second'
    _unit_frequency = 'Hz'
    _title_period = 'Period'
    _title_frequency = 'Frequency'

    def __init__(self, parent, use_period=False):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.use_period = use_period
        self.ui = Ui_groupBoxFrequency_pereiod_single()
        self.ui.setupUi(self)
        self._histogram = FrequencySingle.FrequencyHistogram()
        self.set_use_period(self.use_period)
        # add matplotlib canvas
        self.ui.verticalLayoutFrequencyPeriod.addWidget(self._histogram)
        # connect components
        # self.ui.horizontalSliderPeriod.valueChanged.connect(lambda value: self.update_period_text(value))
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_histogram)
        self.ui.comboBoxPeriod.editTextChanged.connect(self.update_histogram)
        self._histogram.mpl_connect('button_release_event', self._mouse_pick)

    def toggle_time_scale(self, *args):
        self.use_period = not self.use_period
        self.set_use_period(self.use_period)

    def set_use_period(self, use_period=False):
        if use_period:
            self._histogram.set_unit(self._unit_period)
            self._histogram.set_title(self._title_period)
            title = '%s (%s)' % (self._title_period, self._unit_period)
        else:
            self._histogram.set_unit(self._unit_frequency)
            self._histogram.set_title(self._unit_frequency)
            title = '%s (%s)' % (self._title_frequency, self._unit_frequency)
        self.setTitle(title)
        self._update_frequency()

    def _mouse_pick(self, event):
        if not event.inaxes:
            return
        x = event.xdata
        self.ui.comboBoxPeriod.setEditText("%.5f" % x)

    def get_frequency(self):
        return float(self.ui.comboBoxPeriod.currentText())

    def update_histogram(self):
        value = float(self.ui.comboBoxPeriod.currentText())
        self._histogram.set_current_frequency(value)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._update_frequency()

    def _update_frequency(self, ):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))

            if self.use_period:
                all_periods = 1.0 / np.array(all_freqs)
                self._histogram.set_data(all_periods)
                all_unique = sorted(list(set(all_periods)))
            else:
                self._histogram.set_data(all_freqs)
                all_unique = sorted(list(set(all_freqs)))
            self._histogram.update_figure()
            # sort all frequencies in ascending order
            for period in all_unique:
                self.ui.comboBoxPeriod.addItem("%.5f" % period)
            self.ui.comboBoxPeriod.setCurrentIndex(0)
            self.update_histogram()

    class FrequencyHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, hight=3, dpi=100):
            self.artists = dict()
            self._frequency = None
            self._current_period = None
            self._title = None
            self._unit = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self._lx = None
            self.cursor = None

            # self.mpl_connect('motion_notify_event', self.cursor)
            self.mpl_connect('button_release_event', self.mouse_pick)
            self.setMinimumSize(200, 150)

        # def mouse_move(self, event):
        #     if not event.inaxes:
        #         return
        #     x = event.xdata
        #     y = event.ydata
        #     if self._cursor_x is None:
        #         self._cursor_x = self._axes.axvline(linewidth=1, color="green")
        #     if self._cursor_text is None:
        #         self._cursor_text = self._axes.text(0.0, 0.0, '', fontsize=8)
        #     self._cursor_x.set_xdata(x)
        #     self._cursor_text.set_text('period=%.2f' % x)
        #     self._cursor_text.set_position((x, y))
        #     self.draw()

        def set_title(self, title):
            self._title = title

        def set_unit(self, unit):
            if unit != self._unit:
                self._unit = unit
                self.cursor = Cursor(self._axes, track_y=False, text_format="%f " + self._unit, useblit=True)

        def mouse_pick(self, event):
            if not event.inaxes:
                return
            x = event.xdata
            self.set_current_frequency(x)

        def compute_initial_figure(self):
            if self._frequency is not None:
                self._axes.tick_params(axis='both', which='major', labelsize=6)
                self._axes.tick_params(axis='both', which='minor', labelsize=4)
                self._axes.hist(self._frequency, 50, normed=1)
                if self._title and self._unit:
                    self._axes.set_xlabel("%s (%s)" % (self._title, self._unit), fontsize=8)
                    self.figure.suptitle('%s Distribution in Selected Stations' % self._title, fontsize=8)

        def set_data(self, frequency):
            self._frequency = frequency

        def set_current_frequency(self, freq):
            self._current_period = freq
            if self._lx is None:
                self._lx = self._axes.axvline(linewidth=2, color="red")
            self._lx.set_xdata(self._current_period)
            # if self._fig.canvas.supports_blit:
            #     self._axes.draw_artist(self._lx)
            #     self._fig.canvas.blit(self._axes.bbox)
            # else:
            #     self._fig.canvas.draw_idle()
            self._fig.canvas.draw_idle()

        def update_figure(self):
            # clear figure
            self._axes.cla()
            self.compute_initial_figure()
            self.draw()


class Ellipse(QtGui.QGroupBox):
    """
    ellipse_dict defined for mtpy.imagining.phase_tensor_maps.PlogPhaseTensorMaps
    """
    _colorby = ['phimin', 'phimax', 'skew', 'skew_seg', 'normalized_skew', 'normalized_skew_seg', 'phidet', 'ellipticity']
    _cmap = ['mt_yl2rd', 'mt_bl2yl2rd', 'mt_wh2bl', 'mt_rd2bl', 'mt_bl2wh2rd', 'mt_seg_bl2wh2rd', 'mt_rd2gr2bl']

    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxEllipse()
        self.ui.setupUi(self)
        # set tooltips for color by and cmap
        for i in range(self.ui.comboBoxColor_by.count()):
            self.ui.comboBoxColor_by.setItemData(i, self.ui.comboBoxColor_by.itemText(i), QtCore.Qt.ToolTipRole)
        for i in range(self.ui.comboBox_cmap.count()):
            self.ui.comboBox_cmap.setItemData(i, self.ui.comboBox_cmap.itemText(i), QtCore.Qt.ToolTipRole)

        self.ui.doubleSpinBox_min.editingFinished.connect(self._increase_max)
        self.ui.doubleSpinBox_max.editingFinished.connect(self._decrease_min)

    def _increase_max(self):
        value = self.ui.doubleSpinBox_min.value()
        if value > self.ui.doubleSpinBox_max.value():
            self.ui.doubleSpinBox_max.setValue(value)

    def _decrease_min(self):
        value = self.ui.doubleSpinBox_max.value()
        if value < self.ui.doubleSpinBox_min.value():
            self.ui.doubleSpinBox_min.setValue(value)

    def get_ellipse_dict(self):
        ellipse_dict = {
            'size': self.ui.doubleSpinBox_size.value(),
            'colorby': self._colorby[self.ui.comboBoxColor_by.currentIndex()],
            'range': (
                self.ui.doubleSpinBox_min.value(), self.ui.doubleSpinBox_max.value(),
                self.ui.doubleSpinBox_step.value()),
            'cmap': self._cmap[self.ui.comboBox_cmap.currentIndex()]
        }
        return ellipse_dict


class FrequencyTolerance(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxTolerance()
        self.ui.setupUi(self)

    def get_tolerance_in_float(self):
        return self.ui.doubleSpinBox.value() / 100.0


class ColorBar(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_ColorBar()
        self.ui.setupUi(self)

        # connect event
        self.ui.horizontalSlider_x.valueChanged.connect(lambda value: self.ui.doubleSpinBox_x.setValue(value / 100.0))
        self.ui.horizontalSlider_y.valueChanged.connect(lambda value: self.ui.doubleSpinBox_y.setValue(value / 100.0))
        self.ui.horizontalSlider_width.valueChanged.connect(
            lambda value: self.ui.doubleSpinBox_width.setValue(value / 100.0))
        self.ui.horizontalSlider_height.valueChanged.connect(
            lambda value: self.ui.doubleSpinBox_height.setValue(value / 100.0))
        self.ui.comboBox_orientation.currentIndexChanged.connect(self._orientation_changed)
        self.ui.doubleSpinBox_x.editingFinished.connect(self._update_slider_x)
        self.ui.doubleSpinBox_y.editingFinished.connect(self._update_slider_y)
        self.ui.doubleSpinBox_width.editingFinished.connect(self._update_slider_width)
        self.ui.doubleSpinBox_height.editingFinished.connect(self._update_slider_height)

    def _orientation_changed(self, *args):
        x = self.ui.doubleSpinBox_x.value()
        y = self.ui.doubleSpinBox_y.value()
        width = self.ui.doubleSpinBox_width.value()
        height = self.ui.doubleSpinBox_height.value()
        self.ui.horizontalSlider_x.setValue(100 - x * 100)
        self.ui.horizontalSlider_y.setValue(100 - y * 100)
        self.ui.horizontalSlider_width.setValue(height * 100)
        self.ui.horizontalSlider_height.setValue(width * 100)

    def _update_slider_x(self):
        value = int(self.ui.doubleSpinBox_x.value() * 100)
        self.ui.horizontalSlider_x.setValue(value)

    def _update_slider_y(self):
        value = int(self.ui.doubleSpinBox_y.value() * 100)
        self.ui.horizontalSlider_y.setValue(value)

    def _update_slider_width(self):
        value = int(self.ui.doubleSpinBox_width.value() * 100)
        self.ui.horizontalSlider_width.setValue(value)

    def _update_slider_height(self):
        value = int(self.ui.doubleSpinBox_height.value() * 100)
        self.ui.horizontalSlider_height.setValue(value)

    _cb_orientation = ['vertical', 'horizontal']

    def get_colorbar_dict(self):
        if self.isChecked():
            cb_dict = {
                'orientataion': self._cb_orientation[self.ui.comboBox_orientation.currentIndex()],
                'position': (
                    self.ui.doubleSpinBox_x.value(),
                    self.ui.doubleSpinBox_y.value(),
                    self.ui.doubleSpinBox_width.value(),
                    self.ui.doubleSpinBox_height.value()
                )
            }
            return cb_dict
        else:
            return None


class Arrow(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Arrow()
        self.ui.setupUi(self)
        self._simple_color = simple_color
        if not self._simple_color:
            # use all colors available to matplot
            self.ui.comboBox_color_imaginary.clear()
            self.ui.comboBox_color_real.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color_imaginary.addItems(cnames)
            self.ui.comboBox_color_real.addItems(cnames)

    _direction = [0, 1]

    def get_arrow_dict(self):
        if self.ui.groupBox_advanced_options.isChecked():
            arrow_dict = {
                'size': self.ui.doubleSpinBox_size.value(),
                'head_length': self.ui.doubleSpinBox_head_length.value(),
                'head_width': self.ui.doubleSpinBox_head_width.value(),
                'lw': self.ui.doubleSpinBox_line_width.value(),
                'threshold': self.ui.doubleSpinBox_threshold.value(),
                'direction': self._direction[self.ui.comboBox_direction.currentIndex()]
            }
            if self._simple_color:
                arrow_dict['color'] = (SIMPLE_COLORS[self.ui.comboBox_color_real.currentIndex()],
                                       SIMPLE_COLORS[self.ui.comboBox_color_imaginary.currentIndex()])
            else:
                arrow_dict['color'] = (COLORS[self.ui.comboBox_color_real.currentIndex()][1],
                                       COLORS[self.ui.comboBox_color_imaginary.currentIndex()][1])
            return arrow_dict
        else:
            return None

    def get_plot_tipper(self):
        if self.ui.checkBox_real.isChecked() and self.ui.checkBox_imaginary.isChecked():
            return 'yri'
        elif self.ui.checkBox_real.isChecked():
            return 'yr'
        elif self.ui.checkBox_imaginary.isChecked():
            return 'ri'
        else:
            return 'n'


class Padding(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Padding()
        self.ui.setupUi(self)

    def get_x_pad(self):
        return self.ui.doubleSpinBox_x.value()

    def get_y_pad(self):
        return self.ui.doubleSpinBox_y.value()


class Scale(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Scale()
        self.ui.setupUi(self)

    _tscale = ['period', 'freq']
    _mapscale = ['deg', 'm', 'km']

    def get_tscale(self):
        return self._tscale[self.ui.comboBox_time.currentIndex()]

    def hide_mapscale(self):
        self.ui.label_map.hide()
        self.ui.comboBox_map.hide()

    def get_mapscale(self):
        return self._mapscale[self.ui.comboBox_map.currentIndex()]


class Font(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Font()
        self.ui.setupUi(self)
        self.ui.checkBox_size.stateChanged.connect(self.size_state_changed)
        self.ui.checkBox_weight.stateChanged.connect(self.weight_state_changed)
        self.ui.checkBox_color.stateChanged.connect(self.color_state_changed)
        self._simple_color = simple_color
        if not self._simple_color:
            self.ui.comboBox_color.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color.addItems(cnames)

    def size_state_changed(self, p_int):
        if p_int == 0:
            self.ui.spinBox_size.setEnabled(False)
        else:
            self.ui.spinBox_size.setEnabled(True)

    def weight_state_changed(self, p_int):
        if p_int == 0:
            self.ui.comboBox_weight.setEnabled(False)
        else:
            self.ui.comboBox_weight.setEnabled(True)

    def color_state_changed(self, p_int):
        if p_int == 0:
            self.ui.comboBox_color.setEnabled(False)
        else:
            self.ui.comboBox_color.setEnabled(True)

    def hide_size(self):
        self.ui.spinBox_size.hide()
        self.ui.checkBox_size.hide()

    def hide_weight(self):
        self.ui.comboBox_weight.hide()
        self.ui.checkBox_weight.hide()

    def hide_color(self):
        self.ui.comboBox_color.hide()
        self.ui.checkBox_color.hide()

    def get_size(self):
        return self.ui.spinBox_size.value()

    def get_weight(self):
        return str(self.ui.comboBox_weight.currentText())

    def get_color(self):
        if self._simple_color:
            return SIMPLE_COLORS[self.ui.comboBox_color.currentIndex()]
        else:
            return COLORS[self.ui.comboBox_color.currentIndex()][1]


class Stretch(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Stretch()
        self.ui.setupUi(self)
        self.ui.checkBox_x_range.stateChanged.connect(self._x_range_state_change)
        self.ui.checkBox_y_range.stateChanged.connect(self._y_range_state_change)

    def _x_range_state_change(self, p_int):
        if p_int == 0:
            self.ui.doubleSpinBox_x_min.setEnabled(False)
            self.ui.doubleSpinBox_x_max.setEnabled(False)
        else:
            self.ui.doubleSpinBox_x_min.setEnabled(True)
            self.ui.doubleSpinBox_x_max.setEnabled(True)

    def _y_range_state_change(self, p_int):
        if p_int == 0:
            self.ui.doubleSpinBox_y_min.setEnabled(False)
            self.ui.doubleSpinBox_y_max.setEnabled(False)
        else:
            self.ui.doubleSpinBox_y_min.setEnabled(True)
            self.ui.doubleSpinBox_y_max.setEnabled(True)

    def get_stretch(self):
        return self.ui.doubleSpinBox_x.value(), self.ui.doubleSpinBox_y.value()

    def get_x_limits(self):
        return self.ui.doubleSpinBox_x_min.value(), self.ui.doubleSpinBox_x_max.value()

    def get_y_limits(self):
        return self.ui.doubleSpinBox_y_min.value(), self.ui.doubleSpinBox_y_max.value()


class LineDir(QtGui.QGroupBox):
    def __init__(self, parent, simple_color=True):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Linedir()
        self.ui.setupUi(self)

    def get_linedir(self):
        if self.ui.radioButton_ns.isChecked():
            return 'ns'
        elif self.ui.radioButton_ew.isChecked():
            return 'ew'
        else:
            return None
