# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
import numpy as np
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.ui_asset.groupbox_color_bar import Ui_GroupBox_ColorBar
from mtpy.gui.SmartMT.ui_asset.groupbox_ellipse import Ui_GroupBoxEllipse
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_single import Ui_groupBoxFrequency_pereiod_single
from mtpy.gui.SmartMT.ui_asset.groupbox_tolerance import Ui_GroupBoxTolerance
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_multiple import Ui_groupBoxZ_Component_Multiple
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_single import Ui_groupBoxZ_Component_Single
from mtpy.gui.SmartMT.ui_asset.plot_parameters import Ui_GroupBoxParameters


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

    def __init__(self, parent, unit="Hz", distribution='Frequency', inverse=False):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self._invert_frequency = inverse
        self._unit = unit
        self._distribuction = distribution
        self.ui = Ui_groupBoxFrequency_pereiod_single()
        self.ui.setupUi(self)
        self._histogram = FrequencySingle.FrequencyHistogram(unit=self._unit, distribution=self._distribuction)
        # add matplotlib canvas
        self.ui.verticalLayoutFrequencyPeriod.addWidget(self._histogram)
        # connect components
        # self.ui.horizontalSliderPeriod.valueChanged.connect(lambda value: self.update_period_text(value))
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_histogram)
        self.ui.comboBoxPeriod.editTextChanged.connect(self.update_histogram)
        self._histogram.mpl_connect('button_release_event', self._mouse_pick)

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

    def _update_frequency(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))

            if self._invert_frequency:
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
        def __init__(self, parent=None, width=5, hight=3, dpi=100, unit="Hz", distribution='Frequency'):
            self.artists = dict()
            self._frequency = None
            self._current_period = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self._lx = None
            self._unit = unit
            self._distribution = distribution
            self.cursor = Cursor(self._axes, track_y=False, text_format="%f " + self._unit, useblit=True)
            # self.cursor = Cursor(self._axes, useblit=True, color='green', linewidth=1)
            # self._cursor_x = None
            # self._cursor_text = None
            # self.mpl_connect('motion_notify_event', self.mouse_move)

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
                self._axes.set_xlabel("%s (%s)" % (self._distribution, self._unit), fontsize=8)
                self.figure.suptitle('%s Distribution in Selected Stations' % self._distribution, fontsize=8)

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
    _colorby = ['phimin', 'phimax', 'skew', 'skew_seg', 'phidet', 'ellipticity']
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
