# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
import numpy as np
from PyQt4 import QtGui

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas, Cursor
from mtpy.gui.SmartMT.ui_asset.groupbox_frequency_period_single import Ui_groupBoxFrequency_pereiod_single
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


class ZComponentMultiple(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Multiple()
        self.ui.setupUi(self)
        # z-component checkbox logic
        self.ui.checkBox_zyx.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_zxy.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_det.stateChanged.connect(self._multiple_zcomponent_logic)

    def add_parameter_groubox(self, groupbox):
        self.ui.horizontalLayout.addWidget(groupbox)
        self.resize(self.sizeHint())

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
        if self.parameter_ui.ui.radioButton_det.isChecked():
            return 'det'
        elif self.parameter_ui.ui.radioButton_zxy.isChecked():
            return 'zxy'
        elif self.parameter_ui.ui.radioButton_zyx.isChecked():
            return 'zyx'


class FrequencyPeriodSingle(QtGui.QGroupBox):
    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.ui = Ui_groupBoxFrequency_pereiod_single()
        self.ui.setupUi(self)
        self._period_histogram = PlotParameter.PeriodHistogram()
        # add matplotlib canvas
        self.ui.verticalLayoutFrequencyPeriod.addWidget(self._period_histogram)
        # connect components
        # self.ui.horizontalSliderPeriod.valueChanged.connect(lambda value: self.update_period_text(value))
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_period_histogram)
        self.ui.comboBoxPeriod.editTextChanged.connect(self.update_period_histogram)
        # self.ui.doubleSpinBoxPeriod.editingFinished.connect(self.update_period_slider)
        self._period_histogram.mpl_connect('button_release_event', self._mouse_pick)

    def _mouse_pick(self, event):
        if not event.inaxes:
            return
        x = event.xdata
        self.ui.comboBoxPeriod.setEditText("%.5f" % x)

    def update_period_histogram(self):
        # value = self.ui.doubleSpinBoxPeriod.value()
        value = float(self.ui.comboBoxPeriod.currentText())
        self._period_histogram.set_current_period(value)
        # self.ui.horizontalSliderPeriod.setValue((value - self._slider_min) / self._slider_tick_size)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._update_period()

    def _update_period(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))
            all_periods = 1.0 / np.array(all_freqs)
            # self._parameter_ui.slider_min = np.min(all_unique_periods)
            # self._parameter_ui.slider_max = np.max(all_unique_periods)
            # self._parameter_ui.slider_tick_size = (self._parameter_ui.slider_max -
            #                                        self._parameter_ui.slider_min) / 100.0  # 100 is the default ticks of the slidere
            self._period_histogram.set_data(all_periods)
            self._period_histogram.update_figure()
            # sort all frequencies in ascending order
            all_unique_periods = sorted(list(set(all_periods)))
            for period in all_unique_periods:
                self.ui.comboBoxPeriod.addItem("%.5f" % period)
            self.ui.comboBoxPeriod.setCurrentIndex(0)
            self.update_period_histogram()

    class PeriodHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, hight=4, dpi=100):
            self.artists = dict()
            self._periods = None
            self._current_period = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self._lx = None
            self.cursor = Cursor(self._axes, track_y=False, text_format="period=%f", useblit=True)
            # self.cursor = Cursor(self._axes, useblit=True, color='green', linewidth=1)
            # self._cursor_x = None
            # self._cursor_text = None
            # self.mpl_connect('motion_notify_event', self.mouse_move)

            # self.mpl_connect('motion_notify_event', self.cursor)
            self.mpl_connect('button_release_event', self.mouse_pick)
            self.setMinimumSize(200, 50)
            # self.resize(self.minimumSizeHint())

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
            self.set_current_period(x)

        def compute_initial_figure(self):
            if self._periods is not None:
                self._axes.tick_params(axis='both', which='major', labelsize=6)
                self._axes.tick_params(axis='both', which='minor', labelsize=4)
                self._axes.hist(self._periods, 50, normed=1)
                self._axes.set_xlabel('Period', fontsize=8)
                self.figure.suptitle('Period Distribution in Selected Stations', fontsize=8)

        def set_data(self, periods):
            self._periods = periods

        def set_current_period(self, period):
            self._current_period = period
            if self._lx is None:
                self._lx = self._axes.axvline(linewidth=2, color="red")
            self._lx.set_xdata(self._current_period)
            self.draw()

        def update_figure(self):
            # clear figure
            self._axes.cla()
            self.compute_initial_figure()
            self.draw()
