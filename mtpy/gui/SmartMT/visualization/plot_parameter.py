# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui
import numpy as np

from mtpy.gui.SmartMT.ui_asset.plot_parameters import Ui_GroupBoxParameters
from mtpy.gui.SmartMT.visualization.matplotlib_imabedding import MPLCanvas


class PlotParameter(QtGui.QGroupBox):
    _slider_tick_size = 10.0
    _slider_min = 0.0
    _slider_max = 100.0

    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self._mt_objs = None
        self.ui = Ui_GroupBoxParameters()
        self.ui.setupUi(self)
        self._period_histogram = PlotParameter.PeriodHistogram()
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

    def update_period_text(self, value):
        self.ui.comboBoxPeriod.setEditText("%.5f" % ((value * self._slider_tick_size) + self._slider_min))

    def update_period_histogram(self):
        # value = self.ui.doubleSpinBoxPeriod.value()
        value = float(self.ui.comboBoxPeriod.currentText())
        self._period_histogram.set_current_period(value)
        # self.ui.horizontalSliderPeriod.setValue((value - self._slider_min) / self._slider_tick_size)

    def update_ui(self):
        self._update_period()

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs

    def _set_period_tick_size(self, size):
        self._slider_tick_size = size

    def _get_period_tick_size(self):
        return self._slider_tick_size

    def _set_period_min(self, mini):
        self._slider_min = mini
        # self.ui.label_period_min.setText("%.5f" % mini)

    def _get_period_min(self):
        return self._slider_min

    def _set_period_max(self, maxi):
        self._slider_max = maxi
        # self.ui.label_period_max.setText("%.5f" % maxi)

    def _get_period_max(self):
        return self._slider_max

    def _update_period(self):
        if self._mt_objs is not None:
            all_freqs = []
            for mt_obj in self._mt_objs:
                all_freqs.extend(list(mt_obj.Z.freq))
            all_periods = 1.0/np.array(all_freqs)
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

    slider_tick_size = property(_get_period_tick_size, _set_period_tick_size, doc="tick size of period slider")
    slider_min = property(_get_period_min, _set_period_min, doc="minimum on slider")
    slider_max = property(_get_period_max, _set_period_max, doc="maximum on slider")

    class PeriodHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, hight=4, dpi=100):
            self.artists = dict()
            self._periods = None
            self._current_period = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self._lx = None
            self._cursor_x = None
            self._cursor_text = None
            self.mpl_connect('motion_notify_event', self.mouse_move)
            self.mpl_connect('button_release_event', self.mouse_pick)

        def mouse_move(self, event):
            if not event.inaxes:
                return
            x = event.xdata
            y = event.ydata
            if self._cursor_x is None:
                self._cursor_x = self._axes.axvline(linewidth=1, color="green")
            if self._cursor_text is None:
                self._cursor_text = self._axes.text(0.0, 0.0, '', fontsize=8)
            self._cursor_x.set_xdata(x)
            self._cursor_text.set_text('period=%.2f' % x)
            self._cursor_text.set_position((x, y))
            self.draw()

        def mouse_pick(self, event):
            if not event.inaxes:
                return
            x = event.xdata
            self.set_current_period(x)

        def compute_initial_figure(self):
            if self._periods is not None:
                self._axes.hist(self._periods, 50, normed=1)

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
