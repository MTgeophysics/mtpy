# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui

from mtpy.gui.SmartMT.ui_asset.plot_parameters import Ui_GroupBoxParameters
from mtpy.gui.SmartMT.visualization.matplotlib_imabedding import MPLCanvas


class PlotParameter(QtGui.QGroupBox):
    _slider_tick_size = 10.0
    _slider_min = 0.0
    _slider_max = 100.0

    def __init__(self, parent):
        QtGui.QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxParameters()
        self.ui.setupUi(self)
        # self.period_histogram = PlotParameter.PeriodHistogram()
        # self.ui.verticalLayoutFrequencyPeriod.addWidget(self.period_histogram)
        # connect components
        self.ui.horizontalSliderPeriod.valueChanged.connect(lambda value: self.update_period_text(value))
        self.ui.comboBoxPeriod.currentIndexChanged.connect(self.update_period_slider)
        # self.ui.doubleSpinBoxPeriod.editingFinished.connect(self.update_period_slider)

    def update_period_text(self, value):
        self.ui.comboBoxPeriod.setEditText("%.5f" % ((value * self._slider_tick_size) + self._slider_min))

    def update_period_slider(self):
        # value = self.ui.doubleSpinBoxPeriod.value()
        value = float(self.ui.comboBoxPeriod.currentText())
        self.ui.horizontalSliderPeriod.setValue((value - self._slider_min) / self._slider_tick_size)

    def _set_period_tick_size(self, size):
        self._slider_tick_size = size

    def _get_period_tick_size(self):
        return self._slider_tick_size

    def _set_period_min(self, mini):
        self._slider_min = mini
        self.ui.label_period_min.setText("%.5f" % mini)

    def _get_period_min(self):
        return self._slider_min

    def _set_period_max(self, maxi):
        self._slider_max = maxi
        self.ui.label_period_max.setText("%.5f" % maxi)

    def _get_period_max(self):
        return self._slider_max

    slider_tick_size = property(_get_period_tick_size, _set_period_tick_size, doc="tick size of period slider")
    slider_min = property(_get_period_min, _set_period_min, doc="minimum on slider")
    slider_max = property(_get_period_max, _set_period_max, doc="maximum on slider")

    class PeriodHistogram(MPLCanvas):
        def __init__(self, parent=None, width=5, hight=4, dpi=100):
            self.artists = dict()
            self._periods = None
            MPLCanvas.__init__(self, parent, width, hight, dpi)

        def compute_initial_figure(self):
            if self._periods is not None:
                self._axes.hist(self._periods, 50, normed=1)

        def set_data(self, periods):
            self._periods = periods
            self.update_figure()

        def update_figure(self):
            # clear figure
            self._axes.cla()
            self.compute_initial_figure()
            self.draw()
