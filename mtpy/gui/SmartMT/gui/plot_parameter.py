# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
import six
from qtpy import QtCore
from qtpy.QtWidgets import QGroupBox
from matplotlib import colors as mcolors

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


class PlotParameter(QGroupBox):
    _slider_tick_size = 10.0
    _slider_min = 0.0
    _slider_max = 100.0

    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxParameters()
        self.ui.setupUi(self)

        # dpi
        self.ui.spinBox_dpi.valueChanged.connect(self._dpi_changed)
        # inches
        self.ui.doubleSpinBox_width_inches.valueChanged.connect(self._width_inches_changed)
        self.ui.doubleSpinBox_height_inches.valueChanged.connect(self._height_inches_changed)
        # pixels
        self.ui.spinBox_width_pixels.valueChanged.connect(self._width_pixels_changed)
        self.ui.spinBox_height_pixels.valueChanged.connect(self._height_pixels_changed)
        # x
        self.ui.horizontalSlider_x.valueChanged.connect(self._x_slider_changed)
        self.ui.doubleSpinBox_x.valueChanged.connect(self._x_spinbox_changed)
        # y
        self.ui.horizontalSlider_y.valueChanged.connect(self._y_slider_changed)
        self.ui.doubleSpinBox_y.valueChanged.connect(self._y_spinbox_changed)

    def add_parameter_groupbox(self, groupbox):
        self.ui.verticalLayout_plot.addWidget(groupbox, QtCore.Qt.AlignLeft)
        self.resize(self.sizeHint())

    def end_of_parameter_components(self):
        self.ui.verticalLayout_plot.addStretch()
        self.ui.verticalLayout_figure.addStretch()

    def add_figure_groupbox(self, groupbox):
        self.ui.verticalLayout_figure.addWidget(groupbox, QtCore.Qt.AlignLeft)
        self.resize(self.sizeHint())

    def get_title(self):
        return str(self.ui.lineEdit_title.text())

    def set_title(self, title):
        self.ui.lineEdit_title.setText(title)

    def _x_slider_changed(self, value):
        self.ui.doubleSpinBox_x.blockSignals(True)
        self.ui.doubleSpinBox_x.setValue(value / 100.0)
        self.ui.doubleSpinBox_x.blockSignals(False)

    def _y_slider_changed(self, value):
        self.ui.doubleSpinBox_y.blockSignals(True)
        self.ui.doubleSpinBox_y.setValue(value / 100.0)
        self.ui.doubleSpinBox_y.blockSignals(False)

    def _x_spinbox_changed(self, value):
        self.ui.horizontalSlider_x.blockSignals(True)
        self.ui.horizontalSlider_x.setValue(value * 100)
        self.ui.horizontalSlider_x.blockSignals(False)

    def _y_spinbox_changed(self, value):
        self.ui.horizontalSlider_y.blockSignals(True)
        self.ui.horizontalSlider_y.setValue(value * 100)
        self.ui.horizontalSlider_y.blockSignals(False)

    def _dpi_changed(self, dpi):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(
            self.ui.doubleSpinBox_height_inches.value() * dpi
        )
        self.ui.spinBox_width_pixels.setValue(
            self.ui.doubleSpinBox_width_inches.value() * dpi
        )
        self.ui.doubleSpinBox_height_inches.blockSignals(False)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _width_pixels_changed(self, width):
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        new_width_inches = width / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_width_inches.setValue(new_width_inches)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _height_pixels_changed(self, height):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        new_height_inches = height / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_height_inches.setValue(new_height_inches)
        self.ui.doubleSpinBox_height_inches.blockSignals(False)

    def _width_inches_changed(self, width):
        self.ui.spinBox_width_pixels.blockSignals(True)
        self.ui.spinBox_width_pixels.setValue(
            width * self.ui.spinBox_dpi.value()
        )
        self.ui.spinBox_width_pixels.blockSignals(False)

    def _height_inches_changed(self, height):
        self.ui.spinBox_height_pixels.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(
            height * self.ui.spinBox_dpi.value()
        )
        self.ui.spinBox_height_pixels.blockSignals(False)

    def customized_figure_title(self):
        return self.ui.groupBox_title.isChecked()

    def customized_figure_size(self):
        return self.ui.groupBox_figure_size.isChecked()

    def get_size_inches_width(self):
        return self.ui.doubleSpinBox_width_inches.value()

    def get_size_inches_height(self):
        return self.ui.doubleSpinBox_height_inches.value()

    def get_dpi(self):
        return self.ui.spinBox_dpi.value()

    def get_layout(self):
        if self.ui.checkBox_tight_layout.isChecked():
            return True
        else:
            return False

    def get_title_font_dict(self):
        font_properties = {
            'x': self.ui.doubleSpinBox_x.value(),
            'y': self.ui.doubleSpinBox_y.value(),
            'horizontalalignment': self._horizontalalignment[
                self.ui.comboBox_horizontal_alignment.currentIndex()],
            'verticalalignment': self._verticalalignment[
                self.ui.comboBox_vertical_alignment.currentIndex()
            ],
            'fontsize': self.ui.spinBox_fontsize.value()
        }
        return font_properties

    _horizontalalignment = ['right', 'center', 'left']
    _verticalalignment = ['top', 'center', 'bottom', 'baseline']


