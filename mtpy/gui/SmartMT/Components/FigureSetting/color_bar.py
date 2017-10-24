# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_color_bar import Ui_GroupBox_ColorBar


class ColorBar(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_ColorBar()
        self.ui.setupUi(self)

        # connect event
        self.ui.horizontalSlider_x.valueChanged.connect(self._x_slider_value_changed)
        self.ui.horizontalSlider_y.valueChanged.connect(self._y_slider_value_changed)
        self.ui.horizontalSlider_width.valueChanged.connect(self._width_slider_value_changed)
        self.ui.horizontalSlider_height.valueChanged.connect(self._height_slider_value_changed)
        self.ui.comboBox_orientation.currentIndexChanged.connect(self._orientation_changed)
        self.ui.doubleSpinBox_x.editingFinished.connect(self._update_slider_x)
        self.ui.doubleSpinBox_y.editingFinished.connect(self._update_slider_y)
        self.ui.doubleSpinBox_width.editingFinished.connect(self._update_slider_width)
        self.ui.doubleSpinBox_height.editingFinished.connect(self._update_slider_height)

    def _width_slider_value_changed(self, value):
        self.ui.doubleSpinBox_width.setValue(value / 100.0)

    def _height_slider_value_changed(self, value):
        self.ui.doubleSpinBox_height.setValue(value / 100.0)

    def _x_slider_value_changed(self, value):
        self.ui.doubleSpinBox_x.setValue(value / 100.0)

    def _y_slider_value_changed(self, value):
        self.ui.doubleSpinBox_y.setValue(value / 100.0)

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
