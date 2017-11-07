# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_stretch import Ui_GroupBox_Stretch


class Stretch(QGroupBox):
    def __init__(self, parent, simple_color=True):
        QGroupBox.__init__(self, parent)
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
