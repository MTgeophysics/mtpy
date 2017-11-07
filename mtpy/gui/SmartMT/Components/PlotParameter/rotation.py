# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_rotation import Ui_GroupBox_Rotation


class Rotation(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Rotation()
        self.ui.setupUi(self)
        self.ui.dial_rotation.valueChanged.connect(self._dial_value_changed)
        self.ui.doubleSpinBox_rotation.valueChanged.connect(self._text_value_changed)

    def _dial_value_changed(self, p_int):
        degree = (p_int - 180) % 360
        self.ui.doubleSpinBox_rotation.setValue(degree)

    def _text_value_changed(self):
        degree = (int(self.ui.doubleSpinBox_rotation.value()) + 180) % 360
        if degree != self.ui.dial_rotation.value():
            self.ui.dial_rotation.setValue(degree)

    def get_rotation_in_degree(self):
        return self.ui.doubleSpinBox_rotation.value()
