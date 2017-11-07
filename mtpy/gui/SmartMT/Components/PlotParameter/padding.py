# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_padding import Ui_GroupBox_Padding


class Padding(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Padding()
        self.ui.setupUi(self)

    def get_x_pad(self):
        return self.ui.doubleSpinBox_x.value()

    def get_y_pad(self):
        return self.ui.doubleSpinBox_y.value()
