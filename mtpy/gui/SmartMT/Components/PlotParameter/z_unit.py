# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_z_unit import Ui_GroupBox_z_unit


class ZUnit(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_z_unit()
        self.ui.setupUi(self)

    def get_unit(self):
        return 'km' if self.ui.radioButton_km.isChecked() else 'm'
