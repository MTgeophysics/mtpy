# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_linedir import Ui_GroupBox_Linedir


class LineDir(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Linedir()
        self.ui.setupUi(self)

    def get_linedir(self):
        if self.ui.radioButton_ns.isChecked():
            return 'ns'
        elif self.ui.radioButton_ew.isChecked():
            return 'ew'
        else:
            return None
