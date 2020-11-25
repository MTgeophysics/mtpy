from qtpy.QtWidgets import QGroupBox

# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 20/06/2017
"""
from mtpy.gui.SmartMT.ui_asset.groupbox_aspect import Ui_GroupBox_aspect


class AspectRatio(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_aspect()
        self.ui.setupUi(self)

        # connect signal
        self.ui.radioButton_aspect_float.toggled.connect(self._aspect_float_toggled)

    def _aspect_float_toggled(self, checked):
        self.ui.doubleSpinBox_aspect_float.setEnabled(checked)

    def get_aspect(self):
        if self.ui.radioButton_aspect_auto.isChecked():
            return "auto"
        elif self.ui.radioButton_aspect_equal.isChecked():
            return "equal"
        else:
            return self.ui.doubleSpinBox_aspect_float.value()
