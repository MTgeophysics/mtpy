# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""

from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_multiple import (
    Ui_groupBoxZ_Component_Multiple,
)
from mtpy.gui.SmartMT.ui_asset.groupbox_z_component_single import (
    Ui_groupBoxZ_Component_Single,
)


class ZComponentMultiple(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Multiple()
        self.ui.setupUi(self)
        # z-component checkbox logic
        self.ui.checkBox_zyx.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_zxy.stateChanged.connect(self._multiple_zcomponent_logic)
        self.ui.checkBox_det.stateChanged.connect(self._multiple_zcomponent_logic)

    def _multiple_zcomponent_logic(self, int):
        """
        set up at least one component selected
        :return:
        """
        # counter = 0
        if (
            self.ui.checkBox_det.isChecked()
            + self.ui.checkBox_zxy.isChecked()
            + self.ui.checkBox_zyx.isChecked()
            == 1
        ):
            # only one checkbox is checked, lock the checked box
            if self.ui.checkBox_det.isChecked():
                self.ui.checkBox_det.setEnabled(False)
            elif self.ui.checkBox_zxy.isChecked():
                self.ui.checkBox_zxy.setEnabled(False)
            else:
                self.ui.checkBox_zyx.setEnabled(False)
        else:
            self.ui.checkBox_det.setEnabled(True)
            self.ui.checkBox_zxy.setEnabled(True)
            self.ui.checkBox_zyx.setEnabled(True)

    def get_selection(self):
        zcomponent = []
        if self.ui.checkBox_det.isChecked():
            zcomponent.append("det")
        if self.ui.checkBox_zxy.isChecked():
            zcomponent.append("zxy")
        if self.ui.checkBox_zyx.isChecked():
            zcomponent.append("zyx")
        return zcomponent


class ZComponentSingle(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_groupBoxZ_Component_Single()
        self.ui.setupUi(self)

    def get_selection(self):
        if self.ui.radioButton_det.isChecked():
            return "det"
        elif self.ui.radioButton_zxy.isChecked():
            return "zxy"
        elif self.ui.radioButton_zyx.isChecked():
            return "zyx"
