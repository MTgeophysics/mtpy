# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_scale import Ui_GroupBox_Scale


class Scale(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Scale()
        self.ui.setupUi(self)

    _tscale = ["period", "freq"]
    _mapscale = ["deg", "m", "km"]

    def get_tscale(self):
        return self._tscale[self.ui.comboBox_time.currentIndex()]

    def hide_mapscale(self):
        self.ui.label_map.hide()
        self.ui.comboBox_map.hide()

    def get_mapscale(self):
        return self._mapscale[self.ui.comboBox_map.currentIndex()]
