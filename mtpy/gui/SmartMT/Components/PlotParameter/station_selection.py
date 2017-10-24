# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtCore import Signal
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_station_select import Ui_GroupBox_Station_Select


class StationSelection(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Station_Select()
        self.ui.setupUi(self)
        self.mt_objs = None

        self.ui.comboBox_station.currentIndexChanged.connect(self._current_station_changed)

    def _current_station_changed(self):
        self.station_changed.emit()

    station_changed = Signal()

    def set_data(self, mt_objs):
        self.ui.comboBox_station.clear()
        self.mt_objs = []
        for mt_obj in mt_objs:
            self.mt_objs.append(mt_obj)
            self.ui.comboBox_station.addItem(mt_obj.station)

    def get_station(self):
        index = self.ui.comboBox_station.currentIndex()
        return self.mt_objs[index]
