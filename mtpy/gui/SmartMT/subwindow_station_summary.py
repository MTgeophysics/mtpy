# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the station summary subwindow

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.PyQt4.station_status import Ui_StationStatus


class StationSummary(QtGui.QWidget):
    def __init__(self, parent, file_handler):
        QtGui.QWidget.__init__(self, parent)
        self.file_handler = file_handler
        self.ui = Ui_StationStatus()
        self.ui.setupUi(self)
        self.subwindow, _ = parent.create_subwindow(self, self.windowTitle())
        self.subwindow.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)

