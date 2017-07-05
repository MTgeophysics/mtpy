# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\station_status.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_StationStatus(object):
    def setupUi(self, StationStatus):
        StationStatus.setObjectName(_fromUtf8("StationStatus"))
        StationStatus.resize(400, 300)
        self.horizontalLayout = QtGui.QHBoxLayout(StationStatus)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.tabWidget = QtGui.QTabWidget(StationStatus)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab_summary = QtGui.QWidget()
        self.tab_summary.setObjectName(_fromUtf8("tab_summary"))
        self.tabWidget.addTab(self.tab_summary, _fromUtf8(""))
        self.tab_detial = QtGui.QWidget()
        self.tab_detial.setObjectName(_fromUtf8("tab_detial"))
        self.tabWidget.addTab(self.tab_detial, _fromUtf8(""))
        self.horizontalLayout.addWidget(self.tabWidget)

        self.retranslateUi(StationStatus)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(StationStatus)

    def retranslateUi(self, StationStatus):
        StationStatus.setWindowTitle(_translate("StationStatus", "Status", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_summary), _translate("StationStatus", "Station Summary", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_detial), _translate("StationStatus", "Station Detial", None))

