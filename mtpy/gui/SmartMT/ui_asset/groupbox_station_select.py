# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_station_select.ui'
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


class Ui_GroupBox_Station_Select(object):
    def setupUi(self, GroupBox_Station_Select):
        GroupBox_Station_Select.setObjectName(_fromUtf8("GroupBox_Station_Select"))
        GroupBox_Station_Select.resize(300, 53)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Station_Select.sizePolicy().hasHeightForWidth())
        GroupBox_Station_Select.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_Station_Select)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.comboBox_station = QtGui.QComboBox(GroupBox_Station_Select)
        self.comboBox_station.setObjectName(_fromUtf8("comboBox_station"))
        self.verticalLayout.addWidget(self.comboBox_station)

        self.retranslateUi(GroupBox_Station_Select)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Station_Select)

    def retranslateUi(self, GroupBox_Station_Select):
        GroupBox_Station_Select.setWindowTitle(_translate("GroupBox_Station_Select", "GroupBox", None))
        GroupBox_Station_Select.setToolTip(
            _translate("GroupBox_Station_Select", "<html><head/><body><p>Select a station to plot</p></body></html>",
                       None))
        GroupBox_Station_Select.setTitle(_translate("GroupBox_Station_Select", "Stations", None))
