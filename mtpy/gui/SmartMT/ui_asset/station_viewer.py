# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\station_viewer.ui'
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

class Ui_StationViewer(object):
    def setupUi(self, StationViewer):
        StationViewer.setObjectName(_fromUtf8("StationViewer"))
        StationViewer.resize(400, 600)
        StationViewer.setFocusPolicy(QtCore.Qt.NoFocus)
        self.verticalLayout = QtGui.QVBoxLayout(StationViewer)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.treeWidget_stations = QtGui.QTreeWidget(StationViewer)
        self.treeWidget_stations.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.treeWidget_stations.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.treeWidget_stations.setAnimated(False)
        self.treeWidget_stations.setHeaderHidden(True)
        self.treeWidget_stations.setObjectName(_fromUtf8("treeWidget_stations"))
        self.verticalLayout.addWidget(self.treeWidget_stations)
        self.pushButton_showMap = QtGui.QPushButton(StationViewer)
        self.pushButton_showMap.setObjectName(_fromUtf8("pushButton_showMap"))
        self.verticalLayout.addWidget(self.pushButton_showMap)

        self.retranslateUi(StationViewer)
        QtCore.QMetaObject.connectSlotsByName(StationViewer)

    def retranslateUi(self, StationViewer):
        StationViewer.setWindowTitle(_translate("StationViewer", "Stations", None))
        self.treeWidget_stations.headerItem().setText(0, _translate("StationViewer", "Station", None))
        self.treeWidget_stations.headerItem().setText(1, _translate("StationViewer", "File", None))
        self.pushButton_showMap.setText(_translate("StationViewer", "Show Map", None))

