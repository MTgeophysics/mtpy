# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\plot_parameters.ui'
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

class Ui_GroupBoxParameters(object):
    def setupUi(self, GroupBoxParameters):
        GroupBoxParameters.setObjectName(_fromUtf8("GroupBoxParameters"))
        GroupBoxParameters.resize(340, 200)
        self.gridLayout = QtGui.QGridLayout(GroupBoxParameters)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.scrollArea = QtGui.QScrollArea(GroupBoxParameters)
        self.scrollArea.setFrameShape(QtGui.QFrame.Panel)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName(_fromUtf8("scrollArea"))
        self.scrollAreaWidgetContents = QtGui.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 318, 136))
        self.scrollAreaWidgetContents.setObjectName(_fromUtf8("scrollAreaWidgetContents"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_2.setMargin(0)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.gridLayout.addWidget(self.scrollArea, 0, 0, 1, 1)
        self.pushButtonPlot = QtGui.QPushButton(GroupBoxParameters)
        self.pushButtonPlot.setObjectName(_fromUtf8("pushButtonPlot"))
        self.gridLayout.addWidget(self.pushButtonPlot, 1, 0, 1, 1)

        self.retranslateUi(GroupBoxParameters)
        QtCore.QMetaObject.connectSlotsByName(GroupBoxParameters)

    def retranslateUi(self, GroupBoxParameters):
        GroupBoxParameters.setTitle(_translate("GroupBoxParameters", "Parameters", None))
        self.pushButtonPlot.setText(_translate("GroupBoxParameters", "Create Figure", None))

