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
        GroupBoxParameters.resize(400, 300)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBoxParameters)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.pushButtonPlot = QtGui.QPushButton(GroupBoxParameters)
        self.pushButtonPlot.setObjectName(_fromUtf8("pushButtonPlot"))
        self.verticalLayout.addWidget(self.pushButtonPlot)

        self.retranslateUi(GroupBoxParameters)
        QtCore.QMetaObject.connectSlotsByName(GroupBoxParameters)

    def retranslateUi(self, GroupBoxParameters):
        GroupBoxParameters.setTitle(_translate("GroupBoxParameters", "Parameters", None))
        self.pushButtonPlot.setText(_translate("GroupBoxParameters", "Create Figure", None))

