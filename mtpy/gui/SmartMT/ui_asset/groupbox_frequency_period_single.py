# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_frequency_period_single.ui'
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

class Ui_groupBoxFrequency_pereiod_single(object):
    def setupUi(self, groupBoxFrequency_pereiod_single):
        groupBoxFrequency_pereiod_single.setObjectName(_fromUtf8("groupBoxFrequency_pereiod_single"))
        groupBoxFrequency_pereiod_single.resize(400, 53)
        self.verticalLayoutFrequencyPeriod = QtGui.QVBoxLayout(groupBoxFrequency_pereiod_single)
        self.verticalLayoutFrequencyPeriod.setObjectName(_fromUtf8("verticalLayoutFrequencyPeriod"))
        self.comboBoxPeriod = QtGui.QComboBox(groupBoxFrequency_pereiod_single)
        self.comboBoxPeriod.setObjectName(_fromUtf8("comboBoxPeriod"))
        self.verticalLayoutFrequencyPeriod.addWidget(self.comboBoxPeriod)

        self.retranslateUi(groupBoxFrequency_pereiod_single)
        QtCore.QMetaObject.connectSlotsByName(groupBoxFrequency_pereiod_single)

    def retranslateUi(self, groupBoxFrequency_pereiod_single):
        groupBoxFrequency_pereiod_single.setWindowTitle(_translate("groupBoxFrequency_pereiod_single", "GroupBox", None))
        groupBoxFrequency_pereiod_single.setTitle(_translate("groupBoxFrequency_pereiod_single", "Frequency Period", None))

