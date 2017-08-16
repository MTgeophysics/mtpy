# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_frequency_period_index.ui'
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

class Ui_GroupBox_Frequency_Period_Index(object):
    def setupUi(self, GroupBox_Frequency_Period_Index):
        GroupBox_Frequency_Period_Index.setObjectName(_fromUtf8("GroupBox_Frequency_Period_Index"))
        GroupBox_Frequency_Period_Index.resize(300, 225)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Frequency_Period_Index.sizePolicy().hasHeightForWidth())
        GroupBox_Frequency_Period_Index.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_Frequency_Period_Index)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.listWidget_frequency_period = QtGui.QListWidget(GroupBox_Frequency_Period_Index)
        self.listWidget_frequency_period.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.listWidget_frequency_period.setObjectName(_fromUtf8("listWidget_frequency_period"))
        self.verticalLayout.addWidget(self.listWidget_frequency_period)

        self.retranslateUi(GroupBox_Frequency_Period_Index)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Frequency_Period_Index)

    def retranslateUi(self, GroupBox_Frequency_Period_Index):
        GroupBox_Frequency_Period_Index.setWindowTitle(_translate("GroupBox_Frequency_Period_Index", "GroupBox", None))
        GroupBox_Frequency_Period_Index.setTitle(_translate("GroupBox_Frequency_Period_Index", "Frequency/Period", None))

