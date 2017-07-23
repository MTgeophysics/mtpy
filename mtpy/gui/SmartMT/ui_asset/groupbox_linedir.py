# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_linedir.ui'
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

class Ui_GroupBox_Linedir(object):
    def setupUi(self, GroupBox_Linedir):
        GroupBox_Linedir.setObjectName(_fromUtf8("GroupBox_Linedir"))
        GroupBox_Linedir.resize(300, 73)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Linedir.sizePolicy().hasHeightForWidth())
        GroupBox_Linedir.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_Linedir)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.radioButton_ns = QtGui.QRadioButton(GroupBox_Linedir)
        self.radioButton_ns.setChecked(True)
        self.radioButton_ns.setObjectName(_fromUtf8("radioButton_ns"))
        self.verticalLayout.addWidget(self.radioButton_ns)
        self.radioButton_ew = QtGui.QRadioButton(GroupBox_Linedir)
        self.radioButton_ew.setObjectName(_fromUtf8("radioButton_ew"))
        self.verticalLayout.addWidget(self.radioButton_ew)

        self.retranslateUi(GroupBox_Linedir)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Linedir)

    def retranslateUi(self, GroupBox_Linedir):
        GroupBox_Linedir.setWindowTitle(_translate("GroupBox_Linedir", "GroupBox", None))
        GroupBox_Linedir.setToolTip(_translate("GroupBox_Linedir", "<html><head/><body><p>Predominant dirrection of profile line</p></body></html>", None))
        GroupBox_Linedir.setTitle(_translate("GroupBox_Linedir", "Line Direction", None))
        self.radioButton_ns.setText(_translate("GroupBox_Linedir", "North-South Line", None))
        self.radioButton_ew.setText(_translate("GroupBox_Linedir", "East-West Line", None))

