# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_z_unit.ui'
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


class Ui_GroupBox_z_unit(object):
    def setupUi(self, GroupBox_z_unit):
        GroupBox_z_unit.setObjectName(_fromUtf8("GroupBox_z_unit"))
        GroupBox_z_unit.resize(300, 50)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_z_unit.sizePolicy().hasHeightForWidth())
        GroupBox_z_unit.setSizePolicy(sizePolicy)
        self.horizontalLayout = QtGui.QHBoxLayout(GroupBox_z_unit)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.radioButton_m = QtGui.QRadioButton(GroupBox_z_unit)
        self.radioButton_m.setObjectName(_fromUtf8("radioButton_m"))
        self.horizontalLayout.addWidget(self.radioButton_m)
        self.radioButton_km = QtGui.QRadioButton(GroupBox_z_unit)
        self.radioButton_km.setChecked(True)
        self.radioButton_km.setObjectName(_fromUtf8("radioButton_km"))
        self.horizontalLayout.addWidget(self.radioButton_km)

        self.retranslateUi(GroupBox_z_unit)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_z_unit)

    def retranslateUi(self, GroupBox_z_unit):
        GroupBox_z_unit.setWindowTitle(_translate("GroupBox_z_unit", "GroupBox", None))
        GroupBox_z_unit.setTitle(_translate("GroupBox_z_unit", "Z Unit", None))
        self.radioButton_m.setText(_translate("GroupBox_z_unit", "meter (m)", None))
        self.radioButton_km.setText(_translate("GroupBox_z_unit", "kilometre (km)", None))
