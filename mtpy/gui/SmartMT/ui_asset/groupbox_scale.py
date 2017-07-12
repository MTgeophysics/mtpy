# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_scale.ui'
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

class Ui_GroupBox_Scale(object):
    def setupUi(self, GroupBox_Scale):
        GroupBox_Scale.setObjectName(_fromUtf8("GroupBox_Scale"))
        GroupBox_Scale.resize(300, 79)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Scale.sizePolicy().hasHeightForWidth())
        GroupBox_Scale.setSizePolicy(sizePolicy)
        self.formLayout = QtGui.QFormLayout(GroupBox_Scale)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(GroupBox_Scale)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.comboBox_time = QtGui.QComboBox(GroupBox_Scale)
        self.comboBox_time.setObjectName(_fromUtf8("comboBox_time"))
        self.comboBox_time.addItem(_fromUtf8(""))
        self.comboBox_time.addItem(_fromUtf8(""))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.comboBox_time)
        self.label_2 = QtGui.QLabel(GroupBox_Scale)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.comboBox_map = QtGui.QComboBox(GroupBox_Scale)
        self.comboBox_map.setObjectName(_fromUtf8("comboBox_map"))
        self.comboBox_map.addItem(_fromUtf8(""))
        self.comboBox_map.addItem(_fromUtf8(""))
        self.comboBox_map.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.comboBox_map)

        self.retranslateUi(GroupBox_Scale)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Scale)

    def retranslateUi(self, GroupBox_Scale):
        GroupBox_Scale.setTitle(_translate("GroupBox_Scale", "Scale", None))
        self.label.setToolTip(_translate("GroupBox_Scale", "<html><head/><body><p><span style=\" font-weight:600;\">Period:</span> plot vertical scale in period</p><p><span style=\" font-weight:600;\">Frequency:</span> plot vertical scale in frequency</p></body></html>", None))
        self.label.setText(_translate("GroupBox_Scale", "Time Scale", None))
        self.comboBox_time.setToolTip(_translate("GroupBox_Scale", "<html><head/><body><p><span style=\" font-weight:600;\">Period:</span> plot vertical scale in period</p><p><span style=\" font-weight:600;\">Frequency:</span> plot vertical scale in frequency</p></body></html>", None))
        self.comboBox_time.setItemText(0, _translate("GroupBox_Scale", "Period", None))
        self.comboBox_time.setItemText(1, _translate("GroupBox_Scale", "Frequency", None))
        self.label_2.setToolTip(_translate("GroupBox_Scale", "<html><head/><body><p>Scale of the map coordinates.</p></body></html>", None))
        self.label_2.setText(_translate("GroupBox_Scale", "Map Scale", None))
        self.comboBox_map.setToolTip(_translate("GroupBox_Scale", "<html><head/><body><p>Scale of the map coordinates.</p></body></html>", None))
        self.comboBox_map.setItemText(0, _translate("GroupBox_Scale", "Degress in latitude and longitude", None))
        self.comboBox_map.setItemText(1, _translate("GroupBox_Scale", "Meters for easting and northing", None))
        self.comboBox_map.setItemText(2, _translate("GroupBox_Scale", "Kilometers for easting and northing", None))

