# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_padding.ui'
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

class Ui_GroupBox_Padding(object):
    def setupUi(self, GroupBox_Padding):
        GroupBox_Padding.setObjectName(_fromUtf8("GroupBox_Padding"))
        GroupBox_Padding.resize(300, 79)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Padding.sizePolicy().hasHeightForWidth())
        GroupBox_Padding.setSizePolicy(sizePolicy)
        GroupBox_Padding.setCheckable(True)
        GroupBox_Padding.setChecked(False)
        self.formLayout = QtGui.QFormLayout(GroupBox_Padding)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(GroupBox_Padding)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.doubleSpinBox_x = QtGui.QDoubleSpinBox(GroupBox_Padding)
        self.doubleSpinBox_x.setProperty("value", 0.2)
        self.doubleSpinBox_x.setObjectName(_fromUtf8("doubleSpinBox_x"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_x)
        self.label_2 = QtGui.QLabel(GroupBox_Padding)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.doubleSpinBox_y = QtGui.QDoubleSpinBox(GroupBox_Padding)
        self.doubleSpinBox_y.setProperty("value", 0.2)
        self.doubleSpinBox_y.setObjectName(_fromUtf8("doubleSpinBox_y"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_y)

        self.retranslateUi(GroupBox_Padding)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Padding)

    def retranslateUi(self, GroupBox_Padding):
        GroupBox_Padding.setWindowTitle(_translate("GroupBox_Padding", "GroupBox", None))
        GroupBox_Padding.setTitle(_translate("GroupBox_Padding", "Padding", None))
        self.label.setToolTip(_translate("GroupBox_Padding", "<html><head/><body><p>Padding in the east-west direction of plot boundaries.</p><p><span style=\" font-weight:600;\">Note:</span> This is by default set to lat and long units.</p></body></html>", None))
        self.label.setText(_translate("GroupBox_Padding", "X (east-west)", None))
        self.doubleSpinBox_x.setToolTip(_translate("GroupBox_Padding", "<html><head/><body><p>Padding in the east-west direction of plot boundaries.</p><p><span style=\" font-weight:600;\">Note:</span> This is by default set to lat and long units.</p></body></html>", None))
        self.label_2.setToolTip(_translate("GroupBox_Padding", "<html><head/><body><p>Padding in the north-south direction of plot boundaries.</p><p><span style=\" font-weight:600;\">Note:</span> This is by default set to lat and long unites.</p></body></html>", None))
        self.label_2.setText(_translate("GroupBox_Padding", "Y (north-south)", None))
        self.doubleSpinBox_y.setToolTip(_translate("GroupBox_Padding", "<html><head/><body><p>Padding in the north-south direction of plot boundaries.</p><p><span style=\" font-weight:600;\">Note:</span> This is by default set to lat and long unites.</p></body></html>", None))

