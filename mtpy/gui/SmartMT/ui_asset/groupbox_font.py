# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_font.ui'
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

class Ui_GroupBox_Font(object):
    def setupUi(self, GroupBox_Font):
        GroupBox_Font.setObjectName(_fromUtf8("GroupBox_Font"))
        GroupBox_Font.resize(300, 105)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Font.sizePolicy().hasHeightForWidth())
        GroupBox_Font.setSizePolicy(sizePolicy)
        GroupBox_Font.setCheckable(False)
        GroupBox_Font.setChecked(False)
        self.formLayout = QtGui.QFormLayout(GroupBox_Font)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.checkBox_size = QtGui.QCheckBox(GroupBox_Font)
        self.checkBox_size.setObjectName(_fromUtf8("checkBox_size"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.checkBox_size)
        self.spinBox_size = QtGui.QSpinBox(GroupBox_Font)
        self.spinBox_size.setEnabled(False)
        self.spinBox_size.setMinimum(4)
        self.spinBox_size.setProperty("value", 10)
        self.spinBox_size.setObjectName(_fromUtf8("spinBox_size"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.spinBox_size)
        self.checkBox_weight = QtGui.QCheckBox(GroupBox_Font)
        self.checkBox_weight.setObjectName(_fromUtf8("checkBox_weight"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.checkBox_weight)
        self.comboBox_weight = QtGui.QComboBox(GroupBox_Font)
        self.comboBox_weight.setEnabled(False)
        self.comboBox_weight.setObjectName(_fromUtf8("comboBox_weight"))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.comboBox_weight)
        self.checkBox_color = QtGui.QCheckBox(GroupBox_Font)
        self.checkBox_color.setObjectName(_fromUtf8("checkBox_color"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.checkBox_color)
        self.comboBox_color = QtGui.QComboBox(GroupBox_Font)
        self.comboBox_color.setEnabled(False)
        self.comboBox_color.setObjectName(_fromUtf8("comboBox_color"))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.comboBox_color.addItem(_fromUtf8(""))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.comboBox_color)

        self.retranslateUi(GroupBox_Font)
        self.comboBox_weight.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Font)

    def retranslateUi(self, GroupBox_Font):
        GroupBox_Font.setTitle(_translate("GroupBox_Font", "Font", None))
        self.checkBox_size.setText(_translate("GroupBox_Font", "Size", None))
        self.spinBox_size.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font size</p></body></html>", None))
        self.checkBox_weight.setText(_translate("GroupBox_Font", "Weight", None))
        self.comboBox_weight.setItemText(0, _translate("GroupBox_Font", "light", None))
        self.comboBox_weight.setItemText(1, _translate("GroupBox_Font", "normal", None))
        self.comboBox_weight.setItemText(2, _translate("GroupBox_Font", "medium", None))
        self.comboBox_weight.setItemText(3, _translate("GroupBox_Font", "semibold", None))
        self.comboBox_weight.setItemText(4, _translate("GroupBox_Font", "bold", None))
        self.comboBox_weight.setItemText(5, _translate("GroupBox_Font", "heavy", None))
        self.comboBox_weight.setItemText(6, _translate("GroupBox_Font", "black", None))
        self.checkBox_color.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font color</p></body></html>", None))
        self.checkBox_color.setText(_translate("GroupBox_Font", "Color", None))
        self.comboBox_color.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font color</p></body></html>", None))
        self.comboBox_color.setItemText(0, _translate("GroupBox_Font", "Blue", None))
        self.comboBox_color.setItemText(1, _translate("GroupBox_Font", "Green", None))
        self.comboBox_color.setItemText(2, _translate("GroupBox_Font", "Red", None))
        self.comboBox_color.setItemText(3, _translate("GroupBox_Font", "Cyan", None))
        self.comboBox_color.setItemText(4, _translate("GroupBox_Font", "Magenta", None))
        self.comboBox_color.setItemText(5, _translate("GroupBox_Font", "Yellow", None))
        self.comboBox_color.setItemText(6, _translate("GroupBox_Font", "Black", None))
        self.comboBox_color.setItemText(7, _translate("GroupBox_Font", "White", None))

