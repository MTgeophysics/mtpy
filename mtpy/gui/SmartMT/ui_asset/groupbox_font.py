# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_font.ui'
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
        self.gridLayout = QtGui.QGridLayout(GroupBox_Font)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.checkBox_weight = QtGui.QCheckBox(GroupBox_Font)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_weight.sizePolicy().hasHeightForWidth())
        self.checkBox_weight.setSizePolicy(sizePolicy)
        self.checkBox_weight.setObjectName(_fromUtf8("checkBox_weight"))
        self.gridLayout.addWidget(self.checkBox_weight, 1, 0, 1, 1)
        self.checkBox_size = QtGui.QCheckBox(GroupBox_Font)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_size.sizePolicy().hasHeightForWidth())
        self.checkBox_size.setSizePolicy(sizePolicy)
        self.checkBox_size.setObjectName(_fromUtf8("checkBox_size"))
        self.gridLayout.addWidget(self.checkBox_size, 0, 0, 1, 1)
        self.spinBox_size = QtGui.QSpinBox(GroupBox_Font)
        self.spinBox_size.setEnabled(False)
        self.spinBox_size.setMinimum(4)
        self.spinBox_size.setProperty("value", 10)
        self.spinBox_size.setObjectName(_fromUtf8("spinBox_size"))
        self.gridLayout.addWidget(self.spinBox_size, 0, 2, 1, 1)
        self.checkBox_color = QtGui.QCheckBox(GroupBox_Font)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_color.sizePolicy().hasHeightForWidth())
        self.checkBox_color.setSizePolicy(sizePolicy)
        self.checkBox_color.setObjectName(_fromUtf8("checkBox_color"))
        self.gridLayout.addWidget(self.checkBox_color, 2, 0, 1, 1)
        self.comboBox_size = QtGui.QComboBox(GroupBox_Font)
        self.comboBox_size.setEnabled(False)
        self.comboBox_size.setObjectName(_fromUtf8("comboBox_size"))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.comboBox_size, 0, 1, 1, 1)
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
        self.gridLayout.addWidget(self.comboBox_color, 2, 1, 1, 2)
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
        self.gridLayout.addWidget(self.comboBox_weight, 1, 1, 1, 2)

        self.retranslateUi(GroupBox_Font)
        self.comboBox_size.setCurrentIndex(2)
        self.comboBox_weight.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Font)
        GroupBox_Font.setTabOrder(self.checkBox_size, self.comboBox_size)
        GroupBox_Font.setTabOrder(self.comboBox_size, self.spinBox_size)
        GroupBox_Font.setTabOrder(self.spinBox_size, self.checkBox_weight)
        GroupBox_Font.setTabOrder(self.checkBox_weight, self.comboBox_weight)
        GroupBox_Font.setTabOrder(self.comboBox_weight, self.checkBox_color)
        GroupBox_Font.setTabOrder(self.checkBox_color, self.comboBox_color)

    def retranslateUi(self, GroupBox_Font):
        GroupBox_Font.setTitle(_translate("GroupBox_Font", "Font", None))
        self.checkBox_weight.setText(_translate("GroupBox_Font", "Weight", None))
        self.checkBox_size.setText(_translate("GroupBox_Font", "Size", None))
        self.spinBox_size.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font size</p></body></html>", None))
        self.spinBox_size.setSuffix(_translate("GroupBox_Font", "points", None))
        self.checkBox_color.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font color</p></body></html>", None))
        self.checkBox_color.setText(_translate("GroupBox_Font", "Color", None))
        self.comboBox_size.setItemText(0, _translate("GroupBox_Font", "xx-small", None))
        self.comboBox_size.setItemText(1, _translate("GroupBox_Font", "x-small", None))
        self.comboBox_size.setItemText(2, _translate("GroupBox_Font", "small", None))
        self.comboBox_size.setItemText(3, _translate("GroupBox_Font", "medium", None))
        self.comboBox_size.setItemText(4, _translate("GroupBox_Font", "large", None))
        self.comboBox_size.setItemText(5, _translate("GroupBox_Font", "x-large", None))
        self.comboBox_size.setItemText(6, _translate("GroupBox_Font", "xx-large", None))
        self.comboBox_size.setItemText(7, _translate("GroupBox_Font", "size in points", None))
        self.comboBox_color.setToolTip(_translate("GroupBox_Font", "<html><head/><body><p>Font color</p></body></html>", None))
        self.comboBox_color.setItemText(0, _translate("GroupBox_Font", "Blue", None))
        self.comboBox_color.setItemText(1, _translate("GroupBox_Font", "Green", None))
        self.comboBox_color.setItemText(2, _translate("GroupBox_Font", "Red", None))
        self.comboBox_color.setItemText(3, _translate("GroupBox_Font", "Cyan", None))
        self.comboBox_color.setItemText(4, _translate("GroupBox_Font", "Magenta", None))
        self.comboBox_color.setItemText(5, _translate("GroupBox_Font", "Yellow", None))
        self.comboBox_color.setItemText(6, _translate("GroupBox_Font", "Black", None))
        self.comboBox_color.setItemText(7, _translate("GroupBox_Font", "White", None))
        self.comboBox_weight.setItemText(0, _translate("GroupBox_Font", "light", None))
        self.comboBox_weight.setItemText(1, _translate("GroupBox_Font", "normal", None))
        self.comboBox_weight.setItemText(2, _translate("GroupBox_Font", "medium", None))
        self.comboBox_weight.setItemText(3, _translate("GroupBox_Font", "semibold", None))
        self.comboBox_weight.setItemText(4, _translate("GroupBox_Font", "bold", None))
        self.comboBox_weight.setItemText(5, _translate("GroupBox_Font", "heavy", None))
        self.comboBox_weight.setItemText(6, _translate("GroupBox_Font", "black", None))

