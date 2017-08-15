# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_aspect.ui'
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

class Ui_GroupBox_aspect(object):
    def setupUi(self, GroupBox_aspect):
        GroupBox_aspect.setObjectName(_fromUtf8("GroupBox_aspect"))
        GroupBox_aspect.resize(300, 53)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_aspect.sizePolicy().hasHeightForWidth())
        GroupBox_aspect.setSizePolicy(sizePolicy)
        self.gridLayout_2 = QtGui.QGridLayout(GroupBox_aspect)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.radioButton_aspect_float = QtGui.QRadioButton(GroupBox_aspect)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_aspect_float.sizePolicy().hasHeightForWidth())
        self.radioButton_aspect_float.setSizePolicy(sizePolicy)
        self.radioButton_aspect_float.setText(_fromUtf8(""))
        self.radioButton_aspect_float.setObjectName(_fromUtf8("radioButton_aspect_float"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_float, 0, 2, 1, 1)
        self.radioButton_aspect_auto = QtGui.QRadioButton(GroupBox_aspect)
        self.radioButton_aspect_auto.setChecked(True)
        self.radioButton_aspect_auto.setObjectName(_fromUtf8("radioButton_aspect_auto"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_auto, 0, 0, 1, 1)
        self.doubleSpinBox_aspect_float = QtGui.QDoubleSpinBox(GroupBox_aspect)
        self.doubleSpinBox_aspect_float.setEnabled(False)
        self.doubleSpinBox_aspect_float.setProperty("value", 1.0)
        self.doubleSpinBox_aspect_float.setObjectName(_fromUtf8("doubleSpinBox_aspect_float"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_aspect_float, 0, 3, 1, 1)
        self.radioButton_aspect_equal = QtGui.QRadioButton(GroupBox_aspect)
        self.radioButton_aspect_equal.setObjectName(_fromUtf8("radioButton_aspect_equal"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_equal, 0, 1, 1, 1)

        self.retranslateUi(GroupBox_aspect)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_aspect)

    def retranslateUi(self, GroupBox_aspect):
        GroupBox_aspect.setWindowTitle(_translate("GroupBox_aspect", "GroupBox", None))
        GroupBox_aspect.setTitle(_translate("GroupBox_aspect", "Aspect Ratio", None))
        self.radioButton_aspect_auto.setText(_translate("GroupBox_aspect", "Auto", None))
        self.doubleSpinBox_aspect_float.setToolTip(_translate("GroupBox_aspect", "<html><head/><body><p>height/width</p></body></html>", None))
        self.radioButton_aspect_equal.setText(_translate("GroupBox_aspect", "Equal", None))

