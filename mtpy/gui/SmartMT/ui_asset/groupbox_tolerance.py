# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_tolerance.ui'
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

class Ui_GroupBoxTolerance(object):
    def setupUi(self, GroupBoxTolerance):
        GroupBoxTolerance.setObjectName(_fromUtf8("GroupBoxTolerance"))
        GroupBoxTolerance.resize(300, 53)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBoxTolerance.sizePolicy().hasHeightForWidth())
        GroupBoxTolerance.setSizePolicy(sizePolicy)
        self.formLayout = QtGui.QFormLayout(GroupBoxTolerance)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(GroupBoxTolerance)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.doubleSpinBox = QtGui.QDoubleSpinBox(GroupBoxTolerance)
        self.doubleSpinBox.setMaximum(100.0)
        self.doubleSpinBox.setProperty("value", 10.0)
        self.doubleSpinBox.setObjectName(_fromUtf8("doubleSpinBox"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.doubleSpinBox)

        self.retranslateUi(GroupBoxTolerance)
        QtCore.QMetaObject.connectSlotsByName(GroupBoxTolerance)

    def retranslateUi(self, GroupBoxTolerance):
        GroupBoxTolerance.setTitle(_translate("GroupBoxTolerance", "Frequency Tolerance", None))
        self.label.setToolTip(_translate("GroupBoxTolerance", "<html><head/><body><p>Tolerance in frequency range to look for in each file. The frequency will be \'equal\' if it lay in the tolerance range.</p></body></html>", None))
        self.label.setText(_translate("GroupBoxTolerance", "Tolerance (%)", None))
        self.doubleSpinBox.setToolTip(_translate("GroupBoxTolerance", "<html><head/><body><p>Tolerance in frequency range to look for in each file. The frequency will be \'equal\' if it lay in the tolerance range.</p></body></html>", None))

