# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_z_component_single.ui'
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

class Ui_groupBoxZ_Component_Single(object):
    def setupUi(self, groupBoxZ_Component_Single):
        groupBoxZ_Component_Single.setObjectName(_fromUtf8("groupBoxZ_Component_Single"))
        groupBoxZ_Component_Single.resize(300, 50)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(groupBoxZ_Component_Single.sizePolicy().hasHeightForWidth())
        groupBoxZ_Component_Single.setSizePolicy(sizePolicy)
        self.horizontalLayout = QtGui.QHBoxLayout(groupBoxZ_Component_Single)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.radioButton_det = QtGui.QRadioButton(groupBoxZ_Component_Single)
        self.radioButton_det.setChecked(True)
        self.radioButton_det.setObjectName(_fromUtf8("radioButton_det"))
        self.horizontalLayout.addWidget(self.radioButton_det)
        self.radioButton_zxy = QtGui.QRadioButton(groupBoxZ_Component_Single)
        self.radioButton_zxy.setObjectName(_fromUtf8("radioButton_zxy"))
        self.horizontalLayout.addWidget(self.radioButton_zxy)
        self.radioButton_zyx = QtGui.QRadioButton(groupBoxZ_Component_Single)
        self.radioButton_zyx.setObjectName(_fromUtf8("radioButton_zyx"))
        self.horizontalLayout.addWidget(self.radioButton_zyx)

        self.retranslateUi(groupBoxZ_Component_Single)
        QtCore.QMetaObject.connectSlotsByName(groupBoxZ_Component_Single)

    def retranslateUi(self, groupBoxZ_Component_Single):
        groupBoxZ_Component_Single.setWindowTitle(_translate("groupBoxZ_Component_Single", "GroupBox", None))
        groupBoxZ_Component_Single.setTitle(_translate("groupBoxZ_Component_Single", "Z-Component", None))
        self.radioButton_det.setText(_translate("groupBoxZ_Component_Single", "det", None))
        self.radioButton_zxy.setText(_translate("groupBoxZ_Component_Single", "zxy", None))
        self.radioButton_zyx.setText(_translate("groupBoxZ_Component_Single", "zyx", None))

