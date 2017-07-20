# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_z_component_multiple.ui'
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

class Ui_groupBoxZ_Component_Multiple(object):
    def setupUi(self, groupBoxZ_Component_Multiple):
        groupBoxZ_Component_Multiple.setObjectName(_fromUtf8("groupBoxZ_Component_Multiple"))
        groupBoxZ_Component_Multiple.resize(300, 50)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(groupBoxZ_Component_Multiple.sizePolicy().hasHeightForWidth())
        groupBoxZ_Component_Multiple.setSizePolicy(sizePolicy)
        self.horizontalLayout = QtGui.QHBoxLayout(groupBoxZ_Component_Multiple)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.checkBox_det = QtGui.QCheckBox(groupBoxZ_Component_Multiple)
        self.checkBox_det.setObjectName(_fromUtf8("checkBox_det"))
        self.horizontalLayout.addWidget(self.checkBox_det)
        self.checkBox_zxy = QtGui.QCheckBox(groupBoxZ_Component_Multiple)
        self.checkBox_zxy.setObjectName(_fromUtf8("checkBox_zxy"))
        self.horizontalLayout.addWidget(self.checkBox_zxy)
        self.checkBox_zyx = QtGui.QCheckBox(groupBoxZ_Component_Multiple)
        self.checkBox_zyx.setObjectName(_fromUtf8("checkBox_zyx"))
        self.horizontalLayout.addWidget(self.checkBox_zyx)

        self.retranslateUi(groupBoxZ_Component_Multiple)
        QtCore.QMetaObject.connectSlotsByName(groupBoxZ_Component_Multiple)

    def retranslateUi(self, groupBoxZ_Component_Multiple):
        groupBoxZ_Component_Multiple.setWindowTitle(_translate("groupBoxZ_Component_Multiple", "GroupBox", None))
        groupBoxZ_Component_Multiple.setTitle(_translate("groupBoxZ_Component_Multiple", "Z-Component", None))
        self.checkBox_det.setText(_translate("groupBoxZ_Component_Multiple", "det", None))
        self.checkBox_zxy.setText(_translate("groupBoxZ_Component_Multiple", "zxy", None))
        self.checkBox_zyx.setText(_translate("groupBoxZ_Component_Multiple", "zyx", None))

