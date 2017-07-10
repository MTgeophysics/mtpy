# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\plot_parameters.ui'
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

class Ui_GroupBoxParameters(object):
    def setupUi(self, GroupBoxParameters):
        GroupBoxParameters.setObjectName(_fromUtf8("GroupBoxParameters"))
        GroupBoxParameters.resize(400, 300)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBoxParameters)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBoxZ_Component = QtGui.QGroupBox(GroupBoxParameters)
        self.groupBoxZ_Component.setObjectName(_fromUtf8("groupBoxZ_Component"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.groupBoxZ_Component)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.checkBox_zyx = QtGui.QCheckBox(self.groupBoxZ_Component)
        self.checkBox_zyx.setChecked(True)
        self.checkBox_zyx.setObjectName(_fromUtf8("checkBox_zyx"))
        self.horizontalLayout.addWidget(self.checkBox_zyx)
        self.checkBox_zxy = QtGui.QCheckBox(self.groupBoxZ_Component)
        self.checkBox_zxy.setObjectName(_fromUtf8("checkBox_zxy"))
        self.horizontalLayout.addWidget(self.checkBox_zxy)
        self.checkBox_det = QtGui.QCheckBox(self.groupBoxZ_Component)
        self.checkBox_det.setObjectName(_fromUtf8("checkBox_det"))
        self.horizontalLayout.addWidget(self.checkBox_det)
        self.verticalLayout.addWidget(self.groupBoxZ_Component)
        self.groupBoxFrequency_pereiod = QtGui.QGroupBox(GroupBoxParameters)
        self.groupBoxFrequency_pereiod.setObjectName(_fromUtf8("groupBoxFrequency_pereiod"))
        self.verticalLayoutFrequencyPeriod = QtGui.QVBoxLayout(self.groupBoxFrequency_pereiod)
        self.verticalLayoutFrequencyPeriod.setObjectName(_fromUtf8("verticalLayoutFrequencyPeriod"))
        self.comboBoxPeriod = QtGui.QComboBox(self.groupBoxFrequency_pereiod)
        self.comboBoxPeriod.setEditable(True)
        self.comboBoxPeriod.setObjectName(_fromUtf8("comboBoxPeriod"))
        self.verticalLayoutFrequencyPeriod.addWidget(self.comboBoxPeriod)
        self.verticalLayout.addWidget(self.groupBoxFrequency_pereiod)

        self.retranslateUi(GroupBoxParameters)
        QtCore.QMetaObject.connectSlotsByName(GroupBoxParameters)

    def retranslateUi(self, GroupBoxParameters):
        GroupBoxParameters.setTitle(_translate("GroupBoxParameters", "Parameters", None))
        self.groupBoxZ_Component.setTitle(_translate("GroupBoxParameters", "Z-Component", None))
        self.checkBox_zyx.setText(_translate("GroupBoxParameters", "det", None))
        self.checkBox_zxy.setText(_translate("GroupBoxParameters", "zxy", None))
        self.checkBox_det.setText(_translate("GroupBoxParameters", "zyx", None))
        self.groupBoxFrequency_pereiod.setTitle(_translate("GroupBoxParameters", "Frequency Period", None))

