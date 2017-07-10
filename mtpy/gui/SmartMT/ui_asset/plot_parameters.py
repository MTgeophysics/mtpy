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
        self.horizontalSliderPeriod = QtGui.QSlider(self.groupBoxFrequency_pereiod)
        self.horizontalSliderPeriod.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSliderPeriod.setTickPosition(QtGui.QSlider.TicksAbove)
        self.horizontalSliderPeriod.setObjectName(_fromUtf8("horizontalSliderPeriod"))
        self.verticalLayoutFrequencyPeriod.addWidget(self.horizontalSliderPeriod)
        self.widget = QtGui.QWidget(self.groupBoxFrequency_pereiod)
        self.widget.setObjectName(_fromUtf8("widget"))
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout_2.setMargin(0)
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        self.label_period_min = QtGui.QLabel(self.widget)
        self.label_period_min.setObjectName(_fromUtf8("label_period_min"))
        self.horizontalLayout_2.addWidget(self.label_period_min)
        spacerItem = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.label_period_max = QtGui.QLabel(self.widget)
        self.label_period_max.setObjectName(_fromUtf8("label_period_max"))
        self.horizontalLayout_2.addWidget(self.label_period_max)
        self.verticalLayoutFrequencyPeriod.addWidget(self.widget)
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
        self.label_period_min.setText(_translate("GroupBoxParameters", "Min", None))
        self.label_period_max.setText(_translate("GroupBoxParameters", "Max", None))

