# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_plot_control_resistivity_phase_pseudo_section.ui'
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

class Ui_GroupBox_plot_control_resistivity_phase_pseudo_section(object):
    def setupUi(self, GroupBox_plot_control_resistivity_phase_pseudo_section):
        GroupBox_plot_control_resistivity_phase_pseudo_section.setObjectName(_fromUtf8("GroupBox_plot_control_resistivity_phase_pseudo_section"))
        GroupBox_plot_control_resistivity_phase_pseudo_section.resize(300, 300)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_plot_control_resistivity_phase_pseudo_section.sizePolicy().hasHeightForWidth())
        GroupBox_plot_control_resistivity_phase_pseudo_section.setSizePolicy(sizePolicy)
        self.gridLayout = QtGui.QGridLayout(GroupBox_plot_control_resistivity_phase_pseudo_section)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.groupBox = QtGui.QGroupBox(GroupBox_plot_control_resistivity_phase_pseudo_section)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.radioButton_aspect_equal = QtGui.QRadioButton(self.groupBox)
        self.radioButton_aspect_equal.setObjectName(_fromUtf8("radioButton_aspect_equal"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_equal, 0, 1, 1, 1)
        self.radioButton_aspect_float = QtGui.QRadioButton(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.radioButton_aspect_float.sizePolicy().hasHeightForWidth())
        self.radioButton_aspect_float.setSizePolicy(sizePolicy)
        self.radioButton_aspect_float.setText(_fromUtf8(""))
        self.radioButton_aspect_float.setObjectName(_fromUtf8("radioButton_aspect_float"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_float, 0, 2, 1, 1)
        self.radioButton_aspect_auto = QtGui.QRadioButton(self.groupBox)
        self.radioButton_aspect_auto.setChecked(True)
        self.radioButton_aspect_auto.setObjectName(_fromUtf8("radioButton_aspect_auto"))
        self.gridLayout_2.addWidget(self.radioButton_aspect_auto, 0, 0, 1, 1)
        self.doubleSpinBox_aspect_float = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_aspect_float.setEnabled(False)
        self.doubleSpinBox_aspect_float.setProperty("value", 1.0)
        self.doubleSpinBox_aspect_float.setObjectName(_fromUtf8("doubleSpinBox_aspect_float"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_aspect_float, 0, 3, 1, 1)
        self.gridLayout.addWidget(self.groupBox, 0, 1, 1, 1)

        self.retranslateUi(GroupBox_plot_control_resistivity_phase_pseudo_section)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_plot_control_resistivity_phase_pseudo_section)

    def retranslateUi(self, GroupBox_plot_control_resistivity_phase_pseudo_section):
        GroupBox_plot_control_resistivity_phase_pseudo_section.setWindowTitle(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "GroupBox", None))
        GroupBox_plot_control_resistivity_phase_pseudo_section.setTitle(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "Plot Control", None))
        self.groupBox.setToolTip(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "<html><head/><body><p>Aspect ratio of each subplot (height/width)</p></body></html>", None))
        self.groupBox.setTitle(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "Aspect Parameter", None))
        self.radioButton_aspect_equal.setText(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "Equal", None))
        self.radioButton_aspect_auto.setText(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "Auto", None))
        self.doubleSpinBox_aspect_float.setToolTip(_translate("GroupBox_plot_control_resistivity_phase_pseudo_section", "<html><head/><body><p>height/width</p></body></html>", None))

