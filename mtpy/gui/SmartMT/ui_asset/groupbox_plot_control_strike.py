# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_plot_control_strike.ui'
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

class Ui_GroupBox_plot_control_strike(object):
    def setupUi(self, GroupBox_plot_control_strike):
        GroupBox_plot_control_strike.setObjectName(_fromUtf8("GroupBox_plot_control_strike"))
        GroupBox_plot_control_strike.resize(316, 289)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_plot_control_strike.sizePolicy().hasHeightForWidth())
        GroupBox_plot_control_strike.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_plot_control_strike)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox_2 = QtGui.QGroupBox(GroupBox_plot_control_strike)
        self.groupBox_2.setFlat(True)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.checkBox_plot_tipper = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_plot_tipper.setObjectName(_fromUtf8("checkBox_plot_tipper"))
        self.gridLayout_2.addWidget(self.checkBox_plot_tipper, 2, 0, 1, 1)
        self.checkBox_fold = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_fold.setObjectName(_fromUtf8("checkBox_fold"))
        self.gridLayout_2.addWidget(self.checkBox_fold, 2, 1, 1, 1)
        self.radioButton_type_2 = QtGui.QRadioButton(self.groupBox_2)
        self.radioButton_type_2.setChecked(True)
        self.radioButton_type_2.setObjectName(_fromUtf8("radioButton_type_2"))
        self.gridLayout_2.addWidget(self.radioButton_type_2, 1, 0, 1, 2)
        self.radioButton_type_1 = QtGui.QRadioButton(self.groupBox_2)
        self.radioButton_type_1.setObjectName(_fromUtf8("radioButton_type_1"))
        self.gridLayout_2.addWidget(self.radioButton_type_1, 0, 0, 1, 2)
        self.checkBox_max_error = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_max_error.setObjectName(_fromUtf8("checkBox_max_error"))
        self.gridLayout_2.addWidget(self.checkBox_max_error, 3, 0, 1, 1)
        self.doubleSpinBox_max_error = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.doubleSpinBox_max_error.setEnabled(False)
        self.doubleSpinBox_max_error.setMaximum(360.0)
        self.doubleSpinBox_max_error.setObjectName(_fromUtf8("doubleSpinBox_max_error"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_max_error, 3, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.groupBox = QtGui.QGroupBox(GroupBox_plot_control_strike)
        self.groupBox.setFlat(True)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_min = QtGui.QLabel(self.groupBox)
        self.label_min.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_min.sizePolicy().hasHeightForWidth())
        self.label_min.setSizePolicy(sizePolicy)
        self.label_min.setObjectName(_fromUtf8("label_min"))
        self.gridLayout.addWidget(self.label_min, 1, 1, 1, 1)
        self.radioButton_range_data = QtGui.QRadioButton(self.groupBox)
        self.radioButton_range_data.setChecked(True)
        self.radioButton_range_data.setObjectName(_fromUtf8("radioButton_range_data"))
        self.gridLayout.addWidget(self.radioButton_range_data, 0, 0, 1, 1)
        self.radioButton_range_minmax = QtGui.QRadioButton(self.groupBox)
        self.radioButton_range_minmax.setObjectName(_fromUtf8("radioButton_range_minmax"))
        self.gridLayout.addWidget(self.radioButton_range_minmax, 1, 0, 1, 1)
        self.label_max = QtGui.QLabel(self.groupBox)
        self.label_max.setEnabled(False)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_max.sizePolicy().hasHeightForWidth())
        self.label_max.setSizePolicy(sizePolicy)
        self.label_max.setObjectName(_fromUtf8("label_max"))
        self.gridLayout.addWidget(self.label_max, 2, 1, 1, 1)
        self.lineEdit_min = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_min.setEnabled(False)
        self.lineEdit_min.setObjectName(_fromUtf8("lineEdit_min"))
        self.gridLayout.addWidget(self.lineEdit_min, 1, 2, 1, 1)
        self.lineEdit_max = QtGui.QLineEdit(self.groupBox)
        self.lineEdit_max.setEnabled(False)
        self.lineEdit_max.setObjectName(_fromUtf8("lineEdit_max"))
        self.gridLayout.addWidget(self.lineEdit_max, 2, 2, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)

        self.retranslateUi(GroupBox_plot_control_strike)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_plot_control_strike)
        GroupBox_plot_control_strike.setTabOrder(self.radioButton_type_1, self.radioButton_type_2)
        GroupBox_plot_control_strike.setTabOrder(self.radioButton_type_2, self.checkBox_plot_tipper)
        GroupBox_plot_control_strike.setTabOrder(self.checkBox_plot_tipper, self.checkBox_fold)
        GroupBox_plot_control_strike.setTabOrder(self.checkBox_fold, self.checkBox_max_error)
        GroupBox_plot_control_strike.setTabOrder(self.checkBox_max_error, self.doubleSpinBox_max_error)
        GroupBox_plot_control_strike.setTabOrder(self.doubleSpinBox_max_error, self.radioButton_range_data)
        GroupBox_plot_control_strike.setTabOrder(self.radioButton_range_data, self.radioButton_range_minmax)
        GroupBox_plot_control_strike.setTabOrder(self.radioButton_range_minmax, self.lineEdit_min)
        GroupBox_plot_control_strike.setTabOrder(self.lineEdit_min, self.lineEdit_max)

    def retranslateUi(self, GroupBox_plot_control_strike):
        GroupBox_plot_control_strike.setWindowTitle(_translate("GroupBox_plot_control_strike", "GroupBox", None))
        GroupBox_plot_control_strike.setTitle(_translate("GroupBox_plot_control_strike", "Plot Settings", None))
        self.groupBox_2.setTitle(_translate("GroupBox_plot_control_strike", "Plot Type", None))
        self.checkBox_plot_tipper.setText(_translate("GroupBox_plot_control_strike", "Plot Tipper Strike", None))
        self.checkBox_fold.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>check to plot only from 0 to\n"
"                                        180 degree, 0 to 360 degree if otherwise</p></body></html>\n"
"                                    ", None))
        self.checkBox_fold.setText(_translate("GroupBox_plot_control_strike", "Fold", None))
        self.radioButton_type_2.setText(_translate("GroupBox_plot_control_strike", "plot all period ranges into one polar diagram for\n"
"                                        each strike angle estimation\n"
"                                    ", None))
        self.radioButton_type_1.setText(_translate("GroupBox_plot_control_strike", "plot individual decades in one plot", None))
        self.checkBox_max_error.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>Maximum error in degrees that\n"
"                                        is allowed to estimate strike. if not checked/provided, all estimates are\n"
"                                        allowed.</p></body></html>\n"
"                                    ", None))
        self.checkBox_max_error.setText(_translate("GroupBox_plot_control_strike", "Maximum Error", None))
        self.doubleSpinBox_max_error.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>Maximum error in degrees that\n"
"                                        is allowed to estimate strike. if not checked/provided, all estimates are\n"
"                                        allowed.</p></body></html>\n"
"                                    ", None))
        self.doubleSpinBox_max_error.setSuffix(_translate("GroupBox_plot_control_strike", "°", None))
        self.groupBox.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>Period range to estimate the strike angle</p></body></html>", None))
        self.groupBox.setTitle(_translate("GroupBox_plot_control_strike", "Plot Range", None))
        self.label_min.setText(_translate("GroupBox_plot_control_strike", "min=", None))
        self.radioButton_range_data.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>estimating the strike for all\n"
"                                        periods in the data</p></body></html>\n"
"                                    ", None))
        self.radioButton_range_data.setText(_translate("GroupBox_plot_control_strike", "Based on Data", None))
        self.radioButton_range_minmax.setToolTip(_translate("GroupBox_plot_control_strike", "<html><head/><body><p>define minimum and maximum\n"
"                                        periods, will be used as log10(min) and log10(max)</p></body></html>\n"
"                                    ", None))
        self.radioButton_range_minmax.setText(_translate("GroupBox_plot_control_strike", "(min, max)", None))
        self.label_max.setText(_translate("GroupBox_plot_control_strike", "max=", None))
        self.lineEdit_min.setText(_translate("GroupBox_plot_control_strike", "0", None))
        self.lineEdit_max.setText(_translate("GroupBox_plot_control_strike", "inf", None))

