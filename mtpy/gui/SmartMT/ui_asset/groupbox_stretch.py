# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_stretch.ui'
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

class Ui_GroupBox_Stretch(object):
    def setupUi(self, GroupBox_Stretch):
        GroupBox_Stretch.setObjectName(_fromUtf8("GroupBox_Stretch"))
        GroupBox_Stretch.resize(300, 183)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Stretch.sizePolicy().hasHeightForWidth())
        GroupBox_Stretch.setSizePolicy(sizePolicy)
        self.formLayout = QtGui.QFormLayout(GroupBox_Stretch)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(GroupBox_Stretch)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.doubleSpinBox_x = QtGui.QDoubleSpinBox(GroupBox_Stretch)
        self.doubleSpinBox_x.setMaximum(1000000.0)
        self.doubleSpinBox_x.setProperty("value", 200.0)
        self.doubleSpinBox_x.setObjectName(_fromUtf8("doubleSpinBox_x"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_x)
        self.label_2 = QtGui.QLabel(GroupBox_Stretch)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.doubleSpinBox_y = QtGui.QDoubleSpinBox(GroupBox_Stretch)
        self.doubleSpinBox_y.setMaximum(1000000.0)
        self.doubleSpinBox_y.setProperty("value", 5.0)
        self.doubleSpinBox_y.setObjectName(_fromUtf8("doubleSpinBox_y"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_y)
        self.groupBox = QtGui.QGroupBox(GroupBox_Stretch)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.doubleSpinBox_x_min = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_x_min.setEnabled(False)
        self.doubleSpinBox_x_min.setMaximum(1000000.0)
        self.doubleSpinBox_x_min.setObjectName(_fromUtf8("doubleSpinBox_x_min"))
        self.gridLayout.addWidget(self.doubleSpinBox_x_min, 1, 1, 1, 1)
        self.checkBox_x_range = QtGui.QCheckBox(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_x_range.sizePolicy().hasHeightForWidth())
        self.checkBox_x_range.setSizePolicy(sizePolicy)
        self.checkBox_x_range.setObjectName(_fromUtf8("checkBox_x_range"))
        self.gridLayout.addWidget(self.checkBox_x_range, 1, 0, 1, 1)
        self.doubleSpinBox_y_min = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_y_min.setEnabled(False)
        self.doubleSpinBox_y_min.setMaximum(1000000.0)
        self.doubleSpinBox_y_min.setObjectName(_fromUtf8("doubleSpinBox_y_min"))
        self.gridLayout.addWidget(self.doubleSpinBox_y_min, 3, 1, 1, 1)
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout.addWidget(self.label_6, 0, 2, 1, 1)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout.addWidget(self.label_5, 0, 1, 1, 1)
        self.doubleSpinBox_x_max = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_x_max.setEnabled(False)
        self.doubleSpinBox_x_max.setMaximum(1000000.0)
        self.doubleSpinBox_x_max.setProperty("value", 100.0)
        self.doubleSpinBox_x_max.setObjectName(_fromUtf8("doubleSpinBox_x_max"))
        self.gridLayout.addWidget(self.doubleSpinBox_x_max, 1, 2, 1, 1)
        self.doubleSpinBox_y_max = QtGui.QDoubleSpinBox(self.groupBox)
        self.doubleSpinBox_y_max.setEnabled(False)
        self.doubleSpinBox_y_max.setMaximum(1000000.0)
        self.doubleSpinBox_y_max.setProperty("value", 100.0)
        self.doubleSpinBox_y_max.setObjectName(_fromUtf8("doubleSpinBox_y_max"))
        self.gridLayout.addWidget(self.doubleSpinBox_y_max, 3, 2, 1, 1)
        self.checkBox_y_range = QtGui.QCheckBox(self.groupBox)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_y_range.sizePolicy().hasHeightForWidth())
        self.checkBox_y_range.setSizePolicy(sizePolicy)
        self.checkBox_y_range.setObjectName(_fromUtf8("checkBox_y_range"))
        self.gridLayout.addWidget(self.checkBox_y_range, 3, 0, 1, 1)
        self.formLayout.setWidget(2, QtGui.QFormLayout.SpanningRole, self.groupBox)

        self.retranslateUi(GroupBox_Stretch)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Stretch)
        GroupBox_Stretch.setTabOrder(self.doubleSpinBox_x, self.doubleSpinBox_y)
        GroupBox_Stretch.setTabOrder(self.doubleSpinBox_y, self.checkBox_x_range)
        GroupBox_Stretch.setTabOrder(self.checkBox_x_range, self.doubleSpinBox_x_min)
        GroupBox_Stretch.setTabOrder(self.doubleSpinBox_x_min, self.doubleSpinBox_x_max)
        GroupBox_Stretch.setTabOrder(self.doubleSpinBox_x_max, self.checkBox_y_range)
        GroupBox_Stretch.setTabOrder(self.checkBox_y_range, self.doubleSpinBox_y_min)
        GroupBox_Stretch.setTabOrder(self.doubleSpinBox_y_min, self.doubleSpinBox_y_max)

    def retranslateUi(self, GroupBox_Stretch):
        GroupBox_Stretch.setWindowTitle(_translate("GroupBox_Stretch", "GroupBox", None))
        GroupBox_Stretch.setToolTip(_translate("GroupBox_Stretch", "Stretch is a set of factor that scales the distance from one station to the next to make the plot readable", None))
        GroupBox_Stretch.setTitle(_translate("GroupBox_Stretch", "Stretch", None))
        self.label.setText(_translate("GroupBox_Stretch", "X Stretch", None))
        self.label_2.setText(_translate("GroupBox_Stretch", "Y Stretch", None))
        self.groupBox.setTitle(_translate("GroupBox_Stretch", "Range", None))
        self.doubleSpinBox_x_min.setToolTip(_translate("GroupBox_Stretch", "<html><head/><body><p>Minimum and maximum along the x-axis in relative distance of degrees and multiplied by X Stretch.</p></body></html>", None))
        self.checkBox_x_range.setText(_translate("GroupBox_Stretch", "X Range", None))
        self.doubleSpinBox_y_min.setToolTip(_translate("GroupBox_Stretch", "<html><head/><body><p>Minimum and maximum period to plot, note that the scaling will be done in the code. So if you want to plot from (.1s, 100s) </p></body></html>", None))
        self.doubleSpinBox_y_min.setSuffix(_translate("GroupBox_Stretch", "s", None))
        self.label_6.setText(_translate("GroupBox_Stretch", "Maximum", None))
        self.label_5.setText(_translate("GroupBox_Stretch", "Mimimum", None))
        self.doubleSpinBox_x_max.setToolTip(_translate("GroupBox_Stretch", "<html><head/><body><p>Minimum and maximum along the x-axis in relative distance of degrees and multiplied by X Stretch.</p></body></html>", None))
        self.doubleSpinBox_y_max.setToolTip(_translate("GroupBox_Stretch", "<html><head/><body><p>Minimum and maximum period to plot, note that the scaling will be done in the code. So if you want to plot from (.1s, 100s) </p></body></html>", None))
        self.doubleSpinBox_y_max.setSuffix(_translate("GroupBox_Stretch", "s", None))
        self.checkBox_y_range.setToolTip(_translate("GroupBox_Stretch", "<html><head/><body><p>Minimum and maximum period to plot, note that the scaling will be done in the code. So if you want to plot from (.1s, 100s) </p></body></html>", None))
        self.checkBox_y_range.setText(_translate("GroupBox_Stretch", "Y Range", None))

