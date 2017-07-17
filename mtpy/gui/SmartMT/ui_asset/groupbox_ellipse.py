# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_ellipse.ui'
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

class Ui_GroupBoxEllipse(object):
    def setupUi(self, GroupBoxEllipse):
        GroupBoxEllipse.setObjectName(_fromUtf8("GroupBoxEllipse"))
        GroupBoxEllipse.resize(300, 208)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBoxEllipse.sizePolicy().hasHeightForWidth())
        GroupBoxEllipse.setSizePolicy(sizePolicy)
        self.formLayout = QtGui.QFormLayout(GroupBoxEllipse)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.doubleSpinBox_size = QtGui.QDoubleSpinBox(GroupBoxEllipse)
        self.doubleSpinBox_size.setMaximum(10000000.0)
        self.doubleSpinBox_size.setProperty("value", 2.0)
        self.doubleSpinBox_size.setObjectName(_fromUtf8("doubleSpinBox_size"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_size)
        self.label_2 = QtGui.QLabel(GroupBoxEllipse)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.comboBoxColor_by = QtGui.QComboBox(GroupBoxEllipse)
        self.comboBoxColor_by.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToMinimumContentsLength)
        self.comboBoxColor_by.setObjectName(_fromUtf8("comboBoxColor_by"))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.comboBoxColor_by.addItem(_fromUtf8(""))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.comboBoxColor_by)
        self.label_6 = QtGui.QLabel(GroupBoxEllipse)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_6)
        self.comboBox_cmap = QtGui.QComboBox(GroupBoxEllipse)
        self.comboBox_cmap.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToMinimumContentsLength)
        self.comboBox_cmap.setObjectName(_fromUtf8("comboBox_cmap"))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.comboBox_cmap.addItem(_fromUtf8(""))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.comboBox_cmap)
        self.label_7 = QtGui.QLabel(GroupBoxEllipse)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_7)
        self.label_3 = QtGui.QLabel(GroupBoxEllipse)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_3)
        self.doubleSpinBox_min = QtGui.QDoubleSpinBox(GroupBoxEllipse)
        self.doubleSpinBox_min.setMaximum(10000000.0)
        self.doubleSpinBox_min.setObjectName(_fromUtf8("doubleSpinBox_min"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_min)
        self.label_4 = QtGui.QLabel(GroupBoxEllipse)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.LabelRole, self.label_4)
        self.doubleSpinBox_max = QtGui.QDoubleSpinBox(GroupBoxEllipse)
        self.doubleSpinBox_max.setMaximum(10000000.0)
        self.doubleSpinBox_max.setProperty("value", 90.0)
        self.doubleSpinBox_max.setObjectName(_fromUtf8("doubleSpinBox_max"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_max)
        self.label_5 = QtGui.QLabel(GroupBoxEllipse)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout.setWidget(7, QtGui.QFormLayout.LabelRole, self.label_5)
        self.doubleSpinBox_step = QtGui.QDoubleSpinBox(GroupBoxEllipse)
        self.doubleSpinBox_step.setMaximum(10000000.0)
        self.doubleSpinBox_step.setProperty("value", 1.0)
        self.doubleSpinBox_step.setObjectName(_fromUtf8("doubleSpinBox_step"))
        self.formLayout.setWidget(7, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_step)
        self.label = QtGui.QLabel(GroupBoxEllipse)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.line = QtGui.QFrame(GroupBoxEllipse)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.line)

        self.retranslateUi(GroupBoxEllipse)
        self.comboBox_cmap.setCurrentIndex(5)
        QtCore.QMetaObject.connectSlotsByName(GroupBoxEllipse)

    def retranslateUi(self, GroupBoxEllipse):
        GroupBoxEllipse.setWindowTitle(_translate("GroupBoxEllipse", "GroupBox", None))
        GroupBoxEllipse.setToolTip(_translate("GroupBoxEllipse", "<html><head/><body><p>Parameters for the phase tensor ellipses</p></body></html>", None))
        GroupBoxEllipse.setTitle(_translate("GroupBoxEllipse", "Ellipse", None))
        self.doubleSpinBox_size.setToolTip(_translate("GroupBoxEllipse", "<html><head/><body><p>Size of ellipse in points</p><p>Default: 2</p></body></html>", None))
        self.label_2.setToolTip(_translate("GroupBoxEllipse", "<html><head/><body><p>Choose how the phase ellipse is colored</p></body></html>", None))
        self.label_2.setText(_translate("GroupBoxEllipse", "Color by", None))
        self.comboBoxColor_by.setToolTip(_translate("GroupBoxEllipse", "<html><head/><body><p>Choose how the phase ellipse is colored</p></body></html>", None))
        self.comboBoxColor_by.setItemText(0, _translate("GroupBoxEllipse", "colors by minimum phase (phimin)", None))
        self.comboBoxColor_by.setItemText(1, _translate("GroupBoxEllipse", "colors by maximum phase (phimax)", None))
        self.comboBoxColor_by.setItemText(2, _translate("GroupBoxEllipse", "colors by skew (skew)", None))
        self.comboBoxColor_by.setItemText(3, _translate("GroupBoxEllipse", "colors by skew in discrete segments defined by the range (skew_seg)", None))
        self.comboBoxColor_by.setItemText(4, _translate("GroupBoxEllipse", "colors by skew see [Booker, 2014] (normalized_skew)", None))
        self.comboBoxColor_by.setItemText(5, _translate("GroupBoxEllipse", "colors by normalized skew in discrete segments defined by the range (normalized_skew_seg)", None))
        self.comboBoxColor_by.setItemText(6, _translate("GroupBoxEllipse", "colors by determinant of the phase tensor (phidet)", None))
        self.comboBoxColor_by.setItemText(7, _translate("GroupBoxEllipse", "colors by ellipticity (ellipticity)", None))
        self.label_6.setToolTip(_translate("GroupBoxEllipse", "color map of tensor ellipses", None))
        self.label_6.setText(_translate("GroupBoxEllipse", "Color Map", None))
        self.comboBox_cmap.setToolTip(_translate("GroupBoxEllipse", "color map of tensor ellipses", None))
        self.comboBox_cmap.setItemText(0, _translate("GroupBoxEllipse", "yellow to red (mt_yl2rd) ", None))
        self.comboBox_cmap.setItemText(1, _translate("GroupBoxEllipse", "blue to yellow to red (mt_bl2yl2rd)", None))
        self.comboBox_cmap.setItemText(2, _translate("GroupBoxEllipse", "white to blue (mt_wh2bl)", None))
        self.comboBox_cmap.setItemText(3, _translate("GroupBoxEllipse", "red to blue (mt_rd2bl) ", None))
        self.comboBox_cmap.setItemText(4, _translate("GroupBoxEllipse", "blue to white to red (mt_bl2wh2rd)", None))
        self.comboBox_cmap.setItemText(5, _translate("GroupBoxEllipse", "blue to green to red (mt_bl2gr2rd)", None))
        self.comboBox_cmap.setItemText(6, _translate("GroupBoxEllipse", "red to green to blue (mt_rd2gr2bl)", None))
        self.comboBox_cmap.setItemText(7, _translate("GroupBoxEllipse", "discrete blue to white to red (mt_seg_bl2wh2rd)", None))
        self.label_7.setText(_translate("GroupBoxEllipse", "Range:", None))
        self.label_3.setText(_translate("GroupBoxEllipse", "Minimum", None))
        self.label_4.setText(_translate("GroupBoxEllipse", "Maximum", None))
        self.label_5.setText(_translate("GroupBoxEllipse", "Discrete Step", None))
        self.label.setToolTip(_translate("GroupBoxEllipse", "<html><head/><body><p>Size of ellipse in points</p><p>Default: 2</p></body></html>", None))
        self.label.setText(_translate("GroupBoxEllipse", "Size", None))

