# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_arrow.ui'
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

class Ui_GroupBox_Arrow(object):
    def setupUi(self, GroupBox_Arrow):
        GroupBox_Arrow.setObjectName(_fromUtf8("GroupBox_Arrow"))
        GroupBox_Arrow.setEnabled(True)
        GroupBox_Arrow.resize(300, 314)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Arrow.sizePolicy().hasHeightForWidth())
        GroupBox_Arrow.setSizePolicy(sizePolicy)
        GroupBox_Arrow.setCheckable(False)
        GroupBox_Arrow.setChecked(False)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_Arrow)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.checkBox_real = QtGui.QCheckBox(GroupBox_Arrow)
        self.checkBox_real.setObjectName(_fromUtf8("checkBox_real"))
        self.verticalLayout.addWidget(self.checkBox_real)
        self.checkBox_imaginary = QtGui.QCheckBox(GroupBox_Arrow)
        self.checkBox_imaginary.setObjectName(_fromUtf8("checkBox_imaginary"))
        self.verticalLayout.addWidget(self.checkBox_imaginary)
        self.groupBox_advanced_options = QtGui.QGroupBox(GroupBox_Arrow)
        self.groupBox_advanced_options.setCheckable(True)
        self.groupBox_advanced_options.setChecked(False)
        self.groupBox_advanced_options.setObjectName(_fromUtf8("groupBox_advanced_options"))
        self.formLayout_2 = QtGui.QFormLayout(self.groupBox_advanced_options)
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.label_size = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_size.setObjectName(_fromUtf8("label_size"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_size)
        self.doubleSpinBox_size = QtGui.QDoubleSpinBox(self.groupBox_advanced_options)
        self.doubleSpinBox_size.setDecimals(4)
        self.doubleSpinBox_size.setMaximum(10.0)
        self.doubleSpinBox_size.setSingleStep(0.1)
        self.doubleSpinBox_size.setProperty("value", 5.0)
        self.doubleSpinBox_size.setObjectName(_fromUtf8("doubleSpinBox_size"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_size)
        self.label_head_length = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_head_length.setObjectName(_fromUtf8("label_head_length"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_head_length)
        self.doubleSpinBox_head_length = QtGui.QDoubleSpinBox(self.groupBox_advanced_options)
        self.doubleSpinBox_head_length.setDecimals(4)
        self.doubleSpinBox_head_length.setMaximum(10.0)
        self.doubleSpinBox_head_length.setSingleStep(0.1)
        self.doubleSpinBox_head_length.setProperty("value", 1.5)
        self.doubleSpinBox_head_length.setObjectName(_fromUtf8("doubleSpinBox_head_length"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_head_length)
        self.label_head_width = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_head_width.setObjectName(_fromUtf8("label_head_width"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_head_width)
        self.doubleSpinBox_head_width = QtGui.QDoubleSpinBox(self.groupBox_advanced_options)
        self.doubleSpinBox_head_width.setDecimals(4)
        self.doubleSpinBox_head_width.setMaximum(10.0)
        self.doubleSpinBox_head_width.setSingleStep(0.1)
        self.doubleSpinBox_head_width.setProperty("value", 1.5)
        self.doubleSpinBox_head_width.setObjectName(_fromUtf8("doubleSpinBox_head_width"))
        self.formLayout_2.setWidget(2, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_head_width)
        self.label_line_width = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_line_width.setObjectName(_fromUtf8("label_line_width"))
        self.formLayout_2.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_line_width)
        self.doubleSpinBox_line_width = QtGui.QDoubleSpinBox(self.groupBox_advanced_options)
        self.doubleSpinBox_line_width.setDecimals(4)
        self.doubleSpinBox_line_width.setMaximum(10.0)
        self.doubleSpinBox_line_width.setSingleStep(0.1)
        self.doubleSpinBox_line_width.setProperty("value", 0.5)
        self.doubleSpinBox_line_width.setObjectName(_fromUtf8("doubleSpinBox_line_width"))
        self.formLayout_2.setWidget(3, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_line_width)
        self.label_color_real = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_color_real.setObjectName(_fromUtf8("label_color_real"))
        self.formLayout_2.setWidget(4, QtGui.QFormLayout.LabelRole, self.label_color_real)
        self.comboBox_color_real = QtGui.QComboBox(self.groupBox_advanced_options)
        self.comboBox_color_real.setObjectName(_fromUtf8("comboBox_color_real"))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.comboBox_color_real.addItem(_fromUtf8(""))
        self.formLayout_2.setWidget(4, QtGui.QFormLayout.FieldRole, self.comboBox_color_real)
        self.label_color_imaginary = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_color_imaginary.setObjectName(_fromUtf8("label_color_imaginary"))
        self.formLayout_2.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_color_imaginary)
        self.comboBox_color_imaginary = QtGui.QComboBox(self.groupBox_advanced_options)
        self.comboBox_color_imaginary.setObjectName(_fromUtf8("comboBox_color_imaginary"))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.comboBox_color_imaginary.addItem(_fromUtf8(""))
        self.formLayout_2.setWidget(5, QtGui.QFormLayout.FieldRole, self.comboBox_color_imaginary)
        self.label_threshold = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_threshold.setObjectName(_fromUtf8("label_threshold"))
        self.formLayout_2.setWidget(6, QtGui.QFormLayout.LabelRole, self.label_threshold)
        self.doubleSpinBox_threshold = QtGui.QDoubleSpinBox(self.groupBox_advanced_options)
        self.doubleSpinBox_threshold.setDecimals(4)
        self.doubleSpinBox_threshold.setMaximum(10.0)
        self.doubleSpinBox_threshold.setSingleStep(0.1)
        self.doubleSpinBox_threshold.setProperty("value", 1.0)
        self.doubleSpinBox_threshold.setObjectName(_fromUtf8("doubleSpinBox_threshold"))
        self.formLayout_2.setWidget(6, QtGui.QFormLayout.FieldRole, self.doubleSpinBox_threshold)
        self.label_direction = QtGui.QLabel(self.groupBox_advanced_options)
        self.label_direction.setObjectName(_fromUtf8("label_direction"))
        self.formLayout_2.setWidget(7, QtGui.QFormLayout.LabelRole, self.label_direction)
        self.comboBox_direction = QtGui.QComboBox(self.groupBox_advanced_options)
        self.comboBox_direction.setObjectName(_fromUtf8("comboBox_direction"))
        self.comboBox_direction.addItem(_fromUtf8(""))
        self.comboBox_direction.addItem(_fromUtf8(""))
        self.formLayout_2.setWidget(7, QtGui.QFormLayout.FieldRole, self.comboBox_direction)
        self.doubleSpinBox_size.raise_()
        self.label_size.raise_()
        self.label_head_length.raise_()
        self.doubleSpinBox_head_length.raise_()
        self.label_head_width.raise_()
        self.doubleSpinBox_head_width.raise_()
        self.label_line_width.raise_()
        self.doubleSpinBox_line_width.raise_()
        self.label_color_real.raise_()
        self.comboBox_color_real.raise_()
        self.label_color_imaginary.raise_()
        self.comboBox_color_imaginary.raise_()
        self.label_threshold.raise_()
        self.doubleSpinBox_threshold.raise_()
        self.label_direction.raise_()
        self.comboBox_direction.raise_()
        self.verticalLayout.addWidget(self.groupBox_advanced_options)

        self.retranslateUi(GroupBox_Arrow)
        self.comboBox_color_real.setCurrentIndex(6)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Arrow)

    def retranslateUi(self, GroupBox_Arrow):
        GroupBox_Arrow.setWindowTitle(_translate("GroupBox_Arrow", "GroupBox", None))
        GroupBox_Arrow.setTitle(_translate("GroupBox_Arrow", "Arrow", None))
        self.checkBox_real.setText(_translate("GroupBox_Arrow", "Plot the real induction arrows", None))
        self.checkBox_imaginary.setText(_translate("GroupBox_Arrow", "Plot the imaginary induction arrows", None))
        self.groupBox_advanced_options.setTitle(_translate("GroupBox_Arrow", "Advanced Options", None))
        self.label_size.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>multiplier to scale the arrow</p></body></html>", None))
        self.label_size.setText(_translate("GroupBox_Arrow", "Size", None))
        self.doubleSpinBox_size.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>multiplier to scale the arrow</p></body></html>", None))
        self.label_head_length.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>length of the arrwo head</p></body></html>", None))
        self.label_head_length.setText(_translate("GroupBox_Arrow", "Head Length", None))
        self.doubleSpinBox_head_length.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>length of the arrwo head</p></body></html>", None))
        self.label_head_width.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>width of the arrow head</p></body></html>", None))
        self.label_head_width.setText(_translate("GroupBox_Arrow", "Head Width", None))
        self.doubleSpinBox_head_width.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>width of the arrow head</p></body></html>", None))
        self.label_line_width.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>line width of the arrow</p></body></html>", None))
        self.label_line_width.setText(_translate("GroupBox_Arrow", "Line Width", None))
        self.doubleSpinBox_line_width.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>line width of the arrow</p></body></html>", None))
        self.label_color_real.setText(_translate("GroupBox_Arrow", "Color (real)", None))
        self.comboBox_color_real.setItemText(0, _translate("GroupBox_Arrow", "Blue", None))
        self.comboBox_color_real.setItemText(1, _translate("GroupBox_Arrow", "Green", None))
        self.comboBox_color_real.setItemText(2, _translate("GroupBox_Arrow", "Red", None))
        self.comboBox_color_real.setItemText(3, _translate("GroupBox_Arrow", "Cyan", None))
        self.comboBox_color_real.setItemText(4, _translate("GroupBox_Arrow", "Magenta", None))
        self.comboBox_color_real.setItemText(5, _translate("GroupBox_Arrow", "Yellow", None))
        self.comboBox_color_real.setItemText(6, _translate("GroupBox_Arrow", "Black", None))
        self.comboBox_color_real.setItemText(7, _translate("GroupBox_Arrow", "White", None))
        self.label_color_imaginary.setText(_translate("GroupBox_Arrow", "Color (imaginary)", None))
        self.comboBox_color_imaginary.setItemText(0, _translate("GroupBox_Arrow", "Blue", None))
        self.comboBox_color_imaginary.setItemText(1, _translate("GroupBox_Arrow", "Green", None))
        self.comboBox_color_imaginary.setItemText(2, _translate("GroupBox_Arrow", "Red", None))
        self.comboBox_color_imaginary.setItemText(3, _translate("GroupBox_Arrow", "Cyan", None))
        self.comboBox_color_imaginary.setItemText(4, _translate("GroupBox_Arrow", "Magenta", None))
        self.comboBox_color_imaginary.setItemText(5, _translate("GroupBox_Arrow", "Yellow", None))
        self.comboBox_color_imaginary.setItemText(6, _translate("GroupBox_Arrow", "Black", None))
        self.comboBox_color_imaginary.setItemText(7, _translate("GroupBox_Arrow", "White", None))
        self.label_threshold.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>threshold of which any arrow larger than this number will not be plotted, helps clean up if the data is not good.</p><p><span style=\" font-weight:600;\">NOTE:</span> this is applied before scaling by \'size\'</p></body></html>", None))
        self.label_threshold.setText(_translate("GroupBox_Arrow", "Threshold", None))
        self.doubleSpinBox_threshold.setToolTip(_translate("GroupBox_Arrow", "<html><head/><body><p>threshold of which any arrow larger than this number will not be plotted, helps clean up if the data is not good.</p><p><span style=\" font-weight:600;\">NOTE:</span> this is applied before scaling by \'size\'</p></body></html>", None))
        self.label_direction.setText(_translate("GroupBox_Arrow", "Direction", None))
        self.comboBox_direction.setItemText(0, _translate("GroupBox_Arrow", "Point toward a conductor", None))
        self.comboBox_direction.setItemText(1, _translate("GroupBox_Arrow", "Point away from conductor", None))
