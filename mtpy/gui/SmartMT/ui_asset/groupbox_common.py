# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_common.ui'
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

class Ui_GroupBox_common_settings(object):
    def setupUi(self, GroupBox_common_settings):
        GroupBox_common_settings.setObjectName(_fromUtf8("GroupBox_common_settings"))
        GroupBox_common_settings.resize(300, 327)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_common_settings.sizePolicy().hasHeightForWidth())
        GroupBox_common_settings.setSizePolicy(sizePolicy)
        self.gridLayout = QtGui.QGridLayout(GroupBox_common_settings)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.groupBox_figure_size = QtGui.QGroupBox(GroupBox_common_settings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_figure_size.sizePolicy().hasHeightForWidth())
        self.groupBox_figure_size.setSizePolicy(sizePolicy)
        self.groupBox_figure_size.setCheckable(True)
        self.groupBox_figure_size.setChecked(False)
        self.groupBox_figure_size.setObjectName(_fromUtf8("groupBox_figure_size"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_figure_size)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_4 = QtGui.QLabel(self.groupBox_figure_size)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_3.addWidget(self.label_4, 2, 0, 1, 1)
        self.label_3 = QtGui.QLabel(self.groupBox_figure_size)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout_3.addWidget(self.label_3, 1, 0, 1, 1)
        self.doubleSpinBox_height_inches = QtGui.QDoubleSpinBox(self.groupBox_figure_size)
        self.doubleSpinBox_height_inches.setMinimum(1.0)
        self.doubleSpinBox_height_inches.setMaximum(1000.0)
        self.doubleSpinBox_height_inches.setProperty("value", 6.0)
        self.doubleSpinBox_height_inches.setObjectName(_fromUtf8("doubleSpinBox_height_inches"))
        self.gridLayout_3.addWidget(self.doubleSpinBox_height_inches, 1, 1, 1, 1)
        self.spinBox_height_pixels = QtGui.QSpinBox(self.groupBox_figure_size)
        self.spinBox_height_pixels.setMinimum(1)
        self.spinBox_height_pixels.setMaximum(1000000)
        self.spinBox_height_pixels.setProperty("value", 600)
        self.spinBox_height_pixels.setObjectName(_fromUtf8("spinBox_height_pixels"))
        self.gridLayout_3.addWidget(self.spinBox_height_pixels, 1, 2, 1, 1)
        self.doubleSpinBox_width_inches = QtGui.QDoubleSpinBox(self.groupBox_figure_size)
        self.doubleSpinBox_width_inches.setMinimum(1.0)
        self.doubleSpinBox_width_inches.setMaximum(1000.0)
        self.doubleSpinBox_width_inches.setProperty("value", 8.0)
        self.doubleSpinBox_width_inches.setObjectName(_fromUtf8("doubleSpinBox_width_inches"))
        self.gridLayout_3.addWidget(self.doubleSpinBox_width_inches, 0, 1, 1, 1)
        self.spinBox_dpi = QtGui.QSpinBox(self.groupBox_figure_size)
        self.spinBox_dpi.setMinimum(1)
        self.spinBox_dpi.setMaximum(1000)
        self.spinBox_dpi.setProperty("value", 80)
        self.spinBox_dpi.setObjectName(_fromUtf8("spinBox_dpi"))
        self.gridLayout_3.addWidget(self.spinBox_dpi, 2, 1, 1, 1)
        self.spinBox_width_pixels = QtGui.QSpinBox(self.groupBox_figure_size)
        self.spinBox_width_pixels.setMinimum(1)
        self.spinBox_width_pixels.setMaximum(1000000)
        self.spinBox_width_pixels.setProperty("value", 800)
        self.spinBox_width_pixels.setObjectName(_fromUtf8("spinBox_width_pixels"))
        self.gridLayout_3.addWidget(self.spinBox_width_pixels, 0, 2, 1, 1)
        self.checkBox_tight_layout = QtGui.QCheckBox(self.groupBox_figure_size)
        self.checkBox_tight_layout.setChecked(True)
        self.checkBox_tight_layout.setObjectName(_fromUtf8("checkBox_tight_layout"))
        self.gridLayout_3.addWidget(self.checkBox_tight_layout, 2, 2, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox_figure_size)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout_3.addWidget(self.label_2, 0, 0, 1, 1)
        self.gridLayout.addWidget(self.groupBox_figure_size, 2, 0, 1, 3)
        self.groupBox_title = QtGui.QGroupBox(GroupBox_common_settings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_title.sizePolicy().hasHeightForWidth())
        self.groupBox_title.setSizePolicy(sizePolicy)
        self.groupBox_title.setCheckable(True)
        self.groupBox_title.setChecked(False)
        self.groupBox_title.setObjectName(_fromUtf8("groupBox_title"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_title)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.horizontalSlider_y = QtGui.QSlider(self.groupBox_title)
        self.horizontalSlider_y.setMaximum(100)
        self.horizontalSlider_y.setProperty("value", 98)
        self.horizontalSlider_y.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_y.setObjectName(_fromUtf8("horizontalSlider_y"))
        self.gridLayout_2.addWidget(self.horizontalSlider_y, 4, 1, 1, 2)
        self.lineEdit_title = QtGui.QLineEdit(self.groupBox_title)
        self.lineEdit_title.setObjectName(_fromUtf8("lineEdit_title"))
        self.gridLayout_2.addWidget(self.lineEdit_title, 0, 1, 1, 3)
        self.doubleSpinBox_x = QtGui.QDoubleSpinBox(self.groupBox_title)
        self.doubleSpinBox_x.setMaximum(1.0)
        self.doubleSpinBox_x.setSingleStep(0.1)
        self.doubleSpinBox_x.setProperty("value", 0.5)
        self.doubleSpinBox_x.setObjectName(_fromUtf8("doubleSpinBox_x"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_x, 3, 3, 1, 1)
        self.label_5 = QtGui.QLabel(self.groupBox_title)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_2.addWidget(self.label_5, 3, 0, 1, 1)
        self.doubleSpinBox_y = QtGui.QDoubleSpinBox(self.groupBox_title)
        self.doubleSpinBox_y.setMaximum(1.0)
        self.doubleSpinBox_y.setSingleStep(0.1)
        self.doubleSpinBox_y.setProperty("value", 0.98)
        self.doubleSpinBox_y.setObjectName(_fromUtf8("doubleSpinBox_y"))
        self.gridLayout_2.addWidget(self.doubleSpinBox_y, 4, 3, 1, 1)
        self.label = QtGui.QLabel(self.groupBox_title)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.label_6 = QtGui.QLabel(self.groupBox_title)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.gridLayout_2.addWidget(self.label_6, 4, 0, 1, 1)
        self.horizontalSlider_x = QtGui.QSlider(self.groupBox_title)
        self.horizontalSlider_x.setMaximum(100)
        self.horizontalSlider_x.setProperty("value", 50)
        self.horizontalSlider_x.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_x.setObjectName(_fromUtf8("horizontalSlider_x"))
        self.gridLayout_2.addWidget(self.horizontalSlider_x, 3, 1, 1, 2)
        self.label_9 = QtGui.QLabel(self.groupBox_title)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_9.sizePolicy().hasHeightForWidth())
        self.label_9.setSizePolicy(sizePolicy)
        self.label_9.setObjectName(_fromUtf8("label_9"))
        self.gridLayout_2.addWidget(self.label_9, 1, 0, 1, 1)
        self.label_8 = QtGui.QLabel(self.groupBox_title)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.gridLayout_2.addWidget(self.label_8, 6, 0, 1, 2)
        self.comboBox_horizontal_alignment = QtGui.QComboBox(self.groupBox_title)
        self.comboBox_horizontal_alignment.setObjectName(_fromUtf8("comboBox_horizontal_alignment"))
        self.comboBox_horizontal_alignment.addItem(_fromUtf8(""))
        self.comboBox_horizontal_alignment.addItem(_fromUtf8(""))
        self.comboBox_horizontal_alignment.addItem(_fromUtf8(""))
        self.gridLayout_2.addWidget(self.comboBox_horizontal_alignment, 5, 2, 1, 2)
        self.label_7 = QtGui.QLabel(self.groupBox_title)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.gridLayout_2.addWidget(self.label_7, 5, 0, 1, 2)
        self.comboBox_vertical_alignment = QtGui.QComboBox(self.groupBox_title)
        self.comboBox_vertical_alignment.setObjectName(_fromUtf8("comboBox_vertical_alignment"))
        self.comboBox_vertical_alignment.addItem(_fromUtf8(""))
        self.comboBox_vertical_alignment.addItem(_fromUtf8(""))
        self.comboBox_vertical_alignment.addItem(_fromUtf8(""))
        self.comboBox_vertical_alignment.addItem(_fromUtf8(""))
        self.gridLayout_2.addWidget(self.comboBox_vertical_alignment, 6, 2, 1, 2)
        self.spinBox_fontsize = QtGui.QSpinBox(self.groupBox_title)
        self.spinBox_fontsize.setMinimum(4)
        self.spinBox_fontsize.setProperty("value", 12)
        self.spinBox_fontsize.setObjectName(_fromUtf8("spinBox_fontsize"))
        self.gridLayout_2.addWidget(self.spinBox_fontsize, 1, 1, 1, 1)
        self.gridLayout.addWidget(self.groupBox_title, 1, 0, 1, 3)

        self.retranslateUi(GroupBox_common_settings)
        self.comboBox_horizontal_alignment.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_common_settings)
        GroupBox_common_settings.setTabOrder(self.groupBox_title, self.lineEdit_title)
        GroupBox_common_settings.setTabOrder(self.lineEdit_title, self.spinBox_fontsize)
        GroupBox_common_settings.setTabOrder(self.spinBox_fontsize, self.horizontalSlider_x)
        GroupBox_common_settings.setTabOrder(self.horizontalSlider_x, self.doubleSpinBox_x)
        GroupBox_common_settings.setTabOrder(self.doubleSpinBox_x, self.horizontalSlider_y)
        GroupBox_common_settings.setTabOrder(self.horizontalSlider_y, self.doubleSpinBox_y)
        GroupBox_common_settings.setTabOrder(self.doubleSpinBox_y, self.comboBox_horizontal_alignment)
        GroupBox_common_settings.setTabOrder(self.comboBox_horizontal_alignment, self.comboBox_vertical_alignment)
        GroupBox_common_settings.setTabOrder(self.comboBox_vertical_alignment, self.groupBox_figure_size)
        GroupBox_common_settings.setTabOrder(self.groupBox_figure_size, self.doubleSpinBox_width_inches)
        GroupBox_common_settings.setTabOrder(self.doubleSpinBox_width_inches, self.spinBox_width_pixels)
        GroupBox_common_settings.setTabOrder(self.spinBox_width_pixels, self.doubleSpinBox_height_inches)
        GroupBox_common_settings.setTabOrder(self.doubleSpinBox_height_inches, self.spinBox_height_pixels)
        GroupBox_common_settings.setTabOrder(self.spinBox_height_pixels, self.spinBox_dpi)
        GroupBox_common_settings.setTabOrder(self.spinBox_dpi, self.checkBox_tight_layout)

    def retranslateUi(self, GroupBox_common_settings):
        GroupBox_common_settings.setWindowTitle(_translate("GroupBox_common_settings", "GroupBox", None))
        GroupBox_common_settings.setTitle(_translate("GroupBox_common_settings", "General Settings", None))
        self.groupBox_figure_size.setTitle(_translate("GroupBox_common_settings", "Figure Size", None))
        self.label_4.setText(_translate("GroupBox_common_settings", "DPI", None))
        self.label_3.setText(_translate("GroupBox_common_settings", "Height", None))
        self.doubleSpinBox_height_inches.setSuffix(_translate("GroupBox_common_settings", " inches", None))
        self.spinBox_height_pixels.setSuffix(_translate("GroupBox_common_settings", " pixels", None))
        self.doubleSpinBox_width_inches.setSuffix(_translate("GroupBox_common_settings", " inches", None))
        self.spinBox_width_pixels.setSuffix(_translate("GroupBox_common_settings", " pixels", None))
        self.checkBox_tight_layout.setText(_translate("GroupBox_common_settings", "Tight Layout", None))
        self.label_2.setText(_translate("GroupBox_common_settings", "Width", None))
        self.groupBox_title.setTitle(_translate("GroupBox_common_settings", "Title", None))
        self.lineEdit_title.setToolTip(_translate("GroupBox_common_settings", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))
        self.label_5.setText(_translate("GroupBox_common_settings", "X", None))
        self.label.setToolTip(_translate("GroupBox_common_settings", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))
        self.label.setText(_translate("GroupBox_common_settings", "Plot Title", None))
        self.label_6.setText(_translate("GroupBox_common_settings", "Y", None))
        self.label_9.setText(_translate("GroupBox_common_settings", "Font Size", None))
        self.label_8.setText(_translate("GroupBox_common_settings", "Vertical Alignment", None))
        self.comboBox_horizontal_alignment.setItemText(0, _translate("GroupBox_common_settings", "Left", None))
        self.comboBox_horizontal_alignment.setItemText(1, _translate("GroupBox_common_settings", "Center", None))
        self.comboBox_horizontal_alignment.setItemText(2, _translate("GroupBox_common_settings", "Right", None))
        self.label_7.setText(_translate("GroupBox_common_settings", "Horizontal Alignment", None))
        self.comboBox_vertical_alignment.setItemText(0, _translate("GroupBox_common_settings", "Top", None))
        self.comboBox_vertical_alignment.setItemText(1, _translate("GroupBox_common_settings", "Center", None))
        self.comboBox_vertical_alignment.setItemText(2, _translate("GroupBox_common_settings", "Bottom", None))
        self.comboBox_vertical_alignment.setItemText(3, _translate("GroupBox_common_settings", "Baseline", None))
        self.spinBox_fontsize.setSuffix(_translate("GroupBox_common_settings", " points", None))
