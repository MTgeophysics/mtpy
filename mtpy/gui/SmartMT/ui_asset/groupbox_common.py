# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_common.ui'
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
        GroupBox_common_settings.resize(300, 131)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_common_settings.sizePolicy().hasHeightForWidth())
        GroupBox_common_settings.setSizePolicy(sizePolicy)
        self.gridLayout = QtGui.QGridLayout(GroupBox_common_settings)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.lineEdit_title = QtGui.QLineEdit(GroupBox_common_settings)
        self.lineEdit_title.setObjectName(_fromUtf8("lineEdit_title"))
        self.gridLayout.addWidget(self.lineEdit_title, 0, 1, 1, 2)
        self.label = QtGui.QLabel(GroupBox_common_settings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_2 = QtGui.QLabel(GroupBox_common_settings)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.spinBox_width_pixels = QtGui.QSpinBox(GroupBox_common_settings)
        self.spinBox_width_pixels.setMinimum(1)
        self.spinBox_width_pixels.setMaximum(1000000)
        self.spinBox_width_pixels.setProperty("value", 800)
        self.spinBox_width_pixels.setObjectName(_fromUtf8("spinBox_width_pixels"))
        self.gridLayout.addWidget(self.spinBox_width_pixels, 1, 2, 1, 1)
        self.doubleSpinBox_width_inches = QtGui.QDoubleSpinBox(GroupBox_common_settings)
        self.doubleSpinBox_width_inches.setMinimum(1.0)
        self.doubleSpinBox_width_inches.setMaximum(1000.0)
        self.doubleSpinBox_width_inches.setProperty("value", 8.0)
        self.doubleSpinBox_width_inches.setObjectName(_fromUtf8("doubleSpinBox_width_inches"))
        self.gridLayout.addWidget(self.doubleSpinBox_width_inches, 1, 1, 1, 1)
        self.label_3 = QtGui.QLabel(GroupBox_common_settings)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 2, 0, 1, 1)
        self.doubleSpinBox_height_inches = QtGui.QDoubleSpinBox(GroupBox_common_settings)
        self.doubleSpinBox_height_inches.setMinimum(1.0)
        self.doubleSpinBox_height_inches.setMaximum(1000.0)
        self.doubleSpinBox_height_inches.setProperty("value", 6.0)
        self.doubleSpinBox_height_inches.setObjectName(_fromUtf8("doubleSpinBox_height_inches"))
        self.gridLayout.addWidget(self.doubleSpinBox_height_inches, 2, 1, 1, 1)
        self.spinBox_height_pixels = QtGui.QSpinBox(GroupBox_common_settings)
        self.spinBox_height_pixels.setMinimum(1)
        self.spinBox_height_pixels.setMaximum(1000000)
        self.spinBox_height_pixels.setProperty("value", 600)
        self.spinBox_height_pixels.setObjectName(_fromUtf8("spinBox_height_pixels"))
        self.gridLayout.addWidget(self.spinBox_height_pixels, 2, 2, 1, 1)
        self.label_4 = QtGui.QLabel(GroupBox_common_settings)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)
        self.spinBox_dpi = QtGui.QSpinBox(GroupBox_common_settings)
        self.spinBox_dpi.setMinimum(1)
        self.spinBox_dpi.setMaximum(1000)
        self.spinBox_dpi.setProperty("value", 100)
        self.spinBox_dpi.setObjectName(_fromUtf8("spinBox_dpi"))
        self.gridLayout.addWidget(self.spinBox_dpi, 3, 1, 1, 1)

        self.retranslateUi(GroupBox_common_settings)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_common_settings)

    def retranslateUi(self, GroupBox_common_settings):
        GroupBox_common_settings.setWindowTitle(_translate("GroupBox_common_settings", "GroupBox", None))
        GroupBox_common_settings.setTitle(_translate("GroupBox_common_settings", "General Settings", None))
        self.lineEdit_title.setToolTip(_translate("GroupBox_common_settings", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))
        self.label.setToolTip(_translate("GroupBox_common_settings", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))
        self.label.setText(_translate("GroupBox_common_settings", "Plot Title", None))
        self.label_2.setText(_translate("GroupBox_common_settings", "Width", None))
        self.spinBox_width_pixels.setSuffix(_translate("GroupBox_common_settings", " pixels", None))
        self.doubleSpinBox_width_inches.setSuffix(_translate("GroupBox_common_settings", " inches", None))
        self.label_3.setText(_translate("GroupBox_common_settings", "Height", None))
        self.doubleSpinBox_height_inches.setSuffix(_translate("GroupBox_common_settings", " inches", None))
        self.spinBox_height_pixels.setSuffix(_translate("GroupBox_common_settings", " pixels", None))
        self.label_4.setText(_translate("GroupBox_common_settings", "DPI", None))

