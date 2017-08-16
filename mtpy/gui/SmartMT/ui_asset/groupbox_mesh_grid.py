# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_mesh_grid.ui'
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

class Ui_GroupBox_mash_grid(object):
    def setupUi(self, GroupBox_mash_grid):
        GroupBox_mash_grid.setObjectName(_fromUtf8("GroupBox_mash_grid"))
        GroupBox_mash_grid.resize(300, 300)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_mash_grid.sizePolicy().hasHeightForWidth())
        GroupBox_mash_grid.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_mash_grid)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox = QtGui.QGroupBox(GroupBox_mash_grid)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.groupBox)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.radioButton_imshow = QtGui.QRadioButton(self.groupBox)
        self.radioButton_imshow.setChecked(True)
        self.radioButton_imshow.setObjectName(_fromUtf8("radioButton_imshow"))
        self.horizontalLayout.addWidget(self.radioButton_imshow)
        self.radioButton_pcolormesh = QtGui.QRadioButton(self.groupBox)
        self.radioButton_pcolormesh.setObjectName(_fromUtf8("radioButton_pcolormesh"))
        self.horizontalLayout.addWidget(self.radioButton_pcolormesh)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_interpolation_method = QtGui.QGroupBox(GroupBox_mash_grid)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_interpolation_method.sizePolicy().hasHeightForWidth())
        self.groupBox_interpolation_method.setSizePolicy(sizePolicy)
        self.groupBox_interpolation_method.setObjectName(_fromUtf8("groupBox_interpolation_method"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.groupBox_interpolation_method)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.comboBox_interpolation_method = QtGui.QComboBox(self.groupBox_interpolation_method)
        self.comboBox_interpolation_method.setObjectName(_fromUtf8("comboBox_interpolation_method"))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.comboBox_interpolation_method.addItem(_fromUtf8(""))
        self.verticalLayout_2.addWidget(self.comboBox_interpolation_method)
        self.groupBox.raise_()
        self.comboBox_interpolation_method.raise_()
        self.verticalLayout.addWidget(self.groupBox_interpolation_method)

        self.retranslateUi(GroupBox_mash_grid)
        self.comboBox_interpolation_method.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_mash_grid)

    def retranslateUi(self, GroupBox_mash_grid):
        GroupBox_mash_grid.setWindowTitle(_translate("GroupBox_mash_grid", "GroupBox", None))
        GroupBox_mash_grid.setTitle(_translate("GroupBox_mash_grid", "Mesh Grid", None))
        self.groupBox.setToolTip(_translate("GroupBox_mash_grid", "<html><head/><body><p>Type of Grids in the plot</p></body></html>", None))
        self.groupBox.setTitle(_translate("GroupBox_mash_grid", "Grid Type", None))
        self.radioButton_imshow.setToolTip(_translate("GroupBox_mash_grid", "<html><head/><body><p>plots the data as an image and\n"
"                                        can be interpolated, though the image is sstretched to the station spacing and\n"
"                                        plot_period, the cells remain of equal size, so the interpolation might be a\n"
"                                        little skewed.</p></body></html>\n"
"                                    ", None))
        self.radioButton_imshow.setText(_translate("GroupBox_mash_grid", "imshow", None))
        self.radioButton_pcolormesh.setToolTip(_translate("GroupBox_mash_grid", "<html><head/><body><p>plot the data on an irregular\n"
"                                        grid, but with no interpolation, which results an accurate location of\n"
"                                        resistivity values.</p></body></html>\n"
"                                    ", None))
        self.radioButton_pcolormesh.setText(_translate("GroupBox_mash_grid", "pcolormesh", None))
        self.groupBox_interpolation_method.setToolTip(_translate("GroupBox_mash_grid", "<html><head/><body><p>defines the interpolation method if gride\n"
"                            style is \'imshow\'.</p><p><br/></p><p>\'Nearest\' is same as\n"
"                            pcolormesh except the lateral boxes are equal size instead of set in a grid like pcolormesh.\n"
"                            imshow just gives a smoother interpretation of pseudosection.</p></body></html>\n"
"                        ", None))
        self.groupBox_interpolation_method.setTitle(_translate("GroupBox_mash_grid", "Interpolation Method", None))
        self.comboBox_interpolation_method.setItemText(0, _translate("GroupBox_mash_grid", "None", None))
        self.comboBox_interpolation_method.setItemText(1, _translate("GroupBox_mash_grid", "Nearest", None))
        self.comboBox_interpolation_method.setItemText(2, _translate("GroupBox_mash_grid", "Bilinear", None))
        self.comboBox_interpolation_method.setItemText(3, _translate("GroupBox_mash_grid", "Bicubic", None))
        self.comboBox_interpolation_method.setItemText(4, _translate("GroupBox_mash_grid", "Spline16", None))
        self.comboBox_interpolation_method.setItemText(5, _translate("GroupBox_mash_grid", "Spline36", None))
        self.comboBox_interpolation_method.setItemText(6, _translate("GroupBox_mash_grid", "Hanning", None))
        self.comboBox_interpolation_method.setItemText(7, _translate("GroupBox_mash_grid", "Hamming", None))
        self.comboBox_interpolation_method.setItemText(8, _translate("GroupBox_mash_grid", "Hermite", None))
        self.comboBox_interpolation_method.setItemText(9, _translate("GroupBox_mash_grid", "Kaiser", None))
        self.comboBox_interpolation_method.setItemText(10, _translate("GroupBox_mash_grid", "Guadric", None))
        self.comboBox_interpolation_method.setItemText(11, _translate("GroupBox_mash_grid", "Catrom", None))
        self.comboBox_interpolation_method.setItemText(12, _translate("GroupBox_mash_grid", "Gaussian", None))
        self.comboBox_interpolation_method.setItemText(13, _translate("GroupBox_mash_grid", "Bessel", None))
        self.comboBox_interpolation_method.setItemText(14, _translate("GroupBox_mash_grid", "Mitchell", None))
        self.comboBox_interpolation_method.setItemText(15, _translate("GroupBox_mash_grid", "Sinc", None))
        self.comboBox_interpolation_method.setItemText(16, _translate("GroupBox_mash_grid", "Lanczos", None))

