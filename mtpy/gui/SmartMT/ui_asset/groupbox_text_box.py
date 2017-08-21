# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_text_box.ui'
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

class Ui_GroupBox_text_box(object):
    def setupUi(self, GroupBox_text_box):
        GroupBox_text_box.setObjectName(_fromUtf8("GroupBox_text_box"))
        GroupBox_text_box.resize(300, 300)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_text_box.sizePolicy().hasHeightForWidth())
        GroupBox_text_box.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_text_box)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox_location = QtGui.QGroupBox(GroupBox_text_box)
        self.groupBox_location.setFlat(True)
        self.groupBox_location.setCheckable(True)
        self.groupBox_location.setChecked(False)
        self.groupBox_location.setObjectName(_fromUtf8("groupBox_location"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox_location)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label_3 = QtGui.QLabel(self.groupBox_location)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)
        self.horizontalSlider_x = QtGui.QSlider(self.groupBox_location)
        self.horizontalSlider_x.setMaximum(100)
        self.horizontalSlider_x.setProperty("value", 90)
        self.horizontalSlider_x.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_x.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horizontalSlider_x.setTickInterval(10)
        self.horizontalSlider_x.setObjectName(_fromUtf8("horizontalSlider_x"))
        self.gridLayout.addWidget(self.horizontalSlider_x, 0, 1, 1, 1)
        self.doubleSpinBox_x = QtGui.QDoubleSpinBox(self.groupBox_location)
        self.doubleSpinBox_x.setMaximum(1.0)
        self.doubleSpinBox_x.setSingleStep(0.1)
        self.doubleSpinBox_x.setProperty("value", 0.9)
        self.doubleSpinBox_x.setObjectName(_fromUtf8("doubleSpinBox_x"))
        self.gridLayout.addWidget(self.doubleSpinBox_x, 0, 2, 1, 1)
        self.horizontalSlider_y = QtGui.QSlider(self.groupBox_location)
        self.horizontalSlider_y.setMaximum(100)
        self.horizontalSlider_y.setProperty("value", 5)
        self.horizontalSlider_y.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_y.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horizontalSlider_y.setTickInterval(10)
        self.horizontalSlider_y.setObjectName(_fromUtf8("horizontalSlider_y"))
        self.gridLayout.addWidget(self.horizontalSlider_y, 1, 1, 1, 1)
        self.doubleSpinBox_y = QtGui.QDoubleSpinBox(self.groupBox_location)
        self.doubleSpinBox_y.setMaximum(1.0)
        self.doubleSpinBox_y.setSingleStep(0.1)
        self.doubleSpinBox_y.setProperty("value", 0.05)
        self.doubleSpinBox_y.setObjectName(_fromUtf8("doubleSpinBox_y"))
        self.gridLayout.addWidget(self.doubleSpinBox_y, 1, 2, 1, 1)
        self.label_2 = QtGui.QLabel(self.groupBox_location)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_location)
        self.groupBox_2 = QtGui.QGroupBox(GroupBox_text_box)
        self.groupBox_2.setFlat(True)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.gridLayout_2 = QtGui.QGridLayout(self.groupBox_2)
        self.gridLayout_2.setObjectName(_fromUtf8("gridLayout_2"))
        self.checkBox_size = QtGui.QCheckBox(self.groupBox_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_size.sizePolicy().hasHeightForWidth())
        self.checkBox_size.setSizePolicy(sizePolicy)
        self.checkBox_size.setObjectName(_fromUtf8("checkBox_size"))
        self.gridLayout_2.addWidget(self.checkBox_size, 0, 0, 1, 1)
        self.comboBox_size = QtGui.QComboBox(self.groupBox_2)
        self.comboBox_size.setEnabled(False)
        self.comboBox_size.setObjectName(_fromUtf8("comboBox_size"))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.comboBox_size.addItem(_fromUtf8(""))
        self.gridLayout_2.addWidget(self.comboBox_size, 0, 1, 1, 1)
        self.spinBox_size = QtGui.QSpinBox(self.groupBox_2)
        self.spinBox_size.setEnabled(False)
        self.spinBox_size.setMinimum(4)
        self.spinBox_size.setProperty("value", 10)
        self.spinBox_size.setObjectName(_fromUtf8("spinBox_size"))
        self.gridLayout_2.addWidget(self.spinBox_size, 0, 2, 1, 1)
        self.checkBox_weight = QtGui.QCheckBox(self.groupBox_2)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_weight.sizePolicy().hasHeightForWidth())
        self.checkBox_weight.setSizePolicy(sizePolicy)
        self.checkBox_weight.setObjectName(_fromUtf8("checkBox_weight"))
        self.gridLayout_2.addWidget(self.checkBox_weight, 1, 0, 1, 1)
        self.comboBox_weight = QtGui.QComboBox(self.groupBox_2)
        self.comboBox_weight.setEnabled(False)
        self.comboBox_weight.setObjectName(_fromUtf8("comboBox_weight"))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.comboBox_weight.addItem(_fromUtf8(""))
        self.gridLayout_2.addWidget(self.comboBox_weight, 1, 1, 1, 2)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.groupBox_padding = QtGui.QGroupBox(GroupBox_text_box)
        self.groupBox_padding.setFlat(True)
        self.groupBox_padding.setCheckable(True)
        self.groupBox_padding.setChecked(False)
        self.groupBox_padding.setObjectName(_fromUtf8("groupBox_padding"))
        self.gridLayout_3 = QtGui.QGridLayout(self.groupBox_padding)
        self.gridLayout_3.setObjectName(_fromUtf8("gridLayout_3"))
        self.label_5 = QtGui.QLabel(self.groupBox_padding)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.gridLayout_3.addWidget(self.label_5, 0, 0, 1, 1)
        self.horizontalSlider_x_pad = QtGui.QSlider(self.groupBox_padding)
        self.horizontalSlider_x_pad.setMaximum(100)
        self.horizontalSlider_x_pad.setProperty("value", 95)
        self.horizontalSlider_x_pad.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_x_pad.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horizontalSlider_x_pad.setTickInterval(10)
        self.horizontalSlider_x_pad.setObjectName(_fromUtf8("horizontalSlider_x_pad"))
        self.gridLayout_3.addWidget(self.horizontalSlider_x_pad, 0, 1, 1, 1)
        self.doubleSpinBox_x_pad = QtGui.QDoubleSpinBox(self.groupBox_padding)
        self.doubleSpinBox_x_pad.setMaximum(1.0)
        self.doubleSpinBox_x_pad.setSingleStep(0.1)
        self.doubleSpinBox_x_pad.setProperty("value", 0.95)
        self.doubleSpinBox_x_pad.setObjectName(_fromUtf8("doubleSpinBox_x_pad"))
        self.gridLayout_3.addWidget(self.doubleSpinBox_x_pad, 0, 2, 1, 1)
        self.label_4 = QtGui.QLabel(self.groupBox_padding)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.gridLayout_3.addWidget(self.label_4, 1, 0, 1, 1)
        self.horizontalSlider_y_pad = QtGui.QSlider(self.groupBox_padding)
        self.horizontalSlider_y_pad.setMaximum(100)
        self.horizontalSlider_y_pad.setProperty("value", 95)
        self.horizontalSlider_y_pad.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_y_pad.setTickPosition(QtGui.QSlider.TicksBelow)
        self.horizontalSlider_y_pad.setTickInterval(10)
        self.horizontalSlider_y_pad.setObjectName(_fromUtf8("horizontalSlider_y_pad"))
        self.gridLayout_3.addWidget(self.horizontalSlider_y_pad, 1, 1, 1, 1)
        self.doubleSpinBox_y_pad = QtGui.QDoubleSpinBox(self.groupBox_padding)
        self.doubleSpinBox_y_pad.setMaximum(1.0)
        self.doubleSpinBox_y_pad.setSingleStep(0.1)
        self.doubleSpinBox_y_pad.setProperty("value", 0.95)
        self.doubleSpinBox_y_pad.setObjectName(_fromUtf8("doubleSpinBox_y_pad"))
        self.gridLayout_3.addWidget(self.doubleSpinBox_y_pad, 1, 2, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_padding)
        self.label_3.setBuddy(self.horizontalSlider_y)
        self.label_2.setBuddy(self.horizontalSlider_x)
        self.label_5.setBuddy(self.horizontalSlider_x_pad)
        self.label_4.setBuddy(self.horizontalSlider_y_pad)

        self.retranslateUi(GroupBox_text_box)
        self.comboBox_size.setCurrentIndex(2)
        self.comboBox_weight.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_text_box)
        GroupBox_text_box.setTabOrder(self.groupBox_location, self.horizontalSlider_x)
        GroupBox_text_box.setTabOrder(self.horizontalSlider_x, self.doubleSpinBox_x)
        GroupBox_text_box.setTabOrder(self.doubleSpinBox_x, self.horizontalSlider_y)
        GroupBox_text_box.setTabOrder(self.horizontalSlider_y, self.doubleSpinBox_y)
        GroupBox_text_box.setTabOrder(self.doubleSpinBox_y, self.checkBox_size)
        GroupBox_text_box.setTabOrder(self.checkBox_size, self.comboBox_size)
        GroupBox_text_box.setTabOrder(self.comboBox_size, self.spinBox_size)
        GroupBox_text_box.setTabOrder(self.spinBox_size, self.checkBox_weight)
        GroupBox_text_box.setTabOrder(self.checkBox_weight, self.comboBox_weight)
        GroupBox_text_box.setTabOrder(self.comboBox_weight, self.groupBox_padding)
        GroupBox_text_box.setTabOrder(self.groupBox_padding, self.horizontalSlider_x_pad)
        GroupBox_text_box.setTabOrder(self.horizontalSlider_x_pad, self.doubleSpinBox_x_pad)
        GroupBox_text_box.setTabOrder(self.doubleSpinBox_x_pad, self.horizontalSlider_y_pad)
        GroupBox_text_box.setTabOrder(self.horizontalSlider_y_pad, self.doubleSpinBox_y_pad)

    def retranslateUi(self, GroupBox_text_box):
        GroupBox_text_box.setWindowTitle(_translate("GroupBox_text_box", "GroupBox", None))
        GroupBox_text_box.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>Set text box parameters</p><p><br/></p></body></html>", None))
        GroupBox_text_box.setTitle(_translate("GroupBox_text_box", "Text Box", None))
        self.groupBox_location.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>Location for text label for each\n"
"                            resistivity subplot.</p><p>Location is in relative coordinates of the datas.</p><p><br/></p></body></html>\n"
"                        ", None))
        self.groupBox_location.setTitle(_translate("GroupBox_text_box", "Location", None))
        self.label_3.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))
        self.label_3.setText(_translate("GroupBox_text_box", "Y", None))
        self.horizontalSlider_x.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.doubleSpinBox_x.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.horizontalSlider_y.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))
        self.doubleSpinBox_y.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))
        self.label_2.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.label_2.setText(_translate("GroupBox_text_box", "X", None))
        self.groupBox_2.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p><br/></p></body></html>", None))
        self.groupBox_2.setTitle(_translate("GroupBox_text_box", "Font", None))
        self.checkBox_size.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>Size of text for subplot\n"
"                                        labels</p></body></html>\n"
"                                    ", None))
        self.checkBox_size.setText(_translate("GroupBox_text_box", "Size", None))
        self.comboBox_size.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>Size of text for subplot\n"
"                                        labels</p></body></html>\n"
"                                    ", None))
        self.comboBox_size.setItemText(0, _translate("GroupBox_text_box", "xx-small", None))
        self.comboBox_size.setItemText(1, _translate("GroupBox_text_box", "x-small", None))
        self.comboBox_size.setItemText(2, _translate("GroupBox_text_box", "small", None))
        self.comboBox_size.setItemText(3, _translate("GroupBox_text_box", "medium", None))
        self.comboBox_size.setItemText(4, _translate("GroupBox_text_box", "large", None))
        self.comboBox_size.setItemText(5, _translate("GroupBox_text_box", "x-large", None))
        self.comboBox_size.setItemText(6, _translate("GroupBox_text_box", "xx-large", None))
        self.comboBox_size.setItemText(7, _translate("GroupBox_text_box", "size in points", None))
        self.spinBox_size.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>Size of text for subplot\n"
"                                        labels</p></body></html>\n"
"                                    ", None))
        self.spinBox_size.setSuffix(_translate("GroupBox_text_box", " points", None))
        self.checkBox_weight.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>weight of text label font</p></body></html>", None))
        self.checkBox_weight.setText(_translate("GroupBox_text_box", "Weight", None))
        self.comboBox_weight.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>weight of text label font</p></body></html>", None))
        self.comboBox_weight.setItemText(0, _translate("GroupBox_text_box", "light", None))
        self.comboBox_weight.setItemText(1, _translate("GroupBox_text_box", "normal", None))
        self.comboBox_weight.setItemText(2, _translate("GroupBox_text_box", "medium", None))
        self.comboBox_weight.setItemText(3, _translate("GroupBox_text_box", "semibold", None))
        self.comboBox_weight.setItemText(4, _translate("GroupBox_text_box", "bold", None))
        self.comboBox_weight.setItemText(5, _translate("GroupBox_text_box", "heavy", None))
        self.comboBox_weight.setItemText(6, _translate("GroupBox_text_box", "black", None))
        self.groupBox_padding.setTitle(_translate("GroupBox_text_box", "Paddings", None))
        self.label_5.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.label_5.setText(_translate("GroupBox_text_box", "X", None))
        self.horizontalSlider_x_pad.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.doubleSpinBox_x_pad.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>lateral position of the left\n"
"                                        hand corner of the color bar in figure between [0, 1], 0 is the left side.</p></body></html>\n"
"                                    ", None))
        self.label_4.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))
        self.label_4.setText(_translate("GroupBox_text_box", "Y", None))
        self.horizontalSlider_y_pad.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))
        self.doubleSpinBox_y_pad.setToolTip(_translate("GroupBox_text_box", "<html><head/><body><p>vertical position of the\n"
"                                        bottom of the color bar in figure between [0, 1], 0 is bottom side.</p></body></html>\n"
"                                    ", None))

