# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'fifthgui.ui'
#
# Created: Thu Nov 21 11:26:41 2013
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s


class Ui_occamgui2D(object):
    def setupUi(self, occamgui2D):
        occamgui2D.setObjectName(_fromUtf8("occamgui2D"))
        occamgui2D.resize(906, 699)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Maximum
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(occamgui2D.sizePolicy().hasHeightForWidth())
        occamgui2D.setSizePolicy(sizePolicy)
        self.label_29 = QtGui.QLabel(occamgui2D)
        self.label_29.setGeometry(QtCore.QRect(30, 90, 151, 31))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.label_29.setFont(font)
        self.label_29.setObjectName(_fromUtf8("label_29"))
        self.label_4 = QtGui.QLabel(occamgui2D)
        self.label_4.setGeometry(QtCore.QRect(20, 130, 151, 31))
        font = QtGui.QFont()
        font.setFamily(_fromUtf8("Arial"))
        font.setPointSize(18)
        font.setBold(True)
        font.setWeight(75)
        self.label_4.setFont(font)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.pushButton_generateinputfile = QtGui.QPushButton(occamgui2D)
        self.pushButton_generateinputfile.setGeometry(QtCore.QRect(230, 660, 138, 30))
        self.pushButton_generateinputfile.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_generateinputfile.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_generateinputfile.setAutoDefault(False)
        self.pushButton_generateinputfile.setObjectName(
            _fromUtf8("pushButton_generateinputfile")
        )
        self.pushButton_quit = QtGui.QPushButton(occamgui2D)
        self.pushButton_quit.setGeometry(QtCore.QRect(560, 660, 158, 30))
        self.pushButton_quit.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_quit.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_quit.setAutoDefault(False)
        self.pushButton_quit.setObjectName(_fromUtf8("pushButton_quit"))
        self.pushButton_checkparameter = QtGui.QPushButton(occamgui2D)
        self.pushButton_checkparameter.setGeometry(QtCore.QRect(39, 660, 129, 30))
        self.pushButton_checkparameter.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_checkparameter.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_checkparameter.setAutoDefault(False)
        self.pushButton_checkparameter.setObjectName(
            _fromUtf8("pushButton_checkparameter")
        )
        self.pushButton_runoccam = QtGui.QPushButton(occamgui2D)
        self.pushButton_runoccam.setEnabled(False)
        self.pushButton_runoccam.setGeometry(QtCore.QRect(771, 660, 110, 30))
        self.pushButton_runoccam.setMinimumSize(QtCore.QSize(0, 30))
        self.pushButton_runoccam.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_runoccam.setCheckable(False)
        self.pushButton_runoccam.setChecked(False)
        self.pushButton_runoccam.setAutoExclusive(False)
        self.pushButton_runoccam.setAutoDefault(False)
        self.pushButton_runoccam.setDefault(False)
        self.pushButton_runoccam.setFlat(False)
        self.pushButton_runoccam.setObjectName(_fromUtf8("pushButton_runoccam"))
        self.layoutWidget = QtGui.QWidget(occamgui2D)
        self.layoutWidget.setGeometry(QtCore.QRect(177, 66, 691, 263))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.formLayout_9 = QtGui.QFormLayout(self.layoutWidget)
        self.formLayout_9.setFieldGrowthPolicy(QtGui.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout_9.setMargin(0)
        self.formLayout_9.setObjectName(_fromUtf8("formLayout_9"))
        self.horizontalLayout_8 = QtGui.QHBoxLayout()
        self.horizontalLayout_8.setObjectName(_fromUtf8("horizontalLayout_8"))
        self.label_2 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.horizontalLayout_8.addWidget(self.label_2)
        spacerItem = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_8.addItem(spacerItem)
        self.formLayout_9.setLayout(
            0, QtGui.QFormLayout.LabelRole, self.horizontalLayout_8
        )
        self.horizontalLayout = QtGui.QHBoxLayout()
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        spacerItem1 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout.addItem(spacerItem1)
        self.button_browse_wd = QtGui.QToolButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.button_browse_wd.sizePolicy().hasHeightForWidth()
        )
        self.button_browse_wd.setSizePolicy(sizePolicy)
        self.button_browse_wd.setMinimumSize(QtCore.QSize(30, 25))
        self.button_browse_wd.setMaximumSize(QtCore.QSize(30, 25))
        self.button_browse_wd.setObjectName(_fromUtf8("button_browse_wd"))
        self.horizontalLayout.addWidget(self.button_browse_wd)
        self.lineEdit_browse_wd = QtGui.QLineEdit(self.layoutWidget)
        self.lineEdit_browse_wd.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_wd.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_wd.setObjectName(_fromUtf8("lineEdit_browse_wd"))
        self.horizontalLayout.addWidget(self.lineEdit_browse_wd)
        self.formLayout_9.setLayout(
            0, QtGui.QFormLayout.FieldRole, self.horizontalLayout
        )
        self.horizontalLayout_9 = QtGui.QHBoxLayout()
        self.horizontalLayout_9.setObjectName(_fromUtf8("horizontalLayout_9"))
        self.label_5 = QtGui.QLabel(self.layoutWidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.horizontalLayout_9.addWidget(self.label_5)
        spacerItem2 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_9.addItem(spacerItem2)
        self.formLayout_9.setLayout(
            1, QtGui.QFormLayout.LabelRole, self.horizontalLayout_9
        )
        self.horizontalLayout_2 = QtGui.QHBoxLayout()
        self.horizontalLayout_2.setObjectName(_fromUtf8("horizontalLayout_2"))
        spacerItem3 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_2.addItem(spacerItem3)
        self.button_browse_occam = QtGui.QToolButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.button_browse_occam.sizePolicy().hasHeightForWidth()
        )
        self.button_browse_occam.setSizePolicy(sizePolicy)
        self.button_browse_occam.setMinimumSize(QtCore.QSize(30, 25))
        self.button_browse_occam.setMaximumSize(QtCore.QSize(30, 25))
        self.button_browse_occam.setObjectName(_fromUtf8("button_browse_occam"))
        self.horizontalLayout_2.addWidget(self.button_browse_occam)
        self.lineEdit_browse_occam = QtGui.QLineEdit(self.layoutWidget)
        self.lineEdit_browse_occam.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_occam.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_occam.setObjectName(_fromUtf8("lineEdit_browse_occam"))
        self.horizontalLayout_2.addWidget(self.lineEdit_browse_occam)
        self.formLayout_9.setLayout(
            1, QtGui.QFormLayout.FieldRole, self.horizontalLayout_2
        )
        self.horizontalLayout_10 = QtGui.QHBoxLayout()
        self.horizontalLayout_10.setObjectName(_fromUtf8("horizontalLayout_10"))
        self.label_28 = QtGui.QLabel(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_28.sizePolicy().hasHeightForWidth())
        self.label_28.setSizePolicy(sizePolicy)
        self.label_28.setObjectName(_fromUtf8("label_28"))
        self.horizontalLayout_10.addWidget(self.label_28)
        spacerItem4 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_10.addItem(spacerItem4)
        self.formLayout_9.setLayout(
            2, QtGui.QFormLayout.LabelRole, self.horizontalLayout_10
        )
        self.horizontalLayout_3 = QtGui.QHBoxLayout()
        self.horizontalLayout_3.setObjectName(_fromUtf8("horizontalLayout_3"))
        spacerItem5 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_3.addItem(spacerItem5)
        self.button_browse_edis = QtGui.QToolButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.button_browse_edis.sizePolicy().hasHeightForWidth()
        )
        self.button_browse_edis.setSizePolicy(sizePolicy)
        self.button_browse_edis.setMinimumSize(QtCore.QSize(30, 25))
        self.button_browse_edis.setMaximumSize(QtCore.QSize(30, 25))
        self.button_browse_edis.setObjectName(_fromUtf8("button_browse_edis"))
        self.horizontalLayout_3.addWidget(self.button_browse_edis)
        self.lineEdit_browse_edi = QtGui.QLineEdit(self.layoutWidget)
        self.lineEdit_browse_edi.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_edi.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_edi.setObjectName(_fromUtf8("lineEdit_browse_edi"))
        self.horizontalLayout_3.addWidget(self.lineEdit_browse_edi)
        self.formLayout_9.setLayout(
            2, QtGui.QFormLayout.FieldRole, self.horizontalLayout_3
        )
        self.horizontalLayout_11 = QtGui.QHBoxLayout()
        self.horizontalLayout_11.setObjectName(_fromUtf8("horizontalLayout_11"))
        self.checkBox_usestationlist = QtGui.QCheckBox(self.layoutWidget)
        self.checkBox_usestationlist.setObjectName(_fromUtf8("checkBox_usestationlist"))
        self.horizontalLayout_11.addWidget(self.checkBox_usestationlist)
        spacerItem6 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_11.addItem(spacerItem6)
        self.formLayout_9.setLayout(
            3, QtGui.QFormLayout.LabelRole, self.horizontalLayout_11
        )
        self.horizontalLayout_4 = QtGui.QHBoxLayout()
        self.horizontalLayout_4.setObjectName(_fromUtf8("horizontalLayout_4"))
        spacerItem7 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_4.addItem(spacerItem7)
        self.pushButton_loadstations = QtGui.QPushButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.pushButton_loadstations.sizePolicy().hasHeightForWidth()
        )
        self.pushButton_loadstations.setSizePolicy(sizePolicy)
        self.pushButton_loadstations.setMinimumSize(QtCore.QSize(0, 25))
        self.pushButton_loadstations.setMaximumSize(QtCore.QSize(3000, 30))
        self.pushButton_loadstations.setAutoDefault(True)
        self.pushButton_loadstations.setObjectName(_fromUtf8("pushButton_loadstations"))
        self.horizontalLayout_4.addWidget(self.pushButton_loadstations)
        self.lineEdit_browse_stationfile = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_browse_stationfile.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_browse_stationfile.setSizePolicy(sizePolicy)
        self.lineEdit_browse_stationfile.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_stationfile.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_stationfile.setObjectName(
            _fromUtf8("lineEdit_browse_stationfile")
        )
        self.horizontalLayout_4.addWidget(self.lineEdit_browse_stationfile)
        self.formLayout_9.setLayout(
            3, QtGui.QFormLayout.FieldRole, self.horizontalLayout_4
        )
        self.horizontalLayout_12 = QtGui.QHBoxLayout()
        self.horizontalLayout_12.setObjectName(_fromUtf8("horizontalLayout_12"))
        self.checkBox_usedatafile = QtGui.QCheckBox(self.layoutWidget)
        self.checkBox_usedatafile.setEnabled(False)
        self.checkBox_usedatafile.setObjectName(_fromUtf8("checkBox_usedatafile"))
        self.horizontalLayout_12.addWidget(self.checkBox_usedatafile)
        spacerItem8 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_12.addItem(spacerItem8)
        self.formLayout_9.setLayout(
            4, QtGui.QFormLayout.LabelRole, self.horizontalLayout_12
        )
        self.horizontalLayout_5 = QtGui.QHBoxLayout()
        self.horizontalLayout_5.setObjectName(_fromUtf8("horizontalLayout_5"))
        spacerItem9 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_5.addItem(spacerItem9)
        self.pushButton_loaddatafile = QtGui.QPushButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.pushButton_loaddatafile.sizePolicy().hasHeightForWidth()
        )
        self.pushButton_loaddatafile.setSizePolicy(sizePolicy)
        self.pushButton_loaddatafile.setMinimumSize(QtCore.QSize(0, 25))
        self.pushButton_loaddatafile.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_loaddatafile.setAutoDefault(True)
        self.pushButton_loaddatafile.setObjectName(_fromUtf8("pushButton_loaddatafile"))
        self.horizontalLayout_5.addWidget(self.pushButton_loaddatafile)
        self.lineEdit_browse_datafile = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_browse_datafile.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_browse_datafile.setSizePolicy(sizePolicy)
        self.lineEdit_browse_datafile.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_datafile.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_datafile.setObjectName(
            _fromUtf8("lineEdit_browse_datafile")
        )
        self.horizontalLayout_5.addWidget(self.lineEdit_browse_datafile)
        self.formLayout_9.setLayout(
            4, QtGui.QFormLayout.FieldRole, self.horizontalLayout_5
        )
        self.horizontalLayout_13 = QtGui.QHBoxLayout()
        self.horizontalLayout_13.setObjectName(_fromUtf8("horizontalLayout_13"))
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.horizontalLayout_13.addWidget(self.label_3)
        spacerItem10 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_13.addItem(spacerItem10)
        self.formLayout_9.setLayout(
            5, QtGui.QFormLayout.LabelRole, self.horizontalLayout_13
        )
        self.horizontalLayout_6 = QtGui.QHBoxLayout()
        self.horizontalLayout_6.setObjectName(_fromUtf8("horizontalLayout_6"))
        spacerItem11 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_6.addItem(spacerItem11)
        self.lineEdit_datafilename = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_datafilename.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_datafilename.setSizePolicy(sizePolicy)
        self.lineEdit_datafilename.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_datafilename.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_datafilename.setObjectName(_fromUtf8("lineEdit_datafilename"))
        self.horizontalLayout_6.addWidget(self.lineEdit_datafilename)
        self.formLayout_9.setLayout(
            5, QtGui.QFormLayout.FieldRole, self.horizontalLayout_6
        )
        self.horizontalLayout_14 = QtGui.QHBoxLayout()
        self.horizontalLayout_14.setObjectName(_fromUtf8("horizontalLayout_14"))
        self.checkBox_useiterationfile = QtGui.QCheckBox(self.layoutWidget)
        self.checkBox_useiterationfile.setObjectName(
            _fromUtf8("checkBox_useiterationfile")
        )
        self.horizontalLayout_14.addWidget(self.checkBox_useiterationfile)
        spacerItem12 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_14.addItem(spacerItem12)
        self.formLayout_9.setLayout(
            6, QtGui.QFormLayout.LabelRole, self.horizontalLayout_14
        )
        self.horizontalLayout_7 = QtGui.QHBoxLayout()
        self.horizontalLayout_7.setObjectName(_fromUtf8("horizontalLayout_7"))
        spacerItem13 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_7.addItem(spacerItem13)
        self.pushButton_loaditerationfile = QtGui.QPushButton(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.pushButton_loaditerationfile.sizePolicy().hasHeightForWidth()
        )
        self.pushButton_loaditerationfile.setSizePolicy(sizePolicy)
        self.pushButton_loaditerationfile.setMinimumSize(QtCore.QSize(0, 25))
        self.pushButton_loaditerationfile.setMaximumSize(QtCore.QSize(16777215, 30))
        self.pushButton_loaditerationfile.setAutoDefault(True)
        self.pushButton_loaditerationfile.setObjectName(
            _fromUtf8("pushButton_loaditerationfile")
        )
        self.horizontalLayout_7.addWidget(self.pushButton_loaditerationfile)
        self.lineEdit_browse_iterationfile = QtGui.QLineEdit(self.layoutWidget)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_browse_iterationfile.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_browse_iterationfile.setSizePolicy(sizePolicy)
        self.lineEdit_browse_iterationfile.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_iterationfile.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_iterationfile.setObjectName(
            _fromUtf8("lineEdit_browse_iterationfile")
        )
        self.horizontalLayout_7.addWidget(self.lineEdit_browse_iterationfile)
        self.formLayout_9.setLayout(
            6, QtGui.QFormLayout.FieldRole, self.horizontalLayout_7
        )
        self.layoutWidget1 = QtGui.QWidget(occamgui2D)
        self.layoutWidget1.setGeometry(QtCore.QRect(120, 517, 672, 34))
        self.layoutWidget1.setObjectName(_fromUtf8("layoutWidget1"))
        self.horizontalLayout_30 = QtGui.QHBoxLayout(self.layoutWidget1)
        self.horizontalLayout_30.setMargin(0)
        self.horizontalLayout_30.setObjectName(_fromUtf8("horizontalLayout_30"))
        self.horizontalLayout_27 = QtGui.QHBoxLayout()
        self.horizontalLayout_27.setObjectName(_fromUtf8("horizontalLayout_27"))
        self.label_18 = QtGui.QLabel(self.layoutWidget1)
        self.label_18.setObjectName(_fromUtf8("label_18"))
        self.horizontalLayout_27.addWidget(self.label_18)
        self.spinBox_no_layers = QtGui.QSpinBox(self.layoutWidget1)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.spinBox_no_layers.sizePolicy().hasHeightForWidth()
        )
        self.spinBox_no_layers.setSizePolicy(sizePolicy)
        self.spinBox_no_layers.setMinimumSize(QtCore.QSize(0, 30))
        self.spinBox_no_layers.setMaximumSize(QtCore.QSize(600, 30))
        self.spinBox_no_layers.setMinimum(1)
        self.spinBox_no_layers.setMaximum(500)
        self.spinBox_no_layers.setProperty("value", 30)
        self.spinBox_no_layers.setObjectName(_fromUtf8("spinBox_no_layers"))
        self.horizontalLayout_27.addWidget(self.spinBox_no_layers)
        spacerItem14 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_27.addItem(spacerItem14)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_27)
        self.horizontalLayout_26 = QtGui.QHBoxLayout()
        self.horizontalLayout_26.setObjectName(_fromUtf8("horizontalLayout_26"))
        self.label_19 = QtGui.QLabel(self.layoutWidget1)
        self.label_19.setObjectName(_fromUtf8("label_19"))
        self.horizontalLayout_26.addWidget(self.label_19)
        self.doubleSpinBox_model_depth = QtGui.QDoubleSpinBox(self.layoutWidget1)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_model_depth.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_model_depth.setSizePolicy(sizePolicy)
        self.doubleSpinBox_model_depth.setMinimumSize(QtCore.QSize(63, 30))
        self.doubleSpinBox_model_depth.setMaximumSize(QtCore.QSize(50, 30))
        self.doubleSpinBox_model_depth.setDecimals(1)
        self.doubleSpinBox_model_depth.setMinimum(0.1)
        self.doubleSpinBox_model_depth.setMaximum(9999.0)
        self.doubleSpinBox_model_depth.setSingleStep(5.0)
        self.doubleSpinBox_model_depth.setProperty("value", 100.0)
        self.doubleSpinBox_model_depth.setObjectName(
            _fromUtf8("doubleSpinBox_model_depth")
        )
        self.horizontalLayout_26.addWidget(self.doubleSpinBox_model_depth)
        self.label_33 = QtGui.QLabel(self.layoutWidget1)
        self.label_33.setObjectName(_fromUtf8("label_33"))
        self.horizontalLayout_26.addWidget(self.label_33)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_26)
        self.horizontalLayout_23 = QtGui.QHBoxLayout()
        self.horizontalLayout_23.setObjectName(_fromUtf8("horizontalLayout_23"))
        spacerItem15 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_23.addItem(spacerItem15)
        self.label_25 = QtGui.QLabel(self.layoutWidget1)
        self.label_25.setObjectName(_fromUtf8("label_25"))
        self.horizontalLayout_23.addWidget(self.label_25)
        self.spinBox_firstlayer = QtGui.QSpinBox(self.layoutWidget1)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.spinBox_firstlayer.sizePolicy().hasHeightForWidth()
        )
        self.spinBox_firstlayer.setSizePolicy(sizePolicy)
        self.spinBox_firstlayer.setMinimumSize(QtCore.QSize(0, 30))
        self.spinBox_firstlayer.setMaximumSize(QtCore.QSize(500, 30))
        self.spinBox_firstlayer.setMinimum(1)
        self.spinBox_firstlayer.setMaximum(99999)
        self.spinBox_firstlayer.setProperty("value", 100)
        self.spinBox_firstlayer.setObjectName(_fromUtf8("spinBox_firstlayer"))
        self.horizontalLayout_23.addWidget(self.spinBox_firstlayer)
        self.label_27 = QtGui.QLabel(self.layoutWidget1)
        self.label_27.setObjectName(_fromUtf8("label_27"))
        self.horizontalLayout_23.addWidget(self.label_27)
        self.horizontalLayout_30.addLayout(self.horizontalLayout_23)
        self.layoutWidget2 = QtGui.QWidget(occamgui2D)
        self.layoutWidget2.setGeometry(QtCore.QRect(32, 427, 849, 34))
        self.layoutWidget2.setObjectName(_fromUtf8("layoutWidget2"))
        self.horizontalLayout_19 = QtGui.QHBoxLayout(self.layoutWidget2)
        self.horizontalLayout_19.setMargin(0)
        self.horizontalLayout_19.setObjectName(_fromUtf8("horizontalLayout_19"))
        self.horizontalLayout_15 = QtGui.QHBoxLayout()
        self.horizontalLayout_15.setObjectName(_fromUtf8("horizontalLayout_15"))
        self.checkBox_max_no_frequencies = QtGui.QCheckBox(self.layoutWidget2)
        self.checkBox_max_no_frequencies.setObjectName(
            _fromUtf8("checkBox_max_no_frequencies")
        )
        self.horizontalLayout_15.addWidget(self.checkBox_max_no_frequencies)
        self.spinBox_max_no_frequencies = QtGui.QSpinBox(self.layoutWidget2)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.spinBox_max_no_frequencies.sizePolicy().hasHeightForWidth()
        )
        self.spinBox_max_no_frequencies.setSizePolicy(sizePolicy)
        self.spinBox_max_no_frequencies.setMinimumSize(QtCore.QSize(63, 30))
        self.spinBox_max_no_frequencies.setMaximumSize(QtCore.QSize(16777215, 30))
        self.spinBox_max_no_frequencies.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.spinBox_max_no_frequencies.setReadOnly(False)
        self.spinBox_max_no_frequencies.setMinimum(0)
        self.spinBox_max_no_frequencies.setMaximum(999)
        self.spinBox_max_no_frequencies.setObjectName(
            _fromUtf8("spinBox_max_no_frequencies")
        )
        self.horizontalLayout_15.addWidget(self.spinBox_max_no_frequencies)
        spacerItem16 = QtGui.QSpacerItem(
            17, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_15.addItem(spacerItem16)
        self.horizontalLayout_19.addLayout(self.horizontalLayout_15)
        spacerItem17 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_19.addItem(spacerItem17)
        self.horizontalLayout_16 = QtGui.QHBoxLayout()
        self.horizontalLayout_16.setObjectName(_fromUtf8("horizontalLayout_16"))
        self.checkBox_min_frequency = QtGui.QCheckBox(self.layoutWidget2)
        self.checkBox_min_frequency.setObjectName(_fromUtf8("checkBox_min_frequency"))
        self.horizontalLayout_16.addWidget(self.checkBox_min_frequency)
        self.doubleSpinBox_min_frequency = QtGui.QDoubleSpinBox(self.layoutWidget2)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_min_frequency.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_min_frequency.setSizePolicy(sizePolicy)
        self.doubleSpinBox_min_frequency.setMinimumSize(QtCore.QSize(63, 30))
        self.doubleSpinBox_min_frequency.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_min_frequency.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.doubleSpinBox_min_frequency.setReadOnly(False)
        self.doubleSpinBox_min_frequency.setDecimals(4)
        self.doubleSpinBox_min_frequency.setMinimum(0.0)
        self.doubleSpinBox_min_frequency.setMaximum(9999.0)
        self.doubleSpinBox_min_frequency.setSingleStep(1.0)
        self.doubleSpinBox_min_frequency.setObjectName(
            _fromUtf8("doubleSpinBox_min_frequency")
        )
        self.horizontalLayout_16.addWidget(self.doubleSpinBox_min_frequency)
        self.label_31 = QtGui.QLabel(self.layoutWidget2)
        self.label_31.setObjectName(_fromUtf8("label_31"))
        self.horizontalLayout_16.addWidget(self.label_31)
        self.horizontalLayout_19.addLayout(self.horizontalLayout_16)
        spacerItem18 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_19.addItem(spacerItem18)
        self.horizontalLayout_18 = QtGui.QHBoxLayout()
        self.horizontalLayout_18.setObjectName(_fromUtf8("horizontalLayout_18"))
        spacerItem19 = QtGui.QSpacerItem(
            17, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_18.addItem(spacerItem19)
        self.checkBox_max_frequency = QtGui.QCheckBox(self.layoutWidget2)
        self.checkBox_max_frequency.setObjectName(_fromUtf8("checkBox_max_frequency"))
        self.horizontalLayout_18.addWidget(self.checkBox_max_frequency)
        self.doubleSpinBox_max_frequency = QtGui.QDoubleSpinBox(self.layoutWidget2)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_max_frequency.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_max_frequency.setSizePolicy(sizePolicy)
        self.doubleSpinBox_max_frequency.setMinimumSize(QtCore.QSize(63, 30))
        self.doubleSpinBox_max_frequency.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_max_frequency.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.doubleSpinBox_max_frequency.setReadOnly(False)
        self.doubleSpinBox_max_frequency.setDecimals(4)
        self.doubleSpinBox_max_frequency.setMinimum(0.0)
        self.doubleSpinBox_max_frequency.setMaximum(9999.0)
        self.doubleSpinBox_max_frequency.setSingleStep(1.0)
        self.doubleSpinBox_max_frequency.setProperty("value", 1000.0)
        self.doubleSpinBox_max_frequency.setObjectName(
            _fromUtf8("doubleSpinBox_max_frequency")
        )
        self.horizontalLayout_18.addWidget(self.doubleSpinBox_max_frequency)
        self.label_32 = QtGui.QLabel(self.layoutWidget2)
        self.label_32.setObjectName(_fromUtf8("label_32"))
        self.horizontalLayout_18.addWidget(self.label_32)
        self.horizontalLayout_19.addLayout(self.horizontalLayout_18)
        self.layoutWidget3 = QtGui.QWidget(occamgui2D)
        self.layoutWidget3.setGeometry(QtCore.QRect(20, 609, 870, 34))
        self.layoutWidget3.setObjectName(_fromUtf8("layoutWidget3"))
        self.horizontalLayout_25 = QtGui.QHBoxLayout(self.layoutWidget3)
        self.horizontalLayout_25.setMargin(0)
        self.horizontalLayout_25.setObjectName(_fromUtf8("horizontalLayout_25"))
        self.horizontalLayout_20 = QtGui.QHBoxLayout()
        self.horizontalLayout_20.setObjectName(_fromUtf8("horizontalLayout_20"))
        self.checkBox_rho_error = QtGui.QCheckBox(self.layoutWidget3)
        self.checkBox_rho_error.setObjectName(_fromUtf8("checkBox_rho_error"))
        self.horizontalLayout_20.addWidget(self.checkBox_rho_error)
        self.doubleSpinBox_rho_error = QtGui.QDoubleSpinBox(self.layoutWidget3)
        self.doubleSpinBox_rho_error.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_rho_error.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_rho_error.setDecimals(0)
        self.doubleSpinBox_rho_error.setMinimum(0.0)
        self.doubleSpinBox_rho_error.setMaximum(100.0)
        self.doubleSpinBox_rho_error.setSingleStep(1.0)
        self.doubleSpinBox_rho_error.setProperty("value", 10.0)
        self.doubleSpinBox_rho_error.setObjectName(_fromUtf8("doubleSpinBox_rho_error"))
        self.horizontalLayout_20.addWidget(self.doubleSpinBox_rho_error)
        spacerItem20 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_20.addItem(spacerItem20)
        self.horizontalLayout_25.addLayout(self.horizontalLayout_20)
        self.horizontalLayout_22 = QtGui.QHBoxLayout()
        self.horizontalLayout_22.setObjectName(_fromUtf8("horizontalLayout_22"))
        spacerItem21 = QtGui.QSpacerItem(
            13, 27, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_22.addItem(spacerItem21)
        spacerItem22 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_22.addItem(spacerItem22)
        self.checkBox_phase_error = QtGui.QCheckBox(self.layoutWidget3)
        self.checkBox_phase_error.setObjectName(_fromUtf8("checkBox_phase_error"))
        self.horizontalLayout_22.addWidget(self.checkBox_phase_error)
        self.doubleSpinBox_phase_error = QtGui.QDoubleSpinBox(self.layoutWidget3)
        self.doubleSpinBox_phase_error.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_phase_error.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_phase_error.setDecimals(0)
        self.doubleSpinBox_phase_error.setMinimum(0.0)
        self.doubleSpinBox_phase_error.setMaximum(100.0)
        self.doubleSpinBox_phase_error.setSingleStep(1.0)
        self.doubleSpinBox_phase_error.setProperty("value", 15.0)
        self.doubleSpinBox_phase_error.setObjectName(
            _fromUtf8("doubleSpinBox_phase_error")
        )
        self.horizontalLayout_22.addWidget(self.doubleSpinBox_phase_error)
        spacerItem23 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_22.addItem(spacerItem23)
        self.horizontalLayout_25.addLayout(self.horizontalLayout_22)
        spacerItem24 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_25.addItem(spacerItem24)
        self.horizontalLayout_24 = QtGui.QHBoxLayout()
        self.horizontalLayout_24.setObjectName(_fromUtf8("horizontalLayout_24"))
        spacerItem25 = QtGui.QSpacerItem(
            13, 27, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_24.addItem(spacerItem25)
        spacerItem26 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_24.addItem(spacerItem26)
        self.checkBox_tipper_error = QtGui.QCheckBox(self.layoutWidget3)
        self.checkBox_tipper_error.setObjectName(_fromUtf8("checkBox_tipper_error"))
        self.horizontalLayout_24.addWidget(self.checkBox_tipper_error)
        self.doubleSpinBox_tipper_error = QtGui.QDoubleSpinBox(self.layoutWidget3)
        self.doubleSpinBox_tipper_error.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_tipper_error.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_tipper_error.setDecimals(0)
        self.doubleSpinBox_tipper_error.setMinimum(0.0)
        self.doubleSpinBox_tipper_error.setMaximum(100.0)
        self.doubleSpinBox_tipper_error.setSingleStep(1.0)
        self.doubleSpinBox_tipper_error.setProperty("value", 20.0)
        self.doubleSpinBox_tipper_error.setObjectName(
            _fromUtf8("doubleSpinBox_tipper_error")
        )
        self.horizontalLayout_24.addWidget(self.doubleSpinBox_tipper_error)
        spacerItem27 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_24.addItem(spacerItem27)
        self.horizontalLayout_25.addLayout(self.horizontalLayout_24)
        self.layoutWidget4 = QtGui.QWidget(occamgui2D)
        self.layoutWidget4.setGeometry(QtCore.QRect(24, 569, 856, 34))
        self.layoutWidget4.setObjectName(_fromUtf8("layoutWidget4"))
        self.horizontalLayout_36 = QtGui.QHBoxLayout(self.layoutWidget4)
        self.horizontalLayout_36.setMargin(0)
        self.horizontalLayout_36.setObjectName(_fromUtf8("horizontalLayout_36"))
        self.horizontalLayout_31 = QtGui.QHBoxLayout()
        self.horizontalLayout_31.setObjectName(_fromUtf8("horizontalLayout_31"))
        self.label_17 = QtGui.QLabel(self.layoutWidget4)
        self.label_17.setObjectName(_fromUtf8("label_17"))
        self.horizontalLayout_31.addWidget(self.label_17)
        self.doubleSpinBox_rms = QtGui.QDoubleSpinBox(self.layoutWidget4)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_rms.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_rms.setSizePolicy(sizePolicy)
        self.doubleSpinBox_rms.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_rms.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_rms.setMinimum(0.01)
        self.doubleSpinBox_rms.setMaximum(1000.0)
        self.doubleSpinBox_rms.setSingleStep(0.01)
        self.doubleSpinBox_rms.setProperty("value", 1.0)
        self.doubleSpinBox_rms.setObjectName(_fromUtf8("doubleSpinBox_rms"))
        self.horizontalLayout_31.addWidget(self.doubleSpinBox_rms)
        spacerItem28 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_31.addItem(spacerItem28)
        self.horizontalLayout_36.addLayout(self.horizontalLayout_31)
        self.horizontalLayout_35 = QtGui.QHBoxLayout()
        self.horizontalLayout_35.setObjectName(_fromUtf8("horizontalLayout_35"))
        spacerItem29 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_35.addItem(spacerItem29)
        self.label_16 = QtGui.QLabel(self.layoutWidget4)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_16.sizePolicy().hasHeightForWidth())
        self.label_16.setSizePolicy(sizePolicy)
        self.label_16.setObjectName(_fromUtf8("label_16"))
        self.horizontalLayout_35.addWidget(self.label_16)
        self.spinBox_max_no_iterations = QtGui.QSpinBox(self.layoutWidget4)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.spinBox_max_no_iterations.sizePolicy().hasHeightForWidth()
        )
        self.spinBox_max_no_iterations.setSizePolicy(sizePolicy)
        self.spinBox_max_no_iterations.setMinimumSize(QtCore.QSize(0, 30))
        self.spinBox_max_no_iterations.setMaximumSize(QtCore.QSize(50, 30))
        self.spinBox_max_no_iterations.setMinimum(1)
        self.spinBox_max_no_iterations.setMaximum(999)
        self.spinBox_max_no_iterations.setProperty("value", 30)
        self.spinBox_max_no_iterations.setObjectName(
            _fromUtf8("spinBox_max_no_iterations")
        )
        self.horizontalLayout_35.addWidget(self.spinBox_max_no_iterations)
        spacerItem30 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_35.addItem(spacerItem30)
        self.horizontalLayout_36.addLayout(self.horizontalLayout_35)
        self.horizontalLayout_32 = QtGui.QHBoxLayout()
        self.horizontalLayout_32.setObjectName(_fromUtf8("horizontalLayout_32"))
        spacerItem31 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_32.addItem(spacerItem31)
        self.label_24 = QtGui.QLabel(self.layoutWidget4)
        self.label_24.setObjectName(_fromUtf8("label_24"))
        self.horizontalLayout_32.addWidget(self.label_24)
        self.doubleSpinBox_rhostart = QtGui.QDoubleSpinBox(self.layoutWidget4)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_rhostart.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_rhostart.setSizePolicy(sizePolicy)
        self.doubleSpinBox_rhostart.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_rhostart.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_rhostart.setDecimals(1)
        self.doubleSpinBox_rhostart.setMaximum(100000.0)
        self.doubleSpinBox_rhostart.setProperty("value", 100.0)
        self.doubleSpinBox_rhostart.setObjectName(_fromUtf8("doubleSpinBox_rhostart"))
        self.horizontalLayout_32.addWidget(self.doubleSpinBox_rhostart)
        self.label_30 = QtGui.QLabel(self.layoutWidget4)
        self.label_30.setObjectName(_fromUtf8("label_30"))
        self.horizontalLayout_32.addWidget(self.label_30)
        self.horizontalLayout_36.addLayout(self.horizontalLayout_32)
        self.layoutWidget5 = QtGui.QWidget(occamgui2D)
        self.layoutWidget5.setGeometry(QtCore.QRect(140, 471, 628, 34))
        self.layoutWidget5.setObjectName(_fromUtf8("layoutWidget5"))
        self.horizontalLayout_37 = QtGui.QHBoxLayout(self.layoutWidget5)
        self.horizontalLayout_37.setMargin(0)
        self.horizontalLayout_37.setObjectName(_fromUtf8("horizontalLayout_37"))
        self.horizontalLayout_33 = QtGui.QHBoxLayout()
        self.horizontalLayout_33.setObjectName(_fromUtf8("horizontalLayout_33"))
        self.label_21 = QtGui.QLabel(self.layoutWidget5)
        self.label_21.setObjectName(_fromUtf8("label_21"))
        self.horizontalLayout_33.addWidget(self.label_21)
        self.doubleSpinBox_mergethreshold = QtGui.QDoubleSpinBox(self.layoutWidget5)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_mergethreshold.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_mergethreshold.setSizePolicy(sizePolicy)
        self.doubleSpinBox_mergethreshold.setMinimumSize(QtCore.QSize(0, 30))
        self.doubleSpinBox_mergethreshold.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_mergethreshold.setMinimum(0.05)
        self.doubleSpinBox_mergethreshold.setMaximum(2.0)
        self.doubleSpinBox_mergethreshold.setSingleStep(0.05)
        self.doubleSpinBox_mergethreshold.setProperty("value", 0.75)
        self.doubleSpinBox_mergethreshold.setObjectName(
            _fromUtf8("doubleSpinBox_mergethreshold")
        )
        self.horizontalLayout_33.addWidget(self.doubleSpinBox_mergethreshold)
        spacerItem32 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_33.addItem(spacerItem32)
        self.horizontalLayout_37.addLayout(self.horizontalLayout_33)
        spacerItem33 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_37.addItem(spacerItem33)
        self.horizontalLayout_34 = QtGui.QHBoxLayout()
        self.horizontalLayout_34.setObjectName(_fromUtf8("horizontalLayout_34"))
        spacerItem34 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_34.addItem(spacerItem34)
        self.label_26 = QtGui.QLabel(self.layoutWidget5)
        self.label_26.setObjectName(_fromUtf8("label_26"))
        self.horizontalLayout_34.addWidget(self.label_26)
        self.spinBox_maxblockwidth = QtGui.QSpinBox(self.layoutWidget5)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.spinBox_maxblockwidth.sizePolicy().hasHeightForWidth()
        )
        self.spinBox_maxblockwidth.setSizePolicy(sizePolicy)
        self.spinBox_maxblockwidth.setMinimumSize(QtCore.QSize(0, 30))
        self.spinBox_maxblockwidth.setMaximumSize(QtCore.QSize(500, 30))
        self.spinBox_maxblockwidth.setMinimum(5)
        self.spinBox_maxblockwidth.setMaximum(100000)
        self.spinBox_maxblockwidth.setSingleStep(100)
        self.spinBox_maxblockwidth.setProperty("value", 1000)
        self.spinBox_maxblockwidth.setObjectName(_fromUtf8("spinBox_maxblockwidth"))
        self.horizontalLayout_34.addWidget(self.spinBox_maxblockwidth)
        self.label_22 = QtGui.QLabel(self.layoutWidget5)
        self.label_22.setObjectName(_fromUtf8("label_22"))
        self.horizontalLayout_34.addWidget(self.label_22)
        spacerItem35 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_34.addItem(spacerItem35)
        self.horizontalLayout_37.addLayout(self.horizontalLayout_34)
        self.line = QtGui.QFrame(occamgui2D)
        self.line.setGeometry(QtCore.QRect(175, 39, 707, 20))
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.layoutWidget6 = QtGui.QWidget(occamgui2D)
        self.layoutWidget6.setGeometry(QtCore.QRect(190, 10, 679, 32))
        self.layoutWidget6.setObjectName(_fromUtf8("layoutWidget6"))
        self.horizontalLayout_39 = QtGui.QHBoxLayout(self.layoutWidget6)
        self.horizontalLayout_39.setMargin(0)
        self.horizontalLayout_39.setObjectName(_fromUtf8("horizontalLayout_39"))
        self.label = QtGui.QLabel(self.layoutWidget6)
        self.label.setObjectName(_fromUtf8("label"))
        self.horizontalLayout_39.addWidget(self.label)
        spacerItem36 = QtGui.QSpacerItem(
            18, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_39.addItem(spacerItem36)
        self.button_browse_configfile = QtGui.QToolButton(self.layoutWidget6)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.button_browse_configfile.sizePolicy().hasHeightForWidth()
        )
        self.button_browse_configfile.setSizePolicy(sizePolicy)
        self.button_browse_configfile.setMinimumSize(QtCore.QSize(30, 25))
        self.button_browse_configfile.setMaximumSize(QtCore.QSize(30, 25))
        self.button_browse_configfile.setObjectName(
            _fromUtf8("button_browse_configfile")
        )
        self.horizontalLayout_39.addWidget(self.button_browse_configfile)
        self.lineEdit_browse_configfile = QtGui.QLineEdit(self.layoutWidget6)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_browse_configfile.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_browse_configfile.setSizePolicy(sizePolicy)
        self.lineEdit_browse_configfile.setMinimumSize(QtCore.QSize(300, 30))
        self.lineEdit_browse_configfile.setMaximumSize(QtCore.QSize(500, 30))
        self.lineEdit_browse_configfile.setObjectName(
            _fromUtf8("lineEdit_browse_configfile")
        )
        self.horizontalLayout_39.addWidget(self.lineEdit_browse_configfile)
        spacerItem37 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_39.addItem(spacerItem37)
        self.button_load_configfile = QtGui.QToolButton(self.layoutWidget6)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.button_load_configfile.sizePolicy().hasHeightForWidth()
        )
        self.button_load_configfile.setSizePolicy(sizePolicy)
        self.button_load_configfile.setObjectName(_fromUtf8("button_load_configfile"))
        self.horizontalLayout_39.addWidget(self.button_load_configfile)
        self.layoutWidget7 = QtGui.QWidget(occamgui2D)
        self.layoutWidget7.setGeometry(QtCore.QRect(62, 341, 818, 33))
        self.layoutWidget7.setObjectName(_fromUtf8("layoutWidget7"))
        self.horizontalLayout_44 = QtGui.QHBoxLayout(self.layoutWidget7)
        self.horizontalLayout_44.setMargin(0)
        self.horizontalLayout_44.setObjectName(_fromUtf8("horizontalLayout_44"))
        self.horizontalLayout_43 = QtGui.QHBoxLayout()
        self.horizontalLayout_43.setObjectName(_fromUtf8("horizontalLayout_43"))
        self.label_iterationstep = QtGui.QLabel(self.layoutWidget7)
        self.label_iterationstep.setObjectName(_fromUtf8("label_iterationstep"))
        self.horizontalLayout_43.addWidget(self.label_iterationstep)
        self.spinBox_iterationstep = QtGui.QSpinBox(self.layoutWidget7)
        self.spinBox_iterationstep.setMaximum(500)
        self.spinBox_iterationstep.setObjectName(_fromUtf8("spinBox_iterationstep"))
        self.horizontalLayout_43.addWidget(self.spinBox_iterationstep)
        spacerItem38 = QtGui.QSpacerItem(
            17, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_43.addItem(spacerItem38)
        self.horizontalLayout_44.addLayout(self.horizontalLayout_43)
        self.horizontalLayout_42 = QtGui.QHBoxLayout()
        self.horizontalLayout_42.setObjectName(_fromUtf8("horizontalLayout_42"))
        spacerItem39 = QtGui.QSpacerItem(
            13, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_42.addItem(spacerItem39)
        self._label_lagrange = QtGui.QLabel(self.layoutWidget7)
        self._label_lagrange.setObjectName(_fromUtf8("_label_lagrange"))
        self.horizontalLayout_42.addWidget(self._label_lagrange)
        self.doubleSpinBox_lagrange = QtGui.QDoubleSpinBox(self.layoutWidget7)
        self.doubleSpinBox_lagrange.setMaximum(9999.0)
        self.doubleSpinBox_lagrange.setSingleStep(0.5)
        self.doubleSpinBox_lagrange.setProperty("value", 5.0)
        self.doubleSpinBox_lagrange.setObjectName(_fromUtf8("doubleSpinBox_lagrange"))
        self.horizontalLayout_42.addWidget(self.doubleSpinBox_lagrange)
        self.horizontalLayout_44.addLayout(self.horizontalLayout_42)
        self.horizontalLayout_41 = QtGui.QHBoxLayout()
        self.horizontalLayout_41.setObjectName(_fromUtf8("horizontalLayout_41"))
        spacerItem40 = QtGui.QSpacerItem(
            17, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_41.addItem(spacerItem40)
        self.checkBox_misfitreached = QtGui.QCheckBox(self.layoutWidget7)
        self.checkBox_misfitreached.setObjectName(_fromUtf8("checkBox_misfitreached"))
        self.horizontalLayout_41.addWidget(self.checkBox_misfitreached)
        self.horizontalLayout_44.addLayout(self.horizontalLayout_41)
        self.horizontalLayout_40 = QtGui.QHBoxLayout()
        self.horizontalLayout_40.setObjectName(_fromUtf8("horizontalLayout_40"))
        spacerItem41 = QtGui.QSpacerItem(
            17, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_40.addItem(spacerItem41)
        self.label_6 = QtGui.QLabel(self.layoutWidget7)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.horizontalLayout_40.addWidget(self.label_6)
        self.comboBox_debuglevel = QtGui.QComboBox(self.layoutWidget7)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.comboBox_debuglevel.sizePolicy().hasHeightForWidth()
        )
        self.comboBox_debuglevel.setSizePolicy(sizePolicy)
        self.comboBox_debuglevel.setMaxVisibleItems(3)
        self.comboBox_debuglevel.setIconSize(QtCore.QSize(10, 16))
        self.comboBox_debuglevel.setObjectName(_fromUtf8("comboBox_debuglevel"))
        self.comboBox_debuglevel.addItem(_fromUtf8(""))
        self.comboBox_debuglevel.addItem(_fromUtf8(""))
        self.comboBox_debuglevel.addItem(_fromUtf8(""))
        self.horizontalLayout_40.addWidget(self.comboBox_debuglevel)
        self.horizontalLayout_44.addLayout(self.horizontalLayout_40)
        self.layoutWidget8 = QtGui.QWidget(occamgui2D)
        self.layoutWidget8.setGeometry(QtCore.QRect(317, 385, 189, 27))
        self.layoutWidget8.setObjectName(_fromUtf8("layoutWidget8"))
        self.horizontalLayout_17 = QtGui.QHBoxLayout(self.layoutWidget8)
        self.horizontalLayout_17.setMargin(0)
        self.horizontalLayout_17.setObjectName(_fromUtf8("horizontalLayout_17"))
        spacerItem42 = QtGui.QSpacerItem(
            40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_17.addItem(spacerItem42)
        self.label_8 = QtGui.QLabel(self.layoutWidget8)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_8.sizePolicy().hasHeightForWidth())
        self.label_8.setSizePolicy(sizePolicy)
        self.label_8.setMaximumSize(QtCore.QSize(50, 16777215))
        self.label_8.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.horizontalLayout_17.addWidget(self.label_8)
        self.comboBox_mode = QtGui.QComboBox(self.layoutWidget8)
        self.comboBox_mode.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.comboBox_mode.sizePolicy().hasHeightForWidth()
        )
        self.comboBox_mode.setSizePolicy(sizePolicy)
        self.comboBox_mode.setMinimumSize(QtCore.QSize(90, 25))
        self.comboBox_mode.setMaximumSize(QtCore.QSize(16777215, 30))
        self.comboBox_mode.setMaxVisibleItems(5)
        self.comboBox_mode.setFrame(True)
        self.comboBox_mode.setObjectName(_fromUtf8("comboBox_mode"))
        self.comboBox_mode.addItem(_fromUtf8(""))
        self.comboBox_mode.addItem(_fromUtf8(""))
        self.comboBox_mode.addItem(_fromUtf8(""))
        self.comboBox_mode.addItem(_fromUtf8(""))
        self.comboBox_mode.addItem(_fromUtf8(""))
        self.horizontalLayout_17.addWidget(self.comboBox_mode)
        self.layoutWidget9 = QtGui.QWidget(occamgui2D)
        self.layoutWidget9.setGeometry(QtCore.QRect(515, 381, 186, 32))
        self.layoutWidget9.setObjectName(_fromUtf8("layoutWidget9"))
        self.horizontalLayout_28 = QtGui.QHBoxLayout(self.layoutWidget9)
        self.horizontalLayout_28.setMargin(0)
        self.horizontalLayout_28.setObjectName(_fromUtf8("horizontalLayout_28"))
        spacerItem43 = QtGui.QSpacerItem(
            13, 27, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_28.addItem(spacerItem43)
        self.checkBox_strike = QtGui.QCheckBox(self.layoutWidget9)
        self.checkBox_strike.setChecked(True)
        self.checkBox_strike.setObjectName(_fromUtf8("checkBox_strike"))
        self.horizontalLayout_28.addWidget(self.checkBox_strike)
        self.doubleSpinBox_strike = QtGui.QDoubleSpinBox(self.layoutWidget9)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.doubleSpinBox_strike.sizePolicy().hasHeightForWidth()
        )
        self.doubleSpinBox_strike.setSizePolicy(sizePolicy)
        self.doubleSpinBox_strike.setMinimumSize(QtCore.QSize(63, 30))
        self.doubleSpinBox_strike.setMaximumSize(QtCore.QSize(16777215, 30))
        self.doubleSpinBox_strike.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.doubleSpinBox_strike.setDecimals(1)
        self.doubleSpinBox_strike.setMinimum(0.0)
        self.doubleSpinBox_strike.setMaximum(180.0)
        self.doubleSpinBox_strike.setSingleStep(1.0)
        self.doubleSpinBox_strike.setObjectName(_fromUtf8("doubleSpinBox_strike"))
        self.horizontalLayout_28.addWidget(self.doubleSpinBox_strike)
        self.layoutWidget10 = QtGui.QWidget(occamgui2D)
        self.layoutWidget10.setGeometry(QtCore.QRect(712, 384, 164, 29))
        self.layoutWidget10.setObjectName(_fromUtf8("layoutWidget10"))
        self.horizontalLayout_29 = QtGui.QHBoxLayout(self.layoutWidget10)
        self.horizontalLayout_29.setMargin(0)
        self.horizontalLayout_29.setObjectName(_fromUtf8("horizontalLayout_29"))
        spacerItem44 = QtGui.QSpacerItem(
            13, 27, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum
        )
        self.horizontalLayout_29.addItem(spacerItem44)
        self.label_edi_type = QtGui.QLabel(self.layoutWidget10)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Preferred
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.label_edi_type.sizePolicy().hasHeightForWidth()
        )
        self.label_edi_type.setSizePolicy(sizePolicy)
        self.label_edi_type.setMaximumSize(QtCore.QSize(50, 16777215))
        self.label_edi_type.setAlignment(
            QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )
        self.label_edi_type.setObjectName(_fromUtf8("label_edi_type"))
        self.horizontalLayout_29.addWidget(self.label_edi_type)
        self.comboBox_edi_type = QtGui.QComboBox(self.layoutWidget10)
        self.comboBox_edi_type.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.comboBox_edi_type.sizePolicy().hasHeightForWidth()
        )
        self.comboBox_edi_type.setSizePolicy(sizePolicy)
        self.comboBox_edi_type.setMinimumSize(QtCore.QSize(90, 25))
        self.comboBox_edi_type.setMaximumSize(QtCore.QSize(16777215, 30))
        self.comboBox_edi_type.setMaxVisibleItems(3)
        self.comboBox_edi_type.setFrame(True)
        self.comboBox_edi_type.setObjectName(_fromUtf8("comboBox_edi_type"))
        self.comboBox_edi_type.addItem(_fromUtf8(""))
        self.comboBox_edi_type.addItem(_fromUtf8(""))
        self.comboBox_edi_type.addItem(_fromUtf8(""))
        self.horizontalLayout_29.addWidget(self.comboBox_edi_type)
        self.layoutWidget11 = QtGui.QWidget(occamgui2D)
        self.layoutWidget11.setGeometry(QtCore.QRect(59, 385, 252, 26))
        self.layoutWidget11.setObjectName(_fromUtf8("layoutWidget11"))
        self.horizontalLayout_21 = QtGui.QHBoxLayout(self.layoutWidget11)
        self.horizontalLayout_21.setMargin(0)
        self.horizontalLayout_21.setObjectName(_fromUtf8("horizontalLayout_21"))
        self.label_7 = QtGui.QLabel(self.layoutWidget11)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.horizontalLayout_21.addWidget(self.label_7)
        self.lineEdit_modelname = QtGui.QLineEdit(self.layoutWidget11)
        sizePolicy = QtGui.QSizePolicy(
            QtGui.QSizePolicy.MinimumExpanding, QtGui.QSizePolicy.Fixed
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(
            self.lineEdit_modelname.sizePolicy().hasHeightForWidth()
        )
        self.lineEdit_modelname.setSizePolicy(sizePolicy)
        self.lineEdit_modelname.setObjectName(_fromUtf8("lineEdit_modelname"))
        self.horizontalLayout_21.addWidget(self.lineEdit_modelname)
        self.label_18.setBuddy(self.spinBox_no_layers)
        self.label_19.setBuddy(self.doubleSpinBox_model_depth)
        self.label_17.setBuddy(self.doubleSpinBox_rms)
        self.label_16.setBuddy(self.spinBox_max_no_iterations)
        self.label_24.setBuddy(self.doubleSpinBox_rhostart)
        self.label_21.setBuddy(self.doubleSpinBox_mergethreshold)
        self.label_26.setBuddy(self.doubleSpinBox_rms)

        self.retranslateUi(occamgui2D)
        self.comboBox_debuglevel.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(occamgui2D)

    def retranslateUi(self, occamgui2D):
        occamgui2D.setWindowTitle(
            QtGui.QApplication.translate(
                "occamgui2D", "OCCAM 2D - setup", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_29.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Occam2D", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_4.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Parameters", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.pushButton_generateinputfile.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Generate input files",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.pushButton_quit.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Exit", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.pushButton_checkparameter.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Check parameters", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.pushButton_runoccam.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Run OCCAM 2D", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_2.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "working directory", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.button_browse_wd.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "...", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.lineEdit_browse_wd.setText(
            QtGui.QApplication.translate(
                "occamgui2D", ".", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_5.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "OCCAM executable", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.button_browse_occam.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "...", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.lineEdit_browse_occam.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Occam2D", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_28.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "EDI files folder", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.button_browse_edis.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "...", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.lineEdit_browse_edi.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "edi", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_usestationlist.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "use station list file",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.pushButton_loadstations.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Load station list", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_usedatafile.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "use existing data file",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.pushButton_loaddatafile.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Load data file", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_3.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Choose a data file name",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.lineEdit_datafilename.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "OccamInputData.dat", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_useiterationfile.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "use old iteration as startup file",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.pushButton_loaditerationfile.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Load iteration file",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_18.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Number of model layers",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_19.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Model depth", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_33.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "km", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_25.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Thickness of first layer",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_27.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "m", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_max_no_frequencies.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Define max # frequencies",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.checkBox_min_frequency.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Define min frequency",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_31.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Hz", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_max_frequency.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Define max frequency",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_32.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Hz", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_rho_error.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Set Resistivity minimum error (%):",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.checkBox_phase_error.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Set phase minimum error (%):",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.checkBox_tipper_error.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Set Tipper minimum error (%):",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_17.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Target RMS", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_16.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Maximum # iterations",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_24.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Resistivity (homogeneous half space -- starting model)",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_30.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Ohm m", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_21.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "width/thickness ratio for amalgamating blocks",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_26.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Maximum block width",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_22.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "m", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Load old configuration from file:",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.button_browse_configfile.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "...", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.button_load_configfile.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Apply old configuration",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_iterationstep.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Current iteration step",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self._label_lagrange.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Lagrange parameter", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.checkBox_misfitreached.setText(
            QtGui.QApplication.translate(
                "occamgui2D",
                "Misfit reached (smoothing only)",
                None,
                QtGui.QApplication.UnicodeUTF8,
            )
        )
        self.label_6.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Debug level", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.comboBox_debuglevel.setItemText(
            0,
            QtGui.QApplication.translate(
                "occamgui2D", "0", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_debuglevel.setItemText(
            1,
            QtGui.QApplication.translate(
                "occamgui2D", "1", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_debuglevel.setItemText(
            2,
            QtGui.QApplication.translate(
                "occamgui2D", "2", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.label_8.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Mode(s)", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.comboBox_mode.setItemText(
            0,
            QtGui.QApplication.translate(
                "occamgui2D", "TM + TE", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_mode.setItemText(
            1,
            QtGui.QApplication.translate(
                "occamgui2D", "TM", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_mode.setItemText(
            2,
            QtGui.QApplication.translate(
                "occamgui2D", "TE", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_mode.setItemText(
            3,
            QtGui.QApplication.translate(
                "occamgui2D", "Tipper", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_mode.setItemText(
            4,
            QtGui.QApplication.translate(
                "occamgui2D", "All", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.checkBox_strike.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Strike angle", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.label_edi_type.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "EDI type", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.comboBox_edi_type.setItemText(
            0,
            QtGui.QApplication.translate(
                "occamgui2D", "Z", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_edi_type.setItemText(
            1,
            QtGui.QApplication.translate(
                "occamgui2D", "Rho/Phase", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.comboBox_edi_type.setItemText(
            2,
            QtGui.QApplication.translate(
                "occamgui2D", "Spectra", None, QtGui.QApplication.UnicodeUTF8
            ),
        )
        self.label_7.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "Model name", None, QtGui.QApplication.UnicodeUTF8
            )
        )
        self.lineEdit_modelname.setText(
            QtGui.QApplication.translate(
                "occamgui2D", "test1", None, QtGui.QApplication.UnicodeUTF8
            )
        )
