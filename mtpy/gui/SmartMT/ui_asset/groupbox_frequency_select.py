# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_frequency_select.ui'
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

class Ui_GroupBox_frequency_select(object):
    def setupUi(self, GroupBox_frequency_select):
        GroupBox_frequency_select.setObjectName(_fromUtf8("GroupBox_frequency_select"))
        GroupBox_frequency_select.resize(300, 300)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_frequency_select.sizePolicy().hasHeightForWidth())
        GroupBox_frequency_select.setSizePolicy(sizePolicy)
        self.gridLayout = QtGui.QGridLayout(GroupBox_frequency_select)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.checkBox_existing_only = QtGui.QCheckBox(GroupBox_frequency_select)
        self.checkBox_existing_only.setObjectName(_fromUtf8("checkBox_existing_only"))
        self.gridLayout.addWidget(self.checkBox_existing_only, 7, 0, 1, 2)
        self.radioButton_period = QtGui.QRadioButton(GroupBox_frequency_select)
        self.radioButton_period.setObjectName(_fromUtf8("radioButton_period"))
        self.gridLayout.addWidget(self.radioButton_period, 2, 1, 1, 1)
        self.radioButton_frequency = QtGui.QRadioButton(GroupBox_frequency_select)
        self.radioButton_frequency.setChecked(True)
        self.radioButton_frequency.setObjectName(_fromUtf8("radioButton_frequency"))
        self.gridLayout.addWidget(self.radioButton_frequency, 2, 0, 1, 1)
        self.listView_selected = QtGui.QListView(GroupBox_frequency_select)
        self.listView_selected.setObjectName(_fromUtf8("listView_selected"))
        self.gridLayout.addWidget(self.listView_selected, 8, 0, 1, 2)
        self.widget_histgram = QtGui.QWidget(GroupBox_frequency_select)
        self.widget_histgram.setObjectName(_fromUtf8("widget_histgram"))
        self.verticalLayout = QtGui.QVBoxLayout(self.widget_histgram)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(self.widget_histgram)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.gridLayout.addWidget(self.widget_histgram, 0, 0, 1, 2)

        self.retranslateUi(GroupBox_frequency_select)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_frequency_select)

    def retranslateUi(self, GroupBox_frequency_select):
        GroupBox_frequency_select.setWindowTitle(_translate("GroupBox_frequency_select", "GroupBox", None))
        GroupBox_frequency_select.setTitle(_translate("GroupBox_frequency_select", "Frequency/Period", None))
        self.checkBox_existing_only.setToolTip(_translate("GroupBox_frequency_select", "<html><head/><body><p>Check this to select the frequency/period that is exists in the data file. The frequency/period that is closest to the cursor will be selected.</p></body></html>", None))
        self.checkBox_existing_only.setText(_translate("GroupBox_frequency_select", "Selecting Existing Frequency/Period Only", None))
        self.radioButton_period.setText(_translate("GroupBox_frequency_select", "Show Period", None))
        self.radioButton_frequency.setText(_translate("GroupBox_frequency_select", "Show Frequency", None))
        self.label.setText(_translate("GroupBox_frequency_select", "Place Holder for Frequency/Period Diagram", None))

