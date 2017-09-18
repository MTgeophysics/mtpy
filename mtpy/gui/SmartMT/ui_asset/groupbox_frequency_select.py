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
        GroupBox_frequency_select.resize(300, 270)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_frequency_select.sizePolicy().hasHeightForWidth())
        GroupBox_frequency_select.setSizePolicy(sizePolicy)
        self.gridLayout = QtGui.QGridLayout(GroupBox_frequency_select)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.widget_histgram = QtGui.QWidget(GroupBox_frequency_select)
        self.widget_histgram.setObjectName(_fromUtf8("widget_histgram"))
        self.verticalLayout = QtGui.QVBoxLayout(self.widget_histgram)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label_place_holder = QtGui.QLabel(self.widget_histgram)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_place_holder.sizePolicy().hasHeightForWidth())
        self.label_place_holder.setSizePolicy(sizePolicy)
        self.label_place_holder.setObjectName(_fromUtf8("label_place_holder"))
        self.verticalLayout.addWidget(self.label_place_holder)
        self.gridLayout.addWidget(self.widget_histgram, 0, 0, 1, 4)
        self.checkBox_existing_only = QtGui.QCheckBox(GroupBox_frequency_select)
        self.checkBox_existing_only.setObjectName(_fromUtf8("checkBox_existing_only"))
        self.gridLayout.addWidget(self.checkBox_existing_only, 10, 0, 1, 4)
        self.radioButton_frequency = QtGui.QRadioButton(GroupBox_frequency_select)
        self.radioButton_frequency.setChecked(True)
        self.radioButton_frequency.setObjectName(_fromUtf8("radioButton_frequency"))
        self.gridLayout.addWidget(self.radioButton_frequency, 2, 0, 1, 2)
        self.checkBox_show_existing = QtGui.QCheckBox(GroupBox_frequency_select)
        self.checkBox_show_existing.setObjectName(_fromUtf8("checkBox_show_existing"))
        self.gridLayout.addWidget(self.checkBox_show_existing, 9, 0, 1, 4)
        self.radioButton_period = QtGui.QRadioButton(GroupBox_frequency_select)
        self.radioButton_period.setObjectName(_fromUtf8("radioButton_period"))
        self.gridLayout.addWidget(self.radioButton_period, 2, 2, 1, 2)
        self.checkBox_y_log_scale = QtGui.QCheckBox(GroupBox_frequency_select)
        self.checkBox_y_log_scale.setObjectName(_fromUtf8("checkBox_y_log_scale"))
        self.gridLayout.addWidget(self.checkBox_y_log_scale, 8, 0, 1, 3)
        self.pushButton_delete = QtGui.QPushButton(GroupBox_frequency_select)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_delete.sizePolicy().hasHeightForWidth())
        self.pushButton_delete.setSizePolicy(sizePolicy)
        self.pushButton_delete.setObjectName(_fromUtf8("pushButton_delete"))
        self.gridLayout.addWidget(self.pushButton_delete, 11, 3, 1, 1)
        self.pushButton_clear = QtGui.QPushButton(GroupBox_frequency_select)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_clear.sizePolicy().hasHeightForWidth())
        self.pushButton_clear.setSizePolicy(sizePolicy)
        self.pushButton_clear.setObjectName(_fromUtf8("pushButton_clear"))
        self.gridLayout.addWidget(self.pushButton_clear, 12, 3, 1, 1)
        self.listView_selected = QtGui.QListView(GroupBox_frequency_select)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.listView_selected.sizePolicy().hasHeightForWidth())
        self.listView_selected.setSizePolicy(sizePolicy)
        self.listView_selected.setMaximumSize(QtCore.QSize(16777215, 100))
        self.listView_selected.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.listView_selected.setSelectionRectVisible(True)
        self.listView_selected.setObjectName(_fromUtf8("listView_selected"))
        self.gridLayout.addWidget(self.listView_selected, 11, 0, 3, 3)
        self.label_tolerance = QtGui.QLabel(GroupBox_frequency_select)
        self.label_tolerance.setObjectName(_fromUtf8("label_tolerance"))
        self.gridLayout.addWidget(self.label_tolerance, 1, 0, 1, 2)
        self.doubleSpinBox_tolerance = QtGui.QDoubleSpinBox(GroupBox_frequency_select)
        self.doubleSpinBox_tolerance.setMaximum(100.0)
        self.doubleSpinBox_tolerance.setProperty("value", 10.0)
        self.doubleSpinBox_tolerance.setObjectName(_fromUtf8("doubleSpinBox_tolerance"))
        self.gridLayout.addWidget(self.doubleSpinBox_tolerance, 1, 2, 1, 1)
        self.checkBox_x_log_scale = QtGui.QCheckBox(GroupBox_frequency_select)
        self.checkBox_x_log_scale.setObjectName(_fromUtf8("checkBox_x_log_scale"))
        self.gridLayout.addWidget(self.checkBox_x_log_scale, 3, 0, 1, 3)

        self.retranslateUi(GroupBox_frequency_select)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_frequency_select)

    def retranslateUi(self, GroupBox_frequency_select):
        GroupBox_frequency_select.setWindowTitle(_translate("GroupBox_frequency_select", "GroupBox", None))
        GroupBox_frequency_select.setTitle(_translate("GroupBox_frequency_select", "Frequency/Period", None))
        self.label_place_holder.setText(_translate("GroupBox_frequency_select", "Place Holder for Frequency/Period Diagram", None))
        self.checkBox_existing_only.setToolTip(_translate("GroupBox_frequency_select", "<html><head/><body><p>Check this to select the frequency/period\n"
"                            that is exists in the data file. The frequency/period that is closest to the cursor will be\n"
"                            selected.</p></body></html>\n"
"                        ", None))
        self.checkBox_existing_only.setText(_translate("GroupBox_frequency_select", "Selecting Existing Frequencies/Periods Only", None))
        self.radioButton_frequency.setText(_translate("GroupBox_frequency_select", "Show Frequency", None))
        self.checkBox_show_existing.setText(_translate("GroupBox_frequency_select", "Show Existing Frequencies/Periods", None))
        self.radioButton_period.setText(_translate("GroupBox_frequency_select", "Show Period", None))
        self.checkBox_y_log_scale.setText(_translate("GroupBox_frequency_select", "Use log-scale on y-axes", None))
        self.pushButton_delete.setText(_translate("GroupBox_frequency_select", "Delete", None))
        self.pushButton_clear.setText(_translate("GroupBox_frequency_select", "Clear", None))
        self.label_tolerance.setText(_translate("GroupBox_frequency_select", "Selection Tolerance", None))
        self.doubleSpinBox_tolerance.setSuffix(_translate("GroupBox_frequency_select", "%", None))
        self.checkBox_x_log_scale.setText(_translate("GroupBox_frequency_select", "Use log-scale on x-axes", None))

