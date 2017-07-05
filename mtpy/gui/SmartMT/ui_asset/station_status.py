# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\station_status.ui'
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

class Ui_StationStatus(object):
    def setupUi(self, StationStatus):
        StationStatus.setObjectName(_fromUtf8("StationStatus"))
        StationStatus.resize(311, 494)
        self.verticalLayout = QtGui.QVBoxLayout(StationStatus)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.tabWidget = QtGui.QTabWidget(StationStatus)
        self.tabWidget.setObjectName(_fromUtf8("tabWidget"))
        self.tab_summary = QtGui.QWidget()
        self.tab_summary.setObjectName(_fromUtf8("tab_summary"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.tab_summary)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.groupBox_2 = QtGui.QGroupBox(self.tab_summary)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.verticalLayout_2.addWidget(self.groupBox_2)
        self.tabWidget.addTab(self.tab_summary, _fromUtf8(""))
        self.tab_detial = QtGui.QWidget()
        self.tab_detial.setObjectName(_fromUtf8("tab_detial"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.tab_detial)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.groupBox = QtGui.QGroupBox(self.tab_detial)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.formLayout = QtGui.QFormLayout(self.groupBox)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.lineEditStation = QtGui.QLineEdit(self.groupBox)
        self.lineEditStation.setReadOnly(True)
        self.lineEditStation.setObjectName(_fromUtf8("lineEditStation"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.lineEditStation)
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.lineEditLat = QtGui.QLineEdit(self.groupBox)
        self.lineEditLat.setReadOnly(True)
        self.lineEditLat.setObjectName(_fromUtf8("lineEditLat"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.lineEditLat)
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.LabelRole, self.label_3)
        self.lineEditLong = QtGui.QLineEdit(self.groupBox)
        self.lineEditLong.setReadOnly(True)
        self.lineEditLong.setObjectName(_fromUtf8("lineEditLong"))
        self.formLayout.setWidget(2, QtGui.QFormLayout.FieldRole, self.lineEditLong)
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.LabelRole, self.label_4)
        self.lineEditElev = QtGui.QLineEdit(self.groupBox)
        self.lineEditElev.setReadOnly(True)
        self.lineEditElev.setObjectName(_fromUtf8("lineEditElev"))
        self.formLayout.setWidget(3, QtGui.QFormLayout.FieldRole, self.lineEditElev)
        self.label_5 = QtGui.QLabel(self.groupBox)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.LabelRole, self.label_5)
        self.lineEditLocation = QtGui.QLineEdit(self.groupBox)
        self.lineEditLocation.setReadOnly(True)
        self.lineEditLocation.setObjectName(_fromUtf8("lineEditLocation"))
        self.formLayout.setWidget(4, QtGui.QFormLayout.FieldRole, self.lineEditLocation)
        self.label_6 = QtGui.QLabel(self.groupBox)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.LabelRole, self.label_6)
        self.lineEditDate_Acq = QtGui.QLineEdit(self.groupBox)
        self.lineEditDate_Acq.setReadOnly(True)
        self.lineEditDate_Acq.setObjectName(_fromUtf8("lineEditDate_Acq"))
        self.formLayout.setWidget(5, QtGui.QFormLayout.FieldRole, self.lineEditDate_Acq)
        self.label_7 = QtGui.QLabel(self.groupBox)
        self.label_7.setObjectName(_fromUtf8("label_7"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.LabelRole, self.label_7)
        self.lineEditAcquired_by = QtGui.QLineEdit(self.groupBox)
        self.lineEditAcquired_by.setReadOnly(True)
        self.lineEditAcquired_by.setObjectName(_fromUtf8("lineEditAcquired_by"))
        self.formLayout.setWidget(6, QtGui.QFormLayout.FieldRole, self.lineEditAcquired_by)
        self.label_8 = QtGui.QLabel(self.groupBox)
        self.label_8.setObjectName(_fromUtf8("label_8"))
        self.formLayout.setWidget(7, QtGui.QFormLayout.LabelRole, self.label_8)
        self.lineEditFile_Ref = QtGui.QLineEdit(self.groupBox)
        self.lineEditFile_Ref.setReadOnly(True)
        self.lineEditFile_Ref.setObjectName(_fromUtf8("lineEditFile_Ref"))
        self.formLayout.setWidget(7, QtGui.QFormLayout.FieldRole, self.lineEditFile_Ref)
        self.verticalLayout_3.addWidget(self.groupBox)
        self.tabWidget.addTab(self.tab_detial, _fromUtf8(""))
        self.verticalLayout.addWidget(self.tabWidget)

        self.retranslateUi(StationStatus)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(StationStatus)

    def retranslateUi(self, StationStatus):
        StationStatus.setWindowTitle(_translate("StationStatus", "Status", None))
        self.groupBox_2.setTitle(_translate("StationStatus", "Summary", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_summary), _translate("StationStatus", "Station Summary", None))
        self.groupBox.setTitle(_translate("StationStatus", "Metadata", None))
        self.label.setText(_translate("StationStatus", "Station", None))
        self.lineEditStation.setPlaceholderText(_translate("StationStatus", "None", None))
        self.label_2.setText(_translate("StationStatus", "Lat (deg)", None))
        self.lineEditLat.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_3.setText(_translate("StationStatus", "Long (deg)", None))
        self.lineEditLong.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_4.setText(_translate("StationStatus", "Elev (m)", None))
        self.lineEditElev.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_5.setText(_translate("StationStatus", "Location", None))
        self.lineEditLocation.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_6.setText(_translate("StationStatus", "Date Acq", None))
        self.lineEditDate_Acq.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_7.setText(_translate("StationStatus", "Acquired by", None))
        self.lineEditAcquired_by.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.label_8.setText(_translate("StationStatus", "File Location", None))
        self.lineEditFile_Ref.setPlaceholderText(_translate("StationStatus", "N/A", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_detial), _translate("StationStatus", "Station Detial", None))

