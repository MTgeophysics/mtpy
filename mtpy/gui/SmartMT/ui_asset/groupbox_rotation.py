# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_rotation.ui'
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

class Ui_GroupBox_Rotation(object):
    def setupUi(self, GroupBox_Rotation):
        GroupBox_Rotation.setObjectName(_fromUtf8("GroupBox_Rotation"))
        GroupBox_Rotation.resize(300, 233)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_Rotation.sizePolicy().hasHeightForWidth())
        GroupBox_Rotation.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_Rotation)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.label = QtGui.QLabel(GroupBox_Rotation)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(7)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName(_fromUtf8("label"))
        self.verticalLayout.addWidget(self.label)
        self.dial_rotation = QtGui.QDial(GroupBox_Rotation)
        self.dial_rotation.setMaximum(360)
        self.dial_rotation.setSingleStep(1)
        self.dial_rotation.setPageStep(10)
        self.dial_rotation.setProperty("value", 180)
        self.dial_rotation.setWrapping(True)
        self.dial_rotation.setNotchesVisible(True)
        self.dial_rotation.setObjectName(_fromUtf8("dial_rotation"))
        self.verticalLayout.addWidget(self.dial_rotation)
        self.doubleSpinBox_rotation = QtGui.QDoubleSpinBox(GroupBox_Rotation)
        self.doubleSpinBox_rotation.setMaximum(360.0)
        self.doubleSpinBox_rotation.setObjectName(_fromUtf8("doubleSpinBox_rotation"))
        self.verticalLayout.addWidget(self.doubleSpinBox_rotation)

        self.retranslateUi(GroupBox_Rotation)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_Rotation)

    def retranslateUi(self, GroupBox_Rotation):
        GroupBox_Rotation.setWindowTitle(_translate("GroupBox_Rotation", "GroupBox", None))
        GroupBox_Rotation.setToolTip(
            _translate("GroupBox_Rotation", "<html><head/><body><p>Rotation angle (in degree)</p><p><span\n"
                                            "                style=\" font-weight:600;\">Note:</span> All angles are referenced to geographic North,\n"
                                            "                positive in clockwise direction. </p></body></html>\n"
                                            "            ", None))
        GroupBox_Rotation.setTitle(_translate("GroupBox_Rotation", "Rotation", None))
        self.label.setText(_translate("GroupBox_Rotation", "North", None))
        self.doubleSpinBox_rotation.setSuffix(_translate("GroupBox_Rotation", "°", None))

