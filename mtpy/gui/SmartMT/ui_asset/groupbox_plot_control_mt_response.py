# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'groupbox_plot_control_mt_response.ui'
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

class Ui_GroupBox_plot_control_mt_response(object):
    def setupUi(self, GroupBox_plot_control_mt_response):
        GroupBox_plot_control_mt_response.setObjectName(_fromUtf8("GroupBox_plot_control_mt_response"))
        GroupBox_plot_control_mt_response.resize(312, 448)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_plot_control_mt_response.sizePolicy().hasHeightForWidth())
        GroupBox_plot_control_mt_response.setSizePolicy(sizePolicy)
        self.verticalLayout = QtGui.QVBoxLayout(GroupBox_plot_control_mt_response)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBox_plot_style = QtGui.QGroupBox(GroupBox_plot_control_mt_response)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_plot_style.sizePolicy().hasHeightForWidth())
        self.groupBox_plot_style.setSizePolicy(sizePolicy)
        self.groupBox_plot_style.setObjectName(_fromUtf8("groupBox_plot_style"))
        self.horizontalLayout = QtGui.QHBoxLayout(self.groupBox_plot_style)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.radioButton_compare = QtGui.QRadioButton(self.groupBox_plot_style)
        self.radioButton_compare.setChecked(True)
        self.radioButton_compare.setObjectName(_fromUtf8("radioButton_compare"))
        self.horizontalLayout.addWidget(self.radioButton_compare)
        self.radioButton_all = QtGui.QRadioButton(self.groupBox_plot_style)
        self.radioButton_all.setObjectName(_fromUtf8("radioButton_all"))
        self.horizontalLayout.addWidget(self.radioButton_all)
        self.verticalLayout.addWidget(self.groupBox_plot_style)
        self.groupBox = QtGui.QGroupBox(GroupBox_plot_control_mt_response)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName(_fromUtf8("groupBox"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.radioButton_1 = QtGui.QRadioButton(self.groupBox)
        self.radioButton_1.setChecked(True)
        self.radioButton_1.setObjectName(_fromUtf8("radioButton_1"))
        self.verticalLayout_2.addWidget(self.radioButton_1)
        self.radioButton_2 = QtGui.QRadioButton(self.groupBox)
        self.radioButton_2.setObjectName(_fromUtf8("radioButton_2"))
        self.verticalLayout_2.addWidget(self.radioButton_2)
        self.radioButton_3 = QtGui.QRadioButton(self.groupBox)
        self.radioButton_3.setObjectName(_fromUtf8("radioButton_3"))
        self.verticalLayout_2.addWidget(self.radioButton_3)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtGui.QGroupBox(GroupBox_plot_control_mt_response)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setObjectName(_fromUtf8("groupBox_2"))
        self.verticalLayout_3 = QtGui.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.checkBox_strike_i = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_strike_i.setObjectName(_fromUtf8("checkBox_strike_i"))
        self.verticalLayout_3.addWidget(self.checkBox_strike_i)
        self.checkBox_strike_p = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_strike_p.setObjectName(_fromUtf8("checkBox_strike_p"))
        self.verticalLayout_3.addWidget(self.checkBox_strike_p)
        self.checkBox_strike_t = QtGui.QCheckBox(self.groupBox_2)
        self.checkBox_strike_t.setObjectName(_fromUtf8("checkBox_strike_t"))
        self.verticalLayout_3.addWidget(self.checkBox_strike_t)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.groupBox_3 = QtGui.QGroupBox(GroupBox_plot_control_mt_response)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_3.sizePolicy().hasHeightForWidth())
        self.groupBox_3.setSizePolicy(sizePolicy)
        self.groupBox_3.setObjectName(_fromUtf8("groupBox_3"))
        self.gridLayout = QtGui.QGridLayout(self.groupBox_3)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.radioButton_skew_y = QtGui.QRadioButton(self.groupBox_3)
        self.radioButton_skew_y.setObjectName(_fromUtf8("radioButton_skew_y"))
        self.gridLayout.addWidget(self.radioButton_skew_y, 0, 0, 1, 1)
        self.radioButton_skew_n = QtGui.QRadioButton(self.groupBox_3)
        self.radioButton_skew_n.setChecked(True)
        self.radioButton_skew_n.setObjectName(_fromUtf8("radioButton_skew_n"))
        self.gridLayout.addWidget(self.radioButton_skew_n, 0, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_3)
        self.groupBox_4 = QtGui.QGroupBox(GroupBox_plot_control_mt_response)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setObjectName(_fromUtf8("groupBox_4"))
        self.verticalLayout_4 = QtGui.QVBoxLayout(self.groupBox_4)
        self.verticalLayout_4.setObjectName(_fromUtf8("verticalLayout_4"))
        self.radioButton_ellipses_n = QtGui.QRadioButton(self.groupBox_4)
        self.radioButton_ellipses_n.setChecked(True)
        self.radioButton_ellipses_n.setObjectName(_fromUtf8("radioButton_ellipses_n"))
        self.verticalLayout_4.addWidget(self.radioButton_ellipses_n)
        self.radioButton_ellipses_y = QtGui.QRadioButton(self.groupBox_4)
        self.radioButton_ellipses_y.setObjectName(_fromUtf8("radioButton_ellipses_y"))
        self.verticalLayout_4.addWidget(self.radioButton_ellipses_y)
        self.verticalLayout.addWidget(self.groupBox_4)

        self.retranslateUi(GroupBox_plot_control_mt_response)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_plot_control_mt_response)

    def retranslateUi(self, GroupBox_plot_control_mt_response):
        GroupBox_plot_control_mt_response.setWindowTitle(_translate("GroupBox_plot_control_mt_response", "GroupBox", None))
        GroupBox_plot_control_mt_response.setTitle(_translate("GroupBox_plot_control_mt_response", "Plot Control", None))
        self.groupBox_plot_style.setTitle(_translate("GroupBox_plot_control_mt_response", "Plot Style", None))
        self.radioButton_compare.setText(_translate("GroupBox_plot_control_mt_response", "Compare", None))
        self.radioButton_all.setText(_translate("GroupBox_plot_control_mt_response", "Side-by-Side", None))
        self.groupBox.setTitle(_translate("GroupBox_plot_control_mt_response", "Plot Type", None))
        self.radioButton_1.setText(_translate("GroupBox_plot_control_mt_response", "plot just Ex/By and Ey/Bx", None))
        self.radioButton_2.setText(_translate("GroupBox_plot_control_mt_response", "plot all 4 components", None))
        self.radioButton_3.setText(_translate("GroupBox_plot_control_mt_response", "plot off diagonal and the deeterminant", None))
        self.groupBox_2.setToolTip(_translate("GroupBox_plot_control_mt_response", "<html><head/><body><p>Plots thee strike angle from different\n"
"                            parameters</p></body></html>\n"
"                        ", None))
        self.groupBox_2.setTitle(_translate("GroupBox_plot_control_mt_response", "Strike", None))
        self.checkBox_strike_i.setText(_translate("GroupBox_plot_control_mt_response", "plot strike angle determined from the invariants of\n"
"Weaver et al. [2000]", None))
        self.checkBox_strike_p.setText(_translate("GroupBox_plot_control_mt_response", "plot strike angle determined from the phase tensor\n"
"of Caldwell et al. [2004]", None))
        self.checkBox_strike_t.setText(_translate("GroupBox_plot_control_mt_response", "plot strike angle determined from the tipper", None))
        self.groupBox_3.setToolTip(_translate("GroupBox_plot_control_mt_response", "<html><head/><body><p>plot the skew angle calculated from the\n"
"                            phase tensor.</p></body></html>\n"
"                        ", None))
        self.groupBox_3.setTitle(_translate("GroupBox_plot_control_mt_response", "Skew", None))
        self.radioButton_skew_y.setText(_translate("GroupBox_plot_control_mt_response", "plot skew angle", None))
        self.radioButton_skew_n.setText(_translate("GroupBox_plot_control_mt_response", "do not plot skew angle", None))
        self.groupBox_4.setTitle(_translate("GroupBox_plot_control_mt_response", "Ellipses", None))
        self.radioButton_ellipses_n.setText(_translate("GroupBox_plot_control_mt_response", "do not plot phase tensor pllipses", None))
        self.radioButton_ellipses_y.setText(_translate("GroupBox_plot_control_mt_response", "plot phase tensor ellipses", None))

