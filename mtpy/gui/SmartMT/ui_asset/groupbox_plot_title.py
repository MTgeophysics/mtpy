# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\groupbox_plot_title.ui'
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

class Ui_GroupBox_plot_title(object):
    def setupUi(self, GroupBox_plot_title):
        GroupBox_plot_title.setObjectName(_fromUtf8("GroupBox_plot_title"))
        GroupBox_plot_title.resize(300, 53)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(GroupBox_plot_title.sizePolicy().hasHeightForWidth())
        GroupBox_plot_title.setSizePolicy(sizePolicy)
        self.formLayout = QtGui.QFormLayout(GroupBox_plot_title)
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(GroupBox_plot_title)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.lineEdit_title = QtGui.QLineEdit(GroupBox_plot_title)
        self.lineEdit_title.setObjectName(_fromUtf8("lineEdit_title"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.lineEdit_title)

        self.retranslateUi(GroupBox_plot_title)
        QtCore.QMetaObject.connectSlotsByName(GroupBox_plot_title)

    def retranslateUi(self, GroupBox_plot_title):
        GroupBox_plot_title.setWindowTitle(_translate("GroupBox_plot_title", "GroupBox", None))
        GroupBox_plot_title.setTitle(_translate("GroupBox_plot_title", "Plot Title", None))
        self.label.setToolTip(_translate("GroupBox_plot_title", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))
        self.label.setText(_translate("GroupBox_plot_title", "Title", None))
        self.lineEdit_title.setToolTip(_translate("GroupBox_plot_title", "<html><head/><body><p>The title of the plot, leave blank to\n"
"                            remove title.</p><p><span style=\" font-weight:600;\">Note:</span>\n"
"                            The size and location of the title may depends on the selected plot type.</p></body></html>\n"
"                        ", None))

