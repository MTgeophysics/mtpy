# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\u64132\Documents\mtpy2\mtpy\gui\SmartMT\ui_asset\plot_options.ui'
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

class Ui_PlotOption(object):
    def setupUi(self, PlotOption):
        PlotOption.setObjectName(_fromUtf8("PlotOption"))
        PlotOption.setWindowModality(QtCore.Qt.WindowModal)
        PlotOption.resize(400, 700)
        PlotOption.setMinimumSize(QtCore.QSize(340, 600))
        self.verticalLayout_3 = QtGui.QVBoxLayout(PlotOption)
        self.verticalLayout_3.setMargin(0)
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName(_fromUtf8("verticalLayout_3"))
        self.widget = QtGui.QWidget(PlotOption)
        self.widget.setObjectName(_fromUtf8("widget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.widget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.groupBoxPlot_Type = QtGui.QGroupBox(self.widget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBoxPlot_Type.sizePolicy().hasHeightForWidth())
        self.groupBoxPlot_Type.setSizePolicy(sizePolicy)
        self.groupBoxPlot_Type.setMinimumSize(QtCore.QSize(0, 200))
        self.groupBoxPlot_Type.setMaximumSize(QtCore.QSize(16777215, 200))
        self.groupBoxPlot_Type.setObjectName(_fromUtf8("groupBoxPlot_Type"))
        self.verticalLayout_2 = QtGui.QVBoxLayout(self.groupBoxPlot_Type)
        self.verticalLayout_2.setObjectName(_fromUtf8("verticalLayout_2"))
        self.comboBoxSelect_Plot = QtGui.QComboBox(self.groupBoxPlot_Type)
        self.comboBoxSelect_Plot.setObjectName(_fromUtf8("comboBoxSelect_Plot"))
        self.verticalLayout_2.addWidget(self.comboBoxSelect_Plot)
        self.textEditPlot_Description = QtGui.QTextEdit(self.groupBoxPlot_Type)
        self.textEditPlot_Description.setReadOnly(True)
        self.textEditPlot_Description.setObjectName(_fromUtf8("textEditPlot_Description"))
        self.verticalLayout_2.addWidget(self.textEditPlot_Description)
        self.verticalLayout.addWidget(self.groupBoxPlot_Type)
        self.verticalLayout_3.addWidget(self.widget)
        self.pushButton_plot = QtGui.QPushButton(PlotOption)
        self.pushButton_plot.setObjectName(_fromUtf8("pushButton_plot"))
        self.verticalLayout_3.addWidget(self.pushButton_plot)

        self.retranslateUi(PlotOption)
        QtCore.QMetaObject.connectSlotsByName(PlotOption)

    def retranslateUi(self, PlotOption):
        PlotOption.setWindowTitle(_translate("PlotOption", "Plot Option", None))
        self.groupBoxPlot_Type.setTitle(_translate("PlotOption", "Plot Type", None))
        self.pushButton_plot.setText(_translate("PlotOption", "Plot", None))

