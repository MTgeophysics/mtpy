# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'main_window.ui'
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

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(800, 600)
        MainWindow.setMinimumSize(QtCore.QSize(800, 600))
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuEdit = QtGui.QMenu(self.menubar)
        self.menuEdit.setObjectName(_fromUtf8("menuEdit"))
        self.menuView = QtGui.QMenu(self.menubar)
        self.menuView.setObjectName(_fromUtf8("menuView"))
        self.menuWindow = QtGui.QMenu(self.menubar)
        self.menuWindow.setObjectName(_fromUtf8("menuWindow"))
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setObjectName(_fromUtf8("menuHelp"))
        self.menuTools = QtGui.QMenu(self.menubar)
        self.menuTools.setObjectName(_fromUtf8("menuTools"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionNew_Project = QtGui.QAction(MainWindow)
        self.actionNew_Project.setObjectName(_fromUtf8("actionNew_Project"))
        self.actionOpen_Project = QtGui.QAction(MainWindow)
        self.actionOpen_Project.setObjectName(_fromUtf8("actionOpen_Project"))
        self.actionClose_Project = QtGui.QAction(MainWindow)
        self.actionClose_Project.setEnabled(False)
        self.actionClose_Project.setObjectName(_fromUtf8("actionClose_Project"))
        self.actionOpen_edi_File = QtGui.QAction(MainWindow)
        self.actionOpen_edi_File.setObjectName(_fromUtf8("actionOpen_edi_File"))
        self.actionOpen_edi_Folder = QtGui.QAction(MainWindow)
        self.actionOpen_edi_Folder.setObjectName(_fromUtf8("actionOpen_edi_Folder"))
        self.actionSave_as_Project = QtGui.QAction(MainWindow)
        self.actionSave_as_Project.setEnabled(False)
        self.actionSave_as_Project.setObjectName(_fromUtf8("actionSave_as_Project"))
        self.actionSave_Project = QtGui.QAction(MainWindow)
        self.actionSave_Project.setEnabled(False)
        self.actionSave_Project.setObjectName(_fromUtf8("actionSave_Project"))
        self.actionExit = QtGui.QAction(MainWindow)
        self.actionExit.setStatusTip(_fromUtf8(""))
        self.actionExit.setMenuRole(QtGui.QAction.QuitRole)
        self.actionExit.setObjectName(_fromUtf8("actionExit"))
        self.actionExport = QtGui.QAction(MainWindow)
        self.actionExport.setObjectName(_fromUtf8("actionExport"))
        self.actionShow_Data_Collection = QtGui.QAction(MainWindow)
        self.actionShow_Data_Collection.setCheckable(True)
        self.actionShow_Data_Collection.setChecked(True)
        self.actionShow_Data_Collection.setObjectName(_fromUtf8("actionShow_Data_Collection"))
        self.actionFind_Action = QtGui.QAction(MainWindow)
        self.actionFind_Action.setObjectName(_fromUtf8("actionFind_Action"))
        self.actionHelp = QtGui.QAction(MainWindow)
        self.actionHelp.setObjectName(_fromUtf8("actionHelp"))
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setMenuRole(QtGui.QAction.AboutRole)
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.actionOptions = QtGui.QAction(MainWindow)
        self.actionOptions.setMenuRole(QtGui.QAction.PreferencesRole)
        self.actionOptions.setObjectName(_fromUtf8("actionOptions"))
        self.menuFile.addAction(self.actionNew_Project)
        self.menuFile.addAction(self.actionOpen_Project)
        self.menuFile.addAction(self.actionSave_Project)
        self.menuFile.addAction(self.actionClose_Project)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionOpen_edi_File)
        self.menuFile.addAction(self.actionOpen_edi_Folder)
        self.menuFile.addAction(self.actionSave_as_Project)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExport)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionExit)
        self.menuView.addAction(self.actionShow_Data_Collection)
        self.menuHelp.addAction(self.actionFind_Action)
        self.menuHelp.addAction(self.actionHelp)
        self.menuHelp.addSeparator()
        self.menuHelp.addAction(self.actionAbout)
        self.menuTools.addAction(self.actionOptions)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuEdit.menuAction())
        self.menubar.addAction(self.menuView.menuAction())
        self.menubar.addAction(self.menuWindow.menuAction())
        self.menubar.addAction(self.menuTools.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "SmartMT", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuEdit.setTitle(_translate("MainWindow", "Edit", None))
        self.menuView.setTitle(_translate("MainWindow", "View", None))
        self.menuWindow.setTitle(_translate("MainWindow", "Window", None))
        self.menuHelp.setTitle(_translate("MainWindow", "Help", None))
        self.menuTools.setTitle(_translate("MainWindow", "Tools", None))
        self.actionNew_Project.setText(_translate("MainWindow", "New Project...", None))
        self.actionNew_Project.setShortcut(_translate("MainWindow", "Alt+N", None))
        self.actionOpen_Project.setText(_translate("MainWindow", "Open Project...", None))
        self.actionClose_Project.setText(_translate("MainWindow", "Close Project", None))
        self.actionOpen_edi_File.setText(_translate("MainWindow", "Open .edi File...", None))
        self.actionOpen_edi_Folder.setText(_translate("MainWindow", "Open .edi Folder...", None))
        self.actionSave_as_Project.setText(_translate("MainWindow", "Save as Project...", None))
        self.actionSave_Project.setText(_translate("MainWindow", "Save Project...", None))
        self.actionExit.setText(_translate("MainWindow", "Exit", None))
        self.actionExport.setText(_translate("MainWindow", "Export...", None))
        self.actionShow_Data_Collection.setText(_translate("MainWindow", "Show Data Collection", None))
        self.actionFind_Action.setText(_translate("MainWindow", "Find Action...", None))
        self.actionHelp.setText(_translate("MainWindow", "Help", None))
        self.actionAbout.setText(_translate("MainWindow", "About", None))
        self.actionOptions.setText(_translate("MainWindow", "Options...", None))

