# -*- coding: utf-8 -*-
"""
    Description:
        This is the main driver of the SmartMT gui

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""

import os
import sys

from PyQt4 import QtCore, QtGui

from main_window import Ui_SmartMT_MainWindow, _fromUtf8, _translate
from file_handler import FileHandler, FileHandlingException
from station_viewer import Ui_StationViewer

from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog


class StartQt4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self._is_file_dialog_opened = False
        self.ui = Ui_SmartMT_MainWindow()
        self.ui.setupUi(self)
        self.setup_menu()
        self._file_handler = FileHandler()
        self._station_viewer = None
        self._subwindows = {}

    def setup_menu(self):
        # connect exit menu
        self.ui.actionExit.triggered.connect(QtGui.qApp.quit)
        self.ui.actionOpen_edi_File.triggered.connect(self.file_dialog)
        self.ui.actionOpen_edi_Folder.triggered.connect(self.folder_dialog)
        # not yet impleneted
        self.ui.actionAbout.triggered.connect(self.dummy_action)
        self.ui.actionClose_Project.triggered.connect(self.dummy_action)
        self.ui.actionExport.triggered.connect(self.dummy_action)
        self.ui.actionFind_Action.triggered.connect(self.dummy_action)
        self.ui.actionHelp.triggered.connect(self.dummy_action)
        self.ui.actionNew_Project.triggered.connect(self.dummy_action)
        self.ui.actionOpen_Project.triggered.connect(self.dummy_action)
        self.ui.actionOptions.triggered.connect(self.dummy_action)
        self.ui.actionSave_as_Project.triggered.connect(self.dummy_action)
        self.ui.actionSave_Project.triggered.connect(self.dummy_action)
        self.ui.actionShow_Data_Collection.triggered.connect(self.dummy_action)

    def file_dialog(self, *args, **kwargs):
        dialog = QtGui.QFileDialog(self)
        if not self._is_file_dialog_opened:
            # set the initial directory to HOME
            dialog.setDirectory(os.path.expanduser("~"))
            self._is_file_dialog_opened = True
        dialog.setWindowTitle('Open .edi Files...')
        dialog.setNameFilter('.edi files (*.edi)')
        dialog.setFileMode(QtGui.QFileDialog.ExistingFiles)
        if dialog.exec_() == QtGui.QDialog.Accepted:
            file_list = dialog.selectedFiles()
            self._add_files(file_list)
            self._update_tree_view()

    def _add_files(self, file_list):
        for file_ref in file_list:
            try:
                self._file_handler.add_file(str(file_ref), group_id='Default Group')
            except FileHandlingException as exp:
                self._logger.warning(exp.message)
            except Exception as exp:
                self._logger.critical(exp.message)

    def folder_dialog(self, *args, **kwargs):
        dialog = QtGui.QFileDialog(self)
        if not self._is_file_dialog_opened:
            # set the initial directory to HOME
            dialog.setDirectory(os.path.expanduser("~"))
            self._is_file_dialog_opened = True
        dir_name = None
        dialog.setWindowTitle("Open .edi Directory...")
        dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        while dir_name is None:
            if dialog.exec_() == QtGui.QDialog.Accepted:
                dir_name = dialog.selectedFiles()[0]
                dir_name = str(dir_name)
                file_list = [os.path.join(dir_name, edi) for edi in os.listdir(dir_name) if edi.endswith("edi")]
                if not file_list:
                    # empty list
                    QtGui.QMessageBox.information(self, "NOTE",
                                                  "Directory does not contain any .edi file, please select again.")
                    dir_name = None  # will read again
                else:
                    self._add_files(dir_name)
                    self._update_tree_view()
            else:
                break

    def _update_tree_view(self):
        if not self._station_viewer:
            widget = QtGui.QWidget()
            self._station_viewer = Ui_StationViewer()
            self._station_viewer.setupUi(widget)
            self.create_subwindow(widget, "Stations Stats")
        for group in self._file_handler.get_groups():
            node = QtGui.QTreeWidgetItem(self._station_viewer.treeWidget_stations)
            node.setText(0, group)
            for ref in self._file_handler.get_group_members(group):
                mt_obj = self._file_handler.get_MT_obj(ref)
                child_node = QtGui.QTreeWidgetItem(node)
                child_node.setText(1, mt_obj.station)
                child_node.setText(2, ref)

    def create_subwindow(self, widget, title):
        if title not in self._subwindows:
            subwindow = QtGui.QMdiSubWindow()
            subwindow.setWindowTitle(title)
            subwindow.setWidget(widget)
            # todo connect subwindow close event to remove menu item
            self.ui.mdiArea.addSubWindow(subwindow)

            # create menu action
            new_window_action = QtGui.QAction(self)
            new_window_action.setObjectName(_fromUtf8("actionSubwindow%d" % (len(self._subwindows) + 1)))
            new_window_action.setText(_translate("SmartMT_MainWindow", title, None))
            new_window_action.triggered.connect(subwindow.setFocus)

            # add to window menu
            self.ui.menuWindow.addAction(new_window_action)

            # add all references to self._subwindow
            self._subwindows[title]=(subwindow, widget, new_window_action)

            subwindow.show()

    @deprecated(
        "This is an dummy action that should be unsed only when the required function hasn't been implemented")
    def dummy_action(self, *args, **kwargs):
        pass


if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    smartMT = StartQt4()
    smartMT.show()
    sys.exit(app.exec_())
