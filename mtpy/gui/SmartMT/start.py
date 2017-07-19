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

from mtpy.gui.SmartMT.gui.plot_option import PlotOption
from mtpy.gui.SmartMT.gui.progress_bar import ProgressBar
from mtpy.gui.SmartMT.gui.station_summary import StationSummary
from mtpy.gui.SmartMT.gui.station_viewer import StationViewer
from mtpy.gui.SmartMT.utils.file_handler import FileHandler, FileHandlingException
from mtpy.gui.SmartMT.visualization.visualization_base import MPLCanvasWidget
from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog
from ui_asset.main_window import Ui_SmartMT_MainWindow, _fromUtf8, _translate

DEFAULT_GROUP_NAME = str(_translate("SmartMT_MainWindow", "Default Group", None))


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
        self._subwindow_counter = 0
        self._station_summary = None
        self._progress_bar = ProgressBar(title='Loading files...')
        self.subwindows = {}
        # enable export if the activated subwindow is a image window
        self.ui.mdiArea.subWindowActivated.connect(
            lambda subwindow: self.ui.actionExport.setEnabled(True)
            if subwindow and isinstance(subwindow.widget(), MPLCanvasWidget)
            else self.ui.actionExport.setEnabled(False))

    def setup_menu(self):
        # connect exit menu
        self.ui.actionExit.triggered.connect(QtGui.qApp.quit)
        self.ui.actionOpen_edi_File.triggered.connect(self.file_dialog)
        self.ui.actionOpen_edi_Folder.triggered.connect(self.folder_dialog)
        self.ui.actionShow_Data_Collection.triggered.connect(self._toggle_tree_view)
        self.ui.actionShow_Station_Summary.triggered.connect(self._toggle_station_summary)
        self.ui.actionWindowed_View.triggered.connect(self._toggle_windowed_tabbed_view)
        self.ui.actionTabbed_View.triggered.connect(self._toggle_windowed_tabbed_view)
        self.ui.actionTile_Windows.triggered.connect(self._tile_windows)
        self.ui.actionCascade_Windows.triggered.connect(self._cascade_windows)
        self.ui.actionPlot.triggered.connect(self.plot_selected_station)
        self.ui.actionClose_All_Images.triggered.connect(self._close_all_images)
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

    def _tile_windows(self, *args, **kwargs):
        self.ui.mdiArea.tileSubWindows()

    def _cascade_windows(self, *args, **kwargs):
        self.ui.mdiArea.cascadeSubWindows()

    def _close_all_images(self, *args, **kwargs):
        close_later = []
        for title, (subwindow, action) in self.subwindows.iteritems():
            if title != self._station_viewer.windowTitle() and\
                title != self._station_summary.windowTitle() and\
                    not isinstance(subwindow.widget(), PlotOption):
                    close_later.append(subwindow)
        for subwindow in close_later:
            subwindow.close()

    def _toggle_windowed_tabbed_view(self, *args, **kwargs):
        if self.ui.actionTabbed_View.isEnabled() and self.ui.actionTabbed_View.isChecked():
            self.ui.actionTabbed_View.setEnabled(False)
            self.ui.actionWindowed_View.setEnabled(True)
            self.ui.actionWindowed_View.setChecked(False)
            self.ui.actionTile_Windows.setEnabled(False)
            self.ui.actionCascade_Windows.setEnabled(False)
            self.ui.mdiArea.setViewMode(QtGui.QMdiArea.TabbedView)
        elif self.ui.actionWindowed_View.isEnabled() and self.ui.actionWindowed_View.isChecked():
            self.ui.actionWindowed_View.setEnabled(False)
            self.ui.actionTabbed_View.setEnabled(True)
            self.ui.actionTabbed_View.setChecked(False)
            self.ui.actionTile_Windows.setEnabled(True)
            self.ui.actionCascade_Windows.setEnabled(True)
            self.ui.mdiArea.setViewMode(QtGui.QMdiArea.SubWindowView)

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
            self._progress_bar.setMaximumValue(len(file_list))
            self._progress_bar.onStart()
            self._add_files(file_list, DEFAULT_GROUP_NAME)
            self._update_tree_view()
            self._progress_bar.onFinished()

    def _add_files(self, file_list, group_id=DEFAULT_GROUP_NAME):
        for file_ref in file_list:
            try:
                self._file_handler.add_file(os.path.abspath(str(file_ref)), group_id=group_id)
            except FileHandlingException as exp:
                self._logger.warning(exp.message)
            except Exception as exp:
                self._logger.critical(exp.message)
            self._progress_bar.incrementValue()

    def plot_selected_station(self, *args, **kwargs):
        if self._station_viewer and self._station_viewer.fig_canvas.selected_stations:
            plot_option = PlotOption(self, self._file_handler, self._station_viewer.fig_canvas.selected_stations)
            subwindow, _ = self.create_subwindow(plot_option,
                                                 plot_option.windowTitle(),
                                                 tooltip="plot stations:\n%s" % "\n".join(
                                                     self._station_viewer.fig_canvas.selected_stations))
        else:
            self._logger.info("nothing to plot")

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
                    self._progress_bar.setMaximumValue(len(file_list))
                    self._progress_bar.onStart()
                    self._add_files(file_list, os.path.basename(dir_name))
                    self._update_tree_view()
                    self._progress_bar.onFinished()
            else:
                break

    def _toggle_tree_view(self):
        if self._station_viewer:
            subwindow = self.subwindows[self._station_viewer.windowTitle()][0]
            if self.ui.actionShow_Data_Collection.isEnabled():
                if self.ui.actionShow_Data_Collection.isChecked():
                    subwindow.show()
                else:
                    subwindow.hide()

    def _toggle_station_summary(self):
        if self._station_summary:
            subwindow = self.subwindows[self._station_summary.windowTitle()][0]
            if self.ui.actionShow_Station_Summary.isEnabled():
                if self.ui.actionShow_Station_Summary.isChecked():
                    subwindow.show()
                else:
                    subwindow.hide()

    def _update_tree_view(self):
        if not self._station_viewer:
            self._station_viewer = StationViewer(self, self._file_handler)
            self.ui.actionShow_Data_Collection.setEnabled(True)
        if not self._station_summary:
            self._station_summary = StationSummary(self, self._file_handler, self._station_viewer.fig_canvas.selected_stations)
            self.ui.actionShow_Station_Summary.setEnabled(True)
            # connect to tree view to update summary
            self._station_viewer.ui.treeWidget_stations.selectionModel().selectionChanged.connect(self._station_summary.update_view)
            self._station_viewer.setFocus()
        self._station_viewer.update_view()
        self._station_summary.update_view()

    def create_subwindow(self, widget, title, overide=True, tooltip=None):
        subwindow = None
        self._subwindow_counter += 1
        if title in self.subwindows:
            if overide:
                subwindow, window_action = self.subwindows[title]
                subwindow.close()
                self.ui.menuWindow.removeAction(window_action)
                del self.subwindows[title]
            else:
                # find a new window title by adding a number after the title
                counter = 1
                new_title = "%s %d" % (title, counter)
                while new_title in self.subwindows:
                    counter += 1
                    new_title = "%s %d" % (title, counter)
                title = new_title

        widget.setAttribute(QtCore.Qt.WA_DeleteOnClose, True)
        subwindow = StartQt4.MDISubWindow(self)
        subwindow.setWindowTitle(title)
        if tooltip:
            subwindow.setToolTip(tooltip)
        subwindow.setWidget(widget)
        subwindow.resize(widget.size())
        self.ui.mdiArea.addSubWindow(subwindow)

        # create menu action
        new_window_action = QtGui.QAction(self)
        new_window_action.setObjectName(_fromUtf8("actionSubwindow%d" % self._subwindow_counter))
        new_window_action.setText(_translate("SmartMT_MainWindow", title, None))
        new_window_action.triggered.connect(subwindow.show_and_focus)

        # add to window menu
        self.ui.menuWindow.addAction(new_window_action)
        # add all references to self._subwindow
        self.subwindows[title] = (subwindow, new_window_action)
        subwindow.show()

        return subwindow, new_window_action

    class MDISubWindow(QtGui.QMdiSubWindow):
        def __init__(self, main_ui, parent=None, flags=0):
            """

            :param main_ui:
            :type main_ui: StartQt4
            :param parent:
            :param flags:
            """
            super(QtGui.QMdiSubWindow, self).__init__(parent)
            self._main_ui = main_ui

        def show_and_focus(self):
            self.show()
            self.setFocus()

        def closeEvent(self, QCloseEvent):
            """
            override the close event to remove the window menu action
            :param QCloseEvent:
            :return:
            """
            if self._main_ui.subwindows:
                title = str(self.windowTitle())
                if title in self._main_ui.subwindows:
                    (subwindow, window_action) = self._main_ui.subwindows[title]
                    if subwindow.widget().testAttribute(QtCore.Qt.WA_DeleteOnClose):
                        # remove menu action
                        self._main_ui.ui.menuWindow.removeAction(window_action)
                        # remove the window
                        del self._main_ui.subwindows[title]
                    else:
                        subwindow.hide()

    @deprecated(
        "This is an dummy action that should be used only when the required function hasn't been implemented")
    def dummy_action(self, *args, **kwargs):
        pass


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    smartMT = StartQt4()
    smartMT.show()
    sys.exit(app.exec_())
