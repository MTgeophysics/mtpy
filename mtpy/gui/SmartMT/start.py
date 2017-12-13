# -*- coding: utf-8 -*-
"""
    Description:
        This is the main driver of the SmartMT gui

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
import inspect
import os
import sys
import webbrowser

import sip
from qtpy import QtCore
from qtpy.QtWidgets import QMainWindow, QWidget, qApp, QMessageBox, QFileDialog, QDialog, QWizard, QMdiArea, QAction, \
    QMdiSubWindow, QApplication

from qtpy import QT_VERSION
import matplotlib

if QT_VERSION.startswith("4"):
    matplotlib.use("Qt4Agg")
elif QT_VERSION.startswith("5"):
    matplotlib.use("Qt5Agg")

from mtpy.core.edi_collection import EdiCollection
from mtpy.gui.SmartMT.gui.busy_indicators import ProgressBar
from mtpy.gui.SmartMT.gui.export_dialog import ExportDialog
from mtpy.gui.SmartMT.gui.export_dialog_modem import ExportDialogModEm
from mtpy.gui.SmartMT.gui.plot_option import PlotOption
from mtpy.gui.SmartMT.gui.station_summary import StationSummary
from mtpy.gui.SmartMT.gui.station_viewer import StationViewer
from mtpy.gui.SmartMT.ui_asset.main_window import Ui_SmartMT_MainWindow
from mtpy.gui.SmartMT.utils.file_handler import FileHandler, FileHandlingException
from mtpy.gui.SmartMT.visualization.visualization_base import MPLCanvasWidget
from mtpy.utils.decorator import deprecated
from mtpy.utils.mtpylog import MtPyLog

_translate = QtCore.QCoreApplication.translate
DEFAULT_GROUP_NAME = str(_translate("SmartMT_MainWindow", "Default Group", None))


class StartGUI(QMainWindow):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        self._file_handler = FileHandler()
        self._is_file_dialog_opened = False

        self.ui = Ui_SmartMT_MainWindow()

        # init gui
        self.ui.setupUi(self)
        self.setup_menu()

        # export dialogs
        self._export_dialog = ExportDialog(self)
        self._export_dialog_modem = ExportDialogModEm(self)

        # set station viewer
        self._station_viewer = StationViewer(self, file_handler=self._file_handler)
        self._station_viewer.ui.pushButton_plot.clicked.connect(self.plot_selected_station)
        self._station_viewer.selection_changed.connect(self._selected_station_changed)
        self.ui.stackedWidget.addWidget(self._station_viewer)

        # set plot option widget
        self._plot_option = PlotOption(self, self._file_handler, self._station_viewer.selected_stations)
        self._plot_option.ui.pushButton_back.clicked.connect(
            lambda x: self.ui.stackedWidget.setCurrentWidget(self._station_viewer))
        self._station_viewer.selection_changed.connect(self._plot_option.data_changed)
        self.ui.stackedWidget.addWidget(self._plot_option)

        self.ui.stackedWidget.setCurrentWidget(self._station_viewer)

        self._subwindow_counter = 0

        self._station_summary = None

        self._progress_bar = ProgressBar(title='Loading files...')
        self.subwindows = {}
        # enable export if the activated subwindow is a image window
        self.ui.mdiArea.subWindowActivated.connect(self._subwindow_activated)
        self.setWindowState(QtCore.Qt.WindowMaximized)

    def _subwindow_activated(self, subwindow):
        if subwindow and isinstance(subwindow.widget(), MPLCanvasWidget):
            self.ui.actionExport.setEnabled(True)
        else:
            self.ui.actionExport.setEnabled(False)

    def setup_menu(self):
        # connect exit menu
        self.ui.actionExit.triggered.connect(qApp.quit)
        self.ui.actionOpen_edi_File.triggered.connect(self.file_dialog)
        self.ui.actionOpen_edi_Folder.triggered.connect(self.folder_dialog)
        self.ui.actionShow_Station_Summary.triggered.connect(self._toggle_station_summary)
        self.ui.actionWindowed_View.triggered.connect(self._toggle_windowed_tabbed_view)
        self.ui.actionTabbed_View.triggered.connect(self._toggle_windowed_tabbed_view)
        self.ui.actionTile_Windows.triggered.connect(self._tile_windows)
        self.ui.actionCascade_Windows.triggered.connect(self._cascade_windows)
        self.ui.actionPlot.triggered.connect(self.plot_selected_station)
        self.ui.actionClose_All_Images.triggered.connect(self._close_all_images)
        self.ui.actionExport.triggered.connect(self._export_image)
        self.ui.actionExport_ModEM_Data.triggered.connect(self._export_modem)
        self.ui.actionCreate_Shape_File_From_Stations.triggered.connect(self._export_shape_file)
        self.ui.actionCreate_Phase_Tensor_csv_file.triggered.connect(self._export_phase_tensor_csv)
        self.ui.actionCreate_Measurement_csv_file.triggered.connect(self._export_measurement_csv)
        # not yet impleneted
        self.ui.actionAbout.triggered.connect(self.dummy_action)
        self.ui.actionClose_Project.triggered.connect(self.dummy_action)
        self.ui.actionFind_Action.triggered.connect(self.dummy_action)
        self.ui.actionHelp.triggered.connect(self.dummy_action)
        self.ui.actionNew_Project.triggered.connect(self.dummy_action)
        self.ui.actionOpen_Project.triggered.connect(self.dummy_action)
        self.ui.actionOptions.triggered.connect(self.dummy_action)
        self.ui.actionSave_as_Project.triggered.connect(self.dummy_action)
        self.ui.actionSave_Project.triggered.connect(self.dummy_action)

    def _export_measurement_csv(self):
        # show files
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("You are about to create measurement .csv files.")
        msg.setInformativeText("Please select an output directory after click \"OK\"\n"
                               "For the list of .edi files (stations) included in the creation, please click \"Show Details\"")
        msg.setWindowTitle("Note")
        msg.setDetailedText(
            "\n".join(["{station} ({fn})".format(
                station=station, fn=self._file_handler.station2ref(station)
            ) for station in self._station_viewer.selected_stations])
        )
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        if msg.exec_() == QMessageBox.Ok:
            dialog = QFileDialog(self)
            dir_name = None
            dialog.setWindowTitle("Selecting Output Directory ...")
            dialog.setFileMode(QFileDialog.DirectoryOnly)
            while dir_name is None:
                if dialog.exec_() == QDialog.Accepted:
                    dir_name = dialog.selectedFiles()[0]
                    dir_name = str(dir_name)
                    if not os.path.isdir(dir_name):
                        QMessageBox.information(self, "NOTE",
                                                "Please select a directory to save the created .csv files.")
                        dir_name = None  # will read again
                else:
                    break
            if dir_name is not None:
                collect = EdiCollection(
                    mt_objs=[
                        self._file_handler.get_MT_obj(self._file_handler.station2ref(station))
                        for station in self._station_viewer.selected_stations
                    ]
                )
                collect.create_measurement_csv(dir_name)
                QMessageBox.information(self, "Creation Completed", "Output written to %s" % dir_name)
                webbrowser.open(dir_name)

    def _export_phase_tensor_csv(self):
        # show files
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("You are about to create phase tensor .csv files.")
        msg.setInformativeText("Please select an output directory after click \"OK\"\n"
                               "For the list of .edi files (stations) included in the creation, please click \"Show Details\"")
        msg.setWindowTitle("Note")
        msg.setDetailedText(
            "\n".join(["{station} ({fn})".format(
                station=station, fn=self._file_handler.station2ref(station)
            ) for station in self._station_viewer.selected_stations])
        )
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        if msg.exec_() == QMessageBox.Ok:
            dialog = QFileDialog(self)
            dir_name = None
            dialog.setWindowTitle("Selecting Output Directory ...")
            dialog.setFileMode(QFileDialog.DirectoryOnly)
            while dir_name is None:
                if dialog.exec_() == QDialog.Accepted:
                    dir_name = dialog.selectedFiles()[0]
                    dir_name = str(dir_name)
                    if not os.path.isdir(dir_name):
                        QMessageBox.information(self, "NOTE",
                                                "Please select a directory to save the created .csv files.")
                        dir_name = None  # will read again
                else:
                    break
            if dir_name is not None:
                collect = EdiCollection(
                    mt_objs=[
                        self._file_handler.get_MT_obj(self._file_handler.station2ref(station))
                        for station in self._station_viewer.selected_stations
                    ]
                )
                collect.create_phase_tensor_csv(dir_name)
                QMessageBox.information(self, "Creation Completed", "Output written to %s" % dir_name)
                webbrowser.open(dir_name)

    def _export_shape_file(self, *args, **kwargs):
        # show files
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Information)
        msg.setText("You are about to create shape files.")
        msg.setInformativeText("Please select an output directory after click \"OK\"\n"
                               "For the list of .edi files (stations) included in the creation, please click \"Show Details\"")
        msg.setWindowTitle("Note")
        msg.setDetailedText(
            "\n".join(["{station} ({fn})".format(
                station=station, fn=self._file_handler.station2ref(station)
            ) for station in self._station_viewer.selected_stations])
        )
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

        if msg.exec_() == QMessageBox.Ok:
            dialog = QFileDialog(self)
            dir_name = None
            dialog.setWindowTitle("Selecting Output Directory ...")
            dialog.setFileMode(QFileDialog.DirectoryOnly)
            while dir_name is None:
                if dialog.exec_() == QDialog.Accepted:
                    dir_name = dialog.selectedFiles()[0]
                    dir_name = str(dir_name)
                    if not os.path.isdir(dir_name):
                        QMessageBox.information(self, "NOTE",
                                                "Please select a directory to save the created shape files.")
                        dir_name = None  # will read again
                else:
                    break
            if dir_name is not None:
                collect = EdiCollection(
                    mt_objs=[
                        self._file_handler.get_MT_obj(self._file_handler.station2ref(station))
                        for station in self._station_viewer.selected_stations
                    ]
                )
                collect.create_mt_station_gdf(dir_name)
                QMessageBox.information(self, "Creation Completed", "Output written to %s" % dir_name)
                webbrowser.open(dir_name)

    def _export_image(self, *args, **kwargs):
        subwindow = self.ui.mdiArea.activeSubWindow()
        widget = subwindow.widget()
        if isinstance(widget, MPLCanvasWidget):
            try:
                self._export_dialog.export_to_file(widget.get_fig())
            except Exception as e:
                frm = inspect.trace()[-1]
                mod = inspect.getmodule(frm[0])
                QMessageBox.critical(self,
                                     'Exporting Error',
                                     "{}: {}".format(mod.__name__, e.message),
                                     QMessageBox.Close)

    def _export_modem(self, *args, **kwargs):
        mt_objs = []
        for selected_station in self._station_viewer.selected_stations:
            ref = self._file_handler.station2ref(selected_station)
            mt_obj = self._file_handler.get_MT_obj(ref)
            mt_objs.append(mt_obj)
        self._export_dialog_modem.set_data(mt_objs)
        self._export_dialog_modem.restart()
        if self._export_dialog_modem.exec_() == QWizard.Accepted:
            self._export_dialog_modem.export_data()

    def _tile_windows(self, *args, **kwargs):
        self.ui.mdiArea.tileSubWindows()

    def _cascade_windows(self, *args, **kwargs):
        self.ui.mdiArea.cascadeSubWindows()

    def _close_all_images(self, *args, **kwargs):
        close_later = []
        for title, (subwindow, action) in self.subwindows.iteritems():
            if not isinstance(subwindow.widget(), StationSummary):
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
            self.ui.mdiArea.setViewMode(QMdiArea.TabbedView)
        elif self.ui.actionWindowed_View.isEnabled() and self.ui.actionWindowed_View.isChecked():
            self.ui.actionWindowed_View.setEnabled(False)
            self.ui.actionTabbed_View.setEnabled(True)
            self.ui.actionTabbed_View.setChecked(False)
            self.ui.actionTile_Windows.setEnabled(True)
            self.ui.actionCascade_Windows.setEnabled(True)
            self.ui.mdiArea.setViewMode(QMdiArea.SubWindowView)

    def file_dialog(self, *args, **kwargs):
        dialog = QFileDialog(self)
        if not self._is_file_dialog_opened:
            # set the initial directory to HOME
            dialog.setDirectory(os.path.expanduser("~"))
            self._is_file_dialog_opened = True
        dialog.setWindowTitle('Open .edi Files...')
        dialog.setNameFilter('.edi files (*.edi)')
        dialog.setFileMode(QFileDialog.ExistingFiles)
        if dialog.exec_() == QDialog.Accepted:
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
        if self._station_viewer and self._station_viewer.selected_stations:
            self.ui.stackedWidget.setCurrentWidget(self._plot_option)
        else:
            self._logger.info("nothing to plot")
            QMessageBox.information(self, "NOTE",
                                    "Please load and select station data for visualization.")

    def folder_dialog(self, *args, **kwargs):
        dialog = QFileDialog(self)
        if not self._is_file_dialog_opened:
            # set the initial directory to HOME
            dialog.setDirectory(os.path.expanduser("~"))
            self._is_file_dialog_opened = True
        dir_name = None
        dialog.setWindowTitle("Open .edi Directory...")
        dialog.setFileMode(QFileDialog.DirectoryOnly)
        while dir_name is None:
            if dialog.exec_() == QDialog.Accepted:
                dir_name = dialog.selectedFiles()[0]
                dir_name = str(dir_name)
                file_list = [os.path.join(dir_name, edi) for edi in os.listdir(dir_name) if edi.endswith("edi")]
                if not file_list:
                    # empty list
                    QMessageBox.information(self, "NOTE",
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

    def _toggle_station_summary(self):
        self._update_station_summary()
        subwindow = self.subwindows[self._station_summary.windowTitle()][0]
        if self.ui.actionShow_Station_Summary.isChecked():
            subwindow.show_and_focus()
        else:
            subwindow.hide()

    def _update_tree_view(self):
        self._station_viewer.update_view()

    def _update_station_summary(self):
        if not self._station_summary:
            self._station_summary = StationSummary(self, self._file_handler,
                                                   self._station_viewer.selected_stations)
            # connect to tree view to update summary
            self._station_viewer.ui.treeWidget_stations.selectionModel().selectionChanged.connect(
                self._station_summary.update_view)
            # connect to handle selection_changed signal from station_viewer
            self._station_viewer.selection_changed.connect(self._station_summary.update_view)

    def _selected_station_changed(self):
        enable = bool(self._station_viewer.selected_stations)
        self.ui.actionShow_Station_Summary.setEnabled(enable)
        self._station_viewer.ui.pushButton_plot.setEnabled(enable)
        self.ui.actionPlot.setEnabled(enable)
        self.ui.actionExport_ModEM_Data.setEnabled(enable)
        self.ui.actionCreate_Shape_File_From_Stations.setEnabled(enable)
        self.ui.actionCreate_Phase_Tensor_csv_file.setEnabled(enable)
        self.ui.actionCreate_Measurement_csv_file.setEnabled(enable)

    def create_subwindow(self, widget, title, override=True, tooltip=None):
        subwindow = None
        self._subwindow_counter += 1
        if title in self.subwindows:
            if override:
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
        subwindow = StartGUI.MDISubWindow(self)
        subwindow.setWindowTitle(title)
        if tooltip:
            subwindow.setToolTip("<p>" + tooltip + "</p>")
        subwindow.setWidget(widget)
        subwindow.resize(widget.size())
        self.ui.mdiArea.addSubWindow(subwindow)

        # create menu action
        new_window_action = QAction(self)
        new_window_action.setObjectName("actionSubwindow%d" % self._subwindow_counter)
        new_window_action.setText(_translate("SmartMT_MainWindow", title, None))
        new_window_action.triggered.connect(subwindow.show_and_focus)

        # add to window menu
        self.ui.menuWindow.addAction(new_window_action)
        # add all references to self._subwindow
        self.subwindows[title] = (subwindow, new_window_action)
        subwindow.show()

        return subwindow, new_window_action

    class MDISubWindow(QMdiSubWindow):
        def __init__(self, main_ui, parent=None, flags=0):
            """

            :param main_ui:
            :type main_ui: StartGUI
            :param parent:
            :param flags:
            """
            super(QMdiSubWindow, self).__init__(parent)
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
    if __debug__:
        MtPyLog.load_configure(os.path.join(os.path.abspath("tests"), "logging.yml"))  # debug only
        # handle uncaught exceptions to log as since PYQT5.5 will not display any uncaught exceptions
        # ref: http://pyqt.sourceforge.net/Docs/PyQt5/incompatibilities.html#unhandled-python-exceptions
        logger = MtPyLog.get_mtpy_logger(__name__)
        sys.excepthook = lambda exc_type, exc_value, exc_trace: logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_trace))

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    smartMT = StartGUI()
    smartMT.show()

    # hack to fix the "python has stopped working" error,
    # the possible cause is the QtGui4.dll crashes, need to test it on linux environment
    # ref of the issue: http://pyqt.sourceforge.net/Docs/sip4/python_api.html
    #   "When the Python interpreter exits it garbage collects those objects that it can.
    #    This means that any corresponding C++ instances and C structures owned by Python
    #    are destroyed. Unfortunately this happens in an unpredictable order and so can
    #    cause memory faults within the wrapped library. Calling this function with a
    #    value of False disables the automatic destruction of C++ instances and C
    #    structures."
    sip.setdestroyonexit(False)
    # end of hack

    sys.exit(app.exec_())
