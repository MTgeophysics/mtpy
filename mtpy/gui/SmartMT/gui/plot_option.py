# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the plot option subwindow

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui

from mtpy.gui.SmartMT.gui.busy_indicators import BusyOverlay
from mtpy.gui.SmartMT.ui_asset.plot_options import Ui_PlotOption
from mtpy.gui.SmartMT.visualization import *
from mtpy.gui.SmartMT.visualization.visualization_base import MPLCanvasWidget
from mtpy.utils.mtpylog import MtPyLog


class PlotOption(QtGui.QWidget):
    def __init__(self, parent, file_handler, selected_files):
        """

        :param parent:
        :type parent: StartQt4
        :param file_handler:
        :type file_handler: FileHandler
        :param selected_files:
        :type selected_files: set
        """
        QtGui.QWidget.__init__(self, parent)
        self._parent = parent
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self.file_handler = file_handler
        self.selected_stations = selected_files
        self._current_plot = None
        self.ui = Ui_PlotOption()
        self.ui.setupUi(self)

        # hide cancel button
        self.ui.pushButton_cancel.hide()

        # populate dropdown menu
        self.plotOptions = []

        # print VisualizationBase.__subclasses__()

        for child in VisualizationBase.__subclasses__():
            name = child.plot_name()
            if name not in self.plotOptions:
                self.plotOptions.append(child)
                self.ui.comboBoxSelect_Plot.addItem(name)
            else:
                raise Exception("Duplicated Plot Name: %s in class %s" % (name, child.__name__))

        # busy overlay
        self._busy_overlay = BusyOverlay(self)
        self._busy_overlay.hide()

        # connect signals
        self.ui.comboBoxSelect_Plot.currentIndexChanged.connect(self._selection_changed)
        self.ui.pushButton_plot.clicked.connect(self._create_plot)
        self.ui.pushButton_cancel.clicked.connect(self._cancel_plot)

        if VisualizationBase.__subclasses__():
            self.ui.comboBoxSelect_Plot.setEnabled(True)
            self.ui.comboBoxSelect_Plot.setCurrentIndex(0)
            self.ui.comboBoxSelect_Plot.currentIndexChanged.emit(0)
        else:
            self.ui.comboBoxSelect_Plot.setEnabled(False)

    def resizeEvent(self, event):
        size = event.size()
        size.setHeight(size.height() - self.ui.pushButton_plot.height())  # give space to the buttons
        self._busy_overlay.resize(size)
        # self._busy_overlay.resize(event.size())
        event.accept()

    def _selection_changed(self, *args, **kwargs):
        # print "selection changed"
        index = self.ui.comboBoxSelect_Plot.currentIndex()
        plot_option = self.plotOptions[index]
        description = plot_option.plot_description()
        # print description
        self.ui.textEditPlot_Description.setText(description)

        # set parameter ui
        if self._current_plot is not None:
            # self.ui.verticalLayout.removeWidget(self._current_plot.parameter_ui)
            self._current_plot.parameter_ui.deleteLater()

        self._current_plot = plot_option(self)
        # connect signal
        self._current_plot.started.connect(self._busy_overlay.show)
        self._current_plot.started.connect(self._tuggle_plot_cancel)
        self._current_plot.finished.connect(self._busy_overlay.hide)
        self._current_plot.finished.connect(self._tuggle_plot_cancel)
        self._current_plot.plotting_completed.connect(self._show_plot)
        self._current_plot.plotting_error.connect(self._plotting_error)

        self.ui.verticalLayout.addWidget(self._current_plot.parameter_ui)

        # self.resize(self.width(), self.sizeHint().height())
        self.update_ui()

    def _create_plot(self):
        self._current_plot.start(QtCore.QThread.HighPriority)

    def _cancel_plot(self):
        # self._current_plot.terminate()  # this does not work
        self._current_plot.wait()

    def _tuggle_plot_cancel(self):
        self.ui.pushButton_cancel.setHidden(not self.ui.pushButton_cancel.isHidden())
        self.ui.pushButton_plot.setHidden(not self.ui.pushButton_plot.isHidden())

    def _plotting_error(self, msg):
        QtGui.QMessageBox.critical(self,
                                   'Plotting Error', msg,
                                   QtGui.QMessageBox.Close)

    def _show_plot(self):
        fig = self._current_plot.get_fig()
        if fig:
            # self._fig.show()
            widget = MPLCanvasWidget(fig)
            self._parent.create_subwindow(widget, "%s" % self._current_plot.plot_name(), overide=False,
                                          tooltip=self._current_plot.get_plot_tooltip())

    def update_ui(self):
        if self._current_plot is not None:
            self._current_plot.set_data(self.get_mt_objs())

    def get_mt_objs(self):
        mt_objs = []
        for selected_station in self.selected_stations:
            ref = self.file_handler.station2ref(selected_station)
            mt_obj = self.file_handler.get_MT_obj(ref)
            mt_objs.append(mt_obj)
        return mt_objs
