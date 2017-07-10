# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the plot option subwindow

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.ui_asset.plot_options import Ui_PlotOption
from mtpy.gui.SmartMT.visualization.visualization_base import VisualizationBase
from mtpy.utils.mtpylog import MtPyLog
# import all VisualizationBase subclasses here
# todo may need a better way of searching sublasses from unloaded files
from mtpy.gui.SmartMT.visualization.penetration_depth3d import PenetrationDepth3D


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
        self._logger = MtPyLog().get_mtpy_logger(__name__)
        self.file_handler = file_handler
        self.selected_stations = selected_files
        self._current_plot = None
        self.ui = Ui_PlotOption()
        self.ui.setupUi(self)
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

        # self.ui.comboBoxSelect_Plot.addItem("test")

        if VisualizationBase.__subclasses__():
            self.ui.comboBoxSelect_Plot.setEnabled(True)
            self.ui.comboBoxSelect_Plot.setCurrentIndex(0)
            self._selection_changed()
        else:
            self.ui.comboBoxSelect_Plot.setEnabled(False)

        self.ui.comboBoxSelect_Plot.currentIndexChanged.connect(self._selection_changed)

    def _selection_changed(self, *args, **kwargs):
        # print "selection changed"
        index = self.ui.comboBoxSelect_Plot.currentIndex()
        plot_option = self.plotOptions[index]
        if self._current_plot is not None:
            self.ui.verticalLayout.removeWidget(self._current_plot.parameter_ui)
        self._current_plot = plot_option(self)
        description = plot_option.plot_description()
        # print description
        self.ui.textEditPlot_Description.setText(description)
        # set parameter ui
        self.ui.verticalLayout.addWidget(self._current_plot.parameter_ui)
        self.update_ui()

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

