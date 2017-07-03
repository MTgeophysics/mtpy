# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the station viewer subwindow

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui, QtCore
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib_imabedding import MPLCanvas
from station_viewer import Ui_StationViewer, _translate, _fromUtf8

WINDOW_TITLE = _translate("StationViewer","Stations Stats", None)
MAP_DISABLE_TEXT = _translate("StationViewer", "Show Map", None)
MAP_ENABLE_TEXT = _translate("StationViewer", "Hide Map", None)

class StationViewer(QtGui.QWidget):
    def __init__(self, parent, file_handler):
        """

        :param parent:
        :type parent: start.StartQt4
        :param file_handler
        :type file_handler.FileHandler
        """
        QtGui.QWidget.__init__(self, parent)
        self.file_handler = file_handler
        self.ui = Ui_StationViewer()
        self.ui.setupUi(self)
        self.subwindow = parent.create_subwindow(self, WINDOW_TITLE)
        # make station_viewer never been deleted
        self.subwindow = parent.subwindows[WINDOW_TITLE][0]
        self.subwindow.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
        # connect map button
        QtCore.QObject.connect(self.ui.pushButton_showMap, QtCore.SIGNAL("clicked()"), self.toggle_map)
        # add map area and hide by default
        self.fig_canvas = StationViewer.StationMap(self, self.file_handler)
        self.fig_canvas.setHidden(True)
        self.ui.verticalLayout.addWidget(self.fig_canvas)

    def toggle_map(self, *args, **kwargs):
        if self.ui.pushButton_showMap.text() == MAP_DISABLE_TEXT:
            self.fig_canvas.setMinimumHeight(self.height()/2)
            self.fig_canvas.update_figure()
            self.fig_canvas.setHidden(False)
            self.ui.pushButton_showMap.setText(MAP_ENABLE_TEXT)
        else:
            self.fig_canvas.setHidden(True)
            self.ui.pushButton_showMap.setText(MAP_DISABLE_TEXT)

    def update_view(self):
        self.ui.treeWidget_stations.clear()
        for group in self.file_handler.get_groups():
            node = QtGui.QTreeWidgetItem(self.ui.treeWidget_stations)
            node.setText(0, group)
            for ref in self.file_handler.get_group_members(group):
                mt_obj = self.file_handler.get_MT_obj(ref)
                child_node = QtGui.QTreeWidgetItem(node)
                child_node.setText(1, mt_obj.station)
                child_node.setText(2, ref)
        if not self.fig_canvas.isHidden():
            self.fig_canvas.update_figure()

    class StationMap(MPLCanvas):
        def __init__(self, parent=None, file_handler=None, width=5, hight=4, dpi=100):
            """

            :param parent:
            :param file_handler:
            :type file_handler: file_handler.FileHandler
            :param width:
            :param hight:
            :param dpi:
            """
            self.file_handler = file_handler
            MPLCanvas.__init__(self, parent, width, hight, dpi)

        def compute_initial_figure(self):
            self.update_figure()

        def update_figure(self):
            # clear figure
            self.axes.cla()
            # prepare data
            groups = []
            lat = []
            lon = []
            stations = []
            for group in self.file_handler.get_groups():
                for ref in self.file_handler.get_group_members(group):
                    mt_obj = self.file_handler.get_MT_obj(ref)
                    groups.append(group)
                    lat.append(mt_obj.lat)
                    lon.append(mt_obj.lon)
                    stations.append(mt_obj.station)
            df = pd.DataFrame(dict(x=lon, y=lat, group=groups, station=stations))
            print df.head()
            groups = df.groupby("group")

            # plot
            # self.axes = plt.subplots()
            self.axes.margins(0.05)   # 5% padding
            for name, group in groups:
                self.axes.plot(group.x, group.y, marker='o', linestyle='', ms=12, label=name)
            self.axes.legend(numpoints=1, loc="best")





