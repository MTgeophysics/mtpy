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
from itertools import cycle

from matplotlib_imabedding import MPLCanvas
from station_viewer import Ui_StationViewer, _translate, _fromUtf8

WINDOW_TITLE = _translate("StationViewer", "Stations Stats", None)
MAP_DISABLE_TEXT = _translate("StationViewer", "Show Map", None)
MAP_ENABLE_TEXT = _translate("StationViewer", "Hide Map", None)

MARKERS = ['*', 'D', 'H', '^']
COLORS = ['g', 'r', 'c', 'm', 'y', 'k', 'b']

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
        # handle selection changed event
        self.ui.treeWidget_stations.connect(self, QtCore.SIGNAL("itemSelectionChanged()"), self.item_selection_changed)
        # add context menu for tree view
        self.tree_menu = StationViewer.TreeViewMenu()
        self.tree_menu.actionCreate_New_Group.triggered.connect(self.create_new_group)
        self.tree_menu.actionDismiss_Group.triggered.connect(self.dismiss_group)
        self.tree_menu.actionAdd_To_Group.triggered.connect(self.add_selected_to_group)
        self.tree_menu.actionRemove_Station.triggered.connect(self.remove_stations)
        self.ui.treeWidget_stations.customContextMenuRequested.connect(self.open_menu)
        # connect map button
        QtCore.QObject.connect(self.ui.pushButton_showMap, QtCore.SIGNAL("clicked()"), self.toggle_map)
        # add map area and hide by default
        self.fig_canvas = StationViewer.StationMap(self, self.file_handler)
        self.fig_canvas.setHidden(True)
        self.ui.verticalLayout.addWidget(self.fig_canvas)

    def remove_stations(self, *args, **kwargs):
        selected = self.ui.treeWidget_stations.selectedItems()
        selected_stations = set()
        for item in selected:
            if not item.parent():
                # selected a group
                selected_group_id = str(item.text(0))
                members = self.file_handler.get_group_members(selected_group_id)
                if members:
                    for member in members:
                        selected_stations.add((self.file_handler.get_MT_obj(member), member))
            else:
                ref = str(item.text(1))
                selected_stations.add((self.file_handler.get_MT_obj(ref), ref))
        if selected_stations:
            reply = QtGui.QMessageBox.question(self, "Unload Selected Stations",
                                               "Are you sure you want to unload/remove the selected stations?\n(You can always load them back again.)",
                                               QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                for station, ref in selected_stations:
                    self.file_handler.unload(ref)
                self.update_view()

    def add_selected_to_group(self, *args, **kwargs):
        selected = self.ui.treeWidget_stations.selectedItems()
        if selected:
            groups = self.file_handler.get_groups()
            group_id, ok = QtGui.QInputDialog.getItem(self, "Add Selected Items to Group", "Please select one group:",
                                                      groups, 0, False)
            if ok and group_id:
                group_id = str(group_id)
                for item in selected:
                    if not item.parent():
                        # selected a group
                        selected_group_id = str(item.text(0))
                        members = self.file_handler.get_group_members(selected_group_id)
                        if members:
                            for member in members:
                                self.file_handler.add_to_group(group_id, member)
                    else:
                        # selected an item
                        ref = str(item.text(1))
                        self.file_handler.add_to_group(group_id, ref)
                self.update_view()

    def item_selection_changed(self, *args):
        # todo handle selection change
        selected = self.ui.treeWidget_stations.selectedItems()
        if selected:
            for item in selected:
                pass
        print "selection changed!!"
        pass

    def toggle_map(self, *args, **kwargs):
        if self.ui.pushButton_showMap.text() == MAP_DISABLE_TEXT:
            self.fig_canvas.setMinimumHeight(self.height() / 2)
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
                child_node.setText(0, mt_obj.station)
                child_node.setText(1, ref)
        if not self.fig_canvas.isHidden():
            self.fig_canvas.update_figure()

    def open_menu(self, position):
        items = self.ui.treeWidget_stations.selectedItems()
        if items:
            self.tree_menu.actionAdd_To_Group.setEnabled(True)
            self.tree_menu.actionRemove_Station.setEnabled(True)
            # if selected only a group
            if len(items) == 1 and not items[0].parent() and items[0].text(0) != "Default Group":
                # selected only one group
                self.tree_menu.actionDismiss_Group.setEnabled(True)
            else:
                self.tree_menu.actionDismiss_Group.setEnabled(False)
        else:
            self.tree_menu.actionAdd_To_Group.setEnabled(False)
            self.tree_menu.actionRemove_Station.setEnabled(False)
        self.tree_menu.exec_(self.ui.treeWidget_stations.viewport().mapToGlobal(position))

    def create_new_group(self, *args, **kwargs):
        ok = False
        while not ok:
            text, ok = QtGui.QInputDialog.getText(self, 'New Group', 'Enter a New Group Name:')
            if ok:
                text = str(text)
                ok = self.file_handler.create_group(text)
                if ok:
                    self.update_view()
                else:
                    QtGui.QMessageBox.information(self, "NOTE",
                                                  "Group %s already exits" % text)
            else:
                # cancelled
                break

    def dismiss_group(self, *args, **kwargs):
        """
        remove group and more ungrouped data to Default Group
        :return:
        """
        group_id = self.ui.treeWidget_stations.selectedItems()[0].text(0)
        group_id = str(group_id)
        self.file_handler.remove_group(group_id)
        self.update_view()

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
            df.station.astype(str)
            groups = df.groupby("group")
            # plot
            # self.axes = plt.subplots()
            self.axes.margins(0.05)  # 5% padding
            annotated_stations = set()
            for (name, group), marker, color in zip(list(groups), cycle(MARKERS), cycle(COLORS)):
                self.axes.plot(group.x, group.y, marker=marker, linestyle='', ms=10, alpha=.5, label=name,
                               picker=3)  # picker = 5 3 point tolerance
                for index, row in df.iterrows():
                    if row['station'] not in annotated_stations:
                        self.axes.annotate(row['station'], xy=(row['x'], row['y']), size=10)
                        annotated_stations.add(row['station'])
            self.axes.legend(numpoints=1, loc="best")
            self.axes.grid()

        def update_figure(self):
            self.compute_initial_figure()
            self.draw()

    class TreeViewMenu(QtGui.QMenu):
        def __init__(self, *args):
            QtGui.QMenu.__init__(self, *args)
            # create new group action
            self.actionCreate_New_Group = QtGui.QAction(self)
            self.actionCreate_New_Group.setEnabled(True)
            self.actionCreate_New_Group.setObjectName(_fromUtf8("actionCreate_New_Group"))
            self.addAction(self.actionCreate_New_Group)

            self.actionDismiss_Group = QtGui.QAction(self)
            self.actionDismiss_Group.setEnabled(True)
            self.actionDismiss_Group.setObjectName(_fromUtf8("actionDismiss_Group"))
            self.addAction(self.actionDismiss_Group)

            self.actionAdd_To_Group = QtGui.QAction(self)
            self.actionAdd_To_Group.setEnabled(True)
            self.actionAdd_To_Group.setObjectName(_fromUtf8("actionAdd_To_Group"))
            self.addAction(self.actionAdd_To_Group)

            self.addSeparator()
            # item operations
            self.actionRemove_Station = QtGui.QAction(self)
            self.actionRemove_Station.setEnabled(True)
            self.actionRemove_Station.setObjectName(_fromUtf8("actionRemove_Station"))
            self.addAction(self.actionRemove_Station)

            self.addSeparator()
            # plot menu action
            self.actionPlot = QtGui.QAction(self)
            self.actionPlot.setEnabled(False)
            self.actionPlot.setObjectName(_fromUtf8("actionPlot"))
            self.addAction(self.actionPlot)

            self.retranslateUi()

        def retranslateUi(self):
            self.actionCreate_New_Group.setText(_translate("StationViewer", "Create New Group...", None))
            self.actionDismiss_Group.setText(_translate("StationViewer", "Dismiss Selected Group", None))
            self.actionAdd_To_Group.setText(_translate("StationViewer", "Add Selected to Group...", None))
            self.actionRemove_Station.setText(_translate("StationViewer", "Unload Selected Stations...", None))
            self.actionPlot.setText(_translate("StationViewer", "Plot...", None))
