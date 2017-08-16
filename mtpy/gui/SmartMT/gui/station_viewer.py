# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the station viewer subwindow

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
from itertools import cycle

import pandas as pd
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import pyqtSignal
from matplotlib import artist

from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MPLCanvas
from mtpy.gui.SmartMT.ui_asset.station_viewer import Ui_StationViewer, _translate, _fromUtf8

MARKERS = ['*', 'D', 'H', '^']
COLORS = ['g', 'r', 'c', 'm', 'y', 'k', 'b']


class StationViewer(QtGui.QWidget):
    """
    signal when the selected station changed
    """
    selection_changed = pyqtSignal()

    def __init__(self, parent, file_handler):
        """

        :param parent:
        :type parent: start.StartQt4
        :param file_handler
        :type file_handler.FileHandler
        """
        QtGui.QWidget.__init__(self, parent)
        # self._ignore_selection_change = False  # guard that used to ignore the selection change processing
        self.file_handler = file_handler
        self.ui = Ui_StationViewer()
        self.ui.setupUi(self)
        self.subwindow, _ = parent.create_subwindow(self, self.windowTitle())
        # self.subwindow.setMaximumWidth(600)
        # self.subwindow.setMinimumWidth(400)
        # self.subwindow.resize(parent.width()/3, self.height())
        # set the item to be ordered by ascending order
        self.ui.treeWidget_stations.sortByColumn(0, QtCore.Qt.AscendingOrder)
        # make station_viewer never been deleted
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
        # handle selection changed event
        self.ui.treeWidget_stations.selectionModel().selectionChanged.connect(self.item_selection_changed)
        # add context menu for tree view
        self.tree_menu = StationViewer.TreeViewMenu()
        self.tree_menu.actionCreate_New_Group.triggered.connect(self.create_new_group)
        self.tree_menu.actionDismiss_Group.triggered.connect(self.dismiss_group)
        self.tree_menu.actionAdd_To_Group.triggered.connect(self.add_selected_to_group)
        self.tree_menu.actionRemove_Station.triggered.connect(self.remove_stations)
        self.tree_menu.actionPlot.triggered.connect(parent.plot_selected_station)
        self.ui.treeWidget_stations.customContextMenuRequested.connect(self.open_menu_in_tree_view)
        # setup and connect map button
        self.ui.pushButton_hideMap.hide()
        self.ui.pushButton_hideMap.clicked.connect(self.toggle_map)
        self.ui.pushButton_showMap.clicked.connect(self.toggle_map)
        # add map area and hide by default
        self.fig_canvas = StationViewer.StationMap(self, self.file_handler)
        self.fig_canvas.setHidden(True)
        self.fig_canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.fig_canvas.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        self.fig_canvas.customContextMenuRequested.connect(self.open_menu_in_map_view)
        self.ui.verticalLayout.addWidget(self.fig_canvas)
        # share the selected station collection from fig_canvas
        self.selected_stations = self.fig_canvas.selected_stations
        # connect signals between map and tree view
        self.fig_canvas.selection_changed.connect(self.update_selection)
        self.selection_changed.connect(self.fig_canvas.update_figure)

    _selection_update_in_progress = False

    def remove_stations(self, *args, **kwargs):
        selected = self.ui.treeWidget_stations.selectedItems()
        # selected_stations = set()
        # for item in selected:
        #     if not item.parent():
        #         # selected a group
        #         selected_group_id = str(item.text(0))
        #         members = self.file_handler.get_group_members(selected_group_id)
        #         if members:
        #             for member in members:
        #                 selected_stations.add((self.file_handler.get_MT_obj(member), member))
        #     else:
        #         ref = str(item.text(1))
        #         selected_stations.add((self.file_handler.get_MT_obj(ref), ref))
        # just use the container from the figure
        selected_stations = self.selected_stations.copy()
        if selected_stations:
            reply = QtGui.QMessageBox.question(self, "Unload Selected Stations",
                                               "Are you sure you want to unload/remove the selected stations?\n(You can always load them back again.)",
                                               QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if reply == QtGui.QMessageBox.Yes:
                # for station, ref in selected_stations:
                #     self.file_handler.unload(ref)
                for station in selected_stations:
                    self.file_handler.unload(self.file_handler.station2ref(station))
                    self.selected_stations.remove(station)
                self.fig_canvas.selected_stations.clear()
                self.update_view()

                # update station summary/status
                self.selection_changed.emit()

    def add_selected_to_group(self, *args, **kwargs):
        # selected = self.ui.treeWidget_stations.selectedItems()
        # if selected:
        if self.selected_stations:
            groups = self.file_handler.get_groups()
            group_id, ok = QtGui.QInputDialog.getItem(self, "Add Selected Items to Group", "Please select one group:",
                                                      groups, 0, False)
            if ok and group_id:
                group_id = str(group_id)
                # for item in selected:
                for station in self.selected_stations:
                    # if not item.parent():
                    # selected a group
                    # selected_group_id = str(item.text(0))
                    # members = self.file_handler.get_group_members(selected_group_id)
                    # if members:
                    #     for member in members:
                    #         self.file_handler.add_to_group(group_id, member)
                    # else:
                    # selected an item
                    # ref = str(item.text(1))
                    # self.file_handler.add_to_group(group_id, ref)
                    self.file_handler.add_to_group(group_id, self.file_handler.station2ref(station))
                self.update_view()

    def update_selection(self):
        if self._selection_update_in_progress:
            return
        else:
            self._selection_update_in_progress = True
        root = self.ui.treeWidget_stations.invisibleRootItem()
        for i in range(root.childCount()):
            item = root.child(i)
            for j in range(item.childCount()):
                child = item.child(j)
                # print(str(child.text(0)))
                if str(child.text(0)) in self.selected_stations:
                    child.setSelected(True)
                else:
                    child.setSelected(False)
        self.ui.treeWidget_stations.updateGeometry()
        self.selection_changed.emit()
        self._selection_update_in_progress = False

    def item_selection_changed(self, *args):
        # print "selection changed"
        if self._selection_update_in_progress:
            return
        else:
            self._selection_update_in_progress = True

        self.selected_stations.clear()
        selected = self.ui.treeWidget_stations.selectedItems()
        if selected:
            for item in selected:
                if not item.parent():
                    # selected a group
                    selected_group_id = str(item.text(0))
                    members = self.file_handler.get_group_members(selected_group_id)
                    if members:
                        for member in members:
                            self.selected_stations.add(self.file_handler.get_MT_obj(member).station)
                else:
                    self.selected_stations.add(str(item.text(0)))
                    # root = self.ui.treeWidget_stations.invisibleRootItem()
                    # for i in range(root.childCount()):
                    #     group_item = root.child(i)
                    #     group_selected = group_item.isSelected()
                    #     for j in range(group_item.childCount()):
                    #         station_item = group_item.child(j)
                    #         if group_selected or station_item.isSelected():
                    #             # item selected or the group is selected
                    #             self.fig_canvas.selected_stations.add(str(station_item.text(0)))
                    #         else:
                    #             station = str(station_item.text(0))
                    #             if station in self.fig_canvas.selected_stations:
                    #                 self.fig_canvas.selected_stations.remove(station)
        self.selection_changed.emit()
        self._selection_update_in_progress = False

    def toggle_map(self, *args, **kwargs):
        if self.ui.pushButton_showMap.isEnabled():
            self.fig_canvas.setMinimumHeight(self.height() / 2)
            self.fig_canvas.setHidden(False)
            self.fig_canvas.update_figure()
            self.ui.pushButton_showMap.setEnabled(False)
            self.ui.pushButton_showMap.setHidden(True)
            self.ui.pushButton_hideMap.setEnabled(True)
            self.ui.pushButton_hideMap.setHidden(False)
        else:
            self.fig_canvas.setHidden(True)
            self.ui.pushButton_showMap.setEnabled(True)
            self.ui.pushButton_showMap.setHidden(False)
            self.ui.pushButton_hideMap.setEnabled(False)
            self.ui.pushButton_hideMap.setHidden(True)

    def update_view(self):
        # this part of the code is a hack - re-create tree every time
        # self.ui.treeWidget_stations.clear()
        # for group in self.file_handler.get_groups():
        #     node = QtGui.QTreeWidgetItem(self.ui.treeWidget_stations)
        #     node.setText(0, group)
        #     for ref in self.file_handler.get_group_members(group):
        #         mt_obj = self.file_handler.get_MT_obj(ref)
        #         child_node = QtGui.QTreeWidgetItem(node)
        #         child_node.setText(0, mt_obj.station)
        #         child_node.setText(1, ref)
        root = self.ui.treeWidget_stations.invisibleRootItem()
        groups = set()
        to_delete = []
        existing_groups = set(self.file_handler.get_groups())
        for i in range(root.childCount()):
            group_item = root.child(i)
            group_id = str(group_item.text(0))
            if group_id in existing_groups:
                groups.add(group_id)
                refs = set()
                # check station items
                for j in range(group_item.childCount()):
                    station_item = group_item.child(j)
                    station_ref = str(station_item.text(1))
                    mt_obj = self.file_handler.get_MT_obj(station_ref)
                    if mt_obj:
                        refs.add(station_ref)
                    else:
                        to_delete.append(station_item)
                # add new children
                for ref in self.file_handler.get_group_members(group_id):
                    if ref not in refs:
                        child_node = QtGui.QTreeWidgetItem(group_item)
                        child_node.setText(0, self.file_handler.get_MT_obj(ref).station)
                        child_node.setText(1, ref)
            else:
                to_delete.append(group_item)
        # delete all non existed items
        for item in to_delete:
            (item.parent() or root).removeChild(item)
        # add now groups
        for group_id in existing_groups:
            if group_id not in groups:
                # create new node
                node = QtGui.QTreeWidgetItem(self.ui.treeWidget_stations)
                node.setText(0, group_id)
                for ref in self.file_handler.get_group_members(group_id):
                    child_node = QtGui.QTreeWidgetItem(node)
                    child_node.setText(0, self.file_handler.get_MT_obj(ref).station)
                    child_node.setText(1, ref)
        # self.ui.treeWidget_stations.updateGeometry()
        self.selection_changed.emit()

    def open_menu_in_tree_view(self, position):
        self._update_menu_context()
        self.tree_menu.exec_(self.ui.treeWidget_stations.viewport().mapToGlobal(position))

    def _update_menu_context(self):
        items = self.ui.treeWidget_stations.selectedItems()
        if items:
            self.tree_menu.actionAdd_To_Group.setEnabled(True)
            self.tree_menu.actionRemove_Station.setEnabled(True)
            self.tree_menu.actionPlot.setEnabled(True)
            # if selected only a group
            if len(items) == 1 and not items[0].parent() and items[0].text(0) != "Default Group":
                # selected only one group
                self.tree_menu.actionDismiss_Group.setEnabled(True)
            else:
                self.tree_menu.actionDismiss_Group.setEnabled(False)
        else:
            self.tree_menu.actionAdd_To_Group.setEnabled(False)
            self.tree_menu.actionRemove_Station.setEnabled(False)
            self.tree_menu.actionPlot.setEnabled(False)

    def open_menu_in_map_view(self, position):
        self._update_menu_context()
        self.tree_menu.exec_(self.fig_canvas.mapToGlobal(position))

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
        """
        signal when the selected station changed
        """
        selection_changed = pyqtSignal()

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
            self.artists = dict()
            self._annotation_artists = dict()
            self.selected_stations = set()
            MPLCanvas.__init__(self, parent, width, hight, dpi)
            self.useblit = self.supports_blit
            self.mpl_connect('pick_event', self.map_pick)

        def compute_initial_figure(self):
            # clear figure
            self._axes.cla()
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
            self._axes.tick_params(axis='both', which='major', labelsize=8)
            self._axes.tick_params(axis='both', which='minor', labelsize=6)
            self._axes.margins(0.05)  # 5% padding
            annotated_stations = set()
            self.artists.clear()
            for (name, group), marker, color in zip(list(groups), cycle(MARKERS), cycle(COLORS)):
                artist, = self._axes.plot(group.x, group.y, marker=marker, linestyle='', ms=10, alpha=.5, label=name,
                                          picker=5)  # picker = 5  point tolerance
                self.artists[artist] = group.station.tolist()
                for index, row in df.iterrows():
                    station = row['station']
                    if station not in annotated_stations:
                        annotation_artist = self._axes.annotate(station, xy=(row['x'], row['y']), size=6, picker=3,
                                                                color='red' if station in self.selected_stations else 'black')
                        self._annotation_artists[station] = annotation_artist
                        annotated_stations.add(station)
            self._axes.legend(numpoints=1, loc="best", fontsize=8)
            self._axes.grid()

        def update_figure(self, redraw=True):
            if not self.isHidden():
                self.compute_initial_figure()
                self.draw()

        def map_pick(self, event):
            """
            handles the pick event on the station map
            :param event:
            :return:
            """
            # print event.mouseevent.key, event.mouseevent.button
            if event.artist in self.artists:
                if not event.mouseevent.key or event.mouseevent.key != 'control':
                    self.selected_stations.clear()

                stations = self.artists[event.artist]
                if len(event.ind) != 0:
                    for subplotnum, dataind in enumerate(event.ind):
                        self.selected_stations.add(stations[dataind])
                        # update tree view accordingly
                        # self._ignore_selection_change = True
                        # self._ignore_selection_change = False
                        # emit selection changed signal
                        # self.ui.treeWidget_stations.emit(QtCore.SIGNAL("selectionChanged()"))
                if self.useblit:
                    for station, annotation_artist in self._annotation_artists.iteritems():
                        artist.setp(annotation_artist, color='red' if station in self.selected_stations else 'black')
                        self._axes.draw_artist(annotation_artist)
                    self.blit(self._axes.bbox)
                else:
                    self.update_figure()
                # tell others that the map selection is changed
                self.selection_changed.emit()
            else:
                pass
            return True

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
