# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the station summary subwindow

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.ui_asset.station_status import Ui_StationStatus


class StationSummary(QtGui.QWidget):
    def __init__(self, parent, file_handler, selected_files):
        """

        :param parent:
        :param file_handler:
        :type file_handler: file_handler.FileHandler
        :param selected_files:
        :type selected_files: set
        """
        QtGui.QWidget.__init__(self, parent)
        self.file_handler = file_handler
        self.selected_stations = selected_files
        self.ui = Ui_StationStatus()
        self.ui.setupUi(self)
        self.subwindow, _ = parent.create_subwindow(self, self.windowTitle())
        self.subwindow.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)

    def update_view(self):
        mt_objs = []
        # create a list of mt_obj
        for station in self.selected_stations:
            ref = self.file_handler.station2ref(station)
            mt_obj = self.file_handler.get_MT_obj(ref)
            if mt_obj:
                mt_objs.append(mt_obj)
        if mt_objs:
            self._update_station_detail(mt_objs[0])
            self._update_station_summary(mt_objs)
        else:
            self._clear_station_detail()
            self._clear_station_summary()

    def _clear_station_summary(self):
        pass

    def _update_station_detail(self, mt_obj):
        self._clear_station_detail()
        # set text
        self.ui.lineEditStation.setText(mt_obj.station)
        self.ui.lineEditFile_Ref.setText(mt_obj.fn)
        self.ui.lineEditAcquired_by.setText(mt_obj.edi_object.Header.acqby)
        self.ui.lineEditLocation.setText(mt_obj.edi_object.Header.loc)
        if mt_obj.elev:
            self.ui.lineEditElev.setText('%.5f' % mt_obj.elev)
        if mt_obj.lat:
            self.ui.lineEditLat.setText('%.5f' % mt_obj.lat)
        if mt_obj.lon:
            self.ui.lineEditLong.setText('%.5f' % mt_obj.lon)
        if mt_obj.edi_object.Header.filedate:
            self.ui.lineEditDate_Acq.setText(mt_obj.edi_object.Header.filedate)

    def _clear_station_detail(self):
        # clear text
        self.ui.lineEditStation.clear()
        self.ui.lineEditFile_Ref.clear()
        self.ui.lineEditAcquired_by.clear()
        self.ui.lineEditLong.clear()
        self.ui.lineEditLat.clear()
        self.ui.lineEditElev.clear()
        self.ui.lineEditLocation.clear()
        self.ui.lineEditDate_Acq.clear()

    def _update_station_summary(self, mt_objs):
        pass
