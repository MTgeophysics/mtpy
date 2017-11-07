# -*- coding: utf-8 -*-
"""
    Description:
        define the object to initialize the station summary subwindow

    Usage:
        python start.py

    Author: YingzhiGou
    Date: 20/06/2017
"""
from qtpy import QtCore
from qtpy.QtWidgets import QWidget

from mtpy.gui.SmartMT.ui_asset.station_status import Ui_StationStatus


class StationSummary(QWidget):
    def __init__(self, parent, file_handler, selected_files):
        """

        :param parent:
        :type parent: StartQt4
        :param file_handler:
        :type file_handler: FileHandler
        :param selected_files:
        :type selected_files: set
        """
        QWidget.__init__(self, parent)
        self.file_handler = file_handler
        self.selected_stations = selected_files
        self.ui = Ui_StationStatus()
        self.ui.setupUi(self)
        self.subwindow, _ = parent.create_subwindow(self, self.windowTitle())
        # self.subwindow.setMaximumWidth(600)
        # self.subwindow.setMinimumWidth(400)
        # self.subwindow.resize(parent.width()/3, self.height())
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)

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
        if mt_obj.station:
            self.ui.lineEditStation.setText(mt_obj.station)
        if mt_obj.fn:
            self.ui.lineEditFile_Ref.setText(mt_obj.fn)
        if mt_obj.edi_object.Header.acqby:
            self.ui.lineEditAcquired_by.setText(mt_obj.edi_object.Header.acqby)
        if mt_obj.edi_object.Header.loc:
            self.ui.lineEditLocation.setText(mt_obj.edi_object.Header.loc)
        if mt_obj.elev:
            self.ui.lineEditElev.setText('%.5f' % mt_obj.elev)
        if mt_obj.lat:
            self.ui.lineEditLat.setText('%.5f' % mt_obj.lat)
        if mt_obj.lon:
            self.ui.lineEditLong.setText('%.5f' % mt_obj.lon)
        if mt_obj.edi_object.Header.filedate:
            self.ui.lineEditDate_Acq.setText(mt_obj.edi_object.Header.filedate)
        with open(mt_obj.fn, 'r') as edi_file:
            self.ui.plainTextEdit_edi_text.setPlainText(edi_file.read())

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
        self.ui.plainTextEdit_edi_text.clear()

    def _update_station_summary(self, mt_objs):
        pass
