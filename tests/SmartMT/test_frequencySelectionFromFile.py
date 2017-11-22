import glob
from unittest import TestCase

import os

from qtpy import QtCore
from qtpy.QtTest import QTest
from qtpy.QtWidgets import QMainWindow, QWidget, QVBoxLayout

import numpy as np

from mtpy.core import mt
from mtpy.gui.SmartMT.Components.PlotParameter import FrequencySelectionFromFile
from tests import make_temp_dir
from tests.SmartMT import _click_area


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)

        widget = QWidget(self)
        layout = QVBoxLayout(widget)
        self.freq_slct_from_file = FrequencySelectionFromFile(self)
        layout.addWidget(self.freq_slct_from_file)
        widget.setLayout(layout)
        self.setCentralWidget(widget)


def _get_mt_objs(edi_path):
    edi_files = glob.glob(os.path.join(edi_path, '*.edi'))
    mt_objs = [mt.MT(os.path.abspath(file_name)) for file_name in edi_files]
    return mt_objs


class TestFrequencySelectionFromFile(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):
        # create gui
        self.app = MainWindow()
        self.app.show()
        self.app.freq_slct_from_file.set_data(_get_mt_objs("data/edifiles"))
        QTest.qWaitForWindowActive(self.app)

    def tearDown(self):
        self.app.close()

    # def test_manuel(self):
    #     QTest.qWait(10000)

    def test_selection(self):
        self.assertTrue(self.app.freq_slct_from_file.ui.tableWidget_selected.rowCount() == 0)
        _click_area(self.app.freq_slct_from_file.ui.listView_stations.viewport())
        QTest.qWait(1000)
        _click_area(self.app.freq_slct_from_file.ui.listView_stations.viewport(), modifier=QtCore.Qt.ControlModifier)
        QTest.qWait(1000)
        _click_area(self.app.freq_slct_from_file.ui.listView_stations.viewport(), modifier=QtCore.Qt.ControlModifier)
        QTest.qWait(1000)
        self.assertTrue(self.app.freq_slct_from_file.ui.tableWidget_selected.rowCount() > 0)

    def test_get_frequencies(self):
        _click_area(self.app.freq_slct_from_file.ui.listView_stations.viewport())
        QTest.qWait(1000)
        freqs = self.app.freq_slct_from_file.get_selected_frequencies()
        self.assertTrue(len(freqs) == self.app.freq_slct_from_file.ui.tableWidget_selected.rowCount())
        #  check frequency
        index = self.app.freq_slct_from_file.ui.listView_stations.selectedIndexes()[0]
        station = self.app.freq_slct_from_file.model_stations.item(index.row()).data(QtCore.Qt.DisplayRole)
        mtobj = self.app.freq_slct_from_file._mt_obj_dict[station]
        mt_freqs = sorted(list(set(mtobj.Z.freq)))
        self.assertTrue(all([freq == mtfreq for freq, mtfreq in zip(freqs, mt_freqs)]))

    def test_get_periods(self):
        _click_area(self.app.freq_slct_from_file.ui.listView_stations.viewport())
        QTest.qWait(1000)
        periods = self.app.freq_slct_from_file.get_selected_periods()
        self.assertTrue(len(periods) == self.app.freq_slct_from_file.ui.tableWidget_selected.rowCount())
        #  check periods
        index = self.app.freq_slct_from_file.ui.listView_stations.selectedIndexes()[0]
        station = self.app.freq_slct_from_file.model_stations.item(index.row()).data(QtCore.Qt.DisplayRole)
        mtobj = self.app.freq_slct_from_file._mt_obj_dict[station]
        mt_freqs = sorted(list(set(mtobj.Z.freq)))
        mt_periods = list(1./np.array(mt_freqs))
        self.assertTrue(all([freq == mtfreq for freq, mtfreq in zip(periods, mt_periods)]))
