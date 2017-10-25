import glob
from unittest import TestCase

import os

from qtpy.QtTest import QTest
from qtpy.QtWidgets import QMainWindow, QWidget, QVBoxLayout

from mtpy.core import mt
from mtpy.gui.SmartMT.Components.PlotParameter import FrequencySelectionFromFile


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
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    def setUp(self):
        # create gui
        self.app = MainWindow()
        self.app.show()
        self.app.freq_slct_from_file.set_data(_get_mt_objs("tests/data/edifiles"))
        QTest.qWaitForWindowActive(self.app)

    def tearDown(self):
        self.app.close()

    def test_manuel(self):
        QTest.qWait(100000)


