import glob
import os
import sys
from unittest import TestCase

from PyQt4.QtGui import QApplication, QMainWindow, QWidget, QVBoxLayout
from PyQt4.QtTest import QTest

from mtpy.core import mt
from mtpy.gui.SmartMT.gui.plot_parameter_guis import FrequencySelect

app = QApplication(sys.argv)

edi_paths = [
    "",
    "tests\\data\\edifiles",
    "examples\\data\\edi2",
    "examples\\data\\edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests\\data\\edifiles2"
]


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)

        widget = QWidget(self)
        layout = QVBoxLayout(widget)
        self.frequency_select = FrequencySelect(self)
        layout.addWidget(self.frequency_select)
        widget.setLayout(layout)
        self.setCentralWidget(widget)


class TestFrequencySelect(TestCase):
    def setUp(self):
        # create gui
        self.app = MainWindow()
        self.app.show()
        QTest.qWaitForWindowShown(self.app)

    def test_default(self):
        edi_files = glob.glob(os.path.join(edi_paths[1], '*.edi'))
        mt_objs = [mt.MT(os.path.abspath(file_name)) for file_name in edi_files]

        self.app.frequency_select.set_data(mt_objs)

        app.exec_()
