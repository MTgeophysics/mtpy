from __future__ import print_function

import glob
import os
import pprint
import sys
from unittest import TestCase

from PyQt4 import QtGui, QtCore
from PyQt4.QtGui import QApplication
from PyQt4.QtTest import QTest

from mtpy.core import mt
from mtpy.gui.SmartMT.gui.export_dialog_modem import ExportDialogModEm

app = QApplication(sys.argv)


def _fake_exec_accept():
    return QtGui.QFileDialog.Accepted


def _fake_exec_reject():
    return QtGui.QFileDialog.Rejected


def _rewrite_text(widget, text):
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_A, QtCore.Qt.ControlModifier)
    QTest.keyClicks(widget, text)
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_Enter)


edi_paths = [
    "",
    "tests\\data\\edifiles",
    "examples\\data\\edi2",
    "examples\\data\\edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests\\data\\edifiles2"
]


class TestExportDialogModEm(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = "tests/temp"
        if not os.path.isdir(cls._temp_dir):
            os.mkdir(cls._temp_dir)

    def setUp(self):
        # create gui
        self.dialog = ExportDialogModEm()
        self.dialog.show()
        QTest.qWaitForWindowShown(self.dialog)

    def tearDown(self):
        self.dialog.close()

    def test_defaults(self):
        edi_files = glob.glob(os.path.join(edi_paths[1], '*.edi'))
        mt_objs = [mt.MT(os.path.abspath(file_name)) for file_name in edi_files]
        self.dialog.set_data(mt_objs)
        if self.dialog.exec_() == QtGui.QWizard.Accepted:
            print(self.dialog.get_save_file_path())
            pprint.pprint(self.dialog.get_data_kwargs())
            pprint.pprint(self.dialog.get_model_kwargs())

            self.dialog.export_data()
