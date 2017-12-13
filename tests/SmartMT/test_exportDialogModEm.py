from __future__ import print_function

import glob
import os
import pprint
from unittest import TestCase

import pytest
from qtpy import QtCore
from qtpy.QtWidgets import QFileDialog, QWizard
from qtpy.QtTest import QTest

from mtpy.core.mt import MT
from mtpy.gui.SmartMT.gui.export_dialog_modem import ExportDialogModEm
from tests import make_temp_dir, AUS_TOPO_FILE


def _fake_exec_accept():
    return QFileDialog.Accepted


def _fake_exec_reject():
    return QFileDialog.Rejected


def _rewrite_text(widget, text):
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_A, QtCore.Qt.ControlModifier)
    QTest.keyClicks(widget, text)
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_Enter)


edi_paths = [
    "data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM",
    "../MT_Datasets/GA_UA_edited_10s-10000s",
    "data/edifiles2"
]


@pytest.mark.skip("Not yet implemented")
class TestExportDialogModEm(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = make_temp_dir(cls.__name__)

    def setUp(self):
        # create gui
        self.dialog = ExportDialogModEm()
        self.dialog.show()
        QTest.qWaitForWindowActive(self.dialog)

    def tearDown(self):
        self.dialog.close()

    def test_defaults(self):
        edi_files = glob.glob(os.path.join(edi_paths[0], '*.edi'))
        mt_objs = [MT(os.path.abspath(file_name)) for file_name in edi_files]
        self.dialog.set_data(mt_objs)
        _rewrite_text(self.dialog.ui.comboBox_topography_file, AUS_TOPO_FILE)
        self.dialog.exec_ = _fake_exec_accept
        if self.dialog.exec_() == QWizard.Accepted:
            print(self.dialog.get_save_file_path())
            pprint.pprint(self.dialog.get_data_kwargs())
            pprint.pprint(self.dialog.get_model_kwargs())

            self.dialog.export_data()

        self.dialog.close()
