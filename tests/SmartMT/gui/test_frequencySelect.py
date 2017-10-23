import glob
import os
from unittest import TestCase

from qtpy.QtCore import QPoint
from qtpy.QtWidgets import QMainWindow, QWidget, QVBoxLayout
from PyQt4.QtTest import QTest

from mtpy.core import mt
from mtpy.gui.SmartMT.gui.plot_parameter_guis import FrequencySelect
from tests.SmartMT.gui import _click_area

edi_paths = [
    "tests/data/edifiles",
    "examples/data/edi2",
    "examples/data/edi_files",
    "../MT_Datasets/3D_MT_data_edited_fromDuanJM/",
    "../MT_Datasets/GA_UA_edited_10s-10000s/",
    "tests/data/edifiles2"
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


def _get_mt_objs(edi_path):
    edi_files = glob.glob(os.path.join(edi_path, '*.edi'))
    mt_objs = [mt.MT(os.path.abspath(file_name)) for file_name in edi_files]
    return mt_objs


class TestFrequencySelect(TestCase):
    def setUp(self):
        # create gui
        self.app = MainWindow()
        self.app.show()
        QTest.qWaitForWindowShown(self.app)


    def _std_function_tests(self):
        pos_check_box = QPoint(8, 8)
        _click_area(self.app.frequency_select.ui.radioButton_period, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.radioButton_period.isChecked())

        _click_area(self.app.frequency_select.ui.radioButton_frequency, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.radioButton_frequency.isChecked())

        # test frequency selection
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        QTest.qWait(1000)
        self.assertTrue(self.app.frequency_select.model_selected.rowCount() > 0)

        _click_area(self.app.frequency_select.ui.checkBox_existing_only, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.checkBox_existing_only.isChecked())
        self.assertTrue(self.app.frequency_select.histogram._select_existing_only)
        self.assertTrue(self.app.frequency_select.model_selected.rowCount() == 0)

        _click_area(self.app.frequency_select.ui.checkBox_show_existing, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.checkBox_show_existing.isChecked())
        self.assertTrue(self.app.frequency_select.histogram._show_existing)

        _click_area(self.app.frequency_select.ui.checkBox_y_log_scale, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.checkBox_y_log_scale.isChecked())
        self.assertTrue(self.app.frequency_select.histogram._y_log_scale)

        _click_area(self.app.frequency_select.ui.checkBox_x_log_scale, pos_check_box)
        QTest.qWait(2000)
        self.assertTrue(self.app.frequency_select.ui.checkBox_x_log_scale.isChecked())
        self.assertTrue(self.app.frequency_select.histogram._x_log_scale)

        # test clear
        _click_area(self.app.frequency_select.ui.pushButton_clear)
        QTest.qWait(1000)
        self.assertTrue(self.app.frequency_select.model_selected.rowCount() == 0)
        # test delete
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        _click_area(self.app.frequency_select.histogram, offset=self.app.frequency_select.histogram.geometry().topLeft())
        QTest.qWait(1000)
        self.assertTrue(self.app.frequency_select.model_selected.rowCount() > 0)
        self.app.frequency_select.ui.listView_selected.selectAll()
        _click_area(self.app.frequency_select.ui.pushButton_delete)
        QTest.qWait(1000)
        self.assertTrue(self.app.frequency_select.model_selected.rowCount() == 0)


def _generate_tests(edi_path):
    def _test_case(self):
        if os.path.exists(edi_path):
            mt_objs = _get_mt_objs(edi_path)
            self.app.frequency_select.set_data(mt_objs)
            QTest.qWait(2000)
            self._std_function_tests()
        else:
            self.app.close()
            self.skipTest("edi path not found")

    return _test_case


for edi_path in edi_paths:
    _test_case = _generate_tests(edi_path)
    _test_case.__name__ = "test_case_{}".format(os.path.basename(edi_path))
    setattr(TestFrequencySelect, _test_case.__name__, _test_case)
