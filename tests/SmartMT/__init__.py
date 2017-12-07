from __future__ import print_function

import os
from unittest import TestCase

import pytest

from mtpy.gui.SmartMT.start import StartGUI
from mtpy.gui.SmartMT.visualization import VisualizationBase
from tests import EDI_DATA_DIR

qtpy = pytest.importorskip("qtpy")
pytestmark = pytest.mark.skipif(os.name == "posix" and 'DISPLAY' not in os.environ)

import random
import sys

import matplotlib
import sip
from qtpy import QtCore
from qtpy.QtTest import QTest
from qtpy.QtWidgets import QApplication

from mtpy.utils.mtpylog import MtPyLog

sip.setdestroyonexit(False)

app = QApplication(sys.argv)

import matplotlib.pyplot as plt
plt.ion()

MtPyLog.get_mtpy_logger(__name__).info("Testing using matplotlib backend {}".format(matplotlib.rcParams['backend']))

# handle uncaught exceptions to log as since PYQT5.5 will not display any uncaught exceptions
# ref: http://pyqt.sourceforge.net/Docs/PyQt5/incompatibilities.html#unhandled-python-exceptions
logger = MtPyLog.get_mtpy_logger(__name__)
sys.excepthook = lambda exc_type, exc_value, exc_trace: logger.error("Uncaught exception", exc_info=(exc_type, exc_value, exc_trace))

def _click_area(qobj, pos=None, offset=None, modifier=QtCore.Qt.NoModifier, timeout=100):
    geom = qobj.geometry()
    if pos is None:
        x = int(geom.width() * random.uniform(0.2, 0.8))  # avid to click on the edge of widgets
        y = int(geom.height() * random.uniform(0.2, 0.8))
        pos = QtCore.QPoint(x, y)
    if offset is not None:
        pos += offset
    QTest.mouseClick(qobj, QtCore.Qt.LeftButton, modifier=modifier, pos=pos)
    # print(pos.x(), pos.y())
    QTest.qWait(timeout)  # wait for things to happen on gui


class SmartMTGUITestCase(TestCase):
    _pos_check_box = QtCore.QPoint(8, 8)

    @classmethod
    def setUpClass(cls):
        import matplotlib.pyplot as plt
        plt.interactive(False)


    def setUp(self):
        self.smartMT = StartGUI()
        self.smartMT.show()
        QTest.qWaitForWindowActive(self.smartMT)
        print(matplotlib.get_backend())

    def _switch_to_plot(self, plot_type=VisualizationBase):
        name = plot_type.plot_name()
        # load some data
        self._load_data()
        # set the loaded data to be plotted
        self.smartMT._station_viewer.ui.treeWidget_stations.selectAll()
        self.smartMT._station_viewer.item_selection_changed()
        _click_area(self.smartMT._station_viewer.ui.pushButton_plot)  # trigger plot widget
        self.assertTrue(self.smartMT.ui.stackedWidget.currentIndex() == 1)
        # switch to the plot
        index = [i for i in range(self.smartMT._plot_option.ui.comboBoxSelect_Plot.count())
                 if self.smartMT._plot_option.ui.comboBoxSelect_Plot.itemText(i) == name]
        self.assertFalse(len(index) == 0, "plot type not found")
        self.assertFalse(len(index) > 1, "plot type name is not unique")
        self.smartMT._plot_option.ui.comboBoxSelect_Plot.setCurrentIndex(index[0])

        plot_config = self.smartMT._plot_option._current_plot
        self.assertTrue(isinstance(plot_config, plot_type))
        return plot_config

    def tearDown(self):
        self.smartMT.close()

    def _load_data(self):
        file_list = [os.path.join(EDI_DATA_DIR, edi) for edi in os.listdir(EDI_DATA_DIR) if edi.endswith("edi")]
        self.smartMT._progress_bar.setMaximumValue(len(file_list))
        self.smartMT._progress_bar.onStart()
        self.smartMT._add_files(file_list, os.path.basename(EDI_DATA_DIR))
        self.smartMT._update_tree_view()
        self.smartMT._progress_bar.onFinished()

    def _plot(self, timeout=3000):
        subwindow_counter = self.smartMT._subwindow_counter
        loop = QtCore.QEventLoop()

        self.smartMT._plot_option._current_plot.plotting_completed.connect(loop.quit)
        self.smartMT._plot_option._current_plot.plotting_error.connect(loop.quit)

        def handleTimeout():
            # timed out, stop loop
            if loop.isRunning():
                loop.quit()
                self.fail("GUI plotting timed out, maybe consider increasing timeout to wait longer")
        _click_area(self.smartMT._plot_option.ui.pushButton_plot)

        if timeout is not None:
            QtCore.QTimer.singleShot(timeout, handleTimeout)
        loop.exec_()  # wait for plotting

        self.assertTrue(self.smartMT._subwindow_counter == subwindow_counter + 1,
                        "no image created")  # test if the image is created
                        # "no image created, maybe consider to increase the timeout value")  # test if the image is created
