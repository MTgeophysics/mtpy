import os
from unittest import TestCase

import matplotlib
import pytest
from qtpy import QtCore
from qtpy.QtCore import QPoint
from qtpy.QtTest import QTest

from mtpy.gui.SmartMT.start import StartGUI
from tests import EDI_DATA_DIR
from tests.SmartMT import _click_area

_pos_check_box = QPoint(8, 8)


@pytest.mark.last
class TestStartGUI(TestCase):
    """
    only testing the loading of the gui at the moment. to make sure all the necessory packages and dependencies are
    correctly loaded
    """
    @classmethod
    def setUpClass(cls):
        import matplotlib.pyplot as plt
        plt.interactive(False)

    def setUp(self):
        self.smartMT = StartGUI()
        self.smartMT.show()
        QTest.qWaitForWindowActive(self.smartMT)
        print(matplotlib.get_backend())

    def test_main(self):
        self.assertTrue(self.smartMT.isVisible())
        self.assertTrue(self.smartMT.isMaximized())

    def test_plot_mt_response_default(self):
        self._switch_to_plot("MT Response")
        self._plot()

    def test_plot_mt_response_enable_all(self):
        self._switch_to_plot("MT Response")

        # config plot
        plot_config = self.smartMT._plot_option._current_plot
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=_pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=_pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=_pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation, pos=_pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.radioButton_2, pos=_pos_check_box)

        self._plot()
        self.assertTrue(self.smartMT._subwindow_counter == 1, "no image creataed")  # test if the image is created

    def test_multiple_mt_response(self):
        self._switch_to_plot("Multiple MT responses")
        self._plot()

    def _switch_to_plot(self, name):
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
        self.assertTrue(len(index) == 1, "plot type name is not unique")
        self.smartMT._plot_option.ui.comboBoxSelect_Plot.setCurrentIndex(index[0])

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
        if timeout is not None:
            QtCore.QTimer.singleShot(timeout, loop.quit)
        _click_area(self.smartMT._plot_option.ui.pushButton_plot)  # wait for plotting
        loop.exec_()

        self.assertTrue(self.smartMT._subwindow_counter == subwindow_counter + 1, "no image created")  # test if the image is created

