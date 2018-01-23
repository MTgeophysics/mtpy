import os

from mtpy.gui.SmartMT.visualization import MTResponse
from tests import TEST_MTPY_ROOT
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotMTResponse(SmartMTGUITestCase):
    def test_plot_mt_response_enable_type_1(self):
        plot_config = self._switch_to_plot(MTResponse)  # type:MTResponse
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.radioButton_1, pos=self._pos_check_box)

        self._plot()

    def test_plot_mt_response_enable_type_2(self):
        plot_config = self._switch_to_plot(MTResponse)  # type:MTResponse

        # config plot
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.radioButton_2, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot()

    def test_plot_mt_response_enable_type_3(self):
        plot_config = self._switch_to_plot(MTResponse)  # type:MTResponse

        # config plot
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.radioButton_3, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot()

    def test_github_issue_24(self):
        # this test recreate the bug reported in issue #24 on github
        # https://github.com/MTgeophysics/mtpy/issues/24
        data_path = os.path.join(TEST_MTPY_ROOT, "examples/data/edi2")
        # load data under mtpy/examples/data/edi2
        self._load_data([os.path.join(data_path, edi) for edi in os.listdir(data_path) if edi.endswith("edi")])
        # select top 5 stations
        self._select_data(5)
        # switch to visualise
        _click_area(self.smartMT._station_viewer.ui.pushButton_plot)  # trigger plot widget
        self.assertTrue(self.smartMT.ui.stackedWidget.currentIndex() == 1)
        # select MT response visualisation
        index = [i for i in range(self.smartMT._plot_option.ui.comboBoxSelect_Plot.count())
                 if self.smartMT._plot_option.ui.comboBoxSelect_Plot.itemText(i) == MTResponse.plot_name()]
        self.assertFalse(len(index) == 0, "plot type not found")
        self.assertFalse(len(index) > 1, "plot type name is not unique")
        self.smartMT._plot_option.ui.comboBoxSelect_Plot.setCurrentIndex(index[0])
        plot_config = self.smartMT._plot_option._current_plot
        self.assertTrue(isinstance(plot_config, MTResponse))

        # click visualise to start
        self._plot()