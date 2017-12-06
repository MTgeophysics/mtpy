from mtpy.gui.SmartMT.visualization import MultipleMTResponses
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotMultipleMTResponse(SmartMTGUITestCase):
    def test_multiple_mt_response_compare_type_1(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_1, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot

    def test_multiple_mt_response_compare_type_2(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_2, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot

    def test_multiple_mt_response_compare_type_3(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_3, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot


    def test_multiple_mt_response_all_type_1(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_all)
        _click_area(plot_config._plot_control_ui.ui.radioButton_1, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot

    def test_multiple_mt_response_all_type_2(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_all)
        _click_area(plot_config._plot_control_ui.ui.radioButton_2, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot

    def test_multiple_mt_response_all_type_3(self):
        self._switch_to_plot("Multiple MT responses")
        plot_config = self.smartMT._plot_option._current_plot  # type: MultipleMTResponses
        _click_area(plot_config._plot_control_ui.ui.radioButton_all)
        _click_area(plot_config._plot_control_ui.ui.radioButton_3, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_skew, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_pt, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_i, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_p, pos=self._pos_check_box)
        _click_area(plot_config._plot_control_ui.ui.checkBox_strike_t, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_config._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_config._rotation_ui.ui.dial_rotation)

        self._plot(60000)  # this complete image could take very long time to plot
