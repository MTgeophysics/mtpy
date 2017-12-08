from mtpy.gui.SmartMT.visualization import Strike
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotPhaseTensorPseudoSection(SmartMTGUITestCase):
    def test_plot_type_1(self):
        plot_gui = self._switch_to_plot(Strike)  # type: Strike

        _click_area(plot_gui._plot_control_ui.ui.radioButton_type_1, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control_ui.ui.checkBox_plot_tipper, pos=self._pos_check_box)
        _click_area(plot_gui._rotation_ui.ui.dial_rotation)

        self._plot()

        # fold
        _click_area(plot_gui._plot_control_ui.ui.checkBox_fold, pos=self._pos_check_box)
        self._plot()

    def test_plot_type_2(self):
        plot_gui = self._switch_to_plot(Strike)  # type: Strike

        # type 2
        _click_area(plot_gui._plot_control_ui.ui.radioButton_type_2, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control_ui.ui.checkBox_plot_tipper, pos=self._pos_check_box)
        _click_area(plot_gui._rotation_ui.ui.dial_rotation)

        self._plot()

        # fold
        _click_area(plot_gui._plot_control_ui.ui.checkBox_fold, pos=self._pos_check_box)
        self._plot()
