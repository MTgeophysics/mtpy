from mtpy.gui.SmartMT.visualization import Depth1D, Depth2D, Depth3D
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotPenetrationDepth1D(SmartMTGUITestCase):
    def test_all(self):
        plot_gui = self._switch_to_plot(Depth1D)  # type: Depth1D
        self._plot()

        _click_area(plot_gui._z_component_ui.ui.checkBox_det, self._pos_check_box)
        self._plot()
        _click_area(plot_gui._z_component_ui.ui.checkBox_zxy, self._pos_check_box)
        self._plot()

        # last checked checkbox should be disabled and not unchecked
        _click_area(plot_gui._z_component_ui.ui.checkBox_zyx, self._pos_check_box)
        self.assertFalse(plot_gui._z_component_ui.ui.checkBox_zyx.isEnabled())
        self.assertTrue(plot_gui._z_component_ui.ui.checkBox_zyx.isChecked())


class TestGUIPlotPenetrationDepth2D(SmartMTGUITestCase):
    def test_all(self):
        plot_gui = self._switch_to_plot(Depth2D)  # type: Depth2D
        _click_area(plot_gui._z_component_ui.ui.radioButton_det, self._pos_check_box)
        self._plot()
        _click_area(plot_gui._z_component_ui.ui.radioButton_zxy, self._pos_check_box)
        self._plot()
        _click_area(plot_gui._z_component_ui.ui.radioButton_zyx, self._pos_check_box)
        self._plot()

class TestGUIPlotPenetrationDepth2D(SmartMTGUITestCase):
    def test_all(self):
        plot_gui = self._switch_to_plot(Depth3D)  # type: Depth3D

        # select random frequency by randomly click on the frequency selection gui
        # repeat a few time in case some of the random clicks are not valid
        for i in range(3):
            _click_area(plot_gui._frequency_period_ui.histogram, offset=plot_gui._frequency_period_ui.histogram.geometry().topLeft())

        _click_area(plot_gui._z_component_ui.ui.radioButton_det, self._pos_check_box)
        _click_area(plot_gui._z_unit_ui.ui.radioButton_m)
        self._plot()
        _click_area(plot_gui._z_unit_ui.ui.radioButton_km)
        _click_area(plot_gui._z_component_ui.ui.radioButton_zxy, self._pos_check_box)
        self._plot()
        _click_area(plot_gui._z_component_ui.ui.radioButton_zyx, self._pos_check_box)
        self._plot()
