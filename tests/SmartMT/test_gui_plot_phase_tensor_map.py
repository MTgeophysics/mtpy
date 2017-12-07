from mtpy.gui.SmartMT.visualization import PhaseTensorMap
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotPenetrationDepth1D(SmartMTGUITestCase):
    def test_default(self):
        plot_gui = self._switch_to_plot(PhaseTensorMap)  # type: PhaseTensorMap

        # select random frequency by randomly click on the frequency selection gui
        # repeat a few time in case some of the random clicks are not valid
        for i in range(3):
            _click_area(plot_gui._frequency_ui.histogram,
                        offset=plot_gui._frequency_ui.histogram.geometry().topLeft())

        _click_area(plot_gui._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_gui._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        self._plot()
