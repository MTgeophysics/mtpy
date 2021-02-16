from mtpy.gui.SmartMT.visualization import (
    PhaseTensorMap,
    PhaseTensorPseudoSection,
    ResistivityPhasePseudoSection,
)
from tests.SmartMT import SmartMTGUITestCase, _click_area


class TestGUIPlotPhaseTensorMap(SmartMTGUITestCase):
    def test_default(self):
        print("TestGUIPlotPhaseTensorMap.....")
        plot_gui = self._switch_to_plot(PhaseTensorMap)  # type: PhaseTensorMap

        # select random frequency by randomly click on the frequency selection gui
        # repeat a few time in case some of the random clicks are not valid

        for i in range(3):
            _click_area(
                plot_gui._frequency_ui.histogram,
                offset=plot_gui._frequency_ui.histogram.geometry().topLeft(),
            )

        _click_area(plot_gui._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        _click_area(plot_gui._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        self._plot()


class TestGUIPlotPhaseTensorPseudoSection(SmartMTGUITestCase):
    def test_default(self):
        print("TestGUIPlotPhaseTensorPseudoSection.....")

        plot_gui = self._switch_to_plot(
            PhaseTensorPseudoSection
        )  # type: PhaseTensorPseudoSection

        _click_area(plot_gui._plot_control.ui.checkBox_zxx, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_zyy, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_period, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_phase, pos=self._pos_check_box)
        _click_area(
            plot_gui._plot_control.ui.checkBox_resistivity, pos=self._pos_check_box
        )
        _click_area(plot_gui._arrow_ui.ui.checkBox_real, pos=self._pos_check_box)
        _click_area(plot_gui._arrow_ui.ui.checkBox_imaginary, pos=self._pos_check_box)
        self._plot(20000)


class TestGUIPlotResistivityPhasePseudoSection(SmartMTGUITestCase):
    def test_grid_imshow(self):
        print("TestGUIPlotResistivityPhasePseudoSection.....")
        plot_gui = self._switch_to_plot(
            ResistivityPhasePseudoSection
        )  # type: ResistivityPhasePseudoSection

        _click_area(plot_gui._plot_control.ui.checkBox_zxx, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_zyy, pos=self._pos_check_box)
        # _click_area(plot_gui._plot_control.ui.checkBox_period, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_phase, pos=self._pos_check_box)
        _click_area(
            plot_gui._plot_control.ui.checkBox_resistivity, pos=self._pos_check_box
        )

        _click_area(
            plot_gui._mesh_grid_ui.ui.radioButton_imshow, pos=self._pos_check_box
        )

        self._plot(10000)

    def test_grid_pcolormesh(self):
        print("test_grid_pcolormesh.....")
        plot_gui = self._switch_to_plot(
            ResistivityPhasePseudoSection
        )  # type: ResistivityPhasePseudoSection

        _click_area(plot_gui._plot_control.ui.checkBox_zxx, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_zyy, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_period, pos=self._pos_check_box)
        _click_area(plot_gui._plot_control.ui.checkBox_phase, pos=self._pos_check_box)
        _click_area(
            plot_gui._plot_control.ui.checkBox_resistivity, pos=self._pos_check_box
        )

        _click_area(
            plot_gui._mesh_grid_ui.ui.radioButton_pcolormesh, pos=self._pos_check_box
        )

        self._plot(10000)
