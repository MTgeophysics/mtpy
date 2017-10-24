# -*- coding: utf-8 -*-
"""
Description:
This file contains GUI parts for plot parameter configuration that cannot be implemented to a generalized modules (in z_unit.py)
"""
import numpy as np
from qtpy.QtWidgets import QGroupBox
from qtpy.QtGui import QDoubleValidator


from mtpy.gui.SmartMT.ui_asset.groupbox_plot_control_mt_response import Ui_GroupBox_plot_control_mt_response
from mtpy.gui.SmartMT.ui_asset.groupbox_plot_control_resistivity_phase_pseudo_section import \
    Ui_GroupBox_plot_control_resistivity_phase_pseudo_section
from mtpy.gui.SmartMT.ui_asset.groupbox_plot_control_strike import Ui_GroupBox_plot_control_strike
from mtpy.imaging.mtcolors import cmapdict


class PlotControlMTResponse(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_plot_control_mt_response()
        self.ui.setupUi(self)

    def get_plot_num(self):
        if self.ui.radioButton_1.isChecked():
            return 1
        elif self.ui.radioButton_2.isChecked():
            return 2
        elif self.ui.radioButton_3.isChecked():
            return 3
        else:
            return 0  # should never reach here

    def get_strike(self):
        strike = "y"
        if self.ui.checkBox_strike_t.isChecked():
            strike += 't'
        if self.ui.checkBox_strike_p.isChecked():
            strike += 'p'
        if self.ui.checkBox_strike_i.isChecked():
            strike += 'i'
        if len(strike) > 1:
            return strike
        else:
            return 'n'

    def get_skew(self):
        if self.ui.radioButton_skew_y.isChecked():
            return 'y'
        else:
            return 'n'

    def get_ellipses(self):
        if self.ui.radioButton_ellipses_y.isChecked():
            return 'y'
        else:
            return 'n'

    def get_style(self):
        if self.ui.radioButton_compare.isChecked():
            return "compare"
        else:
            return "all"

    def hide_plot_style(self):
        self.ui.groupBox_plot_style.hide()


class PlotControlResistivityPhasePseudoSection(QGroupBox):
    """
    plot settings for resistivity phase pseudo section that cannot be standardized
    """

    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_plot_control_resistivity_phase_pseudo_section()
        # set up gui
        self.ui.setupUi(self)

        self.ui.comboBox_res_cmap.addItems(cmapdict.keys())
        self.ui.comboBox_res_cmap.setCurrentIndex(self.ui.comboBox_res_cmap.findText('mt_rd2gr2bl'))
        self.ui.comboBox_phase_cmap.addItems(cmapdict.keys())
        self.ui.comboBox_phase_cmap.setCurrentIndex(self.ui.comboBox_phase_cmap.findText('mt_bl2gr2rd'))

        self.ui.doubleSpinBox_res_max.setMaximum(np.inf)
        self.ui.doubleSpinBox_res_min.setMaximum(np.inf)
        self.ui.doubleSpinBox_period_min.setMaximum(np.inf)
        self.ui.doubleSpinBox_period_max.setMaximum(np.inf)

        # connect signals
        self.ui.doubleSpinBox_period_min.valueChanged.connect(self._period_min_changed)
        self.ui.doubleSpinBox_period_max.valueChanged.connect(self._period_max_changed)
        self.ui.doubleSpinBox_phase_min.valueChanged.connect(self._phase_min_changed)
        self.ui.doubleSpinBox_phase_max.valueChanged.connect(self._phase_max_changed)
        self.ui.doubleSpinBox_res_min.valueChanged.connect(self._res_min_changed)
        self.ui.doubleSpinBox_res_max.valueChanged.connect(self._res_max_changed)

        self.ui.checkBox_phase.stateChanged.connect(self._phase_state_changed)
        self.ui.checkBox_period.stateChanged.connect(self._period_state_changed)
        self.ui.checkBox_resistivity.stateChanged.connect(self._res_state_changed)

    def _res_min_changed(self, value):
        if value > self.ui.doubleSpinBox_res_max.value():
            self.ui.doubleSpinBox_res_max.blockSignals(True)
            self.ui.doubleSpinBox_res_max.setValue(value)
            self.ui.doubleSpinBox_res_max.blockSignals(False)

    def _res_max_changed(self, value):
        if value < self.ui.doubleSpinBox_res_min.value():
            self.ui.doubleSpinBox_res_min.blockSignals(True)
            self.ui.doubleSpinBox_res_min.setValue(value)
            self.ui.doubleSpinBox_res_min.blockSignals(False)

    def _res_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_res_min.setEnabled(state)
        self.ui.doubleSpinBox_res_max.setEnabled(state)

    def _phase_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_phase_min.setEnabled(state)
        self.ui.doubleSpinBox_phase_max.setEnabled(state)

    def _period_state_changed(self, p_int):
        state = bool(p_int != 0)
        self.ui.doubleSpinBox_period_min.setEnabled(state)
        self.ui.doubleSpinBox_period_max.setEnabled(state)

    def _period_min_changed(self, value):
        if value > self.ui.doubleSpinBox_period_max.value():
            self.ui.doubleSpinBox_period_max.blockSignals(True)
            self.ui.doubleSpinBox_period_max.setValue(value)
            self.ui.doubleSpinBox_period_max.blockSignals(False)

    def _period_max_changed(self, value):
        if value < self.ui.doubleSpinBox_period_min.value():
            self.ui.doubleSpinBox_period_min.blockSignals(True)
            self.ui.doubleSpinBox_period_min.setValue(value)
            self.ui.doubleSpinBox_period_min.blockSignals(False)

    def _phase_min_changed(self, value):
        if value > self.ui.doubleSpinBox_phase_max.value():
            self.ui.doubleSpinBox_phase_max.blockSignals(True)
            self.ui.doubleSpinBox_phase_max.setValue(value)
            self.ui.doubleSpinBox_phase_max.blockSignals(False)

    def _phase_max_changed(self, value):
        if value < self.ui.doubleSpinBox_phase_min.value():
            self.ui.doubleSpinBox_phase_min.blockSignals(True)
            self.ui.doubleSpinBox_phase_min.setValue(value)
            self.ui.doubleSpinBox_phase_min.blockSignals(False)

    def get_phase_limit(self):
        if self.ui.checkBox_phase.isChecked():
            return self.ui.doubleSpinBox_phase_min.value(), self.ui.doubleSpinBox_phase_max.value()
        else:
            return None

    def get_period_limit(self):
        if self.ui.checkBox_period.isChecked():
            return self.ui.doubleSpinBox_period_min, self.ui.doubleSpinBox_period_max
        else:
            return None

    def get_resistivity_limits(self):
        if self.ui.checkBox_resistivity.isChecked():
            return self.ui.doubleSpinBox_res_min, self.ui.doubleSpinBox_res_max
        else:
            return None

    def get_plot_xx(self):
        return 'y' if self.ui.checkBox_zxx.isChecked() else 'n'

    def get_plot_xy(self):
        return 'y' if self.ui.checkBox_zxy.isChecked() else 'n'

    def get_plot_yx(self):
        return 'y' if self.ui.checkBox_zyx.isChecked() else 'n'

    def get_plot_yy(self):
        return 'y' if self.ui.checkBox_zyy.isChecked() else 'n'

    def get_res_cmap(self):
        return cmapdict[str(self.ui.comboBox_res_cmap.currentText())]

    def get_phase_cmap(self):
        return cmapdict[str(self.ui.comboBox_phase_cmap.currentText())]

    def get_tickspace(self):
        return self.ui.spinBox_tickspace.value()


class PlotControlStrike(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_plot_control_strike()
        self.ui.setupUi(self)

        # setup ui
        # format validator
        self._double_validator = QDoubleValidator(-np.inf, np.inf, 1000)
        self.ui.lineEdit_min.setValidator(self._double_validator)
        self.ui.lineEdit_max.setValidator(self._double_validator)

        # connect signals
        self.ui.radioButton_range_minmax.toggled.connect(self._range_minmax_toggled)
        self.ui.lineEdit_min.editingFinished.connect(self._range_min_value_changed)
        self.ui.lineEdit_max.editingFinished.connect(self._range_max_value_changed)
        self.ui.checkBox_max_error.toggled.connect(self._max_error_toggled)

    def _max_error_toggled(self, is_checked):
        self.ui.doubleSpinBox_max_error.setEnabled(is_checked)

    def _range_minmax_toggled(self, is_checked):
        self.ui.lineEdit_min.setEnabled(is_checked)
        self.ui.lineEdit_max.setEnabled(is_checked)
        self.ui.label_min.setEnabled(is_checked)
        self.ui.label_max.setEnabled(is_checked)

    def _range_min_value_changed(self):
        value_min = np.float(str(self.ui.lineEdit_min.text()))
        value_max = np.float(str(self.ui.lineEdit_max.text()))
        if value_min > value_max:
            self.ui.lineEdit_max.blockSignals(True)
            self.ui.lineEdit_max.setText(str(value_min))
            self.ui.lineEdit_max.blockSignals(False)

    def _range_max_value_changed(self):
        value_min = np.float(str(self.ui.lineEdit_min.text()))
        value_max = np.float(str(self.ui.lineEdit_max.text()))
        if value_min > value_max:
            self.ui.lineEdit_min.blockSignals(True)
            self.ui.lineEdit_min.setText(str(value_max))
            self.ui.lineEdit_min.blockSignals(False)

    def get_plot_range(self):
        if self.ui.radioButton_range_data.isChecked():
            return 'data'
        else:
            value_min = np.float(str(self.ui.lineEdit_min.text()))
            value_max = np.float(str(self.ui.lineEdit_max.text()))
            return value_min, value_max

    def get_plot_type(self):
        return 1 if self.ui.radioButton_type_1.isChecked() else 2

    def get_plot_tipper(self):
        return 'y' if self.ui.checkBox_plot_tipper.isChecked() else 'n'

    def get_error_floor(self):
        return self.ui.doubleSpinBox_max_error.value() if self.ui.checkBox_max_error.isChecked() else None

    def get_fold(self):
        return self.ui.checkBox_fold.isChecked()
