# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.Components import COLORS, SIMPLE_COLORS
from mtpy.gui.SmartMT.ui_asset.groupbox_arrow import Ui_GroupBox_Arrow


class Arrow(QGroupBox):
    def __init__(self, parent, simple_color=True):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_Arrow()
        self.ui.setupUi(self)
        self._simple_color = simple_color
        if not self._simple_color:
            # use all colors available to matplot
            self.ui.comboBox_color_imaginary.clear()
            self.ui.comboBox_color_real.clear()
            cnames = [name for name, hex in COLORS]
            self.ui.comboBox_color_imaginary.addItems(cnames)
            self.ui.comboBox_color_real.addItems(cnames)

    _direction = [0, 1]

    def hide_size(self):
        self.ui.label_size.hide()
        self.ui.doubleSpinBox_size.hide()

    def hide_head_length(self):
        self.ui.label_head_length.hide()
        self.ui.doubleSpinBox_head_length.hide()

    def hide_head_width(self):
        self.ui.label_head_width.hide()
        self.ui.doubleSpinBox_head_width.hide()

    def hide_color_real(self):
        self.ui.label_color_real.hide()
        self.ui.comboBox_color_real.hide()

    def hide_color_imaginary(self):
        self.ui.label_color_imaginary.hide()
        self.ui.comboBox_color_imaginary.hide()

    def hide_threshold(self):
        self.ui.label_threshold.hide()
        self.ui.doubleSpinBox_threshold.hide()

    def hide_direction(self):
        self.ui.label_direction.hide()
        self.ui.comboBox_direction.hide()

    def get_arrow_dict(self):
        if self.ui.groupBox_advanced_options.isChecked():
            arrow_dict = {
                'arrow_size': self.ui.doubleSpinBox_size.value(),
                'arrow_head_length': self.ui.doubleSpinBox_head_length.value(),
                'arrow_head_width': self.ui.doubleSpinBox_head_width.value(),
                'arrow_lw': self.ui.doubleSpinBox_line_width.value(),
                'arrow_threshold': self.ui.doubleSpinBox_threshold.value(),
                'arrow_direction': self._direction[self.ui.comboBox_direction.currentIndex()]
            }
            if self._simple_color:
                arrow_dict['arrow_color'] = (SIMPLE_COLORS[self.ui.comboBox_color_real.currentIndex()],
                                             SIMPLE_COLORS[self.ui.comboBox_color_imaginary.currentIndex()])
            else:
                arrow_dict['arrow_color'] = (COLORS[self.ui.comboBox_color_real.currentIndex()][1],
                                             COLORS[self.ui.comboBox_color_imaginary.currentIndex()][1])
            return arrow_dict
        else:
            return None

    def get_plot_tipper(self):
        if self.ui.checkBox_real.isChecked() and self.ui.checkBox_imaginary.isChecked():
            return 'yri'
        elif self.ui.checkBox_real.isChecked():
            return 'yr'
        elif self.ui.checkBox_imaginary.isChecked():
            return 'yi'
        else:
            return 'n'
