# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox
from qtpy import QtCore

from mtpy.gui.SmartMT.ui_asset.groupbox_ellipse import Ui_GroupBoxEllipse


class Ellipse(QGroupBox):
    """
    ellipse_dict defined for mtpy.imagining.phase_tensor_maps.PlogPhaseTensorMaps
    """
    _colorby = ['phimin', 'phimax', 'skew', 'skew_seg', 'normalized_skew', 'normalized_skew_seg', 'phidet',
                'ellipticity']
    _cmap = ['mt_yl2rd', 'mt_bl2yl2rd', 'mt_wh2bl', 'mt_rd2bl', 'mt_bl2wh2rd', 'mt_seg_bl2wh2rd', 'mt_rd2gr2bl']

    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBoxEllipse()
        self.ui.setupUi(self)
        # set tooltips for color by and cmap
        for i in range(self.ui.comboBoxColor_by.count()):
            self.ui.comboBoxColor_by.setItemData(i, self.ui.comboBoxColor_by.itemText(i), QtCore.Qt.ToolTipRole)
        for i in range(self.ui.comboBox_cmap.count()):
            self.ui.comboBox_cmap.setItemData(i, self.ui.comboBox_cmap.itemText(i), QtCore.Qt.ToolTipRole)

        self.ui.doubleSpinBox_min.editingFinished.connect(self._increase_max)
        self.ui.doubleSpinBox_max.editingFinished.connect(self._decrease_min)

    def _increase_max(self):
        value = self.ui.doubleSpinBox_min.value()
        if value > self.ui.doubleSpinBox_max.value():
            self.ui.doubleSpinBox_max.setValue(value)

    def _decrease_min(self):
        value = self.ui.doubleSpinBox_max.value()
        if value < self.ui.doubleSpinBox_min.value():
            self.ui.doubleSpinBox_min.setValue(value)

    def get_ellipse_dict(self, prefix='ellipse_'):
        ellipse_dict = {
            prefix+'size': self.ui.doubleSpinBox_size.value(),
            prefix+'colorby': self._colorby[self.ui.comboBoxColor_by.currentIndex()],
            prefix+'range': (
                self.ui.doubleSpinBox_min.value(), self.ui.doubleSpinBox_max.value(),
                self.ui.doubleSpinBox_step.value()),
            prefix+'cmap': self._cmap[self.ui.comboBox_cmap.currentIndex()]
        }
        return ellipse_dict
