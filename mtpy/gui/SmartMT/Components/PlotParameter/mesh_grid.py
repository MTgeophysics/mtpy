# -*- coding: utf-8 -*-
"""
    Description:


    Usage:

    Author: YingzhiGou
    Date: 24/10/2017
"""
from qtpy.QtWidgets import QGroupBox

from mtpy.gui.SmartMT.ui_asset.groupbox_mesh_grid import Ui_GroupBox_mash_grid


class MeshGrid(QGroupBox):
    def __init__(self, parent):
        QGroupBox.__init__(self, parent)
        self.ui = Ui_GroupBox_mash_grid()
        self.ui.setupUi(self)

        # connect signal
        self.ui.radioButton_imshow.toggled.connect(self._imshow_toggled)

    _grid_types = ["imshow", "pcolormesh"]
    _interpolation_methods = [
        "none",
        "nearest",
        "bilinear",
        "bicubic",
        "spline16",
        "spline36",
        "hanning",
        "hamming",
        "hermite",
        "kaiser",
        "quadric",
        "catrom",
        "gaussian",
        "bessel",
        "mitchell",
        "sinc",
        "lanczos",
    ]

    def _imshow_toggled(self, checked):
        if checked:
            self.ui.groupBox_interpolation_method.setHidden(False)
        else:
            self.ui.groupBox_interpolation_method.setHidden(True)

    def get_grid_type(self):
        if self.ui.radioButton_imshow.isChecked():
            return self._grid_types[0]
        elif self.ui.radioButton_pcolormesh.isChecked():
            return self._grid_types[1]
        else:
            return None  # should never reach here

    def get_interpolation_method(self):
        if not self.ui.groupBox_interpolation_method.isHidden():
            return self._interpolation_methods[
                self.ui.comboBox_interpolation_method.currentIndex()
            ]
        else:
            return None
