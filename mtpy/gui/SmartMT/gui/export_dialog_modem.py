# -*- coding: utf-8 -*-
"""
    Description:
        dialog for exporting the selected .edi files to data files required by ModEm
    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.ui_asset.wizard_export_modem import Ui_Wizard_esport_modem


class ExportDialogModEm(QtGui.QWizard):
    def __init__(self, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = Ui_Wizard_esport_modem()
        self.ui.setupUi(self)

        # setup gui

        # tooltip for error types
        for index, tooltip in enumerate(self._error_type_tool_tip):
            self.ui.comboBox_error_type.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zxx.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zxy.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zyx.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zyy.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)

    _error_type = [
        'floor',
        'value',
        'egbert',
        'floor_egbert',
        'stddev',
        'sqr',
        'meansqr'
    ]

    _error_type_tool_tip = [
        'sets the error floor to the value of Error Floor',
        'sets error to error value',
        'sets error to the value of Error Egbert * sqrt(abs(zxy*zyx), see Egbert & Kelbert',
        'set error floor to error_egbert * sqrt(abs(zxy*zyx))',
        'use the standard deviation of the errors of a component across all frequencies for one station',
        'use square error of a frequency of component of a station',
        'mean square error of a frequency of a component across all frequencies for one station'
    ]
