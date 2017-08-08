# -*- coding: utf-8 -*-
"""
    Description:
        dialog for exporting the selected .edi files to data files required by ModEm
    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import os

from PyQt4 import QtGui, QtCore

from mtpy.gui.SmartMT.ui_asset.wizard_export_modem import Ui_Wizard_esport_modem


class ExportDialogModEm(QtGui.QWizard):
    def __init__(self, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = Ui_Wizard_esport_modem()
        self.ui.setupUi(self)

        # setup gui

        # setup directory and dir dialog
        self._dir_dialog = QtGui.QFileDialog(self)
        # self._dir_dialog.setDirectory(os.path.expanduser("~"))
        self._dir_dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        self._dir_dialog.setWindowTitle("Save to ...")
        self.ui.comboBox_directory.addItem(os.path.expanduser("~"))
        self._update_full_output()

        # tooltip for error types
        for index, tooltip in enumerate(self._error_type_tool_tip):
            self.ui.comboBox_error_type.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zxx.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zxy.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zyx.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)
            self.ui.comboBox_error_type_zyy.setItemData(index, tooltip, QtCore.Qt.ToolTipRole)

        # connect signals
        self.ui.comboBox_output_name.currentIndexChanged.connect(self._update_full_output)
        self.ui.comboBox_output_name.lineEdit().editingFinished.connect(self._output_name_changed)
        self.ui.comboBox_directory.currentIndexChanged.connect(self._update_full_output)
        self.ui.comboBox_directory.lineEdit().editingFinished.connect(self._output_dir_changed)
        self.ui.pushButton_browse.clicked.connect(self._browse)

        self.ui.comboBox_error_type.currentIndexChanged.connect(self._error_type_changed)
        for component in self._impedance_components:
            combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
            checkbox = getattr(self.ui, 'checkBox_{}'.format(component))
            combobox.currentIndexChanged.connect(self._component_error_type_changed)
            checkbox.toggled.connect(self._error_component_checkbox_toggled(combobox))

        # register fields
        self.ui.wizardPage_intro.registerField('output_name', self.ui.comboBox_output_name)
        self.ui.wizardPage_intro.registerField('output_directory', self.ui.comboBox_directory)

    _impedance_components = ['zxx', 'zxy', 'zyx', 'zyy']
    _vertical_components = ['tx', 'ty']

    _error_type = [
        'floor',            # 0
        'value',            # 1
        'egbert',           # 2
        'floor_egbert',     # 3
        'stddev',           # 4
        'sqr',              # 5
        'meansqr'           # 6
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

    def _error_component_checkbox_toggled(self, combobox):
        def _checkbox_toggled(checked):
            combobox.setEnabled(checked)
            self._component_error_type_changed()
        return _checkbox_toggled

    def _component_error_type_changed(self, error_type_index=-1):
        types = set([self.ui.comboBox_error_type.currentIndex()])
        for component in self._impedance_components:
            combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
            if combobox.isEnabled():
                types.add(combobox.currentIndex())
        self.ui.doubleSpinBox_error_floor.setEnabled(0 in types)
        self.ui.doubleSpinBox_error_egbert.setEnabled(2 in types or 3 in types)
        self.ui.doubleSpinBox_error_value.setEnabled(1 in types)

    def _error_type_changed(self, error_type_index):
        # sync the component error types with default
        for component in self._impedance_components:
            combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
            if not combobox.isEnabled():
                combobox.blockSignals(True)
                combobox.setCurrentIndex(error_type_index)
                combobox.blockSignals(False)
        self._component_error_type_changed()

    def _output_name_changed(self, *args, **kwargs):
        output_name = str(self.ui.comboBox_output_name.currentText())
        index = self.ui.comboBox_output_name.findText(output_name)
        if index == -1:
            self.ui.comboBox_output_name.addItem(output_name)
        self.ui.comboBox_output_name.setCurrentIndex(
            index if index >= 0 else self.ui.comboBox_output_name.findText(output_name)
        )

    def _output_dir_changed(self, *args, **kwargs):
        directory = str(self.ui.comboBox_directory.currentText())
        directory = os.path.normpath(directory)
        # update directory
        index = self.ui.comboBox_directory.findText(directory)
        if index == -1:
            self.ui.comboBox_directory.addItem(directory)
        self.ui.comboBox_directory.setCurrentIndex(index
                                                   if index >= 0
                                                   else self.ui.comboBox_directory.findText(directory))

    def _update_full_output(self, *args, **kwargs):
        directory = str(self.ui.comboBox_directory.currentText())
        output_name = str(self.ui.comboBox_output_name.currentText())
        full_output = os.path.normpath(os.path.join(directory, output_name))
        self.ui.lineEdit_full_output.setText(full_output)

    def _browse(self, *args, **kwargs):
        if self._dir_dialog.exec_() == QtGui.QDialog.Accepted:
            directory = str(self._dir_dialog.selectedFiles()[0])
            self.ui.comboBox_directory.setEditText(directory)
            self._output_dir_changed()
