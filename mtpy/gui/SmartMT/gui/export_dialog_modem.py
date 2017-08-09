# -*- coding: utf-8 -*-
"""
    Description:
        dialog for exporting the selected .edi files to data files required by ModEm
    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import os

import numpy as np
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import pyqtSignal
from matplotlib.figure import Figure

from examples.create_modem_input import select_periods
from mtpy.gui.SmartMT.gui.plot_parameter_guis import Rotation
from mtpy.gui.SmartMT.ui_asset.wizard_export_modem import Ui_Wizard_esport_modem
from mtpy.utils.mtpylog import MtPyLog


class ExportDialogModEm(QtGui.QWizard):
    def __init__(self, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = Ui_Wizard_esport_modem()
        self.ui.setupUi(self)

        # setup gui
        # add rotation
        self._rotation_ui = Rotation(self.ui.wizardPage_data)
        self._rotation_ui.setTitle('Data Rotation Angle')
        self.ui.horizontalLayout_data.addWidget(self._rotation_ui)

        self._mesh_rotation_ui = Rotation(self.ui.wizardPage_mesh)
        self._mesh_rotation_ui.setTitle('Mesh Rotation Angle')
        self.ui.gridLayout_mesh.addWidget(self._mesh_rotation_ui)

        # patch some double spinBox
        # patch_QDoubleSpinBox(self.ui.doubleSpinBox_resistivity_air)
        # patch_QDoubleSpinBox(self.ui.doubleSpinBox_resistivity_sea)
        self._double_validator = QtGui.QDoubleValidator(-np.inf, np.inf, 1000)
        self._double_validator.setNotation(QtGui.QDoubleValidator.ScientificNotation)
        self.ui.lineEdit_resistivity_air.setValidator(self._double_validator)
        self.ui.lineEdit_resistivity_sea.setValidator(self._double_validator)

        # setup directory and dir dialog
        self._dir_dialog = QtGui.QFileDialog(self)
        # self._dir_dialog.setDirectory(os.path.expanduser("~"))
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
        self.ui.radioButton_impedance_full.toggled.connect(self._impedance_full_toggled)
        self.ui.radioButton_impedance_off_diagonal.toggled.connect(self._impedance_off_diagonal_toggled)
        self.ui.radioButton_impedance_none.toggled.connect(self._impedance_none_toggled)
        self.ui.radioButton_vertical_full.toggled.connect(self._vertical_full_toggled)

        self.ui.comboBox_error_type.currentIndexChanged.connect(self._error_type_changed)
        for component in self._impedance_components:
            combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
            checkbox = getattr(self.ui, 'checkBox_{}'.format(component))
            combobox.currentIndexChanged.connect(self._component_error_type_changed)
            checkbox.toggled.connect(self._error_component_checkbox_toggled(combobox))

        self.ui.comboBox_output_name.currentIndexChanged.connect(self._update_full_output)
        self.ui.comboBox_output_name.lineEdit().editingFinished.connect(self._output_name_changed)
        self.ui.comboBox_directory.currentIndexChanged.connect(self._update_full_output)
        self.ui.comboBox_directory.lineEdit().editingFinished.connect(self._output_dir_changed)
        self.ui.pushButton_browse.clicked.connect(self._browse)

        # self.ui.comboBox_topography_file.currentIndexChanged.connect()
        self.ui.comboBox_topography_file.lineEdit().editingFinished.connect(
            self._topography_file_changed
        )
        self.ui.pushButton_browse_topography_file.clicked.connect(self._browse_topography_file)

        self.ui.pushButton_test.clicked.connect(self._test_button_clicked)

        # register fields
        self.ui.wizardPage_output.registerField('output_path*', self.ui.lineEdit_full_output)
        self.ui.wizardPage_topography.registerField('topography_file*', self.ui.comboBox_topography_file)
        self.ui.wizardPage_topography.registerField('sea_resistivity',
                                                    self.ui.lineEdit_resistivity_sea)
        self.ui.wizardPage_topography.registerField('air_resistivity',
                                                    self.ui.lineEdit_resistivity_air)

        # attribute
        self._mt_objs = None

    _impedance_components = ['zxx', 'zxy', 'zyx', 'zyy']
    _vertical_components = ['tx', 'ty']

    _error_type = [
        'floor',  # 0
        'value',  # 1
        'egbert',  # 2
        'floor_egbert',  # 3
        'stddev',  # 4
        'sqr',  # 5
        'meansqr'  # 6
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

    def _impedance_full_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, 'checkBox_{}'.format(comp))
                checkbox.setEnabled(True)
            self.ui.radioButton_vertical_none.setEnabled(True)
            self.ui.groupBox_sign_impedance.setEnabled(True)

    def _impedance_off_diagonal_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, 'checkBox_{}'.format(comp))
                checkbox.setEnabled(comp == 'zxy' or comp == 'zyx')
            self.ui.radioButton_vertical_none.setEnabled(True)
            self.ui.groupBox_sign_impedance.setEnabled(True)

    def _impedance_none_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, 'checkBox_{}'.format(comp))
                checkbox.setEnabled(False)
            self.ui.radioButton_vertical_none.setEnabled(False)
            self.ui.groupBox_sign_impedance.setEnabled(False)

    def _vertical_full_toggled(self, checked):
        self.ui.doubleSpinBox_error_tipper.setEnabled(checked)
        self.ui.radioButton_impedance_none.setEnabled(checked)
        self.ui.groupBox_sign_vertical.setEnabled(checked)

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
        self._dir_dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        self._dir_dialog.setWindowTitle("Save to ...")
        if self._dir_dialog.exec_() == QtGui.QDialog.Accepted:
            directory = str(self._dir_dialog.selectedFiles()[0])
            self.ui.comboBox_directory.setEditText(directory)
            self._output_dir_changed()

    def _topography_file_changed(self, *args, **kwargs):
        topo_file = str(self.ui.comboBox_topography_file.currentText())
        topo_file = os.path.normpath(topo_file)
        # update
        index = self.ui.comboBox_topography_file.findText(topo_file)
        if index == -1:
            self.ui.comboBox_topography_file.addItem(topo_file)
        self.ui.comboBox_topography_file.setCurrentIndex(
            index if index >= 0 else self.ui.comboBox_topography_file.findText(topo_file)
        )

    def _browse_topography_file(self, *args, **kwargs):
        self._dir_dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        self._dir_dialog.setWindowTitle("Find Topography File...")
        if self._dir_dialog.exec_() == QtGui.QDialog.Accepted:
            file_path = str(self._dir_dialog.selectedFiles()[0])
            self.ui.comboBox_topography_file.setEditText(file_path)
            self._topography_file_changed()

    def _test_button_clicked(self, *args, **kwargs):
        self.export_data(True)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs

    def get_data_kwargs(self):
        kwargs = {
            'error_type': self._error_type[self.ui.comboBox_error_type.currentIndex()],
            'save_path': self.get_save_file_path(),
            'format': '1' if self.ui.radioButton_format_1.isChecked() else '2',
            'rotation_angle': self._rotation_ui.get_rotation_in_degree()
        }

        # comp_error_type
        if any(
                [getattr(self.ui, 'comboBox_error_type_{}'.format(component)).isEnabled()
                 for component in self._impedance_components]
        ):
            kwargs['comp_error_type'] = dict(
                [
                    (
                        component,
                        self._error_type[
                            getattr(self.ui, 'comboBox_error_type_{}'.format(component)).currentIndex()
                        ]
                    )
                    for component in self._impedance_components
                    if getattr(self.ui, 'comboBox_error_type_{}'.format(component)).isEnabled()
                ]
            )

        # error_floor
        if self.ui.doubleSpinBox_error_floor.isEnabled():
            kwargs['error_floor'] = self.ui.doubleSpinBox_error_floor.value()
        # error_value
        if self.ui.doubleSpinBox_error_value.isEnabled():
            kwargs['error_value'] = self.ui.doubleSpinBox_error_value.value()
        # error_egbert
        if self.ui.doubleSpinBox_error_egbert.isEnabled():
            kwargs['error_egbert'] = self.ui.doubleSpinBox_error_egbert.value()
        # error_tipper
        if self.ui.doubleSpinBox_error_tipper.isEnabled():
            kwargs['error_tipper'] = self.ui.doubleSpinBox_error_tipper.value() / 100.

        # wave signs
        if self.ui.groupBox_sign_impedance.isEnabled():
            kwargs['wave_sign_impedance'] = '+' if self.ui.radioButton_impedance_sign_plus.isChecked() \
                else '-'
        if self.ui.groupBox_sign_vertical.isEnabled():
            kwargs['wave_sign_tipper'] = '+' if self.ui.radioButton_vertical_sign_plus.isChecked() \
                else '-'

        # units
        kwargs['units'] = '[mV/km]/[nT]' if self.ui.radioButton_unit_mvkmnt.isChecked() \
            else '[V/m]/[T]' if self.ui.radioButton_unit_vmt.isChecked() \
            else 'Ohm'

        return kwargs

    def get_model_kwargs(self):
        kwargs = {
            'cell_size_east': self.ui.doubleSpinBox_cell_size_east.value(),
            'cell_size_north': self.ui.doubleSpinBox_cell_szie_north.value(),
            'pad_east': self.ui.spinBox_pad_east.value(),
            'pad_north': self.ui.spinBox_pad_north.value(),
            'pad_z': self.ui.spinBox_pad_z.value(),
            'pad_stretch_h': self.ui.doubleSpinBox_pad_stretch_h.value(),
            'pad_stretch_v': self.ui.doubleSpinBox_pad_stretch_v.value(),
            'z1_layer': self.ui.doubleSpinBox_z1_thickness.value(),
            'z_target_depth': self.ui.doubleSpinBox_target_depth.value(),
            'z_bottom': self.ui.doubleSpinBox_bottom.value(),
            'n_layers': self.ui.spinBox_num_layers.value(),
            'n_airlayers': self.ui.spinBox_num_air_layers.value(),
            'mesh_rotation_angle': self._mesh_rotation_ui.get_rotation_in_degree()
        }
        return kwargs

    def get_save_file_path(self):
        return str(self.ui.lineEdit_full_output.text())

    def get_air_resistivity(self):
        return float(str(self.ui.lineEdit_resistivity_air.text()))

    def get_sea_resistivity(self):
        return float(str(self.ui.lineEdit_resistivity_sea.text()))

    def get_topography_2_mesh_args(self):
        topo_file = str(self.ui.comboBox_topography_file.currentText())
        topo_array = None
        interp_method = 'nearest' if self.ui.radioButton_interpo_method_nearest.isChecked() \
            else 'linear' if self.ui.radioButton_interpo_method_linear.isChecked() \
            else 'cubic'
        air_resis = self.get_air_resistivity()
        sea_resis = self.get_sea_resistivity()
        return topo_file, topo_array, interp_method, air_resis, sea_resis

    def get_covariance_kwargs(self):
        kwargs = {
            'smoothing_east': self.ui.doubleSpinBox_smoothing_east.value(),
            'smoothing_north': self.ui.doubleSpinBox_smoothing_north.value(),
            'smoothing_z': self.ui.doubleSpinBox_smoothing_z.value(),
            'smoothing_num': self.ui.spinBox_smoothing_number.value(),
            'save_path': self.get_save_file_path()
        }
        return kwargs

    def export_data(self, test=False):
        if self._mt_objs is None:
            return

        worker = ModEMWorker(
            self,
            edi_list=[mt_obj.fn for mt_obj in self._mt_objs]
        )

        worker.start()
        worker.wait()


class ModEMWorker(QtCore.QThread):
    def __init__(self, parent, edi_list):
        QtCore.QThread.__init__(self, parent)
        self._logger = MtPyLog().get_mtpy_logger(__name__)

        self._edi_list = edi_list

    status_updated = pyqtSignal(str)
    figure_updated = pyqtSignal(Figure)

    def run(self):
        period_list = select_periods(self._edi_list)
