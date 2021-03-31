# -*- coding: utf-8 -*-
"""
    Description:
        dialog for exporting the selected .edi files to data files required by ModEm
    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import datetime
import inspect
import os
import pprint
import webbrowser

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from qtpy import QtCore
from qtpy.QtGui import QDoubleValidator
from qtpy.QtWidgets import QWizard, QFileDialog, QDialog, QMessageBox

from mtpy.core.edi_collection import EdiCollection
from mtpy.gui.SmartMT.Components.PlotParameter import (
    FrequencySelection,
    Rotation,
    FrequencySelectionFromFile,
)
from mtpy.gui.SmartMT.gui.busy_indicators import ProgressBar
from mtpy.gui.SmartMT.gui.export_dialog import PreviewDialog
from mtpy.gui.SmartMT.gui.matplotlib_imabedding import MathTextLabel
from mtpy.gui.SmartMT.ui_asset.wizard_export_modem import Ui_Wizard_esport_modem
from mtpy.gui.SmartMT.utils.validator import FileValidator, DirectoryValidator
from mtpy.modeling.modem import Data, Covariance
from mtpy.modeling.modem import Model
from mtpy.mtpy_globals import epsg_dict
from mtpy.utils.mtpylog import MtPyLog


class ExportDialogModEm(QWizard):
    def __init__(self, parent=None):
        QWizard.__init__(self, parent)
        self.ui = Ui_Wizard_esport_modem()
        self.ui.setupUi(self)

        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)
        # setup gui
        # self.setWindowTitle("ModEM input file generator")

        # add math label
        self._math_label_sign_impedances = MathTextLabel(
            self,
            self._math_label_sign_text.format(
                "+" if self.ui.radioButton_impedance_sign_plus.isChecked() else "-"
            ),
        )
        self.ui.horizontalLayout_sign_impedance.addWidget(
            self._math_label_sign_impedances
        )
        self._math_label_sign_vertical = MathTextLabel(
            self,
            self._math_label_sign_text.format(
                "+" if self.ui.radioButton_vertical_sign_plus.isChecked() else "-"
            ),
        )
        self.ui.horizontalLayout_sign_vertical.addWidget(self._math_label_sign_vertical)

        # add math formulae of each error type
        self._math_label_elbert = MathTextLabel(
            self, "$E_{egbert}=e_Z\\times |Z_{xy}\\times Z_{yx}|^\\frac{1}{2}$",
        )
        self.ui.verticalLayout_error_types.addWidget(self._math_label_elbert)

        self._math_label_mean = MathTextLabel(
            self, "$E_{mean} = e_Z\\times mean(|Z_{xy}, Z_{yx}|)$"
        )
        self._math_label_mean.setHidden(True)
        self.ui.verticalLayout_error_types.addWidget(self._math_label_mean)

        self._math_label_eigen = MathTextLabel(
            self, "$E_{eigen} = e_Z\\times eigenvalues(Z)$"
        )
        self._math_label_eigen.setHidden(True)
        self.ui.verticalLayout_error_types.addWidget(self._math_label_eigen)

        self._math_label_median = MathTextLabel(
            self, "$E_{median} = e_Z\\times median(|Z_{xx},Z_{xy},Z_{yx},Z_{yy}|)$"
        )
        self._math_label_median.setHidden(True)
        self.ui.verticalLayout_error_types.addWidget(self._math_label_median)

        self._math_label_error_type_note = MathTextLabel(
            self, "where $E$ is the error and $e$ is the error value (z)."
        )
        self.ui.verticalLayout_error_types.addWidget(self._math_label_error_type_note)

        # add period selection
        self._period_select_ui = FrequencySelection(
            self.ui.wizardPage_period,
            show_period=True,
            show_frequency=False,
            allow_range_select=True,
            select_multiple=True,
        )
        self._period_select_ui.ui.checkBox_show_existing.setChecked(True)
        self._period_select_ui.setEnabled(False)
        self._period_select_ui.setHidden(True)
        self.ui.wizardPage_period.layout().addWidget(self._period_select_ui)

        self._period_select_from_file_ui = FrequencySelectionFromFile(
            self.ui.wizardPage_period
        )
        self._period_select_from_file_ui.setEnabled(False)
        self._period_select_from_file_ui.setHidden(True)
        self.ui.wizardPage_period.layout().addWidget(self._period_select_from_file_ui)

        # add rotation
        self._rotation_ui = Rotation(self.ui.wizardPage_data)
        self._rotation_ui.setTitle("Data Rotation Angle")
        self.ui.horizontalLayout_data.addWidget(self._rotation_ui)

        self._mesh_rotation_ui = Rotation(self.ui.wizardPage_mesh)
        self._mesh_rotation_ui.setTitle("Mesh Rotation Angle")
        self.ui.gridLayout_mesh.addWidget(self._mesh_rotation_ui)

        # hide error percents
        self._component_error_type_z_changed()

        # hide bottom in vertical mesh as it is not used in mesh gen
        self.ui.doubleSpinBox_bottom.hide()
        self.ui.label_bottom.hide()

        # epsg
        self.ui.comboBox_epsg.addItems([str(epsg) for epsg in sorted(epsg_dict.keys())])

        # set validators
        self._double_validator = QDoubleValidator(-np.inf, np.inf, 1000)
        self._double_validator.setNotation(QDoubleValidator.ScientificNotation)
        self.ui.lineEdit_resistivity_init.setValidator(self._double_validator)
        self.ui.lineEdit_resistivity_air.setValidator(self._double_validator)
        self.ui.lineEdit_resistivity_sea.setValidator(self._double_validator)

        self._file_validator = FileValidator()
        self.ui.comboBox_topography_file.lineEdit().setValidator(self._file_validator)
        self.ui.comboBox_topography_file.lineEdit().setMaxLength(256)

        self._dir_validator = DirectoryValidator()
        self.ui.comboBox_directory.lineEdit().setValidator(self._dir_validator)

        # setup directory and dir dialog
        self._dir_dialog = QFileDialog(self)
        # self._dir_dialog.setDirectory(os.path.expanduser("~"))
        self.ui.comboBox_directory.addItem(os.path.expanduser("~"))
        self._update_full_output()

        # set maximum
        # self.ui.spinBox_cell_num_ew.setMaximum(0xFFFFFFFF)
        # self.ui.spinBox_cell_num_ns.setMaximum(0xFFFFFFFF)

        # tooltip for error types
        for index, tooltip in enumerate(self._error_type_z_tool_tip):
            self.ui.comboBox_error_type_z.setItemData(
                index, tooltip, QtCore.Qt.ToolTipRole
            )
            self.ui.comboBox_error_type_zxx.setItemData(
                index, tooltip, QtCore.Qt.ToolTipRole
            )
            self.ui.comboBox_error_type_zxy.setItemData(
                index, tooltip, QtCore.Qt.ToolTipRole
            )
            self.ui.comboBox_error_type_zyx.setItemData(
                index, tooltip, QtCore.Qt.ToolTipRole
            )
            self.ui.comboBox_error_type_zyy.setItemData(
                index, tooltip, QtCore.Qt.ToolTipRole
            )

        # hide parts
        self.ui.groupBox_component_error_types.setHidden(True)
        self.ui.label_component_error_types.setHidden(True)

        # connect signals
        self.ui.radioButton_impedance_sign_plus.toggled.connect(
            lambda is_checked: self._math_label_sign_impedances.set_math_text(
                self._math_label_sign_text.format("+" if is_checked else "-")
            )
        )
        self.ui.radioButton_vertical_sign_plus.toggled.connect(
            lambda is_checked: self._math_label_sign_vertical.set_math_text(
                self._math_label_sign_text.format("+" if is_checked else "-")
            )
        )

        self.ui.radioButton_impedance_full.toggled.connect(self._impedance_full_toggled)
        self.ui.radioButton_impedance_off_diagonal.toggled.connect(
            self._impedance_off_diagonal_toggled
        )
        self.ui.radioButton_impedance_none.toggled.connect(self._impedance_none_toggled)
        self.ui.radioButton_vertical_full.toggled.connect(self._vertical_full_toggled)

        self.ui.comboBox_error_type_z.currentIndexChanged.connect(
            self._error_type_z_changed
        )
        for component in self._impedance_components:
            combobox = getattr(self.ui, "comboBox_error_type_{}".format(component))
            checkbox = getattr(self.ui, "checkBox_{}".format(component))
            combobox.currentIndexChanged.connect(self._component_error_type_z_changed)
            checkbox.toggled.connect(self._error_component_checkbox_toggled(combobox))

        self.ui.comboBox_output_name.currentIndexChanged.connect(
            self._update_full_output
        )
        self.ui.comboBox_output_name.lineEdit().editingFinished.connect(
            self._output_name_changed
        )
        self.ui.comboBox_output_name.editTextChanged.connect(self._update_full_output)
        self.ui.comboBox_directory.currentIndexChanged.connect(self._update_full_output)
        self.ui.comboBox_directory.lineEdit().editingFinished.connect(
            self._output_dir_changed
        )
        self.ui.comboBox_directory.editTextChanged.connect(self._update_full_output)
        self.ui.pushButton_browse.clicked.connect(self._browse)

        # self.ui.doubleSpinBox_target_depth.valueChanged.connect(self._target_depth_changed)
        self.ui.doubleSpinBox_target_depth.lineEdit().editingFinished.connect(
            self._target_depth_changed
        )
        self.ui.doubleSpinBox_bottom.lineEdit().editingFinished.connect(
            self._bottom_changed
        )
        # self.ui.doubleSpinBox_bottom.valueChanged.connect(self._bottom_changed)

        # self.ui.comboBox_topography_file.currentIndexChanged.connect()
        self.ui.comboBox_topography_file.lineEdit().editingFinished.connect(
            self._topography_file_changed
        )
        self.ui.pushButton_browse_topography_file.clicked.connect(
            self._browse_topography_file
        )

        self.ui.pushButton_test.clicked.connect(self._test_button_clicked)

        self.ui.checkBox_cell_num_ew.stateChanged.connect(
            lambda p_int: self.ui.spinBox_cell_num_ew.setEnabled(p_int != 0)
        )
        self.ui.checkBox_cell_num_ns.stateChanged.connect(
            lambda p_int: self.ui.spinBox_cell_num_ns.setEnabled(p_int != 0)
        )

        self.ui.checkBox_component_error_types.setEnabled(
            False
        )  # disabled due to missing implementations todo implement this option in modem.Data
        self.ui.checkBox_component_error_types.toggled.connect(
            self._component_error_type_toggled
        )
        self.ui.radioButton_select_period_percent.toggled.connect(
            lambda checked: self.ui.doubleSpinBox_select_period_percent.setEnabled(
                checked
            )
        )
        self.ui.radioButton_select_period.toggled.connect(
            lambda checked: self._period_select_ui.setEnabled(checked)
        )
        self.ui.radioButton_select_period.toggled.connect(
            lambda checked: self._period_select_ui.setHidden(not checked)
        )
        self.ui.radioButton_select_by_file.toggled.connect(
            lambda checked: self._period_select_from_file_ui.setEnabled(checked)
        )
        self.ui.radioButton_select_by_file.toggled.connect(
            lambda checked: self._period_select_from_file_ui.setHidden(not checked)
        )

        # register fields
        self.ui.wizardPage_output.registerField(
            "output_path*", self.ui.lineEdit_full_output
        )
        self.ui.wizardPage_topography.registerField(
            "topography_file*", self.ui.comboBox_topography_file
        )
        self.ui.wizardPage_topography.registerField(
            "sea_resistivity", self.ui.lineEdit_resistivity_sea
        )
        self.ui.wizardPage_topography.registerField(
            "air_resistivity", self.ui.lineEdit_resistivity_air
        )

        # attribute
        self._mt_objs = None

        self._progress_bar = ProgressBar()

    _impedance_components = ["zxx", "zxy", "zyx", "zyy"]
    _vertical_components = ["tx", "ty"]
    _math_label_sign_text = "$exp({}i\\omega t)$"

    _error_type_z = [
        "egbert",  # 0
        "mean_od",  # 1
        "eigen",  # 2
        "median",  # 3
    ]

    _error_type_z_tool_tip = [
        "sets error to the value of Error Egbert * sqrt(abs(zxy*zyx), see Egbert & Kelbert",
        "sets error to error_value_z * mean([Zxy, Zyx]) (non zeros)",
        "sets error to error_value_z * eigenvalues(Z[ii])",
        "sets error to error_value_z * median([Zxx, Zxy, Zyx, Zyy]) (non zeros)",
    ]

    _error_type_tipper = ["abs", "floor"]

    def _component_error_type_toggled(self, is_checked):
        self.ui.label_component_error_types.setHidden(not is_checked)
        self.ui.groupBox_component_error_types.setHidden(not is_checked)

    def _target_depth_changed(self, *args):
        value = self.ui.doubleSpinBox_target_depth.value()
        # target depth has too be at least as deep as the bottom
        if self.ui.doubleSpinBox_bottom.value() < value:
            self.ui.doubleSpinBox_bottom.setValue(value)

    def _bottom_changed(self, *args):
        value = self.ui.doubleSpinBox_bottom.value()
        # bottom as to be at least at least as deep as the target depth
        if self.ui.doubleSpinBox_target_depth.value() > value:
            self.ui.doubleSpinBox_target_depth.setValue(value)

    def _impedance_full_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, "checkBox_{}".format(comp))
                checkbox.setEnabled(True)
            self.ui.radioButton_vertical_none.setEnabled(True)
            self.ui.groupBox_sign_impedance.setEnabled(True)

    def _impedance_off_diagonal_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, "checkBox_{}".format(comp))
                checkbox.setEnabled(comp == "zxy" or comp == "zyx")
            self.ui.radioButton_vertical_none.setEnabled(True)
            self.ui.groupBox_sign_impedance.setEnabled(True)

    def _impedance_none_toggled(self, checked):
        if checked:
            for comp in self._impedance_components:
                checkbox = getattr(self.ui, "checkBox_{}".format(comp))
                checkbox.setEnabled(False)
            self.ui.radioButton_vertical_none.setEnabled(False)
            self.ui.groupBox_sign_impedance.setEnabled(False)

    def _vertical_full_toggled(self, checked):
        self.ui.doubleSpinBox_error_value_tipper.setEnabled(checked)
        self.ui.radioButton_impedance_none.setEnabled(checked)
        self.ui.groupBox_sign_vertical.setEnabled(checked)

    def _error_component_checkbox_toggled(self, combobox):
        def _checkbox_toggled(checked):
            combobox.setEnabled(checked)
            self._component_error_type_z_changed()

        return _checkbox_toggled

    def _component_error_type_z_changed(self, error_type_index=-1):
        # types = {self.ui.comboBox_error_type_z.currentIndex()}
        # for component in self._impedance_components:
        #     combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
        #     if combobox.isEnabled():
        #         types.add(combobox.currentIndex())
        #
        # hidden = bool(0 not in types)
        # self.ui.label_error_floor.setHidden(hidden)
        # self.ui.doubleSpinBox_error_floor.setHidden(hidden)
        # hidden = bool(2 not in types and 3 not in types)
        # self.ui.label_error_egbert.setHidden(hidden)
        # self.ui.doubleSpinBox_error_egbert.setHidden(hidden)
        # self._math_label_elbert.setHidden(hidden)
        # hidden = bool(1 not in types)
        # self.ui.label_error_value.setHidden(hidden)
        # self.ui.doubleSpinBox_error_value.setHidden(hidden)
        pass  # todo re-enable after the sub-type for each component is implemented

    def _error_type_z_changed(self, error_type_index):
        # sync the component error types with default
        # for component in self._impedance_components:
        #     combobox = getattr(self.ui, 'comboBox_error_type_{}'.format(component))
        #     if not combobox.isEnabled():
        #         combobox.blockSignals(True)
        #         combobox.setCurrentIndex(error_type_index)
        #         combobox.blockSignals(False)
        # self._component_error_type_z_changed()
        pass  # todo re-enable after the sub-type for each component is implemented

        self._math_label_elbert.setHidden(error_type_index != 0)
        self._math_label_mean.setHidden(error_type_index != 1)
        self._math_label_eigen.setHidden(error_type_index != 2)
        self._math_label_median.setHidden(error_type_index != 3)

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
        self.ui.comboBox_directory.setCurrentIndex(
            index if index >= 0 else self.ui.comboBox_directory.findText(directory)
        )

    def _update_full_output(self, *args, **kwargs):
        directory = str(self.ui.comboBox_directory.currentText())
        output_name = str(self.ui.comboBox_output_name.currentText())
        full_output = os.path.normpath(os.path.join(directory, output_name))
        self.ui.lineEdit_full_output.setText(full_output)

    def _browse(self, *args, **kwargs):
        self._dir_dialog.setFileMode(QFileDialog.DirectoryOnly)
        self._dir_dialog.setWindowTitle("Save to ...")
        if self._dir_dialog.exec_() == QDialog.Accepted:
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
            index
            if index >= 0
            else self.ui.comboBox_topography_file.findText(topo_file)
        )

    def _browse_topography_file(self, *args, **kwargs):
        self._dir_dialog.setFileMode(QFileDialog.AnyFile)
        self._dir_dialog.setWindowTitle("Find Topography File...")
        if self._dir_dialog.exec_() == QDialog.Accepted:
            file_path = str(self._dir_dialog.selectedFiles()[0])
            self.ui.comboBox_topography_file.setEditText(file_path)
            self._topography_file_changed()

    def _test_button_clicked(self, *args, **kwargs):
        self.export_data(True)

    def set_data(self, mt_objs):
        self._mt_objs = mt_objs
        self._period_select_ui.set_data(self._mt_objs)
        self._period_select_from_file_ui.set_data(self._mt_objs)
        self.ui.listWidget_edi_files.clear()
        for mt_obj in mt_objs:
            self.ui.listWidget_edi_files.addItem(
                "{mt.station} ({mt.fn})".format(mt=mt_obj)
            )

    def get_inversion_mode(self):
        if self.ui.radioButton_impedance_full.isChecked():
            return "1" if self.ui.radioButton_vertical_full.isChecked() else "2"
        elif self.ui.radioButton_impedance_off_diagonal.isChecked():
            return "3" if self.ui.radioButton_vertical_full.isChecked() else "4"
        else:
            return "5"

    def _get_error_type_z(self):
        type = self._error_type_z[self.ui.comboBox_error_type_z.currentIndex()]
        if self.ui.checkBox_error_type_z_floor.isChecked():
            type += "_floor"
        return type

    def get_data_kwargs(self):
        kwargs = {
            "error_type_z": self._get_error_type_z(),
            "error_value_z": self.ui.doubleSpinBox_error_value_z.value(),
            "error_type_tipper": self._error_type_tipper[
                self.ui.comboBox_error_type_tipper.currentIndex()
            ],
            "error_value_tipper": self.ui.doubleSpinBox_error_value_tipper.value()
            / 100.0,
            "save_path": self.get_save_file_path(),
            "format": "1" if self.ui.radioButton_format_1.isChecked() else "2",
            "rotation_angle": self._rotation_ui.get_rotation_in_degree(),
            "model_epsg": self.get_epsg(),
            "inv_mode": self.get_inversion_mode(),
        }

        # comp_error_type
        if any(
            [
                getattr(self.ui, "comboBox_error_type_{}".format(component)).isEnabled()
                for component in self._impedance_components
            ]
        ):
            kwargs["comp_error_type"] = dict(
                [
                    (
                        component,
                        self._error_type_z[
                            getattr(
                                self.ui, "comboBox_error_type_{}".format(component)
                            ).currentIndex()
                        ],
                    )
                    for component in self._impedance_components
                    if getattr(
                        self.ui, "comboBox_error_type_{}".format(component)
                    ).isEnabled()
                ]
            )

        # wave signs
        if self.ui.groupBox_sign_impedance.isEnabled():
            kwargs["wave_sign_impedance"] = (
                "+" if self.ui.radioButton_impedance_sign_plus.isChecked() else "-"
            )
        if self.ui.groupBox_sign_vertical.isEnabled():
            kwargs["wave_sign_tipper"] = (
                "+" if self.ui.radioButton_vertical_sign_plus.isChecked() else "-"
            )

        # units
        kwargs["units"] = (
            "[mV/km]/[nT]"
            if self.ui.radioButton_unit_mvkmnt.isChecked()
            else "[V/m]/[T]"
            if self.ui.radioButton_unit_vmt.isChecked()
            else "Ohm"
        )

        return kwargs

    def get_model_kwargs(self):
        kwargs = {
            "save_path": self.get_save_file_path(),
            "cell_size_east": self.ui.doubleSpinBox_cell_size_east.value(),
            "cell_size_north": self.ui.doubleSpinBox_cell_szie_north.value(),
            "pad_east": self.ui.spinBox_pad_east.value(),
            "pad_north": self.ui.spinBox_pad_north.value(),
            "pad_z": self.ui.spinBox_pad_z.value(),
            "pad_stretch_h": self.ui.doubleSpinBox_pad_stretch_h.value(),
            "pad_stretch_v": self.ui.doubleSpinBox_pad_stretch_v.value(),
            "z1_layer": self.ui.doubleSpinBox_z1_thickness.value(),
            "z_target_depth": self.ui.doubleSpinBox_target_depth.value(),
            "z_bottom": self.ui.doubleSpinBox_bottom.value(),
            "n_layers": self.ui.spinBox_num_layers.value(),
            "n_air_layers": self.ui.spinBox_num_air_layers.value(),
            "res_model": self.get_initial_resistivity(),
            "mesh_rotation_angle": self._mesh_rotation_ui.get_rotation_in_degree(),
            "cell_number_ew": self.ui.spinBox_cell_num_ew.value()
            if self.ui.checkBox_cell_num_ew.isChecked()
            else None,
            "cell_number_ns": self.ui.spinBox_cell_num_ns.value()
            if self.ui.checkBox_cell_num_ns.isChecked()
            else None,
        }
        return kwargs

    def get_save_file_path(self):
        return str(self.ui.lineEdit_full_output.text())

    def get_initial_resistivity(self):
        return float(str(self.ui.lineEdit_resistivity_init.text()))

    def get_air_resistivity(self):
        return float(str(self.ui.lineEdit_resistivity_air.text()))

    def get_sea_resistivity(self):
        return float(str(self.ui.lineEdit_resistivity_sea.text()))

    def get_topography_2_mesh_args(self):
        kwargs = {
            "topography_file": str(self.ui.comboBox_topography_file.currentText()),
            "topography_array": None,
            "interp_method": "nearest"
            if self.ui.radioButton_interpo_method_nearest.isChecked()
            else "linear"
            if self.ui.radioButton_interpo_method_linear.isChecked()
            else "cubic",
            "air_resistivity": self.get_air_resistivity(),
            "sea_resistivity": self.get_sea_resistivity(),
        }

        return kwargs

    def get_covariance_kwargs(self):
        kwargs = {
            "smoothing_east": self.ui.doubleSpinBox_smoothing_east.value(),
            "smoothing_north": self.ui.doubleSpinBox_smoothing_north.value(),
            "smoothing_z": self.ui.doubleSpinBox_smoothing_z.value(),
            "smoothing_num": self.ui.spinBox_smoothing_number.value(),
            "save_path": self.get_save_file_path(),
        }
        return kwargs

    def get_select_period_kwargs(self):
        period_list = None
        if self.ui.radioButton_select_period.isChecked():
            period_list = self._period_select_ui.get_frequencies()
        elif self.ui.radioButton_select_by_file.isChecked():
            period_list = self._period_select_from_file_ui.get_selected_periods()
        elif self.ui.radioButton_select_by_decade.isChecked():
            period_list = []
            unique_periods = set()
            for mt in self._mt_objs:
                unique_periods.update(1.0 / mt.Z.freq)
            unique_periods = np.array(sorted(list(unique_periods), reverse=True))
            num = self.ui.spinBox_num_per_decade.value()
            decade_min = 10.0 ** int(np.log(unique_periods[0]))
            decade_max = decade_min * 10.0
            while (
                unique_periods[-1] < decade_min
            ):  # todo this block of code can be more efficient
                indexes = np.where(
                    (unique_periods <= decade_max) & (unique_periods > decade_min)
                )[0]
                if len(indexes) < num:
                    self._logger.warn(
                        "Selecting {} periods per decade, but there is only {} in [{},{}]".format(
                            num, len(indexes), decade_min, decade_max
                        )
                    )
                else:
                    # sample n
                    np.random.shuffle(indexes)
                    indexes = indexes[0:num]
                period_list += list(unique_periods[indexes])
                decade_max = decade_min
                decade_min /= 10.0

        kwargs = {
            "period_list": period_list,
            "percentage": self.ui.doubleSpinBox_select_period_percent.value(),
        }
        return kwargs

    def get_epsg(self):
        return int(str(self.ui.comboBox_epsg.currentText()))

    def export_data(self, test=False):
        if self._mt_objs is None:
            return
        # self._progress_bar.progressbar.setRange(0, 1)
        self._progress_bar.progressbar.setRange(0, 0)
        self._progress_bar.onStart()

        self._update_full_output()

        worker = ModEMWorker(
            self,
            show=test,
            edi_list=[mt_obj.fn for mt_obj in self._mt_objs],
            select_period_kwargs=self.get_select_period_kwargs(),
            data_kwargs=self.get_data_kwargs(),
            mesh_kwargs=self.get_model_kwargs(),
            topo_args=self.get_topography_2_mesh_args(),
            covariance_kwargs=self.get_covariance_kwargs(),
        )

        if test:
            worker.figure_updated.connect(self._show_figure)
        self._progress_bar.setWindowTitle(
            "Testing ModEM Data..." if test else "Generating ModEM Data..."
        )
        # worker.started.connect(self._progress_bar.onStart)
        # worker.finished.connect(self._progress_bar.onFinished)
        worker.status_updated.connect(self._update_progress_bar_text)
        worker.export_error.connect(self._export_error)
        # self._plot_opened.connect(worker.pause)

        # worker.start()
        # worker.wait()
        worker.run()

        # clean up
        worker.deleteLater()
        if self.ui.checkBox_open_output_dir.isChecked():
            webbrowser.open(self.get_save_file_path())

        self._progress_bar.onFinished()
        self._progress_bar.progressbar.setRange(0, 1)

    # _plot_opened = pyqtSignal(bool)

    def _update_progress_bar_text(self, text):
        self._progress_bar.incrementValue()
        self._progress_bar.updateIndicatorText(text)

    def _show_figure(self, string, fig):
        # self._plot_opened.emit(True)
        # print "plot_opened"
        preview_dialog = PreviewDialog(None, fig)
        preview_dialog.setWindowTitle(string)
        preview_dialog.exec_()
        preview_dialog.deleteLater()
        # self._plot_opened.emit(False)

    def _export_error(self, message):
        QMessageBox.critical(self, "Export Error", message, QMessageBox.Close)


class ModEMWorker(QtCore.QThread):
    def __init__(
        self,
        parent,
        edi_list,
        select_period_kwargs,
        data_kwargs,
        mesh_kwargs,
        topo_args,
        covariance_kwargs,
        show=False,
    ):
        QtCore.QThread.__init__(self, parent)
        self._logger = MtPyLog.get_mtpy_logger(self.__class__.__name__)

        self._edi_list = edi_list
        self._select_period_kwargs = select_period_kwargs
        self._data_kwargs = data_kwargs
        self._mesh_kwagrs = mesh_kwargs
        self._topo_args = topo_args
        self._covariance_kwargs = covariance_kwargs

        self.output_dir = self._data_kwargs["save_path"]

        self.show = show

    status_updated = QtCore.Signal(str)
    figure_updated = QtCore.Signal(str, Figure)
    export_error = QtCore.Signal(str)

    _period_image_name = "periods.png"
    _mesh_image_name = "mesh.png"
    _topo_image_name = "topography.png"
    _readme_name = "README.txt"

    # def pause(self, paused):
    #     print "paused"
    #     self.is_paused = paused

    def run(self):
        # monkey patch plt.show() so the plot is not displayed in worker thread
        true_plt_show = plt.show
        plt.show = _fake_plt_show
        self._figures = []
        try:
            if not os.path.exists(self.output_dir):
                os.mkdir(self.output_dir)

            # get period_list list
            self.status_updated.emit("Selecting Periods...")
            period_list = EdiCollection(self._edi_list).select_periods(
                **self._select_period_kwargs
            )
            # save period plot for reference
            figure = plt.gcf()
            figure.savefig(os.path.join(self.output_dir, self._period_image_name))
            if self.show:
                self.figure_updated.emit("Period Distribution", figure)

            # data object
            self.status_updated.emit("Creating ModEM Data Object...")
            self._data_kwargs["period_list"] = period_list
            data = Data(edi_list=self._edi_list, **self._data_kwargs)
            # write data file
            self.status_updated.emit("Writing Initial ModEM Data File...")
            data.write_data_file()

            # create model
            self.status_updated.emit("Creating Mesh Model...")
            model = Model(data_object=data, **self._mesh_kwagrs)
            model.make_mesh()
            # plot mesh
            model.plot_mesh(fig_num=plt.gcf().number + 1)
            figure = plt.gcf()
            figure.savefig(os.path.join(self.output_dir, self._mesh_image_name))
            if self.show:
                self.figure_updated.emit("Mesh", figure)
                # model.plot_mesh_xy()
                # self.figure_updated.emit("Mesh XY", plt.gcf())
                # model.plot_mesh_xz()
                # self.figure_updated.emit("Mesh XZ", plt.gcf())
            model.write_model_file()

            # add topography
            self.status_updated.emit("Adding Topography...")
            model.add_topography_to_mesh(**self._topo_args)
            if self.show:
                model.plot_topograph()  # this is too slow so only plot and save image when asked
                figure = plt.gcf()
                figure.savefig(os.path.join(self.output_dir, self._topo_image_name))
                self.figure_updated.emit("Topography", figure)
            self.status_updated.emit("Updating Mesh Model...")
            model.write_model_file()

            # covariance
            self.status_updated.emit("Creating Covariance File...")
            self._covariance_kwargs["mask_arr"] = model.covariance_mask
            cov = Covariance(**self._covariance_kwargs)
            self.status_updated.emit("Writing Covariance File...")
            cov.write_covariance_file(
                model_fn=model.model_fn,
                sea_water=self._topo_args["sea_resistivity"],
                air=self._topo_args["air_resistivity"],
            )

            self.status_updated.emit("Creating README File...")
            self.write_readme(os.path.join(self.output_dir, self._readme_name))
            # done
            self.status_updated.emit("Finishing...")

        except Exception as e:
            frm = inspect.trace()[-1]
            mod = inspect.getmodule(frm[0])
            self.export_error.emit("{}: {}".format(mod.__name__, e.message))

        # restore plt.show()
        plt.show = true_plt_show

    def write_readme(self, fname):
        with open(fname, "w") as outf:
            outf.write(
                "{} Data File Generated by MTPy2\n".format(
                    "ModEM"
                    if self._data_kwargs["format"] == "1"
                    else "WinMT (wrong name)"
                )
            )
            outf.write("====================\n")
            outf.write("Created on: {}\n".format(str(datetime.datetime.now())))
            outf.write("\n")
            outf.write(
                '**NOTE** All paths are modified to hide identity. "[*]" indicates the hidden path.\n'
            )
            outf.write("## EDI Files:\n")
            # hide path that may review user's identity
            outf.write(
                "\n".join(
                    [
                        os.path.normpath(
                            os.path.join(
                                "[*]/" + os.path.basename(os.path.dirname(edi_file)),
                                os.path.basename(edi_file),
                            )
                        )
                        for edi_file in self._edi_list
                    ]
                )
            )
            outf.write("\n\n")
            outf.write("## Parameters Used:\n\n")
            outf.write("### Period Included\n")
            pprint.pprint(self._select_period_kwargs, stream=outf)
            outf.write("\n")
            outf.write("\nSelected Period List:\n")
            pprint.pprint(
                EdiCollection(self._edi_list).select_periods(
                    **self._select_period_kwargs
                ),
                stream=outf,
            )
            outf.write("\n")
            outf.write("### Data:\n")
            kwargs = self._data_kwargs.copy()
            # hide path that may review user's identity
            kwargs["save_path"] = os.path.normpath(
                "[*]/" + os.path.basename(kwargs["save_path"])
            )
            pprint.pprint(kwargs, stream=outf)
            outf.write("\n")
            outf.write("### Mesh Model:\n")
            kwargs = self._mesh_kwagrs.copy()
            # hide path that may review user's identity
            kwargs["save_path"] = os.path.normpath(
                "[*]/" + os.path.basename(kwargs["save_path"])
            )
            pprint.pprint(kwargs, stream=outf)
            outf.write("\n")
            outf.write("### Topography:\n")
            kwargs = self._topo_args.copy()
            # hide path that may review user's identity
            kwargs["topography_file"] = os.path.normpath(
                "[*]/" + os.path.basename(kwargs["topography_file"])
            )
            pprint.pprint(kwargs, stream=outf)
            outf.write("\n")
            outf.write("### Covariance:\n")
            kwargs = self._covariance_kwargs.copy()
            # hide path that may review user's identity
            kwargs["save_path"] = os.path.normpath(
                "[*]/" + os.path.basename(kwargs["save_path"])
            )
            pprint.pprint(kwargs, stream=outf)
            outf.write("\n")
            outf.write("## Appendix I - Data Parameter Definition\n")
            outf.write(Data.__doc__)
            outf.write("\n")
            outf.write("## Appendix II - Mesh Model Parameter Definition\n")
            outf.write(Model.__doc__)
            outf.write("\n")
            outf.write("## Appendix III - Covariance Parameter Definition\n")
            outf.write(Covariance.__doc__)  # todo update Covariance class document
            outf.write("\n")


def _fake_plt_show(*args, **kwargs):
    # print "suppressed show"
    pass
