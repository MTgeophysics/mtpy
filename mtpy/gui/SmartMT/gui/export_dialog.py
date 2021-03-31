# -*- coding: utf-8 -*-
"""
    Description:

    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import os
import tempfile
import webbrowser

import matplotlib.pyplot as plt
from PIL import Image
from qtpy import QtCore, QT_VERSION
from qtpy.QtWidgets import QDialog, QFileDialog, QMessageBox

if QT_VERSION.startswith("4"):
    from matplotlib.backends.backend_qt4agg import FigureCanvas
else:
    from matplotlib.backends.backend_qt5agg import FigureCanvas

from mtpy.gui.SmartMT.ui_asset.dialog_export import Ui_Dialog_Export
from mtpy.gui.SmartMT.ui_asset.dialog_preview import Ui_Dialog_preview
from mtpy.gui.SmartMT.utils.validator import DirectoryValidator

IMAGE_FORMATS = []
_temp_fig = plt.figure()
_file_types = _temp_fig.canvas.get_supported_filetypes()
for type, description in _file_types.items():
    IMAGE_FORMATS.append((type, description))
plt.close(_temp_fig)
del _temp_fig, _file_types


class ExportDialog(QDialog):
    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.ui = Ui_Dialog_Export()
        self.ui.setupUi(self)
        self._fig = None

        # setup file types
        for frmt in IMAGE_FORMATS:
            self.ui.comboBox_fileType.addItem("{0[1]} (.{0[0]})".format(frmt))

        self._file_name_changed()  # select the default format

        # setup directory and dir dialog
        self._dir_dialog = QFileDialog(self)
        # self._dir_dialog.setDirectory(os.path.expanduser("~"))
        self._dir_dialog.setFileMode(QFileDialog.DirectoryOnly)
        self._dir_dialog.setWindowTitle("Save to ...")
        self.ui.comboBox_directory.addItem(os.path.expanduser("~"))
        self.ui.pushButton_browse.clicked.connect(self._browse)

        self._dir_validator = DirectoryValidator()
        self.ui.comboBox_directory.lineEdit().setValidator(self._dir_validator)

        # file name
        self.ui.comboBox_fileName.currentIndexChanged.connect(self._file_name_changed)

        # file type
        self.ui.comboBox_fileType.currentIndexChanged.connect(self._file_type_changed)

        # cancel button
        self.ui.pushButton_cancel.clicked.connect(self._cancel_button_clicked)
        # export button
        self.ui.pushButton_export.clicked.connect(self._export_button_clicked)
        # preview button
        self.ui.pushButton_preview.clicked.connect(self._preview_button_clicked)

        # dpi
        self.ui.spinBox_dpi.valueChanged.connect(self._dpi_changed)
        # inches
        self.ui.doubleSpinBox_width_inches.valueChanged.connect(
            self._width_inches_changed
        )
        self.ui.doubleSpinBox_height_inches.valueChanged.connect(
            self._height_inches_changed
        )
        # pixels
        self.ui.spinBox_width_pixels.valueChanged.connect(self._width_pixels_changed)
        self.ui.spinBox_height_pixels.valueChanged.connect(self._height_pixels_changed)

        # message box for when the file already exist
        self._msg_box = QMessageBox(self)
        self._msg_box.setWindowTitle("Export...")
        self._msg_box_button_overwrite = self._msg_box.addButton(
            self.tr("Overwrite"), QMessageBox.AcceptRole
        )
        self._msg_box_button_save_as = self._msg_box.addButton(
            self.tr("Save As"), QMessageBox.ActionRole
        )
        self._msg_box_button_cancel = self._msg_box.addButton(QMessageBox.Cancel)
        self._msg_box.setDefaultButton(self._msg_box_button_save_as)

    _orientation = ["portrait", "landscape"]

    _no_alpha_channel_formats = ("jpg", "jpeg")  # ("png", "gif", "psd")

    def _dpi_changed(self, dpi):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(
            self.ui.doubleSpinBox_height_inches.value() * dpi
        )
        self.ui.spinBox_width_pixels.setValue(
            self.ui.doubleSpinBox_width_inches.value() * dpi
        )
        self.ui.doubleSpinBox_height_inches.blockSignals(False)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _width_pixels_changed(self, width):
        self.ui.doubleSpinBox_width_inches.blockSignals(True)
        new_width_inches = width / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_width_inches.setValue(new_width_inches)
        self.ui.doubleSpinBox_width_inches.blockSignals(False)

    def _height_pixels_changed(self, height):
        self.ui.doubleSpinBox_height_inches.blockSignals(True)
        new_height_inches = height / float(self.ui.spinBox_dpi.value())
        self.ui.doubleSpinBox_height_inches.setValue(new_height_inches)
        self.ui.doubleSpinBox_height_inches.blockSignals(False)

    def _width_inches_changed(self, width):
        self.ui.spinBox_width_pixels.blockSignals(True)
        self.ui.spinBox_width_pixels.setValue(width * self.ui.spinBox_dpi.value())
        self.ui.spinBox_width_pixels.blockSignals(False)

    def _height_inches_changed(self, height):
        self.ui.spinBox_height_pixels.blockSignals(True)
        self.ui.spinBox_height_pixels.setValue(height * self.ui.spinBox_dpi.value())
        self.ui.spinBox_height_pixels.blockSignals(False)

    def _cancel_button_clicked(self, b):
        self.reject()

    def _export_button_clicked(self, b):
        self.accept()

    def _preview_button_clicked(self):
        if self._fig:
            # set figures
            self._fig.set_size_inches(
                self.get_size_inches_width(), self.get_size_inches_height()
            )
            params = self.get_savefig_params()
            self._fig.set_dpi(params["dpi"])
            self._fig.set_tight_layout(
                True if params["bbox_inches"] == "tight" else False
            )

            canvas = FigureCanvas(self._fig)
            canvas.show()

            # dialog
            preview_dialog = PreviewDialog(self, self._fig)
            preview_dialog.exec_()

    def _file_type_changed(self, *args, **kwargs):
        index = self.ui.comboBox_fileType.currentIndex()
        ext, _ = IMAGE_FORMATS[index]
        filename = str(self.ui.comboBox_fileName.currentText())
        filename, _ = os.path.splitext(filename)

        if (
            ext in self._no_alpha_channel_formats
        ):  # enable transparent if the format supports
            self.ui.checkBox_transparent.setEnabled(False)
        else:
            self.ui.checkBox_transparent.setEnabled(True)

        # update file name
        filename = "%s.%s" % (filename, ext)
        index = self.ui.comboBox_fileName.findText(filename)
        if index == -1:
            self.ui.comboBox_fileName.addItem(filename)
        self.ui.comboBox_fileName.setCurrentIndex(
            index if index >= 0 else self.ui.comboBox_fileName.findText(filename)
        )

    def _file_name_changed(self, *args, **kwargs):
        filename = str(self.ui.comboBox_fileName.currentText())
        filename, extension = os.path.splitext(filename)
        extension = extension[1:]  # get ride of .
        # check if the extension is supported
        index = [i for i, (ext, dsc) in enumerate(IMAGE_FORMATS) if ext == extension]
        if index:
            self.ui.comboBox_fileType.setCurrentIndex(index[0])
        elif extension:
            # no extension
            pass
        else:
            # extension not supported:
            pass

    def _browse(self, *args, **kwargs):
        if self._dir_dialog.exec_() == QDialog.Accepted:
            dirs = (
                self._dir_dialog.selectedFiles()
            )  # behave differently in pyqt4 and pyqt5
            directory = str(
                dirs[0] if dirs else self._dir_dialog.directory().absolutePath()
            )  # this makes the behave the same
            # update directory
            index = self.ui.comboBox_directory.findText(directory)
            if index == -1:
                self.ui.comboBox_directory.addItem(directory)
            self.ui.comboBox_directory.setCurrentIndex(
                index if index >= 0 else self.ui.comboBox_directory.findText(directory)
            )

    def export_to_file(self, fig):
        self._fig = fig
        respawn = True
        while respawn:
            respawn = False
            self.ui.spinBox_dpi.setValue(fig.get_dpi())
            self.ui.doubleSpinBox_width_inches.setValue(fig.get_figwidth())
            self.ui.doubleSpinBox_height_inches.setValue(fig.get_figheight())
            response = self.exec_()
            if response == QDialog.Accepted:
                # saving files
                fname = self.get_save_file_name()

                if os.path.exists(fname):
                    new_name = generate_unique_file_name(fname)
                    self._show_file_exist_message(fname, new_name)
                    if self._msg_box.clickedButton() == self._msg_box_button_cancel:
                        respawn = True
                        continue
                    elif self._msg_box.clickedButton() == self._msg_box_button_save_as:
                        fname = new_name  # save_as
                    else:
                        pass  # use the original name to overwrite

                params = self.get_savefig_params()
                # change size
                fig.set_size_inches(
                    self.get_size_inches_width(), self.get_size_inches_height()
                )
                try:
                    fig.savefig(fname, **params)
                except IOError as err:
                    if "RGBA" in err.message:
                        # if the problem is RGBA as the alpha channel is not supported in the selected format
                        # save to png then save as
                        basename = os.path.basename(fname)
                        tmp_dir = tempfile.gettempdir()
                        filename, ext = os.path.splitext(basename)
                        png_file = filename + ".png"
                        final_format = params["format"]
                        params["format"] = "png"
                        new_fname = os.path.join(tmp_dir, png_file)
                        fig.savefig(new_fname, **params)
                        with Image.open(new_fname) as im:
                            rgb_im = im.convert("RGB")
                            # make sure the fname is ended with the right extension
                            fname, _ = os.path.splitext(fname)
                            fname += "." + final_format
                            rgb_im.save(fname)
                    else:
                        raise err

                if self.ui.checkBox_open_after_export.isChecked():
                    # open with the system default application, this should work on all platforms
                    webbrowser.open(fname)
                return fname
                # elif response == QtGu

    def get_size_inches_width(self):
        return self.ui.doubleSpinBox_width_inches.value()

    def get_size_inches_height(self):
        return self.ui.doubleSpinBox_height_inches.value()

    def _show_file_exist_message(self, fname, new_name):
        self._msg_box.setText(
            '<p>File "{0}" already exists. Do you want to overwrite the existing, or save to "{1}" instead?<\p>'.format(
                fname, new_name
            )
        )
        self._msg_box.exec_()

    def get_save_file_name(self):
        name = os.path.join(
            str(self.ui.comboBox_directory.currentText()),
            str(self.ui.comboBox_fileName.currentText()),
        )
        return os.path.normpath(name)

    def get_savefig_params(self):
        params = {
            "dpi": self.ui.spinBox_dpi.value(),
            "orientation": self.get_orientation(),
            "format": self.get_file_format()[0],
            "transparent": self.get_transparent(),
            "bbox_inches": self.get_bbox_inches(),
        }
        return params

    def get_transparent(self):
        return (
            self.ui.checkBox_transparent.isEnabled()
            and self.ui.checkBox_transparent.isChecked()
        )

    def get_file_format(self):
        return IMAGE_FORMATS[self.ui.comboBox_fileType.currentIndex()]

    def get_bbox_inches(self):
        return "tight" if self.ui.checkBox_tightBbox.isChecked() else None

    def get_orientation(self):
        return self._orientation[self.ui.comboBox_orientation.currentIndex()]

    def keyPressEvent(self, event):
        """Capture and ignore all key press events.

        This is used so that return key event does not trigger any button
        from the dialog. We need to allow the return key to be used in filters
        in the widget."""
        if event.key() == QtCore.Qt.Key_Escape:
            # call reject if Escape is pressed.
            self.reject()
        pass

        # def closeEvent(self, event):
        #     self._msg_box.deleteLater()
        #     super(ExportDialog, self).closeEvent(event)


class PreviewDialog(QDialog):
    def __init__(self, parent, fig):
        QDialog.__init__(self, parent)
        self.ui = Ui_Dialog_preview()
        self.ui.setupUi(self)
        self._canvas = FigureCanvas(fig)
        self.ui.verticalLayout_2.addWidget(self._canvas)
        self.resize(self.sizeHint())


def generate_unique_file_name(basename):
    name, ext = os.path.splitext(basename)
    counter = 1  # start from 2 as the first extension
    new_name = basename
    while os.path.exists(new_name):
        counter += 1
        new_name = "{0}_{1}{2}".format(name, counter, ext)
    return new_name
