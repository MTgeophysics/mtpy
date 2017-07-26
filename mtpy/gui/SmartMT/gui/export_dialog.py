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

from PIL import Image
from PyQt4 import QtGui, QtCore
import matplotlib.pyplot as plt

from mtpy.gui.SmartMT.ui_asset.dialog_export import Ui_Dialog_Export

IMAGE_FORMATS = []
filetypes = plt.gcf().canvas.get_supported_filetypes()
for type, description in filetypes.iteritems():
    IMAGE_FORMATS.append((type, description))


class ExportDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = Ui_Dialog_Export()
        self.ui.setupUi(self)

        # setup file types
        for frmt in IMAGE_FORMATS:
            self.ui.comboBox_fileType.addItem("{0[1]} (.{0[0]})".format(frmt))

        self._file_name_changed()  # select the default format

        # setup directory and dir dialog
        self._dir_dialog = QtGui.QFileDialog(self)
        # self._dir_dialog.setDirectory(os.path.expanduser("~"))
        self._dir_dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        self._dir_dialog.setWindowTitle("Save to ...")
        self.ui.comboBox_directory.addItem(os.path.expanduser("~"))
        self.ui.pushButton_browse.clicked.connect(self._browse)

        # file name
        self.ui.comboBox_fileName.currentIndexChanged.connect(self._file_name_changed)

        # file type
        self.ui.comboBox_fileType.currentIndexChanged.connect(self._file_type_changed)

        # cancel button
        self.ui.pushButton_cancel.clicked.connect(self._cancel_button_clicked)
        # export button
        self.ui.pushButton_export.clicked.connect(self._export_button_clicked)

        # message box for when the file already exist
        self._msg_box = QtGui.QMessageBox(self)
        self._msg_box.setWindowTitle("Export...")
        self._msg_box_button_overwrite = self._msg_box.addButton(self.tr("Overwrite"), QtGui.QMessageBox.AcceptRole)
        self._msg_box_button_save_as = self._msg_box.addButton(self.tr("Save As"), QtGui.QMessageBox.ActionRole)
        self._msg_box_button_cancel = self._msg_box.addButton(QtGui.QMessageBox.Cancel)
        self._msg_box.setDefaultButton(self._msg_box_button_save_as)

    _orientation = ['portrait', 'landscape']

    _no_alpha_channel_formats = ('jpg', 'jpeg')  # ("png", "gif", "psd")

    def _cancel_button_clicked(self, b):
        self.reject()

    def _export_button_clicked(self, b):
        self.accept()

    def _file_type_changed(self, *args, **kwargs):
        index = self.ui.comboBox_fileType.currentIndex()
        ext, _ = IMAGE_FORMATS[index]
        filename = str(self.ui.comboBox_fileName.currentText())
        filename, _ = os.path.splitext(filename)

        if ext in self._no_alpha_channel_formats:  # enable transparent if the format supports
            self.ui.checkBox_transparent.setEnabled(False)
        else:
            self.ui.checkBox_transparent.setEnabled(True)

        # update file name
        filename = "%s.%s" % (filename, ext)
        index = self.ui.comboBox_fileName.findText(filename)
        if index == -1:
            self.ui.comboBox_fileName.addItem(filename)
        self.ui.comboBox_fileName.setCurrentIndex(index if index >= 0
                                                  else self.ui.comboBox_fileName.findText(filename))

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
        if self._dir_dialog.exec_() == QtGui.QDialog.Accepted:
            directory = str(self._dir_dialog.selectedFiles()[0])
            # update directory
            index = self.ui.comboBox_directory.findText(directory)
            if index == -1:
                self.ui.comboBox_directory.addItem(directory)
            self.ui.comboBox_directory.setCurrentIndex(index
                                                       if index >= 0
                                                       else self.ui.comboBox_directory.findText(directory))

    def export_to_file(self, fig):
        respawn = True
        while respawn:
            respawn = False
            response = self.exec_()
            if response == QtGui.QDialog.Accepted:
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
                try:
                    fig.savefig(fname, **params)
                except IOError as err:
                    if 'RGBA' in err.message:
                        # if the problem is RGBA as the alpha channel is not supported in the selected format
                        # save to png then save as
                        basename = os.path.basename(fname)
                        tmp_dir = tempfile.gettempdir()
                        filename, ext = os.path.splitext(basename)
                        png_file = filename + ".png"
                        final_format = params['format']
                        params['format'] = 'png'
                        new_fname = os.path.join(tmp_dir, png_file)
                        fig.savefig(new_fname, **params)
                        with Image.open(new_fname) as im:
                            rgb_im = im.convert('RGB')
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

    def _show_file_exist_message(self, fname, new_name):
        self._msg_box.setText(
            "<p>File \"{0}\" already exists. Do you want to overwrite the existing, or save to \"{1}\" instead?<\p>".format(
                fname, new_name))
        self._msg_box.exec_()

    def get_save_file_name(self):
        name = os.path.join(
            str(self.ui.comboBox_directory.currentText()),
            str(self.ui.comboBox_fileName.currentText())
        )
        return name

    def get_savefig_params(self):
        params = {
            'dpi': self.ui.spinBox_dpi.value(),
            'orientation': self.get_orientation(),
            'format': self.get_file_format()[0],
            'transparent': self.get_transparent(),
            'bbox_inches': self.get_bbox_inches()
        }
        return params

    def get_transparent(self):
        return self.ui.checkBox_transparent.isEnabled() and self.ui.checkBox_transparent.isChecked()

    def get_file_format(self):
        return IMAGE_FORMATS[self.ui.comboBox_fileType.currentIndex()]

    def get_bbox_inches(self):
        return 'tight' if self.ui.checkBox_tightBbox.isChecked() else None

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


def generate_unique_file_name(basename):
    name, ext = os.path.splitext(basename)
    counter = 1  # start from 2 as the first extension
    new_name = basename
    while os.path.exists(new_name):
        counter += 1
        new_name = "{0}_{1}{2}".format(name, counter, ext)
    return new_name
