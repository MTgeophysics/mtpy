# -*- coding: utf-8 -*-
"""
    Description:

    Usage:

    Author: YingzhiGou
    Date: 24/07/2017
"""
import os

from PyQt4 import QtGui, QtCore
import matplotlib.pyplot as plt

from mtpy.gui.SmartMT.ui_asset.dialog_export import Ui_Dialog_Export


class ExportDialog(QtGui.QDialog):
    def __init__(self, parent=None):
        QtGui.QDialog.__init__(self, parent)
        self.ui = Ui_Dialog_Export()
        self.ui.setupUi(self)

        # setup file types
        self._formats = []
        filetypes = plt.gcf().canvas.get_supported_filetypes()  # this may need to be set everytime
        for type, description in filetypes.iteritems():
            self.ui.comboBox_fileType.addItem("%s (.%s)" % (description, type))
            self._formats.append((type, description))

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

    def _file_type_changed(self, *args, **kwargs):
        index = self.ui.comboBox_fileType.currentIndex()
        ext, _ = self._formats[index]
        filename = str(self.ui.comboBox_fileName.currentText())
        filename, _ = os.path.splitext(filename)

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
        index = [i for i, (ext, dsc) in enumerate(self._formats) if ext == extension]
        if index:
            self.ui.comboBox_fileType.setCurrentIndex(index[0])

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

    def show_dialog(self, fig):
        self.show()

    def keyPressEvent(self, event):
        """Capture and ignore all key press events.

        This is used so that return key event does not trigger any button
        from the dialog. We need to allow the return key to be used in filters
        in the widget."""
        if event.key() == QtCore.Qt.Key_Escape:
            # call reject if Escape is pressed.
            self.reject()
        pass
