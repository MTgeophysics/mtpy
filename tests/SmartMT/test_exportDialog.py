import os
from unittest import TestCase

import matplotlib.pyplot as plt
import numpy as np
from qtpy import QtCore
from qtpy.QtWidgets import QFileDialog, QMessageBox, QDialog
from qtpy.QtTest import QTest

from mtpy.gui.SmartMT.gui.export_dialog import ExportDialog, IMAGE_FORMATS
from tests import TEST_TEMP_DIR, plt_wait
from tests.imaging import reset_matplotlib


def _fake_exec_accept():
    return QFileDialog.Accepted


def _fake_exec_reject():
    return QFileDialog.Rejected


def _rewrite_text(widget, text, modifier=QtCore.Qt.NoModifier):
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_A, QtCore.Qt.ControlModifier)
    QTest.keyClicks(widget, text, modifier=modifier)
    QTest.keyEvent(QTest.Click, widget, QtCore.Qt.Key_Enter)


def _create_fig():
    t = np.arange(0.0, 2.0, 0.01)
    s = 1 + np.sin(2 * np.pi * t)
    plt.plot(t, s)
    plt.xlabel('time (s)')
    plt.ylabel('voltage (mV)')
    plt.title('About as simple as it gets, folks')
    plt.grid(True)
    # plt.savefig("test.png")
    # plt.show()
    plt_wait(1)
    return plt.gcf()  # get access to the current fig


class TestExportDialog(TestCase):
    @classmethod
    def setUpClass(cls):
        # setup temp dir
        cls._temp_dir = TEST_TEMP_DIR
        reset_matplotlib()

    def setUp(self):
        # create figure
        self._fig = _create_fig()
        # create GUI
        self.dialog = ExportDialog()
        self.dialog.show()
        QTest.qWaitForWindowActive(self.dialog)

    def tearDown(self):
        self.dialog.close()
        plt.close(self._fig)

    def test_defaults(self):
        """ test gui default state"""
        # check row states
        self.assertTrue(self.dialog.ui.comboBox_fileName.currentText() == "figure.png", "Default File Name")
        self.assertTrue(self.dialog.ui.comboBox_directory.currentText() == os.path.expanduser("~"), "Default Path")
        # file type
        self.assertTrue(set(["{} (.{})".format(desc, ext)
                             for ext, desc in self._fig.canvas.get_supported_filetypes().iteritems()]) ==
                        set([str(self.dialog.ui.comboBox_fileType.itemText(i))
                             for i in range(self.dialog.ui.comboBox_fileType.count())]),
                        "Supported Formats")

        self.assertTrue(self.dialog.ui.checkBox_tightBbox.isChecked(), "Tight Layout Default")
        self.assertFalse(self.dialog.get_transparent(), "Transparent Default")
        self.assertTrue(self.dialog.ui.comboBox_orientation.currentText() == "Landscape", "Orientation Default")
        self.assertTrue(self.dialog.ui.spinBox_dpi.value() == 80)
        self.assertTrue(self.dialog.ui.doubleSpinBox_height_inches.value() == 6.)
        self.assertTrue(self.dialog.ui.doubleSpinBox_width_inches.value() == 8.)
        self.assertTrue(self.dialog.ui.spinBox_height_pixels.value() == 480)
        self.assertTrue(self.dialog.ui.spinBox_width_pixels.value() == 640)
        self.assertTrue(self.dialog.ui.checkBox_open_after_export.isChecked())

        # check states from the getters
        self.assertTrue(self.dialog.get_bbox_inches() == 'tight', "Tight Layout Value")
        self.assertTrue(self.dialog.get_file_format()[0] == 'png', "Format Value")
        self.assertTrue(self.dialog.get_orientation() == 'landscape', "Orientation Value")
        self.assertTrue(os.path.normpath(self.dialog.get_save_file_name()) ==
                        os.path.normpath(os.path.join(os.path.expanduser("~"),
                                                      str(self.dialog.ui.comboBox_fileName.currentText()))
                                         ),
                        "Save File Path Value")

    def test_file_name_change(self):
        # select all existing tests
        _rewrite_text(self.dialog.ui.comboBox_fileName, "test_file.jpg")
        # current text should have changed
        self.assertTrue(self.dialog.ui.comboBox_fileName.currentText() == "test_file.jpg",
                        "Changed file name")
        # format should have changed
        self.assertTrue(self.dialog.get_file_format()[0] == "jpg")
        # transparent should be false
        self.assertFalse(self.dialog.get_transparent(), "transparent")

        # change to file with unsupported format
        _rewrite_text(self.dialog.ui.comboBox_fileName, "test_file_2.abcd")
        # current text should have changed
        self.assertTrue(self.dialog.ui.comboBox_fileName.currentText() == "test_file_2.abcd",
                        "Changed file name")
        # current format should not been changed
        self.assertTrue(self.dialog.get_file_format()[0] == "jpg")

    def test_file_type_change(self):
        for i in range(self.dialog.ui.comboBox_fileType.count()):
            self.dialog.ui.comboBox_fileType.setCurrentIndex(i)
            extenion = self.dialog.get_file_format()[0]
            self.assertTrue(self.dialog.ui.comboBox_fileName.currentText() == "figure.{}".format(extenion))

    def test_directory_change(self):
        _rewrite_text(self.dialog.ui.comboBox_directory, os.path.abspath(self._temp_dir))
        self.assertTrue(os.path.dirname(self.dialog.get_save_file_name()) == os.path.abspath(self._temp_dir))
        # print self.dialog.get_save_file_name()

        # select from the browse

        self.dialog._dir_dialog.setDirectory(os.path.normpath(os.path.expanduser("~")))
        self.dialog._dir_dialog.exec_ = _fake_exec_reject  # path should not change
        QTest.mouseClick(self.dialog.ui.pushButton_browse, QtCore.Qt.LeftButton)
        self.assertTrue(os.path.dirname(self.dialog.get_save_file_name()) == os.path.abspath(self._temp_dir))

        self.dialog._dir_dialog.exec_ = _fake_exec_accept
        QTest.mouseClick(self.dialog.ui.pushButton_browse, QtCore.Qt.LeftButton)
        # QTest.qWaitForWindowShown(self.dialog._dir_dialog)
        # self.dialog._dir_dialog.accept()
        self.assertTrue(
            os.path.dirname(
                os.path.normpath(self.dialog.get_save_file_name())
            ) == os.path.normpath(os.path.expanduser("~")))

    def test_export(self):
        # set export dir
        _rewrite_text(self.dialog.ui.comboBox_directory,
                      os.path.abspath(self._temp_dir))
        fname = self.dialog.get_save_file_name()
        if os.path.isfile(fname):
            # if file exist, remove
            os.remove(fname)
        self.assertFalse(os.path.exists(fname), "File exists")

        # set open after to false
        self.dialog.ui.checkBox_open_after_export.setChecked(False)

        self.dialog.exec_ = self._fake_export_dialog_exec_cancel  # should not create file
        self.dialog._msg_box.exec_ = self._fake_msg_dialog_exec_cancel
        fname = self.dialog.export_to_file(self._fig)
        print self._fig.get_dpi(), self.dialog.ui.spinBox_dpi.value()
        self.assertTrue(self.dialog.ui.spinBox_dpi.value() == self._fig.get_dpi())
        self.assertTrue(fname is None)
        self.assertFalse(os.path.exists(self.dialog.get_save_file_name()), "File exists")

        # save the new file now
        self.dialog.exec_ = self._fake_export_dialog_exec_export
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        self.assertTrue(os.path.isfile(fname))

        file_count = len([name for name in os.listdir(self._temp_dir)
                          if os.path.isfile(os.path.join(self._temp_dir, name))])
        # save to the same file and overwrite
        self.dialog._msg_box.exec_ = self._fake_msg_dialog_exec_overwrite
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        new_file_count = len([name for name in os.listdir(self._temp_dir)
                              if os.path.isfile(os.path.join(self._temp_dir, name))])
        self.assertTrue(file_count == new_file_count)  # no file should be created
        # save to the same file and save as new name
        self.dialog._msg_box.exec_ = self._fake_msg_dialog_exec_save_as
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        new_file_count = len([name for name in os.listdir(self._temp_dir)
                              if os.path.isfile(os.path.join(self._temp_dir, name))])
        self.assertTrue(file_count + 1 == new_file_count)  # one extra file should be created
        file_count = new_file_count

    def test_dpi(self):
        # save to higher dpi
        # set export dir
        _rewrite_text(self.dialog.ui.comboBox_directory,
                      os.path.abspath(self._temp_dir))
        self.dialog.exec_ = self._fake_export_dialog_exec_export
        self.dialog._msg_box.exec_ = self._fake_msg_dialog_exec_overwrite
        # set open after to false
        self.dialog.ui.checkBox_open_after_export.setChecked(False)

        QTest.keyClicks(self.dialog.ui.spinBox_dpi, '400')
        _rewrite_text(self.dialog.ui.comboBox_fileName, "400dpi.jpg")
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        new_file_count = len([name for name in os.listdir(self._temp_dir)
                              if os.path.isfile(os.path.join(self._temp_dir, name))])

        QTest.keyClicks(self.dialog.ui.spinBox_dpi, '600')
        _rewrite_text(self.dialog.ui.comboBox_fileName, "600dpi.jpg")
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        new_file_count = len([name for name in os.listdir(self._temp_dir)
                              if os.path.isfile(os.path.join(self._temp_dir, name))])

        QTest.keyClicks(self.dialog.ui.spinBox_dpi, '1000')
        _rewrite_text(self.dialog.ui.comboBox_fileName, "1000dpi.jpg")
        fname = self.dialog.export_to_file(self._fig)
        self.assertTrue(os.path.exists(fname), "File exists")
        new_file_count = len([name for name in os.listdir(self._temp_dir)
                              if os.path.isfile(os.path.join(self._temp_dir, name))])

    def _fake_msg_dialog_exec_overwrite(self):
        self.dialog._msg_box.show()
        QTest.qWaitForWindowActive(self.dialog._msg_box)
        QTest.mouseClick(self.dialog._msg_box_button_overwrite, QtCore.Qt.LeftButton)
        return QMessageBox.Accepted

    def _fake_msg_dialog_exec_save_as(self):
        self.dialog._msg_box.show()
        QTest.qWaitForWindowActive(self.dialog._msg_box)
        QTest.mouseClick(self.dialog._msg_box_button_save_as, QtCore.Qt.LeftButton)
        return QMessageBox.Accepted

    def _fake_msg_dialog_exec_cancel(self):
        self.dialog._msg_box.show()
        QTest.qWaitForWindowActive(self.dialog._msg_box)
        QTest.mouseClick(self.dialog._msg_box_button_cancel, QtCore.Qt.LeftButton)
        return QMessageBox.Cancel

    def _fake_export_dialog_exec_cancel(self):
        QTest.mouseClick(self.dialog.ui.pushButton_cancel, QtCore.Qt.LeftButton)
        return QDialog.Rejected

    def _fake_export_dialog_exec_export(self):
        QTest.mouseClick(self.dialog.ui.pushButton_export, QtCore.Qt.LeftButton)
        return QDialog.Accepted


def _transparent_test_gen(index, ext, description):
    def _test_transparent(self):
        # set to save to tmp dir
        _rewrite_text(self.dialog.ui.comboBox_directory,
                      os.path.abspath(self._temp_dir))

        self.dialog.exec_ = self._fake_export_dialog_exec_export
        self.dialog._msg_box.exec_ = self._fake_msg_dialog_exec_overwrite
        # set open after to false
        self.dialog.ui.checkBox_open_after_export.setChecked(False)

        # print "testing save to {0[1]} (.{0[0]})".format(self.dialog.get_file_format())
        for isTrans in [True, False]:
            _rewrite_text(self.dialog.ui.comboBox_fileName, "transparent_{}.{}".format(isTrans, ext))
            self.dialog.ui.comboBox_fileType.setCurrentIndex(index)
            self.assertTrue((ext, description) == self.dialog.get_file_format(), "sanity check")
            self.dialog.ui.checkBox_transparent.setChecked(isTrans)
            try:
                fname = self.dialog.export_to_file(self._fig)
            except RuntimeError as e:
                self.skipTest(e.message)
            self.assertTrue(os.path.exists(fname),
                            "testing save to {0[1]} (.{0[0]}) without transparent".format(
                                self.dialog.get_file_format()))

    return _test_transparent


# generate tests
for index, (ext, description) in enumerate(IMAGE_FORMATS):
    _test = _transparent_test_gen(index, ext, description)
    _test.__name__ = "test_transparent_{}".format(ext)
    setattr(TestExportDialog, _test.__name__, _test)
