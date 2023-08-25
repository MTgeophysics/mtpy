# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:57:21 2016

@author: jpeacock

Gui to rotate edi files
"""

# ==============================================================================
# Imports
# ==============================================================================
# standard packages
from pathlib import Path
import sys

# 3rd party packages
try:
    from PyQt5 import QtCore, QtGui, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")

from mtpy.gui.my_stream import MyStream
from mtpy import MT

# ==============================================================================
class Rotate_EDI_Files(QtWidgets.QWidget):
    """ """

    def __init__(self):
        super(Rotate_EDI_Files, self).__init__()

        self.edi_list = None
        self.rotation_angle = 0
        self.cwd = Path().cwd()
        self.save_dir = None

        self.setup_ui()

    def setup_ui(self):
        self.setWindowTitle("Rotate EDI Files")

        self.cwd_button = QtWidgets.QPushButton("EDI Folder")
        self.cwd_button.pressed.connect(self.get_cwd)
        self.cwd_edit = QtWidgets.QLineEdit()
        self.cwd_edit.editingFinished.connect(self.set_cwd)

        self.save_dir_button = QtWidgets.QPushButton("Save EDI Folder")
        self.save_dir_button.pressed.connect(self.get_save_path)
        self.save_dir_edit = QtWidgets.QLineEdit()
        self.save_dir_edit.editingFinished.connect(self.set_save_path)

        self.edi_list_box = QtWidgets.QListWidget()
        self.edi_list_box.setMaximumWidth(100)

        self.rotation_label = QtWidgets.QLabel("Rotation Angle (deg)")
        self.rotation_help = QtWidgets.QLabel("Angle is assuming N=0, E=90")
        self.rotation_edit = QtWidgets.QLineEdit(
            "{0:.2f}".format(self.rotation_angle)
        )
        self.rotation_edit.editingFinished.connect(self.set_rotation_angle)

        self.rotate_button = QtWidgets.QPushButton("Rotate")
        self.rotate_button.pressed.connect(self.rotate_edi_files)

        self.output_box = QtWidgets.QTextEdit()
        self.output_box.setMinimumWidth(700)

        self.my_stream = MyStream()
        self.my_stream.message.connect(self.normal_output)

        sys.stdout = self.my_stream

        # --> layout
        cwd_layout = QtWidgets.QHBoxLayout()
        cwd_layout.addWidget(self.cwd_button)
        cwd_layout.addWidget(self.cwd_edit)

        save_layout = QtWidgets.QHBoxLayout()
        save_layout.addWidget(self.save_dir_button)
        save_layout.addWidget(self.save_dir_edit)

        rot_layout = QtWidgets.QHBoxLayout()
        rot_layout.addWidget(self.rotation_label)
        rot_layout.addWidget(self.rotation_edit)
        rot_layout.addWidget(self.rotation_help)

        rot_box = QtWidgets.QVBoxLayout()
        rot_box.addLayout(rot_layout)
        rot_box.addWidget(self.rotate_button)

        right_layout = QtWidgets.QVBoxLayout()
        right_layout.addLayout(cwd_layout)
        right_layout.addLayout(save_layout)
        right_layout.addLayout(rot_box)
        right_layout.addWidget(self.output_box)

        final_layout = QtWidgets.QHBoxLayout()
        final_layout.addWidget(self.edi_list_box)
        final_layout.addLayout(right_layout)

        self.setLayout(final_layout)

    def make_save_path(self):
        rot_str = f"Rotated_{self.rotation_angle:.0f}_deg"
        return self.cwd.joinpath(rot_str.replace("-", "m"))

    def get_cwd(self):
        dir_dialog = QtWidgets.QFileDialog()
        self.cwd = Path(
            str(
                dir_dialog.getExistingDirectory(caption="Choose EDI Directory")
            )
        )
        self.cwd_edit.setText(self.cwd.as_posix())
        self.set_edi_list()
        self.save_dir = self.make_save_path()
        self.save_dir_edit.setText(self.save_dir.as_posix())

    def set_cwd(self):
        self.cwd = Path(str(self.cwd_edit.text()))
        if not self.cwd.exists():
            print(f"Path does not exist {self.cwd}")

    def get_save_path(self):
        dir_dialog = QtWidgets.QFileDialog()
        self.save_dir = Path(
            str(
                dir_dialog.getExistingDirectory(
                    caption="Choose EDI Save Directory"
                )
            )
        )
        self.save_dir_edit.setText(self.save_dir.as_posix())
        if not self.save_dir.exists():
            self.save_dir.mkdir()

    def set_save_path(self):
        self.save_dir = Path(str(self.save_dir_edit.text()))
        if not self.save_dir.exists():
            self.save_dir.mkdir()

    def set_edi_list(self):
        self.edi_list_box.clear()
        self.edi_list = [edi.stem for edi in self.cwd.glob("*.edi")]
        self.edi_list_box.addItems(self.edi_list)

    def set_rotation_angle(self):
        self.rotation_angle = float(str(self.rotation_edit.text()))
        self.rotation_edit.setText(f"{self.rotation_angle:.2f}")
        self.save_dir = self.make_save_path()
        self.save_dir_edit.setText(self.save_dir.as_posix())

    def rotate_edi_files(self):
        if not self.save_dir.exists():
            self.save_dir.mkdir()
            print(f"Made directory {self.save_dir}")

        for edi in self.edi_list:
            print("=" * 40)
            edi_fn = self.cwd.joinpath(f"{edi}.edi")
            mt_obj = MT(fn=edi_fn)
            mt_obj.read()
            mt_obj.rotate(self.rotation_angle, inplace=True)
            mt_obj.write(fn=self.save_dir.joinpath(f"{edi}.edi"))

    @QtCore.pyqtSlot(str)
    def normal_output(self, message):
        self.output_box.moveCursor(QtGui.QTextCursor.End)
        self.output_box.insertPlainText(message)


# ==============================================================================
# create main
# ==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = Rotate_EDI_Files()
    ui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
