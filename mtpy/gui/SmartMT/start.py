import os
import sys
from PyQt4 import QtCore, QtGui

from main_window import Ui_MainWindow

from mtpy.utils.decorator import deprecated


class StartQt4(QtGui.QMainWindow):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        self.setup_menu()

    def setup_menu(self):
        # connect exit menu
        self.ui.actionExit.triggered.connect(QtGui.qApp.quit)
        self.ui.actionOpen_edi_File.triggered.connect(self.file_dialog)
        self.ui.actionOpen_edi_Folder.triggered.connect(self.folder_dialog)
        # not yet impleneted
        self.ui.actionAbout.triggered.connect(self.dummy_action)
        self.ui.actionClose_Project.triggered.connect(self.dummy_action)
        self.ui.actionExport.triggered.connect(self.dummy_action)
        self.ui.actionFind_Action.triggered.connect(self.dummy_action)
        self.ui.actionHelp.triggered.connect(self.dummy_action)
        self.ui.actionNew_Project.triggered.connect(self.dummy_action)
        self.ui.actionOpen_Project.triggered.connect(self.dummy_action)
        self.ui.actionOptions.triggered.connect(self.dummy_action)
        self.ui.actionSave_as_Project.triggered.connect(self.dummy_action)
        self.ui.actionSave_Project.triggered.connect(self.dummy_action)
        self.ui.actionShow_Data_Collection.triggered.connect(self.dummy_action)

    def file_dialog(self, *args, **kwargs):
        dialog = QtGui.QFileDialog(self)
        dialog.setWindowTitle('Open .edi Files...')
        dialog.setNameFilter('.edi files (*.edi)')
        dialog.setFileMode(QtGui.QFileDialog.ExistingFiles)
        if dialog.exec_() == QtGui.QDialog.Accepted:
            file_list = dialog.selectedFiles()
            # todo load in files
            print(" ".join([str(file_name) for file_name in file_list]))

    def folder_dialog(self, *args, **kwargs):
        dir_name = None
        dialog = QtGui.QFileDialog(self)
        dialog.setWindowTitle("Open .edi Directory...")
        dialog.setFileMode(QtGui.QFileDialog.DirectoryOnly)
        while dir_name is None:
            if dialog.exec_() == QtGui.QDialog.Accepted:
                dir_name = dialog.selectedFiles()[0]
                dir_name = str(dir_name)
                file_list = [os.path.join(dir_name, edi) for edi in os.listdir(dir_name) if edi.endswith("edi")]
                if not file_list:
                    # empty list
                    QtGui.QMessageBox.information(self, "NOTE", "Directory does not contain any .edi file, please select again.")
                    dir_name = None   # will read again
                else:
                    # todo load in files
                    for edi in file_list:
                        print edi
            else:
                break
        # dir_name = QtGui.QFileDialog.getExistingDirectory(self, "Open .edi Directory",
        #                                                   '',
        #                                                   QtGui.QFileDialog.DontUseNativeDialog | QtGui.QFileDialog.ShowDirsOnly)



    @deprecated(
        "This is an dommany action that should be unsed only when the required function hasn't been implemented")
    def dummy_action(self, *args, **kwargs):
        pass


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    smartMT = StartQt4()
    smartMT.show()
    sys.exit(app.exec_())
