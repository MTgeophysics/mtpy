# -*- coding: utf-8 -*-
"""
Created on Sun Nov 02 14:31:02 2014

@author: jrpeacock

this module is used by mtpy/gui/modem_mesh_builder.py
"""

from PyQt4 import QtGui
import sys


class Get_EDI_Files(QtGui.QWidget):

    def __init__(self):
        super(Get_EDI_Files, self).__init__()
        self.dialog_box = QtGui.QFileDialog()
        fn_list = self.dialog_box.getOpenFileNames(
            caption='Choose EDI Files',
            filter='*.edi')
        self.edi_list = []
        for fn in fn_list:
            self.edi_list.append(str(fn))


def main():
    app = QtGui.QApplication(sys.argv)
    window = Get_EDI_Files()
    sys.exit(app.exec_())

# if __name__ == '__main__':

#    import sys
#    app = QtGui.QApplication(sys.argv)
#    window = Get_EDI_Files()
#    sys.exit(app.exec_())
