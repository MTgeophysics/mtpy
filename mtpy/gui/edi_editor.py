# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 15:04:27 2016

@author: jpeacock
"""

#==============================================================================
#  Imports
#==============================================================================

from PyQt4 import QtCore, QtGui
import mtpy.modeling.modem_new as modem
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
import os
import copy
import sys

#==============================================================================
# UI
#==============================================================================

class MyStream(QtCore.QObject):
    """
    this class will emit a signal to write standard output to the UI
    """
    message = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(MyStream, self).__init__(parent)

    def write(self, message):
        self.message.emit(str(message))

class EDI_Editor_Window(QtGui.QMainWindow):
    """
    This is the main window for editing an edi file.
    
    Includes:
        * Editing points by removing them
        * Adjusting for static shift
        * Removing distortion
        * Rotations
        * Strike estimations
    """
    
    def __init__(self):
        super(EDI_Editor_Window, self).__init__()
        
        
    def ui_setup(self):
        """
        set up the user interface
        """
        
        self.setWindowTitle('EDI Editor')
        self.resize(1920, 1080)
        
        self.plot_widget = PlotWidget()
        self.centralWidget = self.setCentralWidget(self.plot_widget)
        
        self.menu_file = self.menuBar().addMenu(self.tr("&File"))
        
        self.action_open_file = self.menu_file.addAction(self.tr("&Open"))
        self.action_open_file.triggered.connect(self.get_edi_file)
        
        self.action_close_file = self.menu_file.addAction(self.tr("C&lose"))
        self.action_close_file.triggered.connect(self.close_edi_file)
        
        self.action_save_file = self.menu_file.addAction(self.tr("&Save"))
        self.action_save_file.triggered.connect(self.save_edi_file)
        
        #self.my_stream = MyStream()
        #self.my_stream.message.connect(self.plot_widget.normal_output)
        
        #sys.stdout = self.my_stream
        
        #QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_edi_file(self):
        """
        get edi file
        """
        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose EDI file', 
                                           filter='*.edi'))
        
        return fn 

    def close_edi_file(self):
        pass

    def save_edi_file(self):
        pass
                       
#==============================================================================
# Plot Widget     
#==============================================================================
class PlotWidget(QtGui.QWidget):
    """
    matplotlib plot of the data
    
    """
    
    def __init__(self):
        super(PlotWidget, self).__init__()

        self.mt_obj = None
        
    def setup_ui(self):
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        
        self.mpl_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)
                                      
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)
        
        mpl_vbox = QtGui.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)
        
        self.setLayout(mpl_vbox)
        self.mpl_widget.updateGeometry()
        
    #@QtCore.pyqtSlot(str)
    #def normal_output(self, message):
    #    self.output_box.moveCursor(QtGui.QTextCursor.End)
    #    self.output_box.insertPlainText(message)
    
        
#==============================================================================
# Def Main
#==============================================================================
def main():
    app = QtGui.QApplication(sys.argv)
    ui = EDI_Editor_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()