# -*- coding: utf-8 -*-
"""
Created on Sun Nov 02 13:47:10 2014

@author: jrpeacock
"""

from PyQt4 import QtCore, QtGui
import mtpy.modeling.modem_new as modem
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import numpy as np
import matplotlib.pyplot as plt
import os
import mtpy.analysis.pt as mtpt
import mtpy.utils.exceptions as mtex
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import mtpy.imaging.mtcolors as mtcl
import mtpy.gui.get_edi_files as mt_get_edi_files

class Ui_Mesh_Window(QtGui.QMainWindow):
    """
    main window for building a mesh for ModEM
    
    """
    
    def __init__(self):
        super(Ui_Mesh_Window, self).__init__()
        
        self.ui_setup()
        
    def ui_setup(self):
        """
        set up the user interface
        """
        self.setWindowTitle("Build a mesh for ModEM")
        self.resize(1920, 1080)
        
        self.mesh_widget = MeshWidget()
        self.central_widget = self.setCentralWidget(self.mesh_widget)
        #self.central_widget.setLayout(self.v_layout)
        
        QtCore.QMetaObject.connectSlotsByName(self)

#==============================================================================
# Mesh widget        
#==============================================================================
class MeshWidget(QtGui.QWidget):
    """
    """
    
    def __init__(self):
        super(MeshWidget, self).__init__()
        self.model_obj = modem.Model()
        
        self.setup_ui()
        
    def setup_ui(self):
        
        edi_button = QtGui.QPushButton('Get EDI Files')
        edi_button.clicked.connect(self.get_edi_files)
        
        parameters_label = QtGui.QLabel('Mesh Parameters')
        header_font = QtGui.QFont()
        header_font.setBold = True
        header_font.setPointSize (16)
        parameters_label.setFont(header_font)

        # cell size        
        self.cell_size_label = QtGui.QLabel('Cell Size [E, N] (m)')
        
        self.cell_size_edit_east = QtGui.QLineEdit()
        self.cell_size_edit_east.setText('{0:.2f}'.format(self.model_obj.cell_size_east))
        self.cell_size_edit_east.editingFinished.connect(self.set_cell_size_east_return)
        
        self.cell_size_edit_north = QtGui.QLineEdit()
        self.cell_size_edit_north.setText('{0:.2f}'.format(self.model_obj.cell_size_north))
        self.cell_size_edit_north.editingFinished.connect(self.set_cell_size_north_return)
            
        # cell padding
        self.cell_pad_label = QtGui.QLabel('# of Pad cells [E, N, V]')
        self.cell_pad_east_edit = QtGui.QLineEdit()
        self.cell_pad_east_edit.setText('{0:.0f}'.format(self.model_obj.pad_east))
        self.cell_pad_east_edit.editingFinished.connect(self.set_cell_pad_east)
        
        self.cell_pad_north_edit = QtGui.QLineEdit()
        self.cell_pad_north_edit.setText('{0:.0f}'.format(self.model_obj.pad_north))
        self.cell_pad_north_edit.editingFinished.connect(self.set_cell_pad_north)
        
        self.cell_pad_z_edit = QtGui.QLineEdit()
        self.cell_pad_z_edit.setText('{0:.0f}'.format(self.model_obj.pad_z))
        self.cell_pad_z_edit.editingFinished.connect(self.set_cell_pad_z)
        
        self.grid_layout = QtGui.QGridLayout()
        
        self.grid_layout.addWidget(edi_button, 1, 0)
        self.grid_layout.addWidget(parameters_label, 2, 0)
        
        self.grid_layout.addWidget(self.cell_size_label, 3, 0)
        self.grid_layout.addWidget(self.cell_size_edit_east, 3, 1)
        self.grid_layout.addWidget(self.cell_size_edit_north, 3, 2)
        
        self.grid_layout.addWidget(self.cell_pad_label, 4, 0)
        self.grid_layout.addWidget(self.cell_pad_east_edit, 4, 1)
        self.grid_layout.addWidget(self.cell_pad_north_edit, 4, 2)
        self.grid_layout.addWidget(self.cell_pad_z_edit, 4, 3)
        
        
        self.mpl_widget = MeshPlot()
        self.h_layout = QtGui.QHBoxLayout()
        self.h_layout.addLayout(self.grid_layout)
        self.h_layout.addWidget(self.mpl_widget)
        
        self.setLayout(self.h_layout)
        
    def get_edi_files(self):
        """
        get a list of edi files
        """
        self.dialog_box = QtGui.QFileDialog()
        fn_list = self.dialog_box.getOpenFileNames(
                                    caption='Choose EDI Files',
                                    filter='*.edi')
        self.edi_list = []                                    
        for fn in fn_list:
            self.edi_list.append(str(fn))
            print os.path.basename(str(fn))
         
    def set_cell_size_east_return(self):
        self.model_obj.cell_size_east = float(str(self.cell_size_edit_east.text()))
        self.cell_size_edit_east.setText('{0:.2f}'.format(self.model_obj.cell_size_east))
        print self.model_obj.cell_size_east
            
    def set_cell_size_north_return(self):
        self.model_obj.cell_size_north = float(str(self.cell_size_edit_north.text()))
        self.cell_size_edit_north.setText('{0:.2f}'.format(self.model_obj.cell_size_north))
        print self.model_obj.cell_size_north
        
    def set_cell_pad_east(self):
        self.model_obj.pad_east = float(str(self.cell_pad_east_edit.text()))
        self.cell_pad_east_edit.setText('{0:.2f}'.format(self.model_obj.cell_pad_east))
        print self.model_obj.pad_east
        
    def set_cell_pad_north(self):
        self.model_obj.pad_north = float(str(self.cell_pad_north_edit.text()))
        self.cell_pad_north_edit.setText('{0:.2f}'.format(self.model_obj.cell_pad_north))
        print self.model_obj.pad_north
        
    def set_cell_pad_z(self):
        self.model_obj.pad_z = float(str(self.cell_pad_z_edit.text()))
        self.cell_pad_z_edit.setText('{0:.2f}'.format(self.model_obj.cell_pad_z))
        print self.model_obj.pad_z
#==============================================================================
# Mesh Plot
#==============================================================================
class MeshPlot(QtGui.QWidget):
    """
    plotting the mesh
    """    
    
    def __init__(self):
        super(MeshPlot, self).__init__()
        
        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05        
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        
        self.figure.subplots_adjust(left=self.subplot_left,
                                    right=self.subplot_right,
                                    bottom=self.subplot_bottom,
                                    top=self.subplot_top,
                                    hspace=self.subplot_hspace,
                                    wspace=self.subplot_wspace)
                                    
        #make sure the figure takes up the entire plottable space
        self.mpl_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                     QtGui.QSizePolicy.Expanding)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)
        # set the layout for the plot
        mpl_vbox = QtGui.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)
        
        self.setLayout(mpl_vbox)
        
        self.mpl_widget.updateGeometry()
    
def main():
#if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = Ui_Mesh_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()       
    