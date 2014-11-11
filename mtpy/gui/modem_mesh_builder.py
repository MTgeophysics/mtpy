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
from mtpy.gui.get_edi_files import Get_EDI_Files

class ModEM_Mesh_Window(QtGui.QMainWindow):
    """
    main window for building a mesh for ModEM
    
    """
    
    def __init__(self):
        super(ModEM_Mesh_Window, self).__init__()
        
        self.period_list = []
        self.period_dict = {}
        
        self.modem_model = modem.Model()
        
        self.ui_setup()
        
    def ui_setup(self):
        
        self.setWindowTitle('Make ModEM Mesh')
        self.resize(1920, 1080)

        #-------------- MENU BAR ---------------------------------
        # add a menu bar to the top of the window        
        menu_bar = QtGui.QMenuBar(self)
        menu_bar.setGeometry(QtCore.QRect(0, 0, 1920, 38))
        menu_bar.setObjectName("menu_bar")
        
        # --> Data
        menu_data = QtGui.QMenu(menu_bar)
        menu_data.setTitle("Data File")
        
        # create an action to the data tab
        action_data_open = QtGui.QAction(self)
        action_data_open.setText('Open')
        action_data_open.triggered.connect(self.get_data_fn)

        # add the action to the data tap
        menu_data.addAction(action_data_open)
        menu_bar.addAction(menu_data.menuAction())
        
        #----------> Model -------------------------------------
        menu_model = QtGui.QMenu(menu_bar)
        menu_model.setTitle("Model File")
        
        # create an action to the model tab
        action_model_open = QtGui.QAction(self)
        action_model_open.setText('Open')
        action_model_open.triggered.connect(self.get_model_fn)
        
        action_model_save = QtGui.QAction(self)
        action_model_save.setText('Save')
        action_model_save.triggered.connect(self.save_model_fn)

        # add the action to the model tap
        menu_model.addAction(action_model_open)
        menu_bar.addAction(menu_model.menuAction())
        
        self.setMenuBar(menu_bar)
        
        #------------------Layout of main window------------------------------
        #need to have a central widget 
        self.central_widget = QtGui.QWidget()
        self.setCentralWidget(self.central_widget)
        
        #--> list widget
        self.list_qmodel = QtGui.QStandardItemModel()
        # make a widget that will be the station list
        self.list_widget = QtGui.QListWidget()
        self.list_widget.itemSelectionChanged
        self.list_widget.itemClicked.connect(self.locate_station)
        self.list_widget.setMaximumWidth(150)  
        
        #Load edi button
        load_edi_button = QtGui.QPushButton('Load .edi Files')
        load_edi_button.clicked.connect(self.load_edi_files)
        
        edi_layout = QtGui.QVBoxLayout()
        edi_layout.addWidget(self.list_widget)
        edi_layout.addWidget(load_edi_button)
        
        #add edit fields for mesh properties
        cell_size_label = QtGui.QLabel('Cell Size [East, North] (m)')
        cell_size_east_edit = QtGui.QLineEdit()
        cell_size_east_edit.setText('{0:.2f}'.format(self.modem_model.cell_size_east))
        cell_size_east_edit.textChanged[str].connect(self.set_cell_size_east)
        
        cell_size_north_edit = QtGui.QLineEdit()
        cell_size_north_edit.setText('{0:.2f}'.format(self.modem_model.cell_size_north))
        cell_size_north_edit.textChanged[str].connect(self.set_cell_size_north)
        
        cell_size_grid = QtGui.QHBoxLayout()
        cell_size_grid.addWidget(cell_size_label)
        cell_size_grid.addWidget(cell_size_east_edit)
        cell_size_grid.addWidget(cell_size_north_edit)
        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        
#        self.figure.subplots_adjust(left=self.subplot_left,
#                                    right=self.subplot_right,
#                                    bottom=self.subplot_bottom,
#                                    top=self.subplot_top,
#                                    hspace=self.subplot_hspace,
#                                    wspace=self.subplot_wspace)
        
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
        mpl_vbox.addLayout(cell_size_grid)
        
        # set the layout the main window
        layout = QtGui.QHBoxLayout()
        layout.addLayout(edi_layout)
        layout.addLayout(mpl_vbox)
        self.central_widget.setLayout(layout)

        #set the geometry of each widget        
        self.list_widget.setObjectName("listWidget")
        self.mpl_widget.setObjectName("mpl_widget")
        self.mpl_widget.updateGeometry()
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_data_fn(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)'))
                                       
        self.modem_data = modem.Data()
        self.modem_data.read_data_file(fn)
        self.modem_data_fn = fn
        
        self.dir_path = os.path.dirname(fn)
        
        self.period_list = sorted(self.modem_data.period_list)
        self.period_dict = dict([('{0:.5f}'.format(key), value) for value, key
                                 in enumerate(self.period_list)])
        
        self.list_widget.clear()
        
        #this will add the station name for each station to the qwidget list
        for period in self.period_list:
            self.list_widget.addItem('{0:.5f}'.format(period))
            
        self.plot_period = self.period_list[0]
        
    
        
    def get_model_fn(self):
        """ 
        read in an existing model file
        """

        fn_dialog = QtGui.QFileDialog() 
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                       filter='*.rho'))
                                       
        self.modem_model = modem.Model()
        self.modem_model.read_model_file(fn)

        self.dir_path = os.path.dirname(fn)
        
    def save_model_fn(self):
        """
        save the current mesh settings to a file
        """
        
        fn_dialog = QtGui.QFileDialog()
        save_fn = str(fn_dialog.getSaveFileName(
                                    caption='Choose ModEM model file',
                                    filter='*.rho'))
                                    
        sv_path = os.path.dirname(save_fn)
        sv_basename = os.path.basename(save_fn)
        self.modem_model.write_model_file(save_path=sv_path,
                                          model_fn_basename=sv_basename)
        
        
    def locate_station(self):
        pass
        
    def get_period(self, widget_item):
        """
        get the station name from the clicked station 
        """
        self.plot_period = str(widget_item.text())
        
    def load_edi_files(self):
        edi_obj = Get_EDI_Files()
        self.edi_list = edi_obj.edi_list
        
        self.modem_model = modem.Model(edi_list=self.edi_list)
        
        self.modem_model.get_station_locations()
        
        #self.make_checkable_list(self.modem_model.station_locations['station'])
        
        for station in self.modem_model.station_locations['station']:
            item = QtGui.QListWidgetItem()
            item.setText(station)
            item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.Checked)
            self.list_widget.addItem(item)
            
    def set_cell_size_east(self, text):
        try:
            self.modem_model.cell_size_east = float(text)
        except ValueError:
            print 'Enter floating point number'
            
    def set_cell_size_north(self, text):
        try:
            self.modem_model.cell_size_north = float(text)
        except ValueError:
            print 'Enter floating point number'
        
        
#def main():
    
def main():
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = ModEM_Mesh_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()
        
                
        

        
        
        
        
        
    