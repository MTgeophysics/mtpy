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
        
        make_mesh_button = QtGui.QPushButton('Plot Mesh')
        make_mesh_button.clicked.connect(self.plot_mesh)
        
        edi_layout = QtGui.QVBoxLayout()
        edi_layout.addWidget(self.list_widget)
        edi_layout.addWidget(load_edi_button)
        edi_layout.addWidget(make_mesh_button)
        
        #add edit fields for mesh properties
        cell_size_label = QtGui.QLabel('Cell Size [East, North] (m)')
        cell_size_east_edit = QtGui.QLineEdit()
        cell_size_east_edit.setText('{0:.2f}'.format(self.modem_model.cell_size_east))
        cell_size_east_edit.textEdited[str].connect(self.set_cell_size_east)
        
        cell_size_north_edit = QtGui.QLineEdit()
        cell_size_north_edit.setText('{0:.2f}'.format(self.modem_model.cell_size_north))
        cell_size_north_edit.textEdited[str].connect(self.set_cell_size_north)
        
        cell_size_grid = QtGui.QHBoxLayout()
        cell_size_grid.addWidget(cell_size_label)
        cell_size_grid.addWidget(cell_size_east_edit)
        cell_size_grid.addWidget(cell_size_north_edit)
        
        pad_label = QtGui.QLabel('# of padding cells [East, North, Z]')
        pad_east_edit = QtGui.QLineEdit()
        pad_east_edit.setText('{0:.0f}'.format(self.modem_model.pad_east))
        pad_east_edit.textEdited[str].connect(self.set_pad_east)
        
        pad_north_edit = QtGui.QLineEdit()
        pad_north_edit.setText('{0:.0f}'.format(self.modem_model.pad_north))
        pad_north_edit.textEdited[str].connect(self.set_pad_north)
        
        pad_z_edit = QtGui.QLineEdit()
        pad_z_edit.setText('{0:.0f}'.format(self.modem_model.pad_z))
        pad_z_edit.textEdited[str].connect(self.set_pad_z)
        
        #add to cell size grid
        cell_size_grid.addWidget(pad_label)
        cell_size_grid.addWidget(pad_east_edit)
        cell_size_grid.addWidget(pad_north_edit)
        cell_size_grid.addWidget(pad_z_edit)
        
        pad_stretch_label = QtGui.QLabel('Stretch factor [H, V]')
        pad_stretch_h_edit = QtGui.QLineEdit()
        pad_stretch_h_edit.setText('{0:.2f}'.format(self.modem_model.pad_stretch_h))
        pad_stretch_h_edit.textEdited[str].connect(self.set_pad_stretch_h)
               
        pad_stretch_v_edit = QtGui.QLineEdit()
        pad_stretch_v_edit.setText('{0:.2f}'.format(self.modem_model.pad_stretch_v))
        pad_stretch_v_edit.textEdited[str].connect(self.set_pad_stretch_v)

        #add to cell size grid        
        cell_size_grid.addWidget(pad_stretch_label)
        cell_size_grid.addWidget(pad_stretch_h_edit)
        cell_size_grid.addWidget(pad_stretch_v_edit)
        
        num_z_layers_label = QtGui.QLabel('No. Z-layers')
        num_z_layers_edit = QtGui.QLineEdit()
        num_z_layers_edit.setText('{0:.0f}'.format(self.modem_model.n_layers))
        num_z_layers_edit.textEdited[str].connect(self.set_num_z_layers)
        
        z1_layer_label = QtGui.QLabel('Z1 Layer (m)')
        z1_layer_edit = QtGui.QLineEdit()
        z1_layer_edit.setText('{0:.2f}'.format(self.modem_model.z1_layer))
        z1_layer_edit.textEdited[str].connect(self.set_z1_layer)
        
        z_target_label = QtGui.QLabel('Target Layer (m)')
        z_target_edit = QtGui.QLineEdit()
        z_target_edit.setText('{0:.2f}'.format(self.modem_model.z_target_depth))
        z_target_edit.textEdited[str].connect(self.set_z_target)
        
        z_bottom_label = QtGui.QLabel('Bottom Layer (m)')
        z_bottom_edit = QtGui.QLineEdit()
        z_bottom_edit.setText('{0:.2f}'.format(self.modem_model.z_bottom))
        z_bottom_edit.textEdited[str].connect(self.set_z_bottom)
        
        z_grid = QtGui.QHBoxLayout()
        z_grid.addWidget(num_z_layers_label)
        z_grid.addWidget(num_z_layers_edit)
        z_grid.addWidget(z1_layer_label)
        z_grid.addWidget(z1_layer_edit)
        z_grid.addWidget(z_target_label)
        z_grid.addWidget(z_target_edit)
        z_grid.addWidget(z_bottom_label)
        z_grid.addWidget(z_bottom_edit)
        
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
        mpl_vbox.addLayout(z_grid)
        
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
        
    def plot_mesh(self):
        """
        plot mesh after a parameter is changed
        """
        
        try: 
            self.modem_model.make_mesh()
        except AttributeError:
            print 'No stations exist to create a mesh'
            
        except modem.ModEMError:
            print 'No stations exist to create a mesh'
        
        station_marker = 'v'
        marker_color = 'b'
        marker_size = 2
        
        line_color = 'k'
        line_width = .5
        
        plt.rcParams['figure.subplot.hspace'] = .3
        plt.rcParams['figure.subplot.wspace'] = .3
        plt.rcParams['figure.subplot.left'] = .12
        plt.rcParams['font.size'] = 10
        
        self.figure.clf()
        self.figure.subplots_adjust(left=.12, right=.98, top=.95, bottom=.12)
        
        #make a rotation matrix to rotate data
        #cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        #sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))
        
        #turns out ModEM has not accomodated rotation of the grid, so for
        #now we will not rotate anything.
        cos_ang = 1
        sin_ang = 0
        
        #--->plot map view    
        ax1 = self.figure.add_subplot(1, 2, 1, aspect='equal')
        
        
        #plot station locations
        plot_east = self.modem_model.station_locations['rel_east']
        plot_north = self.modem_model.station_locations['rel_north']
        
        ax1.scatter(plot_east,
                    plot_north, 
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)
                                
        
        east_line_xlist = []
        east_line_ylist = []   
        north_min = self.modem_model.grid_north.min()         
        north_max = self.modem_model.grid_north.max()         
        for xx in self.modem_model.grid_east:
            east_line_xlist.extend([xx*cos_ang+north_min*sin_ang, 
                                    xx*cos_ang+north_max*sin_ang])
            east_line_xlist.append(None)
            east_line_ylist.extend([-xx*sin_ang+north_min*cos_ang, 
                                    -xx*sin_ang+north_max*cos_ang])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                      east_line_ylist,
                      lw=line_width,
                      color=line_color)

        north_line_xlist = []
        north_line_ylist = [] 
        east_max = self.modem_model.grid_east.max()
        east_min = self.modem_model.grid_east.min()
        for yy in self.modem_model.grid_north:
            north_line_xlist.extend([east_min*cos_ang+yy*sin_ang,
                                     east_max*cos_ang+yy*sin_ang])
            north_line_xlist.append(None)
            north_line_ylist.extend([-east_min*sin_ang+yy*cos_ang, 
                                     -east_max*sin_ang+yy*cos_ang])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                      north_line_ylist,
                      lw=line_width,
                      color=line_color)
        
#        ax1.set_xlim(plot_east.min()-10*self.modem_model.cell_size_east,
#                     plot_east.max()+10*self.modem_model.cell_size_east)
#        
#        ax1.set_ylim(plot_north.min()-10*self.modem_model.cell_size_north,
#                     plot_north.max()+ 10*self.modem_model.cell_size_east)
            
        ax1.set_ylabel('Northing (m)', fontdict={'size':9,'weight':'bold'})
        ax1.set_xlabel('Easting (m)', fontdict={'size':9,'weight':'bold'})
        ax1.axis('tight')
        
        ##----plot depth view
        ax2 = self.figure.add_subplot(1, 2, 2, aspect='auto', sharex=ax1)
        

        #plot the grid 
        east_line_xlist = []
        east_line_ylist = []            
        for xx in self.modem_model.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([0, 
                                    self.modem_model.grid_z.max()])
            east_line_ylist.append(None)
        ax2.plot(east_line_xlist,
                 east_line_ylist,
                 lw=line_width,
                 color=line_color)

        z_line_xlist = []
        z_line_ylist = [] 
        for zz in self.modem_model.grid_z:
            z_line_xlist.extend([self.modem_model.grid_east.min(),
                                     self.modem_model.grid_east.max()])
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        ax2.plot(z_line_xlist,
                 z_line_ylist,
                 lw=line_width,
                 color=line_color)
                      
        
        #--> plot stations
        ax2.scatter(plot_east,
                    [0]*self.modem_model.station_locations.shape[0],
                    marker=station_marker,
                    c=marker_color,
                    s=marker_size)

        

        ax2.set_ylim(self.modem_model.grid_z.max(), -200)

            
#        ax1.set_xlim(plot_east.min()-10*self.modem_model.cell_size_east,
#                     plot_east.max()+10*self.modem_model.cell_size_east)
            
        ax2.set_ylabel('Depth (m)', fontdict={'size':11, 'weight':'bold'})
        ax2.set_xlabel('Easting (m)', fontdict={'size':11, 'weight':'bold'}) 
        
        ax2.axis('tight')
        
        # --> be sure to draw the plot
        self.mpl_widget.draw()
        
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
                                          
        self.ui_setup()
        
        
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
            #item.setFlags(item.flags() | QtCore.Qt.ItemIsUserCheckable)
            #item.setCheckState(QtCore.Qt.Checked)
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
            
            
    def set_pad_east(self, text):
        self.modem_model.pad_east = int(text)
        
    def set_pad_north(self, text):
        self.modem_model.pad_north = int(text)
        
    def set_pad_z(self, text):
        self.modem_model.pad_z = int(text)
    
    def set_pad_stretch_h(self, text):
        self.modem_model.pad_stretch_h = float(text)
        
    def set_pad_stretch_v(self, text):
        self.modem_model.pad_stretch_v = float(text)
        
    def set_num_z_layers(self, text):
        self.modem_model.n_layers = int(text)
                
    def set_z1_layer(self, text):
        self.modem_model.z1_layer = float(text)
                
    def set_z_target(self, text):
        self.modem_model.z_target_depth = float(text)
        
    def set_z_bottom(self, text):
        self.modem_model.z_bottom = float(text)
       
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
        
                
        

        
        
        
        
        
    