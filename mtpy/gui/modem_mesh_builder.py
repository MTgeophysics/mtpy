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
import sys


class MyStream(QtCore.QObject):
    """
    this class will emit a signal
    """
    message = QtCore.pyqtSignal(str)
    def __init__(self, parent=None):
        super(MyStream, self).__init__(parent)

    def write(self, message):
        self.message.emit(str(message))

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
        """
        set up the user interface
        """
        self.setWindowTitle("Build a mesh for ModEM")
        self.resize(1920, 1080)
        
        self.mesh_widget = MeshWidget()
        self.central_widget = self.setCentralWidget(self.mesh_widget)
        #self.central_widget.setLayout(self.v_layout)
        self.my_stream = MyStream()
        self.my_stream.message.connect(self.mesh_widget.normal_output)
        
        sys.stdout = self.my_stream
        
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
        self.mpl_widget = MeshPlot()
        
        #sys.stdout = MyStream()
        
        
        self.setup_ui()
        
    def setup_ui(self):
        # get edi files
        self.edi_button = QtGui.QPushButton('Get EDI Files')
        self.edi_button.clicked.connect(self.get_edi_files)
        
        self.list_widget = QtGui.QListWidget()
        self.list_widget.itemClicked.connect(self.select_station)
        self.list_widget.setMaximumWidth(200)
        
        self.plot_mesh_button = QtGui.QPushButton('Make Mesh')
        self.plot_mesh_button.clicked.connect(self.plot_mesh)
        
        self.save_mesh_button = QtGui.QPushButton('Save Mesh')
        self.save_mesh_button.clicked.connect(self.save_mesh)
        
        
        self.output_box = QtGui.QTextEdit()
        
        #make a label for the mesh parameters
        self.parameters_label = QtGui.QLabel('Mesh Parameters')
        header_font = QtGui.QFont()
        header_font.setBold = True
        header_font.setPointSize (16)
        self.parameters_label.setFont(header_font)

        # cell size        
        self.cell_size_label_east = QtGui.QLabel('Cell Size East (m)')
        self.cell_size_label_north = QtGui.QLabel('Cell Size North (m)')
        
        self.cell_size_edit_east = QtGui.QLineEdit()
        self.cell_size_edit_east.setText('{0:.2f}'.format(self.model_obj.cell_size_east))
        self.cell_size_edit_east.editingFinished.connect(self.set_cell_size_east)
        
        self.cell_size_edit_north = QtGui.QLineEdit()
        self.cell_size_edit_north.setText('{0:.2f}'.format(self.model_obj.cell_size_north))
        self.cell_size_edit_north.editingFinished.connect(self.set_cell_size_north)
            
        # cell padding
        self.cell_pad_label_east = QtGui.QLabel('# of Pad cells E')
        self.cell_pad_label_north = QtGui.QLabel('# of Pad cells N')
        self.cell_pad_label_z = QtGui.QLabel('# of Pad cells Z')
        
        self.cell_pad_east_edit = QtGui.QLineEdit()
        self.cell_pad_east_edit.setText('{0:.0f}'.format(self.model_obj.pad_east))
        self.cell_pad_east_edit.editingFinished.connect(self.set_cell_pad_east)
        
        self.cell_pad_north_edit = QtGui.QLineEdit()
        self.cell_pad_north_edit.setText('{0:.0f}'.format(self.model_obj.pad_north))
        self.cell_pad_north_edit.editingFinished.connect(self.set_cell_pad_north)
        
        self.cell_pad_z_edit = QtGui.QLineEdit()
        self.cell_pad_z_edit.setText('{0:.0f}'.format(self.model_obj.pad_z))
        self.cell_pad_z_edit.editingFinished.connect(self.set_cell_pad_z)
        
        self.pad_h_label = QtGui.QLabel('Horiz. Padding Factor')
        self.pad_v_label = QtGui.QLabel('Vert. Padding Factor')
        
        self.pad_h_edit = QtGui.QLineEdit()
        self.pad_h_edit.setText('{0:.2f}'.format(self.model_obj.pad_stretch_h))
        self.pad_h_edit.editingFinished.connect(self.set_pad_h)
        
        self.pad_v_edit = QtGui.QLineEdit()
        self.pad_v_edit.setText('{0:.2f}'.format(self.model_obj.pad_stretch_v))
        self.pad_v_edit.editingFinished.connect(self.set_pad_v)
        
        # vertical layer parameters
        self.n_layers_label = QtGui.QLabel('Number of Vertical Layers')
        self.n_layers_edit = QtGui.QLineEdit()
        self.n_layers_edit.setText('{0:.0f}'.format(self.model_obj.n_layers))
        self.n_layers_edit.editingFinished.connect(self.set_n_layers)
        
        self.z1_layer_label = QtGui.QLabel('Thicknes of 1st layer (m)')
        self.z1_layer_edit = QtGui.QLineEdit()
        self.z1_layer_edit.setText('{0:.2f}'.format(self.model_obj.z1_layer))
        self.z1_layer_edit.editingFinished.connect(self.set_z1_layer)
        
        self.z_target_label = QtGui.QLabel('Target Depth (m)')
        self.z_target_edit = QtGui.QLineEdit()
        self.z_target_edit.setText('{0:.2f}'.format(self.model_obj.z_target_depth))
        self.z_target_edit.editingFinished.connect(self.set_z_target)
        
        self.z_bottom_label = QtGui.QLabel('Bottom of the Model (m)')
        self.z_bottom_edit = QtGui.QLineEdit()
        self.z_bottom_edit.setText('{0:.2f}'.format(self.model_obj.z_bottom))
        self.z_bottom_edit.editingFinished.connect(self.set_z_bottom)
        
        # rotation angle
        self.rot_ang_label = QtGui.QLabel('Mesh Rotation (deg)')
        self.rot_ang_hint = QtGui.QLabel('[N=0, E=90]')
        self.rot_ang_edit = QtGui.QLineEdit()
        self.rot_ang_edit.setText('{0:.2f}'.format(self.model_obj.mesh_rotation_angle))
        self.rot_ang_edit.editingFinished.connect(self.set_rotation_angle)
        
        # starting resistivity
        self.rho_start_label = QtGui.QLabel('Starting rho (Ohmm)')
        self.rho_start_edit = QtGui.QLineEdit()
        self.rho_start_edit.setText('{0:.2f}'.format(100))
        self.rho_start_edit.editingFinished.connect(self.set_rho)
        
        #--- Set the layout ----------
        self.edi_grid = QtGui.QVBoxLayout()
        self.edi_grid.addWidget(self.edi_button)
        self.edi_grid.addWidget(self.list_widget)
        
        self.grid_layout = QtGui.QGridLayout()
        
        
        self.grid_layout.addWidget(self.plot_mesh_button, 0, 0)
        self.grid_layout.addWidget(self.save_mesh_button, 0, 1)
        
        self.grid_layout.addWidget(self.parameters_label, 2, 0, 1, 2)

        self.grid_layout.addWidget(self.cell_size_label_east, 3, 0)
        self.grid_layout.addWidget(self.cell_size_edit_east, 3, 1)
        
        self.grid_layout.addWidget(self.cell_size_label_north, 4, 0)
        self.grid_layout.addWidget(self.cell_size_edit_north, 4, 1)
        
        self.grid_layout.addWidget(self.cell_pad_label_east, 5, 0)
        self.grid_layout.addWidget(self.cell_pad_east_edit, 5, 1)
        
        self.grid_layout.addWidget(self.cell_pad_label_north, 6, 0)
        self.grid_layout.addWidget(self.cell_pad_north_edit, 6, 1)
        
        self.grid_layout.addWidget(self.cell_pad_label_z, 7, 0)
        self.grid_layout.addWidget(self.cell_pad_z_edit, 7, 1)
        
        self.grid_layout.addWidget(self.pad_h_label, 8, 0)
        self.grid_layout.addWidget(self.pad_h_edit, 8, 1)
        
        self.grid_layout.addWidget(self.pad_v_label, 9, 0 )
        self.grid_layout.addWidget(self.pad_v_edit, 9, 1)
        
        self.grid_layout.addWidget(self.n_layers_label, 10, 0)
        self.grid_layout.addWidget(self.n_layers_edit, 10, 1)
        
        self.grid_layout.addWidget(self.z1_layer_label, 11, 0)
        self.grid_layout.addWidget(self.z1_layer_edit, 11, 1)
        
        self.grid_layout.addWidget(self.z_target_label, 12, 0)
        self.grid_layout.addWidget(self.z_target_edit, 12, 1)
        
        self.grid_layout.addWidget(self.z_bottom_label, 13, 0)
        self.grid_layout.addWidget(self.z_bottom_edit, 13, 1)
        
        self.grid_layout.addWidget(self.rot_ang_label, 14, 0)
        self.grid_layout.addWidget(self.rot_ang_edit, 14, 1)
        self.grid_layout.addWidget(self.rot_ang_hint, 15, 1)
        
        self.grid_layout.addWidget(self.rho_start_label, 16, 0)
        self.grid_layout.addWidget(self.rho_start_edit, 16, 1)
#        self.grid_layout.addWidget(self.rot_ang_hint, 11, 2)
        
#        self.grid_layout.addWidget(self.edi_button, 1, 0)
#        self.grid_layout.addWidget(self.plot_mesh_button, 1, 1, 1, 2)
#        self.grid_layout.addWidget(self.save_mesh_button, 1, 3)
#        
#        self.grid_layout.addWidget(self.list_widget, 2, 0)
#        self.grid_layout.addWidget(self.output_box, 2, 1, 1, 3)
#        
#        self.grid_layout.addWidget(self.parameters_label, 3, 0)
#        
#        self.grid_layout.addWidget(self.cell_size_label, 4, 0)
#        self.grid_layout.addWidget(self.cell_size_edit_east, 4, 1)
#        self.grid_layout.addWidget(self.cell_size_edit_north, 4, 2)
#        
#        self.grid_layout.addWidget(self.cell_pad_label, 5, 0)
#        self.grid_layout.addWidget(self.cell_pad_east_edit, 5, 1)
#        self.grid_layout.addWidget(self.cell_pad_north_edit, 5, 2)
#        self.grid_layout.addWidget(self.cell_pad_z_edit, 5, 3)
#        
#        self.grid_layout.addWidget(self.pad_h_label, 6, 0)
#        self.grid_layout.addWidget(self.pad_h_edit, 6, 1)
#        self.grid_layout.addWidget(self.pad_v_edit, 6, 2)
#        
#        self.grid_layout.addWidget(self.n_layers_label, 7, 0)
#        self.grid_layout.addWidget(self.n_layers_edit, 7, 1)
#        
#        self.grid_layout.addWidget(self.z1_layer_label, 8, 0)
#        self.grid_layout.addWidget(self.z1_layer_edit, 8, 1)
#        
#        self.grid_layout.addWidget(self.z_target_label, 9, 0)
#        self.grid_layout.addWidget(self.z_target_edit, 9, 1)
#        
#        self.grid_layout.addWidget(self.z_bottom_label, 10, 0)
#        self.grid_layout.addWidget(self.z_bottom_edit, 10, 1)
#        
#        self.grid_layout.addWidget(self.rot_ang_label, 11, 0)
#        self.grid_layout.addWidget(self.rot_ang_edit, 11, 1)
#        self.grid_layout.addWidget(self.rot_ang_hint, 11, 2)

        self.h_layout = QtGui.QHBoxLayout()
        self.h_layout.addLayout(self.edi_grid)
        self.h_layout.addLayout(self.grid_layout)
        
        self.v_layout = QtGui.QVBoxLayout()
        self.v_layout.addLayout(self.h_layout)
        self.v_layout.addWidget(self.output_box)
        
        self.full_layout = QtGui.QHBoxLayout()
        self.full_layout.addLayout(self.v_layout)
        self.full_layout.addWidget(self.mpl_widget)
        
        self.setLayout(self.full_layout)
        
    def get_edi_files(self):
        """
        get a list of edi files
        """
        self.dialog_box = QtGui.QFileDialog()
        fn_list = self.dialog_box.getOpenFileNames(
                                    caption='Choose EDI Files',
                                    filter='*.edi')
        edi_list = []  
        self.list_widget.clear()                                  
        for fn in fn_list:
            self.list_widget.addItem(os.path.basename(str(fn))[:-4])
            edi_list.append(str(fn))
            
        self.model_obj.edi_list = edi_list
         
    def set_cell_size_east(self):
        self.model_obj.cell_size_east = float(str(self.cell_size_edit_east.text()))
        self.cell_size_edit_east.setText('{0:.2f}'.format(self.model_obj.cell_size_east))
            
    def set_cell_size_north(self):
        self.model_obj.cell_size_north = float(str(self.cell_size_edit_north.text()))
        self.cell_size_edit_north.setText('{0:.2f}'.format(self.model_obj.cell_size_north))

        
    def set_cell_pad_east(self):
        self.model_obj.pad_east = int(str(self.cell_pad_east_edit.text()))
        self.cell_pad_east_edit.setText('{0:.0f}'.format(self.model_obj.pad_east))
        
    def set_cell_pad_north(self):
        self.model_obj.pad_north = int(str(self.cell_pad_north_edit.text()))
        self.cell_pad_north_edit.setText('{0:.0f}'.format(self.model_obj.pad_north))
        
    def set_cell_pad_z(self):
        self.model_obj.pad_z = int(str(self.cell_pad_z_edit.text()))
        self.cell_pad_z_edit.setText('{0:.0f}'.format(self.model_obj.pad_z))
        
    def set_pad_h(self):
        self.model_obj.pad_stretch_h = float(str(self.pad_h_edit.text()))
        self.pad_h_edit.setText('{0:.2f}'.format(self.model_obj.pad_stretch_h))
        
    def set_pad_v(self):
        self.model_obj.pad_stretch_v = float(str(self.pad_v_edit.text()))
        self.pad_v_edit.setText('{0:.2f}'.format(self.model_obj.pad_stretch_v))
        
    def set_n_layers(self):
        self.model_obj.n_layers = int(str(self.n_layers_edit.text()))
        self.n_layers_edit.setText('{0:.0f}'.format(self.model_obj.n_layers))
     
    def set_z1_layer(self):
         self.model_obj.z1_layer = float(str(self.z1_layer_edit.text()))
         self.z1_layer_edit.setText('{0:.2f}'.format(self.model_obj.z1_layer))
         
    def set_z_target(self):
        self.model_obj.z_target_depth = float(str(self.z_target_edit.text()))
        self.z_target_edit.setText('{0:.2f}'.format(self.model_obj.z_target_depth))
        
    def set_z_bottom(self):
        self.model_obj.z_bottom = float(str(self.z_bottom_edit.text()))
        self.z_bottom_edit.setText('{0:.2f}'.format(self.model_obj.z_bottom))
        
    def set_rotation_angle(self):
        self.model_obj.mesh_rotation_angle = float(str(self.rot_ang_edit.text()))
        self.rot_ang_edit.setText('{0:.2f}'.format(self.model_obj.mesh_rotation_angle))
        
    def select_station(self):
        pass
    
    def plot_mesh(self):
        self.mpl_widget.plot_mesh(self.model_obj)

    def save_mesh(self):
        fn = str(QtGui.QFileDialog.getSaveFileName(self,
                                    caption='Choose Model File',
                                    directory=os.getcwd()))
                                    
        self.model_obj.write_model_file(save_path=os.path.dirname(fn),
                                        model_fn_basename=os.path.basename(fn))
                                        
    def set_rho(self):
        if self.model_obj.res_model is None:
            self.model_obj.res_model = np.zeros(self.model_obj.grid_north.shape[0],
                                                self.model_obj.grid_east.shape[0],
                                                self.model_obj.grid_z.shape[0])
                                                
        self.model_obj.res_model[:, :, :] = float(str(self.rho_start_edit.text()))
        self.rho_start_edit.setText('{0:.2f}'.format(self.rho_start_edit.text()))
    
    @QtCore.pyqtSlot(str)
    def normal_output(self, message):
        self.output_box.moveCursor(QtGui.QTextCursor.End)
        self.output_box.insertPlainText(message)
        
        

        
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
        self.subplot_left = .12
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .3
        self.subplot_wspace = .3
        
        self.station_marker = 'v'
        self.marker_color = 'b'
        self.marker_size =  2
        
        self.line_color =  'k'
        self.line_width =  .5
        
        self.fs = 10
        
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
    
    def plot_mesh(self, model_obj, east_limits=None, north_limits=None, 
                  z_limits=None):
                      
        try:
            model_obj.make_mesh()
        except AttributeError:
            QtGui.QMessageBox.warning(self,
                                      'Cannot Make Mesh -- Need EDI Files',
                                      "Please press the 'Get EDI Files' button to get files", 
                                      QtGui.QMessageBox.Cancel,
                                      )
            return
            
        
        self.figure.clf()
        
        #make a rotation matrix to rotate data
        #cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        #sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))
        
        #turns out ModEM has not accomodated rotation of the grid, so for
        #now we will not rotate anything.
        cos_ang = 1
        sin_ang = 0
        
        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1])
        
        #--->plot map view    
        ax1 = self.figure.add_subplot(gs[0], aspect='equal')
        
        
        #plot station locations
        plot_east = model_obj.station_locations['rel_east']/1000.
        plot_north = model_obj.station_locations['rel_north']/1000.
        
        plot_grid_east = model_obj.grid_east.copy()/1000.
        plot_grid_north = model_obj.grid_north.copy()/1000.
        plot_grid_z = model_obj.grid_z.copy()/1000.
        
        ax1.scatter(plot_east,
                    plot_north, 
                    marker=self.station_marker,
                    c=self.marker_color,
                    s=self.marker_size)
                                
        
        east_line_xlist = []
        east_line_ylist = []   
        north_min = plot_grid_north.min()         
        north_max = plot_grid_north.max()         
        for xx in plot_grid_east:
            east_line_xlist.extend([xx*cos_ang+north_min*sin_ang, 
                                    xx*cos_ang+north_max*sin_ang])
            east_line_xlist.append(None)
            east_line_ylist.extend([-xx*sin_ang+north_min*cos_ang, 
                                    -xx*sin_ang+north_max*cos_ang])
            east_line_ylist.append(None)
        ax1.plot(east_line_xlist,
                      east_line_ylist,
                      lw=self.line_width,
                      color=self.line_color)

        north_line_xlist = []
        north_line_ylist = [] 
        east_max = plot_grid_east.max()
        east_min = plot_grid_east.min()
        for yy in plot_grid_north:
            north_line_xlist.extend([east_min*cos_ang+yy*sin_ang,
                                     east_max*cos_ang+yy*sin_ang])
            north_line_xlist.append(None)
            north_line_ylist.extend([-east_min*sin_ang+yy*cos_ang, 
                                     -east_max*sin_ang+yy*cos_ang])
            north_line_ylist.append(None)
        ax1.plot(north_line_xlist,
                      north_line_ylist,
                      lw=self.line_width,
                      color=self.line_color)
        
        if east_limits == None:
            ax1.set_xlim(plot_grid_east.min(),
                         plot_grid_east.max())
            pass
        else:
            ax1.set_xlim(east_limits)
        
        if north_limits == None:
            ax1.set_ylim(plot_grid_north.min(),
                         plot_grid_north.max())
            pass
        else:
            ax1.set_ylim(north_limits)
            
        ax1.set_ylabel('Northing (km)', fontdict={'size':12, 'weight':'bold'})
        ax1.set_xlabel('Easting (km)', fontdict={'size':12, 'weight':'bold'})
        
        ##----plot depth view
        ax2 = self.figure.add_subplot(gs[1], aspect='auto')
        

        #plot the grid 
        east_line_xlist = []
        east_line_ylist = []            
        for xx in plot_grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend([0, 
                                    plot_grid_z.max()])
            east_line_ylist.append(None)
        ax2.plot(east_line_xlist,
                 east_line_ylist,
                 lw=self.line_width,
                 color=self.line_color)

        z_line_xlist = []
        z_line_ylist = [] 
        for zz in plot_grid_z:
            z_line_xlist.extend([plot_grid_east.min(),
                                     plot_grid_east.max()])
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        ax2.plot(z_line_xlist,
                 z_line_ylist,
                 lw=self.line_width,
                 color=self.line_color)
                      
        
        #--> plot stations
        ax2.scatter(plot_east,
                    [0]*model_obj.station_locations.shape[0],
                    marker=self.station_marker,
                    c=self.marker_color,
                    s=self.marker_size)

        
        if z_limits == None:
            ax2.set_ylim(model_obj.z_target_depth/1000., -1)
        else:
            ax2.set_ylim(z_limits)
            
        ax2.set_xlim(-model_obj.cell_size_east/6000., 
                     model_obj.cell_size_east/6000.)
        plt.setp(ax2.xaxis.get_ticklabels(), visible=False)
            
        ax2.set_ylabel('Depth (km)', fontdict={'size':9, 'weight':'bold'})
#        ax2.set_xlabel('Easting (m)', fontdict={'size':9, 'weight':'bold'})  
        
        self.mpl_widget.draw()
        
#==============================================================================
#  DEFINE MAIN   
#==============================================================================
def main():
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = ModEM_Mesh_Window()
    ui.ui_setup()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()      
    