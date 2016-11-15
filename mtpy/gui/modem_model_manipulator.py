# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 12:20:42 2016

@author: jpeacock
"""


#==============================================================================
# Imports
#==============================================================================
import os
import sys

from PyQt4 import QtCore, QtGui

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import numpy as np

from mtpy.gui.my_stream import MyStream
import mtpy.modeling.modem_new as modem

#==============================================================================
# Main Window
#==============================================================================
class ModEM_Model_Manipulator(QtGui.QMainWindow):
    """
    main window for manipulating a model
    """
    
    def __init__(self):
        super(ModEM_Model_Manipulator, self).__init__()
        
        self.model_widget = ModelWidget()
        
        self.ui_setup()
        
    def ui_setup(self):
        # window basics
        self.setWindowTitle("Manipulate ModEM Model")
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        self.central_widget = self.setCentralWidget(self.model_widget)
        
        #-------------- MENU BAR ---------------------------------
        # add a menu bar to the top of the window        
        self.menu_data_file = self.menuBar().addMenu("Data &File")
        self.menu_data_open_action = self.menu_data_file.addAction("Open")
        self.menu_data_open_action.triggered.connect(self.get_data_fn)
        
        self.menu_model_file = self.menuBar().addMenu("&Model File")
        self.menu_model_open_action = self.menu_model_file.addAction("Open")
        self.menu_model_open_action.triggered.connect(self.get_model_fn)
        
        self.menu_model_save_action = self.menu_model_file.addAction("Save")
        self.menu_model_save_action.triggered.connect(self.save_model_fn)
        
        #------------Output stream box-------------------------------------
#        self.my_stream = MyStream()
#        self.my_stream.message.connect(self.mesh_widget.normal_output)
#        
#        sys.stdout = self.my_stream
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_data_fn(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)'))
                                       
        self.model_widget.data_fn = fn        
        
    def get_model_fn(self):
        """ 
        read in an existing model file
        """

        fn_dialog = QtGui.QFileDialog() 
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                       filter='*.rho'))
                                       
        self.model_widget.model_fn = fn
        
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
        self.model_widget.model_obj.write_model_file(save_path=sv_path,
                                                    model_fn_basename=sv_basename)
                                                    
#==============================================================================
# Model Widget
#==============================================================================
class ModelWidget(QtGui.QWidget):
    """
    make the model plot its own widget
    """
    
    def __init__(self):
        super(ModelWidget, self).__init__()
        
        self.model_obj = None
        self.data_obj = None
        self.cov_obj = None
        
        self._data_fn = None
        self._model_fn = None
        
        self.map_index = 0
        self.east_index = 0
        self.north_index = 0 
        
        self.plot_east_map = None
        self.plot_north_map = None
        self.plot_east_z = None
        self.plot_z_east = None
        self.plot_north_z = None
        self.plot_z_north = None
        
        self.units = 'km'
        self.scale = 1000.
        
        self.cmap = 'jet_r'
        self.res_limits = (0, 4)

        self.ui_setup()
 
    def ui_setup(self):
#        
#        self.colorbar_widget = mcb.ColorbarBase(self.ax2,cmap=self.cmap,
#                                   norm=colors.Normalize(vmin=self.cmin,
#                                                         vmax=self.cmax),
#                                    orientation='horizontal')
#                                                         
#                            
#        self.cb.set_label('Resistivity ($\Omega \cdot$m)',
#                          fontdict={'size':self.font_size})
#        self.cb.set_ticks(np.arange(self.cmin, self.cmax+1))
#        self.cb.set_ticklabels([mtplottools.labeldict[cc] 
#                                for cc in np.arange(self.cmin, self.cmax+1)])
    
        self.map_figure = Figure()
        self.map_canvas = FigureCanvas(self.map_figure)
        self.map_canvas.mpl_connect('pick_event', self.map_on_pick)
        self.map_canvas.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)
        self.map_toolbar = NavigationToolbar(self.map_canvas, self)
                                      
        self.map_depth_label = QtGui.QLabel('Depth {0:>10.2f} {1}'.format(0, 
                                                                    self.units))
        self.map_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.map_slider.valueChanged.connect(self.set_map_index)
        self.map_slider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.map_slider.setMinimum(0)
        self.map_slider.setMaximum(0)
        self.map_slider.setTickInterval(1)
        

        self.east_figure= Figure()
        self.east_canvas = FigureCanvas(self.east_figure)
        self.east_canvas.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                       QtGui.QSizePolicy.Expanding)
        self.east_toolbar = NavigationToolbar(self.east_canvas, self)
                                       
        self.east_label = QtGui.QLabel('Easting {0:>10.2f} {1}'.format(0, 
                                                                   self.units))
        self.east_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.east_slider.valueChanged.connect(self.set_east_index)
        self.east_slider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.east_slider.setMinimum(0)
        self.east_slider.setMaximum(0)
        self.east_slider.setTickInterval(1)
        
        self.north_figure = Figure()
        self.north_canvas = FigureCanvas(self.north_figure)
        self.north_canvas.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)
        self.north_toolbar = NavigationToolbar(self.north_canvas, self)
                                      
        self.north_label = QtGui.QLabel('Northing {0:>10.2f} m'.format(0,
                                                                    self.units))
        self.north_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.north_slider.valueChanged.connect(self.set_north_index)
        self.north_slider.setTickPosition(QtGui.QSlider.TicksBelow)
        self.north_slider.setMinimum(0)
        self.north_slider.setMaximum(0)
        self.north_slider.setTickInterval(1)
        
        
        self.figure_3d = Figure()
        self.canvas_3d = FigureCanvas(self.figure_3d)
        self.canvas_3d.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                       QtGui.QSizePolicy.Expanding)

        ##------------------------------------------------
        ## Layout
        
        map_bottom_layout = QtGui.QHBoxLayout()
        map_bottom_layout.addWidget(self.map_depth_label)
        map_bottom_layout.addWidget(self.map_slider)
        map_layout = QtGui.QVBoxLayout()
        map_layout.addWidget(self.map_toolbar)
        map_layout.addWidget(self.map_canvas)
        map_layout.addLayout(map_bottom_layout)
        
        east_bottom_layout = QtGui.QHBoxLayout()
        east_bottom_layout.addWidget(self.east_label)
        east_bottom_layout.addWidget(self.east_slider)
        east_layout = QtGui.QVBoxLayout()
        east_layout.addWidget(self.east_toolbar)
        east_layout.addWidget(self.east_canvas)
        east_layout.addLayout(east_bottom_layout)
        
        north_bottom_layout = QtGui.QHBoxLayout()
        north_bottom_layout.addWidget(self.north_label)
        north_bottom_layout.addWidget(self.north_slider)
        north_layout = QtGui.QVBoxLayout()
        north_layout.addWidget(self.north_toolbar)
        north_layout.addWidget(self.north_canvas)
        north_layout.addLayout(north_bottom_layout)
        
        
        self.grid_layout = QtGui.QGridLayout()
        self.grid_layout.addLayout(map_layout, 1, 1)
        self.grid_layout.addLayout(east_layout, 1, 2)
        self.grid_layout.addLayout(north_layout, 2, 1)
        self.grid_layout.addWidget(self.canvas_3d, 2, 2)
        
        self.setLayout(self.grid_layout)
      
    @property
    def data_fn(self):
        self._data_fn
    
    @data_fn.getter
    def data_fn(self):
        return self._data_fn
    @data_fn.setter
    def data_fn(self, data_fn):
        self._data_fn = data_fn
        
        self.data_obj = modem.Data()
        self.data_obj.read_data_file(self._data_fn)
        
        
    @property
    def model_fn(self):
        self._model_fn
    
    @model_fn.getter
    def model_fn(self):
        return self._model_fn
        
    @model_fn.setter
    def model_fn(self, model_fn):
        self._model_fn = model_fn
        self.model_obj = modem.Model()
        self.model_obj.read_model_file(self._model_fn)
        
        # set slider bar intervals
        # need the minus 1 cause we are using the value of the slider as
        # the index.
        self.map_slider.setMaximum(self.model_obj.grid_z.size-1)
        self.east_slider.setMaximum(self.model_obj.grid_north.size-1)
        self.north_slider.setMaximum(self.model_obj.grid_east.size-1)
        
        ## plot the model
        plot_east = np.append(self.model_obj.grid_east,
                              self.model_obj.grid_east[-1]*1.2)/self.scale
        plot_north = np.append(self.model_obj.grid_north, 
                               self.model_obj.grid_north[-1]*1.2)/self.scale
        plot_z = np.append(self.model_obj.grid_z, 
                               self.model_obj.grid_z[-1]*1.2)/self.scale
                               
        self.plot_east_map, self.plot_north_map = np.meshgrid(plot_east, 
                                                              plot_north,
                                                              indexing='ij')
        self.plot_east_z, self.plot_z_east = np.meshgrid(plot_east, 
                                                         plot_z,
                                                         indexing='ij')
        self.plot_north_z, self.plot_z_north = np.meshgrid(plot_north, 
                                                           plot_z,
                                                           indexing='ij')
                                                           
                                                           
        self.map_ax = self.map_figure.add_subplot(1, 1, 1, aspect='equal')
        self.map_ax.pcolormesh(self.plot_east_map, 
                               self.plot_north_map,
                               np.log10(self.model_obj.res_model[:, :, self.map_index].T),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.map_ax.set_xlabel('Easting {0}'.format(self.units))
        self.map_ax.set_ylabel('Northing {0}'.format(self.units))
        self.map_ax.axis('tight')
        self.map_canvas.draw()
        
        self.north_ax = self.north_figure.add_subplot(1, 1, 1, aspect='equal')
        self.north_ax.pcolormesh(self.plot_east_z, 
                               self.plot_z_east,
                               np.log10(self.model_obj.res_model[self.north_index, :, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.north_ax.set_xlabel('Easting {0}'.format(self.units))
        self.north_ax.set_ylabel('Depth {0}'.format(self.units))
        z_lim = self.north_ax.get_ylim()
        self.north_ax.set_ylim(z_lim[1], z_lim[0])
        self.north_ax.axis('tight')
        self.north_canvas.draw()
                               
        self.east_ax = self.east_figure.add_subplot(1, 1, 1, aspect='equal')
        self.east_ax.pcolormesh(self.plot_north_z, 
                               self.plot_z_north,
                               np.log10(self.model_obj.res_model[:, self.east_index, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
                               
        self.east_ax.set_xlabel('Northing {0}'.format(self.units))
        self.east_ax.set_ylabel('Depth {0}'.format(self.units))
        z_lim = self.east_ax.get_ylim()
        self.east_ax.set_ylim(z_lim[1], z_lim[0])
        self.east_ax.axis('tight')
        self.east_canvas.draw()
        
        
        
        
    def set_map_index(self):
        self.map_index = int(self.map_slider.value())
        depth = self.model_obj.grid_z[self.map_index]/self.scale
        self.map_depth_label.setText('Depth {0:>10.2f} {1}'.format(depth,
                                                                self.units))
                                                                
        self.map_ax.pcolormesh(self.plot_east_map, 
                               self.plot_north_map,
                               np.log10(self.model_obj.res_model[:, :, self.map_index].T),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.map_canvas.draw()
                                        
    def set_east_index(self):
        self.east_index = int(self.east_slider.value())
        easting = self.model_obj.grid_north[self.east_index]/self.scale

        self.east_label.setText('Easting {0:>10.2f} {1}'.format(easting,
                                                                self.units))
                                                                
        self.east_ax.pcolormesh(self.plot_north_z, 
                               self.plot_z_north,
                               np.log10(self.model_obj.res_model[:, self.east_index, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.east_canvas.draw()
        
    def set_north_index(self):
        self.north_index = int(self.north_slider.value())
        northing = self.model_obj.grid_north[self.north_index]/self.scale
        self.north_label.setText('Northing {0:>10.2f} {1}'.format(northing,
                                                                self.units))
                                                                
        self.north_ax.pcolormesh(self.plot_east_z, 
                               self.plot_z_east,
                               np.log10(self.model_obj.res_model[self.north_index, :, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.north_canvas.draw()
        
    def map_on_pick(self):
        pass
    
#==============================================================================
#  generic plot widget
#==============================================================================
class ModelPlotWidget(QtGui.QWidget):
    """
    generic plot widget that has a plot and slider bar and navigation bar
    """
    
    def __init__(self):
        super(ModelPlotWidget, self).__init__()
        
        self.plot_index = 0
        
        self.res = 2.0
        self.res_model = None
        self.nodes_x = None
        self.nodes_y = None
        
        self.subplot_right = .99
        self.subplot_left = .12
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .3
        self.subplot_wspace = .3
        
        self.station_marker = 'v'
        self.marker_color = 'k'
        self.marker_size = 2
        self.fig_dpi = 150
        
        self.ui_setup()
    
    def ui_setup(self):
        
        self.figure = Figure(dpi=self.fig_dpi)
        self.figure.subplots_adjust(left=self.subplot_left,
                                    right=self.subplot_right,
                                    bottom=self.subplot_bottom,
                                    top=self.subplot_top,
                                    hspace=self.subplot_hspace,
                                    wspace=self.subplot_wspace)
        
        self.mpl_widget = FigureCanvas(self.figure)

        self.mpl_widget.mpl_connect('pick_event', self.on_pick)
        self.mpl_widget.mpl_connect('axes_enter_event', self.in_axes)
                                     
        self.mpl_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                     QtGui.QSizePolicy.Expanding)
        
        self.mpl_navigation_bar = NavigationToolbar(self.mpl_widget, self)
        
        self.slider_bar = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider_bar.setMinimum(0)
        self.slider_bar.setMaximum(10)
        self.slider_bar.setValue(0)
        self.slider_bar.setTickPosition(QtGui.QSlider.TicksBelow)
        self.slider_bar.setTickInterval(1)
        
        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.mpl_navigation_bar)
        self.layout.addWidget(self.mpl_widget)
        self.layout.addWidget(self.slider_bar)
        
        self.setLayout(self.layout)
        
        
    def on_pick(self):
        pass
    
    def in_axes(self):
        pass
    
    def plot(self):
        pass
    

#==============================================================================
#  DEFINE MAIN   
#==============================================================================
def main():
    app = QtGui.QApplication(sys.argv)
    ui = ModEM_Model_Manipulator()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()   
        

