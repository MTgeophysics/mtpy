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
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

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
        
        self.ui_setup()
        
    def ui_setup(self):
        # window basics
        self.setWindowTitle("Manipulate ModEM Model")
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        self.model_plot_widget = ModelWidget()
        self.central_widget = self.setCentralWidget(self.model_plot_widget)
        
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
                                       
        self.mesh_widget.modem_data = modem.Data()
        self.mesh_widget.modem_data.read_data_file(fn)
        self.mesh_widget.modem_data_fn = fn
        
        self.mesh_widget.dir_path = os.path.dirname(fn)
        
        
    def get_model_fn(self):
        """ 
        read in an existing model file
        """

        fn_dialog = QtGui.QFileDialog() 
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                       filter='*.rho'))
                                       
        self.mesh_widget.model_obj = modem.Model()
        self.mesh_widget.model_obj.read_model_file(fn)

        self.mesh_widget.dir_path = os.path.dirname(fn)
        
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
        self.mesh_widget.model_obj.write_model_file(save_path=sv_path,
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
        
        self.model_obj = modem.Model()
        self.data_obj = modem.Data()
        self.cov_obj = modem.Covariance()

        self.plot_map = ModelPlotWidget()
        self.plot_north = ModelPlotWidget()
        self.plot_east = ModelPlotWidget()
        self.plot_3d = ModelPlotWidget()

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
        self.map_canvas.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)
                                      
        self.map_input_label = QtGui.QLabel('Index Value')
        self.map_input_edit = QtGui.QLineEdit()
        self.map_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        
        map_bottom_layout = QtGui.QGridLayout()
        map_bottom_layout.addWidget(self.map_input_label, 1, 1, 1, 1)
        map_bottom_layout.addWidget(self.map_input_edit, 1, 2, 1, 1)
        map_bottom_layout.addWidget(self.map_slider, 1, 3, 1, 150)
        map_layout = QtGui.QVBoxLayout()
        map_layout.addWidget(self.map_canvas)
        map_layout.addLayout(map_bottom_layout)
        
        
        self.figure_east = Figure()
        self.canvas_east = FigureCanvas(self.figure_east)
        self.canvas_east.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                       QtGui.QSizePolicy.Expanding)
        self.figure_north = Figure()
        self.canvas_north = FigureCanvas(self.figure_north)
        self.canvas_north.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)
        self.figure_3d = Figure()
        self.canvas_3d = FigureCanvas(self.figure_3d)
        self.canvas_3d.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                       QtGui.QSizePolicy.Expanding)
                                       
        
#        self.plot_map = ModelPlotWidget()
#        self.plot_east = ModelPlotWidget()
#        self.plot_north = ModelPlotWidget()
#        self.plot_3d = ModelPlotWidget()
        
        self.grid_layout = QtGui.QGridLayout()
        self.grid_layout.addLayout(map_layout, 1, 1)
        self.grid_layout.addWidget(self.canvas_east, 1, 2)
        self.grid_layout.addWidget(self.canvas_north, 2, 1)
        self.grid_layout.addWidget(self.canvas_3d, 2, 2)
        
        self.setLayout(self.grid_layout)
    
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
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = ModEM_Model_Manipulator()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()   
        

