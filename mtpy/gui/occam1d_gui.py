# -*- coding: utf-8 -*-
"""
Occam 1D 
-----------------

All encompassing ploting data, model and model response.


JP 2015
"""
# 

from PyQt4 import QtCore, QtGui
import mtpy.modeling.occam1d as occam1d
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
import mtpy.imaging.mtplottools as mtplottools
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt
import os
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

class Occam1D_GUI(QtGui.QMainWindow):
    def __init__(self):
        super(Occam1D_GUI, self).__init__()
        self.ms  = 5
        self.lw = 1.5
        self.data_marker_color = (0, 0, 0)
        self.data_marker = 's'
        self.model_marker_color = (1, 0, 0)
        self.model_marker = 'x'
        self.e_capthick = 1
        self.e_capsize =  5
         
        self.res_limits = None 
        self.phase_limits = None 
 
        self.subplot_wspace = .25
        self.subplot_hspace = .0
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .93
        self.subplot_bottom = .08
        
        self.legend_loc = 'upper center'
        self.legend_pos = (.5, 1.15)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.fs = 11
        self.ylabel_pad = 1.25
        self.station = None
        
        self.occam_data = None
        self.occam_resp = None
        self.dir_path = os.getcwd()
        
        self.ui_setup()
        
    def ui_setup(self):
        """
        setup the window layout
        """
        
        self.setWindowTitle("Run Occam 1D")
        self.resize(1920, 1080)
        
        self.occam_widget = OccamWidget()
        self.central_widget = self.setCentralWidget(self.occam_widget)

        #--------------------------------------------------------
        # stream the output of occam 1D        
        self.my_stream = MyStream()
        self.my_stream.message.connect(self.occam_widget.normal_output)
        
        sys.stdout = self.my_stream
        
        QtCore.QMetaObject.connectSlotsByName(self)
            
#==============================================================================
# Occam 1D widget
#==============================================================================
class OccamWidget(QtGui.QWidget):
    """
    occam 1D widget
    """
        
    def __init__(self):
        super(OccamWidget, self).__init__()
        
        self.occam_data = occam1d.Data()
        self.occam_model = occam1d.Model()
        self.occam_startup = occam1d.Startup()
        self.mpl_widget = OccamPlot()        
        
        self.res_err = 10
        self.phase_err = 5
        self.data_mode = 'TE'

        self.setup_ui()
        
    def setup_ui(self):
        """
        setup the user interface
        """
        # font type to use for labels        
        label_font = QtGui.QFont()
        label_font.setBold = True
        label_font.setPointSize (16)
        
        #---------------------------------------------------
        # menu bar
#        self.menu_bar = QtGui.QMenuBar(self)
#        self.menu_bar.setGeometry(QtCore.QRect(0, 0, 1920, 40))
#        
#        self.menu_data = QtGui.QMenu(self.menu_bar)
#        self.menu_data.setTitle("Data")
#        
#        self.menu_model = QtGui.QMenu(self.menu_bar)
#        self.menu_model.setTitle("Model")
#        
#        self.menu_startup = QtGui.QMenu(self.menu_bar)
#        self.menu_startup.setTitle('Startup')
#        
#        self.setMenuBar(self.menu_bar)
        
        self.get_edi_button = QtGui.QPushButton('Get EDI File')
        self.get_edi_button.clicked.connect(self.get_edi_file)
        
        self.data_label = QtGui.QLabel('Data Parameters')
        self.data_label.setFont(label_font)         
        
        self.data_res_err_label = QtGui.QLabel('Res. Error (%)')
        self.data_res_err_edit = QtGui.QLineEdit()
        self.data_res_err_edit.setText('{0:.2f}'.format(self.res_err))
        self.data_res_err_edit.editingFinished.connect(self.set_res_err)
        
        self.data_phase_err_label = QtGui.QLabel('Phase Error (%)')
        self.data_phase_err_edit = QtGui.QLineEdit()
        self.data_phase_err_edit.setText('{0:.2f}'.format(self.phase_err))
        self.data_phase_err_edit.editingFinished.connect(self.set_phase_err)
        
        self.data_mode_label = QtGui.QLabel('Mode')
        self.data_mode_combo = QtGui.QComboBox()
        self.data_mode_combo.addItem('TE')
        self.data_mode_combo.addItem('TM')
        self.data_mode_combo.addItem('Det')
        self.data_mode_combo.activated[str].connect(self.set_data_mode)
        
        # vertical layer parameters
        self.model_label = QtGui.QLabel('Model Parameters')
        self.model_label.setFont(label_font)
        
        self.n_layers_label = QtGui.QLabel('Number of Vertical Layers')
        self.n_layers_edit = QtGui.QLineEdit()
        self.n_layers_edit.setText('{0:.0f}'.format(self.occam_model.n_layers))
        self.n_layers_edit.editingFinished.connect(self.set_n_layers)
        
        self.z1_layer_label = QtGui.QLabel('Thicknes of 1st layer (m)')
        self.z1_layer_edit = QtGui.QLineEdit()
        self.z1_layer_edit.setText('{0:.2f}'.format(self.occam_model.z1_layer))
        self.z1_layer_edit.editingFinished.connect(self.set_z1_layer)
        
        self.z_target_label = QtGui.QLabel('Target Depth (m)')
        self.z_target_edit = QtGui.QLineEdit()
        self.z_target_edit.setText('{0:.2f}'.format(self.occam_model.target_depth))
        self.z_target_edit.editingFinished.connect(self.set_z_target)
        
        self.z_bottom_label = QtGui.QLabel('Bottom of the Model (m)')
        self.z_bottom_edit = QtGui.QLineEdit()
        self.z_bottom_edit.setText('{0:.2f}'.format(self.occam_model.bottom_layer))
        self.z_bottom_edit.editingFinished.connect(self.set_z_bottom)
        
        # starting resistivity
        self.rho_start_label = QtGui.QLabel('Starting rho (Ohmm)')
        self.rho_start_edit = QtGui.QLineEdit()
        self.rho_start_edit.setText('{0:.2f}'.format(100))
        self.rho_start_edit.editingFinished.connect(self.set_rho)

        #---set the layout---------------
        data_grid = QtGui.QGridLayout()
        data_grid.addWidget(self.data_label, 0, 0)
        
        data_grid.addWidget(self.data_res_err_label, 1, 0)
        data_grid.addWidget(self.data_res_err_edit, 1, 1)
        
        data_grid.addWidget(self.data_phase_err_label, 2, 0)
        data_grid.addWidget(self.data_phase_err_edit, 2, 1)
        
        data_grid.addWidget(self.data_mode_label, 3, 0)
        data_grid.addWidget(self.data_mode_combo, 3, 1)
        
        model_grid = QtGui.QGridLayout()
        model_grid.addWidget(self.model_label, 0, 0)
        
        model_grid.addWidget(self.n_layers_label, 1, 0)
        model_grid.addWidget(self.n_layers_edit, 1, 1)
        
        model_grid.addWidget(self.z1_layer_label, 2, 0)
        model_grid.addWidget(self.z1_layer_edit, 2, 1)
        
        model_grid.addWidget(self.z_target_label, 3, 0)
        model_grid.addWidget(self.z_target_edit, 3, 1)
        
        model_grid.addWidget(self.z_bottom_label, 4, 0)
        model_grid.addWidget(self.z_bottom_edit, 4, 1)
        
        
        model_grid.addWidget(self.rho_start_label, 5, 0)
        model_grid.addWidget(self.rho_start_edit, 5, 1)
        
        edit_layout = QtGui.QVBoxLayout()
        edit_layout.addWidget(self.get_edi_button)
        edit_layout.addLayout(data_grid)
        edit_layout.addLayout(model_grid)
        
        self.output_box = QtGui.QTextEdit()
        
        plot_layout = QtGui.QVBoxLayout()
        plot_layout.addWidget(self.mpl_widget)
        plot_layout.addWidget(self.output_box)
        
        window_layout = QtGui.QHBoxLayout()
        window_layout.addLayout(edit_layout)
        window_layout.addLayout(plot_layout)
        
        self.setLayout(window_layout)
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_edi_file(self):
        """
        get edi file to invert
        """
        
        edi_dialog = QtGui.QFileDialog()
        fn = str(edi_dialog.getOpenFileName(caption='Pick .edi file',
                                            filter='*.edi'))
        self.edi_fn = fn
        
        self.mpl_widget.plot_mesh()
            
    def set_res_err(self):
        self.res_err = float(str(self.data_res_err_edit.text()))
        
    def set_phase_err(self):
        self.phase_err = float(str(self.data_phase_err_edit.text()))
        
    def set_data_mode(self, text):
        self.data_mode = text
        
    def set_n_layers(self):
        self.occam_model.n_layers = int(str(self.n_layers_edit.text()))
        self.n_layers_edit.setText('{0:.0f}'.format(self.occam_model.n_layers))
     
    def set_z1_layer(self):
         self.occam_model.z1_layer = float(str(self.z1_layer_edit.text()))
         self.z1_layer_edit.setText('{0:.2f}'.format(self.occam_model.z1_layer))
         
    def set_z_target(self):
        self.occam_model.target_depth = float(str(self.z_target_edit.text()))
        self.z_target_edit.setText('{0:.2f}'.format(self.occam_model.target_depth))
        
    def set_z_bottom(self):
        self.occam_model.bottom_layer = float(str(self.z_bottom_edit.text()))
        self.z_bottom_edit.setText('{0:.2f}'.format(self.occam_model.bottom_layer))
        
    def set_rho(self):
        self.occam_startup.start_rho = float(str(self.rho_start_edit.text()))
        self.rho_start_edit.setText('{0:.2f}'.format(self.occam_startup.start_rho))
        
    @QtCore.pyqtSlot(str)
    def normal_output(self, message):
        self.output_box.moveCursor(QtGui.QTextCursor.End)
        self.output_box.insertPlainText(message)
        
#==============================================================================
# Mesh Plot
#==============================================================================
class OccamPlot(QtGui.QWidget):
    """
    plotting the mesh
    """    
    
    def __init__(self):
        super(OccamPlot, self).__init__()
        
        self.subplot_wspace = .25
        self.subplot_hspace = .15
        self.subplot_right = .92
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1
        
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
    
    def plot_mesh(self):
        """
        plot response and depth model
        """
        self.figure.clf()
        
        #make a grid of subplots
        gs=gridspec.GridSpec(6, 5, hspace=self.subplot_hspace, 
                             wspace=self.subplot_wspace)
        #subplot resistivity
        self.axr = self.figure.add_subplot(gs[:4, :4])
        
        #subplot for phase
        self.axp = self.figure.add_subplot(gs[4:,:4], sharex=self.axr)

        #subplot for model
        self.axm = self.figure.add_subplot(gs[:, 4])
        
        self.mpl_widget.draw()

#==============================================================================
# Main execution        
#==============================================================================
def main():
#if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    ui = Occam1D_GUI()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()