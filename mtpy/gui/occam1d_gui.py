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
from matplotlib.ticker import MultipleLocator
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

        ## menu bar
        self.menu_bar = self.menuBar()
        #self.menu_bar.setGeometry(QtCore.QRect(0, 0, 1920, 40))
        
        self.menu_help = self.menu_bar.addMenu('Help')
        help_action = QtGui.QAction('help doc', self)
        help_action.triggered.connect(self.display_help)
        
        
        self.menu_help.addAction(help_action)
        
        #self.menu_model = QtGui.QMenu(self.menu_bar)
        #self.menu_model.setTitle("Model")
        
        #self.menu_startup = QtGui.QMenu(self.menu_bar)
        #self.menu_startup.setTitle('Startup')
        
        #self.setMenuBar(self.menu_bar)
        #self.show()
        #--------------------------------------------------------
        # stream the output of occam 1D        
        self.my_stream = MyStream()
        self.my_stream.message.connect(self.occam_widget.normal_output)
        
        sys.stdout = self.my_stream
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def display_help(self):
        ll = ['***Be sure you have a working executable of Occam1D first***\n',
              'To begin: ',
              '\t* select an edi file to invert by clicking the ',
              '\t  "Get EDI File" button at the top left.',
              '\t  The TE mode will be plotted meanint the file was read',
              '\t* Change the parameters in the Data, Model, and Startup fields',
              '\t* Locate Occam1D on your system by clicking "Occam1D Path"',
              '\t* Hit the "Run" button to run an inversion.',
              '\t  The first iteration will be plotted once it is finished.',
              '',
              'An L2 curve will be shown in the lower plot to give you an',
              'idea of which iteration is the optimum. To change iterations',
              'pick a number on the combination box labeled "Iteration".',
              '',
              'Change the parameters and try again.',
              '',
              'The output will be shown on the left handside.',
              '',
              'To save an image of the model and response click the disk icon']
            
        help_string = '\n'.join(ll)        
        
        QtGui.QMessageBox.information(self.central_widget, 'Help', help_string)
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
        self.occam_exec = ''
        self.mpl_widget = OccamPlot() 
        
        self.l2_widget = PlotL2()
        self.l2_widget.l2_widget.mpl_connect('pick event', self.on_click)
        self.l2_widget.l2_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.l2_widget.l2_widget.setFocus()        
        
        self.res_err = 10.
        self.phase_err = 5.
        self.data_mode = 'TE'
        self.edi_fn = ''
        
        self.save_dir = None
        self.station_dir = None

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
        
        self.get_occam_path_button = QtGui.QPushButton('Occam1D Path')
        self.get_occam_path_button.clicked.connect(self.get_occam_path)
        
        self.get_occam_path_edit = QtGui.QLineEdit()
        self.get_occam_path_edit.setText(self.occam_exec)
        self.get_occam_path_edit.editingFinished.connect(self.get_occam_path)
        
        self.get_edi_button = QtGui.QPushButton('Get EDI File')
        self.get_edi_button.clicked.connect(self.get_edi_file)
        
        self.get_edi_edit = QtGui.QLineEdit()
        self.get_edi_edit.setText(self.edi_fn)
        self.get_edi_edit.editingFinished.connect(self.get_edi_file)
        
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
        #self.data_mode_combo.addItem('Both')
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
        self.startup_label = QtGui.QLabel('Startup Parameters')
        self.startup_label.setFont(label_font)
        
        self.start_rho_label = QtGui.QLabel('Starting rho (Ohmm)')
        self.start_rho_edit = QtGui.QLineEdit()
        self.start_rho_edit.setText('{0:.2f}'.format(self.occam_startup.start_rho))
        self.start_rho_edit.editingFinished.connect(self.set_rho)
        
        self.max_iter_label = QtGui.QLabel('Num of Iterations')
        self.max_iter_edit = QtGui.QLineEdit()
        self.max_iter_edit.setText('{0:.0f}'.format(self.occam_startup.max_iter))
        self.max_iter_edit.editingFinished.connect(self.set_max_iter)
        
        self.target_rms_label = QtGui.QLabel('Target RMS')
        self.target_rms_edit = QtGui.QLineEdit()
        self.target_rms_edit.setText('{0:.2f}'.format(self.occam_startup.target_rms))
        self.target_rms_edit.editingFinished.connect(self.set_target_rms)
        
        self.start_roughness_label = QtGui.QLabel('Starting Roughness')
        self.start_roughness_edit = QtGui.QLineEdit()
        self.start_roughness_edit.setText('{0:.2f}'.format(self.occam_startup.start_rough))
        self.start_roughness_edit.editingFinished.connect(self.set_start_rough)

        self.start_lagrange_label = QtGui.QLabel('Starting Lagrange')
        self.start_lagrange_edit = QtGui.QLineEdit()
        self.start_lagrange_edit.setText('{0:.2f}'.format(self.occam_startup.start_lagrange))
        self.start_lagrange_edit.editingFinished.connect(self.set_start_lagrange)
        
        self.iter_combo_label = QtGui.QLabel('Iteration')
        self.iter_combo_edit = QtGui.QComboBox()
        self.iter_combo_edit.addItem('1')
        self.iter_combo_edit.activated[str].connect(self.set_iteration)
        self.iter_combo_edit.setSizeAdjustPolicy(QtGui.QComboBox.AdjustToContents)
        self.iter_combo_edit.setMinimumWidth(50)
        
        self.output_box = QtGui.QTextEdit()
        #---set the layout---------------
        path_layout = QtGui.QHBoxLayout()
        path_layout.addWidget(self.get_occam_path_button)
        path_layout.addWidget(self.get_occam_path_edit)        
        
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

        
        startup_grid = QtGui.QGridLayout()
        startup_grid.addWidget(self.startup_label, 0, 0)
        
        startup_grid.addWidget(self.target_rms_label, 1, 0)
        startup_grid.addWidget(self.target_rms_edit, 1, 1)
        
        startup_grid.addWidget(self.max_iter_label, 2, 0)
        startup_grid.addWidget(self.max_iter_edit, 2, 1)
        
        startup_grid.addWidget(self.start_rho_label, 3, 0)
        startup_grid.addWidget(self.start_rho_edit, 3, 1)
        
        startup_grid.addWidget(self.start_lagrange_label, 4, 0)
        startup_grid.addWidget(self.start_lagrange_edit, 4, 1)
        
        startup_grid.addWidget(self.start_roughness_label, 5, 0)
        startup_grid.addWidget(self.start_roughness_edit, 5, 1)
        
        run_button = QtGui.QPushButton()
        run_button.setText('Run')
        run_button.clicked.connect(self.run_occam)
        
        edi_layout = QtGui.QHBoxLayout()
        edi_layout.addWidget(self.get_edi_button)
        edi_layout.addWidget(self.get_edi_edit)
        
        edit_layout = QtGui.QVBoxLayout()
        edit_layout.addLayout(edi_layout)
        edit_layout.addLayout(data_grid)
        edit_layout.addLayout(model_grid)
        edit_layout.addLayout(startup_grid)
        edit_layout.addWidget(self.output_box)
        edit_layout.addLayout(path_layout)
        edit_layout.addWidget(run_button)
        
        bottom_plot_layout = QtGui.QHBoxLayout()
        bottom_plot_layout.addWidget(self.iter_combo_label)
        bottom_plot_layout.addWidget(self.iter_combo_edit)
        bottom_plot_layout.addWidget(self.l2_widget)
        
        plot_layout = QtGui.QGridLayout()
        plot_layout.addWidget(self.mpl_widget, 0, 0, 1, 1)
        plot_layout.addLayout(bottom_plot_layout, 2, 0, 2, 1)
        
#        window_layout = QtGui.QHBoxLayout()
#        window_layout.addLayout(edit_layout)
#        window_layout.addLayout(plot_layout)
        
        window_grid = QtGui.QGridLayout()
        window_grid.addLayout(edit_layout, 0, 0, 1, 5)
        window_grid.addLayout(plot_layout, 0, 5, 1, 1)
        
        self.setLayout(window_grid)
        
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def get_occam_path(self):
        """
        get occam path
        """        
        
        occam_path_dialog = QtGui.QFileDialog()
        fn = str(occam_path_dialog.getOpenFileName(
                                        caption='Locate Occam1D executable'))
                                        
        self.occam_exec = os.path.abspath(fn) 
        self.get_occam_path_edit.setText(self.occam_exec)
        
    def get_edi_file(self):
        """
        get edi file to invert
        """
        if self.edi_fn is not '':
            edi_path = os.path.dirname(self.edi_fn)
            edi_dialog = QtGui.QFileDialog()
            fn = str(edi_dialog.getOpenFileName(caption='Pick .edi file',
                                                filter='*.edi',
                                                directory=edi_path))
        else:
            edi_dialog = QtGui.QFileDialog()
            fn = str(edi_dialog.getOpenFileName(caption='Pick .edi file',
                                                filter='*.edi'))
        self.edi_fn = fn
        self.get_edi_edit.setText(self.edi_fn)
        
        station = os.path.basename(self.edi_fn)[:-4]
        self.station_dir = os.path.join(os.path.dirname(self.edi_fn), 
                                        station)
        if not os.path.isdir(self.station_dir):
            os.mkdir(self.station_dir)
            print 'Made director {0}'.format(self.station_dir)
        
        self.save_dir = os.path.join(self.station_dir)
        
        # make an initial data file
        self.occam_data.write_data_file(edi_file=self.edi_fn,
                                        save_path=self.save_dir,
                                        mode=self.data_mode,
                                        res_err=self.res_err,
                                        phase_err=self.phase_err,
                                        thetar=0)
        
        self.mpl_widget.plot_data(data_fn=self.occam_data.data_fn)
            
    def set_res_err(self):
        self.res_err = float(str(self.data_res_err_edit.text()))
        self.data_res_err_edit.setText('{0:.2f}'.format(self.res_err))
        
    def set_phase_err(self):
        self.phase_err = float(str(self.data_phase_err_edit.text()))
        self.data_phase_err_edit.setText('{0:.2f}'.format(self.phase_err))
        
    def set_data_mode(self, text):
        self.data_mode = str(text)
        
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
        self.occam_startup.start_rho = float(str(self.start_rho_edit.text()))
        self.start_rho_edit.setText('{0:.2f}'.format(self.occam_startup.start_rho))
    
    def set_max_iter(self):
        self.occam_startup.max_iter = int(str(self.max_iter_edit.text()))
        self.max_iter_edit.setText('{0:.0f}'.format(self.occam_startup.max_iter))

    def set_target_rms(self):
        self.occam_startup.target_rms = float(str(self.target_rms_edit.text()))
        self.target_rms_edit.setText('{0:.2f}'.format(self.occam_startup.target_rms))

    def set_start_rough(self):
        self.occam_startup.start_rough = float(str(self.start_roughness_edit.text()))
        self.start_rough_edit.setText('{0:.2f}'.format(self.occam_startup.start_rough))

    def set_start_lagrange(self):
        self.occam_startup.start_lagrange = float(str(self.start_lagrange_edit.text()))
        self.start_lagrange_edit.setText('{0:.2f}'.format(self.occam_startup.start_lagrange))

    def _get_inv_folder(self):
        """
        create an inversion folder for each run
        """
        dir_path = os.path.join(self.station_dir, self.data_mode)
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)
            print 'Made directory {0}'.format(dir_path)
        
        dir_list = []
        for roots, dirs, files in os.walk(dir_path):
            dir_list.append(dirs)
            
        inv_num = len(dir_list[0])+1
        
        self.save_dir = os.path.join(self.station_dir, self.data_mode,
                                     'Inv_{0:02}'.format(inv_num))
        
    def run_occam(self):
        """
        write all the needed files and run occam then plot
        """
        
        self._get_inv_folder()
        
        if not os.path.isdir(self.save_dir):
            os.mkdir(self.save_dir)
            print 'Made directory {0}'.format(self.save_dir)
        
        # write data file
        self.occam_data.write_data_file(edi_file=self.edi_fn,
                                        save_path=self.save_dir,
                                        mode=self.data_mode,
                                        res_err=self.res_err,
                                        phase_err=self.phase_err,
                                        thetar=0)
                                        
        # write model file
        self.occam_model.write_model_file(save_path=self.save_dir)
        
        # write startup file
        self.occam_startup.data_fn = self.occam_data.data_fn
        self.occam_startup.model_fn = self.occam_model.model_fn
        self.occam_startup.write_startup_file(save_path=self.save_dir)
        
        warning_txt = '\n'.join(['Cannot find Occam1D executable. ',
                                 'Looked for {0}'.format(self.occam_exec),
                                 'Click Occam1D button and rerun.'])
        if not os.path.isfile(self.occam_exec):
            QtGui.QMessageBox.warning(self, "Warning", warning_txt)
            return 
                                   
        # run occam
        occam_run = occam1d.Run(startup_fn=self.occam_startup.startup_fn,
                                occam_path=self.occam_exec,
                                mode=self.data_mode) 
                                
        ini_resp_fn = os.path.join(self.save_dir, 
                                   '{0}_{1}.resp'.format(self.data_mode, 1))
        ini_model_fn = os.path.join(self.save_dir, 
                                    '{0}_{1}.iter'.format(self.data_mode, 1))

        ini_resp_fn = os.path.abspath(ini_resp_fn)
        ini_model_fn = os.path.abspath(ini_model_fn)                                            
        self.mpl_widget.plot_data(data_fn=self.occam_data.data_fn,
                                  resp_fn=ini_resp_fn,
                                  iter_fn=ini_model_fn,
                                  model_fn=self.occam_model.model_fn)
                                  
                                  
        self.l2_widget.plot_l2(dir_path=self.save_dir, 
                               model_fn=self.occam_model.model_fn)
                               
        # add iteration values to combo box 
        for ii in range(self.iter_combo_edit.count()):
            self.iter_combo_edit.removeItem(0)
        for ii in range(1, self.l2_widget.rms_arr.shape[0]):
            self.iter_combo_edit.addItem(str(ii))
            
        # resize the combo box to have width of max iteration
            
        self.iter_combo_edit.resize(self.iter_combo_edit.size())
        self.iter_combo_edit.update()
        self.iter_combo_edit.repaint()
        
                               
    def on_click(self, event):
        data_point = event.artist
        iteration = data_point.get_xdata()[event.ind]
        print 'Picked iteration {0}'.format(iteration)
        ini_resp_fn = os.path.join(self.save_dir, 
                                   '{0}_{1}.resp'.format(self.data_mode,
                                                         iteration))
        ini_model_fn = os.path.join(self.save_dir, 
                                    '{0}_{1}.iter'.format(self.data_mode,
                                                         iteration))

        ini_resp_fn = os.path.abspath(ini_resp_fn)
        ini_model_fn = os.path.abspath(ini_model_fn)                                            
        self.mpl_widget.plot_data(data_fn=self.occam_data.data_fn,
                                  resp_fn=ini_resp_fn,
                                  iter_fn=ini_model_fn,
                                  model_fn=self.occam_model.model_fn)
                                  
    def set_iteration(self, text):
        iteration = text
        rms = self.l2_widget.rms_arr['rms'][int(iteration)-1]
        roughness = self.l2_widget.rms_arr['roughness'][int(iteration)-1]
        print 'Iteration {0}, RMS={1:.2f}, Roughnes={2:.2f}'.format(
                iteration, rms, roughness)
        
        ini_resp_fn = os.path.join(self.save_dir, 
                                   '{0}_{1}.resp'.format(self.data_mode,
                                                         iteration))
        ini_model_fn = os.path.join(self.save_dir, 
                                    '{0}_{1}.iter'.format(self.data_mode,
                                                         iteration))

        ini_resp_fn = os.path.abspath(ini_resp_fn)
        ini_model_fn = os.path.abspath(ini_model_fn)                                            
        self.mpl_widget.plot_data(data_fn=self.occam_data.data_fn,
                                  resp_fn=ini_resp_fn,
                                  iter_fn=ini_model_fn,
                                  model_fn=self.occam_model.model_fn)
        
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
        
        self.subplot_wspace = .15
        self.subplot_hspace = .2
        self.subplot_right = .90
        self.subplot_left = .085
        self.subplot_top = .93
        self.subplot_bottom = .1
        
        self.fig = None
        self.axr = None
        self.axp = None
        self.axm = None
        
        self.res_limits = None
        self.phase_limits = (-5, 95)
        self.depth_scale = 'log'
        self.depth_units = 'km'
        self.depth_limits = None
        
        self.marker_data = 's'
        self.marker_data_color = 'k'
        self.marker_resp = 'h'
        self.marker_resp_color = 'b'
        self.marker_size =  2
        
        self.lw = .75
        self.ls = ':'
        self.e_capthick = .75
        self.e_capsize = 3
        
        self.font_size = 8
        
        self.setup_ui()
        
    def setup_ui(self):
        
        self.figure = Figure(dpi=200)
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
    
    def plot_data(self, data_fn=None, resp_fn=None, model_fn=None, 
                  iter_fn=None):
        """
        plot response and depth model
        """
        self.figure.clf()
        
        d_kwargs = {'ls':self.ls,
                    'marker':self.marker_data,
                    'ms':self.marker_size,
                    'mfc':self.marker_data_color,
                    'mec':self.marker_data_color,
                    'color':self.marker_data_color,
                    'ecolor':self.marker_data_color,
                    'picker':2,
                    'lw':self.lw,
                    'elinewidth':self.lw,
                    'capsize':self.e_capsize,
                    'capthick':self.e_capthick}
                    
        r_kwargs = {'ls':self.ls,
                    'marker':self.marker_resp,
                    'ms':self.marker_size,
                    'mfc':self.marker_resp_color,
                    'mec':self.marker_resp_color,
                    'color':self.marker_resp_color,
                    'ecolor':self.marker_resp_color,
                    'picker':2,
                    'lw':self.lw,
                    'elinewidth':self.lw,
                    'capsize':self.e_capsize,
                    'capthick':self.e_capthick}
        
        #make a grid of subplots
        gs=gridspec.GridSpec(6, 5, hspace=self.subplot_hspace, 
                             wspace=self.subplot_wspace)
        #subplot resistivity
        self.axr = self.figure.add_subplot(gs[:4, :4])
        
        #subplot for phase
        self.axp = self.figure.add_subplot(gs[4:,:4], sharex=self.axr)

        #subplot for model
        self.axm = self.figure.add_subplot(gs[:, 4])
        
        #-----------------------------------------------------------------
        #--> plot data apparent resistivity and phase-------------------------
        if data_fn is not None:
            d1 = occam1d.Data()
            d1.read_data_file(data_fn)
            
            #--> cut out missing data
            rxy = np.where(d1.res_te[0] != 0)[0]
            
            #--> TE mode Data 
            if len(rxy) > 0:
                rte = self.axr.errorbar(1./d1.freq[rxy],
                                        d1.res_te[0][rxy],
                                        yerr=d1.res_te[1][rxy],
                                        **d_kwargs)
                #legend_marker_list_te.append(rte[0])
                #legend_label_list_te.append('$Obs_{TE}$')
            else:
                pass
            
            #--> cut out missing data
            ryx = np.where(d1.res_tm[0] != 0)[0]
            
            #--> TE mode Data 
            if len(ryx) > 0:
                rtm = self.axr.errorbar(1./d1.freq[ryx],
                                        d1.res_tm[0][ryx],
                                        yerr=d1.res_tm[1][ryx],
                                        **d_kwargs)
                #legend_marker_list_te.append(rte[0])
                #legend_label_list_te.append('$Obs_{TE}$')
            else:
                pass
            #--------------------plot phase--------------------------------
            #cut out missing data points first
            pxy = np.where(d1.phase_te[0]!=0)[0]
            
            #--> TE mode data
            if len(pxy) > 0:
                self.axp.errorbar(1./d1.freq[pxy],
                                   d1.phase_te[0][pxy],
                                   yerr=d1.phase_te[1][pxy],
                                   **d_kwargs)
            else:
                pass
            
            #cut out missing data points first
            pyx = np.where(d1.phase_tm[0]!=0)[0]
            
            #--> TE mode data
            if len(pyx) > 0:
                self.axp.errorbar(1./d1.freq[pyx],
                                   d1.phase_tm[0][pyx],
                                   yerr=d1.phase_tm[1][pyx],
                                  **d_kwargs)
            else:
                pass
            
        #-----------------------------------------------------------------
        #--> plot data apparent resistivity and phase-------------------------
        if resp_fn is not None:
            r1 = occam1d.Data()
            r1.read_resp_file(resp_fn, data_fn=data_fn)
            
            #--> cut out missing data
            rxy = np.where(r1.res_te[2] != 0)[0]
            
            #--> TE mode Data 
            if len(rxy) > 0:
                rte = self.axr.errorbar(1./r1.freq[rxy],
                                        r1.res_te[2][rxy],
                                        yerr=None,
                                        **r_kwargs)
                #legend_marker_list_te.append(rte[0])
                #legend_label_list_te.append('$Obs_{TE}$')
            else:
                pass
            
            #--> cut out missing data
            ryx = np.where(r1.res_tm[2] != 0)[0]
            
            #--> TE mode Data 
            if len(ryx) > 0:
                rtm = self.axr.errorbar(1./r1.freq[ryx],
                                        r1.res_tm[2][ryx],
                                        yerr=None,
                                        **r_kwargs)
                #legend_marker_list_te.append(rte[0])
                #legend_label_list_te.append('$Obs_{TE}$')
            else:
                pass
            #--------------------plot phase--------------------------------
            #cut out missing data points first
            pxy = np.where(r1.phase_te[2]!=0)[0]
            
            #--> TE mode data
            if len(pxy) > 0:
                self.axp.errorbar(1./r1.freq[pxy],
                                   r1.phase_te[2][pxy],
                                   yerr=None,
                                   **r_kwargs)
            else:
                pass
            
            #cut out missing data points first
            pyx = np.where(r1.phase_tm[2]!=0)[0]
            
            #--> TE mode data
            if len(pyx) > 0:
                self.axp.errorbar(1./r1.freq[pyx],
                                   r1.phase_tm[2][pyx],
                                   yerr=None,
                                  **r_kwargs)
            else:
                pass
        
        
        #--> set axis properties-----------------------------------------------
        self.axr.set_xscale('log')
        self.axp.set_xscale('log')
        self.axr.set_yscale('log')        
        self.axr.grid(True, alpha=.75, which='both', 
                      color=(.75, .75, .75))
        plt.setp(self.axr.xaxis.get_ticklabels(),visible=False)
        self.axp.grid(True, alpha=.75, which='both', 
                      color=(.75, .75, .75))
        self.axp.yaxis.set_major_locator(MultipleLocator(15))
        self.axp.yaxis.set_minor_locator(MultipleLocator(3))
        
        if self.res_limits is not None:
            self.axr.set_ylim(self.res_limits)
            
        if self.phase_limits is not None:
            self.axp.set_ylim(self.phase_limits)
            
        self.axr.set_ylabel('App. Res. ($\Omega \cdot m$)',
                       fontdict={'size':self.font_size,'weight':'bold'})
        self.axp.set_ylabel('Phase (deg)',
                       fontdict={'size':self.font_size,'weight':'bold'})
        self.axp.set_xlabel('Period (s)',
                            fontdict={'size':self.font_size,'weight':'bold'})
        #plt.suptitle(self.title_str,fontsize=self.font_size+2,fontweight='bold')
        for ax in [self.axr, self.axp, self.axm]:
            ax.tick_params(axis='both', which='major', 
                           labelsize=self.font_size-2)
        
        #--> plot depth model--------------------------------------------------     
        if model_fn is not None:
            if self.depth_units == 'km':
                dscale = 1000.
            else:
                dscale = 1.
            
            #--> plot te models
            m1 = occam1d.Model()
            m1.read_iter_file(iter_fn, model_fn)
            plot_depth = m1.model_depth[1:]/dscale
            plot_model = abs(10**m1.model_res[1:,1])
            self.axm.semilogx(plot_model[::-1],
                              plot_depth[::-1],
                              ls='steps-',
                              color='b', 
                              lw=self.lw)
                           
            #if self.depth_limits == None:
            dmin = min(plot_depth)
            if dmin == 0:
                dmin = 1
            dmax = max(plot_depth)
            self.depth_limits = (dmin, dmax)
            
            self.axm.set_ylim(ymin=max(self.depth_limits), 
                              ymax=min(self.depth_limits))
            
        if self.depth_scale == 'log':
            self.axm.set_yscale('log')
        self.axm.set_ylabel('Depth ({0})'.format(self.depth_units),
                            fontdict={'size':self.font_size,'weight':'bold'})
        self.axm.set_xlabel('Resistivity ($\Omega \cdot m$)',
                       fontdict={'size':self.font_size,'weight':'bold'})
        self.axm.grid(True, which='both', alpha=.75, color=(.75, .75, .75))
        self.axm.yaxis.set_label_position('right')
        self.axm.yaxis.tick_right()
            
        self.mpl_widget.draw()
        
#==============================================================================
# plot L2
#==============================================================================
class PlotL2(QtGui.QWidget):
    """
    plot the l2 curve and it will be pickable for each iteration
    """
    
    def __init__(self, dir_path=None, model_fn=None):
        super(PlotL2, self).__init__()
        self.dir_path = dir_path 
        self.model_fn = model_fn

        self.fig_dpi = 200
        self.font_size = 8
        
        self.subplot_right = .90
        self.subplot_left = .085
        self.subplot_top = .86
        self.subplot_bottom = .15
        
        self.rms_lw = 1
        self.rms_marker = 'd'
        self.rms_color = 'k'
        self.rms_marker_size = 5
        self.rms_median_color = 'red'
        self.rms_mean_color = 'orange'
        
        self.rough_lw = .75
        self.rough_marker = 's'
        self.rough_color = 'b'
        self.rough_marker_size = 3
        self.rough_font_size =  8
        
        self.int = 1
        
        self.setup_ui()
        
        if self.dir_path is not None:
            self.plot_l2()
        
    def setup_ui(self):
        
        self.figure = Figure(dpi=200)
        self.l2_widget = FigureCanvas(self.figure)
        
        #self.l2_widget.mpl_connect('pick event', self.on_click)
        
        self.figure.subplots_adjust(left=self.subplot_left,
                                    right=self.subplot_right,
                                    bottom=self.subplot_bottom,
                                    top=self.subplot_top)
                                    
        #make sure the figure takes up the entire plottable space
        self.l2_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                     QtGui.QSizePolicy.Expanding)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.l2_widget, self)
        # set the layout for the plot
        mpl_vbox = QtGui.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.l2_widget)
        
        self.setLayout(mpl_vbox)
        
        self.l2_widget.updateGeometry()
        
    def _get_iter_list(self):
        """
        get all iteration files in dir_path        
        """
        
        if os.path.isdir(self.dir_path) == False:
            raise IOError('Could not find {0}'.format(self.dir_path))
            
        iter_list = [os.path.join(self.dir_path, fn) 
                     for fn in os.listdir(self.dir_path) 
                     if fn.find('.iter')>0]
            
        self.rms_arr = np.zeros(len(iter_list), 
                                dtype=np.dtype([('iteration', np.int),
                                                ('rms', np.float),
                                                ('roughness', np.float)]))             
        for ii, fn in enumerate(iter_list):
            m1 = occam1d.Model()
            m1.read_iter_file(fn, self.model_fn)
            self.rms_arr[ii]['iteration'] = int(m1.itdict['Iteration'])
            self.rms_arr[ii]['rms'] = float(m1.itdict['Misfit Value'])
            self.rms_arr[ii]['roughness'] = float(m1.itdict['Roughness Value'])
            
        self.rms_arr.sort(order='iteration')
        
    def plot_l2(self, dir_path=None, model_fn=None):
        """
        plot l2 curve rms vs iteration
        """
        self.figure.clf()
        
        if dir_path is not None:
            self.dir_path = dir_path
            
        if model_fn is not None:
            self.model_fn = model_fn
            
        self._get_iter_list()
        
        nr = self.rms_arr.shape[0]
        med_rms = np.median(self.rms_arr['rms'][1:])
        mean_rms = np.mean(self.rms_arr['rms'][1:])
        
        #set the dimesions of the figure
        plt.rcParams['font.size'] = self.font_size
        plt.rcParams['figure.subplot.left'] = self.subplot_left
        plt.rcParams['figure.subplot.right'] = self.subplot_right
        plt.rcParams['figure.subplot.bottom'] = self.subplot_bottom
        plt.rcParams['figure.subplot.top'] = self.subplot_top
        
        #make figure instance
        
        #make a subplot for RMS vs Iteration
        self.ax1 = self.figure.add_subplot(1, 1, 1)
        
        #plot the rms vs iteration
        l1, = self.ax1.plot(self.rms_arr['iteration'],
                            self.rms_arr['rms'],
                            '-k', 
                            lw=1,
                            marker='d',
                            ms=5,
                            picker=3)
        
        #plot the median of the RMS
        m1, = self.ax1.plot(self.rms_arr['iteration'],
                            np.repeat(med_rms, nr),
                            ls='--',
                            color=self.rms_median_color,
                            lw=self.rms_lw*.75)
        
        #plot the mean of the RMS
        m2, = self.ax1.plot(self.rms_arr['iteration'],
                            np.repeat(mean_rms, nr),
                            ls='--',
                            color=self.rms_mean_color,
                            lw=self.rms_lw*.75)
                            
        self.ax2 = self.ax1.twinx()
        l2, = self.ax2.plot(self.rms_arr['iteration'],
                            self.rms_arr['roughness'],
                            ls='-',
                            color=self.rough_color,
                            lw=self.rough_lw,
                            marker=self.rough_marker,
                            ms=self.rough_marker_size,
                            mfc=self.rough_color) 
    
        
        #make a legend
        self.figure.legend([l1, l2, m1, m2],
                        ['RMS', 'Roughness',
                         'Median_RMS={0:.2f}'.format(med_rms),
                         'Mean_RMS={0:.2f}'.format(mean_rms)],
                         ncol=4,
                         loc='upper center',
                         columnspacing=.25,
                         markerscale=.75,
                         handletextpad=.15,
                         borderaxespad=.02,
                         prop={'size':self.font_size})
                    
        #set the axis properties for RMS vs iteration
#        self.ax1.yaxis.set_minor_locator(MultipleLocator(.1))
        self.ax1.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax1.set_ylabel('RMS', 
                            fontdict={'size':self.font_size+2,
                                      'weight':'bold'})                                   
        self.ax1.set_xlabel('Iteration',
                            fontdict={'size':self.font_size+2,
                                      'weight':'bold'})
        self.ax1.grid(alpha=.25, which='both', lw=self.rough_lw)
        self.ax2.set_ylabel('Roughness',
                            fontdict={'size':self.font_size+2,
                                      'weight':'bold',
                                      'color':self.rough_color})
                                      
        self.ax1.set_ylim(np.floor(self.rms_arr['rms'][1:].min()), 
                          np.ceil(self.rms_arr['rms'][1:].max()))
        self.ax2.set_ylim(np.floor(self.rms_arr['roughness'][1:].min()),
                          np.ceil(self.rms_arr['roughness'][1:].max()))


        
        for t2 in self.ax2.get_yticklabels():
            t2.set_color(self.rough_color)
            
        self.l2_widget.draw()
        
    def on_click(self, event):
        data_point = event.artist
        iteration = data_point.get_xdata()[event.ind]        
        print iteration
        
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