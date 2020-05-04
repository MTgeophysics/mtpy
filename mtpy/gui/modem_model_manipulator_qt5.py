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

#from PyQt5 import QtCore, QtWidgets
try:
    from PyQt5 import QtCore, QtGui, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")


from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.widgets as widgets
from matplotlib.figure import Figure

import numpy as np
from scipy import signal

import mtpy.modeling.modem as modem

#==============================================================================
# Main Window
#==============================================================================
class ModEM_Model_Manipulator(QtWidgets.QMainWindow):
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

        self.menu_properties = self.menuBar().addMenu("Properties")
        self.menu_properties_cb_action = self.menu_properties.addAction("Resistivity Limits")
        self.menu_properties_cb_action.triggered.connect(self.set_res_limits)

        self.menu_tools = self.menuBar().addMenu("Tools")
        self.menu_tools_pad_action = self.menu_tools.addAction("Pad Fill")
        self.menu_tools_pad_action.triggered.connect(self.pad_fill)
        self.menu_tools_smooth_action = self.menu_tools.addAction("Smooth")
        self.menu_tools_smooth_action.triggered.connect(self.smooth)

        QtCore.QMetaObject.connectSlotsByName(self)

    def get_data_fn(self):
        """
        get the filename from a file dialogue

        """

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)')[0])

        self.model_widget.data_fn = fn

    def get_model_fn(self):
        """
        read in an existing model file
        """

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                       filter='*.rho')[0])

        self.model_widget.model_fn = fn

    def save_model_fn(self):
        """
        save the current mesh settings to a file
        """

        fn_dialog = QtWidgets.QFileDialog()
        save_fn = str(fn_dialog.getSaveFileName(
                                    caption='Choose ModEM model file',
                                    filter='*.rho')[0])

        sv_path = os.path.dirname(save_fn)
        sv_basename = os.path.basename(save_fn)
        self.model_widget.model_obj.write_model_file(save_path=sv_path,
                                                     model_fn_basename=sv_basename,
                                                     res_model=self.model_widget.new_res_model)

    def set_res_limits(self):
        """
        set the resistivity limits
        """

        self.res_popup = ResLimits(self.model_widget.res_limits[0],
                               self.model_widget.res_limits[1])
        #self.popup.show()
        self.res_popup.res_changed.connect(self.set_res)

    def set_res(self):
        self.model_widget.res_limits = (self.res_popup.res_min,
                                        self.res_popup.res_max)

    def pad_fill(self):
        self.model_widget.set_fill_params()

    def smooth(self):
        self.model_widget.set_smooth_params()

# =============================================================================
# Resistivity limits widget
# =============================================================================
class ResLimits(QtWidgets.QWidget):
    res_changed = QtCore.pyqtSignal()
    def __init__(self, res_limit_min, res_limit_max):
        super(ResLimits, self).__init__()
        self.res_min = res_limit_min
        self.res_max = res_limit_max

        self.setup_ui()

    def setup_ui(self):

        self.label = QtWidgets.QLabel("Resistivty (Log10)")
        self.res_min_label = QtWidgets.QLabel('min')
        self.res_min_edit = QtWidgets.QLineEdit()
        self.res_min_edit.setText('{0:0.5g}'.format(self.res_min))
        self.res_min_edit.editingFinished.connect(self.set_res_min)

        self.res_max_label = QtWidgets.QLabel('max')
        self.res_max_edit = QtWidgets.QLineEdit()
        self.res_max_edit.setText('{0:0.5g}'.format(self.res_max))
        self.res_max_edit.editingFinished.connect(self.set_res_max)

        grid_layout = QtWidgets.QGridLayout()
        grid_layout.addWidget(self.label, 1, 2, 2, 2)
        grid_layout.addWidget(self.res_min_label, 2, 1)
        grid_layout.addWidget(self.res_min_edit, 2, 2)
        grid_layout.addWidget(self.res_max_label, 3, 1)
        grid_layout.addWidget(self.res_max_edit, 3, 2)

        self.setLayout(grid_layout)
        self.show()

    def set_res_min(self):
        self.res_min = float(str(self.res_min_edit.text()))
        self.res_min_edit.setText('{0:.5g}'.format(self.res_min))
        self.res_changed.emit()

    def set_res_max(self):
        self.res_max = float(str(self.res_max_edit.text()))
        self.res_max_edit.setText('{0:.5g}'.format(self.res_max))
        self.res_changed.emit()

# =============================================================================
# Fill the padding cells widget
# =============================================================================
class PadFill(QtWidgets.QWidget):
    """
    widget to get pad filling parameters
    """
    apply_button_pushed = QtCore.pyqtSignal()
    undo_button_pushed = QtCore.pyqtSignal()
    def __init__(self):
        super(PadFill, self).__init__()

        self.n_pad = 3
        self.avg_range = 7

        self.setup_ui()

    def setup_ui(self):
        self.pad_edit = QtWidgets.QLineEdit()
        self.pad_edit.setText('{0:.0f}'.format(self.n_pad))
        self.pad_edit.editingFinished.connect(self.set_pad_num)
        self.pad_label = QtWidgets.QLabel('Number of Cells')

        self.avg_edit = QtWidgets.QLineEdit()
        self.avg_edit.setText('{0:.0f}'.format(self.avg_range))
        self.avg_edit.editingFinished.connect(self.set_avg_range)
        self.avg_label = QtWidgets.QLabel('Number of cells to find average')

        self.apply_button = QtWidgets.QPushButton()
        self.apply_button.setText('Apply')
        self.apply_button.clicked.connect(self.emit_apply_signal)

        self.undo_button = QtWidgets.QPushButton()
        self.undo_button.setText('Undo')
        self.undo_button.clicked.connect(self.emit_undo_signal)

        pad_layout = QtWidgets.QHBoxLayout()
        pad_layout.addWidget(self.pad_label)
        pad_layout.addWidget(self.pad_edit)

        avg_layout = QtWidgets.QHBoxLayout()
        avg_layout.addWidget(self.avg_label)
        avg_layout.addWidget(self.avg_edit)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(pad_layout)
        layout.addLayout(avg_layout)
        layout.addWidget(self.apply_button)
        layout.addWidget(self.undo_button)

        self.setLayout(layout)
        self.show()

    def set_pad_num(self):
        self.n_pad = int(str(self.pad_edit.text()))
        self.pad_edit.setText('{0:.0f}'.format(self.n_pad))

    def set_avg_range(self):
        self.avg_range = int(str(self.avg_edit.text()))
        self.avg_edit.setText('{0:.0f}'.format(self.avg_range))

    def emit_apply_signal(self):
        self.apply_button_pushed.emit()

    def emit_undo_signal(self):
        self.undo_button_pushed.emit()

# =============================================================================
# Smoothing widget
# =============================================================================
class Smooth(QtWidgets.QWidget):
    """
    smoothing widget
    """
    apply_button_pushed = QtCore.pyqtSignal()
    undo_button_pushed = QtCore.pyqtSignal()
    def __init__(self):
        super(Smooth, self).__init__()

        self.radius = 5
        self.sigma = 1

        self.setup_ui()

    def setup_ui(self):
        self.radius_edit = QtWidgets.QLineEdit()
        self.radius_edit.setText('{0:.0f}'.format(self.radius))
        self.radius_edit.editingFinished.connect(self.set_radius)
        self.radius_label = QtWidgets.QLabel('Gaussian Radius')

        self.sigma_edit = QtWidgets.QLineEdit()
        self.sigma_edit.setText('{0:.3f}'.format(self.sigma))
        self.sigma_edit.editingFinished.connect(self.set_sigma)
        self.sigma_label = QtWidgets.QLabel('Gaussian Full Width Half Max')

        self.apply_button = QtWidgets.QPushButton()
        self.apply_button.setText('Apply')
        self.apply_button.clicked.connect(self.emit_apply_signal)

        self.undo_button = QtWidgets.QPushButton()
        self.undo_button.setText('Undo')
        self.undo_button.clicked.connect(self.emit_undo_signal)


        radius_layout = QtWidgets.QHBoxLayout()
        radius_layout.addWidget(self.radius_label)
        radius_layout.addWidget(self.radius_edit)

        sigma_layout = QtWidgets.QHBoxLayout()
        sigma_layout.addWidget(self.sigma_label)
        sigma_layout.addWidget(self.sigma_edit)

        layout = QtWidgets.QVBoxLayout()
        layout.addLayout(radius_layout)
        layout.addLayout(sigma_layout)
        layout.addWidget(self.apply_button)
        layout.addWidget(self.undo_button)

        self.setLayout(layout)
        self.show()

    def set_radius(self):
        self.radius = int(str(self.radius_edit.text()))
        self.radius_edit.setText('{0:.0f}'.format(self.radius))

    def set_sigma(self):
        self.sigma = float(str(self.sigma_edit.text()))
        self.sigma_edit.setText('{0:.3f}'.format(self.sigma))

    def emit_apply_signal(self):
        self.apply_button_pushed.emit()

    def emit_undo_signal(self):
        self.undo_button_pushed.emit()

#==============================================================================
# Model Widget
#==============================================================================
class ModelWidget(QtWidgets.QWidget):
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
        self.north_line = None
        self.east_line = None
        self.north_line_xlist = None
        self.north_line_ylist = None
        self.east_line_xlist = None
        self.east_line_ylist = None

        self.map_ax = None
        self.north_ax = None
        self.east_ax = None
        self.cb_ax = None
        self.location_ax = None
        self.new_res_model = None

        self.units = 'km'
        self.scale = 1000.
        self.res_value = 100

        self.npad = 10
        self.avg_pad = 12

        self.cmap = 'jet_r'
        self._res_limits = (0, 4)
        self.map_copy_num = 1
        self.east_copy_num = 1
        self.north_copy_num = 1

        self.make_cb()

        self.ui_setup()

    def ui_setup(self):

        self.screen_size = QtWidgets.QDesktopWidget().screenGeometry()

        ## --> map view of the model
        self.map_figure = Figure()
        self.map_canvas = FigureCanvas(self.map_figure)
        self.map_canvas.mpl_connect('pick_event', self.map_on_pick)
        self.map_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Expanding)
        self.map_toolbar = NavigationToolbar(self.map_canvas, self)

        self.map_copy_down_button = QtWidgets.QPushButton('Copy Down (N layers)')
        self.map_copy_down_button.pressed.connect(self.map_copy_down)

        self.map_copy_up_button = QtWidgets.QPushButton('Copy Up (N layers)')
        self.map_copy_up_button.pressed.connect(self.map_copy_up)

        self.map_copy_number_edit = QtWidgets.QLineEdit()
        self.map_copy_number_edit.setText('{0:0.0f}'.format(self.map_copy_num))
        self.map_copy_number_edit.setMaximumWidth(35)
        self.map_copy_number_edit.editingFinished.connect(self.set_map_copy_num)
        self.map_copy_number_label = QtWidgets.QLabel('N')

        self.map_depth_label = QtWidgets.QLabel('Depth {0:>10.2f} {1}'.format(0,
                                                                    self.units))
        self.map_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.map_slider.valueChanged.connect(self.set_map_index)
        self.map_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.map_slider.setMinimum(0)
        self.map_slider.setMaximum(0)
        self.map_slider.setTickInterval(1)

        ## --> a N-S cross section that moves east to west
        self.east_figure= Figure()
        self.east_canvas = FigureCanvas(self.east_figure)
        self.east_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                       QtWidgets.QSizePolicy.Expanding)
        self.east_toolbar = NavigationToolbar(self.east_canvas, self)

        self.east_copy_west_button = QtWidgets.QPushButton('Copy West (N layers)')
        self.east_copy_west_button.pressed.connect(self.east_copy_west)

        self.east_copy_east_button = QtWidgets.QPushButton('Copy East (N layers)')
        self.east_copy_east_button.pressed.connect(self.east_copy_east)

        self.east_copy_number_edit = QtWidgets.QLineEdit()
        self.east_copy_number_edit.setText('{0:0.0f}'.format(self.east_copy_num))
        self.east_copy_number_edit.setMaximumWidth(35)
        self.east_copy_number_edit.editingFinished.connect(self.set_east_copy_num)
        self.east_copy_number_label = QtWidgets.QLabel('N')

        self.east_label = QtWidgets.QLabel('Easting {0:>10.2f} {1}'.format(0,
                                                                   self.units))
        self.east_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.east_slider.valueChanged.connect(self.set_east_index)
        self.east_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.east_slider.setMinimum(0)
        self.east_slider.setMaximum(0)
        self.east_slider.setTickInterval(1)

        ## --> a E-W cross section that moves N-S
        self.north_figure = Figure()
        self.north_canvas = FigureCanvas(self.north_figure)
        self.north_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                        QtWidgets.QSizePolicy.Expanding)
        self.north_toolbar = NavigationToolbar(self.north_canvas, self)

        self.north_copy_south_button = QtWidgets.QPushButton('Copy South (N layers)')
        self.north_copy_south_button.pressed.connect(self.north_copy_south)

        self.north_copy_north_button = QtWidgets.QPushButton('Copy North (N layers)')
        self.north_copy_north_button.pressed.connect(self.north_copy_north)

        self.north_copy_number_edit = QtWidgets.QLineEdit()
        self.north_copy_number_edit.setText('{0:0.0f}'.format(self.north_copy_num))
        self.north_copy_number_edit.setMaximumWidth(35)
        self.north_copy_number_edit.editingFinished.connect(self.set_north_copy_num)
        self.north_copy_number_label = QtWidgets.QLabel('N')

        self.north_label = QtWidgets.QLabel('Northing {0:>10.2f} m'.format(0,
                                                                    self.units))
        self.north_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.north_slider.valueChanged.connect(self.set_north_index)
        self.north_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.north_slider.setMinimum(0)
        self.north_slider.setMaximum(0)
        self.north_slider.setTickInterval(1)

        self.location_figure = Figure()
        self.location_canvas = FigureCanvas(self.location_figure)
        self.location_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                           QtWidgets.QSizePolicy.Expanding)
        self.location_toolbar = NavigationToolbar(self.location_canvas, self)
        self.location_canvas.mpl_connect('pick_event', self.location_pick)

        self.cb_figure = Figure()
        self.cb_canvas = FigureCanvas(self.cb_figure)
        self.cb_canvas.setMaximumWidth(int(self.screen_size.width()*.05))
        self.cb_canvas.setMinimumHeight(int(self.screen_size.height()*.85))
        self.cb_ax = self.cb_figure.add_axes([0.45, 0.005, 1.0, 1.0])
        self.cb_ax.pcolormesh(self.cb_x, self.cb_y, self.cb_bar,
                              vmin=self.res_limits[0],
                              vmax=self.res_limits[1],
                              cmap=self.cmap,
                              picker=5)
        self.cb_ax.set_yticks(np.arange(self._res_limits[0],
                                        self._res_limits[1],
                                        (self.res_limits[1]-self._res_limits[0])/10.))
        self.cb_ax.set_yticklabels(['{0:.4g}'.format(np.round(ii, 0)) for ii in
                                    np.logspace(self._res_limits[0],
                                                self._res_limits[1],
                                                num=11)])

        self.res_line, = self.cb_ax.plot([0, 1],
                                         [np.log10(self.res_value),
                                         np.log10(self.res_value)],
                                         lw=3,
                                         color='k',
                                         picker=5)
        self.cb_canvas.mpl_connect('button_press_event', self.on_res_pick)
        self.cb_ax.set_xticks([0, 1])
        self.cb_ax.set_xticklabels(['', ''])
        self.cb_ax.axis('tight')
        self.cb_canvas.draw()

        self.cb_line_edit = QtWidgets.QLineEdit()
        self.cb_line_edit.setMaximumWidth(140)
        self.cb_line_edit.setText('{0:.2f}'.format(self.res_value))
        self.cb_line_edit.editingFinished.connect(self.set_res_value)

        self.cb_label = QtWidgets.QLabel('Ohm-m')

        ##------------------------------------------------
        ## Layout

#        button_layout = QtWidgets.QHBoxLayout()
#        button_layout.addWidget(self.fill_outside_npad_label)
#        button_layout.addWidget(self.fill_outside_npad_edit)
#        button_layout.addWidget(self.fill_outside_avg_pad_label)
#        button_layout.addWidget(self.fill_outside_avg_pad_edit)
#        button_layout.addWidget(self.fill_outside_button)

        map_bottom_layout = QtWidgets.QHBoxLayout()
        map_bottom_layout.addWidget(self.map_depth_label)
        map_bottom_layout.addWidget(self.map_slider)
        map_top_layout = QtWidgets.QHBoxLayout()
        map_top_layout.addWidget(self.map_toolbar)
        map_top_layout.addWidget(self.map_copy_down_button)
        map_top_layout.addWidget(self.map_copy_up_button)
        map_top_layout.addWidget(self.map_copy_number_label)
        map_top_layout.addWidget(self.map_copy_number_edit)
        map_layout = QtWidgets.QVBoxLayout()
        map_layout.addLayout(map_top_layout)
        map_layout.addWidget(self.map_canvas)
        map_layout.addLayout(map_bottom_layout)

        east_bottom_layout = QtWidgets.QHBoxLayout()
        east_bottom_layout.addWidget(self.east_label)
        east_bottom_layout.addWidget(self.east_slider)
        east_top_layout = QtWidgets.QHBoxLayout()
        east_top_layout.addWidget(self.east_toolbar)
        east_top_layout.addWidget(self.east_copy_west_button)
        east_top_layout.addWidget(self.east_copy_east_button)
        east_top_layout.addWidget(self.east_copy_number_label)
        east_top_layout.addWidget(self.east_copy_number_edit)
        east_layout = QtWidgets.QVBoxLayout()
        east_layout.addLayout(east_top_layout)
        east_layout.addWidget(self.east_canvas)
        east_layout.addLayout(east_bottom_layout)

        north_bottom_layout = QtWidgets.QHBoxLayout()
        north_bottom_layout.addWidget(self.north_label)
        north_bottom_layout.addWidget(self.north_slider)
        north_top_layout = QtWidgets.QHBoxLayout()
        north_top_layout.addWidget(self.north_toolbar)
        north_top_layout.addWidget(self.north_copy_south_button)
        north_top_layout.addWidget(self.north_copy_north_button)
        north_top_layout.addWidget(self.north_copy_number_label)
        north_top_layout.addWidget(self.north_copy_number_edit)
        north_layout = QtWidgets.QVBoxLayout()
        north_layout.addLayout(north_top_layout)
        north_layout.addWidget(self.north_canvas)
        north_layout.addLayout(north_bottom_layout)

        location_layout = QtWidgets.QVBoxLayout()
        location_layout.addWidget(self.location_toolbar)
        location_layout.addWidget(self.location_canvas)

        grid_layout = QtWidgets.QGridLayout()
        grid_layout.addLayout(map_layout, 1, 1)
        grid_layout.addLayout(east_layout, 1, 2)
        grid_layout.addLayout(north_layout, 2, 1)
        grid_layout.addLayout(location_layout, 2, 2)

        cb_edit = QtWidgets.QVBoxLayout()
        cb_edit.addWidget(self.cb_label)
        cb_edit.addWidget(self.cb_line_edit)
        cb_layout = QtWidgets.QVBoxLayout()
        cb_layout.addWidget(self.cb_canvas)
        cb_layout.addLayout(cb_edit)

        final_layout = QtWidgets.QHBoxLayout()
        final_layout.addLayout(cb_layout)
        final_layout.addLayout(grid_layout)

        self.setLayout(final_layout)

    @property
    def res_limits(self):
        return self._res_limits

    @res_limits.setter
    def res_limits(self, res_limits):
        self._res_limits = res_limits

        self.redraw_cb()

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
        if self.map_ax is not None:
            self.redraw_plots()

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
        ## make a copy of the resistivity model to manipulate
        self.new_res_model = self.model_obj.res_model.copy()

        # set slider bar intervals
        # need the minus 1 cause we are using the value of the slider as
        # the index.
        self.map_slider.setMaximum(self.model_obj.plot_z.size-1)
        self.east_slider.setMaximum(self.model_obj.plot_east.size-1)
        self.north_slider.setMaximum(self.model_obj.plot_north.size-1)

        self.east_label.setText('{0:.2f}'.format(self.model_obj.grid_east[0]))
        self.north_label.setText('{0:.2f}'.format(self.model_obj.grid_north[0]))

        ##--------------plot the model-----------------------------------------
        ## get the grid coordinates first
        self.initialize_vectors()

        ## --> make map axes
        self.map_ax = self.map_figure.add_subplot(1, 1, 1)
        self.map_ax.set_xlabel('Easting {0}'.format(self.units))
        self.map_ax.set_ylabel('Northing {0}'.format(self.units))
        self.map_ax.set_aspect('equal')
        self.map_ax.plot(self.map_east_line_xlist,
                         self.map_east_line_ylist,
                         lw=.25,
                         color='k')
        self.map_ax.plot(self.map_north_line_xlist,
                         self.map_north_line_ylist,
                         lw=.25,
                         color='k')

        self.redraw_map()

        ## --> make EW cross section axes
        self.north_ax = self.north_figure.add_subplot(1, 1, 1,
                                                      sharex=self.map_ax,)
                                                      #aspect='equal')
        self.north_ax.plot(self.north_east_line_xlist,
                           self.north_east_line_ylist,
                           lw=.25,
                           color='k')
        self.north_ax.plot(self.north_z_line_xlist,
                           self.north_z_line_ylist,
                           lw=.25,
                           color='k')
        self.north_ax.set_xlabel('Easting {0}'.format(self.units))
        self.north_ax.set_ylabel('Depth {0}'.format(self.units))
        #self.north_ax.set_aspect('equal')
        self.redraw_north()
        # need to reverse the depth limits to plot properly
        z_lim = self.north_ax.get_ylim()
        self.north_ax.set_ylim([z_lim[-1], z_lim[0]])
        self.north_canvas.draw()

        ## --> make NS cross section axes
        self.east_ax = self.east_figure.add_subplot(1, 1, 1,
                                                    #aspect='equal',
                                                    sharex=self.map_ax,
                                                    sharey=self.north_ax)
        ## --> plot the mesh lines, this way only do it once
        self.east_ax.plot(self.east_north_line_xlist,
                          self.east_north_line_ylist,
                          lw=.25,
                          color='k')
        self.east_ax.plot(self.east_z_line_xlist,
                          self.east_z_line_ylist,
                          lw=.25,
                          color='k')
        self.east_ax.set_xlabel('Northing {0}'.format(self.units))
        self.east_ax.set_ylabel('Depth {0}'.format(self.units))
        #self.east_ax.set_aspect('equal')
        self.redraw_east()

        ## plot the location grid
        self.location_ax = self.location_figure.add_subplot(1, 1, 1,
                                                            aspect='equal')
        self.location_ax.set_xlabel('Easting {0}'.format(self.units))
        self.location_ax.set_ylabel('Northing {0}'.format(self.units))
        self.location_ax.set_aspect('equal')
        self.location_ax.plot(self.map_east_line_xlist,
                              self.map_east_line_ylist,
                              lw=.25,
                              color='k',
                              picker=3)
        self.location_ax.plot(self.map_north_line_xlist,
                              self.map_north_line_ylist,
                              lw=.25,
                              color='k',
                              picker=3)

        self.location_ax.set_xlim((self.model_obj.grid_east[self.model_obj.pad_east]/self.scale,
                                   self.model_obj.grid_east[-self.model_obj.pad_east]/self.scale))
        self.location_ax.set_ylim((self.model_obj.grid_north[self.model_obj.pad_north]/self.scale,
                                   self.model_obj.grid_north[-self.model_obj.pad_north]/self.scale))

        # make lines that can move around
        self.east_line = self.location_ax.plot([self.model_obj.grid_east[self.east_index]/self.scale,
                                                self.model_obj.grid_east[self.east_index]/self.scale],
                                                [self.model_obj.grid_north.min()/self.scale,
                                                 self.model_obj.grid_north.max()/self.scale],
                                                 'g',
                                                 lw=2)[0]
        self.north_line = self.location_ax.plot([self.model_obj.grid_east.min()/self.scale,
                                                self.model_obj.grid_east.max()/self.scale],
                                                [self.model_obj.grid_north[self.north_index]/self.scale,
                                                 self.model_obj.grid_north[self.north_index]/self.scale],
                                                 'b',
                                                 lw=2)[0]
        if self.data_fn is not None:
            self.location_ax.scatter(self.data_obj.station_locations.rel_east/self.scale,
                                     self.data_obj.station_locations.rel_north/self.scale,
                                     marker='v',
                                     c='k',
                                     s=10)
            self.location_ax.set_xlim((self.data_obj.station_locations.rel_east.min()/self.scale - 1,
                                       self.data_obj.station_locations.rel_east.max()/self.scale + 1))
            self.location_ax.set_ylim((self.data_obj.station_locations.rel_north.min()/self.scale - 1,
                                       self.data_obj.station_locations.rel_north.max()/self.scale+ 1))

        self.location_canvas.draw()

        #make a rectangular selector
        self.map_selector = widgets.RectangleSelector(self.map_ax,
                                                      self.map_on_pick,
                                                      drawtype='box',
                                                      useblit=True)
        self.east_selector = widgets.RectangleSelector(self.east_ax,
                                                       self.east_on_pick,
                                                       drawtype='box',
                                                       useblit=True)
        self.north_selector = widgets.RectangleSelector(self.north_ax,
                                                        self.north_on_pick,
                                                        drawtype='box',
                                                        useblit=True)


    def undo(self):
        """
        reset the resistivity model to its original
        """
        self.new_res_model = self.model_obj.res_model.copy()

    def initialize_vectors(self):
        """
        get all the plotting vectors
        """

        ### --> get mesh grids for plotting pcolormeshes
        self.plot_east_map, self.plot_north_map = np.meshgrid(self.model_obj.grid_east/self.scale,
                                                              self.model_obj.grid_north/self.scale,
                                                              indexing='ij')
        self.plot_east_z, self.plot_z_east = np.meshgrid(self.model_obj.grid_east/self.scale,
                                                         self.model_obj.grid_z/self.scale,
                                                         indexing='ij')
        self.plot_north_z, self.plot_z_north = np.meshgrid(self.model_obj.grid_north/self.scale,
                                                           self.model_obj.grid_z/self.scale,
                                                           indexing='ij')

        # get line lists for plotting grid lines
        ## --> map view
        self.map_east_line_xlist = []
        self.map_east_line_ylist = []
        for xx in self.model_obj.grid_east:
            self.map_east_line_xlist.extend([xx / self.scale, xx / self.scale])
            self.map_east_line_xlist.append(None)
            self.map_east_line_ylist.extend([self.model_obj.grid_north.min() / self.scale,
                                             self.model_obj.grid_north.max() / self.scale])
            self.map_east_line_ylist.append(None)

        self.map_north_line_xlist = []
        self.map_north_line_ylist = []
        for yy in self.model_obj.grid_north:
            self.map_north_line_xlist.extend([self.model_obj.grid_east.min() / self.scale,
                                              self.model_obj.grid_east.max() / self.scale])
            self.map_north_line_xlist.append(None)
            self.map_north_line_ylist.extend([yy / self.scale, yy / self.scale])
            self.map_north_line_ylist.append(None)

        ##--> NS cross section that move E-W
        self.east_north_line_xlist = []
        self.east_north_line_ylist = []
        for xx in self.model_obj.grid_north:
            self.east_north_line_xlist.extend([xx/self.scale, xx/self.scale])
            self.east_north_line_xlist.append(None)
            self.east_north_line_ylist.extend([self.model_obj.grid_z.min()/self.scale,
                                               self.model_obj.grid_z.max()/self.scale])
            self.east_north_line_ylist.append(None)

        self.east_z_line_xlist = []
        self.east_z_line_ylist = []
        for yy in self.model_obj.grid_z:
            self.east_z_line_xlist.extend([self.model_obj.grid_north.min() / self.scale,
                                           self.model_obj.grid_north.max() / self.scale])
            self.east_z_line_xlist.append(None)
            self.east_z_line_ylist.extend([yy / self.scale, yy / self.scale])
            self.east_z_line_ylist.append(None)

        ##--> EW cross section that move N-S
        self.north_east_line_xlist = []
        self.north_east_line_ylist = []
        for xx in self.model_obj.grid_east:
            self.north_east_line_xlist.extend([xx/self.scale, xx/self.scale])
            self.north_east_line_xlist.append(None)
            self.north_east_line_ylist.extend([self.model_obj.grid_z.min()/self.scale,
                                               self.model_obj.grid_z.max()/self.scale])
            self.north_east_line_ylist.append(None)

        self.north_z_line_xlist = []
        self.north_z_line_ylist = []
        for yy in self.model_obj.grid_z:
            self.north_z_line_xlist.extend([self.model_obj.grid_east.min() / self.scale,
                                           self.model_obj.grid_east.max() / self.scale])
            self.north_z_line_xlist.append(None)
            self.north_z_line_ylist.extend([yy / self.scale, yy / self.scale])
            self.north_z_line_ylist.append(None)


    def make_cb(self):
        res = np.arange(self.res_limits[0],
                        self.res_limits[1],
                       (self.res_limits[1]-self.res_limits[0])/256.)
        self.cb_x, self.cb_y = np.meshgrid(np.array([0, 1]), res, indexing='ij')
        self.cb_bar = np.zeros((2, 256))
        self.cb_bar[:, :] = res

    def on_res_pick(self, event):
        try:
            y_data = 10**event.ydata
            y_data = np.log10(np.round(y_data,
                                       -int(np.floor(np.log10(y_data)))))

            self.res_line.set_xdata([0, 1])
            self.res_line.set_ydata([y_data, y_data])
            self.cb_canvas.draw()
            self.res_value = 10**y_data

            self.cb_line_edit.setText('{0:.2f}'.format(self.res_value))

        except TypeError:
            return

    def set_res_value(self):
        self.res_value = float(self.cb_line_edit.text())
        self.res_line.set_ydata([np.log10(self.res_value),
                                 np.log10(self.res_value)])
        self.cb_canvas.draw()
        self.cb_line_edit.setText('{0:.2f}'.format(self.res_value))

    def set_map_index(self):
        self.map_index = int(self.map_slider.value())
        depth = self.model_obj.grid_z[self.map_index]/self.scale
        self.map_depth_label.setText('Depth {0:>10.2f} {1}'.format(depth,
                                                                self.units))

        self.redraw_map()

    def redraw_map(self):
        """
        redraw map view
        """
        self.map_ax.pcolormesh(self.plot_east_map,
                               self.plot_north_map,
                               np.log10(self.new_res_model[:, :, self.map_index].T),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        if self.data_fn is not None:
            self.map_ax.scatter(self.data_obj.station_locations.rel_east/self.scale,
                                self.data_obj.station_locations.rel_north/self.scale,
                                marker='v',
                                c='k',
                                s=10)
        self.map_canvas.draw()

    def set_east_index(self):
        self.east_index = int(self.east_slider.value())
        easting = self.model_obj.grid_east[self.east_index]/self.scale

        self.east_label.setText('Easting {0:>10.2f} {1}'.format(easting,
                                                                self.units))
        self.redraw_east()
        self.redraw_location()

    def redraw_east(self):
        """
        redraw east view
        """
        xlim = self.east_ax.get_xlim()
        ylim = self.east_ax.get_ylim()

        self.east_ax.cla()
        self.east_ax.plot(self.east_north_line_xlist,
                          self.east_north_line_ylist,
                          lw=.25,
                          color='k')
        self.east_ax.plot(self.east_z_line_xlist,
                          self.east_z_line_ylist,
                          lw=.25,
                          color='k')
        self.east_ax.pcolormesh(self.plot_north_z,
                               self.plot_z_north,
                               np.log10(self.new_res_model[:, self.east_index, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])

        line = self.get_stations_east()
        if line is not None:
            self.east_ax.scatter(line['rel_north']/self.scale,
                                 line['rel_elev']/self.scale,
                                 marker='v',  c='cyan', s=50,
                                 edgecolors='k')

            self.location_ax.scatter(line['rel_east']/self.scale,
                                     line['rel_north']/self.scale,
                                     marker='v',  c='k', s=30,
                                     edgecolors='cyan')
            self.location_canvas.draw()

            for ss in line:
                self.east_ax.text(ss['rel_north']/self.scale,
                                  ss['rel_elev']/self.scale - .2,
                                  ss['station'],
                                  va='bottom', ha='center',
                                  fontdict={'weight':'bold',
                                            'size':10},
                                  clip_on=True,
                                  bbox={'boxstyle':"square",
                                        'ec':'k',
                                        'fc':'w'})

        self.east_ax.set_xlabel('Northing {0}'.format(self.units))
        self.east_ax.set_ylabel('Depth {0}'.format(self.units))
        self.east_ax.set_ylim(ylim)
        self.east_ax.set_xlim(xlim)
        self.east_canvas.draw()

    def redraw_location(self):
        """
        redraw the location map with the indication lines on it
        """

        self.location_ax.cla()
        self.location_ax.plot(self.map_east_line_xlist,
                              self.map_east_line_ylist,
                              lw=.25,
                              color='k',
                              picker=3)
        self.location_ax.plot(self.map_north_line_xlist,
                              self.map_north_line_ylist,
                              lw=.25,
                              color='k',
                              picker=3)
        # make lines that can move around
        self.east_line = self.location_ax.plot([self.model_obj.grid_east[self.east_index]/self.scale,
                                                self.model_obj.grid_east[self.east_index]/self.scale],
                                                [self.model_obj.grid_north.min()/self.scale,
                                                 self.model_obj.grid_north.max()/self.scale],
                                                 'g',
                                                 lw=2)[0]
        self.north_line = self.location_ax.plot([self.model_obj.grid_east.min()/self.scale,
                                                self.model_obj.grid_east.max()/self.scale],
                                                [self.model_obj.grid_north[self.north_index]/self.scale,
                                                 self.model_obj.grid_north[self.north_index]/self.scale],
                                                 'b',
                                                 lw=2)[0]
        self.location_ax.set_xlim((self.model_obj.grid_east[self.model_obj.pad_east]/self.scale,
                                   self.model_obj.grid_east[-self.model_obj.pad_east]/self.scale))
        self.location_ax.set_ylim((self.model_obj.grid_north[self.model_obj.pad_north]/self.scale,
                                   self.model_obj.grid_north[-self.model_obj.pad_north]/self.scale))

        if self.data_fn is not None:
            self.location_ax.scatter(self.data_obj.station_locations.rel_east/self.scale,
                                     self.data_obj.station_locations.rel_north/self.scale,
                                     marker='v',
                                     c='k',
                                     edgecolors='k',
                                     s=30)
            self.location_ax.set_xlim((self.data_obj.station_locations.rel_east.min()/self.scale - 1,
                                       self.data_obj.station_locations.rel_east.max()/self.scale + 1))
            self.location_ax.set_ylim((self.data_obj.station_locations.rel_north.min()/self.scale - 1,
                                       self.data_obj.station_locations.rel_north.max()/self.scale + 1))

        self.location_ax.set_xlabel('Easting {0}'.format(self.units))
        self.location_ax.set_ylabel('Northing {0}'.format(self.units))

        self.location_canvas.draw()

    def set_north_index(self):
        self.north_index = int(self.north_slider.value())
        northing = self.model_obj.grid_north[self.north_index]/self.scale
        self.north_label.setText('Northing {0:>10.2f} {1}'.format(northing,
                                                                self.units))

        self.redraw_north()
        self.redraw_location()

    def redraw_north(self):
        """
        redraw north view
        """
        ylim = self.north_ax.get_ylim()
        xlim = self.north_ax.get_xlim()
        self.north_ax.cla()
        self.north_ax.plot(self.north_east_line_xlist,
                           self.north_east_line_ylist,
                           lw=.25,
                           color='k')
        self.north_ax.plot(self.north_z_line_xlist,
                           self.north_z_line_ylist,
                           lw=.25,
                           color='k')
        self.north_ax.pcolormesh(self.plot_east_z,
                                 self.plot_z_east,
                                 np.log10(self.new_res_model[self.north_index, :, :]),
                                 cmap=self.cmap,
                                 vmin=self.res_limits[0],
                                 vmax=self.res_limits[1])

        line = self.get_stations_north()

        if line is not None:
            self.north_ax.scatter(line['rel_east']/self.scale,
                                  line['rel_elev']/self.scale,
                                  marker='v',  c='cyan', s=50,
                                  edgecolors='k')
            self.location_ax.scatter(line['rel_east']/self.scale,
                                     line['rel_north']/self.scale,
                                     marker='v',  c='k', s=30,
                                     edgecolors='cyan')
            self.location_canvas.draw()

            for ss in line:
                self.north_ax.text(ss['rel_east']/self.scale,
                                   ss['rel_elev']/self.scale - .2,
                                   ss['station'],
                                   va='bottom', ha='center',
                                   fontdict={'weight':'bold',
                                            'size':10},
                                  clip_on=True,
                                  bbox={'boxstyle':"square",
                                        'ec':'k',
                                        'fc':'w'})

        self.north_ax.set_xlabel('Easting {0}'.format(self.units))
        self.north_ax.set_ylabel('Elevation {0}'.format(self.units))
        self.north_ax.set_ylim(ylim)
        self.north_ax.set_xlim(xlim)
        self.north_canvas.draw()

    def get_stations_north(self):
        """
        get stations close to the line north
        Returns
        -------
        None.

        """
        ymin = self.model_obj.grid_north[self.north_index] - \
                self.model_obj.cell_size_north
        ymax = self.model_obj.grid_north[self.north_index] + \
                self.model_obj.cell_size_north
        if self.data_fn is not None:
            s_find = np.where((self.data_obj.data_array['rel_north'] >= ymin) &
                              (self.data_obj.data_array['rel_north'] <= ymax))
            return self.data_obj.data_array[s_find[0]]

        else:
            return None

    def get_stations_east(self):
        """
        get stations close to the line east
        Returns
        -------
        None.

        """
        ymin = self.model_obj.grid_east[self.east_index] - \
                self.model_obj.cell_size_east
        ymax = self.model_obj.grid_east[self.east_index] + \
                self.model_obj.cell_size_east
        if self.data_fn is not None:
            s_find = np.where((self.data_obj.data_array['rel_east'] >= ymin) &
                              (self.data_obj.data_array['rel_east'] <= ymax))
            return self.data_obj.data_array[s_find[0]]

        else:
            return None

    def redraw_cb(self):
        """
        redraw the colorbar
        """
        self.cb_ax.cla()
        self.make_cb()
        self.cb_ax.pcolormesh(self.cb_x, self.cb_y, self.cb_bar,
                              vmin=self.res_limits[0],
                              vmax=self.res_limits[1],
                              cmap=self.cmap,
                              picker=5)
        self.cb_ax.set_yticks(np.arange(self._res_limits[0],
                                        self._res_limits[1],
                                        (self.res_limits[1]-self._res_limits[0])/10))
        self.cb_ax.set_yticklabels(['{0:.4g}'.format(np.round(ii, 0)) for ii in
                                    np.logspace(self._res_limits[0],
                                                self._res_limits[1],
                                                num=11)])
        self.res_line, = self.cb_ax.plot([0, 1],
                                         [np.log10(self.res_value),
                                          np.log10(self.res_value)],
                                          lw=3,
                                          color='k',
                                          picker=5)
        self.cb_canvas.draw()
        self.redraw_plots()

    def map_on_pick(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        x_change = self._get_change_index(x1, x2, self.model_obj.grid_east)
        y_change = self._get_change_index(y1, y2, self.model_obj.grid_north)

        #reset values of resistivity
        for xx in x_change:
            for yy in y_change:
                if self.model_obj.res_model[yy, xx, self.map_index] < 1E10:
                    self.new_res_model[yy, xx, self.map_index] = self.res_value
        self.redraw_plots()

    def east_on_pick(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        x_change = self._get_change_index(x1, x2, self.model_obj.grid_north)
        y_change = self._get_change_index(y1, y2, self.model_obj.grid_z)

        #reset values of resistivity
        for xx in x_change:
            for yy in y_change:
                if self.model_obj.res_model[xx, self.east_index, yy] < 1E10:
                    self.new_res_model[xx, self.east_index, yy] = self.res_value

        self.redraw_plots()

    def north_on_pick(self, eclick, erelease):
        """
        on selecting a rectangle change the colors to the resistivity values
        """
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        x_change = self._get_change_index(x1, x2, self.model_obj.grid_east)
        y_change = self._get_change_index(y1, y2, self.model_obj.grid_z)

        #reset values of resistivity
        for xx in x_change:
            for yy in y_change:
                if self.model_obj.res_model[self.north_index, xx, yy] < 1E10:
                    self.new_res_model[self.north_index, xx, yy] = self.res_value

        self.redraw_plots()

    def _get_change_index(self, y1, y2, grid_dir):
        """
        get the index value of the points to be changed in north direction

        need to flip the index because the plot is flipped

        """

        if y1 < y2:
            ychange = np.where((grid_dir/self.scale > y1) & \
                               (grid_dir/self.scale < y2))[0]
            if len(ychange) == 0:
                ychange = np.where(grid_dir/self.scale >= y1)[0][0]-1
                return [ychange]

        elif y1 > y2:
            ychange = np.where((grid_dir/self.scale < y1) & \
                               (grid_dir/self.scale > y2))[0]
            if len(ychange) == 0:
                ychange = np.where(grid_dir/self.scale >= y2)[0][0]-1
                return [ychange]

        ychange -= 1
        ychange = np.append(ychange, ychange[-1]+1)

        return ychange

    def redraw_plots(self):
        """
        redraw all plots
        """
        self.redraw_map()
        self.redraw_east()
        self.redraw_north()
        self.redraw_location()

    def location_pick(self, event):
        """
        change the index of the location line either e-w or n-s.
        """
        if event.mouseevent.button == 1:
            data_point = event.mouseevent

            # figure out if data point is close to x or y
            try:
                x_index = np.where(self.model_obj.grid_east/self.scale >= data_point.xdata)[0][0]
            except IndexError:
                return
            try:
                y_index = np.where(self.model_obj.grid_north/self.scale >= data_point.ydata)[0][0]
            except IndexError:
                return
            dx = np.abs(data_point.xdata-self.model_obj.grid_east[x_index]/self.scale)
            dy = np.abs(data_point.ydata-self.model_obj.grid_north[y_index]/self.scale)

            if dx < dy:
                self.east_index = x_index
                self.east_line.set_xdata([self.model_obj.grid_east[self.east_index]/self.scale,
                                          self.model_obj.grid_east[self.east_index]/self.scale])
                self.east_slider.setValue(self.east_index)
                self.east_slider.triggerAction(QtWidgets.QAbstractSlider.SliderMove)

            elif dx > dy:
                self.north_index = y_index
                self.north_line.set_ydata([self.model_obj.grid_north[self.north_index]/self.scale,
                                          self.model_obj.grid_north[self.north_index]/self.scale])
                self.north_slider.setValue(self.north_index)
                self.north_slider.triggerAction(QtWidgets.QAbstractSlider.SliderMove)

    def set_fill_params(self):
        """
        fill the area outside the main station area
        """

        self.avg_widget = PadFill()
        self.avg_widget.apply_button_pushed.connect(self.fill_outside_area)
        self.avg_widget.undo_button_pushed.connect(self.undo_tools)


    def fill_outside_area(self):
        """
        fill areas outside given area
        """
        avg_range = self.avg_widget.avg_range
        n_pad = self.avg_widget.n_pad

        x_range = np.append(np.arange(avg_range), np.arange(-avg_range, 0, 1))
        y_range = np.append(np.arange(avg_range), np.arange(-avg_range, 0, 1))

        x_index, y_index = np.meshgrid(x_range, y_range)
        for zz in range(self.new_res_model.shape[2]):
            self.new_res_model[:, :, zz] = self.mask_elevation_cells(self.new_res_model[:, :, zz])
            avg_res_value = np.mean([np.median(self.new_res_model[x_index, y_index, zz]),
                                     np.median(self.new_res_model[avg_range:-avg_range, 0:avg_range, zz]),
                                     np.median(self.new_res_model[avg_range:-avg_range, -avg_range:, zz]),
                                     np.median(self.new_res_model[0:avg_range, avg_range:-avg_range, zz]),
                                     np.median(self.new_res_model[-avg_range:, avg_range:-avg_range, zz])])

            #self.new_res_model[x_index, y_index, zz] = avg_res_value
            self.new_res_model[:, -n_pad:, zz] = avg_res_value
            self.new_res_model[:, 0:n_pad, zz] = avg_res_value
            self.new_res_model[0:n_pad, n_pad:-n_pad, zz] = avg_res_value
            self.new_res_model[-n_pad:, n_pad:-n_pad, zz] = avg_res_value

        ### need to elevation
        elev_index = np.where(self.model_obj.res_model > 1E10)
        self.new_res_model[elev_index] = 1E12

        self.redraw_plots()

    def undo_tools(self):
        """
        undo fill outside area
        """

        self.undo()
        self.redraw_plots()

    def set_smooth_params(self):
        """
        set smoothing parameters
        """
        self.smooth_widget = Smooth()
        self.smooth_widget.apply_button_pushed.connect(self.smooth_model)
        self.smooth_widget.undo_button_pushed.connect(self.undo_tools)

    def smooth_model(self):
        """
        smooth model with 2D gaussian filter
        """
        radius = self.smooth_widget.radius
        sigma = self.smooth_widget.sigma**2

        gx, gy = np.mgrid[-radius:radius+1,
                          -radius:radius+1]

        gauss = (1./(2*np.pi*sigma))*np.exp(-((gx**2)+(gy**2))/(2*sigma))



        for zz in range(self.new_res_model.shape[2]):
            ### need to take into account elevation cells
            self.new_res_model[:, :, zz] = self.mask_elevation_cells(self.new_res_model[:, :, zz])
            self.new_res_model[:, :, zz] = signal.convolve(self.new_res_model[:, :, zz],
                                                           gauss,
                                                           mode='same')
        ### need to elevation
        elev_index = np.where(self.model_obj.res_model > 1E10)
        self.new_res_model[elev_index] = 1E12

        self.redraw_plots()

    def mask_elevation_cells(self, res_array):
        """
        remove the effect of elevation cells
        """
        mean_value = res_array[np.where(res_array <1E10)].mean()
        res_array[np.where(res_array > 1E10)] = mean_value

        return res_array

    def map_copy_down(self):
        """
        copy the current map down the number of layers given
        """
        o_shape = (self.new_res_model.shape[0],
                   self.new_res_model.shape[1],
                   1)
        # need to add 1 to the index to make sure that copy number is observed
        copy_index = self.map_index+(self.map_copy_num+1)
        if copy_index > self.new_res_model.shape[2]:
            copy_index = self.new_res_model.shape[2]
#        for m_index in range(self.map_index, copy_index, 1):
#            nax = np.where(self.new_res_model[:, :, m_index] <1E12)
#            self.new_res_model[nax] =
#        na_index = np.where(self.new_res_model[:, :, self.map_index] < 1E10)
#        print(na_index)
        self.new_res_model[:, : , self.map_index:copy_index] = \
            self.new_res_model[:, :, self.map_index].reshape(o_shape)

        #self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_map()

    def map_copy_up(self):
        """
        copy the current map up the number of layers given
        """
        o_shape = (self.new_res_model.shape[0],
                   self.new_res_model.shape[1],
                   1)

        copy_index = self.map_index-(self.map_copy_num+1)
        if copy_index < 0:
            copy_index = 0
        self.new_res_model[:, :, copy_index:self.map_index] = \
            self.new_res_model[:, :, self.map_index].reshape(o_shape)
        self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_map()

    def set_map_copy_num(self):
        """
        set number of layers to copy
        """
        self.map_copy_num = int(round(float(str(self.map_copy_number_edit.text()))))
        self.map_copy_number_edit.setText('{0:.0f}'.format(self.map_copy_num))

    def east_copy_east(self):
        """
        copy the current cross section east by east_copy_num
        """
        o_shape = (self.new_res_model.shape[0],
                   1,
                   self.new_res_model.shape[2])

        copy_index = self.east_index+(self.east_copy_num+1)
        if copy_index > self.new_res_model.shape[1]:
            copy_index = self.new_res_model.shape[1]

        self.new_res_model[:, self.east_index:copy_index, :] = \
            self.new_res_model[:, self.east_index, :].reshape(o_shape)

        self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_east()

    def east_copy_west(self):
        """
        copy the current cross section west by east_copy_num
        """
        o_shape = (self.new_res_model.shape[0],
                   1,
                   self.new_res_model.shape[2])

        copy_index = self.east_index-(self.east_copy_num+1)
        if copy_index < 0:
            copy_index = 0

        self.new_res_model[:, copy_index:self.east_index, :] = \
            self.new_res_model[:, self.east_index, :].reshape(o_shape)

        self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_east()

    def set_east_copy_num(self):
        """
        set the number of layers to copy in the east direction
        """
        self.east_copy_num = int(round(float(str(self.east_copy_number_edit.text()))))
        self.east_copy_number_edit.setText('{0:.0f}'.format(self.east_copy_num))

    def north_copy_south(self):
        """
        copy the current cross section south by north_copy_num
        """
        o_shape = (1,
                   self.new_res_model.shape[1],
                   self.new_res_model.shape[2])

        copy_index = self.north_index-(self.north_copy_num+1)
        if copy_index > self.new_res_model.shape[0]:
            copy_index = self.new_res_model.shape[0]

        self.new_res_model[copy_index:self.north_index, :, :] = \
            self.new_res_model[self.north_index, :, :].reshape(o_shape)
        self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_north()

    def north_copy_north(self):
        """
        copy the current cross section north by north_copy_num
        """
        o_shape = (1,
                   self.new_res_model.shape[1],
                   self.new_res_model.shape[2])

        copy_index = self.north_index+(self.north_copy_num+1)
        if copy_index > self.new_res_model.shape[0]:
            copy_index = self.new_res_model.shape[0]

        self.new_res_model[self.north_index:copy_index:, :, :] = \
            self.new_res_model[self.north_index, :, :].reshape(o_shape)
        self.new_res_model[np.where(self.model_obj.res_model > 1E10)] = 1E12

        self.redraw_north()

    def set_north_copy_num(self):
        """
        set the number of layers to copy in the north direction
        """
        self.north_copy_num = int(round(float(str(self.north_copy_number_edit.text()))))
        self.north_copy_number_edit.setText('{0:.0f}'.format(self.north_copy_num))



#==============================================================================
#  DEFINE MAIN
#==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = ModEM_Model_Manipulator()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()


