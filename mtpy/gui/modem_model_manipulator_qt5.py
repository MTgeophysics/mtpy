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

from mtpy.gui.my_stream import MyStream
import mtpy.modeling.modem_new as modem

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

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)'))

        self.model_widget.data_fn = fn

    def get_model_fn(self):
        """
        read in an existing model file
        """

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                       filter='*.rho'))

        self.model_widget.model_fn = fn

    def save_model_fn(self):
        """
        save the current mesh settings to a file
        """

        fn_dialog = QtWidgets.QFileDialog()
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

        self.units = 'km'
        self.scale = 1000.
        self.res_value = 100
        
        self.npad = 10
        self.avg_pad = 12

        self.cmap = 'jet_r'
        self.res_limits = (0, 4)

        self.make_cb()
        
        self.ui_setup()

    def ui_setup(self):
        
        self.screen_size = QtWidgets.QDesktopWidget().screenGeometry()
        
        self.fill_outside_avg_pad_label = QtWidgets.QLabel('Avg Cells Max')
        self.fill_outside_avg_pad_edit = QtWidgets.QLineEdit('{0:.0f}'.format(self.avg_pad))
        self.fill_outside_avg_pad_edit.editingFinished.connect(self.set_avg_pad)
        
        self.fill_outside_npad_label = QtWidgets.QLabel('Num Cells to Pad')
        self.fill_outside_npad_edit = QtWidgets.QLineEdit('{0:.0f}'.format(self.npad))
        self.fill_outside_npad_edit.editingFinished.connect(self.set_npad)
        
        self.fill_outside_button = QtWidgets.QPushButton('Filter Outside Area')
        self.fill_outside_button.pressed.connect(self.fill_outside_area)

        self.map_figure = Figure()
        self.map_canvas = FigureCanvas(self.map_figure)
        self.map_canvas.mpl_connect('pick_event', self.map_on_pick)
        self.map_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Expanding)
        self.map_toolbar = NavigationToolbar(self.map_canvas, self)

        self.map_depth_label = QtWidgets.QLabel('Depth {0:>10.2f} {1}'.format(0,
                                                                    self.units))
        self.map_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.map_slider.valueChanged.connect(self.set_map_index)
        self.map_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.map_slider.setMinimum(0)
        self.map_slider.setMaximum(0)
        self.map_slider.setTickInterval(1)

        self.east_figure= Figure()
        self.east_canvas = FigureCanvas(self.east_figure)
        self.east_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                       QtWidgets.QSizePolicy.Expanding)
        self.east_toolbar = NavigationToolbar(self.east_canvas, self)

        self.east_label = QtWidgets.QLabel('Easting {0:>10.2f} {1}'.format(0,
                                                                   self.units))
        self.east_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.east_slider.valueChanged.connect(self.set_east_index)
        self.east_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.east_slider.setMinimum(0)
        self.east_slider.setMaximum(0)
        self.east_slider.setTickInterval(1)

        self.north_figure = Figure()
        self.north_canvas = FigureCanvas(self.north_figure)
        self.north_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Expanding)
        self.north_toolbar = NavigationToolbar(self.north_canvas, self)

        self.north_label = QtWidgets.QLabel('Northing {0:>10.2f} m'.format(0,
                                                                    self.units))
        self.north_slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.north_slider.valueChanged.connect(self.set_north_index)
        self.north_slider.setTickPosition(QtWidgets.QSlider.TicksBelow)
        self.north_slider.setMinimum(0)
        self.north_slider.setMaximum(0)
        self.north_slider.setTickInterval(1)


        self.figure_3d = Figure()
        self.canvas_3d = FigureCanvas(self.figure_3d)
        self.canvas_3d.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                       QtWidgets.QSizePolicy.Expanding)

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
        self.cb_ax.set_yticks(np.arange(self.res_limits[0], self.res_limits[1]))
        self.cb_ax.set_yticklabels(['10$^{0}$'.format(ii) for ii in
                                    np.arange(self.res_limits[0], self.res_limits[1])])

        self.res_line, = self.cb_ax.plot([0, 1], 
                                         [np.log10(self.res_value),
                                         np.log10(self.res_value)],
                                         lw=3, 
                                         color='k',
                                         picker=5)
        self.cb_canvas.mpl_connect('button_press_event', self.on_res_pick)
#        self.cb_canvas.mpl_connect('pick_event', self.on_res_pick)
        self.cb_ax.set_xticks([0, 1])
        self.cb_ax.set_xticklabels(['', ''])
        self.cb_ax.axis('tight')
        self.cb_canvas.draw()

        self.cb_line_edit = QtWidgets.QLineEdit()
        self.cb_line_edit.setMaximumWidth(80)
        self.cb_line_edit.setText('{0:.2f}'.format(self.res_value))
        self.cb_line_edit.editingFinished.connect(self.set_res_value)

        self.cb_label = QtWidgets.QLabel('Ohm-m')
        
        ##------------------------------------------------
        ## Layout
        
        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addWidget(self.fill_outside_npad_label)
        button_layout.addWidget(self.fill_outside_npad_edit)
        button_layout.addWidget(self.fill_outside_avg_pad_label)
        button_layout.addWidget(self.fill_outside_avg_pad_edit)
        button_layout.addWidget(self.fill_outside_button)

        map_bottom_layout = QtWidgets.QHBoxLayout()
        map_bottom_layout.addWidget(self.map_depth_label)
        map_bottom_layout.addWidget(self.map_slider)
        map_layout = QtWidgets.QVBoxLayout()
        map_layout.addWidget(self.map_toolbar)
        map_layout.addWidget(self.map_canvas)
        map_layout.addLayout(map_bottom_layout)

        east_bottom_layout = QtWidgets.QHBoxLayout()
        east_bottom_layout.addWidget(self.east_label)
        east_bottom_layout.addWidget(self.east_slider)
        east_layout = QtWidgets.QVBoxLayout()
        east_layout.addWidget(self.east_toolbar)
        east_layout.addWidget(self.east_canvas)
        east_layout.addLayout(east_bottom_layout)

        north_bottom_layout = QtWidgets.QHBoxLayout()
        north_bottom_layout.addWidget(self.north_label)
        north_bottom_layout.addWidget(self.north_slider)
        north_layout = QtWidgets.QVBoxLayout()
        north_layout.addWidget(self.north_toolbar)
        north_layout.addWidget(self.north_canvas)
        north_layout.addLayout(north_bottom_layout)


        grid_layout = QtWidgets.QGridLayout()
        grid_layout.addLayout(map_layout, 1, 1)
        grid_layout.addLayout(east_layout, 1, 2)
        grid_layout.addLayout(north_layout, 2, 1)
        grid_layout.addWidget(self.canvas_3d, 2, 2)

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

        # set slider bar intervals
        # need the minus 1 cause we are using the value of the slider as
        # the index.
        self.map_slider.setMaximum(self.model_obj.grid_z.size-1)
        self.east_slider.setMaximum(self.model_obj.grid_north.size-1)
        self.north_slider.setMaximum(self.model_obj.grid_east.size-1)
        
        self.east_label.setText('{0:.2f}'.format(self.model_obj.grid_east[0]))
        self.north_label.setText('{0:.2f}'.format(self.model_obj.grid_north[0]))

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
        self.map_ax.set_aspect('equal', 'box-forced')
        

        self.north_ax = self.north_figure.add_subplot(1, 1, 1, sharex=self.map_ax,
                                                      aspect='equal')
        self.north_ax.pcolormesh(self.plot_east_z,
                               self.plot_z_east,
                               np.log10(self.model_obj.res_model[self.north_index, :, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        self.north_ax.set_xlabel('Easting {0}'.format(self.units))
        self.north_ax.set_ylabel('Depth {0}'.format(self.units))
        self.north_ax.set_aspect('equal', 'box-forced')
        z_lim = self.north_ax.get_ylim()
        self.north_ax.set_ylim(z_lim[1], z_lim[0])
        

        self.east_ax = self.east_figure.add_subplot(1, 1, 1,
                                                    aspect='equal',
                                                    sharex=self.map_ax,
                                                    sharey=self.north_ax)
        self.east_ax.pcolormesh(self.plot_north_z,
                               self.plot_z_north,
                               np.log10(self.model_obj.res_model[:, self.east_index, :]),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])

        self.east_ax.set_xlabel('Northing {0}'.format(self.units))
        self.east_ax.set_ylabel('Depth {0}'.format(self.units))
        self.east_ax.set_aspect('equal', 'box-forced')
        
        if self.data_fn is not None:
            self.map_ax.scatter(self.data_obj.station_locations['rel_east']/self.scale,
                                self.data_obj.station_locations['rel_north']/self.scale,
                                marker='v',
                                c='k',
                                s=10)
        
        self.north_canvas.draw()
        self.map_canvas.draw()
        self.east_canvas.draw()

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
            y_data = np.log10(np.round(y_data, -int(np.floor(np.log10(y_data)))))
            
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

        self.map_ax.pcolormesh(self.plot_east_map,
                               self.plot_north_map,
                               np.log10(self.model_obj.res_model[:, :, self.map_index].T),
                               cmap=self.cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
        if self.data_fn is not None:
            self.map_ax.scatter(self.data_obj.station_locations['rel_east']/self.scale,
                                self.data_obj.station_locations['rel_north']/self.scale,
                                marker='v',
                                c='k',
                                s=10)
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
                self.model_obj.res_model[yy, xx, self.map_index] = self.res_value
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
                self.model_obj.res_model[xx, self.east_index, yy] = self.res_value

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
                    self.model_obj.res_model[self.north_index, xx, yy] = self.res_value

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
        redraw all the plots after things have changed
        """
        self.map_ax.pcolormesh(self.plot_east_map,
                       self.plot_north_map,
                       np.log10(self.model_obj.res_model[:, :, self.map_index].T),
                       cmap=self.cmap,
                       vmin=self.res_limits[0],
                       vmax=self.res_limits[1])
                       
        if self.data_fn is not None:
            self.map_ax.scatter(self.data_obj.station_locations['rel_east']/self.scale,
                                self.data_obj.station_locations['rel_north']/self.scale,
                                marker='v',
                                c='k',
                                s=10)

        self.east_ax.pcolormesh(self.plot_north_z,
                        self.plot_z_north,
                        np.log10(self.model_obj.res_model[:, self.east_index, :]),
                        cmap=self.cmap,
                        vmin=self.res_limits[0],
                        vmax=self.res_limits[1])

        self.north_ax.pcolormesh(self.plot_east_z,
                         self.plot_z_east,
                         np.log10(self.model_obj.res_model[self.north_index, :, :]),
                         cmap=self.cmap,
                         vmin=self.res_limits[0],
                         vmax=self.res_limits[1])

        self.north_canvas.draw()
        self.east_canvas.draw()
        self.map_canvas.draw()

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


