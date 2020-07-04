# -*- coding: utf-8 -*-
"""
ModEM data and response visualization with a gui.

The user will be able to choose from stations within the data to look at
in either impedance or apparent resistivity and phase.

The functionality is quite simple at the moment

JP 2016
"""
# 
#==============================================================================
# Imports
#==============================================================================
# standard imports
import os
import sys

try:
    from PyQt5 import QtCore, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

import mtpy.imaging.mtplottools as mtplottools
import mtpy.modeling.modem as modem

#==============================================================================
class ModEMPlotResponse(QtWidgets.QMainWindow):
    """
    main window
    """
    
    def __init__(self):
        super(ModEMPlotResponse, self).__init__()
        
        self.data_fn = None
        self.resp_fn = None
        self.modem_data = None
        
        self.setup_ui()
        
    def setup_ui(self):
        self.setWindowTitle("Plot ModEM Responses")
        self.setWindowState(QtCore.Qt.WindowMaximized)
        screen_shape = QtWidgets.QDesktopWidget().screenGeometry()         
        
        #create a menu bar on the window with 4 different items
        self.menubar = QtWidgets.QMenuBar(self)
        self.menubar.setGeometry(QtCore.QRect(0, 0, screen_shape.width(), 38))

        # add a tab for File --> open, close, save
        self.menu_data_file = QtWidgets.QMenu(self.menubar)
        self.menu_data_file.setTitle("Data File")
        
        self.menu_resp_file = QtWidgets.QMenu(self.menubar)
        self.menu_resp_file.setTitle("Response File")
        
        # add a tab for chaning the display
        self.menu_display = QtWidgets.QMenu(self.menubar)
        self.menu_display.setTitle("Display")
        
        # add a tab for help
        self.menu_help = QtWidgets.QMenu(self.menubar)
        self.menu_help.setTitle("Help")

        self.setMenuBar(self.menubar)
     
        # set the actions for the data file menu item 
        # set an open option that on click opens a modem file
        self.action_open_data = QtWidgets.QAction(self)
        self.action_open_data.setText("&Open")
        self.action_open_data.setShortcut("Ctrl+o")
        self.action_open_data.triggered.connect(self.get_data_file)

        # set a close that closes the main window
        self.action_close = QtWidgets.QAction(self)
        self.action_close.setText("Close")
        self.action_close.setShortcut("Ctrl+x")
        self.action_close.triggered.connect(self.close)

        # set a save option that will eventually save the masked data
        self.action_save_data = QtWidgets.QAction(self)
        self.action_save_data.setText("&Save Edits")
        self.action_save_data.setShortcut("Ctrl+s")
        self.action_save_data.triggered.connect(self.save_edits)

        # add the action on the menu tab
        self.menu_data_file.addAction(self.action_open_data)
        self.menu_data_file.addAction(self.action_close)
        self.menu_data_file.addAction(self.action_save_data)
        self.menubar.addAction(self.menu_data_file.menuAction())
        
        # set the action items for the response file
        self.action_resp_open = QtWidgets.QAction(self)
        self.action_resp_open.setText("Open")
        self.action_resp_open.triggered.connect(self.get_resp_fn)
        self.menu_resp_file.addAction(self.action_resp_open)
        self.menubar.addAction(self.menu_resp_file.menuAction())
#        
        #adding options for display plot type        
        self.menu_plot_type = QtWidgets.QMenu(self)
        self.menu_plot_type.setTitle("Plot Type")
        self.menu_display.addMenu(self.menu_plot_type)
        self.menubar.addAction(self.menu_display.menuAction())
        
        #set plot impedance or resistivity and phase
        self.action_plot_z = QtWidgets.QAction(self)
        self.action_plot_z.setText('Impedance')
        self.action_plot_z.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_z)
        self.action_plot_z.toggled.connect(self.status_checked_ptz)
        
        self.action_plot_rp = QtWidgets.QAction(self)
        self.action_plot_rp.setText('Resistivity-Phase')
        self.action_plot_rp.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_rp)
        self.action_plot_rp.toggled.connect(self.status_checked_ptrp)
        
        self.action_plot_settings = QtWidgets.QAction(self)
        self.action_plot_settings.setText('Settings')
        self.action_plot_settings.triggered.connect(self.show_settings)
        self.menu_display.addAction(self.action_plot_settings)
        self.menubar.addAction(self.menu_display.menuAction())

        self.menu_display.addAction(self.menu_plot_type.menuAction())
        
        self.action_help = QtWidgets.QAction(self)
        self.action_help.setText('Help Documentation')
        self.action_help.triggered.connect(self.disp_help)
        self.menu_help.addAction(self.action_help)
        self.menubar.addAction(self.menu_help.menuAction())
        
        self.plot_response = PlotResponses(self.data_fn, self.resp_fn)
        self.setCentralWidget(self.plot_response)
        
    
        #self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(self)
        
    def status_checked_ptz(self, toggled):
        """
        be sure that only one plot style is checked
        """
        
        self.plot_response.plot_z = toggled
        if toggled == True:
            untoggled = False
            
        elif toggled == False:
            untoggled = True
        
        self.action_plot_z.setChecked(toggled)
        self.action_plot_rp.setChecked(untoggled)
        
    def status_checked_ptrp(self, toggled):
        """
        be sure that only one plot style is checked
        """
        
        if toggled == True:
            untoggled = False
            self.plot_response.plot_z = False
        elif toggled == False:
            untoggled = True
            self.plot_response.plot_z = True
        
        self.action_plot_z.setChecked(untoggled)
        self.action_plot_rp.setChecked(toggled)
                        
    def get_data_file(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                           filter='(*.dat);; (*.data)')[0])
        
                                       
        self.plot_response.data_fn = os.path.abspath(fn)
        self.dir_path = os.path.dirname(fn)
                                       
        
    def save_edits(self):
        """
        save edits to another file
        """
        fn_dialog = QtWidgets.QFileDialog()
        save_fn = str(fn_dialog.getSaveFileName(caption='Choose File to save',
                                                filter='*.dat')[0])
        
        self.plot_response.modem_data.write_data_file(save_path=os.path.dirname(save_fn),
                                                      fn_basename=os.path.basename(save_fn),
                                                      compute_error=False,
                                                      fill=False,
                                                      elevation=True)
        
    def get_resp_fn(self):
        """
        get response file name
        """
        
        fn_dialog = QtWidgets.QFileDialog(directory=self.dir_path)
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM response file',
                                       filter='(*.dat);; (*.data)')[0])
                                       
        self.plot_response.resp_fn = os.path.abspath(fn)
        
    def show_settings(self):
        self.settings_window = PlotSettings(**self.__dict__)
        self.settings_window.show()
        self.settings_window.settings_updated.connect(self.update_settings)
        
    def update_settings(self):
        
        for attr in sorted(self.settings_window.__dict__.keys()):
            setattr(self, attr, self.settings_window.__dict__[attr])
            
        self.plot()
        
    def disp_help(self):
        """
        display a help dialogue
        """
        ll = ['This GUI will allow you to edit your data by masking points',
              'and adding error bars to dodgy data points.  Only the top row',
              'row is editable for now. However, all edits to the top row ',
              'are applied to the bottom row (real and imaginary parts).\n',
              '   * Left-Click the mouse to mask a point this will mask both',
              '     the real and imaginary part of that component.\n',
              '   * Right-Click the mouse to add error bars to the data point',
              '     again this will apply to both real and imaginary parts of',
              '     the selected component. Current it goes up by 5%\n',
              '   * To save your masking, go to Data File -> Save Edits' ]
        
        help_string = '\n'.join(ll)        
        
        help_popup = QtWidgets.QMessageBox.information(self.central_widget, 'Help', 
                                                   help_string)
        
        
                                  
#==============================================================================
# plot part
#==============================================================================
class PlotResponses(QtWidgets.QWidget):
    """
    the plot and list of stations
    """
    
    def __init__(self, data_fn=None, resp_fn=None):
        super(PlotResponses, self).__init__()
        
        self.file_watcher_dfn = QtCore.QFileSystemWatcher()
        self.file_watcher_dfn.fileChanged.connect(self.file_changed_dfn)
    
        
        self.modem_data = None
        self.modem_resp = None
        
        self.station = None
        
        self._modem_data_copy = None
        
        self._plot_z = False
        self.plot_settings = PlotSettings()
        
        self._ax = None
        self._ax2 = None
        self._key = 'z'
        self._ax_index = 0
        self.ax_list = None
        
        self.setup_ui()
        
        self._data_fn = data_fn
        self._resp_fn = resp_fn
    
    #------------------------------------------------
    # make the data_fn and resp_fn properties so that if they are reset
    # they will read in the data to a new modem.Data object
    # trying to use decorators for syntactical sugar    
    @property
    def data_fn(self):
        return self._data_fn
        
    @data_fn.setter
    def data_fn(self, data_fn):
        self._data_fn = os.path.abspath(data_fn)
        self.file_watcher_dfn.addPath(self._data_fn)
        
        # create new modem data object
        self.modem_data = modem.Data()
        self.modem_data.read_data_file(self._data_fn)
        
        # make a back up copy that will be unchanged
        # that way we can revert back
        self._modem_data_copy = modem.Data()
        self._modem_data_copy.read_data_file(self._data_fn)
        
        self.dirpath = os.path.dirname(self._data_fn)
        
        # fill list of stations
        station_list = np.array(sorted(self.modem_data.mt_dict.keys()))
        self.list_widget.clear()
        for station in station_list:
            self.list_widget.addItem(station)
            
        if self.station is None:
            self.station = station_list[0]
            
        self.plot()
        
    @property
    def resp_fn(self):
        return self._resp_fn
        
    @resp_fn.setter
    def resp_fn(self, resp_fn):
        self._resp_fn = os.path.abspath(resp_fn)
        self.modem_resp = modem.Data()

        self.modem_resp.read_data_file(self._resp_fn)
        self.plot() 
        
    @property
    def plot_z(self):
        return self._plot_z
        
    @plot_z.setter
    def plot_z(self, value):
        self._plot_z = value
        self.plot()
        
    #----------------------------
    def setup_ui(self):
        """
        setup the user interface with list of stations on the left and the 
        plot on the right.  There will be a button for save edits.
        """
        
        #make a widget that will be the station list
        self.list_widget = QtWidgets.QListWidget()
        self.list_widget.itemClicked.connect(self.get_station)
        self.list_widget.currentItemChanged.connect(self.get_station)
        self.list_widget.setMaximumWidth(150)
        
        self.save_edits_button = QtWidgets.QPushButton()
        self.save_edits_button.setText("Save Edits")
        self.save_edits_button.setStyleSheet("background-color: #42f489")
        self.save_edits_button.pressed.connect(self.save_edits)
        
        self.apply_edits_button = QtWidgets.QPushButton()
        self.apply_edits_button.setText('Apply Edits')
        self.apply_edits_button.setStyleSheet("background-color: #c6dcff")
        self.apply_edits_button.pressed.connect(self.apply_edits)
        
#        self.undo_edit_button = QtWidgets.QPushButton()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        self.mpl_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.mpl_widget.setFocus()
        
        # be able to edit the data
        self.mpl_widget.mpl_connect('pick_event', self.on_pick)
        self.mpl_widget.mpl_connect('axes_enter_event', self.in_axes)
        
        #make sure the figure takes up the entire plottable space
        self.mpl_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                     QtWidgets.QSizePolicy.Expanding)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, self)
         
        # set the layout for the plot
        mpl_vbox = QtWidgets.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)
        
        left_layout = QtWidgets.QVBoxLayout()
        left_layout.addWidget(self.list_widget)
        left_layout.addWidget(self.apply_edits_button)
        left_layout.addWidget(self.save_edits_button)
        
        # set the layout the main window
        layout = QtWidgets.QHBoxLayout()
        layout.addLayout(left_layout)
        layout.addLayout(mpl_vbox)
        
        self.setLayout(layout)
    
    def get_station(self, widget_item):
        """
        get the station name from the clicked station 
        """
        try:
            self.station = str(widget_item.text()) 
        except AttributeError:
            print('Error: Widget is None')
            
        self.plot()
        
    def file_changed_dfn(self):
        """
        data file changed outside the program reload it
        """
        
        print('{0} changed'.format(self.data_fn))
        self.data_fn = os.path.abspath(self._data_fn)
        
    def save_edits(self):
        """
        save edits to another file
        """
        fn_dialog = QtWidgets.QFileDialog()
        save_fn = str(fn_dialog.getSaveFileName(caption='Choose File to save',
                                                filter='*.dat')[0])
        
        self.modem_data.write_data_file(save_path=os.path.dirname(save_fn),
                                        fn_basename=os.path.basename(save_fn),
                                        compute_error=False,
                                        fill=False,
                                        elevation=True)
        
    def apply_edits(self):
        self.plot()
        
    def plot(self):
        """
        plot the data
        """

        if self.station is None:
            return
            
        z_obj = self.modem_data.mt_dict[self.station].Z
        t_obj = self.modem_data.mt_dict[self.station].Tipper
        period = self.modem_data.period_list
        
        # need to make sure that resistivity and phase is computed
        z_obj.compute_resistivity_phase()

        plt.rcParams['font.size'] = self.plot_settings.fs
        fontdict = {'size':self.plot_settings.fs+2, 'weight':'bold'} 

        #--> make key word dictionaries for plotting
        kw_xx = {'color':self.plot_settings.cted,
                 'marker':self.plot_settings.mted,
                 'ms':self.plot_settings.ms,
                 'ls':':',
                 'lw':self.plot_settings.lw,
                 'e_capsize':self.plot_settings.e_capsize,
                 'e_capthick':self.plot_settings.e_capthick,
                 'picker':3}        
       
        kw_yy = {'color':self.plot_settings.ctmd,
                 'marker':self.plot_settings.mtmd,
                 'ms':self.plot_settings.ms,
                 'ls':':',
                 'lw':self.plot_settings.lw,
                 'e_capsize':self.plot_settings.e_capsize,
                 'e_capthick':self.plot_settings.e_capthick,
                 'picker':3} 

        #convert to apparent resistivity and phase
        if self.plot_z == True:
            scaling = np.zeros_like(z_obj.z)
            for ii in range(2):
                for jj in range(2):
                    scaling[:, ii, jj] = 1./np.sqrt(z_obj.freq)
            plot_res = abs(z_obj.z.real*scaling)
            plot_res_err = abs(z_obj.z_err*scaling)
            plot_phase = abs(z_obj.z.imag*scaling)
            plot_phase_err = abs(z_obj.z_err*scaling)
            h_ratio = [1, 1, .5]
            
        elif self.plot_z == False:
            plot_res = z_obj.resistivity
            plot_res_err = z_obj.resistivity_err
            plot_phase = z_obj.phase
            plot_phase_err = z_obj.phase_err
            h_ratio = [1.5, 1, .5]
        
        #find locations where points have been masked
        nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
        nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
        nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
        nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
        ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
        nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]
        

        self.figure.clf()
        self.figure.suptitle(str(self.station), fontdict=fontdict)
        
        #set the grid of subplots
        if np.all(t_obj.tipper == 0.0) == True:
            self.plot_tipper = False
        else:
            self.plot_tipper = True
            
        gs = gridspec.GridSpec(3, 4, height_ratios=h_ratio)
        gs.update(wspace=self.plot_settings.subplot_wspace,
                   left=self.plot_settings.subplot_left,
                   top=self.plot_settings.subplot_top,
                   bottom=self.plot_settings.subplot_bottom, 
                   right=self.plot_settings.subplot_right, 
                   hspace=self.plot_settings.subplot_hspace)

        axrxx = self.figure.add_subplot(gs[0, 0])
        axrxy = self.figure.add_subplot(gs[0, 1], sharex=axrxx)
        axryx = self.figure.add_subplot(gs[0, 2], sharex=axrxx)
        axryy = self.figure.add_subplot(gs[0, 3], sharex=axrxx)
        
        axpxx = self.figure.add_subplot(gs[1, 0])
        axpxy = self.figure.add_subplot(gs[1, 1], sharex=axrxx)
        axpyx = self.figure.add_subplot(gs[1, 2], sharex=axrxx)
        axpyy = self.figure.add_subplot(gs[1, 3], sharex=axrxx)
        
        axtxr = self.figure.add_subplot(gs[2, 0], sharex=axrxx)
        axtxi = self.figure.add_subplot(gs[2, 1], sharex=axrxx)
        axtyr = self.figure.add_subplot(gs[2, 2], sharex=axrxx)
        axtyi = self.figure.add_subplot(gs[2, 3], sharex=axrxx)
        
        self.ax_list = [axrxx, axrxy, axryx, axryy,
                        axpxx, axpxy, axpyx, axpyy,
                        axtxr, axtxi, axtyr, axtyi]
        

        # plot data response
        erxx = mtplottools.plot_errorbar(axrxx, 
                                         period[nzxx], 
                                         plot_res[nzxx, 0, 0], 
                                         plot_res_err[nzxx, 0, 0],
                                         **kw_xx)
        erxy = mtplottools.plot_errorbar(axrxy, 
                                         period[nzxy], 
                                         plot_res[nzxy, 0, 1], 
                                         plot_res_err[nzxy, 0, 1],
                                         **kw_xx)
        eryx = mtplottools.plot_errorbar(axryx, 
                                         period[nzyx], 
                                         plot_res[nzyx, 1, 0], 
                                         plot_res_err[nzyx, 1, 0],
                                         **kw_yy)
        eryy = mtplottools.plot_errorbar(axryy, 
                                         period[nzyy], 
                                         plot_res[nzyy, 1, 1], 
                                         plot_res_err[nzyy, 1, 1],
                                         **kw_yy)
        #plot phase                         
        epxx = mtplottools.plot_errorbar(axpxx, 
                                         period[nzxx], 
                                         plot_phase[nzxx, 0, 0], 
                                         plot_phase_err[nzxx, 0, 0],
                                         **kw_xx)
        epxy = mtplottools.plot_errorbar(axpxy, 
                                         period[nzxy], 
                                         plot_phase[nzxy, 0, 1], 
                                         plot_phase_err[nzxy, 0, 1],
                                         **kw_xx)
        epyx = mtplottools.plot_errorbar(axpyx, 
                                         period[nzyx], 
                                         plot_phase[nzyx, 1, 0], 
                                         plot_phase_err[nzyx, 1, 0],
                                         **kw_yy)
        epyy = mtplottools.plot_errorbar(axpyy, 
                                         period[nzyy], 
                                         plot_phase[nzyy, 1, 1], 
                                         plot_phase_err[nzyy, 1, 1],
                                         **kw_yy)
                                      
        #plot tipper
        if self.plot_tipper == True:
            ertx = mtplottools.plot_errorbar(axtxr, 
                                             period[ntx],
                                             t_obj.tipper[ntx, 0, 0].real,
                                             t_obj.tipper_err[ntx, 0, 0],
                                             **kw_xx)
            erty = mtplottools.plot_errorbar(axtyr, 
                                             period[nty],
                                             t_obj.tipper[nty, 0, 1].real,
                                             t_obj.tipper_err[nty, 0, 1],
                                             **kw_yy)
                                     
            eptx = mtplottools.plot_errorbar(axtxi, 
                                             period[ntx],
                                             t_obj.tipper[ntx, 0, 0].imag,
                                             t_obj.tipper_err[ntx, 0, 0],
                                             **kw_xx)
            epty = mtplottools.plot_errorbar(axtyi, 
                                             period[nty],
                                             t_obj.tipper[nty, 0, 1].imag,
                                             t_obj.tipper_err[nty, 0, 1],
                                             **kw_yy)
        
        #----------------------------------------------
        # get error bar list for editing later        
        if self.plot_tipper == False: 
            try:                   
                self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                  [erxy[1][0], erxy[1][1], erxy[2][0]],
                                  [eryx[1][0], eryx[1][1], eryx[2][0]],
                                  [eryy[1][0], eryy[1][1], eryy[2][0]], 
                                  [epxx[1][0], epxx[1][1], epxx[2][0]],
                                  [epxy[1][0], epxy[1][1], epxy[2][0]],
                                  [epyx[1][0], epyx[1][1], epyx[2][0]],
                                  [epyy[1][0], epyy[1][1], epyy[2][0]]]
                line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
            except IndexError:
                print('Found no Z components for {0}'.format(self.station))
                line_list = [[None], [None], 
                             [None], [None]]
                                     
                self._err_list = [[None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None]]

        else:
            try:                    
                line_list = [[erxx[0]], [erxy[0]], 
                             [eryx[0]], [eryy[0]],
                             [ertx[0]], [erty[0]]]
                                     
                self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                  [erxy[1][0], erxy[1][1], erxy[2][0]],
                                  [eryx[1][0], eryx[1][1], eryx[2][0]],
                                  [eryy[1][0], eryy[1][1], eryy[2][0]],
                                  [epxx[1][0], epxx[1][1], epxx[2][0]],
                                  [epxy[1][0], epxy[1][1], epxy[2][0]],
                                  [epyx[1][0], epyx[1][1], epyx[2][0]],
                                  [epyy[1][0], epyy[1][1], epyy[2][0]],
                                  [ertx[1][0], ertx[1][1], ertx[2][0]],
                                  [eptx[1][0], eptx[1][1], eptx[2][0]],
                                  [erty[1][0], erty[1][1], erty[2][0]],
                                  [epty[1][0], epty[1][1], epty[2][0]]]
            except IndexError:
                print('Found no Z components for {0}'.format(self.station))
                line_list = [[None], [None], 
                             [None], [None],
                             [None], [None]]
                                     
                self._err_list = [[None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None],
                                  [None, None, None]]
        
        #------------------------------------------
        # make things look nice        
        # set titles of the Z components
        label_list = [['$Z_{xx}$'], ['$Z_{xy}$'], 
                       ['$Z_{yx}$'], ['$Z_{yy}$']] 
        for ax, label in zip(self.ax_list[0:4], label_list):
            ax.set_title(label[0],fontdict={'size':self.plot_settings.fs+2, 
                                          'weight':'bold'}) 
                                          
        # set legends for tipper components
        # fake a line
        l1 = plt.Line2D([0], [0], linewidth=0, color='w', linestyle='None', 
                        marker='.')
        t_label_list = ['Re{$T_x$}', 'Im{$T_x$}', 'Re{$T_y$}', 'Im{$T_y$}']
        label_list += [['$T_{x}$'], ['$T_{y}$']]
        for ax, label in zip(self.ax_list[-4:], t_label_list):
            ax.legend([l1], [label], loc='upper left',
                      markerscale=.01,
                      borderaxespad=.05,
                      labelspacing=.01,
                      handletextpad=.05,
                      borderpad=.05,
                      prop={'size':max([self.plot_settings.fs, 5])})
        
        #--> set limits if input
        if self.plot_settings.res_xx_limits is not None:
            axrxx.set_ylim(self.plot_settings.res_xx_limits)    
        if self.plot_settings.res_xy_limits is not None:
            axrxy.set_ylim(self.plot_settings.res_xy_limits)    
        if self.plot_settings.res_yx_limits is not None:
            axryx.set_ylim(self.plot_settings.res_yx_limits)    
        if self.plot_settings.res_yy_limits is not None:
            axryy.set_ylim(self.plot_settings.res_yy_limits) 
            
        if self.plot_settings.phase_xx_limits is not None:
            axpxx.set_ylim(self.plot_settings.phase_xx_limits)    
        if self.plot_settings.phase_xy_limits is not None:
            axpxy.set_ylim(self.plot_settings.phase_xy_limits)    
        if self.plot_settings.phase_yx_limits is not None:
            axpyx.set_ylim(self.plot_settings.phase_yx_limits)    
        if self.plot_settings.phase_yy_limits is not None:
            axpyy.set_ylim(self.plot_settings.phase_yy_limits) 
    
        #set axis properties
        for aa, ax in enumerate(self.ax_list):
            ax.tick_params(axis='y', pad=self.plot_settings.ylabel_pad)
            ylabels = ax.get_yticks().tolist()
            if aa < 8:
                ylabels[-1] = ''
                ylabels[0] = ''
                ax.set_yticklabels(ylabels)
                plt.setp(ax.get_xticklabels(), visible=False)
                if self.plot_z == True:
                    ax.set_yscale('log', nonposy='clip')

            else:
                ax.set_xlabel('Period (s)', fontdict=fontdict)
                
            if aa < 4 and self.plot_z is False:
                ax.set_yscale('log', nonposy='clip')
                    
            #set axes labels
            if aa == 0:
                if self.plot_z == False:
                    ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                  fontdict=fontdict)
                elif self.plot_z == True:
                    ax.set_ylabel('Re[Z (mV/km nT)]',
                                  fontdict=fontdict)
            elif aa == 4:
                if self.plot_z == False:
                    ax.set_ylabel('Phase (deg)',
                                  fontdict=fontdict)
                elif self.plot_z == True:
                    ax.set_ylabel('Im[Z (mV/km nT)]',
                                  fontdict=fontdict)
            elif aa == 8:
                ax.set_ylabel('Tipper',
                              fontdict=fontdict)
                    
            if aa > 7:
                if self.plot_settings.tipper_limits is not None:
                    ax.set_ylim(self.plot_settings.tipper_limits)
                else:
                    pass

            ax.set_xscale('log', nonposx='clip')
            ax.set_xlim(xmin=10**(np.floor(np.log10(period[0])))*1.01,
                     xmax=10**(np.ceil(np.log10(period[-1])))*.99)
            ax.grid(True, alpha=.25)
            
        ##----------------------------------------------
        #plot model response
        if self.modem_resp is not None:
            try:
                resp_z_obj = self.modem_resp.mt_dict[self.station].Z
                resp_z_err = np.nan_to_num((z_obj.z-resp_z_obj.z)/z_obj.z_err)
                resp_z_obj.compute_resistivity_phase()
                
                resp_t_obj = self.modem_resp.mt_dict[self.station].Tipper
                resp_t_err = np.nan_to_num((t_obj.tipper-resp_t_obj.tipper)/t_obj.tipper_err)
            except KeyError:
                print('Could not find {0} in .resp file'.format(self.station))
                self.mpl_widget.draw()
                return
            
            try:
                resp_period = 1./resp_z_obj.freq
            except TypeError:
                resp_period = 1./resp_t_obj.freq
                
            #convert to apparent resistivity and phase
            if self.plot_z == True:
                scaling = np.zeros_like(resp_z_obj.z)
                for ii in range(2):
                    for jj in range(2):
                        scaling[:, ii, jj] = 1./np.sqrt(resp_z_obj.freq)
                r_plot_res = abs(resp_z_obj.z.real*scaling)
                r_plot_phase = abs(resp_z_obj.z.imag*scaling)
                
            elif self.plot_z == False:
                r_plot_res = resp_z_obj.resistivity
                r_plot_phase = resp_z_obj.phase
                
            #find locations where points have been masked
            nzxx_r = np.nonzero(resp_z_obj.z[:, 0, 0])[0]
            nzxy_r = np.nonzero(resp_z_obj.z[:, 0, 1])[0]
            nzyx_r = np.nonzero(resp_z_obj.z[:, 1, 0])[0]
            nzyy_r = np.nonzero(resp_z_obj.z[:, 1, 1])[0]
            ntx_r = np.nonzero(resp_t_obj.tipper[:, 0, 0])[0]
            nty_r = np.nonzero(resp_t_obj.tipper[:, 0, 1])[0]

            rms_xx = resp_z_err[nzxx_r, 0, 0].std()
            rms_xy = resp_z_err[nzxy_r, 0, 1].std()
            rms_yx = resp_z_err[nzyx_r, 1, 0].std()
            rms_yy = resp_z_err[nzyy_r, 1, 1].std()
            
            #--> make key word dictionaries for plotting
            kw_xx = {'color':self.plot_settings.ctem,
                     'marker':self.plot_settings.mtem,
                     'ms':self.plot_settings.ms,
                     'ls':':',
                     'lw':self.plot_settings.lw,
                     'e_capsize':self.plot_settings.e_capsize,
                     'e_capthick':self.plot_settings.e_capthick}        
           
            kw_yy = {'color':self.plot_settings.ctmm,
                     'marker':self.plot_settings.mtmm,
                     'ms':self.plot_settings.ms,
                     'ls':':',
                     'lw':self.plot_settings.lw,
                     'e_capsize':self.plot_settings.e_capsize,
                     'e_capthick':self.plot_settings.e_capthick}
            
            # plot data response
            rerxx = mtplottools.plot_errorbar(axrxx, 
                                             resp_period[nzxx_r], 
                                             r_plot_res[nzxx_r, 0, 0], 
                                             None,
                                             **kw_xx)
            rerxy = mtplottools.plot_errorbar(axrxy, 
                                             resp_period[nzxy_r], 
                                             r_plot_res[nzxy_r, 0, 1], 
                                             None,
                                             **kw_xx)
            reryx = mtplottools.plot_errorbar(axryx, 
                                             resp_period[nzyx_r], 
                                             r_plot_res[nzyx_r, 1, 0], 
                                             None,
                                             **kw_yy)
            reryy = mtplottools.plot_errorbar(axryy, 
                                             resp_period[nzyy_r], 
                                             r_plot_res[nzyy_r, 1, 1], 
                                             None,
                                             **kw_yy)
            #plot phase                         
            repxx = mtplottools.plot_errorbar(axpxx, 
                                             resp_period[nzxx_r], 
                                             r_plot_phase[nzxx_r, 0, 0], 
                                             None,
                                             **kw_xx)
            repxy = mtplottools.plot_errorbar(axpxy, 
                                             resp_period[nzxy_r], 
                                             r_plot_phase[nzxy_r, 0, 1], 
                                             None,
                                             **kw_xx)
            repyx = mtplottools.plot_errorbar(axpyx, 
                                             resp_period[nzyx_r], 
                                             r_plot_phase[nzyx_r, 1, 0], 
                                             None,
                                             **kw_yy)
            repyy = mtplottools.plot_errorbar(axpyy, 
                                             resp_period[nzyy_r], 
                                             r_plot_phase[nzyy_r, 1, 1], 
                                             None,
                                             **kw_yy)
                                          
            #plot tipper
            if self.plot_tipper == True:
                rertx = mtplottools.plot_errorbar(axtxr, 
                                                 resp_period[ntx_r],
                                                 resp_t_obj.tipper[ntx_r, 0, 0].real,
                                                 None,
                                                 **kw_xx)
                rerty = mtplottools.plot_errorbar(axtyr, 
                                                 resp_period[nty_r],
                                                 resp_t_obj.tipper[nty_r, 0, 1].real,
                                                 None,
                                                 **kw_yy)
                                         
                reptx = mtplottools.plot_errorbar(axtxi, 
                                                 resp_period[ntx_r],
                                                 resp_t_obj.tipper[ntx_r, 0, 0].imag,
                                                 None,
                                                 **kw_xx)
                repty = mtplottools.plot_errorbar(axtyi, 
                                                 resp_period[nty_r],
                                                 resp_t_obj.tipper[nty_r, 0, 1].imag,
                                                 None,
                                                 **kw_yy)
                             
            if self.plot_tipper == False:
                line_list[0] += [rerxx[0]]
                line_list[1] += [rerxy[0]]
                line_list[2] += [reryx[0]]
                line_list[3] += [reryy[0]]
                label_list[0] += ['$Z^m_{xx}$ '+
                                   'rms={0:.2f}'.format(rms_xx)]
                label_list[1] += ['$Z^m_{xy}$ '+
                               'rms={0:.2f}'.format(rms_xy)]
                label_list[2] += ['$Z^m_{yx}$ '+
                               'rms={0:.2f}'.format(rms_yx)]
                label_list[3] += ['$Z^m_{yy}$ '+
                               'rms={0:.2f}'.format(rms_yy)]
            else:
                line_list[0] += [rerxx[0]]
                line_list[1] += [rerxy[0]]
                line_list[2] += [reryx[0]]
                line_list[3] += [reryy[0]]
                line_list[4] += [rertx[0]]
                line_list[5] += [rerty[0]]
                label_list[0] += ['$Z^m_{xx}$ '+
                                   'rms={0:.2f}'.format(rms_xx)]
                label_list[1] += ['$Z^m_{xy}$ '+
                               'rms={0:.2f}'.format(rms_xy)]
                label_list[2] += ['$Z^m_{yx}$ '+
                               'rms={0:.2f}'.format(rms_yx)]
                label_list[3] += ['$Z^m_{yy}$ '+
                               'rms={0:.2f}'.format(rms_yy)]
                label_list[4] += ['$T^m_{x}$ '+
                                'rms={0:.2f}'.format(resp_t_err[ntx, 0, 0].std())]
                label_list[5] += ['$T^m_{y}$'+
                                'rms={0:.2f}'.format(resp_t_err[nty, 0, 1].std())]
            
            legend_ax_list = self.ax_list[0:4]
            if self.plot_tipper == True:
                legend_ax_list += [self.ax_list[-4], self.ax_list[-2]]
                
            for aa, ax in enumerate(legend_ax_list):
                ax.legend(line_list[aa],
                          label_list[aa],
                          loc=self.plot_settings.legend_loc,
                          bbox_to_anchor=self.plot_settings.legend_pos,
                          markerscale=self.plot_settings.legend_marker_scale,
                          borderaxespad=self.plot_settings.legend_border_axes_pad,
                          labelspacing=self.plot_settings.legend_label_spacing,
                          handletextpad=self.plot_settings.legend_handle_text_pad,
                          borderpad=self.plot_settings.legend_border_pad,
                          prop={'size':max([self.plot_settings.fs, 5])})
        
        self.mpl_widget.draw()
        
    def on_pick(self, event):
        """
        mask a data point when it is clicked on.  
        """         
        data_point = event.artist
        data_period = data_point.get_xdata()[event.ind]
        data_value = data_point.get_ydata()[event.ind]
        
        # get the indicies where the data point has been edited
        try:
            p_index = np.where(self.modem_data.period_list==data_period)[0][0]
            s_index = np.where(self.modem_data.data_array['station']==self.station)[0][0]
        except IndexError:
            return

        if self._key == 'tip':
            data_value_2 = self.modem_data.mt_dict[self.station].Tipper.tipper[p_index,
                                self._comp_index_x, self._comp_index_y]
            if self._ax_index%2 == 0:
                data_value_2 = data_value_2.imag
            else:
                data_value_2 = data_value_2.real
                
        elif self._key == 'z':
            if self.plot_z == True:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.z[p_index,
                                self._comp_index_x, self._comp_index_y]
                                
                if self._ax_index % 2 == 0:
                    data_value_2 = data_value_2.imag
                else:
                    data_value_2 = data_value_2.real
            elif self.plot_z == False and self._ax_index < 4:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.phase[p_index,
                                self._comp_index_x, self._comp_index_y]
            elif self.plot_z == False and self._ax_index >= 4:
                data_value_2 = self.modem_data.mt_dict[self.station].Z.resistivity[p_index,
                                self._comp_index_x, self._comp_index_y]
                        
        if event.mouseevent.button == 1:
            # mask the point in the data mt_dict

            self.modem_data.data_array[s_index][self._key][p_index, 
                        self._comp_index_x, self._comp_index_y] = 0+0j
            
            if self._key == 'tip':
                self.modem_data.mt_dict[self.station].Tipper.tipper[p_index, 
                            self._comp_index_x, self._comp_index_y] = 0+0j
            elif self._key == 'z':
                self.modem_data.mt_dict[self.station].Z.z[p_index, 
                            self._comp_index_x, self._comp_index_y] = 0+0j
            
            # plot the points as masked
            self._ax.plot(data_period, data_value, color=(0, 0, 0),
                          marker='x', 
                          ms=self.plot_settings.ms*2,
                          mew=4)
                          
            self._ax2.plot(data_period, data_value_2, color=(0, 0, 0),
                          marker='x', 
                          ms=self.plot_settings.ms*2,
                          mew=4)
                          
            self._ax.figure.canvas.draw()
            self._ax2.figure.canvas.draw()
        
        # Increase error bars
        if event.mouseevent.button == 3:
            # make sure just checking the top plots 
            
            #put the new error into the error array
            if self._key == 'tip':
                err = self.modem_data.mt_dict[self.station].Tipper.tipper_err[p_index, 
                        self._comp_index_x, self._comp_index_y]
                err = err+abs(err)*self.plot_settings.t_err_increase
                self.modem_data.mt_dict[self.station].Tipper.tipper_err[p_index, 
                        self._comp_index_x, self._comp_index_y] = err
                
                    
            if self._key == 'z':
                err = self.modem_data.mt_dict[self.station].Z.z_err[p_index, 
                        self._comp_index_x, self._comp_index_y]
                err = err+abs(err)*self.plot_settings.z_err_increase

                self.modem_data.mt_dict[self.station].Z.z_err[p_index, 
                    self._comp_index_x, self._comp_index_y] = err
                    
            self.modem_data.data_array[s_index][self._key+'_err'][p_index, 
                        self._comp_index_x, self._comp_index_y] = err

            # make error bar array
            try:
                e_index = event.ind[0]
                eb = self._err_list[self._ax_index][2].get_paths()[e_index].vertices
            except IndexError:
                return
            
            # make ecap array
            ecap_l = self._err_list[self._ax_index][0].get_data()[1][e_index]
            ecap_u = self._err_list[self._ax_index][1].get_data()[1][e_index]
            
            # change apparent resistivity error
            if self._key == 'tip':
                neb_u = eb[0, 1]-self.plot_settings.t_err_increase*abs(eb[0,1])
                neb_l = eb[1, 1]+self.plot_settings.t_err_increase*abs(eb[1,1])
                ecap_l = ecap_l-self.plot_settings.t_err_increase*abs(ecap_l)
                ecap_u = ecap_u+self.plot_settings.t_err_increase*abs(ecap_u)
            elif self._key == 'z':
                neb_u = eb[0, 1]-self.plot_settings.z_err_increase*abs(eb[0,1])
                neb_l = eb[1, 1]+self.plot_settings.z_err_increase*abs(eb[1,1])
                ecap_l = ecap_l-self.plot_settings.z_err_increase*abs(ecap_l)
                ecap_u = ecap_u+self.plot_settings.z_err_increase*abs(ecap_u)
                
            #set the new error bar values
            eb[0,1] = neb_u
            eb[1,1] = neb_l
            
            #reset the error bars and caps
            ncap_l = self._err_list[self._ax_index][0].get_data()
            ncap_u = self._err_list[self._ax_index][1].get_data()
            ncap_l[1][e_index] = ecap_l
            ncap_u[1][e_index] = ecap_u
            
            #set the values 
            self._err_list[self._ax_index][0].set_data(ncap_l)
            self._err_list[self._ax_index][1].set_data(ncap_u)
            self._err_list[self._ax_index][2].get_paths()[e_index].vertices = eb
                                       
            # need to redraw the figure
            self._ax.figure.canvas.draw()
                          
    def in_axes(self, event):
        """
        figure out which axes you just chose the point from
        """

        ax_index_dict = {0:(0, 0),
                         1:(0, 1),
                         2:(1, 0),
                         3:(1, 1),
                         4:(0, 0),
                         5:(0, 1),
                         6:(1, 0),
                         7:(1, 1),
                         8:(0, 0),
                         9:(0, 0),
                         10:(0, 1),
                         11:(0, 1)}
        
        ax_pairs = {0:4,
                    1:5,
                    2:6,
                    3:7,
                    4:0,
                    5:1,
                    6:2,
                    7:3,
                    8:9,
                    9:8,
                    10:11,
                    11:10}
        # make the axis an attribute
        self._ax = event.inaxes
        
        # find the component index so that it can be masked
        for ax_index, ax in enumerate(self.ax_list):
            if ax == event.inaxes:
                self._comp_index_x, self._comp_index_y = ax_index_dict[ax_index]
                self._ax_index = ax_index
                self._ax2 = self.ax_list[ax_pairs[ax_index]]
                if ax_index < 8:
                    self._key = 'z'
                    
                else:
                    self._key = 'tip'


#==============================================================================
#  Plot setting        
#==============================================================================
class PlotSettings(object):
    def __init__(self, **kwargs):
        
        self.fs = kwargs.pop('fs', 10)
        self.lw = kwargs.pop('lw', 1.5)
        self.ms = kwargs.pop('ms', 5)
        
        self.e_capthick = kwargs.pop('e_capthick', 1)
        self.e_capsize =  kwargs.pop('e_capsize', 5)

        #color mode
        self.cted = kwargs.pop('cted', (0, 0, .75))
        self.ctmd = kwargs.pop('ctmd', (.75, 0, 0))
        self.mted = kwargs.pop('mted', 's')
        self.mtmd = kwargs.pop('mtmd', 'o')
        
        #color for occam2d model
        self.ctem = kwargs.pop('ctem', (0, .6, .3))
        self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
        self.mtem = kwargs.pop('mtem', '+')
        self.mtmm = kwargs.pop('mtmm', '+')
         
        self.res_xx_limits = kwargs.pop('res_xx_limits', None)   
        self.res_xy_limits = kwargs.pop('res_xy_limits', None)   
        self.res_yx_limits = kwargs.pop('res_yx_limits', None)   
        self.res_yy_limits = kwargs.pop('res_yy_limits', None) 
        
        self.phase_xx_limits = kwargs.pop('phase_xx_limits', None)   
        self.phase_xy_limits = kwargs.pop('phase_xy_limits', None)   
        self.phase_yx_limits = kwargs.pop('phase_yx_limits', None)   
        self.phase_yy_limits = kwargs.pop('phase_yy_limits', None) 
        
        self.tipper_limits = kwargs.pop('tipper_limits', None)
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .2)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .08)
        
        self.z_err_increase = kwargs.pop('z_err_increase', 0.40)
        self.t_err_increase = kwargs.pop('t_err_increase', 0.10)
        
        self.legend_loc = kwargs.pop('legend_loc', 'upper left')
        self.legend_pos = kwargs.pop('legend_pos', None)
        self.legend_marker_scale = kwargs.pop('legend_marker_scale', 1)
        self.legend_border_axes_pad = kwargs.pop('legend_border_axes_pad', .01)
        self.legend_label_spacing = kwargs.pop('legend_label_spacing', 0.07)
        self.legend_handle_text_pad = kwargs.pop('legend_handle_text_pad', .05)
        self.legend_border_pad = kwargs.pop('legend_border_pad', .05)

        self.ylabel_pad = 1.25
#==============================================================================
# Def Main
#==============================================================================
def main():
    app = QtWidgets.QApplication(sys.argv)
    ui = ModEMPlotResponse()
    ui.show()
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    main()  
