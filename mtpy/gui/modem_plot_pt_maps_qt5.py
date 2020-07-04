# -*- coding: utf-8 -*-
"""
ModEM data and response visualization with a gui.

The user will be able to choose from stations within the data to look at
in either impedance or apparent resistivity and phase.

The functionality is quite simple at the moment

JP 2014
"""
# 
#==============================================================================
import os
import numpy as np

try:
    from PyQt5 import QtCore, QtGui, QtWidgets
except ImportError:
    raise ImportError("This version needs PyQt5")

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse
from matplotlib.colors import Normalize
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from matplotlib import cm

import mtpy.modeling.modem as modem
import mtpy.imaging.mtplottools as mtplottools
import mtpy.analysis.pt as mtpt
import mtpy.utils.exceptions as mtex

import mtpy.imaging.mtcolors as mtcl
import mtpy.analysis.niblettbostick as mtnb

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtWidgets.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtWidgets.QApplication.translate(context, text, disambig)
# 
#==============================================================================


class ModEMPlotPTMap(QtWidgets.QMainWindow, mtplottools.MTArrows,
                     mtplottools.MTEllipse):
    def __init__(self):
        
        super(ModEMPlotPTMap, self).__init__()
        
        self.modem_model_fn = None
        self.modem_data_fn = None
        self.modem_resp_fn = None
                
        self.save_plots = 'y'
        self.plot_period_index = None
        self.plot_period_list = None
        self.period_dict = None
        
        self.map_scale = 'km'
        #make map scale
        if self.map_scale == 'km':
            self.dscale = 1000.
            self.ellipse_size = 3
            self.arrow_size = 2
            self.arrow_head_length = .35
            self.arrow_head_width = .35
            self.arrow_lw = 1.
                                
        elif self.map_scale == 'm':
            self.dscale = 1.
            self.ellipse_size = 500
            self.arrow_size = 500
            self.arrow_head_length = 50
            self.arrow_head_width = 50
            self.arrow_lw = .75
        
        self.ew_limits = None
        self.ns_limits = None
        
        self.pad_east = 2*self.ellipse_size
        self.pad_north = 2*self.ellipse_size
        
        self.plot_grid = 'n'
        self.plot_stations = False
        
        
        self.xminorticks = 1000/self.dscale
        self.yminorticks = 1000/self.dscale
        
        self.residual_cmap = 'mt_wh2or'
        self.font_size = 9
        
        self.cb_tick_step = 45
        self.cb_residual_tick_step = 3
        self.cb_pt_pad = 1.25
        self.cb_res_pad = .65
        
        
        self.res_limits = (-1,4)
        self.res_cmap = 'jet_r'
        
        #--> set the ellipse properties -------------------
        
        self.subplot_right = .99
        self.subplot_left = .085
        self.subplot_top = .92
        self.subplot_bottom = .1
        self.subplot_hspace = .2
        self.subplot_wspace = .05
        
        # arrays to put data into        
        self.pt_data_arr = None
        self.pt_resp_arr = None
        self.pt_resid_arr = None

        self.dir_path = os.getcwd()
        
        self.setup_ui()
        
    def setup_ui(self):
        self.setWindowTitle("Plot ModEM MT Response as PT Maps")
        self.setWindowState(QtCore.Qt.WindowMaximized)
        
        #make a central widget that everything is tied to.
        self.central_widget = QtWidgets.QWidget(self)
        self.central_widget.setWindowTitle("Plot MT Response")
        
        #make a widget that will be the period list
        self.list_widget = QtWidgets.QListWidget()
        self.list_widget.itemClicked.connect(self.get_period)
        self.list_widget.currentItemChanged.connect(self.get_period)
        self.list_widget.setMaximumWidth(150)
        
        # make a depth text bar
        self.depth_label = QtWidgets.QLabel('Depth (m):')
        depth_font = QtGui.QFont()
        depth_font.setBold = True
        depth_font.setPointSize (16)
        self.depth_label.setFont(depth_font)
        
        self.depth_text = QtWidgets.QLabel('0.0')
        self.depth_text.setFont(depth_font)
        self.depth_text.setAlignment(QtCore.Qt.AlignCenter)
        depth_vbox = QtWidgets.QVBoxLayout()
        depth_vbox.addWidget(self.depth_label)
        depth_vbox.addWidget(self.depth_text)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        
        self.figure.subplots_adjust(left=self.subplot_left,
                                    right=self.subplot_right,
                                    bottom=self.subplot_bottom,
                                    top=self.subplot_top,
                                    hspace=self.subplot_hspace,
                                    wspace=self.subplot_wspace)
        
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
        left_layout.addLayout(depth_vbox)
        
        # set the layout the main window
        layout = QtWidgets.QHBoxLayout()
        layout.addLayout(left_layout)
        layout.addLayout(mpl_vbox)
        self.central_widget.setLayout(layout)

        #set the geometry of each widget        
        self.list_widget.setObjectName(_fromUtf8("listWidget"))
        self.mpl_widget.setObjectName(_fromUtf8("mpl_widget"))
        self.mpl_widget.updateGeometry()

        #set the central widget
        self.setCentralWidget(self.central_widget)

        #create a menu bar on the window
        self.menubar = QtWidgets.QMenuBar()
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1920, 38))
        self.menubar.setObjectName(_fromUtf8("menubar"))

        # add a tab for File --> open, close, save
        self.menu_data_file = QtWidgets.QMenu(self.menubar)
        self.menu_data_file.setTitle("Data File")
        
        self.menu_resp_file = QtWidgets.QMenu(self.menubar)
        self.menu_resp_file.setTitle("Response File")
        
        self.menu_model_file = QtWidgets.QMenu(self.menubar)
        self.menu_model_file.setTitle("Model File")
        
        # add a tab for chaning the display
        self.menu_display = QtWidgets.QMenu(self.menubar)
        self.menu_display.setTitle("Display")

        self.setMenuBar(self.menubar)

        # add a status bar on the bottom of the main window
        self.statusbar = QtWidgets.QStatusBar()
        self.statusbar.setObjectName(_fromUtf8("statusbar"))

        self.setStatusBar(self.statusbar)
        
        # set an open option that on click opens a modem file
        self.action_data_open = QtWidgets.QAction(self)
        self.action_data_open.setText("Open")
        self.action_data_open.triggered.connect(self.get_data_fn)

        # set a close that closes the main window
        self.action_close = QtWidgets.QAction(self)
        self.action_close.setText("Close")
        self.action_close.triggered.connect(self.close)

        # set a save option that will eventually save the masked data
        self.action_save = QtWidgets.QAction(self)
        self.action_save.setText("Save")

        # add the action on the menu tab
        self.menu_data_file.addAction(self.action_data_open)
        self.menu_data_file.addAction(self.action_close)
        self.menu_data_file.addAction(self.action_save)
        self.menubar.addAction(self.menu_data_file.menuAction())
        
        self.action_resp_open = QtWidgets.QAction(self)
        self.action_resp_open.setText("Open")
        self.action_resp_open.triggered.connect(self.get_resp_fn)
        self.menu_resp_file.addAction(self.action_resp_open)
        self.menubar.addAction(self.menu_resp_file.menuAction())
        
        self.action_model_open = QtWidgets.QAction(self)
        self.action_model_open.setText("Open")
        self.action_model_open.triggered.connect(self.get_model_fn)
        self.menu_model_file.addAction(self.action_model_open)
        self.menubar.addAction(self.menu_model_file.menuAction())

        self.action_plot_settings = QtWidgets.QAction(self)
        self.action_plot_settings.setText('Settings')
        self.action_plot_settings.triggered.connect(self.show_settings)
        self.menu_display.addAction(self.action_plot_settings)
        
        self.action_plot_stations = QtWidgets.QAction(self)
        self.action_plot_stations.setText('Plot Stations')
        self.action_plot_stations.setCheckable(True)
        self.action_plot_stations.toggled.connect(self.set_plot_stations)
        self.menu_display.addAction(self.action_plot_stations)
        
        self.menubar.addAction(self.menu_display.menuAction())
        
        # be sure to connnect all slots first
        QtCore.QMetaObject.connectSlotsByName(self)
        
                    
    def get_data_fn(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)')[0])
        
        fn = os.path.abspath(fn)
                                       
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
            
        self.plot_period = '{0:.5f}'.format(self.period_list[0])
        
        self._get_pt()
        
        self.get_depth_array()
            
    def get_model_fn(self):
        """
        get the filename from a file dialogue
        
        """        

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM model file',
                                           filter='(*.rho);; (*.ws)',
                                           directory=self.dir_path)[0])
        fn = os.path.abspath(fn)
        self.modem_model = modem.Model()
        self.modem_model.read_model_file(fn)
        self.modem_model_fn = fn
        self.get_depth_array()
        self.plot()
        
        
    def get_period(self, widget_item):
        """
        get the station name from the clicked station 
        """
        self.plot_period = '{0:.5f}'.format(float(str(widget_item.text()))) 
        self.plot()
        
    def get_resp_fn(self):
        """
        get response file name
        """

        fn_dialog = QtWidgets.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM response file',
                                           filter='*.dat', 
                                           directory=self.dir_path)[0])
                                       
        self.modem_resp = modem.Data()
        self.modem_resp.read_data_file(fn)
        self.modem_resp_fn = fn
        self._get_pt()
        self.plot()
        
    def show_settings(self):
        """
        show setting window
        """
        self.settings_window = PlotSettings(None, **self.__dict__)
        self.settings_window.show()
        self.settings_window.settings_updated.connect(self.update_settings)
        
    def update_settings(self):
        """
        update all the new plot settings
        """
        
        for attr in sorted(self.settings_window.__dict__.keys()):
            setattr(self, attr, self.settings_window.__dict__[attr])
            
        self.plot()
        
    def set_plot_stations(self, toggled):
        """
        plot station names if desired
        """
        
        self.plot_stations = toggled
        
        self.plot()
            
        
    def _get_pt(self):
        """
        put pt parameters into something useful for plotting
        """
        ns = len(list(self.modem_data.mt_dict.keys()))
        nf = len(self.modem_data.period_list)
        
        data_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                ('phimax', np.float),
                                                ('skew', np.float),
                                                ('azimuth', np.float),
                                                ('east', np.float),
                                                ('north', np.float),
                                                ('txr', np.float),
                                                ('tyr', np.float),
                                                ('txi', np.float),
                                                ('tyi', np.float)])
        if self.modem_resp_fn is not None:
            model_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                    ('phimax', np.float),
                                                    ('skew', np.float),
                                                    ('azimuth', np.float),
                                                    ('east', np.float),
                                                    ('north', np.float),
                                                    ('txr', np.float),
                                                    ('tyr', np.float),
                                                    ('txi', np.float),
                                                    ('tyi', np.float)])
            
            res_pt_arr = np.zeros((nf, ns), dtype=[('phimin', np.float),
                                                    ('phimax', np.float),
                                                    ('skew', np.float),
                                                    ('azimuth', np.float),
                                                    ('east', np.float),
                                                    ('north', np.float),
                                                    ('geometric_mean', np.float),
                                                    ('txr', np.float),
                                                    ('tyr', np.float),
                                                    ('txi', np.float),
                                                    ('tyi', np.float)])
                                                
        for ii, key in enumerate(self.modem_data.mt_dict.keys()):
            east = self.modem_data.mt_dict[key].grid_east/self.dscale
            north = self.modem_data.mt_dict[key].grid_north/self.dscale            
            dpt = self.modem_data.mt_dict[key].pt
            data_pt_arr[:, ii]['east'] = east
            data_pt_arr[:, ii]['north'] = north
            data_pt_arr[:, ii]['phimin'] = dpt.phimin
            data_pt_arr[:, ii]['phimax'] = dpt.phimax
            data_pt_arr[:, ii]['azimuth'] = dpt.azimuth
            data_pt_arr[:, ii]['skew'] = dpt.beta

            # compute tipper data
            tip = self.modem_data.mt_dict[key].Tipper
            tip.compute_mag_direction()
                        
            data_pt_arr[:, ii]['txr'] = tip.mag_real*\
                                        np.sin(np.deg2rad(tip.angle_real))
            data_pt_arr[:, ii]['tyr'] = tip.mag_real*\
                                        np.cos(np.deg2rad(tip.angle_real))
            data_pt_arr[:, ii]['txi'] = tip.mag_imag*\
                                        np.sin(np.deg2rad(tip.angle_imag))
            data_pt_arr[:, ii]['tyi'] = tip.mag_imag*\
                                        np.cos(np.deg2rad(tip.angle_imag))
            if self.modem_resp_fn is not None:
                mpt = self.modem_resp.mt_dict[key].pt
                
                model_pt_arr[:, ii]['east'] = east
                model_pt_arr[:, ii]['north'] = north
                model_pt_arr[:, ii]['phimin'] = mpt.phimin
                model_pt_arr[:, ii]['phimax'] = mpt.phimax
                model_pt_arr[:, ii]['azimuth'] = mpt.azimuth
                model_pt_arr[:, ii]['skew'] = mpt.beta
                
                mtip = self.modem_resp.mt_dict[key].Tipper
                mtip.compute_mag_direction()
                            
                model_pt_arr[:, ii]['txr'] = mtip.mag_real*\
                                            np.sin(np.deg2rad(mtip.angle_real))
                model_pt_arr[:, ii]['tyr'] = mtip.mag_real*\
                                            np.cos(np.deg2rad(mtip.angle_real))
                model_pt_arr[:, ii]['txi'] = mtip.mag_imag*\
                                            np.sin(np.deg2rad(mtip.angle_imag))
                model_pt_arr[:, ii]['tyi'] = mtip.mag_imag*\
                                            np.cos(np.deg2rad(mtip.angle_imag))
                try:
                    rpt = mtpt.ResidualPhaseTensor(pt_object1=dpt, 
                                                   pt_object2=mpt)
                    rpt = rpt.residual_pt
                    
                    res_pt_arr[:, ii]['phimin'] = rpt.phimin
                    res_pt_arr[:, ii]['phimax'] = rpt.phimax
                    res_pt_arr[:, ii]['azimuth'] = rpt.azimuth
                    res_pt_arr[:, ii]['skew'] = rpt.beta
                    res_pt_arr[:, ii]['geometric_mean'] = np.sqrt(abs(rpt.phimin*\
                                                                  rpt.phimax))
                                                                  
            
                except mtex.MTpyError_PT:
                    print('Could not calculate residual PT for {0}'.format(key))
                    
                res_pt_arr[:, ii]['east'] = east
                res_pt_arr[:, ii]['north'] = north
                    
                res_pt_arr[:, ii]['txr'] = data_pt_arr[:, ii]['txr']-\
                                            model_pt_arr[:, ii]['txr']
                res_pt_arr[:, ii]['tyr'] = data_pt_arr[:, ii]['tyr']-\
                                            model_pt_arr[:, ii]['tyr']
                res_pt_arr[:, ii]['txi'] = data_pt_arr[:, ii]['txi']-\
                                            model_pt_arr[:, ii]['txi']
                res_pt_arr[:, ii]['tyi'] = data_pt_arr[:, ii]['tyi']-\
                                            model_pt_arr[:, ii]['tyi']
                
                
        #make these attributes   
        self.pt_data_arr = data_pt_arr
        
        if self.modem_resp_fn is not None:
            self.pt_resp_arr = model_pt_arr
            self.pt_resid_arr = res_pt_arr
            
    def get_depth_array(self):
        """
        estimate a niblett-bostick depth from the impedance tensors
        
        find the average depth for each station at each period
        """
        if self.modem_data.mt_dict is None:
            return
            
        d_arr_min = np.zeros((self.modem_data.period_list.shape[0],
                              len(list(self.modem_data.mt_dict.keys()))))
        d_arr_max = np.zeros((self.modem_data.period_list.shape[0],
                              len(list(self.modem_data.mt_dict.keys()))))
#        print self.modem_data.mt_dict[self.modem_data.mt_dict.keys()[0]].Z.z                      
        for ii, mt_key in enumerate(sorted(self.modem_data.mt_dict.keys())):
            mt_obj = self.modem_data.mt_dict[mt_key]
            d_arr = mtnb.calculate_depth_nb(z_object=mt_obj.Z)
            
            d_arr_min[:, ii] = d_arr['depth_min']
            d_arr_max[:, ii] = d_arr['depth_max']
        
        # average only the non zero terms
        d_avg_min = np.array([d_arr_min[kk, np.nonzero(d_arr_min[kk, :])].mean()
                              for kk in range(d_arr_min.shape[0])])
        d_avg_max = np.array([d_arr_max[kk, np.nonzero(d_arr_max[kk, :])].mean()
                              for kk in range(d_arr_min.shape[0])])

        # find the average and leave in meters cause grid_z is in meters
        d_avg = np.nan_to_num(((d_avg_min+d_avg_max)/2.))
        
        self.depth_array = d_avg
                
    def plot(self):
        """
        plot phase tensor maps for data and or response, each figure is of a
        different period.  If response is input a third column is added which is 
        the residual phase tensor showing where the model is not fitting the data 
        well.  The data is plotted in km.
        
        """
        
        plt.rcParams['font.size'] = self.font_size
        
        #make sure there is PT data
        if self.pt_data_arr is None:
            self._get_pt()
            
        if self.modem_resp_fn is not None:
            if self.pt_resp_arr is None:
                self._get_pt()
        
        # make a grid of subplots 
        gs = gridspec.GridSpec(1, 3, hspace=self.subplot_hspace,
                               wspace=self.subplot_wspace)
                               
        font_dict = {'size':self.font_size+2, 'weight':'bold'}
        
        #set some parameters for the colorbar
        ckmin = float(self.ellipse_range[0])
        ckmax = float(self.ellipse_range[1])
        try:
            ckstep = float(self.ellipse_range[2])
        except IndexError:
            if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                raise ValueError('Need to input range as (min, max, step)')
            else:
                ckstep = 3
        bounds = np.arange(ckmin, ckmax+ckstep, ckstep)
        nseg = float((ckmax-ckmin)/(2*ckstep))
        
        # set plot limits to be the station area
        if self.ew_limits == None:
            east_min = self.pt_data_arr['east'].min()-self.pad_east
            east_max = self.pt_data_arr['east'].max()+self.pad_east

            self.ew_limits = (east_min, east_max)
            
        if self.ns_limits == None:
            north_min = self.pt_data_arr['north'].min()-self.pad_north
            north_max = self.pt_data_arr['north'].max()+self.pad_north

            self.ns_limits = (north_min, north_max)
        #-------------plot phase tensors------------------------------------                    
        data_ii = self.period_dict[self.plot_period]
        print('Ploting period {0}'.format(data_ii))
        
        self.figure.clf()
                         
        if self.modem_resp_fn is not None:
            
            axd = self.figure.add_subplot(gs[0, 0], aspect='equal')
            axm = self.figure.add_subplot(gs[0, 1], 
                                          aspect='equal',
                                          sharex=axd, 
                                          sharey=axd)
            axr = self.figure.add_subplot(gs[0, 2], 
                                          aspect='equal',
                                          sharex=axd,
                                          sharey=axd)
            ax_list = [axd, axm, axr]
#        
        else:
            axd = self.figure.add_subplot(gs[0, :], aspect='equal')
            ax_list = [axd]
            
        arr_dir = (-1)**self.arrow_direction
        
        #plot model below the phase tensors
        if self.modem_model_fn is not None:
            #self.get_depth_array()
            if self.depth_array[data_ii] == 0:
                print('Could not estimate depth for period {0:.5g}'.format(
                    float(self.plot_period)))
                d_index = 0
            else:
                try:
                    d_index = np.where(self.modem_model.grid_z >= 
                                        self.depth_array[data_ii])[0][0]
                                        
                    print('Estimated depth for period {0:.5g} is {1:.2f} m'.format(
                            float(self.plot_period), self.depth_array[data_ii]))
                            
                    self.depth_text.setText('{0:.5g}'.format(self.depth_array[data_ii]))
                    
                except IndexError:
                    print('Could not estimate depth for period {0:.2f}'.format(
                            float(self.plot_period)))
                    d_index = 0
            
            #need to add an extra row and column to east and north to make sure 
            #all is plotted see pcolor for details.
            plot_east = self.modem_model.grid_east/self.dscale
            plot_north = self.modem_model.grid_north/self.dscale

            
            #make a mesh grid for plotting
            #the 'ij' makes sure the resulting grid is in east, north
            self.mesh_east, self.mesh_north = np.meshgrid(plot_east, 
                                                          plot_north,
                                                          indexing='ij')
            
            for ax in ax_list:
                plot_res = np.log10(self.modem_model.res_model[:, :, d_index].T)
                ax.pcolormesh(self.mesh_east,
                               self.mesh_north, 
                               plot_res,
                               cmap=self.res_cmap,
                               vmin=self.res_limits[0],
                               vmax=self.res_limits[1])
                
            
        #--> plot data phase tensors
        for pt in self.pt_data_arr[data_ii]:
            if pt['phimin'] == 0 and pt['phimax'] == 0:
                pass
            else:
                eheight = pt['phimin']/\
                          self.pt_data_arr[data_ii]['phimax'].max()*\
                          self.ellipse_size
                ewidth = pt['phimax']/\
                          self.pt_data_arr[data_ii]['phimax'].max()*\
                          self.ellipse_size
                          
                ellipse = Ellipse((pt['east'],
                                   pt['north']),
                                   width=ewidth,
                                   height=eheight,
                                   angle=90-pt['azimuth'])
                
                #get ellipse color
                if self.ellipse_cmap.find('seg')>0:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax,
                                                         bounds=bounds))
                else:
                    ellipse.set_facecolor(mtcl.get_plot_color(pt[self.ellipse_colorby],
                                                         self.ellipse_colorby,
                                                         self.ellipse_cmap,
                                                         ckmin,
                                                         ckmax))
                
                axd.add_artist(ellipse)
            
            #-----------Plot Induction Arrows---------------------------
            if pt['txr'] != 0.0:
                real_mag = np.sqrt(abs(pt['txr'])**2+abs(pt['tyr'])**2)
                imag_mag = np.sqrt(abs(pt['txi'])**2+abs(pt['tyi'])**2)
                #plot real tipper
                if real_mag <= self.arrow_threshold:
                    axd.arrow(pt['east'],
                              pt['north'],
                              self.arrow_size*pt['txr']*arr_dir,
                              self.arrow_size*pt['tyr']*arr_dir,
                              lw=self.arrow_lw,
                              facecolor=self.arrow_color_real,
                              edgecolor=self.arrow_color_real,
                              length_includes_head=False,
                              head_width=self.arrow_head_width,
                              head_length=self.arrow_head_length)
                else:
                    pass
                    
                #plot imaginary tipper
                if imag_mag <= self.arrow_threshold:
                    axd.arrow(pt['east'],
                              pt['north'],
                              self.arrow_size*pt['txi']*arr_dir,
                              self.arrow_size*pt['tyi']*arr_dir,
                              lw=self.arrow_lw,
                              facecolor=self.arrow_color_imag,
                              edgecolor=self.arrow_color_imag,
                              length_includes_head=False,
                              head_width=self.arrow_head_width,
                              head_length=self.arrow_head_length)
                else:
                    pass
                
        #-----------plot response phase tensors---------------
        if self.modem_resp_fn is not None:
            rcmin = np.floor(self.pt_resid_arr['geometric_mean'].min())
            rcmax = np.floor(self.pt_resid_arr['geometric_mean'].max())
            for mpt, rpt in zip(self.pt_resp_arr[data_ii], 
                                self.pt_resid_arr[data_ii]):
                if mpt['phimin'] == 0 and mpt['phimax'] == 0:
                    pass
                else:
                    eheight = mpt['phimin']/\
                              self.pt_data_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                    ewidth = mpt['phimax']/\
                              self.pt_data_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                              
                    ellipsem = Ellipse((mpt['east'],
                                       mpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=90-mpt['azimuth'])
                    
                    #get ellipse color
                    if self.ellipse_cmap.find('seg')>0:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                             self.ellipse_colorby,
                                                             self.ellipse_cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                    else:
                        ellipsem.set_facecolor(mtcl.get_plot_color(mpt[self.ellipse_colorby],
                                                             self.ellipse_colorby,
                                                             self.ellipse_cmap,
                                                             ckmin,
                                                             ckmax))
                
                    axm.add_artist(ellipsem)
                
                #-----------Plot Induction Arrows---------------------------
                if mpt['txr'] != 0.0:
                    real_mag = np.sqrt(abs(mpt['txr'])**2+abs(mpt['tyr'])**2)
                    imag_mag = np.sqrt(abs(mpt['txi'])**2+abs(mpt['tyi'])**2)
                    #plot real tipper
                    if real_mag <= self.arrow_threshold:
                        axm.arrow(mpt['east'],
                                  mpt['north'],
                                  self.arrow_size*mpt['txr']*arr_dir,
                                  self.arrow_size*mpt['tyr']*arr_dir,
                                  lw=self.arrow_lw,
                                  facecolor=self.arrow_color_real,
                                  edgecolor=self.arrow_color_real,
                                  length_includes_head=False,
                                  head_width=self.arrow_head_width,
                                  head_length=self.arrow_head_length)
                    else:
                        pass
                        
                    #plot imaginary tipper
                    if imag_mag <= self.arrow_threshold:
                        axm.arrow(mpt['east'],
                                  mpt['north'],
                                  self.arrow_size*mpt['txi']*arr_dir,
                                  self.arrow_size*mpt['tyi']*arr_dir,
                                  lw=self.arrow_lw,
                                  facecolor=self.arrow_color_imag,
                                  edgecolor=self.arrow_color_imag,
                                  length_includes_head=False,
                                  head_width=self.arrow_head_width,
                                  head_length=self.arrow_head_length)
                    else:
                        pass
                
                #-----------plot residual phase tensors---------------
                if mpt['phimin'] == 0 and mpt['phimax'] == 0:
                    pass
                else:
                    eheight = rpt['phimin']/\
                              self.pt_data_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                    ewidth = rpt['phimax']/\
                              self.pt_data_arr[data_ii]['phimax'].max()*\
                              self.ellipse_size
                              
                    ellipser = Ellipse((rpt['east'],
                                       rpt['north']),
                                       width=ewidth,
                                       height=eheight,
                                       angle=rpt['azimuth'])
                    
                    #get ellipse color
                    rpt_color = np.sqrt(abs(rpt['phimin']*rpt['phimax']))
                    if self.ellipse_cmap.find('seg')>0:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                             'geometric_mean',
                                                             self.residual_cmap,
                                                             ckmin,
                                                             ckmax,
                                                             bounds=bounds))
                    else:
                        ellipser.set_facecolor(mtcl.get_plot_color(rpt_color,
                                                             'geometric_mean',
                                                             self.residual_cmap,
                                                             ckmin,
                                                             ckmax))
                    
                    
                    axr.add_artist(ellipser)
                
                #-----------Plot Induction Arrows---------------------------
                if rpt['txr'] != 0.0:
                    real_mag = np.sqrt(abs(rpt['txr'])**2+abs(rpt['tyr'])**2)
                    imag_mag = np.sqrt(abs(rpt['txi'])**2+abs(rpt['tyi'])**2)
                    #plot real tipper
                    if real_mag <= self.arrow_threshold:
                        axr.arrow(rpt['east'],
                                  rpt['north'],
                                  self.arrow_size*rpt['txr']*arr_dir,
                                  self.arrow_size*rpt['tyr']*arr_dir,
                                  lw=self.arrow_lw,
                                  facecolor=self.arrow_color_real,
                                  edgecolor=self.arrow_color_real,
                                  length_includes_head=False,
                                  head_width=self.arrow_head_width,
                                  head_length=self.arrow_head_length)
                    else:
                        pass
                        
                    #plot imaginary tipper
                    if imag_mag <= self.arrow_threshold:
                        axr.arrow(rpt['east'],
                                  rpt['north'],
                                  self.arrow_size*rpt['txi']*arr_dir,
                                  self.arrow_size*rpt['tyi']*arr_dir,
                                  lw=self.arrow_lw,
                                  facecolor=self.arrow_color_imag,
                                  edgecolor=self.arrow_color_imag,
                                  length_includes_head=False,
                                  head_width=self.arrow_head_width,
                                  head_length=self.arrow_head_length)
                    else:
                        pass
                
        #--> set axes properties
        # data
        axd.set_xlim(self.ew_limits)
        axd.set_ylim(self.ns_limits)

        axd.set_xlabel('Easting ({0})'.format(self.map_scale), 
                       fontdict=font_dict)
        axd.set_ylabel('Northing ({0})'.format(self.map_scale),
                       fontdict=font_dict)
        #make a colorbar for phase tensors
        #bb = axd.axes.get_position().bounds
        bb = axd.get_position().bounds
        y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                 (self.ew_limits[1]-self.ew_limits[0]))
        cb_location = (3.35*bb[2]/5+bb[0], 
                        y1*self.cb_pt_pad, .295*bb[2], .02)
        cbaxd = self.figure.add_axes(cb_location)
        
        if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
            #make a color list
            clist = [(cc, cc ,1) 
                         for cc in np.arange(0, 1+1./(nseg), 1./(nseg))]+\
                       [(1, cc, cc) 
                         for cc in np.arange(1, -1./(nseg), -1./(nseg))]
            
            #make segmented colormap
            mt_seg_bl2wh2rd = colors.ListedColormap(clist)

            #make bounds so that the middle is white
            bounds = np.arange(ckmin-ckstep, ckmax+2*ckstep, ckstep)
            
            #normalize the colors
            norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)
            
            #make the colorbar
            cb_pt = mcb.ColorbarBase(cbaxd,
                                   cmap=mt_seg_bl2wh2rd,
                                   norm=norms,
                                   orientation='horizontal',
                                   ticks=bounds[1:-1])
        else:
            
            cb_pt = mcb.ColorbarBase(cbaxd, 
                                   cmap=mtcl.cmapdict[self.ellipse_cmap],
                                   norm=Normalize(vmin=ckmin,
                                                  vmax=ckmax),
                                   orientation='horizontal')
        cb_pt.ax.xaxis.set_label_position('top')
        cb_pt.ax.xaxis.set_label_coords(.5, 1.75)
        cb_pt.set_label(mtplottools.ckdict[self.ellipse_colorby])
        cb_pt.set_ticks([ckmin, (ckmax-ckmin)/2, ckmax])
                                
        axd.text(self.ew_limits[0]*.95,
                 self.ns_limits[1]*.95,
                 'Data',
                 horizontalalignment='left',
                 verticalalignment='top',
                 bbox={'facecolor':'white'},
                 fontdict={'size':self.font_size+1})
                
        #Model and residual
        if self.modem_resp_fn is not None:
            for aa, ax in enumerate([axm, axr]):
                ax.set_xlabel('Easting ({0})'.format(self.map_scale), 
                               fontdict=font_dict)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                #make a colorbar ontop of axis
                bb = ax.axes.get_position().bounds
                y1 = .25*(2+(self.ns_limits[1]-self.ns_limits[0])/
                         (self.ew_limits[1]-self.ew_limits[0]))
                cb_location = (3.35*bb[2]/5+bb[0], 
                               y1*self.cb_pt_pad, .295*bb[2], .02)
                cbax = self.figure.add_axes(cb_location)
                if aa == 0:
                    if self.ellipse_cmap == 'mt_seg_bl2wh2rd':
                        #make a color list
                        clist = [(cc, cc ,1) 
                                     for cc in np.arange(0, 1+1./(nseg), 1./(nseg))]+\
                                   [(1, cc, cc) 
                                     for cc in np.arange(1, -1./(nseg), -1./(nseg))]
                        
                        #make segmented colormap
                        mt_seg_bl2wh2rd = colors.ListedColormap(clist)
            
                        #make bounds so that the middle is white
                        bounds = np.arange(ckmin-ckstep, ckmax+2*ckstep, ckstep)
                        
                        #normalize the colors
                        norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)
                        
                        #make the colorbar
                        cb_ptr = mcb.ColorbarBase(cbax,
                                               cmap=mt_seg_bl2wh2rd,
                                               norm=norms,
                                               orientation='horizontal',
                                               ticks=bounds[1:-1])
                    else:
                        cb_ptr = mcb.ColorbarBase(cbax, 
                                              cmap=mtcl.cmapdict[self.ellipse_cmap],
                                              norm=Normalize(vmin=ckmin,
                                                             vmax=ckmax),
                                               orientation='horizontal')
                    cb_ptr.ax.xaxis.set_label_position('top')
                    cb_ptr.ax.xaxis.set_label_coords(.5, 1.75)
                    cb_ptr.set_label(mtplottools.ckdict[self.ellipse_colorby])
                    cb_ptr.set_ticks(np.arange(ckmin, ckmax+self.cb_tick_step, 
                                self.cb_tick_step))
                    ax.text(self.ew_limits[0]*.95,
                            self.ns_limits[1]*.95,
                            'Model',
                            horizontalalignment='left',
                            verticalalignment='top',
                            bbox={'facecolor':'white'},
                             fontdict={'size':self.font_size+1})
                else:
                    cb_ptr = mcb.ColorbarBase(cbax, 
                                          cmap=mtcl.cmapdict[self.residual_cmap],
                                           norm=Normalize(vmin=rcmin,
                                                          vmax=rcmax),
                                           orientation='horizontal')
                    cb_ptr.ax.xaxis.set_label_position('top')
                    cb_ptr.ax.xaxis.set_label_coords(.5, 1.75)
                    cb_ptr.set_label(r"$\sqrt{\Phi_{min} \Phi_{max}}$")
                    cb_ptr.set_ticks([rcmin, (rcmax-rcmin)/2, rcmax])
                    ax.text(self.ew_limits[0]*.95,
                            self.ns_limits[1]*.95,
                            'Residual',
                            horizontalalignment='left',
                            verticalalignment='top',
                            bbox={'facecolor':'white'},
                            fontdict={'size':self.font_size+1})
        
        if self.modem_model_fn is not None:
            for ax in ax_list:
                ax.tick_params(direction='out')
                bb = ax.axes.get_position().bounds
                y1 = .25*(2-(self.ns_limits[1]-self.ns_limits[0])/
                         (self.ew_limits[1]-self.ew_limits[0]))
                cb_position = (3.0*bb[2]/5+bb[0], 
                               y1*self.cb_res_pad, .35*bb[2], .02)
                cbax = self.figure.add_axes(cb_position)
                cb_res = mcb.ColorbarBase(cbax, 
                                      cmap=cm.get_cmap(self.res_cmap),
                                      norm=Normalize(vmin=self.res_limits[0],
                                                     vmax=self.res_limits[1]),
                                      orientation='horizontal')
                cb_res.ax.xaxis.set_label_position('top')
                cb_res.ax.xaxis.set_label_coords(.5, 1.5)
                cb_res.set_label('Resistivity ($\Omega \cdot$m)')
                cb_ticks = np.arange(np.floor(self.res_limits[0]), 
                                     np.ceil(self.res_limits[1]+1), 1)
                cb_res.set_ticks(cb_ticks)
                cb_res.set_ticklabels([mtplottools.labeldict[ctk] for ctk in cb_ticks])

        if self.plot_stations == True:
            for ax in ax_list:
                for s_arr in self.modem_data.station_locations.station_locations:
                    ax.text(s_arr['rel_east']/self.dscale,
                            s_arr['rel_north']/self.dscale,
                            s_arr['station'],
                            horizontalalignment='center',
                            verticalalignment='baseline',
                            fontdict={'size':self.font_size},
                            clip_on=True) 
                
        # draw plot
        self.mpl_widget.draw()
        
class PlotSettings(QtWidgets.QWidget):
    settings_updated = QtCore.pyqtSignal()
    def __init__(self, parent, **kwargs):
        super(PlotSettings, self).__init__(parent)
        
        self.font_size = kwargs.pop('font_size', 10)
        
        self.map_scale = kwargs.pop('map_scale', 'km')
        
        if self.map_scale == 'km': 
            self.ellipse_size = kwargs.pop('ellipse_size', 1)
            self.arrow_head_length = kwargs.pop('arrow_head_length', .025)
            self.arrow_head_width = kwargs.pop('arrow_head_width', .025)
            self.arrow_size = kwargs.pop('arrow_size', 1)
            self.ew_limits = kwargs.pop('ew_limits', [-10, 10])
            self.ns_limits = kwargs.pop('ns_limits', [-10, 10])
            
        if self.map_scale == 'm': 
            self.ellipse_size = kwargs.pop('ellipse_size', 500)
            self.arrow_head_length = kwargs.pop('arrow_head_length', 50)
            self.arrow_head_width = kwargs.pop('arrow_head_width', 50)
            self.arrow_size = kwargs.pop('arrow_size', 500)
            self.ew_limits = kwargs.pop('ew_limits', [-10000, 10000])
            self.ns_limits = kwargs.pop('ns_limits', [-10000, 10000])
            
        if type(self.ns_limits) is tuple:
            self.ns_limits = list(self.ns_limits)
        if type(self.ew_limits) is tuple:
            self.ew_limits = list(self.ew_limits)
            
        self.ellipse_cmap = kwargs.pop('ellipse_cmap', 'mt_bl2wh2rd')
        self.ellipse_range = kwargs.pop('ellipse_range', [0, 90, 5])
        self.ellipse_colorby = kwargs.pop('ellipse_colorby', 'phimin')
        
        if type(self.ellipse_range) == tuple:
            self.ellipse_range = list(self.ellipse_range)
        
        self.arrow_threshold = kwargs.pop('arrow_threshold', 2)
        self.arrow_color_imag = kwargs.pop('arrow_color_imag', 'b')
        self.arrow_color_real = kwargs.pop('arrow_color_real', 'k')
        self.arrow_direction = kwargs.pop('arrow_direction', 0)
        self.arrow_lw = kwargs.pop('arrow_lw', .75)
        
        self.cb_pt_pad = kwargs.pop('cb_pt_pad', 0.5)
        self.cb_res_pad = kwargs.pop('cb_res_pad', 1.2)
        
        self.res_limits = kwargs.pop('res_limits', [0, 4])
        if type(self.res_limits) is tuple:
            self.res_limits = list(self.res_limits)
        
        self.subplot_wspace = kwargs.pop('subplot_wspace', .2)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .08)
        
        #--> run the init function to make the gui
        self.initUI()

    def initUI(self):
        #--> line properties
        fs_label = QtWidgets.QLabel('Font Size')
        fs_edit = QtWidgets.QLineEdit()
        fs_edit.setText('{0:.1f}'.format(self.font_size))
        fs_edit.textChanged[str].connect(self.set_text_fs)
        
        #--> Map properties
        mapscale_label = QtWidgets.QLabel('Map Scale')
        mapscale_combo = QtWidgets.QComboBox()
        mapscale_combo.addItem('km')
        mapscale_combo.addItem('m')
        mapscale_combo.activated[str].connect(self.set_mapscale) 
        
        ew_limits_label = QtWidgets.QLabel('E-W Limits (min, max)')
        ew_limits_min_edit = QtWidgets.QLineEdit()
        ew_limits_min_edit.setText('{0:.3f}'.format(self.ew_limits[0]))
        ew_limits_min_edit.textChanged[str].connect(self.set_ew_limits_min)
        
        ew_limits_max_edit = QtWidgets.QLineEdit()
        ew_limits_max_edit.setText('{0:.3f}'.format(self.ew_limits[1]))
        ew_limits_max_edit.textChanged[str].connect(self.set_ew_limits_max)
        
        ew_limits_grid = QtWidgets.QGridLayout()
        ew_limits_grid.setSpacing(5)
        ew_limits_grid.addWidget(ew_limits_label, 1, 0)
        ew_limits_grid.addWidget(ew_limits_min_edit, 1, 1)
        ew_limits_grid.addWidget(ew_limits_max_edit, 1, 2)
        
        ns_limits_label = QtWidgets.QLabel('N-S Limits (min, max)')
        ns_limits_min_edit = QtWidgets.QLineEdit()
        ns_limits_min_edit.setText('{0:.3f}'.format(self.ns_limits[0]))
        ns_limits_min_edit.textChanged[str].connect(self.set_ns_limits_min)
        
        ns_limits_max_edit = QtWidgets.QLineEdit()
        ns_limits_max_edit.setText('{0:.3f}'.format(self.ns_limits[1]))
        ns_limits_max_edit.textChanged[str].connect(self.set_ns_limits_max)
        
        ns_limits_grid = QtWidgets.QGridLayout()
        ns_limits_grid.setSpacing(5)
        ns_limits_grid.addWidget(ns_limits_label, 1, 0)
        ns_limits_grid.addWidget(ns_limits_min_edit, 1, 1)
        ns_limits_grid.addWidget(ns_limits_max_edit, 1, 2)
        
        grid_line = QtWidgets.QGridLayout()
        grid_line.setSpacing(10)
        
        grid_line.addWidget(fs_label, 1, 0)
        grid_line.addWidget(fs_edit, 1, 1)
        
        grid_line.addWidget(mapscale_label, 1, 2)
        grid_line.addWidget(mapscale_combo, 1, 3)
        
        grid_line.addLayout(ew_limits_grid, 1, 4)
        
        grid_line.addLayout(ns_limits_grid, 1, 5)
        
        #--> ellipse properties
        ellipse_size_label = QtWidgets.QLabel('Ellipse Size')
        ellipse_size_edit = QtWidgets.QLineEdit()
        ellipse_size_edit.setText('{0:.2f}'.format(self.ellipse_size))
        ellipse_size_edit.textChanged[str].connect(self.set_ellipse_size)
        
        ellipse_range_label = QtWidgets.QLabel('Ellipse Range (min, max, step)')
        
        ellipse_range_edit_min = QtWidgets.QLineEdit()
        try:
            ellipse_range_edit_min.setText('{0:.2f}'.format(self.ellipse_range[0]))
        except IndexError:
            if self.ellipse_colorby == 'skew':
                ellipse_range_edit_min.setText('{0:.2f}'.format(-9))
                self.ellipse_range = [-9.0]
            else:
                ellipse_range_edit_min.setText('{0:.2f}'.format(0))
                self.ellipse_range = [0]
        ellipse_range_edit_min.textChanged[str].connect(self.set_ellipse_range_min)
        
        ellipse_range_edit_max = QtWidgets.QLineEdit()
        try:
            ellipse_range_edit_max.setText('{0:.2f}'.format(self.ellipse_range[1]))
        except IndexError:
            if self.ellipse_colorby == 'skew':
                ellipse_range_edit_max.setText('{0:.2f}'.format(9))
                self.ellipse_range.append(9.0)
            else:
                ellipse_range_edit_max.setText('{0:.2f}'.format(90))
                self.ellipse_range.append(90.0)
        ellipse_range_edit_max.textChanged[str].connect(self.set_ellipse_range_max)
        
        ellipse_range_edit_step = QtWidgets.QLineEdit()
        try:
            ellipse_range_edit_step.setText('{0:.2f}'.format(self.ellipse_range[2]))
        except IndexError:
            if self.ellipse_colorby == 'skew':
                ellipse_range_edit_step.setText('{0:.2f}'.format(3))
                self.ellipse_range.append(3.0)
            else:
                ellipse_range_edit_step.setText('{0:.2f}'.format(5))
                self.ellipse_range.append(5)
        ellipse_range_edit_step.textChanged[str].connect(self.set_ellipse_range_step)

        range_grid = QtWidgets.QGridLayout()
        range_grid.setSpacing(5)
        range_grid.addWidget(ellipse_range_edit_min, 1, 0)
        range_grid.addWidget(ellipse_range_edit_max, 1, 1)
        range_grid.addWidget(ellipse_range_edit_step, 1, 2)

        ellipse_colorby_label = QtWidgets.QLabel('Ellipse Color By')
        ellipse_colorby_combo = QtWidgets.QComboBox()
        ellipse_colorby_combo.addItem('phimin')
        ellipse_colorby_combo.addItem('phimax')
        ellipse_colorby_combo.addItem('ellipticty')
        ellipse_colorby_combo.addItem('skew')
        ellipse_colorby_combo.activated[str].connect(self.set_ellipse_colorby)
        
        ellipse_cmap_label = QtWidgets.QLabel('Ellipse Color Map')
        ellipse_cmap_combo = QtWidgets.QComboBox()
        ellipse_cmap_combo.addItem('mt_bl2wh2rd')
        ellipse_cmap_combo.addItem('mt_yl2rd')
        ellipse_cmap_combo.addItem('mt_wh2bl')
        ellipse_cmap_combo.addItem('mt_bl2gr2rd')
        ellipse_cmap_combo.addItem('mt_rd2gr2bl')
        ellipse_cmap_combo.addItem('mt_seg_bl2wh2rd')
        ellipse_cmap_combo.activated[str].connect(self.set_ellipse_cmap)
        
        ellipse_grid = QtWidgets.QGridLayout()
        ellipse_grid.setSpacing(10)
        
        ellipse_grid.addWidget(ellipse_size_label, 1, 0)
        ellipse_grid.addWidget(ellipse_size_edit, 1, 1)
        
        ellipse_grid.addWidget(ellipse_range_label, 1, 2)
        ellipse_grid.addLayout(range_grid, 1, 3)
        
        ellipse_grid.addWidget(ellipse_colorby_label, 1, 4)
        ellipse_grid.addWidget(ellipse_colorby_combo, 1, 5)
        
        ellipse_grid.addWidget(ellipse_cmap_label, 1, 6)
        ellipse_grid.addWidget(ellipse_cmap_combo, 1, 7)
        
        #--> arrow settings
        arrow_size_label = QtWidgets.QLabel('Induction Arrow Size')
        arrow_size_edit = QtWidgets.QLineEdit()
        arrow_size_edit.setText('{0:.2f}'.format(self.arrow_size))
        arrow_size_edit.textChanged[str].connect(self.set_arrow_size)
        
        arrow_lw_label = QtWidgets.QLabel('Arrow Line Width')
        arrow_lw_edit = QtWidgets.QLineEdit()
        arrow_lw_edit.setText('{0:.2f}'.format(self.arrow_lw))
        arrow_lw_edit.textChanged[str].connect(self.set_arrow_lw)
        
        arrow_head_length_label = QtWidgets.QLabel('Arrow Head Size')
        arrow_head_length_edit = QtWidgets.QLineEdit()
        arrow_head_length_edit.setText('{0:.2f}'.format(self.arrow_head_length))
        arrow_head_length_edit.textChanged[str].connect(self.set_arrow_head_length)
        
        arrow_head_width_label = QtWidgets.QLabel('Arrow Head Width')
        arrow_head_width_edit = QtWidgets.QLineEdit()
        arrow_head_width_edit.setText('{0:.2f}'.format(self.arrow_head_width))
        arrow_head_width_edit.textChanged[str].connect(self.set_arrow_head_width)
        
        arrow_direction_label = QtWidgets.QLabel('Arrow Direction')
        arrow_direction_combo = QtWidgets.QComboBox()
        arrow_direction_combo.addItem('Parkinson')
        arrow_direction_combo.addItem('Weise')
        arrow_direction_combo.activated[str].connect(self.set_arrow_direction)
        
        arrow_threshold_label = QtWidgets.QLabel('Arrow Threshold')
        arrow_threshold_edit = QtWidgets.QLineEdit()
        arrow_threshold_edit.setText('{0:.2f}'.format(self.arrow_threshold))
        arrow_threshold_edit.textChanged[str].connect(self.set_arrow_threshold)
        
        arrow_color_real = QtWidgets.QPushButton('Arrow Color Real', self)
        arrow_color_real.clicked.connect(self.get_arrow_color_real)    
        
        arrow_color_imag = QtWidgets.QPushButton('Arrow Color Imaginary', self)
        arrow_color_imag.clicked.connect(self.get_arrow_color_imag)        
        
        arrow_grid = QtWidgets.QGridLayout()
        arrow_grid.setSpacing(10)
        
        arrow_grid.addWidget(arrow_size_label, 1, 0)
        arrow_grid.addWidget(arrow_size_edit, 1, 1)
        
        arrow_grid.addWidget(arrow_lw_label, 1, 2)
        arrow_grid.addWidget(arrow_lw_edit, 1, 3)
        
        arrow_grid.addWidget(arrow_direction_label, 1, 4)
        arrow_grid.addWidget(arrow_direction_combo, 1, 5)
        
        arrow_grid.addWidget(arrow_head_length_label, 2, 0)
        arrow_grid.addWidget(arrow_head_length_edit, 2, 1)
        
        arrow_grid.addWidget(arrow_head_width_label, 2, 2)
        arrow_grid.addWidget(arrow_head_width_edit, 2, 3)
        
        arrow_grid.addWidget(arrow_threshold_label, 2, 4)
        arrow_grid.addWidget(arrow_threshold_edit, 2, 5)
        
        arrow_grid.addWidget(arrow_color_real, 3, 1)
        arrow_grid.addWidget(arrow_color_imag, 3, 3)
        
        #--> colorbar properties
        cb_pt_label = QtWidgets.QLabel('PT Colorbar Pad')
        cb_pt_edit = QtWidgets.QLineEdit()
        cb_pt_edit.setText('{0:.2f}'.format(self.cb_pt_pad))
        cb_pt_edit.textChanged[str].connect(self.set_cb_pt_pad)
        
        cb_res_label = QtWidgets.QLabel('Resistivity Colorbar Pad')
        cb_res_edit = QtWidgets.QLineEdit()
        cb_res_edit.setText('{0:.2f}'.format(self.cb_res_pad))
        cb_res_edit.textChanged[str].connect(self.set_cb_res_pad)
        
        res_limits_label = QtWidgets.QLabel('Resistivity Limits (log scale)')
        res_limits_min_edit = QtWidgets.QLineEdit()
        res_limits_min_edit.setText('{0:.1f}'.format(self.res_limits[0]))
        res_limits_min_edit.textChanged[str].connect(self.set_res_limits_min)
        
        res_limits_max_edit = QtWidgets.QLineEdit()
        res_limits_max_edit.setText('{0:.1f}'.format(self.res_limits[1]))
        res_limits_max_edit.textChanged[str].connect(self.set_res_limits_max)
    
        res_grid = QtWidgets.QGridLayout()
        res_grid.addWidget(res_limits_min_edit, 1, 0)        
        res_grid.addWidget(res_limits_max_edit, 1, 1)        
        
        cb_grid = QtWidgets.QGridLayout()
        cb_grid.setSpacing(5)
        
        cb_grid.addWidget(cb_pt_label, 1, 0)
        cb_grid.addWidget(cb_pt_edit, 1, 1)
        
        cb_grid.addWidget(cb_res_label, 1, 2)
        cb_grid.addWidget(cb_res_edit, 1, 3)
        
        cb_grid.addWidget(res_limits_label, 1, 4)
        cb_grid.addLayout(res_grid, 1, 5)
        
        #--> subplot parameters
        subplot_left_label = QtWidgets.QLabel('Subplot Left')
        subplot_left_edit = QtWidgets.QLineEdit()
        subplot_left_edit.setText('{0:.2f}'.format(self.subplot_left))
        subplot_left_edit.textChanged[str].connect(self.set_subplot_left)
        
        subplot_right_label = QtWidgets.QLabel('Subplot right')
        subplot_right_edit = QtWidgets.QLineEdit()
        subplot_right_edit.setText('{0:.2f}'.format(self.subplot_right))
        subplot_right_edit.textChanged[str].connect(self.set_subplot_right)
        
        subplot_bottom_label = QtWidgets.QLabel('Subplot bottom')
        subplot_bottom_edit = QtWidgets.QLineEdit()
        subplot_bottom_edit.setText('{0:.2f}'.format(self.subplot_bottom))
        subplot_bottom_edit.textChanged[str].connect(self.set_subplot_bottom)
        
        subplot_top_label = QtWidgets.QLabel('Subplot top')
        subplot_top_edit = QtWidgets.QLineEdit()
        subplot_top_edit.setText('{0:.2f}'.format(self.subplot_top))
        subplot_top_edit.textChanged[str].connect(self.set_subplot_top)
        
        subplot_hspace_label = QtWidgets.QLabel('Subplot Horizontal Spacing')
        subplot_hspace_edit = QtWidgets.QLineEdit()
        subplot_hspace_edit.setText('{0:.2f}'.format(self.subplot_wspace))
        subplot_hspace_edit.textChanged[str].connect(self.set_subplot_hspace)
        
        subplot_vspace_label = QtWidgets.QLabel('Subplot Vertical Spacing')
        subplot_vspace_edit = QtWidgets.QLineEdit()
        subplot_vspace_edit.setText('{0:.2f}'.format(self.subplot_wspace))
        subplot_vspace_edit.textChanged[str].connect(self.set_subplot_vspace)
        
        subplot_grid = QtWidgets.QGridLayout()
        subplot_grid.setSpacing(5)
        
        subplot_grid.addWidget(subplot_left_label, 1, 0)
        subplot_grid.addWidget(subplot_left_edit, 1, 1)
        
        subplot_grid.addWidget(subplot_right_label, 1, 2)
        subplot_grid.addWidget(subplot_right_edit, 1, 3)
        
        subplot_grid.addWidget(subplot_bottom_label, 1, 4)
        subplot_grid.addWidget(subplot_bottom_edit, 1, 5)
        
        subplot_grid.addWidget(subplot_top_label, 1, 6)
        subplot_grid.addWidget(subplot_top_edit, 1, 7)
        
        subplot_grid.addWidget(subplot_hspace_label, 2, 0)
        subplot_grid.addWidget(subplot_hspace_edit, 2, 1)
        
        subplot_grid.addWidget(subplot_vspace_label, 2, 1)
        subplot_grid.addWidget(subplot_vspace_edit, 2, 2)
        
        #--> update button        
        update_button = QtWidgets.QPushButton('Update')
        update_button.clicked.connect(self.update_settings) 
        
        #--> set the final layout as a vertical box
        vbox = QtWidgets.QVBoxLayout()
        vbox.addLayout(grid_line)
        vbox.addLayout(ellipse_grid)
        vbox.addLayout(arrow_grid)
        vbox.addLayout(cb_grid)
        #vbox.addLayout(subplot_grid)
        vbox.addWidget(update_button)
        
        self.setLayout(vbox) 
        
        self.setGeometry(300, 300, 350, 300)
        self.resize(1050, 500)
        self.setWindowTitle('Plot Settings')    
        self.show()

    def set_text_fs(self, text):
        try:
            self.font_size = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_mapscale(self, text):
        self.map_scale = str(text)
            
    def set_ew_limits_min(self, text):
        try:
            self.ew_limits[0] = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_ew_limits_max(self, text):
        try:
            self.ew_limits[1] = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_ns_limits_min(self, text):
        try:
            self.ns_limits[0] = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_ns_limits_max(self, text):
        try:
            self.ns_limits[1] = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_ellipse_size(self, text):
        try:
            self.ellipse_size = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_ellipse_range_min(self, text):
        try:
            self.ellipse_range[0] = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_ellipse_range_max(self, text):
        try:
            self.ellipse_range[1] = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_ellipse_range_step(self, text):
        try:
            self.ellipse_range[2] = float(text)
        except IndexError:
            self.ellipse_range.append(float(text))
        except ValueError:
            print("Enter a floating point number")
            
    def set_ellipse_cmap(self, text):
        self.ellipse_cmap = str(text)
        
    def set_ellipse_colorby(self, text):
        self.ellipse_colorby = str(text)
        
    def set_arrow_size(self, text):
        try:
            self.arrow_size = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_arrow_lw(self, text):
        try:
            self.arrow_lw = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_arrow_threshold(self, text):
        try:
            self.arrow_threshold = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_arrow_head_length(self, text):
        try:
            self.arrow_head_length = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_arrow_head_width(self, text):
        try:
            self.arrow_head_width = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_arrow_direction(self, text):
        text = str(text)
        
        if text.lower() == 'parkinson':
            self.arrow_direction = 0
        elif text.lower() == 'weise':
            self.arrow_direction = 1
            
    def get_arrow_color_real(self):
        real_color = QtWidgets.QColorDialog().getColor()
        
        if real_color.isValid():
            self.arrow_color_real = real_color.getRgbF()
        else:
            print('Not a valid color')
            
    def get_arrow_color_imag(self):
        imag_color = QtWidgets.QColorDialog().getColor()
        
        if imag_color.isValid():
            self.arrow_color_imag = imag_color.getRgbF()
        else:
            print('Not a valid color')
            
    def set_cb_pt_pad(self, text):
        try:
            self.cb_pt_pad = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_cb_res_pad(self, text):
        try:
            self.cb_res_pad = float(text)
        except ValueError:
            print("Enter a floating point number")
            
    def set_res_limits_min(self, text):
        try:
            self.res_limits[0] = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_res_limits_max(self, text):
        try:
            self.res_limits[1] = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_left(self, text):
        try:
            self.subplot_left = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_right(self, text):
        try:
            self.subplot_right = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_bottom(self, text):
        try:
            self.subplot_bottom = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_top(self, text):
        try:
            self.subplot_top = float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_hspace(self, text):
        try:
            self.subplot_wspace= float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def set_subplot_vspace(self, text):
        try:
            self.subplot_hspace= float(text) 
        except ValueError:
            print("Enter a floating point number")
            
    def update_settings(self):
        self.settings_updated.emit()

#def main():
    
def main():
#if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ui = ModEMPlotPTMap()
    ui.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()
