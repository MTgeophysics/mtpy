# -*- coding: utf-8 -*-
"""
ModEM data and response visualization with a gui.

The user will be able to choose from stations within the data to look at
in either impedance or apparent resistivity and phase.

The functionality is quite simple at the moment

JP 2014
"""
#

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from PyQt4 import QtCore, QtGui
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

import mtpy.modeling.modem as modem
import mtpy.imaging.mtplottools as mtplottools

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8

    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)


class Ui_MainWindow(object):

    def __init__(self):

        self.ms = 5
        self.lw = 1.1
        self.e_capthick = 1
        self.e_capsize = 5

        # color mode
        self.cted = (0, 0, 1)
        self.ctmd = (1, 0, 0)
        self.ctrd = (.5, .5, 0)
        self.ctid = (0, .5, .5)
        self.mted = 's'
        self.mtmd = 'o'

        # color for occam2d model
        self.ctem = (0, .6, .3)
        self.ctmm = (.9, 0, .8)
        self.ctrm = (.9, .4, .2)
        self.ctim = (.2, .4, .9)

        self.mtem = self.mted
        self.mtmm = self.mtmd

        self.res_xx_limits = None
        self.res_xy_limits = None
        self.res_yx_limits = None
        self.res_yy_limits = None

        self.phase_xx_limits = None
        self.phase_xy_limits = None
        self.phase_yx_limits = None
        self.phase_yy_limits = None

        self.subplot_wspace = .25
        self.subplot_hspace = .0
        self.subplot_right = .98
        self.subplot_left = .08
        self.subplot_top = .90
        self.subplot_bottom = .08

        self.legend_loc = 'upper center'
        self.legend_pos = (.5, 1.22)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.fs = 11

        self.plot_component = 4
        self.plot_tipper = True
        self.plot_z = True
        self.ylabel_pad = 1.25
        self.station = None

        self.modem_data = None
        self.modem_resp = None
        self.dirpath = os.getcwd()

    def setupUi(self, MainWindow):
        MainWindow.setWindowTitle("Plot ModEM MT Response")
        MainWindow.resize(1920, 1080)

        # make a central widget that everything is tied to.
        self.central_widget = QtGui.QWidget(MainWindow)
        self.central_widget.setWindowTitle("Plot MT Response")

        # make a widget that will be the station list
        self.list_widget = QtGui.QListWidget()
        self.list_widget.itemClicked.connect(self.get_station)
        self.list_widget.setMaximumWidth(150)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)

        # make sure the figure takes up the entire plottable space
        self.mpl_widget.setSizePolicy(QtGui.QSizePolicy.Expanding,
                                      QtGui.QSizePolicy.Expanding)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.mpl_toolbar = NavigationToolbar(self.mpl_widget, MainWindow)

        # set the layout for the plot
        mpl_vbox = QtGui.QVBoxLayout()
        mpl_vbox.addWidget(self.mpl_toolbar)
        mpl_vbox.addWidget(self.mpl_widget)

        # set the layout the main window
        layout = QtGui.QHBoxLayout()
        layout.addWidget(self.list_widget)
        layout.addLayout(mpl_vbox)
        self.central_widget.setLayout(layout)

        # set the geometry of each widget
        self.list_widget.setObjectName(_fromUtf8("listWidget"))
        self.mpl_widget.setObjectName(_fromUtf8("mpl_widget"))
        self.mpl_widget.updateGeometry()

        # set the central widget
        MainWindow.setCentralWidget(self.central_widget)

        # create a menu bar on the window
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1920, 38))
        self.menubar.setObjectName(_fromUtf8("menubar"))

        # add a tab for File --> open, close, save
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setTitle("Data File")

        self.menuModelFile = QtGui.QMenu(self.menubar)
        self.menuModelFile.setTitle("Response File")

        # add a tab for chaning the display
        self.menuDisplay = QtGui.QMenu(self.menubar)
        self.menuDisplay.setTitle("Display")

        MainWindow.setMenuBar(self.menubar)

        # add a status bar on the bottom of the main window
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))

        MainWindow.setStatusBar(self.statusbar)

        # set an open option that on click opens a modem file
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setText("Open")
        self.actionOpen.triggered.connect(self.get_filename)

        # set a close that closes the main window
        self.actionClose = QtGui.QAction(MainWindow)
        self.actionClose.setText("Close")
        self.actionClose.triggered.connect(MainWindow.close)

        # set a save option that will eventually save the masked data
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setText("Save")

        # add the action on the menu tab
        self.menuFile.addAction(self.actionOpen)
        self.menuFile.addAction(self.actionClose)
        self.menuFile.addAction(self.actionSave)
        self.menubar.addAction(self.menuFile.menuAction())

        self.actionModelOpen = QtGui.QAction(MainWindow)
        self.actionModelOpen.setText("Open")
        self.actionModelOpen.triggered.connect(self.get_resp_fn)
        self.menuModelFile.addAction(self.actionModelOpen)
        self.menubar.addAction(self.menuModelFile.menuAction())
#
        # adding options for display plot type
        self.menu_plot_type = QtGui.QMenu(MainWindow)
        self.menu_plot_type.setTitle("Plot Type")
        self.menuDisplay.addMenu(self.menu_plot_type)
        self.menubar.addAction(self.menuDisplay.menuAction())

        # set plot impedance or resistivity and phase
        self.action_plot_z = QtGui.QAction(MainWindow)
        self.action_plot_z.setText('Impedance')
        self.action_plot_z.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_z)
        self.action_plot_z.toggled.connect(self.status_checked_ptz)

        self.action_plot_rp = QtGui.QAction(MainWindow)
        self.action_plot_rp.setText('Resistivity-Phase')
        self.action_plot_rp.setCheckable(True)
        self.menu_plot_type.addAction(self.action_plot_rp)
        self.action_plot_rp.toggled.connect(self.status_checked_ptrp)

        self.action_plot_settings = QtGui.QAction(MainWindow)
        self.action_plot_settings.setText('Settings')
        self.action_plot_settings.triggered.connect(self.show_settings)
        self.menuDisplay.addAction(self.action_plot_settings)
        self.menubar.addAction(self.menuDisplay.menuAction())

#        self.menuDisplay.addAction(self.menu_plot_style.menuAction())
        self.menuDisplay.addAction(self.menu_plot_type.menuAction())

        # self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def status_checked_ptz(self, toggled):
        """
        be sure that only one plot style is checked
        """

        self.plot_z = toggled
        if toggled == True:
            untoggled = False
        elif toggled == False:
            untoggled = True

        self.action_plot_z.setChecked(toggled)
        self.action_plot_rp.setChecked(untoggled)

        self.plot()

    def status_checked_ptrp(self, toggled):
        """
        be sure that only one plot style is checked
        """

        if toggled == True:
            untoggled = False
            self.plot_z = False
        elif toggled == False:
            untoggled = True
            self.plot_z = True

        self.action_plot_z.setChecked(untoggled)
        self.action_plot_rp.setChecked(toggled)

        self.plot()

    def get_filename(self):
        """
        get the filename from a file dialogue

        """

        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                           filter='(*.dat);; (*.data)'))

        self.modem_data = modem.Data()
        self.modem_data.read_data_file(fn)

        self.dirpath = os.path.dirname(fn)

        print(self.dirpath)

        station_list = sorted(self.modem_data.mt_dict.keys())

        self.list_widget.clear()

        # this will add the station name for each station to the qwidget list
        for station in station_list:
            self.list_widget.addItem(station)

        self.station = str(station_list[0])
        self.plot()

    def get_station(self, widget_item):
        """
        get the station name from the clicked station
        """
        self.station = str(widget_item.text())
        self.plot()

    def get_resp_fn(self):
        """
        get response file name
        """
        print(self.dirpath)
        fn_dialog = QtGui.QFileDialog(directory=self.dirpath)
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM response file',
                                           filter='*.dat'))

        if fn is None:
            return

        self.modem_resp = modem.Data()
        self.modem_resp.read_data_file(fn)
        self.plot()

    def show_settings(self):

        self.settings_window = PlotSettings(None, **self.__dict__)
        self.settings_window.show()
        self.settings_window.settings_updated.connect(self.update_settings)

    def update_settings(self):

        for attr in sorted(self.settings_window.__dict__.keys()):
            setattr(self, attr, self.settings_window.__dict__[attr])

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

        plt.rcParams['font.size'] = self.fs

        fontdict = {'size': self.fs + 2, 'weight': 'bold'}

        #--> make key word dictionaries for plotting
        kw_xx = {'color': self.cted,
                 'marker': self.mted,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        kw_yy = {'color': self.ctmd,
                 'marker': self.mtmd,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        kw_tr = {'color': self.ctrd,
                 'marker': self.mtmd,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        kw_ti = {'color': self.ctid,
                 'marker': self.mtmd,
                 'ms': self.ms,
                 'ls': ':',
                 'lw': self.lw,
                 'e_capsize': self.e_capsize,
                 'e_capthick': self.e_capthick}

        # convert to apparent resistivity and phase
        rp = mtplottools.ResPhase(z_object=z_obj)

        # find locations where points have been masked
        nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
        nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
        nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
        nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
        ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
        nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]

        self.figure.clf()
        self.figure.suptitle(str(self.station), fontdict=fontdict)

        # set the grid of subplots
        gs = gridspec.GridSpec(3, 2,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom,
                               right=self.subplot_right,
                               hspace=self.subplot_hspace,
                               height_ratios=[1, 1, .4])

        # od for off diagonal components, d for diagonal
        axrod = self.figure.add_subplot(gs[0, 0])
        axrd = self.figure.add_subplot(gs[0, 1], sharex=axrod)

        axpod = self.figure.add_subplot(gs[1, 0])
        axpd = self.figure.add_subplot(gs[1, 1], sharex=axrod)

        axtx = self.figure.add_subplot(gs[2, 0], sharex=axrod)
        axty = self.figure.add_subplot(gs[2, 1], sharex=axrod)

        if self.plot_z == False:
            # plot resistivity
            erxx = mtplottools.plot_errorbar(axrd,
                                             period[nzxx],
                                             rp.resxx[nzxx],
                                             rp.resxx_err[nzxx],
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axrod,
                                             period[nzxy],
                                             rp.resxy[nzxy],
                                             rp.resxy_err[nzxy],
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axrod,
                                             period[nzyx],
                                             rp.resyx[nzyx],
                                             rp.resyx_err[nzyx],
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axrd,
                                             period[nzyy],
                                             rp.resyy[nzyy],
                                             rp.resyy_err[nzyy],
                                             **kw_yy)
            # plot phase
            erxx = mtplottools.plot_errorbar(axpd,
                                             period[nzxx],
                                             rp.phasexx[nzxx],
                                             rp.phasexx_err[nzxx],
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axpod,
                                             period[nzxy],
                                             rp.phasexy[nzxy],
                                             rp.phasexy_err[nzxy],
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axpod,
                                             period[nzyx],
                                             rp.phaseyx[nzyx],
                                             rp.phaseyx_err[nzyx],
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axpd,
                                             period[nzyy],
                                             rp.phaseyy[nzyy],
                                             rp.phaseyy_err[nzyy],
                                             **kw_yy)
        elif self.plot_z == True:
            # plot real
            erxx = mtplottools.plot_errorbar(axrd,
                                             period[nzxx],
                                             z_obj.z[nzxx, 0, 0].real,
                                             z_obj.z_err[nzxx, 0, 0].real,
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axrod,
                                             period[nzxy],
                                             z_obj.z[nzxy, 0, 1].real,
                                             z_obj.z_err[nzxy, 0, 1].real,
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axrod,
                                             period[nzyx],
                                             z_obj.z[nzyx, 1, 0].real,
                                             z_obj.z_err[nzyx, 1, 0].real,
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axrd,
                                             period[nzyy],
                                             z_obj.z[nzyy, 1, 1].real,
                                             z_obj.z_err[nzyy, 1, 1].real,
                                             **kw_yy)
            # plot phase
            erxx = mtplottools.plot_errorbar(axpd,
                                             period[nzxx],
                                             z_obj.z[nzxx, 0, 0].imag,
                                             z_obj.z_err[nzxx, 0, 0].real,
                                             **kw_xx)
            erxy = mtplottools.plot_errorbar(axpod,
                                             period[nzxy],
                                             z_obj.z[nzxy, 0, 1].imag,
                                             z_obj.z_err[nzxy, 0, 1].real,
                                             **kw_xx)
            eryx = mtplottools.plot_errorbar(axpod,
                                             period[nzyx],
                                             z_obj.z[nzyx, 1, 0].imag,
                                             z_obj.z_err[nzyx, 1, 0].real,
                                             **kw_yy)
            eryy = mtplottools.plot_errorbar(axpd,
                                             period[nzyy],
                                             z_obj.z[nzyy, 1, 1].imag,
                                             z_obj.z_err[nzyy, 1, 1].real,
                                             **kw_yy)

        # plot tipper
        ertxr = mtplottools.plot_errorbar(axtx,
                                          period[ntx],
                                          t_obj.tipper[ntx, 0, 0].real,
                                          t_obj.tipper_err[ntx, 0, 0],
                                          **kw_tr)
        ertxi = mtplottools.plot_errorbar(axtx,
                                          period[ntx],
                                          t_obj.tipper[ntx, 0, 0].imag,
                                          t_obj.tipper_err[ntx, 0, 1],
                                          **kw_ti)
        ertyr = mtplottools.plot_errorbar(axty,
                                          period[nty],
                                          t_obj.tipper[nty, 0, 1].real,
                                          t_obj.tipper_err[nty, 0, 0],
                                          **kw_tr)

        ertyi = mtplottools.plot_errorbar(axty,
                                          period[nty],
                                          t_obj.tipper[nty, 0, 1].imag,
                                          t_obj.tipper_err[nty, 0, 1],
                                          **kw_ti)
        ax_list = [axrod, axrd,
                   axpod, axpd,
                   axtx, axty]

        od_line_list = [erxy[0], eryx[0], ertxr[0], ertxi[0]]
        d_line_list = [erxx[0], eryy[0], ertyr[0], ertyi[0]]

        od_label_list = ['$Z_{xy}$', '$Z_{yx}$', 'Re{$T_{x}$}', 'Im{$T_{x}$}']
        d_label_list = ['$Z_{xx}$', '$Z_{yy}$', 'Re{$T_{y}$}', 'Im{$T_{y}$}']

        #--> set limits if input
        if self.res_xx_limits is not None:
            axrd.set_ylim(self.res_xx_limits)
        if self.res_xy_limits is not None:
            axrod.set_ylim(self.res_xy_limits)

        if self.phase_xx_limits is not None:
            axpd.set_ylim(self.phase_xx_limits)
        if self.phase_xy_limits is not None:
            axpod.set_ylim(self.phase_xy_limits)

#        axtx.set_ylim((-1.1, 1.1))
#        axty.set_ylim((-1.1, 1.1))

        # set axis properties
        for aa, ax in enumerate(ax_list):
            ax.tick_params(axis='y', pad=self.ylabel_pad)
            ylabels = ax.get_yticks().tolist()
            if aa == 2 or aa == 3 or aa == 4 or aa == 5:
                ylabels[-1] = ''
                ylabels[0] = ''
                ax.set_yticklabels(ylabels)

            if aa < 2:
                if self.plot_z == False:
                    ax.set_yscale('log', nonposy='clip')

            if aa < 4:
                plt.setp(ax.get_xticklabels(), visible=False)

            else:
                ax.set_xlabel('Period (s)', fontdict=fontdict)

            # set axes labels
            if aa == 0:
                if self.plot_z == False:
                    ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                  fontdict=fontdict)
                elif self.plot_z == True:
                    ax.set_ylabel('Re[Z (mV/km nT)]',
                                  fontdict=fontdict)

                ax.legend(od_line_list, od_label_list,
                          loc=self.legend_loc,
                          ncol=2,
                          bbox_to_anchor=self.legend_pos,
                          markerscale=self.legend_marker_scale,
                          borderaxespad=self.legend_border_axes_pad,
                          labelspacing=self.legend_label_spacing,
                          handletextpad=self.legend_handle_text_pad,
                          borderpad=self.legend_border_pad,
                          prop={'size': max([self.fs - 1, 5])})

            elif aa == 1:
                ax.legend(d_line_list, d_label_list,
                          loc=self.legend_loc,
                          ncol=2,
                          bbox_to_anchor=self.legend_pos,
                          markerscale=self.legend_marker_scale,
                          borderaxespad=self.legend_border_axes_pad,
                          labelspacing=self.legend_label_spacing,
                          handletextpad=self.legend_handle_text_pad,
                          borderpad=self.legend_border_pad,
                          prop={'size': max([self.fs - 1, 5])})
            elif aa == 2:
                if self.plot_z == False:
                    ax.set_ylabel('Phase (deg)',
                                  fontdict=fontdict)
                elif self.plot_z == True:
                    ax.set_ylabel('Im[Z (mV/km nT)]',
                                  fontdict=fontdict)

            elif aa == 4:
                ax.set_ylabel(r'W$_x$', fontdict=fontdict)

            elif aa == 5:
                ax.set_ylabel(r'W$_y$', fontdict=fontdict)

            ax.set_xscale('log', nonposx='clip')
            ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))) * 1.01,
                        xmax=10**(np.ceil(np.log10(period[-1]))) * .99)
            ax.grid(True, alpha=.25)

        # plot model response
        if self.modem_resp is not None:
            resp_z_obj = self.modem_resp.mt_dict[self.station].Z
            resp_z_err = np.nan_to_num((z_obj.z - resp_z_obj.z) / z_obj.z_err)

            resp_t_obj = self.modem_resp.mt_dict[self.station].Tipper

            rrp = mtplottools.ResPhase(resp_z_obj)

            rms = resp_z_err.std()
            rms_xx = resp_z_err[:, 0, 0].std()
            rms_xy = resp_z_err[:, 0, 1].std()
            rms_yx = resp_z_err[:, 1, 0].std()
            rms_yy = resp_z_err[:, 1, 1].std()

            #--> make key word dictionaries for plotting
            kw_xx = {'color': self.ctem,
                     'marker': self.mtem,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}

            kw_yy = {'color': self.ctmm,
                     'marker': self.mtmm,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}

            kw_tr = {'color': self.ctrm,
                     'marker': self.mtmd,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}

            kw_ti = {'color': self.ctim,
                     'marker': self.mtmd,
                     'ms': self.ms,
                     'ls': ':',
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}

            if self.plot_z == False:
                # plot resistivity
                rerxx = mtplottools.plot_errorbar(axrd,
                                                  period[nzxx],
                                                  rrp.resxx[nzxx],
                                                  **kw_xx)
                rerxy = mtplottools.plot_errorbar(axrod,
                                                  period[nzxy],
                                                  rrp.resxy[nzxy],
                                                  **kw_xx)
                reryx = mtplottools.plot_errorbar(axrod,
                                                  period[nzyx],
                                                  rrp.resyx[nzyx],
                                                  **kw_yy)
                reryy = mtplottools.plot_errorbar(axrd,
                                                  period[nzyy],
                                                  rrp.resyy[nzyy],
                                                  **kw_yy)
                # plot phase
                rerxx = mtplottools.plot_errorbar(axpd,
                                                  period[nzxx],
                                                  rrp.phasexx[nzxx],
                                                  **kw_xx)
                rerxy = mtplottools.plot_errorbar(axpod,
                                                  period[nzxy],
                                                  rrp.phasexy[nzxy],
                                                  **kw_xx)
                reryx = mtplottools.plot_errorbar(axpod,
                                                  period[nzyx],
                                                  rrp.phaseyx[nzyx],
                                                  **kw_yy)
                reryy = mtplottools.plot_errorbar(axpd,
                                                  period[nzyy],
                                                  rrp.phaseyy[nzyy],
                                                  **kw_yy)
            elif self.plot_z == True:
                # plot real
                rerxx = mtplottools.plot_errorbar(axrd,
                                                  period[nzxx],
                                                  resp_z_obj.z[
                                                      nzxx, 0, 0].real,
                                                  **kw_xx)
                rerxy = mtplottools.plot_errorbar(axrod,
                                                  period[nzxy],
                                                  resp_z_obj.z[
                                                      nzxy, 0, 1].real,
                                                  **kw_xx)
                reryx = mtplottools.plot_errorbar(axrod,
                                                  period[nzyx],
                                                  resp_z_obj.z[
                                                      nzyx, 1, 0].real,
                                                  **kw_yy)
                reryy = mtplottools.plot_errorbar(axrd,
                                                  period[nzyy],
                                                  resp_z_obj.z[
                                                      nzyy, 1, 1].real,
                                                  **kw_yy)
                # plot phase
                rerxx = mtplottools.plot_errorbar(axpd,
                                                  period[nzxx],
                                                  resp_z_obj.z[
                                                      nzxx, 0, 0].imag,
                                                  **kw_xx)
                rerxy = mtplottools.plot_errorbar(axpod,
                                                  period[nzxy],
                                                  resp_z_obj.z[
                                                      nzxy, 0, 1].imag,
                                                  **kw_xx)
                reryx = mtplottools.plot_errorbar(axpod,
                                                  period[nzyx],
                                                  resp_z_obj.z[
                                                      nzyx, 1, 0].imag,
                                                  **kw_yy)
                reryy = mtplottools.plot_errorbar(axpd,
                                                  period[nzyy],
                                                  resp_z_obj.z[
                                                      nzyy, 1, 1].imag,
                                                  **kw_yy)

            rertx = mtplottools.plot_errorbar(axtx,
                                              period[ntx],
                                              resp_t_obj.tipper[
                                                  ntx, 0, 0].real,
                                              **kw_tr)
            rerty = mtplottools.plot_errorbar(axtx,
                                              period[ntx],
                                              resp_t_obj.tipper[
                                                  ntx, 0, 1].imag,
                                              **kw_ti)

            rertx = mtplottools.plot_errorbar(axty,
                                              period[nty],
                                              resp_t_obj.tipper[
                                                  nty, 0, 0].real,
                                              **kw_tr)
            rerty = mtplottools.plot_errorbar(axty,
                                              period[nty],
                                              resp_t_obj.tipper[
                                                  nty, 0, 1].imag,
                                              **kw_ti)

            # add laines and labels
            od_line_list += [rerxy[0], reryx[0], rertx[0], rerty[0]]
            d_line_list += [rerxx[0], reryy[0], rertx[0], rerty[0]]

            od_label_list += ['$Z^m_{xy}$' + 'rms={0:.2f}'.format(rms_xy),
                              '$Z^m_{yx}$' + 'rms={0:.2f}'.format(rms_yx),
                              'Re{$T^m_{x}$}', 'Im{$T^m_{x}$}']
            d_label_list += ['$Z^m_{xx}$' + 'rms={0:.2f}'.format(rms_xx),
                             '$Z^m_{yy}$' + 'rms={0:.2f}'.format(rms_yy),
                             'Re{$T^m_{x}$}', 'Im{$T^m_{x}$}']

            legend_ax_list = [axrod, axrd]
            od_label_list = od_label_list[0:2] + od_label_list[4:6] +\
                od_label_list[2:4] + od_label_list[6:8]
            d_label_list = d_label_list[0:2] + d_label_list[4:6] +\
                d_label_list[2:4] + d_label_list[6:8]

            od_line_list = od_line_list[0:2] + od_line_list[4:6] +\
                od_line_list[2:4] + od_line_list[6:8]
            d_line_list = d_line_list[0:2] + d_line_list[4:6] +\
                d_line_list[2:4] + d_line_list[6:8]

            legend_ll = [od_line_list, d_line_list]
            legend_tl = [od_label_list, d_label_list]

            for ax, ll, tl in zip(legend_ax_list, legend_ll, legend_tl):
                ax.legend(ll,
                          tl,
                          loc=self.legend_loc,
                          ncol=4,
                          bbox_to_anchor=self.legend_pos,
                          markerscale=self.legend_marker_scale,
                          borderaxespad=self.legend_border_axes_pad,
                          labelspacing=self.legend_label_spacing,
                          handletextpad=self.legend_handle_text_pad,
                          borderpad=self.legend_border_pad,
                          prop={'size': max([self.fs - 1, 5])})

        self.mpl_widget.draw()


class PlotSettings(QtGui.QWidget):
    settings_updated = QtCore.pyqtSignal()

    def __init__(self, parent, **kwargs):
        super(PlotSettings, self).__init__(parent)

        self.fs = kwargs.pop('fs', 10)
        self.lw = kwargs.pop('lw', 1.5)
        self.ms = kwargs.pop('ms', 5)

        self.e_capthick = kwargs.pop('e_capthick', 1)
        self.e_capsize = kwargs.pop('e_capsize', 5)

        # color mode
        self.cted = kwargs.pop('cted', (0, 0, 1))
        self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
        self.ctrd = kwargs.pop('ctrd', (.5, .5, 0))
        self.ctid = kwargs.pop('ctid', (0, .5, .5))
        self.mted = kwargs.pop('mted', 's')
        self.mtmd = kwargs.pop('mtmd', 'o')

        # color for occam2d model
        self.ctem = kwargs.pop('ctem', (0, .6, .3))
        self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
        self.ctrm = kwargs.pop('ctrm', (.5, .5, .1))
        self.ctim = kwargs.pop('ctim', (.1, .5, .5))
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

        self.subplot_wspace = kwargs.pop('subplot_wspace', .2)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .93)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .08)

        self.legend_loc = kwargs.pop('legend_loc', 'upper center')
        self.legend_pos = kwargs.pop('legend_pos', (.5, 1.11))
        self.legend_marker_scale = kwargs.pop('legend_marker_scale', 1)
        self.legend_border_axes_pad = kwargs.pop('legend_border_axes_pad', .01)
        self.legend_label_spacing = kwargs.pop('legend_label_spacing', 0.07)
        self.legend_handle_text_pad = kwargs.pop('legend_handle_text_pad', .2)
        self.legend_border_pad = kwargs.pop('legend_border_pad', .15)

        self.initUI()

    def initUI(self):
        #--> line properties
        fs_label = QtGui.QLabel('Font Size')
        fs_edit = QtGui.QLineEdit()
        fs_edit.setText('{0:.1f}'.format(self.fs))
        fs_edit.textChanged[str].connect(self.set_text_fs)

        lw_label = QtGui.QLabel('Line Width')
        lw_edit = QtGui.QLineEdit()
        lw_edit.setText('{0:.1f}'.format(self.lw))
        lw_edit.textChanged[str].connect(self.set_text_lw)

        e_capthick_label = QtGui.QLabel('Error cap thickness')
        e_capthick_edit = QtGui.QLineEdit()
        e_capthick_edit.setText('{0:.1f}'.format(self.e_capthick))
        e_capthick_edit.textChanged[str].connect(self.set_text_e_capthick)

        e_capsize_label = QtGui.QLabel('Error cap size')
        e_capsize_edit = QtGui.QLineEdit()
        e_capsize_edit.setText('{0:.1f}'.format(self.e_capsize))
        e_capsize_edit.textChanged[str].connect(self.set_text_e_capsize)

        grid_line = QtGui.QGridLayout()
        grid_line.setSpacing(10)

        grid_line.addWidget(fs_label, 1, 0)
        grid_line.addWidget(fs_edit, 1, 1)

        grid_line.addWidget(lw_label, 1, 2)
        grid_line.addWidget(lw_edit, 1, 3)

        grid_line.addWidget(e_capthick_label, 1, 4)
        grid_line.addWidget(e_capthick_edit, 1, 5)

        grid_line.addWidget(e_capsize_label, 1, 6)
        grid_line.addWidget(e_capsize_edit, 1, 7)

        #--> marker properties
        ms_label = QtGui.QLabel('Marker Size')
        ms_edit = QtGui.QLineEdit()
        ms_edit.setText('{0:.1f}'.format(self.ms))
        ms_edit.textChanged[str].connect(self.set_text_ms)

        dcxy_button = QtGui.QPushButton('Data Color xy', self)
        dcxy_button.clicked.connect(self.set_cted)

        dcyx_button = QtGui.QPushButton('Data Color yx', self)
        dcyx_button.clicked.connect(self.set_ctmd)

        dctr_button = QtGui.QPushButton('Data Color Re{T}', self)
        dctr_button.clicked.connect(self.set_ctrd)

        dcti_button = QtGui.QPushButton('Data Color Im{T}', self)
        dcti_button.clicked.connect(self.set_ctid)

        mcxy_button = QtGui.QPushButton('Model Color xy', self)
        mcxy_button.clicked.connect(self.set_ctem)

        mcyx_button = QtGui.QPushButton('Model Color yx', self)
        mcyx_button.clicked.connect(self.set_ctmm)

        mctr_button = QtGui.QPushButton('Model Color Re{T}', self)
        mctr_button.clicked.connect(self.set_ctrm)

        mcti_button = QtGui.QPushButton('Model Color Im{T}', self)
        mcti_button.clicked.connect(self.set_ctim)

        dmxy_label = QtGui.QLabel('Data Marker xy')
        dmxy_edit = QtGui.QLineEdit()
        dmxy_edit.setText('{0}'.format(self.mted))
        dmxy_edit.textChanged[str].connect(self.set_text_mted)

        dmyx_label = QtGui.QLabel('Data Marker yx')
        dmyx_edit = QtGui.QLineEdit()
        dmyx_edit.setText('{0}'.format(self.mtmd))
        dmyx_edit.textChanged[str].connect(self.set_text_mtmd)

        mmxy_label = QtGui.QLabel('Model Marker xy')
        mmxy_edit = QtGui.QLineEdit()
        mmxy_edit.setText('{0}'.format(self.mtem))
        mmxy_edit.textChanged[str].connect(self.set_text_mtem)

        mmyx_label = QtGui.QLabel('Model Marker yx')
        mmyx_edit = QtGui.QLineEdit()
        mmyx_edit.setText('{0}'.format(self.mtmm))
        mmyx_edit.textChanged[str].connect(self.set_text_mtmm)

        marker_label = QtGui.QLabel('Maker Properties:')

        marker_grid = QtGui.QGridLayout()
        marker_grid.setSpacing(10)

        marker_grid.addWidget(marker_label, 1, 0)
        marker_grid.addWidget(ms_label, 1, 2)
        marker_grid.addWidget(ms_edit, 1, 3)

        marker_grid.addWidget(dcxy_button, 2, 0)
        marker_grid.addWidget(dcyx_button, 2, 1)
        marker_grid.addWidget(dctr_button, 2, 2)
        marker_grid.addWidget(dcti_button, 2, 3)

        marker_grid.addWidget(dmxy_label, 2, 4)
        marker_grid.addWidget(dmxy_edit, 2, 5)

        marker_grid.addWidget(dmyx_label, 2, 6)
        marker_grid.addWidget(dmyx_edit, 2, 7)

        marker_grid.addWidget(mcxy_button, 3, 0)
        marker_grid.addWidget(mcyx_button, 3, 1)
        marker_grid.addWidget(mctr_button, 3, 2)
        marker_grid.addWidget(mcti_button, 3, 3)

        marker_grid.addWidget(mmxy_label, 3, 4)
        marker_grid.addWidget(mmxy_edit, 3, 5)

        marker_grid.addWidget(mmyx_label, 3, 6)
        marker_grid.addWidget(mmyx_edit, 3, 7)

        #--> plot limits
        ylimr_xx_label = QtGui.QLabel('Res_xx')
        ylimr_xx_edit = QtGui.QLineEdit()
        ylimr_xx_edit.setText('{0}'.format(self.res_xx_limits))
        ylimr_xx_edit.textChanged[str].connect(self.set_text_res_xx)

        ylimr_xy_label = QtGui.QLabel('Res_xy')
        ylimr_xy_edit = QtGui.QLineEdit()
        ylimr_xy_edit.setText('{0}'.format(self.res_xy_limits))
        ylimr_xy_edit.textChanged[str].connect(self.set_text_res_xy)

        ylimr_yx_label = QtGui.QLabel('Res_yx')
        ylimr_yx_edit = QtGui.QLineEdit()
        ylimr_yx_edit.setText('{0}'.format(self.res_yx_limits))
        ylimr_yx_edit.textChanged[str].connect(self.set_text_res_yx)

        ylimr_yy_label = QtGui.QLabel('Res_yy')
        ylimr_yy_edit = QtGui.QLineEdit()
        ylimr_yy_edit.setText('{0}'.format(self.res_yy_limits))
        ylimr_yy_edit.textChanged[str].connect(self.set_text_res_yy)

        ylimp_xx_label = QtGui.QLabel('phase_xx')
        ylimp_xx_edit = QtGui.QLineEdit()
        ylimp_xx_edit.setText('{0}'.format(self.phase_xx_limits))
        ylimp_xx_edit.textChanged[str].connect(self.set_text_phase_xx)

        ylimp_xy_label = QtGui.QLabel('phase_xy')
        ylimp_xy_edit = QtGui.QLineEdit()
        ylimp_xy_edit.setText('{0}'.format(self.phase_xy_limits))
        ylimp_xy_edit.textChanged[str].connect(self.set_text_phase_xy)

        ylimp_yx_label = QtGui.QLabel('phase_yx')
        ylimp_yx_edit = QtGui.QLineEdit()
        ylimp_yx_edit.setText('{0}'.format(self.phase_yx_limits))
        ylimp_yx_edit.textChanged[str].connect(self.set_text_phase_yx)

        ylimp_yy_label = QtGui.QLabel('phase_yy')
        ylimp_yy_edit = QtGui.QLineEdit()
        ylimp_yy_edit.setText('{0}'.format(self.phase_yy_limits))
        ylimp_yy_edit.textChanged[str].connect(self.set_text_phase_yy)

        limits_grid = QtGui.QGridLayout()
        limits_grid.setSpacing(10)

        limits_label = QtGui.QLabel('Plot Limits: (Res=Real, Phase=Imaginary)'
                                    ' --> input on a linear scale')

        limits_grid.addWidget(limits_label, 1, 0, 1, 7)

        limits_grid.addWidget(ylimr_xx_label, 2, 0)
        limits_grid.addWidget(ylimr_xx_edit, 2, 1)
        limits_grid.addWidget(ylimr_xy_label, 2, 2)
        limits_grid.addWidget(ylimr_xy_edit, 2, 3)
        limits_grid.addWidget(ylimr_yx_label, 2, 4)
        limits_grid.addWidget(ylimr_yx_edit, 2, 5)
        limits_grid.addWidget(ylimr_yy_label, 2, 6)
        limits_grid.addWidget(ylimr_yy_edit, 2, 7)

        limits_grid.addWidget(ylimp_xx_label, 3, 0)
        limits_grid.addWidget(ylimp_xx_edit, 3, 1)
        limits_grid.addWidget(ylimp_xy_label, 3, 2)
        limits_grid.addWidget(ylimp_xy_edit, 3, 3)
        limits_grid.addWidget(ylimp_yx_label, 3, 4)
        limits_grid.addWidget(ylimp_yx_edit, 3, 5)
        limits_grid.addWidget(ylimp_yy_label, 3, 6)
        limits_grid.addWidget(ylimp_yy_edit, 3, 7)

        #--> legend properties
        legend_pos_label = QtGui.QLabel('Legend Position')
        legend_pos_edit = QtGui.QLineEdit()
        legend_pos_edit.setText('{0}'.format(self.legend_pos))
        legend_pos_edit.textChanged[str].connect(self.set_text_legend_pos)

        legend_grid = QtGui.QGridLayout()
        legend_grid.setSpacing(10)

        legend_grid.addWidget(QtGui.QLabel('Legend Properties:'), 1, 0)
        legend_grid.addWidget(legend_pos_label, 1, 2,)
        legend_grid.addWidget(legend_pos_edit, 1, 3)

        update_button = QtGui.QPushButton('Update')
        update_button.clicked.connect(self.update_settings)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid_line)
        vbox.addLayout(marker_grid)
        vbox.addLayout(limits_grid)
        vbox.addLayout(legend_grid)
        vbox.addWidget(update_button)

        self.setLayout(vbox)

        self.setGeometry(300, 300, 350, 300)
        self.resize(1350, 500)
        self.setWindowTitle('Plot Settings')
        self.show()

    def set_text_fs(self, text):
        try:
            self.fs = float(text)
        except ValueError:
            print("Enter a float point number")

    def set_text_e_capthick(self, text):
        try:
            self.e_capthick = float(text)
        except ValueError:
            print("Enter a float point number")

    def set_text_e_capsize(self, text):
        try:
            self.e_capsize = float(text)
        except ValueError:
            print("Enter a float point number")

    def set_text_lw(self, text):
        try:
            self.lw = float(text)
        except ValueError:
            print("Enter a float point number")

    def set_text_ms(self, text):
        try:
            self.ms = float(text)
        except ValueError:
            print("Enter a float point number")

    def set_cted(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.cted = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctmd(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctmd = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctrd(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctrd = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctid(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctid = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_text_mted(self, text):
        try:
            self.mted = str(text)
        except ValueError:
            print("Enter a string")

    def set_text_mtmd(self, text):
        try:
            self.mtmd = str(text)
        except ValueError:
            print("Enter a string")

    def set_ctem(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctem = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctmm(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctmm = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctrm(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctrm = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_ctim(self):
        imag_color = QtGui.QColorDialog().getColor()

        if imag_color.isValid():
            self.ctim = imag_color.getRgbF()
        else:
            print('Not a valid color')

    def set_text_mtem(self, text):
        try:
            self.mtem = str(text)
        except ValueError:
            print("Enter a string")

    def set_text_mtmm(self, text):
        try:
            self.mtmm = str(text)
        except ValueError:
            print("Enter a string")

    def set_text_res_xx(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_xx_limits = tuple(l_list)

    def set_text_res_xy(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_xy_limits = tuple(l_list)

    def set_text_res_yx(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_yx_limits = tuple(l_list)

    def set_text_res_yy(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.res_yy_limits = tuple(l_list)

    def set_text_phase_xx(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_xx_limits = tuple(l_list)

    def set_text_phase_xy(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_xy_limits = tuple(l_list)

    def set_text_phase_yx(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_yx_limits = tuple(l_list)

    def set_text_phase_yy(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.phase_yy_limits = tuple(l_list)

    def set_text_legend_pos(self, text):
        if text == 'None':
            return
        text = text.replace('(', '').replace(')', '')
        t_list = text.split(',')
        if len(t_list) != 2:
            print('enter as (min, max)')
        l_list = []
        for txt in t_list:
            try:
                l_list.append(float(txt))
            except ValueError:
                pass
        if len(l_list) == 2:
            self.legend_pos = tuple(l_list)

    def update_settings(self):
        self.settings_updated.emit()

# def main():


def main():
    # if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()
