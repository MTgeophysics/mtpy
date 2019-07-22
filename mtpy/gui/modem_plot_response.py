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
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
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

        self.ms  = 5
        self.lw = 1.5
        self.e_capthick = 1
        self.e_capsize =  5

        #color mode
        self.cted = (0, 0, .75)
        self.ctmd = (.75, 0, 0)
        self.mted = 's'
        self.mtmd = 'o'

        #color for occam2d model
        self.ctem = (0, .6, .3)
        self.ctmm = (.9, 0, .8)
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
        self.tipper_limits = None


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

        self.z_err_increase = .05
        self.t_err_increase = .05

        self.fs = 11

        self.plot_component = 4
        self.plot_tipper = True
        self.plot_z = False
        self.ylabel_pad = 1.25
        self.station = None

        self.modem_data = None
        self.modem_resp = None
        self.dirpath = os.getcwd()

    def setupUi(self, MainWindow):
        MainWindow.setWindowTitle("Plot ModEM MT Response")
        MainWindow.setWindowState(QtCore.Qt.WindowMaximized)

        #make a central widget that everything is tied to.
        self.central_widget = QtGui.QWidget(MainWindow)
        self.central_widget.setWindowTitle("Plot MT Response")

        #make a widget that will be the station list
        self.list_widget = QtGui.QListWidget()
        self.list_widget.itemClicked.connect(self.get_station)
        self.list_widget.setMaximumWidth(150)

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.figure = Figure(dpi=150)
        self.mpl_widget = FigureCanvas(self.figure)
        self.mpl_widget.setParent(self.central_widget)
        self.mpl_widget.setFocusPolicy(QtCore.Qt.ClickFocus)
        self.mpl_widget.setFocus()

        # be able to edit the data
        self.mpl_widget.mpl_connect('pick_event', self.on_pick)
        self.mpl_widget.mpl_connect('axes_enter_event', self.in_axes)

        #make sure the figure takes up the entire plottable space
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

        #set the geometry of each widget
        self.list_widget.setObjectName(_fromUtf8("listWidget"))
        self.mpl_widget.setObjectName(_fromUtf8("mpl_widget"))
        self.mpl_widget.updateGeometry()

        #set the central widget
        MainWindow.setCentralWidget(self.central_widget)

        #create a menu bar on the window
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

        # add a tab for help
        self.menuHelp = QtGui.QMenu(self.menubar)
        self.menuHelp.setTitle("Help")

        MainWindow.setMenuBar(self.menubar)

        # add a status bar on the bottom of the main window
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))

        MainWindow.setStatusBar(self.statusbar)

        # set an open option that on click opens a modem file
        self.actionOpen = QtGui.QAction(MainWindow)
        self.actionOpen.setText("Open")
        self.actionOpen.triggered.connect(self.get_data_file)

        # set a close that closes the main window
        self.actionClose = QtGui.QAction(MainWindow)
        self.actionClose.setText("Close")
        self.actionClose.triggered.connect(MainWindow.close)

        # set a save option that will eventually save the masked data
        self.actionSave = QtGui.QAction(MainWindow)
        self.actionSave.setText("Save Edits")
        self.actionSave.triggered.connect(self.save_edits)

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
        #adding options for display plot type
        self.menu_plot_type = QtGui.QMenu(MainWindow)
        self.menu_plot_type.setTitle("Plot Type")
        self.menuDisplay.addMenu(self.menu_plot_type)
        self.menubar.addAction(self.menuDisplay.menuAction())

        #set plot impedance or resistivity and phase
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

#        self.action_edit_settings = QtGui.QAction(MainWindow)
#        self.action_edit_settings.setText('Editing Settings')
#        self.action_edit_settings.triggered.connect(self.show_edit_settings)


#        self.menuDisplay.addAction(self.menu_plot_style.menuAction())
        self.menuDisplay.addAction(self.menu_plot_type.menuAction())

        self.action_help = QtGui.QAction(MainWindow)
        self.action_help.setText('help doc')
        self.action_help.triggered.connect(self.disp_help)
        self.menuHelp.addAction(self.action_help)
        self.menubar.addAction(self.menuHelp.menuAction())


        #self.retranslateUi(MainWindow)
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


    def get_data_file(self):
        """
        get the filename from a file dialogue

        """

        fn_dialog = QtGui.QFileDialog()
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM data file',
                                       filter='(*.dat);; (*.data)'))

        self.modem_data = modem.Data()
        self.modem_data.read_data_file(fn)

        self.dirpath = os.path.dirname(fn)

        station_list = sorted(self.modem_data.mt_dict.keys())

        self.list_widget.clear()

        #this will add the station name for each station to the qwidget list
        for station in station_list:
#            try:
#                self.list_widget.addItem('{0:.0f}'.format(station))
#            except ValueError:
            self.list_widget.addItem(str(station))

        self.station = station_list[0]
        self.plot()

    def save_edits(self):
        """
        save edits to another file
        """
        fn_dialog = QtGui.QFileDialog()
        save_fn = str(fn_dialog.getSaveFileName(caption='Choose File to save',
                                                filter='*.dat'))

        self.modem_data.write_data_file(save_path=os.path.dirname(save_fn),
                                        fn_basename=os.path.basename(save_fn),
                                        compute_error=False,
                                        fill=False)

    def get_station(self, widget_item):
        """
        get the station name from the clicked station
        """
        self.station = str(widget_item.text())
        try:
            self.station = float(self.station)
        except ValueError:
            pass
        self.plot()

    def get_resp_fn(self):
        """
        get response file name
        """

        fn_dialog = QtGui.QFileDialog(directory=self.dirpath)
        fn = str(fn_dialog.getOpenFileName(caption='Choose ModEM response file',
                                       filter='*.dat'))

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

        help_popup = QtGui.QMessageBox.information(self.central_widget, 'Help',
                                                   help_string)

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

        fontdict = {'size':self.fs+2, 'weight':'bold'}
        if self.plot_z == True:
            h_ratio = [1,1]
        elif self.plot_z == False:
            h_ratio = [2, 1.5]

        #--> make key word dictionaries for plotting
        kw_xx = {'color':self.cted,
                 'marker':self.mted,
                 'ms':self.ms,
                 'ls':':',
                 'lw':self.lw,
                 'e_capsize':self.e_capsize,
                 'e_capthick':self.e_capthick,
                 'picker':3}

        kw_yy = {'color':self.ctmd,
                 'marker':self.mtmd,
                 'ms':self.ms,
                 'ls':':',
                 'lw':self.lw,
                 'e_capsize':self.e_capsize,
                 'e_capthick':self.e_capthick,
                 'picker':3}



        #convert to apparent resistivity and phase
        rp = mtplottools.ResPhase(z_object=z_obj)

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

        if self.plot_tipper == True:
            gs = gridspec.GridSpec(2, 6,
                               wspace=self.subplot_wspace,
                               left=self.subplot_left,
                               top=self.subplot_top,
                               bottom=self.subplot_bottom,
                               right=self.subplot_right,
                               hspace=self.subplot_hspace,
                               height_ratios=h_ratio)
        else:
            gs = gridspec.GridSpec(2, 4,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio)

        if self.plot_tipper == False:
            axrxx = self.figure.add_subplot(gs[0, 0])
            axrxy = self.figure.add_subplot(gs[0, 1], sharex=axrxx)
            axryx = self.figure.add_subplot(gs[0, 2], sharex=axrxx)
            axryy = self.figure.add_subplot(gs[0, 3], sharex=axrxx)

            axpxx = self.figure.add_subplot(gs[1, 0])
            axpxy = self.figure.add_subplot(gs[1, 1], sharex=axrxx)
            axpyx = self.figure.add_subplot(gs[1, 2], sharex=axrxx)
            axpyy = self.figure.add_subplot(gs[1, 3], sharex=axrxx)

            self.ax_list = [axrxx, axrxy, axryx, axryy,
                            axpxx, axpxy, axpyx, axpyy]
        else:
            axrxx = self.figure.add_subplot(gs[0, 0])
            axrxy = self.figure.add_subplot(gs[0, 1], sharex=axrxx)
            axryx = self.figure.add_subplot(gs[0, 2], sharex=axrxx)
            axryy = self.figure.add_subplot(gs[0, 3], sharex=axrxx)

            axpxx = self.figure.add_subplot(gs[1, 0])
            axpxy = self.figure.add_subplot(gs[1, 1], sharex=axrxx)
            axpyx = self.figure.add_subplot(gs[1, 2], sharex=axrxx)
            axpyy = self.figure.add_subplot(gs[1, 3], sharex=axrxx)

            axtxr = self.figure.add_subplot(gs[0, 4], sharex=axrxx)
            axtxi = self.figure.add_subplot(gs[1, 4], sharex=axrxx)
            axtyr = self.figure.add_subplot(gs[0, 5], sharex=axrxx)
            axtyi = self.figure.add_subplot(gs[1, 5], sharex=axrxx)

            self.ax_list = [axrxx, axrxy, axryx, axryy,
                            axpxx, axpxy, axpyx, axpyy,
                            axtxr, axtxi, axtyr, axtyi]

        if self.plot_z == False:
            #plot resistivity
            erxx = mtplottools.plot_errorbar(axrxx,
                                      period[nzxx],
                                      rp.resxx[nzxx],
                                      rp.resxx_err[nzxx],
                                      **kw_xx)
            erxy = mtplottools.plot_errorbar(axrxy,
                                      period[nzxy],
                                      rp.resxy[nzxy],
                                      rp.resxy_err[nzxy],
                                      **kw_xx)
            eryx = mtplottools.plot_errorbar(axryx,
                                      period[nzyx],
                                      rp.resyx[nzyx],
                                      rp.resyx_err[nzyx],
                                      **kw_yy)
            eryy = mtplottools.plot_errorbar(axryy,
                                      period[nzyy],
                                      rp.resyy[nzyy],
                                      rp.resyy_err[nzyy],
                                      **kw_yy)
            #plot phase
            epxx = mtplottools.plot_errorbar(axpxx,
                                      period[nzxx],
                                      rp.phasexx[nzxx],
                                      rp.phasexx_err[nzxx],
                                      **kw_xx)
            epxy = mtplottools.plot_errorbar(axpxy,
                                      period[nzxy],
                                      rp.phasexy[nzxy],
                                      rp.phasexy_err[nzxy],
                                      **kw_xx)
            epyx = mtplottools.plot_errorbar(axpyx,
                                      period[nzyx],
                                      rp.phaseyx[nzyx],
                                      rp.phaseyx_err[nzyx],
                                      **kw_yy)
            epyy = mtplottools.plot_errorbar(axpyy,
                                      period[nzyy],
                                      rp.phaseyy[nzyy],
                                      rp.phaseyy_err[nzyy],
                                      **kw_yy)
        elif self.plot_z == True:
            #plot real
            erxx = mtplottools.plot_errorbar(axrxx,
                                      period[nzxx],
                                      z_obj.z[nzxx,0,0].real,
                                      z_obj.z_err[nzxx,0,0].real,
                                      **kw_xx)
            erxy = mtplottools.plot_errorbar(axrxy,
                                      period[nzxy],
                                      z_obj.z[nzxy,0,1].real,
                                      z_obj.z_err[nzxy,0,1].real,
                                      **kw_xx)
            eryx = mtplottools.plot_errorbar(axryx,
                                      period[nzyx],
                                      z_obj.z[nzyx,1,0].real,
                                      z_obj.z_err[nzyx,1,0].real,
                                      **kw_yy)
            eryy = mtplottools.plot_errorbar(axryy,
                                      period[nzyy],
                                      z_obj.z[nzyy,1,1].real,
                                      z_obj.z_err[nzyy,1,1].real,
                                      **kw_yy)
            #plot phase
            epxx = mtplottools.plot_errorbar(axpxx,
                                      period[nzxx],
                                      z_obj.z[nzxx,0,0].imag,
                                      z_obj.z_err[nzxx,0,0].imag,
                                      **kw_xx)
            epxy = mtplottools.plot_errorbar(axpxy,
                                      period[nzxy],
                                      z_obj.z[nzxy,0,1].imag,
                                      z_obj.z_err[nzxy,0,1].imag,
                                      **kw_xx)
            epyx = mtplottools.plot_errorbar(axpyx,
                                      period[nzyx],
                                      z_obj.z[nzyx,1,0].imag,
                                      z_obj.z_err[nzyx,1,0].imag,
                                      **kw_yy)
            epyy = mtplottools.plot_errorbar(axpyy,
                                      period[nzyy],
                                      z_obj.z[nzyy,1,1].imag,
                                      z_obj.z_err[nzyy,1,1].imag,
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
                                     t_obj.tipper_err[nty, 0, 0],
                                     **kw_yy)

            eptx = mtplottools.plot_errorbar(axtxi,
                                     period[ntx],
                                     t_obj.tipper[ntx, 0, 0].imag,
                                     t_obj.tipper_err[ntx, 0, 1],
                                     **kw_xx)
            epty = mtplottools.plot_errorbar(axtyi,
                                     period[nty],
                                     t_obj.tipper[nty, 0, 1].imag,
                                     t_obj.tipper_err[nty, 0, 1],
                                     **kw_yy)
        if self.plot_tipper == False:
            line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
            label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                          ['$Z_{yx}$'], ['$Z_{yy}$']]
            for ax, label in zip(self.ax_list, label_list):
                ax.set_title(label[0],fontdict={'size':self.fs+2,
                                              'weight':'bold'})
            self._err_list = [[erxx[1][0],erxx[1][1],erxx[2][0]],
                              [erxy[1][0],erxy[1][1],erxy[2][0]],
                              [eryx[1][0],eryx[1][1],eryx[2][0]],
                              [eryy[1][0],eryy[1][1],eryy[2][0]]]
        else:
            line_list = [[erxx[0]], [erxy[0]],
                         [eryx[0]], [eryy[0]],
                         [ertx[0]], [erty[0]]]
            label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                          ['$Z_{yx}$'], ['$Z_{yy}$'],
                          ['$T_{x}$'], ['$T_{y}$']]
            for ax, label in zip([axrxx, axrxy, axryx, axryy, axtxr, axtyr],
                                 label_list):
                ax.set_title(label[0], fontdict={'size':self.fs+2,
                                              'weight':'bold'})

            self._err_list = [[erxx[1][0],erxx[1][1],erxx[2][0]],
                              [erxy[1][0],erxy[1][1],erxy[2][0]],
                              [eryx[1][0],eryx[1][1],eryx[2][0]],
                              [eryy[1][0],eryy[1][1],eryy[2][0]],
                              [ertx[1][0],ertx[1][1],ertx[2][0]],
                              [erty[1][0],erty[1][1],erty[2][0]]]

        #--> set limits if input
        if self.res_xx_limits is not None:
            axrxx.set_ylim(self.res_xx_limits)
        if self.res_xy_limits is not None:
            axrxy.set_ylim(self.res_xy_limits)
        if self.res_yx_limits is not None:
            axryx.set_ylim(self.res_yx_limits)
        if self.res_yy_limits is not None:
            axryy.set_ylim(self.res_yy_limits)

        if self.phase_xx_limits is not None:
            axpxx.set_ylim(self.phase_xx_limits)
        if self.phase_xy_limits is not None:
            axpxy.set_ylim(self.phase_xy_limits)
        if self.phase_yx_limits is not None:
            axpyx.set_ylim(self.phase_yx_limits)
        if self.phase_yy_limits is not None:
            axpyy.set_ylim(self.phase_yy_limits)

        #set axis properties
        for aa, ax in enumerate(self.ax_list):
            ax.tick_params(axis='y', pad=self.ylabel_pad)
            ylabels = ax.get_yticks().tolist()
            if aa < 8:
                ylabels[-1] = ''
                ylabels[0] = ''
                ax.set_yticklabels(ylabels)

            if len(self.ax_list) == 4 or len(self.ax_list) == 6:
                if aa < 2:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    if self.plot_z == False:
                        ax.set_yscale('log', nonposy='clip')
                else:
                    ax.set_xlabel('Period (s)', fontdict=fontdict)

                #set axes labels
                if aa == 0:
                    if self.plot_z == False:
                        ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                      fontdict=fontdict)
                    elif self.plot_z == True:
                        ax.set_ylabel('Re[Z (mV/km nT)]',
                                      fontdict=fontdict)
                elif aa == 2:
                    if self.plot_z == False:
                        ax.set_ylabel('Phase (deg)',
                                      fontdict=fontdict)
                    elif self.plot_z == True:
                        ax.set_ylabel('Im[Z (mV/km nT)]',
                                      fontdict=fontdict)

            elif len(self.ax_list) == 8 or len(self.ax_list) == 12:
                if aa < 4:
                    plt.setp(ax.get_xticklabels(), visible=False)
                    if self.plot_z == False:
                        ax.set_yscale('log', nonposy='clip')
                else:
                    if aa == 8 or aa == 10:
                        plt.setp(ax.get_xticklabels(), visible=False)
                    else:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)

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
            if aa > 7:
                if self.tipper_limits is not None:
                    ax.set_ylim(self.tipper_limits)
                else:
                    ax.set_ylim(-1.1, 1.1)

            ax.set_xscale('log', nonposx='clip')
            ax.set_xlim(xmin=10**(np.floor(np.log10(period[0]))) * 1.01,
                        xmax=10**(np.ceil(np.log10(period[-1]))) * .99)
            ax.grid(True, alpha=.25)

        #plot model response
        if self.modem_resp is not None:
            resp_z_obj = self.modem_resp.mt_dict[self.station].Z
            resp_z_err = np.nan_to_num((z_obj.z-resp_z_obj.z)/z_obj.z_err)

            resp_t_obj = self.modem_resp.mt_dict[self.station].Tipper

            rrp = mtplottools.ResPhase(resp_z_obj)

            rms_xx = resp_z_err[:, 0, 0].std()
            rms_xy = resp_z_err[:, 0, 1].std()
            rms_yx = resp_z_err[:, 1, 0].std()
            rms_yy = resp_z_err[:, 1, 1].std()

            #--> make key word dictionaries for plotting
            kw_xx = {'color':self.ctem,
                     'marker':self.mtem,
                     'ms':self.ms,
                     'ls':':',
                     'lw':self.lw,
                     'e_capsize':self.e_capsize,
                     'e_capthick':self.e_capthick}

            kw_yy = {'color':self.ctmm,
                     'marker':self.mtmm,
                     'ms':self.ms,
                     'ls':':',
                     'lw':self.lw,
                     'e_capsize':self.e_capsize,
                     'e_capthick':self.e_capthick}
            if self.plot_z == False:
                #plot resistivity
                rerxx= mtplottools.plot_errorbar(axrxx,
                                          period[nzxx],
                                          rrp.resxx[nzxx],
                                          **kw_xx)
                rerxy = mtplottools.plot_errorbar(axrxy,
                                          period[nzxy],
                                          rrp.resxy[nzxy],
                                          **kw_xx)
                reryx = mtplottools.plot_errorbar(axryx,
                                          period[nzyx],
                                          rrp.resyx[nzyx],
                                          **kw_yy)
                reryy = mtplottools.plot_errorbar(axryy,
                                          period[nzyy],
                                          rrp.resyy[nzyy],
                                          **kw_yy)
                #plot phase
                rerxx= mtplottools.plot_errorbar(axpxx,
                                          period[nzxx],
                                          rrp.phasexx[nzxx],
                                          **kw_xx)
                rerxy = mtplottools.plot_errorbar(axpxy,
                                          period[nzxy],
                                          rrp.phasexy[nzxy],
                                          **kw_xx)
                reryx = mtplottools.plot_errorbar(axpyx,
                                          period[nzyx],
                                          rrp.phaseyx[nzyx],
                                          **kw_yy)
                reryy = mtplottools.plot_errorbar(axpyy,
                                          period[nzyy],
                                          rrp.phaseyy[nzyy],
                                          **kw_yy)
            elif self.plot_z == True:
                #plot real
                rerxx = mtplottools.plot_errorbar(axrxx,
                                          period[nzxx],
                                          resp_z_obj.z[nzxx,0,0].real,
                                          **kw_xx)
                rerxy = mtplottools.plot_errorbar(axrxy,
                                          period[nzxy],
                                          resp_z_obj.z[nzxy,0,1].real,
                                          **kw_xx)
                reryx = mtplottools.plot_errorbar(axryx,
                                          period[nzyx],
                                          resp_z_obj.z[nzyx,1,0].real,
                                          **kw_yy)
                reryy = mtplottools.plot_errorbar(axryy,
                                          period[nzyy],
                                          resp_z_obj.z[nzyy,1,1].real,
                                          **kw_yy)
                #plot phase
                rerxx = mtplottools.plot_errorbar(axpxx,
                                          period[nzxx],
                                          resp_z_obj.z[nzxx,0,0].imag,
                                          **kw_xx)
                rerxy = mtplottools.plot_errorbar(axpxy,
                                          period[nzxy],
                                          resp_z_obj.z[nzxy,0,1].imag,
                                          **kw_xx)
                reryx = mtplottools.plot_errorbar(axpyx,
                                          period[nzyx],
                                          resp_z_obj.z[nzyx,1,0].imag,
                                          **kw_yy)
                reryy = mtplottools.plot_errorbar(axpyy,
                                          period[nzyy],
                                          resp_z_obj.z[nzyy,1,1].imag,
                                          **kw_yy)
            if self.plot_tipper == True:
                rertx = mtplottools.plot_errorbar(axtxr,
                             period[ntx],
                             resp_t_obj.tipper[ntx, 0, 0].real,
                             **kw_xx)
                rerty = mtplottools.plot_errorbar(axtyr,
                             period[nty],
                             resp_t_obj.tipper[nty, 0, 1].real,
                             **kw_yy)

                rertx = mtplottools.plot_errorbar(axtxi,
                             period[ntx],
                             resp_t_obj.tipper[ntx, 0, 0].imag,
                             **kw_xx)
                rerty = mtplottools.plot_errorbar(axtyi,
                             period[nty],
                             resp_t_obj.tipper[nty, 0, 1].imag,
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
                label_list[4] += ['$T^m_{x}$ ']
                label_list[5] += ['$T^m_{y}$']

            legend_ax_list = self.ax_list[0:self.plot_component]
            if self.plot_tipper == True:
                legend_ax_list.append(self.ax_list[8])
                legend_ax_list.append(self.ax_list[10])

            for aa, ax in enumerate(legend_ax_list):
                ax.legend(line_list[aa],
                          label_list[aa],
                          loc=self.legend_loc,
                          bbox_to_anchor=self.legend_pos,
                          markerscale=self.legend_marker_scale,
                          borderaxespad=self.legend_border_axes_pad,
                          labelspacing=self.legend_label_spacing,
                          handletextpad=self.legend_handle_text_pad,
                          borderpad=self.legend_border_pad,
                          prop={'size':max([self.fs, 5])})

        self.mpl_widget.draw()

    def on_pick(self, event):
        """
        mask a data point when it is clicked on.
        """
        data_point = event.artist
        data_period = data_point.get_xdata()[event.ind]
        data_value = data_point.get_ydata()[event.ind]
        p_index = np.where(self.modem_data.period_list==data_period)[0][0]
        s_index = np.where(self.modem_data.data_array['station']==self.station)[0][0]

        if event.mouseevent.button == 1:
            # mask the point in the data mt_dict
            if len(self.ax_list) == 8:
                self.modem_data.data_array[s_index]['z'][p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
                self.modem_data.mt_dict[self.station].Z.z[p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
            elif len(self.ax_list) == 12:
                if self._ax_index == 4 or self._ax_index == 5 or \
                   self._ax_index == 10 or self._ax_index == 11:
                    self.modem_data.data_array[s_index]['tip'][p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
                    self.modem_data.mt_dict[self.station].Tipper.tipper[p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
                else:
                    self.modem_data.data_array[s_index]['z'][p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
                    self.modem_data.mt_dict[self.station].Z.z[p_index,
                                self._comp_index_x, self._comp_index_y] = 0+0j
            # plot the point as masked
            self._ax.plot(data_period, data_value, color=(.7, .7, .7),
                          marker=self.mted, ms=self.ms)

        # Increase error bars
        if event.mouseevent.button == 3:
            # make sure just checking the top plots
            ax_index = self._ax_index%len(self._err_list)

            #put the new error into the error array
            if len(self.ax_list) == 8:
                err = self.modem_data.mt_dict[self.station].Z.z_err[p_index,
                        self._comp_index_x, self._comp_index_y]
                err = err+abs(err)*self.z_err_increase
                self.modem_data.data_array[s_index]['z_err'][p_index,
                            self._comp_index_x, self._comp_index_y] = err
                self.modem_data.mt_dict[self.station].Z.z_err[p_index,
                        self._comp_index_x, self._comp_index_y] = err
            elif len(self.ax_list) == 12:
                if self._ax_index == 4 or self._ax_index == 5 or \
                   self._ax_index == 10 or self._ax_index == 11:
                    err = self.modem_data.mt_dict[self.station].Tipper.tipper_err[p_index,
                                self._comp_index_x, self._comp_index_y]
                    self.modem_data.data_array[s_index]['tip_err'][p_index,
                                self._comp_index_x, self._comp_index_y] += abs(err)*self.t_err_increase
                    self.modem_data.mt_dict[self.station].Tipper.tipper_err[p_index,
                                self._comp_index_x, self._comp_index_y] += abs(err)*self.t_err_increase
                else:
                    err = self.modem_data.mt_dict[self.station].Z.z_err[p_index,
                        self._comp_index_x, self._comp_index_y]
                    self.modem_data.data_array[s_index]['z_err'][p_index,
                        self._comp_index_x, self._comp_index_y] += abs(err)*self.z_err_increase
                    self.modem_data.mt_dict[self.station].Z.z_err[p_index,
                        self._comp_index_x, self._comp_index_y] += abs(err)*self.z_err_increase

            # make error bar array
            eb = self._err_list[ax_index][2].get_paths()[p_index].vertices

            # make ecap array
            ecap_l = self._err_list[ax_index][0].get_data()[1][p_index]
            ecap_u = self._err_list[ax_index][1].get_data()[1][p_index]

            # change apparent resistivity error
            neb_u = eb[0,1]-.025*abs(eb[0,1])
            neb_l = eb[1,1]+.025*abs(eb[1,1])
            ecap_l = ecap_l-.025*abs(ecap_l)
            ecap_u = ecap_u+.025*abs(ecap_u)

            #set the new error bar values
            eb[0,1] = neb_u
            eb[1,1] = neb_l

            #reset the error bars and caps
            ncap_l = self._err_list[ax_index][0].get_data()
            ncap_u = self._err_list[ax_index][1].get_data()
            ncap_l[1][p_index] = ecap_l
            ncap_u[1][p_index] = ecap_u

            #set the values
            self._err_list[ax_index][0].set_data(ncap_l)
            self._err_list[ax_index][1].set_data(ncap_u)
            self._err_list[ax_index][2].get_paths()[p_index].vertices = eb

        # need to redraw the figure
        self._ax.figure.canvas.draw()

    def in_axes(self, event):
        """
        figure out which axes you just chose the point from
        """
        if len(self.ax_list) == 8:
            ax_index_dict = {0:(0, 0),
                             1:(0, 1),
                             2:(1, 0),
                             3:(1, 1),
                             4:(0, 0),
                             5:(0, 1),
                             6:(1, 0),
                             7:(1, 1)}
        if len(self.ax_list) == 12:
            ax_index_dict = {0:(0, 0),
                             1:(0, 1),
                             2:(1, 0),
                             3:(1, 1),
                             4:(0, 0),
                             5:(0, 1),
                             6:(1, 0),
                             7:(1, 1),
                             8:(0, 0),
                             9:(0, 1),
                             10:(0, 0),
                             11:(0, 1)}
        # make the axis an attribute
        self._ax = event.inaxes

        # find the component index so that it can be masked
        for ax_index, ax in enumerate(self.ax_list):
            if ax == event.inaxes:
                self._comp_index_x, self._comp_index_y = ax_index_dict[ax_index]
                self._ax_index = ax_index

#==============================================================================
#  Plot setting
#==============================================================================
class PlotSettings(QtGui.QWidget):
    settings_updated = QtCore.pyqtSignal()
    def __init__(self, parent, **kwargs):
        super(PlotSettings, self).__init__(parent)

        self.fs = kwargs.pop('fs', 10)
        self.lw = kwargs.pop('lw', 1.5)
        self.ms = kwargs.pop('ms', 5)

        self.e_capthick = kwargs.pop('e_capthick', 1)
        self.e_capsize =  kwargs.pop('e_capsize', 5)

        #color mode
        self.cted = kwargs.pop('cted', (0, 0, 1))
        self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
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

        self.z_err_increase = kwargs.pop('z_err_increase', 0.05)
        self.t_err_increase = kwargs.pop('t_err_increase', 0.05)

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
        self.fs_label = QtGui.QLabel('Font Size')
        self.fs_edit = QtGui.QLineEdit('{0:.1f}'.format(self.fs))
        self.fs_edit.editingFinished.connect(self.set_text_fs)

        self.lw_label = QtGui.QLabel('Line Width')
        self.lw_edit = QtGui.QLineEdit('{0:.1f}'.format(self.lw))
        self.lw_edit.editingFinished.connect(self.set_text_lw)

        self.e_capthick_label = QtGui.QLabel('Error cap thickness')
        self.e_capthick_edit = QtGui.QLineEdit('{0:.1f}'.format(self.e_capthick))
        self.e_capthick_edit.editingFinished.connect(self.set_text_e_capthick)

        self.e_capsize_label = QtGui.QLabel('Error cap size')
        self.e_capsize_edit = QtGui.QLineEdit('{0:.1f}'.format(self.e_capsize))
        self.e_capsize_edit.editingFinished.connect(self.set_text_e_capsize)

        grid_line = QtGui.QGridLayout()

        grid_line.addWidget(self.fs_label, 1, 0)
        grid_line.addWidget(self.fs_edit, 1, 1)

        grid_line.addWidget(self.lw_label, 1, 2)
        grid_line.addWidget(self.lw_edit, 1, 3)

        grid_line.addWidget(self.e_capthick_label, 1, 4)
        grid_line.addWidget(self.e_capthick_edit, 1, 5)

        grid_line.addWidget(self.e_capsize_label, 1, 6)
        grid_line.addWidget(self.e_capsize_edit, 1, 7)

        #--> marker properties
        self.ms_label = QtGui.QLabel('Marker Size')
        self.ms_edit = QtGui.QLineEdit()
        self.ms_edit.setText('{0:.1f}'.format(self.ms))
        self.ms_edit.editingFinished.connect(self.set_text_ms)

        self.mted_label = QtGui.QLabel("Marker x components (data)")
        self.mted_combo = QtGui.QComboBox()
        self.mted_combo.addItem(self.mted)
        self.mted_combo.addItem('.')
        self.mted_combo.addItem(',')
        self.mted_combo.addItem('o')
        self.mted_combo.addItem('v')
        self.mted_combo.addItem('^')
        self.mted_combo.addItem('<')
        self.mted_combo.addItem('>')
        self.mted_combo.addItem('s')
        self.mted_combo.addItem('p')
        self.mted_combo.addItem('*')
        self.mted_combo.addItem('h')
        self.mted_combo.addItem('H')
        self.mted_combo.addItem('+')
        self.mted_combo.addItem('x')
        self.mted_combo.addItem('D')
        self.mted_combo.addItem('d')
        self.mted_combo.addItem('|')
        self.mted_combo.addItem('_')
        self.mted_combo.activated[str].connect(self.set_mted)

        self.mtmd_label = QtGui.QLabel("Marker y components (data)")
        self.mtmd_combo = QtGui.QComboBox()
        self.mtmd_combo.addItem(self.mtmd)
        self.mtmd_combo.addItem('.')
        self.mtmd_combo.addItem(',')
        self.mtmd_combo.addItem('o')
        self.mtmd_combo.addItem('v')
        self.mtmd_combo.addItem('^')
        self.mtmd_combo.addItem('<')
        self.mtmd_combo.addItem('>')
        self.mtmd_combo.addItem('s')
        self.mtmd_combo.addItem('p')
        self.mtmd_combo.addItem('*')
        self.mtmd_combo.addItem('h')
        self.mtmd_combo.addItem('H')
        self.mtmd_combo.addItem('+')
        self.mtmd_combo.addItem('x')
        self.mtmd_combo.addItem('D')
        self.mtmd_combo.addItem('d')
        self.mtmd_combo.addItem('|')
        self.mtmd_combo.addItem('_')
        self.mtmd_combo.activated[str].connect(self.set_mtmd)

        self.mtem_label = QtGui.QLabel("Marker x components (model)")
        self.mtem_combo = QtGui.QComboBox()
        self.mtem_combo.addItem(self.mtem)
        self.mtem_combo.addItem('.')
        self.mtem_combo.addItem(',')
        self.mtem_combo.addItem('o')
        self.mtem_combo.addItem('v')
        self.mtem_combo.addItem('^')
        self.mtem_combo.addItem('<')
        self.mtem_combo.addItem('>')
        self.mtem_combo.addItem('s')
        self.mtem_combo.addItem('p')
        self.mtem_combo.addItem('*')
        self.mtem_combo.addItem('h')
        self.mtem_combo.addItem('H')
        self.mtem_combo.addItem('+')
        self.mtem_combo.addItem('x')
        self.mtem_combo.addItem('D')
        self.mtem_combo.addItem('d')
        self.mtem_combo.addItem('|')
        self.mtem_combo.addItem('_')
        self.mtem_combo.activated[str].connect(self.set_mtem)

        self.mtmm_label = QtGui.QLabel("Marker y components (model)")
        self.mtmm_combo = QtGui.QComboBox()
        self.mtmm_combo.addItem(self.mtmm)
        self.mtmm_combo.addItem('.')
        self.mtmm_combo.addItem(',')
        self.mtmm_combo.addItem('o')
        self.mtmm_combo.addItem('v')
        self.mtmm_combo.addItem('^')
        self.mtmm_combo.addItem('<')
        self.mtmm_combo.addItem('>')
        self.mtmm_combo.addItem('s')
        self.mtmm_combo.addItem('p')
        self.mtmm_combo.addItem('*')
        self.mtmm_combo.addItem('h')
        self.mtmm_combo.addItem('H')
        self.mtmm_combo.addItem('+')
        self.mtmm_combo.addItem('x')
        self.mtmm_combo.addItem('D')
        self.mtmm_combo.addItem('d')
        self.mtmm_combo.addItem('|')
        self.mtmm_combo.addItem('_')
        self.mtmm_combo.activated[str].connect(self.set_mtmm)

        self.cted_button = QtGui.QPushButton("Set Z_xi Color (data)")
        self.cted_button.pressed.connect(self.set_cted)

        self.ctmd_button = QtGui.QPushButton("Set Z_yi Color (data)")
        self.ctmd_button.pressed.connect(self.set_ctmd)


        self.ctem_button = QtGui.QPushButton("Set Z_xi Color (model)")
        self.ctem_button.pressed.connect(self.set_ctem)

        self.ctmm_button = QtGui.QPushButton("Set Z_yi Color (model)")
        self.ctmm_button.pressed.connect(self.set_ctmm)


        marker_size_grid = QtGui.QGridLayout()
        marker_size_grid.addWidget(self.ms_label, 0, 0)
        marker_size_grid.addWidget(self.ms_edit, 0, 1)

        marker_grid = QtGui.QGridLayout()

        marker_grid.addWidget(self.mted_label, 0, 0)
        marker_grid.addWidget(self.mted_combo, 0, 1)

        marker_grid.addWidget(self.mtmd_label, 0, 2)
        marker_grid.addWidget(self.mtmd_combo, 0, 3)

        marker_grid.addWidget(self.mtem_label, 0, 4)
        marker_grid.addWidget(self.mtem_combo, 0, 5)

        marker_grid.addWidget(self.mtmm_label, 0, 6)
        marker_grid.addWidget(self.mtmm_combo, 0, 7)

        marker_grid.addWidget(self.cted_button, 1, 0, 2, 2)
        marker_grid.addWidget(self.ctmd_button, 1, 2, 2, 2)
        marker_grid.addWidget(self.ctem_button, 1, 4, 2, 2)
        marker_grid.addWidget(self.ctmm_button, 1, 6, 2, 2)

        marker_layout = QtGui.QVBoxLayout()
        marker_layout.addLayout(marker_size_grid)
        marker_layout.addLayout(marker_grid)


        #-> increase error by this percent
        self.z_err_inc_label = QtGui.QLabel('Z error increase (%)')
        self.z_err_inc_edit = QtGui.QLineEdit('{0:.2f}'.format(self.z_err_increase*100))
        self.z_err_inc_edit.editingFinished.connect(self.set_z_err_inc)

        self.t_err_inc_label = QtGui.QLabel('Tipper error increase (%)')
        self.t_err_inc_edit = QtGui.QLineEdit('{0:.2f}'.format(self.z_err_increase*100))
        self.t_err_inc_edit.editingFinished.connect(self.set_z_err_inc)

        err_grid = QtGui.QGridLayout()
        err_grid.addWidget(self.z_err_inc_label, 0, 0)
        err_grid.addWidget(self.z_err_inc_edit, 0, 1)
        err_grid.addWidget(self.t_err_inc_label, 0, 2)
        err_grid.addWidget(self.t_err_inc_edit, 0, 3)

        #--> plot limits
        self.ylimr_xx_label = QtGui.QLabel('Res_xx')
        self.ylimr_xx_edit = QtGui.QLineEdit()
        self.ylimr_xx_edit.setText('{0}'.format(self.res_xx_limits))
        self.ylimr_xx_edit.editingFinished.connect(self.set_text_res_xx)

        self.ylimr_xy_label = QtGui.QLabel('Res_xy')
        self.ylimr_xy_edit = QtGui.QLineEdit()
        self.ylimr_xy_edit.setText('{0}'.format(self.res_xy_limits))
        self.ylimr_xy_edit.editingFinished.connect(self.set_text_res_xy)

        self.ylimr_yx_label = QtGui.QLabel('Res_yx')
        self.ylimr_yx_edit = QtGui.QLineEdit()
        self.ylimr_yx_edit.setText('{0}'.format(self.res_yx_limits))
        self.ylimr_yx_edit.editingFinished.connect(self.set_text_res_yx)

        self.ylimr_yy_label = QtGui.QLabel('Res_yy')
        self.ylimr_yy_edit = QtGui.QLineEdit()
        self.ylimr_yy_edit.setText('{0}'.format(self.res_yy_limits))
        self.ylimr_yy_edit.editingFinished.connect(self.set_text_res_yy)

        self.ylimp_xx_label = QtGui.QLabel('phase_xx')
        self.ylimp_xx_edit = QtGui.QLineEdit()
        self.ylimp_xx_edit.setText('{0}'.format(self.phase_xx_limits))
        self.ylimp_xx_edit.editingFinished.connect(self.set_text_phase_xx)

        self.ylimp_xy_label = QtGui.QLabel('phase_xy')
        self.ylimp_xy_edit = QtGui.QLineEdit()
        self.ylimp_xy_edit.setText('{0}'.format(self.phase_xy_limits))
        self.ylimp_xy_edit.editingFinished.connect(self.set_text_phase_xy)

        self.ylimp_yx_label = QtGui.QLabel('phase_yx')
        self.ylimp_yx_edit = QtGui.QLineEdit()
        self.ylimp_yx_edit.setText('{0}'.format(self.phase_yx_limits))
        self.ylimp_yx_edit.editingFinished.connect(self.set_text_phase_yx)

        self.ylimp_yy_label = QtGui.QLabel('phase_yy')
        self.ylimp_yy_edit = QtGui.QLineEdit()
        self.ylimp_yy_edit.setText('{0}'.format(self.phase_yy_limits))
        self.ylimp_yy_edit.editingFinished.connect(self.set_text_phase_yy)

        self.ytlim_label = QtGui.QLabel('tipper')
        self.ytlim_edit = QtGui.QLineEdit()
        self.ytlim_edit.setText('{0}'.format(self.tipper_limits))
        self.ytlim_edit.editingFinished.connect(self.set_text_tipper_lim)

        limits_grid = QtGui.QGridLayout()

        limits_label = QtGui.QLabel('Plot Limits: (Res=Real, Phase=Imaginary)'
                                    ' --> input on a linear scale')

        limits_grid.addWidget(limits_label, 0, 0, 7, 1)

        limits_grid.addWidget(self.ylimr_xx_label, 2, 0)
        limits_grid.addWidget(self.ylimr_xx_edit, 2, 1)
        limits_grid.addWidget(self.ylimr_xy_label, 2, 2)
        limits_grid.addWidget(self.ylimr_xy_edit, 2, 3)
        limits_grid.addWidget(self.ylimr_yx_label, 2, 4)
        limits_grid.addWidget(self.ylimr_yx_edit, 2, 5)
        limits_grid.addWidget(self.ylimr_yy_label, 2, 6)
        limits_grid.addWidget(self.ylimr_yy_edit, 2, 7)

        limits_grid.addWidget(self.ylimp_xx_label, 3, 0)
        limits_grid.addWidget(self.ylimp_xx_edit, 3, 1)
        limits_grid.addWidget(self.ylimp_xy_label, 3, 2)
        limits_grid.addWidget(self.ylimp_xy_edit, 3, 3)
        limits_grid.addWidget(self.ylimp_yx_label, 3, 4)
        limits_grid.addWidget(self.ylimp_yx_edit, 3, 5)
        limits_grid.addWidget(self.ylimp_yy_label, 3, 6)
        limits_grid.addWidget(self.ylimp_yy_edit, 3, 7)

        limits_grid.addWidget(self.ytlim_label, 4, 0)
        limits_grid.addWidget(self.ytlim_edit, 4, 1)

        #--> legend properties
        self.legend_pos_label = QtGui.QLabel('Legend Position')
        self.legend_pos_edit = QtGui.QLineEdit()
        self.legend_pos_edit.setText('{0}'.format(self.legend_pos))
        self.legend_pos_edit.editingFinished.connect(self.set_text_legend_pos)

        legend_grid = QtGui.QGridLayout()
        legend_grid.setSpacing(10)

        legend_grid.addWidget(QtGui.QLabel('Legend Properties:'), 1, 0)
        legend_grid.addWidget(self.legend_pos_label, 1, 2,)
        legend_grid.addWidget(self.legend_pos_edit, 1, 3)

        update_button = QtGui.QPushButton('Update')
        update_button.clicked.connect(self.update_settings)

        vbox = QtGui.QVBoxLayout()
        vbox.addLayout(grid_line)
        vbox.addLayout(marker_layout)
        vbox.addLayout(err_grid)
        vbox.addLayout(limits_grid)
        vbox.addLayout(legend_grid)
        vbox.addWidget(update_button)

        self.setLayout(vbox)

        #self.setGeometry(300, 300, 350, 300)
        #self.resize(1350, 500)
        self.setWindowTitle('Plot Settings')
        self.show()

    def _test_float(self, text):
        try:
            return_float = float(str(text))
        except ValueError:
            print('Enter a floating point number')
            return_float = 1.0

        return return_float

    def set_text_fs(self):
        self.fs = self._test_float(self.fs_edit.text())
        self.fs_edit.setText('{0:.1f}'.format(self.fs))

    def set_text_e_capthick(self):
        self.e_capthick_edit = self._test_float(self.e_capthick_edit.text())
        self.e_capthick_edit.setText('{0:.1f}'.format(self.e_capthick))

    def set_text_e_capsize(self):
        self.e_capsize_edit = self._test_float(self.e_capsize_edit.text())
        self.e_capsize_edit.setText('{0:.1f}'.format(self.e_capsize))

    def set_text_lw(self):
        self.lw = self._test_float(self.lw_edit.text())
        self.lw_edit.setText('{0:.2f}'.format(self.lw))

    def set_text_ms(self, text):
        self.ms = self._test_float(self.ms_edit.text())
        self.ms_edit.setText('{0:.1f}'.format(self.ms))

    def set_mted(self, text):
        self.mted = text

    def set_mtmd(self, text):
        self.mtmd = text

    def set_mtem(self, text):
        self.mtem = text

    def set_mtmm(self, text):
        self.mtmm = text

    def convert_color_to_qt(self, color):
        """
        convert decimal tuple to QColor object
        """
        r = int(color[0]*255)
        g = int(color[1]*255)
        b = int(color[2]*255)

        return QtGui.QColor(r, g, b)

    def set_cted(self):
        initial_color = self.convert_color_to_qt(self.cted)
        new_color = QtGui.QColorDialog.getColor(initial_color)

        r,g,b,a = new_color.getRgbF()

        self.cted = (r, g, b)

    def set_ctmd(self):
        initial_color = self.convert_color_to_qt(self.ctmd)
        new_color = QtGui.QColorDialog.getColor(initial_color)

        r,g,b,a = new_color.getRgbF()

        self.ctmd = (r, g, b)

    def set_ctem(self):
        initial_color = self.convert_color_to_qt(self.ctem)
        new_color = QtGui.QColorDialog.getColor(initial_color)

        r,g,b,a = new_color.getRgbF()

        self.ctem = (r, g, b)

    def set_ctmm(self):
        initial_color = self.convert_color_to_qt(self.ctmm)
        new_color = QtGui.QColorDialog.getColor(initial_color)

        r,g,b,a = new_color.getRgbF()

        self.ctmm = (r, g, b)

    def set_text_res_xx(self):
        text = self.ylimr_xx_edit.text()
        if text =='None':
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
        text = self.ylimr_xy_edit.text()
        if text =='None':
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
        text = self.ylimr_yx_edit.text()
        if text =='None':
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
        text = self.ylimr_yy_edit.text()
        if text =='None':
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
        text = self.ylimp_xx_edit.text()
        if text =='None':
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
        text = self.ylimp_xy_edit.text()
        if text =='None':
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
        text = self.ylimp_yx_edit.text()
        if text =='None':
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
        text = self.ylimp_yy_edit.text()
        if text =='None':
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

    def set_text_tipper_lim(self, text):
        text = self.ytlim_yy_edit.text()
        if text =='None':
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
            self.tipper_limits = tuple(l_list)

    def set_text_legend_pos(self):
        text = self.legend_pos_edit.text()
        if text =='None':
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

    def set_z_err_inc(self):
        err = float(str(self.z_err_inc_edit.text()))/100.
        self.z_err_increase = err
        self.z_err_inc_edit.setText('{0:.2f}'.format(err*100.))

    def set_t_err_inc(self):
        err = float(str(self.t_err_inc_edit.text()))/100.
        self.t_err_increase = err
        self.t_err_inc_edit.setText('{0:.2f}'.format(err*100.))


    def update_settings(self):
        self.settings_updated.emit()

#def main():

def main():
#if __name__ == "__main__":
    import sys
    app = QtGui.QApplication(sys.argv)
    MainWindow = QtGui.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

if __name__ == '__main__':

    main()
