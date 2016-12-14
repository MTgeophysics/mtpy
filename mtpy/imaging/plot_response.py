# ==============================================================================
# plot response
# ==============================================================================

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

import mtpy.imaging.mtplottools as mtplottools
from mtpy.modeling.modem import Data

try:
    from pyevtk.hl import gridToVTK, pointsToVTK
except ImportError:
    print ('If you want to write a vtk file for 3d viewing, you need to pip install PyEVTK:'
           ' https://bitbucket.org/pauloh/pyevtk')

    print ('Note: if you are using Windows you should build evtk first with'
           'either MinGW or cygwin using the command: \n'
           '    python setup.py build -compiler=mingw32  or \n'
           '    python setup.py build -compiler=cygwin')


class PlotResponse(object):
    """
    plot data and response

    Plots the real and imaginary impedance and induction vector if present.

    :Example: ::

        >>> import mtpy.modeling.new_modem as modem
        >>> dfn = r"/home/MT/ModEM/Inv1/DataFile.dat"
        >>> rfn = r"/home/MT/ModEM/Inv1/Test_resp_000.dat"
        >>> mrp = modem.PlotResponse(data_fn=dfn, resp_fn=rfn)
        >>> # plot only the TE and TM modes
        >>> mrp.plot_component = 2
        >>> mrp.redraw_plot()

    ======================== ==================================================
    Attributes               Description
    ======================== ==================================================
    color_mode               [ 'color' | 'bw' ] color or black and white plots
    cted                     color for data TE mode
    ctem                     color for data TM mode
    ctmd                     color for model TE mode
    ctmm                     color for model TM mode
    data_fn                  full path to data file
    data_object              WSResponse instance
    e_capsize                cap size of error bars in points (*default* is .5)
    e_capthick               cap thickness of error bars in points (*default*
                             is 1)
    fig_dpi                  resolution of figure in dots-per-inch (300)
    fig_list                 list of matplotlib.figure instances for plots
    fig_size                 size of figure in inches (*default* is [6, 6])
    font_size                size of font for tick labels, axes labels are
                             font_size+2 (*default* is 7)
    legend_border_axes_pad   padding between legend box and axes
    legend_border_pad        padding between border of legend and symbols
    legend_handle_text_pad   padding between text labels and symbols of legend
    legend_label_spacing     padding between labels
    legend_loc               location of legend
    legend_marker_scale      scale of symbols in legend
    lw                       line width response curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    mted                     marker for data TE mode
    mtem                     marker for data TM mode
    mtmd                     marker for model TE mode
    mtmm                     marker for model TM mode
    phase_limits             limits of phase
    plot_component           [ 2 | 4 ] 2 for TE and TM or 4 for all components
    plot_style               [ 1 | 2 ] 1 to plot each mode in a seperate
                             subplot and 2 to plot xx, xy and yx, yy in same
                             plots
    plot_type                [ '1' | list of station name ] '1' to plot all
                             stations in data file or input a list of station
                             names to plot if station_fn is input, otherwise
                             input a list of integers associated with the
                             index with in the data file, ie 2 for 2nd station
    plot_z                   [ True | False ] *default* is True to plot
                             impedance, False for plotting resistivity and
                             phase
    plot_yn                  [ 'n' | 'y' ] to plot on instantiation
    res_limits               limits of resistivity in linear scale
    resp_fn                  full path to response file
    resp_object              WSResponse object for resp_fn, or list of
                             WSResponse objects if resp_fn is a list of
                             response files
    station_fn               full path to station file written by WSStation
    subplot_bottom           space between axes and bottom of figure
    subplot_hspace           space between subplots in vertical direction
    subplot_left             space between axes and left of figure
    subplot_right            space between axes and right of figure
    subplot_top              space between axes and top of figure
    subplot_wspace           space between subplots in horizontal direction
    ======================== ==================================================
    """

    def __init__(self, data_fn=None, resp_fn=None, **kwargs):
        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.data_object = None
        self.resp_object = []

        self.color_mode = kwargs.pop('color_mode', 'color')

        self.ms = kwargs.pop('ms', 1.5)
        self.lw = kwargs.pop('lw', .5)
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0, .6, .3))
            self.ctmm = kwargs.pop('ctmm', (.9, 0, .8))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', '+')

        # black and white mode
        elif self.color_mode == 'bw':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 0))
            self.ctmd = kwargs.pop('ctmd', (0, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            self.ctem = kwargs.pop('ctem', (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop('ctmm', (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop('mtem', '+')
            self.mtmm = kwargs.pop('mtmm', 'x')

        self.phase_limits = kwargs.pop('phase_limits', None)
        self.res_limits = kwargs.pop('res_limits', None)

        self.fig_num = kwargs.pop('fig_num', 1)
        self.fig_size = kwargs.pop('fig_size', [6, 6])
        self.fig_dpi = kwargs.pop('dpi', 300)

        self.subplot_wspace = kwargs.pop('subplot_wspace', .3)
        self.subplot_hspace = kwargs.pop('subplot_hspace', .0)
        self.subplot_right = kwargs.pop('subplot_right', .98)
        self.subplot_left = kwargs.pop('subplot_left', .08)
        self.subplot_top = kwargs.pop('subplot_top', .85)
        self.subplot_bottom = kwargs.pop('subplot_bottom', .1)

        self.legend_loc = 'upper center'
        self.legend_pos = (.5, 1.21)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')
        self.plot_style = kwargs.pop('plot_style', 1)
        self.plot_component = kwargs.pop('plot_component', 4)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.plot_z = kwargs.pop('plot_z', True)
        self.ylabel_pad = kwargs.pop('ylabel_pad', 1.25)

        self.fig_list = []

        # this __init__ is a constructor which creates an object pObj, call pObj.plot() method will do
        # if self.plot_yn == 'y':
        #     self.plot()

        return

    def plot(self, save2file=None):
        """
        plot show figure and optionally save to a file named save2file
        """

        self.data_object = Data()
        self.data_object.read_data_file(self.data_fn)

        # get shape of impedance tensors
        ns = len(self.data_object.mt_dict.keys())

        # read in response files
        if self.resp_fn != None:
            self.resp_object = []
            if type(self.resp_fn) is not list:
                resp_obj = Data()
                resp_obj.read_data_file(self.resp_fn)
                self.resp_object = [resp_obj]
            else:
                for rfile in self.resp_fn:
                    resp_obj = Data()
                    resp_obj.read_data_file(rfile)
                    self.resp_object.append(resp_obj)

        # get number of response files
        nr = len(self.resp_object)

        if type(self.plot_type) is list:
            ns = len(self.plot_type)

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}
        if self.plot_z == True:
            h_ratio = [1, 1]
        elif self.plot_z == False:
            h_ratio = [2, 1.5]

        ax_list = []
        line_list = []
        label_list = []

        # --> make key word dictionaries for plotting
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

        if self.plot_type != '1':
            pstation_list = []
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]
            for ii, station in enumerate(self.data_object.mt_dict.keys()):
                if type(station) is not int:
                    for pstation in self.plot_type:
                        if station.find(str(pstation)) >= 0:
                            pstation_list.append(station)
                else:
                    for pstation in self.plot_type:
                        if station == int(pstation):
                            pstation_list.append(ii)
        else:
            pstation_list = self.data_object.mt_dict.keys()

        for jj, station in enumerate(pstation_list):
            z_obj = self.data_object.mt_dict[station].Z
            t_obj = self.data_object.mt_dict[station].Tipper
            period = self.data_object.period_list
            print 'Plotting: {0}'.format(station)

            # convert to apparent resistivity and phase
            rp = mtplottools.ResPhase(z_object=z_obj)

            # find locations where points have been masked
            nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
            nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
            nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
            nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
            ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
            nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]

            if self.resp_fn != None:
                plotr = True
            else:
                plotr = False

            # make figure
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)

            # set the grid of subplots
            tipper_zero = (np.round(abs(t_obj.tipper.mean()), 4) == 0.0)
            if tipper_zero == False:
                # makes more sense if plot_tipper is True to plot tipper
                plot_tipper = True
            else:
                plot_tipper = False

            if plot_tipper == True:

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
            # ---------plot the apparent resistivity-----------------------------------
            # plot each component in its own subplot
            if self.plot_style == 1:
                # plot xy and yx
                if self.plot_component == 2:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axryx = fig.add_subplot(gs[0, 2:], sharex=axrxy)

                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        axpyx = fig.add_subplot(gs[1, 2:], sharex=axrxy)
                    else:

                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axryx = fig.add_subplot(gs[0, 2:4], sharex=axrxy)

                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        axpyx = fig.add_subplot(gs[1, 2:4], sharex=axrxy)

                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)
                        axtr.set_ylim(-1.2, 1.2)
                        axti.set_ylim(-1.2, 1.2)

                    if self.plot_z == False:
                        # plot resistivity
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period,
                                                         rp.resxy[nzxy],
                                                         rp.resxy_err[nzxy],
                                                         **kw_xx)

                        eryx = mtplottools.plot_errorbar(axryx,
                                                         period[nzyx],
                                                         rp.resyx[nzyx],
                                                         rp.resyx_err[nzyx],
                                                         **kw_yy)
                        # plot phase
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         rp.phasexy[nzxy],
                                                         rp.phasexy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx,
                                                         period[nzyx],
                                                         rp.phaseyx[nzyx],
                                                         rp.phaseyx_err[nzyx],
                                                         **kw_yy)

                    elif self.plot_z == True:
                        # plot real
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axryx,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].real),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                        # plot phase
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                    # plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr,
                                                         period[ntx],
                                                         t_obj.tipper[ntx, 0, 0].real,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr,
                                                         period[nty],
                                                         t_obj.tipper[nty, 0, 1].real,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                        ertx = mtplottools.plot_errorbar(axti,
                                                         period[ntx],
                                                         t_obj.tipper[ntx, 0, 0].imag,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axti,
                                                         period[nty],
                                                         t_obj.tipper[nty, 0, 1].imag,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                    if plot_tipper == False:
                        ax_list = [axrxy, axryx, axpxy, axpyx]
                        line_list = [[erxy[0]], [eryx[0]]]
                        label_list = [['$Z_{xy}$'], ['$Z_{yx}$']]
                    else:
                        ax_list = [axrxy, axryx, axpxy, axpyx, axtr, axti]
                        line_list = [[erxy[0]], [eryx[0]],
                                     [ertx[0], erty[0]]]
                        label_list = [['$Z_{xy}$'], ['$Z_{yx}$'],
                                      ['$T_{x}$', '$T_{y}$']]

                elif self.plot_component == 4:
                    if plot_tipper == False:
                        axrxx = fig.add_subplot(gs[0, 0])
                        axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
                        axryx = fig.add_subplot(gs[0, 2], sharex=axrxx)
                        axryy = fig.add_subplot(gs[0, 3], sharex=axrxx)

                        axpxx = fig.add_subplot(gs[1, 0])
                        axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
                        axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
                        axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)
                    else:
                        axrxx = fig.add_subplot(gs[0, 0])
                        axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
                        axryx = fig.add_subplot(gs[0, 2], sharex=axrxx)
                        axryy = fig.add_subplot(gs[0, 3], sharex=axrxx)

                        axpxx = fig.add_subplot(gs[1, 0])
                        axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
                        axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
                        axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)

                        axtxr = fig.add_subplot(gs[0, 4], sharex=axrxx)
                        axtxi = fig.add_subplot(gs[1, 4], sharex=axrxx)
                        axtyr = fig.add_subplot(gs[0, 5], sharex=axrxx)
                        axtyi = fig.add_subplot(gs[1, 5], sharex=axrxx)

                        axtxr.set_ylim(-1.2, 1.2)
                        axtxi.set_ylim(-1.2, 1.2)
                        axtyr.set_ylim(-1.2, 1.2)
                        axtyi.set_ylim(-1.2, 1.2)

                    if self.plot_z == False:
                        # plot resistivity
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
                        # plot phase
                        erxx = mtplottools.plot_errorbar(axpxx,
                                                         period[nzxx],
                                                         rp.phasexx[nzxx],
                                                         rp.phasexx_err[nzxx],
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         rp.phasexy[nzxy],
                                                         rp.phasexy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx,
                                                         period[nzyx],
                                                         rp.phaseyx[nzyx],
                                                         rp.phaseyx_err[nzyx],
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpyy,
                                                         period[nzyy],
                                                         rp.phaseyy[nzyy],
                                                         rp.phaseyy_err[nzyy],
                                                         **kw_yy)
                    elif self.plot_z == True:
                        # plot real
                        erxx = mtplottools.plot_errorbar(axrxx,
                                                         period[nzxx],
                                                         abs(z_obj.z[nzxx, 0, 0].real),
                                                         abs(z_obj.z_err[nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axryx,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].real),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axryy,
                                                         period[nzyy],
                                                         abs(z_obj.z[nzyy, 1, 1].real),
                                                         abs(z_obj.z_err[nzyy, 1, 1].real),
                                                         **kw_yy)
                        # plot phase
                        erxx = mtplottools.plot_errorbar(axpxx,
                                                         period[nzxx],
                                                         abs(z_obj.z[nzxx, 0, 0].imag),
                                                         abs(z_obj.z_err[nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpyx,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpyy,
                                                         period[nzyy],
                                                         abs(z_obj.z[nzyy, 1, 1].imag),
                                                         abs(z_obj.z_err[nzyy, 1, 1].real),
                                                         **kw_yy)

                    # plot tipper
                    if plot_tipper == True:
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

                        ertx = mtplottools.plot_errorbar(axtxi,
                                                         period[ntx],
                                                         t_obj.tipper[ntx, 0, 0].imag,
                                                         t_obj.tipper_err[ntx, 0, 1],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtyi,
                                                         period[nty],
                                                         t_obj.tipper[nty, 0, 1].imag,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)
                    if plot_tipper == False:
                        ax_list = [axrxx, axrxy, axryx, axryy,
                                   axpxx, axpxy, axpyx, axpyy]
                        line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                        label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                                      ['$Z_{yx}$'], ['$Z_{yy}$']]
                    else:
                        ax_list = [axrxx, axrxy, axryx, axryy,
                                   axpxx, axpxy, axpyx, axpyy,
                                   axtxr, axtxi, axtyr, axtyi]
                        line_list = [[erxx[0]], [erxy[0]],
                                     [eryx[0]], [eryy[0]],
                                     [ertx[0]], [erty[0]]]
                        label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                                      ['$Z_{yx}$'], ['$Z_{yy}$'],
                                      ['$T_{x}$'], ['$T_{y}$']]

                # set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    #                    ylabels = ax.get_yticks().tolist()
                    #                    ylabels[-1] = ''
                    #                    ylabels[0] = ''
                    #                    ax.set_yticklabels(ylabels)
                    #                    print ylabels

                    #                    dy = abs(ax.yaxis.get_ticklocs()[1]-
                    #                             ax.yaxis.get_ticklocs()[0])
                    #                    ylim = ax.get_ylim()
                    #                    ax.set_ylim(ylim[0]-.25*dy, ylim[1]+1.25*dy)
                    #                    ax.yaxis.set_major_locator(MultipleLocator(dy))

                    if len(ax_list) == 4:
                        #                        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                        if self.plot_z == True:
                            ax.set_yscale('log', nonposy='clip')
                            ylim = ax.get_ylim()
                            ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                       10 ** np.ceil(np.log10(ylim[1])))
                            ax.set_ylim(ylimits)
                            ylabels = [' '] + \
                                      [mtplottools.labeldict[ii] for ii
                                       in np.arange(np.log10(ylimits[0]),
                                                    np.log10(ylimits[1]), 1)] + \
                                      [' ']
                            ax.set_yticklabels(ylabels)
                    if len(ax_list) == 6:
                        if aa < 4:
                            #                            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                            if self.plot_z == True:
                                ax.set_yscale('log', nonposy='clip')
                                ylim = ax.get_ylim()
                                ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                           10 ** np.ceil(np.log10(ylim[1])))
                                ax.set_ylim(ylimits)
                                ylabels = [' '] + \
                                          [mtplottools.labeldict[ii] for ii
                                           in np.arange(np.log10(ylimits[0]),
                                                        np.log10(ylimits[1]), 1)] + \
                                          [' ']
                                ax.set_yticklabels(ylabels)
                    if len(ax_list) == 8:
                        #                        ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                        if self.plot_z == True:
                            ax.set_yscale('log', nonposy='clip')
                            ylim = ax.get_ylim()
                            ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                       10 ** np.ceil(np.log10(ylim[1])))
                            ax.set_ylim(ylimits)
                            ylabels = [' '] + \
                                      [mtplottools.labeldict[ii] for ii
                                       in np.arange(np.log10(ylimits[0]),
                                                    np.log10(ylimits[1]), 1)] + \
                                      [' ']
                            ax.set_yticklabels(ylabels)
                    if len(ax_list) == 12:
                        if aa < 4:
                            ylabels = ax.get_yticks().tolist()
                            ylabels[0] = ''
                            ax.set_yticklabels(ylabels)
                        if aa < 8:
                            #                            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                            if self.plot_z == True:
                                ax.set_yscale('log', nonposy='clip')
                                ylim = ax.get_ylim()
                                ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                           10 ** np.ceil(np.log10(ylim[1])))
                                ax.set_ylim(ylimits)
                                ylabels = [' '] + \
                                          [mtplottools.labeldict[ii] for ii
                                           in np.arange(np.log10(ylimits[0]),
                                                        np.log10(ylimits[1]), 1)] + \
                                          [' ']
                                ax.set_yticklabels(ylabels)
                    if len(ax_list) == 4 or len(ax_list) == 6:
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)

                        # set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Re[Z]| (mV/km nT)',
                                              fontdict=fontdict)
                        elif aa == 2:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Im[Z]| (mV/km nT)',
                                              fontdict=fontdict)

                    elif len(ax_list) == 8 or len(ax_list) == 12:
                        if aa < 4:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                                ylim = ax.get_ylim()
                                ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                           10 ** np.ceil(np.log10(ylim[1])))
                                ax.set_ylim(ylimits)
                                ylabels = [' ', ' '] + \
                                          [mtplottools.labeldict[ii] for ii
                                           in np.arange(np.log10(ylimits[0]) + 1,
                                                        np.log10(ylimits[1]) + 1, 1)]
                                ax.set_yticklabels(ylabels)
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            if aa == 8 or aa == 10:
                                plt.setp(ax.get_xticklabels(), visible=False)
                            else:
                                ax.set_ylim(self.phase_limits)
                                ax.set_xlabel('Period (s)', fontdict=fontdict)

                        # set axes labels
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Re[Z]| (mV/km nT)',
                                              fontdict=fontdict)
                        elif aa == 4:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Im[Z]| (mV/km nT)',
                                              fontdict=fontdict)

                    ax.set_xscale('log')
                    ax.set_xlim(xmin=10 ** (np.floor(np.log10(period[0]))) * 1.01,
                                xmax=10 ** (np.ceil(np.log10(period[-1]))) * .99)
                    ax.grid(True, alpha=.25)

            # plot xy and yx together and xx, yy together
            elif self.plot_style == 2:
                if self.plot_component == 2:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:])
                        axpxy = fig.add_subplot(gs[1, 0:], sharex=axrxy)
                    else:
                        axrxy = fig.add_subplot(gs[0, 0:4])
                        axpxy = fig.add_subplot(gs[1, 0:4], sharex=axrxy)
                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)

                    if self.plot_z == False:
                        # plot resistivity
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         rp.resxy[nzxy],
                                                         rp.resxy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzyx],
                                                         rp.resyx[nzyx],
                                                         rp.resyx_err[nzyx],
                                                         **kw_yy)
                        # plot phase
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         rp.phasexy[nzxy],
                                                         rp.phasexy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         rp.phaseyx[nzyx],
                                                         rp.phaseyx_err[nzyx],
                                                         **kw_yy)
                    elif self.plot_z == True:
                        # plot real
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 1, 0].real),
                                                         abs(z_obj.z_err[nzxy, 1, 0].real),
                                                         **kw_yy)
                        # plot phase
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                    # plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr,
                                                         period,
                                                         t_obj.tipper[ntx, 0, 0].real,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr,
                                                         period,
                                                         t_obj.tipper[nty, 0, 1].real,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                        ertx = mtplottools.plot_errorbar(axti,
                                                         period,
                                                         t_obj.tipper[ntx, 0, 0].imag,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axti,
                                                         period,
                                                         t_obj.tipper[nty, 0, 1].imag,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                    if plot_tipper == False:
                        ax_list = [axrxy, axpxy]
                        line_list = [erxy[0], eryx[0]]
                        label_list = ['$Z_{xy}$', '$Z_{yx}$']
                    else:
                        ax_list = [axrxy, axpxy, axtr, axti]
                        line_list = [[erxy[0], eryx[0]],
                                     [ertx[0], erty[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$T_{x}$', '$T_{y}$']]

                elif self.plot_component == 4:
                    if plot_tipper == False:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)

                        axrxx = fig.add_subplot(gs[0, 2:], sharex=axrxy)
                        axpxx = fig.add_subplot(gs[1, 2:], sharex=axrxy)
                    else:
                        axrxy = fig.add_subplot(gs[0, 0:2])
                        axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)

                        axrxx = fig.add_subplot(gs[0, 2:4], sharex=axrxy)
                        axpxx = fig.add_subplot(gs[1, 2:4], sharex=axrxy)

                        axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)

                    if self.plot_z == False:
                        # plot resistivity
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
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzyx],
                                                         rp.resyx[nzyx],
                                                         rp.resyx_err[nzyx],
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axrxx,
                                                         period[nzyy],
                                                         rp.resyy[nzyy],
                                                         rp.resyy_err[nzyy],
                                                         **kw_yy)
                        # plot phase
                        erxx = mtplottools.plot_errorbar(axpxx,
                                                         period[nzxx],
                                                         rp.phasexx[nzxx],
                                                         rp.phasexx_err[nzxx],
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         rp.phasexy[nzxy],
                                                         rp.phasexy_err[nzxy],
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         rp.phaseyx[nzyx],
                                                         rp.phaseyx_err[nzyx],
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpxx,
                                                         period[nzyy],
                                                         rp.phaseyy[nzyy],
                                                         rp.phaseyy_err[nzyy],
                                                         **kw_yy)
                    elif self.plot_z == True:
                        # plot real
                        erxx = mtplottools.plot_errorbar(axrxx,
                                                         period[nzxx],
                                                         abs(z_obj.z[nzxx, 0, 0].real),
                                                         abs(z_obj.z_err[nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].real),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axrxx,
                                                         period[nzyy],
                                                         abs(z_obj.z[nzyy, 1, 1].real),
                                                         abs(z_obj.z_err[nzyy, 1, 1].real),
                                                         **kw_yy)
                        # plot phase
                        erxx = mtplottools.plot_errorbar(axpxx,
                                                         period[nzxx],
                                                         abs(z_obj.z[nzxx, 0, 0].imag),
                                                         abs(z_obj.z_err[nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpxx,
                                                         period[nzyy],
                                                         abs(z_obj.z[nzyy, 1, 1].imag),
                                                         abs(z_obj.z_err[nzyy, 1, 1].real),
                                                         **kw_yy)
                    # plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr,
                                                         period[ntx],
                                                         t_obj.tipper[ntx, 0, 0].real,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr,
                                                         period[nty],
                                                         t_obj.tipper[nty, 0, 1].real,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                        ertx = mtplottools.plot_errorbar(axti,
                                                         period[ntx],
                                                         t_obj.tipper[ntx, 0, 0].imag,
                                                         t_obj.tipper_err[ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axti,
                                                         period[nty],
                                                         t_obj.tipper[nty, 0, 1].imag,
                                                         t_obj.tipper_err[nty, 0, 1],
                                                         **kw_yy)

                    if plot_tipper == False:
                        ax_list = [axrxy, axrxx, axpxy, axpxx]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$Z_{xx}$', '$Z_{yy}$']]
                    else:
                        ax_list = [axrxy, axrxx, axpxy, axpxx, axtr, axti]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]],
                                     [ertx[0]], erty[0]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$Z_{xx}$', '$Z_{yy}$'],
                                      ['$T_x$', '$T_y$']]

                # set axis properties
                for aa, ax in enumerate(ax_list):
                    ax.tick_params(axis='y', pad=self.ylabel_pad)
                    #                    ylabels = ax.get_yticks().tolist()
                    #                    ylabels[-1] = ''
                    #                    ylabels[0] = ''
                    #                    ax.set_yticklabels(ylabels)
                    if len(ax_list) == 2:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if self.plot_z == True:
                            ax.set_yscale('log')
                            ylim = ax.get_ylim()
                            ylimits = (10 ** np.floor(np.log10(ylim[0])),
                                       10 ** np.ceil(np.log10(ylim[1])))
                            ax.set_ylim(ylimits)
                            ylabels = [' '] + \
                                      [mtplottools.labeldict[ii] for ii
                                       in np.arange(np.log10(ylimits[0]),
                                                    np.log10(ylimits[1]), 1)] + \
                                      [' ']
                            ax.set_yticklabels(ylabels)
                        if aa == 0:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Re[Z (mV/km nT)]|',
                                              fontdict=fontdict)
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('|Im[Z (mV/km nT)]|',
                                              fontdict=fontdict)
                    elif len(ax_list) == 4 and plot_tipper == False:
                        if self.plot_z == True:
                            ax.set_yscale('log')
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            if self.plot_z == False:
                                ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
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

                    elif len(ax_list) == 4 and plot_tipper == True:
                        if aa == 0 or aa == 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        else:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res. ($\mathbf{\Omega \cdot m}$)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Re[Z (mV/km nT)]',
                                              fontdict=fontdict)
                        elif aa == 1:
                            if self.plot_z == False:
                                ax.set_ylabel('Phase (deg)',
                                              fontdict=fontdict)
                            elif self.plot_z == True:
                                ax.set_ylabel('Im[Z (mV/km nT)]',
                                              fontdict=fontdict)
                        if aa <= 2:
                            ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
                            if self.plot_z == True:
                                ax.set_yscale('log')
                                #                        else:
                                #                            plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    ax.set_xscale('log')
                    ax.set_xlim(xmin=10 ** (np.floor(np.log10(period[0]))) * 1.01,
                                xmax=10 ** (np.ceil(np.log10(period[-1]))) * .99)
                    ax.grid(True, alpha=.25)

            if plotr == True:
                for rr in range(nr):
                    if self.color_mode == 'color':
                        cxy = (0, .4 + float(rr) / (3 * nr), 0)
                        cyx = (.7 + float(rr) / (4 * nr), .13, .63 - float(rr) / (4 * nr))
                    elif self.color_mode == 'bw':
                        cxy = tuple(3 * [1 - .5 / (rr + 1)])
                        cyx = tuple(3 * [1 - .5 / (rr + 1)])

                    resp_z_obj = self.resp_object[rr].mt_dict[station].Z
                    resp_z_err = np.nan_to_num((z_obj.z - resp_z_obj.z) / z_obj.z_err)

                    resp_t_obj = self.resp_object[rr].mt_dict[station].Tipper
                    resp_t_err = np.nan_to_num((t_obj.tipper - resp_t_obj.tipper) /
                                               t_obj.tipper_err)

                    rrp = mtplottools.ResPhase(resp_z_obj)

                    rms = resp_z_err.std()
                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()
                    rms_tx = resp_t_err[:, 0, 0].std()
                    rms_ty = resp_t_err[:, 0, 1].std()
                    print ' --- response {0} ---'.format(rr)
                    print '  RMS = {:.2f}'.format(rms)
                    print '      RMS_xx = {:.2f}'.format(rms_xx)
                    print '      RMS_xy = {:.2f}'.format(rms_xy)
                    print '      RMS_yx = {:.2f}'.format(rms_yx)
                    print '      RMS_yy = {:.2f}'.format(rms_yy)
                    print '      RMS_Tx = {:.2f}'.format(rms_tx)
                    print '      RMS_Ty = {:.2f}'.format(rms_ty)

                    # --> make key word dictionaries for plotting
                    kw_xx = {'color': cxy,
                             'marker': self.mtem,
                             'ms': self.ms,
                             'ls': ':',
                             'lw': self.lw,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    kw_yy = {'color': cyx,
                             'marker': self.mtmm,
                             'ms': self.ms,
                             'ls': ':',
                             'lw': self.lw,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    if self.plot_style == 1:
                        if self.plot_component == 2:
                            if self.plot_z == False:
                                # plot resistivity
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  rrp.resxy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx,
                                                                  period[nzyx],
                                                                  rrp.resyx[nzyx],
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  rrp.phasexy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx,
                                                                  period[nzyx],
                                                                  rrp.phaseyx[nzyx],
                                                                  **kw_yy)
                            elif self.plot_z == True:
                                # plot real
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].imag),
                                                                  **kw_yy)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].real,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].real,
                                                                  **kw_yy)

                                rertx = mtplottools.plot_errorbar(axti,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].imag,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].imag,
                                                                  **kw_yy)
                            if plot_tipper == False:
                                line_list[0] += [rerxy[0]]
                                line_list[1] += [reryx[0]]
                                label_list[0] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy)]
                                label_list[1] += ['$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                            else:
                                line_list[0] += [rerxy[0]]
                                line_list[1] += [reryx[0]]
                                line_list[2] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy)]
                                label_list[1] += ['$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[2] += ['$T^m_{x}$' +
                                                  'rms={0:.2f}'.format(rms_tx),
                                                  '$T^m_{y}$' +
                                                  'rms={0:.2f}'.format(rms_ty)]
                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                # plot resistivity
                                rerxx = mtplottools.plot_errorbar(axrxx,
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
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
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
                                # plot real
                                rerxx = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[nzxx, 0, 0].real),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axryx,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].real),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axryy,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[nzyy, 1, 1].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[nzxx, 0, 0].imag),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpyx,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].imag),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpyy,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[nzyy, 1, 1].imag),
                                                                  **kw_yy)
                            if plot_tipper == True:
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

                            if plot_tipper == False:
                                line_list[0] += [rerxx[0]]
                                line_list[1] += [rerxy[0]]
                                line_list[2] += [reryx[0]]
                                line_list[3] += [reryy[0]]
                                label_list[0] += ['$Z^m_{xx}$ ' +
                                                  'rms={0:.2f}'.format(rms_xx)]
                                label_list[1] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy)]
                                label_list[2] += ['$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[3] += ['$Z^m_{yy}$ ' +
                                                  'rms={0:.2f}'.format(rms_yy)]
                            else:
                                line_list[0] += [rerxx[0]]
                                line_list[1] += [rerxy[0]]
                                line_list[2] += [reryx[0]]
                                line_list[3] += [reryy[0]]
                                line_list[4] += [rertx[0]]
                                line_list[5] += [rerty[0]]
                                label_list[0] += ['$Z^m_{xx}$ ' +
                                                  'rms={0:.2f}'.format(rms_xx)]
                                label_list[1] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy)]
                                label_list[2] += ['$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[3] += ['$Z^m_{yy}$ ' +
                                                  'rms={0:.2f}'.format(rms_yy)]
                                label_list[4] += ['$T^m_{x}$' +
                                                  'rms={0:.2f}'.format(rms_tx)]
                                label_list[5] += ['$T^m_{y}$' +
                                                  'rms={0:.2f}'.format(rms_ty)]

                    elif self.plot_style == 2:
                        if self.plot_component == 2:
                            if self.plot_z == False:
                                # plot resistivity
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  rrp.resxy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  rrp.resyx[nzyx],
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  rrp.phasexy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  rrp.phaseyx[nzyx],
                                                                  **kw_yy)
                            elif self.plot_z == True:
                                # plot real
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].imag),
                                                                  **kw_xx)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].real,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].real,
                                                                  **kw_yy)

                                rertx = mtplottools.plot_errorbar(axti,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].imag,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].imag,
                                                                  **kw_yy)

                            if plot_tipper == False:
                                line_list += [rerxy[0], reryx[0]]
                                label_list += ['$Z^m_{xy}$ ' +
                                               'rms={0:.2f}'.format(rms_xy),
                                               '$Z^m_{yx}$ ' +
                                               'rms={0:.2f}'.format(rms_yx)]
                            else:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy),
                                                  '$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$T^m_{x}$' +
                                                  'rms={0:.2f}'.format(rms_tx),
                                                  '$T^m_{y}$' +
                                                  'rms={0:.2f}'.format(rms_ty)]

                        elif self.plot_component == 4:
                            if self.plot_z == False:
                                # plot resistivity
                                rerxx = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzxx],
                                                                  rrp.resxx[nzxx],
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  rrp.resxy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  rrp.resyx[nzyx],
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzyy],
                                                                  rrp.resyy[nzyy],
                                                                  **kw_yy)
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzxx],
                                                                  rrp.phasexx[nzxx],
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  rrp.phasexy[nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  rrp.phaseyx[nzyx],
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzyy],
                                                                  rrp.phaseyy[nzyy],
                                                                  **kw_yy)
                            elif self.plot_z == True:
                                # plot real
                                rerxx = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[nzxx, 0, 0].real),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].real),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[nzyy, 1, 1].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[nzxx, 0, 0].imag),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[nzyx, 1, 0].imag),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[nzyy, 1, 1].imag),
                                                                  **kw_yy)

                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].real,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].real,
                                                                  **kw_yy)

                                rertx = mtplottools.plot_errorbar(axti,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[ntx, 0, 0].imag,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[nty, 0, 1].imag,
                                                                  **kw_yy)

                            if plot_tipper == False:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rerxx[0], reryy[0]]
                                label_list[0] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy),
                                                  '$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$Z^m_{xx}$ ' +
                                                  'rms={0:.2f}'.format(rms_xx),
                                                  '$Z^m_{yy}$ ' +
                                                  'rms={0:.2f}'.format(rms_yy)]
                            else:
                                line_list[0] += [rerxy[0], reryx[0]]
                                line_list[1] += [rerxx[0], reryy[0]]
                                line_list[2] += [rertx[0], rerty[0]]
                                label_list[0] += ['$Z^m_{xy}$ ' +
                                                  'rms={0:.2f}'.format(rms_xy),
                                                  '$Z^m_{yx}$ ' +
                                                  'rms={0:.2f}'.format(rms_yx)]
                                label_list[1] += ['$Z^m_{xx}$ ' +
                                                  'rms={0:.2f}'.format(rms_xx),
                                                  '$Z^m_{yy}$ ' +
                                                  'rms={0:.2f}'.format(rms_yy)]
                                label_list[2] += ['$T^m_{x}$' +
                                                  'rms={0:.2f}'.format(rms_tx),
                                                  '$T^m_{y}$' +
                                                  'rms={0:.2f}'.format(rms_ty)]

            # make legends
            if self.plot_style == 1:
                legend_ax_list = ax_list[0:self.plot_component]
                if plot_tipper == True:
                    if self.plot_component == 2:
                        legend_ax_list.append(ax_list[4])
                    elif self.plot_component == 4:
                        legend_ax_list.append(ax_list[8])
                        legend_ax_list.append(ax_list[10])
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
                              prop={'size': max([self.font_size / (nr + 1), 5])})
            if self.plot_style == 2:
                if self.plot_component == 2:
                    legend_ax_list = [ax_list[0]]
                    if plot_tipper == True:
                        legend_ax_list.append(ax_list[2])
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
                                  prop={'size': max([self.font_size / (nr + 1), 5])})
                else:
                    legend_ax_list = ax_list[0:self.plot_component / 2]
                    if plot_tipper == True:
                        if self.plot_component == 2:
                            legend_ax_list.append(ax_list[2])
                        elif self.plot_component == 4:
                            legend_ax_list.append(ax_list[4])
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
                                  prop={'size': max([self.font_size / (nr + 1), 5])})

        ##--> BE SURE TO SHOW THE PLOT
        plt.show()

        if save2file is not None:
            fig.savefig(save2file, dpi=self.fig_dpi, bbox_inches='tight')
            # self.save_figure(save2file)
        return plt

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotAllResponses()
            >>> #change line width
            >>> p1.lw = 2
            >>> p1.redraw_plot()
        """
        for fig in self.fig_list:
            plt.close(fig)
        self.plot()

    def save_figure(self, save_fn, file_format='png', orientation='portrait',
                    fig_dpi=None, close_fig='y'):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_PhaseTensor.file_format

                          * full path -> file will be save to the given
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps...

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at
                          which the figure was made.  I don't think that
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotPseudoSection()
            >>> ps1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        fig = plt.gcf()
        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        else:
            save_fn = os.path.join(save_fn, '_L2.' +
                                   file_format)
            fig.savefig(save_fn, dpi=fig_dpi, format=file_format,
                        orientation=orientation, bbox_inches='tight')

        if close_fig == 'y':
            plt.clf()
            plt.close(fig)
        else:
            pass

        self.fig_fn = save_fn
        print 'Saved figure to: ' + self.fig_fn

        return self.fig_fn

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.modeling.occam2d as occam2d
            >>> dfn = r"/home/occam2d/Inv1/data.dat"
            >>> ocd = occam2d.Occam2DData(dfn)
            >>> ps1 = ocd.plotAllResponses()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return ("Plots data vs model response computed by WS3DINV")
