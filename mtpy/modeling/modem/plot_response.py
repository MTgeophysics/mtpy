"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""

import numpy as np
import os
from matplotlib import pyplot as plt, gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter,LogFormatterSciNotation
from mtpy.imaging import mtplottools as mtplottools
from mtpy.modeling.modem import Data, Residual
from mtpy.core.z import Z, Tipper
import sys

__all__ = ['PlotResponse']


class PlotResponse(object):
    """
    plot data and response

    Plots the real and imaginary impedance and induction vector if present.

    :Example: ::

        >>> import mtpy.modeling.modem as modem
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
    cted                     color for data Z_XX and Z_XY mode
    ctem                     color for model Z_XX and Z_XY mode
    ctmd                     color for data Z_YX and Z_YY mode
    ctmm                     color for model Z_YX and Z_YY mode
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
    lw                       line width data curves (*default* is .5)
    ms                       size of markers (*default* is 1.5)
    lw_r                     line width response curves (*default* is .5)
    ms_r                     size of markers response curves (*default* is 1.5)
    mted                     marker for data Z_XX and Z_XY mode
    mtem                     marker for model Z_XX and Z_XY mode
    mtmd                     marker for data Z_YX and Z_YY mode
    mtmm                     marker for model Z_YX and Z_YY mode
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
        self.ms_r = kwargs.pop('ms_r', 3)
        self.lw = kwargs.pop('lw', .5)
        self.lw_r = kwargs.pop('lw_r', 1.0)
        self.ls = kwargs.pop('ls',':')
        self.e_capthick = kwargs.pop('e_capthick', .5)
        self.e_capsize = kwargs.pop('e_capsize', 2)

        self.plot_style = kwargs.pop('plot_style', 1)
        if self.plot_style not in [1, 2, 3]:
            print(("self.plot_style = %s. It MUST be either 1 (default; 4 column figures) or 2 (2 column figures) or 3 (1 column figures)" % str(self.plot_style)))
            self.plot_style = 1

        # color mode
        if self.color_mode == 'color':
            # color for data
            self.cted = kwargs.pop('cted', (0, 0, 1))
            self.ctmd = kwargs.pop('ctmd', (1, 0, 0))
            self.mted = kwargs.pop('mted', 's')
            self.mtmd = kwargs.pop('mtmd', 'o')

            # color for occam2d model
            if self.plot_style == 3:
                # if plot_style is 3, set default color for model response to same as data
                self.ctem = kwargs.pop('ctem',self.cted)
                self.ctmm = kwargs.pop('ctmm',self.ctmd)
            else:
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
        self.phase_limits_d = kwargs.pop('phase_limits_d', None)
        self.res_limits_d = kwargs.pop('res_limits_d', None)
        self.res_limits_od = kwargs.pop('res_limits_od', None)
        self.tipper_limits = kwargs.pop('tipper_limits', None)
        self.period_limits = kwargs.pop('period_limits', None)

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
        self.legend_pos = (.5, 1.18)
        self.legend_pos_tipper = (.5, 1.18)
        self.legend_marker_scale = 1
        self.legend_border_axes_pad = .01
        self.legend_label_spacing = 0.07
        self.legend_handle_text_pad = .2
        self.legend_border_pad = .15

        self.font_size = kwargs.pop('font_size', 6)

        self.plot_type = kwargs.pop('plot_type', '1')

        self.plot_component = kwargs.pop('plot_component', 4)
        self.plot_yn = kwargs.pop('plot_yn', 'y')
        self.save_plots = kwargs.pop('plot_yn', False)
        self.plot_z = kwargs.pop('plot_z', True)
        self.ylabel_pad = kwargs.pop('ylabel_pad', 1.25)
        self.label_axes = kwargs.pop('label_axes',True)
        self.shift_yx_phase = kwargs.pop('shift_yx_phase',False)

        self.fig_list = []


        # if self.plot_yn == 'y':
        #     self.plot()

    def plot(self):
        if self.plot_style == 1: # and has tipper data
            self._plot()
        if self.plot_style == 2:
            self._plot_2col()
        if self.plot_style == 3:
            self._plot_1col()


    def _read_files(self):
        
        self.data_object = Data()
        self.data_object.read_data_file(self.data_fn)


        # read in response files
        if self.resp_fn is not None:
            self.resp_object = []
            if not isinstance(self.resp_fn, list):
                resp_obj = Data()
                resp_obj.read_data_file(self.resp_fn)
                self.resp_object = [resp_obj]
            else:
                for rfile in self.resp_fn:
                    resp_obj = Data()
                    resp_obj.read_data_file(rfile)
                    self.resp_object.append(resp_obj)



    def _plot(self):
        """
        plot as an internal function of this class
        """

        self._read_files()
        
        # get shape of impedance tensors
        ns = len(list(self.data_object.mt_dict.keys()))

        # get number of response files
        nr = len(self.resp_object)

        if type(self.plot_type) is list:
            ns = len(self.plot_type)

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}
        if self.plot_z:
            h_ratio = [1, 1, .5]
        elif not self.plot_z:
            h_ratio = [1.5, 1, .5]

        ax_list = []
        line_list = []
        label_list = []

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
            pstation_list = list(self.data_object.mt_dict.keys())

        for jj, station in enumerate(pstation_list):
            z_obj = self.data_object.mt_dict[station].Z
            t_obj = self.data_object.mt_dict[station].Tipper
            period = self.data_object.period_list
            print('Plotting: {0}'.format(station))


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


            # convert to apparent resistivity and phase
            z_obj.compute_resistivity_phase()

            # find locations where points have been masked
            nzxx = np.nonzero(z_obj.z[:, 0, 0])[0]
            nzxy = np.nonzero(z_obj.z[:, 0, 1])[0]
            nzyx = np.nonzero(z_obj.z[:, 1, 0])[0]
            nzyy = np.nonzero(z_obj.z[:, 1, 1])[0]
            ntx = np.nonzero(t_obj.tipper[:, 0, 0])[0]
            nty = np.nonzero(t_obj.tipper[:, 0, 1])[0]

            # convert to apparent resistivity and phase
            if self.plot_z:
                scaling = np.zeros_like(z_obj.z)
                for ii in range(2):
                    for jj in range(2):
                        scaling[:, ii, jj] = 1. / np.sqrt(z_obj.freq)
                plot_res = abs(z_obj.z.real * scaling)
                plot_res_err = abs(z_obj.z_err * scaling)
                plot_phase = abs(z_obj.z.imag * scaling)
                plot_phase_err = abs(z_obj.z_err * scaling)
                h_ratio = [1, 1, .5]

            elif not self.plot_z:
                plot_res = z_obj.resistivity
                plot_res_err = z_obj.resistivity_err
                plot_phase = z_obj.phase
                plot_phase_err = z_obj.phase_err
                h_ratio = [1.5, 1, .5]

                try:
                    self.res_limits_d = (10 ** (np.floor(np.log10(min([plot_res[nzxx, 0, 0].min(),
                                                                       plot_res[nzyy, 1, 1].min()])))),
                                         10 ** (np.ceil(np.log10(max([plot_res[nzxx, 0, 0].max(),
                                                                      plot_res[nzyy, 1, 1].max()])))))
                except ValueError:
                    self.res_limits_d = None
                try:
                    self.res_limits_od = (10 ** (np.floor(np.log10(min([plot_res[nzxy, 0, 1].min(),
                                                                        plot_res[nzyx, 1, 0].min()])))),
                                          10 ** (np.ceil(np.log10(max([plot_res[nzxy, 0, 1].max(),
                                                                       plot_res[nzyx, 1, 0].max()])))))
                except ValueError:
                    self.res_limits_od = None

            # make figure
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            plt.clf()
            fig.suptitle(str(station), fontdict=fontdict)

            # set the grid of subplots
            if np.all(t_obj.tipper == 0.0) == True:
                self.plot_tipper = False
                gs = gridspec.GridSpec(2, 4,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio[:2])

            else:
                self.plot_tipper = True
                self.tipper_limits = (np.round(min([t_obj.tipper[ntx, 0, 0].real.min(),
                                                    t_obj.tipper[nty, 0, 1].real.min(),
                                                    t_obj.tipper[ntx, 0, 0].imag.min(),
                                                    t_obj.tipper[nty, 0, 1].imag.min()]),
                                               1),
                                      np.round(max([t_obj.tipper[ntx, 0, 0].real.max(),
                                                    t_obj.tipper[nty, 0, 1].real.max(),
                                                    t_obj.tipper[ntx, 0, 0].imag.max(),
                                                    t_obj.tipper[nty, 0, 1].imag.max()]),
                                               1))

                gs = gridspec.GridSpec(3, 4,
                                   wspace=self.subplot_wspace,
                                   left=self.subplot_left,
                                   top=self.subplot_top,
                                   bottom=self.subplot_bottom,
                                   right=self.subplot_right,
                                   hspace=self.subplot_hspace,
                                   height_ratios=h_ratio)

            axrxx = fig.add_subplot(gs[0, 0])
            axrxy = fig.add_subplot(gs[0, 1], sharex=axrxx)
            axryx = fig.add_subplot(gs[0, 2], sharex=axrxx, sharey=axrxy)
            axryy = fig.add_subplot(gs[0, 3], sharex=axrxx, sharey=axrxx)

            axpxx = fig.add_subplot(gs[1, 0])
            axpxy = fig.add_subplot(gs[1, 1], sharex=axrxx)
            axpyx = fig.add_subplot(gs[1, 2], sharex=axrxx)
            axpyy = fig.add_subplot(gs[1, 3], sharex=axrxx)

            if self.plot_tipper == True:
                axtxr = fig.add_subplot(gs[2, 0], sharex=axrxx)
                axtxi = fig.add_subplot(gs[2, 1], sharex=axrxx, sharey=axtxr)
                axtyr = fig.add_subplot(gs[2, 2], sharex=axrxx)
                axtyi = fig.add_subplot(gs[2, 3], sharex=axrxx, sharey=axtyr)

                self.ax_list = [axrxx, axrxy, axryx, axryy,
                                axpxx, axpxy, axpyx, axpyy,
                                axtxr, axtxi, axtyr, axtyi]
            else:

                self.ax_list = [axrxx, axrxy, axryx, axryy,
                                axpxx, axpxy, axpyx, axpyy]


            # ---------plot the apparent resistivity-----------------------------------
            # plot each component in its own subplot
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
            # plot phase
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

            # plot tipper
            if self.plot_tipper:
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


            print(("self.plot_tipper = {}".format(self.plot_tipper)))
            # ----------------------------------------------
            # get error bar list for editing later
            if not self.plot_tipper:
                try:
                    self._err_list = [[erxx[1][0], erxx[1][1], erxx[2][0]],
                                      [erxy[1][0], erxy[1][1], erxy[2][0]],
                                      [eryx[1][0], eryx[1][1], eryx[2][0]],
                                      [eryy[1][0], eryy[1][1], eryy[2][0]]]
                    line_list = [[erxx[0]], [erxy[0]], [eryx[0]], [eryy[0]]]
                except IndexError:
                    print('Found no Z components for {0}'.format(self.station))
                    line_list = [[None], [None],
                                 [None], [None]]

                    self._err_list = [[None, None, None],
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
                                      [ertx[1][0], ertx[1][1], ertx[2][0]],
                                      [erty[1][0], erty[1][1], erty[2][0]]]
                except IndexError:
                    print('Found no Z components for {0}'.format(station))
                    line_list = [[None], [None],
                                 [None], [None],
                                 [None], [None]]

                    self._err_list = [[None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None],
                                      [None, None, None]]
            # ------------------------------------------
            # make things look nice
            # set titles of the Z components
            label_list = [['$Z_{xx}$'], ['$Z_{xy}$'],
                          ['$Z_{yx}$'], ['$Z_{yy}$']]
            for ax, label in zip(self.ax_list[0:4], label_list):
                ax.set_title(label[0], fontdict={'size': self.font_size + 2,
                                                 'weight': 'bold'})

                # set legends for tipper components
            # fake a line
            l1 = plt.Line2D([0], [0], linewidth=0, color='w', linestyle='None',
                            marker='.')
            t_label_list = ['Re{$T_x$}', 'Im{$T_x$}', 'Re{$T_y$}', 'Im{$T_y$}']
            label_list += [['$T_{x}$'], ['$T_{y}$']]
            if self.plot_tipper:
                for ax, label in zip(self.ax_list[-4:], t_label_list):
                    ax.legend([l1], [label], loc='upper left',
                              markerscale=.01,
                              borderaxespad=.05,
                              labelspacing=.01,
                              handletextpad=.05,
                              borderpad=.05,
                              prop={'size': max([self.font_size, 6])})



                # set axis properties
            for aa, ax in enumerate(self.ax_list):
                ax.tick_params(axis='y', pad=self.ylabel_pad)
                if self.plot_tipper==False:
                    if aa < 4:
                        if self.plot_z == True:
                            ax.set_yscale('log', nonposy='clip')
                    else:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)

                if aa < 8:
                    #                    ylabels[-1] = ''
                    #                    ylabels[0] = ''
                    #                    ax.set_yticklabels(ylabels)
                    #                    plt.setp(ax.get_xticklabels(), visible=False)
                    if self.plot_z == True:
                        ax.set_yscale('log', nonposy='clip')

                else:
                    ax.set_xlabel('Period (s)', fontdict=fontdict)

                if aa < 4 and self.plot_z is False:
                    ylabels = ax.get_yticklabels()
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    ax.set_yscale('log', nonposy='clip')
                    if aa == 0 or aa == 3:
                        ax.set_ylim(self.res_limits_d)
                    elif aa == 1 or aa == 2:
                        ax.set_ylim(self.res_limits_od)

                if aa > 3 and aa < 8 and self.plot_z is False:
                    #ax.yaxis.set_major_locator(MultipleLocator(10.0))
                    if self.phase_limits_d is not None:
                        ax.set_ylim(self.phase_limits_d)
                # set axes labels
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
                    ax.yaxis.set_major_locator(MultipleLocator(.1))
                    if self.tipper_limits is not None:
                        ax.set_ylim(self.tipper_limits)
                    else:
                        pass

                ax.set_xscale('log', nonposx='clip')
                # set period limits
                if self.period_limits is None:
                    self.period_limits = (10 ** (np.floor(np.log10(period[0]))) * 1.01,
                                          10 ** (np.ceil(np.log10(period[-1]))) * .99)
                ax.set_xlim(xmin=self.period_limits[0],
                            xmax=self.period_limits[1])
                ax.grid(True, alpha=.25)

                ylabels = ax.get_yticks().tolist()
                if aa < 8:
                    ylabels[-1] = ''
                    ylabels[0] = ''
                    ax.set_yticklabels(ylabels)
                    plt.setp(ax.get_xticklabels(), visible=False)

            # plot model response
            if self.resp_object is not None:
                for resp_obj in self.resp_object:
                    resp_z_obj = resp_obj.mt_dict[station].Z
                    resp_z_err = np.nan_to_num((z_obj.z - resp_z_obj.z) / z_obj.z_err)
                    resp_z_obj.compute_resistivity_phase()

                    resp_t_obj = resp_obj.mt_dict[station].Tipper
                    resp_t_err = np.nan_to_num((t_obj.tipper - resp_t_obj.tipper) / t_obj.tipper_err)

                    # convert to apparent resistivity and phase
                    if self.plot_z == True:
                        scaling = np.zeros_like(resp_z_obj.z)
                        for ii in range(2):
                            for jj in range(2):
                                scaling[:, ii, jj] = 1. / np.sqrt(resp_z_obj.freq)
                        r_plot_res = abs(resp_z_obj.z.real * scaling)
                        r_plot_phase = abs(resp_z_obj.z.imag * scaling)

                    elif self.plot_z == False:
                        r_plot_res = resp_z_obj.resistivity
                        r_plot_phase = resp_z_obj.phase

                    rms_xx = resp_z_err[:, 0, 0].std()
                    rms_xy = resp_z_err[:, 0, 1].std()
                    rms_yx = resp_z_err[:, 1, 0].std()
                    rms_yy = resp_z_err[:, 1, 1].std()

                    # --> make key word dictionaries for plotting
                    kw_xx = {'color': self.ctem,
                             'marker': self.mtem,
                             'ms': self.ms_r,
                             'ls': ':',
                             'lw': self.lw_r,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    kw_yy = {'color': self.ctmm,
                             'marker': self.mtmm,
                             'ms': self.ms_r,
                             'ls': ':',
                             'lw': self.lw_r,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    # plot data response
                    rerxx = mtplottools.plot_errorbar(axrxx,
                                                      period[nzxx],
                                                      r_plot_res[nzxx, 0, 0],
                                                      None,
                                                      **kw_xx)
                    rerxy = mtplottools.plot_errorbar(axrxy,
                                                      period[nzxy],
                                                      r_plot_res[nzxy, 0, 1],
                                                      None,
                                                      **kw_xx)
                    reryx = mtplottools.plot_errorbar(axryx,
                                                      period[nzyx],
                                                      r_plot_res[nzyx, 1, 0],
                                                      None,
                                                      **kw_yy)
                    reryy = mtplottools.plot_errorbar(axryy,
                                                      period[nzyy],
                                                      r_plot_res[nzyy, 1, 1],
                                                      None,
                                                      **kw_yy)
                    # plot phase
                    repxx = mtplottools.plot_errorbar(axpxx,
                                                      period[nzxx],
                                                      r_plot_phase[nzxx, 0, 0],
                                                      None,
                                                      **kw_xx)
                    repxy = mtplottools.plot_errorbar(axpxy,
                                                      period[nzxy],
                                                      r_plot_phase[nzxy, 0, 1],
                                                      None,
                                                      **kw_xx)
                    repyx = mtplottools.plot_errorbar(axpyx,
                                                      period[nzyx],
                                                      r_plot_phase[nzyx, 1, 0],
                                                      None,
                                                      **kw_yy)
                    repyy = mtplottools.plot_errorbar(axpyy,
                                                      period[nzyy],
                                                      r_plot_phase[nzyy, 1, 1],
                                                      None,
                                                      **kw_yy)

                    # plot tipper
                    if self.plot_tipper == True:
                        rertx = mtplottools.plot_errorbar(axtxr,
                                                          period[ntx],
                                                          resp_t_obj.tipper[ntx, 0, 0].real,
                                                          None,
                                                          **kw_xx)
                        rerty = mtplottools.plot_errorbar(axtyr,
                                                          period[nty],
                                                          resp_t_obj.tipper[nty, 0, 1].real,
                                                          None,
                                                          **kw_yy)

                        reptx = mtplottools.plot_errorbar(axtxi,
                                                          period[ntx],
                                                          resp_t_obj.tipper[ntx, 0, 0].imag,
                                                          None,
                                                          **kw_xx)
                        repty = mtplottools.plot_errorbar(axtyi,
                                                          period[nty],
                                                          resp_t_obj.tipper[nty, 0, 1].imag,
                                                          None,
                                                          **kw_yy)

                    if self.plot_tipper == False:
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
                        label_list[4] += ['$T^m_{x}$ ' +
                                          'rms={0:.2f}'.format(resp_t_err[:, 0, 0].std())]
                        label_list[5] += ['$T^m_{y}$' +
                                          'rms={0:.2f}'.format(resp_t_err[:, 0, 1].std())]

                legend_ax_list = self.ax_list[0:4]
                #                if self.plot_tipper == True:
                #                    legend_ax_list += [self.ax_list[-4], self.ax_list[-2]]

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
                              framealpha=1,
                              prop={'size': max([self.font_size, 5])})

            plt.show()
            
            if self.save_plots:
                save_filename = os.path.join(os.path.dirname(self.data_fn),station+'.png')
                self.save_figure(save_filename,fig_dpi=self.fig_dpi)
            

    def _plot_2col(self):

        """
        Internal method to plot responses in 2 columns of subplots.
        Show the figure and optionally save to a file named save2file
        """

        self._read_files()

        # get shape of impedance tensors
        ns = len(list(self.data_object.mt_dict.keys()))

        # get number of response files
        nr = len(self.resp_object)

        if isinstance(self.plot_type, list):
            ns = len(self.plot_type)

        # --> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}
        if self.plot_z == True:
            h_ratio = [1, 1]
        elif self.plot_z == False:
            # h_ratio = [2, 1.5]
            h_ratio = [2.0, 1.5, 0.75]

        self.ax_list = []
        line_list = []
        label_list = []


        if self.plot_type != '1':
            pstation_list = []
            if not isinstance(self.plot_type, list):
                self.plot_type = [self.plot_type]
            for ii, station in enumerate(self.data_object.mt_dict.keys()):
                if not isinstance(station, int):
                    for pstation in self.plot_type:
                        if station.find(str(pstation)) >= 0:
                            pstation_list.append(station)
                else:
                    for pstation in self.plot_type:
                        if station == int(pstation):
                            pstation_list.append(ii)
        else:
            pstation_list = list(self.data_object.mt_dict.keys())

        for jj, station in enumerate(pstation_list):
            z_obj = self.data_object.mt_dict[station].Z
            t_obj = self.data_object.mt_dict[station].Tipper
            period = self.data_object.period_list


            # --> make key word dictionaries for plotting
            kw_xx = {'color': self.cted,
                     'marker': self.mted,
                     'ms': self.ms,
                     'ls': self.ls,
                     'lw': self.lw,
                     'e_capsize': self.e_capsize,
                     'e_capthick': self.e_capthick}
    
            kw_yy = {'color': self.ctmd,
                     'marker': self.mtmd,
                     'ms': self.ms,
                     'ls': self.ls,
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

            if self.resp_fn is not None:
                plotr = True
            else:
                plotr = False

            # make figure
            fig = plt.figure(station, self.fig_size, dpi=self.fig_dpi)
            self.fig_list.append(fig)
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

                # gs = gridspec.GridSpec(2, 6,
                #                        wspace=self.subplot_wspace,
                #                        left=self.subplot_left,
                #                        top=self.subplot_top,
                #                        bottom=self.subplot_bottom,
                #                        right=self.subplot_right,
                #                        hspace=self.subplot_hspace,
                #                        height_ratios=h_ratio)
                if len(h_ratio) < 3 :
                    h_ratio = [2.0, 1.5, 0.75]
                gs = gridspec.GridSpec(3, 2,
                                       wspace=self.subplot_wspace,
                                       left=self.subplot_left,
                                       top=self.subplot_top,
                                       bottom=self.subplot_bottom,
                                       right=self.subplot_right,
                                       hspace=self.subplot_hspace,
                                       height_ratios=h_ratio)

            else:
                if len(h_ratio) >= 3 :
                    h_ratio = [1,1]
                gs = gridspec.GridSpec(2, 4,
                                       wspace=self.subplot_wspace,
                                       left=self.subplot_left,
                                       top=self.subplot_top,
                                       bottom=self.subplot_bottom,
                                       right=self.subplot_right,
                                       hspace=self.subplot_hspace,
                                       height_ratios=h_ratio)
            # ---------plot the apparent resistivity---------------------------
            # plot each component in its own subplot




            # plot xy and yx together and xx, yy together
            if self.plot_style == 2:
#                rp.phasexy[rp.phasexy > 180] -= 360
#                rp.phaseyx[rp.phaseyx > 180] -= 360
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
                                                         abs(z_obj.z[
                                                             nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[
                                                             nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[
                                                             nzxy, 1, 0].real),
                                                         abs(z_obj.z_err[
                                                             nzxy, 1, 0].real),
                                                         **kw_yy)
                        # plot phase
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[
                                                             nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[
                                                             nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[
                                                             nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[
                                                             nzyx, 1, 0].real),
                                                         **kw_yy)
                    # plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr,
                                                         period,
                                                         t_obj.tipper[
                                                             ntx, 0, 0].real,
                                                         t_obj.tipper_err[
                                                             ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr,
                                                         period,
                                                         t_obj.tipper[
                                                             nty, 0, 1].real,
                                                         t_obj.tipper_err[
                                                             nty, 0, 1],
                                                         **kw_yy)

                        ertx = mtplottools.plot_errorbar(axti,
                                                         period,
                                                         t_obj.tipper[
                                                             ntx, 0, 0].imag,
                                                         t_obj.tipper_err[
                                                             ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axti,
                                                         period,
                                                         t_obj.tipper[
                                                             nty, 0, 1].imag,
                                                         t_obj.tipper_err[
                                                             nty, 0, 1],
                                                         **kw_yy)
                    if plot_tipper == False:
                        self.ax_list = [axrxy, axpxy]
                        line_list = [erxy[0], eryx[0]]
                        label_list = ['$Z_{xy}$', '$Z_{yx}$']
                    else:
                        self.ax_list = [axrxy, axpxy, axtr, axti]
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
                        # axrxy = fig.add_subplot(gs[0, 0:2])
                        # axpxy = fig.add_subplot(gs[1, 0:2], sharex=axrxy)
                        #
                        #
                        # axrxx = fig.add_subplot(gs[0, 2:4], sharex=axrxy)
                        # axpxx = fig.add_subplot(gs[1, 2:4], sharex=axrxy)
                        #
                        # axtr = fig.add_subplot(gs[0, 4:], sharex=axrxy)
                        # axti = fig.add_subplot(gs[1, 4:], sharex=axrxy)

                        axrxy = fig.add_subplot(gs[0, 0])
                        axpxy = fig.add_subplot(gs[1, 0], sharex=axrxy)


                        axrxx = fig.add_subplot(gs[0, 1], sharex=axrxy)
                        axpxx = fig.add_subplot(gs[1, 1], sharex=axrxy)

                        axtr = fig.add_subplot(gs[2, 0], sharex=axrxy)
                        axti = fig.add_subplot(gs[2, 1], sharex=axrxy)



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
                                                         abs(z_obj.z[
                                                             nzxx, 0, 0].real),
                                                         abs(z_obj.z_err[
                                                             nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axrxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[
                                                             nzxy, 0, 1].real),
                                                         abs(z_obj.z_err[
                                                             nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axrxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[
                                                             nzyx, 1, 0].real),
                                                         abs(z_obj.z_err[
                                                             nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axrxx,
                                                         period[nzyy],
                                                         abs(z_obj.z[
                                                             nzyy, 1, 1].real),
                                                         abs(z_obj.z_err[
                                                             nzyy, 1, 1].real),
                                                         **kw_yy)
                        # plot phase
                        erxx = mtplottools.plot_errorbar(axpxx,
                                                         period[nzxx],
                                                         abs(z_obj.z[
                                                             nzxx, 0, 0].imag),
                                                         abs(z_obj.z_err[
                                                             nzxx, 0, 0].real),
                                                         **kw_xx)
                        erxy = mtplottools.plot_errorbar(axpxy,
                                                         period[nzxy],
                                                         abs(z_obj.z[
                                                             nzxy, 0, 1].imag),
                                                         abs(z_obj.z_err[
                                                             nzxy, 0, 1].real),
                                                         **kw_xx)
                        eryx = mtplottools.plot_errorbar(axpxy,
                                                         period[nzyx],
                                                         abs(z_obj.z[
                                                             nzyx, 1, 0].imag),
                                                         abs(z_obj.z_err[
                                                             nzyx, 1, 0].real),
                                                         **kw_yy)
                        eryy = mtplottools.plot_errorbar(axpxx,
                                                         period[nzyy],
                                                         abs(z_obj.z[
                                                             nzyy, 1, 1].imag),
                                                         abs(z_obj.z_err[
                                                             nzyy, 1, 1].real),
                                                         **kw_yy)
                    # plot tipper
                    if plot_tipper == True:
                        ertx = mtplottools.plot_errorbar(axtr,
                                                         period[ntx],
                                                         t_obj.tipper[
                                                             ntx, 0, 0].real,
                                                         t_obj.tipper_err[
                                                             ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axtr,
                                                         period[nty],
                                                         t_obj.tipper[
                                                             nty, 0, 1].real,
                                                         t_obj.tipper_err[
                                                             nty, 0, 1],
                                                         **kw_yy)

                        ertx = mtplottools.plot_errorbar(axti,
                                                         period[ntx],
                                                         t_obj.tipper[
                                                             ntx, 0, 0].imag,
                                                         t_obj.tipper_err[
                                                             ntx, 0, 0],
                                                         **kw_xx)
                        erty = mtplottools.plot_errorbar(axti,
                                                         period[nty],
                                                         t_obj.tipper[
                                                             nty, 0, 1].imag,
                                                         t_obj.tipper_err[
                                                             nty, 0, 1],
                                                         **kw_yy)



                    if plot_tipper == False:
                        self.ax_list = [axrxy, axrxx, axpxy, axpxx]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$Z_{xx}$', '$Z_{yy}$']]
                    else:
                        self.ax_list = [axrxy, axrxx, axpxy, axpxx, axtr, axti]
                        line_list = [[erxy[0], eryx[0]], [erxx[0], eryy[0]],
                                     [ertx[0], erty[0]]]
                        label_list = [['$Z_{xy}$', '$Z_{yx}$'],
                                      ['$Z_{xx}$', '$Z_{yy}$'],
                                      ['$T_x$', '$T_y$']]


                # set axis properties
                for aa, ax in enumerate(self.ax_list):
                    ax.set_xscale('log', nonposx='clip')
                    #                    ylabels = ax.get_yticks().tolist()
                    #                    ylabels[-1] = ''
                    #                    ylabels[0] = ''
                    #                    ax.set_yticklabels(ylabels)
                    if len(self.ax_list) == 2:
                        ax.set_xlabel('Period (s)', fontdict=fontdict)
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
                        if aa == 0:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
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
                    elif len(self.ax_list) == 4 and plot_tipper == False:
                        if self.plot_z == True:
                            ax.set_yscale('log', nonposy='clip')
                        if aa < 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
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

                    elif len(self.ax_list) == 4 and plot_tipper == True:
                        if aa == 0 or aa == 2:
                            plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                ax.set_yscale('log', nonposy='clip')
                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        elif aa < 4:
                            ax.set_ylim(self.phase_limits)
                            
                        else:
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                            if self.tipper_limits is not None:
                                ax.set_ylim(self.tipper_limits)
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

                    elif len(self.ax_list) == 6 and plot_tipper == True:

                        if aa < 2: # Changes applied
                            # plt.setp(ax.get_xticklabels(), visible=False)
                            if self.plot_z == False:
                                if aa == 0 or aa == 1:
                                    ax.set_yscale('log', nonposy='clip')
                                    ylim = ax.get_ylim()
                                    ylimits = (10 ** (np.floor(np.log10(ylim[0]))),
                                               10 ** (np.ceil(np.log10(ylim[1]))))
                                    ax.set_ylim(ylimits)

                            if self.res_limits is not None:
                                ax.set_ylim(self.res_limits)
                        elif aa < 4:
                            ax.set_ylim(self.phase_limits)
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                        else:
                            ax.set_xlabel('Period (s)', fontdict=fontdict)
                            if self.tipper_limits is not None:
                                ax.set_ylim(self.tipper_limits)
                        if aa == 0:
                            if self.plot_z == False:
                                ax.set_ylabel('App. Res . ($\mathbf{\Omega \cdot m}$)',
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

                        if aa == 0:
                            ax.yaxis.set_major_formatter(LogFormatterSciNotation())
                        if aa == 2: # Setting the decimal places
                            ax.yaxis.set_major_formatter(
                                FormatStrFormatter('%.0f'))
                            pass
                            if self.plot_z == True:
                                ax.set_yscale('log', nonposy='clip')
                                # else:
                                #     plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                        if aa == 4:
                            if plot_tipper:
                                ax.set_ylabel('Tipper', fontdict=fontdict)
                    # writing x axis ticks and making it visible

                    if ((aa == 4 or aa == 5) and plot_tipper == True) or \
                        ((aa == 2 or aa == 3) and plot_tipper == False):
                        plt.setp(ax.get_xticklabels(), visible=True)
                    else:
                        plt.setp(ax.get_xticklabels(), visible=False)

                    ax.set_xscale('log', nonposx='clip')
                    # set period limits
                    if self.period_limits is None:
                        self.period_limits = (10 ** (np.floor(np.log10(period[0]))) * 1.01,
                                              10 ** (np.ceil(np.log10(period[-1]))) * .99)
                    ax.set_xlim(xmin=self.period_limits[0],
                                xmax=self.period_limits[1])
                    ax.grid(True, alpha=.25)

            if plotr == True:
                for rr in range(nr):
                    if self.color_mode == 'color':
                        cxy = (0, .4 + float(rr) / (3 * nr), 0)
                        cyx = (.7 + float(rr) / (4 * nr), .13, .63 -
                               float(rr) / (4 * nr))
                    elif self.color_mode == 'bw':
                        cxy = tuple(3 * [1 - .5 / (rr + 1)])
                        cyx = tuple(3 * [1 - .5 / (rr + 1)])

                    resp_z_obj = self.resp_object[rr].mt_dict[station].Z
                    resp_z_err = np.nan_to_num(
                        (z_obj.z - resp_z_obj.z) / z_obj.z_err)

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
                    # print ' --- response {0} ---'.format(rr)
                    # print '  RMS = {:.2f}'.format(rms)
                    # print '      RMS_xx = {:.2f}'.format(rms_xx)
                    # print '      RMS_xy = {:.2f}'.format(rms_xy)
                    # print '      RMS_yx = {:.2f}'.format(rms_yx)
                    # print '      RMS_yy = {:.2f}'.format(rms_yy)
                    # print '      RMS_Tx = {:.2f}'.format(rms_tx)
                    # print '      RMS_Ty = {:.2f}'.format(rms_ty)

                    # --> make key word dictionaries for plotting
                    kw_xx = {'color': self.ctem,#cxy,
                             'marker': self.mtem,
                             'ms': self.ms,
                             'ls': self.ls,
                             'lw': self.lw,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

                    kw_yy = {'color': self.ctmm,#cyx,
                             'marker': self.mtmm,
                             'ms': self.ms,
                             'ls': self.ls,
                             'lw': self.lw,
                             'e_capsize': self.e_capsize,
                             'e_capthick': self.e_capthick}

#                    if self.plot_style == 1:
#                        if self.plot_component == 2:
#                            if self.plot_z == False:
#                                # plot resistivity
#                                rerxy = mtplottools.plot_errorbar(axrxy,
#                                                                  period[nzxy],
#                                                                  rrp.resxy[
#                                                                      nzxy],
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axryx,
#                                                                  period[nzyx],
#                                                                  rrp.resyx[
#                                                                      nzyx],
#                                                                  **kw_yy)
#                                # plot phase
#                                rerxy = mtplottools.plot_errorbar(axpxy,
#                                                                  period[nzxy],
#                                                                  rrp.phasexy[
#                                                                      nzxy],
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axpyx,
#                                                                  period[nzyx],
#                                                                  rrp.phaseyx[
#                                                                      nzyx],
#                                                                  **kw_yy)
#                            elif self.plot_z == True:
#                                # plot real
#                                rerxy = mtplottools.plot_errorbar(axrxy,
#                                                                  period[nzxy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxy, 0, 1].real),
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axryx,
#                                                                  period[nzyx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyx, 1, 0].real),
#                                                                  **kw_yy)
#                                # plot phase
#                                rerxy = mtplottools.plot_errorbar(axpxy,
#                                                                  period[nzxy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxy, 0, 1].imag),
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axpyx,
#                                                                  period[nzyx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyx, 1, 0].imag),
#                                                                  **kw_yy)
#                            if plot_tipper == True:
#                                rertx = mtplottools.plot_errorbar(axtr,
#                                                                  period[ntx],
#                                                                  resp_t_obj.tipper[
#                                                                      ntx, 0, 0].real,
#                                                                  **kw_xx)
#                                rerty = mtplottools.plot_errorbar(axtr,
#                                                                  period[nty],
#                                                                  resp_t_obj.tipper[
#                                                                      nty, 0, 1].real,
#                                                                  **kw_yy)
#
#                                rertx = mtplottools.plot_errorbar(axti,
#                                                                  period[ntx],
#                                                                  resp_t_obj.tipper[
#                                                                      ntx, 0, 0].imag,
#                                                                  **kw_xx)
#                                rerty = mtplottools.plot_errorbar(axti,
#                                                                  period[nty],
#                                                                  resp_t_obj.tipper[
#                                                                      nty, 0, 1].imag,
#                                                                  **kw_yy)
#
#                            if plot_tipper == False:
#                                line_list[0] += [rerxy[0]]
#                                line_list[1] += [reryx[0]]
#                                label_list[0] += ['$Z^m_{xy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xy)]
#                                label_list[1] += ['$Z^m_{yx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yx)]
#                            else:
#                                line_list[0] += [rerxy[0]]
#                                line_list[1] += [reryx[0]]
#                                line_list[2] += [rertx[0], rerty[0]]
#                                label_list[0] += ['$Z^m_{xy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xy)]
#                                label_list[1] += ['$Z^m_{yx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yx)]
#                                label_list[2] += ['$T^m_{x}$' +
#                                                  'rms={0:.2f}'.format(rms_tx),
#                                                  '$T^m_{y}$' +
#                                                  'rms={0:.2f}'.format(rms_ty)]
#                        elif self.plot_component == 4:
#                            if self.plot_z == False:
#                                # plot resistivity
#                                rerxx = mtplottools.plot_errorbar(axrxx,
#                                                                  period[nzxx],
#                                                                  rrp.resxx[
#                                                                      nzxx],
#                                                                  **kw_xx)
#                                rerxy = mtplottools.plot_errorbar(axrxy,
#                                                                  period[nzxy],
#                                                                  rrp.resxy[
#                                                                      nzxy],
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axryx,
#                                                                  period[nzyx],
#                                                                  rrp.resyx[
#                                                                      nzyx],
#                                                                  **kw_yy)
#                                reryy = mtplottools.plot_errorbar(axryy,
#                                                                  period[nzyy],
#                                                                  rrp.resyy[
#                                                                      nzyy],
#                                                                  **kw_yy)
#                                # plot phase
#                                rerxx = mtplottools.plot_errorbar(axpxx,
#                                                                  period[nzxx],
#                                                                  rrp.phasexx[
#                                                                      nzxx],
#                                                                  **kw_xx)
#                                rerxy = mtplottools.plot_errorbar(axpxy,
#                                                                  period[nzxy],
#                                                                  rrp.phasexy[
#                                                                      nzxy],
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axpyx,
#                                                                  period[nzyx],
#                                                                  rrp.phaseyx[
#                                                                      nzyx],
#                                                                  **kw_yy)
#                                reryy = mtplottools.plot_errorbar(axpyy,
#                                                                  period[nzyy],
#                                                                  rrp.phaseyy[
#                                                                      nzyy],
#                                                                  **kw_yy)
#                            elif self.plot_z == True:
#                                # plot real
#                                rerxx = mtplottools.plot_errorbar(axrxx,
#                                                                  period[nzxx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxx, 0, 0].real),
#                                                                  **kw_xx)
#                                rerxy = mtplottools.plot_errorbar(axrxy,
#                                                                  period[nzxy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxy, 0, 1].real),
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axryx,
#                                                                  period[nzyx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyx, 1, 0].real),
#                                                                  **kw_yy)
#                                reryy = mtplottools.plot_errorbar(axryy,
#                                                                  period[nzyy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyy, 1, 1].real),
#                                                                  **kw_yy)
#                                # plot phase
#                                rerxx = mtplottools.plot_errorbar(axpxx,
#                                                                  period[nzxx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxx, 0, 0].imag),
#                                                                  **kw_xx)
#                                rerxy = mtplottools.plot_errorbar(axpxy,
#                                                                  period[nzxy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzxy, 0, 1].imag),
#                                                                  **kw_xx)
#                                reryx = mtplottools.plot_errorbar(axpyx,
#                                                                  period[nzyx],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyx, 1, 0].imag),
#                                                                  **kw_yy)
#                                reryy = mtplottools.plot_errorbar(axpyy,
#                                                                  period[nzyy],
#                                                                  abs(resp_z_obj.z[
#                                                                      nzyy, 1, 1].imag),
#                                                                  **kw_yy)
#                            if plot_tipper == True:
#                                rertx = mtplottools.plot_errorbar(axtxr,
#                                                                  period[ntx],
#                                                                  resp_t_obj.tipper[
#                                                                      ntx, 0, 0].real,
#                                                                  **kw_xx)
#                                rerty = mtplottools.plot_errorbar(axtyr,
#                                                                  period[nty],
#                                                                  resp_t_obj.tipper[
#                                                                      nty, 0, 1].real,
#                                                                  **kw_yy)
#
#                                rertx = mtplottools.plot_errorbar(axtxi,
#                                                                  period[ntx],
#                                                                  resp_t_obj.tipper[
#                                                                      ntx, 0, 0].imag,
#                                                                  **kw_xx)
#                                rerty = mtplottools.plot_errorbar(axtyi,
#                                                                  period[nty],
#                                                                  resp_t_obj.tipper[
#                                                                      nty, 0, 1].imag,
#                                                                  **kw_yy)
#
#                            if plot_tipper == False:
#                                line_list[0] += [rerxx[0]]
#                                line_list[1] += [rerxy[0]]
#                                line_list[2] += [reryx[0]]
#                                line_list[3] += [reryy[0]]
#                                label_list[0] += ['$Z^m_{xx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xx)]
#                                label_list[1] += ['$Z^m_{xy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xy)]
#                                label_list[2] += ['$Z^m_{yx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yx)]
#                                label_list[3] += ['$Z^m_{yy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yy)]
#                            else:
#                                line_list[0] += [rerxx[0]]
#                                line_list[1] += [rerxy[0]]
#                                line_list[2] += [reryx[0]]
#                                line_list[3] += [reryy[0]]
#                                line_list[4] += [rertx[0]]
#                                line_list[5] += [rerty[0]]
#                                label_list[0] += ['$Z^m_{xx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xx)]
#                                label_list[1] += ['$Z^m_{xy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_xy)]
#                                label_list[2] += ['$Z^m_{yx}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yx)]
#                                label_list[3] += ['$Z^m_{yy}$ ' +
#                                                  'rms={0:.2f}'.format(rms_yy)]
#                                label_list[4] += ['$T^m_{x}$' +
#                                                  'rms={0:.2f}'.format(rms_tx)]
#                                label_list[5] += ['$T^m_{y}$' +
#                                                  'rms={0:.2f}'.format(rms_ty)]

                    if self.plot_style == 2:
#                        rrp.phasexy[rrp.phasexy > 180] -= 360
#                        rrp.phaseyx[rrp.phaseyx > 180] -= 360

                        if self.plot_component == 2:
                            if self.plot_z == False:
                                # plot resistivity
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  rrp.resxy[
                                                                      nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  rrp.resyx[
                                                                      nzyx],
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  rrp.phasexy[
                                                                      nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  rrp.phaseyx[
                                                                      nzyx],
                                                                  **kw_yy)
                            elif self.plot_z == True:
                                # plot real
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[
                                                                      nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[
                                                                      nzyx, 1, 0].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[
                                                                      nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[
                                                                      nzyx, 1, 0].imag),
                                                                  **kw_xx)
                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[
                                                                      ntx, 0, 0].real,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[
                                                                      nty, 0, 1].real,
                                                                  **kw_yy)

                                rertx = mtplottools.plot_errorbar(axti,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[
                                                                      ntx, 0, 0].imag,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[
                                                                      nty, 0, 1].imag,
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
                                                                  rrp.resxx[
                                                                      nzxx],
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  rrp.resxy[
                                                                      nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  rrp.resyx[
                                                                      nzyx],
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzyy],
                                                                  rrp.resyy[
                                                                      nzyy],
                                                                  **kw_yy)
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzxx],
                                                                  rrp.phasexx[
                                                                      nzxx],
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  rrp.phasexy[
                                                                      nzxy],
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  rrp.phaseyx[
                                                                      nzyx],
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzyy],
                                                                  rrp.phaseyy[
                                                                      nzyy],
                                                                  **kw_yy)
                            elif self.plot_z == True:
                                # plot real
                                rerxx = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[
                                                                      nzxx, 0, 0].real),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[
                                                                      nzxy, 0, 1].real),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axrxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[
                                                                      nzyx, 1, 0].real),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axrxx,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[
                                                                      nzyy, 1, 1].real),
                                                                  **kw_yy)
                                # plot phase
                                rerxx = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzxx],
                                                                  abs(resp_z_obj.z[
                                                                      nzxx, 0, 0].imag),
                                                                  **kw_xx)
                                rerxy = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzxy],
                                                                  abs(resp_z_obj.z[
                                                                      nzxy, 0, 1].imag),
                                                                  **kw_xx)
                                reryx = mtplottools.plot_errorbar(axpxy,
                                                                  period[nzyx],
                                                                  abs(resp_z_obj.z[
                                                                      nzyx, 1, 0].imag),
                                                                  **kw_yy)
                                reryy = mtplottools.plot_errorbar(axpxx,
                                                                  period[nzyy],
                                                                  abs(resp_z_obj.z[
                                                                      nzyy, 1, 1].imag),
                                                                  **kw_yy)

                            if plot_tipper == True:
                                rertx = mtplottools.plot_errorbar(axtr,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[
                                                                      ntx, 0, 0].real,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axtr,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[
                                                                      nty, 0, 1].real,
                                                                  **kw_yy)

                                rertx = mtplottools.plot_errorbar(axti,
                                                                  period[ntx],
                                                                  resp_t_obj.tipper[
                                                                      ntx, 0, 0].imag,
                                                                  **kw_xx)
                                rerty = mtplottools.plot_errorbar(axti,
                                                                  period[nty],
                                                                  resp_t_obj.tipper[
                                                                      nty, 0, 1].imag,
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
#            if self.plot_style == 1:
#                legend_ax_list = self.ax_list[0:self.plot_component]
#                if plot_tipper == True:
#                    if self.plot_component == 2:
#                        legend_ax_list.append(self.ax_list[4])
#                    elif self.plot_component == 4:
#                        legend_ax_list.append(self.ax_list[8])
#                        legend_ax_list.append(self.ax_list[10])
#                for aa, ax in enumerate(legend_ax_list):
#                    ax.legend(line_list[aa],
#                              label_list[aa],
#                              loc=self.legend_loc,
#                              bbox_to_anchor=self.legend_pos,
#                              markerscale=self.legend_marker_scale,
#                              borderaxespad=self.legend_border_axes_pad,
#                              labelspacing=self.legend_label_spacing,
#                              handletextpad=self.legend_handle_text_pad,
#                              borderpad=self.legend_border_pad,
#                              prop={'size': max([self.font_size / (nr + 1), 5])})
            if self.plot_style == 2:
                if self.plot_component == 2:
                    legend_ax_list = [self.ax_list[0]]
                    if plot_tipper == True:
                        legend_ax_list.append(self.ax_list[2])

                    for aa, ax in enumerate(legend_ax_list):
                        # work out if tipper or z legend and update position accordingly
                        if aa < 1:
                            legend_pos = self.legend_pos
                        else:
                            legend_pos = self.legend_pos_tipper
                        ax.legend(line_list[aa],
                                  label_list[aa],
                                  loc=self.legend_loc,
                                  bbox_to_anchor=legend_pos,
                                  markerscale=self.legend_marker_scale,
                                  borderaxespad=self.legend_border_axes_pad,
                                  labelspacing=self.legend_label_spacing,
                                  handletextpad=self.legend_handle_text_pad,
                                  borderpad=self.legend_border_pad,
                                  prop={'size': max([self.font_size / (nr + 1), 4])})

                        
                else:
                    legend_ax_list = self.ax_list[0:int(self.plot_component / 2)]
                    if plot_tipper == True:
                        if self.plot_component == 2:
                            legend_ax_list.append(self.ax_list[2])
                        elif self.plot_component == 4:
                            legend_ax_list.append(self.ax_list[4])
                            
#                    # add text to distinguish real and imaginary tipper
#                    for aa, ax in enumerate(self.ax_list[4:]):
#                        ax.text(0.5,0.8,['Real','Imaginary'][aa],
#                                ha='center',va='center',
#                                transform=ax.transAxes)
                            

                    for aa, ax in enumerate(legend_ax_list):
                        if aa < 2:
                            legend_pos = self.legend_pos
                        else:
                            legend_pos = self.legend_pos_tipper

                        ax.legend(line_list[aa],
                                  label_list[aa],
                                  loc=self.legend_loc,
                                  bbox_to_anchor=legend_pos,
                                  markerscale=self.legend_marker_scale,
                                  borderaxespad=self.legend_border_axes_pad,
                                  labelspacing=self.legend_label_spacing,
                                  handletextpad=self.legend_handle_text_pad,
                                  borderpad=self.legend_border_pad,
                                  prop={'size': max([self.font_size / (nr + 1), 4])})
                        
            if self.save_plots:
                save_filename = os.path.join(os.path.dirname(self.data_fn),station+'.png')
                self.save_figure(save_filename,fig_dpi=self.fig_dpi)
        else:
            pass

        plt.show()   # --> BE SURE TO SHOW THE PLOT


    def _get_station_index_list(self):
        """
        get list of station indices to plot
        
        """
        
        allstations = self.data_object.station_locations.station
        allstations_up = [str.upper(val) for val in allstations]
        if self.plot_type == '1':
            # plot all stations
            s_index_list = range(len(allstations))
        else:
            s_index_list = []
            for st in self.plot_type:
                # if integer, assume it is an index
                if type(st) == int:
                    s_index_list.append(st)
                # else, search for station name (case-sensitive)
                elif st in allstations:
                    s_index_list.append(list(allstations).index(st))
                # else, search for station name (not case-sensitive) 
                elif str.upper(st) in allstations_up:
                    s_index_list.append(list(allstations_up).index(str.upper(st)))
                    
        return s_index_list
                    

    def _plot_1col(self):
        """
        plot data and response in a 1-column plot
        All elements of the impedance tensor are overlaid, with diagonals semi-
        transparent (XX has same color as ctem/cted, YY has same color as ctmm/ctmd)
        
        """
        # read files
        self._read_files()
        
        # get station indices to plot
        s_index_list = self._get_station_index_list()
        
        # collate a list of data objects to plot
        data_objects = [self.data_object]
        if type(self.resp_object) in [list,np.ndarray]:
            data_objects += self.resp_object
        else:
            data_objects.append(self.resp_object)
            
        # get residual object if data and response are provided
        if ((self.data_fn is not None) and (self.resp_fn is not None)):
            rsObj = Residual()
            rsObj.calculate_residual_from_data(data_fn=self.data_fn,
                                               resp_fn=self.resp_fn,
                                               save=False)
            rms = rsObj.rms
        else:
            rms = np.nan

        # get period limits
        if self.period_limits is None:
            self.period_limits = (10 ** (np.floor(np.log10(self.data_object.period_list[0]))) * 1.01,
                                  10 ** (np.ceil(np.log10(self.data_object.period_list[-1]))) * .99)
        
        # initialise color/marker/linestyle/transparency lists for plotting
        color_z = np.array([[[self.cted,self.cted],[self.ctmd,self.ctmd]],
                            [[self.ctem,self.ctem],[self.ctmm,self.ctmm]]])
        markers = np.array([[[self.mted,self.mted],[self.mtmd,self.mtmd]],
                            [[self.mtem,self.mtem],[self.mtmm,self.mtmm]]])
        linestyles=['','--']
        alpha = 1.-np.eye(2)*0.8
        color_tip = np.array([[self.cted,self.ctmd],[self.ctem,self.ctmd]])
        
        
        #--> set default font size
        plt.rcParams['font.size'] = self.font_size

        fontdict = {'size': self.font_size + 2, 'weight': 'bold'}
        
        for ii,si in enumerate(s_index_list):
            fig = plt.figure(figsize=self.fig_size)
        
            # make some axes for the plot
            if self.label_axes:
                left = 0.25
            else:
                left = 0.15
            
            
            if self.plot_z:
                gs = gridspec.GridSpec(4,1,height_ratios=[0.85,0.85,0.3,0.3],hspace=0.1,left=left)
            else:
                gs = gridspec.GridSpec(4,1,height_ratios=[1,0.7,0.3,0.3],hspace=0.1,left=left)
            axr = fig.add_subplot(gs[0,0],xscale='log',yscale='log',xlim=self.period_limits)
            axp = fig.add_subplot(gs[1,0],xscale='log',sharex=axr)
            axt = [fig.add_subplot(gs[2,0],xscale='log',sharex=axr),
                   fig.add_subplot(gs[3,0],xscale='log',sharex=axr)]
            
            
            for di in range(2):
                dObj = data_objects[di]
                
                zObj = Z(z_array=dObj.data_array['z'][si],
                         z_err_array=dObj.data_array['z_err'][si],
                         freq=1./dObj.period_list)
                
                tObj = Tipper(tipper_array=dObj.data_array['tip'][si],
                              tipper_err_array=dObj.data_array['tip_err'][si],
                              freq=1./dObj.period_list)
                
                if self.plot_z:
                    data1 = np.abs(zObj.z.real)
                    data2 = np.abs(zObj.z.imag)
                    data1label = 'Z (real)'
                    data2label = 'Z (imag)'
                    data1err = zObj.z_err
                    data2err = zObj.z_err
                    axp.set_yscale('log')
                else:
                    data1 = zObj.resistivity
                    data2 = zObj.phase
                    if self.shift_yx_phase:
                        data2[:,1,0] += 180.
                    data1label = 'Resistivity, $\Omega$m'
                    data2label = 'Phase, degree'
                    data1err = zObj.resistivity_err
                    data2err = zObj.phase_err
                    
                if di >= 1:
                    data1err = np.zeros_like(data1err)
                    data2err = np.zeros_like(data2err)
                
                for i in range(2):
                    for j in range(2):
                        kwargs={'color':color_z[di,i,j],
                                'ls':linestyles[di],
                                'marker':markers[di,i,j],
                                'alpha':alpha[i,j]}
                        nonzero = np.nonzero(data1[:,i,j])[0]
                        axr.errorbar(1./zObj.freq[nonzero],data1[nonzero][:,i,j],yerr=data1err[nonzero][:,i,j],label='Res'+'XY'[i]+'XY'[j],**kwargs)
                        axp.errorbar(1./zObj.freq[nonzero],data2[nonzero][:,i,j],yerr=data2err[nonzero][:,i,j],label='Phs'+'XY'[i]+'XY'[j],**kwargs)
                     
                    nonzerot = np.nonzero(tObj.tipper[:,0,i])
                    tipper_err = tObj.tipper_err[nonzerot][:,0,i]
                    if di >= 1:
                        tipper_err = np.zeros_like(tipper_err)
                    
                    kwargs['alpha'] = 1.
                    kwargs['color'] = color_tip[di,0]
                    axt[i].errorbar(1./tObj.freq[nonzerot],tObj.tipper.real[nonzerot][:,0,i],tipper_err,label='Tip'+'XY'[i]+'R',**kwargs)
                    kwargs['color'] = color_tip[di,1]
                    axt[i].errorbar(1./tObj.freq[nonzerot],tObj.tipper.imag[nonzerot][:,0,i],tipper_err,label='Tip'+'XY'[i]+'I',**kwargs)
                    axt[i].set_ylim(self.tipper_limits)

                    
            sname = self.data_object.station_locations.station[si]
            
            axp.set_ylim(self.phase_limits)
            axr.set_ylim(self.res_limits)
            axr.set_title('%s, RMS=%.1f'%(sname,rms),fontdict=fontdict)
            
            axt[1].set_xlabel('Period, s',fontsize=self.font_size)
            axr.set_xlim(self.period_limits)
        
            for ax in [axr,axp,axt[0]]:
                plt.setp(ax.xaxis.get_ticklabels(),visible=False)
            
            if self.tipper_limits is not None:
                yticks = [np.round(self.tipper_limits[0]*0.75,1),0,np.round(self.tipper_limits[1]*0.75,1)]
                for axti in axt:
                    axti.set_yticks(yticks)
                    
            if self.label_axes:
                axr.set_ylabel(data1label,labelpad=0)
                axp.set_ylabel(data2label,labelpad=0)
                
                axt[0].set_ylabel('Tipper, X',labelpad=0)
                axt[1].set_ylabel('Tipper, Y',labelpad=0)
        
        

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

    def save_figure(self, save_fn, file_format='pdf', orientation='portrait',
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
        print('Saved figure to: ' + self.fig_fn)
