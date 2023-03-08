# -*- coding: utf-8 -*-
"""
Spin-off from 'occamtools'
(Created August 2011, re-written August 2013)

Tools for Occam2D

authors: JP/LK


Classes:
    - Data
    - Model
    - Setup
    - Run
    - Plot
    - Mask


Functions:
    - getdatetime
    - makestartfiles
    - writemeshfile
    - writemodelfile
    - writestartupfile
    - read_datafile
    - get_model_setup
    - blocks_elements_setup


"""
# ==============================================================================
import os
import os.path as op
import sys
import time

import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.interpolate as spi
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator
from scipy.stats import mode

import mtpy.analysis.geometry as MTgy
import mtpy.core.mt as mt
import mtpy.modeling.winglinktools as MTwl
from mtpy.imaging.mtplottools import plot_errorbar
from mtpy.utils.calculator import centre_point
from mtpy.utils.gis_tools import get_epsg, project_point_ll2utm


# ==============================================================================


# ==============================================================================
# plot the MT and model responses
# ==============================================================================
class PlotResponse:
    """
    Helper class to deal with plotting the MT response and occam2d model.

    Arguments:
    -------------
        **data_fn** : string
                      full path to data file

        **resp_fn** : string or list
                      full path(s) to response file(s)


    ==================== ======================================================
    Attributes/key words            description
    ==================== ======================================================
    ax_list              list of matplotlib.axes instances for use with
                         OccamPointPicker
    color_mode           [ 'color' | 'bw' ] plot figures in color or
                         black and white ('bw')
    cted                 color of Data TE marker and line
    ctem                 color of Model TE marker and line
    ctewl                color of Winglink Model TE marker and line
    ctmd                 color of Data TM marker and line
    ctmm                 color of Model TM marker and line
    ctmwl                color of Winglink Model TM marker and line
    e_capsize            size of error bar caps in points
    e_capthick           line thickness of error bar caps in points
    err_list             list of line properties of error bars for use with
                         OccamPointPicker
    fig_dpi              figure resolution in dots-per-inch
    fig_list             list of dictionaries with key words
                         station --> station name
                         fig --> matplotlib.figure instance
                         axrte --> matplotlib.axes instance for TE app.res
                         axrtm --> matplotlib.axes instance for TM app.res
                         axpte --> matplotlib.axes instance for TE phase
                         axptm --> matplotlib.axes instance for TM phase

    fig_num              starting number of figure
    fig_size             size of figure in inches (width, height)
    font_size            size of axes ticklabel font in points
    line_list            list of matplotlib.Line instances for use with
                         OccamPointPicker
    lw                   line width of lines in points
    ms                   marker size in points
    mted                 marker for Data TE mode
    mtem                 marker for Model TE mode
    mtewl                marker for Winglink Model TE
    mtmd                 marker for Data TM mode
    mtmm                 marker for Model TM mode
    mtmwl                marker for Winglink TM mode
    period               np.ndarray of periods to plot
    phase_limits         limits on phase plots in degrees (min, max)
    plot_num             [ 1 | 2 ]
                         1 to plot both modes in a single plot
                         2 to plot modes in separate plots (default)
    plot_tipper          [ 'y' | 'n' ] plot tipper data if desired
    plot_type            [ '1' | station_list]
                         '1' --> to plot all stations in different figures
                         station_list --> to plot a few stations, give names
                         of stations ex. ['mt01', 'mt07']
    plot_yn              [ 'y' | 'n']
                         'y' --> to plot on instantiation
                         'n' --> to not plot on instantiation
    res_limits           limits on resistivity plot in log scale (min, max)
    rp_list               list of dictionaries from read2Ddata
    station_list          station_list list of stations in rp_list
    subplot_bottom       subplot spacing from bottom (relative coordinates)
    subplot_hspace       vertical spacing between subplots
    subplot_left         subplot spacing from left
    subplot_right        subplot spacing from right
    subplot_top          subplot spacing from top
    subplot_wspace       horizontal spacing between subplots
    wl_fn                Winglink file name (full path)
    ==================== ======================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots the apparent resistiviy and phase of data and
                        model if given.  called on instantiation if plot_yn
                        is 'y'.
    redraw_plot         call redraw_plot to redraw the figures,
                        if one of the attributes has been changed
    save_figures        save all the matplotlib.figure instances in fig_list
    =================== =======================================================


    :Example: ::
        >>> data_fn = r"/home/occam/line1/inv1/OccamDataFile.dat"
        >>> resp_list = [r"/home/occam/line1/inv1/test_{0:02}".format(ii)
                         for ii in range(2, 8, 2)]
        >>> pr_obj = occam2d.PlotResponse(data_fn, resp_list, plot_tipper='y')

    """

    def __init__(self, data_fn, resp_fn=None, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn
        if self.resp_fn is not None:
            if type(self.resp_fn) != list:
                self.resp_fn = [self.resp_fn]

        self.wl_fn = kwargs.pop("wl_fn", None)

        self.color_mode = kwargs.pop("color_mode", "color")

        self.ms = kwargs.pop("ms", 1.5)
        self.lw = kwargs.pop("lw", 0.5)
        self.e_capthick = kwargs.pop("e_capthick", 0.5)
        self.e_capsize = kwargs.pop("e_capsize", 2)

        self.ax_list = []
        self.line_list = []
        self.err_list = []

        # color mode
        if self.color_mode == "color":
            # color for data
            self.cted = kwargs.pop("cted", (0, 0, 1))
            self.ctmd = kwargs.pop("ctmd", (1, 0, 0))
            self.mted = kwargs.pop("mted", "s")
            self.mtmd = kwargs.pop("mtmd", "o")

            # color for occam2d model
            self.ctem = kwargs.pop("ctem", (0, 0.6, 0.3))
            self.ctmm = kwargs.pop("ctmm", (0.9, 0, 0.8))
            self.mtem = kwargs.pop("mtem", "+")
            self.mtmm = kwargs.pop("mtmm", "+")

            # color for Winglink model
            self.ctewl = kwargs.pop("ctewl", (0, 0.6, 0.8))
            self.ctmwl = kwargs.pop("ctmwl", (0.8, 0.7, 0))
            self.mtewl = kwargs.pop("mtewl", "x")
            self.mtmwl = kwargs.pop("mtmwl", "x")

            # color of tipper
            self.ctipr = kwargs.pop("ctipr", self.cted)
            self.ctipi = kwargs.pop("ctipi", self.ctmd)

        # black and white mode
        elif self.color_mode == "bw":
            # color for data
            self.cted = kwargs.pop("cted", (0, 0, 0))
            self.ctmd = kwargs.pop("ctmd", (0, 0, 0))
            self.mted = kwargs.pop("mted", "*")
            self.mtmd = kwargs.pop("mtmd", "v")

            # color for occam2d model
            self.ctem = kwargs.pop("ctem", (0.6, 0.6, 0.6))
            self.ctmm = kwargs.pop("ctmm", (0.6, 0.6, 0.6))
            self.mtem = kwargs.pop("mtem", "+")
            self.mtmm = kwargs.pop("mtmm", "x")

            # color for Winglink model
            self.ctewl = kwargs.pop("ctewl", (0.3, 0.3, 0.3))
            self.ctmwl = kwargs.pop("ctmwl", (0.3, 0.3, 0.3))
            self.mtewl = kwargs.pop("mtewl", "|")
            self.mtmwl = kwargs.pop("mtmwl", "_")

            self.ctipr = kwargs.pop("ctipr", self.cted)
            self.ctipi = kwargs.pop("ctipi", self.ctmd)

        self.phase_limits = kwargs.pop("phase_limits", (-5, 95))
        self.res_limits = kwargs.pop("res_limits", None)
        self.tip_limits = kwargs.pop("tip_limits", (-0.5, 0.5))

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)

        self.subplot_wspace = 0.1
        self.subplot_hspace = 0.15
        self.subplot_right = 0.98
        self.subplot_left = 0.085
        self.subplot_top = 0.93
        self.subplot_bottom = 0.1

        self.font_size = kwargs.pop("font_size", 6)

        self.plot_type = kwargs.pop("plot_type", "1")
        self.plot_num = kwargs.pop("plot_num", 2)
        self.plot_tipper = kwargs.pop("plot_tipper", "n")
        self.plot_model_error = kwargs.pop("plot_model_err", "y")
        self.plot_yn = kwargs.pop("plot_yn", "y")

        if self.plot_num == 1:
            self.ylabel_coord = kwargs.pop("ylabel_coords", (-0.055, 0.5))
        elif self.plot_num == 2:
            self.ylabel_coord = kwargs.pop("ylabel_coords", (-0.12, 0.5))

        self.fig_list = []

        if self.plot_yn == "y":
            self.plot()

    def plot(self):
        """
        plot the data and model response, if given, in individual plots.

        """

        data_obj = Data()
        data_obj.read_data_file(self.data_fn)

        rp_list = data_obj.data
        nr = len(rp_list)

        # create station list
        self.station_list = [rp["station"] for rp in rp_list]

        # boolean for adding winglink output to the plots 0 for no, 1 for yes
        addwl = 0
        # read in winglink data file
        if self.wl_fn != None:
            addwl = 1
            self.subplot_hspace + 0.1
            wld, wlrp_list, wlplist, wlslist, wltlist = MTwl.readOutputFile(
                self.wl_fn
            )
            sdict = dict(
                [
                    (ostation, wlistation)
                    for wlistation in wlslist
                    for ostation in self.station_list
                    if wlistation.find(ostation) >= 0
                ]
            )

        # set a local parameter period for less typing
        period = data_obj.period

        # ---------------plot each respones in a different figure--------------
        if self.plot_type == "1":
            pstation_list = list(range(len(self.station_list)))

        else:
            if type(self.plot_type) is not list:
                self.plot_type = [self.plot_type]

            pstation_list = []
            for ii, station in enumerate(self.station_list):
                for pstation in self.plot_type:
                    if station.find(pstation) >= 0:
                        pstation_list.append(ii)

                        # set the grid of subplots
        if self.plot_tipper == "y":
            gs = gridspec.GridSpec(
                3,
                2,
                wspace=self.subplot_wspace,
                left=self.subplot_left,
                top=self.subplot_top,
                bottom=self.subplot_bottom,
                right=self.subplot_right,
                hspace=self.subplot_hspace,
                height_ratios=[2, 1.5, 1],
            )
        else:
            gs = gridspec.GridSpec(
                2,
                2,
                wspace=self.subplot_wspace,
                left=self.subplot_left,
                top=self.subplot_top,
                bottom=self.subplot_bottom,
                right=self.subplot_right,
                hspace=self.subplot_hspace,
                height_ratios=[2, 1.5],
            )

        # --> set default font size
        plt.rcParams["font.size"] = self.font_size

        # loop over each station to plot
        for ii, jj in enumerate(pstation_list):
            fig = plt.figure(
                self.station_list[jj], self.fig_size, dpi=self.fig_dpi
            )
            plt.clf()

            # --> set subplot instances
            # ---plot both TE and TM in same subplot---
            if self.plot_num == 1:
                axrte = fig.add_subplot(gs[0, :])
                axrtm = axrte
                axpte = fig.add_subplot(gs[1, :], sharex=axrte)
                axptm = axpte
                if self.plot_tipper == "y":
                    axtipre = fig.add_subplot(gs[2, :], sharex=axrte)
                    axtipim = axtipre

            # ---plot TE and TM in separate subplots---
            elif self.plot_num == 2:
                axrte = fig.add_subplot(gs[0, 0])
                axrtm = fig.add_subplot(gs[0, 1])
                axpte = fig.add_subplot(gs[1, 0], sharex=axrte)
                axptm = fig.add_subplot(gs[1, 1], sharex=axrtm)
                if self.plot_tipper == "y":
                    axtipre = fig.add_subplot(gs[2, 0], sharex=axrte)
                    axtipim = fig.add_subplot(gs[2, 1], sharex=axrtm)

            # plot the data, it should be the same for all response files
            # empty lists for legend marker and label
            rlistte = []
            llistte = []
            rlisttm = []
            llisttm = []
            # ------------Plot Resistivity----------------------------------
            # cut out missing data points first
            # --> data
            rxy = np.where(rp_list[jj]["te_res"][0] != 0)[0]
            ryx = np.where(rp_list[jj]["tm_res"][0] != 0)[0]

            # --> TE mode Data
            if len(rxy) > 0:
                rte_err = (
                    rp_list[jj]["te_res"][1, rxy]
                    * rp_list[jj]["te_res"][0, rxy]
                )
                rte = plot_errorbar(
                    axrte,
                    period[rxy],
                    rp_list[jj]["te_res"][0, rxy],
                    ls=":",
                    marker=self.mted,
                    ms=self.ms,
                    color=self.cted,
                    y_error=rte_err,
                    lw=self.lw,
                    e_capsize=self.e_capsize,
                    e_capthick=self.e_capthick,
                )

                rlistte.append(rte[0])
                llistte.append("$Obs_{TE}$")
            else:
                rte = [None, [None, None, None], [None, None, None]]

                # --> TM mode data
            if len(ryx) > 0:
                rtm_err = (
                    rp_list[jj]["tm_res"][1, ryx]
                    * rp_list[jj]["tm_res"][0, ryx]
                )
                rtm = plot_errorbar(
                    axrtm,
                    period[ryx],
                    rp_list[jj]["tm_res"][0, ryx],
                    ls=":",
                    marker=self.mtmd,
                    ms=self.ms,
                    color=self.ctmd,
                    y_error=rtm_err,
                    lw=self.lw,
                    e_capsize=self.e_capsize,
                    e_capthick=self.e_capthick,
                )

                rlisttm.append(rtm[0])
                llisttm.append("$Obs_{TM}$")
            else:
                rtm = [None, [None, None, None], [None, None, None]]
            # --------------------plot phase--------------------------------
            # cut out missing data points first
            # --> data
            pxy = np.where(rp_list[jj]["te_phase"][0] != 0)[0]
            pyx = np.where(rp_list[jj]["tm_phase"][0] != 0)[0]

            # --> TE mode data
            if len(pxy) > 0:
                pte = plot_errorbar(
                    axpte,
                    period[pxy],
                    rp_list[jj]["te_phase"][0, pxy],
                    ls=":",
                    marker=self.mted,
                    ms=self.ms,
                    color=self.cted,
                    y_error=rp_list[jj]["te_phase"][1, pxy],
                    lw=self.lw,
                    e_capsize=self.e_capsize,
                    e_capthick=self.e_capthick,
                )

            else:
                pte = [None, [None, None, None], [None, None, None]]

            # --> TM mode data
            if len(pyx) > 0:
                ptm = plot_errorbar(
                    axptm,
                    period[pyx],
                    rp_list[jj]["tm_phase"][0, pyx],
                    ls=":",
                    marker=self.mtmd,
                    ms=self.ms,
                    color=self.ctmd,
                    y_error=rp_list[jj]["tm_phase"][1, pyx],
                    lw=self.lw,
                    e_capsize=self.e_capsize,
                    e_capthick=self.e_capthick,
                )

            else:
                ptm = [None, [None, None, None], [None, None, None]]

            # append axis properties to lists that can be used by
            # OccamPointPicker
            self.ax_list.append([axrte, axrtm, axpte, axptm])
            self.line_list.append([rte[0], rtm[0], pte[0], ptm[0]])
            self.err_list.append(
                [
                    [rte[1][0], rte[1][1], rte[2][0]],
                    [rtm[1][0], rtm[1][1], rtm[2][0]],
                    [pte[1][0], pte[1][1], pte[2][0]],
                    [ptm[1][0], ptm[1][1], ptm[2][0]],
                ]
            )

            # ---------------------plot tipper---------------------------------
            if self.plot_tipper == "y":
                t_list = []
                t_label = []
                txy = np.where(rp_list[jj]["re_tip"][0] != 0)[0]
                tyx = np.where(rp_list[jj]["im_tip"][0] != 0)[0]
                # --> real tipper  data
                if len(txy) > 0:
                    per_list_p = []
                    tpr_list_p = []
                    per_list_n = []
                    tpr_list_n = []
                    for per, tpr in zip(
                        period[txy], rp_list[jj]["re_tip"][0, txy]
                    ):
                        if tpr >= 0:
                            per_list_p.append(per)
                            tpr_list_p.append(tpr)
                        else:
                            per_list_n.append(per)
                            tpr_list_n.append(tpr)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipre.stem(
                            per_list_p, tpr_list_p, markerfmt="^", basefmt="k"
                        )
                        plt.setp(m_line, "markerfacecolor", self.ctipr)
                        plt.setp(m_line, "markeredgecolor", self.ctipr)
                        plt.setp(m_line, "markersize", self.ms)
                        plt.setp(s_line, "linewidth", self.lw)
                        plt.setp(s_line, "color", self.ctipr)
                        plt.setp(b_line, "linewidth", 0.01)
                        t_list.append(m_line)
                        t_label.append("Real")
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipre.stem(
                            per_list_n, tpr_list_n, markerfmt="v", basefmt="k"
                        )
                        plt.setp(m_line, "markerfacecolor", self.ctipr)
                        plt.setp(m_line, "markeredgecolor", self.ctipr)
                        plt.setp(m_line, "markersize", self.ms)
                        plt.setp(s_line, "linewidth", self.lw)
                        plt.setp(s_line, "color", self.ctipr)
                        plt.setp(b_line, "linewidth", 0.01)
                        if len(t_list) == 0:
                            t_list.append(m_line)
                            t_label.append("Real")

                else:
                    pass
                if len(tyx) > 0:
                    per_list_p = []
                    tpi_list_p = []
                    per_list_n = []
                    tpi_list_n = []
                    for per, tpi in zip(
                        period[tyx], rp_list[jj]["im_tip"][0, tyx]
                    ):
                        if tpi >= 0:
                            per_list_p.append(per)
                            tpi_list_p.append(tpi)
                        else:
                            per_list_n.append(per)
                            tpi_list_n.append(tpi)
                    if len(per_list_p) > 0:
                        m_line, s_line, b_line = axtipim.stem(
                            per_list_p, tpi_list_p, markerfmt="^", basefmt="k"
                        )
                        plt.setp(m_line, "markerfacecolor", self.ctipi)
                        plt.setp(m_line, "markeredgecolor", self.ctipi)
                        plt.setp(m_line, "markersize", self.ms)
                        plt.setp(s_line, "linewidth", self.lw)
                        plt.setp(s_line, "color", self.ctipi)
                        plt.setp(b_line, "linewidth", 0.01)
                        t_list.append(m_line)
                        t_label.append("Imag")
                    if len(per_list_n) > 0:
                        m_line, s_line, b_line = axtipim.stem(
                            per_list_n, tpi_list_n, markerfmt="v", basefmt="k"
                        )
                        plt.setp(m_line, "markerfacecolor", self.ctipi)
                        plt.setp(m_line, "markeredgecolor", self.ctipi)
                        plt.setp(m_line, "markersize", self.ms)
                        plt.setp(s_line, "linewidth", self.lw)
                        plt.setp(s_line, "color", self.ctipi)
                        plt.setp(b_line, "linewidth", 0.01)
                        if len(t_list) <= 1:
                            t_list.append(m_line)
                            t_label.append("Imag")

                else:
                    pass

            # ------------------- plot model response -------------------------
            if self.resp_fn is not None:
                num_resp = len(self.resp_fn)
                for rr, rfn in enumerate(self.resp_fn):
                    resp_obj = Response()
                    resp_obj.read_response_file(rfn)

                    rp = resp_obj.resp
                    # create colors for different responses
                    if self.color_mode == "color":
                        cxy = (0, 0.4 + float(rr) / (3 * num_resp), 0)
                        cyx = (
                            0.7 + float(rr) / (4 * num_resp),
                            0.13,
                            0.63 - float(rr) / (4 * num_resp),
                        )
                    elif self.color_mode == "bw":
                        cxy = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )
                        cyx = (
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                            1 - 1.25 / (rr + 2.0),
                        )

                    # calculate rms's
                    rmslistte = np.hstack(
                        (rp[jj]["te_res"][1], rp[jj]["te_phase"][1])
                    )
                    rmslisttm = np.hstack(
                        (rp[jj]["tm_res"][1], rp[jj]["tm_phase"][1])
                    )
                    rmste = np.sqrt(
                        np.sum([rms**2 for rms in rmslistte])
                        / len(rmslistte)
                    )
                    rmstm = np.sqrt(
                        np.sum([rms**2 for rms in rmslisttm])
                        / len(rmslisttm)
                    )

                    # ------------Plot Resistivity-----------------------------
                    # cut out missing data points first
                    # --> response
                    mrxy = np.where(rp[jj]["te_res"][0] != 0)[0]
                    mryx = np.where(rp[jj]["tm_res"][0] != 0)[0]

                    # --> TE mode Model Response
                    if len(mrxy) > 0:
                        r3 = plot_errorbar(
                            axrte,
                            period[mrxy],
                            rp[jj]["te_res"][0, mrxy],
                            ls="--",
                            marker=self.mtem,
                            ms=self.ms,
                            color=cxy,
                            y_error=None,
                            lw=self.lw,
                            e_capsize=self.e_capsize,
                            e_capthick=self.e_capthick,
                        )

                        rlistte.append(r3[0])
                        llistte.append("$Mod_{TE}$ " + "{0:.2f}".format(rmste))
                    else:
                        pass

                    # --> TM mode model response
                    if len(mryx) > 0:
                        r4 = plot_errorbar(
                            axrtm,
                            period[mryx],
                            rp[jj]["tm_res"][0, mryx],
                            ls="--",
                            marker=self.mtmm,
                            ms=self.ms,
                            color=cyx,
                            y_error=None,
                            lw=self.lw,
                            e_capsize=self.e_capsize,
                            e_capthick=self.e_capthick,
                        )
                        rlisttm.append(r4[0])
                        llisttm.append("$Mod_{TM}$ " + "{0:.2f}".format(rmstm))
                    else:
                        pass

                    # --------------------plot phase---------------------------
                    # cut out missing data points first
                    # --> reponse
                    mpxy = np.where(rp[jj]["te_phase"][0] != 0)[0]
                    mpyx = np.where(rp[jj]["tm_phase"][0] != 0)[0]

                    # --> TE mode response
                    if len(mpxy) > 0:
                        p3 = plot_errorbar(
                            axpte,
                            period[mpxy],
                            rp[jj]["te_phase"][0, mpxy],
                            ls="--",
                            ms=self.ms,
                            color=cxy,
                            y_error=None,
                            lw=self.lw,
                            e_capsize=self.e_capsize,
                            e_capthick=self.e_capthick,
                        )

                    else:
                        pass

                    # --> TM mode response
                    if len(mpyx) > 0:
                        p4 = plot_errorbar(
                            axptm,
                            period[mpyx],
                            rp[jj]["tm_phase"][0, mpyx],
                            ls="--",
                            marker=self.mtmm,
                            ms=self.ms,
                            color=cyx,
                            y_error=None,
                            lw=self.lw,
                            e_capsize=self.e_capsize,
                            e_capthick=self.e_capthick,
                        )
                    else:
                        pass

                    # ---------------------plot tipper-------------------------
                    if self.plot_tipper == "y":
                        txy = np.where(rp[jj]["re_tip"][0] != 0)[0]
                        tyx = np.where(rp[jj]["im_tip"][0] != 0)[0]
                        # --> real tipper  data
                        if len(txy) > 0:
                            per_list_p = []
                            tpr_list_p = []
                            per_list_n = []
                            tpr_list_n = []
                            for per, tpr in zip(
                                period[txy], rp[jj]["re_tip"][0, txy]
                            ):
                                if tpr >= 0:
                                    per_list_p.append(per)
                                    tpr_list_p.append(tpr)
                                else:
                                    per_list_n.append(per)
                                    tpr_list_n.append(tpr)
                            if len(per_list_p) > 0:
                                m_line, s_line, b_line = axtipre.stem(
                                    per_list_p,
                                    tpr_list_p,
                                    markerfmt="^",
                                    basefmt="k",
                                )
                                plt.setp(m_line, "markerfacecolor", cxy)
                                plt.setp(m_line, "markeredgecolor", cxy)
                                plt.setp(m_line, "markersize", self.ms)
                                plt.setp(s_line, "linewidth", self.lw)
                                plt.setp(s_line, "color", cxy)
                                plt.setp(b_line, "linewidth", 0.01)
                            if len(per_list_n) > 0:
                                m_line, s_line, b_line = axtipre.stem(
                                    per_list_n,
                                    tpr_list_n,
                                    markerfmt="v",
                                    basefmt="k",
                                )
                                plt.setp(m_line, "markerfacecolor", cxy)
                                plt.setp(m_line, "markeredgecolor", cxy)
                                plt.setp(m_line, "markersize", self.ms)
                                plt.setp(s_line, "linewidth", self.lw)
                                plt.setp(s_line, "color", cxy)
                                plt.setp(b_line, "linewidth", 0.01)

                        else:
                            pass
                        if len(tyx) > 0:
                            per_list_p = []
                            tpi_list_p = []
                            per_list_n = []
                            tpi_list_n = []
                            for per, tpi in zip(
                                period[tyx], rp[jj]["im_tip"][0, tyx]
                            ):
                                if tpi >= 0:
                                    per_list_p.append(per)
                                    tpi_list_p.append(tpi)
                                else:
                                    per_list_n.append(per)
                                    tpi_list_n.append(tpi)
                            if len(per_list_p) > 0:
                                m_line, s_line, b_line = axtipim.stem(
                                    per_list_p,
                                    tpi_list_p,
                                    markerfmt="^",
                                    basefmt="k",
                                )
                                plt.setp(m_line, "markerfacecolor", cyx)
                                plt.setp(m_line, "markeredgecolor", cyx)
                                plt.setp(m_line, "markersize", self.ms)
                                plt.setp(s_line, "linewidth", self.lw)
                                plt.setp(s_line, "color", cyx)
                                plt.setp(b_line, "linewidth", 0.01)
                            if len(per_list_n) > 0:
                                m_line, s_line, b_line = axtipim.stem(
                                    per_list_n,
                                    tpi_list_n,
                                    markerfmt="v",
                                    basefmt="k",
                                )
                                plt.setp(m_line, "markerfacecolor", cyx)
                                plt.setp(m_line, "markeredgecolor", cyx)
                                plt.setp(m_line, "markersize", self.ms)
                                plt.setp(s_line, "linewidth", self.lw)
                                plt.setp(s_line, "color", cyx)
                                plt.setp(b_line, "linewidth", 0.01)

                        else:
                            pass
            # --------------add in winglink responses------------------------
            if addwl == 1:
                try:
                    wlrms = wld[sdict[self.station_list[jj]]]["rms"]
                    axrte.set_title(
                        self.station_list[jj]
                        + "\n rms_occ_TE={0:.2f}".format(rmste)
                        + "rms_occ_TM={0:.2f}".format(rmstm)
                        + "rms_wl={0:.2f}".format(wlrms),
                        fontdict={"size": self.font_size, "weight": "bold"},
                    )
                    for ww, wlistation in enumerate(wlslist):
                        if wlistation.find(self.station_list[jj]) == 0:
                            print(
                                "{0} was Found {0} in winglink file".format(
                                    self.station_list[jj], wlistation
                                )
                            )
                            wlrpdict = wlrp_list[ww]

                    zrxy = [np.where(wlrpdict["te_res"][0] != 0)[0]]
                    zryx = [np.where(wlrpdict["tm_res"][0] != 0)[0]]

                    # plot winglink resistivity
                    r5 = axrte.loglog(
                        wlplist[zrxy],
                        wlrpdict["te_res"][1][zrxy],
                        ls="-.",
                        marker=self.mtewl,
                        ms=self.ms,
                        color=self.ctewl,
                        mfc=self.ctewl,
                        lw=self.lw,
                    )
                    r6 = axrtm.loglog(
                        wlplist[zryx],
                        wlrpdict["tm_res"][1][zryx],
                        ls="-.",
                        marker=self.mtmwl,
                        ms=self.ms,
                        color=self.ctmwl,
                        mfc=self.ctmwl,
                        lw=self.lw,
                    )

                    # plot winglink phase
                    axpte.semilogx(
                        wlplist[zrxy],
                        wlrpdict["te_phase"][1][zrxy],
                        ls="-.",
                        marker=self.mtewl,
                        ms=self.ms,
                        color=self.ctewl,
                        mfc=self.ctewl,
                        lw=self.lw,
                    )

                    axptm.semilogx(
                        wlplist[zryx],
                        wlrpdict["tm_phase"][1][zryx],
                        ls="-.",
                        marker=self.mtmwl,
                        ms=self.ms,
                        color=self.ctmwl,
                        mfc=self.ctmwl,
                        lw=self.lw,
                    )

                    rlistte.append(r5[0])
                    rlisttm.append(r6[0])
                    llistte.append("$WLMod_{TE}$ " + "{0:.2f}".format(wlrms))
                    llisttm.append("$WLMod_{TM}$ " + "{0:.2f}".format(wlrms))
                except (IndexError, KeyError):
                    print("Station not present")
            else:
                if self.plot_num == 1:
                    axrte.set_title(
                        self.station_list[jj],
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )
                elif self.plot_num == 2:
                    fig.suptitle(
                        self.station_list[jj],
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

            # set the axis properties
            ax_list = [axrte, axrtm]
            for aa, axr in enumerate(ax_list):
                # set both axes to logarithmic scale
                axr.set_xscale("log", nonposx="clip")

                try:
                    axr.set_yscale("log", nonposy="clip")
                except ValueError:
                    pass

                # put on a grid
                axr.grid(True, alpha=0.3, which="both", lw=0.5 * self.lw)
                axr.yaxis.set_label_coords(
                    self.ylabel_coord[0], self.ylabel_coord[1]
                )

                # set resistivity limits if desired
                if self.res_limits != None:
                    axr.set_ylim(
                        10 ** self.res_limits[0], 10 ** self.res_limits[1]
                    )

                # set the tick labels to invisible
                plt.setp(axr.xaxis.get_ticklabels(), visible=False)
                if aa == 0:
                    axr.set_ylabel(
                        "App. Res. ($\Omega \cdot m$)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

                # set legend based on the plot type
                if self.plot_num == 1:
                    if aa == 0:
                        axr.legend(
                            rlistte + rlisttm,
                            llistte + llisttm,
                            loc=2,
                            markerscale=1,
                            borderaxespad=0.05,
                            labelspacing=0.08,
                            handletextpad=0.15,
                            borderpad=0.05,
                            prop={"size": self.font_size + 1},
                        )
                elif self.plot_num == 2:
                    if aa == 0:
                        axr.legend(
                            rlistte,
                            llistte,
                            loc=2,
                            markerscale=1,
                            borderaxespad=0.05,
                            labelspacing=0.08,
                            handletextpad=0.15,
                            borderpad=0.05,
                            prop={"size": self.font_size + 1},
                        )

                    if aa == 1:
                        axr.legend(
                            rlisttm,
                            llisttm,
                            loc=2,
                            markerscale=1,
                            borderaxespad=0.05,
                            labelspacing=0.08,
                            handletextpad=0.15,
                            borderpad=0.05,
                            prop={"size": self.font_size + 1},
                        )

            # set Properties for the phase axes
            for aa, axp in enumerate([axpte, axptm]):
                # set the x-axis to log scale
                axp.set_xscale("log", nonposx="clip")

                # set the phase limits
                axp.set_ylim(self.phase_limits)

                # put a grid on the subplot
                axp.grid(True, alpha=0.3, which="both", lw=0.5 * self.lw)

                # set the tick locations
                axp.yaxis.set_major_locator(MultipleLocator(10))
                axp.yaxis.set_minor_locator(MultipleLocator(2))

                # set the x axis label
                if self.plot_tipper == "y":
                    plt.setp(axp.get_xticklabels(), visible=False)
                else:
                    axp.set_xlabel(
                        "Period (s)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

                # put the y label on the far left plot
                axp.yaxis.set_label_coords(
                    self.ylabel_coord[0], self.ylabel_coord[1]
                )
                if aa == 0:
                    axp.set_ylabel(
                        "Phase (deg)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

            # set axes properties of tipper axis
            if self.plot_tipper == "y":
                for aa, axt in enumerate([axtipre, axtipim]):
                    axt.set_xscale("log", nonposx="clip")

                    # set tipper limits
                    axt.set_ylim(self.tip_limits)

                    # put a grid on the subplot
                    axt.grid(True, alpha=0.3, which="both", lw=0.5 * self.lw)

                    # set the tick locations
                    axt.yaxis.set_major_locator(MultipleLocator(0.2))
                    axt.yaxis.set_minor_locator(MultipleLocator(0.1))

                    # set the x axis label
                    axt.set_xlabel(
                        "Period (s)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

                    axt.set_xlim(
                        10 ** np.floor(np.log10(data_obj.period.min())),
                        10 ** np.ceil(np.log10(data_obj.period.max())),
                    )

                    # put the y label on the far left plot
                    axt.yaxis.set_label_coords(
                        self.ylabel_coord[0], self.ylabel_coord[1]
                    )
                    if aa == 0:
                        axt.set_ylabel(
                            "Tipper",
                            fontdict={
                                "size": self.font_size + 2,
                                "weight": "bold",
                            },
                        )
                        if self.plot_num == 2:
                            axt.text(
                                axt.get_xlim()[0] * 1.25,
                                self.tip_limits[1] * 0.9,
                                "Real",
                                horizontalalignment="left",
                                verticalalignment="top",
                                bbox={"facecolor": "white"},
                                fontdict={"size": self.font_size + 1},
                            )
                        else:
                            axt.legend(
                                t_list,
                                t_label,
                                loc=2,
                                markerscale=1,
                                borderaxespad=0.05,
                                labelspacing=0.08,
                                handletextpad=0.15,
                                borderpad=0.05,
                                prop={"size": self.font_size + 1},
                            )
                    if aa == 1:
                        if self.plot_num == 2:
                            axt.text(
                                axt.get_xlim()[0] * 1.25,
                                self.tip_limits[1] * 0.9,
                                "Imag",
                                horizontalalignment="left",
                                verticalalignment="top",
                                bbox={"facecolor": "white"},
                                fontdict={"size": self.font_size + 1},
                            )

            # make sure the axis and figure are accessible to the user
            self.fig_list.append(
                {
                    "station": self.station_list[jj],
                    "fig": fig,
                    "axrte": axrte,
                    "axrtm": axrtm,
                    "axpte": axpte,
                    "axptm": axptm,
                }
            )

        # set the plot to be full screen well at least try
        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> #change color of te markers to a gray-blue
            >>> p1.cted = (.5, .5, .7)
            >>> p1.redraw_plot()
        """

        plt.close("all")
        self.plot()

    def save_figures(
        self, save_path, fig_fmt="pdf", fig_dpi=None, close_fig="y"
    ):
        """
        save all the figure that are in self.fig_list

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plot2DResponses()
            >>> p1.save_figures(r"/home/occam2d/Figures", fig_fmt='jpg')
        """

        if not os.path.exists(save_path):
            os.mkdir(save_path)

        for fdict in self.fig_list:
            svfn = "{0}_resp.{1}".format(fdict["station"], fig_fmt)
            fdict["fig"].savefig(
                os.path.join(save_path, svfn), dpi=self.fig_dpi
            )
            if close_fig == "y":
                plt.close(fdict["fig"])

            print("saved figure to {0}".format(os.path.join(save_path, svfn)))

            # ==============================================================================


# plot model
# ==============================================================================
class PlotModel(Model):
    """
    plot the 2D model found by Occam2D.  The model is displayed as a meshgrid
    instead of model bricks.  This speeds things up considerably.

    Inherets the Model class to take advantage of the attributes and methods
    already coded.

    Arguments:
    -----------
        **iter_fn** : string
                      full path to iteration file.  From here all the
                      necessary files can be found assuming they are in the
                      same directory.  If they are not then need to input
                      manually.


    ======================= ===============================================
    keywords                description
    ======================= ===============================================
    block_font_size         font size of block number is blocknum == 'on'
    blocknum                [ 'on' | 'off' ] to plot regulariztion block
                            numbers.
    cb_pad                  padding between axes edge and color bar
    cb_shrink               percentage to shrink the color bar
    climits                 limits of the color scale for resistivity
                            in log scale (min, max)
    cmap                    name of color map for resistivity values
    femesh                  plot the finite element mesh
    femesh_triangles        plot the finite element mesh with each block
                            divided into four triangles
    fig_aspect              aspect ratio between width and height of
                            resistivity image. 1 for equal axes
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    grid                    [ 'both' | 'major' |'minor' | None ] string
                            to tell the program to make a grid on the
                            specified axes.
    meshnum                 [ 'on' | 'off' ] 'on' will plot finite element
                            mesh numbers
    meshnum_font_size       font size of mesh numbers if meshnum == 'on'
    ms                      size of station marker
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    regmesh                 [ 'on' | 'off' ] plot the regularization mesh
                            plots as blue lines
    station_color           color of station marker
    station_font_color      color station label
    station_font_pad        padding between station label and marker
    station_font_rotation   angle of station label in degrees 0 is
                            horizontal
    station_font_size       font size of station label
    station_font_weight     font weight of station label
    station_id              index to take station label from station name
    station_marker          station marker.  if inputing a LaTex marker
                            be sure to input as r"LaTexMarker" otherwise
                            might not plot properly
    subplot_bottom          subplot spacing from bottom
    subplot_left            subplot spacing from left
    subplot_right           subplot spacing from right
    subplot_top             subplot spacing from top
    title                   title of plot.  If None then the name of the
                            iteration file and containing folder will be
                            the title with RMS and Roughness.
    xlimits                 limits of plot in x-direction in (km)
    xminorticks             increment of minor ticks in x direction
    xpad                    padding in x-direction in km
    ylimits                 depth limits of plot positive down (km)
    yminorticks             increment of minor ticks in y-direction
    ypad                    padding in negative y-direction (km)
    yscale                  [ 'km' | 'm' ] scale of plot, if 'm' everything
                            will be scaled accordingly.
    ======================= ===============================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots resistivity model.
    redraw_plot         call redraw_plot to redraw the figures,
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired
                        location and format
    =================== ====================================================

    :Example:
    ---------------
        >>> import mtpy.modeling.occam2d as occam2d
        >>> model_plot = occam2d.PlotModel(r"/home/occam/Inv1/mt_01.iter")
        >>> # change the color limits
        >>> model_plot.climits = (1, 4)
        >>> model_plot.redraw_plot()
        >>> #change len of station name
        >>> model_plot.station_id = [2, 5]
        >>> model_plot.redraw_plot()


    """

    def __init__(self, iter_fn=None, data_fn=None, **kwargs):
        Model.__init__(self, iter_fn, **kwargs)

        self.yscale = kwargs.pop("yscale", "km")

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)
        self.fig_aspect = kwargs.pop("fig_aspect", 1)
        self.title = kwargs.pop("title", "on")

        self.xpad = kwargs.pop("xpad", 1.0)
        self.ypad = kwargs.pop("ypad", 1.0)

        self.ms = kwargs.pop("ms", 10)

        self.station_locations = None
        self.station_list = None
        self.station_id = kwargs.pop("station_id", None)
        self.station_font_size = kwargs.pop("station_font_size", 8)
        self.station_font_pad = kwargs.pop("station_font_pad", 1.0)
        self.station_font_weight = kwargs.pop("station_font_weight", "bold")
        self.station_font_rotation = kwargs.pop("station_font_rotation", 60)
        self.station_font_color = kwargs.pop("station_font_color", "k")
        self.station_marker = kwargs.pop(
            "station_marker", r"$\blacktriangledown$"
        )
        self.station_color = kwargs.pop("station_color", "k")

        self.ylimits = kwargs.pop("ylimits", None)
        self.xlimits = kwargs.pop("xlimits", None)

        self.xminorticks = kwargs.pop("xminorticks", 5)
        self.yminorticks = kwargs.pop("yminorticks", 1)

        self.climits = kwargs.pop("climits", (0, 4))
        self.cmap = kwargs.pop("cmap", "jet_r")
        self.font_size = kwargs.pop("font_size", 8)

        self.femesh = kwargs.pop("femesh", "off")
        self.femesh_triangles = kwargs.pop("femesh_triangles", "off")
        self.femesh_lw = kwargs.pop("femesh_lw", 0.4)
        self.femesh_color = kwargs.pop("femesh_color", "k")
        self.meshnum = kwargs.pop("meshnum", "off")
        self.meshnum_font_size = kwargs.pop("meshnum_font_size", 3)

        self.regmesh = kwargs.pop("regmesh", "off")
        self.regmesh_lw = kwargs.pop("regmesh_lw", 0.4)
        self.regmesh_color = kwargs.pop("regmesh_color", "b")
        self.blocknum = kwargs.pop("blocknum", "off")
        self.block_font_size = kwargs.pop("block_font_size", 3)
        self.grid = kwargs.pop("grid", None)

        self.cb_shrink = kwargs.pop("cb_shrink", 0.8)
        self.cb_pad = kwargs.pop("cb_pad", 0.01)

        self.subplot_right = 0.99
        self.subplot_left = 0.085
        self.subplot_top = 0.92
        self.subplot_bottom = 0.1

        self.plot_yn = kwargs.pop("plot_yn", "y")
        if self.plot_yn == "y":
            self.plot()

    def plot(self):
        """
        plotModel will plot the model output by occam2d in the iteration file.


        :Example: ::

            >>> import mtpy.modeling.occam2d as occam2d
            >>> itfn = r"/home/Occam2D/Line1/Inv1/Test_15.iter"
            >>> model_plot = occam2d.PlotModel(itfn)
            >>> model_plot.ms = 20
            >>> model_plot.ylimits = (0,.350)
            >>> model_plot.yscale = 'm'
            >>> model_plot.spad = .10
            >>> model_plot.ypad = .125
            >>> model_plot.xpad = .025
            >>> model_plot.climits = (0,2.5)
            >>> model_plot.aspect = 'equal'
            >>> model_plot.redraw_plot()

        """
        # --> read in iteration file and build the model
        self.read_iter_file()
        self.build_model()

        # --> get station locations and names from data file
        d_object = Data()
        d_object.read_data_file(self.data_fn)
        setattr(self, "station_locations", d_object.station_locations.copy())
        setattr(self, "station_list", d_object.station_list.copy())

        # set the scale of the plot
        if self.yscale == "km":
            df = 1000.0
            pf = 1.0
        elif self.yscale == "m":
            df = 1.0
            pf = 1000.0
        else:
            df = 1000.0
            pf = 1.0

        # set some figure properties to use the maiximum space
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top

        # station font dictionary
        fdict = {
            "size": self.station_font_size,
            "weight": self.station_font_weight,
            "rotation": self.station_font_rotation,
            "color": self.station_font_color,
        }

        # plot the model as a mesh
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # add a subplot to the figure with the specified aspect ratio
        ax = self.fig.add_subplot(1, 1, 1, aspect=self.fig_aspect)

        # plot the model as a pcolormesh so the extents are constrained to
        # the model coordinates
        ax.pcolormesh(
            self.mesh_x / df,
            self.mesh_z / df,
            self.res_model,
            cmap=self.cmap,
            vmin=self.climits[0],
            vmax=self.climits[1],
        )

        # make a colorbar for the resistivity
        cbx = mcb.make_axes(ax, shrink=self.cb_shrink, pad=self.cb_pad)
        cb = mcb.ColorbarBase(
            cbx[0],
            cmap=self.cmap,
            norm=Normalize(vmin=self.climits[0], vmax=self.climits[1]),
        )

        cb.set_label(
            "Resistivity ($\Omega \cdot$m)",
            fontdict={"size": self.font_size + 1, "weight": "bold"},
        )
        cb.set_ticks(np.arange(int(self.climits[0]), int(self.climits[1]) + 1))
        cb.set_ticklabels(
            [
                "10$^{0}$".format("{" + str(nn) + "}")
                for nn in np.arange(
                    int(self.climits[0]), int(self.climits[1]) + 1
                )
            ]
        )

        # set the offsets of the stations and plot the stations
        # need to figure out a way to set the marker at the surface in all
        # views.
        for offset, name in zip(self.station_locations, self.station_list):
            # plot the station marker
            # plots a V for the station cause when you use scatter the spacing
            # is variable if you change the limits of the y axis, this way it
            # always plots at the surface.
            ax.text(
                offset / df,
                self.plot_z.min(),
                self.station_marker,
                horizontalalignment="center",
                verticalalignment="baseline",
                fontdict={"size": self.ms, "color": self.station_color},
            )

            # put station id onto station marker
            # if there is a station id index
            if self.station_id != None:
                ax.text(
                    offset / df,
                    -self.station_font_pad * pf,
                    name[self.station_id[0] : self.station_id[1]],
                    horizontalalignment="center",
                    verticalalignment="baseline",
                    fontdict=fdict,
                )
            # otherwise put on the full station name found form data file
            else:
                ax.text(
                    offset / df,
                    -self.station_font_pad * pf,
                    name,
                    horizontalalignment="center",
                    verticalalignment="baseline",
                    fontdict=fdict,
                )

        # set the initial limits of the plot to be square about the profile
        # line
        if self.ylimits == None:
            ax.set_ylim(
                abs(
                    self.station_locations.max() - self.station_locations.min()
                )
                / df,
                -self.ypad * pf,
            )
        else:
            ax.set_ylim(
                self.ylimits[1] * pf, (self.ylimits[0] - self.ypad) * pf
            )
        if self.xlimits == None:
            ax.set_xlim(
                self.station_locations.min() / df - (self.xpad * pf),
                self.station_locations.max() / df + (self.xpad * pf),
            )
        else:
            ax.set_xlim(self.xlimits[0] * pf, self.xlimits[1] * pf)

        # set the axis properties
        ax.xaxis.set_minor_locator(MultipleLocator(self.xminorticks * pf))
        ax.yaxis.set_minor_locator(MultipleLocator(self.yminorticks * pf))

        # set axes labels
        ax.set_xlabel(
            "Horizontal Distance ({0})".format(self.yscale),
            fontdict={"size": self.font_size + 2, "weight": "bold"},
        )
        ax.set_ylabel(
            "Depth ({0})".format(self.yscale),
            fontdict={"size": self.font_size + 2, "weight": "bold"},
        )

        # put a grid on if one is desired
        if self.grid is not None:
            ax.grid(alpha=0.3, which=self.grid, lw=0.35)

        # set title as rms and roughness
        if type(self.title) is str:
            if self.title == "on":
                titlestr = os.path.join(
                    os.path.basename(os.path.dirname(self.iter_fn)),
                    os.path.basename(self.iter_fn),
                )
                ax.set_title(
                    "{0}: RMS={1:.2f}, Roughness={2:.0f}".format(
                        titlestr, self.misfit_value, self.roughness_value
                    ),
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )
            else:
                ax.set_title(
                    "{0}; RMS={1:.2f}, Roughness={2:.0f}".format(
                        self.title, self.misfit_value, self.roughness_value
                    ),
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )
        else:
            print(
                "RMS {0:.2f}, Roughness={1:.0f}".format(
                    self.misfit_value, self.roughness_value
                )
            )

            # plot forward model mesh
        # making an extended list seperated by None's speeds up the plotting
        # by as much as 99 percent, handy
        if self.femesh == "on":
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.plot_x / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.plot_zy[0] / df])
                row_line_ylist.append(None)

            # plot column lines (variables are a little bit of a misnomer)
            ax.plot(row_line_xlist, row_line_ylist, color="k", lw=0.5)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.plot_z / df:
                col_line_xlist.extend(
                    [self.plot_x[0] / df, self.plot_x[-1] / df]
                )
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot row lines (variables are a little bit of a misnomer)
            ax.plot(col_line_xlist, col_line_ylist, color="k", lw=0.5)

        if self.femesh_triangles == "on":
            row_line_xlist = []
            row_line_ylist = []
            for xx in self.plot_x / df:
                row_line_xlist.extend([xx, xx])
                row_line_xlist.append(None)
                row_line_ylist.extend([0, self.plot_z[0] / df])
                row_line_ylist.append(None)

            # plot columns
            ax.plot(row_line_xlist, row_line_ylist, color="k", lw=0.5)

            col_line_xlist = []
            col_line_ylist = []
            for yy in self.plot_z / df:
                col_line_xlist.extend(
                    [self.plot_x[0] / df, self.plot_x[-1] / df]
                )
                col_line_xlist.append(None)
                col_line_ylist.extend([yy, yy])
                col_line_ylist.append(None)

            # plot rows
            ax.plot(col_line_xlist, col_line_ylist, color="k", lw=0.5)

            diag_line_xlist = []
            diag_line_ylist = []
            for xi, xx in enumerate(self.plot_x[:-1] / df):
                for yi, yy in enumerate(self.plot_z[:-1] / df):
                    diag_line_xlist.extend([xx, self.plot_x[xi + 1] / df])
                    diag_line_xlist.append(None)
                    diag_line_xlist.extend([xx, self.plot_x[xi + 1] / df])
                    diag_line_xlist.append(None)

                    diag_line_ylist.extend([yy, self.plot_z[yi + 1] / df])
                    diag_line_ylist.append(None)
                    diag_line_ylist.extend([self.plot_z[yi + 1] / df, yy])
                    diag_line_ylist.append(None)

            # plot diagonal lines.
            ax.plot(diag_line_xlist, diag_line_ylist, color="k", lw=0.5)

        # plot the regularization mesh
        if self.regmesh == "on":
            line_list = []
            for ii in range(len(self.model_rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.model_rows[:ii, 0].sum()

                # the second index  in the vertical direction
                ny2 = ny1 + self.model_rows[ii][0]

                # make the list of amalgamated columns an array for ease
                lc = np.array(self.model_cols[ii])
                yline = ax.plot(
                    [self.plot_x[0] / df, self.plot_x[-1] / df],
                    [self.plot_z[-ny1] / df, self.plot_z[-ny1] / df],
                    color="b",
                    lw=0.5,
                )

                line_list.append(yline)

                # loop over the number of amalgamated blocks
                for jj in range(len(self.model_cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()

                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        xline = ax.plot(
                            [self.plot_x[nx1] / df, self.plot_x[nx1] / df],
                            [self.plot_z[-ny1] / df, self.plot_z[-ny2] / df],
                            color="b",
                            lw=0.5,
                        )
                        line_list.append(xline)
                    except IndexError:
                        pass

        # plot the mesh block numbers
        if self.meshnum == "on":
            kk = 1
            for yy in self.plot_z[::-1] / df:
                for xx in self.plot_x / df:
                    ax.text(
                        xx,
                        yy,
                        "{0}".format(kk),
                        fontdict={"size": self.meshnum_font_size},
                    )
                    kk += 1

        # plot regularization block numbers
        if self.blocknum == "on":
            kk = 1
            for ii in range(len(self.model_rows)):
                # get the number of layers to combine
                # this index will be the first index in the vertical direction
                ny1 = self.model_rows[:ii, 0].sum()

                # the second index  in the vertical direction
                ny2 = ny1 + self.model_rows[ii][0]
                # make the list of amalgamated columns an array for ease
                lc = np.array(self.model_cols[ii])
                # loop over the number of amalgamated blocks
                for jj in range(len(self.model_cols[ii])):
                    # get first in index in the horizontal direction
                    nx1 = lc[:jj].sum()
                    # get second index in horizontal direction
                    nx2 = nx1 + lc[jj]
                    try:
                        if ny1 == 0:
                            ny1 = 1
                        # get center points of the blocks
                        yy = (
                            self.plot_z[-ny1]
                            - (self.plot_z[-ny1] - self.plot_z[-ny2]) / 2
                        )
                        xx = (
                            self.plot_x[nx1]
                            - (self.plot_x[nx1] - self.plot_x[nx2]) / 2
                        )
                        # put the number
                        ax.text(
                            xx / df,
                            yy / df,
                            "{0}".format(kk),
                            fontdict={"size": self.block_font_size},
                            horizontalalignment="center",
                            verticalalignment="center",
                        )
                        kk += 1
                    except IndexError:
                        pass

        plt.show()

        # make attributes that can be manipulated
        self.ax = ax
        self.cbax = cb

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

        plt.close(self.fig)
        self.plot()

    def save_figure(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_fig="y",
    ):
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
            >>> model_plot.save_figure(r"/home/occam/figures",
                                       file_format='jpg')

        """

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        else:
            save_fn = os.path.join(save_fn, "OccamModel." + file_format)
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        if close_fig == "y":
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print("Saved figure to: " + self.fig_fn)

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

        return "Plots the resistivity model found by Occam2D."


# ==============================================================================
# plot L2 curve of iteration vs rms
# ==============================================================================
class PlotL2:
    """
    Plot L2 curve of iteration vs rms and rms vs roughness.

    Need to only input an .iter file, will read all similar .iter files
    to get the rms, iteration number and roughness of all similar .iter files.

    Arguments:
    ----------
        **iter_fn** : string
                      full path to an iteration file output by Occam2D.

    ======================= ===================================================
    Keywords/attributes     Description
    ======================= ===================================================
    ax1                     matplotlib.axes instance for rms vs iteration
    ax2                     matplotlib.axes instance for roughness vs rms
    fig                     matplotlib.figure instance
    fig_dpi                 resolution of figure in dots-per-inch
    fig_num                 number of figure instance
    fig_size                size of figure in inches (width, height)
    font_size               size of axes tick labels, axes labels is +2
    plot_yn                 [ 'y' | 'n']
                            'y' --> to plot on instantiation
                            'n' --> to not plot on instantiation
    rms_arr                 structure np.array as described above
    rms_color               color of rms marker and line
    rms_lw                  line width of rms line
    rms_marker              marker for rms values
    rms_marker_size         size of marker for rms values
    rms_mean_color          color of mean line
    rms_median_color        color of median line
    rough_color             color of roughness line and marker
    rough_font_size         font size for iteration number inside roughness
                            marker
    rough_lw                line width for roughness line
    rough_marker            marker for roughness
    rough_marker_size       size of marker for roughness
    subplot_bottom          subplot spacing from bottom
    subplot_left            subplot spacing from left
    subplot_right           subplot spacing from right
    subplot_top             subplot spacing from top
    ======================= ===================================================

    =================== =======================================================
    Methods             Description
    =================== =======================================================
    plot                plots L2 curve.
    redraw_plot         call redraw_plot to redraw the figures,
                        if one of the attributes has been changed
    save_figure         saves the matplotlib.figure instance to desired
                        location and format
    =================== ======================================================

    """

    def __init__(self, iter_fn, **kwargs):
        self.iter_path = os.path.dirname(iter_fn)
        self.iter_basename = os.path.basename(iter_fn)[:-7]
        self.iter_fn_list = []
        self.rms_arr = None
        self.rough_arr = None

        self.subplot_right = 0.98
        self.subplot_left = 0.085
        self.subplot_top = 0.91
        self.subplot_bottom = 0.1

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)
        self.font_size = kwargs.pop("font_size", 8)

        self.rms_lw = kwargs.pop("rms_lw", 1)
        self.rms_marker = kwargs.pop("rms_marker", "d")
        self.rms_color = kwargs.pop("rms_color", "k")
        self.rms_marker_size = kwargs.pop("rms_marker_size", 5)
        self.rms_median_color = kwargs.pop("rms_median_color", "red")
        self.rms_mean_color = kwargs.pop("rms_mean_color", "orange")

        self.rough_lw = kwargs.pop("rough_lw", 0.75)
        self.rough_marker = kwargs.pop("rough_marker", "o")
        self.rough_color = kwargs.pop("rough_color", "b")
        self.rough_marker_size = kwargs.pop("rough_marker_size", 7)
        self.rough_font_size = kwargs.pop("rough_font_size", 6)

        self.plot_yn = kwargs.pop("plot_yn", "y")
        if self.plot_yn == "y":
            self.plot()

    def _get_iterfn_list(self):
        """
        get all iteration files for a given inversion

        """

        self.iter_fn_list = [
            os.path.join(self.iter_path, fn)
            for fn in os.listdir(self.iter_path)
            if fn.find(self.iter_basename) == 0 and fn.find(".iter") > 0
        ]

    def _get_values(self):
        """
        get rms and roughness values from iteration files
        """
        self._get_iterfn_list()
        self.rms_arr = np.zeros((len(self.iter_fn_list), 2))
        self.rough_arr = np.zeros((len(self.iter_fn_list), 2))

        for ii, itfn in enumerate(self.iter_fn_list):
            m_object = Model(itfn)
            m_object.read_iter_file()
            m_index = int(m_object.iteration)
            self.rms_arr[ii, 1] = float(m_object.misfit_value)
            self.rms_arr[ii, 0] = m_index
            self.rough_arr[ii, 1] = float(m_object.roughness_value)
            self.rough_arr[ii, 0] = m_index

            # sort by iteration number
            #        self.rms_arr = np.sort(self.rms_arr, axis=1)
            #        self.rough_arr = np.sort(self.rough_arr, axis=1)

    def plot(self):
        """
        plot L2 curve
        """

        self._get_values()
        nr = self.rms_arr.shape[0]
        med_rms = np.median(self.rms_arr[1:, 1])
        mean_rms = np.mean(self.rms_arr[1:, 1])

        # set the dimesions of the figure
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top

        # make figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        # make a subplot for RMS vs Iteration
        self.ax1 = self.fig.add_subplot(1, 1, 1)

        # plot the rms vs iteration
        (l1,) = self.ax1.plot(
            self.rms_arr[:, 0],
            self.rms_arr[:, 1],
            "-k",
            lw=1,
            marker="d",
            ms=5,
        )

        # plot the median of the RMS
        (m1,) = self.ax1.plot(
            self.rms_arr[:, 0],
            np.repeat(med_rms, nr),
            ls="--",
            color=self.rms_median_color,
            lw=self.rms_lw * 0.75,
        )

        # plot the mean of the RMS
        (m2,) = self.ax1.plot(
            self.rms_arr[:, 0],
            np.repeat(mean_rms, nr),
            ls="--",
            color=self.rms_mean_color,
            lw=self.rms_lw * 0.75,
        )

        # make subplot for RMS vs Roughness Plot
        self.ax2 = self.ax1.twiny()

        self.ax2.set_xlim(
            self.rough_arr[1:, 1].min(), self.rough_arr[1:, 1].max()
        )

        self.ax1.set_ylim(
            np.floor(self.rms_arr[1:, 1].min()), self.rms_arr[1:, 1].max()
        )

        # plot the rms vs roughness
        (l2,) = self.ax2.plot(
            self.rough_arr[:, 1],
            self.rms_arr[:, 1],
            ls="--",
            color=self.rough_color,
            lw=self.rough_lw,
            marker=self.rough_marker,
            ms=self.rough_marker_size,
            mfc="white",
        )

        # plot the iteration number inside the roughness marker
        for rms, ii, rough in zip(
            self.rms_arr[:, 1], self.rms_arr[:, 0], self.rough_arr[:, 1]
        ):
            # need this because if the roughness is larger than this number
            # matplotlib puts the text out of bounds and a draw_text_image
            # error is raised and file cannot be saved, also the other
            # numbers are not put in.
            if rough > 1e8:
                pass
            else:
                self.ax2.text(
                    rough,
                    rms,
                    "{0:.0f}".format(ii),
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontdict={
                        "size": self.rough_font_size,
                        "weight": "bold",
                        "color": self.rough_color,
                    },
                )

        # make a legend
        self.ax1.legend(
            [l1, l2, m1, m2],
            [
                "RMS",
                "Roughness",
                "Median_RMS={0:.2f}".format(med_rms),
                "Mean_RMS={0:.2f}".format(mean_rms),
            ],
            ncol=1,
            loc="upper right",
            columnspacing=0.25,
            markerscale=0.75,
            handletextpad=0.15,
        )

        # set the axis properties for RMS vs iteration
        self.ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
        self.ax1.xaxis.set_minor_locator(MultipleLocator(1))
        self.ax1.xaxis.set_major_locator(MultipleLocator(1))
        self.ax1.set_ylabel(
            "RMS", fontdict={"size": self.font_size + 2, "weight": "bold"}
        )
        self.ax1.set_xlabel(
            "Iteration",
            fontdict={"size": self.font_size + 2, "weight": "bold"},
        )
        self.ax1.grid(alpha=0.25, which="both", lw=self.rough_lw)
        self.ax2.set_xlabel(
            "Roughness",
            fontdict={
                "size": self.font_size + 2,
                "weight": "bold",
                "color": self.rough_color,
            },
        )

        for t2 in self.ax2.get_xticklabels():
            t2.set_color(self.rough_color)

        plt.show()

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

        plt.close(self.fig)
        self.plot()

    def save_figure(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_fig="y",
    ):
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

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        else:
            save_fn = os.path.join(save_fn, "_L2." + file_format)
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        if close_fig == "y":
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print("Saved figure to: " + self.fig_fn)

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

        return "Plots RMS vs Iteration computed by Occam2D"


# ==============================================================================
# plot pseudo section of data and model response
# ==============================================================================
class PlotPseudoSection(object):
    """
     plot a pseudo section of the data and response if given


     Arguments:
     -------------
         **data_fn** : string
                       full path to data file.

         **resp_fn** : string
                       full path to response file.

     ==================== ======================================================
     key words            description
     ==================== ======================================================
     axmpte               matplotlib.axes instance for TE model phase
     axmptm               matplotlib.axes instance for TM model phase
     axmrte               matplotlib.axes instance for TE model app. res
     axmrtm               matplotlib.axes instance for TM model app. res
     axpte                matplotlib.axes instance for TE data phase
     axptm                matplotlib.axes instance for TM data phase
     axrte                matplotlib.axes instance for TE data app. res.
     axrtm                matplotlib.axes instance for TM data app. res.
     cb_pad               padding between colorbar and axes
     cb_shrink            percentage to shrink the colorbar to
     fig                  matplotlib.figure instance
     fig_dpi              resolution of figure in dots per inch
     fig_num              number of figure instance
     fig_size             size of figure in inches (width, height)
     font_size            size of font in points
     label_list            list to label plots
     ml                   factor to label stations if 2 every other station
                          is labeled on the x-axis
     period               np.array of periods to plot
     phase_cmap           color map name of phase
     phase_limits_te      limits for te phase in degrees (min, max)
     phase_limits_tm      limits for tm phase in degrees (min, max)
     plot_resp            [ 'y' | 'n' ] to plot response
     plot_tipper          [ 'y' | 'n' ] to plot tipper
     plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation
     res_cmap             color map name for resistivity
     res_limits_te        limits for te resistivity in log scale (min, max)
     res_limits_tm        limits for tm resistivity in log scale (min, max)
     rp_list               list of dictionaries as made from read2Dresp
     station_id           index to get station name (min, max)
     station_list          station list got from rp_list
     subplot_bottom       subplot spacing from bottom (relative coordinates)
     subplot_hspace       vertical spacing between subplots
     subplot_left         subplot spacing from left
     subplot_right        subplot spacing from right
     subplot_top          subplot spacing from top
     subplot_wspace       horizontal spacing between subplots
     ==================== ======================================================

     =================== =======================================================
     Methods             Description
     =================== =======================================================
     plot                plots a pseudo-section of apparent resistiviy and phase
                         of data and model if given.  called on instantiation
                         if plot_yn is 'y'.
     redraw_plot         call redraw_plot to redraw the figures,
                         if one of the attributes has been changed
     save_figure         saves the matplotlib.figure instance to desired
                         location and format
     =================== =======================================================

    :Example: ::

         >>> import mtpy.modeling.occam2d as occam2d
         >>> r_fn = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
         >>> d_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
         >>> ps_plot = occam2d.PlotPseudoSection(d_fn, r_fn)

    """

    def __init__(self, data_fn, resp_fn=None, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.plot_resp = kwargs.pop("plot_resp", "y")
        if self.resp_fn is None:
            self.plot_resp = "n"

        self.label_list = [
            r"$\rho_{TE-Data}$",
            r"$\rho_{TE-Model}$",
            r"$\rho_{TM-Data}$",
            r"$\rho_{TM-Model}$",
            "$\phi_{TE-Data}$",
            "$\phi_{TE-Model}$",
            "$\phi_{TM-Data}$",
            "$\phi_{TM-Model}$",
            "$\Re e\{T_{Data}\}$",
            "$\Re e\{T_{Model}\}$",
            "$\Im m\{T_{Data}\}$",
            "$\Im m\{T_{Model}\}$",
        ]

        self.phase_limits_te = kwargs.pop("phase_limits_te", (-5, 95))
        self.phase_limits_tm = kwargs.pop("phase_limits_tm", (-5, 95))
        self.res_limits_te = kwargs.pop("res_limits_te", (0, 3))
        self.res_limits_tm = kwargs.pop("res_limits_tm", (0, 3))
        self.tip_limits_re = kwargs.pop("tip_limits_re", (-1, 1))
        self.tip_limits_im = kwargs.pop("tip_limits_im", (-1, 1))

        self.phase_cmap = kwargs.pop("phase_cmap", "jet")
        self.res_cmap = kwargs.pop("res_cmap", "jet_r")
        self.tip_cmap = kwargs.pop("res_cmap", "Spectral")

        self.ml = kwargs.pop("ml", 2)
        self.station_id = kwargs.pop("station_id", [0, 4])

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)

        self.subplot_wspace = 0.025
        self.subplot_hspace = 0.0
        self.subplot_right = 0.95
        self.subplot_left = 0.085
        self.subplot_top = 0.97
        self.subplot_bottom = 0.1

        self.font_size = kwargs.pop("font_size", 6)

        self.plot_type = kwargs.pop("plot_type", "1")
        self.plot_num = kwargs.pop("plot_num", 2)
        self.plot_tipper = kwargs.pop("plot_tipper", "n")
        self.plot_yn = kwargs.pop("plot_yn", "y")

        self.cb_shrink = 0.7
        self.cb_pad = 0.015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axmrte = None
        self.axmrtm = None
        self.axmpte = None
        self.axmptm = None
        self.axtpr = None
        self.axtpi = None
        self.axmtpr = None
        self.axmtpi = None

        self.te_res_arr = None
        self.tm_res_arr = None
        self.te_phase_arr = None
        self.tm_phase_arr = None
        self.tip_real_arr = None
        self.tip_imag_arr = None

        self.fig = None

        if self.plot_yn == "y":
            self.plot()

    def plot(self):
        """
        plot pseudo section of data and response if given

        """
        if self.plot_resp == "y":
            nr = 2
        else:
            nr = 1

        data_obj = Data()
        data_obj.read_data_file(self.data_fn)

        if self.resp_fn is not None:
            resp_obj = Response()
            resp_obj.read_response_file(self.resp_fn)

        ns = len(data_obj.station_list)
        nf = len(data_obj.period)
        ylimits = (data_obj.period.max(), data_obj.period.min())

        # make a grid for pcolormesh so you can have a log scale
        # get things into arrays for plotting
        offset_list = np.zeros(ns + 1)
        te_res_arr = np.ones((nf, ns, nr))
        tm_res_arr = np.ones((nf, ns, nr))
        te_phase_arr = np.zeros((nf, ns, nr))
        tm_phase_arr = np.zeros((nf, ns, nr))
        tip_real_arr = np.zeros((nf, ns, nr))
        tip_imag_arr = np.zeros((nf, ns, nr))

        for ii, d_dict in enumerate(data_obj.data):
            offset_list[ii] = d_dict["offset"]
            te_res_arr[:, ii, 0] = d_dict["te_res"][0]
            tm_res_arr[:, ii, 0] = d_dict["tm_res"][0]
            te_phase_arr[:, ii, 0] = d_dict["te_phase"][0]
            tm_phase_arr[:, ii, 0] = d_dict["tm_phase"][0]
            tip_real_arr[:, ii, 0] = d_dict["re_tip"][0]
            tip_imag_arr[:, ii, 0] = d_dict["im_tip"][0]

        # read in response data
        if self.plot_resp == "y":
            for ii, r_dict in enumerate(resp_obj.resp):
                te_res_arr[:, ii, 1] = r_dict["te_res"][0]
                tm_res_arr[:, ii, 1] = r_dict["tm_res"][0]
                te_phase_arr[:, ii, 1] = r_dict["te_phase"][0]
                tm_phase_arr[:, ii, 1] = r_dict["tm_phase"][0]
                tip_real_arr[:, ii, 1] = r_dict["re_tip"][0]
                tip_imag_arr[:, ii, 1] = r_dict["im_tip"][0]

        # need to make any zeros 1 for taking log10
        te_res_arr[np.where(te_res_arr == 0)] = 1.0
        tm_res_arr[np.where(tm_res_arr == 0)] = 1.0

        self.te_res_arr = te_res_arr
        self.tm_res_arr = tm_res_arr
        self.te_phase_arr = te_phase_arr
        self.tm_phase_arr = tm_phase_arr
        self.tip_real_arr = tip_real_arr
        self.tip_imag_arr = tip_imag_arr

        # need to extend the last grid cell because meshgrid expects n+1 cells
        offset_list[-1] = offset_list[-2] * 1.15
        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, data_obj.period[::-1])

        # make list for station labels
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [
            data_obj.station_list[ss][sindex_1:sindex_2]
            for ss in range(0, ns, self.ml)
        ]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * data_obj.period[1]

        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.left"] = self.subplot_left

        log_labels_te = [
            "10$^{0}$".format("{" + str(nn) + "}")
            for nn in np.arange(
                int(self.res_limits_te[0]), int(self.res_limits_te[1]) + 1
            )
        ]
        log_labels_tm = [
            "10$^{0}$".format("{" + str(nn) + "}")
            for nn in np.arange(
                int(self.res_limits_tm[0]), int(self.res_limits_tm[1]) + 1
            )
        ]

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_resp == "y":
            if self.plot_tipper == "y":
                gs1 = gridspec.GridSpec(
                    1,
                    3,
                    left=self.subplot_left,
                    right=self.subplot_right,
                    wspace=self.subplot_wspace,
                )
                gs4 = gridspec.GridSpecFromSubplotSpec(
                    2,
                    2,
                    hspace=self.subplot_hspace,
                    wspace=0,
                    subplot_spec=gs1[2],
                )
            else:
                gs1 = gridspec.GridSpec(
                    1,
                    2,
                    left=self.subplot_left,
                    right=self.subplot_right,
                    wspace=self.subplot_wspace,
                )
            gs2 = gridspec.GridSpecFromSubplotSpec(
                2, 2, hspace=self.subplot_hspace, wspace=0, subplot_spec=gs1[0]
            )
            gs3 = gridspec.GridSpecFromSubplotSpec(
                2, 2, hspace=self.subplot_hspace, wspace=0, subplot_spec=gs1[1]
            )

            # plot TE resistivity data
            self.axrte = plt.Subplot(self.fig, gs2[0, 0])
            self.fig.add_subplot(self.axrte)
            self.axrte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(te_res_arr[:, :, 0])),
                cmap=self.res_cmap,
                vmin=self.res_limits_te[0],
                vmax=self.res_limits_te[1],
            )

            # plot TE resistivity model
            self.axmrte = plt.Subplot(self.fig, gs2[0, 1])
            self.fig.add_subplot(self.axmrte)
            self.axmrte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(te_res_arr[:, :, 1])),
                cmap=self.res_cmap,
                vmin=self.res_limits_te[0],
                vmax=self.res_limits_te[1],
            )

            # plot TM resistivity data
            self.axrtm = plt.Subplot(self.fig, gs3[0, 0])
            self.fig.add_subplot(self.axrtm)
            self.axrtm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(tm_res_arr[:, :, 0])),
                cmap=self.res_cmap,
                vmin=self.res_limits_tm[0],
                vmax=self.res_limits_tm[1],
            )

            # plot TM resistivity model
            self.axmrtm = plt.Subplot(self.fig, gs3[0, 1])
            self.fig.add_subplot(self.axmrtm)
            self.axmrtm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(tm_res_arr[:, :, 1])),
                cmap=self.res_cmap,
                vmin=self.res_limits_tm[0],
                vmax=self.res_limits_tm[1],
            )

            # plot TE phase data
            self.axpte = plt.Subplot(self.fig, gs2[1, 0])
            self.fig.add_subplot(self.axpte)
            self.axpte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(te_phase_arr[:, :, 0]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_te[0],
                vmax=self.phase_limits_te[1],
            )

            # plot TE phase model
            self.axmpte = plt.Subplot(self.fig, gs2[1, 1])
            self.fig.add_subplot(self.axmpte)
            self.axmpte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(te_phase_arr[:, :, 1]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_te[0],
                vmax=self.phase_limits_te[1],
            )

            # plot TM phase data
            self.axptm = plt.Subplot(self.fig, gs3[1, 0])
            self.fig.add_subplot(self.axptm)
            self.axptm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(tm_phase_arr[:, :, 0]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_tm[0],
                vmax=self.phase_limits_tm[1],
            )

            # plot TM phase model
            self.axmptm = plt.Subplot(self.fig, gs3[1, 1])
            self.fig.add_subplot(self.axmptm)
            self.axmptm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(tm_phase_arr[:, :, 1]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_tm[0],
                vmax=self.phase_limits_tm[1],
            )

            ax_list = [
                self.axrte,
                self.axmrte,
                self.axrtm,
                self.axmrtm,
                self.axpte,
                self.axmpte,
                self.axptm,
                self.axmptm,
            ]

            if self.plot_tipper == "y":
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs4[0, 0])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_real_arr[:, :, 0]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )
                # plot real tipper  model
                self.axmtpr = plt.Subplot(self.fig, gs4[0, 1])
                self.fig.add_subplot(self.axmtpr)
                self.axmtpr.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_real_arr[:, :, 1]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )

                # plot imag tipper  data
                self.axtpi = plt.Subplot(self.fig, gs4[1, 0])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_imag_arr[:, :, 0]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )
                # plot imag tipper  model
                self.axmtpi = plt.Subplot(self.fig, gs4[1, 1])
                self.fig.add_subplot(self.axmtpi)
                self.axmtpi.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_imag_arr[:, :, 1]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )

                ax_list.append(self.axtpr)
                ax_list.append(self.axmtpr)
                ax_list.append(self.axtpi)
                ax_list.append(self.axmtpi)

            # make everthing look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.set_xlim(offset_list.min(), offset_list.max())
                if np.remainder(xx, 2.0) == 1:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cbx = mcb.make_axes(
                        ax, shrink=self.cb_shrink, pad=self.cb_pad
                    )
                if xx == 2 or xx == 6 or xx == 8 or xx == 10:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                if xx < 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    if xx == 1:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.res_cmap,
                            norm=Normalize(
                                vmin=self.res_limits_te[0],
                                vmax=self.res_limits_te[1],
                            ),
                        )
                        cb.set_ticks(
                            np.arange(
                                int(self.res_limits_te[0]),
                                int(self.res_limits_te[1]) + 1,
                            )
                        )
                        cb.set_ticklabels(log_labels_te)
                    if xx == 3:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.res_cmap,
                            norm=Normalize(
                                vmin=self.res_limits_tm[0],
                                vmax=self.res_limits_tm[1],
                            ),
                        )
                        cb.set_label(
                            "App. Res. ($\Omega \cdot$m)",
                            fontdict={
                                "size": self.font_size + 1,
                                "weight": "bold",
                            },
                        )
                        cb.set_label(
                            "Resistivity ($\Omega \cdot$m)",
                            fontdict={
                                "size": self.font_size + 1,
                                "weight": "bold",
                            },
                        )
                        cb.set_ticks(
                            np.arange(
                                int(self.res_limits_tm[0]),
                                int(self.res_limits_tm[1]) + 1,
                            )
                        )
                        cb.set_ticklabels(log_labels_tm)
                else:
                    # color bar TE phase
                    if xx == 5:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.phase_cmap,
                            norm=Normalize(
                                vmin=self.phase_limits_te[0],
                                vmax=self.phase_limits_te[1],
                            ),
                        )
                    # color bar TM phase
                    if xx == 7:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.phase_cmap,
                            norm=Normalize(
                                vmin=self.phase_limits_tm[0],
                                vmax=self.phase_limits_tm[1],
                            ),
                        )
                        cb.set_label(
                            "Phase (deg)",
                            fontdict={
                                "size": self.font_size + 1,
                                "weight": "bold",
                            },
                        )
                    # color bar tipper Imag
                    if xx == 9:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.tip_cmap,
                            norm=Normalize(
                                vmin=self.tip_limits_re[0],
                                vmax=self.tip_limits_re[1],
                            ),
                        )
                        cb.set_label(
                            "Re{T}",
                            fontdict={
                                "size": self.font_size + 1,
                                "weight": "bold",
                            },
                        )
                    if xx == 11:
                        cb = mcb.ColorbarBase(
                            cbx[0],
                            cmap=self.tip_cmap,
                            norm=Normalize(
                                vmin=self.tip_limits_im[0],
                                vmax=self.tip_limits_im[1],
                            ),
                        )
                        cb.set_label(
                            "Im{T}",
                            fontdict={
                                "size": self.font_size + 1,
                                "weight": "bold",
                            },
                        )

                ax.text(
                    xloc,
                    yloc,
                    self.label_list[xx],
                    fontdict={"size": self.font_size + 1},
                    bbox={"facecolor": "white"},
                    horizontalalignment="left",
                    verticalalignment="top",
                )
                if xx == 0 or xx == 4:
                    ax.set_ylabel(
                        "Period (s)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )
                if xx > 3:
                    ax.set_xlabel(
                        "Station",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

            plt.show()

        else:
            if self.plot_tipper == "y":
                gs1 = gridspec.GridSpec(
                    2,
                    3,
                    left=self.subplot_left,
                    right=self.subplot_right,
                    hspace=self.subplot_hspace,
                    wspace=self.subplot_wspace,
                )

            else:
                gs1 = gridspec.GridSpec(
                    2,
                    2,
                    left=self.subplot_left,
                    right=self.subplot_right,
                    hspace=self.subplot_hspace,
                    wspace=self.subplot_wspace,
                )

            # plot TE resistivity data
            self.axrte = self.fig.add_subplot(gs1[0, 0])
            self.axrte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(te_res_arr[:, :, 0])),
                cmap=self.res_cmap,
                vmin=self.res_limits_te[0],
                vmax=self.res_limits_te[1],
            )

            # plot TM resistivity data
            self.axrtm = self.fig.add_subplot(gs1[0, 1])
            self.axrtm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(np.log10(tm_res_arr[:, :, 0])),
                cmap=self.res_cmap,
                vmin=self.res_limits_tm[0],
                vmax=self.res_limits_tm[1],
            )

            # plot TE phase data
            self.axpte = self.fig.add_subplot(gs1[1, 0])
            self.axpte.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(te_phase_arr[:, :, 0]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_te[0],
                vmax=self.phase_limits_te[1],
            )

            # plot TM phase data
            self.axptm = self.fig.add_subplot(gs1[1, 1])
            self.axptm.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(tm_phase_arr[:, :, 0]),
                cmap=self.phase_cmap,
                vmin=self.phase_limits_tm[0],
                vmax=self.phase_limits_tm[1],
            )
            ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]
            if self.plot_tipper == "y":
                # plot real tipper  data
                self.axtpr = plt.Subplot(self.fig, gs1[0, 2])
                self.fig.add_subplot(self.axtpr)
                self.axtpr.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_real_arr[:, :, 0]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )
                # plot real tipper  data
                self.axtpi = plt.Subplot(self.fig, gs1[1, 2])
                self.fig.add_subplot(self.axtpi)
                self.axtpi.pcolormesh(
                    dgrid,
                    fgrid,
                    np.flipud(tip_imag_arr[:, :, 0]),
                    cmap=self.tip_cmap,
                    vmin=self.tip_limits_re[0],
                    vmax=self.tip_limits_re[1],
                )
                ax_list.append(self.axtpr)
                ax_list.append(self.axtpi)

            # make everything look tidy
            for xx, ax in enumerate(ax_list):
                ax.semilogy()
                ax.set_ylim(ylimits)
                ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
                ax.xaxis.set_ticks(offset_list, minor=True)
                ax.xaxis.set_ticklabels(slabel)
                ax.grid(True, alpha=0.25)
                ax.set_xlim(offset_list.min(), offset_list.max())
                cbx = mcb.make_axes(ax, shrink=self.cb_shrink, pad=self.cb_pad)
                if xx == 0:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.res_cmap,
                        norm=Normalize(
                            vmin=self.res_limits_te[0],
                            vmax=self.res_limits_te[1],
                        ),
                    )
                    cb.set_ticks(
                        np.arange(
                            self.res_limits_te[0], self.res_limits_te[1] + 1
                        )
                    )
                    cb.set_ticklabels(log_labels_te)
                elif xx == 1:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)

                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.res_cmap,
                        norm=Normalize(
                            vmin=self.res_limits_tm[0],
                            vmax=self.res_limits_tm[1],
                        ),
                    )
                    cb.set_label(
                        "App. Res. ($\Omega \cdot$m)",
                        fontdict={
                            "size": self.font_size + 1,
                            "weight": "bold",
                        },
                    )
                    cb.set_ticks(
                        np.arange(
                            self.res_limits_tm[0], self.res_limits_tm[1] + 1
                        )
                    )
                    cb.set_ticklabels(log_labels_tm)
                elif xx == 2:
                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.phase_cmap,
                        norm=Normalize(
                            vmin=self.phase_limits_te[0],
                            vmax=self.phase_limits_te[1],
                        ),
                    )
                    cb.set_ticks(
                        np.arange(
                            self.phase_limits_te[0],
                            self.phase_limits_te[1] + 1,
                            15,
                        )
                    )
                elif xx == 3:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.phase_cmap,
                        norm=Normalize(
                            vmin=self.phase_limits_tm[0],
                            vmax=self.phase_limits_tm[1],
                        ),
                    )
                    cb.set_label(
                        "Phase (deg)",
                        fontdict={
                            "size": self.font_size + 1,
                            "weight": "bold",
                        },
                    )
                    cb.set_ticks(
                        np.arange(
                            self.phase_limits_te[0],
                            self.phase_limits_te[1] + 1,
                            15,
                        )
                    )

                # real tipper
                elif xx == 4:
                    plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.tip_cmap,
                        norm=Normalize(
                            vmin=self.tip_limits_re[0],
                            vmax=self.tip_limits_re[1],
                        ),
                    )
                    cb.set_label(
                        "Re{T}",
                        fontdict={
                            "size": self.font_size + 1,
                            "weight": "bold",
                        },
                    )
                    # imag tipper
                elif xx == 5:
                    plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                    cb = mcb.ColorbarBase(
                        cbx[0],
                        cmap=self.tip_cmap,
                        norm=Normalize(
                            vmin=self.tip_limits_im[0],
                            vmax=self.tip_limits_im[1],
                        ),
                    )
                    cb.set_label(
                        "Im{T}",
                        fontdict={
                            "size": self.font_size + 1,
                            "weight": "bold",
                        },
                    )

                ax.text(
                    xloc,
                    yloc,
                    self.label_list[2 * xx],
                    fontdict={"size": self.font_size + 1},
                    bbox={"facecolor": "white"},
                    horizontalalignment="left",
                    verticalalignment="top",
                )
                if xx == 0 or xx == 2:
                    ax.set_ylabel(
                        "Period (s)",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )
                if xx > 1:
                    ax.set_xlabel(
                        "Station",
                        fontdict={
                            "size": self.font_size + 2,
                            "weight": "bold",
                        },
                    )

            plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotPseudoSection()
            >>> #change color of te markers to a gray-blue
            >>> p1.res_cmap = 'seismic_r'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_plot="y",
    ):
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

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        else:
            save_fn = os.path.join(
                save_fn, "OccamPseudoSection." + file_format
            )
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        if close_plot == "y":
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print("Saved figure to: " + self.fig_fn)

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
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return (
            "Plots a pseudo section of TE and TM modes for data and "
            "response if given."
        )

        # ==============================================================================


# plot misfits as a pseudo-section
# ==============================================================================
class PlotMisfitPseudoSection(object):
    """
     plot a pseudo section of the data and response if given


     Arguments:
     -------------
         **rp_list** : list of dictionaries for each station with keywords:

                 * *station* : string
                              station name

                 * *offset* : float
                              relative offset

                 * *resxy* : np.array(nf,4)
                             TE resistivity and error as row 0 and 1 respectively

                 * *resyx* : np.array(fn,4)
                             TM resistivity and error as row 0 and 1 respectively

                 * *phasexy* : np.array(nf,4)
                               TE phase and error as row 0 and 1 respectively

                 * *phaseyx* : np.array(nf,4)
                               Tm phase and error as row 0 and 1 respectively

                 * *realtip* : np.array(nf,4)
                               Real Tipper and error as row 0 and 1 respectively

                 * *imagtip* : np.array(nf,4)
                               Imaginary Tipper and error as row 0 and 1
                               respectively

                 Note: that the resistivity will be in log10 space.  Also, there
                 are 2 extra rows in the data arrays, this is to put the
                 response from the inversion.

         **period** : np.array of periods to plot that correspond to the index
                      values of each rp_list entry ie. resxy.

     ==================== ==================================================
     key words            description
     ==================== ==================================================
     axmpte               matplotlib.axes instance for TE model phase
     axmptm               matplotlib.axes instance for TM model phase
     axmrte               matplotlib.axes instance for TE model app. res
     axmrtm               matplotlib.axes instance for TM model app. res
     axpte                matplotlib.axes instance for TE data phase
     axptm                matplotlib.axes instance for TM data phase
     axrte                matplotlib.axes instance for TE data app. res.
     axrtm                matplotlib.axes instance for TM data app. res.
     cb_pad               padding between colorbar and axes
     cb_shrink            percentage to shrink the colorbar to
     fig                  matplotlib.figure instance
     fig_dpi              resolution of figure in dots per inch
     fig_num              number of figure instance
     fig_size             size of figure in inches (width, height)
     font_size            size of font in points
     label_list            list to label plots
     ml                   factor to label stations if 2 every other station
                          is labeled on the x-axis
     period               np.array of periods to plot
     phase_cmap           color map name of phase
     phase_limits_te      limits for te phase in degrees (min, max)
     phase_limits_tm      limits for tm phase in degrees (min, max)
     plot_resp            [ 'y' | 'n' ] to plot response
     plot_yn              [ 'y' | 'n' ] 'y' to plot on instantiation

     res_cmap             color map name for resistivity
     res_limits_te        limits for te resistivity in log scale (min, max)
     res_limits_tm        limits for tm resistivity in log scale (min, max)
     rp_list               list of dictionaries as made from read2Dresp
     station_id           index to get station name (min, max)
     station_list          station list got from rp_list
     subplot_bottom       subplot spacing from bottom (relative coordinates)
     subplot_hspace       vertical spacing between subplots
     subplot_left         subplot spacing from left
     subplot_right        subplot spacing from right
     subplot_top          subplot spacing from top
     subplot_wspace       horizontal spacing between subplots
     ==================== ==================================================

     =================== =======================================================
     Methods             Description
     =================== =======================================================
     plot                plots a pseudo-section of apparent resistiviy and phase
                         of data and model if given.  called on instantiation
                         if plot_yn is 'y'.
     redraw_plot         call redraw_plot to redraw the figures,
                         if one of the attributes has been changed
     save_figure         saves the matplotlib.figure instance to desired
                         location and format
     =================== =======================================================

    :Example: ::

         >>> import mtpy.modeling.occam2d as occam2d
         >>> ocd = occam2d.Occam2DData()
         >>> rfile = r"/home/Occam2D/Line1/Inv1/Test_15.resp"
         >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/DataRW.dat"
         >>> ps1 = ocd.plot2PseudoSection(resp_fn=rfile)

    """

    def __init__(self, data_fn, resp_fn, **kwargs):

        self.data_fn = data_fn
        self.resp_fn = resp_fn

        self.label_list = [
            r"$\rho_{TE}$",
            r"$\rho_{TM}$",
            "$\phi_{TE}$",
            "$\phi_{TM}$",
            "$\Re e\{T\}$",
            "$\Im m\{T\}$",
        ]

        self.phase_limits_te = kwargs.pop("phase_limits_te", (-10, 10))
        self.phase_limits_tm = kwargs.pop("phase_limits_tm", (-10, 10))
        self.res_limits_te = kwargs.pop("res_limits_te", (-2, 2))
        self.res_limits_tm = kwargs.pop("res_limits_tm", (-2, 2))
        self.tip_limits_re = kwargs.pop("tip_limits_re", (-0.2, 0.2))
        self.tip_limits_im = kwargs.pop("tip_limits_im", (-0.2, 0.2))

        self.phase_cmap = kwargs.pop("phase_cmap", "BrBG")
        self.res_cmap = kwargs.pop("res_cmap", "BrBG_r")
        self.tip_cmap = kwargs.pop("tip_cmap", "PuOr")
        self.plot_tipper = kwargs.pop("plot_tipper", "n")

        self.ml = kwargs.pop("ml", 2)
        self.station_id = kwargs.pop("station_id", [0, 4])

        self.fig_num = kwargs.pop("fig_num", 1)
        self.fig_size = kwargs.pop("fig_size", [6, 6])
        self.fig_dpi = kwargs.pop("dpi", 300)

        self.subplot_wspace = 0.0025
        self.subplot_hspace = 0.0
        self.subplot_right = 0.95
        self.subplot_left = 0.085
        self.subplot_top = 0.97
        self.subplot_bottom = 0.1

        self.font_size = kwargs.pop("font_size", 6)
        self.plot_yn = kwargs.pop("plot_yn", "y")

        self.cb_shrink = 0.7
        self.cb_pad = 0.015

        self.axrte = None
        self.axrtm = None
        self.axpte = None
        self.axptm = None
        self.axtpr = None
        self.axtpi = None

        self.misfit_te_res = None
        self.misfit_te_phase = None
        self.misfit_tm_res = None
        self.misfit_tm_phase = None
        self.misfit_tip_real = None
        self.misfit_tip_imag = None

        self.fig = None
        self._data_obj = None

        if self.plot_yn == "y":
            self.plot()

    def get_misfit(self):
        """
        compute misfit of MT response found from the model and the data.

        Need to normalize correctly
        """
        data_obj = Data()
        data_obj.read_data_file(self.data_fn)
        self._data_obj = data_obj

        resp_obj = Response()
        resp_obj.read_response_file(self.resp_fn)

        n_stations = len(data_obj.station_list)
        n_periods = len(data_obj.freq)

        self.misfit_te_res = np.zeros((n_periods, n_stations))
        self.misfit_te_phase = np.zeros((n_periods, n_stations))
        self.misfit_tm_res = np.zeros((n_periods, n_stations))
        self.misfit_tm_phase = np.zeros((n_periods, n_stations))
        self.misfit_tip_real = np.zeros((n_periods, n_stations))
        self.misfit_tip_imag = np.zeros((n_periods, n_stations))

        for rr, r_dict in zip(list(range(n_stations)), resp_obj.resp):
            self.misfit_te_res[:, rr] = r_dict["te_res"][1]
            self.misfit_tm_res[:, rr] = r_dict["tm_res"][1]
            self.misfit_te_phase[:, rr] = r_dict["te_phase"][1]
            self.misfit_tm_phase[:, rr] = r_dict["tm_phase"][1]
            self.misfit_tip_real[:, rr] = r_dict["re_tip"][1]
            self.misfit_tip_imag[:, rr] = r_dict["im_tip"][1]

        self.misfit_te_res = np.nan_to_num(self.misfit_te_res)
        self.misfit_te_phase = np.nan_to_num(self.misfit_te_phase)
        self.misfit_tm_res = np.nan_to_num(self.misfit_tm_res)
        self.misfit_tm_phase = np.nan_to_num(self.misfit_tm_phase)
        self.misfit_tip_real = np.nan_to_num(self.misfit_tip_real)
        self.misfit_tip_imag = np.nan_to_num(self.misfit_tip_imag)

    def plot(self):
        """
        plot pseudo section of data and response if given

        """

        self.get_misfit()

        ylimits = (self._data_obj.period.max(), self._data_obj.period.min())

        offset_list = np.append(
            self._data_obj.station_locations,
            self._data_obj.station_locations[-1] * 1.15,
        )

        # make a meshgrid for plotting
        # flip frequency so bottom corner is long period
        dgrid, fgrid = np.meshgrid(offset_list, self._data_obj.period[::-1])

        # make list for station labels
        ns = len(self._data_obj.station_list)
        sindex_1 = self.station_id[0]
        sindex_2 = self.station_id[1]
        slabel = [
            self._data_obj.station_list[ss][sindex_1:sindex_2]
            for ss in range(0, ns, self.ml)
        ]

        xloc = offset_list[0] + abs(offset_list[0] - offset_list[1]) / 5
        yloc = 1.10 * self._data_obj.period[1]

        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.hspace"] = self.subplot_hspace
        plt.rcParams["figure.subplot.wspace"] = self.subplot_wspace

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        if self.plot_tipper != "y":
            self.axrte = self.fig.add_subplot(2, 2, 1)
            self.axrtm = self.fig.add_subplot(2, 2, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 2, 3, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 2, 4, sharex=self.axrte)

        else:
            self.axrte = self.fig.add_subplot(2, 3, 1)
            self.axrtm = self.fig.add_subplot(2, 3, 2, sharex=self.axrte)
            self.axpte = self.fig.add_subplot(2, 3, 4, sharex=self.axrte)
            self.axptm = self.fig.add_subplot(2, 3, 5, sharex=self.axrte)
            self.axtpr = self.fig.add_subplot(2, 3, 3, sharex=self.axrte)
            self.axtpi = self.fig.add_subplot(2, 3, 6, sharex=self.axrte)

        # --> TE Resistivity
        self.axrte.pcolormesh(
            dgrid,
            fgrid,
            np.flipud(self.misfit_te_res),
            cmap=self.res_cmap,
            vmin=self.res_limits_te[0],
            vmax=self.res_limits_te[1],
        )
        # --> TM Resistivity
        self.axrtm.pcolormesh(
            dgrid,
            fgrid,
            np.flipud(self.misfit_tm_res),
            cmap=self.res_cmap,
            vmin=self.res_limits_tm[0],
            vmax=self.res_limits_tm[1],
        )
        # --> TE Phase
        self.axpte.pcolormesh(
            dgrid,
            fgrid,
            np.flipud(self.misfit_te_phase),
            cmap=self.phase_cmap,
            vmin=self.phase_limits_te[0],
            vmax=self.phase_limits_te[1],
        )
        # --> TM Phase
        self.axptm.pcolormesh(
            dgrid,
            fgrid,
            np.flipud(self.misfit_tm_phase),
            cmap=self.phase_cmap,
            vmin=self.phase_limits_tm[0],
            vmax=self.phase_limits_tm[1],
        )

        ax_list = [self.axrte, self.axrtm, self.axpte, self.axptm]

        if self.plot_tipper == "y":
            self.axtpr.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(self.misfit_tip_real),
                cmap=self.tip_cmap,
                vmin=self.tip_limits_re[0],
                vmax=self.tip_limits_re[1],
            )
            self.axtpi.pcolormesh(
                dgrid,
                fgrid,
                np.flipud(self.misfit_tip_imag),
                cmap=self.tip_cmap,
                vmin=self.tip_limits_im[0],
                vmax=self.tip_limits_im[1],
            )

            ax_list.append(self.axtpr)
            ax_list.append(self.axtpi)
            # make everthing look tidy
        for xx, ax in enumerate(ax_list):
            ax.semilogy()
            ax.set_ylim(ylimits)
            ax.xaxis.set_ticks(offset_list[np.arange(0, ns, self.ml)])
            ax.xaxis.set_ticks(offset_list, minor=True)
            ax.xaxis.set_ticklabels(slabel)
            ax.set_xlim(offset_list.min(), offset_list.max())
            cbx = mcb.make_axes(ax, shrink=self.cb_shrink, pad=self.cb_pad)

            # te res
            if xx == 0:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.res_cmap,
                    norm=Normalize(
                        vmin=self.res_limits_te[0], vmax=self.res_limits_te[1]
                    ),
                )
            # tm res
            elif xx == 1:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.res_cmap,
                    norm=Normalize(
                        vmin=self.res_limits_tm[0], vmax=self.res_limits_tm[1]
                    ),
                )
                cb.set_label(
                    "Log$_{10}$ App. Res. ($\Omega \cdot$m)",
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )
            # te phase
            elif xx == 2:
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.phase_cmap,
                    norm=Normalize(
                        vmin=self.phase_limits_te[0],
                        vmax=self.phase_limits_te[1],
                    ),
                )
            # tm phase
            elif xx == 3:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.phase_cmap,
                    norm=Normalize(
                        vmin=self.phase_limits_tm[0],
                        vmax=self.phase_limits_tm[1],
                    ),
                )
                cb.set_label(
                    "Phase (deg)",
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )

            # real tipper
            elif xx == 4:
                plt.setp(ax.xaxis.get_ticklabels(), visible=False)
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.tip_cmap,
                    norm=Normalize(
                        vmin=self.tip_limits_re[0], vmax=self.tip_limits_re[1]
                    ),
                )
                cb.set_label(
                    "Re{Tip}",
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )
                # imag tipper
            elif xx == 5:
                plt.setp(ax.yaxis.get_ticklabels(), visible=False)
                cb = mcb.ColorbarBase(
                    cbx[0],
                    cmap=self.tip_cmap,
                    norm=Normalize(
                        vmin=self.tip_limits_im[0], vmax=self.tip_limits_im[1]
                    ),
                )
                cb.set_label(
                    "Im{Tip}",
                    fontdict={"size": self.font_size + 1, "weight": "bold"},
                )

            # make label for plot
            ax.text(
                xloc,
                yloc,
                self.label_list[xx],
                fontdict={"size": self.font_size + 2},
                bbox={"facecolor": "white"},
                horizontalalignment="left",
                verticalalignment="top",
            )

            if xx == 0 or xx == 2:
                ax.set_ylabel(
                    "Period (s)",
                    fontdict={"size": self.font_size + 2, "weight": "bold"},
                )
            if xx > 1:
                ax.set_xlabel(
                    "Station",
                    fontdict={"size": self.font_size + 2, "weight": "bold"},
                )

        plt.show()

    def redraw_plot(self):
        """
        redraw plot if parameters were changed

        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.modeling.occam2d as occam2d
            >>> ocd = occam2d.Occam2DData(r"/home/occam2d/Data.dat")
            >>> p1 = ocd.plotPseudoSection()
            >>> #change color of te markers to a gray-blue
            >>> p1.res_cmap = 'seismic_r'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def save_figure(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_plot="y",
    ):
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

        if fig_dpi == None:
            fig_dpi = self.fig_dpi

        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        else:
            save_fn = os.path.join(
                save_fn, "OccamMisfitPseudoSection." + file_format
            )
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
                bbox_inches="tight",
            )

        if close_plot == "y":
            plt.clf()
            plt.close(self.fig)

        else:
            pass

        self.fig_fn = save_fn
        print("Saved figure to: " + self.fig_fn)

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
            >>> ps1 = ocd.plotPseudoSection()
            >>> [ax.grid(True, which='major') for ax in [ps1.axrte,ps1.axtep]]
            >>> ps1.update_plot()

        """

        self.fig.canvas.draw()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return (
            "Plots a pseudo section of TE and TM modes for data and "
            "response if given."
        )


class OccamPointPicker(object):
    """
    This class helps the user interactively pick points to mask and add
    error bars.

    Useage:
    -------
    To mask just a single point right click over the point and a gray point
    will appear indicating it has been masked

    To mask both the apparent resistivity and phase left click over the point.
    Gray points will appear over both the apparent resistivity and phase.
    Sometimes the points don't exactly matchup, haven't quite worked that bug
    out yet, but not to worry it picks out the correct points

    To add error bars to a point click the middle or scroll bar button.  This
    only adds error bars to the point and does not reduce them so start out
    with reasonable errorbars.  You can change the increment that the error
    bars are increased with res_err_inc and phase_err_inc.

    .. note:: There is a bug when only plotting TE or TM that you cannot mask
              points in the phase.  I'm not sure where it comes from, but
              works with all modes.  So my suggestion is to make a data file
              with all modes, mask data points and then rewrite that data file
              if you want to use just one of the modes.  That's the work
              around for the moment.



    Arguments:
    ----------
        **ax_list** : list of the resistivity and phase axis that have been
                    plotted as [axr_te,axr_tm,axp_te,axp_tm]

        **line_list** : list of lines used to plot the responses, not the error
                      bars as [res_te,res_tm,phase_te,phase_tm]

        **err_list** : list of the errorcaps and errorbar lines as
                   [[cap1,cap2,bar],...]

        **res_err_inc** : increment to increase the errorbars for resistivity.
                        put .20 for 20 percent change. *Default* is .05

        **phase_err_inc** : increment to increase the errorbars for the phase
                          put .10 for 10 percent change. *Defualt* is .02

        **marker** : marker type for masked points.  See matplotlib.pyplot.plot
                    for options of markers.  *Default* is h for hexagon.

    Attributes:
    -----------

        **ax_list** : axes list used to plot the data

        **line_list** : line list used to plot the data

        **err_list** : error list used to plot the data

        **data** : list of data points that were not masked for each plot.

        **fdict** : dictionary of frequency arrays for each plot and data set.

        **fndict** : dictionary of figure numbers to corresponed with data.

        **cid_list** : list of event ids.

        **res_err_inc** : increment to increase resistivity error bars

        **phase_inc** : increment to increase phase error bars

        **marker** : marker of masked points

        **fig_num** : figure numbers

        **data_list** : list of lines to write into the occam2d data file.

    :Example: ::

        >>> ocd = occam2d.Occam2DData()
        >>> ocd.data_fn = r"/home/Occam2D/Line1/Inv1/Data.dat"
        >>> ocd.plotMaskPoints()
    """

    def __init__(
        self,
        ax_list,
        line_list,
        err_list,
        res_err_inc=0.05,
        phase_err_inc=0.02,
        marker="h",
    ):

        # give the class some attributes
        self.ax_list = ax_list
        self.line_list = line_list
        self.err_list = err_list
        self.data = []
        self.error = []
        self.fdict = []
        self.fndict = {}
        # see if just one figure is plotted or multiple figures are plotted
        self.ax = ax_list[0][0]
        self.line = line_list[0][0]
        self.cidlist = []
        self.ax_num = None
        self.res_index = None
        self.phase_index = None
        self.fig_num = 0
        for nn in range(len(ax_list)):
            self.data.append([])
            self.error.append([])
            self.fdict.append([])

            # get data from lines and make a dictionary of frequency points for
            # easy indexing
            # line_find = False
            for ii, line in enumerate(line_list[nn]):
                try:
                    self.data[nn].append(line.get_data()[1])
                    self.fdict[nn].append(
                        dict(
                            [
                                ("{0:.5g}".format(kk), ff)
                                for ff, kk in enumerate(line.get_data()[0])
                            ]
                        )
                    )
                    self.fndict["{0}".format(line.figure.number)] = nn

                    # set some events
                    # if line_find == False:
                    cid1 = line.figure.canvas.mpl_connect("pick_event", self)
                    cid2 = line.figure.canvas.mpl_connect(
                        "axes_enter_event", self.inAxes
                    )
                    cid3 = line.figure.canvas.mpl_connect(
                        "key_press_event", self.on_close
                    )
                    cid4 = line.figure.canvas.mpl_connect(
                        "figure_enter_event", self.inFigure
                    )
                    self.cidlist.append([cid1, cid2, cid3, cid4])

                    # set the figure number
                    self.fig_num = self.line.figure.number

                    # line_find = True

                except AttributeError:
                    self.data[nn].append([])
                    self.fdict[nn].append([])

            # read in the error in a useful way so that it can be translated to
            # the data file.  Make the error into an array
            for ee, err in enumerate(err_list[nn]):
                try:
                    errpath = err[2].get_paths()
                    errarr = np.zeros(len(list(self.fdict[nn][ee].keys())))
                    for ff, epath in enumerate(errpath):
                        errv = epath.vertices
                        errarr[ff] = abs(errv[0, 1] - self.data[nn][ee][ff])
                    self.error[nn].append(errarr)
                except AttributeError:
                    self.error[nn].append([])

        # set the error bar increment values
        self.res_err_inc = res_err_inc
        self.phase_err_inc = phase_err_inc

        # set the marker
        self.marker = marker

        # make a list of occam2d lines to write later
        self.data_list = []

    def __call__(self, event):
        """
        When the function is called the mouse events will be recorder for
        picking points to mask or change error bars.  The axes is redrawn with
        a gray marker to indicate a masked point and/or increased size in
        errorbars.

        Arguments:
        ----------
            **event** : type mouse_click_event

        Useage:
        -------

            **Left mouse button** will mask both resistivity and phase point

            **Right mouse button** will mask just the point selected

            **Middle mouse button** will increase the error bars

            **q** will close the figure.
        """
        self.event = event
        # make a new point that is an PickEvent type
        npoint = event.artist
        # if the right button is clicked mask the point
        if event.mouseevent.button == 3:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fig_num][self.ax_num]["{0:.5g}".format(xd[0])]

            # change the data to be a zero
            self.data[self.fig_num][self.ax_num][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(
                xd,
                yd,
                ls="None",
                color=(0.7, 0.7, 0.7),
                marker=self.marker,
                ms=4,
            )

            # redraw the canvas
            self.ax.figure.canvas.draw()

        # if the left button is clicked change both resistivity and phase
        # points
        elif event.mouseevent.button == 1:
            # get the point that was clicked on
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]

            # set the x index from the frequency dictionary
            ll = self.fdict[self.fig_num][self.ax_num]["{0:.5g}".format(xd[0])]

            # set the data point to zero
            # print self.data[self.fig_num][self.ax_num][ll]
            self.data[self.fig_num][self.ax_num][ll] = 0

            # reset the point to be a gray x
            self.ax.plot(
                xd,
                yd,
                ls="None",
                color=(0.7, 0.7, 0.7),
                marker=self.marker,
                ms=4,
            )

            self.ax.figure.canvas.draw()

            # check to make sure there is a corresponding res/phase point
            try:
                kk = (self.ax_num + 2) % 4
                print(kk)
                # get the corresponding y-value
                yd2 = self.data[self.fig_num][kk][ll]

                # set that data point to 0 as well
                self.data[self.fig_num][kk][ll] = 0

                # make that data point a gray x
                self.ax_list[self.fig_num][kk].plot(
                    xd,
                    yd2,
                    ls="None",
                    color=(0.7, 0.7, 0.7),
                    marker=self.marker,
                    ms=4,
                )
                # redraw the canvas
                self.ax.figure.canvas.draw()
            except KeyError:
                print("Axis does not contain res/phase point")

        # if click the scroll button or middle button change increase the
        # errorbars by the given amount
        elif event.mouseevent.button == 2:
            ii = event.ind
            xd = npoint.get_xdata()[ii]
            yd = npoint.get_ydata()[ii]
            jj = self.ax_num

            # get x index
            ll = self.fdict[self.fig_num][jj]["{0:.5g}".format(xd[0])]

            # make error bar array
            eb = self.err_list[self.fig_num][jj][2].get_paths()[ll].vertices

            # make ecap array
            ecapl = self.err_list[self.fig_num][jj][0].get_data()[1][ll]
            ecapu = self.err_list[self.fig_num][jj][1].get_data()[1][ll]

            # change apparent resistivity error
            if jj == 0 or jj == 1:
                nebu = eb[0, 1] - self.res_err_inc * eb[0, 1]
                nebl = eb[1, 1] + self.res_err_inc * eb[1, 1]
                ecapl = ecapl - self.res_err_inc * ecapl
                ecapu = ecapu + self.res_err_inc * ecapu

            # change phase error
            elif jj == 2 or jj == 3:
                nebu = eb[0, 1] - eb[0, 1] * self.phase_err_inc
                nebl = eb[1, 1] + eb[1, 1] * self.phase_err_inc
                ecapl = ecapl - ecapl * self.phase_err_inc
                ecapu = ecapu + ecapu * self.phase_err_inc

            # put the new error into the error array
            self.error[self.fig_num][jj][ll] = abs(
                nebu - self.data[self.fig_num][jj][ll]
            )

            # set the new error bar values
            eb[0, 1] = nebu
            eb[1, 1] = nebl

            # reset the error bars and caps
            ncapl = self.err_list[self.fig_num][jj][0].get_data()
            ncapu = self.err_list[self.fig_num][jj][1].get_data()
            ncapl[1][ll] = ecapl
            ncapu[1][ll] = ecapu

            # set the values
            self.err_list[self.fig_num][jj][0].set_data(ncapl)
            self.err_list[self.fig_num][jj][1].set_data(ncapu)
            self.err_list[self.fig_num][jj][2].get_paths()[ll].vertices = eb

            # redraw the canvas
            self.ax.figure.canvas.draw()

    # get the axis number that the mouse is in and change to that axis
    def inAxes(self, event):
        """
        gets the axes that the mouse is currently in.

        Arguments:
        ---------
            **event**: is a type axes_enter_event

        Returns:
        --------

            **OccamPointPicker.jj** : index of resistivity axes for ax_list

            **OccamPointPicker.kk** : index of phase axes for ax_list

        """

        self.event2 = event
        self.ax = event.inaxes
        for jj, axj in enumerate(self.ax_list):
            for ll, axl in enumerate(axj):
                if self.ax == axl:
                    self.ax_num = ll
        self.line = self.line_list[self.fig_num][self.ax_num]

    # get the figure number that the mouse is in
    def inFigure(self, event):
        """
        gets the figure number that the mouse is in

        Arguments:
        ----------
            **event** : figure_enter_event

        Returns:
        --------
            **OccamPointPicker.fig_num** : figure number that corresponds to the
                                          index in the ax_list, datalist, errorlist
                                          and line_list.

        """
        self.event3 = event
        self.fig_num = self.fndict["{0}".format(event.canvas.figure.number)]

    # type the q key to quit the figure and disconnect event handling
    def on_close(self, event):
        """
        close the figure with a 'q' key event and disconnect the event ids

        Arguments:
        ----------
            **event** : key_press_event

        Returns:
        --------
            print statement saying the figure is closed
        """
        self.event3 = event
        if self.event3.key == "q":
            for cid in self.cidlist[self.fig_num]:
                event.canvas.mpl_disconnect(cid)
            plt.close(event.canvas.figure)
            print("Closed figure ", self.fig_num)


class Run:
    """
    Run Occam2D by system call.

    Future plan: implement Occam in Python and call it from here directly.
    """


class Mask(Data):
    """
    Allow masking of points from data file (effectively commenting them out,
    so the process is reversable). Inheriting from Data class.
    """


class OccamInputError(Exception):
    pass


# ======================================
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print(
            (
                "\n please provide path to edi files\n USAGE:  %s path2edifiles"
                % sys.argv[0]
            )
        )
        sys.exit(1)
    else:
        edi_dir = sys.argv[1]

    # stations = ['151{0:02}A'.format(s) for s in xrange(24, 31)]
    # print (stations)
    # pr = Profile(edi_path=edi_dir, station_list=stations)
    # OR pr = Profile(edi_path=edi_dir,
    # station_list=['16125A','16124A','16123A','16127A','16126A', '16122A'])

    pr = Profile(edi_path=edi_dir)

    # pr.geoelectric_strike = 45

    pr.generate_profile()

    print((pr.profile_angle))

    # set station labels to only be from 1st to 4th index of station name
    pr.plot_profile(station_id=[0, 4])
