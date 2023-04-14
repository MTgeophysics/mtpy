# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:28:24 2013

@author: jpeacock-pr
"""

# ==============================================================================

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from mtpy.imaging.mtplot_tools import PlotBase
from mtpy.core import Tipper

# ==============================================================================


class PlotStrike(PlotBase):
    """
    PlotStrike will plot the strike estimated from the invariants, phase tensor
    and the tipper in either a rose diagram of xy plot


    plots the strike angle as determined by phase tensor azimuth (Caldwell et
    al. [2004]) and invariants of the impedance tensor (Weaver et al. [2003]).

    The data is split into decades where the histogram for each is plotted in
    the form of a rose diagram with a range of 0 to 180 degrees.
    Where 0 is North and 90 is East.   The median angle of the period band is
    set in polar diagram.  The top row is the strike estimated from
    the invariants of the impedance tensor.  The bottom row is the azimuth
    estimated from the phase tensor.  If tipper is 'y' then the 3rd row is the
    strike determined from the tipper, which is orthogonal to the induction
    arrow direction.

    Arguments
    ----------


        :param fs: font size for labels of plotting. *Default* is 10

        :param rot_z: angle of rotation clockwise positive. *Default* is 0

        :param period_tolerance: float
                   Tolerance level to match periods from different edi files.
                   *Default* is 0.05

        :param text_pad: padding of the angle label at the bottom of each
                         polar diagram.  *Default* is 1.65

        :param text_size: font size


        :param plot_range: [ 'data' | (period_min,period_max) ]
                    period range to estimate the strike angle. Options are:
                        * *'data'* for estimating the strike for all periods
                            in the data.
                        * (pmin,pmax) for period min and period max, input as
                          (log10(pmin),log10(pmax))

        :param plot_type: [ 1 | 2 ]
                        - *1* to plot individual decades in one plot
                        - *2* to plot all period ranges into one polar diagram
                              for each strike angle estimation

        :param plot_tipper: [ True | False ]
                          - True to plot the tipper strike
                          - False to not plot tipper strike

        :param pt_error_floor: Maximum error in degrees that is allowed to
                               estimate strike. *Default* is None allowing all
                               estimates to be used.

        :param fold: [ True | False ]
                    * True to plot only from 0 to 180
                    * False to plot from 0 to 360

        :param plot_orthogonal: [ True | False]
                                * True to plot the orthogonal strike directions
                                * False to not

        :param color: [ True | False ]
                      * True to plot shade colors
                      * False to plot all in one color

        :param color_inv: color of invariants plots

        :param color_pt: color of phase tensor plots

        :param color_tip: color of tipper plots

        :param ring_spacing: spacing of rings in polar plots

        :param ring_limits: (min count, max count) set each plot have these
                            limits

        :param plot_orientation: [ 'h' | 'v' ] horizontal or vertical plots

    :Example: ::

        >>> import glob
        >>> import mtpy.imaging.mtplot as mtplot
        >>> edi_dir = r"/home/EDIFiles"
        >>> edi_list = glob.glob("{0}\*.edi".format(edi_dir)
        >>> #---plot rose plots in decades
        >>> strike = mtplot.plot_strike(fn_list=edilist, plot_type=1)
        >>> #---Turn on Tipper
        >>> strike.plot_tipper = True
        >>> #---Plot only main directions
        >>> strike.plot_orthogonal = False
        >>> # Redraw plot
        >>> strike.redraw_plot()
        >>> # plot only from 0-180
        >>> strike.fold = True
        >>> strike.redraw_plot()
        >>> #---save the plot---
        >>> strike.save_plot(r"/home/Figures")

    """

    def __init__(self, mt_data, **kwargs):

        self._rotation_angle = 0
        self.mt_data = mt_data

        super().__init__(**kwargs)

        # ------Set attributes of the class-----------------
        # --> set plot properties
        self.plot_range = "data"
        self.plot_orientation = "h"
        self.plot_type = 2

        self.period_tolerance = 0.05
        self.pt_error_floor = None
        self.fold = False
        self.bin_width = 5
        self.color = True
        self.color_inv = (0.7, 0, 0.2)
        self.color_pt = (0.2, 0, 0.7)
        self.color_tip = (0.2, 0.65, 0.2)
        self.ring_spacing = 10
        self.ring_limits = None
        self.plot_orthogonal = False
        self.plot_pt = True
        self.plot_tipper = True
        self.plot_invariant = True
        self.print_stats = False

        self.polar_limits = (np.deg2rad(-180), np.deg2rad(180))

        # make a dictionary for plotting titles
        self.title_dict = {
            -5: "10$^{-5}$ - 10$^{-4}$ s",
            -4: "10$^{-4}$ - 10$^{-3}$ s",
            -3: "10$^{-3}$ - 10$^{-2}$ s",
            -2: "10$^{-2}$ - 10$^{-1}$ s",
            -1: "10$^{-1}$ - 10$^{0}$ s",
            0: "10$^{0}$ - 10$^{1}$ s",
            1: "10$^{1}$ - 10$^{2}$ s",
            2: "10$^{2}$ - 10$^{3}$ s",
            3: "10$^{3}$ - 10$^{4}$ s",
            4: "10$^{4}$ - 10$^{5}$ s",
            5: "10$^{5}$ - 10$^{6}$ s",
            6: "10$^{6}$ - 10$^{7}$ s",
        }

        self.strike_df = None
        self.subplot_hspace = 0.3
        self.subplot_wspace = 0.2
        self.font_weight = "normal"
        self.text_y_pad = 1.5

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    # ---need to rotate data on setting rotz
    @property
    def rotation_angle(self):
        return self._rotation_angle

    @rotation_angle.setter
    def rotation_angle(self, value):
        """
        only a single value is allowed
        """
        for mt in self.mt_data.values():
            mt.rotation_angle = value
        self._rotation_angle = value

        self.make_strike_df()

    def make_strike_df(self):
        """
        make strike array

        .. note:: Polar plots assume the azimuth is an angle measured
                counterclockwise positive from x = 0.  Therefore all angles
                are calculated as 90 - angle to make them conform to the
                polar plot convention.

        """

        entries = []

        for mt in self.mt_data.values():
            # -----------get strike angle from invariants----------------------
            if mt.has_impedance():
                zinv = mt.Z.invariants

                # subtract 90 because polar plot assumes 0 is on the x an 90 is
                # on the y
                zs = 90 - zinv.strike
                # make a dictionary of strikes with keys as period
                for period, plot_strike, strike in zip(
                    mt.period, zs, zinv.strike
                ):
                    entry = {
                        "estimate": "invariant",
                        "period": period,
                        "plot_strike": plot_strike,
                        "measured_strike": strike,
                    }
                    entries.append(entry)

                # ------------get strike from phase tensor strike angle------------
                # subtract 90 because polar plot assumes 0 is on the x an 90 is
                # on the y
                pt = mt.pt
                az = 90 - pt.azimuth
                az_err = pt.azimuth_error
                az[pt.phimax == 0] = np.nan

                # put an error max on the estimation of strike angle
                if self.pt_error_floor:
                    az[np.where(az_err > self.pt_error_floor)] = 0.0

                # make a dictionary of strikes with keys as period
                for period, plot_strike, strike in zip(
                    mt.period, az, pt.azimuth
                ):
                    entry = {
                        "estimate": "pt",
                        "period": period,
                        "plot_strike": plot_strike,
                        "measured_strike": strike,
                    }
                    entries.append(entry)

            # -----------get tipper strike------------------------------------
            if mt.has_tipper():
                tip = mt.Tipper
                if isinstance(tip, Tipper):
                    if tip.tipper is None:
                        tip = Tipper(
                            np.zeros((len(mt.period), 1, 2), dtype="complex"),
                            frequency=[1],
                        )

                # # subtract 90 because polar plot assumes 0 is on the x an 90 is
                # on the y
                tipr = 270 - tip.angle_real

                tipr[np.where(abs(tipr) == 180.0)] = np.nan
                tipr[np.where(abs(tipr) == 0)] = np.nan

                # make a dictionary of strikes with keys as period
                for period, plot_strike, strike in zip(
                    mt.period, tipr, tip.angle_real
                ):
                    entry = {
                        "estimate": "tipper",
                        "period": period,
                        "plot_strike": plot_strike,
                        "measured_strike": strike,
                    }
                    entries.append(entry)

        self.strike_df = pd.DataFrame(entries)

    def get_mean(self, estimate_df):
        """
        get mean value
        """
        s_mean = estimate_df.measured_strike.mean(skipna=True)
        s_mean %= 360

        return s_mean

    def get_median(self, estimate_df):
        """
        get median value
        """
        s_median = estimate_df.measured_strike.median(skipna=True)
        s_median %= 360

        return s_median

    def get_mode(self, estimate_df):
        """
        get mode from a historgram
        """

        bins = np.linspace(-360, 360, 146)

        binned = pd.cut(estimate_df["measured_strike"], bins).value_counts()
        s_mode = binned.index[np.argmax(binned)].mid
        s_mode %= 360

        return s_mode

    def get_estimate(self, estimate, period_range=None):
        if period_range is None:
            return self.strike_df.loc[self.strike_df.estimate == estimate]

        else:
            estimate_df = self.strike_df.loc[
                self.strike_df.estimate == estimate
            ]

            return estimate_df.loc[
                (self.strike_df.period >= period_range[0])
                & ((self.strike_df.period < period_range[1]))
            ]

    def get_stats(self, estimate, period_range=None):
        """
        print stats nicely
        """
        estimate_df = self.get_estimate(estimate, period_range)
        # print out the statistics of the strike angles
        s_mean = self.get_mean(estimate_df)
        s_median = self.get_median(estimate_df)
        s_mode = self.get_mode(estimate_df)

        msg = f"Strike statistics for {estimate} "
        if period_range is None:
            msg += "in all periods "
        else:
            msg += f"period range {period_range[0]:.3g} to {period_range[1]:.3g} (s) "
        msg += f"median={s_median:.1f} mode={s_mode:.1f} mean={s_mean:.1f}"
        self.logger.debug(msg)
        if self.print_stats:
            print(msg)

        return s_median, s_mode, s_mean

    def get_plot_array(self, estimate, period_range=None):
        """
        get a plot array that has the min and max angles
        """
        estimate_df = self.get_estimate(estimate, period_range)
        # array goes from
        st_array = estimate_df.plot_strike.to_numpy().flatten()
        st_array = st_array[np.isfinite(st_array)]
        plot_array = np.hstack([st_array, (st_array + 180) % 360])

        if self.plot_orthogonal:
            plot_array = np.hstack([plot_array, (plot_array + 90) % 360])
        if self.fold:
            plot_array %= 180

        return plot_array

    def _get_histogram_range(self):
        """
        get histogram range based on fold
        """
        if self.fold == True:
            return (0, 180)
        elif self.fold == False:
            return (0, 360)

    def _get_bin_range(self):
        """
        get the bin range
        """
        ### get the range in periods to plot
        if self.plot_range == "data":
            return np.arange(
                np.floor(np.log10(self.strike_df.period.min())),
                np.ceil(np.log10(self.strike_df.period.max())),
                1,
            )
        else:
            return np.arange(
                np.floor(self.plot_range[0]), np.ceil(self.plot_range[1]), 1
            )

    def _get_n_subplots(self):
        """
        get number of subplots
        """

        n_subplots = 0
        if self.plot_pt:
            n_subplots += 1
        if self.plot_invariant:
            n_subplots += 1
        if self.plot_tipper:
            n_subplots += 1

        return n_subplots

    def _get_subplots(self, index=None):

        ax_inv = None
        ax_pt = None
        ax_tip = None
        if self.plot_type == 1:
            n_subplots = self._get_n_subplots()
            nb = len(self._get_bin_range())

            if "h" in self.plot_orientation:
                if self.plot_invariant:
                    ax_inv = self.fig.add_subplot(
                        n_subplots, nb, index, polar=True
                    )
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(
                            n_subplots, nb, index + nb, polar=True
                        )
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                n_subplots, nb, index + 2 * nb, polar=True
                            )
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                n_subplots, nb, index + nb, polar=True
                            )
                else:
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(
                            n_subplots, nb, index, polar=True
                        )
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                n_subplots, nb, index + nb, polar=True
                            )
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                n_subplots, nb, index, polar=True
                            )

            elif "v" in self.plot_orientation:
                if self.plot_invariant:
                    ax_inv = self.fig.add_subplot(
                        nb,
                        n_subplots,
                        index * n_subplots - (n_subplots - 1),
                        polar=True,
                    )
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(
                            nb,
                            n_subplots,
                            index * n_subplots - (n_subplots - 2),
                            polar=True,
                        )
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                nb,
                                n_subplots,
                                index * n_subplots(n_subplots - 3),
                                polar=True,
                            )
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                nb,
                                n_subplots,
                                index * n_subplots(n_subplots - 2),
                                polar=True,
                            )
                else:
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(
                            nb,
                            n_subplots,
                            index * n_subplots - (n_subplots - 1),
                            polar=True,
                        )
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                nb,
                                n_subplots,
                                index * n_subplots(n_subplots - 2),
                                polar=True,
                            )
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(
                                nb,
                                n_subplots,
                                index * n_subplots(n_subplots - 1),
                                polar=True,
                            )

        elif self.plot_type == 2:

            nb = self._get_n_subplots()

            if "h" in self.plot_orientation:

                if self.plot_invariant:
                    ax_inv = self.fig.add_subplot(1, nb, 1, polar=True)
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(1, nb, 2, polar=True)
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(1, nb, 3, polar=True)
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(1, nb, 2, polar=True)
                else:
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(1, nb, 1, polar=True)
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(1, nb, 2, polar=True)
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(1, nb, 1, polar=True)
            elif "v" in self.plot_orientation:

                if self.plot_invariant:
                    ax_inv = self.fig.add_subplot(nb, 1, 1, polar=True)
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(nb, 1, 2, polar=True)
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(nb, 1, 3, polar=True)
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(nb, 1, 2, polar=True)
                else:
                    if self.plot_pt:
                        ax_pt = self.fig.add_subplot(nb, 1, 1, polar=True)
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(nb, 1, 2, polar=True)
                    else:
                        if self.plot_tipper:
                            ax_tip = self.fig.add_subplot(nb, 1, 1, polar=True)

        return ax_inv, ax_pt, ax_tip

    def _plot_bars(self, ax, strike, hist_range):
        """
        plot rose diagram of given strike array
        """
        hist = np.histogram(
            strike[np.nonzero(strike)].flatten(),
            bins=int(360 / self.bin_width),
            range=hist_range,
        )

        # make a bar graph with each bar being width of bw degrees
        return hist, ax.bar(
            np.deg2rad(hist[1][:-1]),
            hist[0],
            width=np.deg2rad(self.bin_width),
        )

    def _set_bar_color(self, bars, hist, estimate):
        """
        set the bar colors according to estimate

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if estimate == "tipper":
            # set color of the bars according to the number in that bin
            # tipper goes from dark blue (low) to light blue (high)
            for cc, bar in enumerate(bars):
                if self.color:
                    try:
                        fc = float(hist[0][cc]) / hist[0].max() * 0.9
                        bar.set_facecolor(
                            (self.color_tip[0], 1 - fc, self.color_tip[-1])
                        )
                    except ZeroDivisionError:
                        pass
                else:
                    bar.set_facecolor(self.color_tip)

        elif estimate == "invariant":
            # set the color of the bars according to the number in that bin
            # invariants go from purple (low) to red (high)
            for cc, bar in enumerate(bars):
                if self.color:
                    try:
                        fc = float(hist[0][cc]) / hist[0].max() * 0.8
                        bar.set_facecolor((0.75, 1 - fc, 0))
                    except ZeroDivisionError:
                        pass
                else:
                    bar.set_facecolor(self.color_inv)

        elif estimate == "pt":
            # pt goes from green (low) to orange (high)
            for cc, bar in enumerate(bars):
                if self.color:
                    try:
                        fc = float(hist[0][cc]) / hist[0].max() * 0.8
                        bar.set_facecolor((1 - fc, 0, 1 - fc))
                    except ZeroDivisionError:
                        pass
                else:
                    bar.set_facecolor(self.color_pt)

    def _set_ax_label(self, ax, label, box_color):
        """

        :param ax: DESCRIPTION
        :type ax: TYPE
        :param label: DESCRIPTION
        :type label: TYPE
        :param box_color: DESCRIPTION
        :type box_color: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if "h" in self.plot_orientation:
            ax.set_ylabel(
                label,
                fontdict=self.font_dict,
                labelpad=self.font_size,
                bbox={
                    "facecolor": box_color,
                    "alpha": 0.25,
                },
            )
            # --> set the title offset
            ax.titleOffsetTrans._t = (0, 0.1)
        elif "v" in self.plot_orientation:
            ax.set_ylabel(
                label,
                fontdict=self.font_dict,
                bbox={"facecolor": "white", "alpha": 0.25},
                rotation=0,
                labelpad=50,
            )
            ax.yaxis.set_label_position("right")

    def _plot_per_period(self):
        """
        plot per period range
        """

        # plot specs
        self.subplot_hspace = 0.3
        self.subplot_wspace = 0.3
        self._set_subplot_params()

        hist_range = self._get_histogram_range()
        bin_range = self._get_bin_range()

        self.fig = plt.figure(self.fig_num, dpi=self.fig_dpi)
        plt.clf()

        for jj, bb in enumerate(bin_range, 1):
            # make subplots for invariants and phase tensor azimuths
            # dependent on vertical or horizontal orientation
            ax_inv, ax_pt, ax_tip = self._get_subplots(jj)

            # make a list of indicies for each decades
            period_range = [10**bb, 10 ** (bb + 1)]

            ax_dict = {}
            if self.plot_invariant:
                ### plot invariants azimuths
                inv = self.get_plot_array("invariant", period_range)
                inv_hist, inv_bar = self._plot_bars(ax_inv, inv, hist_range)
                self._set_bar_color(inv_bar, inv_hist, "invariant")
                ax_dict["invariant"] = ax_inv

                # set plot labels
                if jj == 1:
                    if "h" in self.plot_orientation:
                        self._set_ax_label(
                            ax_inv, "Strike (Z)", self.color_inv
                        )
                    elif "v" in self.plot_orientation:
                        self._set_ax_label(ax_inv, self.title_dict[bb], None)

            if self.plot_pt:
                ### plot phase tensor azimuths
                pt = self.get_plot_array("pt", period_range)
                pt_hist, pt_bar = self._plot_bars(ax_pt, pt, hist_range)
                self._set_bar_color(pt_bar, pt_hist, "pt")
                ax_dict["pt"] = ax_pt

                # set plot labels
                if jj == 1:
                    if "h" in self.plot_orientation:
                        self._set_ax_label(ax_pt, "PT Azimuth", self.color_pt)
                    elif "v" in self.plot_orientation:
                        self._set_ax_label(ax_pt, self.title_dict[bb], None)

            if self.plot_tipper:
                tr = self.get_plot_array("tipper", period_range)
                tr_hist, tr_bar = self._plot_bars(ax_tip, tr, hist_range)
                self._set_bar_color(tr_bar, tr_hist, "tipper")
                ax_dict["tipper"] = ax_tip
                # set plot labels
                if jj == 1:
                    if "h" in self.plot_orientation:
                        self._set_ax_label(ax_tip, "Tipper", self.color_tip)
                    elif "v" in self.plot_orientation:
                        self._set_ax_label(ax_tip, self.title_dict[bb], None)

            # make axis look correct with N to the top at 90.
            for estimate, axh in ax_dict.items():
                # set multiple locator to be every 15 degrees
                axh.xaxis.set_major_locator(MultipleLocator(30 * np.pi / 180))

                ### set labels on the correct axis
                # axh.xaxis.set_ticklabels(
                #     ["", "", "", "", "", "", "", "", "", "", "", ""]
                # )
                ### set y limits if asked for
                if self.ring_limits is not None:
                    axh.set_ylim(self.ring_limits)
                axh.yaxis.set_major_locator(MultipleLocator(self.ring_spacing))

                # make a light grid
                axh.grid(alpha=0.25, zorder=0)
                axh.set_xlim(self.polar_limits)

                # properties for the invariants
                # limits need to be rotate 90 counter clockwise because
                # we already rotated by 90 degrees so the range is
                # from -90 to 270 with -90 being east

                # label the plot with the mode value of strike
                # need to subtract 90 again because the histogram is
                # for ploting 0 east, 90 north measuring
                # counter-clockwise
                st_median, st_mode, st_mean = self.get_stats(
                    estimate, period_range
                )
                if estimate == "invariant":
                    box_color = self.color_inv
                elif estimate == "pt":
                    box_color = self.color_pt
                elif estimate == "tipper":
                    box_color = self.color_tip

                ### place the estimated strike
                axh.text(
                    -np.pi / 2,
                    axh.get_ylim()[1] * self.text_y_pad,
                    f"{st_mode:.1f}$^o$",
                    horizontalalignment="center",
                    verticalalignment="baseline",
                    fontdict={"size": self.text_size},
                    bbox={"facecolor": box_color, "alpha": 0.25},
                )

                # --> set title of subplot
                if "h" in self.plot_orientation:
                    axh.set_title(
                        self.title_dict[bb],
                        fontdict=self.font_dict,
                        bbox={"facecolor": "white", "alpha": 0.25},
                    )

                    # --> set the title offset
                    axh.titleOffsetTrans._t = (0, 0.1)
                elif "v" in self.plot_orientation:
                    axh.set_ylabel(
                        self.title_dict[bb],
                        fontdict=self.font_dict,
                        bbox={"facecolor": "white", "alpha": 0.25},
                        rotation=0,
                        labelpad=50,
                    )
                    axh.yaxis.set_label_position("right")

                plt.setp(axh.yaxis.get_ticklabels(), visible=False)
                plt.setp(axh.xaxis.get_ticklabels(), visible=False)
        self.logger.info(
            "Note: North is assumed to be 0 and the strike angle is measured"
            + "clockwise positive."
        )

        plt.show()

    def _plot_all_periods(self):
        """
        plot all periods into one rose diagram

        :return: DESCRIPTION
        :rtype: TYPE

        """
        hist_range = self._get_histogram_range()

        # plot specs
        self._set_subplot_params()

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        ax_inv, ax_pt, ax_tip = self._get_subplots()

        ax_dict = {}
        if self.plot_invariant:
            ### plot invariants azimuths
            inv = self.get_plot_array("invariant")
            inv_hist, inv_bar = self._plot_bars(ax_inv, inv, hist_range)
            self._set_bar_color(inv_bar, inv_hist, "invariant")
            ax_dict["invariant"] = ax_inv

        if self.plot_pt:
            ### plot phase tensor azimuths
            pt = self.get_plot_array("pt")
            pt_hist, pt_bar = self._plot_bars(ax_pt, pt, hist_range)
            self._set_bar_color(pt_bar, pt_hist, "pt")
            ax_dict["pt"] = ax_pt

        if self.plot_tipper:
            tr = self.get_plot_array("tipper")
            tr_hist, tr_bar = self._plot_bars(ax_tip, tr, hist_range)
            self._set_bar_color(tr_bar, tr_hist, "tipper")
            ax_dict["tipper"] = ax_tip
            # set plot labels

        # make axis look correct with N to the top at 90.
        for estimate, axh in ax_dict.items():
            # set major ticks to be every 30 degrees
            axh.xaxis.set_major_locator(MultipleLocator(2 * np.pi / 12))

            # set a light grid
            axh.grid(alpha=0.25)  # works in 2.0.2 not 2.1.0

            st_median, st_mode, st_mean = self.get_stats(estimate)
            if estimate == "invariant":
                box_color = self.color_inv
            elif estimate == "pt":
                box_color = self.color_pt
            elif estimate == "tipper":
                box_color = self.color_tip

            axh.text(
                170 * np.pi / 180,
                axh.get_ylim()[1] * 0.65,
                f"{st_mode:.1f}$^o$",
                horizontalalignment="center",
                verticalalignment="baseline",
                fontdict={"size": self.font_size},
                bbox={"facecolor": box_color, "alpha": 0.25},
            )

            if estimate == "invariant":
                axh.set_title(
                    "Strike (Z)",
                    fontdict=self.font_dict,
                    bbox={"facecolor": self.color_inv, "alpha": 0.25},
                )
            # set pt axes properties
            elif estimate == "pt":
                axh.set_title(
                    "PT Azimuth",
                    fontdict=self.font_dict,
                    bbox={"facecolor": self.color_pt, "alpha": 0.25},
                )
            # set tipper axes properties
            elif estimate == "tipper":

                axh.set_title(
                    "Tipper Strike",
                    fontdict=self.font_dict,
                    bbox={"facecolor": self.color_tip, "alpha": 0.25},
                )
            # move title up a little to make room for labels
            axh.titleOffsetTrans._t = (0, 0.15)
            plt.setp(axh.yaxis.get_ticklabels(), visible=False)
            plt.setp(axh.xaxis.get_ticklabels(), visible=False)
        # remind the user what the assumptions of the strike angle are
        self.logger.info(
            "Note: North is assumed to be 0 and the strike angle is "
            + "measured clockwise positive."
        )

        plt.show()

    def plot(self, show=True):
        """
        plot Strike angles as rose plots

        """
        if self.strike_df is None:
            self.make_strike_df()

        # -----Plot Histograms of the strike angles-----------------------------

        # ------------------plot indivdual decades------------------------------
        if self.plot_type == 1:
            self._plot_per_period()

        # ------------------Plot strike angles for all period ranges------------
        elif self.plot_type == 2:
            self._plot_all_periods()
