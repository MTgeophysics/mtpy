# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as mcb
import matplotlib.colors as colors

from mtpy.imaging.mtplot_tools import PlotBaseProfile
import mtpy.imaging.mtcolors as mtcl

# ==============================================================================


class PlotResPhasePseudoSection(PlotBaseProfile):
    """
    plot a resistivity and phase pseudo section for different components

    Need to input one of the following lists:

    Arguments:
    ----------
        **fn_list** : list of strings
                     full paths to .edi files to plot. *default* is None

        **z_object** : list of class mtpy.core.z.Z
                      object of mtpy.core.z.  If this is input be sure the
                      attribute z.frequency is filled.  *default* is None

        **res_object_list** : list of class mtpy.mtplottools.ResPhase
                             list of ResPhase objects. *default* is None

        **mt_object** : class mtpy.imaging.mtplot.MTplot
                        object of mtpy.imaging.mtplot.MTplot
                        *default* is None
    Optional Key Words:
    -------------------

        *aspect*: [ 'equal' | 'auto' | float ]
                  aspect ratio of each subplot (height/width),
                  *default* is 'auto'

        *cb_orientation*: [ 'vertical' | 'horizontal' ]
                          orientation of colorbars, *default* is 'vertical'

        *cb_pad*: float
                  padding between edge of plot and edge of colorbar,
                  *default* is .0375

        *cb_position*: (res_position (x, y, ds, dy),
                        phase_position(x, y, dx, dy))
                        position of colorbars on the plot, need to input both
                        resistivity and phase positions.  *default* is None
                        and will automatically place the colorbars

        *cb_shrink*: float
                     factor to shrink the colorbar relative to the height of
                     the y-axis. *default* is 0.75

        *fig_dpi*: float
                   dots-per-inch resolution of figure, *default* is 300

        *fig_num*: int
                  number of figure instance, *default* is 1

        *fig_size*: (x, y) in inches
                    figure size in inches, *default* is (8, 4)

        *font_size*: float
                     size of font for axes tick labels, note label are +1.
                     *default* is 7

        *ftol*: float
                tolerance to extract periods relative to plot_period.
                *default* is 0.1

        *imshow_interp*: [ 'none' | 'nearest' | 'bilinear' | 'bicubic' |
                           'spline16' | 'spline36' | 'hanning' | 'hamming' |
                           'hermite' | 'kaiser' | 'quadric' | 'catrom' |
                           'gaussian' | 'bessel' | 'mitchell' | 'sinc' |
                           'lanczos']
                         defines the interpolation method if plot_style is
                         'imshow', if you use 'nearest' same as pcolormesh
                         except the lateral boxes are equal size instead of
                         set in a grid like pcolormesh.  Imshow just gives
                         a smoother interpretation of the pseudosection.
                         *default* is 'bicubic'

        *linedir*: [ 'ew' | 'ns' ]
                   predominant direction of profile line. *default* is 'ew'

        *period_limits*: (min_period, max_period)
                         limits on the plotting period in a linear scale.
                         *default* is None and extracts limits from data

        *phase_cmap*: [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 'mt_wh2bl' |
                        'mt_rd2bl' | 'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                        'mt_rd2gr2bl' ]
                      color map for phase plots, at the moment only supports
                      colormaps from mtpy.imaging.mtcolors, which are:

                           - 'mt_yl2rd' -> yellow to red
                           - 'mt_bl2yl2rd' -> blue to yellow to red
                           - 'mt_wh2bl' -> white to blue
                           - 'mt_rd2bl' -> red to blue
                           - 'mt_bl2wh2rd' -> blue to white to red
                           - 'mt_bl2gr2rd' -> blue to green to red *default*
                           - 'mt_rd2gr2bl' -> red to green to blue
                           - 'mt_seg_bl2wh2rd' -> discrete blue to
                                                 white to red

        *phase_limits*: (min_phase, max_phase)
                        minimum and maximum phase in degrees for coloring
                        *default* is (0, 90)

        *plot_period*: np.ndarray(periods)
                       array of periods to plot.  *default* is None, which
                       extracts the longest period from the input data
                       assuming that the longest is the most complete, if it
                       is not input manually.  If there are same lengths, it
                       picks the first one it finds.

        *plot_style*: [ 'imshow' | 'pcolormesh' ]
                      type of gridding for the plot. 'imshow' plots the data
                      as an image and can be interpolated, though the image
                      is stretched to the station spacing and plot_period, the
                      cells remain of equal size, so the interpolation might be
                      a little skewed.  For an accurate location of resistivity
                      values use pcolormesh, which can plot the data on an
                      irregular grid, but with no interpolation.
                      *default* is 'imshow'

        *plot_xx*: [ 'y' | 'n' ]
                  boolean to plot Z_xx, *default* is 'n'

        *plot_xy*: [ 'y' | 'n' ]
                  boolean to plot Z_xy, *default* is 'y'

        *plot_yx*: [ 'y' | 'n' ]
                  boolean to plot Z_yx, *default* is 'y'

        *plot_yy*: [ 'y' | 'n' ]
                  boolean to plot Z_yy, *default* is 'n'

        *plot_yn*: [ 'y' | 'n' ]
                  boolean to plot on instance creation, *default* is 'y'

        *res_cmap*: [ 'mt_yl2rd' | 'mt_bl2yl2rd' | 'mt_wh2bl' |
                        'mt_rd2bl' | 'mt_bl2wh2rd' | 'mt_seg_bl2wh2rd' |
                        'mt_rd2gr2bl' ]
                  color map for phase plots, at the moment only supports
                  colormaps from mtpy.imaging.mtcolors, which are:

                       - 'mt_yl2rd' -> yellow to red
                       - 'mt_bl2yl2rd' -> blue to yellow to red
                       - 'mt_wh2bl' -> white to blue
                       - 'mt_rd2bl' -> red to blue
                       - 'mt_bl2wh2rd' -> blue to white to red
                       - 'mt_bl2gr2rd' -> blue to green to red
                       - 'mt_rd2gr2bl' -> red to green to blue *default*
                       - 'mt_seg_bl2wh2rd' -> discrete blue to
                                             white to red
        *res_limits*: (min_resistivity, max_resistivity)
                      limits on resistivity in log scale, *default* is (0,3)

        *stationid*: (min, max)
                     min and max indicies to extract from each station name.
                     *default* is (0,4)

        *text_location*: (x, y)
                         location for text label for each resistivity subplot.
                         location is in relative coordinates of the data.
                         *default* is None, which locates the label in the
                         upper left hand corner

        *text_size*: [ size in points | 'xx-small' | 'x-small' | 'small' |
                      'medium' | 'large' | 'x-large' | 'xx-large' ]
                     size of text for subplot label, *default* is font_size

        *text_weight*: [ a numeric value in range 0-1000 | 'ultralight' |
                        'light' | 'normal' | 'regular' | 'book' | 'medium' |
                        'roman' | 'semibold' | 'demibold' | 'demi' | 'bold' |
                        'heavy' | 'extra bold' | 'black' ]
                       weight of text label font

        *text_xpad*: float
                     padding from x-axis, as a percentage of the xmin,
                     *default* is 0.95

        *text_ypad*: float
                     padding from y-axis, as a percentage of the ymin,
                     *default* is 0.95

        *xtickspace*: int
                      integer telling at what interval to place station names
                      as the tick labels, in case they are closely spaced.
                      *default* is 1

    """

    def __init__(self, tf_list, **kwargs):
        """
        Initialize parameters
        """

        # read in key word arguments and set defaults if none given
        self.tf_list = tf_list

        # --> set figure parameters
        self.aspect = kwargs.pop("aspect", "auto")

        self.xtickspace = kwargs.pop("xtickspace", 1)
        self.stationid = kwargs.pop("stationid", [0, 4])
        self.linedir = kwargs.pop("linedir", "ew")

        # --> set plots to plot and how to plot them
        self.plot_xx = False
        self.plot_xy = True
        self.plot_yx = True
        self.plot_yy = False
        self.plot_det = False
        self.plot_resistivity = True
        self.plot_phase = True

        # --> set plot limits
        self.cmap_limits = {
            "res_xx": (-1, 2),
            "res_xy": (0, 3),
            "res_yx": (0, 3),
            "res_yy": (-1, 2),
            "res_det": (0, 3),
            "phase_xx": (-180, 180),
            "phase_xy": (0, 100),
            "phase_yx": (0, 100),
            "phase_yy": (-180, 180),
            "phase_det": (0, 100),
        }

        self.label_dict = {
            "res_xx": "$\\rho_{xx}  \\mathrm{[\Omega m]}$",
            "res_xy": "$\\rho_{xy}  \\mathrm{[\Omega m]}$",
            "res_yx": "$\\rho_{yx}  \\mathrm{[\Omega m]}$",
            "res_yy": "$\\rho_{yy}  \\mathrm{[\Omega m]}$",
            "res_det": "$\\rho_{det}  \\mathrm{[\Omega m]}$",
            "phase_xx": "$\\phi_{xx}$",
            "phase_xy": "$\\phi_{xy}$",
            "phase_yx": "$\\phi_{yx}$",
            "phase_yy": "$\\phi_{yy}$",
            "phase_det": "$\\phi_{det}$",
        }

        # --> set colormaps Note only mtcolors is supported
        self.res_cmap = mtcl.cmapdict["mt_rd2gr2bl"]
        self.phase_cmap = mtcl.cmapdict["mt_bl2gr2rd"]

        for key, value in kwargs.items():
            setattr(self, key, value)

    def _get_period_array(self):
        """
        Get the period array to interpolate on to
        """

        tf_periods = np.array([tf.period for tf in self.tf_list]).flatten()

        p_min = np.log10(tf_periods.min())
        p_max = np.log10(tf_periods.max())

        n_periods = (np.ceil(p_max) - np.floor(p_min)) * 10

        return np.logspace(p_min, p_max, n_periods)

    def _get_n_rows(self):
        """
        Get the number of rows in the subplot

        :return: DESCRIPTION
        :rtype: TYPE

        """
        n = 0
        if self.plot_resistivity:
            n += 1
        if self.plot_phase:
            n += 1
        return n

    def _get_n_columns(self):
        """get the number of columns in the subplot"""
        n = 0

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if getattr(self, f"plot_{cc}"):
                n += 1

        return n

    def _get_n_subplots(self):
        """
        Get the subplot indices
        """
        nr = self._get_n_rows()
        nc = self._get_n_columns()

        subplot_dict = {
            "res_xx": None,
            "res_xy": None,
            "res_yx": None,
            "res_yy": None,
            "res_det": None,
            "phase_xx": None,
            "phase_xy": None,
            "phase_yx": None,
            "phase_yy": None,
            "phase_det": None,
        }

        plot_num = 0
        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_resistivity:
                if getattr(self, f"plot_{cc}"):
                    plot_num += 1
                    subplot_dict[f"res_{cc}"] = (nr, nc, plot_num)

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_phase:
                if getattr(self, f"plot_{cc}"):
                    plot_num += 1
                    subplot_dict[f"phase_{cc}"] = (nr, nc, plot_num)

        return subplot_dict

    def _get_subplots(self):
        """
        get the subplots

        :return: DESCRIPTION
        :rtype: TYPE

        """
        subplot_dict = self._get_n_subplots()
        ax_dict = {}

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_resistivity:
                comp = f"res_{cc}"
                if getattr(self, f"plot_{cc}"):
                    ax_dict[comp] = self.fig.add_subplot(
                        *subplot_dict[comp], aspect="equal"
                    )

            if self.plot_phase:
                comp = f"phase_{cc}"
                if getattr(self, f"plot_{cc}"):
                    ax_dict[comp] = self.fig.add_subplot(
                        *subplot_dict[comp], aspect="equal"
                    )

        share = [ax for comp, ax in ax_dict.items() if ax is not None]

        # share x and y across all subplots for easier zooming
        for ax in share[1:]:
            ax.sharex(share[0])
            ax.sharey(share[0])

        return ax_dict

    def _get_data_array(self):
        """
        get resistivity and phase values in the correct order according to
        offsets and periods.

        """

        self._get_profile_line()

        entries = []

        for ii, tf in enumerate(self.tf_list):
            offset = self._get_offset(tf)
            rp = tf.Z

            for ii, period in enumerate(tf.period):
                entry = {
                    "x": offset,
                    "y": period,
                    "res_xx": rp.res_xx[ii],
                    "res_xy": rp.res_xy[ii],
                    "res_yx": rp.res_yx[ii],
                    "res_yy": rp.res_yy[ii],
                    "res_det": rp.res_det[ii],
                    "phase_xx": rp.phase_xx[ii],
                    "phase_xy": rp.phase_xy[ii],
                    "phase_yx": rp.phase_yx[ii] + 180,
                    "phase_yy": rp.phase_yy[ii],
                    "phase_det": rp.phase_det[ii],
                }
                entries.append(entry)

        return pd.DataFrame(entries)

    def plot(self, show=True, get_rp_arrays=True):

        # --> set subplot spacing
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = 0.12
        plt.rcParams["figure.subplot.right"] = 0.90
        plt.rcParams["figure.subplot.bottom"] = 0.09
        plt.rcParams["figure.subplot.top"] = 0.98

        # get apparent resistivity and phase
        if get_rp_arrays:
            self.get_rp_arrays()
            if self.shift_yx_phase:
                self.phaseyx = self.phaseyx + 180
        # make a list of tuples to see how many subplots are needed
        ynlist = [
            self.plot_xx + "xx",
            self.plot_xy + "xy",
            self.plot_yx + "yx",
            self.plot_yy + "yy",
        ]
        reslist = [self.resxx, self.resxy, self.resyx, self.resyy]
        phaselist = [self.phasexx, self.phasexy, self.phaseyx, self.phaseyy]
        plist = [
            (yn[1:], res, phase)
            for yn, res, phase in zip(ynlist, reslist, phaselist)
            if yn[0] == "y"
        ]

        # make a general subplot array
        gs = gridspec.GridSpec(
            2, len(plist), height_ratios=[1, 1], hspace=0.00, wspace=0.025
        )

        # get ylimits for plot
        if self.period_limits is None:
            self.period_limits = (
                self.plot_period.min(),
                self.plot_period.max(),
            )
        font_dict = {"size": self.font_size + 2, "weight": "bold"}
        ns = len(self.station_list)
        # --> plot data
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # plot as a mesh where the data are left as blocks
        if self.plot_style == "pcolormesh":
            # need to add another element at the end of the array so pcolor
            # will plot the full array
            # first, get median station spacing
            mss = np.median(
                np.abs(self.offset_list[1:] - self.offset_list[:-1])
            )
            xgrid_edges = np.mean(
                [self.offset_list[1:], self.offset_list[:-1]], axis=0
            )
            xgrid = np.hstack(
                [self.offset_list[:1], xgrid_edges, self.offset_list[-1:]]
            )

            ygrid_edges = 10 ** np.mean(
                [
                    np.log10(self.plot_period[1:]),
                    np.log10(self.plot_period[:-1]),
                ],
                axis=0,
            )
            ygrid = np.hstack(
                [self.plot_period[:1], ygrid_edges, self.plot_period[-1:]]
            )
            xgrid, ygrid = np.meshgrid(xgrid, ygrid)

            #            xgrid, ygrid = np.meshgrid(np.append(self.offset_list,
            #                                                 self.offset_list[-1]*1.1),
            #                                       np.append(self.plot_period,
            #                                                 self.plot_period[-1]*1.1))

            for ii, tt in enumerate(plist):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])

                # plot apparent resistivity
                axr.pcolormesh(
                    xgrid,
                    ygrid,
                    tt[1],  # np.flipud(tt[1]),
                    cmap=self.res_cmap,
                    vmin=self.res_limits[0],
                    vmax=self.res_limits[1],
                )

                axr.set_aspect(self.aspect)
                axp.set_aspect(self.aspect)

                axr.set_xticks(
                    self.offset_list[list(range(0, ns, self.xtickspace))]
                )
                if self.xtickspace != 1:
                    axr.set_xticks(self.offset_list)  # , minor=True)
                plt.setp(axr.get_xticklabels(), visible=False)
                if self.show_grid:
                    axr.grid(which="major", alpha=0.25)
                axr.set_yscale("log", nonposy="clip")
                axr.set_xlim(
                    self.offset_list.min(), self.offset_list.max() * 1.1
                )
                axr.set_ylim(self.period_limits)

                # label the plot with a text box
                if self.text_location is None:
                    txloc = self.offset_list.min() * self.text_xpad
                    tyloc = self.period_limits[1] * self.text_ypad
                else:
                    txloc = self.text_location[0]
                    tyloc = self.text_location[1]
                self.text = axr.text(
                    txloc,
                    tyloc,
                    "$Z_{" + tt[0] + "}$",
                    fontdict={
                        "size": self.text_size,
                        "weight": self.text_weight,
                    },
                    verticalalignment="top",
                    horizontalalignment="left",
                    bbox={"facecolor": "white", "alpha": 0.5},
                )

                # plot phase
                axp.pcolormesh(
                    xgrid,
                    ygrid,
                    tt[2],  # np.flipud(tt[2]),
                    cmap=self.phase_cmap,
                    vmin=self.phase_limits[0],
                    vmax=self.phase_limits[1],
                )
                if self.show_grid:
                    axp.grid(which="major", alpha=0.25)
                axp.set_xticks(
                    self.offset_list[list(range(0, ns, self.xtickspace))]
                )
                axp.set_xticklabels(
                    [
                        self.station_list[st]
                        for st in range(0, ns, self.xtickspace)
                    ],
                    rotation=self.station_label_rotation,
                    fontsize=self.text_size,
                )
                if self.xtickspace != 1:
                    axp.set_xticks(self.offset_list, minor=True)
                axp.set_yscale("log", nonposy="clip")
                axp.set_xlim(
                    self.offset_list.min(), self.offset_list.max() * 1.1
                )
                axp.set_ylim(self.period_limits)
                if ii == 0:
                    axp.set_ylabel("Period (s)", font_dict)
                    axr.set_ylabel("Period (s)", font_dict)
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                # add colorbars
                if ii == len(plist) - 1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    # add colorbar for res
                    axrpos = axr.get_position()

                    # set position just to the right of the figure
                    if self.cb_position is None:
                        cbr_position = (
                            axrpos.bounds[0] + axrpos.bounds[2] + self.cb_pad,
                            axrpos.bounds[1] + 0.05,
                            0.015,
                            axrpos.bounds[3] * self.cb_shrink,
                        )
                    else:
                        cbr_position = self.cb_position[0]
                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(
                        self.cbaxr,
                        cmap=self.res_cmap,
                        norm=colors.Normalize(vmin=cminr, vmax=cmaxr),
                        orientation=self.cb_orientation,
                    )
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)

                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax + 1))
                    cbr_ticklabels = [
                        mtpl.labeldict[ll]
                        for ll in np.arange(tkrmin, tkrmax + 1)
                    ]

                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position("right")
                    self.cbr.ax.yaxis.set_label_coords(1.35, 0.5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis="y", direction="in", pad=1)
                    self.cbr.set_label(
                        "App. Res ($\Omega \cdot$m)",
                        fontdict={"size": self.font_size},
                    )

                    # --> add colorbar for phase
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]

                    axppos = axp.get_position()

                    # set position just to the right of the figure
                    if self.cb_position is None:
                        cbp_position = (
                            axppos.bounds[0] + axppos.bounds[2] + self.cb_pad,
                            axppos.bounds[1] + 0.05,
                            0.015,
                            axppos.bounds[3] * self.cb_shrink,
                        )
                    else:
                        cbp_position = self.cb_position[1]
                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(
                        self.cbaxp,
                        cmap=self.phase_cmap,
                        norm=colors.Normalize(vmin=cminp, vmax=cmaxp),
                        orientation=self.cb_orientation,
                    )
                    self.cbp.set_ticks([cminp, (cmaxp - cminp) / 2, cmaxp])
                    self.cbp.set_ticklabels(
                        [
                            "{0:.0f}".format(cminp),
                            "{0:.0f}".format((cmaxp - cminp) / 2),
                            "{0:.0f}".format(cmaxp),
                        ]
                    )
                    self.cbp.ax.yaxis.set_label_position("right")
                    self.cbp.ax.yaxis.set_label_coords(1.35, 0.5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis="y", direction="in", pad=0.5)
                    self.cbp.set_label(
                        "Phase (deg)", fontdict={"size": self.font_size}
                    )
                # make axes attributes for user editing
                if tt == "xx":
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt == "xy":
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt == "yx":
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt == "yy":
                    self.ax_ryy = axr
                    self.ax_pyy = axp
            if show:
                plt.show()
        # plot data as an image which can have interpolation
        elif self.plot_style == "imshow":
            # make ticks simulate a log scale in the y direction
            # --> set major and minor ticks with appropriate labels
            major_yticks = np.arange(
                np.ceil(
                    np.log10(self.period_limits[0])
                    if self.period_limits[0] != 0
                    else 0
                ),
                np.floor(
                    np.log10(self.period_limits[1])
                    if self.period_limits[0] != 0
                    else 0
                )
                + 1,
            )

            # make minor ticks look like they are on a log scale
            minor_yticks = []
            for ll in major_yticks:
                minor_yticks += [np.arange(1, 10) * 10**ll]
            minor_yticks = np.array(minor_yticks)
            minor_yticks = np.log10(minor_yticks.flatten())

            # set ticklabels as 10**
            yticklabels = [mtpl.labeldict[ll] for ll in major_yticks]

            for ii, tt in enumerate(plist):
                axr = self.fig.add_subplot(gs[0, ii])
                axp = self.fig.add_subplot(gs[1, ii])

                # plot apparent resistivity
                axr.imshow(
                    tt[1],
                    cmap=self.res_cmap,
                    vmin=self.res_limits[0],
                    vmax=self.res_limits[1],
                    aspect=self.aspect,
                    interpolation=self.imshow_interp,
                    extent=(
                        self.offset_list.min(),
                        self.offset_list.max(),
                        np.log10(self.plot_period.min()),
                        np.log10(self.plot_period.max()),
                    ),
                )

                # set x ticks but remove labels
                axr.set_xticks(
                    self.offset_list[list(range(0, ns, self.xtickspace))]
                )
                if self.xtickspace != 1:
                    axr.set_xticks(self.offset_list, minor=True)
                plt.setp(axr.get_xticklabels(), visible=False)

                # set y-axis ticks
                axr.yaxis.set_ticks(major_yticks)
                axr.yaxis.set_ticks(minor_yticks, minor=True)
                axr.set_yticklabels(yticklabels[::-1])

                axr.grid(which="major", alpha=0.25)
                axr.set_xlim(self.offset_list.min(), self.offset_list.max())
                axr.set_ylim(
                    np.log10(self.period_limits[0])
                    if self.period_limits[0] != 0
                    else 0,
                    np.log10(self.period_limits[1])
                    if self.period_limits[0] != 0
                    else 0,
                )

                # label the plot with a text box
                if self.text_location is None:
                    txloc = self.offset_list.min() * self.text_xpad
                    tyloc = (
                        np.log10(self.period_limits[1] * self.text_ypad)
                        if self.period_limits[0] != 0
                        else 0
                    )
                else:
                    txloc = self.text_location[0]
                    tyloc = self.text_location[1]
                self.text = axr.text(
                    txloc,
                    tyloc,
                    "$Z_{" + tt[0] + "}$",
                    fontdict={
                        "size": self.text_size,
                        "weight": self.text_weight,
                    },
                    verticalalignment="top",
                    horizontalalignment="left",
                    bbox={"facecolor": "white", "alpha": 0.5},
                )

                if ii == 0:
                    axr.set_ylabel("Period (s)", font_dict)
                # plot phase
                axp.imshow(
                    tt[2],
                    cmap=self.phase_cmap,
                    vmin=self.phase_limits[0],
                    vmax=self.phase_limits[1],
                    aspect=self.aspect,
                    interpolation=self.imshow_interp,
                    extent=(
                        self.offset_list.min(),
                        self.offset_list.max(),
                        np.log10(self.plot_period.min()),
                        np.log10(self.plot_period.max()),
                    ),
                )

                axp.grid(which="major", alpha=0.25)
                axp.set_xticks(
                    self.offset_list[list(range(0, ns, self.xtickspace))]
                )
                axp.set_xticklabels(
                    [
                        self.station_list[st]
                        for st in range(0, ns, self.xtickspace)
                    ]
                )
                if self.xtickspace != 1:
                    axp.set_xticks(self.offset_list, minor=True)
                # remove tick labels if not the first subplot
                if ii != 0:
                    plt.setp(axr.get_yticklabels(), visible=False)
                    plt.setp(axp.get_yticklabels(), visible=False)
                # set y-axis ticks
                axp.yaxis.set_ticks(major_yticks)
                axp.yaxis.set_ticks(minor_yticks, minor=True)
                axp.set_yticklabels(yticklabels[::-1])

                axp.set_xlim(self.offset_list.min(), self.offset_list.max())
                axp.set_ylim(
                    np.log10(
                        self.period_limits[0]
                        if self.period_limits[0] != 0
                        else 0
                    ),
                    np.log10(self.period_limits[1])
                    if self.period_limits[0] != 0
                    else 0,
                )

                if ii == 0:
                    axp.set_ylabel("Period (s)", font_dict)
                # add colorbars
                if ii == len(plist) - 1:
                    cminr = self.res_limits[0]
                    cmaxr = self.res_limits[1]
                    # add colorbar for res
                    axrpos = axr.get_position()

                    # set position just to the right of the figure
                    if self.cb_position is None:
                        cbr_position = (
                            axrpos.bounds[0] + axrpos.bounds[2] + self.cb_pad,
                            axrpos.bounds[1] + 0.05,
                            0.015,
                            axrpos.bounds[3] * self.cb_shrink,
                        )
                    else:
                        cbr_position = self.cb_position[0]
                    self.cbaxr = self.fig.add_axes(cbr_position)
                    self.cbr = mcb.ColorbarBase(
                        self.cbaxr,
                        cmap=self.res_cmap,
                        norm=colors.Normalize(vmin=cminr, vmax=cmaxr),
                        orientation=self.cb_orientation,
                    )
                    tkrmin = np.ceil(cminr)
                    tkrmax = np.floor(cmaxr)

                    self.cbr.set_ticks(np.arange(tkrmin, tkrmax + 1))
                    cbr_ticklabels = [
                        mtpl.labeldict[ll]
                        for ll in np.arange(tkrmin, tkrmax + 1)
                    ]

                    self.cbr.set_ticklabels(cbr_ticklabels)
                    self.cbr.ax.yaxis.set_label_position("right")
                    self.cbr.ax.yaxis.set_label_coords(1.35, 0.5)
                    self.cbr.ax.yaxis.tick_left()
                    self.cbr.ax.tick_params(axis="y", direction="in", pad=1)
                    self.cbr.set_label(
                        "App. Res ($\Omega \cdot$m)",
                        fontdict={"size": self.font_size},
                    )

                    # --> add colorbar for phase
                    cminp = self.phase_limits[0]
                    cmaxp = self.phase_limits[1]

                    axppos = axp.get_position()

                    # set position just to the right of the figure
                    if self.cb_position is None:
                        cbp_position = (
                            axppos.bounds[0] + axppos.bounds[2] + self.cb_pad,
                            axppos.bounds[1] + 0.05,
                            0.015,
                            axppos.bounds[3] * self.cb_shrink,
                        )
                    else:
                        cbp_position = self.cb_position[1]
                    self.cbaxp = self.fig.add_axes(cbp_position)
                    self.cbp = mcb.ColorbarBase(
                        self.cbaxp,
                        cmap=self.phase_cmap,
                        norm=colors.Normalize(vmin=cminp, vmax=cmaxp),
                        orientation=self.cb_orientation,
                    )
                    self.cbp.set_ticks([cminp, (cmaxp - cminp) / 2, cmaxp])
                    self.cbp.set_ticklabels(
                        [
                            "{0:.0f}".format(cminp),
                            "{0:.0f}".format((cmaxp - cminp) / 2),
                            "{0:.0f}".format(cmaxp),
                        ]
                    )
                    self.cbp.ax.yaxis.set_label_position("right")
                    self.cbp.ax.yaxis.set_label_coords(1.35, 0.5)
                    self.cbp.ax.yaxis.tick_left()
                    self.cbp.ax.tick_params(axis="y", direction="in", pad=0.5)
                    self.cbp.set_label(
                        "Phase (deg)", fontdict={"size": self.font_size}
                    )
                if tt[0] == "xx":
                    self.ax_rxx = axr
                    self.ax_pxx = axp
                elif tt[0] == "xy":
                    self.ax_rxy = axr
                    self.ax_pxy = axp
                elif tt[0] == "yx":
                    self.ax_ryx = axr
                    self.ax_pyx = axp
                elif tt[0] == "yy":
                    self.ax_ryy = axr
                    self.ax_pyy = axp
            if show:
                plt.show()

    def save_plot(
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
                            save_fn/station_name_ResPhase.file_format

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
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.save_plot(r'/home/MT/figures', file_format='jpg')

        """

        if fig_dpi is None:
            fig_dpi = self.fig_dpi
        if os.path.isdir(save_fn) == False:
            file_format = save_fn[-3:]
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
            )
            # plt.clf()
            # plt.close(self.fig)
        else:
            save_fn = os.path.join(
                save_fn,
                self._mt.station + "_ResPhasePseudoSection." + file_format,
            )
            self.fig.savefig(
                save_fn,
                dpi=fig_dpi,
                format=file_format,
                orientation=orientation,
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
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.xy_color = (.5,.5,.9)
            >>> p1.xy_marker = '*'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def writeTextFiles(self, save_path=None, ptol=0.10):
        """
        This will write text files for all the phase tensor parameters
        """

        if save_path is None:
            try:
                svpath = os.path.dirname(self.mt_list[0].fn)
            except TypeError:
                raise IOError("Need to input save_path, could not find a path")
        else:
            svpath = save_path
        if self.resxy.mean() == 0:
            self.get_rp_arrays()
        header_list = (
            ["{0:^10}".format("period(s)")]
            + ["{0:^8}".format(ss) for ss in self.station_list]
            + ["\n"]
        )

        fn_dict = {
            "resxx": self.resxx,
            "resxy": self.resxy,
            "resyx": self.resyx,
            "resyy": self.resyy,
            "phasexx": self.phasexx,
            "phasexy": self.phasexy,
            "phaseyx": self.phaseyx,
            "phaseyy": self.phaseyy,
        }

        # write the arrays into lines properly formatted
        t1_kwargs = {
            "spacing": "{0:^10} ",
            "value_format": "{0:.2e}",
            "append": False,
            "add": False,
        }

        tr_kwargs = {
            "spacing": "{0:^8}",
            "value_format": "{0: .2f}",
            "append": False,
            "add": False,
        }

        tp_kwargs = {
            "spacing": "{0:^8}",
            "value_format": "{0: .2f}",
            "append": False,
            "add": False,
        }

        for key in list(fn_dict.keys()):
            fid = file(os.path.join(svpath, "PseudoSection." + key), "w")
            fid.write("".join(header_list))
            for ii, per in enumerate(self.plot_period):
                if key[0] == "r":
                    line = (
                        [mtpl.make_value_str(per, **t1_kwargs)]
                        + [
                            mtpl.make_value_str(rr, **tr_kwargs)
                            for rr in fn_dict[key][ii]
                        ]
                        + ["\n"]
                    )
                elif key[0] == "p":
                    line = (
                        [mtpl.make_value_str(per, **t1_kwargs)]
                        + [
                            mtpl.make_value_str(rr, **tp_kwargs)
                            for rr in fn_dict[key][ii]
                        ]
                        + ["\n"]
                    )
                fid.write("".join(line))
            fid.close()
        print(
            "Wrote files to: "
            + os.path.join(svpath, "PseudoSection.component")
        )

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return "Plots Resistivity and phase as a pseudo section."
