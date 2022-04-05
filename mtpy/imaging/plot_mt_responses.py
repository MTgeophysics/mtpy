# -*- coding: utf-8 -*-
"""
plots multiple MT responses simultaneously

Created on Thu May 30 17:02:39 2013
@author: jpeacock-pr

YG: the code there is massey, todo may need to rewrite it sometime

"""

# ============================================================================

import matplotlib.colorbar as mcb
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

import mtpy.imaging.mtcolors as mtcl
from mtpy.imaging.mtplot_tools import (
    PlotBase,
    plot_pt_lateral,
    get_log_tick_labels,
    plot_resistivity,
    plot_phase,
    plot_tipper_lateral,
)


# ============================================================================


class PlotMultipleResponses(PlotBase):
    """
    plots multiple MT responses simultaneously either in single plots or in
    one plot of sub-figures or in a single plot with subfigures for each
    component.


    Arguments:
    ----------
        **fn_list** : list of filenames to plot
                     ie. [fn_1, fn_2, ...], *default* is None

         **z_object_list** : list of mtpy.core.z.Z objects
                            *default* is None

         **res_object_list** : list of mtpy.imaging.mtplot.ResPhase objects
                              *default* is None

         **tipper_object_list** : list of mtpy.imaging.mtplot.Tipper objects
                                 *default* is None

         **mt_object_list** : list of mtpy.imaging.mtplot.MTplot objects
                             *default* is None

        **fig_num** : int
                     figure number
                     *default* is 1
        **fig_size** : [width, height] of figure size in inches

        **rot_z** : float or np.ndarray
                   rotation angle of impedance tensor (deg or radians),
                   *Note* : rotaion is clockwise positive
                   *default* is 0
                   Can input so each station is rotated at a constant angle or
                   each period is rotated differently, or both.

        **plot_num** : [ 1 | 2 | 3 ]
                        * 1 for just Ex/By and Ey/Bx *default*
                        * 2 for all 4 components
                        * 3 for off diagonal plus the determinant

        **plot_style** : [ '1' | 'all' | 'compare' ]
                        determines the plotting style:
                            * '1' for plotting each station in a different
                                  figure. *default*

                            * 'all' for plotting each station in a subplot
                                    all in the same figure

                            * 'compare' for comparing the responses all in
                                        one plot.  Here the responses are
                                        colored from dark to light.  This
                                        plot can get messy if too many stations
                                        are plotted.




    """

    def __init__(self, tf_list, **kwargs):
        """
        Initialize parameters
        """
        self.plot_num = 1
        self.plot_style = "1"

        super().__init__(**kwargs)

        self.tf_list = tf_list

        self.plot_dict = dict(
            [
                (kk, vv)
                for kk, vv in zip(
                    ["tip", "pt", "strike", "skew"],
                    [self.plot_tipper, self.plot_pt, self.plot_strike, self.plot_skew,],
                )
            ]
        )

        # set arrow properties
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = 0.5

        # ellipse_properties
        self.ellipse_size = 0.25

        # plot on initializing
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
        for tf in self.tf_list:
            tf.rotation_angle = value
        self._rotation_angle = value

    def _has_tipper(self, t_obj):
        if self.plot_tipper.find("y") >= 0:
            if t_obj is None or (t_obj.tipper == 0 + 0j).all():
                self._logger.info(f"No Tipper data for station {self.station}")
                self.plot_tipper = "n"

    def _has_pt(self, pt):
        if self.plot_pt:
            # if np.all(self.Z.z == 0 + 0j) or self.Z is None:
            if pt is None:  # no phase tensor object provided
                self._logger.info(f"No PT data for station {self.station}")
                self.plot_pt = False

    def _plot_resistivity(self, axr, period, z_obj, mode="od", index=0):

        if mode == "od":
            comps = ["xy", "yx"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
        elif mode == "d":
            comps = ["xx", "yy"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
        elif mode == "det":
            comps = ["det"]
            props = [self.det_error_bar_properties]
        res_limits = self.set_resistivity_limits(z_obj.resistivity, mode=mode)
        x_limits = self.set_period_limits(period)

        eb_list = []
        label_list = []
        for comp, prop in zip(comps, props):
            ebax = plot_resistivity(
                axr,
                period,
                getattr(z_obj, f"res_{comp}"),
                getattr(z_obj, f"res_err_{comp}"),
                **prop,
            )
            eb_list.append(ebax[0])
            label_list.append(f"$Z_{comp}$")
        # --> set axes properties
        plt.setp(axr.get_xticklabels(), visible=False)

        axr.set_yscale("log", nonpositive="clip")
        axr.set_xscale("log", nonpositive="clip")
        axr.set_xlim(x_limits)
        axr.set_ylim(res_limits)
        axr.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)

        if index == 0:
            axr.set_ylabel(
                "App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=self.font_dict
            )
        axr.legend(
            eb_list,
            label_list,
            loc=3,
            markerscale=1,
            borderaxespad=0.01,
            labelspacing=0.07,
            handletextpad=0.2,
            borderpad=0.02,
        )
        return axr, eb_list, label_list

    def _plot_phase(self, axp, period, z_obj, mode="od", index=0):
        if mode == "od":
            comps = ["xy", "yx"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
        elif mode == "d":
            comps = ["xx", "yy"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
        elif mode == "det":
            comps = ["det"]
            props = [self.det_error_bar_properties]
        phase_limits = self.set_phase_limits(z_obj.phase, mode=mode)
        x_limits = self.set_period_limits(period)

        for comp, prop in zip(comps, props):
            phase = getattr(z_obj, f"phase_{comp}")
            if comp == "yx":
                phase += 180
            ebax = plot_resistivity(
                axr, period, phase, getattr(z_obj, f"phase_err_{comp}"), **prop,
            )
        # --> set axes properties
        if index == 0:
            axp.set_ylabel("Phase (deg)", self.font_dict)
        axp.set_xscale("log", nonpositive="clip")

        axp.set_ylim(phase_limits)
        if phase_limits[0] < -10 or phase_limits[1] > 100:
            axp.yaxis.set_major_locator(MultipleLocator(30))
            axp.yaxis.set_minor_locator(MultipleLocator(10))
        else:
            axp.yaxis.set_major_locator(MultipleLocator(15))
            axp.yaxis.set_minor_locator(MultipleLocator(5))
        axp.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)

        return axp

    def _plot_determinant(self):

        self.axr.legend(
            (self.ebxyr[0], self.ebyxr[0], self.ebdetr[0]),
            ("$Z_{xy}$", "$Z_{yx}$", "$\det(\mathbf{\hat{Z}})$"),
            loc=3,
            markerscale=1,
            borderaxespad=0.01,
            labelspacing=0.07,
            handletextpad=0.2,
            borderpad=0.02,
        )

    def _plot_tipper(self, axt, period, t_obj, index=0):
        axt = plot_tipper_lateral(
            axt,
            t_obj,
            self.plot_tipper,
            self.arrow_real_properties,
            self.arrow_imag_properties,
            self.font_size,
        )

        axt.set_xlabel("Period (s)", fontdict=self.font_dict)

        axt.yaxis.set_major_locator(MultipleLocator(0.2))
        axt.yaxis.set_minor_locator(MultipleLocator(0.1))
        axt.set_xlabel("Period (s)", fontdict=self.font_dict)
        axt.set_ylabel("Tipper", fontdict=self.font_dict)

        # set th xaxis tick labels to invisible
        if self.plot_pt:
            plt.setp(axt.xaxis.get_ticklabels(), visible=False)
            axt.set_xlabel("")

    def _plot_pt(self):
        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt:

            color_array = self.get_pt_color_array(self.pt)

            # -------------plot ellipses-----------------------------------
            self.cbax, self.cbpt, = plot_pt_lateral(
                self.axpt, self.pt, color_array, self.ellipse_properties, self.fig,
            )

            # ----set axes properties-----------------------------------------------
            # --> set tick labels and limits
            self.axpt.set_xlim(np.log10(self.x_limits[0]), np.log10(self.x_limits[1]))

            tklabels, xticks = get_log_tick_labels(self.axpt)

            self.axpt.set_xticks(xticks)
            self.axpt.set_xticklabels(tklabels, fontdict={"size": self.font_size})
            self.axpt.set_xlabel("Period (s)", fontdict=self.font_dict)

            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            self.axpt.set_xlim(np.log10(self.x_limits[0]), np.log10(self.x_limits[1]))
            self.axpt.grid(
                True, alpha=0.25, which="major", color=(0.25, 0.25, 0.25), lw=0.25
            )

            plt.setp(self.axpt.get_yticklabels(), visible=False)

            self.cbpt.set_label(
                self.cb_label_dict[self.ellipse_colorby],
                fontdict={"size": self.font_size},
            )

    def _setup_subplots(self, n_stations=1, index=0):
        # create a dictionary for the number of subplots needed
        pdict = {"res": 0, "phase": 1}
        # start the index at 2 because resistivity and phase is permanent for
        # now
        axr = None
        axp = None
        axr2 = None
        axp2 = None
        axt = None
        axpt = None

        index = 0
        nrows = 1
        if self.plot_tipper.find("y") >= 0:
            pdict["tip"] = index
            index += 1
            nrows = 2
        if self.plot_pt:
            pdict["pt"] = index
            nrows = 2
            index += 1
        if nrows == 1:
            hr = [1]
        elif nrows == 2:
            hr = [2, 1]
        gs_master = gridspec.GridSpec(nrows, n_stations, hspace=0.15, height_ratios=hr)
        gs_rp = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=gs_master[index],
            height_ratios=[2, 1.5],
            hspace=0.05,
            wspace=0.15,
        )
        if nrows == 2:
            gs_aux = gridspec.GridSpecFromSubplotSpec(
                index, 1, subplot_spec=gs_master[1], hspace=0.05
            )
        # --> make figure for xy,yx components
        if self.plot_num == 1 or self.plot_num == 3:
            # set label coordinates
            label_coords = (-0.075, 0.5)

            # --> create the axes instances
            # apparent resistivity axis
            axr = self.fig.add_subplot(gs_rp[0, :])

            # phase axis that shares period axis with resistivity
            axp = self.fig.add_subplot(gs_rp[1, :], sharex=axr)
        # --> make figure for all 4 components
        elif self.plot_num == 2:
            # set label coordinates
            label_coords = (-0.095, 0.5)

            # --> create the axes instances
            # apparent resistivity axis
            axr = self.fig.add_subplot(gs_rp[0, 0])
            axr2 = self.fig.add_subplot(gs_rp[0, 1], sharex=axr)
            axr2.yaxis.set_label_coords(-0.1, 0.5)

            # phase axis that shares period axis with resistivity
            axp = self.fig.add_subplot(gs_rp[1, 0], sharex=axr)
            axp2 = self.fig.add_subplot(gs_rp[1, 1], sharex=axr)
            axp2.yaxis.set_label_coords(-0.1, 0.5)
        # set albel coordinates
        axr.yaxis.set_label_coords(label_coords[0], label_coords[1])
        axp.yaxis.set_label_coords(label_coords[0], label_coords[1])

        # --> plot tipper
        if self.plot_tipper.find("y") >= 0:
            axt = self.fig.add_subplot(gs_aux[pdict["tip"], :],)
            axt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        # --> plot phase tensors
        if self.plot_pt:
            # can't share axis because not on the same scale
            # Removed aspect = "equal" for now, it flows better, if you want
            # a detailed analysis look at plot pt
            axpt = self.fig.add_subplot(gs_aux[pdict["pt"], :])
            axpt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        return axr, axp, axr2, axp2, axt, axpt, label_coords

    def _plot_all(self):
        ns = len(self.mt_list)

        # set figure size according to what the plot will be.
        if self.fig_size is None:
            if self.plot_num == 1 or self.plot_num == 3:
                self.fig_size = [ns * 4, 6]
            elif self.plot_num == 2:
                self.fig_size = [ns * 8, 6]
        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # make subplots as columns for all stations that need to be plotted
        gs0 = gridspec.GridSpec(1, ns)

        # space out the subplots
        gs0.update(hspace=0.025, wspace=0.025, left=0.085)

        for ii, mt in enumerate(self.tf_list):
            axr, axp, axr2, axp2, axt, axpt, label_coords = self._setup_subplots(
                n_stations=ns, index=ii
            )

            x_limits = self.set_period_limits(mt.period)
            res_limits = self.set_resistivity_limits(mt.Z.resistivity)

            # ---------plot the apparent resistivity----------------------
            # --> plot as error bars and just as points xy-blue, yx-red
            # res_xy
            ebxyr = plot_resistivity(
                axr,
                mt.period,
                mt.Z.res_xy,
                mt.Z.res_err_xy,
                **self.xy_error_bar_properties,
            )

            # res_yx
            ebyxr = plot_resistivity(
                axr,
                mt.period,
                mt.Z.res_yx,
                mt.Z.res_err_yx,
                **self.yx_error_bar_properties,
            )

            # --> set axes properties
            plt.setp(axr.get_xticklabels(), visible=False)
            axr.set_yscale("log", nonpositive="clip")
            axr.set_xscale("log", nonpositive="clip")
            axr.set_xlim(x_limits)
            axr.set_ylim(res_limits)
            axr.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)
            if ii == 0:
                axr.set_ylabel(
                    "App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=self.font_dict
                )
                axr.legend(
                    (ebxyr[0], ebyxr[0]),
                    ("$Z_{xy}$", "$Z_{yx}$"),
                    loc=3,
                    markerscale=1,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.02,
                )
            else:
                plt.setp(axr.get_yticklabels(), visible=False)
            # -----Plot the phase----------------------------------------
            # phase_xy
            phase_limits = self.set_period_limits(mt.Z.phase)
            ebxyp = plot_phase(
                axp,
                mt.period,
                mt.Z.phase_xy,
                mt.Z.phase_err_xy,
                **self.xy_error_bar_properties,
            )

            # phase_yx:
            ebyxp = plot_phase(
                axp,
                mt.period,
                mt.Z.phase_yx + 180,
                mt.Z.phase_err_yx,
                **self.yx_error_bar_properties,
            )

            # --> set axes properties
            if ii == 0:
                axp.set_ylabel("Phase (deg)", self.font_dict)
            else:
                plt.setp(axp.get_yticklabels(), visible=False)
                axp.set_xlabel("Period (s)", self.font_dict)
            axp.set_xscale("log", nonpositive="clip")
            if self.phase_limits is None:
                self.phase_limits = (-179.9, 179.9)
            axp.set_ylim(self.phase_limits)
            axp.set_xlim(self.x_limits)
            axp.yaxis.set_major_locator(MultipleLocator(15))
            axp.yaxis.set_minor_locator(MultipleLocator(5))
            axp.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)

            tklabels = [
                mtpl.labeldict[tt]
                for tt in np.arange(
                    np.log10(self.xlimits[0]), np.log10(self.xlimits[1]) + 1
                )
            ]
            tklabels[0] = ""
            tklabels[-1] = ""

            axp.set_xticklabels(tklabels, fontdict={"size": self.font_size})

            # -----plot tipper--------------------------------------------
            if self._plot_tipper.find("y") == 0:
                plt.setp(axp.xaxis.get_ticklabels(), visible=False)

                tp = mt.Tipper

                txr = tp.mag_real * np.sin(
                    tp.angle_real * np.pi / 180 + np.pi * self.arrow_direction
                )
                tyr = tp.mag_real * np.cos(
                    tp.angle_real * np.pi / 180 + np.pi * self.arrow_direction
                )

                txi = tp.mag_imag * np.sin(
                    tp.angle_imag * np.pi / 180 + np.pi * self.arrow_direction
                )
                tyi = tp.mag_imag * np.cos(
                    tp.angle_imag * np.pi / 180 + np.pi * self.arrow_direction
                )

                nt = len(txr)

                tiplist = []
                tiplabel = []

                for aa in range(nt):
                    xlenr = txr[aa] * mt.period[aa]
                    xleni = txi[aa] * mt.period[aa]

                    # --> plot real arrows
                    if self._plot_tipper.find("r") > 0:
                        axt.arrow(
                            np.log10(mt.period[aa]),
                            0,
                            xlenr,
                            tyr[aa],
                            lw=self.arrow_lw,
                            facecolor=self.arrow_color_real,
                            edgecolor=self.arrow_color_real,
                            head_width=self.arrow_head_width,
                            head_length=self.arrow_head_length,
                            length_includes_head=False,
                        )

                        if aa == 0:
                            line1 = axt.plot(0, 0, self.arrow_color_real)
                            tiplist.append(line1[0])
                            tiplabel.append("real")
                    # --> plot imaginary arrows
                    if self.plot_tipper.find("i") > 0:
                        axt.arrow(
                            np.log10(mt.period[aa]),
                            0,
                            xleni,
                            tyi[aa],
                            facecolor=self.arrow_color_imag,
                            edgecolor=self.arrow_color_imag,
                            lw=self.arrow_lw,
                            head_width=self.arrow_head_width,
                            head_length=self.arrow_head_length,
                            length_includes_head=False,
                        )
                        if aa == 0:
                            line2 = axt.plot(0, 0, self.arrow_color_imag)
                            tiplist.append(line2[0])
                            tiplabel.append("imag")
                # make a line at 0 for reference
                axt.plot(mt.period, [0] * nt, "k", lw=0.5)

                if ii == 0:
                    axt.legend(
                        tiplist,
                        tiplabel,
                        loc="upper left",
                        markerscale=1,
                        borderaxespad=0.01,
                        labelspacing=0.07,
                        handletextpad=0.2,
                        borderpad=0.1,
                        prop={"size": self.font_size},
                    )

                    axt.set_ylabel("Tipper", fontdict=self.font_dict)
                else:
                    plt.setp(axt.get_yticklabels(), visible=False)
                # set axis properties
                axt.yaxis.set_major_locator(MultipleLocator(0.2))
                axt.yaxis.set_minor_locator(MultipleLocator(0.1))
                axt.set_xlabel("Period (s)", fontdict=self.font_dict)

                axt.set_xscale("log", nonpositive="clip")
                if self.tipper_limits is None:
                    tmax = max([tyr.max(), tyi.max()])
                    if tmax > 1:
                        tmax = 0.899
                    tmin = min([tyr.min(), tyi.min()])
                    if tmin < -1:
                        tmin = -0.899
                    self.tipper_limits = (tmin - 0.1, tmax + 0.1)
                axt.set_ylim(self.tipper_limits)
                axt.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )

                tklabels = []
                xticks = []
                for tk in axt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk])
                        xticks.append(tk)
                    except KeyError:
                        pass
                axt.set_xticks(xticks)
                axt.set_xticklabels(tklabels, fontdict={"size": self.font_size})

                if pdict["tip"] != nrows - 1:
                    plt.setp(axt.get_yticklabels(), visible=False)
                # need to reset the xlimits caouse they get reset when calling
                # set_ticks for some reason
                axt.set_xlim(np.log10(self.xlimits[0]), np.log10(self.xlimits[1]))
            # ------plot strike angles----------------------------------------------
            if self._plot_strike.find("y") == 0:

                if self._plot_strike.find("i") > 0:
                    # strike from invariants
                    zinv = Zinvariants(mt.Z)
                    s1 = zinv.strike

                    # fold angles so go from -90 to 90
                    s1[np.where(s1 > 90)] -= -180
                    s1[np.where(s1 < -90)] += 180

                    # plot strike with error bars
                    ps1 = axst.errorbar(
                        mt.period,
                        s1,
                        marker=self.strike_inv_marker,
                        ms=self.marker_size,
                        mfc=self.strike_inv_color,
                        mec=self.strike_inv_color,
                        mew=self.marker_lw,
                        ls="none",
                        yerr=zinv.strike_err,
                        ecolor=self.strike_inv_color,
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    stlist.append(ps1[0])
                    stlabel.append("Z_inv")
                    st_maxlist.append(s1.max())
                    st_minlist.append(s1.min())
                if self._plot_strike.find("p") > 0:
                    # strike from phase tensor
                    pt = mt.pt  # type: PhaseTensor
                    s2, s2_err = pt.azimuth, pt.azimuth_err

                    # fold angles to go from -90 to 90
                    s2[np.where(s2 > 90)] -= 180
                    s2[np.where(s2 < -90)] += 180

                    # plot strike with error bars
                    ps2 = axst.errorbar(
                        mt.period,
                        s2,
                        marker=self.strike_pt_marker,
                        ms=self.marker_size,
                        mfc=self.strike_pt_color,
                        mec=self.strike_pt_color,
                        mew=self.marker_lw,
                        ls="none",
                        yerr=s2_err,
                        ecolor=self.strike_pt_color,
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    stlist.append(ps2[0])
                    stlabel.append("PT")
                    st_maxlist.append(s2.max())
                    st_minlist.append(s2.min())
                if self._plot_strike.find("t") > 0:
                    # strike from tipper
                    tp = mt.Tipper
                    s3 = tp.angle_real + 90

                    # fold to go from -90 to 90
                    s3[np.where(s3 > 90)] -= 180
                    s3[np.where(s3 < -90)] += 180

                    # plot strike with error bars
                    ps3 = axst.errorbar(
                        mt.period,
                        s3,
                        marker=self.strike_tip_marker,
                        ms=self.marker_size,
                        mfc=self.strike_tip_color,
                        mec=self.strike_tip_color,
                        mew=self.marker_lw,
                        ls="none",
                        yerr=np.zeros_like(s3),
                        ecolor=self.strike_tip_color,
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    stlist.append(ps3[0])
                    stlabel.append("Tip")
                    st_maxlist.append(s3.max())
                    st_minlist.append(s3.min())
                # --> set axes properties
                if self.strike_limits is None:
                    stmin = min(st_minlist)
                    if stmin - 3 < -90:
                        stmin -= 3
                    else:
                        stmin = -89.99
                    stmax = max(st_maxlist)
                    if stmin + 3 < 90:
                        stmin += 3
                    else:
                        stmin = 89.99
                    self.strike_limits = (
                        -max([abs(stmin), abs(stmax)]),
                        max([abs(stmin), abs(stmax)]),
                    )
                axst.plot(axr.get_xlim(), [0, 0], color="k", lw=0.5)
                if ii == 0:
                    axst.set_ylabel("Strike", fontdict=self.font_dict)
                else:
                    plt.setp(axst.get_yticklabels(), visible=False)
                axst.set_xlabel("Period (s)", fontdict=self.font_dict)
                axst.set_ylim(self.strike_limits)
                axst.yaxis.set_major_locator(MultipleLocator(30))
                axst.yaxis.set_minor_locator(MultipleLocator(5))
                axst.set_xscale("log", nonpositive="clip")
                axst.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )
                if ii == 0:
                    try:
                        axst.legend(
                            stlist,
                            stlabel,
                            loc=3,
                            markerscale=1,
                            borderaxespad=0.01,
                            labelspacing=0.07,
                            handletextpad=0.2,
                            borderpad=0.02,
                            prop={"size": self.font_size - 1},
                        )
                    except:
                        pass
                # set th xaxis tick labels to invisible
                if pdict["strike"] != nrows - 1:
                    plt.setp(axst.xaxis.get_ticklabels(), visible=False)
            # ------plot skew angle---------------------------------------------
            if self._plot_skew == "y":
                # strike from phase tensor
                pt = mt.pt
                sk, sk_err = pt.beta, pt.beta_err

                ps4 = axsk.errorbar(
                    mt.period,
                    sk,
                    marker=self.skew_marker,
                    ms=self.marker_size,
                    mfc=self.skew_color,
                    mec=self.skew_color,
                    mew=self.marker_lw,
                    ls="none",
                    yerr=sk_err,
                    ecolor=self.skew_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )
                stlist.append(ps4[0])
                stlabel.append("Skew")
                if self.skew_limits is None:
                    self.skew_limits = (-9, 9)
                axsk.set_ylim(self.skew_limits)
                axsk.yaxis.set_major_locator(MultipleLocator(3))
                axsk.yaxis.set_minor_locator(MultipleLocator(1))
                if ii == 0:
                    axsk.set_ylabel("Skew", fontdict)
                else:
                    plt.setp(axsk.get_yticklabels(), visible=False)
                axsk.set_xlabel("Period (s)", fontdict)
                axsk.set_xscale("log", nonpositive="clip")

                # set th xaxis tick labels to invisible
                if pdict["skew"] != nrows - 1:
                    plt.setp(axsk.xaxis.get_ticklabels(), visible=False)
            # ----plot phase tensor ellipse---------------------------------------
            if self._plot_pt == "y":
                # get phase tensor instance
                pt = mt.pt

                cmap = self.ellipse_cmap
                ckmin = self.ellipse_range[0]
                ckmax = self.ellipse_range[1]
                try:
                    ckstep = float(self.ellipse_range[2])
                except IndexError:
                    ckstep = 3
                if cmap == "mt_seg_bl2wh2rd":
                    bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                    nseg = float((ckmax - ckmin) / (2 * ckstep))
                # get the properties to color the ellipses by
                if (
                    self.ellipse_colorby == "phiminang"
                    or self.ellipse_colorby == "phimin"
                ):
                    colorarray = pt.phimin
                elif self.ellipse_colorby == "phidet":
                    colorarray = np.sqrt(abs(pt.det)) * (180 / np.pi)
                elif (
                    self.ellipse_colorby == "skew" or self.ellipse_colorby == "skew_seg"
                ):
                    colorarray = pt.beta
                elif self.ellipse_colorby == "ellipticity":
                    colorarray = pt.ellipticity
                elif self.ellipse_colorby in ["strike", "azimuth"]:
                    colorarray = self.fold_strike(pt.azimuth)
                    self.ellipse_range = (-90, 90)
                    ckmin = self.ellipse_range[0]
                    ckmax = self.ellipse_range[1]
                else:
                    raise NameError(self.ellipse_colorby + " is not supported")
                # -------------plot ellipses-----------------------------------
                for kk, ff in enumerate(mt.period):
                    # make sure the ellipses will be visable
                    eheight = pt.phimin[kk] / pt.phimax[kk] * self.ellipse_size
                    ewidth = pt.phimax[kk] / pt.phimax[kk] * self.ellipse_size

                    # create an ellipse scaled by phimin and phimax and
                    # oriented along the azimuth which is calculated as
                    # clockwise but needs to be plotted counter-clockwise
                    # hence the negative sign.
                    ellipd = patches.Ellipse(
                        (np.log10(ff) * self.ellipse_spacing, 0),
                        width=ewidth,
                        height=eheight,
                        angle=90 - pt.azimuth[kk],
                    )

                    axpt.add_patch(ellipd)

                    # get ellipse color
                    if cmap.find("seg") > 0:
                        ellipd.set_facecolor(
                            mtcl.get_plot_color(
                                colorarray[kk],
                                self.ellipse_colorby,
                                cmap,
                                ckmin,
                                ckmax,
                                bounds=bounds,
                            )
                        )
                    else:
                        ellipd.set_facecolor(
                            mtcl.get_plot_color(
                                colorarray[kk],
                                self.ellipse_colorby,
                                cmap,
                                ckmin,
                                ckmax,
                            )
                        )
                # ----set axes properties-----------------------------------------------
                # --> set tick labels and limits
                axpt.set_xlim(
                    np.floor(np.log10(self.xlimits[0])),
                    np.ceil(np.log10(self.xlimits[1])),
                )

                tklabels = []
                xticks = []
                for tk in axpt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk])
                        xticks.append(tk)
                    except KeyError:
                        pass
                axpt.set_xticks(xticks)
                axpt.set_xticklabels(tklabels, fontdict={"size": self.font_size})
                axpt.set_xlabel("Period (s)", fontdict=self.font_dict)
                axpt.set_ylim(
                    ymin=-1.5 * self.ellipse_size, ymax=1.5 * self.ellipse_size
                )

                axpt.grid(
                    True, alpha=0.25, which="major", color=(0.25, 0.25, 0.25), lw=0.25,
                )

                plt.setp(axpt.get_yticklabels(), visible=False)
                if pdict["pt"] != nrows - 1:
                    plt.setp(axpt.get_xticklabels(), visible=False)
                # add colorbar for PT only for first plot
                if ii == 0:
                    axpos = axpt.get_position()
                    cb_position = (
                        axpos.bounds[0] - 0.0575,
                        axpos.bounds[1] + 0.02,
                        0.01,
                        axpos.bounds[3] * 0.75,
                    )
                    cbax = self.fig.add_axes(cb_position)
                    if cmap == "mt_seg_bl2wh2rd":
                        # make a color list
                        clist = [
                            (cc, cc, 1)
                            for cc in np.arange(0, 1 + 1.0 / (nseg), 1.0 / (nseg))
                        ] + [
                            (1, cc, cc)
                            for cc in np.arange(1, -1.0 / (nseg), -1.0 / (nseg))
                        ]

                        # make segmented colormap
                        mt_seg_bl2wh2rd = colors.ListedColormap(clist)

                        # make bounds so that the middle is white
                        bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

                        # normalize the colors
                        norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

                        # make the colorbar
                        cbpt = mcb.ColorbarBase(
                            cbax,
                            cmap=mt_seg_bl2wh2rd,
                            norm=norms,
                            orientation="vertical",
                            ticks=bounds[1:-1],
                        )
                    else:
                        cbpt = mcb.ColorbarBase(
                            cbax,
                            cmap=mtcl.cmapdict[cmap],
                            norm=colors.Normalize(vmin=ckmin, vmax=ckmax),
                            orientation="vertical",
                        )
                    cbpt.set_ticks([ckmin, (ckmax - ckmin) / 2, ckmax])
                    cbpt.set_ticklabels(
                        [
                            "{0:.0f}".format(ckmin),
                            "{0:.0f}".format((ckmax - ckmin) / 2),
                            "{0:.0f}".format(ckmax),
                        ]
                    )
                    cbpt.ax.yaxis.set_label_position("left")
                    cbpt.ax.yaxis.set_label_coords(-1.05, 0.5)
                    cbpt.ax.yaxis.tick_right()
                    cbpt.ax.tick_params(axis="y", direction="in")
                    cbpt.set_label(
                        mtpl.ckdict[self.ellipse_colorby],
                        fontdict={"size": self.font_size},
                    )
            # ==  == Plot the Z_xx, Z_yy components if desired ==
            if self.plot_num == 2:
                # ---------plot the apparent resistivity----------------
                axr2 = self.fig.add_subplot(gs[0, 1], sharex=axr, sharey=axr)
                axr2.yaxis.set_label_coords(-0.1, 0.5)

                # res_xx
                ebxxr = axr2.errorbar(
                    mt.period,
                    mt.Z.res_xx,
                    marker=self.xy_marker,
                    ms=self.marker_size,
                    mfc=self.xy_mfc,
                    mec=self.xy_color,
                    mew=self.marker_lw,
                    ls=self.xy_ls,
                    yerr=mt.Z.res_err_xx,
                    ecolor=self.xy_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # res_yy
                ebyyr = axr2.errorbar(
                    mt.period,
                    mt.Z.res_yy,
                    marker=self.yx_marker,
                    ms=self.marker_size,
                    mfc=self.yx_mfc,
                    mec=self.yx_color,
                    mew=self.marker_lw,
                    ls=self.yx_ls,
                    yerr=mt.Z.res_err_yy,
                    ecolor=self.yx_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # --> set axes properties
                plt.setp(axr2.get_xticklabels(), visible=False)
                plt.setp(axr2.get_yticklabels(), visible=False)
                axr2.set_yscale("log", nonpositive="clip")
                axr2.set_xscale("log", nonpositive="clip")
                axr2.set_xlim(self.x_limits)
                axr2.set_ylim(self.res_limits)
                axr2.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )
                if ii == 0:
                    axr2.legend(
                        (ebxxr[0], ebyyr[0]),
                        ("$Z_{xx}$", "$Z_{yy}$"),
                        loc=3,
                        markerscale=1,
                        borderaxespad=0.01,
                        labelspacing=0.07,
                        handletextpad=0.2,
                        borderpad=0.02,
                    )
                # -----Plot the phase-----------------------------------
                axp2 = self.fig.add_subplot(gs[1, 1], sharex=axr, sharey=axp)
                axp2.yaxis.set_label_coords(-0.1, 0.5)

                # phase_xx
                ebxxp = axp2.errorbar(
                    mt.period,
                    mt.Z.phase_xx,
                    marker=self.xy_marker,
                    ms=self.marker_size,
                    mfc=self.xy_mfc,
                    mec=self.xy_color,
                    mew=self.marker_lw,
                    ls=self.xy_ls,
                    yerr=mt.Z.phase_err_xx,
                    ecolor=self.xy_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # phase_yy
                ebyyp = axp2.errorbar(
                    mt.period,
                    mt.Z.phase_yy,
                    marker=self.yx_marker,
                    ms=self.marker_size,
                    mfc=self.yx_mfc,
                    mec=self.yx_color,
                    mew=self.marker_lw,
                    ls=self.yx_ls,
                    yerr=mt.Z.phase_err_yy,
                    ecolor=self.yx_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # --> set axes properties
                plt.setp(axp2.get_xticklabels(), visible=False)
                plt.setp(axp2.get_yticklabels(), visible=False)
                axp2.set_xlabel("Period (s)", fontdict)
                axp2.set_xscale("log", nonpositive="clip")
                if self.phase_limits is None:
                    self.phase_limits = (-179.9, 179.9)
                axp2.set_ylim(self.phase_limits)
                axp2.set_xlim(self.x_limits)
                axp2.yaxis.set_major_locator(MultipleLocator(30))
                axp2.yaxis.set_minor_locator(MultipleLocator(5))
                # axp2.set_xticklabels(tklabels,
                #                      fontdict={'size': self.font_size})
                axp2.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )

                if len(list(pdict.keys())) > 2:
                    plt.setp(axp2.xaxis.get_ticklabels(), visible=False)
                    plt.setp(axp2.xaxis.get_label(), visible=False)
            # == =Plot the Determinant if desired ==  ==  ==  ==
            if self.plot_num == 3:

                # res_det
                ebdetr = axr.errorbar(
                    mt.period,
                    rp.res_det,
                    marker=self.det_marker,
                    ms=self.marker_size,
                    mfc=self.det_mfc,
                    mec=self.det_color,
                    mew=self.marker_lw,
                    ls=self.det_ls,
                    yerr=rp.res_det_err,
                    ecolor=self.det_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # phase_det
                ebdetp = axp.errorbar(
                    mt.period,
                    rp.phase_det,
                    marker=self.det_marker,
                    ms=self.marker_size,
                    mfc=self.det_mfc,
                    mec=self.det_color,
                    mew=self.marker_lw,
                    ls=self.det_ls,
                    yerr=rp.phase_det_err,
                    ecolor=self.det_color,
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # --> set axes properties
                plt.setp(axr.get_xticklabels(), visible=False)
                if ii == 0:
                    axr.set_ylabel(
                        "App. Res. ($\mathbf{\Omega \cdot m}$)",
                        fontdict=self.font_dict,
                    )
                else:
                    plt.setp(axr.get_yticklabels(), visible=False)
                axr.set_yscale("log", nonpositive="clip")
                axr.set_xscale("log", nonpositive="clip")
                axr.set_ylim(self.res_limits)
                axr.set_xlim(self.xlimits)
                axr.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )

                # --> set axes properties
                axp.set_xlabel("Period (s)", fontdict)

                if ii == 0:
                    axp.set_ylabel("Phase (deg)", fontdict)
                else:
                    plt.setp(axp.get_yticklabels(), visible=False)
                axp.set_xscale("log", nonpositive="clip")
                axp.set_ylim(self.phase_limits)
                axp.yaxis.set_major_locator(MultipleLocator(15))
                axp.yaxis.set_minor_locator(MultipleLocator(5))
                tklabels = [
                    mtpl.labeldict[tt]
                    for tt in np.arange(
                        np.log10(self.xlimits[0]), np.log10(self.xlimits[1]) + 1
                    )
                ]
                tklabels[0] = ""
                tklabels[-1] = ""

                axp.set_xticklabels(tklabels, fontdict={"size": self.font_size})
                axp.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25,
                )
            # make title and show

            axr.set_title(mt.station, fontsize=self.font_size, fontweight="bold")
        if show:
            plt.show()

    # ---plot the resistivity and phase
    def plot(self):
        """
        plot the apparent resistivity and phase
        """
        self._has_pt()
        self._has_tipper()
        self._set_subplot_params()

        # -----Plot All in one figure with each plot as a subfigure------------
        if self.plot_style == "all":
            pass
        # ===Plot all responses into one plot to compare changes ==
        if self.plot_style == "compare":
            ns = len(self.mt_list)

            # make color lists for the plots going light to dark
            cxy = [(0, 0 + float(cc) / ns, 1 - float(cc) / ns) for cc in range(ns)]
            cyx = [(1, float(cc) / ns, 0) for cc in range(ns)]
            cdet = [(0, 1 - float(cc) / ns, 0) for cc in range(ns)]
            ctipr = [
                (0.75 * cc / ns, 0.75 * cc / ns, 0.75 * cc / ns) for cc in range(ns)
            ]
            ctipi = [(float(cc) / ns, 1 - float(cc) / ns, 0.25) for cc in range(ns)]
            cst = [(0.5 * cc / ns, 0, 0.5 * cc / ns) for cc in range(ns)]

            # make marker lists for the different components
            mxy = ["s", "D", "x", "+", "*", "1", "3", "4"] * 5
            myx = ["o", "h", "8", "p", "H", 7, 4, 6] * 5

            legendlistxy = []
            legendlistyx = []
            stationlist = []
            tiplist = []
            stlist = []
            sklist = []

            # set some parameters of the figure and subplot spacing
            plt.rcParams["font.size"] = self.font_size
            plt.rcParams["figure.subplot.bottom"] = 0.1
            plt.rcParams["figure.subplot.top"] = 0.97
            plt.rcParams["figure.subplot.left"] = 0.08
            plt.rcParams["figure.subplot.right"] = 0.98

            # set the font properties for the axis labels
            fontdict = {"size": self.font_size + 1, "weight": "bold"}

            # set figure size according to what the plot will be.
            if self.fig_size is None:
                if self.plot_num == 1 or self.plot_num == 3:
                    self.fig_size = [8, 6]
                    pass
                elif self.plot_num == 2:
                    self.fig_size = [8, 6]
                    nrows += 1
            # make a figure instance
            self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

            # make a grid as usual, but put xy and yx in different plots
            # otherwise the plot is too busy to see what's going on.
            hr = [2, 1.5] + [1] * (nrows - 2)
            gs = gridspec.GridSpec(nrows, 2, height_ratios=hr, hspace=0.05)

            # --> make figure for xy,yx components
            if self.plot_num == 1 or self.plot_num == 3:
                # set label coordinates
                labelcoords = (-0.125, 0.5)

                # space out the subplots
                gs.update(hspace=0.05, wspace=0.02, left=0.1)
            # --> make figure for all 4 components
            elif self.plot_num == 2:
                # set label coordinates
                labelcoords = (-0.125, 0.5)

                # space out the subplots
                gs.update(hspace=0.05, wspace=0.02, left=0.07)

                for key in pdict:
                    if key != "res" and key != "phase":
                        pdict[key] += 1
            # --> create the axes instances
            # apparent resistivity axis
            self.axrxy = self.fig.add_subplot(gs[0, 0])
            self.axryx = self.fig.add_subplot(
                gs[0, 1], sharex=self.axrxy, sharey=self.axrxy
            )

            # phase axis that shares period axis with resistivity
            self.axpxy = self.fig.add_subplot(gs[1, 0], sharex=self.axrxy)
            self.axpyx = self.fig.add_subplot(
                gs[1, 1], sharex=self.axrxy, sharey=self.axpxy
            )

            # place y coordinate labels in the same location
            # self.axrxy.yaxis.set_label_coords(labelcoords[0], labelcoords[1])
            # self.axpxy.yaxis.set_label_coords(labelcoords[0], labelcoords[1])

            # --> plot tipper
            try:
                self.axt = self.fig.add_subplot(gs[pdict["tip"], :])
                # self.axt.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                 labelcoords[1])
            except KeyError:
                pass
            # --> plot phase tensors
            try:
                # can't share axis because not on the same scale
                self.axpt = self.fig.add_subplot(gs[pdict["pt"], :], aspect="equal")
                # self.axpt.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass
            # --> plot strike
            try:
                self.axst = self.fig.add_subplot(
                    gs[pdict["strike"], :], sharex=self.axrxy
                )
                # self.axst.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass
            # --> plot skew
            try:
                self.axsk = self.fig.add_subplot(
                    gs[pdict["skew"], :], sharex=self.axrxy
                )
                # self.axsk.yaxis.set_label_coords(labelcoords[0] * .5,
                #                                  labelcoords[1])
            except KeyError:
                pass
            for ii, mt in enumerate(self.mt_list):
                # get the reistivity and phase object

                # set x-axis limits from short period to long period
                if self.xlimits is None:
                    self.xlimits = (
                        10 ** (np.floor(np.log10(mt.period.min()))),
                        10 ** (np.ceil(np.log10(mt.period.max()))),
                    )
                else:
                    self.xlimits = (
                        10
                        ** min(
                            [
                                np.floor(np.log10(self.xlimits[0])),
                                np.floor(np.log10(mt.period.min())),
                            ]
                        ),
                        10
                        ** max(
                            [
                                np.ceil(np.log10(self.xlimits[1])),
                                np.ceil(np.log10(mt.period.max())),
                            ]
                        ),
                    )
                if self.phase_limits is None:
                    self.phase_limits = (0, 89.9)
                stationlist.append(mt.station)

                # ==  ==  ==  == =Plot Z_xy and Z_yx ==
                # if self.plot_num == 1 or self.plot_num == 2:
                # ---------plot the apparent resistivity--------------------
                # --> plot as error bars and just as points xy-blue, yx-red
                # res_xy
                ebxyr = self.axrxy.errorbar(
                    mt.period,
                    mt.Z.res_xy,
                    color=cxy[ii],
                    marker=mxy[ii % len(mxy)],
                    ms=self.marker_size,
                    mfc="None",
                    mec=cxy[ii],
                    mew=self.marker_lw,
                    ls=self.xy_ls,
                    yerr=mt.Z.res_err_xy,
                    ecolor=cxy[ii],
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # res_yx
                ebyxr = self.axryx.errorbar(
                    mt.period,
                    mt.Z.res_yx,
                    color=cyx[ii],
                    marker=myx[ii % len(myx)],
                    ms=self.marker_size,
                    mfc="None",
                    mec=cyx[ii],
                    mew=self.marker_lw,
                    ls=self.yx_ls,
                    yerr=mt.Z.res_err_yx,
                    ecolor=cyx[ii],
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # -----Plot the phase---------------------------------------
                # phase_xy
                self.axpxy.errorbar(
                    mt.period,
                    mt.Z.phase_xy,
                    color=cxy[ii],
                    marker=mxy[ii % len(mxy)],
                    ms=self.marker_size,
                    mfc="None",
                    mec=cxy[ii],
                    mew=self.marker_lw,
                    ls=self.xy_ls,
                    yerr=mt.Z.phase_err_xy,
                    ecolor=cxy[ii],
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                # phase_yx: Note add 180 to place it in same quadrant as
                # phase_xy
                self.axpyx.errorbar(
                    mt.period,
                    mt.Z.phase_yx + 180,
                    color=cyx[ii],
                    marker=myx[ii % len(myx)],
                    ms=self.marker_size,
                    mfc="None",
                    mec=cyx[ii],
                    mew=self.marker_lw,
                    ls=self.yx_ls,
                    yerr=mt.Z.phase_err_yx,
                    ecolor=cyx[ii],
                    capsize=self.marker_size,
                    elinewidth=self.marker_lw,
                )

                legendlistxy.append(ebxyr)
                legendlistyx.append(ebyxr)

                # ==== Plot the Z_xx, Z_yy components if desired ==
                if self.plot_num == 2:
                    # ---------plot the apparent resistivity----------------
                    self.axr2xx = self.fig.add_subplot(gs[2, 0], sharex=self.axrxy)
                    self.axr2xx.yaxis.set_label_coords(-0.095, 0.5)
                    self.axr2yy = self.fig.add_subplot(gs[2, 1], sharex=self.axrxy)

                    # res_xx
                    ebxxr = self.axr2xx.errorbar(
                        mt.period,
                        mt.Z.res_xx,
                        color=cxy[ii],
                        marker=mxy[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cxy[ii],
                        mew=self.marker_lw,
                        ls=self.xy_ls,
                        yerr=mt.Z.res_err_xx,
                        ecolor=cxy[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    # res_yy
                    ebyyr = self.axr2yy.errorbar(
                        mt.period,
                        mt.Z.res_yy,
                        color=cyx[ii],
                        marker=myx[ii % len(myx)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cyx[ii],
                        mew=self.marker_lw,
                        ls=self.yx_ls,
                        yerr=mt.Z.res_err_yy,
                        ecolor=cyx[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    # -----Plot the phase-----------------------------------
                    self.axp2xx = self.fig.add_subplot(gs[2, 0], sharex=self.axrxy)
                    self.axp2xx.yaxis.set_label_coords(-0.095, 0.5)
                    self.axp2yy = self.fig.add_subplot(gs[2, 1], sharex=self.axrxy)

                    # phase_xx
                    ebxxp = self.axp2xx.errorbar(
                        mt.period,
                        mt.Z.phase_xx,
                        color=cxy[ii],
                        marker=mxy[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cxy[ii],
                        mew=self.marker_lw,
                        ls=self.xy_ls,
                        yerr=mt.Z.phase_err_xx,
                        ecolor=cxy[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    # phase_yy
                    ebyyp = self.axp2yy.errorbar(
                        mt.period,
                        mt.Z.phase_yy,
                        color=cyx[ii],
                        marker=myx[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cyx[ii],
                        mew=self.marker_lw,
                        ls=self.yx_ls,
                        yerr=mt.Z.phase_err_yy,
                        ecolor=cyx[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )
                # ===Plot the Determinant if desired ==
                if self.plot_num == 3:
                    # res_det
                    ebdetr = self.axrxy.errorbar(
                        mt.period,
                        mt.Z.res_det,
                        color=cxy[ii],
                        marker=mxy[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cdet[ii],
                        mew=self.marker_lw,
                        ls=self.det_ls,
                        yerr=mt.Z.res_det_err,
                        ecolor=cdet[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    # phase_det
                    ebdetp = self.axpxy.errorbar(
                        mt.period,
                        mt.Z.phase_det,
                        color=cyx[ii],
                        marker=mxy[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc="None",
                        mec=cdet[ii],
                        mew=self.marker_lw,
                        ls=self.det_ls,
                        yerr=mt.Z.phase_det_err,
                        ecolor=cdet[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )

                    legendlistxy.append(ebdetr)
                # -----plot tipper----------------------------------------------
                if self._plot_tipper.find("y") == 0:

                    txr = mt.Tipper.mag_real * np.sin(
                        mt.Tipper.angle_real * np.pi / 180
                        + np.pi * self.arrow_direction
                    )
                    tyr = mt.Tipper.mag_real * np.cos(
                        mt.Tipper.angle_real * np.pi / 180
                        + np.pi * self.arrow_direction
                    )

                    txi = mt.Tipper.mag_imag * np.sin(
                        mt.Tipper.angle_imag * np.pi / 180
                        + np.pi * self.arrow_direction
                    )
                    tyi = mt.Tipper.mag_imag * np.cos(
                        mt.Tipper.angle_imag * np.pi / 180
                        + np.pi * self.arrow_direction
                    )

                    nt = len(txr)

                    for aa in range(nt):
                        xlenr = txr[aa] * np.log10(mt.period[aa])
                        xleni = txi[aa] * np.log10(mt.period[aa])

                        if self.tipper_limits is None:
                            tmax = max([tyr.max(), tyi.max()])
                            tmin = min([tyr.min(), tyi.min()])
                            if np.isnan(tmax):
                                tmax = 1.0
                            if np.isnan(tmin):
                                tmin = -1.0
                            tmin = max([-1, tmin])
                            tmax = min([1, tmax])
                            self.tipper_limits = (tmin - 0.1, tmax + 0.1)
                        else:
                            tmax = (
                                max([tyr.max(), tyi.max(), self.tipper_limits[1] - 0.1])
                                + 0.1
                            )
                            tmin = (
                                min([tyr.min(), tyi.min(), self.tipper_limits[0] + 0.1])
                                - 0.1
                            )
                            if np.isnan(tmax):
                                tmax = 1.0
                            if np.isnan(tmin):
                                tmin = -1.0
                            tmin = max([-1, tmin])
                            tmax = min([1, tmax])
                            self.tipper_limits = (tmin, tmax)
                        # --> plot real arrows
                        if self._plot_tipper.find("r") > 0:
                            self.axt.arrow(
                                np.log10(mt.period[aa]),
                                0,
                                xlenr,
                                tyr[aa],
                                lw=self.arrow_lw,
                                facecolor=ctipr[ii],
                                edgecolor=ctipr[ii],
                                head_width=self.arrow_head_width,
                                head_length=self.arrow_head_length,
                                length_includes_head=False,
                            )
                        # --> plot imaginary arrows
                        if self._plot_tipper.find("i") > 0:
                            self.axt.arrow(
                                np.log10(mt.period[aa]),
                                0,
                                xleni,
                                tyi[aa],
                                lw=self.arrow_lw,
                                head_width=self.arrow_head_width,
                                head_length=self.arrow_head_length,
                                length_includes_head=False,
                            )
                    lt = self.axt.plot(0, 0, lw=1, color=ctipr[ii])
                    tiplist.append(lt[0])
                # ------plot strike angles----------------------------------------------
                if self._plot_strike.find("y") == 0:
                    if self._plot_strike.find("p") > 0:
                        # strike from phase tensor
                        s2 = mt.pt.azimuth
                        s2_err = mt.pt.azimuth_err

                        # fold angles to go from -90 to 90
                        s2[np.where(s2 > 90)] -= 180
                        s2[np.where(s2 < -90)] += 180

                        # plot strike with error bars
                        ps2 = self.axst.errorbar(
                            mt.period,
                            s2,
                            marker=myx[ii % len(myx)],
                            ms=self.marker_size,
                            mfc=cxy[ii],
                            mec=cxy[ii],
                            mew=self.marker_lw,
                            ls="none",
                            yerr=s2_err,
                            ecolor=cxy[ii],
                            capsize=self.marker_size,
                            elinewidth=self.marker_lw,
                        )

                        stlist.append(ps2[0])
                    if self._plot_strike.find("t") > 0:
                        # strike from tipper
                        s3 = mt.Tipper.angle_real + 90

                        # fold to go from -90 to 90
                        s3[np.where(s3 > 90)] -= 180
                        s3[np.where(s3 < -90)] += 180

                        # plot strike with error bars
                        ps3 = self.axst.errorbar(
                            mt.period,
                            s3,
                            marker=mxy[ii % len(mxy)],
                            ms=self.marker_size,
                            mfc=ctipr[ii],
                            mec=ctipr[ii],
                            mew=self.marker_lw,
                            ls="none",
                            yerr=np.zeros_like(s3),
                            ecolor=ctipr[ii],
                            capsize=self.marker_size,
                            elinewidth=self.marker_lw,
                        )

                        stlist.append(ps3[0])
                # ------plot skew angle---------------------------------------------
                if self._plot_skew == "y":
                    # strike from phase tensor
                    sk = mt.pt.beta
                    sk_err = mt.pt.beta_err

                    ps4 = self.axsk.errorbar(
                        mt.period,
                        sk,
                        marker=mxy[ii % len(mxy)],
                        ms=self.marker_size,
                        mfc=cxy[ii],
                        mec=cxy[ii],
                        mew=self.marker_lw,
                        ls="none",
                        yerr=sk_err,
                        ecolor=cxy[ii],
                        capsize=self.marker_size,
                        elinewidth=self.marker_lw,
                    )
                    stlist.append(ps4[0])
                # ----plot phase tensor ellipse---------------------------------------
                if self._plot_pt == "y":
                    # get phase tensor instance
                    pt = mt.pt

                    cmap = self.ellipse_cmap
                    ckmin = self.ellipse_range[0]
                    ckmax = self.ellipse_range[1]
                    try:
                        ckstep = float(self.ellipse_range[2])
                    except IndexError:
                        ckstep = 3
                    if cmap == "mt_seg_bl2wh2rd":
                        bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                        nseg = float((ckmax - ckmin) / (2 * ckstep))
                    # get the properties to color the ellipses by
                    if (
                        self.ellipse_colorby == "phiminang"
                        or self.ellipse_colorby == "phimin"
                    ):
                        colorarray = mt.pt.phimin
                    elif self.ellipse_colorby == "phidet":
                        colorarray = np.sqrt(abs(mt.pt.det)) * (180 / np.pi)
                    elif (
                        self.ellipse_colorby == "skew"
                        or self.ellipse_colorby == "skew_seg"
                    ):
                        colorarray = mt.pt.beta
                    elif self.ellipse_colorby == "ellipticity":
                        colorarray = mt.pt.ellipticity
                    else:
                        raise NameError(self.ellipse_colorby + " is not supported")
                    # -------------plot ellipses-----------------------------------
                    for kk, ff in enumerate(mt.period):
                        # make sure the ellipses will be visable
                        eheight = (
                            mt.pt.phimin[kk] / mt.pt.phimax[kk] * self.ellipse_size
                        )
                        ewidth = mt.pt.phimax[kk] / mt.pt.phimax[kk] * self.ellipse_size

                        # create an ellipse scaled by phimin and phimax and oriented
                        # along the azimuth which is calculated as clockwise but needs
                        # to be plotted counter-clockwise hence the negative sign.
                        ellipd = patches.Ellipse(
                            (
                                np.log10(ff) * self.ellipse_spacing,
                                ii * self.ellipse_size * 1.5,
                            ),
                            width=ewidth,
                            height=eheight,
                            angle=90 - pt.azimuth[kk],
                        )

                        self.axpt.add_patch(ellipd)

                        # get ellipse color
                        if cmap.find("seg") > 0:
                            ellipd.set_facecolor(
                                mtcl.get_plot_color(
                                    colorarray[kk],
                                    self.ellipse_colorby,
                                    cmap,
                                    ckmin,
                                    ckmax,
                                    bounds=bounds,
                                )
                            )
                        else:
                            ellipd.set_facecolor(
                                mtcl.get_plot_color(
                                    colorarray[kk],
                                    self.ellipse_colorby,
                                    cmap,
                                    ckmin,
                                    ckmax,
                                )
                            )
                        ellipd.set_edgecolor(cxy[ii])
            # -------set axis properties---------------------------------------
            self.axrxy.set_yscale("log", nonpositive="clip")
            self.axrxy.set_xscale("log", nonpositive="clip")
            self.axrxy.set_ylim(self.res_limits)
            self.axrxy.set_xlim(self.x_limits)
            self.axrxy.grid(
                True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
            )

            # make a text label in upper left hand corner
            # label the plot with a text box
            if self.text_location is None:
                txloc = self.xlimits[0] * self.text_xpad
                tyloc = self.axrxy.get_ylim()[1] * self.text_ypad
            else:
                txloc = self.text_location[0]
                tyloc = self.text_location[1]
            self.text = self.axrxy.text(
                txloc,
                tyloc,
                "$Z_{xy}$",
                fontdict={"size": self.text_size, "weight": self.text_weight},
                verticalalignment="top",
                horizontalalignment="left",
                bbox={"facecolor": "white", "alpha": 1},
            )

            plt.setp(self.axrxy.get_xticklabels(), visible=False)

            self.axrxy.set_ylabel(
                "App. Resistivity($\Omega \cdot$m)", fontdict=self.font_dict
            )

            self.axryx.set_yscale("log", nonpositive="clip")
            self.axryx.set_xscale("log", nonpositive="clip")
            self.axryx.set_ylim(self.res_limits)
            self.axryx.set_xlim(self.x_limits)
            self.axryx.grid(
                True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
            )

            self.text = self.axryx.text(
                txloc,
                tyloc,
                "$Z_{yx}$",
                fontdict={"size": self.text_size, "weight": self.text_weight},
                verticalalignment="top",
                horizontalalignment="left",
                bbox={"facecolor": "white", "alpha": 1},
            )

            plt.setp(self.axryx.get_xticklabels(), visible=False)
            plt.setp(self.axryx.get_yticklabels(), visible=False)

            # check the phase to see if any point are outside of [0:90]
            if self.phase_limits is None:
                self.phase_limits = (0, 89.99)
            # --> set axes properties
            self.axpxy.set_xlabel("Period(s)", fontdict=self.font_dict)
            self.axpxy.set_ylabel("Phase(deg)", fontdict=self.font_dict)
            self.axpxy.set_xscale("log", nonpositive="clip")
            self.axpxy.set_ylim(self.phase_limits)
            self.axpxy.yaxis.set_major_locator(MultipleLocator(15))
            self.axpxy.yaxis.set_minor_locator(MultipleLocator(5))
            self.axpxy.grid(
                True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
            )
            if len(list(pdict.keys())) > 2:
                plt.setp(self.axpxy.xaxis.get_ticklabels(), visible=False)
                self.axpxy.set_xlabel("")
            self.axpyx.set_xlabel("Period(s)", fontdict=self.font_dict)
            self.axpyx.set_xscale("log", nonpositive="clip")
            self.axpyx.set_ylim(self.phase_limits)
            self.axpyx.yaxis.set_major_locator(MultipleLocator(15))
            self.axpyx.yaxis.set_minor_locator(MultipleLocator(5))
            self.axpyx.grid(
                True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
            )
            plt.setp(self.axpyx.yaxis.get_ticklabels(), visible=False)

            if len(list(pdict.keys())) > 2:
                plt.setp(self.axpyx.xaxis.get_ticklabels(), visible=False)
                self.axpyx.set_xlabel("")
            # make legend
            if self.plot_num == 1:
                self.axrxy.legend(
                    legendlistxy,
                    stationlist,
                    loc=3,
                    ncol=2,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )

                self.axryx.legend(
                    legendlistyx,
                    stationlist,
                    loc=3,
                    ncol=2,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )
            elif self.plot_num == 3:
                llist = [ll[0] for ll in legendlistxy]
                slist = [ss + "_det" for ss in stationlist]

                self.axrxy.legend(
                    llist,
                    slist,
                    loc=3,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )
                self.axryx.legend(
                    llist,
                    slist,
                    loc=3,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )
            if self.plot_num == 2:
                # --> set axes properties for resxx
                self.axrxy.set_yscale("log", nonpositive="clip")
                self.axrxy.set_xscale("log", nonpositive="clip")
                self.axrxy.set_xlim(self.x_limits)
                self.axrxy.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )
                plt.setp(self.axrxy.get_xticklabels(), visible=False)

                # --> set axes properties for resyy
                self.axryx.set_yscale("log", nonpositive="clip")
                self.axryx.set_xscale("log", nonpositive="clip")
                self.axryx.set_xlim(self.x_limits)
                self.axryx.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )
                plt.setp(self.axryx.get_xticklabels(), visible=False)

                # --> set axes properties Phasexx
                self.axpxy.set_xlabel("Period(s)", fontdict)
                self.axpxy.set_xscale("log", nonpositive="clip")
                self.axpxy.set_ylim(ymin=-179.9, ymax=179.9)
                self.axpxy.yaxis.set_major_locator(MultipleLocator(30))
                self.axpxy.yaxis.set_minor_locator(MultipleLocator(5))
                self.axpxy.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )

                # --> set axes properties Phaseyy
                self.axpyx.set_xlabel("Period(s)", fontdict)
                self.axpyx.set_xscale("log", nonpositive="clip")
                self.axpyx.set_ylim(ymin=-179.9, ymax=179.9)
                self.axpyx.yaxis.set_major_locator(MultipleLocator(30))
                self.axpyx.yaxis.set_minor_locator(MultipleLocator(5))
                self.axpyx.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )
                if len(list(pdict.keys())) > 3:
                    plt.setp(self.axpxy.xaxis.get_ticklabels(), visible=False)
                    self.axpxy.set_xlabel("")
                    plt.setp(self.axpyx.xaxis.get_ticklabels(), visible=False)
                    self.axpyx.set_xlabel("")
            if self._plot_tipper.find("y") == 0:
                self.axt.plot(self.axt.get_xlim(), [0, 0], color="k", lw=0.5)
                # --> set axis properties Tipper
                if self.plot_num == 2:
                    plt.setp(self.axpxy.get_xticklabels(), visible=False)
                    self.axpxy.set_xlabel("")
                    plt.setp(self.axpyx.get_xticklabels(), visible=False)
                    self.axpyx.set_xlabel("")
                self.axt.yaxis.set_major_locator(MultipleLocator(0.2))
                self.axt.yaxis.set_minor_locator(MultipleLocator(0.1))
                self.axt.set_xlabel("Period(s)", fontdict=self.font_dict)
                self.axt.set_ylabel("Tipper", fontdict=self.font_dict)
                self.axt.set_xlim(np.log10(self.xlimits[0]), np.log10(self.xlimits[1]))
                tklabels = []
                xticks = []
                for tk in self.axt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk])
                        xticks.append(tk)
                    except KeyError:
                        pass
                self.axt.set_xticks(xticks)
                self.axt.set_xticklabels(tklabels, fontdict={"size": self.font_size})

                self.axt.set_ylim(self.tipper_limits)
                self.axt.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )

                self.axt.legend(
                    tiplist,
                    stationlist,
                    loc=3,
                    ncol=2,
                    markerscale=1,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.02,
                )

                # need to reset the xlimits caouse they get reset when calling
                # set_ticks for some reason
                self.axt.set_xlim(np.log10(self.xlimits[0]), np.log10(self.xlimits[1]))

                if pdict["tip"] != nrows - 1:
                    plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
                    self.axt.set_xlabel(" ")
            # --> set axes properties for strike and skew
            if self._plot_strike[0] == "y":

                if self.strike_limits is None:
                    self.strike_limits = (-89.99, 89.99)
                self.axst.plot(self.axrxy.get_xlim(), [0, 0], color="k", lw=0.5)

                self.axst.set_ylabel("Strike(deg)", fontdict=self.font_dict)
                self.axst.set_xlabel("Period(s)", fontdict=self.font_dict)
                self.axst.set_ylim(self.strike_limits)
                self.axst.yaxis.set_major_locator(MultipleLocator(30))
                self.axst.yaxis.set_minor_locator(MultipleLocator(5))
                self.axst.set_xscale("log", nonpositive="clip")
                self.axst.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )
                # self.axst.legend(stlist,
                #                stationlist,
                #                loc=3,
                #                ncol=2,
                #                markerscale=1,
                #                borderaxespad=.01,
                #                labelspacing=.07,
                #                handletextpad=.2,
                #                borderpad=.02)
                if pdict["strike"] != nrows - 1:
                    plt.setp(self.axst.xaxis.get_ticklabels(), visible=False)
                    self.axst.set_xlabel(" ")
            # --> set axes properties for skew
            if self._plot_skew == "y":
                self.axsk.set_ylim(self.skew_limits)
                self.axsk.yaxis.set_major_locator(MultipleLocator(3))
                self.axsk.yaxis.set_minor_locator(MultipleLocator(1))
                self.axsk.set_ylabel("Skew(deg)", fontdict)
                self.axsk.set_xlabel("Period(s)", fontdict=self.font_dict)
                self.axsk.set_xscale("log", nonpositive="clip")
                self.axsk.grid(
                    True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
                )

                # self.axsk.legend(sklist,
                #                 stationlist,
                #                 loc=4,
                #                 ncol=2,
                #                 markerscale=1,
                #                 borderaxespad=.01,
                #                 labelspacing=.07,
                #                 handletextpad=.2,
                #                 borderpad=.02)
                if pdict["skew"] != nrows - 1:
                    plt.setp(self.axsk.xaxis.get_ticklabels(), visible=False)
                    self.axsk.set_xlabel(" ")
            # ----set axes properties for pt-----------------------------------
            if self._plot_pt == "y":
                self.axpt.set_xlim(
                    np.floor(np.log10(self.xlimits[0])) * self.ellipse_spacing,
                    np.ceil(np.log10(self.xlimits[1])) * self.ellipse_spacing,
                )

                tklabels = []
                xticks = []
                for tk in self.axpt.get_xticks():
                    try:
                        tklabels.append(mtpl.labeldict[tk / self.ellipse_spacing])
                        xticks.append(tk)
                    except KeyError:
                        pass
                self.axpt.set_xticks(xticks)
                self.axpt.set_xticklabels(tklabels, fontdict={"size": self.font_size})
                self.axpt.set_xlabel("Period (s)", fontdict=self.font_dict)
                self.axpt.set_ylim(
                    ymin=-1.5 * self.ellipse_size,
                    ymax=1.5 * self.ellipse_size * (ii + 1),
                )

                self.axpt.grid(
                    True, alpha=0.25, which="major", color=(0.25, 0.25, 0.25), lw=0.25
                )

                plt.setp(self.axpt.get_yticklabels(), visible=False)
                if pdict["pt"] != nrows - 1:
                    plt.setp(self.axpt.get_xticklabels(), visible=False)
                # add colorbar for PT
                axpos = self.axpt.get_position()
                cb_position = (
                    axpos.bounds[0] - 0.0575,
                    axpos.bounds[1] + 0.02,
                    0.01,
                    axpos.bounds[3] * 0.75,
                )
                self.cbax = self.fig.add_axes(cb_position)
                if self.ellipse_cmap == "mt_seg_bl2wh2rd":
                    # make a color list
                    clist = [
                        (cc, cc, 1)
                        for cc in np.arange(0, 1 + 1.0 / (nseg), 1.0 / (nseg))
                    ] + [
                        (1, cc, cc) for cc in np.arange(1, -1.0 / (nseg), -1.0 / (nseg))
                    ]

                    # make segmented colormap
                    mt_seg_bl2wh2rd = colors.ListedColormap(clist)

                    # make bounds so that the middle is white
                    bounds = np.arange(ckmin - ckstep, ckmax + 2 * ckstep, ckstep)

                    # normalize the colors
                    norms = colors.BoundaryNorm(bounds, mt_seg_bl2wh2rd.N)

                    # make the colorbar
                    self.cbpt = mcb.ColorbarBase(
                        self.cbax,
                        cmap=mt_seg_bl2wh2rd,
                        norm=norms,
                        orientation="vertical",
                        ticks=bounds[1:-1],
                    )
                else:
                    self.cbpt = mcb.ColorbarBase(
                        self.cbax,
                        cmap=mtcl.cmapdict[cmap],
                        norm=colors.Normalize(vmin=ckmin, vmax=ckmax),
                        orientation="vertical",
                    )
                self.cbpt.set_ticks([ckmin, (ckmax - ckmin) / 2, ckmax])
                self.cbpt.set_ticklabels(
                    [
                        "{0:.0f}".format(ckmin),
                        "{0:.0f}".format((ckmax - ckmin) / 2),
                        "{0:.0f}".format(ckmax),
                    ]
                )
                self.cbpt.ax.yaxis.set_label_position("left")
                self.cbpt.ax.yaxis.set_label_coords(-1.05, 0.5)
                self.cbpt.ax.yaxis.tick_right()
                self.cbpt.ax.tick_params(axis="y", direction="in")
                self.cbpt.set_label(
                    mtpl.ckdict[self.ellipse_colorby], fontdict={"size": self.font_size}
                )

                if pdict["pt"] != nrows - 1:
                    plt.setp(self.axpt.xaxis.get_ticklabels(), visible=False)
                    self.axpt.set_xlabel(" ")
            if show:
                plt.show()

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

        plt.close("all")
        self.plot()

    def __str__(self):
        """
        rewrite the string builtin to give a useful message
        """

        return (
            "Plots resistivity and phase for the different modes of the MT \n"
            + "response for multiple sites. At the moment it supports the \n"
            + "input of an .edi file. Other formats that will be supported\n"
            + "are the impedance tensor and errors with an array of periods\n"
            + "and .j format.\n"
        )
