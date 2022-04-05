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
        self.tf_list = tf_list

        super().__init__(**kwargs)

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

    def _plot_resistivity(self, axr, period, z_obj, mode="od", index=0, axr2=None):

        if mode == "od":
            comps = ["xy", "yx"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
            if axr2 is not None:
                ax_list = [axr, axr2]
            else:
                ax_list = [axr, axr]
        elif mode == "d":
            comps = ["xx", "yy"]
            props = [self.xy_error_bar_properties, self.yx_error_bar_properties]
            if axr2 is not None:
                ax_list = [axr, axr2]
            else:
                ax_list = [axr, axr]
        elif mode == "det":
            comps = ["xy", "yx", "det"]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
                self.det_error_bar_properties,
            ]
            if axr2 is not None:
                ax_list = [axr, axr2, axr]
            else:
                ax_list = [axr, axr, axr]
        res_limits = self.set_resistivity_limits(z_obj.resistivity, mode=mode)
        x_limits = self.set_period_limits(period)

        eb_list = []
        label_list = []
        for comp, prop, ax in zip(comps, props, ax_list):
            ebax = plot_resistivity(
                ax,
                period,
                getattr(z_obj, f"res_{comp}"),
                getattr(z_obj, f"res_err_{comp}"),
                **prop,
            )
            eb_list.append(ebax[0])
            label_list.append(f"$Z_{comp}$")
            # --> set axes properties
            plt.setp(ax.get_xticklabels(), visible=False)

            ax.set_yscale("log", nonpositive="clip")
            ax.set_xscale("log", nonpositive="clip")
            ax.set_xlim(x_limits)
            ax.set_ylim(res_limits)
            ax.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)
        if index == 0:
            axr.set_ylabel(
                "App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=self.font_dict
            )
        else:
            plt.setp(axr.get_yticklabels(), visible=False)
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

        return eb_list, label_list

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

        for comp, prop in zip(comps, props):
            if comp == "yx":
                plot_phase(
                    axp,
                    period,
                    getattr(z_obj, f"phase_{comp}") + 180,
                    getattr(z_obj, f"phase_err_{comp}"),
                    **prop,
                )
            else:
                plot_phase(
                    axp,
                    period,
                    getattr(z_obj, f"phase_{comp}"),
                    getattr(z_obj, f"phase_err_{comp}"),
                    **prop,
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

    def _plot_pt(self, axpt, period, pt_obj, index=0):
        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt:

            color_array = self.get_pt_color_array(pt_obj)
            x_limits = self.set_period_limits(period)

            # -------------plot ellipses-----------------------------------
            self.cbax, self.cbpt, = plot_pt_lateral(
                axpt, pt_obj, color_array, self.ellipse_properties, self.fig,
            )

            # ----set axes properties-----------------------------------------------
            # --> set tick labels and limits
            axpt.set_xlim(np.log10(x_limits[0]), np.log10(x_limits[1]))

            tklabels, xticks = get_log_tick_labels(axpt)

            axpt.set_xticks(xticks)
            axpt.set_xticklabels(tklabels, fontdict={"size": self.font_size})
            axpt.set_xlabel("Period (s)", fontdict=self.font_dict)

            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            axpt.set_xlim(np.log10(x_limits[0]), np.log10(x_limits[1]))
            axpt.grid(
                True, alpha=0.25, which="major", color=(0.25, 0.25, 0.25), lw=0.25
            )

            plt.setp(axpt.get_yticklabels(), visible=False)

            if index == 0:
                self.cbpt.set_label(
                    self.cb_label_dict[self.ellipse_colorby],
                    fontdict={"size": self.font_size},
                )

    def _get_nrows(self):
        pdict = {"res": 0, "phase": 1}
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
        return nrows, index, hr

    def _setup_subplots(self, gs_master, n_stations=1, n_index=0):
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

        nrows, index, hr = self._get_nrows()

        gs_rp = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=gs_master[0, n_index],
            height_ratios=[2, 1.5],
            hspace=0.05,
            wspace=0.15,
        )
        if nrows == 2:
            gs_aux = gridspec.GridSpecFromSubplotSpec(
                index, 1, subplot_spec=gs_master[1, n_index], hspace=0.05
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
        nrows, n_index, hr = self._get_nrows()
        gs_master = gridspec.GridSpec(nrows, ns, hspace=0.05, height_ratios=hr)
        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        for ii, mt in enumerate(self.tf_list):
            axr, axp, axr2, axp2, axt, axpt, label_coords = self._setup_subplots(
                gs_master, n_stations=ns, index=ii
            )

            # plot apparent resistivity od
            if self.plot_num == 1:
                self._plot_resistivity(axr, mt.period, mt.Z, mode="od", index=ii)

                # plot phase od
                self._plot_phase(axp, mt.period, mt.Z, mode="od", index=ii)
            # Plot Determinant
            elif self.plot_num == 3:
                # plot apparent resistivity od
                self._plot_resistivity(axr, mt.period, mt.Z, mode="det", index=ii)

                # plot phase od
                self._plot_phase(axp, mt.period, mt.Z, mode="det", index=ii)
            # plot diagonal components
            if self.plot_num == 2:
                # plot apparent resistivity od
                self._plot_resistivity(axr2, mt.period, mt.Z, mode="d", index=ii)

                # plot phase od
                self._plot_phase(axp2, mt.period, mt.Z, mode="d", index=ii)
            # plot tipper
            self._plot_tipper(axt, mt.period, mt.Tipper, index=ii)

            # plot phase tensor
            self._plot_pt(axpt, mt.period, mt.pt, index=ii)

            axr.set_title(mt.station, fontsize=self.font_size, fontweight="bold")

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
            self._plot_all()
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

            # make marker lists for the different components
            mxy = ["s", "D", "x", "+", "*", "1", "3", "4"] * ns
            myx = ["o", "h", "8", "p", "H", 7, 4, 6] * ns

            legend_list_xy = []
            legend_list_yx = []
            station_list = []

            nrows, n_index, hr = self._get_nrows()
            gs_master = gridspec.GridSpec(nrows, 1, hspace=0.15, height_ratios=hr)
            axr, axp, axr2, axp2, axt, axpt, label_coords = self._setup_subplots(
                gs_master
            )

            # make a figure instance
            self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

            for ii, mt in enumerate(self.mt_list):
                self.xy_color = cxy[ii]
                self.xy_marker = mxy[ii]
                self.yx_color = cyx[ii]
                self.yx_marker = myx[ii]
                self.det_color = cdet[ii]
                self.det_marker = mxy[ii]
                self.arrow_color_real = ctipr[ii]
                self.arrow_color_imag = ctipi[ii]

                # plot apparent resistivity od
                if self.plot_num == 1:
                    eb_list, label_list = self._plot_resistivity(
                        axr, mt.period, mt.Z, mode="od"
                    )

                    # plot phase od
                    self._plot_phase(axp, mt.period, mt.Z, mode="od", index=ii)
                # Plot Determinant
                elif self.plot_num == 3:
                    # plot apparent resistivity od
                    eb_list, label_list = self._plot_resistivity(
                        axr, mt.period, mt.Z, mode="det"
                    )

                    # plot phase od
                    self._plot_phase(axp, mt.period, mt.Z, mode="det")
                # plot diagonal components
                if self.plot_num == 2:
                    # plot apparent resistivity od
                    self._plot_resistivity(axr2, mt.period, mt.Z, mode="d")

                    # plot phase od
                    self._plot_phase(axp2, mt.period, mt.Z, mode="d")
                # plot tipper
                self._plot_tipper(axt, mt.period, mt.Tipper)

                # plot phase tensor
                self._plot_pt(axpt, mt.period, mt.pt)

                legend_list_xy += eb_list[0]
                legend_list_yx += eb_list[1]
                station_list.append(mt.station)
            # make legend
            if self.plot_num == 1:
                axr.legend(
                    legend_list_xy,
                    station_list,
                    loc=3,
                    ncol=2,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )

                axr2.legend(
                    legend_list_yx,
                    station_list,
                    loc=3,
                    ncol=2,
                    markerscale=0.75,
                    borderaxespad=0.01,
                    labelspacing=0.07,
                    handletextpad=0.2,
                    borderpad=0.25,
                )
