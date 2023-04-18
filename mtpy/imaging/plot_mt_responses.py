# -*- coding: utf-8 -*-
"""
plots multiple MT responses simultaneously

Created on Thu May 30 17:02:39 2013
@author: jpeacock-pr

YG: the code there is massey, todo may need to rewrite it sometime

"""

# ============================================================================

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

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

    def __init__(self, mt_data, **kwargs):
        """
        Initialize parameters
        """
        self.plot_num = 1
        self.plot_style = "1"
        self.mt_data = mt_data
        self.include_survey = True

        super().__init__(**kwargs)

        self.plot_dict = dict(
            [
                (kk, vv)
                for kk, vv in zip(
                    ["tip", "pt", "strike", "skew"],
                    [
                        self.plot_tipper,
                        self.plot_pt,
                        self.plot_strike,
                        self.plot_skew,
                    ],
                )
            ]
        )

        # set arrow properties
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = 0.5
        self.plot_model_error = None

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
        for tf in self.mt_data:
            tf.rotation_angle = value
        self._rotation_angle = value

    @property
    def plot_model_error(self):
        return self._plot_model_error

    @plot_model_error.setter
    def plot_model_error(self, value):
        if value:
            self._error_str = "model_error"
        else:
            self._error_str = "error"

        self._plot_model_error = value

    def _plot_resistivity(
        self, axr, period, z_obj, mode="od", index=0, axr2=None
    ):

        if mode == "od":
            comps = ["xy", "yx"]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
            ]
            if axr2 is not None:
                ax_list = [axr, axr2]
            else:
                ax_list = [axr, axr]
        elif mode == "d":
            comps = ["xx", "yy"]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
            ]
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
        elif mode == "det_only":
            comps = ["det"]
            props = [self.det_error_bar_properties]
            ax_list = [axr]
        res_limits = self.set_resistivity_limits(z_obj.resistivity, mode=mode)
        x_limits = self.set_period_limits(period)

        eb_list = []
        label_list = []
        for comp, prop, ax in zip(comps, props, ax_list):
            ebax = plot_resistivity(
                ax,
                period,
                getattr(z_obj, f"res_{comp}"),
                getattr(z_obj, f"res_{self._error_str}_{comp}"),
                **prop,
            )
            eb_list.append(ebax[0])
            label_list.append("$Z_{" + comp + "}$")
            # --> set axes properties
            plt.setp(ax.get_xticklabels(), visible=False)

            ax.set_yscale("log", nonpositive="clip")
            ax.set_xscale("log", nonpositive="clip")
            ax.set_xlim(x_limits)
            ax.set_ylim(res_limits)
            ax.grid(
                True,
                alpha=0.25,
                which="both",
                color=(0.25, 0.25, 0.25),
                lw=0.25,
            )
        if index == 0:
            axr.set_ylabel(
                "App. Res. ($\mathbf{\Omega \cdot m}$)",
                fontdict=self.font_dict,
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

    def _plot_phase(self, axp, period, z_obj, mode="od", index=0, axp2=None):
        if mode == "od":
            comps = ["xy", "yx"]
            if axp2 is not None:
                ax_list = [axp, axp2]
            else:
                ax_list = [axp, axp]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
            ]
        elif mode == "d":
            comps = ["xx", "yy"]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
            ]
        elif mode == "det":
            comps = ["xy", "yx", "det"]
            props = [
                self.xy_error_bar_properties,
                self.yx_error_bar_properties,
                self.det_error_bar_properties,
            ]
            if axp2 is not None:
                ax_list = [axp, axp2, axp]
            else:
                ax_list = [axp, axp, axp]
        elif mode == "det_only":
            comps = ["det"]
            props = [self.det_error_bar_properties]
            ax_list = [axp]
        phase_limits = self.set_phase_limits(z_obj.phase, mode=mode)

        for comp, prop, ax in zip(comps, props, ax_list):
            if comp == "yx":
                plot_phase(
                    ax,
                    period,
                    getattr(z_obj, f"phase_{comp}"),
                    getattr(z_obj, f"phase_{self._error_str}_{comp}"),
                    yx=True,
                    **prop,
                )
            else:
                plot_phase(
                    ax,
                    period,
                    getattr(z_obj, f"phase_{comp}"),
                    getattr(z_obj, f"phase_{self._error_str}_{comp}"),
                    yx=False,
                    **prop,
                )
            ax.set_ylim(phase_limits)
            if phase_limits[0] < -10 or phase_limits[1] > 100:
                ax.yaxis.set_major_locator(MultipleLocator(30))
                ax.yaxis.set_minor_locator(MultipleLocator(10))
            else:
                ax.yaxis.set_major_locator(MultipleLocator(15))
                ax.yaxis.set_minor_locator(MultipleLocator(5))
            ax.grid(
                True,
                alpha=0.25,
                which="both",
                color=(0.25, 0.25, 0.25),
                lw=0.25,
            )
            ax.set_xscale("log", nonpositive="clip")
            if "y" not in self.plot_tipper and not self.plot_pt:
                ax.set_xlabel("Period (s)", self.font_dict)
        # --> set axes properties
        if index == 0:
            axp.set_ylabel("Phase (deg)", self.font_dict)

    def _plot_tipper(
        self, axt, period, t_obj, index=0, legend=False, zero_reference=False
    ):
        if t_obj is None:
            return None, None

        axt, tip_list, tip_label = plot_tipper_lateral(
            axt,
            t_obj,
            self.plot_tipper,
            self.arrow_real_properties,
            self.arrow_imag_properties,
            self.font_size,
            legend=legend,
            zero_reference=zero_reference,
        )
        if axt is None:
            return None, None

        axt.set_xlabel("Period (s)", fontdict=self.font_dict)

        axt.yaxis.set_major_locator(MultipleLocator(0.2))
        axt.yaxis.set_minor_locator(MultipleLocator(0.1))
        axt.set_xlabel("Period (s)", fontdict=self.font_dict)
        if index == 0:
            axt.set_ylabel("Tipper", fontdict=self.font_dict)
        # set th xaxis tick labels to invisible
        if self.plot_pt:
            plt.setp(axt.xaxis.get_ticklabels(), visible=False)
            axt.set_xlabel("")
        return tip_list, tip_label

    def _plot_pt(
        self, axpt, period, pt_obj, index=0, y_shift=0, edge_color=None
    ):
        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt:

            color_array = self.get_pt_color_array(pt_obj)
            x_limits = self.set_period_limits(period)

            # -------------plot ellipses-----------------------------------
            self.cbax, self.cbpt, = plot_pt_lateral(
                axpt,
                pt_obj,
                color_array,
                self.ellipse_properties,
                y_shift,
                self.fig,
                edge_color,
                index,
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
                True,
                alpha=0.25,
                which="major",
                color=(0.25, 0.25, 0.25),
                lw=0.25,
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
        return nrows, index, hr, pdict

    def _setup_subplots(
        self,
        gs_master,
        n_stations=1,
        n_index=0,
        plot_num=1,
        hspace=0.05,
        wspace=0.15,
    ):
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

        nrows, index, hr, pdict = self._get_nrows()

        gs_rp = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=gs_master[0, n_index],
            height_ratios=[2, 1.5],
            hspace=hspace,
            wspace=wspace,
        )
        if nrows == 2:
            gs_aux = gridspec.GridSpecFromSubplotSpec(
                index, 1, subplot_spec=gs_master[1, n_index], hspace=hspace
            )
        # --> make figure for xy,yx components
        if plot_num == 1 or plot_num == 3:
            # set label coordinates
            if self.plot_style == "compare":
                label_coords = (-0.075, 0.5)
            elif self.plot_style == "compare":
                label_coords = (-0.095, 0.5)
            # --> create the axes instances
            # apparent resistivity axis
            axr = self.fig.add_subplot(gs_rp[0, :])

            # phase axis that shares period axis with resistivity
            axp = self.fig.add_subplot(gs_rp[1, :], sharex=axr)
        # --> make figure for all 4 components
        elif plot_num == 2:
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
            axt = self.fig.add_subplot(
                gs_aux[pdict["tip"], :],
            )
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
        ns = self.mt_data.n_stations

        # set figure size according to what the plot will be.
        if self.fig_size is None:
            if self.plot_num == 1 or self.plot_num == 3:
                self.fig_size = [ns * 4, 6]
            elif self.plot_num == 2:
                self.fig_size = [ns * 8, 6]
        nrows, n_index, hr, pdict = self._get_nrows()
        gs_master = gridspec.GridSpec(nrows, ns, hspace=0.05, height_ratios=hr)
        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        for ii, mt in enumerate(self.mt_data.values()):
            (
                axr,
                axp,
                axr2,
                axp2,
                axt,
                axpt,
                label_coords,
            ) = self._setup_subplots(
                gs_master,
                n_stations=ns,
                n_index=ii,
                plot_num=self.plot_num,
                hspace=0.075,
                wspace=0.09,
            )

            # plot apparent resistivity od
            if self.plot_num == 1:
                self._plot_resistivity(
                    axr, mt.period, mt.Z, mode="od", index=ii
                )
                if self.res_limits is not None:
                    axr.set_ylim(self.res_limits)
                # plot phase od
                self._plot_phase(axp, mt.period, mt.Z, mode="od", index=ii)
                if self.phase_limits is not None:
                    axp.set_ylim(self.phase_limits)
            # Plot Determinant
            elif self.plot_num == 3:
                # plot apparent resistivity od
                self._plot_resistivity(
                    axr, mt.period, mt.Z, mode="det", index=ii
                )
                if self.res_limits is not None:
                    axr.set_ylim(self.res_limits)
                # plot phase od
                self._plot_phase(axp, mt.period, mt.Z, mode="det", index=ii)
                if self.phase_limits is not None:
                    axp.set_ylim(self.phase_limits)
            # plot diagonal components
            if self.plot_num == 2:
                # plot apparent resistivity od
                self._plot_resistivity(
                    axr2, mt.period, mt.Z, mode="d", index=ii
                )

                # plot phase od
                self._plot_phase(axp2, mt.period, mt.Z, mode="d", index=ii)
            # plot tipper
            self._plot_tipper(axt, mt.period, mt.Tipper, index=ii)
            if self.tipper_limits is not None:
                axt.set_ylim(self.tipper_limits)
            # plot phase tensor
            self._plot_pt(axpt, mt.period, mt.pt, index=ii)

            axr.set_title(
                mt.station, fontsize=self.font_size, fontweight="bold"
            )

    def _plot_compare(self):
        # plot diagonal components
        if self.plot_num == 2:
            raise ValueError(
                "Compare mode does not support plotting diagonal components yet"
            )
        ns = self.mt_data.n_stations

        # make color lists for the plots going light to dark
        cxy = [(0, 0 + float(cc) / ns, 1 - float(cc) / ns) for cc in range(ns)]
        cyx = [(1, float(cc) / ns, 0) for cc in range(ns)]
        cdet = [(0, 1 - float(cc) / ns, 0) for cc in range(ns)]
        ctipr = [
            (0.75 * cc / ns, 0.75 * cc / ns, 0.75 * cc / ns)
            for cc in range(ns)
        ]
        ctipi = [
            (float(cc) / ns, 1 - float(cc) / ns, 0.25) for cc in range(ns)
        ]

        # make marker lists for the different components
        mxy = ["s", "D", "x", "+", "*", "1", "3", "4"] * ns
        myx = ["o", "h", "8", "p", "H", 7, 4, 6] * ns

        legend_list_xy = []
        legend_list_yx = []
        legend_list_tip = []
        station_list = []
        station_list_t = []

        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        nrows, n_index, hr, pdict = self._get_nrows()
        gs_master = gridspec.GridSpec(nrows, 1, hspace=0.15, height_ratios=hr)
        if self.plot_num == 1:
            (
                axr,
                axp,
                axr2,
                axp2,
                axt,
                axpt,
                label_coords,
            ) = self._setup_subplots(gs_master, plot_num=2)
        elif self.plot_num == 3:
            (
                axr,
                axp,
                axr2,
                axp2,
                axt,
                axpt,
                label_coords,
            ) = self._setup_subplots(gs_master, plot_num=1)

        period = []

        for ii, mt in enumerate(self.mt_data.values()):
            period.append(mt.period.min())
            period.append(mt.period.max())
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
                    axr,
                    mt.period,
                    mt.Z,
                    mode="od",
                    axr2=axr2,
                )

                # plot phase od
                self._plot_phase(
                    axp, mt.period, mt.Z, mode="od", index=ii, axp2=axp2
                )
            # Plot Determinant
            elif self.plot_num == 3:
                # plot apparent resistivity od
                eb_list, label_list = self._plot_resistivity(
                    axr, mt.period, mt.Z, mode="det_only"
                )

                # plot phase od
                self._plot_phase(axp, mt.period, mt.Z, mode="det_only")
            # plot tipper
            tip_list, tip_label = self._plot_tipper(axt, mt.period, mt.Tipper)
            if tip_list is not None:
                if self.plot_tipper.find("r") > 0:
                    legend_list_tip.append(tip_list[0])
                    station_list_t.append(f"{mt.station}_{tip_label[0]}")
                    if self.plot_tipper.find("i") > 0:
                        legend_list_tip.append(tip_list[1])
                        station_list_t.append(f"{mt.station}_{tip_label[1]}")
                elif self.plot_tipper.find("i") > 0:
                    legend_list_tip.append(tip_list[0])
                    station_list_t.append(f"{mt.station}_{tip_label[0]}")
            # plot phase tensor
            self._plot_pt(
                axpt,
                mt.period,
                mt.pt,
                y_shift=ii * self.ellipse_size,
                edge_color=cxy[ii],
            )

            legend_list_xy += [eb_list[0]]
            if self.plot_num in [1, 2]:
                legend_list_yx += [eb_list[1]]
            if self.include_survey:
                station_list.append(f"{mt.station}_{mt.survey_metadata.id}")
            else:
                station_list.append(f"{mt.station}")

        # set limits
        if self.res_limits is not None:
            axr.set_ylim(self.res_limits)
            if axr2 is not None:
                axr2.set_ylim(self.res_limits)
        if self.phase_limits is not None:
            axp.set_ylim(self.phase_limits)
            if axp2 is not None:
                axp2.set_ylim(self.phase_limits)
        if self.tipper_limits is not None:
            axt.set_ylim(self.tipper_limits)

        period_limits = [
            10 ** np.floor(np.log10(min(period))),
            10 ** np.ceil(np.log10(max(period))),
        ]
        for ax in [axr, axp]:
            if ax is not None:
                ax.set_xlim(period_limits)
        if ax in [axt, axpt]:
            if ax is not None:
                ax.set_xlim(
                    [np.log10(period_limits[0]), np.log10(period_limits[1])]
                )

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
                prop={"size": self.font_size - 2},
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
                prop={"size": self.font_size - 2},
            )
        elif self.plot_num == 3:
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
        if self.plot_tipper.find("y") >= 0:
            axt.legend(
                legend_list_tip,
                station_list_t,
                loc=3,
                ncol=2,
                markerscale=0.75,
                borderaxespad=0.01,
                labelspacing=0.07,
                handletextpad=0.2,
                borderpad=0.25,
            )
        self.axr = axr
        self.axp = axp
        self.axr2 = axr2
        self.axp2 = axp2
        self.axt = axt
        self.axpt = axpt

    def _plot_single(self):
        p_dict = {}
        for ii, tf in enumerate(self.mt_data.values(), 1):
            p = tf.plot_mt_response(
                **{
                    "fig_num": ii,
                    "plot_tipper": self.plot_tipper,
                    "plot_pt": self.plot_pt,
                    "plot_num": self.plot_num,
                }
            )
            p_dict[tf.station] = p
        return p_dict

    # ---plot the resistivity and phase
    def plot(self):
        """
        plot the apparent resistivity and phase
        """

        plt.clf()
        self.subplot_right = 0.98
        self._set_subplot_params()

        # Plot all in one figure as subplots
        if self.plot_style == "all":
            self.subplot_left = 0.04
            self.subplot_top = 0.96
            self._set_subplot_params()
            self._plot_all()
        #  Plot all responses into one plot to compare changes
        elif self.plot_style == "compare":
            self._plot_compare()
        elif self.plot_style in [1, "1", "single"]:
            return self._plot_single()
