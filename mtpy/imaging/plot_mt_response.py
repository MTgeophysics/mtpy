# -*- coding: utf-8 -*-
"""
=================
plot_mt_response
=================

Plots the resistivity and phase for different modes and components

Created 2017

@author: jpeacock
"""
# ==============================================================================
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.colors as colors
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.gridspec as gridspec

import mtpy.imaging.mtcolors as mtcl
from mtpy.imaging.mtplot_tools import (
    PlotBase,
    plot_errorbar,
    plot_pt_lateral,
    get_log_tick_labels,
    make_color_list,
)

# ==============================================================================
#  Plot apparent resistivity and phase
# ==============================================================================


class PlotMTResponse(PlotBase):
    """
    Plots Resistivity and phase for the different modes of the MT response.  At
    the moment it supports the input of an .edi file. Other formats that will
    be supported are the impedance tensor and errors with an array of periods
    and .j format.

    The normal use is to input an .edi file, however it would seem that not
    everyone uses this format, so you can input the data and put it into 
    arrays or objects like class mtpy.core.z.Z.  Or if the data is in 
    resistivity and phase format they can be input as arrays or a class
    mtpy.imaging.mtplot.ResPhase.  Or you can put it into a class
    mtpy.imaging.mtplot.MTplot.

    The plot places the apparent resistivity in log scale in the top panel(s), 
    depending on the plot_num.  The phase is below this, note that 180 degrees
    has been added to the yx phase so the xy and yx phases plot in the same
    quadrant.  Both the resistivity and phase share the same x-axis which is in
    log period, short periods on the left to long periods on the right.  So
    if you zoom in on the plot both plots will zoom in to the same 
    x-coordinates.  If there is tipper information, you can plot the tipper
    as a third panel at the bottom, and also shares the x-axis.  The arrows are
    in the convention of pointing towards a conductor.  The xx and yy 
    components can be plotted as well, this adds two panels on the right.  
    Here the phase is left unwrapped.  Other parameters can be added as 
    subplots such as strike, skew and phase tensor ellipses.

    To manipulate the plot you can change any of the attributes listed below
    and call redraw_plot().  If you know more aout matplotlib and want to 
    change axes parameters, that can be done by changing the parameters in the
    axes attributes and then call update_plot(), note the plot must be open.


    """

    def __init__(
        self, z_object=None, t_object=None, pt_obj=None, station="MT Response", **kwargs
    ):
        super().__init__(**kwargs)

        self.Z = z_object
        self.Tipper = t_object
        self.pt = pt_obj
        self.station = station
        self._basename = f"{self.station}_mt_response"

        self.phase_quadrant = 1

        self.plot_num = kwargs.pop("plot_num", 1)
        self.rotation_angle = kwargs.pop("rotation_angle", 0)

        if self.Tipper is not None:
            self.plot_tipper = "yri"
        if self.pt is not None:
            self.plot_pt = True
        # set arrow properties
        self.arrow_size = 1
        self.arrow_head_length = 0.03
        self.arrow_head_width = 0.03
        self.arrow_lw = 0.5

        # ellipse_properties
        self.ellipse_size = 0.25
        self.ellipse_spacing = 1
        if self.ellipse_size == 2 and self.ellipse_spacing == 1:
            self.ellipse_size = 0.25
        # layout params
        self.show_resphase_xticklabels = False

        # plot on initializing
        if self.show_plot:
            self.plot()

    @property
    def period(self):
        """
        plot period
        """
        if self.Z is not None:
            return 1.0 / self.Z.freq
        elif self.Tipper is not None:
            return 1.0 / self.Tipper.freq
        else:
            return None

    def _has_tipper(self):
        if self.plot_tipper.find("y") >= 0:
            if self.Tipper is None or (self.Tipper.tipper == 0 + 0j).all():
                self._logger.info(f"No Tipper data for station {self.station}")
                self.plot_tipper = "n"

    def _has_pt(self):
        if self.plot_pt:
            # if np.all(self.Z.z == 0 + 0j) or self.Z is None:
            if self.pt is None:  # no phase tensor object provided
                self._logger.info(f"No PT data for station {self.station}")
                self.plot_pt = False

    def _set_subplot_params(self):
        # set some parameters of the figure and subplot spacing
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.bottom"] = 0.1
        plt.rcParams["figure.subplot.top"] = 0.93
        plt.rcParams["figure.subplot.left"] = 0.12
        plt.rcParams["figure.subplot.right"] = 0.98

    def _setup_subplots(self):
        # create a dictionary for the number of subplots needed
        pdict = {"res": 0, "phase": 1}
        # start the index at 2 because resistivity and phase is permanent for
        # now
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
        gs_master = gridspec.GridSpec(nrows, 1, hspace=0.15, height_ratios=hr)
        gs_rp = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=gs_master[0],
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
            self.axr = self.fig.add_subplot(gs_rp[0, :])

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs_rp[1, :], sharex=self.axr)
        # --> make figure for all 4 components
        elif self.plot_num == 2:
            # set label coordinates
            label_coords = (-0.095, 0.5)

            # --> create the axes instances
            # apparent resistivity axis
            self.axr = self.fig.add_subplot(gs_rp[0, 0])
            self.axr2 = self.fig.add_subplot(gs_rp[0, 1], sharex=self.axr)
            self.axr2.yaxis.set_label_coords(-0.1, 0.5)

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs_rp[1, 0], sharex=self.axr)
            self.axp2 = self.fig.add_subplot(gs_rp[1, 1], sharex=self.axr)
            self.axp2.yaxis.set_label_coords(-0.1, 0.5)
        # set albel coordinates
        self.axr.yaxis.set_label_coords(label_coords[0], label_coords[1])
        self.axp.yaxis.set_label_coords(label_coords[0], label_coords[1])

        # --> plot tipper
        if self.plot_tipper.find("y") >= 0:
            self.axt = self.fig.add_subplot(gs_aux[pdict["tip"], :],)
            self.axt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        # --> plot phase tensors
        if self.plot_pt:
            # can't share axis because not on the same scale
            # Removed aspect = "equal" for now, it flows better, if you want
            # a detailed analysis look at plot pt
            self.axpt = self.fig.add_subplot(gs_aux[pdict["pt"], :])
            self.axpt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        return label_coords

    def _get_nonzero_indices(self):
        self._nz_xx = np.nonzero(self.Z.z[:, 0, 0])
        self._nz_xy = np.nonzero(self.Z.z[:, 0, 1])
        self._nz_yx = np.nonzero(self.Z.z[:, 1, 0])
        self._nz_yy = np.nonzero(self.Z.z[:, 1, 1])

        self._nz_tx = None
        self._nz_ty = None

        if self.Tipper is not None:  # fix github issue #24.
            # NOTE the following lines seems not have any effect anyway
            self._nz_tx = np.nonzero(self.Tipper.tipper[:, 0, 0])
            self._nz_ty = np.nonzero(self.Tipper.tipper[:, 0, 1])

    def _plot_resistivity_od(self):

        res_limits = self.set_resistivity_limits(self.Z.resistivity, mode="od")

        # res_xy
        self.ebxyr = plot_errorbar(
            self.axr,
            self.period[self._nz_xy],
            self.Z.res_xy[self._nz_xy],
            y_error=self.Z.res_err_xy[self._nz_xy],
            **self.xy_error_bar_properties,
        )

        # res_yx
        self.ebyxr = plot_errorbar(
            self.axr,
            self.period[self._nz_yx],
            self.Z.res_yx[self._nz_yx],
            y_error=self.Z.res_err_yx[self._nz_yx],
            **self.yx_error_bar_properties,
        )

        # --> set axes properties
        plt.setp(self.axr.get_xticklabels(), visible=False)
        self.axr.set_ylabel(
            "App. Res. ($\mathbf{\Omega \cdot m}$)", fontdict=self.font_dict
        )
        self.axr.set_yscale("log", nonpositive="clip")
        self.axr.set_xscale("log", nonpositive="clip")
        self.axr.set_xlim(self.x_limits)
        self.axr.set_ylim(res_limits)
        self.axr.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)

        self.axr.legend(
            (self.ebxyr[0], self.ebyxr[0]),
            ("$Z_{xy}$", "$Z_{yx}$"),
            loc=3,
            markerscale=1,
            borderaxespad=0.01,
            labelspacing=0.07,
            handletextpad=0.2,
            borderpad=0.02,
        )

    def _plot_resistivity_d(self):

        res_limits = self.set_resistivity_limits(self.Z.resistivity, mode="d")
        # res_xx
        self.ebxxr = plot_errorbar(
            self.axr2,
            self.period[self._nz_xx],
            self.Z.res_xx[self._nz_xx],
            y_error=self.Z.res_err_xx[self._nz_xx],
            **self.xy_error_bar_properties,
        )

        # res_yy
        self.ebyyr = plot_errorbar(
            self.axr2,
            self.period[self._nz_yy],
            self.Z.res_yy[self._nz_yy],
            y_error=self.Z.res_err_yy[self._nz_yy],
            **self.yx_error_bar_properties,
        )

        # --> set axes properties
        plt.setp(self.axr2.get_xticklabels(), visible=False)
        self.axr2.set_yscale("log", nonpositive="clip")
        self.axr2.set_xscale("log", nonpositive="clip")
        self.axr2.set_xlim(self.x_limits)
        self.axr2.set_ylim(res_limits)
        self.axr2.grid(
            True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
        )

        self.axr2.legend(
            (self.ebxxr[0], self.ebyyr[0]),
            ("$Z_{xx}$", "$Z_{yy}$"),
            loc=3,
            markerscale=1,
            borderaxespad=0.01,
            labelspacing=0.07,
            handletextpad=0.2,
            borderpad=0.02,
        )

    def _plot_phase_od(self):
        # phase_xy
        self.ebxyp = plot_errorbar(
            self.axp,
            self.period[self._nz_xy],
            self.Z.phase_xy[self._nz_xy],
            y_error=self.Z.phase_err_xy[self._nz_xy],
            **self.xy_error_bar_properties,
        )

        # phase_yx: Note add 180 to place it in same quadrant as phase_xy
        self.ebyxp = plot_errorbar(
            self.axp,
            self.period[self._nz_yx],
            self.Z.phase_yx[self._nz_yx] + 180,
            y_error=self.Z.phase_err_yx[self._nz_yx],
            **self.yx_error_bar_properties,
        )

        # check the phase to see if any point are outside of [0:90]
        phase_limits = self.set_phase_limits(self.Z.phase)
        # --> set axes properties
        if self.plot_tipper.find("y") < 0 or not self.plot_pt:
            self.axp.set_xlabel("Period (s)", self.font_dict)
        self.axp.set_ylabel("Phase (deg)", self.font_dict)
        self.axp.set_xscale("log", nonpositive="clip")
        self.axp.set_ylim(phase_limits)
        self.axp.yaxis.set_major_locator(MultipleLocator(15))
        self.axp.yaxis.set_minor_locator(MultipleLocator(5))
        self.axp.grid(True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25)
        # # set th xaxis tick labels to invisible
        # if self.plot_tipper.find("y") >= 0 or self.plot_pt == "y":
        #     plt.setp(self.axp.xaxis.get_ticklabels(), visible=False)
        #     self.axp.set_xlabel("")

    def _plot_phase_d(self):
        # phase_xx
        self.ebxyp = plot_errorbar(
            self.axp2,
            self.period[self._nz_xx],
            self.Z.phase_xx[self._nz_xx],
            y_error=self.Z.phase_err_xx[self._nz_xx],
            **self.xy_error_bar_properties,
        )

        # phase_yy
        self.ebyxp = plot_errorbar(
            self.axp2,
            self.period[self._nz_yy],
            self.Z.phase_yy[self._nz_yy],
            y_error=self.Z.phase_err_yy[self._nz_yy],
            **self.yx_error_bar_properties,
        )

        # --> set axes properties
        if self.plot_tipper.find("y") < 0 or not self.plot_pt:
            self.axp2.set_xlabel("Period (s)", self.font_dict)
        self.axp2.set_xscale("log", nonpositive="clip")
        self.axp2.set_ylim(ymin=-179.9, ymax=179.9)
        self.axp2.yaxis.set_major_locator(MultipleLocator(30))
        self.axp2.yaxis.set_minor_locator(MultipleLocator(5))
        self.axp2.grid(
            True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
        )

    def _plot_determinant(self):
        # res_det
        self.ebdetr = self.axr.errorbar(
            self.period,
            self.Z.res_det,
            yerr=self.Z.res_det_err,
            **self.det_error_bar_properties,
        )

        # phase_det
        self.ebdetp = self.axp.errorbar(
            self.period,
            self.Z.phase_det,
            yerr=self.Z.phase_det_err,
            **self.det_error_bar_properties,
        )

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

    def _plot_tipper(self):
        # -----plot tipper----------------------------------------------------
        if self.plot_tipper.find("y") == 0 or self.plot_tipper:

            txr = self.Tipper.mag_real * np.cos(np.deg2rad(self.Tipper.angle_real))
            tyr = self.Tipper.mag_real * np.sin(np.deg2rad(self.Tipper.angle_real))

            txi = self.Tipper.mag_imag * np.cos(np.deg2rad(self.Tipper.angle_imag))
            tyi = self.Tipper.mag_imag * np.sin(np.deg2rad(self.Tipper.angle_imag))

            nt = len(txr)

            tiplist = []
            tiplabel = []

            for aa in range(nt):
                xlenr = txr[aa] * np.log10(self.period[aa])
                xleni = txi[aa] * np.log10(self.period[aa])

                # --> plot real arrows
                if self.plot_tipper.find("r") > 0:
                    self.axt.arrow(
                        np.log10(self.period[aa]),
                        0,
                        xlenr,
                        tyr[aa],
                        **self.arrow_real_properties,
                    )

                    if aa == 0:
                        line1 = self.axt.plot(0, 0, self.arrow_color_real)
                        tiplist.append(line1[0])
                        tiplabel.append("real")
                # --> plot imaginary arrows
                if self.plot_tipper.find("i") > 0:
                    self.axt.arrow(
                        np.log10(self.period[aa]),
                        0,
                        xleni,
                        tyi[aa],
                        **self.arrow_imag_properties,
                    )
                    if aa == 0:
                        line2 = self.axt.plot(0, 0, self.arrow_color_imag)
                        tiplist.append(line2[0])
                        tiplabel.append("imag")
            # make a line at 0 for reference
            self.axt.plot(np.log10(self.period), [0] * nt, "k", lw=0.5)

            self.axt.legend(
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

            # set axis properties

            self.axt.set_xlim(np.log10(self.x_limits[0]), np.log10(self.x_limits[1]))

            tklabels = []
            xticks = []

            for tk in self.axt.get_xticks():
                try:
                    tklabels.append(self.period_label_dict[tk])
                    xticks.append(tk)
                except KeyError:
                    pass
            self.axt.set_xticks(xticks)
            self.axt.set_xticklabels(tklabels, fontdict={"size": self.font_size})
            self.axt.set_xlabel("Period (s)", fontdict=self.font_dict)
            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            self.axt.set_xlim(np.log10(self.x_limits[0]), np.log10(self.x_limits[1]))

            self.axt.yaxis.set_major_locator(MultipleLocator(0.2))
            self.axt.yaxis.set_minor_locator(MultipleLocator(0.1))
            self.axt.set_xlabel("Period (s)", fontdict=self.font_dict)
            self.axt.set_ylabel("Tipper", fontdict=self.font_dict)

            # self.axt.set_xscale('log', nonpositive='clip')
            if self.tipper_limits is None:
                tmax = max([np.nanmax(tyr), np.nanmax(tyi)])
                if tmax > 1:
                    tmax = 0.899
                tmin = min([np.nanmin(tyr), np.nanmin(tyi)])
                if tmin < -1:
                    tmin = -0.899
                self.tipper_limits = (tmin - 0.1, tmax + 0.1)
            self.axt.set_ylim(self.tipper_limits)
            self.axt.grid(
                True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
            )

            # set th xaxis tick labels to invisible
            if self.plot_pt:
                plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
                self.axt.set_xlabel("")

    def _plot_pt(self):
        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt:

            cmap = self.ellipse_cmap
            ckmin = self.ellipse_range[0]
            ckmax = self.ellipse_range[1]
            bounds = None
            try:
                ckstep = float(self.ellipse_range[2])
            except IndexError:
                ckstep = 3
            if cmap == "mt_seg_bl2wh2rd":
                bounds = np.arange(ckmin, ckmax + ckstep, ckstep)
                nseg = float((ckmax - ckmin) / (2 * ckstep))
            color_array = self.get_pt_color_array(self.pt)

            # -------------plot ellipses-----------------------------------
            self.cbax, self.cbpt, = plot_pt_lateral(
                self.axpt,
                self.pt,
                color_array,
                self.ellipse_properties,
                bounds,
                self.fig,
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

    def plot(self, show=True):
        """
        plotResPhase(filename,fig_num) will plot the apparent resistivity and 
        phase for a single station. 

        """

        self._has_tipper()
        self._has_pt()
        self._get_nonzero_indices()

        # set x-axis limits from short period to long period
        if self.x_limits is None:
            self.x_limits = self.set_period_limits(self.period)
        if self.res_limits is None:
            self.res_limits = self.set_resistivity_limits(self.Z.resistivity)
        # make figure instance
        self._set_subplot_params()
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        self.fig.clf()

        self._setup_subplots()

        self._plot_resistivity_od()
        self._plot_phase_od()
        self._plot_tipper()
        self._plot_pt()

        # ===Plot the xx, yy components if desired==============================
        if self.plot_num == 2:
            self._plot_resistivity_d()
            self._plot_phase_d()
        # ===Plot the Determinant if desired==================================
        if self.plot_num == 3:
            self._plot_determinant()
        if self.show_resphase_xticklabels:
            if self.plot_num in [1, 3]:
                self.gs.update(hspace=0.2, wspace=0.15, left=0.1)
            else:
                self.gs.update(hspace=0.2, wspace=0.15, left=0.07)
                plt.setp(self.axp2.xaxis.get_ticklabels(), visible=True)
                plt.setp(self.axr2.xaxis.get_ticklabels(), visible=True)
                self.axr2.tick_params(
                    axis="x",
                    pad=2,
                    direction="in",
                    which="both",
                    labelsize=self.font_size - 1,
                )
                self.axp2.tick_params(
                    axis="x",
                    pad=2,
                    direction="in",
                    which="both",
                    labelsize=self.font_size - 1,
                )
                self.axp2.set_xlabel(
                    "Period (s)", fontsize=self.font_size - 1, labelpad=0
                )  #
            plt.setp(self.axr.xaxis.get_ticklabels(), visible=True)
            plt.setp(self.axp.xaxis.get_ticklabels(), visible=True)
            self.axr.tick_params(
                axis="x",
                pad=2,
                direction="in",
                which="both",
                labelsize=self.font_size - 1,
            )
            self.axp.tick_params(
                axis="x",
                pad=2,
                direction="in",
                which="both",
                labelsize=self.font_size - 1,
            )
        #            self.axp.set_xlabel('Period (s)',fontsize=self.font_size-2,labelpad=0)

        # make plot_title and show
        if self.plot_title is None:
            self.plot_title = self.station
        self.fig.suptitle(self.plot_title, fontdict=self.font_dict)

        # be sure to show
        if show:
            plt.show()
