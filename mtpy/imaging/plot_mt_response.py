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
import matplotlib.gridspec as gridspec

from mtpy.imaging.mtplot_tools import (
    PlotBase,
    plot_pt_lateral,
    get_log_tick_labels,
    plot_resistivity,
    plot_phase,
    plot_tipper_lateral,
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
        self,
        z_object=None,
        t_object=None,
        pt_obj=None,
        station="MT Response",
        **kwargs,
    ):
        self.Z = z_object
        self.Tipper = t_object
        self.pt = pt_obj
        self.station = station
        self._basename = f"{self.station}_mt_response"
        self.plot_num = 1
        self.rotation_angle = 0
        super().__init__(**kwargs)

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
        # subplot parameters
        self.subplot_left = 0.12
        self.subplot_right = 0.98
        self.subplot_bottom = 0.1
        self.subplot_top = 0.93

        self.plot_model_error = False

        for key, value in kwargs.items():
            setattr(self, key, value)

        # plot on initializing
        if self.show_plot:
            self.plot()

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

    @property
    def period(self):
        """
        plot period
        """
        if self.Z is not None:
            return 1.0 / self.Z.frequency
        elif self.Tipper is not None:
            return 1.0 / self.Tipper.frequency
        else:
            return None

    # ---need to rotate data on setting rotz
    @property
    def rotation_angle(self):
        return self._rotation_angle

    @rotation_angle.setter
    def rotation_angle(self, theta_r):
        """
        only a single value is allowed
        """
        if not theta_r == 0:

            self.Z.rotate(theta_r)
            self.Tipper.rotate(theta_r)
            self.pt.rotate(theta_r)

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
                index, 1, subplot_spec=gs_master[1], hspace=0.075
            )
        # --> make figure for xy,yx components
        if self.plot_num == 1 or self.plot_num == 3:
            # set label coordinates
            label_coords = (-0.095, 0.5)

            # --> create the axes instances
            # apparent resistivity axis
            self.axr = self.fig.add_subplot(gs_rp[0, :])

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs_rp[1, :], sharex=self.axr)
        # --> make figure for all 4 components
        elif self.plot_num == 2:
            # set label coordinates
            label_coords = (-0.14, 0.5)

            # --> create the axes instances
            # apparent resistivity axis
            self.axr = self.fig.add_subplot(gs_rp[0, 0])
            self.axr2 = self.fig.add_subplot(gs_rp[0, 1], sharex=self.axr)
            self.axr2.yaxis.set_label_coords(label_coords[0], label_coords[1])

            # phase axis that shares period axis with resistivity
            self.axp = self.fig.add_subplot(gs_rp[1, 0], sharex=self.axr)
            self.axp2 = self.fig.add_subplot(gs_rp[1, 1], sharex=self.axr)
            self.axp2.yaxis.set_label_coords(label_coords[0], label_coords[1])
        # set label coordinates
        self.axr.yaxis.set_label_coords(label_coords[0], label_coords[1])
        self.axp.yaxis.set_label_coords(label_coords[0], label_coords[1])

        # --> plot tipper
        if self.plot_tipper.find("y") >= 0:
            self.axt = self.fig.add_subplot(
                gs_aux[pdict["tip"], :],
            )
            self.axt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        # --> plot phase tensors
        if self.plot_pt:
            # can't share axis because not on the same scale
            # Removed aspect = "equal" for now, it flows better, if you want
            # a detailed analysis look at plot pt
            self.axpt = self.fig.add_subplot(gs_aux[pdict["pt"], :])
            self.axpt.yaxis.set_label_coords(label_coords[0], label_coords[1])
        return label_coords

    def _plot_resistivity(self, axr, period, z_obj, mode="od"):

        if mode == "od":
            comps = ["xy", "yx"]
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
        res_limits = self.set_resistivity_limits(z_obj.resistivity, mode=mode)
        x_limits = self.set_period_limits(period)

        eb_list = []
        label_list = []
        for comp, prop in zip(comps, props):
            ebax = plot_resistivity(
                axr,
                period,
                getattr(z_obj, f"res_{comp}"),
                getattr(z_obj, f"res_{self._error_str}_{comp}"),
                **prop,
            )
            eb_list.append(ebax)
            label_list.append(f"$Z_{'{'}{comp}{'}'}$")
        # --> set axes properties
        plt.setp(axr.get_xticklabels(), visible=False)

        axr.set_yscale("log", nonpositive="clip")
        axr.set_xscale("log", nonpositive="clip")
        axr.set_xlim(x_limits)
        axr.set_ylim(res_limits)
        axr.grid(
            True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
        )

        if mode == "od":
            axr.set_ylabel(
                "App. Res. ($\mathbf{\Omega \cdot m}$)",
                fontdict=self.font_dict,
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
        return eb_list, label_list

    def _plot_phase(self, axp, period, z_obj, mode="od", index=0):
        if mode == "od":
            comps = ["xy", "yx"]
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
        phase_limits = self.set_phase_limits(z_obj.phase, mode=mode)

        for comp, prop in zip(comps, props):
            if comp == "yx":
                plot_phase(
                    axp,
                    period,
                    getattr(z_obj, f"phase_{comp}"),
                    getattr(z_obj, f"phase_{self._error_str}_{comp}"),
                    yx=True,
                    **prop,
                )
            else:
                plot_phase(
                    axp,
                    period,
                    getattr(z_obj, f"phase_{comp}"),
                    getattr(z_obj, f"phase_{self._error_str}_{comp}"),
                    yx=False,
                    **prop,
                )
        # --> set axes properties
        if mode == "od":
            axp.set_ylabel("Phase (deg)", self.font_dict)
        axp.set_xscale("log", nonpositive="clip")

        axp.set_ylim(phase_limits)
        if phase_limits[0] < -10 or phase_limits[1] > 100:
            axp.yaxis.set_major_locator(MultipleLocator(30))
            axp.yaxis.set_minor_locator(MultipleLocator(10))
        else:
            axp.yaxis.set_major_locator(MultipleLocator(15))
            axp.yaxis.set_minor_locator(MultipleLocator(5))
        axp.grid(
            True, alpha=0.25, which="both", color=(0.25, 0.25, 0.25), lw=0.25
        )

    def _plot_determinant(self):
        # res_det
        self.ebdetr = plot_resistivity(
            self.axr,
            self.period,
            self.Z.res_det,
            self.Z.res_error_det,
            **self.det_error_bar_properties,
        )

        # phase_det
        self.ebdetp = plot_phase(
            self.axp,
            self.period,
            self.Z.phase_det,
            self.Z.phase_error_det,
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

        self.axr.set_ylim(
            self.set_resistivity_limits(self.Z.resistivity, mode="det")
        )

    def _plot_tipper(self):
        if self.Tipper is not None:
            self.axt, _, _ = plot_tipper_lateral(
                self.axt,
                self.Tipper,
                self.plot_tipper,
                self.arrow_real_properties,
                self.arrow_imag_properties,
                self.font_size,
            )
            if self.plot_tipper.find("y") >= 0:
                self.axt.set_xlabel("Period (s)", fontdict=self.font_dict)
                # need to reset the x_limits caouse they get reset when calling
                # set_ticks for some reason
                self.axt.set_xlim(
                    np.log10(self.x_limits[0]), np.log10(self.x_limits[1])
                )

                self.axt.yaxis.set_major_locator(MultipleLocator(0.2))
                self.axt.yaxis.set_minor_locator(MultipleLocator(0.1))
                self.axt.set_xlabel("Period (s)", fontdict=self.font_dict)
                self.axt.set_ylabel("Tipper", fontdict=self.font_dict)

                # set th xaxis tick labels to invisible
                if self.plot_pt:
                    plt.setp(self.axt.xaxis.get_ticklabels(), visible=False)
                    self.axt.set_xlabel("")

    def _plot_pt(self):
        # ----plot phase tensor ellipse---------------------------------------
        if self.plot_pt:

            color_array = self.get_pt_color_array(self.pt)

            # -------------plot ellipses-----------------------------------
            self.cbax, self.cbpt, = plot_pt_lateral(
                self.axpt,
                self.pt,
                color_array,
                self.ellipse_properties,
                fig=self.fig,
            )

            # ----set axes properties-----------------------------------------------
            # --> set tick labels and limits
            self.axpt.set_xlim(
                np.log10(self.x_limits[0]), np.log10(self.x_limits[1])
            )

            tklabels, xticks = get_log_tick_labels(self.axpt)

            self.axpt.set_xticks(xticks)
            self.axpt.set_xticklabels(
                tklabels, fontdict={"size": self.font_size}
            )
            self.axpt.set_xlabel("Period (s)", fontdict=self.font_dict)

            # need to reset the x_limits caouse they get reset when calling
            # set_ticks for some reason
            self.axpt.set_xlim(
                np.log10(self.x_limits[0]), np.log10(self.x_limits[1])
            )
            self.axpt.grid(
                True,
                alpha=0.25,
                which="major",
                color=(0.25, 0.25, 0.25),
                lw=0.25,
            )

            plt.setp(self.axpt.get_yticklabels(), visible=False)

            self.cbpt.set_label(
                self.cb_label_dict[self.ellipse_colorby],
                fontdict={"size": self.font_size},
            )

    def plot(self):
        """
        plotResPhase(filename,fig_num) will plot the apparent resistivity and
        phase for a single station.

        """

        self._has_tipper()
        self._has_pt()

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

        eb_list, labels = self._plot_resistivity(
            self.axr, self.period, self.Z, mode="od"
        )
        self.ebxyr = eb_list[0]
        self.ebyxr = eb_list[1]
        self._plot_phase(self.axp, self.period, self.Z, mode="od")
        self._plot_tipper()
        self._plot_pt()

        # ===Plot the xx, yy components if desired==============================
        if self.plot_num == 2:
            self._plot_resistivity(self.axr2, self.period, self.Z, mode="d")
            self._plot_phase(self.axp2, self.period, self.Z, mode="d")
        # ===Plot the Determinant if desired==================================
        if self.plot_num == 3:
            self._plot_determinant()
        # make plot_title and show
        if self.plot_title is None:
            self.plot_title = self.station
        self.fig.suptitle(self.plot_title, fontdict=self.font_dict)

        # be sure to show
        if self.show_plot:
            plt.show()
