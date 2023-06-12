# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 15:20:43 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np

import matplotlib.colors as colors
import matplotlib.colorbar as mcb

from . import MTEllipse, MTArrows

import mtpy.imaging.mtcolors as mtcl

# =============================================================================
# ==============================================================================
# Plot settings
# ==============================================================================
class PlotSettings(MTArrows, MTEllipse):
    """
    Hold all the plot settings that one might need
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # figure properties:
        self.fig_num = 1
        self.fig_dpi = 150
        self.fig_size = None
        self.show_plot = True

        self.font_size = 7
        self.font_weight = "bold"
        self.marker_size = 2.5
        self.marker_lw = 0.75
        self.marker_color = "b"
        self.marker = "v"
        self.lw = 1
        self.plot_title = None

        # line styles:
        self.xy_ls = ":"
        self.yx_ls = ":"
        self.det_ls = ":"
        self.skew_ls = ":"
        self.strike_ls = ":"

        # marker styles:
        self.xy_marker = "s"
        self.yx_marker = "o"
        self.det_marker = "v"
        self.skew_marker = "d"
        self.strike_inv_marker = "v"
        self.strike_pt_marker = "^"
        self.strike_tip_marker = ">"

        # marker color styles:
        self.xy_color = (0.25, 0.35, 0.75)
        self.yx_color = (0.75, 0.25, 0.25)
        self.det_color = (0.25, 0.75, 0.25)
        self.skew_color = (0.85, 0.35, 0)
        self.strike_inv_color = (0.2, 0.2, 0.7)
        self.strike_pt_color = (0.7, 0.2, 0.2)
        self.strike_tip_color = (0.2, 0.7, 0.2)

        # marker face color styles:
        self.xy_mfc = (0.25, 0.35, 0.75)
        self.yx_mfc = (0.75, 0.25, 0.25)
        self.det_mfc = (0.25, 0.75, 0.25)
        self.skew_mfc = (0.85, 0.35, 0)
        self.strike_inv_mfc = (0.2, 0.2, 0.7)
        self.strike_pt_mfc = (0.7, 0.2, 0.2)
        self.strike_tip_mfc = (0.2, 0.7, 0.2)

        # plot limits
        self.x_limits = None
        self.y_limits = None
        self.res_limits = None
        self.phase_limits = None
        self.tipper_limits = None
        self.strike_limits = None
        self.skew_limits = None
        self.pt_limits = None

        # Show Plot
        self.show_plot = True
        self.plot_tipper = "n"
        self.plot_pt = False
        self.plot_strike = False
        self.plot_skew = False

        self.text_size = 7
        self.text_weight = "normal"
        self.text_color = "k"
        self.text_ha = "center"
        self.text_va = "baseline"
        self.text_angle = 0
        self.text_x_pad = 0
        self.text_y_pad = 0
        self.text_rotation = 0

        self.subplot_left = 0.09
        self.subplot_right = 0.9
        self.subplot_bottom = 0.09
        self.subplot_top = 0.98
        self.subplot_wspace = None
        self.subplot_hspace = None

        self.cb_orientation = "vertical"
        self.cb_position = None

        # Set class property values from kwargs and pop them
        for key, value in kwargs.items():
            setattr(self, key, value)
        self.cb_label_dict = {
            "phiminang": r"$\Phi_{min}$ (deg)",
            "phimin": r"$\Phi_{min}$ (deg)",
            "phimaxang": r"$\Phi_{max}$ (deg)",
            "phimax": r"$\Phi_{max}$ (deg)",
            "phidet": r"Det{$\Phi$} (deg)",
            "skew": r"Skew (deg)",
            "normalized_skew": r"Normalized Skew (deg)",
            "ellipticity": r"Ellipticity",
            "skew_seg": r"Skew (deg)",
            "normalized_skew_seg": r"Normalized Skew (deg)",
            "geometric_mean": r"$\sqrt{\Phi_{min} \cdot \Phi_{max}}$",
            "strike": r"Azimuth (deg)",
            "azimuth": r"Azimuth (deg)",
        }

    @property
    def period_label_dict(self):
        """
        log 10 labels

        :return: DESCRIPTION
        :rtype: TYPE

        """

        return dict([(ii, "$10^{" + str(ii) + "}$") for ii in range(-20, 21)])

    def set_period_limits(self, period):
        """
        set period limits

        :return: DESCRIPTION
        :rtype: TYPE

        """

        return (
            10 ** (np.floor(np.log10(period.min()))),
            10 ** (np.ceil(np.log10(period.max()))),
        )

    def set_resistivity_limits(self, resistivity, mode="od", scale="log"):
        """
        set resistivity limits

        :param resistivity: DESCRIPTION
        :type resistivity: TYPE
        :param mode: DESCRIPTION, defaults to "od"
        :type mode: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mode in ["od"]:
            try:
                nz_xy = np.nonzero(resistivity[:, 0, 1])
                nz_yx = np.nonzero(resistivity[:, 1, 0])
                limits = [
                    10
                    ** (
                        np.floor(
                            np.log10(
                                min(
                                    [
                                        np.nanmin(resistivity[nz_xy, 0, 1]),
                                        np.nanmin(resistivity[nz_yx, 1, 0]),
                                    ]
                                )
                            )
                        )
                    ),
                    10
                    ** (
                        np.ceil(
                            np.log10(
                                max(
                                    [
                                        np.nanmax(resistivity[nz_xy, 0, 1]),
                                        np.nanmax(resistivity[nz_yx, 1, 0]),
                                    ]
                                )
                            )
                        )
                    ),
                ]
            except (ValueError, TypeError):
                limits = [0.1, 10000]
        elif mode == "d":
            try:
                nz_xx = np.nonzero(resistivity[:, 0, 1])
                nz_yy = np.nonzero(resistivity[:, 1, 0])
                limits = [
                    10
                    ** (
                        np.floor(
                            np.log10(
                                min(
                                    [
                                        np.nanmin(resistivity[nz_xx, 0, 0]),
                                        np.nanmin(resistivity[nz_yy, 1, 1]),
                                    ]
                                )
                            )
                        )
                    ),
                    10
                    ** (
                        np.ceil(
                            np.log10(
                                max(
                                    [
                                        np.nanmax(resistivity[nz_xx, 0, 0]),
                                        np.nanmax(resistivity[nz_yy, 1, 1]),
                                    ]
                                )
                            )
                        )
                    ),
                ]
            except (ValueError, TypeError):
                limits = [0.1, 10000]
        elif mode in ["det", "det_only"]:
            try:
                nz = np.nonzero(resistivity)
                limits = [
                    10 ** (np.floor(np.log10(np.nanmin(resistivity[nz])))),
                    10 ** (np.ceil(np.log10(np.nanmax(resistivity[nz])))),
                ]
            except (ValueError, TypeError):
                limits = [0.1, 10000]
        if scale == "log":
            if limits[0] == 0:
                limits[0] = 0.1
        return limits

    def set_phase_limits(self, phase, mode="od"):
        if mode in ["od"]:
            try:
                nz_xy = np.nonzero(phase[:, 0, 1])
                nz_yx = np.nonzero(phase[:, 1, 0])

                ph_min = min(
                    [
                        np.nanmin(phase[nz_xy, 0, 1]),
                        np.nanmin(phase[nz_yx, 1, 0] + 180),
                    ]
                )
                if ph_min > 0:
                    ph_min = 0
                else:
                    ph_min = round(ph_min / 5) * 5
                ph_max = max(
                    [
                        np.nanmax(phase[:, 0, 1]),
                        np.nanmax(phase[:, 1, 0] + 180),
                    ]
                )
                if ph_max < 91:
                    ph_max = 89.9
                else:
                    ph_max = round(ph_max / 5) * 5
                return (ph_min, ph_max)
            except (ValueError, TypeError):
                return [0, 90]

        elif mode == "d":
            return (-180, 180)

        elif mode in ["det", "det_only"]:
            try:
                phase_det = np.linalg.det(phase)
                nz = np.nonzero(phase_det)
                phase_limits = [
                    np.amin(phase_det[nz]),
                    np.amax(phase_det[nz]),
                ]
                if phase_limits[0] < -180:
                    phase_limits[0] = -180
                if phase_limits[1] > 180:
                    phase_limits[1] = 180
                return phase_limits

            except (ValueError, TypeError):
                return [-180, 180]

    @property
    def xy_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.xy_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.xy_color,
            "color": self.xy_color,
            "ecolor": self.xy_color,
            "ls": self.xy_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def yx_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.yx_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.yx_color,
            "color": self.yx_color,
            "ecolor": self.yx_color,
            "ls": self.yx_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def det_error_bar_properties(self):
        """
        xy error bar properties for xy mode
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {
            "marker": self.det_marker,
            "ms": self.marker_size,
            "mew": self.lw,
            "mec": self.det_color,
            "color": self.det_color,
            "ecolor": self.det_color,
            "ls": self.det_ls,
            "lw": self.lw,
            "capsize": self.marker_size,
            "capthick": self.lw,
        }

    @property
    def font_dict(self):
        return {"size": self.font_size + 2, "weight": self.font_weight}

    def make_pt_cb(self, ax):

        cmap = mtcl.cmapdict[self.ellipse_cmap]
        if "seg" in self.ellipse_cmap:
            # normalize the colors
            norms = colors.BoundaryNorm(self.ellipse_cmap_bounds, cmap.N)

            # make the colorbar
            cb = mcb.Colorbar(
                ax,
                cmap=cmap,
                norm=norms,
                orientation=self.cb_orientation,
                ticks=self.ellipse_cmap_bounds[1:-1],
            )
        else:
            cb = mcb.Colorbar(
                ax,
                cmap=cmap,
                norm=colors.Normalize(
                    vmin=self.ellipse_range[0],
                    vmax=self.ellipse_range[1],
                ),
                orientation=self.cb_orientation,
            )
        # label the color bar accordingly
        cb.set_label(
            self.cb_label_dict[self.ellipse_colorby],
            fontdict=self.font_dict,
        )

        # place the label in the correct location
        if self.cb_orientation == "horizontal":
            cb.ax.xaxis.set_label_position("top")
            cb.ax.xaxis.set_label_coords(0.5, 1.3)
        elif self.cb_orientation == "vertical":
            cb.ax.yaxis.set_label_position("right")
            cb.ax.yaxis.set_label_coords(1.25, 0.5)
            cb.ax.yaxis.tick_left()
            cb.ax.tick_params(axis="y", direction="in")

        return cb

    @property
    def arrow_real_properties(self):
        return {
            "lw": self.arrow_lw,
            "facecolor": self.arrow_color_real,
            "edgecolor": self.arrow_color_real,
            "head_width": self.arrow_head_width,
            "head_length": self.arrow_head_length,
            "length_includes_head": False,
        }

    @property
    def arrow_imag_properties(self):
        return {
            "lw": self.arrow_lw,
            "facecolor": self.arrow_color_imag,
            "edgecolor": self.arrow_color_imag,
            "head_width": self.arrow_head_width,
            "head_length": self.arrow_head_length,
            "length_includes_head": False,
        }

    @property
    def text_dict(self):
        return {
            "size": self.text_size,
            "weight": self.text_weight,
            "rotation": self.text_angle,
            "color": self.text_color,
        }
