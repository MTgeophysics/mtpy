# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:50:39 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""

# ==============================================================================
#  Plot setting
# ==============================================================================
class PlotSettings(object):
    def __init__(self, **kwargs):

        self.fs = kwargs.pop("fs", 10)
        self.lw = kwargs.pop("lw", 1.5)
        self.ms = kwargs.pop("ms", 5)

        self.e_capthick = kwargs.pop("e_capthick", 1)
        self.e_capsize = kwargs.pop("e_capsize", 5)

        # color mode
        self.cted = kwargs.pop("cted", (0, 0, 0.75))
        self.ctmd = kwargs.pop("ctmd", (0.75, 0, 0))
        self.mted = kwargs.pop("mted", "s")
        self.mtmd = kwargs.pop("mtmd", "o")

        # color for occam2d model
        self.ctem = kwargs.pop("ctem", (0, 0.6, 0.3))
        self.ctmm = kwargs.pop("ctmm", (0.9, 0, 0.8))
        self.mtem = kwargs.pop("mtem", "+")
        self.mtmm = kwargs.pop("mtmm", "+")

        self.res_xx_limits = kwargs.pop("res_xx_limits", None)
        self.res_xy_limits = kwargs.pop("res_xy_limits", None)
        self.res_yx_limits = kwargs.pop("res_yx_limits", None)
        self.res_yy_limits = kwargs.pop("res_yy_limits", None)

        self.phase_xx_limits = kwargs.pop("phase_xx_limits", None)
        self.phase_xy_limits = kwargs.pop("phase_xy_limits", None)
        self.phase_yx_limits = kwargs.pop("phase_yx_limits", None)
        self.phase_yy_limits = kwargs.pop("phase_yy_limits", None)

        self.tipper_limits = kwargs.pop("tipper_limits", (-1.1, 1.1))

        self.subplot_wspace = kwargs.pop("subplot_wspace", 0.2)
        self.subplot_hspace = kwargs.pop("subplot_hspace", 0.0)
        self.subplot_right = kwargs.pop("subplot_right", 0.98)
        self.subplot_left = kwargs.pop("subplot_left", 0.08)
        self.subplot_top = kwargs.pop("subplot_top", 0.93)
        self.subplot_bottom = kwargs.pop("subplot_bottom", 0.08)

        self.z_err_increase = kwargs.pop("z_err_increase", 1.5)
        self.t_err_increase = kwargs.pop("t_err_increase", 0.10)

        self.legend_loc = kwargs.pop("legend_loc", "upper left")
        self.legend_pos = kwargs.pop("legend_pos", None)
        self.legend_marker_scale = kwargs.pop("legend_marker_scale", 1)
        self.legend_border_axes_pad = kwargs.pop("legend_border_axes_pad", 0.01)
        self.legend_label_spacing = kwargs.pop("legend_label_spacing", 0.07)
        self.legend_handle_text_pad = kwargs.pop("legend_handle_text_pad", 0.05)
        self.legend_border_pad = kwargs.pop("legend_border_pad", 0.05)

        self.ylabel_pad = 1.25
