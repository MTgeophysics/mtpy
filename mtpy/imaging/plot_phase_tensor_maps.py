# -*- coding: utf-8 -*-
"""
Plot phase tensor map in Lat-Lon Coordinate System

Revision History:
    Created by @author: jpeacock-pr on Thu May 30 18:20:04 2013

    Modified by Fei.Zhang@ga.gov.au 2017-03:

    brenainn.moushall 26-03-2020 15:07:14 AEDT:
        Add plotting of geotiff as basemap background.
        
    Updated 2022 by J. Peacock to work with v2
        
        - Using rasterio to plot geotiffs
        - factorized
        - using interp function for faster plotting.
    
"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colorbar as mcb
import matplotlib.colors as colors

try:
    import contextily as cx

    has_cx = True
except ModuleNotFoundError:
    has_cx = False
from mtpy.imaging.mtplot_tools import PlotBaseMaps, add_raster
from mtpy.core import Tipper
from mtpy.imaging import mtcolors
from mtpy.core.transfer_function import PhaseTensor


# ==============================================================================


class PlotPhaseTensorMaps(PlotBaseMaps):
    """
    Plots phase tensor ellipses in map view from a list of edi files

    """

    def __init__(self, mt_data, **kwargs):
        """
        Initialise the object
        :param kwargs: keyword-value pairs
        """
        super().__init__(**kwargs)

        self._rotation_angle = 0
        self.mt_data = mt_data

        # set the frequency to plot
        self.plot_station = False
        self.plot_period = 1.0

        self.pt_type = "wedges"
        self.skew_cmap = "mt_seg_bl2wh2rd"
        self.phase_limits = (0, 90)
        self.skew_limits = (-9, 9)
        self.skew_step = 3
        self.skew_lw = 2.5
        self.ellipse_alpha = 0.85
        self.wedge_width = 7.0

        # read in map scale
        self.map_scale = "deg"
        self.map_utm_zone = None
        self.map_epsg = None

        self.minorticks_on = True
        self.x_pad = 0.01
        self.y_pad = 0.01

        # arrow legend
        self.arrow_legend_position = "lower right"
        self.arrow_legend_xborderpad = 0.2
        self.arrow_legend_yborderpad = 0.2
        self.arrow_legend_fontpad = 0.05

        self.ref_ax_loc = (0.85, 0.1, 0.1, 0.1)

        # set a central reference point
        self.reference_point = (0, 0)

        self.cx_source = None
        self.cx_zoom = None
        if has_cx:
            self.cx_source = cx.providers.USGS.USTopo
        # station labels
        self.station_id = (0, 2)
        self.station_pad = 0.0005

        self.arrow_legend_fontdict = {"size": self.font_size, "weight": "bold"}
        self.station_font_dict = {"size": self.font_size, "weight": "bold"}

        for key, value in kwargs.items():
            setattr(self, key, value)
        # --> plot if desired ------------------------
        if self.show_plot:
            self.plot()

    @property
    def skew_cmap_bounds(self):
        """skew bounds for a segmented colorbar"""
        return np.arange(
            self.skew_limits[0],
            self.skew_limits[1] + self.skew_step,
            self.skew_step,
        )

    @property
    def map_scale(self):
        return self._map_scale

    @map_scale.setter
    def map_scale(self, map_scale):
        self._map_scale = map_scale

        if self._map_scale == "deg":
            self.xpad = 0.005
            self.ypad = 0.005
            self.ellipse_size = 0.005
            self.arrow_size = 0.005
            self.arrow_head_length = 0.0025
            self.arrow_head_width = 0.0035
            self.arrow_lw = 0.00075

            self.tickstrfmt = "%.3f"
            self.y_label = "Latitude (deg)"
            self.x_label = "Longitude (deg)"
        elif self._map_scale == "m":
            self.xpad = 1000
            self.ypad = 1000
            self.ellipse_size = 500
            self.arrow_size = 500
            self.arrow_head_length = 250
            self.arrow_head_width = 350
            self.arrow_lw = 50
            self.tickstrfmt = "%.0f"
            self.x_label = "Easting (m)"
            self.y_label = "Northing (m)"
        elif self._map_scale == "km":
            self.xpad = 1
            self.ypad = 1
            self.ellipse_size = 0.500
            self.arrow_size = 0.5
            self.arrow_head_length = 0.25
            self.arrow_head_width = 0.35
            self.arrow_lw = 0.075
            self.tickstrfmt = "%.0f"
            self.x_label = "Easting (km)"
            self.y_label = "Northing (km)"
        else:
            raise ValueError(f"map scale {map_scale} is not supported.")

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

    def _get_pt(self, tf):
        """
        Get phase tensor object from TF object

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        try:
            z = self._get_interpolated_z(tf)
            z_error = self._get_interpolated_z_error(tf)

            pt_obj = PhaseTensor(z=z, z_error=z_error)
        except ValueError as error:
            self.logger.warning(
                f"Could not estimate phase tensor for {tf.station} at period "
                f"{self.plot_period} s."
            )
            self.logger.error(error)
            pt_obj = None
        new_t_obj = None
        if tf.has_tipper():
            try:
                t = self._get_interpolated_t(tf)
                t_err = self._get_interpolated_t_err(tf)
                if (t != 0).all():
                    new_t_obj = Tipper(t, t_err, [1.0 / self.plot_period])
            except ValueError:
                self.logger.warning(
                    f"Could not estimate tipper for {tf.station} at period "
                    f"{self.plot_period} s."
                )
        return pt_obj, new_t_obj

    def _get_tick_format(self):
        """

        :return: DESCRIPTION
        :rtype: TYPE

        """

        # set tick parameters depending on the mapscale
        if self.map_scale == "deg":
            self.tickstrfmt = "%.2f"
        elif self.map_scale == "m" or self.map_scale == "km":
            self.tickstrfmt = "%.0f"

    def _set_axis_labels(self):
        # --> set axes properties depending on map scale------------------------
        if self.map_scale == "deg":
            self.ax.set_xlabel(
                "Longitude", fontsize=self.font_size, fontweight="bold"  # +2,
            )
            self.ax.set_ylabel(
                "Latitude", fontsize=self.font_size, fontweight="bold"  # +2,
            )
        elif self.map_scale == "m":
            self.ax.set_xlabel(
                "Easting (m)",
                fontsize=self.font_size,
                fontweight="bold",  # +2,
            )
            self.ax.set_ylabel(
                "Northing (m)",
                fontsize=self.font_size,
                fontweight="bold",  # +2,
            )
        elif self.map_scale == "km":
            self.ax.set_xlabel(
                "Easting (km)",
                fontsize=self.font_size,
                fontweight="bold",  # +2,
            )
            self.ax.set_ylabel(
                "Northing (km)",
                fontsize=self.font_size,
                fontweight="bold",  # +2,
            )

    def _get_location(self, tf):
        """
        Get plot_x and plot_y in appropriate units

        :param tf: DESCRIPTION
        :type tf: TYPE
        :raises NameError: DESCRIPTION
        :return: DESCRIPTION
        :rtype: TYPE

        """
        # if map scale is lat lon set parameters
        if self.map_scale == "deg":
            plot_x = tf.longitude - self.reference_point[0]
            plot_y = tf.latitude - self.reference_point[1]
        # if map scale is in meters easting and northing
        elif self.map_scale in ["m", "km"]:
            tf.project_point_ll2utm(
                epsg=self.map_epsg, utm_zone=self.map_utm_zone
            )

            plot_x = tf.east - self.reference_point[0]
            plot_y = tf.north - self.reference_point[1]
            if self.map_scale in ["km"]:
                plot_x /= 1000.0
                plot_y /= 1000.0
        else:
            raise NameError("mapscale not recognized")
        return plot_x, plot_y

    def _get_patch_ellipse(self, tf):
        """
        Get ellipse patch

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        pt_obj, t_obj = self._get_pt(tf)

        if pt_obj is None and t_obj is None:
            return (0, 0)
        if pt_obj is None:
            has_ellipse = False
        if pt_obj is not None:
            plot_x, plot_y = self._get_location(tf)

            # --> set local variables
            phimin = np.nan_to_num(pt_obj.phimin)[0]
            phimax = np.nan_to_num(pt_obj.phimax)[0]
            eangle = np.nan_to_num(pt_obj.azimuth)[0]

            color_array = self.get_pt_color_array(pt_obj)
            bounds = np.arange(
                self.ellipse_range[0],
                self.ellipse_range[1] + self.ellipse_range[2],
                self.ellipse_range[2],
            )

            has_ellipse = True
            # --> get ellipse properties
            # if the ellipse size is not physically correct make it a dot
            if phimax == 0 or phimax > 100 or phimin == 0 or phimin > 100:
                eheight = 0.0000001
                ewidth = 0.0000001
                has_ellipse = False
            else:
                scaling = self.ellipse_size / phimax
                eheight = phimin * scaling
                ewidth = phimax * scaling
                # make an ellipse
                ellipd = patches.Ellipse(
                    (plot_x, plot_y),
                    width=ewidth,
                    height=eheight,
                    angle=90 - eangle,
                    lw=self.lw,
                )

                # get ellipse color
                ellipd.set_facecolor(
                    mtcolors.get_plot_color(
                        color_array[0],
                        self.ellipse_colorby,
                        self.ellipse_cmap,
                        self.ellipse_range[0],
                        self.ellipse_range[1],
                        bounds=bounds,
                    )
                )

                self.ax.add_artist(ellipd)
        has_tipper = self._get_tipper_patch(plot_x, plot_y, t_obj)

        if has_ellipse or has_tipper:
            return plot_x, plot_y
        else:
            return (0, 0)

    def _get_tipper_patch(self, plot_x, plot_y, t_obj):
        """

        :param t_obj: DESCRIPTION
        :type t_obj: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        has_tipper = False
        if t_obj is not None:
            if "r" in self.plot_tipper == "yri":
                if t_obj.mag_real[0] <= self.arrow_threshold:
                    has_tipper = True
                    txr = (
                        t_obj.mag_real[0]
                        * self.arrow_size
                        * np.sin(
                            np.deg2rad(t_obj.angle_real[0])
                            + self.arrow_direction * np.pi
                        )
                    )
                    tyr = (
                        t_obj.mag_real[0]
                        * self.arrow_size
                        * np.cos(
                            np.deg2rad(t_obj.angle_real[0])
                            + self.arrow_direction * np.pi
                        )
                    )

                    self.ax.arrow(
                        plot_x,
                        plot_y,
                        txr,
                        tyr,
                        width=self.arrow_lw,
                        facecolor=self.arrow_color_real,
                        edgecolor=self.arrow_color_real,
                        length_includes_head=False,
                        head_width=self.arrow_head_width,
                        head_length=self.arrow_head_length,
                    )
                else:
                    pass
            # plot imaginary tipper
            if "i" in self.plot_tipper:
                if t_obj.mag_imag[0] <= self.arrow_threshold:
                    has_tipper = True
                    txi = (
                        t_obj.mag_imag[0]
                        * self.arrow_size
                        * np.sin(
                            np.deg2rad(t_obj.angle_imag[0])
                            + self.arrow_direction * np.pi
                        )
                    )
                    tyi = (
                        t_obj.mag_imag[0]
                        * self.arrow_size
                        * np.cos(
                            np.deg2rad(t_obj.angle_imag[0])
                            + self.arrow_direction * np.pi
                        )
                    )

                    self.ax.arrow(
                        plot_x,
                        plot_y,
                        txi,
                        tyi,
                        width=self.arrow_lw,
                        facecolor=self.arrow_color_imag,
                        edgecolor=self.arrow_color_imag,
                        length_includes_head=False,
                        head_width=self.arrow_head_width,
                        head_length=self.arrow_head_length,
                    )
        return has_tipper

    def _get_patch_wedges(self, tf):
        """
        Add patches to plot
        """

        pt_obj, t_obj = self._get_pt(tf)
        plot_x, plot_y = self._get_location(tf)
        has_ellipse = False
        if pt_obj is not None:

            # --> set local variables
            phimin = np.nan_to_num(pt_obj.phimin)[0]
            phimax = np.nan_to_num(pt_obj.phimax)[0]
            eangle = np.nan_to_num(pt_obj.azimuth)[0]

            has_ellipse = True
            w1 = patches.Wedge(
                (plot_x, plot_y),
                self.ellipse_size,
                90 - eangle - self.wedge_width,
                90 - eangle + self.wedge_width,
                color=mtcolors.get_plot_color(
                    phimax,
                    "phimax",
                    self.ellipse_cmap,
                    self.phase_limits[0],
                    self.phase_limits[1],
                ),
            )
            w2 = patches.Wedge(
                (plot_x, plot_y),
                self.ellipse_size,
                270 - eangle - self.wedge_width,
                270 - eangle + self.wedge_width,
                color=mtcolors.get_plot_color(
                    phimax,
                    "phimax",
                    self.ellipse_cmap,
                    self.phase_limits[0],
                    self.phase_limits[1],
                ),
            )

            w3 = patches.Wedge(
                (plot_x, plot_y),
                self.ellipse_size * phimin / phimax,
                -1 * eangle - self.wedge_width,
                -1 * eangle + self.wedge_width,
                color=mtcolors.get_plot_color(
                    phimin,
                    "phimin",
                    self.ellipse_cmap,
                    self.phase_limits[0],
                    self.phase_limits[1],
                ),
            )
            w4 = patches.Wedge(
                (plot_x, plot_y),
                self.ellipse_size * phimin / phimax,
                180 - eangle - self.wedge_width,
                180 - eangle + self.wedge_width,
                color=mtcolors.get_plot_color(
                    phimin,
                    "phimin",
                    self.ellipse_cmap,
                    self.phase_limits[0],
                    self.phase_limits[1],
                ),
            )

            # make an ellipse
            e1 = patches.Ellipse(
                (plot_x, plot_y),
                width=2 * self.ellipse_size,
                height=2 * self.ellipse_size * phimin / phimax,
                angle=90 - eangle,
            )

            # geometric mean of phimin and phimax
            gm = np.sqrt(abs(phimin) * abs(phimax))
            e1.set_facecolor(
                mtcolors.get_plot_color(
                    gm,
                    "geometric_mean",
                    self.ellipse_cmap,
                    self.phase_limits[0],
                    self.phase_limits[1],
                )
            )
            e1.set_edgecolor(
                mtcolors.get_plot_color(
                    np.nan_to_num(pt_obj.skew),
                    "skew_seg",
                    self.skew_cmap,
                    self.skew_limits[0],
                    self.skew_limits[1],
                    self.skew_cmap_bounds,
                )
            )
            e1.set_linewidth(self.skew_lw)
            e1.set_alpha(self.ellipse_alpha)
            ### add patches
            for patch in [e1, w1, w2, w3, w4]:
                self.ax.add_patch(patch)
        has_tipper = self._get_tipper_patch(plot_x, plot_y, t_obj)

        if has_ellipse or has_tipper:
            return plot_x, plot_y
        else:
            return (0, 0)

    def _add_colorbar_ellipse(self):
        """
        Add phase tensor color bar

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.cb_position is None:
            self.ax2, kw = mcb.make_axes(
                self.ax, orientation=self.cb_orientation, shrink=0.35
            )
        else:
            self.ax2 = self.fig.add_axes(self.cb_position)
        # make the colorbar
        if self.ellipse_cmap in list(mtcolors.cmapdict.keys()):
            cmap_input = mtcolors.cmapdict[self.ellipse_cmap]
        else:
            cmap_input = mtcolors.cm.get_cmap(self.ellipse_cmap)
        if "seg" in self.ellipse_cmap:
            norms = colors.BoundaryNorm(self.ellipse_cmap_bounds, cmap_input.N)
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=cmap_input,
                norm=norms,
                orientation=self.cb_orientation,
                ticks=self.ellipse_cmap_bounds,
            )
        else:
            self.cb = mcb.ColorbarBase(
                self.ax2,
                cmap=cmap_input,
                norm=colors.Normalize(
                    vmin=self.ellipse_range[0], vmax=self.ellipse_range[1]
                ),
                orientation=self.cb_orientation,
            )
        # label the color bar accordingly
        self.cb.set_label(
            self.cb_label_dict[self.ellipse_colorby],
            fontdict={"size": self.font_size, "weight": "bold"},
        )

        # place the label in the correct location
        if self.cb_orientation == "horizontal":
            self.cb.ax.xaxis.set_label_position("top")
            self.cb.ax.xaxis.set_label_coords(0.5, 1.3)
        elif self.cb_orientation == "vertical":
            self.cb.ax.yaxis.set_label_position("right")
            self.cb.ax.yaxis.set_label_coords(1.25, 0.5)
            self.cb.ax.yaxis.tick_left()
            self.cb.ax.tick_params(axis="y", direction="in")

    def _add_colorbar_wedges(self):
        """
        Add phase tensor color bar

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.cb_orientation == "vertical":
            self.ax2, kw = mcb.make_axes(
                self.ax,
                orientation=self.cb_orientation,
                shrink=0.15,
                anchor=(0, 0.60),
            )

            pos = self.ax2.get_position()
            x3_pos = (pos.x0, pos.y0 - 0.175, pos.width, pos.height)

            self.ax3 = self.fig.add_axes(x3_pos)
        elif self.cb_orientation == "horizontal":
            self.ax2, kw = mcb.make_axes(
                self.ax,
                orientation=self.cb_orientation,
                shrink=0.15,
                anchor=(0.4, 0),
            )
            pos = self.ax2.get_position()
            x3_pos = (pos.x0 + 0.175, pos.y0, pos.width, pos.height)

            self.ax3 = self.fig.add_axes(x3_pos)
        # make the colorbar
        for key, ax in zip(["ellipse", "skew"], [self.ax2, self.ax3]):
            dict_key = f"{key}_cmap"
            cmap = getattr(self, dict_key)
            cmap_input = mtcolors.cm.get_cmap(cmap)
            if "seg" in cmap and "ellipse" in key:

                norms = colors.BoundaryNorm(
                    self.ellipse_cmap_bounds, cmap_input.N
                )
                self.cb = mcb.ColorbarBase(
                    ax,
                    cmap=cmap_input,
                    norm=norms,
                    orientation=self.cb_orientation,
                    ticks=self.ellipse_cmap_bounds,
                )
            elif "skew" in key:
                norms = colors.BoundaryNorm(
                    self.skew_cmap_bounds, cmap_input.N
                )
                cb = mcb.ColorbarBase(
                    ax,
                    cmap=cmap_input,
                    norm=norms,
                    orientation=self.cb_orientation,
                    ticks=self.skew_cmap_bounds,
                )
            else:
                cb = mcb.ColorbarBase(
                    ax,
                    cmap=cmap_input,
                    norm=colors.Normalize(
                        vmin=self.ellipse_range[0], vmax=self.ellipse_range[1]
                    ),
                    orientation=self.cb_orientation,
                )
            # label the color bar accordingly
            if key == "ellipse":
                label = "Phase (deg)"
            else:
                label = "Skew (deg)"
            cb.set_label(
                label,
                fontdict={"size": self.font_size, "weight": "bold"},
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

    # -----------------------------------------------
    # The main plot method for this module
    # -----------------------------------------------
    def plot(
        self,
        fig=None,
        save_path=None,
        show=True,
        raster_file=None,
        raster_kwargs={},
    ):
        """
        Plots the phase tensor map.
        :param fig: optional figure object
        :param save_path: path to folder for saving plots
        :param show: show plots if True
        :param raster_dict: dictionary containing information for a raster
         to be plotted below the phase tensors.

        """
        self._set_subplot_params()

        # make figure instance
        self.fig = plt.figure(
            self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi
        )

        # clear the figure if there is already one up
        plt.clf()

        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        self._get_tick_format()

        # make some empty arrays
        self.plot_xarr = np.zeros(len(self.mt_data))
        self.plot_yarr = np.zeros(len(self.mt_data))
        for index, tf in enumerate(self.mt_data.values()):
            if self.pt_type == "ellipses":
                plot_x, plot_y = self._get_patch_ellipse(tf)
            elif self.pt_type == "wedges":
                plot_x, plot_y = self._get_patch_wedges(tf)
            else:
                raise ValueError(
                    f"{self.pt_type} not supported. Use ['ellipses' | 'wedges']"
                )
            self.plot_xarr[index] = plot_x
            self.plot_yarr[index] = plot_y

            # ------------Plot station name------------------------------
            if self.plot_station:
                self.ax.text(
                    plot_x,
                    plot_y + self.station_pad,
                    tf.station[self.station_id[0] : self.station_id[1]],
                    horizontalalignment="center",
                    verticalalignment="baseline",
                    fontdict=self.station_font_dict,
                )
        self._set_axis_labels()
        # --> set plot limits
        #    need to exclude zero values from the calculation of min/max!!!!
        self.ax.set_xlim(
            self.plot_xarr[np.nonzero(self.plot_xarr)].min() - self.x_pad,
            self.plot_xarr[np.nonzero(self.plot_xarr)].max() + self.x_pad,
        )
        self.ax.set_ylim(
            self.plot_yarr[np.nonzero(self.plot_yarr)].min() - self.y_pad,
            self.plot_yarr[np.nonzero(self.plot_xarr)].max() + self.y_pad,
        )

        # --> set tick label format
        self.ax.xaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        self.ax.yaxis.set_major_formatter(FormatStrFormatter(self.tickstrfmt))
        plt.setp(self.ax.get_xticklabels(), rotation=45)

        ## rasterio for plotting geotiffs or other geophysical data.
        if raster_file is not None:
            self.raster_ax, self.raster_cb = add_raster(
                self.ax, raster_file, **raster_kwargs
            )
        else:
            if has_cx:
                try:
                    cx_kwargs = {"source": self.cx_source, "crs": "EPSG:4326"}
                    if self.cx_zoom is not None:
                        cx_kwargs["zoom"] = self.cx_zoom
                    cx.add_basemap(
                        self.ax,
                        **cx_kwargs,
                    )
                except Exception as error:
                    self.logger.warning(
                        f"Could not add base map because {error}"
                    )
        # --> set title in period or frequency
        titlefreq = "{0:.5g} (s)".format(self.plot_period)

        if not self.plot_title:
            self.ax.set_title(
                "Phase Tensor Map for " + titlefreq,
                fontsize=self.font_size + 2,
                fontweight="bold",
            )
        else:
            self.ax.set_title(
                self.plot_title + titlefreq,
                fontsize=self.font_size + 2,
                fontweight="bold",
            )
        # make a grid with color lines
        self.ax.grid(
            True, alpha=0.3, which="major", color=(0.5, 0.5, 0.5), lw=0.75
        )
        self.ax.grid(
            True,
            alpha=0.3,
            which="minor",
            color=(0.5, 0.5, 0.5),
            lw=0.5,
            ls=":",
        )
        if self.minorticks_on:
            plt.minorticks_on()  # turn on minor ticks automatically
        self.ax.set_axisbelow(True)

        if self.pt_type == "ellipses":
            self._add_colorbar_ellipse()
        elif self.pt_type == "wedges":
            self._add_colorbar_wedges()
