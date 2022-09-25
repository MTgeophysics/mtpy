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

import mtpy.imaging.mtcolors as mtcl
import mtpy.analysis.pt as mtpt
import mtpy.core.z as mtz

# ==============================================================================


class PlotPhaseTensorMaps(PlotBaseMaps):
    """
    Plots phase tensor ellipses in map view from a list of edi files

    """

    def __init__(self, tf_list, **kwargs):
        """
        Initialise the object
        :param kwargs: keyword-value pairs
        """
        super().__init__(**kwargs)

        self._rotation_angle = 0
        self.tf_list = tf_list

        # set the freq to plot
        self.plot_station = False
        self.plot_period = 1.0
        self.ftol = 0.1
        self.interpolate = True
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
    def map_scale(self):
        return self._map_scale

    @map_scale.setter
    def map_scale(self, map_scale):
        self._map_scale = map_scale

        if self._map_scale == "deg":
            self.xpad = 0.005
            self.ypad = 0.005
            self.ellipse_size = 0.005
            self.tickstrfmt = "%.3f"
            self.y_label = "Latitude (deg)"
            self.x_label = "Longitude (deg)"

        elif self._map_scale == "m":
            self.xpad = 1000
            self.ypad = 1000
            self.ellipse_size = 500
            self.tickstrfmt = "%.0f"
            self.x_label = "Easting (m)"
            self.y_label = "Northing (m)"

        elif self._map_scale == "km":
            self.xpad = 1
            self.ypad = 1
            self.ellipse_size = 0.500
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
        for tf in self.tf_list:
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
        z = self._get_interpolated_z(tf)
        z_err = self._get_interpolated_z_err(tf)

        new_z_obj = mtz.Z(z, z_err, freq=[1.0 / self.plot_period])
        new_z_obj.compute_resistivity_phase()
        pt_obj = mtpt.PhaseTensor(z_object=new_z_obj)

        new_t_obj = None
        if tf.t_interp_dict is not None:

            t = self._get_interpolated_t(tf)
            t_err = self._get_interpolated_t_err(tf)

            new_t_obj = mtz.Tipper(t, t_err, [1.0 / self.plot_period])

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

    def _get_patch(self, tf):
        """
        Get ellipse patch

        :param tf: DESCRIPTION
        :type tf: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        pt_obj, t_obj = self._get_pt(tf)

        # if map scale is lat lon set parameters
        if self.map_scale == "deg":
            plotx = tf.longitude - self.reference_point[0]
            ploty = tf.latitude - self.reference_point[1]
        # if map scale is in meters easting and northing
        elif self.map_scale in ["m", "km"]:
            tf.project_point_ll2utm(
                epsg=self.map_epsg, utm_zone=self.map_utm_zone
            )

            plotx = tf.east - self.reference_point[0]
            ploty = tf.north - self.reference_point[1]
            if self.map_scale in ["km"]:
                plotx /= 1000.0
                ploty /= 1000.0
        else:
            raise NameError("mapscale not recognized")
        # --> set local variables
        phimin = np.nan_to_num(pt_obj.phimin)
        phimax = np.nan_to_num(pt_obj.phimax)
        eangle = np.nan_to_num(pt_obj.azimuth)

        color_array = self.get_pt_color_array(pt_obj)
        bounds = np.arange(
            self.ellipse_range[0],
            self.ellipse_range[1] + self.ellipse_range[2],
            self.ellipse_range[2],
        )

        # --> get ellipse properties
        # if the ellipse size is not physically correct make it a dot
        if phimax == 0 or phimax > 100 or phimin == 0 or phimin > 100:
            eheight = 0.0000001
            ewidth = 0.0000001
        else:
            scaling = self.ellipse_size / phimax
            eheight = phimin * scaling
            ewidth = phimax * scaling
        # make an ellipse
        ellipd = patches.Ellipse(
            (plotx, ploty),
            width=ewidth,
            height=eheight,
            angle=90 - eangle,
            lw=self.lw,
        )

        # get ellipse color
        ellipd.set_facecolor(
            mtcl.get_plot_color(
                color_array[0],
                self.ellipse_colorby,
                self.ellipse_cmap,
                self.ellipse_range[0],
                self.ellipse_range[1],
                bounds=bounds,
            )
        )

        if t_obj is not None:
            if "r" in self.plot_tipper == "yri":

                if t_obj.mag_real[0] <= self.arrow_threshold:
                    txr = (
                        t_obj.mag_real[0]
                        * self.arrow_size
                        * np.sin(
                            (t_obj.angle_real[0]) * np.pi / 180
                            + self.arrow_direction * np.pi
                        )
                    )
                    tyr = (
                        t_obj.mag_real[0]
                        * self.arrow_size
                        * np.cos(
                            (t_obj.angle_real[0]) * np.pi / 180
                            + self.arrow_direction * np.pi
                        )
                    )

                    self.ax.arrow(
                        plotx,
                        ploty,
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
                    txi = (
                        t_obj.mag_imag[0]
                        * self.arrow_size
                        * np.sin(
                            (t_obj.angle_imag[0]) * np.pi / 180
                            + self.arrow_direction * np.pi
                        )
                    )
                    tyi = (
                        t_obj.mag_imag[0]
                        * self.arrow_size
                        * np.cos(
                            (t_obj.angle_imag[0]) * np.pi / 180
                            + self.arrow_direction * np.pi
                        )
                    )

                    self.ax.arrow(
                        plotx,
                        ploty,
                        txi,
                        tyi,
                        width=self.arrow_lw,
                        facecolor=self.arrow_color_imag,
                        edgecolor=self.arrow_color_imag,
                        length_includes_head=False,
                        head_width=self.arrow_head_width,
                        head_length=self.arrow_head_length,
                    )
        return ellipd, plotx, ploty

    def _add_colorbar(self):
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
        if self.ellipse_cmap in list(mtcl.cmapdict.keys()):
            cmap_input = mtcl.cmapdict[self.ellipse_cmap]
        else:
            cmap_input = mtcl.cm.get_cmap(self.ellipse_cmap)

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
        self.plot_xarr = np.zeros(len(self.tf_list))
        self.plot_yarr = np.zeros(len(self.tf_list))
        for index, tf in enumerate(self.tf_list):
            ellipse_patch, plot_x, plot_y = self._get_patch(tf)
            self.plot_xarr[index] = plot_x
            self.plot_yarr[index] = plot_y

            # ==> add ellipse to the plot
            self.ax.add_artist(ellipse_patch)

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
            self.plot_xarr[self.plot_xarr != 0.0].min() - self.x_pad,
            self.plot_xarr[self.plot_xarr != 0.0].max() + self.x_pad,
        )
        self.ax.set_ylim(
            self.plot_yarr[self.plot_yarr != 0.0].min() - self.y_pad,
            self.plot_yarr[self.plot_xarr != 0.0].max() + self.y_pad,
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

        # --> set title in period or freq
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
        self.ax.grid(True, alpha=0.3, which="both", color=(0.5, 0.5, 0.5))
        if self.minorticks_on:
            plt.minorticks_on()  # turn on minor ticks automatically

        self._add_colorbar()
