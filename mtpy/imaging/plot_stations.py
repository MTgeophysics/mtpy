# -*- coding: utf-8 -*-
"""
===============
PlotStations
===============

Plots station locations in map view.


Created on Fri Jun 07 18:20:00 2013

@author: jpeacock-pr
"""

# =============================================================================
import matplotlib.pyplot as plt
import numpy as np

try:
    import contextily as cx

    has_cx = True
except ModuleNotFoundError:
    has_cx = False
from mtpy.imaging.mtplot_tools import PlotBase
import mtpy.utils.exceptions as mtex

# ==============================================================================


class PlotStations(PlotBase):
    """
    plot station locations in map view.

    Uses contextily to get the basemap.
    See https://contextily.readthedocs.io/en/latest/index.html for more
    information about options.

    """

    def __init__(self, geo_df, **kwargs):

        # --> set plot properties
        self.plot_title = None
        self.station_id = None
        self.ref_point = (0, 0)

        self.map_epsg = 4326
        self.plot_names = True

        self.image_file = None
        self.image_extent = None
        self.pad = 0.005

        super().__init__(**kwargs)

        self._basename = "stations_map"
        self.gdf = geo_df

        self.cx_source = None
        self.cx_zoom = None
        if has_cx:
            self.cx_source = cx.providers.USGS.USTopo
        self._set_subplot_parameters()

        for key, value in kwargs.items():
            setattr(self, key, value)
        if self.image_file is not None:
            if self.image_extent is None:
                raise mtex.MTpyError_inputarguments(
                    "Need to input extents "
                    + "of the image as"
                    + "(x0, y0, x1, y1)"
                )
        # --> plot if desired
        if self.show_plot:
            self.plot()

    def _set_subplot_parameters(self):
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = self.subplot_left
        plt.rcParams["figure.subplot.right"] = self.subplot_right
        plt.rcParams["figure.subplot.bottom"] = self.subplot_bottom
        plt.rcParams["figure.subplot.top"] = self.subplot_top

    def _get_xlimits(self, x):
        if np.sign(x.min()) == -1:
            self.xlimits = (
                x.min() * (1 + self.pad),
                x.max() * (1 - self.pad),
            )
        else:
            self.xlimits = (
                x.min() * (1 - self.pad),
                x.max() * (1 + self.pad),
            )

    def _get_ylimits(self, y):
        if np.sign(y.min()) == -1:
            self.xlimits = (
                y.min() * (1 + self.pad),
                y.max() * (1 - self.pad),
            )
        else:
            self.xlimits = (
                y.min() * (1 - self.pad),
                y.max() * (1 + self.pad),
            )

    def plot(self):
        """

        :param cx_source: DESCRIPTION, defaults to cx.providers.USGS.USTopo
        :type cx_source: TYPE, optional
        :param cx_zoom: DESCRIPTION, defaults to None
        :type cx_zoom: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # add and axes
        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        # --> plot the background image if desired-----------------------
        if self.image_file is not None:
            im = plt.imread(self.image_file)
            self.ax.imshow(
                im, origin="lower", extent=self.image_extent, aspect="auto"
            )
        # plot stations
        gax = self.gdf.plot(
            ax=self.ax,
            marker=self.marker,
            color=self.marker_color,
            markersize=self.marker_size,
        )

        for x, y, label in zip(
            self.gdf.geometry.x, self.gdf.geometry.y, self.gdf.station
        ):
            gax.annotate(
                label,
                xy=(x, y),
                ha="center",
                va="baseline",
                xytext=(x, y + self.text_y_pad),
                rotation=self.text_angle,
                color=self.text_color,
                fontsize=self.text_size,
                fontweight=self.text_weight,
            )
        if self.image_file is None:
            if has_cx:
                try:
                    cx_kwargs = {
                        "crs": self.gdf.crs.to_string(),
                        "source": self.cx_source,
                    }
                    if self.cx_zoom is not None:
                        cx_kwargs["zoom"] = self.cx_zoom
                    cx.add_basemap(
                        gax,
                        **cx_kwargs,
                    )
                except Exception as error:
                    self.logger.warning(
                        f"Could not add base map because {error}"
                    )
        # set axis properties
        self.ax.set_xlabel("latitude", fontdict=self.font_dict)
        self.ax.set_ylabel("longitude", fontdict=self.font_dict)
        self.ax.grid(alpha=0.35, color=(0.35, 0.35, 0.35), lw=0.35)
        self.ax.set_xlim(self._get_xlimits(self.gdf.geometry.x))
        self.ax.set_ylim(self._get_ylimits(self.gdf.geometry.y))
        self.ax.set_axisbelow(True)

        plt.show()
