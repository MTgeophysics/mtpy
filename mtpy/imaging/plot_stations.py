# -*- coding: utf-8 -*-
"""
===============
PlotStations
===============

Plots station locations in map view.


Created on Fri Jun 07 18:20:00 2013

@author: jpeacock-pr
"""

# ==============================================================================
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import contextily as cx
from mtpy.imaging.mtplot_tools import PlotBase
import mtpy.utils.exceptions as mtex

# ==============================================================================


class PlotStations(PlotBase):
    """
    plot station locations in map view.

    Need to input one of the following lists:

   

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

        super().__init__(**kwargs)

        self._basename = "stations_map"
        self.gdf = geo_df

        self._set_subplot_parameters()

        for key, value in kwargs.items():
            setattr(self, key, value)
        if self.image_file is not None:
            if self.image_extent is None:
                raise mtex.MTpyError_inputarguments(
                    "Need to input extents " + "of the image as" + "(x0, y0, x1, y1)"
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
                x.min() * 1.002,
                x.max() * 0.998,
            )
        else:
            self.xlimits = (
                x.min() * 0.998,
                x.max() * 1.002,
            )

    def _get_ylimits(self, y):
        if np.sign(y.min()) == -1:
            self.xlimits = (
                y.min() * 1.002,
                y.max() * 0.998,
            )
        else:
            self.xlimits = (
                y.min() * 0.998,
                y.max() * 1.002,
            )

    def plot(self):
        """
        plots the station locations

        """

        # make a figure instance
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)

        # add and axes
        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        # --> plot the background image if desired-----------------------
        if self.image_file is not None:
            im = plt.imread(self.image_file)
            self.ax.imshow(im, origin="lower", extent=self.image_extent, aspect="auto")
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
            )
        if self.image_file is None:
            try:
                cx.add_basemap(
                    gax, crs=self.gdf.crs.to_string(), source=cx.providers.USGS.USTopo
                )
            except Exception as error:
                self._logger.warning(f"Could not add base map because {error}")
        # set axis properties
        self.ax.set_xlabel("latitude", fontdict=self.font_dict)
        self.ax.set_ylabel("longitude", fontdict=self.font_dict)
        self.ax.grid(alpha=0.35, color=(0.25, 0.25, 0.25))
        self.ax.set_xlim(self._get_xlimits(self.gdf.geometry.x))
        self.ax.set_ylim(self._get_ylimits(self.gdf.geometry.x))

        plt.show()
