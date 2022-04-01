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
from mtpy.imaging.mtplot_tools import PlotSettings
import mtpy.utils.exceptions as mtex

# ==============================================================================


class PlotStations(PlotSettings):
    """
    plot station locations in map view.

    Need to input one of the following lists:

   

    """

    def __init__(self, geo_df, **kwargs):

        super().__init__(**kwargs)
        self.gdf = geo_df

        # --> set plot properties
        self.plot_title = None
        self.station_id = None
        self.ref_point = (0, 0)

        self.map_epsg = 4326
        self.plot_names = True

        self.image_file = kwargs.pop("image_file", None)
        self.image_extent = kwargs.pop("image_extent", None)

        if self.image_file is not None:
            if self.image_extent is None:
                raise mtex.MTpyError_inputarguments(
                    "Need to input extents " + "of the image as" + "(x0, y0, x1, y1)"
                )
        self._set_subplot_parameters()
        # --> plot if desired
        if self.show_plot:
            self.plot()

    def _set_subplot_parameters(self):
        plt.rcParams["font.size"] = self.font_size
        plt.rcParams["figure.subplot.left"] = 0.09
        plt.rcParams["figure.subplot.right"] = 0.98
        plt.rcParams["figure.subplot.bottom"] = 0.09
        plt.rcParams["figure.subplot.top"] = 0.98

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

    def _get_xy(self):
        pass

    def plot(self):
        """
        plots the station locations

        """

        # if self.map_scale == "latlon":
        #     xlabel = "Longitude (deg)"
        #     ylabel = "Latitude (deg)"

        #     if self.xlimits is None:
        #         self._get_xlimits(self.df.longitude)
        #     if self.ylimits is None:
        #         self._get_ylimits(self.df.latitude)
        # elif self.map_scale == "eastnorth":
        #     xlabel = "Easting (m)"
        #     ylabel = "Northing (m)"
        # elif self.map_scale == "eastnorthkm":
        #     xlabel = "Easting (km)"
        #     ylabel = "Northing (km)"
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
            cx.add_basemap(
                gax, crs=self.gdf.crs.to_string(), source=cx.providers.USGS.USTopo
            )
        # set axis properties
        self.ax.set_xlabel("latitude", fontdict=self.font_dict)
        self.ax.set_ylabel("longitude", fontdict=self.font_dict)
        self.ax.grid(alpha=0.35, color=(0.25, 0.25, 0.25))
        # self.ax.set_xlim(self.xlimits)
        # self.ax.set_ylim(self.ylimits)

        plt.show()

    def save_plot(
        self,
        save_fn,
        file_format="pdf",
        orientation="portrait",
        fig_dpi=None,
        close_plot=True,
    ):
        """
        save_plot will save the figure to save_fn.

        Arguments:
        -----------

            **save_fn** : string
                          full path to save figure to, can be input as
                          * directory path -> the directory path to save to
                            in which the file will be saved as
                            save_fn/station_name_ResPhase.file_format

                          * full path -> file will be save to the given
                            path.  If you use this option then the format
                            will be assumed to be provided by the path

            **file_format** : [ pdf | eps | jpg | png | svg ]
                              file type of saved figure pdf,svg,eps...

            **orientation** : [ landscape | portrait ]
                              orientation in which the file will be saved
                              *default* is portrait

            **fig_dpi** : int
                          The resolution in dots-per-inch the file will be
                          saved.  If None then the dpi will be that at
                          which the figure was made.  I don't think that
                          it can be larger than dpi of the figure.

            **close_plot** : [ y | n ]
                             * 'y' will close the plot after saving.
                             * 'n' will leave plot open

        :Example: ::

            >>> # to save plot as jpg
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotPhaseTensorMaps(edilist,freqspot=10)
            >>> p1.save_plot(r'/home/MT', file_format='jpg')
            'Figure saved to /home/MT/PTMaps/PTmap_phimin_10Hz.jpg'

        """

        save_fn = Path(save_fn)
        if fig_dpi is None:
            fig_dpi = self.fig_dpi
        if not save_fn.is_dir():
            file_format = save_fn.suffix
        else:
            save_fn = save_fn.joinpath("station.png")
            self.fig.savefig(
                save_fn, dpi=fig_dpi, format=file_format, orientation=orientation
            )
            plt.close(self.fig)
        if close_plot:
            plt.close(self.fig)
        self.fig_fn = save_fn
        print("Saved figure to: " + self.fig_fn)

    def update_plot(self):
        """
        update any parameters that where changed using the built-in draw from
        canvas.

        Use this if you change an of the .fig or axes properties

        :Example: ::

            >>> # to change the grid lines to only be on the major ticks
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> [ax.grid(True, which='major') for ax in [p1.axr,p1.axp]]
            >>> p1.update_plot()

        """

        self.fig.canvas.draw()

    def redraw_plot(self):
        """
        use this function if you updated some attributes and want to re-plot.

        :Example: ::

            >>> # change the color and marker of the xy components
            >>> import mtpy.imaging.mtplottools as mtplot
            >>> p1 = mtplot.PlotResPhase(r'/home/MT/mt01.edi')
            >>> p1.xy_color = (.5,.5,.9)
            >>> p1.xy_marker = '*'
            >>> p1.redraw_plot()
        """

        plt.close(self.fig)
        self.plot()

    def __str__(self):
        return "Plot station locations"
