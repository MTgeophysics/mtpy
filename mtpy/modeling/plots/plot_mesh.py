# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 08:37:48 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mtpy.imaging.mtplot_tools import PlotBase
from mtpy.imaging.mtcolors import FixPointNormalize

# =============================================================================


class PlotMesh(PlotBase):
    def __init__(self, model_obj, **kwargs):
        self.model_obj = model_obj
        self.grid_lw = 0.5
        self.grid_color = (0.25, 0.25, 0.25)
        self.plot_station_id = False
        self.plot_topography = True
        self.z_limits = None

        super().__init__(**kwargs)

        if self.show_plot:
            self.plot()

    def _plot_topography(self):
        """
        Plot topography if asked

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not "topography" in self.model_obj.surface_dict.keys():

            topo = self.model_obj._get_topography_from_model()
            if topo is not None:
                self.model_obj.surface_dict["topography"] = topo

            else:
                self.logger.warning(
                    "Cannot find topography information, skipping"
                )
                return

        x, y = np.meshgrid(self.model_obj.grid_east, self.model_obj.grid_north)
        norm = FixPointNormalize(
            sealevel=0,
            vmax=np.round(self.model_obj.surface_dict["topography"].max(), -2),
            vmin=min(
                [
                    np.round(
                        self.model_obj.surface_dict["topography"].min(), -2
                    ),
                    0,
                ]
            ),
        )
        imgplot = self.ax1.pcolormesh(
            x,
            y,
            self.model_obj.surface_dict["topography"],
            cmap=plt.get_cmap("cut_terrain"),
            norm=norm,
        )
        divider = make_axes_locatable(self.ax1)
        cax = divider.append_axes("right", size="3%", pad=0.2)
        mycb = plt.colorbar(imgplot, cax=cax, use_gridspec=True)

        mycb.outline.set_linewidth(2)
        mycb.set_label(label="Elevation (m)", size=12)

    def plot(self):
        """
        Plot the mesh to show model grid

        Arguments:
        ----------

            **z_limits** : tuple (zmin,zmax)
                            plot min and max distances in meters for the
                            vertical direction.  If None, the z_limits is
                            set to the number of layers.  Z is positive down
                            *default* is None
        """

        self._set_subplot_params()
        self.fig = plt.figure(
            self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi
        )
        plt.clf()

        # make a rotation matrix to rotate data
        # cos_ang = np.cos(np.deg2rad(self.mesh_rotation_angle))
        # sin_ang = np.sin(np.deg2rad(self.mesh_rotation_angle))

        # turns out ModEM has not accomodated rotation of the grid, so for
        # now we will not rotate anything (angle=0.0)
        cos_ang = 1
        sin_ang = 0

        # --->plot map view
        self.ax1 = self.fig.add_subplot(1, 2, 1, aspect="equal")
        if self.plot_topography:
            self._plot_topography()

        x_min = self.model_obj.grid_east[self.model_obj.pad_east]
        x_max = self.model_obj.grid_east[-self.model_obj.pad_east]
        y_min = self.model_obj.grid_north[self.model_obj.pad_north]
        y_max = self.model_obj.grid_north[-self.model_obj.pad_north]

        # plot station locations
        if self.model_obj.station_locations is not None:
            plot_east = self.model_obj.station_locations.model_east
            plot_north = self.model_obj.station_locations.model_north

            # plot stations
            self.ax1.scatter(
                plot_east,
                plot_north,
                marker=self.marker,
                c=self.marker_color,
                s=self.marker_size,
            )
            if self.plot_station_id:
                for row in self.model_obj.station_locations.itertuples():
                    self.ax1.annotate(
                        row.station,
                        (row.model_east, row.model_north + 0.05),
                        ha="center",
                        va="baseline",
                        clip_on=True,
                    )
            if self.x_limits is None:
                x_min = plot_east.min() * 1.1
                x_max = plot_east.min() * 1.1

            if self.y_limits is None:
                y_min = plot_north.min() * 1.1
                y_max = plot_north.max() * 1.1

        east_line_xlist = []
        east_line_ylist = []
        north_min = self.model_obj.grid_north.min()
        north_max = self.model_obj.grid_north.max()
        for xx in self.model_obj.grid_east:
            east_line_xlist.extend(
                [
                    xx * cos_ang + north_min * sin_ang,
                    xx * cos_ang + north_max * sin_ang,
                ]
            )
            east_line_xlist.append(None)
            east_line_ylist.extend(
                [
                    -xx * sin_ang + north_min * cos_ang,
                    -xx * sin_ang + north_max * cos_ang,
                ]
            )
            east_line_ylist.append(None)
        self.ax1.plot(
            east_line_xlist,
            east_line_ylist,
            lw=self.grid_lw,
            color=self.grid_color,
        )

        north_line_xlist = []
        north_line_ylist = []
        east_max = self.model_obj.grid_east.max()
        east_min = self.model_obj.grid_east.min()
        for yy in self.model_obj.grid_north:
            north_line_xlist.extend(
                [
                    east_min * cos_ang + yy * sin_ang,
                    east_max * cos_ang + yy * sin_ang,
                ]
            )
            north_line_xlist.append(None)
            north_line_ylist.extend(
                [
                    -east_min * sin_ang + yy * cos_ang,
                    -east_max * sin_ang + yy * cos_ang,
                ]
            )
            north_line_ylist.append(None)
        self.ax1.plot(
            north_line_xlist,
            north_line_ylist,
            lw=self.grid_lw,
            color=self.grid_color,
        )

        if self.x_limits is not None:
            self.ax1.set_xlim(self.x_limits)
        else:
            self.ax1.set_xlim((x_min, x_max))

        if self.y_limits is not None:
            self.ax1.set_ylim(self.y_limits)
        else:
            self.ax1.set_ylim((y_min, y_max))

        self.ax1.set_ylabel("Northing (m)", fontdict=self.font_dict)
        self.ax1.set_xlabel("Easting (m)", fontdict=self.font_dict)

        # ---------------------------------------
        # plot depth view along the east direction
        self.ax2 = self.fig.add_subplot(
            1, 2, 2, aspect="auto", sharex=self.ax1
        )

        # plot the grid
        east_line_xlist = []
        east_line_ylist = []
        for xx in self.model_obj.grid_east:
            east_line_xlist.extend([xx, xx])
            east_line_xlist.append(None)
            east_line_ylist.extend(
                [self.model_obj.grid_z.min(), self.model_obj.grid_z.max()]
            )
            east_line_ylist.append(None)
        self.ax2.plot(
            east_line_xlist,
            east_line_ylist,
            lw=self.grid_lw,
            color=self.grid_color,
        )

        z_line_xlist = []
        z_line_ylist = []
        for zz in self.model_obj.grid_z:
            z_line_xlist.extend(
                [
                    self.model_obj.grid_east.min(),
                    self.model_obj.grid_east.max(),
                ]
            )
            z_line_xlist.append(None)
            z_line_ylist.extend([zz, zz])
            z_line_ylist.append(None)
        self.ax2.plot(
            z_line_xlist, z_line_ylist, lw=self.grid_lw, color=self.grid_color
        )

        # --> plot stations
        if self.model_obj.station_locations is not None:
            self.ax2.scatter(
                plot_east,
                self.model_obj.station_locations.model_elevation,
                marker=self.marker,
                c=self.marker_color,
                s=self.marker_size,
            )

        if self.z_limits is None:
            self.ax2.set_ylim(
                self.model_obj.z_target_depth,
                self.model_obj.grid_z.min() - 200,
            )
        else:
            self.ax2.set_ylim(self.z_limits)

        if self.x_limits is not None:
            self.ax2.set_xlim(self.x_limits)
        else:
            self.ax2.set_xlim((x_min, x_max))

        self.ax2.set_ylabel("Depth (m)", fontdict=self.font_dict)
        self.ax2.set_xlabel("Easting (m)", fontdict=self.font_dict)

        self.fig.tight_layout()

        plt.show()
