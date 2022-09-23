#!/bin/env python
"""
Description:
    Plots resistivity and phase maps for a given frequency
   
References:
 
CreationDate:   4/19/18
Developer:      rakib.hassan@ga.gov.au
 
Revision History:
    LastUpdate:     4/19/18   RH
                    2022-09 JP
    

"""
# =============================================================================
# Imports
# =============================================================================

import numpy as np

import matplotlib.pyplot as plt

import matplotlib.colors as colors
from matplotlib import ticker
import matplotlib.tri as tri

from mtpy.core.z import Z
from mtpy.imaging.mtplot_tools import PlotBase, in_hull, period_label_dict

from scipy.spatial import cKDTree
from scipy import interpolate

# =============================================================================


class PlotResPhaseMaps(PlotBase):
    """
    Plots apparent resistivity and phase in map view from a list of edi files

    Arguments:
    -------------

        **fn_list** : list of strings
                          full paths to .edi files to plot

        **fig_size** : tuple or list (x, y) in inches
                      dimensions of the figure box in inches, this is a default
                      unit of matplotlib.  You can use this so make the plot
                      fit the figure box to minimize spaces from the plot axes
                      to the figure box.  *default* is [8, 8]

        **mapscale** : [ 'deg' | 'm' | 'km' ]
                       Scale of the map coordinates.

                       * 'deg' --> degrees in latitude and longitude

                       * 'm' --> meters for easting and northing

                       * 'km' --> kilometers for easting and northing

        **plot_yn** : [ 'y' | 'n' ]
                      *'y' to plot on creating an instance

                      *'n' to not plot on creating an instance

        **title** : string
                    figure title

        **dpi** : int
                  dots per inch of the resolution. *default* is 300

        **font_size** : float
                        size of the font that labels the plot, 2 will be added
                        to this number for the axis labels.
    """

    def __init__(self, tf_list, **kwargs):
        """
        Initialise the object
        :param kwargs: keyword-value pairs
        """
        super().__init__(**kwargs)

        self.tf_list = tf_list

        # read in map scale
        self.map_units = "deg"
        self.scale = 1
        self.res_cmap = "rainbow_r"
        self.phase_cmap = "rainbow"
        self.plot_period = 1

        self.cell_size = 0.002
        self.n_padding_cells = 5

        self.interpolation_method = "cubic"

        self.plot_xx = False
        self.plot_xy = True
        self.plot_yx = True
        self.plot_yy = False
        self.plot_det = False

        self.plot_resistivity = True
        self.plot_phase = True

        self.plot_stations = True

        self.marker_color = "k"
        self.marker_size = 10

        self.cmap_limits = {
            "res_xx": (-1, 2),
            "res_xy": (0, 3),
            "res_yx": (0, 3),
            "res_yy": (-1, 2),
            "res_det": (0, 3),
            "phase_xx": (-180, 180),
            "phase_xy": (0, 100),
            "phase_yx": (0, 100),
            "phase_yy": (-180, 180),
            "phase_det": (0, 100),
        }

        self.label_dict = {
            "res_xx": "$\\rho_{xx}  \\mathrm{[\Omega m]}$",
            "res_xy": "$\\rho_{xy}  \\mathrm{[\Omega m]}$",
            "res_yx": "$\\rho_{yx}  \\mathrm{[\Omega m]}$",
            "res_yy": "$\\rho_{yy}  \\mathrm{[\Omega m]}$",
            "res_det": "$\\rho_{det}  \\mathrm{[\Omega m]}$",
            "phase_xx": "$\\phi_{xx}$",
            "phase_xy": "$\\phi_{xy}$",
            "phase_yx": "$\\phi_{yx}$",
            "phase_yy": "$\\phi_{yy}$",
            "phase_det": "$\\phi_{det}$",
        }

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    @property
    def map_units(self):
        return self._map_units

    @map_units.setter
    def map_units(self, value):
        self._map_units = value
        if value in ["km"]:
            self.scale = 1.0 / 1000
            self.cell_size = 0.2
        if value in ["m"]:
            self.scale = 1
            self.cell_size = 200
        else:
            self.scale = 1.0

    def _get_n_rows(self):
        """
        Get the number of rows in the subplot

        :return: DESCRIPTION
        :rtype: TYPE

        """
        n = 0
        if self.plot_resistivity:
            n += 1
        if self.plot_phase:
            n += 1
        return n

    def _get_n_columns(self):
        """get the number of columns in the subplot"""
        n = 0

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if getattr(self, f"plot_{cc}"):
                n += 1

        return n

    def _get_n_subplots(self):
        """
        Get the subplot indices
        """
        nr = self._get_n_rows()
        nc = self._get_n_columns()

        subplot_dict = {
            "res_xx": None,
            "res_xy": None,
            "res_yx": None,
            "res_yy": None,
            "res_det": None,
            "phase_xx": None,
            "phase_xy": None,
            "phase_yx": None,
            "phase_yy": None,
            "phase_det": None,
        }

        plot_num = 0
        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_resistivity:
                if getattr(self, f"plot_{cc}"):
                    plot_num += 1
                    subplot_dict[f"res_{cc}"] = (nr, nc, plot_num)

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_phase:
                if getattr(self, f"plot_{cc}"):
                    plot_num += 1
                    subplot_dict[f"phase_{cc}"] = (nr, nc, plot_num)

        return subplot_dict

    def _get_subplots(self):
        """
        get the subplots

        :return: DESCRIPTION
        :rtype: TYPE

        """
        subplot_dict = self._get_n_subplots()
        ax_dict = {}

        for cc in ["xx", "xy", "yx", "yy", "det"]:
            if self.plot_resistivity:
                comp = f"res_{cc}"
                if getattr(self, f"plot_{cc}"):
                    ax_dict[comp] = self.fig.add_subplot(
                        *subplot_dict[comp], aspect="equal"
                    )

            if self.plot_phase:
                comp = f"phase_{cc}"
                if getattr(self, f"plot_{cc}"):
                    print()
                    ax_dict[comp] = self.fig.add_subplot(
                        *subplot_dict[comp], aspect="equal"
                    )

        share = [ax for comp, ax in ax_dict.items() if ax is not None]

        # share x and y across all subplots for easier zooming
        for ax in share[1:]:
            ax.sharex(share[0])
            ax.sharey(share[0])

        return ax_dict

    def _get_interpolated_z(self, tf):
        return np.array(
            [
                [
                    tf.z_interp_dict["zxx"]["real"](1 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zxx"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                    tf.z_interp_dict["zxy"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zxy"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                ],
                [
                    tf.z_interp_dict["zyx"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zyx"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                    tf.z_interp_dict["zyy"]["real"](1.0 / self.plot_period)[0]
                    + 1j
                    * tf.z_interp_dict["zyy"]["imag"](1.0 / self.plot_period)[
                        0
                    ],
                ],
            ]
        )

    def _get_data_array(self):
        """
        make a data array to plot
        """

        plot_array = np.zeros(
            len(self.tf_list),
            dtype=[
                ("station", "U20"),
                ("latitude", float),
                ("longitude", float),
                ("elevation", float),
                ("res_xx", float),
                ("res_xy", float),
                ("res_yx", float),
                ("res_yy", float),
                ("res_det", float),
                ("phase_xx", float),
                ("phase_xy", float),
                ("phase_yx", float),
                ("phase_yy", float),
                ("phase_det", float),
            ],
        )

        for ii, tf in enumerate(self.tf_list):
            z = self._get_interpolated_z(tf)
            z_object = Z(z, freq=[1.0 / self.plot_period])

            plot_array["station"][ii] = tf.station
            plot_array["latitude"][ii] = tf.latitude
            plot_array["longitude"][ii] = tf.longitude
            if tf.elevation is not None:
                plot_array["elevation"][ii] = tf.elevation * self.scale

            plot_array["res_xx"][ii] = z_object.res_xx[0]
            plot_array["res_xy"][ii] = z_object.res_xy[0]
            plot_array["res_yx"][ii] = z_object.res_yx[0]
            plot_array["res_yy"][ii] = z_object.res_yy[0]
            plot_array["res_det"][ii] = z_object.res_det[0]

            plot_array["phase_xx"][ii] = z_object.phase_xx[0]
            plot_array["phase_xy"][ii] = z_object.phase_xy[0]
            plot_array["phase_yx"][ii] = z_object.phase_yx[0] + 180
            plot_array["phase_yy"][ii] = z_object.phase_yy[0]
            plot_array["phase_det"][ii] = z_object.phase_det[0]

        return plot_array

    def _interpolate_to_map(self, plot_x, plot_y, plot_array, component):
        """
        Interpolate points onto a map


        :param plot_array: DESCRIPTION
        :type plot_array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        points = np.array([plot_array["longitude"], plot_array["latitude"]]).T
        values = plot_array[component]

        grid_x, grid_y = np.meshgrid(plot_x, plot_y)
        image = interpolate.griddata(
            points,
            values,
            (grid_x, grid_y),
            method=self.interpolation_method,
        )

        if "res" in component:
            image = np.log10(image)

        return grid_x, grid_y, image

    def _fancy_interpolate_to_map(
        self,
        plot_x,
        plot_y,
        plot_array,
        component,
        nearest_neighbors=7,
        interp_pow=4,
    ):
        """
        Using some fancy triangulation and delauny interpolation

        :param plot_x: DESCRIPTION
        :type plot_x: TYPE
        :param plot_y: DESCRIPTION
        :type plot_y: TYPE
        :param plot_array: DESCRIPTION
        :type plot_array: TYPE
        :param component: DESCRIPTION
        :type component: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        x = plot_array["longitude"]
        y = plot_array["latitude"]

        # add padding to the locations
        ds = self.cell_size * self.n_padding_cells

        ex = plot_array["longitude"]
        ey = plot_array["latitude"]

        ex[np.argmin(ex)] -= ds
        ex[np.argmax(ex)] += ds
        ey[np.argmin(ey)] -= ds
        ey[np.argmax(ey)] += ds

        rx, ry = np.meshgrid(plot_x, plot_y)
        rx = rx.flatten()
        ry = ry.flatten()

        triangulation = tri.Triangulation(rx, ry)

        mx = rx[triangulation.triangles].mean(axis=1)
        my = ry[triangulation.triangles].mean(axis=1)

        mxmy = np.array([mx, my]).T
        exey = np.array([ex, ey]).T

        inside_indices = in_hull(mxmy, exey)
        inside_indices = np.bool_(inside_indices)
        triangulation.set_mask(~inside_indices)

        tree = cKDTree(np.array([x, y]).T)

        xy = np.array([rx, ry]).T
        d, l = tree.query(xy, k=nearest_neighbors)

        image = None
        values = plot_array[component]

        if nearest_neighbors == 1:
            # extract nearest neighbour values
            image = values[l]
        else:
            image = np.zeros((xy.shape[0]))

            # field values are directly assigned for coincident locations
            coincident_indices = d[:, 0] == 0
            image[coincident_indices] = values[l[coincident_indices, 0]]

            # perform idw interpolation for non-coincident locations
            idw_indices = d[:, 0] != 0
            w = np.zeros(d.shape)
            w[idw_indices, :] = 1.0 / np.power(d[idw_indices, :], interp_pow)

            image[idw_indices] = np.sum(
                w[idw_indices, :] * values[l[idw_indices, :]], axis=1
            ) / np.sum(w[idw_indices, :], axis=1)

        if "res" in component:
            image = np.log10(image)

        return triangulation, image, inside_indices

    def _get_cmap(self, component):
        """
        get color map with proper limits
        """
        if "res" in component:
            cmap = self.res_cmap
        elif "phase" in component:
            cmap = self.phase_cmap

        # if isinstance(cmap, str):
        #     cmap = plt.get_cmap(cmap)

        # # set cmap values for over and under
        # norm = colors.Normalize(
        #     vmin=self.cmap_limits[component][0],
        #     vmax=self.cmap_limits[component][1],
        # )
        # cmap.set_over(cmap(norm(self.cmap_limits[component][1])))
        # cmap.set_under(cmap(norm(self.cmap_limits[component][0])))

        return cmap

    def _get_colorbar(self, ax, im_mappable, component):
        """

        :param component: DESCRIPTION
        :type component: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if "res" in component:
            cb = plt.colorbar(
                im_mappable,
                ax=ax,
                ticks=ticker.FixedLocator(
                    np.arange(
                        int(np.round(self.cmap_limits[component][0])),
                        int(np.round(self.cmap_limits[component][1])) + 1,
                    )
                ),
                shrink=0.6,
            )
            labels = [
                period_label_dict[dd]
                for dd in np.arange(
                    int(np.round(self.cmap_limits[component][0])),
                    int(np.round(self.cmap_limits[component][1])) + 1,
                )
            ]
            cb.ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))

        elif "phase" in component:
            cb = plt.colorbar(
                im_mappable,
                ax=ax,
                shrink=0.6,
            )

        cb.ax.tick_params(
            axis="both", which="major", labelsize=self.font_size - 1
        )
        cb.ax.tick_params(
            axis="both", which="minor", labelsize=self.font_size - 1
        )

        return cb

    def _get_plot_xy(self, plot_array):
        """
        Get the plot_x and plot_y along a regular grid

        :return: DESCRIPTION
        :rtype: TYPE

        """
        # create uniform x, y to plot on.
        ds = self.cell_size * self.n_padding_cells
        n_points = int(
            abs(
                plot_array["longitude"].max()
                - plot_array["longitude"].min()
                + 2 * ds
            )
            / self.cell_size
        )
        plot_x = np.linspace(
            plot_array["longitude"].min() - ds,
            plot_array["longitude"].max() + ds,
            n_points,
        )

        n_points = int(
            abs(
                plot_array["latitude"].max()
                - plot_array["latitude"].min()
                + 2 * ds
            )
            / self.cell_size
        )
        plot_y = np.linspace(
            plot_array["latitude"].min() - ds,
            plot_array["latitude"].max() + ds,
            n_points,
        )

        return plot_x, plot_y

    # -----------------------------------------------
    # The main plot method for this module
    # -------------------------------------------------
    def plot(self):
        """ """

        # set position properties for the plot
        self._set_subplot_params()

        # make figure instance
        self.fig = plt.figure(
            self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi
        )

        # clear the figure if there is already one up
        plt.clf()

        subplot_dict = self._get_subplots()

        plot_array = self._get_data_array()
        plot_x, plot_y = self._get_plot_xy(plot_array)

        # plot results
        subplot_numbers = self._get_n_subplots()
        for comp, ax in subplot_dict.items():
            cmap = self._get_cmap(comp)

            if self.interpolation_method in ["nearest", "linear", "cubic"]:
                x, y, image = self._interpolate_to_map(
                    plot_x, plot_y, plot_array, comp
                )

                im = ax.pcolormesh(
                    x,
                    y,
                    image,
                    cmap=cmap,
                    vmin=self.cmap_limits[comp][0],
                    vmax=self.cmap_limits[comp][1],
                )
            elif self.interpolation_method in [
                "fancy",
                "delaunay",
                "triangulate",
            ]:
                triangulation, image, indices = self._fancy_interpolate_to_map(
                    plot_x, plot_y, plot_array, comp
                )
                im = ax.tricontourf(
                    triangulation,
                    image,
                    mask=indices,
                    levels=np.linspace(
                        self.cmap_limits[comp][0],
                        self.cmap_limits[comp][1],
                        50,
                    ),
                    extend="both",
                    cmap=cmap,
                )

            cb = self._get_colorbar(ax, im, comp)

            # show stations
            if self.plot_stations:
                if self.plot_stations:
                    ax.scatter(
                        plot_array["longitude"],
                        plot_array["latitude"],
                        marker=self.marker,
                        s=self.marker_size,
                        c=self.marker_color,
                    )
            # Label plots
            ax.text(
                0.01,
                0.9,
                self.label_dict[comp],
                fontdict={"size": self.font_size + 2},
                transform=ax.transAxes,
            )

            if (
                subplot_numbers[comp][2] == 1
                or subplot_numbers[comp][2] == subplot_numbers[comp][1] + 1
            ):
                ax.set_ylabel("Latitude (deg)", fontdict=self.font_dict)
            if subplot_numbers[comp][0] == 1:
                ax.set_xlabel("Longitude (deg)", fontdict=self.font_dict)
            elif subplot_numbers[comp][0] == 2:
                if subplot_numbers[comp][2] > (subplot_numbers[comp][1]):
                    ax.set_xlabel("Longitude (deg)", fontdict=self.font_dict)

        # Plot title
        self.fig.suptitle(f"Plot Period: {self.plot_period:.5g} s", y=0.985)
