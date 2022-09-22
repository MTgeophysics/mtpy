# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:58:58 2022

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate


from mtpy.imaging.mtplot_tools import PlotBase
from mtpy.analysis.niblettbostick import calculate_depth_of_investigation
from mtpy.core.z import Z


# =============================================================================


class PlotPenetrationDepthMap(PlotBase):
    """
    Plot the depth of penetration based on the Niblett-Bostick approximation.


    """

    def __init__(self, tf_list, **kwargs):

        self.tf_list = tf_list

        super().__init__(**kwargs)

        self.depth_units = "km"
        self.plot_period = 1
        self.plot_det = True
        self.plot_te = True
        self.plot_tm = True
        self.grid_element_size = 0.002
        self.grid_pad = 2
        self.interpolation_method = "cubic"
        self.plot_stations = True
        self.depth_cmap = "magma"
        self.marker_color = "k"
        self.marker_size = 10
        self.subplot_title_dict = {
            "det": "Determinant",
            "xy": "TE Mode",
            "yx": "TM Mode",
        }
        self.depth_range = [None, None]

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    @property
    def depth_units(self):
        return self._depth_units

    @depth_units.setter
    def depth_units(self, value):
        self._depth_units = value
        if value in ["km"]:
            self.depth_scale = 1.0 / 1000
        if value in ["m"]:
            self.depth_scale = 1

    def _get_nb_estimation(self, z_object):
        """
        get the depth of investigation estimation

        """

        return calculate_depth_of_investigation(z_object)

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

    def _get_depth_array(self):
        """
        get a depth array with xyz values

        :return: DESCRIPTION
        :rtype: TYPE

        """

        depth_array = np.zeros(
            len(self.tf_list),
            dtype=[
                ("station", "U20"),
                ("latitude", float),
                ("longitude", float),
                ("elevation", float),
                ("det", float),
                ("xy", float),
                ("yx", float),
            ],
        )

        for ii, tf in enumerate(self.tf_list):
            z = self._get_interpolated_z(tf)
            z_object = Z(z, freq=[1.0 / self.plot_period])
            d = self._get_nb_estimation(z_object)

            depth_array["station"][ii] = tf.station
            depth_array["latitude"][ii] = tf.latitude
            depth_array["longitude"][ii] = tf.longitude
            elev = 0
            if tf.elevation is not None:
                depth_array["elevation"][ii] = tf.elevation * self.depth_scale
                elev = tf.elevation
            depth_array["det"][ii] = (
                d["depth_det"][0] - elev
            ) * self.depth_scale
            depth_array["xy"][ii] = (
                d["depth_xy"][0] - elev
            ) * self.depth_scale
            depth_array["yx"][ii] = (
                d["depth_yx"][0] - elev
            ) * self.depth_scale

        return depth_array

    def _interpolate_to_map(self, depth_array, component):
        """
        Interpolate points onto a map


        :param depth_array: DESCRIPTION
        :type depth_array: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        points = np.array(
            [depth_array["longitude"], depth_array["latitude"]]
        ).T
        values = depth_array[component]

        x = np.linspace(
            depth_array["longitude"].min()
            - (self.grid_element_size * self.grid_pad),
            depth_array["longitude"].max()
            + (self.grid_element_size * self.grid_pad),
        )

        y = np.linspace(
            depth_array["latitude"].min()
            - (self.grid_element_size * self.grid_pad),
            depth_array["latitude"].max()
            + (self.grid_element_size * self.grid_pad),
        )

        grid_x, grid_y = np.meshgrid(x, y)

        return (
            grid_x,
            grid_y,
            interpolate.griddata(
                points,
                values,
                (grid_x, grid_y),
                method=self.interpolation_method,
            ),
        )

    def _get_n_subplots(self):
        """
        get the number of subplots

        :return: DESCRIPTION
        :rtype: TYPE

        """
        n = 0
        if self.plot_det:
            n += 1
        if self.plot_te:
            n += 1
        if self.plot_tm:
            n += 1

        return n

    def _get_subplots(self):
        """
        get subplots
        """
        n = self._get_n_subplots()
        ax_det = None
        ax_xy = None
        ax_yx = None
        if self.plot_det:
            ax_det = self.fig.add_subplot(1, n, 1, aspect="equal")
            if self.plot_te:
                ax_xy = self.fig.add_subplot(1, n, 2, aspect="equal")
                if self.plot_tm:
                    ax_yx = self.fig.add_subplot(1, n, 3, aspect="equal")
            else:
                if self.plot_tm:
                    ax_yx = self.fig.add_subplot(1, n, 2, aspect="equal")
        else:
            if self.plot_te:
                ax_xy = self.fig.add_subplot(1, n, 1, aspect="equal")
                if self.plot_tm:
                    ax_yx = self.fig.add_subplot(1, n, 2, aspect="equal")
            else:
                if self.plot_tm:
                    ax_yx = self.fig.add_subplot(1, n, 1, aspect="equal")

        return ax_det, ax_xy, ax_yx

    def _get_plot_component_dict(self):
        """
        Get all the components to plot
        """
        ax_det, ax_xy, ax_yx = self._get_subplots()

        components = {}
        if self.plot_det:
            components["det"] = ax_det
        if self.plot_te:
            components["xy"] = ax_xy
        if self.plot_tm:
            components["yx"] = ax_yx

        if len(components.keys()) == 0:
            raise ValueError(
                "must set at least one of the following to True  'plot_det', "
                "'plot_te', 'plot_tm"
            )

        return components

    def plot(self):
        """
        plot the depth of investigation as a 1d plot with period on the y-axis
        and depth on the x axis

        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        plot_components = self._get_plot_component_dict()

        depth_array = self._get_depth_array()

        for comp, ax in plot_components.items():
            plot_x, plot_y, plot_z = self._interpolate_to_map(
                depth_array, comp
            )

            im = ax.pcolormesh(
                plot_x,
                plot_y,
                plot_z,
                cmap=self.depth_cmap,
                vmin=self.depth_range[0],
                vmax=self.depth_range[1],
            )

            if self.plot_stations:
                ax.scatter(
                    depth_array["longitude"],
                    depth_array["latitude"],
                    marker=self.marker,
                    s=self.marker_size,
                    c=self.marker_color,
                )

            ax.set_xlabel("Longitude (deg)", fontdict=self.font_dict)
            ax.set_ylabel("Latitude (deg)", fontdict=self.font_dict)
            ax.set_title(
                self.subplot_title_dict[comp], fontdict=self.font_dict
            )

            plt.colorbar(
                im,
                ax=ax,
                label=f"Penetration Depth ({self.depth_units})",
                shrink=0.35,
            )

        self.fig.tight_layout()

        self.fig.suptitle(
            f"Depth of investigation for period {self.plot_period:5g} (s)",
            fontproperties=self.font_dict,
        )

        plt.show()
