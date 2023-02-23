# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:58:58 2022

@author: jpeacock
"""
# =============================================================================
# Imports
# =============================================================================
import matplotlib.pyplot as plt

from mtpy.imaging.mtplot_tools import PlotBase

# =============================================================================


class PlotPenetrationDepth1D(PlotBase):
    """
    Plot the depth of penetration based on the Niblett-Bostick approximation.


    """

    def __init__(self, tf, **kwargs):

        self.tf = tf

        super().__init__(**kwargs)

        self.depth_units = "km"

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

    def _get_nb_estimation(self):
        """
        get the depth of investigation estimation

        """

        return self.tf.Z.estimate_depth_of_investigation()

    def plot(self):
        """
        plot the depth of investigation as a 1d plot with period on the y-axis
        and depth on the x axis

        :return: DESCRIPTION
        :rtype: TYPE

        """

        depth_array = self._get_nb_estimation()

        self.fig = plt.figure(self.fig_num, self.fig_size, dpi=self.fig_dpi)
        plt.clf()

        self.ax = self.fig.add_subplot(1, 1, 1, aspect="equal")

        self.ax.fill_betweenx(
            depth_array["period"],
            depth_array["depth_min"] * self.depth_scale,
            depth_array["depth_max"] * self.depth_scale,
            alpha=0.5,
            color=(0.5, 0.5, 0.5),
        )

        line_list = []
        label_list = ["TE", "TM", "DET"]
        for comp in ["xy", "yx", "det"]:
            (line,) = self.ax.loglog(
                depth_array[f"depth_{comp}"] * self.depth_scale,
                depth_array["period"],
                marker=getattr(self, f"{comp}_marker"),
                color=getattr(self, f"{comp}_color"),
                ls=getattr(self, f"{comp}_ls"),
                ms=self.marker_size,
            )

            line_list.append(line)

        self.ax.set_xlabel(
            f"Depth ({self.depth_units})", fontdict=self.font_dict
        )
        self.ax.set_ylabel("Period (s)", fontdict=self.font_dict)

        self.ax.set_ylim(self.set_period_limits(depth_array["period"])[::-1])
        self.fig.suptitle(
            f"Depth of investigation for {self.tf.station}",
            fontproperties=self.font_dict,
        )

        self.ax.grid(which="major", lw=0.75, ls="--", color=(0.25, 0.25, 0.25))
        self.ax.grid(which="minor", lw=0.5, ls=":", color=(0.75, 0.75, 0.75))

        self.ax.set_axisbelow(True)

        self.ax.legend(
            line_list,
            label_list,
            prop={"size": self.font_size},
            loc="upper right",
        )

        plt.show()

        self.fig.tight_layout()
