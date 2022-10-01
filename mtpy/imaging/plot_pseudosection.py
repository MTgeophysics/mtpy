# -*- coding: utf-8 -*-
"""
Created on Thu May 30 18:39:58 2013

@author: jpeacock-pr
"""

# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as mcb
import matplotlib.colors as colors
from matplotlib import ticker

from mtpy.imaging.mtplot_tools import (
    PlotBaseProfile,
    griddata_interpolate,
    triangulate_interpolation,
)
import mtpy.imaging.mtcolors as mtcl

# ==============================================================================


class PlotResPhasePseudoSection(PlotBaseProfile):
    """
    plot a resistivity and phase pseudo section for different components

    Need to input one of the following lists:


    """

    def __init__(self, tf_list, **kwargs):
        """
        Initialize parameters
        """

        super().__init__(**kwargs)
        # read in key word arguments and set defaults if none given
        self.tf_list = tf_list

        # --> set figure parameters
        self.aspect = kwargs.pop("aspect", "auto")

        self.xtickspace = kwargs.pop("xtickspace", 1)
        self.stationid = kwargs.pop("stationid", [0, 4])
        self.linedir = kwargs.pop("linedir", "ew")

        # --> set plots to plot and how to plot them
        self.plot_xx = False
        self.plot_xy = True
        self.plot_yx = True
        self.plot_yy = False
        self.plot_det = False
        self.plot_resistivity = True
        self.plot_phase = True
        self.data_df = None
        self.n_periods = 60
        self.interpolation_method = "nearest"

        # --> set plot limits
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

        # --> set colormaps Note only mtcolors is supported
        self.res_cmap = mtcl.cmapdict["mt_rd2gr2bl"]
        self.phase_cmap = mtcl.cmapdict["mt_bl2gr2rd"]

        for key, value in kwargs.items():
            setattr(self, key, value)

        if self.show_plot:
            self.plot()

    def _get_period_array(self, df):
        """
        Get the period array to interpolate on to
        """

        p_min = np.log10(df.period.min())
        p_max = np.log10(df.period.max())

        return np.logspace(p_min, p_max, self.n_periods)

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
                    ax_dict[comp] = self.fig.add_subplot(
                        *subplot_dict[comp], aspect="equal"
                    )

        share = [ax for comp, ax in ax_dict.items() if ax is not None]

        # share x and y across all subplots for easier zooming
        for ax in share[1:]:
            ax.sharex(share[0])
            ax.sharey(share[0])

        return ax_dict

    def _get_data_df(self):
        """
        get resistivity and phase values in the correct order according to
        offsets and periods.

        """

        self._get_profile_line()

        entries = []

        for ii, tf in enumerate(self.tf_list):
            offset = self._get_offset(tf)
            rp = tf.Z

            for ii, period in enumerate(tf.period):
                entry = {
                    "x": offset,
                    "y": period,
                    "res_xx": rp.res_xx[ii],
                    "res_xy": rp.res_xy[ii],
                    "res_yx": rp.res_yx[ii],
                    "res_yy": rp.res_yy[ii],
                    "res_det": rp.res_det[ii],
                    "phase_xx": rp.phase_xx[ii],
                    "phase_xy": rp.phase_xy[ii],
                    "phase_yx": rp.phase_yx[ii] + 180,
                    "phase_yy": rp.phase_yy[ii],
                    "phase_det": rp.phase_det[ii],
                }
                entries.append(entry)

        return pd.DataFrame(entries)

    def _get_cmap(self, component):
        """
        get color map with proper limits
        """
        if "res" in component:
            cmap = self.res_cmap
        elif "phase" in component:
            cmap = self.phase_cmap

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
                extend="both",
            )
            labels = [
                self.period_label_dict[dd]
                for dd in np.arange(
                    int(np.round(self.cmap_limits[component][0])),
                    int(np.round(self.cmap_limits[component][1])) + 1,
                )
            ]
            cb.ax.yaxis.set_major_formatter(ticker.FixedFormatter(labels))

        elif "phase" in component:
            cb = plt.colorbar(im_mappable, ax=ax, shrink=0.6, extend="both")

        cb.ax.tick_params(
            axis="both", which="major", labelsize=self.font_size - 1
        )
        cb.ax.tick_params(
            axis="both", which="minor", labelsize=self.font_size - 1
        )

        return cb

    def plot(self):
        """

        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.data_df == None:
            self.data_df = self._get_data_df()

        plot_periods = self._get_period_array(self.data_df)

        # set position properties for the plot
        self._set_subplot_params()

        # make figure instance
        self.fig = plt.figure(
            self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi
        )

        # clear the figure if there is already one up
        plt.clf()

        # Get dictionary of subplots
        subplot_dict = self._get_subplots()

        # plot results
        subplot_numbers = self._get_n_subplots()
        for comp, ax in subplot_dict.items():
            cmap = self._get_cmap(comp)

            comp_df = self.data_df.iloc[np.nonzero(self.data_df[comp])]
            if self.interpolation_method in ["nearest", "linear", "cubic"]:

                x, y, image = griddata_interpolate(
                    comp_df.offset,
                    comp_df.period,
                    comp_df[comp],
                    self.data_df.offset,
                    plot_periods,
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
                triangulation, image, indices = triangulate_interpolation(
                    comp_df.offset,
                    comp_df.period,
                    comp_df[comp],
                    comp_df.offset,
                    comp_df.period,
                    self.data_df.offset,
                    plot_periods,
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

            self._get_colorbar(ax, im, comp)

            # Label plots
            ax.title(
                self.label_dict[comp],
                fontdict={"size": self.font_size + 2},
            )

            if (
                subplot_numbers[comp][2] == 1
                or subplot_numbers[comp][2] == subplot_numbers[comp][1] + 1
            ):
                ax.set_ylabel("Period (s)", fontdict=self.font_dict)
            if subplot_numbers[comp][0] == 1:
                ax.set_xlabel("Station", fontdict=self.font_dict)
            elif subplot_numbers[comp][0] == 2:
                if subplot_numbers[comp][2] > (subplot_numbers[comp][1]):
                    ax.set_xlabel("Station", fontdict=self.font_dict)

        self.fig.tight_layout()
        plt.show()
