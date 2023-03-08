# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 10:57:29 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colors as colors
from matplotlib import colorbar as mcb
from matplotlib import cm
from matplotlib import gridspec

try:
    import contextily as cx

    has_cx = True
except ModuleNotFoundError:
    has_cx = False

from mtpy.core.mt_dataframe import MTDataFrame
from mtpy.imaging.mtplot_tools import PlotBaseMaps

# =============================================================================


class PlotRMS(PlotBaseMaps):
    def __init__(self, dataframe, **kwargs):
        super().__init__(**kwargs)

        self.dataframe = dataframe
        self.dx = 0.0075
        self.rms_min = 0
        self.rms_max = 5
        self.rms_step = 0.5
        self.plot_station = True
        self.station_id = None
        self.stack_bottom = False

        self.comp_list = [
            "rms_zxx",
            "rms_zxy",
            "rms_zyx",
            "rms_zyy",
            "rms_tzx",
            "rms_tzy",
        ]
        self.distance_multiplier = [
            (-0.5, 1),
            (0.5, 1),
            (-0.5, 0),
            (0.5, 0),
            (-0.5, -1),
            (0.5, -1),
        ]

        self.color_dict = {
            "rms_z": (0, 162 / 255, 255 / 255),
            "rms_t": (255 / 255, 162 / 255, 0),
            "rms_zxx": (136 / 255, 235 / 255, 193 / 255),
            "rms_zxy": (84 / 255, 189 / 255, 215 / 255),
            "rms_zyx": (136 / 255, 84 / 255, 215 / 255),
            "rms_zyy": (206 / 255, 84 / 255, 215 / 255),
            "rms_tzx": (215 / 255, 210 / 255, 84 / 255),
            "rms_tzy": (215 / 255, 154 / 255, 84 / 255),
        }

        self.label_dict = {
            "rms_z": "Z",
            "rms_t": "Tipper",
            "rms_zxx": "Z_{xx}$",
            "rms_zxy": "Z_{xy}$",
            "rms_zyx": "Z_{yx}$",
            "rms_zyy": "Z_{yy}$",
            "rms_tzx": "T_{zx}$",
            "rms_tzy": "T_{zy}$",
        }

        self.rms_cmap = "jet"

        self.subplot_left = 0.05
        self.subplot_right = 0.99
        self.subplot_bottom = 0.09
        self.subplot_top = 0.99

        self.box_size = 30

        self.cx_source = None
        self.cx_zoom = None
        if has_cx:
            self.cx_source = cx.providers.USGS.USTopo

        for key, value in kwargs.items():
            setattr(self, key, value)

    @property
    def dataframe(self):
        return self._mt_dataframe.dataframe

    @dataframe.setter
    def dataframe(self, df):
        """
        Set dataframe to an MTDataframe
        :param df: DESCRIPTION
        :type df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if df is None:
            self._mt_dataframe = MTDataFrame()

        elif isinstance(df, (pd.DataFrame, MTDataFrame, np.ndarray)):
            self._mt_dataframe = MTDataFrame(df)

        else:
            raise TypeError(
                f"Input must be a dataframe or MTDataFrame object not {type(df)}"
            )

    @property
    def rms_cmap(self):
        return self._rms_cmap

    @rms_cmap.setter
    def rms_cmap(self, value):
        if isinstance(value, str):
            self._rms_cmap = cm.get_cmap(value)

        elif isinstance(value, colors.LinearSegmentedColormap):
            self._rms_cmap = value

        else:
            self._rms_cmap = cm.get_cmap("jet")

    def _plot_rms_map(self):
        """
        plot rms map

        :return: DESCRIPTION
        :rtype: TYPE

        """

        cb_norm = colors.BoundaryNorm(
            np.arange(
                self.rms_min, self.rms_max + self.rms_step, self.rms_step
            ),
            self.rms_cmap.N,
        )

        for dm, comp in zip(self.distance_multiplier, self.comp_list):
            for station in self.dataframe.station.unique():

                sdf = self._mt_dataframe.get_station_df(station)
                rms = sdf[comp].mean()
                self.ax1.scatter(
                    sdf.longitude.iloc[0] + (self.dx / 2) * dm[0],
                    sdf.latitude.iloc[0] + (self.dx / 2) * dm[1],
                    c=rms,
                    marker="s",
                    s=self.box_size,
                    edgecolors=(0, 0, 0),
                    cmap=self.rms_cmap,
                    norm=cb_norm,
                )
                if self.plot_station:
                    self.ax1.text(
                        sdf.longitude.iloc[0],
                        sdf.latitude.iloc[0] + self.dx,
                        station,
                        ha="center",
                        va="baseline",
                    )

        if has_cx:
            if has_cx:
                try:
                    cx_kwargs = {"source": self.cx_source, "crs": "EPSG:4326"}
                    if self.cx_zoom is not None:
                        cx_kwargs["zoom"] = self.cx_zoom
                    cx.add_basemap(
                        self.ax1,
                        **cx_kwargs,
                    )
                except Exception as error:
                    self.logger.warning(
                        f"Could not add base map because {error}"
                    )

        cb_ax, _ = mcb.make_axes(self.ax1, shrink=0.5)
        cb = mcb.ColorbarBase(cb_ax, cmap=self.rms_cmap, norm=cb_norm)

    @property
    def rms_per_period_all(self):
        """
        RMS per period
        """

        if self.dataframe is not None:
            rms_list = []
            for period in self.dataframe.period.unique():
                comp_df = self.dataframe.loc[
                    self.dataframe.period == period,
                    [
                        "rms_zxx",
                        "rms_zxy",
                        "rms_zyx",
                        "rms_zyy",
                        "rms_tzx",
                        "rms_tzy",
                    ],
                ]

                mean_dict = {"period": period}
                for comp in comp_df.columns:
                    mean_dict[comp] = comp_df.loc[:, comp].mean()

                rms_list.append(mean_dict)

            df = pd.DataFrame(rms_list)
            df = df.set_index("period")

            return df

    @property
    def rms_per_station(self):
        """
        RMS per period
        """

        if self.dataframe is not None:
            rms_list = []
            for station in self.dataframe.station.unique():
                z_df = self.dataframe.loc[
                    self.dataframe.station == station,
                    ["rms_zxx", "rms_zxy", "rms_zyx", "rms_zyy"],
                ]
                t_df = self.dataframe.loc[
                    self.dataframe.station == station, ["rms_tzx", "rms_tzy"]
                ]

                rms_list.append(
                    {
                        "station": station,
                        "rms_z": z_df.mean().mean(),
                        "rms_t": t_df.mean().mean(),
                    }
                )

            df = pd.DataFrame(rms_list)
            df = df.set_index("station")

            return df

    def _plot_by_period(self):
        """
        plot by period

        :return: DESCRIPTION
        :rtype: TYPE

        """

        df = self.rms_per_period_all.copy()
        plot_list = []
        color_list = []
        for comp in df.columns:
            if not np.all(np.isnan(df[comp])):
                plot_list.append(comp)
                color_list.append(self.color_dict[comp])

        ax = df.plot.bar(
            y=plot_list,
            color=color_list,
            xlabel="Period (s)",
            ylabel="normalized RMS",
            grid=True,
            ax=self.ax2,
        )
        ax.set_axisbelow(True)

        ax.set_xticklabels(
            [f"{float(x.get_text()):.4g}" for x in ax.get_xticklabels()]
        )

        return ax

    def _plot_by_station(self):
        """
        plot by station

        :return: DESCRIPTION
        :rtype: TYPE

        """

        df = self.rms_per_station.copy()
        plot_list = []
        color_list = []
        for comp in df.columns:
            if not np.all(np.isnan(df[comp])):
                plot_list.append(comp)
                color_list.append(self.color_dict[comp])

        ax = df.plot.bar(
            y=plot_list,
            color=color_list,
            xlabel="Station",
            ylabel="normalized RMS",
            grid=True,
            ax=self.ax3,
        )
        ax.set_axisbelow(True)

        return ax

    def _get_subplots(self, fig):

        if self.stack_bottom:
            gs1 = gridspec.GridSpec(2, 2, hspace=0.25, wspace=0.075)

            self.ax1 = fig.add_subplot(gs1[0, :], aspect="equal")
            self.ax2 = fig.add_subplot(gs1[1, 0])
            self.ax3 = fig.add_subplot(gs1[1, 1], sharey=self.ax2)
        else:
            gs1 = gridspec.GridSpec(2, 2, hspace=0.35, wspace=0.075)

            self.ax1 = fig.add_subplot(gs1[:, 0], aspect="equal")
            self.ax2 = fig.add_subplot(gs1[0, 1])
            self.ax3 = fig.add_subplot(gs1[1, 1], sharey=self.ax2)

    def plot(self, **kwargs):
        """

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self._set_subplot_params()

        self.fig = plt.figure(
            self.fig_num, figsize=self.fig_size, dpi=self.fig_dpi
        )

        plt.clf()

        self._get_subplots(self.fig)

        self._plot_rms_map()
        self._plot_by_period()
        self._plot_by_station()
