"""
==================
ModEM
==================

residuals class to contain RMS information

revised by JP 2017
revised by AK 2017 to bring across functionality from ak branch

"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

import numpy as np
import pandas as pd

from .data import Data
from mtpy.modeling.plots import PlotRMS

# =============================================================================


class Residual(Data):
    """
    class to contain residuals for each data point, and rms values for each
    station

    ====================== ====================================================
    Attributes/Key Words   Description
    ====================== ====================================================
    work_dir
    residual_fn            full path to data file
    residual_array         numpy.ndarray (num_stations) structured to store
                           data.  keys are:
                               * station --> station name
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * zone --> UTM zone
                               * z --> impedance tensor residual (measured - modelled)
                                       (num_freq, 2, 2)
                               * z_err --> impedance tensor error array with
                                       shape (num_freq, 2, 2)
                               * tip --> Tipper residual (measured - modelled)
                                       (num_freq, 1, 2)
                               * tipperr --> Tipper array with shape
                                       (num_freq, 1, 2)
    rms
    rms_array              numpy.ndarray structured to store station
                           location values and rms.  Keys are:
                               * station --> station name
                               * east --> UTM east (m)
                               * north --> UTM north (m)
                               * lat --> latitude in decimal degrees
                               * lon --> longitude in decimal degrees
                               * elev --> elevation (m)
                               * zone --> UTM zone
                               * rel_east -- > relative east location to
                                               center_position (m)
                               * rel_north --> relative north location to
                                               center_position (m)
                               * rms --> root-mean-square residual for each
                                         station
    rms_tip
    rms_z
    ====================== ====================================================
    """

    # todo complete the doc above
    def __init__(self, **kwargs):
        self.work_dir = Path()
        self.residual_fn = None
        self.residual_array = None
        self.rms = None
        self.rms_array = None
        self.rms_tip = None
        self.rms_z = None

        super().__init__(**kwargs)

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

        for key, value in kwargs.items():
            setattr(self, key, value)

    def read_residual_file(self, residual_fn):
        """

        :param residual_fn: DESCRIPTION, defaults to None
        :type residual_fn: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.dataframe = self.read_data_file(residual_fn)
        self.calculate_rms()

    def calculate_rms(self):
        """
        add columns for rms
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.dataframe is None:
            return

        for col in ["zxx", "zxy", "zyx", "zyy", "tzx", "tzy"]:
            with np.errstate(divide="ignore", invalid="ignore"):
                self.dataframe[f"rms_{col}"] = np.abs(self.dataframe[col]) / (
                    np.real(self.dataframe[f"{col}_model_error"]) * np.sqrt(2)
                )

    @property
    def rms_per_period_all(self):
        """
        RMS per period
        """

        if self.dataframe is not None:
            rms_list = []
            for period in self.dataframe.period.unique():
                z_df = self.dataframe.loc[
                    self.dataframe.period == period,
                    ["rms_zxx", "rms_zxy", "rms_zyx", "rms_zyy"],
                ]
                t_df = self.dataframe.loc[
                    self.dataframe.period == period, ["rms_tzx", "rms_tzy"]
                ]

                rms_list.append(
                    {
                        "period": period,
                        "rms_z": z_df.mean().mean(),
                        "rms_t": t_df.mean().mean(),
                    }
                )

            df = pd.DataFrame(rms_list)
            df = df.set_index("period")

            return df

    @property
    def rms_per_period_per_component(self):
        """
        RMS per period by component

        :return: DESCRIPTION
        :rtype: TYPE

        """

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

    def plot_rms_per_period(self, plot_type="all", **kwargs):
        """

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if plot_type == "all":
            df = self.rms_per_period_all.copy()
        elif plot_type == "comp":
            df = self.rms_per_period_per_component.copy()

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
            **kwargs,
        )
        ax.set_axisbelow(True)

        return ax

    def plot_rms(self, **kwargs):
        """
        plot RMS in different views

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        plot_rms = PlotRMS(self.dataframe, **kwargs)
        plot_rms.plot()

        return plot_rms
