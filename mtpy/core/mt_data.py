# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:58:56 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================

from collections import OrderedDict

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from .mt import MT
from .mt_stations import MTStations


# =============================================================================


class MTData(OrderedDict, MTStations):
    def __init__(self, tf_list=None, **kwargs):

        if tf_list is not None:
            for tf in tf_list:
                self.add_station(tf)

        MTStations.__init__(self, tf_list, None, **kwargs)

    def _validate_item(self, tf):
        if not isinstance(tf, MT):
            raise TypeError(
                f"entry must be a mtpy.core.MT object not type({type(tf)})"
            )
        return tf

    def add_station(self, mt_object):
        """
        Add a new station's mt_object

        :param mt_object: DESCRIPTION
        :type mt_object: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        mt_object = self._validate_item(mt_object)
        self.__setitem__(mt_object.station, mt_object)

    def remove_station(self, station_id):
        """
        remove a station from the dictionary

        :param station_id: DESCRIPTION
        :type station_id: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if station_id in self.keys():
            del self[station_id]

    @property
    def n_stations(self):
        if self.mt_list is not None:
            return len(self.mt_list)

    def to_dataframe(self, utm_crs=None, cols=None):
        """

        :param utm_crs: DESCRIPTION, defaults to None
        :type utm_crs: TYPE, optional
        :param cols: DESCRIPTION, defaults to None
        :type cols: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        df_list = [
            tf.to_dataframe(utm_crs=utm_crs, cols=cols) for tf in self.values()
        ]

        return pd.concat(df_list)

    def from_dataframe(self, df):
        """
        Create an dictionary of MT objects from a dataframe

        :param df: dataframe of mt data
        :type df: `pandas.DataFrame`
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for station in df.station.unique():
            sdf = df.loc[df.station == station]
            mt_object = MT()
            mt_object.from_dataframe(sdf)
            self.update(OrderedDict([(mt_object.station, mt_object)]))

    def interpolate(self, new_periods):
        """
        Interpolate onto common period range

        :param new_periods: DESCRIPTION
        :type new_periods: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for mt_obj in self.values():
            interp_periods = new_periods[
                np.where(
                    (new_periods >= mt_obj.period.max())
                    & (new_periods <= mt_obj.period.min())
                )
            ]

            interp_z, interp_t = mt_obj.interpolate(1.0 / interp_periods)

            mt_obj.Z = interp_z
            mt_obj.Tipper = interp_t

    def rotate(self, rotation_angle):
        """
        rotate the data by the given angle assuming positive clockwise with
        north = 0, east = 90.

        :param rotation_angle: DESCRIPTION
        :type rotation_angle: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for mt_obj in self.values():
            mt_obj.rotation_angle = rotation_angle

    def compute_model_errors(
        self,
        z_error_value=5,
        z_error_type="geometric_mean",
        z_floor=True,
        t_error_value=0.02,
        t_error_type="absolute",
        t_floor=True,
    ):

        """
        Compute mode errors based on the error type

        ========================== ===========================================
        key                        definition
        ========================== ===========================================
        egbert                     error_value * sqrt(Zxy * Zyx)
        geometric_mean             error_value * sqrt(Zxy * Zyx)
        arithmetic_mean            error_value * (Zxy + Zyx) / 2
        mean_od                    error_value * (Zxy + Zyx) / 2
        off_diagonals              zxx_err == zxy_err, zyx_err == zyy_err
        median                     error_value * median(z)
        eigen                      error_value * mean(eigen(z))
        percent                    error_value * z
        absolute                   error_value
        ========================== ===========================================

        :param z_error_value: DESCRIPTION, defaults to 5
        :type z_error_value: TYPE, optional
        :param z_error_type: DESCRIPTION, defaults to "geometric_mean"
        :type z_error_type: TYPE, optional
        :param z_floor: DESCRIPTION, defaults to True
        :type z_floor: TYPE, optional
        :param t_error_value: DESCRIPTION, defaults to 0.02
        :type t_error_value: TYPE, optional
        :param t_error_type: DESCRIPTION, defaults to "absolute"
        :type t_error_type: TYPE, optional
        :param t_floor: DESCRIPTION, defaults to True
        :type t_floor: TYPE, optional
        :param : DESCRIPTION
        :type : TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        for mt_obj in self.values():
            mt_obj.compute_model_z_errors(z_error_value, z_error_type, z_floor)
            mt_obj.compute_model_t_errors(t_error_value, t_error_type, t_floor)

    def estimate_starting_rho(self):
        """
        Estimate starting resistivity from the data.
        Creates a plot of the mean and median apparent resistivity values.

        :return: array of the median rho per period
        :rtype: np.ndarray(n_periods)
        :return: array of the mean rho per period
        :rtype: np.ndarray(n_periods)

        >>> d = Data()
        >>> d.read_data_file(r"example/data.dat")
        >>> rho_median, rho_mean = d.estimate_starting_rho()

        """

        entries = []
        for mt_obj in self.values():
            for period, res_det in zip(mt_obj.period, mt_obj.Z.res_det):
                entries.append({"period": period, "res_det": res_det})

        res_df = pd.DataFrame(entries)

        mean_rho = res_df.groupby("period").mean()
        median_rho = res_df.groupby("period").median()

        fig = plt.figure()

        ax = fig.add_subplot(1, 1, 1)
        (l1,) = ax.loglog(
            mean_rho.index, mean_rho.res_det, lw=2, color=(0.75, 0.25, 0)
        )
        (l2,) = ax.loglog(
            median_rho.index, median_rho.res_det, lw=2, color=(0, 0.25, 0.75)
        )

        ax.loglog(
            mean_rho.index,
            np.repeat(mean_rho.res_det.mean(), mean_rho.shape[0]),
            ls="--",
            lw=2,
            color=(0.75, 0.25, 0),
        )
        ax.loglog(
            median_rho.index,
            np.repeat(median_rho.res_det.median(), median_rho.shape[0]),
            ls="--",
            lw=2,
            color=(0, 0.25, 0.75),
        )

        ax.set_xlabel("Period (s)", fontdict={"size": 12, "weight": "bold"})
        ax.set_ylabel(
            "Resistivity (Ohm-m)", fontdict={"size": 12, "weight": "bold"}
        )

        ax.legend(
            [l1, l2],
            [
                f"Mean = {mean_rho.res_det.mean():.1f}",
                f"Median = {median_rho.res_det.median():.1f}",
            ],
            loc="upper left",
        )
        ax.grid(which="both", ls="--", color=(0.75, 0.75, 0.75))
        ax.set_xlim((res_df.period.min(), res_df.period.max()))

        plt.show()
