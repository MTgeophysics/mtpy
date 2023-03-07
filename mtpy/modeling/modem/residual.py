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

        super().__init__(self, **kwargs)

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
    def rms_by_period_all(self):
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
                        "z_rms": np.nanmean(z_df),
                        "t_rms": np.nanmean(t_df),
                    }
                )

            return pd.DataFrame(rms_list)

    def get_rms(self, residual_fn=None):

        if residual_fn is None:
            residual_fn = self.residual_fn

        if self.residual_array is None:
            self.read_residual_file(residual_fn)
        if self.residual_array is None:
            return

        rms_z_comp = np.zeros((len(self.rms_array), 2, 2))
        rms_tip_comp = np.zeros((len(self.rms_array), 2))
        rms_value_list_all = np.zeros(0)
        rms_value_list_z = np.zeros(0)
        rms_value_list_tip = np.zeros(0)

        for station_name in self.rms_array["station"]:
            rms_value_list = []
            rms_value_list_ztip = np.zeros(
                (self.residual_array["z"].shape[1], 6)
            )
            sta_ind = np.where(self.rms_array["station"] == station_name)[0][0]
            sta_indd = np.where(self.residual_array["station"] == station_name)[
                0
            ][0]
            res_vals = self.residual_array[sta_indd]
            z_norm, tip_norm = None, None
            if np.amax(np.abs(res_vals["z"])) > 0:
                # sum over absolute value of z
                # need to divide by sqrt(2) to normalise (ModEM data file has one error for real and imag components)
                z_norm = np.abs(res_vals["z"]) / (
                    np.real(res_vals["z_err"]) * 2.0**0.5
                )
                rms_value_list_ztip[:, :4] = z_norm.reshape(z_norm.shape[0], 4)

                # count number of values for tipper, all
                count_z = np.count_nonzero(
                    np.nan_to_num(z_norm.reshape(z_norm.shape[0], 4)), axis=1
                ).astype(float)
                #                count_z = count_z.reshape(count_z.shape[0],1)

                # normalized error split by period
                self.rms_array["rms_z_period"][sta_ind] = (
                    np.sum(z_norm.reshape(z_norm.shape[0], 4) ** 2, axis=1)
                    / count_z
                ) ** 0.5

                z_norm_nz = z_norm[np.all(np.isfinite(z_norm), axis=(1, 2))]

                # append individual normalised errors to a master list for all stations
                rms_value_list_all = np.append(
                    rms_value_list_all, z_norm_nz.flatten()
                )
                rms_value_list_z = np.append(
                    rms_value_list_z, z_norm_nz.flatten()
                )

                # normalised error for separate components
                rms_z_comp[sta_ind] = (
                    ((z_norm_nz**2.0).sum(axis=0)) / (z_norm_nz.shape[0])
                ) ** 0.5
                rms_value_list.append(rms_z_comp[sta_ind])

            if np.amax(np.abs(res_vals["tip"])) > 0:
                # sum over absolute value of tipper
                # need to divide by sqrt(2) to normalise (code applies same error to real and imag components)
                tip_norm = np.abs(res_vals["tip"]) / (
                    np.real(res_vals["tip_err"]) * 2.0**0.5
                )
                rms_value_list_ztip[:, 4:] = tip_norm.reshape(
                    tip_norm.shape[0], 2
                )

                # count number of values for tipper, all
                count_tip = np.count_nonzero(
                    np.nan_to_num(tip_norm.reshape(tip_norm.shape[0], 2)),
                    axis=1,
                ).astype(float)

                # normalized error split by period
                self.rms_array["rms_tip_period"][sta_ind] = (
                    np.nansum(
                        tip_norm.reshape(tip_norm.shape[0], 2) ** 2, axis=1
                    )
                    / count_tip
                ) ** 0.5

                tip_norm_nz = tip_norm[
                    np.all(np.isfinite(tip_norm), axis=(1, 2))
                ]

                # append individual normalised errors to a master list for all stations
                rms_value_list_all = np.append(
                    rms_value_list_all, tip_norm_nz.flatten()
                )
                rms_value_list_tip = np.append(
                    rms_value_list_tip, tip_norm_nz.flatten()
                )

                # normalised error for separate components
                rms_tip_comp[sta_ind] = (
                    ((tip_norm_nz**2.0).sum(axis=0)) / len(tip_norm_nz)
                ) ** 0.5
                rms_value_list.append(rms_tip_comp[sta_ind])

            # compute overall rms by period
            count_ztip = np.count_nonzero(
                np.nan_to_num(rms_value_list_ztip), axis=1
            ).astype(float)
            self.rms_array["rms_period"][sta_ind] = (
                np.nansum(rms_value_list_ztip**2, axis=1) / count_ztip
            ) ** 0.5

            rms_value_list = np.vstack(rms_value_list).flatten()

            rms_value = (
                (rms_value_list**2.0).sum() / rms_value_list.size
            ) ** 0.5

            self.rms_array[sta_ind]["rms"] = rms_value

            if z_norm is not None:
                self.rms_array[sta_ind]["rms_z"] = (
                    (rms_z_comp[sta_ind] ** 2.0).sum()
                    / rms_z_comp[sta_ind].size
                ) ** 0.5
            if tip_norm is not None:
                self.rms_array[sta_ind]["rms_tip"] = (
                    (rms_tip_comp[sta_ind] ** 2.0).sum()
                    / rms_z_comp[sta_ind].size
                ) ** 0.5

        self.rms = np.mean(rms_value_list_all**2.0) ** 0.5
        self.rms_z = np.mean(rms_value_list_z**2.0) ** 0.5
        self.rms_tip = np.mean(rms_value_list_tip**2.0) ** 0.5

        # by component
        for cpt in ["z", "tip"]:
            # normalised residuals
            res_vals_cpt = np.abs(self.residual_array[cpt]) / (
                np.real(self.residual_array[cpt + "_err"]) * 2.0**0.5
            )
            ijvals = res_vals_cpt.shape[2:]
            for i in range(ijvals[0]):
                for j in range(ijvals[1]):
                    self.rms_array["rms_{}_component".format(cpt)][:, i, j] = (
                        np.nansum(res_vals_cpt[:, :, i, j] ** 2.0, axis=1)
                        / np.nansum(
                            np.isfinite(res_vals_cpt[:, :, i, j]), axis=1
                        )
                    ) ** 0.5
                    self.rms_array["rms_{}_component_period".format(cpt)][
                        :, :, i, j
                    ] = (
                        res_vals_cpt[:, :, i, j] ** 2
                        / np.isfinite(res_vals_cpt[:, :, i, j])
                    ) ** 0.5

    # def write_rms_to_file(self):
    #     """
    #     write rms station data to file
    #     """

    #     fn = op.join(self.work_dir, "rms_values.dat")

    #     if not hasattr(self, "rms"):
    #         self.get_rms()

    #     header_list = [
    #         "station",
    #         "lon",
    #         "lat",
    #         "rel_east",
    #         "rel_north",
    #         "rms",
    #         "rms_z",
    #         "rms_tip",
    #     ]

    #     dtype = []
    #     for val in header_list:
    #         if val == "station":
    #             dtype.append((val, "S10"))
    #         else:
    #             dtype.append((val, np.float))

    #     save_list = np.zeros(len(self.rms_array), dtype=dtype)
    #     for val in header_list:
    #         save_list[val] = self.rms_array[val]

    #     header = " ".join(header_list)

    #     np.savetxt(
    #         fn,
    #         save_list,
    #         header=header,
    #         fmt=["%s", "%.6f", "%.6f", "%.1f", "%.1f", "%.3f", "%.3f", "%.3f"],
    #     )
