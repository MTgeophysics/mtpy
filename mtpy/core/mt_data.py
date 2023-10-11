# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 11:58:56 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from pathlib import Path
from collections import OrderedDict
from copy import deepcopy


import numpy as np
import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt

from .mt import MT
from .mt_stations import MTStations

from mtpy.modeling.errors import ModelErrors
from mtpy.modeling.modem import Data
from mtpy.imaging import (
    PlotStations,
    PlotMultipleResponses,
    PlotPhaseTensorMaps,
    PlotPhaseTensorPseudoSection,
    PlotStrike,
    PlotPenetrationDepthMap,
    PlotResPhaseMaps,
    PlotResPhasePseudoSection,
    PlotResidualPTMaps,
)

from mth5.helpers import validate_name

# =============================================================================


class MTData(OrderedDict, MTStations):
    """
    Collection of MT objects as an OrderedDict where keys are formatted as
    `survey_id.station_id`.  Has all functionality of an OrderedDict for
    example can iterate of `.keys()`, `.values()` or `.items`.  Values are
    a list of MT objects.

    Inherits :class:`mtpyt.core.MTStations` to deal with geographic locations
    of stations.

    Is not optimized yet for speed, works fine for smaller surveys, but for
    large can be slow.  Might try using a dataframe as the base.
    """

    def __init__(self, mt_list=None, **kwargs):

        if mt_list is not None:
            for mt_obj in mt_list:
                self.add_station(mt_obj, compute_relative_location=False)

        MTStations.__init__(self, None, None, **kwargs)

        self.z_model_error = ModelErrors(
            error_value=5,
            error_type="geometric_mean",
            floor=True,
            mode="impedance",
        )
        self.t_model_error = ModelErrors(
            error_value=0.02,
            error_type="absolute",
            floor=True,
            mode="tipper",
        )
        self.data_rotation_angle = 0

        self.model_parameters = {}

        self._copy_attrs = [
            "z_model_error",
            "t_model_error",
            "utm_crs",
            "datum_crs",
            "_center_lat",
            "_center_lon",
            "_center_elev",
            "shift_east",
            "shift_north",
            "rotation_angle",
            "data_rotation_angle",
            "model_parameters",
        ]

    def _validate_item(self, mt_obj):
        """
        Make sure intpu is an MT object.  If the input is a string assume
        its a file path.

        :param mt_obj: Input MT object
        :type mt_obj: :class:`mtpy.MT`, str, :class:`pathlib.Path`
        :raises TypeError: If not any of the input types
        :return: MT object
        :rtype: :class:`mtpy.MT`

        """
        if isinstance(mt_obj, (str, Path)):
            m = MT()
            m.read(mt_obj)
            return m

        elif not isinstance(mt_obj, MT):
            raise TypeError(
                f"entry must be a mtpy.core.MT object not type({type(mt_obj)})"
            )
        return mt_obj

    def __eq__(self, other):
        """
        test for other is equal

        :param other: Other MTData object
        :type other: :class:`mtpy.MTData`
        :return: True if equal, False if not equal
        :rtype: bool

        """

        if not isinstance(other, MTData):
            raise TypeError(f"Can not compare MTData to {type(other)}.")

        for attr in self._copy_attrs:
            value_og = getattr(self, attr)
            value_other = getattr(other, attr)

            if value_og != value_other:
                self.logger.info(
                    f"Attribute {attr}: {value_og} != {value_other}"
                )
                return False
        fail = False
        if len(self) == len(other):
            for key in self.keys():
                mt1 = self[key]
                try:
                    mt2 = other[key]
                    if mt1 != mt2:
                        self.logger.info(f"Station {key} is not equal.")
                        fail = True
                except KeyError:
                    self.logger.info(f"Could not find {key} in other.")
                    fail = True
            if fail:
                return False
        else:
            self.logger.info(
                f"Length of MTData not the same {len(self)} != {len(other)}"
            )
            return False
        return True

    def __deepcopy__(self, memo):
        """
        Deep copy overwrite to make sure that logger is skipped.
        :param memo: DESCRIPTION
        :type memo: TYPE
        :return: Deep copy of original MTData
        :rtype: :class:`mtpy.MTData`

        """
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for key in self._copy_attrs:
            value = getattr(self, key)
            setattr(result, key, deepcopy(value, memo))

        for mt_obj in self.values():
            result.add_station(mt_obj.copy(), compute_relative_location=False)

        return result

    def copy(self):
        """
        Deep copy of original MTData object

        :param memo: DESCRIPTION
        :type memo: TYPE
        :return: Deep copy of original MTData
        :rtype: :class:`mtpy.MTData`

        """
        copied = deepcopy(self)
        copied.logger = self.logger
        return copied

    def clone_empty(self):
        """
        Return a copy of MTData excluding MT objects.

        :return: Copy of MTData object excluding MT objects
        :rtype: :class:`mtpy.MTData`

        """

        md = MTData()
        for attr in self._copy_attrs:
            setattr(self, attr, deepcopy(getattr(self, attr)))

        return md

    @property
    def mt_list(self):
        """

        :return: List of MT objects
        :rtype: list

        """
        return self.values()

    @mt_list.setter
    def mt_list(self, value):
        """
        At this point not implemented, mainly here for inheritance of MTStations
        """
        pass

    @property
    def survey_ids(self):
        """
        Survey IDs for all MT objects

        :return: list of survey IDs
        :rtype: list

        """
        return list(set([key.split(".")[0] for key in self.keys()]))

    def get_survey(self, survey_id):
        """
        Get all MT objects that belong to the 'survey_id' from the data set.

        :param survey_id: survey ID
        :type survey_id: str
        :return: MTData object including only those with the desired 'survey_id'
        :rtype: :class:`mtpy.MTData`

        """

        survey_list = [
            mt_obj for key, mt_obj in self.items() if survey_id in key
        ]
        return MTData(survey_list)

    def add_station(
        self,
        mt_object,
        survey=None,
        compute_relative_location=True,
        interpolate_periods=None,
    ):
        """
        Add a MT object

        :param mt_object: MT object for a single station
        :type mt_object: :class:`mtpy.MT`
        :param survey: new survey name, defaults to None
        :type survey: str, optional
        :param compute_relative_location: Compute relative location,
         can be slow if adding single stations in a loop.  If looping over
         station set to False and compute at the end, defaults to True
        :type compute_relative_location: bool, optional
        :param interpolate_periods: periods to interpolate onto, defaults to None
        :type interpolate_periods: np.array, optional

        """

        if not isinstance(mt_object, (list, tuple)):
            mt_object = [mt_object]

        for m in mt_object:
            m = self._validate_item(m)
            if self.utm_crs is not None:
                m.utm_crs = self.utm_crs
            if survey is not None:
                m.survey = survey

            if interpolate_periods is not None:
                if not isinstance(interpolate_periods, np.ndarray):
                    interpolate_periods = np.array(interpolate_periods)

                m = m.interpolate(interpolate_periods, bounds_error=False)

            self.__setitem__(f"{validate_name(m.survey)}.{m.station}", m)

        if compute_relative_location:
            self.compute_relative_locations()

    def add_tf(self, tf, **kwargs):
        """
        Add a MT object

        :param mt_object: MT object for a single station
        :type mt_object: :class:`mtpy.MT`
        :param survey: new survey name, defaults to None
        :type survey: str, optional
        :param compute_relative_location: Compute relative location,
         can be slow if adding single stations in a loop.  If looping over
         station set to False and compute at the end, defaults to True
        :type compute_relative_location: bool, optional
        :param interpolate_periods: periods to interpolate onto, defaults to None
        :type interpolate_periods: np.array, optional

        """
        self.add_station(tf, **kwargs)

    def remove_station(self, station_id, survey_id=None):
        """
        remove a station from the dictionary based on the key

        :param station_id: station ID
        :type station_id: str
        :param survey_id: survey ID
        :type survey_id: str

        """

        key = self._get_station_key(station_id, survey_id)
        if key in self.keys():
            del self[key]

    def _get_station_key(self, station_id, survey_id):
        """
        get station key from station id and survey id

        :param station_id: station ID
        :type station_id: str
        :param survey_id: survey ID
        :type survey_id: str
        :return: station key
        :rtype: str

        """

        if station_id is not None:
            if survey_id is not None:
                return f"{validate_name(survey_id)}.{station_id}"
            else:
                for key in self.keys():
                    if station_id in key:
                        if key.split(".")[1] == station_id:
                            return key
        raise KeyError(
            f"Could not find station_id = {station_id}, survey_id = {survey_id}"
        )

    def get_station(self, station_id=None, survey_id=None, station_key=None):
        """

        if 'station_key' is None, tries to find key from `station_id` and
        'survey_id' using MTData._get_station_key()

        :param station_key: full station key {survey_id}.{station_id},
         defaults to None
        :type station_key: str, optional
        :param station_id: station ID, defaults to None
        :type station_id: str, optional
        :param survey_id: survey ID, defaults to None
        :type survey_id: str, optional
        :raises KeyError: If cannot find station_key
        :return: MT object
        :rtype: :class:`mtpy.MT`

        """
        if station_key is not None:
            station_key = station_key
        else:
            station_key = self._get_station_key(
                station_id, validate_name(survey_id)
            )

        try:
            return self[station_key]
        except KeyError:
            raise KeyError(f"Could not find {station_key} in MTData.")

    def get_subset(self, station_list):
        """
        get a subset of the data from a list of stations, could be station_id
        or station_keys

        :param station_list: list of station keys as {survey_id}.{station_id}
        :type station_list: list
        :return: Returns just those stations within station_list
        :rtype: :class:`mtpy.MTData`

        """
        mt_data = MTData()
        for station in station_list:
            if station.count(".") > 0:
                mt_data.add_station(
                    self.get_station(station_key=station),
                    compute_relative_location=False,
                )
            else:
                mt_data.add_station(
                    self.get_station(station_id=station),
                    compute_relative_location=False,
                )

        return mt_data

    @property
    def n_stations(self):
        """number of stations in MT data"""

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
            mt_obj.to_dataframe(utm_crs=utm_crs, cols=cols).dataframe
            for mt_obj in self.values()
        ]

        df = pd.concat(df_list)
        df.reset_index(drop=True, inplace=True)
        return df

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
            self.add_station(mt_object, compute_relative_location=False)

    def to_geo_df(self):
        """
        Make a geopandas dataframe for easier GIS manipulation

        """

        df = self.station_locations

        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs=self.datum_crs,
        )

        return gdf

    def interpolate(self, new_periods, f_type="period", inplace=True):
        """
        Interpolate onto common period range

        :param new_periods: DESCRIPTION
        :type new_periods: TYPE
        :param f_type: frequency type can be [ 'frequency' | 'period' ]
        :type f_type: string, defaults to 'period'
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not inplace:
            mt_data = self.clone_empty()

        for mt_obj in self.values():
            interp_periods = new_periods[
                np.where(
                    (new_periods <= mt_obj.period.max())
                    & (new_periods >= mt_obj.period.min())
                )
            ]

            new_mt_obj = mt_obj.interpolate(interp_periods, f_type=f_type)

            if inplace:
                self.update(
                    {
                        f"{new_mt_obj.survey_metadata.id}.{new_mt_obj.station}": new_mt_obj
                    }
                )

            else:
                mt_data.add_station(
                    new_mt_obj, compute_relative_location=False
                )

        if not inplace:
            return mt_data

    def rotate(self, rotation_angle, inplace=True):
        """
        rotate the data by the given angle assuming positive clockwise with
        north = 0, east = 90.

        :param rotation_angle: DESCRIPTION
        :type rotation_angle: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.data_rotation_angle = rotation_angle

        if not inplace:
            mt_data = self.clone_empty
        for mt_obj in self.values():
            if not inplace:
                rot_mt_obj = mt_obj.copy()
                rot_mt_obj.rotation_angle = rotation_angle
                mt_data.add_station(
                    rot_mt_obj, compute_relative_location=False
                )
            else:
                mt_obj.rotation_angle = rotation_angle

        if not inplace:
            return mt_data

    def get_profile(self, x1, y1, x2, y2, radius):
        """
        Get stations along a profile line given the (x1, y1) and (x2, y2)
        coordinates within a given radius (in meters).

        These can be in (longitude, latitude) or (easting, northing).
        The calculation is done in UTM, therefore a UTM CRS must be input

        :param x1: DESCRIPTION
        :type x1: TYPE
        :param y1: DESCRIPTION
        :type y1: TYPE
        :param x2: DESCRIPTION
        :type x2: TYPE
        :param y2: DESCRIPTION
        :type y2: TYPE
        :param radius: DESCRIPTION
        :type radius: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        profile_stations = self._extract_profile(x1, y1, x2, y2, radius)

        mt_data = self.clone_empty()
        for mt_obj in profile_stations:
            mt_data.add_station(mt_obj, compute_relative_location=False)

        return mt_data

    def compute_model_errors(
        self,
        z_error_value=None,
        z_error_type=None,
        z_floor=None,
        t_error_value=None,
        t_error_type=None,
        t_floor=None,
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

        if z_error_value is not None:
            self.z_model_error.error_value = z_error_value
        if z_error_type is not None:
            self.z_model_error.error_type = z_error_type
        if z_floor is not None:
            self.z_model_error.floor = z_floor

        if t_error_value is not None:
            self.t_model_error.error_value = t_error_value
        if t_error_type is not None:
            self.t_model_error.error_type = t_error_type
        if t_floor is not None:
            self.t_model_error.floor = t_floor

        for mt_obj in self.values():
            mt_obj.compute_model_z_errors(
                **self.z_model_error.error_parameters
            )
            mt_obj.compute_model_t_errors(
                **self.t_model_error.error_parameters
            )

    def get_nearby_stations(self, station_key, radius, radius_units="m"):
        """
        get stations close to a given station

        :param station_key: DESCRIPTION
        :type station_key: TYPE
        :param radius: DESCRIPTION
        :type radius: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        # get the local station
        local_station = self.get_station(station_key=station_key)

        sdf = self.station_locations.copy()
        if radius_units in ["m", "meters", "metres"]:
            if self.utm_crs is None:
                raise ValueError(
                    "Cannot estimate distances in meters without a UTM CRS. Set 'utm_crs' first."
                )
            sdf["radius"] = np.sqrt(
                (local_station.east - sdf.east) ** 2
                + (local_station.north - sdf.north)
            )
        elif radius_units in ["deg", "degrees"]:
            sdf["radius"] = np.sqrt(
                (local_station.longitude - sdf.longitude) ** 2
                + (local_station.latitude - sdf.latitude)
            )

        return [
            f"{row.survey}.{row.station}"
            for row in sdf.loc[
                (sdf.radius <= radius) & (sdf.radius > 0)
            ].itertuples()
        ]

    def estimate_spatial_static_shift(
        self,
        station_key,
        radius,
        period_min,
        period_max,
        radius_units="m",
        shift_tolerance=0.15,
    ):
        """
        Estimate static shift for a station by estimating the median resistivity
        values for nearby stations within a radius given.  Can set the period
        range to estimate the resistivity values.

        :param station_key: DESCRIPTION
        :type station_key: TYPE
        :param radius: DESCRIPTION
        :type radius: TYPE
        :param period_min: DESCRIPTION
        :type period_min: TYPE
        :param period_max: DESCRIPTION
        :type period_max: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        md = self.get_subset(
            self.get_nearby_stations(station_key, radius, radius_units)
        )
        if len(md) == 0:
            self.logger.warning(
                f"Could not find any nearby stations for {station_key}."
            )
            return 1.0, 1.0

        local_site = self.get_station(station_key=station_key)

        interp_periods = local_site.period[
            np.where(
                (local_site.period >= period_min)
                & (local_site.period <= period_max)
            )
        ]

        local_site = local_site.interpolate(interp_periods)
        md.interpolate(interp_periods)

        df = md.to_dataframe()

        sx = np.nanmedian(df.res_xy) / np.nanmedian(local_site.Z.res_xy)
        sy = np.nanmedian(df.res_yx) / np.nanmedian(local_site.Z.res_yx)

        # check to see if the estimated static shift is within given tolerance
        if 1 - shift_tolerance < sx and sx < 1 + shift_tolerance:
            sx = 1.0
        # check to see if the estimated static shift is within given tolerance
        if 1 - shift_tolerance < sy and sy < 1 + shift_tolerance:
            sy = 1.0

        return sx, sy

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

    def to_modem_data(self, data_filename=None, **kwargs):
        """
        Create a modem data file

        :param data_filename: DESCRIPTION
        :type data_filename: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        modem_kwargs = dict(self.model_parameters)
        modem_kwargs.update(kwargs)

        if np.all(self.station_locations.model_east == 0):
            if self.utm_crs is None:
                raise ValueError(
                    "Need to input data UTM EPSG or CRS to compute relative "
                    "station locations"
                )
            self.compute_relative_locations()

        modem_data = Data(
            dataframe=self.to_dataframe(),
            center_point=self.center_point,
            **modem_kwargs,
        )
        modem_data.z_model_error = self.z_model_error
        modem_data.t_model_error = self.t_model_error
        if data_filename is not None:
            modem_data.write_data_file(file_name=data_filename)

        return modem_data

    def from_modem_data(self, data_filename, file_type="data", **kwargs):
        """
        read in a modem data file

        :param data_filename: DESCRIPTION
        :type data_filename: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        modem_data = Data(**kwargs)
        mdf = modem_data.read_data_file(data_filename)
        if file_type in ["data"]:
            mdf.dataframe["survey"] = "data"
        elif file_type in ["response", "model"]:
            mdf.dataframe["survey"] = "model"
        self.from_dataframe(mdf.dataframe)
        self.z_model_error = ModelErrors(
            mode="impedance", **modem_data.z_model_error.error_parameters
        )
        self.t_model_error = ModelErrors(
            mode="tipper", **modem_data.t_model_error.error_parameters
        )
        self.data_rotation_angle = modem_data.rotation_angle
        self._center_lat = modem_data.center_point.latitude
        self._center_lon = modem_data.center_point.longitude
        self._center_elev = modem_data.center_point.elevation
        self.utm_epsg = modem_data.center_point.utm_epsg

        self.model_parameters = dict(
            [
                (key, value)
                for key, value in modem_data.model_parameters.items()
                if "." not in key
            ]
        )

    def add_white_noise(self, value, inplace=True):
        """
        Add white noise to the data, useful for synthetic tests.

        :param value: DESCRIPTION
        :type value: TYPE
        :param inplace: DESCRIPTION, defaults to True
        :type inplace: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if value > 1:
            value = value / 100.0

        if not inplace:
            mt_list = []

        for station, mt_obj in self.items():
            if inplace:
                mt_obj.add_white_noise(value)

            else:
                mt_list.append(mt_obj.add_white_noise(value, inplace=False))

        if not inplace:
            return_data = self.clone_empty()
            return_data.add_station(mt_list)
            return return_data

    def plot_mt_response(
        self, station_key=None, station_id=None, survey_id=None, **kwargs
    ):
        """

        :param tf_id: DESCRIPTION
        :type tf_id: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        if input as list, tuple, np.ndarray, pd.series assuming first column
        is tf_id, and if needed the second column should be the survey id
        for that tf.

        """

        if isinstance(station_key, (list, tuple)):
            mt_data = MTData()
            for sk in station_key:
                mt_data.add_station(
                    self.get_station(station_key=sk),
                    compute_relative_location=False,
                )
            return PlotMultipleResponses(mt_data, **kwargs)

        elif isinstance(station_id, (list, tuple)):
            mt_data = MTData()
            if isinstance(survey_id, (list, tuple)):
                if len(survey_id) != len(station_key):
                    raise ValueError(
                        "Number of survey must match number of stations"
                    )
            elif isinstance(survey_id, (str, type(None))):
                survey_id = [survey_id] * len(station_id)
            for survey, station in zip(survey_id, station_id):
                mt_data.add_station(
                    self.get_station(station_id=station, survey_id=survey),
                    compute_relative_location=False,
                )
            return PlotMultipleResponses(mt_data, **kwargs)

        else:
            mt_object = self.get_station(
                station_id=station_id,
                survey_id=survey_id,
                station_key=station_key,
            )
            return mt_object.plot_mt_response(**kwargs)

    def plot_stations(self, map_epsg=4326, bounding_box=None, **kwargs):
        """
        plot stations

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        gdf = self.to_geo_df()
        return PlotStations(gdf, **kwargs)

    def plot_strike(self, **kwargs):
        """
        Plot strike angle

        .. seealso:: :class:`mtpy.imaging.PlotStrike`
        """

        return PlotStrike(self, **kwargs)

    def plot_phase_tensor(
        self, station_key=None, station_id=None, survey_id=None, **kwargs
    ):
        """
        plot phase tensor elements

        :param tf_id: DESCRIPTION
        :type tf_id: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        mt_object = self.get_station(
            station_id=station_id,
            survey_id=survey_id,
            station_key=station_key,
        )
        return mt_object.plot_phase_tensor(**kwargs)

    def plot_phase_tensor_map(self, **kwargs):
        """
        Plot Phase tensor maps for transfer functions in the working_dataframe

        .. seealso:: :class:`mtpy.imaging.PlotPhaseTensorMaps`

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotPhaseTensorMaps(mt_data=self, **kwargs)

    def plot_phase_tensor_pseudosection(self, mt_data=None, **kwargs):
        """
        Plot a pseudo section of  phase tensor ellipses and induction vectors
        if specified

        .. seealso:: :class:`mtpy.imaging.PlotPhaseTensorPseudosection`

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotPhaseTensorPseudoSection(mt_data=self, **kwargs)

    def plot_penetration_depth_1d(
        self, station_key=None, station_id=None, survey_id=None, **kwargs
    ):
        """
        Plot 1D penetration depth based on the Niblett-Bostick transformation

        Note that data is rotated to estimated strike previous to estimation
        and strike angles are interpreted for data points that are 3D.

        .. seealso:: :class:`mtpy.analysis.niblettbostick.calculate_depth_of_investigation`


        :param tf_object: DESCRIPTION
        :type tf_object: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        mt_object = self.get_station(
            station_id=station_id,
            survey_id=survey_id,
            station_key=station_key,
        )

        return mt_object.plot_depth_of_penetration(**kwargs)

    def plot_penetration_depth_map(self, **kwargs):
        """
        Plot Penetration depth in map view for a single period

        .. seealso:: :class:`mtpy.imaging.PlotPenetrationDepthMap`

        :param mt_data: DESCRIPTION
        :type mt_data: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotPenetrationDepthMap(self, **kwargs)

    def plot_resistivity_phase_maps(self, **kwargs):
        """
        Plot apparent resistivity and/or impedance phase maps from the
        working dataframe

        .. seealso:: :class:`mtpy.imaging.PlotResPhaseMaps`
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotResPhaseMaps(self, **kwargs)

    def plot_resistivity_phase_pseudosections(self, **kwargs):
        """
        Plot resistivity and phase in a pseudosection along a profile line

        :param mt_data: DESCRIPTION, defaults to None
        :type mt_data: TYPE, optional
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        return PlotResPhasePseudoSection(self, **kwargs)

    def plot_residual_phase_tensor_maps(self, survey_01, survey_02, **kwargs):
        """

        :param survey_01: DESCRIPTION
        :type survey_01: TYPE
        :param survey_02: DESCRIPTION
        :type survey_02: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        survey_data_01 = self.get_survey(survey_01)
        survey_data_02 = self.get_survey(survey_02)

        return PlotResidualPTMaps(survey_data_01, survey_data_02, **kwargs)
