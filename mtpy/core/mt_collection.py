# -*- coding: utf-8 -*-
"""
Collection of MT stations

Created on Mon Jan 11 15:36:38 2021

:copyright: 
    Jared Peacock (jpeacock@usgs.gov)

:license: MIT

"""
# =============================================================================
# Imports
# =============================================================================
from pathlib import Path

from loguru import logger
import numpy as np
import pandas as pd
import geopandas as gpd

from mtpy import MT
from mtpy.core.mt_data import MTData
from mtpy.imaging import (
    PlotStations,
    PlotMultipleResponses,
    PlotResidualPTMaps,
    PlotResidualPTPseudoSection,
    PlotPhaseTensorMaps,
    PlotPhaseTensorPseudoSection,
    PlotStrike,
    PlotPenetrationDepth1D,
    PlotPenetrationDepthMap,
    PlotResPhaseMaps,
    PlotResPhasePseudoSection,
)

from mth5.mth5 import MTH5

# =============================================================================
#
# =============================================================================


class MTCollection:
    """
    Collection of transfer functions

    The main working variable is `MTCollection.dataframe` which is a property
    that returns either the `master dataframe` that contains all the TF's in
    the MTH5 file, or the `working_dataframe` which is a dataframe that has
    been queried in some way.  Therefore all the user has to do is set
    the working directory as a subset of the master_dataframe

    :Example:

        >>> mc = MTCollection()
        >>> mc.open_collection(filename="path/to/example/mth5.h5")
        >>> mc.working_dataframe = mc.master_dataframe.iloc[0:5]

    """

    def __init__(self, working_directory=None):

        self._cwd = Path().cwd()
        self.mth5_basename = "mt_collection"
        self.working_directory = working_directory
        self.working_dataframe = None

        self.mth5_collection = MTH5()

        self._added = False

        self.logger = logger

    def __str__(self):
        lines = [f"Working Directory: {self.working_directory}"]
        lines.append(f"MTH5 file:         {self.mth5_filename}")
        if self.mth5_collection.h5_is_read():
            lines.append(
                f"\tNumber of Transfer Functions: {len(self.dataframe)}"
            )
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close_collection()
        return False

    @property
    def working_directory(self):
        return self._cwd

    @working_directory.setter
    def working_directory(self, value):
        if value is None:
            return
        value = Path(value)
        if not value.exists():
            raise IOError(f"could not find directory {value}")
        self._cwd = value

    @property
    def mth5_filename(self):
        if self.mth5_basename.find(".h5") > 0:
            return self.working_directory.joinpath(f"{self.mth5_basename}")
        else:
            return self.working_directory.joinpath(f"{self.mth5_basename}.h5")

    @mth5_filename.setter
    def mth5_filename(self, value):
        value = Path(value)
        self.working_directory = value.parent
        self.mth5_basename = value.stem

    @property
    def master_dataframe(self):
        """
        This is the full summary of all transfer functions in the MTH5
        file.  It is a property because if a user adds TF's then the
        master_df will be automatically updated.  the tranformation is quick
        for now.

        """

        if self.mth5_collection.h5_is_read():
            return self.mth5_collection.tf_summary.to_dataframe()
        return None

    @property
    def dataframe(self):
        """
        This property returns the working dataframe or master dataframe if
        the working dataframe is None.

        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.working_dataframe is None:
            return self.master_dataframe
        elif isinstance(self.working_dataframe, pd.DataFrame):
            if self.working_dataframe.empty:
                return None
            else:
                return self.working_dataframe
        return self.working_dataframe

    def has_data(self):
        if self.master_dataframe is not None:
            return True
        return False

    @staticmethod
    def make_file_list(mt_path, file_types=["edi"]):
        """
        Get a list of MT file from a given path

        :param mt_path: full path to where the MT transfer functions are stored
        or a list of paths
        :type mt_path: string or :class:`pathlib.Path` or list

        :param file_types: List of file types to look for given their extension
        :type file_types: list

        Currently available file types are or will be:
            - edi - EDI files
            - zmm - EMTF output file
            - j - BIRRP output file
            - avg - Zonge output file

        """

        def check_path(mt_path):
            if mt_path is None:
                return None
            else:
                mt_path = Path(mt_path)
                if not mt_path.exists():
                    msg = f"{mt_path} does not exists"
                    raise IOError(msg)
                return mt_path

        if isinstance(mt_path, (str, Path)):
            mt_path = [check_path(mt_path)]
        elif isinstance(mt_path, list):
            mt_path = [check_path(path) for path in mt_path]
        else:
            raise TypeError(f"Not sure what to do with {type(mt_path)}")
        fn_list = []
        for path in mt_path:
            for ext in file_types:
                fn_list += list(path.glob(f"*.{ext}"))
        return fn_list

    def open_collection(
        self,
        filename=None,
        basename=None,
        working_directory=None,
        mode="a",
    ):
        """
        Initialize an mth5

        :param basename: DESCRIPTION, defaults to "mt_collection"
        :type basename: TYPE, optional
        :param working_directory: DESCRIPTION, defaults to None
        :type working_directory: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if filename is not None:
            self.mth5_filename = filename

        if basename is not None:
            self.mth5_basename = basename

        if working_directory is not None:
            self.working_directory = working_directory

        self.mth5_collection.open_mth5(self.mth5_filename, mode)

    def close_collection(self):
        """
        close mth5

        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.mth5_collection.close_mth5()

    def add_tf(self, transfer_function):
        """
        transfer_function could be a transfer function object, a file name,
        a list of either.
        g
        :param transfer_function: DESCRIPTION
        :type transfer_function: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if not isinstance(
            transfer_function, (list, tuple, np.ndarray, MTData)
        ):
            transfer_function = [transfer_function]
        for item in transfer_function:
            if isinstance(item, MT):
                self._from_mt_object(item)
            elif isinstance(item, (str, Path)):
                self._from_file(item)
            else:
                raise TypeError(f"Not sure want to do with {type(item)}.")
        self.mth5_collection.tf_summary.summarize()

    def get_tf(self, tf_id, survey=None):
        """

        Get transfer function

        :param tf_id: DESCRIPTION
        :type tf_id: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if survey is None:
            try:
                find_df = self.master_dataframe.loc[
                    self.master_dataframe.tf_id == tf_id
                ]
                find_df = find_df.iloc[0]
                self.logger.warning(
                    f"Found multiple transfer functions with ID {tf_id}. "
                    "Suggest setting survey, otherwise returning the "
                    f"TF from survey {find_df.survey}."
                )
            except IndexError:
                raise ValueError(f"Could not find {tf_id} in collection.")
        else:
            try:
                find_df = self.master_dataframe.loc[
                    (self.master_dataframe.tf_id == tf_id)
                    & (self.master_dataframe.survey == survey)
                ]
                find_df = find_df.iloc[0]
            except IndexError:
                raise ValueError(
                    f"Could not find {survey}.{tf_id} in collection."
                )

        ref = find_df.hdf5_reference

        mt_object = MT()
        tf_object = self.mth5_collection.from_reference(ref)

        mt_object.__dict__.update(tf_object.__dict__)
        mt_object.station_metadata.update_time_period()
        mt_object.survey_metadata.update_time_period()

        return mt_object

    def _from_file(self, filename):
        """
        Add transfer functions for a list of file names

        :param file_list: DESCRIPTION
        :type file_list: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not self.mth5_collection.h5_is_write():
            raise ValueError("Must initiate an MTH5 file first.")
        if not isinstance(filename, (str, Path)):
            raise TypeError(
                f"filename must be a string or Path not {type(filename)}"
            )
        mt_object = MT(filename)
        mt_object.read()

        self._from_mt_object(mt_object)

    def _from_mt_object(self, mt_object):
        """

        :param mt_object: DESCRIPTION
        :type mt_object: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mt_object.survey_metadata.id in [None, ""]:
            mt_object.survey_metadata.id = "unknown_survey"
        self.mth5_collection.add_transfer_function(mt_object)
        self.logger.debug("added %s " % mt_object.station)

    def to_mt_data(self, bounding_box=None):
        """
        Get a list of transfer functions

        :param tf_ids: DESCRIPTION, defaults to None
        :type tf_ids: TYPE, optional
        :param bounding_box: DESCRIPTION, defaults to None
        :type bounding_box: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if bounding_box is not None:
            self.apply_bbox(*bounding_box)

        mt_data = MTData()

        for row in self.dataframe.itertuples():
            tf = self.get_tf(row.tf_id, survey=row.survey)

            mt_data.add_station(tf, compute_relative_location=False)

        # compute locations at the end
        mt_data.compute_relative_locations()

        return mt_data

    def from_mt_data(self, mt_data, new_survey=None, tf_id_extra=None):
        """
        Add data from a MTData object to an MTH5 collection.

        Can use 'new_survey' to create a new survey to load to.

        Can use 'tf_id_extra' to add a string onto the existing 'tf_id',
        useful if data have been edited or manipulated in some way.  For
        example could set 'tf_id_extra' = 'rotated' for rotated data. This will
        help you organize the tf's for each station.

        :param mt_data: MTData object
        :type mt_data: :class:`mtpy.core.mt_data.MTData`
        :param new_survey: new survey name, defaults to None
        :type new_survey: str, optional
        :param tf_id_extra: additional text onto existing 'tf_id',
         defaults to None
        :type tf_id_extra: string, optional
        :raises IOError: If an MTH5 is not writable raises

        """
        if self.mth5_collection.h5_is_write():
            if new_survey is not None or tf_id_extra is not None:
                for mt_obj in mt_data:
                    if new_survey is not None:
                        mt_obj.survey = new_survey
                    if tf_id_extra is not None:
                        mt_obj.tf_id = f"{mt_obj.tf_id}_{tf_id_extra}"
            self.add_tf(mt_data)

        else:
            raise IOError("MTH5 is not writeable, use 'open_mth5()'")

    def check_for_duplicates(self, locate="location", sig_figs=6):
        """
        Check for duplicate station locations in a MT DataFrame

        :param dataframe: DESCRIPTION
        :type dataframe: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.has_data():
            if locate == "location":
                self.dataframe.latitude = np.round(
                    self.dataframe.latitude, sig_figs
                )
                self.dataframe.longitude = np.round(
                    self.dataframe.longitude, sig_figs
                )

                query = ["latitude", "longitude"]
            elif locate not in self.dataframe.columns:
                raise ValueError(f"Not sure what to do with {locate}.")
            else:
                query = [locate]
            return self.dataframe[self.dataframe.duplicated(query)]
        return None

    def apply_bbox(self, lon_min, lon_max, lat_min, lat_max):
        """
        Return :class:`pandas.DataFrame` of station within bounding box

        :param longitude_min: Minimum longitude
        :type longitude_min: float
        :param longitude_max: Maximum longitude
        :type longitude_max: float
        :param latitude_min: Minimum latitude
        :type latitude_min: float
        :param latitude_max: Maximum longitude
        :type latitude_max: float
        :return: Only stations within the given bounding box
        :rtype: :class:`pandas.DataFrame`

        """

        if self.has_data():
            msg = (
                "Applying bounding box: "
                f"lon_min = {lon_min:.6g}, "
                f"lon_max = {lon_max:.6g}, "
                f"lat_min = {lat_min:.6g}, "
                f"lat_max = {lat_max:.6g}"
            )
            self.logger.debug(msg)

            self.working_dataframe = self.master_dataframe.loc[
                (self.master_dataframe.longitude >= lon_min)
                & (self.master_dataframe.longitude <= lon_max)
                & (self.master_dataframe.latitude >= lat_min)
                & (self.master_dataframe.latitude <= lat_max)
            ]

    def to_geo_df(self, bounding_box=None, epsg=4326):
        """
        Make a geopandas dataframe for easier GIS manipulation

        """
        coordinate_system = f"epsg:{epsg}"
        if bounding_box is not None:
            self.apply_bbox(*bounding_box)

        df = self.dataframe

        gdf = gpd.GeoDataFrame(
            df[
                df.columns[
                    ~df.columns.isin(
                        ["hdf5_reference", "station_hdf5_reference"]
                    )
                ]
            ],
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs=coordinate_system,
        )

        return gdf

    def to_shp(self, filename, bounding_box=None, epsg=4326):
        """

        :param filename: DESCRIPTION
        :type filename: TYPE
        :param bounding_box: DESCRIPTION, defaults to None
        :type bounding_box: TYPE, optional
        :param epsg: DESCRIPTION, defaults to 4326
        :type epsg: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.has_data():
            gdf = self.to_geo_df(bounding_box=bounding_box, epsg=epsg)
            gdf.to_file(self.working_directory.joinpath(filename))

            return gdf
        return None

    def average_stations(
        self,
        cell_size_m,
        bounding_box=None,
        count=1,
        n_periods=48,
        new_file=True,
    ):
        """
        Average nearby stations to make it easier to invert

        :param cell_size_m: DESCRIPTION
        :type cell_size_m: TYPE
        :param bounding_box: DESCRIPTION, defaults to None
        :type bounding_box: TYPE, optional
        :param save_dir: DESCRIPTION, defaults to None
        :type save_dir: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        # cell size in degrees (bit of a hack for now)
        r = cell_size_m / 111000.0

        if bounding_box:
            self.apply_bbox(*bounding_box)
            df = self.dataframe
        else:
            df = self.master_dataframe
        new_fn_list = []
        for ee in np.arange(
            df.longitude.min() - r / 2, df.longitude.max() + r, r
        ):
            for nn in np.arange(
                df.latitude.min() - r / 2, df.latitude.max() + r, r
            ):
                bbox = (ee, ee + r, nn, nn + r)
                self.apply_bbox(*bbox)
                if self.dataframe is None:
                    continue

                if len(self.dataframe) > 1:
                    m_list = [
                        self.get_tf(row.tf_id, row.survey)
                        for row in self.dataframe.itertuples()
                    ]
                    # interpolate onto a similar period range
                    f_list = []
                    for m in m_list:
                        f_list += m.period.tolist()
                    f = np.unique(np.array(f_list))
                    f = np.logspace(
                        np.log10(f.min()), np.log10(f.max()), n_periods
                    )

                    m_list_interp = []
                    for m in m_list:
                        m_list_interp.append(
                            m.interpolate(f, bounds_error=False)
                        )
                    avg_z = np.array(
                        [
                            m.impedance.data
                            for m in m_list_interp
                            if m.has_impedance()
                        ]
                    )

                    avg_z_err = np.array(
                        [
                            m.impedance_error.data
                            for m in m_list_interp
                            if m.has_impedance()
                        ]
                    )
                    avg_t = np.array(
                        [
                            m.tipper.data
                            for m in m_list_interp
                            if m.has_tipper()
                        ]
                    )
                    avg_t_err = np.array(
                        [
                            m.tipper_error.data
                            for m in m_list_interp
                            if m.has_tipper()
                        ]
                    )

                    avg_z[np.where(avg_z == 0 + 0j)] = np.nan + 1j * np.nan
                    avg_z_err[np.where(avg_z_err == 0)] = np.nan
                    avg_t[np.where(avg_t == 0 + 0j)] = np.nan + 1j * np.nan
                    avg_t_err[np.where(avg_t_err == 0 + 0j)] = np.nan

                    avg_z = np.nanmean(avg_z, axis=0)
                    avg_z_err = np.nanmean(avg_z_err, axis=0)
                    avg_t = np.nanmean(avg_t, axis=0)
                    avg_t_err = np.nanmean(avg_t_err, axis=0)

                    mt_avg = MT()
                    mt_avg.frequency = f
                    mt_avg.impedance = avg_z
                    mt_avg.impedance_error = avg_z_err

                    mt_avg.tipper = avg_t
                    mt_avg.tipper_error = avg_t_err

                    mt_avg.latitude = np.mean(
                        np.array([m.latitude for m in m_list])
                    )
                    mt_avg.longitude = np.mean(
                        np.array([m.longitude for m in m_list])
                    )
                    mt_avg.elevation = np.mean(
                        np.array([m.elevation for m in m_list])
                    )
                    mt_avg.station = f"AVG{count:03}"
                    mt_avg.station_metadata.comments = (
                        "avgeraged_stations = "
                        + ",".join([m.station for m in m_list])
                    )
                    mt_avg.survey_metadata.id = "averaged"
                    self.add_tf(mt_avg)

                    try:

                        if new_file:
                            edi_obj = mt_avg.write(
                                save_dir=self.working_directory
                            )
                            self.logger.info(
                                f"wrote average file {edi_obj.fn}"
                            )
                        new_fn_list.append(edi_obj.fn)
                        count += 1
                    except Exception as error:
                        self.logger.exception(
                            "Failed to average files %s", error
                        )
                else:
                    continue

    def plot_mt_response(self, tf_id, survey=None, **kwargs):
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
        mt_data = MTData()
        if isinstance(tf_id, str):
            mt_object = self.get_tf(tf_id, survey=survey)
            return mt_object.plot_mt_response(**kwargs)
        elif isinstance(tf_id, (list, tuple, np.ndarray, pd.Series)):
            tf_request = np.array(tf_id)
            if len(tf_request.shape) > 1:
                for row in tf_request:
                    mt_data.add_station(self.get_tf(row[0], survey=row[1]))

            else:
                for row in tf_request:
                    mt_data.add_station(self.get_tf(row, survey=survey))
            return PlotMultipleResponses(mt_data, **kwargs)

        elif isinstance(tf_id, pd.DataFrame):
            for row in tf_id.itertuples():
                mt_data.add_station(self.get_tf(row.tf_id, survey=row.survey))
            return PlotMultipleResponses(mt_data, **kwargs)

    def plot_stations(self, map_epsg=4326, bounding_box=None, **kwargs):
        """
        plot stations

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        if self.dataframe is not None:
            gdf = self.to_geo_df(epsg=map_epsg, bounding_box=bounding_box)
            return PlotStations(gdf, **kwargs)

    def plot_strike(self, mt_data=None, **kwargs):
        """
        Plot strike angle

        .. seealso:: :class:`mtpy.imaging.PlotStrike`
        """
        if mt_data is None:
            mt_data = self.to_mt_data()
        return PlotStrike(mt_data, **kwargs)

    def plot_phase_tensor(self, tf_id, survey=None, **kwargs):
        """
        plot phase tensor elements

        :param tf_id: DESCRIPTION
        :type tf_id: TYPE
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        tf_obj = self.get_tf(tf_id, survey=survey)
        return tf_obj.plot_phase_tensor(**kwargs)

    def plot_phase_tensor_map(self, mt_data=None, **kwargs):
        """
        Plot Phase tensor maps for transfer functions in the working_dataframe

        .. seealso:: :class:`mtpy.imaging.PlotPhaseTensorMaps`

        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mt_data is None:
            mt_data = self.to_mt_data()

        return PlotPhaseTensorMaps(mt_data=mt_data, **kwargs)

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

        if mt_data is None:
            mt_data = self.to_mt_data()

        return PlotPhaseTensorPseudoSection(mt_data=mt_data, **kwargs)

    def plot_residual_phase_tensor(
        self, mt_data_01, mt_data_02, plot_type="map", **kwargs
    ):
        """

        :param mt_data_01: DESCRIPTION
        :type mt_data_01: TYPE
        :param mt_data_02: DESCRIPTION
        :type mt_data_02: TYPE
        :param plot_type: DESCRIPTION, defaults to "map"
        :type plot_type: TYPE, optional
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if plot_type in ["map"]:
            return PlotResidualPTMaps(mt_data_01, mt_data_02, **kwargs)

        if plot_type in ["pseudosection", "ps"]:
            return PlotResidualPTPseudoSection(
                mt_data_01, mt_data_02, **kwargs
            )

    def plot_penetration_depth_1d(self, tf_id, survey=None, **kwargs):
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

        tf_object = self.get_tf(tf_id, survey=survey)

        return PlotPenetrationDepth1D(tf_object, **kwargs)

    def plot_penetration_depth_map(self, mt_data=None, **kwargs):
        """
        Plot Penetration depth in map view for a single period

        .. seealso:: :class:`mtpy.imaging.PlotPenetrationDepthMap`

        :param mt_data: DESCRIPTION
        :type mt_data: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mt_data is None:
            mt_data = self.to_mt_data()

        return PlotPenetrationDepthMap(mt_data, **kwargs)

    def plot_resistivity_phase_maps(self, mt_data=None, **kwargs):
        """
        Plot apparent resistivity and/or impedance phase maps from the
        working dataframe

        .. seealso:: :class:`mtpy.imaging.PlotResPhaseMaps`
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mt_data is None:
            mt_data = self.to_mt_data()

        return PlotResPhaseMaps(mt_data, **kwargs)

    def plot_resistivity_phase_pseudosections(self, mt_data=None, **kwargs):
        """
        Plot resistivity and phase in a pseudosection along a profile line

        :param mt_data: DESCRIPTION, defaults to None
        :type mt_data: TYPE, optional
        :param **kwargs: DESCRIPTION
        :type **kwargs: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if mt_data is None:
            mt_data = self.to_mt_data()

        return PlotResPhasePseudoSection(mt_data, **kwargs)
