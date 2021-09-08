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

import numpy as np
import pandas as pd

import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
from mtpy.core import mt
from mtpy.utils.mtpy_logger import get_mtpy_logger

# =============================================================================
#
# =============================================================================


class MTCollection:
    """
    Collection of mt data
    """

    def __init__(self, mt_path=None):
        self._mt_path = self.check_path(mt_path)
        self.mt_df = None
        self.logger = get_mtpy_logger(
            f"{__name__}.{self.__class__.__name__}", fn="mt_collection"
        )

    def __str__(self):
        lines = [f"MT file path: {self.mt_path}"]
        if self.mt_df is not None:
            lines += [f"\t Number of Stations: {len(self.mt_df)}"]
        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def mt_path(self):
        return self._mt_path

    @mt_path.setter
    def mt_path(self, value):
        self._mt_path = self.check_path(value)

    def check_path(self, mt_path):
        if mt_path is None:
            return None
        else:
            mt_path = Path(mt_path)
            if not mt_path.exists():
                msg = f"{mt_path} does not exists"
                self.logger.error(msg)
                raise IOError(msg)
            return mt_path

    def make_mt_file_list(self, mt_path=None, file_types=["edi"]):
        """
        Get a list of MT file from a given path

        :param mt_path: full path to where the MT transfer functions are stored
        :type mt_path: string or :class:`pathlib.Path`

        :param file_types: List of file types to look for given their extension
        :type file_types: list

        Currently available file types are or will be:
            - edi - EDI files
            - zmm - EMTF output file
            - j - BIRRP output file
            - avg - Zonge output file

        """
        if mt_path is not None:
            self.mt_path = mt_path
        if self.mt_path is None:
            self.logger.info("MT file path is None, returning None.")
            return None

        fn_list = []
        for ext in file_types:
            fn_list += list(self.mt_path.glob(f"*.{ext}"))

        return fn_list

    def make_dataframe_from_file_list(self, mt_file_list=None, move_duplicates=True):
        """
        create a :class:`pandas.DataFrame` from information in the file list

        :param mt_file_list: list of paths to MT tranfer function files
        :type mt_file_list: list of :class:`pathlib.Path` objects
        :return: Data frame with information about each transfer function
        :rtype: :class:`pandas.DataFrame`


        """

        if mt_file_list is None:
            mt_file_list = self.make_mt_file_list()

        station_list = []
        for fn in mt_file_list:
            m = mt.MT(fn)
            entry = {}
            entry["ID"] = m.station
            entry["start"] = m.station_metadata.time_period.start
            entry["end"] = m.station_metadata.time_period.end
            entry["latitude"] = m.latitude
            entry["longitude"] = m.longitude
            entry["elevation"] = m.elevation
            entry["easting"] = m.east
            entry["northing"] = m.north
            entry["utm_zone"] = m.utm_zone
            entry["acquired_by"] = m.station_metadata.acquired_by.author
            entry["period_min"] = 1.0 / m.Z.freq.max()
            entry["period_max"] = 1.0 / m.Z.freq.min()
            entry["n_periods"] = m.Z.freq.size
            entry["survey"] = m.survey_metadata.id
            entry["fn"] = fn
            entry["file_date"] = m.station_metadata.provenance.creation_time

            # add entry to list to put into data frame
            station_list.append(entry)

        mt_df = pd.DataFrame(station_list)
        if move_duplicates:
            mt_df = self._check_for_duplicates(mt_df)

        return mt_df

    def add_stations_from_file_list(self, fn_list, remove_duplicates=True):
        """
        Add stations from a list of files
        
        :param fn_list: DESCRIPTION
        :type fn_list: TYPE
        :param remove_duplicates: DESCRIPTION, defaults to True
        :type remove_duplicates: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        new_df = self.make_dataframe_from_file_list(fn_list)
        self.mt_df.append(new_df, ignore_index=True)
        return self._check_for_duplicates(self.mt_df.append(new_df, ignore_index=True))

    def _check_for_duplicates(self, mt_df, locate="location"):
        """
        Check for duplicate station locations in a MT DataFrame

        :param mt_df: DESCRIPTION
        :type mt_df: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        # write csv file to querry
        duplicates = mt_df[mt_df.duplicated(["latitude", "longitude"])]
        if len(duplicates) > 0:
            self.logger.info(
                f"Found {len(mt_df)} duplicates, moving oldest to 'Duplicates'"
            )
            dup_path = self.mt_path.joinpath("Duplicates")
            if not dup_path.exists():
                dup_path.mkdir()
            for ii, row in duplicates.iterrows():
                fn = Path(row.fn)
                new_fn = dup_path.joinpath(Path(row.fn).name)
                try:
                    self.logger.info("Moved %s to Duplicates", new_fn.name)
                    fn.rename(new_fn)
                except FileNotFoundError:
                    self.logger.debug(f"Could not find {fn} --> skipping")
        mt_df = mt_df.drop_duplicates(subset=["latitude", "longitude"], keep="first")

        return mt_df

    def from_csv(self, csv_fn):
        """
        Read in a csv of 
        :param csv_fn: DESCRIPTION
        :type csv_fn: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        self.mt_df = pd.read_csv(csv_fn)
        self.mt_path = Path(csv_fn).parent

    def to_csv(self, filename):
        """
        Write dataframe to a CSV file

        :param filename: DESCRIPTION
        :type filename: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        filename = Path(filename)
        if filename.parent.as_posix() == ".":
            if self.mt_path is not None:
                filename = self.mt_path.joinpath(filename)
        self.mt_df.to_csv(filename, index=False)
        self.logger.info(f"Wrote CSV file to {filename}")

    def apply_bbox(self, x_min, x_max, y_min, y_max, units="degrees", utm_zone=None):
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
        if units in ["degrees"]:
            msg = (
                "Applying bounding box: "
                f"lon_min = {x_min:.6g}, "
                f"lon_max = {x_max:.6g}, "
                f"lat_min = {y_min:.6g}, "
                f"lat_max = {y_max:.6g}"
            )
            self.logger.debug(msg)
    
            return self.mt_df.loc[
                (self.mt_df.longitude >= x_min)
                & (self.mt_df.longitude <= x_max)
                & (self.mt_df.latitude >= y_min)
                & (self.mt_df.latitude <= y_max)
            ]
        
        if units in ["m"]:
            if utm_zone is None:
                raise ValueError("UTM Zone must be input")
                
            msg = (
                "Applying bounding box: "
                f"east_min = {x_min:.6g}, "
                f"east_max = {x_max:.6g}, "
                f"north_min = {y_min:.6g}, "
                f"north_max = {y_max:.6g}"
                f"utm_zone = {utm_zone}"
            )
            self.logger.debug(msg)
    
            return self.mt_df[self.mt_df.utm_zone == utm_zone].loc[
                (self.mt_df.easting >= x_min)
                & (self.mt_df.easting <= x_max)
                & (self.mt_df.northing >= y_min)
                & (self.mt_df.northing <= y_max)
            ]
        
    def write_shp_file(self, filename, bounding_box=None, epsg=4326):
        """
        
        :param filename: DESCRIPTION
        :type filename: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """
        coordinate_system = {"init": f"epsg:{epsg}"}
        # write shape file
        geometry_list = []
        if self.mt_df is None:
            self.mt_df = self.make_dataframe_from_file_list()
        for ii, row in self.mt_df.iterrows():
            geometry_list.append(Point(row.longitude, row.latitude))
        
        gdf = gpd.GeoDataFrame(self.mt_df, crs=coordinate_system, geometry=geometry_list)
        gdf.fn = gdf.fn.astype("str")
        gdf.to_file(self.mt_path.joinpath(filename))
        
    def average_stations(self, cell_size_m, bounding_box=None, save_dir=None,
                         units="deg", utm_zone=None, count=1):
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
        if save_dir is None:
            save_dir = self.mt_path
        r = cell_size_m
        
        if bounding_box:
            df = self.apply_bbox(*bounding_box, units=units, utm_zone=utm_zone)
        
        else:
            df = self.mt_df
            
        new_fn_list = []
        for ee in np.arange(df.easting.min(), df.easting.max(), r):
            for nn in np.arange(df.northing.min(), df.northing.max(), r):
                bbox = (ee, ee + r, nn, nn + r)
                avg_df = self.apply_bbox(*bbox, units="m", utm_zone=utm_zone)
                if len(avg_df) > 1:
                    m_list = [mt.MT(row.fn) for row in avg_df.itertuples()]
                    # interpolate onto a similar period range
                    print(m_list)
                    f_list = []
                    for m in m_list:
                        f_list += m.Z.freq.tolist()
                    f = np.unique(np.array(f_list))
                    print(f)
                    f = np.logspace(np.log10(f.min()), np.log10(f.max()), 48)
                    for m in m_list:
                        m.Z, m.Tipper = m.interpolate(f, bounds_error=False)
                    avg_z = np.array([m.Z.z for m in m_list])
                    avg_z_err = np.array([m.Z.z_err for m in m_list])
                    avg_t = np.array([m.Tipper.tipper for m in m_list])
                    avg_t_err = np.array([m.Tipper.tipper_err for m in m_list])
                    
                    avg_z[np.where(avg_z==0 + 0j)] = np.nan + 1j * np.nan
                    avg_z_err[np.where(avg_z_err==0)] = np.nan
                    avg_t[np.where(avg_t==0 + 0j)] = np.nan + 1j * np.nan
                    avg_t_err[np.where(avg_t_err==0 + 0j)] = np.nan
                    
                    avg_z = np.nanmean(avg_z, axis=0)
                    avg_z_err = np.nanmean(avg_z_err, axis=0)
                    avg_t = np.nanmean(avg_t, axis=0)
                    avg_t_err = np.nanmean(avg_t_err, axis=0)
                    
                    mt_avg = mt.MT()
                    mt_avg.Z.freq = f
                    mt_avg.Z.z = avg_z
                    mt_avg.Z.z_err = avg_z_err
                    
                    mt_avg.Tipper.freq = f
                    mt_avg.Tipper.tipper = avg_t
                    mt_avg.Tipper.tipper_err = avg_t_err
                    
                    mt_avg.latitude = np.mean(np.array([m.latitude for m in m_list]))
                    mt_avg.longitude = np.mean(np.array([m.longitude for m in m_list]))
                    mt_avg.elevation = np.mean(np.array([m.elevation for m in m_list]))
                    mt_avg.station = f"AVG{count:03}"
                    mt_avg.station_metadata.comments = (
                        "avgeraged_stations = " + ",".join([m.station for m in m_list])
                        )
                    
                    try:
                        edi_obj = mt_avg.write_mt_file(save_dir=save_dir)
                        print(f"wrote average file {edi_obj.fn}")
                        new_fn_list.append(edi_obj.fn)
                        count += 1 
            
  
                    except Exception as error:
                        print(f"{error} ")
        
                else:
                    continue
                
        new_df = self.make_dataframe_from_file_list(new_fn_list)
        
        return new_df
