"""
==================
ModEM
==================

# Generate files for ModEM

# revised by JP 2017
# revised by AK 2017 to bring across functionality from ak branch

"""
from pathlib import Path
import numpy as np
from loguru import logger

import geopandas as gpd
from shapely.geometry import Point

from mtpy.core import mt as mt
from mtpy.utils import gis_tools as gis_tools


# in module imports
from .exception import ModEMError

__all__ = ["Stations"]


class Stations(object):
    """
    station locations class

    """

    def __init__(self, **kwargs):

        self.logger = logger

        self.dtype = [
            ("station", "|U50"),
            ("lat", float),
            ("lon", float),
            ("elev", float),
            ("rel_east", float),
            ("rel_north", float),
            ("rel_elev", float),
            ("east", float),
            ("north", float),
            ("zone", "U4"),
        ]
        self.station_locations = np.zeros(0, dtype=self.dtype)
        self._model_epsg = None
        self._model_utm_zone = None
        self._center_lat = None
        self._center_lon = None
        self._center_elev = 0.0

        for key in list(kwargs.keys()):
            if hasattr(self, key):
                setattr(self, key, kwargs[key])

    def __str__(self):
        fmt_dict = dict(
            [
                ("station", "<8"),
                ("lat", "<10.4f"),
                ("lon", "<10.4f"),
                ("elev", "<8.2f"),
                ("rel_east", "<13.2f"),
                ("rel_north", "<13.2f"),
                ("rel_elev", "<8.2f"),
                ("east", "<12.2f"),
                ("north", "<12.2f"),
                ("zone", "<6"),
            ]
        )
        lines = [
            "".join(
                [
                    f"{n.capitalize():<10}"
                    for n in self.station_locations.dtype.names
                ]
            )
        ]
        lines.append("-" * 72)
        for ss in self.station_locations:
            l = []
            for k in self.station_locations.dtype.names:
                l.append(f"{ss[k]:{fmt_dict[k]}}")
            lines.append("".join(l))

        lines.append("\nModel Center:")
        l = []
        for n in ["lat", "lon", "elev", "east", "north", "zone"]:
            l.append(f"{self.center_point[n][0]:{fmt_dict[n]}}")
        lines.append("".join(l))

        lines.append("\nMean Values:")
        l = []
        for n in ["lat", "lon", "elev", "east", "north"]:
            l.append(f"{self.station_locations[n].mean():{fmt_dict[n]}}")
        lines.append("".join(l) + f"{self.center_point.zone[0]:<6}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    ## --> define properties that can only be returned and not set
    @property
    def lat(self):
        return self.station_locations["lat"]

    @property
    def lon(self):
        return self.station_locations["lon"]

    @property
    def east(self):
        return self.station_locations["east"]

    @property
    def north(self):
        return self.station_locations["north"]

    @property
    def elev(self):
        return self.station_locations["elev"]

    @property
    def rel_east(self):
        return self.station_locations["rel_east"]

    @property
    def rel_north(self):
        return self.station_locations["rel_north"]

    @property
    def rel_elev(self):
        return self.station_locations["rel_elev"]

    @property
    def utm_zone(self):
        return self.station_locations["zone"]

    @property
    def station(self):
        return self.station_locations["station"]

    @property
    def model_epsg(self):
        return self._model_epsg

    @model_epsg.setter
    def model_epsg(self, value):
        """
        set the model epsg number an project east, north
        """
        self._model_epsg = value
        if self.station_locations.size < 2:
            for ss, ii in enumerate(self.station_locations):
                east, north, utm_zone = gis_tools.project_point_ll2utm(
                    ss["lat"],
                    ss["lon"],
                    epsg=self._model_epsg,
                )
                self.station_locations[ii]["east"] = east
                self.station_locations[ii]["north"] = north
                self.station_locations[ii]["zone"] = utm_zone

    @property
    def model_utm_zone(self):
        return self._model_utm_zone

    @model_utm_zone.setter
    def model_utm_zone(self, value):
        """
        set the model epsg number an project east, north
        """
        if value is None:
            return

        self.logger.debug(f"Setting model utm zone to {value}")

        self._model_utm_zone = value
        if self.station_locations.size > 1:
            for ii, ss in enumerate(self.station_locations):
                east, north, utm_zone = gis_tools.project_point_ll2utm(
                    ss["lat"],
                    ss["lon"],
                    utm_zone=self._model_utm_zone,
                )
                self.station_locations[ii]["east"] = east
                self.station_locations[ii]["north"] = north
                self.station_locations[ii]["zone"] = utm_zone

    def _get_mt_objs_from_list(self, input_list):
        """
        get mt_objects from a list of files or mt_objects
        """
        if isinstance(input_list, (list, np.ndarray)):
            if isinstance(input_list[0], mt.MT):
                return input_list

            elif isinstance(input_list[0], (str, Path)):
                return [mt.MT(fn) for fn in input_list]
        else:
            raise ValueError(f"type {type(input_list)} is not supported yet")

    def get_station_locations(self, input_list):
        """
        get station locations from a list of edi files

        Arguments
        -------------
            **input_list** : list
                             list of edi file names, or mt_objects


        Returns
        ------------
            * fills station_locations array

        """
        # self.logger.debug input_list
        mt_obj_list = self._get_mt_objs_from_list(input_list)

        # if station locations are not input read from the edi files
        if mt_obj_list is None:
            raise AttributeError(
                "mt_obj_list is None, need to input a list of "
                "mt objects to read in."
            )

        n_stations = len(mt_obj_list)

        if n_stations == 0:
            raise ModEMError(
                "No .edi files in edi_list, please check " "file locations."
            )

        # make a structured array to put station location information into
        self.station_locations = np.zeros(n_stations, dtype=self.dtype)
        # get station locations in meters
        for ii, mt_obj in enumerate(mt_obj_list):
            self.station_locations[ii]["lat"] = mt_obj.latitude
            self.station_locations[ii]["lon"] = mt_obj.longitude
            self.station_locations[ii]["station"] = mt_obj.station
            self.station_locations[ii]["elev"] = mt_obj.elevation

            if (self.model_epsg is not None) or (
                self.model_utm_zone is not None
            ):
                east, north, utm_zone = gis_tools.project_point_ll2utm(
                    mt_obj.latitude,
                    mt_obj.longitude,
                    utm_zone=self.model_utm_zone,
                    epsg=self.model_epsg,
                )
                self.station_locations[ii]["east"] = east
                self.station_locations[ii]["north"] = north
                self.station_locations[ii]["zone"] = utm_zone
            else:
                self.station_locations[ii]["east"] = mt_obj.east
                self.station_locations[ii]["north"] = mt_obj.north
                self.station_locations[ii]["zone"] = mt_obj.utm_zone

        # get relative station locations
        self.calculate_rel_locations()

    def calculate_rel_locations(self, shift_east=0, shift_north=0):
        """
        put station in a coordinate system relative to
        (shift_east, shift_north)
        (+) shift right or up
        (-) shift left or down

        """

        # translate the stations so they are relative to 0,0
        # east_center = (self.east.max() + self.east.min()) / 2.0
        # north_center = (self.north.max() + self.north.min()) / 2.0

        # self.station_locations["rel_east"] = self.east - east_center
        # self.station_locations["rel_north"] = self.north - north_center
        self.station_locations["rel_east"] = (
            self.east - self.center_point.east[0]
        )
        self.station_locations["rel_north"] = (
            self.north - self.center_point.north[0]
        )

        # BM: Before topograhy is applied to the model, the station
        #  elevation isn't relative to anything (according to
        #  Data.project_stations_on_topography, station elevation is
        #  relevant to topography). So rel_elev and elev are the same.
        #  Once topography has been applied, rel_elev can be calcuated
        #  by calling Data.project_stations_on_topography.
        self.station_locations["rel_elev"] = self.elev

    # make center point a get method, can't set it.
    @property
    def center_point(self):
        """
        calculate the center point from the given station locations

        Returns
        -------------
            **center_location** : np.ndarray
                                  structured array of length 1
                                  dtype includes (east, north, zone, lat, lon)
        """
        dtype = [
            ("lat", float),
            ("lon", float),
            ("east", float),
            ("north", float),
            ("elev", float),
            ("zone", "U4"),
        ]
        center_location = np.recarray(1, dtype=dtype)
        if self._center_lat is not None and self._center_lon is not None:
            self.logger.debug("assigning center from user set values")
            center_location.lat[0] = self._center_lat
            center_location.lon[0] = self._center_lon
            center_location.elev[0] = self._center_elev

            # get the median utm zone
            if self.model_utm_zone is None:
                zone = self.utm_zone.copy()
                zone.sort()
                # get the median zone
                center_utm_zone = zone[int(zone.size / 2)]
                center_location.zone[0] = center_utm_zone
            else:
                center_location.zone[0] = self.model_utm_zone

            # project center
            east, north, zone = gis_tools.project_point_ll2utm(
                center_location.lat[0],
                center_location.lon[0],
                utm_zone=center_location.zone[0],
            )

            center_location.east[0] = east
            center_location.north[0] = north
            return center_location

        # safer to get center from lat and lon if not all zones are the same
        if not np.all(self.utm_zone == self.utm_zone[0]):
            self.logger.debug("Not all stations are in same UTM zone")
            center_location.lat[0] = (self.lat.max() + self.lat.min()) / 2.0
            center_location.lon[0] = (self.lon.max() + self.lon.min()) / 2.0
            # get the median utm zone
            if self.model_utm_zone is None:
                self.logger.info(
                    "Getting median UTM zone of stations for center point"
                )
                zone = self.utm_zone.copy()
                zone.sort()
                center_utm_zone = zone[int(zone.size / 2)]
                center_location.zone[0] = center_utm_zone
            else:
                self.logger.info(
                    f"Using user defined center point UTM zone {self.model_utm_zone}"
                )
                center_location.zone[0] = self.model_utm_zone

            self.logger.info(
                f"Projecting lat, lon to UTM zone {center_location.zone[0]}"
            )
            east, north, zone = gis_tools.project_point_ll2utm(
                center_location.lat[0],
                center_location.lon[0],
                utm_zone=center_location.zone[0],
            )

            center_location.east[0] = east
            center_location.north[0] = north

        else:
            self.logger.debug("locating center from UTM grid")
            center_location.east[0] = (self.east.max() + self.east.min()) / 2
            center_location.north[0] = (
                self.north.max() + self.north.min()
            ) / 2

            # get the median utm zone
            zone = self.utm_zone.copy()
            zone.sort()
            center_utm_zone = zone[int(zone.size / 2)]
            center_location.zone[0] = center_utm_zone

            center_ll = gis_tools.project_point_utm2ll(
                center_location.east[0],
                center_location.north[0],
                center_location.zone[0],
                epsg=self.model_epsg,
            )

            center_location.lat[0] = center_ll[0]
            center_location.lon[0] = center_ll[1]
        # BM: Because we are now writing center_point.elev to ModEm
        #  data file, we need to provide it.
        #  The center point elevation is the highest point of the
        #  model. Before topography is applied, this is the highest
        #  station. After it's applied, it's the highest point
        #  point of the surface model (this will be set by calling
        #  Data.project_stations_on_topography).
        if self._center_elev:
            center_location.elev[0] = self._center_elev
        else:
            center_location.elev[0] = -self.elev.max()

        return center_location

    def rotate_stations(self, rotation_angle):
        """
        Rotate stations assuming N is 0

        Arguments
        -------------
            **rotation_angle** : float
                                 angle in degrees assuming N is 0

        Returns
        -------------
            * refils rel_east and rel_north in station_locations.  Does this
              because you will still need the original locations for plotting
              later.

        """

        cos_ang = np.cos(np.deg2rad(rotation_angle))
        sin_ang = np.sin(np.deg2rad(rotation_angle))
        rot_matrix = np.array([[cos_ang, sin_ang], [-sin_ang, cos_ang]])

        coords = np.array(
            [
                self.station_locations["rel_east"],
                self.station_locations["rel_north"],
            ]
        )

        # rotate the relative station locations
        new_coords = np.array(np.dot(rot_matrix, coords))

        self.station_locations["rel_east"] = new_coords[0, :]
        self.station_locations["rel_north"] = new_coords[1, :]

        self.logger.info(
            f"Rotated stations by {rotation_angle:.1f} deg clockwise from N"
        )

    def to_geopd(self, epsg=None, default_epsg=4326):
        """
        create a geopandas dataframe

        :param epsg: EPSG number to project to
        :type epsg: integer, defaults to None
        :param default_epsg: the default EPSG number that the stations are
        referenced to
        :type default_epsg: integer, defaults to 4326

        """

        default_crs = {"init": f"epsg:{default_epsg}"}
        station_list = []
        geometry_list = []
        for sarr in self.station_locations:
            entry = {
                "station": sarr["station"],
                "latitude": sarr["lat"],
                "longitude": sarr["lon"],
                "elevation": sarr["elev"],
                "easting": sarr["east"],
                "northing": sarr["north"],
                "utm_zone": sarr["zone"],
                "model_east": sarr["rel_east"],
                "model_north": sarr["rel_north"],
                "model_elev": sarr["rel_elev"],
            }
            geometry_list.append(Point(sarr["lon"], sarr["lat"]))
            station_list.append(entry)
        sdf = gpd.GeoDataFrame(
            station_list, crs=default_crs, geometry=geometry_list
        )
        if epsg is not None:
            sdf = sdf.to_crs(epsg=epsg)

        return sdf

    def to_shp(self, shp_fn, epsg=None, default_epsg=4326):
        """
        Write a shape file of the station locations using geopandas which only takes
        in epsg numbers

        :param shp_fn: full path to new shapefile
        :type shp_fn: string
        :param epsg: EPSG number to project to
        :type epsg: integer, defaults to None
        :param default_epsg: the default EPSG number that the stations are
        referenced to
        :type default_epsg: integer, defaults to 4326

        """
        sdf = self.to_geopd(epsg=epsg, default_epsg=default_epsg)

        sdf.to_file(shp_fn)

    def to_csv(self, csv_fn, epsg=None, default_epsg=4326, geometry=False):
        """
        Write a shape file of the station locations using geopandas which only takes
        in epsg numbers

        :param shp_fn: full path to new shapefile
        :type shp_fn: string
        :param epsg: EPSG number to project to
        :type epsg: integer, defaults to None
        :param default_epsg: the default EPSG number that the stations are
        referenced to
        :type default_epsg: integer, defaults to 4326

        """
        sdf = self.to_geopd(epsg=epsg, default_epsg=default_epsg)
        use_columns = list(sdf.columns)
        if not geometry:
            use_columns.remove("geometry")
        sdf.to_csv(csv_fn, index=False, columns=use_columns)
