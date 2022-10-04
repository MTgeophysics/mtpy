# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:04:12 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from pyproj import CRS

from mtpy.utils.mtpy_logger import get_mtpy_logger
from mtpy.utils.gis_tools import (
    assert_lat_value,
    assert_lon_value,
    assert_elevation_value,
    project_point,
)

# =============================================================================


class MTLocation:
    """
    Location for a MT site or point measurement

    """

    def __init__(self, **kwargs):

        self.logger = get_mtpy_logger(f"{__name__}.{self.__class__.__name__}")

        self._latitude = 0
        self._longitude = 0
        self._elevation = 0
        self._east = 0
        self._north = 0
        self._datum_epsg = CRS.from_epsg(4326)
        self._utm_epsg = None
        self._geoid_epsg = None
        self.model_east = 0
        self.model_north = 0
        self.model_elevation = 0

    def __str__(self):
        lines = ["MT Location: ", "-" * 20]
        lines.append(f"  Latitude (deg):   {self.latitude:.6f}")
        lines.append(f"  Longitude (deg):  {self.longitude:.6f}")
        lines.append(f"  Elevation (m):    {self.elevation:.4f}")
        lines.append(f"  Datum EPSG:       {self.datum_epsg.to_epsg()}")
        lines.append(f"  Easting (m):      {self.east:.3f}")
        lines.append(f"  Northing (m):     {self.north:.3f}")
        lines.append(f"  UTM EPSG:         {self.utm_epsg}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    @property
    def datum_epsg(self):
        return self._datum_epsg

    @datum_epsg.setter
    def datum_epsg(self, value):
        new_crs = CRS.from_user_input(value)

        if new_crs != self._datum_epsg:
            if (
                self._datum_epsg is not None
                and self.latitude != 0
                and self.longitude != 0
            ):
                self._longitude, self._latitude = project_point(
                    self.longitude, self.latitude, self._datum_epsg, new_crs
                )

                self._east, self._north = project_point(
                    self.longitude, self.latitude, new_crs, self.utm_epsg
                )

            elif (
                self.datum_epsg is not None
                and self.east != 0
                and self.north != 0
                and self.latitude == 0
                and self.longitude == 0
            ):
                self._longitude, self._latitude = project_point(
                    self.east,
                    self.north,
                    self.utm_epsg,
                    new_crs,
                )
            self._datum_epsg = new_crs

    @property
    def utm_epsg(self):
        return self._utm_epsg

    @utm_epsg.setter
    def utm_epsg(self, value):
        if value != self._utm_epsg:
            # reproject easting, northing to new zone
            if (
                self._utm_epsg is not None
                and self.east != 0
                and self.north != 0
            ):
                self._east, self._north = project_point(
                    self.east, self.north, self._utm_epsg, value
                )

            if (
                self.datum_epsg is not None
                and self.east != 0
                and self.north != 0
            ):
                # reproject lat and lon base on new UTM datum
                self._latitude, self._longitude = project_point(
                    self.east,
                    self.north,
                    value,
                    self.datum_epsg,
                )

            # if east and north == 0 and lat and lon != 0 project to utm
            elif (
                self.datum_epsg is not None
                and self.east == 0
                and self.north == 0
                and self.latitude != 0
                and self.longitude != 0
            ):
                self._east, self._north = project_point(
                    self.longitude,
                    self.latitude,
                    self.datum_epsg,
                    value,
                )

            self._utm_epsg = value

    @property
    def east(self):
        """easting"""
        return self._east

    @east.setter
    def east(self, value):
        """set east"""
        self._east = value
        if self.datum_epsg is not None and self.utm_epsg is not None:
            self._longitude, self._latitude = project_point(
                self._east, self._north, self.utm_epsg, self.datum_epsg
            )

    @property
    def north(self):
        """northing"""
        return self._north

    @north.setter
    def north(self, value):
        """set north"""
        self._north = value
        if self.datum_epsg is not None and self.utm_epsg is not None:
            self._longitude, self._latitude = project_point(
                self._east, self._north, self.utm_epsg, self.datum_epsg
            )

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, lat):
        self._latitude = assert_lat_value(lat)
        if self.utm_epsg is not None and self.datum_epsg is not None:
            self._east, self._north = project_point(
                self._longitude, self._latitude, self.datum_epsg, self.utm_epsg
            )

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, lon):
        self._longitude = assert_lon_value(lon)
        if self.utm_epsg is not None and self.datum_epsg is not None:
            self._east, self._north = project_point(
                self._longitude, self._latitude, self.datum_epsg, self.utm_epsg
            )

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, elev):
        self._elevation = assert_elevation_value(elev)

    @property
    def model_east(self):
        return self._model_east

    @model_east.setter
    def model_east(self, value):
        try:
            self._model_east = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Input should be a float not type {type(value)}")

    @property
    def model_north(self):
        return self._model_north

    @model_north.setter
    def model_north(self, value):
        try:
            self._model_north = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Input should be a float not type {type(value)}")

    @property
    def model_elevation(self):
        return self._model_elevation

    @model_elevation.setter
    def model_elevation(self, value):
        try:
            self._model_elevation = float(value)
        except (TypeError, ValueError):
            raise ValueError(f"Input should be a float not type {type(value)}")
