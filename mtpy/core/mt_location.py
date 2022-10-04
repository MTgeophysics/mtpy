# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 15:04:12 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from pyproj import Transformer

from mtpy.utils.mtpy_logger import get_mtpy_logger
from mtpy.utils.gis_tools import (
    assert_lat_value,
    assert_lon_value,
    assert_elevation_value,
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
        self._datum_epsg = 4326
        self._utm_epsg = None
        self._geoid_epsg = None

    def __str__(self):
        lines = ["MT Location: ", "-" * 20]
        lines.append(f"  Latitude (deg):   {self.latitude:.6f}")
        lines.append(f"  Longitude (deg):  {self.longitude:.6f}")
        lines.append(f"  Elevation (m):    {self.elevation:.4f}")
        lines.append(f"  Datum EPSG:       {self.datum_epsg}")
        lines.append(f"  Easting (m):      {self.east:.3f}")
        lines.append(f"  Northing (m):     {self.north:.3f}")
        lines.append(f"  UTM EPSG:         {self.utm_epsg}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def project_point(self, x, y, old_epsg, new_epsg):
        """
        Transform point to new epsg

        :param x: DESCRIPTION
        :type x: TYPE
        :param y: DESCRIPTION
        :type y: TYPE
        :param old_epsg: DESCRIPTION
        :type old_epsg: TYPE
        :param new_epsg: DESCRIPTION
        :type new_epsg: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if old_epsg is None:
            raise ValueError("Original EPSG must not be None")
        if new_epsg is None:
            raise ValueError("New EPSG must not be None")
        if x == 0 or y == 0:
            raise ValueError("Should not project with 0 value")

        transformer = Transformer.from_crs(
            f"epsg:{old_epsg}", f"epsg:{new_epsg}", always_xy=True
        )

        return transformer.transform(x, y)

    @property
    def datum_epsg(self):
        return self._datum_epsg

    @datum_epsg.setter
    def datum_epsg(self, value):
        if value != self._datum_epsg:
            if (
                self._datum_epsg is not None
                and self.latitude != 0
                and self.longitude != 0
            ):
                self._longitude, self._latitude = self.project_point(
                    self.longitude, self.latitude, self._datum_epsg, value
                )

                self._east, self._north = self.project_point(
                    self.longitude, self.latitude, value, self.utm_epsg
                )

            elif (
                self.datum_epsg is not None
                and self.east != 0
                and self.north != 0
                and self.latitude == 0
                and self.longitude == 0
            ):
                self._longitude, self._latitude = self.project_point(
                    self.east,
                    self.north,
                    self.utm_epsg,
                    value,
                )
            self._datum_epsg = value

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
                self._east, self._north = self.project_point(
                    self.east, self.north, self._utm_epsg, value
                )

                # reproject lat and lon base on new UTM datum
                self._latitude, self._longitude = self.project_point(
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
                self._east, self._north = self.project_point(
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
            self._longitude, self._latitude = self.project_point(
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
            self._longitude, self._latitude = self.project_point(
                self._east, self._north, self.utm_epsg, self.datum_epsg
            )

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, lat):
        self._latitude = assert_lat_value(lat)
        if self.utm_epsg is not None and self.datum_epsg is not None:
            self._east, self._north = self.project_point(
                self._latitude, self.longitude, self.datum_epsg, self.utm_epsg
            )

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, lon):
        self._longitude = assert_lon_value(lon)
        if self.utm_epsg is not None and self.datum_epsg is not None:
            self._east, self._north = self.project_point(
                self._latitude, self.longitude, self.datum_epsg, self.utm_epsg
            )

    @property
    def elevation(self):
        return self._elevation

    @elevation.setter
    def elevation(self, elev):
        self._elevation = assert_elevation_value(elev)
