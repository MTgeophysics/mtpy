# -*- coding: utf-8 -*-
"""
Might think about adding declination

Created on Mon Oct  3 15:04:12 2022

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
from copy import deepcopy
from pyproj import CRS
import numpy as np
from loguru import logger

from mtpy.utils.gis_tools import project_point

from mt_metadata.transfer_functions.tf import Station, Survey, Run
from mt_metadata.transfer_functions.io.tools import get_nm_elev

# =============================================================================


class MTLocation:
    """
    Location for a MT site or point measurement

    """

    def __init__(self, survey_metadata=None, **kwargs):

        self.logger = logger
        if survey_metadata is None:
            self._survey_metadata = self._initiate_metadata()
        else:
            self._survey_metadata = self._validate_metadata(survey_metadata)

        self._east = 0
        self._north = 0
        self._datum_crs = CRS.from_epsg(4326)
        self._utm_crs = None
        self._geoid_crs = None
        self.model_east = 0
        self.model_north = 0
        self.model_elevation = 0
        self.profile_offset = 0

        self._key_attrs = [
            "latitude",
            "longitude",
            "elevation",
            "east",
            "north",
            "model_east",
            "model_north",
            "model_elevation",
            "datum_crs",
            "utm_crs",
            "datum_epsg",
            "utm_epsg",
            "profile_offset",
        ]

        for key, value in kwargs.items():
            if key in self._key_attrs:
                setattr(self, key, value)

        if self.east != 0 and self.north != None:
            if self.utm_crs is None:
                raise ValueError(
                    "Need to input UTM CRS if only setting east and north"
                )

    def _initiate_metadata(self):
        survey_metadata = Survey(id=0)
        survey_metadata.add_station(Station(id=0))
        survey_metadata.stations[0].add_run(Run(id=0))

        return survey_metadata

    def _validate_metadata(self, survey_metadata):
        if not isinstance(survey_metadata, Survey):
            raise TypeError(
                "Input metadata must be type "
                "mt_metadata.transfer_functions.tf.Survey, "
                f"not {type(survey_metadata)}."
            )
        if len(survey_metadata.stations) < 1:
            survey_metadata.add_station(Station(id=0))

        if len(survey_metadata.stations[0].runs) < 1:
            survey_metadata.stations[0].add_run(Run(id=0))

        return survey_metadata

    def __str__(self):
        lines = ["MT Location: ", "-" * 20]
        lines.append(f"  Latitude (deg):   {self.latitude:.6f}")
        lines.append(f"  Longitude (deg):  {self.longitude:.6f}")
        lines.append(f"  Elevation (m):    {self.elevation:.4f}")
        lines.append(f"  Datum crs:        {self.datum_crs}")
        lines.append("")
        lines.append(f"  Easting (m):      {self.east:.3f}")
        lines.append(f"  Northing (m):     {self.north:.3f}")
        lines.append(f"  UTM crs:          {self.utm_crs}")
        lines.append("")
        lines.append(f"  Model Easting (m):      {self.model_east:.3f}")
        lines.append(f"  Model Northing (m):     {self.model_north:.3f}")
        lines.append(f"  Model Elevation (m):    {self.model_elevation:.3f}")
        lines.append(f"  Profile Offset (m):     {self.profile_offset:.3f}")

        return "\n".join(lines)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        """
        equals
        :param other: DESCRIPTION
        :type other: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if not isinstance(other, MTLocation):
            raise TypeError(f"Can not compare MTLocation with {type(other)}")

        for key in self._key_attrs:
            og_value = getattr(self, key)
            other_value = getattr(other, key)

            if isinstance(og_value, float):
                if not np.isclose(og_value, other_value):
                    self.logger.info(
                        f"{key} not equal {og_value} != {other_value}"
                    )
                    return False
            else:
                if not og_value == other_value:
                    self.logger.info(
                        f"{key} not equal {og_value} != {other_value}"
                    )
                    return False
        return True

    def copy(self):
        copied = type(self)()
        copied._survey_metadata = self._survey_metadata.copy()
        # not sure why this is needed, survey metadata copies fine, but here
        # it does not.
        if len(copied._survey_metadata.stations) == 0:
            copied._survey_metadata.add_station(
                self._survey_metadata.stations[0]
            )
        for key in self._key_attrs:
            setattr(copied, key, deepcopy(getattr(self, key)))

        return copied

    @property
    def datum_crs(self):
        if self._datum_crs is not None:
            return self._datum_crs

    @property
    def datum_name(self):
        if self._datum_crs is not None:
            return self._datum_crs.name

    @property
    def datum_epsg(self):
        if self._datum_crs is not None:
            return self._datum_crs.to_epsg()

    @datum_epsg.setter
    def datum_epsg(self, value):
        if value not in ["", None, "None"]:
            self.datum_crs = value

    @datum_crs.setter
    def datum_crs(self, value):
        if value in [None, "None", "none", "null", ""]:
            return

        new_crs = CRS.from_user_input(value)

        if new_crs != self._datum_crs:
            if (
                self._datum_crs is not None
                and self.latitude != 0
                and self.longitude != 0
            ):
                (
                    self._survey_metadata.stations[0].location.longitude,
                    self._survey_metadata.stations[0].location.latitude,
                ) = project_point(
                    self.longitude, self.latitude, self._datum_crs, new_crs
                )

                self._east, self._north = project_point(
                    self.longitude, self.latitude, new_crs, self.utm_crs
                )

            elif (
                self.datum_crs is not None
                and self.east != 0
                and self.north != 0
                and self.latitude == 0
                and self.longitude == 0
            ):
                (
                    self._survey_metadata.stations[0].location.longitude,
                    self._survey_metadata.stations[0].location.latitude,
                ) = project_point(
                    self.east,
                    self.north,
                    self.utm_crs,
                    new_crs,
                )
            self._datum_crs = new_crs

    @property
    def utm_crs(self):
        if self._utm_crs is not None:
            return self._utm_crs

    @property
    def utm_name(self):
        if self._utm_crs is not None:
            return self._utm_crs.name

    @property
    def utm_epsg(self):
        if self._utm_crs is not None:
            return self._utm_crs.to_epsg()

    @utm_epsg.setter
    def utm_epsg(self, value):
        if value not in ["", None, "None"]:
            self.utm_crs = value

    @property
    def utm_zone(self):
        if self._utm_crs is not None:
            return self._utm_crs.utm_zone

    @utm_crs.setter
    def utm_crs(self, value):
        if value in [None, "None", "none", "null", ""]:
            return

        new_crs = CRS.from_user_input(value)
        if value != self._utm_crs:
            # reproject easting, northing to new zone
            if (
                self._utm_crs is not None
                and self.east != 0
                and self.north != 0
            ):
                self._east, self._north = project_point(
                    self.east, self.north, self._utm_crs, new_crs
                )

            if (
                self.datum_crs is not None
                and self.east != 0
                and self.north != 0
            ):
                # reproject lat and lon base on new UTM datum
                (
                    self._survey_metadata.stations[0].location.longitude,
                    self._survey_metadata.stations[0].location.latitude,
                ) = project_point(
                    self.east,
                    self.north,
                    new_crs,
                    self.datum_crs,
                )

            # if east and north == 0 and lat and lon != 0 project to utm
            elif (
                self.datum_crs is not None
                and self.east == 0
                and self.north == 0
                and self.latitude != 0
                and self.longitude != 0
            ):
                self._east, self._north = project_point(
                    self.longitude,
                    self.latitude,
                    self.datum_crs,
                    new_crs,
                )

            self._utm_crs = new_crs

    @property
    def east(self):
        """easting"""
        return self._east

    @east.setter
    def east(self, value):
        """set east"""
        self._east = value
        if (
            self.datum_crs is not None
            and self.utm_crs is not None
            and self._north != 0
        ):
            (
                self._survey_metadata.stations[0].location.longitude,
                self._survey_metadata.stations[0].location.latitude,
            ) = project_point(
                self._east, self._north, self.utm_crs, self.datum_crs
            )

    @property
    def north(self):
        """northing"""
        return self._north

    @north.setter
    def north(self, value):
        """set north"""
        self._north = value
        if (
            self.datum_crs is not None
            and self.utm_crs is not None
            and self._east != 0
        ):
            (
                self._survey_metadata.stations[0].location.longitude,
                self._survey_metadata.stations[0].location.latitude,
            ) = project_point(
                self._east, self._north, self.utm_crs, self.datum_crs
            )

    @property
    def latitude(self):
        return self._survey_metadata.stations[0].location.latitude

    @latitude.setter
    def latitude(self, lat):
        self._survey_metadata.stations[0].location.latitude = lat
        if (
            self.utm_crs is not None
            and self.datum_crs is not None
            and self._survey_metadata.stations[0].location.longitude != 0
        ):
            self._east, self._north = project_point(
                self._survey_metadata.stations[0].location.longitude,
                self._survey_metadata.stations[0].location.latitude,
                self.datum_crs,
                self.utm_crs,
            )

    @property
    def longitude(self):
        return self._survey_metadata.stations[0].location.longitude

    @longitude.setter
    def longitude(self, lon):
        self._survey_metadata.stations[0].location.longitude = lon
        if (
            self.utm_crs is not None
            and self.datum_crs is not None
            and self._survey_metadata.stations[0].location.latitude != 0
        ):
            self._east, self._north = project_point(
                self._survey_metadata.stations[0].location.longitude,
                self._survey_metadata.stations[0].location.latitude,
                self.datum_crs,
                self.utm_crs,
            )

    @property
    def elevation(self):
        return self._survey_metadata.stations[0].location.elevation

    @elevation.setter
    def elevation(self, elev):
        self._survey_metadata.stations[0].location.elevation = elev

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

    def compute_model_location(self, center_location):
        """
        compute model location based on model center and model epsg

        :param model_center: DESCRIPTION
        :type model_center: TYPE
        :return: DESCRIPTION
        :rtype: TYPE

        """

        self.model_east = self.east - center_location.model_east
        self.model_north = self.north - center_location.model_north
        self.model_elevation = self.elevation - center_location.model_elevation

    def project_onto_profile_line(self, profile_slope, profile_intersection):
        """

        :param profile_slope: DESCRIPTION
        :type profile_slope: TYPE
        :param profile_intersection: DESCRIPTION
        :type profile_intersection: TYPE
        :param units: DESCRIPTION, defaults to "deg"
        :type units: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        if self.utm_crs is None:
            raise ValueError(
                "utm_crs is None, cannot project onto profile line."
            )

        profile_vector = np.array([1, profile_slope], dtype=float)
        profile_vector /= np.linalg.norm(profile_vector)

        station_vector = np.array(
            [self.east, self.north - profile_intersection]
        )

        self.profile_offset = np.linalg.norm(
            np.dot(profile_vector, station_vector) * profile_vector
        )

    def get_elevation_from_national_map(self):
        """
        Get elevation from DEM data of the US National Map.  Plan to extend
        this to the globe.

        Pulls data from the USGS national map DEM

        :return: DESCRIPTION
        :rtype: TYPE

        """

        elev = get_nm_elev(self.latitude, self.longitude)
        if elev != 0:
            self.elevation = elev
        else:
            self.logger.warning(
                "Could not get elevation data, not setting elevation"
            )
