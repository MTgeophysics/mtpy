# -*- coding: utf-8 -*-
"""
GIS_TOOLS
==================

This module contains tools to help project between coordinate systems.  The 
module will first use GDAL if installed.  If GDAL is not installed then 
pyproj is used. A test has been made for new versions of GDAL which swap the
input lat and lon when using transferPoint, so the user should not have to 
worry about which version they have. 

Main functions are:
    
    * project_point_ll2utm
    * project_point_utm2ll
    
These can take in a point or an array or list of points to project.

latitude and longitude can be input as:
    * 'DD:mm:ss.ms'
    * 'DD.decimal_degrees'
    * float(DD.decimal_degrees)

Created on Fri Apr 14 14:47:48 2017
Revised: 5/2020 JP 

@author: jrpeacock
"""

# ==============================================================================
# Imports
# ==============================================================================
import numpy as np
from loguru import logger
from pyproj import Transformer, CRS

# =============================================================================
# GIS Error container
# =============================================================================


class GISError(Exception):
    pass


# ==============================================================================
# Make sure lat and lon are in decimal degrees
# ==============================================================================


def assert_minutes(minutes):
    if not 0 <= minutes < 60.0:
        msg = (
            "minutes should be 0 < > 60, currently {0:.0f}".format(minutes)
            + " conversion will account for non-uniform"
            + "timne. Be sure to check accuracy."
        )
        logger.warning(msg)

    return minutes


def assert_seconds(seconds):
    if not 0 <= seconds < 60.0:
        msg = (
            "seconds should be 0 < > 60, currently {0:.0f}".format(seconds)
            + " conversion will account for non-uniform"
            + "timne. Be sure to check accuracy."
        )
        logger.warning(msg)

    return seconds


def convert_position_str2float(position_str):
    """
    Convert a position string in the format of DD:MM:SS to decimal degrees

    :param position: latitude or longitude om DD:MM:SS.ms
    :type position: float

    :returns: latitude or longitude as a float
    """

    if position_str in [None, "None"]:
        return None

    p_list = position_str.split(":")
    if len(p_list) != 3:
        msg = "{0} not correct format, should be DD:MM:SS".format(position_str)
        logger.error(msg)
        raise ValueError(msg)

    deg = float(p_list[0])
    minutes = assert_minutes(float(p_list[1]))
    sec = assert_seconds(float(p_list[2]))

    # get the sign of the position so that when all are added together the
    # position is in the correct place
    sign = 1
    if deg < 0:
        sign = -1

    position_value = sign * (abs(deg) + minutes / 60.0 + sec / 3600.0)

    logger.debug("Converted {0} to {1}".format(position_str, position_value))

    return position_value


def assert_lat_value(latitude):
    """
    Make sure the latitude value is in decimal degrees, if not change it.
    And that the latitude is within -90 < lat > 90.

    :param latitude: latitude in decimal degrees or other format
    :type latitude: float or string
    """
    if latitude in [None, "None", "none", "unknown"]:
        logger.debug("Latitude is None, setting to 0")
        return 0.0
    try:
        lat_value = float(latitude)

    except TypeError:
        logger.debug("Could not convert {0} setting to 0".format(latitude))
        return 0.0

    except ValueError:
        logger.debug("Latitude is a string {0}".format(latitude))
        lat_value = convert_position_str2float(latitude)

    if abs(lat_value) >= 90:
        msg = (
            "latitude value = {0} is unacceptable!".format(lat_value)
            + ".  Must be |Latitude| > 90"
        )
        logger.error(msg)
        raise ValueError(msg)

    return lat_value


def assert_lon_value(longitude):
    """
    Make sure the longitude value is in decimal degrees, if not change it.
    And that the latitude is within -180 < lat > 180.

    :param latitude: longitude in decimal degrees or other format
    :type latitude: float or string
    """
    if longitude in [None, "None", "none", "unknown"]:
        logger.debug("Longitude is None, setting to 0")
        return 0.0
    try:
        lon_value = float(longitude)

    except TypeError:
        logger.debug("Could not convert {0} setting to 0".format(longitude))
        return 0.0

    except ValueError:
        logger.debug("Longitude is a string {0}".format(longitude))
        lon_value = convert_position_str2float(longitude)

    if abs(lon_value) >= 180:
        msg = (
            "longitude value = {0} is unacceptable!".format(lon_value)
            + ".  Must be |longitude| > 180"
        )
        logger.error(msg)
        raise ValueError(msg)

    return lon_value


def assert_elevation_value(elevation):
    """
    make sure elevation is a floating point number

    :param elevation: elevation as a float or string that can convert
    :type elevation: float or str
    """

    try:
        elev_value = float(elevation)
    except (ValueError, TypeError):
        msg = "Could not convert {0} to a number setting to 0".format(
            elevation
        )
        logger.debug(msg)
        elev_value = 0.0

    return elev_value


def convert_position_float2str(position):
    """
    Convert position float to a string in the format of DD:MM:SS.

    :param position: decimal degrees of latitude or longitude
    :type position: float

    :returns: latitude or longitude in format of DD:MM:SS.ms
    """

    assert type(position) is float, "Given value is not a float"

    deg = int(position)
    sign = 1
    if deg < 0:
        sign = -1

    deg = abs(deg)
    minutes = (abs(position) - deg) * 60.0
    # need to round seconds to 4 decimal places otherwise machine precision
    # keeps the 60 second roll over and the string is incorrect.
    sec = np.round((minutes - int(minutes)) * 60.0, 4)
    if sec >= 60.0:
        minutes += 1
        sec = 0

    if int(minutes) == 60:
        deg += 1
        minutes = 0

    position_str = "{0}:{1:02.0f}:{2:05.2f}".format(
        sign * int(deg), int(minutes), sec
    )
    logger.debug("Converted {0} to {1}".format(position, position_str))

    return position_str


# ==============================================================================
# Project a point
# ==============================================================================


def validate_input_values(values, location_type=None):
    """
    make sure the input values for lat, lon, easting, northing will be an
    numpy array with a float data type

    can input a string as a comma separated list

    :param values: values to project, can be given as:
        * float
        * string of a single value or a comma separate string '34.2, 34.5'
        * list of floats or string
        * numpy.ndarray

    :type values: [ float | string | list | numpy.ndarray ]

    :return: array of floats
    :rtype: numpy.ndarray(dtype=float)

    """
    if isinstance(values, (int, float)):
        values = np.array([values], dtype=float)
    elif isinstance(values, (list, tuple)):
        values = np.array(values, dtype=float)
    elif isinstance(values, str):
        values = [ss.strip() for ss in values.strip().split(",")]
        values = np.array(values)
    elif isinstance(values, np.ndarray):
        values = values.astype(float)
    # Flatten to 1D
    values = values.flatten()

    if location_type in ["lat", "latitude"]:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lat_value(value)
            except GISError as error:
                raise GISError(
                    "{0}\n Bad input value at index {1}".format(error, ii)
                )
        values = values.astype(float)
    if location_type in ["lon", "longitude"]:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lon_value(value)
            except GISError as error:
                raise GISError(
                    "{0}\n Bad input value at index {1}".format(error, ii)
                )
        values = values.astype(float)
    return values


def project_point(x, y, old_epsg, new_epsg):
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
    if isinstance(x, (np.ndarray, list, tuple)):
        x = np.array(x)
        if (x == 0).all():
            print(x)
            raise ValueError("Should not project with 0 value")
    if isinstance(y, (np.ndarray, list, tuple)):
        y = np.array(y)
        if (y == 0).all():
            raise ValueError("Should not project with 0 value")

    if isinstance(x, (float, int)):
        x = float(x)
        if x == 0:
            raise ValueError("Should not project with 0 value")

    if isinstance(y, (float, int)):
        y = float(y)
        if y == 0:
            raise ValueError("Should not project with 0 value")

    old_crs = CRS.from_user_input(old_epsg)
    new_crs = CRS.from_user_input(new_epsg)

    transformer = Transformer.from_crs(old_crs, new_crs, always_xy=True)

    return transformer.transform(x, y)


def project_point_ll2utm(lat, lon, datum="WGS84", epsg=None):
    """
    Project a point that is in latitude and longitude to the specified
    UTM coordinate system.

    :param latitude: latitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type latitude: [ string | float ]

    :param longitude: longitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type longitude: [ string | float ]

    :param datum: well known datum
    :type datum: string

    :param epsg: EPSG number defining projection
                (see http://spatialreference.org/ref/ for moreinfo)
                Overrides utm_zone if both are provided
    :type epsg: [ int | string ]

    :return: project point(s)
    :rtype: tuple if a single point, np.recarray if multiple points
        * tuple is (easting, northing,utm_zone)
        * recarray has attributes (easting, northing, utm_zone, elevation)

    :Single Point: ::

        >>> gis_tools.project_point_ll2utm('-34:17:57.99', '149.2010301')
        (702562.6911014864, 6202448.5654573515, '55H')

    :Multiple Points: ::

        >>> lat = np.arange(20, 40, 5)
        >>> lon = np.arange(-110, -90, 5)
        >>> gis_tools.project_point_ll2utm(lat, lon, datum='NAD27')
        rec.array([( -23546.69921068, 2219176.82320353, 0., '13R'),
                   ( 500000.        , 2764789.91224626, 0., '13R'),
                   ( 982556.42985037, 3329149.98905941, 0., '13R'),
                   (1414124.6019547 , 3918877.48599922, 0., '13R')],
                  dtype=[('easting', '<f8'), ('northing', '<f8'),
                         ('elev', '<f8'), ('utm_zone', '<U3')])


    """
    if lat is None or lon is None:
        return None, None, None
    # make sure the lat and lon are in decimal degrees
    lat = validate_input_values(lat, location_type="lat")
    lon = validate_input_values(lon, location_type="lon")

    epsg_crs = CRS.from_epsg(epsg)
    datum_crs = CRS.from_user_input(datum)
    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(
        lat,
        dtype=[
            ("easting", float),
            ("northing", float),
            ("elev", float),
            ("utm_zone", "U4"),
        ],
    )

    easting, northing = project_point(lon, lat, datum_crs, epsg_crs)
    projected_point["easting"][:] = easting
    projected_point["northing"][:] = northing
    projected_point["utm_zone"][:] = epsg_crs.utm_zone

    # if just projecting one point, then return as a tuple so as not to break
    # anything.  In the future we should adapt to just return a record array
    if len(projected_point) == 1:
        return (
            projected_point["easting"][0],
            projected_point["northing"][0],
            projected_point["utm_zone"][0],
        )
    else:
        return np.rec.array(projected_point)


def project_point_utm2ll(easting, northing, utm_epsg, datum_epsg=4326):
    """
    Project a point that is in UTM to the specified geographic coordinate
    system.

    :param easting: easting in meters
    :type easting: float

    :param northing: northing in meters
    :type northing: float

    :param datum: well known datum
    :type datum: string

    :param utm_zone: utm_zone {0-9}{0-9}{C-X} or {+, -}{0-9}{0-9}
    :type utm_zone: [ string | int ]

    :param epsg: EPSG number defining projection
                (see http://spatialreference.org/ref/ for moreinfo)
                Overrides utm_zone if both are provided
    :type epsg: [ int | string ]

    :return: project point(s)
    :rtype: tuple if a single point, np.recarray if multiple points
        * tuple is (easting, northing,utm_zone)
        * recarray has attributes (easting, northing, utm_zone, elevation)

    :Single Point: ::

        >>> gis_tools.project_point_utm2ll(670804.18810336,
        ...                                4429474.30215206,
        ...                                datum='WGS84',
        ...                                utm_zone='11T',
        ...                                epsg=26711)
        (40.000087, -114.999128)

    :Multiple Points: ::

        >>> gis_tools.project_point_utm2ll([670804.18810336, 680200],
        ...                                [4429474.30215206, 4330200],
        ...                                datum='WGS84', utm_zone='11T',
        ...                                epsg=26711)
        rec.array([(40.000087, -114.999128), (39.104208, -114.916058)],
                  dtype=[('latitude', '<f8'), ('longitude', '<f8')])

    """
    easting = validate_input_values(easting)
    northing = validate_input_values(northing)

    utm_crs = CRS.from_epsg(utm_epsg)
    datum_crs = CRS.from_user_input(datum_epsg)

    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(
        easting, dtype=[("latitude", float), ("longitude", float)]
    )

    lon, lat = project_point(easting, northing, utm_crs, datum_crs)
    projected_point["latitude"][:] = np.round(lat, 6)
    projected_point["longitude"][:] = np.round(lon, 6)
    # if just projecting one point, then return as a tuple so as not to break
    # anything.  In the future we should adapt to just return a record array
    if len(projected_point) == 1:
        return (
            projected_point["latitude"][0],
            projected_point["longitude"][0],
        )
    else:
        return np.rec.array(projected_point)
