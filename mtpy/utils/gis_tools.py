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
from mtpy.utils.mtpylog import MtPyLog
from mtpy.utils import HAS_GDAL, EPSG_DICT, NEW_GDAL

if HAS_GDAL:
    from osgeo import osr
    from osgeo.ogr import OGRERR_NONE

else:
    import pyproj

_logger = MtPyLog.get_mtpy_logger(__name__)
if NEW_GDAL:
    _logger.info('INFO: GDAL version 3 detected')

# =============================================================================
# GIS Error container
# =============================================================================


class GISError(Exception):
    pass

# ==============================================================================
# Make sure lat and lon are in decimal degrees
# ==============================================================================


def _assert_minutes(minutes):
    assert 0 <= minutes < 60., \
        'minutes needs to be <60 and >0, currently {0:.0f}'.format(minutes)

    return minutes


def _assert_seconds(seconds):
    assert 0 <= seconds < 60., \
        'seconds needs to be <60 and >0, currently {0:.3f}'.format(seconds)
    return seconds


def convert_position_str2float(position_str):
    """
    Convert a position string in the format of DD:MM:SS to decimal degrees

    :type position_str: string [ 'DD:MM:SS.ms' | 'DD.degrees' ]
    :param position_str: degrees of latitude or longitude

    :rtype: float
    :return: latitude or longitude in decimal degrees

    :Example: ::

        >>> from mtpy.utils import gis_tools
        >>> gis_tools.convert_position_str2float('-118:34:56.3')
        -118.58230555555555

    """

    if position_str in [None, 'None']:
        return None

    if ':' in position_str:
        if position_str.count(':') != 2:
            msg = '{0} not correct format.\n'.format(position_str) +\
                  'Position needs to be DD:MM:SS.ms'
            raise GISError(msg)
        p_list = position_str.split(':')
        deg = float(p_list[0])
        minutes = _assert_minutes(float(p_list[1]))
        sec = _assert_seconds(float(p_list[2]))
        sign = np.sign(deg)

        position_value = sign * (abs(deg) + minutes / 60. + sec / 3600.)
    else:
        try:
            position_value = float(position_str)
        except ValueError:
            msg = '{0} not correct format.\n'.format(position_str) +\
                  'Position needs to be DD.decimal_degrees'
            raise GISError(msg)

    return position_value


def assert_lat_value(latitude):
    """
    make sure latitude is in decimal degrees
    """
    if latitude in [None, 'None']:
        return None
    try:
        lat_value = float(latitude)

    except TypeError:
        return None

    except ValueError:
        lat_value = convert_position_str2float(latitude)

    if abs(lat_value) >= 90:
        raise GISError('|Latitude = {0:.5f}| > 90, unacceptable!'.format(
            lat_value))

    return lat_value


def assert_lon_value(longitude):
    """
    make sure longitude is in decimal degrees
    """
    if longitude in [None, 'None']:
        return None
    try:
        lon_value = float(longitude)

    except TypeError:
        return None

    except ValueError:
        lon_value = convert_position_str2float(longitude)

    if abs(lon_value) >= 180:
        raise GISError('|Longitude = {0:.5f}| > 180, unacceptable!'.format(
                       lon_value))

    return lon_value


def assert_elevation_value(elevation):
    """
    make sure elevation is a floating point number
    """

    try:
        elev_value = float(elevation)
    except (ValueError, TypeError):
        elev_value = 0.0
        _logger.warn('{0} is not a number, setting elevation to 0'.format(elevation))

    return elev_value


def convert_position_float2str(position):
    """
    convert position float to a string in the format of DD:MM:SS

    :type position: float
    :param position: decimal degrees of latitude or longitude

    :rtype: float
    :return: latitude or longitude in DD:MM.SS.ms 

    :Example: ::
        >>> import mtpy.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_float2str(-118.34563)
        '-118:34:56.30'

    """

    if not isinstance(position, float):
        raise GISError('Given value is not a float')

    deg = int(position)
    minutes = (abs(position) - abs(deg)) * 60.

    # need to round seconds to 4 decimal places otherwise machine precision
    # keeps the 60 second roll over and the string is incorrect.
    sec = np.round((minutes - int(minutes)) * 60., 4)
    if sec >= 60.:
        minutes += 1
        sec = 0

    if int(minutes) == 60:
        deg += 1
        minutes = 0

    return '{0:.0f}:{1:02.0f}:{2:05.2f}'.format(deg, int(minutes), sec)


# ==============================================================================
# Project a point
# ==============================================================================
def get_utm_zone(latitude, longitude):
    """
    Get utm zone from a given latitude and longitude

    :param latitude: latitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type latitude: [ string | float ]

    :param longitude: longitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type longitude: [ string | float ]

    :return: zone number
    :rtype: int

    :return: is northern
    :rtype: [ True | False ]

    :return: UTM zone
    :rtype: string

    :Example: ::

        >>> lat, lon = ('-34:17:57.99', 149.2010301)
        >>> zone_number, is_northing, utm_zone = gis_tools.get_utm_zone(lat, lon)
        >>> print(zone_number, is_northing, utm_zone)
        (55, False, '55H')
    """
    latitude = assert_lat_value(latitude)
    longitude = assert_lon_value(longitude)

    zone_number = (int(1 + (longitude + 180.0) / 6.0))
    is_northern = bool(latitude >= 0)
    n_str = utm_letter_designator(latitude)

    return zone_number, is_northern, '{0:02.0f}{1}'.format(zone_number, n_str)


def utm_letter_designator(latitude):
    """
    Get the UTM zone letter designation for a given latitude

    :param latitude: latitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type latitude: [ string | float ]

    :return: UTM zone letter designation
    :rtype: string

    :Example: ::

        >>> gis_utils.utm_letter_designator('-34:17:57.99')
        H
    """
    latitude = assert_lat_value(latitude)

    letter_dict = {'C': (-80, -72),
                   'D': (-72, -64),
                   'E': (-64, -56),
                   'F': (-56, -48),
                   'G': (-48, -40),
                   'H': (-40, -32),
                   'J': (-32, -24),
                   'K': (-24, -16),
                   'L': (-16, -8),
                   'M': (-8, 0),
                   'N': (0, 8),
                   'P': (8, 16),
                   'Q': (16, 24),
                   'R': (24, 32),
                   'S': (32, 40),
                   'T': (40, 48),
                   'U': (48, 56),
                   'V': (56, 64),
                   'W': (64, 72),
                   'X': (72, 84)}

    for key, value in letter_dict.items():
        if value[1] >= latitude >= value[0]:
            return key

    return 'Z'


def split_utm_zone(utm_zone):
    """
    Split utme zone into zone number and is northing

    :param utm_zone: utm zone string as {0-9}{0-9}{C-X} or {+,-}{0-9}{0-9}
    :type utm_zone: [ string | int ]

    :return: utm zone number
    :rtype: int

    :return: is_northern
    :rtype: boolean [ True | False ]

    :Example: ::

        >>> gis_tools.split_utm_zone('11S')
        11, True
    """
    utm_zone = validate_utm_zone(utm_zone)

    if isinstance(utm_zone, int):
        # std UTM code returned by gdal
        is_northern = False if utm_zone < 0 else True
        zone_number = abs(utm_zone)
    elif isinstance(utm_zone, str):
        zone_number = int(utm_zone[0:-1])
        is_northern = True if utm_zone[-1].lower() > 'n' else False

    else:
        msg = "utm_zone type {0}, {1} not supported".format(type(utm_zone),
                                                            str(utm_zone))
        raise NotImplementedError(msg)

    return zone_number, is_northern


def utm_zone_to_epsg(zone_number, is_northern):
    """
    get epsg code (WGS84 datum) for a given utm zone

    :param zone_number: UTM zone number
    :type zone_number: int

    :param is_northing: Boolean of UTM is in northern hemisphere
    :type is_northing: [ True | False ]

    :return: EPSG number
    :rtype: int

    :Example: ::

        >>> gis_tools.utm_zone_to_epsg(55, False)
        32755

    """
    for key in list(EPSG_DICT.keys()):
        val = EPSG_DICT[key]
        if ('+zone={:<2}'.format(zone_number) in val) and \
                ('+datum=WGS84' in val):
            if is_northern:
                if '+south' not in val:
                    return key
            else:
                if '+south' in val:
                    return key


def get_epsg(latitude, longitude):
    """
    get epsg code for the utm projection (WGS84 datum) of a given latitude
    and longitude pair

    :param latitude: latitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type latitude: [ string | float ]

    :param longitude: longitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type longitude: [ string | float ]

    :return: EPSG number
    :rtype: int

    :Example: ::

        >>> gis_tools.get_epsg(-34.299442, '149:12:03.71')
        32755

    """
    zone_number, is_northern, utm_str = get_utm_zone(latitude, longitude)

    return utm_zone_to_epsg(zone_number, is_northern)


def _get_gdal_coordinate_system(datum):
    """
    Get coordinate function from GDAL give a datum or EPSG number

    :param datum: can be a well know datum or an EPSG number
    :type: [ int | string ]

    :return: spatial reference coordinate system
    :rtype: osr.SpatialReference

    """
    # set lat lon coordinate system
    cs = osr.SpatialReference()
    if isinstance(datum, int):
        ogrerr = cs.ImportFromEPSG(datum)
        if ogrerr != OGRERR_NONE:
            raise GISError("GDAL/osgeo ogr error code: {}".format(ogrerr))
    elif isinstance(datum, str):
        ogrerr = cs.SetWellKnownGeogCS(datum)
        if ogrerr != OGRERR_NONE:
            raise GISError("GDAL/osgeo ogr error code: {}".format(ogrerr))
    else:
        raise GISError("""datum {0} not understood, needs to be EPSG as int
                           or a well known datum as a string""".format(datum))

    return cs


def validate_epsg(epsg):
    """
    Make sure epsg is an integer

    :param epsg: EPSG number
    :type epsg: [ int | str ]

    :return: EPSG number
    :rtype: int

    """
    if isinstance(epsg, int):
        return epsg

    else:
        if epsg is None:
            return None

        try:
            epsg = int(epsg)
            return epsg
        except ValueError:
            raise GISError('EPSG must be an integer')


def validate_utm_zone(utm_zone):
    """
    Make sure utm zone is a valid string

    :param utm_zone: UTM zone as {0-9}{0-9}{C-X} or {+, -}{0-9}{0-9}
    :type utm_zone: [ int | string ]
    :return: valid UTM zone
    :rtype: [ int | string ]

    """

    if utm_zone is None:
        return None
    # JP: if its unicode then its already a valid string in python 3
    if isinstance(utm_zone, (np.bytes_, bytes)):
        utm_zone = utm_zone.decode('UTF-8')
    elif isinstance(utm_zone, (float, int)):
        utm_zone = int(utm_zone)
    else:
        utm_zone = str(utm_zone)

    return utm_zone


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
        values = np.array([values], dtype=np.float)
    elif isinstance(values, (list, tuple)):
        values = np.array(values, dtype=np.float)
    elif isinstance(values, str):
        values = [ss.strip() for ss in values.strip().split(',')]
        values = np.array(values)
    elif isinstance(values, np.ndarray):
        values = values.astype(np.float)

    # Flatten to 1D
    values = values.flatten()

    if location_type in ['lat', 'latitude']:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lat_value(value)
            except GISError as error:
                raise GISError('{0}\n Bad input value at index {1}'.format(
                               error, ii))
        values = values.astype(np.float)

    if location_type in ['lon', 'longitude']:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lon_value(value)
            except GISError as error:
                raise GISError('{0}\n Bad input value at index {1}'.format(
                               error, ii))
        values = values.astype(np.float)

    return values


def _get_gdal_projection_ll2utm(datum, utm_zone, epsg):
    """
    Get the GDAL transfrom point function for given datum, utm_zone, epsg to
    transform a latitude and longitude point to UTM coordinates.

    ..note:: Have to input either UTM zone or EPSG number

    :param datum: well known datum
    :type datum: string

    :param utm_zone: utm_zone {0-9}{0-9}{C-X} or {+, -}{0-9}{0-9}
    :type utm_zone: [ string | int ]

    :param epsg: EPSG number
    :type epsg: [ int | string ]

    :return: tranform point function
    :rtype: osr.TransformPoint function

    """
    if utm_zone is None and epsg is None:
        raise GISError('Need to input either UTM zone or EPSG number')

    ll_cs = _get_gdal_coordinate_system(datum)

    # project point on to EPSG coordinate system if given
    if isinstance(epsg, int):
        utm_cs = _get_gdal_coordinate_system(validate_epsg(epsg))
    # otherwise project onto given datum
    elif epsg is None:
        utm_cs = _get_gdal_coordinate_system(datum)

        zone_number, is_northern = split_utm_zone(utm_zone)
        utm_cs.SetUTM(zone_number, is_northern)

    return osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoint


def _get_gdal_projection_utm2ll(datum, utm_zone, epsg):
    """
    Get the GDAL transfrom point function for given datum, utm_zone, epsg to
    transform a UTM point to latitude and longitude.

    ..note:: Have to input either UTM zone or EPSG number

    :param datum: well known datum
    :type datum: string

    :param utm_zone: utm_zone {0-9}{0-9}{C-X} or {+, -}{0-9}{0-9}
    :type utm_zone: [ string | int ]

    :param epsg: EPSG number
    :type epsg: [ int | string ]

    :return: tranform point function
    :rtype: osr.TransformPoint function

    """
    if utm_zone is None and epsg is None:
        raise GISError('Need to input either UTM zone or EPSG number')

    # zone_number, is_northern = split_utm_zone(utm_zone)

    if epsg is not None:
        utm_cs = _get_gdal_coordinate_system(validate_epsg(epsg))
    else:
        zone_number, is_northern = split_utm_zone(utm_zone)
        utm_cs = _get_gdal_coordinate_system(datum)
        utm_cs.SetUTM(zone_number, is_northern)

    ll_cs = utm_cs.CloneGeogCS()
    return osr.CoordinateTransformation(utm_cs, ll_cs).TransformPoint


def _get_pyproj_projection(datum, utm_zone, epsg):
    """

    Get the pyproj transfrom point function for given datum, utm_zone, epsg to
    transform either a UTM point to latitude and longitude, or latitude
    and longitude point to UTM.

    ..note:: Have to input either UTM zone or EPSG number

    :param datum: well known datum
    :type datum: string

    :param utm_zone: utm_zone {0-9}{0-9}{C-X} or {+, -}{0-9}{0-9}
    :type utm_zone: [ string | int ]

    :param epsg: EPSG number
    :type epsg: [ int | string ]

    :return: pyproj transform function
    :rtype: pyproj.Proj function

    """
    if utm_zone is None and epsg is None:
        raise GISError('Need to input either UTM zone or EPSG number')

    if isinstance(epsg, int):
        pp = pyproj.Proj('+init=EPSG:%d' % (epsg))

    elif epsg is None:
        zone_number, is_northern = split_utm_zone(utm_zone)
        zone = 'north' if is_northern else 'south'
        proj_str = '+proj=utm +zone=%d +%s +datum=%s' % (zone_number, zone,
                                                         datum)
        pp = pyproj.Proj(proj_str)

    return pp


def project_point_ll2utm(lat, lon, datum='WGS84', utm_zone=None, epsg=None):
    """
    Project a point that is in latitude and longitude to the specified
    UTM coordinate system.

    :param latitude: latitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type latitude: [ string | float ]

    :param longitude: longitude in [ 'DD:mm:ss.ms' | 'DD.decimal' | float ]
    :type longitude: [ string | float ]

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
    lat = validate_input_values(lat, location_type='lat')
    lon = validate_input_values(lon, location_type='lon')

    if utm_zone in [None, 'none', 'None']:
        # get the UTM zone in the datum coordinate system, otherwise
        zone_number, is_northern, utm_zone = get_utm_zone(lat.mean(),
                                                          lon.mean())
    epsg = validate_epsg(epsg)
    if HAS_GDAL:
        ll2utm = _get_gdal_projection_ll2utm(datum, utm_zone, epsg)
    else:
        ll2utm = _get_pyproj_projection(datum, utm_zone, epsg)

    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(lat, dtype=[('easting', np.float),
                                                ('northing', np.float),
                                                ('elev', np.float),
                                                ('utm_zone', 'U3')])

    for ii in range(lat.size):
        if NEW_GDAL:
            point = ll2utm(lat[ii], lon[ii])
        else:
            point = ll2utm(lon[ii], lat[ii])

        projected_point['easting'][ii] = point[0]
        projected_point['northing'][ii] = point[1]
        if HAS_GDAL:
            projected_point['elev'][ii] = point[2]

        projected_point['utm_zone'][ii] = utm_zone

    # if just projecting one point, then return as a tuple so as not to break
    # anything.  In the future we should adapt to just return a record array
    if len(projected_point) == 1:
        return (projected_point['easting'][0],
                projected_point['northing'][0],
                projected_point['utm_zone'][0])
    else:
        return np.rec.array(projected_point)


def project_point_utm2ll(easting, northing, utm_zone, datum='WGS84', epsg=None):
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
    epsg = validate_epsg(epsg)

    if HAS_GDAL:
        utm2ll = _get_gdal_projection_utm2ll(datum, utm_zone, epsg)
    else:
        utm2ll = _get_pyproj_projection(datum, utm_zone, epsg)

    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(easting,
                                    dtype=[('latitude', np.float),
                                           ('longitude', np.float)])
    for ii in range(easting.size):
        if HAS_GDAL:
            point = utm2ll(easting[ii], northing[ii], 0.0)

            try:
                assert_lat_value(point[0])
                projected_point['latitude'][ii] = round(point[0], 6)
                projected_point['longitude'][ii] = round(point[1], 6)
            except GISError:
                projected_point['latitude'][ii] = round(point[1], 6)
                projected_point['longitude'][ii] = round(point[0], 6)

        else:
            point = utm2ll(easting[ii], northing[ii], inverse=True)
            projected_point['latitude'][ii] = round(point[1], 6)
            projected_point['longitude'][ii] = round(point[0], 6)

    # if just projecting one point, then return as a tuple so as not to break
    # anything.  In the future we should adapt to just return a record array
    if len(projected_point) == 1:
        return (projected_point['latitude'][0],
                projected_point['longitude'][0])
    else:
        return np.rec.array(projected_point)


def epsg_project(x, y, epsg_from, epsg_to):
    """
    project some xy points using the pyproj modules
    """

    try:
        import pyproj
    except ImportError:
        print("please install pyproj")
        return
    if epsg_from is not None:
        try:
            p1 = pyproj.Proj(EPSG_DICT[epsg_from])
            p2 = pyproj.Proj(EPSG_DICT[epsg_to])
        except KeyError:
            print("Surface or data epsg either not in dictionary or None")
            return

    return pyproj.transform(p1, p2, x, y)


def utm_wgs84_conv(lat, lon):
    """
    Bidirectional UTM-WGS84 converter https://github.com/Turbo87/utm/blob/master/utm/conversion.py
    :param lat:
    :param lon:
    :return: tuple(e, n, zone, lett)
    """

    import utm  # pip install utm
    tup = utm.from_latlon(lat, lon)

    (new_lat, new_lon) = utm.to_latlon(tup[0], tup[1], tup[2], tup[3])
    # print (new_lat,new_lon)  # should be same as the input param

    # checking correctess
    if abs(lat - new_lat) > 1.0 * np.e - 10:
        print("Warning: lat and new_lat should be equal!")

    if abs(lon - new_lon) > 1.0 * np.e - 10:
        print("Warning: lon and new_lon should be equal!")

    return tup


#################################################################
# Example usages of this script/module
# python gis_tools.py
# =================================================================
if __name__ == "__main__":

    my_lat = -35.0
    my_lon = 149.5
    utm_point = project_point_ll2utm(my_lat, my_lon)
    print("project_point_ll2utm(mylat, mylon) =:  ", utm_point)
