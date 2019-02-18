# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:47:48 2017

@author: jrpeacock
"""

# ==============================================================================
# Imports
# ==============================================================================

from mtpy.utils import HAS_GDAL, EPSG_DICT
from mtpy.utils.decorator import gdal_data_check, deprecated
import numpy as np
from mtpy.utils.mtpylog import MtPyLog

if HAS_GDAL:
    from osgeo import osr
    from osgeo.ogr import OGRERR_NONE
else:
    import pyproj
# end if

_logger = MtPyLog.get_mtpy_logger(__name__)


class GIS_ERROR(Exception):
    pass
#==============================================================================
# Make sure lat and lon are in decimal degrees
#==============================================================================
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
    
    Arguments
    -------------
        **position_str** : string ('DD:MM:SS.ms')
                           degrees of latitude or longitude
                       
    Returns
    --------------
        **position** : float
                       latitude or longitude in decimal degrees
                          
    Example
    -------------
        >>> import mtpy.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_str2float('-118:34:56.3')
    """

    if position_str in [None, 'None']:
        return None
    
    p_list = position_str.split(':')
    if len(p_list) != 3:
        raise ValueError('{0} not correct format, should be DD:MM:SS'.format(position_str))

    deg = float(p_list[0])
    minutes = _assert_minutes(float(p_list[1]))
    sec = _assert_seconds(float(p_list[2]))

    # get the sign of the position so that when all are added together the
    # position is in the correct place
    sign = 1
    if deg < 0:
        sign = -1

    position_value = sign * (abs(deg) + minutes / 60. + sec / 3600.)

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
        print("==> The lat_value =", lat_value)
        raise ValueError('|Latitude| > 90, unacceptable!')

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
        raise ValueError('|Longitude| > 180, unacceptable!')

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
    
    Arguments
    -------------
        **position** : float
                       decimal degrees of latitude or longitude
                       
    Returns
    --------------
        **position_str** : string
                          latitude or longitude in format of DD:MM:SS.ms
                          
    Example
    -------------
        >>> import mtpy.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_float2str(-118.34563)
        
    """

    assert type(position) is float, 'Given value is not a float'

    deg = int(position)
    sign = 1
    if deg < 0:
        sign = -1

    deg = abs(deg)
    minutes = (abs(position) - deg) * 60.
    sec = (minutes - int(minutes)) * 60.
    if sec == 60:
        minutes += 1
        sec = 0

    if minutes == 60:
        deg += 1
        minutes = 0

    position_str = '{0}:{1:02.0f}:{2:02.2f}'.format(sign * int(deg),
                                                    int(minutes),
                                                    float(sec))

    return position_str


# ==============================================================================
# Project a point
# ==============================================================================
@deprecated("NATO UTM zone is used in other part of mtpy, this function is for Standard UTM")
def get_utm_string_from_sr(spatialreference):
    """
    return utm zone string from spatial reference instance
    """
    zone_number = spatialreference.GetUTMZone()
    if zone_number > 0:
        return str(zone_number) + 'N'
    elif zone_number < 0:
        return str(abs(zone_number)) + 'S'
    else:
        return str(zone_number)


def get_utm_zone(latitude, longitude):
    """
    Get utm zone from a given latitude and longitude
    """
    zone_number = (int(1 + (longitude + 180.0) / 6.0))
    ## why is this needed here, GDAL only needs N-orth or S-outh
    ## n_str = _utm_letter_designator(latitude)
    is_northern = bool(latitude >= 0)
    
    # if latitude < 0.0:
    #     n_str = 'S'
    # else:
    #     n_str = 'N'
    n_str = _utm_letter_designator(latitude)

    return zone_number, is_northern, '{0:02.0f}{1}'.format(zone_number, n_str)

def utm_zone_to_epsg(zone_number, is_northern):
    """
    get epsg code (WGS84 datum) for a given utm zone
    
    """
    for key in list(EPSG_DICT.keys()):
        val = EPSG_DICT[key]
        if ('+zone={:<2}'.format(zone_number) in val) and ('+datum=WGS84' in val):
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
    """  
    zone_no,is_northern,utm_str = get_utm_zone(latitude, longitude)
    
    return utm_zone_to_epsg(zone_no,is_northern)

def project_point_ll2utm(lat, lon, datum='WGS84', utm_zone=None, epsg=None):
    """
    Project a point that is in Lat, Lon (will be converted to decimal degrees)
    into UTM coordinates.

    Arguments:
    ---------------
        **lat** : float or string (DD:MM:SS.ms)
                  latitude of point

        **lon** : float or string (DD:MM:SS.ms)
                  longitude of point

        **datum** : string
                    well known datum ex. WGS84, NAD27, NAD83, etc.

        **utm_zone** : string
                       zone number and 'S' or 'N' e.g. '55S'

        **epsg** : int
                   epsg number defining projection (see
                   http://spatialreference.org/ref/ for moreinfo)
                   Overrides utm_zone if both are provided

    Returns:
    --------------
        **proj_point**: tuple(easting, northing, zone)
                        projected point in UTM in Datum

    """
    if lat is None or lon is None:
        return None, None, None

    # make sure the lat and lon are in decimal degrees
    if np.iterable(lon) and np.iterable(lat):
        lat = np.array([assert_lat_value(lat_value) for lat_value in lat])
        lon = np.array([assert_lon_value(lon_value) for lon_value in lon])
        assert lat.size == lon.size
    else:
        lat = np.array([assert_lat_value(lat)])
        lon = np.array([assert_lon_value(lon)])

    if HAS_GDAL:
        # set lat lon coordinate system
        ll_cs = osr.SpatialReference()
        if isinstance(datum, int):
            ogrerr = ll_cs.ImportFromEPSG(datum)
            if ogrerr != OGRERR_NONE:
                raise GIS_ERROR("GDAL/osgeo ogr error code: {}".format(ogrerr))
        elif isinstance(datum, str):
            ogrerr = ll_cs.SetWellKnownGeogCS(datum)
            if ogrerr != OGRERR_NONE:
                raise GIS_ERROR("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            raise GIS_ERROR("""datum {0} not understood, needs to be EPSG as int
                               or a well known datum as a string""".format(datum))

        # set utm coordinate system
        utm_cs = osr.SpatialReference()
    # end if

    # project point on to EPSG coordinate system if given
    if isinstance(epsg, int):
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise GIS_ERROR("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    # otherwise project onto given datum
    elif epsg is None:
        if HAS_GDAL:
            ogrerr = utm_cs.CopyGeogCSFrom(ll_cs)
            if ogrerr != OGRERR_NONE:
                raise GIS_ERROR("GDAL/osgeo ogr error code: {}".format(ogrerr))
        # end if
        if utm_zone is None or not isinstance(None, str) or utm_zone.lower() == 'none':
            # get the UTM zone in the datum coordinate system, otherwise
            zone_number, is_northern, utm_zone = get_utm_zone(lat.mean(),
                                                              lon.mean())
        else:
            # get zone number and is_northern from utm_zone string
            zone_number = int(utm_zone[0:-1])
            is_northern = True if utm_zone[-1].lower() > 'n' else False

        if(HAS_GDAL):
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if

    # return different results depending on if lat/lon are iterable
    projected_point = np.zeros_like(lat, dtype=[('easting', np.float),
                                                ('northing', np.float),
                                                ('elev', np.float),
                                                ('utm_zone', 'S4')])

    if(HAS_GDAL):
        ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoint
    else:
        ll2utm = pp
    # end if

    for ii in range(lat.size):
        point = ll2utm(lon[ii], lat[ii])
        projected_point['easting'][ii] = point[0]
        projected_point['northing'][ii] = point[1]
        if(HAS_GDAL): projected_point['elev'][ii] = point[2]
        projected_point['utm_zone'][ii] = utm_zone if utm_zone is not None else get_utm_zone(lat[ii], lon[ii])[2]
    # end for

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
    Project a point that is in Lat, Lon (will be converted to decimal degrees)
    into UTM coordinates.
    
    Arguments:
    ---------------
        **easting** : float
                    easting coordinate in meters
                    
        **northing** : float
                    northing coordinate in meters
        
        **utm_zone** : string (##N or ##S)
                      utm zone in the form of number and North or South
                      hemisphere, 10S or 03N
        
        **datum** : string
                    well known datum ex. WGS84, NAD27, etc.
                    
    Returns:
    --------------
        **proj_point**: tuple(lat, lon)
                        projected point in lat and lon in Datum, as decimal
                        degrees.
                    
    """
    try:
        easting = float(easting)
    except ValueError:
        raise GIS_ERROR("easting is not a float")
    try:
        northing = float(northing)
    except ValueError:
        raise GIS_ERROR("northing is not a float")

    if HAS_GDAL:
        # set utm coordinate system
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)
    # end if

    if epsg is not None:
        if HAS_GDAL:
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise Exception("GDAL/osgeo ogr error code: {}".format(ogrerr))
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    elif isinstance(utm_zone, str) or isinstance(utm_zone, np.bytes_):
        # the isinstance(utm_zone, str) could be False in python3 due to numpy datatype change.
        # So FZ added  isinstance(utm_zone, np.bytes_) and convert the utm_zone into string
        if isinstance(utm_zone, np.bytes_):
            utm_zone = utm_zone.decode('UTF-8') # b'54J'
        try:
            zone_number = int(utm_zone[0:-1])  #b'54J'
            zone_letter = utm_zone[-1]
        except ValueError:
            raise ValueError('Zone number {0} is not a number'.format(utm_zone[0:-1]))
        is_northern = True if zone_letter.lower() >= 'n' else False
    elif isinstance(utm_zone, int):
        # std UTM code returned by gdal
        is_northern = False if utm_zone < 0 else True
        zone_number = abs(utm_zone)
    else:
        print("epsg and utm_zone", str(epsg), str(utm_zone))

        raise NotImplementedError("utm_zone type (%s, %s) not supported"%(type(utm_zone), str(utm_zone)))
    
    if epsg is None:
        if HAS_GDAL:
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if

    if HAS_GDAL:
        ll_cs = utm_cs.CloneGeogCS()
        utm2ll = osr.CoordinateTransformation(utm_cs, ll_cs).TransformPoint
        ll_point = list(utm2ll(easting, northing, 0.))
    else:
        ll_point = pp(easting, northing, inverse=True)
    # end if

    # be sure to round out the numbers to remove computing with floats
    return round(ll_point[1], 6), round(ll_point[0], 6)


def project_points_ll2utm(lat, lon, datum='WGS84', utm_zone=None, epsg=None):
    """
    Project a list of points that is in Lat, Lon (will be converted to decimal 
    degrees) into UTM coordinates.
    
    Arguments:
    ---------------
        **lat** : float or string (DD:MM:SS.ms)
                  latitude of point
                  
        **lon** : float or string (DD:MM:SS.ms)
                  longitude of point
        
        **datum** : string
                    well known datum ex. WGS84, NAD27, NAD83, etc.

        **utm_zone** : string
                       zone number and 'S' or 'N' e.g. '55S'. Defaults to the
                       centre point of the provided points
                       
        **epsg** : int
                   epsg number defining projection (see 
                   http://spatialreference.org/ref/ for moreinfo)
                   Overrides utm_zone if both are provided

    Returns:
    --------------
        **proj_point**: tuple(easting, northing, zone)
                        projected point in UTM in Datum
                    
    """

    lat = np.array(lat)
    lon = np.array(lon)

    # check length of arrays
    if np.shape(lat) != np.shape(lon):
        raise ValueError("latitude and longitude arrays are of different lengths")

    # flatten, if necessary
    flattened = False
    llshape = np.shape(lat)
    if llshape[0] > 1:
        flattened = True
        lat = lat.flatten()
        lon = lon.flatten()

    '''
    # check lat/lon values
    # this is incredibly slow; disabling for the time being
    for ii in range(len(lat)):
        lat[ii] = assert_lat_value(lat[ii])
        lon[ii] = assert_lon_value(lon[ii])
    '''
    
    if lat is None or lon is None:
        return None, None, None

    if HAS_GDAL:
        # set utm coordinate system
        utm_cs = osr.SpatialReference()
        utm_cs.SetWellKnownGeogCS(datum)

        # set lat, lon coordinate system
        ll_cs = utm_cs.CloneGeogCS()
        ll_cs.ExportToPrettyWkt()
    # end if

    # get zone number, north and zone name
    if epsg is not None:
        if HAS_GDAL:
            # set projection info
            ogrerr = utm_cs.ImportFromEPSG(epsg)
            if ogrerr != OGRERR_NONE:
                raise Exception("GDAL/osgeo ogr error code: {}".format(ogrerr))
            # get utm zone (for information) if applicable
            utm_zone = utm_cs.GetUTMZone()

            # Whilst some projections e.g. Geoscience Australia Lambert (epsg3112) do
            # not yield UTM zones, they provide eastings and northings for the whole
            # Australian region. We therefore set UTM zones, only when a valid UTM zone
            # is available
            if(utm_zone>0):
                # set projection info
                utm_cs.SetUTM(abs(utm_zone), utm_zone > 0)
        else:
            pp = pyproj.Proj('+init=EPSG:%d'%(epsg))
        # end if
    else:
        if utm_zone is not None:
            # get zone number and is_northern from utm_zone string
            zone_number = int(utm_zone[0:-1])
            is_northern = True if utm_zone[-1].lower() > 'n' else False
        else:
            # get centre point and get zone from that
            latc = (np.nanmax(lat) + np.nanmin(lat)) / 2.
            lonc = (np.nanmax(lon) + np.nanmin(lon)) / 2.
            zone_number, is_northern, utm_zone = get_utm_zone(latc, lonc)
        # set projection info

        if(HAS_GDAL):
            utm_cs.SetUTM(zone_number, is_northern)
        else:
            projstring = '+proj=utm +zone=%d +%s +datum=%s' % \
                         (zone_number, 'north' if is_northern else 'south', datum)
            pp = pyproj.Proj(projstring)
        # end if
    # end if
    
    if HAS_GDAL:
        ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoints
        easting, northing, elev = np.array(ll2utm(np.array([lon, lat]).T)).T

    else:
        ll2utm = pp
        easting, northing = ll2utm(lon, lat)
    # end if

    projected_point = (easting, northing, utm_zone)

    # reshape back into original shape
    if flattened:
        lat = lat.reshape(llshape)
        lon = lon.reshape(llshape)

    return projected_point
# end func

# =================================
# functions from latlon_utm_conversion.py


_deg2rad = np.pi / 180.0
_rad2deg = 180.0 / np.pi
_equatorial_radius = 2
_eccentricity_squared = 3

_ellipsoid = [
    #  id, Ellipsoid name, Equatorial Radius, square of eccentricity
    # first once is a placeholder only, To allow array indices to match id
    # numbers
    [-1, "Placeholder", 0, 0],
    [1, "Airy", 6377563, 0.00667054],
    [2, "Australian National", 6378160, 0.006694542],
    [3, "Bessel 1841", 6377397, 0.006674372],
    [4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
    [5, "Clarke 1866", 6378206, 0.006768658],
    [6, "Clarke 1880", 6378249, 0.006803511],
    [7, "Everest", 6377276, 0.006637847],
    [8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
    [9, "Fischer 1968", 6378150, 0.006693422],
    [10, "GRS 1967", 6378160, 0.006694605],
    [11, "GRS 1980", 6378137, 0.00669438],
    [12, "Helmert 1906", 6378200, 0.006693422],
    [13, "Hough", 6378270, 0.00672267],
    [14, "International", 6378388, 0.00672267],
    [15, "Krassovsky", 6378245, 0.006693422],
    [16, "Modified Airy", 6377340, 0.00667054],
    [17, "Modified Everest", 6377304, 0.006637847],
    [18, "Modified Fischer 1960", 6378155, 0.006693422],
    [19, "South American 1969", 6378160, 0.006694542],
    [20, "WGS 60", 6378165, 0.006693422],
    [21, "WGS 66", 6378145, 0.006694542],
    [22, "WGS-72", 6378135, 0.006694318],
    [23, "WGS-84", 6378137, 0.00669438]
]


# Reference ellipsoids derived from Peter H. Dana's website-
# http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
# Department of Geography, University of Texas at Austin
# Internet: pdana@mail.utexas.edu
# 3/22/95

# Source
# Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
# 1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency

@deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_ll2utm() should be "
            "used instead.")
def ll_to_utm(reference_ellipsoid, lat, lon):
    """
    converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees
    Written by Chuck Gantz- chuck.gantz@globalstar.com

    Outputs:
        UTMzone, easting, northing"""

    a = _ellipsoid[reference_ellipsoid][_equatorial_radius]
    ecc_squared = _ellipsoid[reference_ellipsoid][_eccentricity_squared]
    k0 = 0.9996

    # Make sure the longitude is between -180.00 .. 179.9
    long_temp = (lon + 180) - int((lon + 180) / 360) * 360 - 180  # -180.00 .. 179.9

    lat_rad = lat * _deg2rad
    long_rad = long_temp * _deg2rad

    zone_number = int((long_temp + 180) / 6) + 1

    if 56.0 <= lat < 64.0 and 3.0 <= long_temp < 12.0:
        zone_number = 32

    # Special zones for Svalbard
    if 72.0 <= lat < 84.0:
        if 0.0 <= long_temp < 9.0:
            zone_number = 31
        elif 9.0 <= long_temp < 21.0:
            zone_number = 33
        elif 21.0 <= long_temp < 33.0:
            zone_number = 35
        elif 33.0 <= long_temp < 42.0:
            zone_number = 37

    long_origin = (zone_number - 1) * 6 - 180 + 3  # +3 puts origin in middle of zone
    long_origin_rad = long_origin * _deg2rad

    # compute the UTM Zone from the latitude and longitude
    utm_zone = "%d%c" % (zone_number, _utm_letter_designator(lat))

    ecc_prime_squared = ecc_squared / (1 - ecc_squared)
    N = a / np.sqrt(1 - ecc_squared * np.sin(lat_rad) ** 2)
    T = np.tan(lat_rad) ** 2
    C = ecc_prime_squared * np.cos(lat_rad) ** 2
    A = np.cos(lat_rad) * (long_rad - long_origin_rad)

    M = a * (
        (1
         - ecc_squared / 4
         - 3 * ecc_squared ** 2 / 64
         - 5 * ecc_squared ** 3 / 256) * lat_rad
        - (3 * ecc_squared / 8
           + 3 * ecc_squared ** 2 / 32
           + 45 * ecc_squared ** 3 / 1024) * np.sin(2 * lat_rad)
        + (15 * ecc_squared ** 2 / 256
           + 45 * ecc_squared ** 3 / 1024) * np.sin(4 * lat_rad)
        - (35 * ecc_squared ** 3 / 3072) * np.sin(6 * lat_rad))

    utm_easting = (k0 * N * (A
                             + (1 - T + C) * A ** 3 / 6
                             + (5 - 18 * T
                                + T ** 2
                                + 72 * C
                                - 58 * ecc_prime_squared) * A ** 5 / 120)
                   + 500000.0)

    utm_northing = (k0 * (M
                          + N * np.tan(lat_rad) * (A ** 2 / 2
                                                   + (5
                                                      - T
                                                      + 9 * C
                                                      + 4 * C ** 2) * A ** 4 / 24
                                                   + (61
                                                      - 58 * T
                                                      + T ** 2
                                                      + 600 * C
                                                      - 330 * ecc_prime_squared) * A ** 6 / 720)))

    if lat < 0:
        utm_northing = utm_northing + 10000000.0  # 10000000 meter offset for southern hemisphere
    return utm_zone, utm_easting, utm_northing


def _utm_letter_designator(lat):
    # This routine determines the correct UTM letter designator for the given latitude
    # returns 'Z' if latitude is outside the UTM limits of 84N to 80S
    # Written by Chuck Gantz- chuck.gantz@globalstar.com

    if 84 >= lat >= 72:
        return 'X'
    elif 72 > lat >= 64:
        return 'W'
    elif 64 > lat >= 56:
        return 'V'
    elif 56 > lat >= 48:
        return 'U'
    elif 48 > lat >= 40:
        return 'T'
    elif 40 > lat >= 32:
        return 'S'
    elif 32 > lat >= 24:
        return 'R'
    elif 24 > lat >= 16:
        return 'Q'
    elif 16 > lat >= 8:
        return 'P'
    elif 8 > lat >= 0:
        return 'N'
    elif 0 > lat >= -8:
        return 'M'
    elif -8 > lat >= -16:
        return 'L'
    elif -16 > lat >= -24:
        return 'K'
    elif -24 > lat >= -32:
        return 'J'
    elif -32 > lat >= -40:
        return 'H'
    elif -40 > lat >= -48:
        return 'G'
    elif -48 > lat >= -56:
        return 'F'
    elif -56 > lat >= -64:
        return 'E'
    elif -64 > lat >= -72:
        return 'D'
    elif -72 > lat >= -80:
        return 'C'
    else:
        return 'Z'  # if the Latitude is outside the UTM limits


@deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_utm2ll() should be "
            "used instead.")
def utm_to_ll(reference_ellipsoid, northing, easting, zone):
    """
    converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees.
    Written by Chuck Gantz- chuck.gantz@globalstar.com
    Converted to Python by Russ Nelson <nelson@crynwr.com>

    Outputs:
        Lat,Lon
    """

    k0 = 0.9996
    a = _ellipsoid[reference_ellipsoid][_equatorial_radius]
    ecc_squared = _ellipsoid[reference_ellipsoid][_eccentricity_squared]
    e1 = (1 - np.sqrt(1 - ecc_squared)) / (1 + np.sqrt(1 - ecc_squared))
    # NorthernHemisphere; //1 for northern hemispher, 0 for southern

    x = easting - 500000.0  # remove 500,000 meter offset for longitude
    y = northing

    zone_letter = zone[-1]
    zone_number = int(zone[:-1])
    if zone_letter >= 'N':
        NorthernHemisphere = 1  # point is in northern hemisphere
    else:
        NorthernHemisphere = 0  # point is in southern hemisphere
        y -= 10000000.0  # remove 10,000,000 meter offset used for southern hemisphere

    # +3 puts origin in middle of zone
    long_origin = (zone_number - 1) * 6 - 180 + 3

    ecc_prime_squared = ecc_squared / (1 - ecc_squared)

    M = y / k0
    mu = M / (a * (1 - ecc_squared / 4 - 3 * ecc_squared ** 2 /
                   64 - 5 * ecc_squared ** 3 / 256))

    phi1_rad = (mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * np.sin(2 * mu)
                + (21 * e1 ** 2 / 16 - 55 * e1 ** 4 / 32) * np.sin(4 * mu)
                + (151 * e1 ** 3 / 96) * np.sin(6 * mu))
    phi1 = phi1_rad * _rad2deg

    n1 = a / np.sqrt(1 - ecc_squared * np.sin(phi1_rad) ** 2)
    t1 = np.tan(phi1_rad) ** 2
    c1 = ecc_prime_squared * np.cos(phi1_rad) ** 2
    r1 = a * (1 - ecc_squared) / np.power(1 - ecc_squared *
                                          np.sin(phi1_rad) ** 2, 1.5)
    d = x / (n1 * k0)

    lat = phi1_rad - (n1 * np.tan(phi1_rad) / r1) * (
        d ** 2 / 2 - (5 + 3 * t1 + 10 * c1 - 4 * c1 ** 2 - 9 * ecc_prime_squared) * d ** 4 / 24
        + (
            61 + 90 * t1 + 298 * c1 + 45 * t1 ** 2 - 252 * ecc_prime_squared - 3 * c1 ** 2) * d ** 6 / 720)
    lat = lat * _rad2deg

    lon = (d - (1 + 2 * t1 + c1) * d ** 3 / 6 + (
        5 - 2 * c1 + 28 * t1 - 3 * c1 ** 2 + 8 * ecc_prime_squared + 24 * t1 ** 2)
           * d ** 5 / 120) / np.cos(phi1_rad)
    lon = long_origin + lon * _rad2deg
    return lat, lon

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


@gdal_data_check
@deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_utm2ll() should be "
            "used instead.")
def transform_utm_to_ll(easting, northing, zone,
                        reference_ellipsoid='WGS84'):
    utm_coordinate_system = osr.SpatialReference()
    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid)

    try:
        zone_number = int(zone[0:-1])
        zone_letter = zone[-1]
    except ValueError:
        raise ValueError('Zone number {0} is not a number'.format(zone[0:-1]))
    is_northern = True if zone_letter.lower() >= 'n' else False

    utm_coordinate_system.SetUTM(zone_number, is_northern)

    # Clone ONLY the geographic coordinate system
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS()

    # create transform component
    utm_to_ll_geo_transform = osr.CoordinateTransformation(utm_coordinate_system,
                                                           ll_coordinate_system)
    # returns lon, lat, altitude
    return utm_to_ll_geo_transform.TransformPoint(easting, northing, 0)


@gdal_data_check
# @deprecated("This function may be removed in later release. mtpy.utils.gis_tools.project_point_ll2utm() should be "
#             "used instead.")
def transform_ll_to_utm(lon, lat, reference_ellipsoid='WGS84'):
    """
    transform a (lon,lat) to  a UTM coordinate.
    The UTM zone number will be determined by longitude. South-North will be determined by Lat.
    :param lon: degree
    :param lat: degree
    :param reference_ellipsoid:
    :return: utm_coordinate_system, utm_point
    """

    def get_utm_zone(longitude):
        return (int(1 + (longitude + 180.0) / 6.0))

    def is_northern(latitude):
        """
        Determines if given latitude is a northern for UTM
        """
        # if (latitude < 0.0):
        #     return 0
        # else:
        #     return 1
        return latitude >= 0

    utm_coordinate_system = osr.SpatialReference()
    # Set geographic coordinate system to handle lat/lon
    utm_coordinate_system.SetWellKnownGeogCS(reference_ellipsoid)
    utm_coordinate_system.SetUTM(get_utm_zone(lon), is_northern(lat))

    # Clone ONLY the geographic coordinate system
    ll_coordinate_system = utm_coordinate_system.CloneGeogCS()
    # create transform component
    ll_to_utm_geo_transform = osr.CoordinateTransformation(ll_coordinate_system,
                                                           utm_coordinate_system)

    utm_point = ll_to_utm_geo_transform.TransformPoint(lon, lat, 0)

    # returns easting, northing, altitude
    return utm_coordinate_system, utm_point



#################################################################
# Example usages of this script/module
# python gis_tools.py
#=================================================================
if __name__ == "__main__":

    mylat=-35.0
    mylon=149.5
    utm = project_point_ll2utm(mylat, mylon)
    print ("project_point_ll2utm(mylat, mylon) =:  ", utm)

    if HAS_GDAL:
        utm2 = transform_ll_to_utm(mylon, mylat)
        print ("The transform_ll_to_utm(mylon, mylat) results lat, long, elev =: ", utm2[1])

        spref_obj=utm2[0]
        print("The spatial ref string =:", str(spref_obj))
        print("The spatial ref WKT format =:", spref_obj.ExportToWkt())
        print("*********  Details of the spatial reference object ***************")
        print ("spref.GetAttrValue('projcs')", spref_obj.GetAttrValue('projcs'))
        print( "spref.GetUTMZone()", spref_obj.GetUTMZone())
    # end if


