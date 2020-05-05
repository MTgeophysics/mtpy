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
    
    :type position_str: string ('DD:MM:SS.ms')
    :param position_str: degrees of latitude or longitude
        
    :rtype: float
    :return: latitude or longitude in decimal degrees
                          
    :Example: ::

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
    # need to round seconds to 4 decimal places otherwise machine precision
    # keeps the 60 second roll over and the string is incorrect.
    sec = np.round((minutes - int(minutes)) * 60., 4)
    if sec >= 60.:
        minutes += 1
        sec = 0

    if int(minutes) == 60:
        deg += 1
        minutes = 0
        
    position_str = '{0}:{1:02.0f}:{2:05.2f}'.format(sign * int(deg),
                                                    int(minutes),
                                                    sec)

    return position_str

# ==============================================================================
# Project a point
# ==============================================================================
def get_utm_zone(latitude, longitude):
    """
    Get utm zone from a given latitude and longitude
    """
    zone_number = (int(1 + (longitude + 180.0) / 6.0))
    is_northern = bool(latitude >= 0)
    n_str = _utm_letter_designator(latitude)

    return zone_number, is_northern, '{0:02.0f}{1}'.format(zone_number, n_str)

def utm_zone_to_epsg(zone_number, is_northern):
    """
    get epsg code (WGS84 datum) for a given utm zone
    
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
                
def _utm_letter_designator(lat):
    # This routine determines the correct UTM letter designator for the 
    # given latitude returns 'Z' if latitude is outside the UTM limits of 
    # 84N to 80S Written by Chuck Gantz- chuck.gantz@globalstar.com

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
                
def split_utm_zone(utm_zone):
    """
    
    :param utm_zone: DESCRIPTION
    :type utm_zone: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

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
    
def get_epsg(latitude, longitude):
    """
    get epsg code for the utm projection (WGS84 datum) of a given latitude
    and longitude pair
    """  
    zone_number, is_northern, utm_str = get_utm_zone(latitude, longitude)
    
    return utm_zone_to_epsg(zone_number, is_northern)

def get_gdal_coordinate_system(datum):
    """
    Get coordinate function from GDAL
    """
    # set lat lon coordinate system
    ll_cs = osr.SpatialReference()
    if isinstance(datum, int):
        ogrerr = ll_cs.ImportFromEPSG(datum)
        if ogrerr != OGRERR_NONE:
            raise GISError("GDAL/osgeo ogr error code: {}".format(ogrerr))
    elif isinstance(datum, str):
        ogrerr = ll_cs.SetWellKnownGeogCS(datum)
        if ogrerr != OGRERR_NONE:
            raise GISError("GDAL/osgeo ogr error code: {}".format(ogrerr))
    else:
        raise GISError("""datum {0} not understood, needs to be EPSG as int
                           or a well known datum as a string""".format(datum))

    return ll_cs

def validate_epsg(epsg):
    """
    Make sure epsg is an integer

    :param epsg: DESCRIPTION
    :type epsg: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

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
    
    :param utm_zone: DESCRIPTION
    :type utm_zone: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    #JP: if its unicode then its already a valid string in python 3
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
    
    :param values: DESCRIPTION
    :type values: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    if isinstance(values, (int, float)):
        values = np.array([values], dtype=np.float)
    elif isinstance(values, (list, tuple)):
        values = np.array(values, dtype=np.float)
    elif isinstance(values, str):
        values = [ss.strip() for ss in values.strip().split(',')]
        values = np.array(values, dtype=np.float)
    elif isinstance(values, np.ndarray):
        values = values.astype(np.float)
        
    if location_type in ['lat', 'latitude']:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lat_value(value)
            except GISError as error:
                raise GISError('{0}\n Bad input value at index {1}'.format(
                               error, ii))
            
    if location_type in ['lon', 'longitude']:
        for ii, value in enumerate(values):
            try:
                values[ii] = assert_lon_value(value)
            except GISError as error:
                raise GISError('{0}\n Bad input value at index {1}'.format(
                               error, ii)) 
            
    return values

def _get_gdal_projection_ll2utm(datum, utm_zone, epsg):
    """
    project point using GDAL
    """
    ll_cs = get_gdal_coordinate_system(datum)
    
    # project point on to EPSG coordinate system if given
    if isinstance(epsg, int):
        utm_cs = get_gdal_coordinate_system(epsg)
    # otherwise project onto given datum
    elif epsg is None:
        utm_cs = get_gdal_coordinate_system(datum)

        zone_number, is_northern = split_utm_zone(utm_zone)
        utm_cs.SetUTM(zone_number, is_northern)

    return osr.CoordinateTransformation(ll_cs, utm_cs).TransformPoint

def _get_gdal_projection_utm2ll(datum, utm_zone, epsg):
    """
    
    :param datum: DESCRIPTION
    :type datum: TYPE
    :param utm_zone: DESCRIPTION
    :type utm_zone: TYPE
    :param epsg: DESCRIPTION
    :type epsg: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    zone_number, is_northern = split_utm_zone(utm_zone)
    
    if epsg is not None:
        utm_cs = get_gdal_coordinate_system(epsg)
    else:
        utm_cs = get_gdal_coordinate_system(datum)
        utm_cs.SetUTM(zone_number, is_northern)
    
    ll_cs = utm_cs.CloneGeogCS()
    return osr.CoordinateTransformation(utm_cs, ll_cs).TransformPoint
        

def _get_pyproj_projection(datum, utm_zone, epsg):
    """
    
    :param datum: DESCRIPTION
    :type datum: TYPE
    :param utm_zone: DESCRIPTION
    :type utm_zone: TYPE
    :param epsg: DESCRIPTION, defaults to None
    :type epsg: TYPE, optional
    :return: DESCRIPTION
    :rtype: TYPE

    """

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
#=================================================================
if __name__ == "__main__":

    my_lat= -35.0
    my_lon= 149.5
    utm_point = project_point_ll2utm(my_lat, my_lon)
    print ("project_point_ll2utm(mylat, mylon) =:  ", utm_point)


