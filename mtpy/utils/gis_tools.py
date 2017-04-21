# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:47:48 2017

@author: jrpeacock
"""

from osgeo import ogr
from osgeo import osr
import mtpy.utils.format as MTft

#==============================================================================
# Make sure lat and lon are in decimal degrees
#==============================================================================
def _assert_minutes(minutes):
    assert minutes >= 0 and minutes < 60.,\
           'minutes needs to be <60 and >0, currently {0:.0f}'.format(minutes)
           
    return minutes
           
def _assert_seconds(seconds):
    assert seconds >= 0 and seconds < 60.,\
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
        >>> import mpty.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_str2float('-118:34:56.3')
    """
    
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
    
    position_value = sign*(abs(deg)+minutes/60.+sec/3600.)
    
    return position_value 
    
def assert_lat_value(latitude):
    """
    make sure latitude is in decimal degrees
    """
    try:
        lat_value = float(latitude)
        
    except ValueError:
        lat_value = convert_position_str2float(latitude)
        
    if abs(lat_value) >= 90:
            raise ValueError('|Latitude| > 90, unacceptable!')
            
    return lat_value
    
def assert_lon_value(longitude):
    """
    make sure longitude is in decimal degrees
    """
    try:
        lon_value = float(longitude)
        
    except ValueError:
        lon_value = convert_position_str2float(longitude)
        
    if abs(lon_value) >= 180:
            raise ValueError('|Longitude| > 180, unacceptable!')
            
    return lon_value

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
        >>> import mpty.utils.gis_tools as gis_tools
        >>> gis_tools.convert_position_float2str(-118.34563)
        
    """
    
    assert type(position) is float, 'Given value is not a float'
    
    deg = int(position)
    sign = 1
    if deg < 0: 
        sign = -1
        
    deg = abs(deg)
    minutes = (abs(position)-deg)*60.
    sec = (minutes-int(minutes))*60.
    
    position_str = '{0}:{1:02.0f}:{2:02f}'.format(sign*deg, int(minutes), sec)
    
    return position_str

#==============================================================================
# Project a point
#==============================================================================

def get_utm_zone(latitude, longitude):
    """
    Get utm zone from a given latitude and longitude
    """
    zone_number = (int(1+(longitude+180.0)/6.0))
    if (latitude < 0.0):
        is_northern = 0
        n_str = 'S'
    else:
        is_northern = 1
        n_str = 'N'
    
    return zone_number, is_northern, '{0:.0f}{1}'.format(zone_number, n_str) 

def project_point_ll2utm(point, datum='WGS84'):
    """
    Project a point that is in Lat, Lon (will be converted to decimal degrees)
    into UTM coordinates.
    
    Arguments:
    ---------------
        **point** : tuple (lat, lon)
                    latitude and longitude of point to be projected.
        
        **datum** : string
                    well known datum ex. WGS84, NAD27, etc.
                    
    Returns:
    --------------
        **proj_point**: tuple(easting, northing, zone)
                        projected point in UTM in Datum
                    
    """
    # make sure the lat and lon are in decimal degrees
    lat = assert_lat_value(point[0])
    lon = assert_lon_value(point[1])
    
    # get zone number, north and zone name
    zone_number, is_northern, utm_zone = get_utm_zone(lat, lon)
    
    ## set utm coordinate system
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS(datum)
    utm_cs.SetUTM(zone_number, is_northern);
       
    ## set lat, lon coordinate system
    wgs84_cs = utm_cs.CloneGeogCS()
    wgs84_cs.ExportToPrettyWkt()
       
    ## set the transform wgs84_to_utm and do the transform
    ll2utm = osr.CoordinateTransformation(wgs84_cs, utm_cs)
    easting, northing, elev = list(ll2utm.TransformPoint(lon, lat))
    projected_point = (easting, northing, utm_zone)    

    return projected_point
       
    
p = (40.0, -118.0)
new_p = project_point_ll2utm(p)

print new_p
#    ## set the transform utm_to_wgs84 and do the transform
#    transform_UTM_To_WGS84 = osr.CoordinateTransformation(utm_cs,wgs84_cs)
#    latlon_points = list(transform_UTM_To_WGS84.TransformPoint(utm_points[0],
#                                                               utm_points[1]))

