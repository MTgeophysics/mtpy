# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 14:47:48 2017

@author: jrpeacock
"""

#==============================================================================
# Imports
#==============================================================================
from osgeo import osr
import numpy as np

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
        
    except TypeError:
        return None
        
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
        print 'WARNING -- {0} is not a number, setting elevation to 0'.format(elevation)
        
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
    if sec == 60:
        minutes += 1
        sec = 0
        
    if minutes == 60:
        deg += 1
        minutes = 0
    
    position_str = '{0}:{1:02.0f}:{2:02.2f}'.format(sign*int(deg), 
                                                    int(minutes),
                                                    float(sec))
    
    return position_str

#==============================================================================
# Project a point
#==============================================================================
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
    zone_number = (int(1+(longitude+180.0)/6.0))
    if (latitude < 0.0):
        is_northern = 0
        n_str = 'S'
    else:
        is_northern = 1
        n_str = 'N'
    
    return zone_number, is_northern, '{0:02.0f}{1}'.format(zone_number, n_str) 

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
    # make sure the lat and lon are in decimal degrees
    lat = assert_lat_value(lat)
    lon = assert_lon_value(lon)
    
    if lat is None or lon is None:
        return None, None, None
    
    # get zone number, north and zone name
    if utm_zone is None:
        zone_number, is_northern, utm_zone = get_utm_zone(lat, lon)
    else:
        # get zone number and is_northern from utm_zone string
        zone_number = int(filter(str.isdigit,utm_zone))
        is_northern = min(1,utm_zone.lower().count('s'))
    
    ## set utm coordinate system
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS(datum)
    
    if type(epsg) is not int:
        utm_cs.SetUTM(zone_number, is_northern);
    else:
        utm_cs.ImportFromEPSG(epsg)
        utm_zone = get_utm_string_from_sr(utm_cs)
       
    ## set lat, lon coordinate system
    ll_cs = utm_cs.CloneGeogCS()
    ll_cs.ExportToPrettyWkt()
       
    ## set the transform wgs84_to_utm and do the transform
    ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs)
    
    ## return different results depending on if lat/lon are iterable
    easting, northing, elev = list(ll2utm.TransformPoint(lon, lat))
    projected_point = (easting, northing, utm_zone)    

    return projected_point
    
    
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
        raise ValueError("easting is not a float")
    try:
        northing = float(northing)
    except ValueError:
        raise ValueError("northing is not a float")

    ## set utm coordinate system
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS(datum)

    if ((utm_zone is None) or (len(utm_zone) == 0) or (utm_zone=='0')):
        if epsg is None:
            raise ValueError('Please provide either utm_zone or epsg')
        else:
            utm_cs.ImportFromEPSG(epsg)
            utm_zone = get_utm_string_from_sr(utm_cs)
    else:
        assert len(utm_zone) == 3, 'UTM zone should be imput as ##N or ##S'
        
        try:
            zone_number = int(utm_zone[0:2])
        except ValueError:
            raise ValueError('Zone number {0} is not a number'.format(utm_zone[0:2]))
        is_northern = 1
        if 's' in utm_zone.lower():
            is_northern = 0
        
        utm_cs.SetUTM(zone_number, is_northern);

    ## set lat, lon coordinate system
    ll_cs = utm_cs.CloneGeogCS()
    ll_cs.ExportToPrettyWkt()
       
    ## set the transform utm to lat lon
    transform_utm2ll = osr.CoordinateTransformation(utm_cs, ll_cs)
    ll_point = list(transform_utm2ll.TransformPoint(easting, northing)) 
    
    # be sure to round out the numbers to remove computing with floats
    return (round(ll_point[1], 6), round(ll_point[0], 6))
  

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

    # check length of arrays
    if np.shape(lat) != np.shape(lon):
        raise ValueError("latitude and longitude arrays are of different lengths")
        
    # flatten, if necessary
    flattened = False
    if np.shape(lat) > 1:
        flattened = True
        llshape = np.shape(lat)
        lat = np.array(lat).flatten()
        lon = np.array(lon).flatten()
      
    # check lat/lon values
    for ii in range(len(lat)):
        lat[ii] = assert_lat_value(lat[ii])
        lon[ii] = assert_lon_value(lon[ii])

    if lat is None or lon is None:
        return None, None, None

    ## set utm coordinate system
    utm_cs = osr.SpatialReference()
    utm_cs.SetWellKnownGeogCS(datum)

    # get zone number, north and zone name
    use_epsg = False
    if utm_zone is None:
        if epsg is None:
            # get centre point and get zone from that
            latc = (np.nanmax(lat) + np.nanmin(lat))/2.
            lonc = (np.nanmax(lon) + np.nanmin(lon))/2.
            zone_number, is_northern, utm_zone = get_utm_zone(latc, lonc)
        else:
            use_epsg = True
    else:
        # get zone number and is_northern from utm_zone string
        zone_number = int(filter(str.isdigit,utm_zone))
        is_northern = min(1,utm_zone.count('S'))
    
    # set projection info
    if use_epsg:
        utm_cs.ImportFromEPSG(epsg)
        # get utm zone (for information) if applicable
        utm_zone = get_utm_string_from_sr(utm_cs)
    else:
        utm_cs.SetUTM(zone_number, is_northern);
        
    ## set lat, lon coordinate system
    ll_cs = utm_cs.CloneGeogCS()
    ll_cs.ExportToPrettyWkt()
       
    ## set the transform wgs84_to_utm and do the transform
    ll2utm = osr.CoordinateTransformation(ll_cs, utm_cs)
    
    ## return different results depending on if lat/lon are iterable
    easting, northing, elev = np.array(ll2utm.TransformPoints(np.array([lon, lat]).T)).T
    projected_point = (easting, northing, utm_zone)    
    
    # reshape back into original shape
    if flattened:
        lat = lat.reshape(*llshape)
        lon = lon.reshape(*llshape)

    return projected_point
    