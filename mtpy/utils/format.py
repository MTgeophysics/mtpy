#!/usr/bin/env python

"""
This module contains functions for formatting of files.

The various functions deal with re-formatting e.g. data, file-structure, file-name.
Additionally, it contains (private) functions, which check correct formats, 
e.g. lat/lon or datetimes. 



@UofA, 2013
(LK)

"""
import numpy as np
import  os
import sys
import re

from mtpy.utils.exceptions import *


def _assert_position_format(coordinate, value): 

    """
    Check, if value is valid for the three coordinates of a position: 'latitude, longitude, elevation' - raise 'MTpyError_config_file'-exception, if not.

    If lat/lon are given in deg,min,sec it is converted do degrees. The value is returned, if no exception was raised 

    """
    if coordinate in ['ele','elev','elevation']:

        try:
            elev = float(value)
        except: 
            raise MTpyError_config_file('Config file error: elevation value is not a number')
        
        value = elev


    if coordinate in ['latitude' , 'longitude','lat','lon','long']:
        
        #check, if latlon is given in degrees
        try:
            latlon = float(value)
        #if not, check, if it's in (deg,min,sec) format
        except: 
            try:
                latlon_raw = str(value)
                #allow separation by ':','.' or 'space' for (deg min sec) format
                latlon_list = re.split('[ :,]', latlon_raw)
                try:
                    latlon_list = [float(i) for  i in latlon_list]
                except:
                    raise MTpyError_config_file('Config file error: lat/lon is in invalid format')

                latlon = convert_dms_tuple2degrees(latlon_list)
                
            except:
                raise MTpyError_config_file('Config file error: lat/lon is in invalid format')


        if coordinate in ['latitude','lat'] and ( not -90 <= latlon <= 90):
            raise MTpyError_config_file('Error - Latitude out of range')

        if coordinate in ['longitude','lon','long'] and ( not -180 <= latlon <= 180):
            raise MTpyError_config_file('Error - Longitude out of range')

        value = latlon


    return value


def convert_dms_tuple2degrees(latlon_triple):
    """
    Convert a triple (list, tuple, array) of degrees, minuts, seconds into degrees.

    Validity of the triple is assumed and has to be asserted in advanced.
    """


    try:
        latlon_list = [float(i) for i in latlon_triple]
    except:
        #if triple is given as string:
        latlon_raw = latlon_triple
        #allow separation by :,. or space for (deg min sec) format
        try:
            latlon_list = re.split('[ :,]', latlon_raw)
        except:
            raise MTpyError_config_file('Config file error: lat/lon is in invalid format')
        try:
            latlon_list = [float(i) for  i in latlon_list]
        except:
            raise MTpyError_config_file('Config file error: lat/lon is in invalid format')

    deg = latlon_list[0]

    #take out sign for easier conversion into degrees
    sign = deg/np.abs(deg)
    if deg < 0:
        deg *= -1

    degrees = deg + 1/60. * latlon_list[1] + 1/3600.* latlon_list[2]


    return degrees * sign



def convert_degrees2_dms_tuple(degrees):
    """
    Convert a geographical degree value into a triple (array) of degrees, minutes, seconds.
    """

    deg = float(degrees)

    #take out sign for easier conversion:
    sign = deg/np.abs(deg)
    if deg < 0:
        deg *= -1

    d = int(deg)

    minutes = 60 * (deg - d)

    m = int(minutes)

    seconds = 60 * (minutes - m)

    dms_triple = np.zeros((3))


    dms_triple[0] = sign * d
    dms_triple[1] = m
    dms_triple[2] = seconds


    return dms_triple