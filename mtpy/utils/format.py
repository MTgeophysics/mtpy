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

import mtpy.utils.exceptions as MTex


def _assert_position_format(coordinate, value): 

    """
    Check, if value is valid for the three coordinates of a position: 'latitude, longitude, elevation' - raise 'MTpyError_config_file'-exception, if not.

    If lat/lon are given in deg,min,sec it is converted do degrees. The value is returned, if no exception was raised 

    """
    if coordinate in ['ele','elev','elevation']:

        try:
            elev = float(value)
        except: 
            raise MTex.MTpyError_config_file('Config file error: elevation value is not a number')
        
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
                if len(latlon_list) == 3:
                    try:
                        latlon_list = [float(i) for  i in latlon_list]
                    except:
                        raise 
                    latlon = convert_dms_tuple2degrees(latlon_list)

                elif len(latlon_list) == 2:
                    #TODO:
                    #humbug...just add zero to the pair and then use the standard!
                    try:
                        latlon_list = [float(i) for  i in latlon_list]
                    except:
                        raise 
                    latlon = convert_degmin_tuple2degrees(latlon_list)
                
            except:
                raise MTpyError_config_file('Config file error: lat/lon is in invalid format')


        if coordinate in ['latitude','lat'] and ( not -90 <= latlon <= 90):
            raise MTex.MTpyError_config_file('Error - Latitude out of range')

        if coordinate in ['longitude','lon','long'] and ( not -180 <= latlon <= 180):
            raise MTex.MTpyError_config_file('Error - Longitude out of range')

        value = latlon


    return value


def assert_decimal_coordinates(coord):

    try:
        #if it's in decimal degrees already:
        dec_coord = float(coord)
        return dec_coord

    except:
        #see, if it's tuple

        dms_list = re.split('[ :,]', coord)
        if len(dms_list) == 3:
            dec_coord = convert_dms_tuple2degrees(dms_list)
        elif len(dms_list) == 2:
            dec_coord = convert_degmin_tuple2degrees(dms_list)

    return dec_coord


def convert_dms_tuple2degrees(latlon_triple):
    """
    Convert a triple (list, tuple, array) of degrees, minuts, seconds into degrees.

    Validity of the triple is assumed and has to be asserted in advanced.
    """

     
    sign = 1.
    try:
        latlon_list = [float(i) for i in latlon_triple]

    except:
        #if triple is given as string:
        latlon_raw = latlon_triple
        #allow separation by :,. or space for (deg min sec) format
        try:
            latlon_list = re.split('[ :,]', latlon_raw)
        except:
            raise MTex.MTpyError_config_file('Config file error: lat/lon is in invalid format')
        try:
            latlon_list = [float(i) for  i in latlon_list]
        except:
            raise MTex.MTpyError_config_file('Config file error: lat/lon is in invalid format')
    
    if str(latlon_list[0])[0]=='-':
        sign = -1.
    deg = latlon_list[0]

    #take out sign for easier conversion into degrees
    
    if deg < 0:
        deg *= -1

    degrees = deg + 1/60. * latlon_list[1] + 1/3600.* latlon_list[2]


    return degrees * sign



def convert_degrees2dms_tuple(degrees):
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




def convert_degmin_tuple2degrees(latlon_list):

    """
    Convert a 2-tuple (list, tuple, array) of form "(degrees,minutes)" into degrees.

    Validity of the triple is assumed and has to be asserted in advance.
    """

    sign = 1.

    try:
        latlon_list = [float(i) for i in latlon_list]
    except:
        raise

    deg = latlon_list[0]
    minutes = latlon_list[1]

    #take out sign for easier conversion into degrees
    if str(latlon_list[0])[0]=='-':
        sign = -1.
    
    if deg < 0:
        deg *= -1

    degrees = deg + 1/60. * minutes

    return degrees * sign


def convert_degrees2degmin_tuple(degrees):
    """
    Convert a geographical degree value into a 2-tuple (array) of form "(degrees, minutes)"
    """

    deg = float(degrees)

    #take out sign for easier conversion:
    sign = np.sign(deg)
    if deg < 0:
        deg *= -1

    d = int(deg)

    minutes = 60 * (deg - d)

    dms_tuple = np.zeros((2))


    dms_tuple[0] = sign * d
    dms_tuple[1] = minutes

    return dms_tuple
