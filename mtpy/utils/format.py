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
import os
import sys
import re

import mtpy.utils.exceptions as MTex



def _assert_position_format(coordinate, value): 

    """
    Check, if value is valid for the three coordinates of a position: 
    'latitude, longitude, elevation' - raise 'MTpyError_config_file'-exception if not.

    If lat/lon are given in deg,min,sec it is converted into degrees. 
    The value is returned if no exception was raised 

    """
    if coordinate in ['ele','elev','elevation']:

        try:
            elev = float(value)
        except: 
            print 'WARNING - elevation value is not a number - set to 0 (zero)'
            elev = 0.
        
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
                raise MTex.MTpyError_value('Config file error: lat/lon is in invalid format')


        if coordinate in ['latitude','lat'] and ( not -90 <= latlon <= 90):
            raise MTex.MTpyError_value('Error - Latitude out of range')

        if coordinate in ['longitude','lon','long'] and ( not -180 <= latlon <= 360):
            raise MTex.MTpyError_value('Error - Longitude out of range')

        value = latlon


    return value

def assert_decimal_coordinates(coordinate_string):
    """
    Assert that coordinate component is given in decimal degrees. Converts into it otherwise.

    input:
    - string containing coordinte component (degrees or dms triple)
      * separation symbol for dms triple can be :
        - blank space
        - :
        - , 

    output:
    - input value expressed in decimal degrees

    """
    

    cs = coordinate_string.strip()
    lo_cardinals = re.findall('[neswNESW]', cs)
    if len(lo_cardinals)>0:
        #remove cardinal symbol 
        cs = cs.replace(lo_cardinals[0],'')
        if lo_cardinals[0] in ['s','S','e','E']:
            if cs.startswith('-'):
                cs = cs[1:]
            else:
                cs = '-'+cs

    #try to split by the pre-determined separators
    latlon_list = re.split('[ :,]',cs)
    try:
        latlon_list = [float(i) for i in latlon_list]
    except:
        print 'coordinate value format invalid - non numerical values: {0}'.format(cs)

    if len(latlon_list) == 1:
        #only degrees
        try:
            value = float(cs)
            return value
        except:
            print 'coordinate value format invalid - not float: {0}'.format(cs)
            return str(cs)

    elif len(latlon_list) > 3 :
        print 'coordinate value format invalid - too many components: {0}'.format(cs)
        return str(cs)

    if len(latlon_list) == 3:
        return convert_dms_tuple2degrees(latlon_list)
    else:
        return convert_degmin_tuple2degrees(latlon_list)

def convert_dms_string2tuple(coordinatestring):

    dms = [0,0,0]
    try:
        dms_raw = coordinatestring.split(':')
        dms[0] = float(dms_raw[0])
        dms[1] = float(dms_raw[1])
        dms[2] = float(dms_raw[2])
        return dms

    except:
        print 'coordinate string invalid - must be of form "deg:min:sec.ss"'
        return


def convert_dms_tuple2string(dmstuple):

    d = dmstuple[0]
    m = dmstuple[1]
    s = int(round(dmstuple[2],5))
    try:
        dms_string = "%i:%02i:%02i"%(d,m,int(s))
        remainder = s - int(s)
        remainder = round(remainder,5)
        if remainder !=0:
            print remainder
            print round(remainder,5)*1e4
            dms_string += '.%i'%(int(round(remainder,5)*1e4))
        return dms_string

    except:
        print 'coordinate tuple invalid - must be of form "(deg,min,sec)"'
        return


def convert_dms_tuple2degrees(latlon_list):

    """
    Convert a triple (list, tuple, array) of degrees, minuts, seconds into degrees.

    Validity of the triple is assumed and has to be asserted in advance.
    """

    sign = 1.
    try:
            latlon_list = [float(i) for i in latlon_list]
            if str(latlon_list[0]).startswith('-'):
                sign = -1.

    except:

        #if triple is given as string:
        latlon_raw = latlon_list
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
    minutes = latlon_list[1]
    seconds = latlon_list[2]
    if not (-180<=deg<=360 and  0<=minutes<60 and 0<=seconds<60):
        print 'Value for degrees, minutes, or seconds invalid'
        raise MTex.MTpyError_inputarguments()

    #take out sign for easier conversion into degrees

    if deg < 0:
        deg *= -1

    degrees = deg + 1/60. * minutes + 1/3600.* seconds

    return degrees * sign



def convert_degrees2dms_tuple(degrees):
    """
    Convert a geographical degree value into a triple (array) of degrees, minutes, seconds.
    """

    deg = float(degrees)
    dms_triple = np.zeros((3))

    if deg == 0 :
        return dms_triple
        
    #take out sign for easier conversion:
    sign = deg/np.abs(deg)
    if deg < 0:
        deg *= -1

    d = int(deg)

    minutes = 60 * (deg - d)

    m = int(minutes)

    seconds = 60. * (minutes - m)


    dms_triple[0] = sign * d
    dms_triple[1] = m
    dms_triple[2] = seconds


    return dms_triple




def convert_degmin_tuple2degrees(latlon_list):

    """
    Convert a 2tuple (list, tuple, array) of form "(degrees, minutes)" into degrees.

    Validity of the triple is assumed and has to be asserted in advance.
    """

    sign = 1.

    try:
        latlon_list = [float(i) for i in latlon_list]
        if str(latlon_list[0]).startswith('-'):
            sign = -1.
    except:
        raise

    deg = latlon_list[0]
    minutes = latlon_list[1]
    if not (0<=minutes<60):
        raise MTex.MTpyError_inputarguments('Minutes value invalid')
    #take out sign for easier conversion into degrees

    if deg < 0:
        deg *= -1.


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
