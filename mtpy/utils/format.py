#!/usr/bin/env python

"""
This module contains functions for formatting of files.

The various functions deal with re-formatting e.g. data, file-structure, file-name.
Additionally, it contains (private) functions, which check correct formats, 
e.g. lat/lon or datetimes. 



@UofA, 2013
(LK)

"""
import numpy as numpy
import  os
import sys
import re

def _assert_position_format(coordinate,value): 

    """
    Check, if value is valid for the three coordinates of a position: latitude, longitude, elevation raise MTpyError_config_file - exception, if not.

    If lat/lon are given in deg,min,sec it is converted do degrees. The value is returned, if no exception was raised 

    """

    if coordinate == 'elevation':
        try:
            elev = float(value)
        except: 
            raise MTpyError_config_file('Config file error: elevation value is not a number')
        
        value = elev


    if coordinate in ['latitude' , 'longitude']:
        
        #check, if latlon is given in degrees
        try:
            latlon = float(value)
        #if not, check, if it's in (deg,min,sec) format
        except: 
            try:
                latlon_raw = value
                #allow separation by :,. or space for (deg min sec) format
                latlon_list = re.split('[ :,]', latlon_raw)
                try:
                    latlon_list = [float(i) for  i in latlon_list]
                except:
                    raise MTpyError_config_file('Config file error: lat/lon is in invalid format')

                latlon = convert_degminsec_tuple2degrees(latlon_list)
                

 

            raise MTpyError_config_file('Config file error: lat/lon is in invalid format')


        if coordinate == 'latitude' and ( not -90 <= latlon <= 90):
            raise MTpyError_config_file('Error - Latitude out of range')

        if coordinate == 'longitude' and ( not -180 <= latlon <= 180):
            raise MTpyError_config_file('Error - Longitude out of range')

        value = latlon


    return value


def convert_degminsec_tuple2degrees(latlon_triple):


    latlon_list = list(latlon_triple)

    #take out sign for easier conversion into degrees
    sign = latlon_list[0]/np.abs(latlon_list[0])
    if latlon_list[0] < 0:
        latlon_list[0] *= -1

.
.
.

    return degrees