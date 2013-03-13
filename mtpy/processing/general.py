#!/usr/bin/env python

"""
mtpy/processing/general.py

Helper functions for standard time series processing. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np

import math, cmath

import mtpy.utils.exceptions as MTexceptions


#=================================================================



def correct4sensor_orientation(x_values, y_values, x_sensor_angle = 0 , y_sensor_angle = 90):
    """
        Correct time series of a sensor pair, which has not been in default (x=0, y=90) orientation.

        Input:
        - x-values - Numpy array
        - y-values - Numpy array
            (same length for both!)

        Optional:
        - Angle of the x-sensor - measured in degrees, clockwise from North (0) 
        - Angle of the y-sensor - measured in degrees, clockwise from North (0) 

        Output:
        - corrected x-values 
        - corrected y-values 
    """

    x_values = np.array(x_values)
    y_values = np.array(y_values)

    try:
        if x_values.dtype not in ['complex', 'float', 'int']:
            raise
        if len(x_values) != len(y_values):
            raise
    except:
        raise MTexceptions.MTpyError_inputarguments('ERROR - both input arrays must be of same length')

    in_array = np.zeros((len(x_values), 2), x_values.dtype)

    in_array[:,0] = x_values
    in_array[:,1] = y_values

    try:
        x_angle = math.radians(x_sensor_angle)
        y_angle = math.radians(y_sensor_angle)
    except:
        raise MTexceptions.MTpyError_inputarguments('ERROR - both angles must be of type int or float')
       

    T = np.matrix( [[ np.real(cmath.rect(1,x_angle)), np.imag(cmath.rect(1,x_angle))],[np.real(cmath.rect(1,y_angle)), np.imag(cmath.rect(1,y_angle))]])

    try:
        new_array = np.dot(in_array, T.I)
    except:
        raise MTexceptions.MTpyError_inputarguments('ERROR - angles must define independent axes to span 2D')


    return new_array[0,:], new_array[1,:]