#!/usr/bin/env python

"""
mtpy/processing/filter.py

Functions for the frequency filtering of raw time series. 



For calling a batch filtering rather than just one file, use the appropriate scripts from the mtpy.utils subpackage. 


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys, os
import os.path as op
import time
import copy
import scipy.signal as SS
from scipy.signal import butter, lfilter, buttord



import  mtpy.utils.exceptions as MTex

#=================================================================


def butter_bandpass(lowcut, highcut, samplingrate, order=4):
    nyq = 0.5 * samplingrate
    low = lowcut / nyq
    high = highcut / nyq
    print high, low
    if high >=1. and low == 0.:
        b = np.array([1.])
        a = np.array([1.])

    elif high < 0.95 and low > 0. :
        wp = [1.05*low,high-0.05]
        ws = [0.95*low,high+0.05]
 
        print wp,ws
        order,wn = buttord(wp,ws,0., 30.)
        b, a = butter(order, wn, btype='band')
    
    elif high>= 0.95:
        print 'highpass',low,1.2*low,0.8*low
        order,wn = buttord( 15*low,0.05*low,gpass=0.0, gstop=10.0)
        print order,wn
        b, a = butter(order, wn, btype='high')
    elif low <= 0.05:
        print 'lowpass',high
        order,wn = buttord( high-0.05,high+0.05,gpass=0.0, gstop=10.0)
        b, a = butter(order, wn, btype='low')

    return b, a



def butter_bandpass_filter(data, lowcut, highcut, samplingrate, order=4):
    b, a = butter_bandpass(lowcut, highcut, samplingrate, order=order)
    print b,a
    y = lfilter(b, a, data)
    return y


def tukey(window_length, alpha=0.2):
    '''The Tukey window, also known as the tapered cosine window, can be regarded as a cosine lobe of width alpha * N / 2  that is convolved with a rectangle window of width (1 - alpha / 2). At alpha = 0 it becomes rectangular, and at alpha = 1 it becomes a Hann window.
 
    output
 
    Reference
    ---------
 
    http://www.mathworks.com/access/helpdesk/help/toolbox/signal/tukeywin.html
 
    '''
    # Special cases
    if alpha <= 0:
        return np.ones(window_length) #rectangular window
    elif alpha >= 1:
        return np.hanning(window_length)
 
    # Normal case
    x = np.linspace(0, 1, window_length)
    w = np.ones(x.shape)
 
    # first condition 0 <= x < alpha/2
    first_condition = x<alpha/2
    w[first_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[first_condition] - alpha/2) ))
 
    # second condition already taken care of
 
    # third condition 1 - alpha / 2 <= x <= 1
    third_condition = x>=(1 - alpha/2)
    w[third_condition] = 0.5 * (1 + np.cos(2*np.pi/alpha * (x[third_condition] - 1 + alpha/2)))
 
    return w    

