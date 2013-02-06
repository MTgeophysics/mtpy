#!/usr/bin/env python

"""
This module contains helper functions for miscellaneous stuff. 

Calculations, optimisations, ....


@UofA, 2013
(LK)

"""
import sys
import numpy as np

def find_longest_common_time_window_from_list(lo_time_windows, sampling_rate):
    """
    Determine the longest consecutive time window, in which all given input traces are present.

    input:
    list of lists - 4(5) lists, each containing 2-tuples with start and end values of time series.

    output:
    3-tuple (start time, end time, n_samples) - end time refers to the start time of the last sample
    """

    #find absolute min and max:
    mins = []
    maxs = []

    for ch_list in lo_time_windows:
        for minmaxtuple in ch_list:
            mins.append(minmaxtuple[0])
            maxs.append(minmaxtuple[1])

    totalmin = np.min(mins)
    totalmax = np.max(maxs)

    #do not correct for last sample, since 'end' in the MTpy handling is defined as the start time of the last sample already
    totallength = (totalmax - totalmin ) * sampling_rate
    

    #define time axis:
    ta = np.arange(totallength)/sampling_rate + totalmin

    #set up array for 4/5 components
    d = np.zeros((totallength,len(lo_time_windows)))

    for idx_ch, ch_list in enumerate(lo_time_windows):
        for minmaxtuple in ch_list:
            s = np.argmin(np.abs(ta-minmaxtuple[0]))
            e = np.argmin(np.abs(ta-minmaxtuple[1]))+1

            #fill array with value 1 where data are available
            d[s:e,idx_ch] = 1

    start_idx = None
    end_idx = None

    t1 = 0
    ts_tmp = None
    te_tmp = None
    longest_window = 0 

    window_idx = 0

    while t1 < totallength:
        if np.prod(d[t1,:]) == 0 :
            #check, if it's been a data window before
            if ts_tmp != None:
                #define sample before as end sample
                te_tmp = t1-1
                #get window length
                window_length = te_tmp - ts_tmp
                #check if it's been the longest window so far
                if window_length > longest_window:
                    longest_window = window_length
                    start_idx = ts_tmp
                    end_idx = te_tmp 
                #re-initialise temporary variables
                ts_tmp = None
                te_tmp = None
            #otherwise just continue to next step in outer loop
            t1 += 1
            continue
        #check if it's the first sample of a data window:
        if ts_tmp == None:
            #if yes, initialise temporary variable
            ts_tmp = t1
            window_idx += 1
        t1 += 1

    #after the loop, check, if last sample belogs to a data window:
    if ts_tmp != None:
        te_tmp = t1 -1
        window_length = te_tmp - ts_tmp
        if window_length > longest_window:
            longest_window = window_length
            start_idx = ts_tmp
            end_idx = te_tmp  

    #rounding limits of the time window to precision defined by the sampling rate
    precision = -int(np.log10(1./sampling_rate))
    return (round(ta[start_idx], precision), round(ta[end_idx], precision), window_length)

def add_birrp_simple_parameters_to_dictionary(birrp_dictionary):

    birrp_dictionary['ilev'] = 0 
    birrp_dictionary['ninp'] = 2
    birrp_dictionary['tbw'] = 2
    birrp_dictionary['uin'] = 0 
    birrp_dictionary['ainuin'] = 0.999
    birrp_dictionary['nlev'] = 0 
    birrp_dictionary['npcs'] = 1 
    birrp_dictionary['nar'] = 5
    birrp_dictionary['imode'] = 2
    birrp_dictionary['jmode'] = 0 
    birrp_dictionary['nfil'] = 0
    birrp_dictionary['nskip'] = 0 
    birrp_dictionary['theta1'] = 0 
    birrp_dictionary['theta2'] = 90
    birrp_dictionary['phi'] = 0
 
    return birrp_dictionary
