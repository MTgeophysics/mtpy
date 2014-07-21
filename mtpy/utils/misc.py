#!/usr/bin/env python

"""
This module contains helper functions for miscellaneous stuff. 

Calculations, optimisations, ....


@UofA, 2013
(LK)

"""
import sys
import os 
import numpy as np
import ctypes

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
    totallength = int((totalmax - totalmin ) * sampling_rate + 1)
    

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
    print '\t\tMaximum time window covered by data files:', totalmax-totalmin 

    print '\t\tCheck data availablility - while-loop until "maximum time window" is reached...'
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
        if t1%(int(totallength/100.)) == 0:
            sys.stdout.write('\t\t{0:3} %\r'.format(int(np.round(t1/float(totallength) *100))))
            sys.stdout.flush()

    
    print '\n\t\t...Done!'#' \n\t\tChecking for last sample (include/exclude)...'
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
    #print 'return time window parameters:'
    #print ta[start_idx],ta[end_idx] ,window_length, len(ta)
    #print '\t\tStart time, end time, samples: ',round(ta[start_idx], precision),\
    #                     round(ta[end_idx], precision), window_length#, len(ta))
    

    window_length = (round(ta[end_idx], precision) - round(ta[start_idx], precision)) *sampling_rate
    
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



class MemoryCheck():
    """Checks memory of a given system"""
 
    def __init__(self):
 
        if os.name == "posix":
            self.value = self.linuxRam()
        elif os.name == "nt":
            self.value = self.windowsRam()
        else:
            print "I only work with Win or Linux "
 
    def windowsRam(self):
        """Uses Windows API to check RAM in this OS"""
        kernel32 = ctypes.windll.kernel32
        c_ulong = ctypes.c_ulong
        class MEMORYSTATUS(ctypes.Structure):
            _fields_ = [
                ("dwLength", c_ulong),
                ("dwMemoryLoad", c_ulong),
                ("dwTotalPhys", c_ulong),
                ("dwAvailPhys", c_ulong),
                ("dwTotalPageFile", c_ulong),
                ("dwAvailPageFile", c_ulong),
                ("dwTotalVirtual", c_ulong),
                ("dwAvailVirtual", c_ulong)
            ]
        memoryStatus = MEMORYSTATUS()
        memoryStatus.dwLength = ctypes.sizeof(MEMORYSTATUS)
        kernel32.GlobalMemoryStatus(ctypes.byref(memoryStatus))
 
        return int(memoryStatus.dwTotalPhys/1024**2),int(memoryStatus.dwAvailPhys/1024**2)
 
    def linuxRam(self):
        """Returns the RAM of a linux system"""
        totalmemory = os.popen("free -m").readlines()[1].split()[1]
        freememory = os.popen("free -m").readlines()[1].split()[3]
        return int(totalmemory),int(freememory)

def show_memory():

    M = MemoryCheck()
    return M.value
