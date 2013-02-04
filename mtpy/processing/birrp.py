#!/usr/bin/env python

"""
This module contains functions for handling BIRRP software. 
 
- validate_data -- validation of BIRRP input data files
- run -- call the external function BIRRP as system process
- sortoutput -- sort output functions
- validate output -- validate the output functions
- renamecoherencefiles -- rename the coherence output files
- setup_arguments -- calculating and validating the input arguments for the call of BIRRP
- convert -- convert the output of BIRRP into "Z" and "Tipper" values



@UofA, 2013
(LK)

"""

#=================================================================



import numpy as np
import re
import sys, os
import glob
import os.path as op
import subprocess

from mtpy.utils.exceptions import *
import mtpy.utils.format as MTformat
import mtpy.utils.filehandling as FH

#=================================================================


def runbirrp2in2out_simple(birrp_exe, stationname, ts_directory, coherence_threshold = 0.5):
    """
    Call BIRRP for 2 input and 2 output channels with the simplest setup. 

    Provide stationname and directory containing the data folders. Data must be in 1-column ascii format, including one header line for identification of station, channel and sampling rate. For this function, the files must have aligned time series!
    All parameters are automatically determined, the threshold for coherence is set to 0.5 by default.
    Output files are put into a subfolder of the source directory, named 'birrp_processed'.

    """

    if not op.isfile(birrp_exe):
        raise MTpyError_inputarguments('birrp executable not found: %s'%birrp_exe)
    if not op.isdir(ts_directory):
        raise MTpyError_inputarguments('time series files directory: %s'%ts_directory)

    current_dir = op.abspath(os.curdir)

    wd = op.abspath(op.realpath(ts_directory))
    if not op.isdir(wd):
        try:
            os.makedirs(wd) 
        except:
            raise MTpyError_file_handling('cannot create working directory:%s'%(wd))

    os.chdir(wd)

    inputstring = generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold)


    #correct inputstring for potential errorneous line endings:
    tempstring = inputstring.split()
    tempstring = [i.strip() for i in tempstring]
    inputstring = '\n'.join(tempstring)

    logfile = open('birrp_logfile.log','w')

    birrpprocess = subprocess.Popen(birrp_exe, stdin=subprocess.PIPE, stdout=logfile,stderr=logfile)

    out,err = birrpprocess.communicate(inputstring)

    logfile.close()

    os.chdir(current_dir)



def generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold, output_channels=2):

    if not output_channels in [2,3]:
        raise MTpyError_inputarguments( 'Output channels must be 2 or 3' )


    input_filename, length, sampling_rate = set_birrp_input_file_simple(stationname, ts_directory, output_channels)


    longest_section, number_of_bisections = get_optimal_window_decimation(length, sampling_rate)

    
    if output_channels == 2:
        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n\i\n4,3,4\n%s\n0\n0,90,180\n0,90,180\n0,90,180'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,stationname)

    elif output_channels == 3:
        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n2\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n\i\n4,3,4\n%s\n0\n0,90,180\n0,90,180\n0,90,180'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,stationname)


    return inputstring



def set_birrp_input_file_simple(stationname, ts_directory, output_channels):
    """
    File handling: collect longest possible input for BIRRP from different files for the given station name. Generate a new input file in the working directory and return the name of this file, as well as time series properties. 

    Scan all files in the directory by their headers: if the stationname matches, add the data file to the list.
    Additionally read out channels and start-/end times. Find longest consecutive time available on all channels.
    Then generate nx4 array for n consecutive time steps. Array in order Ex, Ey, Bx, By (,Bz) saved into file 'birrp_input_data.txt'
    """

    from scipy.stats import mode


    lo_files = []
    lo_channels = []
    lo_starttimes = []
    lo_endtimes = []
    lo_sampling_rates = []

    channels =  ['ex', 'ey', 'bx', 'by']
    if output_channels == 3:
        channels.append('bz')

    for entry in os.listdir(ts_directory):
        fn = op.join(ts_directory,entry)

        try:
            header = FH.read_data_header(fn)
        except:
            continue


        stationname_read = header[0]
        if not stationname_read == stationname.upper():
            continue

        lo_channels.append(header[1])
        lo_sampling_rates.append(header[2])
        lo_starttimes.append(header[3])
        lo_endtimes.append(header[4]+1./float(header[2]))

    #take the most common sampling rate, if there are more than one:
    sampling_rate = float(mode(lo_sampling_rates)[0])

    #get a list with all existing time windows of consecutive data for all the channels
    lo_time_windows = []
    #find sorting of the files by their start time:
    starttime_sorting = np.argsort(lo_starttimes)
    tmp_starttime = None

    #loop over the components:
    for ch in channels:
        tmp_timewindows_list_per_channel = []
        for st in starttime_sorting:
            ch_read = lo_channels[st]
            if not ch == ch_read:
                continue

            if tmp_starttime != None:
                if np.abs(lo_starttimes[st] - tmp_endtime) > 0.5*1./sampling_rate:
                    tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))
                    tmp_starttime = lo_starttimes[st] 
            else:
                tmp_starttime = lo_starttimes[st] 

            tmp_endtime = lo_endtimes[st]

        lo_time_windows.append(tmp_timewindows_list_per_channel)

    longest_common_time_window = MISC.find_longest_common_time_window_from_list(lo_time_windows, sampling_rate)








def get_optimal_window_decimation(length, sampling_rate):
    """

    input:
    - data series length in samples
    - sampling rate in Hz

    """

    # longest time window cannot exceed 1/30 of the total length in order to obtain statistically sound results (maybe correcting this to 1/100 later); value given in samples
    longest_window = int(length/30.)

    #shortest_window cannot be smaller than 2^4=16 times the sampling interval in order to suit Nyquist and other criteria:
    shortest_window = 16 #int(16.*1./sampling_rate)


    #find maximal number of bisections so that the current window is still longer than the shortest window allowed
    
    number_of_bisections = int(np.ceil(np.log(16./longest_window) / np.log(0.5) ))


    return longest_window, number_of_bisections


def run():

    pass


def validate_data():

    pass

def sortoutput():

    pass

def validate_outputfiles():

    pass

def renamecoherencefiles():


    pass


def setup_arguments():

    pass


def convert():

    pass

    


