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
import time 

from mtpy.utils.exceptions import *
import mtpy.utils.format as MTformat
import mtpy.utils.filehandling as FH
reload(FH)
import mtpy.utils.misc as MISC
reload(MISC)


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

    wd = op.abspath(op.realpath(op.join(ts_directory,'birrp_processed')))
    if not op.isdir(wd):
        try:
            os.makedirs(wd) 
        except:
            raise MTpyError_file_handling('cannot create working directory:%s'%(wd))

    os.chdir(wd)

    inputstring, birrp_stationdict = generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold)



    #print inputstring
    #sys.exit()
    #correct inputstring for potential errorneous line endings:
    tempstring = inputstring.split()
    tempstring = [i.strip() for i in tempstring]
    inputstring = '\n'.join(tempstring)

    logfile = open('birrp_logfile.log','w')

    birrpprocess = subprocess.Popen(birrp_exe, stdin=subprocess.PIPE, stdout=logfile,stderr=logfile)

    out,err = birrpprocess.communicate(inputstring)

    logfile.close()

    #generate a local configuration file, containing information about all BIRRP and station parameters
    #required for the header of the EDI file 
    FH.write_dict_to_configfile(birrp_stationdict, 'stationname_birrpconfig.cfg')

    #go back to initial directory
    os.chdir(current_dir)



def generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold, output_channels=2):

    if not output_channels in [2,3]:
        raise MTpyError_inputarguments( 'Output channels must be 2 or 3' )


    (input_filename, length, sampling_rate), birrp_stationdict = set_birrp_input_file_simple(stationname, ts_directory, output_channels, op.join(ts_directory,'birrp_wd'))


    longest_section, number_of_bisections = get_optimal_window_decimation(length, sampling_rate)

    birrp_stationdict['max_window_length'] = longest_section
    birrp_stationdict['n_bisections'] = number_of_bisections

    if output_channels == 2:
        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n%i\n4,3,4\n%s\n0\n0,90,0\n0,90,0\n0,90,0'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,input_filename)

    elif output_channels == 3:
        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n2\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n\i\n4,3,4\n%s\n0\n0,90,0\n0,90,0\n0,90,0'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,stationname)

    string_file = op.join(ts_directory,'birrp_wd','birrp_input_string.txt')
    F = open(string_file,'w')
    F.write(inputstring)
    F.close()

    return inputstring, birrp_stationdict



def set_birrp_input_file_simple(stationname, ts_directory, output_channels, w_directory = '.'):
    """
    File handling: collect longest possible input for BIRRP from different files for the given station name. Generate a new input file in the working directory and return the name of this file, as well as time series properties. 

    Scan all files in the directory by their headers: if the stationname matches, add the data file to the list.
    Additionally read out channels and start-/end times. Find longest consecutive time available on all channels.
    Then generate nx4/5 array for n consecutive time steps. Array in order Ex, Ey, Bx, By (,Bz) saved into file 'birrp_input_data.txt'

    Output_channels determine the number of output channels: 2 for Ex, Ey - 3 for additional Bz
    """


    lo_files = []
    lo_channels = []
    lo_starttimes = []
    lo_endtimes = []
    lo_sampling_rates = []
    birrp_stationdict = {}


    channels =  ['ex', 'ey', 'bx', 'by']
    if output_channels == 3:
        channels.append('bz')

    for entry in os.listdir(ts_directory):
        fn = op.join(ts_directory,entry)
        lo_files.append(fn)

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
        endtime = (np.arange(header[4]+1)/header[2] + header[3])[-1]
        lo_endtimes.append(endtime)

    #take the most common sampling rate, if there are more than one:
    from collections import Counter
    tmp_dummy1 = lo_sampling_rates
    tmp_dummy2 = Counter(tmp_dummy1)
    sampling_rate = tmp_dummy2.most_common(1)[0][0]


    if not len(set(lo_channels)) in [4,5]:
        raise MTpyError_ts_data( 'Missing data files in directory %s - not all channels found'%ts_directory )

    #get a list with all existing time windows of consecutive data for all the channels
    lo_time_windows = []
    #find sorting of the files by their start time:
    starttime_sorting = np.argsort(lo_starttimes)

    #loop over the components:
    for ch in channels:
        tmp_starttime = None
        tmp_endtime = None

        tmp_timewindows_list_per_channel = []
        
        for st in starttime_sorting:
            ch_read = lo_channels[st]
            if not ch == ch_read:
                continue

            if tmp_starttime != None:
                if (tmp_endtime != None) and (np.abs(lo_starttimes[st] - tmp_endtime) > 0.5*1./sampling_rate):
                    tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))
                    tmp_starttime = lo_starttimes[st] 
            
            else:
                tmp_starttime = lo_starttimes[st] 
            tmp_endtime = lo_endtimes[st]
        if tmp_starttime != None:
            tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))


        lo_time_windows.append(tmp_timewindows_list_per_channel)

    longest_common_time_window = MISC.find_longest_common_time_window_from_list(lo_time_windows, sampling_rate)


    # print longest_common_time_window
    # print time.gmtime(longest_common_time_window[0]), time.gmtime(longest_common_time_window[1])
    # os.chdir('/data/TS')
    # sys.exit()


    #data array to hold time series for longest possible time window for the files given 
    #order Ex, Ey, Bx, By (,Bz)
    data = np.zeros((longest_common_time_window[2],output_channels+2))

    #define time axis for referencing time and position in the output array
    #correct by rounding for internal floating point errors
    ta = np.array([ np.round( i , -int(np.log10(1./sampling_rate))) for i in  np.linspace(*longest_common_time_window)])


    for idx_ch, ch in enumerate(channels):
        for st in starttime_sorting:
            ch_read = lo_channels[st]
            if not ch == ch_read:
                continue

            sampling_rate_read = lo_sampling_rates[st]
            if not sampling_rate_read == sampling_rate:
                continue

            #read in data
            data_in = np.loadtxt(lo_files[st])
            #define time axis for read in data
            ta_file = np.arange(len(data_in))/sampling_rate + lo_starttimes[st]
            #find overlap of overall time axis and the ta of current data set
            overlap = np.sort(list(set(ta_file) & set(ta)))

            #find starting index of overlap for current data file time axis
            idx_ta_file = np.argmin(np.abs(ta_file - overlap[0]))

            #find starting index of overlap for overall time axis
            idx_overall_ta = np.argmin(np.abs(ta - overlap[0])) 

            #set data entries
            data[idx_overall_ta:idx_overall_ta+len(overlap), idx_ch] = data_in[idx_ta_file:idx_ta_file+len(overlap)]

            #print 'wrote %i data to array'%len(overlap)

    #define output file for storing the output data array to:
    if not op.isdir(w_directory):
        os.makedirs(w_directory)
        print 'created directory:%s'%w_directory

    print 'size of data arry',data.shape

    try:
        outfn = op.join(w_directory, 'birrp_input_data.txt') 
        np.savetxt(outfn, data)
    except:
        raise MTpyError_file_handling('Error - cannot write data to file:%s'%outfn)

    print 'Wrote input data to file:%s'%outfn


    return (op.abspath(outfn), len(data), sampling_rate), birrp_stationdict


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

def rename_coherencefiles():


    pass


def setup_arguments():

    pass


def convert(stationname, in_dir, configfile, out_dir = None):
    """
    Convert BIRRP output files into .edi and .coh file.

    List of BIRRP output files is searched for in the in_dir directory for the given stationname (base part of filenames). All meta-data must be provided in a config file. If MTpy standard processing has been applied, this is the same file as used from the beginning of the processing. 
    
    If the overall config file is missing, a temporary config file must been created from the header information of the time series files that have been used as BIRRP input.

    The outputfiles 'stationname.edi' and 'stationname.coh' are stored in the out_dir directory. If out_dir is not given, the files are stored in the in_dir.
    """ 

    input_dir = op.abspath(op.realpath(in_dir))
    if not op.isdir(input_dir):
        raise MTpyError_inputarguments('Directory not existing:%s'%(input_dir))

    if out_dir == None:
        output_dir = input_dir
    else:
        output_dir = op.abspath(op.realpath(out_dir))
        if not op.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except:
                print 'output directory could not be created - using input directory instead'
                output_dir = input_dir

    if not op.isfile(configfile):
        raise MTpyError_inputarguments('Config file not existing:%s'%(configfile))
    
    #read the config file:
    try:
        config_dir = FH.read_configfile(configfile)
    except:
        raise EX.MTpyError_config_file( 'Config file cannot be read: %s' % (configfile) )

    if not stationname.upper() in config_dir:
        raise EX.MTpyError_config_file( 'No information about station %s found in config file: %s' % (stationname, configfile) )





    pass

    


