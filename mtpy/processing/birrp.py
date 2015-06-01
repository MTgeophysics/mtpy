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

import gc
import numpy as np
import sys
import os
import os.path as op
import subprocess
import time
import datetime 
import fnmatch
import math
import scipy.signal as SS

import mtpy.utils.exceptions as MTex
import mtpy.utils.format as MTft
import mtpy.utils.filehandling as MTfh
import mtpy.utils.configfile as MTcf
import mtpy.utils.misc as MTmc
import mtpy.utils.interpolation as MTip

#=================================================================
#for time stamp differences:
epsilon = 1e-5
#=================================================================

def runbirrp_Nin2out_simple(birrp_exe, stationname, ts_directory, 
                           coherence_threshold=0.0, rr_station=None,
                           output_channels=2, output_dir=None, 
                           starttime=None, endtime=None):


    """
    Call BIRRP for N input and 2 output channels with the simplemost setup. 

    Provide stationname and directory containing the data folders. Data must 
    be in 1-column ascii format, including one header line for identification 
    of station, channel and sampling rate. For this function, the files must 
    have aligned time series!
    
    All parameters are automatically determined, the threshold for coherence 
    is set to 0.5 by default.
    
    If no output directory is specified, output files are put into a subfolder 
    of the source directory, named 'birrp_processed'.

    Optionally, starttime and endtime can be given (in epoch seconds). If that
    yields an interval which is smaller than the provided data, only that 
    section is used.

    Additionally, a configuration file is created. It contains information 
    about the processing paramters for the station. Keys are generic for the 
    common parameters and named after BIRRP input keywords for the processing 
    parameters.

    """

    if not op.isfile(birrp_exe):
        print '\n Error - Birrp executable not found: {0}'.format(birrp_exe)
        raise MTex.MTpyError_inputarguments()

    if not op.isdir(ts_directory):
        print '\n Error - time series files directory not existing: {0}'.format(ts_directory)
        raise MTex.MTpyError_inputarguments()

    current_dir = op.abspath(os.curdir)

    wd = op.abspath(op.realpath(op.join(current_dir,'birrp_wd')))
    if rr_station is not None:
        wd = op.abspath(op.realpath(op.join(current_dir,'birrp_wd_rr')))

    if output_dir != None:
        output_dir = op.abspath(op.join(current_dir,output_dir))
        if not op.isdir(output_dir) :
            try:
                os.makedirs(output_dir)
                wd = output_dir
            except:
                print '\nWarning - Could not find or generate specified output '+\
                      'directory {0} - using default instead!'.format(output_dir)
        else:
            wd = output_dir 

    if not op.isdir(wd):
        try:
            os.makedirs(wd) 
        except:
            print '\nError - cannot create working directory: {0}'.format(wd)
            raise MTex.MTpyError_file_handling()
    
    print "\nFound/generated output directory: {0}".format(wd)

    os.chdir(wd)
    inputstring, birrp_stationdict, inputfilename = generate_birrp_inputstring_simple(
                                            stationname, rr_station, ts_directory, 
                                            coherence_threshold,output_channels, starttime, endtime)

    print "Inputstring and configuration dictionary generated for station {0}\n".format(stationname)
    #print inputstring
    #sys.exit()
    #correct inputstring for potential errorneous line endings due to strange operating systems:
    tempstring = inputstring.split()
    tempstring = [i.strip() for i in tempstring]
    inputstring = '\n'.join(tempstring)
    inputstring += '\n'

    #print 'opening logfile...'
    
    # dummy1 = sys.stdout
    # dummy2 = sys.stderr
    logfile = open('birrp_logfile.log','w')
    # sys.stdout = logfile
    # sys.stderr =  logfile

    print 'Start Birrp processing...'
    #os.system("{0} < {1}".format(birrp_exe,inputfilename))

    birrpprocess = subprocess.Popen(birrp_exe, stdin=subprocess.PIPE, stdout=logfile,stderr=logfile)
    #instringhandler = StringIO.StringIO(inputstring)
    #birrpprocess = subprocess.Popen(birrp_exe, stdin=instringhandler, stdout=logfile,stderr=logfile)

    out,err = birrpprocess.communicate(inputstring)
    
    #sys.stdout = dummy1
    #sys.stderr = dummy2
    logfile.close()

    print '...Done!\nLogfile closed: {0}\n'.format(op.abspath(logfile.name))

    #generate a local configuration file, containing information about all BIRRP and station parameters
    #required for the header of the EDI file 
    # print out
    # print err
    print 'Generating configuration file containing the applied processing parameters...'
    station_config_file = '{0}_birrpconfig.cfg'.format(stationname)
    MTcf.write_dict_to_configfile(birrp_stationdict, station_config_file)
    print '...Done!'
    print 'Wrote BIRRP and time series configurations to file: {0}'.format(op.abspath(station_config_file))

    #go back to initial directory
    os.chdir(current_dir)

    print '\n \t\tDONE !!!\n'



def generate_birrp_inputstring_simple(stationname, rr_station, ts_directory, 
                                        coherence_threshold, output_channels=2,
                                        starttime=None, endtime=None):

    if not output_channels in [2,3]:
        raise MTex.MTpyError_inputarguments( 'Output channels must be 2 or 3' )


    print '\nSetting basic input components,e.g. filenames, samplingrate,...\n'
    input_filename, length, sampling_rate, birrp_stationdict = set_birrp_input_file_simple(
                                    stationname, rr_station, ts_directory, output_channels, 
                                    os.curdir, starttime, endtime)

    print '...Done!\n\nCalculating optimal time window bisection parameters...'
    longest_section, number_of_bisections = get_optimal_window_bisection(length, sampling_rate)
    #BIRRP automatically divides by 4 ... for what ever reason....
    longest_section *= 4
    print '...Done!\n'

    birrp_stationdict['max_window_length'] = longest_section
    birrp_stationdict['n_bisections'] = number_of_bisections
    birrp_stationdict['coherence_threshold'] = coherence_threshold

    #for self referencing:
    birrp_stationdict['rr_station'] = birrp_stationdict['station']
    #otherwise:
    if rr_station is not None:
        birrp_stationdict['rr_station'] = rr_station.upper()

    birrp_stationdict = MTmc.add_birrp_simple_parameters_to_dictionary(birrp_stationdict)


    if output_channels == 2:
        birrp_stationdict['nout'] = 2

        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n%s\n0\n1\n3\n2\n0'\
                    '\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n%i\n4,3,4\n%s\n0'\
                    '\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,longest_section,
                        number_of_bisections,coherence_threshold,stationname,
                        input_filename,length,input_filename)
        if rr_station is not None:
            inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n%s\n0\n1\n3\n'\
                        '2\n0\n0\n0\n0\n0\n0\n0\n6,1,2,3,4\n%s\n0\n%i\n6,5,6\n%s'\
                        '\n0\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,
                            longest_section,number_of_bisections,
                            coherence_threshold,stationname,input_filename,
                            length,input_filename)

    #TODO: check order of columns!!!
    elif output_channels == 3:
        birrp_stationdict['nout'] = 3
        birrp_stationdict['nz'] = 2
        inputstring = '0\n3\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%.2f\n2\n%s\n0\n1\n0\n2\n0'\
                    '\n0\n0\n0\n0\n0\n0\n0\n5,3,4,5,3,4\n%s\n0\n%i\n5,3,4\n%s\n0'\
                    '\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,longest_section,
                        number_of_bisections,coherence_threshold,stationname,
                        input_filename,length,input_filename)
        if rr_station is not None:
            inputstring = '0\n3\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%.2f\n2\n%s\n0\n1\n0'\
                        '\n2\n0\n0\n0\n0\n0\n0\n0\n0\n7,1,2,5,3,4\n%s\n0\n%i\n'\
                        '7,6,7\n%s\n0\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,
                            longest_section,number_of_bisections,
                            coherence_threshold,stationname,input_filename,
                            length,input_filename)



    string_file = 'birrp_string.txt'
    string_file = MTfh.make_unique_filename(string_file)
    with open(string_file,'w') as F:
        F.write(inputstring)
        F.write('\n')
    

    return inputstring, birrp_stationdict,string_file



def set_birrp_input_file_simple(stationname, rr_station, ts_directory, 
                                output_channels = 2, w_directory = '.', 
                                starttime=None, endtime=None):
    """
    File handling: collect longest possible input for BIRRP from different files 
    for the given station name. Cut out section, if start and/or end times are given.
    Generate a new input file in the working directory and return the name of 
    this file, as well as time series and processing properties in form of a 
    dictionary. 

    Scan all files in the directory by their headers: if the stationname matches, 
    add the data file to the list. Additionally read out channels and start-/end 
    times. Find longest consecutive time available on all channels.
    Include remote reference station data, if rr_station is given.
    Then generate n x 4/5 array for n consecutive time steps. Array in order 
    Ex, Ey, Bx, By (,Bz) saved into file 'birrp_data.txt'

    If remote reference is include, generate n x 6/7 array in order 
    Ex, Ey, Bx, By (,Bz), remoteBx, remoteBy

    Output_channels determine the number of output channels: 
    "2" for Ex, Ey - "3" for additional Bz

    Output:
    - filename for birrp input data_in
    - length of data (samples)
    - sampling_rate (in Hz)
    - configuration dictionary for processed station

    """


    lo_station_files = []
    lo_station_channels = []
    lo_station_starttimes = []
    lo_station_endtimes = []
    
    lo_sampling_rates = []

    lo_rr_files = []
    lo_rr_channels = []
    lo_rr_starttimes = []
    lo_rr_endtimes = []
   

    #check for data from (primary) station:

    station_channels =  ['ex', 'ey', 'bx', 'by']
    rr_channels = ['bx', 'by']

    if output_channels == 3:
        station_channels.append('bz')

    for entry in os.listdir(ts_directory):
        fn = op.join(ts_directory,entry)
        rr_station_flag = False 
        
        if not op.isfile(fn):
            continue
        try:
            header = MTfh.read_ts_header(fn)
        except:
            continue
        stationname_read = header['station'].upper()

        if not stationname_read == stationname.upper():
            if rr_station is None:
                continue
            else:
                if not stationname_read == rr_station.upper():
                    continue

        if rr_station is not None:
            if stationname_read == rr_station.upper():
                rr_station_flag = True


        if rr_station_flag is True:
            if header['channel'].lower() in rr_channels:
                
                lo_rr_channels.append(header['channel'].lower())
                lo_rr_starttimes.append(np.float64(header['t_min']))
                ta_rr = np.arange(int(float(header['nsamples']))+1)/float(
                                    header['samplingrate']) + np.float64(header['t_min'])
                ta_rr_endtime = ta_rr[-1]
                lo_rr_endtimes.append(ta_rr_endtime)
                lo_rr_files.append(fn)
                
                if stationname.upper() != rr_station.upper():
                    continue
            

        if not header['channel'].lower() in station_channels:
            continue
        if not stationname_read == stationname.upper():
            continue


        lo_station_files.append(fn)
                
        lo_station_channels.append(header['channel'].lower())
        lo_sampling_rates.append(float(header['samplingrate']))
        lo_station_starttimes.append(np.float64(header['t_min']))
        ta_station = np.arange(int(float(header['nsamples']))+1)/float(header['samplingrate']) + float(header['t_min'])
        ta_station_endtime = ta_station[-1]
        lo_station_endtimes.append(ta_station_endtime)

 
    if (len(lo_sampling_rates) == 0) or (len(lo_station_files) == 0):
        print '\n\tERROR - no MTpy data files for station'\
            ' {1} found in directory {0} !!\n'.format(ts_directory,stationname)
        raise MTex.MTpyError_ts_data()
    
    if rr_station is not None:
        if len(lo_rr_files) == 0:
            print '\n\tERROR - no MTpy data files for '\
                        'remote reference station {1} found in directory {0} '\
                                        '!!\n'.format(ts_directory,rr_station)
            raise MTex.MTpyError_ts_data()

    try:
    #take the most common sampling rate, if there are more than one
    #assuming typo in this case!!
        from collections import Counter
        tmp_dummy1 = lo_sampling_rates
        tmp_dummy2 = Counter(tmp_dummy1)
        sampling_rate = tmp_dummy2.most_common(1)[0][0]
        del Counter
    except:
        #for older Python versions :-( ... 
        sampling_rate = lo_sampling_rates[0]

    if not len(set(lo_station_channels)) in [4,5]:
        print 'Error - Missing data files in directory {0} - not all channels found\n'.format(ts_directory)
        sys.exit( )
    if not len(set(lo_rr_channels)) in [2]:
        if rr_station is not None:
            print 'Error - Missing data files in directory {0} - not all remote channels'\
                        ' found\n'.format(ts_directory)
            sys.exit( )


    #get a list with all existing time windows of consecutive data for all the channels
    lo_time_windows = []

    #find sorting of the files by their start time:
    station_starttime_sorting = np.argsort(lo_station_starttimes)
    rr_starttime_sorting = np.argsort(lo_rr_starttimes)

    print '\tlooping over all components to find time windows with data...'
    #loop over the components:
    for ch in station_channels:
        tmp_starttime = None
        tmp_endtime = None

        tmp_timewindows_list_per_channel = []
        
        for sst in station_starttime_sorting:
            ch_read = lo_station_channels[sst]
            if not ch == ch_read:
                continue

            if tmp_starttime != None:
                if (tmp_endtime != None) and (np.abs(lo_station_starttimes[sst] - tmp_endtime) > 0.5*1./sampling_rate):
                    tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))
                    tmp_starttime = lo_station_starttimes[sst] 
            
            else:
                tmp_starttime = lo_station_starttimes[sst] 
            tmp_endtime = lo_station_endtimes[sst]
        if tmp_starttime != None:
            tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))

        lo_time_windows.append(tmp_timewindows_list_per_channel)

    if rr_station is not None:
        #loop over the remote reference time windows as well 
        for ch in rr_channels:
            tmp_starttime = None
            tmp_endtime = None

            tmp_timewindows_list_per_channel = []
            
            for rst in rr_starttime_sorting:
                ch_read = lo_rr_channels[rst]
                if not ch == ch_read:
                    continue

                if tmp_starttime != None:
                    if (tmp_endtime != None) and (np.abs(lo_rr_starttimes[rst] - tmp_endtime) > 0.5*1./sampling_rate):
                        tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))
                        tmp_starttime = lo_rr_starttimes[rst] 
                
                else:
                    tmp_starttime = lo_rr_starttimes[rst] 
                tmp_endtime = lo_rr_endtimes[rst]
            if tmp_starttime != None:
                tmp_timewindows_list_per_channel.append((tmp_starttime, tmp_endtime))

            lo_time_windows.append(tmp_timewindows_list_per_channel)

    print '\t...Done!\n\tFind longest common time window for all channels...'
    longest_common_time_window = MTmc.find_longest_common_time_window_from_list(lo_time_windows, sampling_rate)
    print '\t...Done:{0}'.format(longest_common_time_window)

    sampling_interval = (longest_common_time_window[1] - longest_common_time_window[0]) / longest_common_time_window[2]

    
    
    try:
        if starttime is None:
            raise
        t_start_given = np.float64(starttime)
        if (t_start_given <  longest_common_time_window[0]) and np.abs(t_start_given - longest_common_time_window[0]) > epsilon:
            print 'Warning - given start time is too small - using data start time'
        if t_start_given >=  longest_common_time_window[1]  and np.abs(t_start_given - longest_common_time_window[1]) > epsilon:
            print 'Warning - provided starttime {0} too large'\
                ' - data set ends already at time {1}'.format(t_start_given,
                                                longest_common_time_window[1])
            raise
        else:
            t_start = max(longest_common_time_window[0],t_start_given)
    except:
        print 'Warning - given start time {0} could not be processed - '\
                'using start time of data instead: {1} '.format(starttime,
                                                longest_common_time_window[0])
        t_start = longest_common_time_window[0]


    try:
        if endtime is None:
            raise
        t_end_given = np.float64(endtime)
        if t_end_given >  longest_common_time_window[1] and np.abs(t_end_given - longest_common_time_window[1]) > epsilon:
            print 'Warning - given end time is too large - using data end time' 
        if t_end_given <=  longest_common_time_window[0] and np.abs(t_end_given - longest_common_time_window[0]) > epsilon:
            print 'Warning - provided end time {0} too small'\
                ' - data set does not start until {1}'.format(endtime,
                                                longest_common_time_window[0])
            raise
        else:
            t_end = min(longest_common_time_window[1],t_end_given)
    except:
        print 'Warning - given end time {0} could not be processed - '\
                'using end time of data instead: {1} '.format(endtime,
                                                longest_common_time_window[1])

        t_end = longest_common_time_window[1]
    

    t_start  = np.round(t_start,5)
    t_end    = np.round(t_end,5)


    #define time axis for referencing time and position in the output array
    #correct by rounding for internal floating point errors
    #ta = np.array([ np.round( i , -int(np.log10(1./sampling_rate))) for i in  np.linspace(*longest_common_time_window)])
    #alternative : no rounding
    ta_full_dataset = np.linspace(*longest_common_time_window, endpoint=False)
 
    idx_start = np.abs(ta_full_dataset-t_start).argmin()
    if t_start > ta_full_dataset[idx_start]:
        idx_start += 1

    #careful to EXCLUDE last sample for t_end...time series ends before a sample!:
    
    #find index for value closest to end time on the time axis
    #reduce by 1 for the sampling step (excluding last sample)  
    idx_end = np.abs(ta_full_dataset - t_end).argmin()  

    #check, if index yields a time value, which is at least one sampling before the end time
    if (t_end - sampling_interval) < ta_full_dataset[idx_end]:
        #check, if variation is outside numerical uncertainty
        if np.abs((t_end - sampling_interval) - ta_full_dataset[idx_end]) > epsilon:
            #if so, reduce index one more to be safely inside the given interval
            idx_end -= 1

    #new time axis - potentially just a sub-section of the original,
    #maximal the same size though
    ta =  ta_full_dataset[idx_start: idx_end + 1]
    
    #print ta[-1]-ta[0],sampling_interval
    
    print '\n\tTime section set to {0} - {1} ({2} samples)'.format(ta[0],
                                            ta[-1]+sampling_interval,len(ta))     
    #sys.exit()

    #print ta[0]-ta_full_dataset[0], ta[-1]-ta_full_dataset[0], len(ta)
    
    #data array to hold time series for longest possible time window for the files given 
    #order Ex, Ey, Bx, By (,Bz)
    #if remote reference is set: Ex, Ey, Bx, By (,Bz), remoteBx, remoteBy
    totalmemory, freememory = MTmc.show_memory()
    print '\n\tTotal memory available: {0} MB - free: {1} MB'.format(totalmemory,freememory)

    data = np.zeros((len(ta),output_channels+2))

    if rr_station is not None:
        data = np.zeros((len(ta),output_channels+4))

    print '\tSize of data array: {0} MB\n'.format(np.round(data.nbytes/1024.**2,2))

    print '\tData array (size: {0}) and time axis (length: {1}) initialised - '\
            'looping over channels to read in data...\n'.format(data.shape, len(ta))
    for idx_ch, ch in enumerate(station_channels):

        #arrayindex = 0
        
        for st in station_starttime_sorting:
            ch_read = lo_station_channels[st]
            if not ch == ch_read:
                continue

            #sampling_rate_read = lo_sampling_rates[st]
            #if not sampling_rate_read == sampling_rate:
            #    continue

            header = MTfh.read_ts_header(lo_station_files[st])
            
            #determine the time axis of the file:
            ta_file = np.arange(int(float(header['nsamples']))) / \
                                    float(header['samplingrate']) + \
                                    float(header['t_min'])
            
            #skip file, if file ends before the chosen section starts
            if ta[0] > ta_file[-1]:
                continue
            #finished station, if file starts after end of chosen section
            if ta_file[0] >= ta[-1]:
                break

            #read in data
            print '\t...reading station data from file {0}'.format(lo_station_files[st])
            data_in = np.loadtxt(lo_station_files[st])
            
            
            if ta_file[0] >= ta[0]: 
                startindex = np.abs(ta - ta_file[0]).argmin()
                if ta_file[-1] <= ta[-1]:
                    #find starting index w.r.t. overall time axis:
                    data[startindex:startindex+len(data_in),idx_ch] = data_in
                else:
                    endindex = np.abs(ta - ta_file[-1]).argmin()
                    data_section_length = endindex - startindex + 1
                    data[startindex:,idx_ch] = data_in[:data_section_length]
            else:
                startindex =  np.abs(ta[0] - ta_file).argmin()
                if ta_file[-1] >= ta[-1]:
                    data_section_length = len(ta) 
                    data[:,idx_ch] = data_in[startindex:startindex+data_section_length]
                else:
                    endindex = np.abs(ta - ta_file[-1]).argmin() +1 
                    data_section_length = len(data[:endindex])
                    data[:endindex,idx_ch] = data_in[startindex:]



            #print len(data_in), sampling_rate , lo_starttimes[st]

            
            print 'file time section: ',ta_file[0],ta_file[-1], '(overall window: ', ta[0],ta[-1],')'

            gc.collect()

    if rr_station is not None:
        #same as before, just with the remote reference data files 
        for idx_ch, ch in enumerate(rr_channels):
            for st in rr_starttime_sorting:
                ch_read = lo_rr_channels[st]
                if not ch == ch_read:
                    continue

                
                header = MTfh.read_ts_header(lo_rr_files[st])
                #determine the time axis of the file:
                ta_file = np.arange(int(float(header['nsamples']))) / \
                                    float(header['samplingrate']) + \
                                    float(header['t_min'])

                #skip file, if file ends before the chosen section starts
                if ta[0] > ta_file[-1]:
                    continue
                #finished station, if file starts after end of chosen section
                if ta_file[0] >= ta[-1]:
                    break
            
                #read in data
                print '\t...reading remote data from file {0}'.format(lo_rr_files[st])
                data_in = np.loadtxt(lo_rr_files[st])
                
                if ta_file[0] >= ta[0]: 
                    startindex = np.abs(ta - ta_file[0]).argmin()
                    if ta_file[-1] <= ta[-1]:
                        #find starting index w.r.t. overall time axis:
                        data[startindex:startindex+len(data_in),idx_ch+2+output_channels] = data_in
                    else:
                        endindex = np.abs(ta - ta_file[-1]).argmin()
                        data_section_length = endindex - startindex + 1
                        data[startindex:,idx_ch+2+output_channels] = data_in[:data_section_length]
                else:
                    startindex =  np.abs(ta[0] - ta_file).argmin()
                    if ta_file[-1] >= ta[-1]:
                        data_section_length = len(ta) 
                        data[:,idx_ch+2+output_channels] = data_in[startindex:startindex+data_section_length]
                    else:
                        endindex = np.abs(ta - ta_file[-1]).argmin() +1 
                        data_section_length = len(data[:endindex])
                        data[:endindex,idx_ch+2+output_channels] = data_in[startindex:]
     
                print 'file time section: ',ta_file[0],ta_file[-1], '(overall window: ', ta[0],ta[-1],')'



    #define output file for storing the output data array to:
    w_directory = op.abspath(op.join(os.curdir, w_directory))
    if not op.isdir(w_directory):
        os.makedirs(w_directory)
        print '\t(Created temporary working directory: {0})'.format(w_directory)

    #print '\tSize of usable data arry: ',data.shape
    
    # for all channels wihtin the 'data'-array:
    # linear detrending as first order pre-whitening filter for BIRRP

    for i in range(data.shape[1]):
        data[:,i] = SS.detrend(data[:,i])

    try:
        outfn = op.join(w_directory, 'birrp_data.txt') 
        outfn = MTfh.make_unique_filename(outfn)
        print '\n\tSave input data array to file: {0} ...'.format(outfn)
        np.savetxt(outfn, data)
    except:
        raise MTex.MTpyError_file_handling('Error - cannot write data to file:{0}'.format(outfn))

    print '\t...Done!\n'

    birrp_stationdict = {}
    birrp_stationdict['station'] =  stationname.upper()
    birrp_stationdict['n_output_channels'] = output_channels
    birrp_stationdict['sampling_rate'] = sampling_rate
    birrp_stationdict['n_samples'] =  len(data)
    birrp_stationdict['processing_window_start'] = ta[0]
    birrp_stationdict['processing_window_end'] = ta[-1] + sampling_interval
    birrp_stationdict['processing_window'] = ta[-1] + sampling_interval - ta[0]
    birrp_stationdict['recording_starttime'] =  lo_station_starttimes[0]

    gc.collect()

    return op.abspath(outfn), birrp_stationdict['n_samples'] , sampling_rate, birrp_stationdict


def get_optimal_window_bisection(length, sampling_rate):
    """

    input:
    - data series length in samples
    - sampling rate in Hz

    """

    # longest time window cannot exceed 1/30 of the total length in order to obtain 
    #statistically sound results (maybe correcting this to 1/100 later); value given in samples
    longest_window = int(length/30.)

    #restrict lentght to power of 2...optimisinf the FFT and memory allocation:
    longest_window = 2**(int(np.log2(longest_window)))

    #shortest_window cannot be smaller than 2^4=16 times the sampling interval 
    #in order to suit Nyquist and other criteria; value given in samples
    shortest_window = 16 


    #find maximal number of bisections so that the current window is still longer than the shortest window allowed
    number_of_bisections = int(np.ceil(np.log(float(shortest_window)/longest_window) / np.log(0.5) ))

    return longest_window, number_of_bisections



def write_script_file(processing_dict, save_path=None):
    """
    writeScriptfile(processingdict will write a script file for BIRRP using 
    info in processingdict which is a dictionary with keys:
    
    ================== ========================================================
    parameter          description
    ================== ======================================================== 
    station            station name
    fn_list             list of file names to be processed, this must be in 
                       the correct order [EX, EY, HZ, HX, HY] and if multiple
                       sections are to be processed at the same time then 
                       must be input as a nested loop 
                       [[EX1, EY1, HZ1, HX1, HY1], 
                       [EX2, EY2, HZ2, HX2, HY2], ...]
    rrfn_list           list of remote reference file names, similar to the 
                       fn_list [[HX1, HY1], [HX2, HY2], ...]
    ilev               processing mode 0 for basic and 1 for advanced RR-2 
                       stage
    nout               Number of Output time series (2 or 3-> for BZ)
    ninp               Number of input time series for E-field (1,2,3) 
    nref               Number of reference channels (2 for MT)
    nrr                bounded remote reference (0) or 2 stage bounded 
                       influence (1)
    tbw                Time bandwidth for Sepian sequence
    deltat             Sampling rate (+) for (s), (-) for (Hz)
    nfft               Length of FFT (should be even)
    nsctinc            section increment divisor (2 to divide by half)
    nsctmax            Number of windows used in FFT
    nf1                1st frequency to extract from FFT window (>=3)
    nfinc              frequency extraction increment 
    nfsect             number of frequencies to extract
    mfft               AR filter factor, window divisor (2 for half)
    uin                Quantile factor determination
    ainlin             Residual rejection factor low end (usually 0)
    ainuin             Residual rejection factor high end (.95-.99)
    c2threshb          Coherence threshold for magnetics (0 if undesired)
    c2threshe          Coherence threshold for electrics (0 if undesired)
    nz                 Threshold for Bz (0=separate from E, 1=E threshold, 
                                         2=E and B) 
                       Input if 3 B components else None
    c2thresh1          Squared coherence for Bz, input if NZ=0, Nout=3
    perlo              longest period to apply coherence threshold over
    perhi              shortes period to apply coherence threshold over
    ofil               Output file root(usually three letters, can add full
                                        path)
    nlev               Output files (0=Z; 1=Z,qq; 2=Z,qq,w; 3=Z,qq,w,d)
    nprej              number of frequencies to reject
    prej               frequencies to reject (+) for period, (-) for frequency
    npcs               Number of independent data to be processed (1 for one 
                       segement)
    nar                Prewhitening Filter (3< >15) or 0 if not desired',
    imode              Output file mode (0=ascii; 1=binary; 2=headerless ascii; 
                       3=ascii in TS mode',
    jmode              input file mode (0=user defined; 1=sconvert2tart time 
                                        YYYY-MM-DD HH:MM:SS)',
    nread              Number of points to be read for each data set  
                       (if segments>1 -> npts1,npts2...)',
    nfil               Filter parameters (0=none; >0=input parameters; 
                                          <0=filename)
    nskip              Skip number of points in time series (0) if no skip, 
                        (if segements >1 -> input1,input2...)',
    nskipr             Number of points to skip over (0) if none,
                       (if segements >1 -> input1,input2...)',
    thetae             Rotation angles for electrics (relative to geomagnetic 
                       North)(N,E,rot)',
    thetab             Rotation angles for magnetics (relative to geomagnetic 
                       North)(N,E,rot)',
    thetar             Rotation angles for calculation (relative to geomagnetic 
                       North)(N,E,rot)'
    hx_cal             full path to calibration file for hx 
    hy_cal             full path to calibration file for hy 
    hz_cal             full path to calibration file for hz 
    rrhx_cal           full path to calibration file for remote reference hx 
    rrhy_cal           full path to calibration file for remote reference hy 
                       Note that BIRRP assumes the calibrations are the same
                       as in the first segment, so the remote reference must
                       be the same for all segments, or at least the same
                       instruments have to be used so the calibration is 
                       the same.
    ================== ========================================================              
    
    .. see also::                
                
        => see BIRRP Manual for more details on the parameters
        => see A. D. Chave and D. J. Thomson [1989,2003,2004] for more
            information on Bounded influence and robust processing.
            
    Arguments:
    -----------
        **processing_dict** : dictionary with keys as above
        
        **save_path** : string (full path to directory to save script file)
                        if none saves as:
                            os.path.join(os.path.dirname(fn_list[0]),'BF')
        
    Outputs:
    --------
        **script_file** : full path to script file to guide birrp
        
        **birrp_dict** : dictionary of birrp parameters input into script file
        
    
    """
     
    #===================================================================
    # Write a script file for BIRRP, Chave et al. [2004]            
    #===================================================================
    #print processingdict
    #compute how many timeseries and days there are 
    #ndf = # of files per day
    #nds = # of day
    pdict = dict(processing_dict)
    
    try:
        fn_array = np.array(pdict['fn_list'])
    except KeyError:
        raise KeyError('fn_list --> Need to input a list of files to process')
    
    try:
        nds, ndf = fn_array.shape
    except ValueError:
        ndf = fn_array.shape[0]
        nds = 0
    if save_path is None:
        if nds == 0:
            bfpath = os.path.join(os.path.dirname(pdict['fn_list'][0]),
                                 'BF')
        else:
            bfpath = os.path.join(os.path.dirname(pdict['fn_list'][0][0]),
                                  'BF')
    else:
        bfpath = save_path
    
    
    if nds == 0:
        npcs=1
    elif nds==1:
        nds=0
        npcs=1
        pdict['fn_list'] = pdict['fn_list'][0]
        try:
            pdict['rrfn_list'] = pdict['rrfn_list'][0]  
        except KeyError:
            pass

    else:
        npcs=int(nds)
        
    #make a directory to put BIRRP Files (BF)
    if not os.path.exists(bfpath):
        os.mkdir(bfpath)
        print 'Made directory: ', bfpath

    #output file stem, full path
    ofil = os.path.join(bfpath,pdict['station']) 
    
    #mode to process: default is basic
    ilev = int(pdict.pop('ilev', 0))
    
    #number of output channels
    try:
        nout = int(pdict['nout'])
    except KeyError:
        if nds!=0:
            if ndf==5:
                nout = 3
            else:
                nout = 2
        elif ndf==5:
            nout = 3
        else:
            nout = 2
    
    #number of input channels default is 2
    ninp = int(pdict.pop('ninp', 2))
        
    #time bandwidth window size            
    tbw = int(pdict.pop('tbw', 2))
        
    #------------Options for Advanced mode-------------------------
    if ilev == 1:
        #Advanced: number of remote reference channels
        nref = int(pdict.pop('nref', 2))
            
        #Advanced: remote reference type processing
        nrr = int(pdict.pop('nrr', 1))          
            
        #Advanced: magnetic coherence threshold
        c2threshb = float(pdict.pop('c2threshb', 0))
            
        #Advanced: window increment divisor
        nsctinc = int(pdict.pop('nsctinc', 2))
            
        #Advanced: first frequency to extract
        nf1 = int(pdict.pop('nf1', tbw+2))
        
        #Advanced: frequency increment
        nfinc = int(pdict.pop('nfinc', tbw))
            
        #number of frequencies to extract
        nfsect = int(pdict.pop('nfsec', 2))
            
        #number AR filter is divided by
        mfft = int(pdict.pop('mfft', 2))
        
        #Advanced: lower bound of leverage point rejection
        ainlin = float(pdict.pop('ainlin', .0001))
        
        #Advanced: coherence threshold low period
        perlo = int(pdict.pop('perlo', 1000))
            
        #Advanced: coherenct threshold high period
        perhi = float(pdict.pop('perhi', .0001))
            
        #Advanced:  number of frequencies to reject
        nprej = int(pdict.pop('nprej', 0))
        
        #Advanced
        try:
            prej = pdict['prej'].split(',')
            if type(prej) is list:
                prej = [float(ff) for ff in prej]
            if nprej != len(prej):
                nprej = len(prej)
        except KeyError:
            prej = []
            
    #---------------------Options for Basic Mode--------------------------
    
    #time series sampling rate
    deltat = int(pdict.pop('deltat', -100))
    
    #max length of fft window 
    nfft = int(pdict.pop('nfft', 2**16))

    #maximum number of sections
    nsctmax = int(pdict.pop('nsctmax', 12))
    
    #quantile factor
    uin = int(pdict.pop('uin', 0))            
    
    #upper bound of leverage point rejection
    ainuin = float(pdict.pop('ainuin', .9999))
 
    #electric channel coherence threshold
    c2threshe = int(pdict.pop('c2threshe', 0))

    #Bz coherency threshold mode
    try:
        nz = int(pdict['nz'])
    except KeyError:
        if nout == 3:
            nz = 0
        else:
            nz = None
    
    #Bz coherence threshold
    try:
        c2threshe1 = float(pdict['c2threshe1'])
    except KeyError:
        if nout == 3:
            c2threshe1 = 0
        else:
            c2threshe1 = None
    
    #output level
    nlev = int(pdict.pop('nlev', 0))
    
    #order of prewhitening auto regressive filter
    nar = int(pdict.pop('nar', 5))
        
    #input mode
    imode = int(pdict.pop('imode', 0))
    
    #output mode
    jmode = int(pdict.pop('jmode', 0))       
    
    #name of filter file
    nfil = int(pdict.pop('nfil', 0))
    
    #calibration for coils
    hx_cal = pdict.pop('hx_cal', None) 
    hy_cal = pdict.pop('hy_cal', None) 
    hz_cal = pdict.pop('hz_cal', None)
    
    rrhx_cal = pdict.pop('rrhx_cal', None) 
    rrhy_cal = pdict.pop('rrhy_cal', None) 
    
    if jmode == 0:
        #number of points to read
        nread = pdict.pop('nread', 1440000)
        if type(nread) is not list and type(nread) is not np.ndarray:
            nread = int(nread)
        #number of points to skip in time series
        try:
            nskip = pdict['nskip']
            if nds != 0 and type(nskip) is not list and \
               type(nskip) is not np.ndarray:
                nskip = [nskip for ii in range(nds)]
        except KeyError:
            if nds != 0:
                nskip = [0 for ii in range(nds)]
            else:
                nskip = 0
        
        #number of point to skip from remote reference time series    
        try:
            nskipr = pdict['nskipr']
            if nds != 0 and type(nskipr) is not list and \
               type(nskipr) is not np.ndarray:
                nskipr = [nskipr for ii in range(nds)]
        except KeyError:
            if nds==0:
                nskipr = 0
            else:
                nskipr = [0 for ii in range(nds)]
    
    if jmode == 1:
        #start time of data
        dstim = pdict.pop('dstim', '1970-01-01 00:00:00')
        
        #window start time
        wstim = pdict.pop('wstim', '1970-01-01 00:00:00')
        
        #window end time
        wetim = pdict.pop('wetim', '1970-01-02 00:00:00')
        
    
    #rotation angle of electric channels
    thetae = pdict.pop('thetae', '0,90,0')
    
    #rotation angle of magnetic channels
    thetab = pdict.pop('thetab', '0,90,0')
    
    #rotation angle of final impedance tensor
    thetaf = pdict.pop('thetaf', '0,90,0')

    #===================================================================
    # Write values to a .script file
    #===================================================================
    #print '+++ ',nskipr
    #print ndf,nds
    #write to a file
    scriptfile=ofil+'.script'
    fid=file(scriptfile,'w')
    if ilev==0: 
        fid.write('{0:d} \n'.format(ilev))
        fid.write('{0:d} \n'.format(nout))
        fid.write('{0:d} \n'.format(ninp))
        fid.write('{0:.3f} \n'.format(tbw))
        fid.write('{0:.3f} \n'.format(deltat))
        fid.write('{0:d},{1:d} \n'.format(nfft,nsctmax))
        fid.write('y \n')
        fid.write('{0:.5f},{1:.5f} \n'.format(uin,ainuin))
        fid.write('{0:.3f} \n'.format(c2threshe))
        #parameters for bz component if ninp=3
        if nout==3:
            if c2threshe==0:
                fid.write('{0:d} \n'.format(0))
                fid.write('{0:.3f} \n'.format(c2threshe1))
            else:
                fid.write('{0:d} \n'.format(nz))
                fid.write('{0:.3f} \n'.format(c2threshe1))
        else:
            pass
        fid.write(ofil+'\n')
        fid.write('{0:d} \n'.format(nlev))
        
    elif ilev == 1:
        print 'Writing Advanced mode'
        fid.write('{0:d} \n'.format(ilev))
        fid.write('{0:d} \n'.format(nout))
        fid.write('{0:d} \n'.format(ninp))
        fid.write('{0:d} \n'.format(nref))
        if nref>3:
            nrrlist=np.array([len(rrlist) 
                            for rrlist in pdict['rrfn_list']])
            nr3=len(np.where(nrrlist==3)[0])
            nr2=len(np.where(nrrlist==2)[0])
            fid.write('{0:d},{1:d} \n'.format(nr3,nr2))
        fid.write('{0:d} \n'.format(nrr))
        #if remote referencing
        if int(nrr) == 0:
            fid.write('{0:.3f} \n'.format(tbw))
            fid.write('{0:.3f} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshe))
            #parameters for bz component if ninp=3
            if nout == 3:
                if c2threshe != 0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                else:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                if c2threshe1 != 0.0 or c2threshe != 0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            else:
                if c2threshe != 0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej = [prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        #if 2 stage processing
        elif int(nrr) == 1:
            fid.write('{0:.5g} \n'.format(tbw))
            fid.write('{0:.5g} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))        
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshb))        
            fid.write('{0:.3f} \n'.format(c2threshe))
            if nout == 3:
                if c2threshb != 0 or c2threshe != 0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                elif c2threshb == 0 and c2threshe == 0:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(0))
            if c2threshb != 0.0 or c2threshe != 0.0:
                fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej = [prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        
    fid.write('{0:d} \n'.format(npcs))    
    fid.write('{0:d} \n'.format(nar))    
    fid.write('{0:d} \n'.format(imode))    
    fid.write('{0:d} \n'.format(jmode)) 
    
    #!!!NEED TO SORT FILE NAMES SUCH THAT EX, EY, HZ, HX, HY or EX, EY, HX, HY
    
    #write in filenames 
    if npcs != 1:
        if jmode == 0:
            fid.write(str(nread[0])+'\n')
            #--> write filenames to process with other information for first
            #    time section
            for tt, tfile in enumerate(pdict['fn_list'][0]):
                #write in calibration files if given
                if tt == 2:
                    if hx_cal is not None:
                        fid.write('-2\n')
                        fid.write(hx_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                elif tt == 3:
                    if hy_cal is not None:
                        fid.write('-2\n')
                        fid.write(hy_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                elif tt == 4:
                    if hz_cal is not None:
                        fid.write('-2\n')
                        fid.write(hz_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                else:
                    fid.write(str(nfil)+'\n')
                fid.write(tfile+'\n')
                fid.write(str(nskip[0])+'\n')
                
            #--> write remote reference time series
            for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                if rr == 0:
                    if rrhx_cal is not None:
                        fid.write('-2\n')
                        fid.write(rrhx_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                elif rr == 1:
                    if rrhy_cal is not None:
                        fid.write('-2\n')
                        fid.write(rrhy_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                fid.write(rfile+'\n')
                fid.write(str(nskipr[0])+'\n')
            
            #--> write in other pieces if there are more, note calibrations 
            #    are only given for the first time block so it is assumed
            #    that the same remote referenc is used for all time blocks
            for nn in range(1,npcs):
                fid.write(str(nread[nn])+'\n')            
                #write filenames
                for tfile in pdict['fn_list'][nn]:
                    fid.write(tfile+'\n')
                    fid.write(str(nskip[0])+'\n')
                for rfile in pdict['rrfn_list'][nn]:
                    fid.write(rfile+'\n')
                    fid.write(str(nskipr[nn])+'\n')
                    
        #--> if start and end time are give write in those
        elif jmode == 1:
            #write filenames
            for tt, tfile in enumerate(pdict['fn_list'][0]):
                #write in calibration files if given
                if tt == 2:
                    if hx_cal is not None:
                        fid.write('-2\n')
                        fid.write(hx_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                elif tt == 3:
                    if hy_cal is not None:
                        fid.write('-2\n')
                        fid.write(hy_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                elif tt == 4:
                    if hz_cal is not None:
                        fid.write('-2\n')
                        fid.write(hz_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                else:
                    fid.write(str(nfil)+'\n')
                fid.write(tfile+'\n')
                fid.write(dstim+'\n')
                fid.write(wstim+'\n')
                fid.write(wetim+'\n')
                
            #--> write remote referenc information
            for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                if rr == 0:
                    if rrhx_cal is not None:
                        fid.write('-2\n')
                        fid.write(rrhx_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                if rr == 1:
                    if rrhy_cal is not None:
                        fid.write('-2\n')
                        fid.write(rrhy_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                fid.write(rfile+'\n')
                fid.write(dstim+'\n')
                fid.write(wstim+'\n')
                fid.write(wetim+'\n')
                
            #--> write other time blocks
            for nn in range(1,npcs):
                fid.write(str(nread[nn])+'\n')            
                #write filenames
                for tfile in pdict['fn_list'][nn]:
                    fid.write(tfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
                for rfile in pdict['rrfn_list'][nn]:
                    fid.write(rfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
    else:
        if jmode == 0:
            if type(nread) is list:
                fid.write(str(nread[0])+'\n')
            else:
                fid.write(str(nread)+'\n')
            #--> write filenames for first block
            if nds==0:
                for tt, tfile in enumerate(pdict['fn_list']):
                    if tt == 2:
                        if hx_cal is not None:
                            fid.write('-2\n')
                            fid.write(hx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 3:
                        if hy_cal is not None:
                            fid.write('-2\n')
                            fid.write(hy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 4:
                        if hz_cal is not None:
                            fid.write('-2\n')
                            fid.write(hz_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                    fid.write(tfile+'\n')
                    if type(nskip) is list:
                        fid.write(str(nskip[0])+'\n')
                    else:
                        fid.write(str(nskip)+'\n')
                for rr, rfile in enumerate(pdict['rrfn_list']):
                    if rr == 0:
                        if rrhx_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    if rr == 1:
                        if rrhy_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    fid.write(rfile+'\n')
                    if type(nskipr) is list:
                        fid.write(str(nskipr[0])+'\n')
                    else:
                        fid.write(str(nskipr)+'\n')
            else:
                for tfile in pdict['fn_list'][0]:
                    fid.write(str(nfil)+'\n')
                    fid.write(tfile+'\n')
                    if type(nskip) is list:
                        fid.write(str(nskip[0])+'\n')
                    else:
                        fid.write(str(nskip)+'\n')
                for rfile in pdict['rrfn_list'][0]:
                    fid.write(str(nfil)+'\n')
                    fid.write(rfile+'\n')
                    if type(nskipr) is list:
                        fid.write(str(nskipr[0])+'\n')
                    else:
                        fid.write(str(nskipr)+'\n')
                        
        elif jmode == 1:
            #write filenames
            if nds==0:
                for tt, tfile in enumerate(pdict['fn_list']):
                    if tt == 2:
                        if hx_cal is not None:
                            fid.write('-2\n')
                            fid.write(hx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 3:
                        if hy_cal is not None:
                            fid.write('-2\n')
                            fid.write(hy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 4:
                        if hz_cal is not None:
                            fid.write('-2\n')
                            fid.write(hz_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                    fid.write(tfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
                for rr, rfile in enumerate(pdict['rrfn_list']):
                    if rr == 0:
                        if rrhx_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    if rr == 1:
                        if rrhy_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    fid.write(rfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
            else:
                for tt, tfile in enumerate(pdict['fn_list'][0]):
                    if tt == 2:
                        if hx_cal is not None:
                            fid.write('-2\n')
                            fid.write(hx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 3:
                        if hy_cal is not None:
                            fid.write('-2\n')
                            fid.write(hy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    elif tt == 4:
                        if hz_cal is not None:
                            fid.write('-2\n')
                            fid.write(hz_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                    fid.write(tfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
                for rr, rfile in enumerate(pdict['rrfn_list'][0]):
                    if rr == 0:
                        if rrhx_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhx_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    if rr == 1:
                        if rrhy_cal is not None:
                            fid.write('-2\n')
                            fid.write(rrhy_cal+'\n')
                        else:
                            fid.write(str(nfil)+'\n')
                    fid.write(rfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
                    
    #write rotation angles
    fid.write(thetae.replace(',',' ')+'\n')
    fid.write(thetab.replace(',',' ')+'\n')
    fid.write(thetaf.replace(',',' ')+'\n')    
    fid.close()
    
    birrp_dict = {}
    
    if ilev == 0:
        birrp_dict['ilev'] = ilev
        birrp_dict['nout'] = nout
        birrp_dict['ninp'] = ninp
        birrp_dict['tbw'] = tbw
        birrp_dict['sampling_rate'] = deltat
        birrp_dict['nfft'] = nfft
        birrp_dict['nsctmax'] = nsctmax
        birrp_dict['uin'] = uin
        birrp_dict['ainuin'] = ainuin
        birrp_dict['c2threshe'] = c2threshe
        birrp_dict['nz'] = nz
        birrp_dict['c2threshe1'] = c2threshe1
        birrp_dict['ofil'] = ofil
        birrp_dict['nlev'] = nlev
        birrp_dict['npcs'] = npcs
        birrp_dict['n_samples'] = nread
        birrp_dict['nar'] = nar
        birrp_dict['imode'] = imode
        birrp_dict['jmode'] = jmode
        birrp_dict['nfil'] = nfil
        birrp_dict['nskip'] = nskip
        birrp_dict['nskipr'] = nskipr
        birrp_dict['thetae'] = thetae
        birrp_dict['thetab'] = thetab
        birrp_dict['thetaf'] = thetaf
    elif ilev == 1:
        birrp_dict['ilev'] = ilev
        birrp_dict['nout'] = nout
        birrp_dict['ninp'] = ninp
        birrp_dict['nref'] = nref
        birrp_dict['nrr'] = nrr
        birrp_dict['tbw'] = tbw
        birrp_dict['sampling_rate'] = deltat
        birrp_dict['nfft'] = nfft
        birrp_dict['nsctinc'] = nsctinc
        birrp_dict['nsctmax'] = nsctmax
        birrp_dict['nf1'] = nf1
        birrp_dict['nfinc'] = nfinc
        birrp_dict['nfsect'] = nfsect
        birrp_dict['uin'] = uin
        birrp_dict['ainlin'] = ainlin
        birrp_dict['ainuin'] = ainuin
        if nrr == 1:
            birrp_dict['c2threshb'] = c2threshb
            birrp_dict['c2threshe'] = c2threshe
            if c2threshe == 0 and c2threshb == 0:
                birrp_dict['nz'] = 0
            else:
                birrp_dict['nz'] = 0
                birrp_dict['perlo'] = perlo
                birrp_dict['perhi'] = perhi
        elif nrr == 0:
            birrp_dict['c2threshb'] = 0
            birrp_dict['c2threshe'] = c2threshe
        birrp_dict['nprej'] = nprej
        birrp_dict['prej'] = prej
        birrp_dict['c2threshe1'] = c2threshe1
        birrp_dict['ofil'] = ofil
        birrp_dict['npcs'] = npcs
        birrp_dict['n_samples'] = nread
        birrp_dict['nlev'] = nlev
        birrp_dict['nar'] = nar
        birrp_dict['imode'] = imode
        birrp_dict['jmode'] = jmode
        birrp_dict['nfil'] = nfil
        if jmode == 0:
            birrp_dict['nskip'] = nskip
            birrp_dict['nskipr'] = nskipr
        elif jmode == 1:
            birrp_dict['dstim'] = dstim
            birrp_dict['wstim'] = wstim
            birrp_dict['wetim'] = wetim
        birrp_dict['thetae'] = thetae
        birrp_dict['thetab'] = thetab
        birrp_dict['thetaf'] = thetaf

    print 'Wrote BIRRP script file: {0}.script'.format(ofil)
    
    return scriptfile,birrp_dict

def run(birrp_exe, script_file):
    """
    run a birrp script file
    
    """
    if not op.isfile(birrp_exe):
        raise MTex.MTpyError_inputarguments('birrp executable not found:'+
                                            '{0}'.format(birrp_exe))

    current_dir = op.abspath(os.curdir)

    #change directory to directory of the script file
    os.chdir(os.path.dirname(script_file))

    sfid = file(script_file,'r')
    inputstring = ''.join(sfid.readlines())
    sfid.close()
    #correct inputstring for potential errorneous line endings due to strange
    #operating systems:
    tempstring = inputstring.split()
    tempstring = [i.strip() for i in tempstring]
    inputstring = '\n'.join(tempstring)
    inputstring += '\n'

    #open a log file to catch process and errors of BIRRP executable
    logfile = open('birrp_logfile.log','w')

    print 'starting Birrp processing at {0}...'.format(time.ctime())

    birrpprocess = subprocess.Popen(birrp_exe, 
                                    stdin=subprocess.PIPE, 
                                    stdout=logfile,
                                    stderr=logfile)

    out, err = birrpprocess.communicate(inputstring)
    
    logfile.close()

    print 'logfile closed: {0} at {1}'.format(logfile.name, time.ctime())
 
    print 'Outputs: {0}'.format(out)
    print 'Errors: {0}'.format(err)
    
    #go back to initial directory
    os.chdir(current_dir)

    print '\n{0} DONE !!! {0}\n'.format('='*20)


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


def convert2edi(stationname, in_dir, survey_configfile, birrp_configfile, 
                out_dir = None):
    """
    Convert BIRRP output files into EDI file.

    The list of BIRRP output files is searched for in the in_dir directory 
    for the given stationname (base part of filenames). Meta-data must be 
    provided in two config files.  If MTpy standard processing has been 
    applied, the first one is the same file as used from the beginning of
    the processing. If this survey-config-file is missing, a temporary 
    config file must been created from the header information of the time 
    series files that have been used as BIRRP input.
    
    A second config file contains information about the BIRRP processing 
    parameters. It's generated when BIRRP is called with MTpy.

    The outputfile 'stationname.edi' and is stored in the out_dir directory.
    If out_dir is not given, the files are stored in the in_dir.

    Input:
    - name of the station
    - directory, which contains the BIRRP output files
    - configuration file of the survey, containing all station setup information
    - configuration file for the processing of the station, containing all 
      BIRRP and other processing parameters
    [- location to store the EDI file]
    """ 

    stationname = stationname.upper()

    input_dir = op.abspath(op.realpath(in_dir))
    if not op.isdir(input_dir):
        raise MTex.MTpyError_inputarguments('Directory not existing:%s'%(input_dir))

    if out_dir == None:
        output_dir = input_dir
    else:
        output_dir = op.abspath(op.realpath(op.join(os.getcwd(),out_dir)))
        if not op.isdir(output_dir):
            try:
                os.makedirs(output_dir)
            except:
                print 'output directory could not be created - using input directory instead'
                output_dir = input_dir

    out_fn = op.join(output_dir,'{0}.edi'.format(stationname))

    if not op.isfile(survey_configfile):
        raise MTex.MTpyError_inputarguments('Survey - configfile not existing: "{0}"'.format(survey_configfile))
    if birrp_configfile is not None:
        if not op.isfile(birrp_configfile):
            raise MTex.MTpyError_inputarguments('BIRRP - Configfile not existing: "{0}"'.format(birrp_configfile))
     
    #read the survey config file:
    #try:
    survey_config_dict = MTcf.read_survey_configfile(survey_configfile)

    # except:
    #     raise EX.MTpyError_config_file( 'Config file cannot be read: %s' % (survey_configfile) )

    if not stationname in survey_config_dict:
        raise MTex.MTpyError_config_file( 'No information about station {0} found in configuration file: {1}'.format(stationname, survey_configfile) )

    station_config_dict = survey_config_dict[stationname]


    #read the BIRRP/processing config file:
    birrp_config_dict = {}
    if birrp_configfile is not None:
        try:
            birrp_config_dict = MTcf.read_configfile(birrp_configfile)
        except:
            print 'Config file with BIRRP processing parameters could not'\
            ' be read: {0} - using generic values'.format(birrp_configfile) 
            birrp_config_dict = {}
    

    #find the birrp-output j-file for the current station 
    #j_filename_list = [i for i in os.listdir(input_dir) if op.basename(i).upper() == ('%s.j'%stationname).upper() ]
    #find the birrp-output j-file for the current station 
    j_filename_list = [i for i in os.listdir(input_dir) if i.lower().endswith('.j') ]
    j_filename_list = [i for i in  j_filename_list if '{0}'.format(stationname.upper()) in op.basename(i).upper() ]
    j_filename_list = [op.join(input_dir,i) for i in j_filename_list]
    try:
        j_filename = j_filename_list[0]
    except:
        print 'j-file for station %s not found in directory %s'%(stationname, input_dir)
        raise MTex.MTpyError_file_handling
    
    if len(j_filename_list) > 1:
        print 'Warning - more than one j-file found - taking the first one only: {0}'.format(j_filename)



    #Having now:
    # station_config_dict - contains information about station setup
    # birrp_config_dict - contains information about the processing (BIRRP parameters, selected time window, Rem.Ref.,..)
    # directory - contains BIRRP output files, coded by stationname

    # To be converted into .EDI
    # Dictionaries information goes into EDI header: HEAD and INFO section - check for other sections though
    # output EDI file is out_fn
      
    periods, Z_array, tipper_array,processing_dict,sorting_dict = read_j_file(j_filename)


    HEAD = _set_edi_head(station_config_dict,birrp_config_dict)

    INFO = _set_edi_info(station_config_dict,birrp_config_dict, sorting_dict)

    DATA = _set_edi_data(periods, Z_array, tipper_array)

    DEFINEMEAS = _set_edi_defmeas(station_config_dict)
    MTSECT = _set_edi_mtsect(station_config_dict,periods)

    
    out_fn = MTfh.make_unique_filename(out_fn)

    F_out = open(out_fn,'w') 
    
    F_out.write(HEAD)
    F_out.write(INFO)
    F_out.write(DEFINEMEAS)
    F_out.write(MTSECT)
    F_out.write(DATA)
    F_out.write('>END\n')

    F_out.close()

    return out_fn

    
def convert2edi_incl_instrument_correction(stationname, in_dir, 
                                        survey_configfile, birrp_configfile, 
                                        instr_response_file, out_dir = None, instr_type='lemi'):
    """
    Convert BIRRP output files into EDI file.

    The list of BIRRP output files is searched for in the in_dir directory for the given stationname (base part of filenames). Meta-data must be provided in two config files.  If MTpy standard processing has been applied, the first one is the same file as used from the beginning of the processing. If this survey-config-file is missing, a temporary config file must been created from the header information of the time series files that have been used as BIRRP input.
    A second config file contains information about the BIRRP processing parameters. It's generated when BIRRP is called with MTpy.

    The outputfile 'stationname.edi' and is stored in the out_dir directory. If out_dir is not given, the files are stored in the in_dir.

    Input:
    - name of the station
    - directory, which contains the BIRRP output files
    - configuration file of the survey, containing all station setup information
    - configuration file for the processing of the station, containing all BIRRP and other processing parameters
    - instrument response file (3 column data: frequencies, real, imaginary)
    [- location to store the EDI file]
    """ 

    stationname = stationname.upper()

    input_dir = op.abspath(op.realpath(in_dir))
    if not op.isdir(input_dir):
        raise MTex.MTpyError_inputarguments('Directory not existing:%s'%(input_dir))

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
    out_fn = op.join(output_dir,'{0}.edi'.format(stationname))

    if not op.isfile(survey_configfile):
        raise MTex.MTpyError_inputarguments('Survey - configfile not existing: "{0}"'.format(survey_configfile))
    
    if birrp_configfile is not None:
        if not op.isfile(birrp_configfile):
            raise MTex.MTpyError_inputarguments('BIRRP - Configfile not existing: "{0}"'.format(birrp_configfile))

     
    #read the survey config file:
    #try:
    survey_config_dict = MTcf.read_survey_configfile(survey_configfile)

    # except:
    #     raise EX.MTpyError_config_file( 'Config file cannot be read: %s' % (survey_configfile) )

    if not stationname in survey_config_dict:
        print  'No information about station {0} found in configuration file: {1}'.format(stationname, survey_configfile)  
        raise MTex.MTpyError_config_file()

    station_config_dict = survey_config_dict[stationname]

    instr_response_file = op.abspath(instr_response_file)
    try: 
        instr_resp = np.loadtxt(instr_response_file)
    except:
        sys.exit('ERROR - cannot read instrument response file {0}'.format(instr_response_file))


    try:
        if not instr_resp.shape[1]==3:
            raise
    except:
        sys.exit('ERROR - instrument response file {0} has wrong format - need 3 columns'.format(instr_response_file))

    
    #read the BIRRP/processing config file:
    birrp_config_dict = {}
    if birrp_configfile is not None:
        try:
            birrp_config_dict = MTcf.read_configfile(birrp_configfile)
        except:
            raise MTex.MTpyError_config_file( 'Config file with BIRRP processing parameters could not be read: %s' % (birrp_configfile) )

    #find the birrp-output j-file for the current station 
    j_filename_list = [i for i in os.listdir(input_dir) if i.lower().endswith('.j') ]
    j_filename_list = [i for i in  j_filename_list if '{0}'.format(stationname.upper()) in op.basename(i).upper() ]
    j_filename_list = [op.join(input_dir,i) for i in j_filename_list]
    try:
        j_filename = j_filename_list[0]
    except:
        print 'j-file for station %s not found in directory %s'%(stationname, input_dir)
        raise MTex.MTpyError_file_handling
    
    if len(j_filename_list) > 1:
        print 'Warning - more than one j-file found - taking the first one only: {0}'.format(j_filename)
        #raise MTex.MTpyError_file_handling('More than one j-file for station %s found in directory %s'%(stationname, input_dir))

    #Having now:
    # station_config_dict - contains information about station setup
    # birrp_config_dict - contains information about the processing (BIRRP parameters, selected time window, Rem.Ref.,..)
    # directory - contains BIRRP output files, coded by stationname

    # To be converted into .EDI
    # Dictionaries information goes into EDI header: HEAD and INFO section - check for other sections though
    # output EDI file is out_fn
    
    periods, Z_array, tipper_array,processing_dict,sorting_dict = read_j_file(j_filename)

    #add meta data from j-file to the processing dictionary:
    birrp_config_dict.update(processing_dict)

    frequencies = 1./periods
    def correct_z_for_instrument_response(Z_array, instr_resp, frequencies, instr_type):

        for idx_f, freq in enumerate(frequencies):
            if not (instr_resp[0,0] <= np.abs(freq) <= instr_resp[-1,0]):
                print 'no instrument response in this frequency range - array values set to zero here: ', freq
                correction_factor = 0.
                continue
            #find the appropriate frequencies ( since the current freq-value is 
            #most likely inbetween two values on the instr_freqs-axis) 
            #- get the respective value by interpolation !

            #find the value closest to the current freq, assume it's lower
            closest_lower = np.abs(freq-instr_resp[:,0]).argmin()

            #if it coincides with the highest frequency/last entry:
            if freq == instr_resp[-1,0]:
                correction_factor = np.complex(instr_resp[-1,1],instr_resp[-1,2])
    
            #if it coincides with the lowest frequency/first entry:
            elif freq == instr_resp[0,0]:
                correction_factor = np.complex(instr_resp[0,1],instr_resp[0,2])
    
            else:

                correction_factor = MTip.interpolate_instrumentresponse(freq,instr_resp ,instr_type)

            
            #finally correct Z for the instrument influence by multiplying with the instrument response value: 
            for i in range(4):
                zentry = np.complex(Z_array[idx_f,0,i], Z_array[idx_f,1,i] )
                corrected = zentry * correction_factor 
                Z_array[idx_f,0,i] = np.real(corrected)
                Z_array[idx_f,1,i] = np.imag(corrected)
                #correct the error: stretch by the ratio of amplitudes of original and correct Z value
                Z_array[idx_f,2,i] = Z_array[idx_f,2,i]/np.abs(zentry)*np.abs(corrected)


        return Z_array

    Z_array = correct_z_for_instrument_response(Z_array, instr_resp, frequencies,instr_type)



    HEAD = _set_edi_head(station_config_dict,birrp_config_dict)
    INFO = _set_edi_info(station_config_dict,birrp_config_dict,sorting_dict)
    DATA = _set_edi_data(periods, Z_array, tipper_array)
    DEFINEMEAS = _set_edi_defmeas(station_config_dict)
    MTSECT = _set_edi_mtsect(station_config_dict,periods)


    out_fn = MTfh.make_unique_filename(out_fn)

    F_out = open(out_fn,'w') 
    
    F_out.write(HEAD)
    F_out.write(INFO)
    F_out.write(DEFINEMEAS)
    F_out.write(MTSECT)
    F_out.write(DATA)
    F_out.write('>END\n')

    F_out.close()

    return out_fn



def _set_edi_data(lo_periods, Z_array, tipper_array):
    

    periods = lo_periods

    datastring = ''
    
    datastring += '>ZROT // %i\n'%(len(periods))
    for i,period in enumerate(periods):
        freq = 1./period
        datastring += '\t%E'%(0.)
        if (i+1)%5 == 0 and (i != len(periods) - 1) and i > 0:
            datastring += '\n'

    datastring += '\n'


    #datastring += '>!****FREQUENCIES****!\n'
    datastring += '>FREQ // %i\n'%(len(periods))
    for i,period in enumerate(periods):
        freq = 1./period
        datastring += '\t%E'%(freq)
        if (i+1)%5 == 0 and (i != len(periods) - 1) and i > 0:
            datastring += '\n'

    datastring += '\n'

    #datastring += '>!****IMPEDANCES****!\n'
    compstrings = ['ZXX','ZXY','ZYX','ZYY']
    Z_entries = ['R','I','.VAR']
    
    for Z_comp in range(4):
        for entry in range(3):
            datastring += '>%s%s ROT=ZROT // %i\n'%(compstrings[Z_comp], Z_entries[entry], len(periods))
            for i,period in enumerate(periods):
                data = Z_array[i,entry,Z_comp]
                #EDI files carries variances, not standard deviations:
                if entry == 2:
                    data = data**2
                datastring += '\t%E'%(data)
                if (i+1)%5 == 0 and (i != len(periods) - 1) and i > 0:
                    datastring += '\n'
                
            datastring += '\n'

        
        
    #datastring += '\n'
    #datastring += '>!****TIPPER****!\n'

    compstrings = ['TX','TY']
    T_entries = ['R.EXP','I.EXP','VAR.EXP']
    
    for T_comp in range(2):
        for entry in range(3):
            datastring += '>%s%s ROT=ZROT // %i\n'%(compstrings[T_comp], T_entries[entry], len(periods))
            for i,period in enumerate(periods):
                if tipper_array != None :
                    data = tipper_array[i,entry,T_comp]
                else:
                    data = 0.
                datastring += '\t%E'%(data)
                if (i+1)%5 == 0 and (i != len(periods) - 1) and i > 0:
                    datastring += '\n'
                

            datastring += '\n'


    datastring += '\n'

    return datastring.expandtabs(4)


def _set_edi_info(station_config_dict,birrp_config_dict,sorting_dict):

    infostring = ''
    infostring += '>INFO\t max lines=1000\n'
    infostring += '\t\tedifile_generated_with: MTpy\n'
    infostring += '\t\tTF_processing: BIRRP \n'
    infostring += '\t\tZ_unit: km/s \n\n'
   
  
    infostring += '\tStation parameters\n'

    for key in sorted(station_config_dict.iterkeys()):
        infostring += '\t\t{0}: {1}  \n'.format(str(key),
                                                str(station_config_dict[key]))   
    

    infostring += '\n'
    infostring += '\tBIRRP processing parameters\n'

    temp_dict = birrp_config_dict.copy()
    #first write out all imtems, which are vcontained in the sorting dict, then 
    #append the ones that are left
    for key in sorted(sorting_dict):
        newkey = sorting_dict[key]
        try:
            infostring += '\t\t{0}: {1}  \n'.format(str(newkey),
                                                temp_dict.pop(newkey))
        except:
            continue


    for key in sorted(temp_dict.iterkeys()):
        infostring += '\t\t{0}: {1}  \n'.format(str(key),
                                                str(temp_dict[key]))   

    # for key in sorted(birrp_config_dict.iterkeys()):
    #     infostring += '\t\t{0}: {1}  \n'.format(str(key),
    #                                             str(birrp_config_dict[key]))   

    infostring += '\n'    
    if len(birrp_config_dict) == 0 :
        infostring += '\t\tunknown\n\n'


    infostring += '\n'

    return infostring.expandtabs(4)


def _set_edi_head(station_config_dict,birrp_config_dict):
    """
    set header string
    
    set date to format YYYY/MM/DD HH:MM:SS UTC
    """
    frmt = '%Y/%m/%d %H:%M:%S UTC'

    headstring = ''
    headstring += '>HEAD\n'
    if len((station_config_dict['station'].split())) != 1: 
        headstring += '\tdataid="%s"\n'%(station_config_dict['station'])
    else:
        headstring += '\tdataid={0}\n'.format(station_config_dict['station'])


    if station_config_dict.has_key('company'):
        acqby = station_config_dict.has_key('company')
    elif station_config_dict.has_key('acqby'):
        acqby = station_config_dict.has_key('acqby')
    else:
        acqby = ''

    if len(acqby.split()) != 1:
        headstring += '\tacqby="%s"\n'%(acqby)
    else:
        headstring += '\tacqby={0}\n'.format(acqby)


    if len(birrp_config_dict) !=0 :
        try:
            sampling_rate = float(birrp_config_dict['sampling_rate'])
        except:
            sampling_rate = float(birrp_config_dict['sampling'])

        try:
            n_samples = int(birrp_config_dict['n_samples'])
        except ValueError:
            n_samples = sum([int(ii) for ii in 
                             birrp_config_dict['n_samples'][1:-1].strip().split(',')])
        #new:
        try:
            acq_starttime = float(birrp_config_dict['processing_window_start'])
            dt_start = datetime.datetime.fromtimestamp(acq_starttime)
            acq_start = dt_start.combine(dt_start.date(),dt_start.time()).strftime(frmt)

            dt_end_time =  acq_starttime + 1./sampling_rate*(n_samples) 
            dt_end = datetime.datetime.fromtimestamp(dt_end_time)
            acq_end = dt_end.combine(dt_end.date(),dt_end.time()).strftime(frmt)
            
            headstring +='\tacqdate=%s \n'%(acq_start)
            headstring +='\tenddate=%s \n'%(acq_end)
        except KeyError:
            try:
                acq_start = station_config_dict.has_key('acq_date')
                headstring += '\tacqdate={0} \n'.format(acq_start)
            except KeyError:
                 headstring += '\tacqdate={0} \n'.format('1970/01/01 00:00:00 UTC')
                

    todaystring = datetime.datetime.utcnow().strftime(frmt)
    headstring += '\tfiledate="%s"\n'%(todaystring)


    network = ''
    if station_config_dict.has_key('network'):
        location = station_config_dict.has_key('network')
    headstring += '\tprospect="%s"\n'%(network)

    location = ''
    if station_config_dict.has_key('location'):
        location = station_config_dict.has_key('location')
    headstring += '\tloc="%s"\n'%(location)

    #headstring += '\tlat=%.5f\n'%float(station_config_dict['latitude'])
    #headstring += '\tlong=%.5f\n'%float(station_config_dict['longitude'])

    lattuple = MTft.convert_degrees2dms_tuple(float(station_config_dict['latitude']))
    lontuple = MTft.convert_degrees2dms_tuple(float(station_config_dict['longitude']))
    headstring += '\tlat=%s\n'%(MTft.convert_dms_tuple2string(lattuple))
    headstring += '\tlong=%s\n'%(MTft.convert_dms_tuple2string(lontuple))
    headstring += '\telev=%.1f\n'%float(station_config_dict['elevation'])

    headstring += '\n'

    return headstring.upper().expandtabs(4)


def _set_edi_defmeas(station_config_dict):

    dmeasstring = ''
    dmeasstring += '>=DEFINEMEAS\n'
    dmeasstring += '\n'

    dmeasstring += '\tmaxchan=7\n'
    dmeasstring += '\tmaxrun=999\n'
    dmeasstring += '\tmaxmeas=9999\n'

    #NOT necessary:
    #dmeasstring += '\tunits=m\n'
    #dmeasstring += '\treftype="WGS 84"\n'
    lattuple = MTft.convert_degrees2dms_tuple(float(station_config_dict['latitude']))
    lontuple = MTft.convert_degrees2dms_tuple(float(station_config_dict['longitude']))    
    dmeasstring += '\treflat=%s\n'%(MTft.convert_dms_tuple2string(lattuple))
    dmeasstring += '\treflong=%s\n'%(MTft.convert_dms_tuple2string(lontuple))
    dmeasstring += '\trefelev=%.1f\n'%float(station_config_dict['elevation'])
    
    dmeasstring += '\n'
    dmeasstring += '>HMEAS id=1001.001 chtype=hx x=0. y=0. azm=0.\n'
    dmeasstring += '>HMEAS id=1002.001 chtype=hy x=0. y=0. azm=90.\n'
    dmeasstring += '>HMEAS id=1003.001 chtype=hz x=0. y=0. azm=0.\n'

    try:
        dmeasstring += '>EMEAS id=1004.001 chtype=ex x=0. y=0. x2=%.1f y2=0\n'%float(station_config_dict['e_xaxis_length'])
    except:
        dmeasstring += '>EMEAS id=1004.001 chtype=ex x=0. y=0. x2=0. y2=0.\n'
        
    try:
        dmeasstring += '>EMEAS id=1005.001 chtype=ey x=0. y=0. x2=0. y2=%.1f\n'%float(station_config_dict['e_yaxis_length'])
    except:
        dmeasstring += '>EMEAS id=1005.001 chtype=ey x=0. y=0. x2=0. y2=0.\n'

    dmeasstring += '>HMEAS id=1006.001 chtype=rx x=0. y=0. azm=0.\n'
    dmeasstring += '>HMEAS id=1007.001 chtype=ry x=0. y=0. azm=90.\n'


    dmeasstring += '\n'

    return dmeasstring.expandtabs(4)


def _set_edi_mtsect(station_config_dict,periods):
    mtsectstring = ''
    mtsectstring += '>=MTSECT\n' 
    mtsectstring += '\tsectid=%s\n'%station_config_dict['station']
    mtsectstring += '\tnfreq=%i\n'%(len(periods))
    mtsectstring += '\thx=1001.001\n'
    mtsectstring += '\thy=1002.001\n'
    mtsectstring += '\thz=1003.001\n'
    mtsectstring += '\tex=1004.001\n'
    mtsectstring += '\tey=1005.001\n'
    mtsectstring += '\trx=1006.001\n'
    mtsectstring += '\try=1007.001\n'

    mtsectstring += '\n'

    return mtsectstring.expandtabs(4)




def read_j_file(fn):
    """
    read_j_file will read in a *.j file output by BIRRP (better than reading lots of *.<k>r<l>.rf files)

    Input:
    j-filename

    Output: 4-tuple
    - periods : N-array
    - Z_array : 2-tuple - values and errors
    - tipper_array : 2-tuple - values and errors
    - processing_dict : parsed processing parameters from j-file header

    """   

    j_fn = op.abspath(fn)
    if not op.isfile(j_fn):
        raise MTex.MTpyError_inputarguments('Cannot read j-file %s - file is not existing'%(j_fn))

    
    
    with open(j_fn,'r') as F_in:
        j_lines = F_in.readlines()

    processing_dict,sorting_dict = parse_jfile_header(j_lines)

    
    Z_start_row = None
    tipper_start_row = None
    tipper = None

    for idx_jline,j_line in enumerate(j_lines):

        if 'ZXX' == j_line.upper().strip()[:3]:
            Z_start_row = idx_jline
        if 'TZX' == j_line.upper().strip()[:3]:
            tipper_start_row = idx_jline 

    try:
        n_periods = int(float(j_lines[Z_start_row + 1] ))
    except:
        raise MTex.MTpyError_inputarguments('File is not a proper j-file: %s'%(j_fn))

    Z = np.zeros((n_periods,3,4))
    periods = np.zeros((n_periods,4))
    if not tipper_start_row == None:
        tipper = np.zeros((n_periods,3,2))
        periods = np.zeros((n_periods,6))


    for idx_comp in range(4):
        starting_row = Z_start_row + 2 + ((n_periods +2)* idx_comp)
        for idx_per in range(n_periods):
            idx_row = starting_row + idx_per
            cur_row = j_lines[idx_row]
            #print idx_row, cur_row
            row_entries = cur_row.strip().split()
            try:
                periods[idx_per,idx_comp] = float(row_entries[0])
            except:
                periods[idx_per,idx_comp] = np.nan

            if periods[idx_per,idx_comp] == -999:
                periods[idx_per,idx_comp] = np.nan

            for idx_z_entry in range(3):
                raw_value = row_entries[idx_z_entry + 1]
                try:
                    value = float(raw_value)
                except:
                    value = np.nan
                if value == -999:
                    value = np.nan


                Z[idx_per,idx_z_entry,idx_comp] = value

    if tipper != None :
            for idx_comp in range(2):
                starting_row = tipper_start_row+2+((n_periods +2)*idx_comp)
                for idx_per in range(n_periods):
                    idx_row = starting_row + idx_per
                    cur_row = j_lines[idx_row]
                    row_entries = cur_row.strip().split()
                    try:
                        periods[idx_per,idx_comp+4] = float(row_entries[0])
                    except:
                        periods[idx_per,idx_comp+4] = np.nan
                    if periods[idx_per,idx_comp+4] == -999:
                        periods[idx_per,idx_comp+4] = np.nan

                    for idx_z_entry in range(3):
                        raw_value = row_entries[idx_z_entry + 1]
                        try:
                            value = float(raw_value)
                        except:
                            value = np.nan
                        if value == -999:
                            value = np.nan
                        tipper[idx_per,idx_z_entry,idx_comp] = value

    
    #NOTE: j files can contain periods that are NOT sorted increasingly, but random
    indexorder = np.array([iii[0] for iii in periods]).argsort()
    periods = periods[indexorder]
    Z = Z[indexorder]
    if tipper is not None:
        tipper = tipper[indexorder]
    
    periods,Z,tipper = _check_j_file_content(periods, Z, tipper)

    return periods, Z, tipper, processing_dict,sorting_dict

def parse_jfile_header(j_lines):

    """
    Parsing the header lines of a j-file to extract processing information.

    Input:
    - j-file as list of lines (output of readlines())

    Output:
    - Dictionary with all parameters found

    !TODO!
    
    """
    header_dict = {}
    sorting_dict = {}
    tuples = []

    for line in j_lines:
        if not '=' in line: continue
        if not '#' in line: continue
        line = line.strip().replace('#','')
        line=[i.strip() for i in line.split('=')]
        no_keys = len(line)-1
        elements = [line[0]]
        for i in np.arange(no_keys-1)+1:
            elements.extend(line[i].split())
        elements.append(line[-1])

        for i in np.arange(no_keys)*2:
            tuples.append([elements[i],elements[i+1]])

        #print elements
    #print tuples 

    for dict_idx,pair in enumerate(tuples):
        k= pair[0]
        v= pair[1]
        if len(v) == 0:
            continue
        try:
            v = float(v)
            try:
                if v%1 ==0:
                    v = int(v)
            except:
                pass
        except:
            pass
        if k=='deltat':
            header_dict['sampling_rate'] = 1./v
        if k=='nread':
            header_dict['n_samples']=v

        idx = 2
        if k in header_dict.keys():
            knew = k
            while knew in header_dict.keys():
                knew = '{0}_{1}'.format(k,idx)
                idx += 1
            k = knew
        header_dict[k] = v

        sorting_dict[dict_idx+1]=k


    return header_dict,sorting_dict
 

def _check_j_file_content( periods_array, Z_array, tipper_array):
    """ 
    Check the content of j file.
    
    If 'nan' appears at any part for some period, the respective period must be 
    deleted together with all respective entries of the Z_array and tipper_array.
    Additionally, check the entries of the period array. This should have fully 
    redundant entries. If this is not the case for at least one period for at 
    least one component, the period and all respective entries of the arrays 
    have to be deleted.
    """
    period_epsilon = 1E-7
    lo_periods = []

    lo_all_periods_raw = list(set(periods_array.flatten()))
    lo_all_periods_raw = [i for i in lo_all_periods_raw if not np.isnan(i)]
    #print lo_all_periods_raw
    #lo_all_periods_raw.sort()
    lo_all_periods = np.array(sorted(lo_all_periods_raw))


    n_period_entries = periods_array.shape[1]

    for idx_period, period in enumerate(lo_all_periods):
        tmp_lo_period_idxs = []
        foundnan = 0
        for i in range(n_period_entries):
            #loop over all 4/6 components of Z and tipper
            
            #check, where the current period appears for the current component
            coinc = np.where(period == periods_array[:,i])[0]

            #check, if period is found exactly once for this component
            if len(coinc) == 1:
                #check all components for NaN:
                for j in range(3):
                    #only Z:
                    if i < 4:
                        if math.isnan(Z_array[coinc[0], j, i]):
                            foundnan = 1
                    else:
                        if math.isnan(tipper_array[coinc[0], j, i - 4]):
                            foundnan = 1
                if foundnan == 0:
                    tmp_lo_period_idxs.append(coinc[0])

        if len(tmp_lo_period_idxs) == n_period_entries:
            lo_periods.append( (period,tuple(tmp_lo_period_idxs)) )

    Z_array_out = np.zeros((len(lo_periods),3,4))
    tipper_array_out = None
    if n_period_entries == 6:
        tipper_array_out = np.zeros((len(lo_periods),3,2))

    lo_periods_out = []

    for idx in range(len(lo_periods)):
        lo_periods_out.append(lo_periods[idx][0])
        idx_tuple = lo_periods[idx][1]
        for j in range(4):
            Z_array_out[idx,:,j] = Z_array[idx_tuple[j],:,j]

        if n_period_entries == 6:
            for k in range(2):
                tipper_array_out[idx,:,k] = tipper_array[idx_tuple[k+4],:,k]

 

    return np.array(lo_periods_out), Z_array_out, tipper_array_out  


    

def convert2coh(stationname, birrp_output_directory):
    """
        Convert BIRRP output coherence files into just one *.coh file.
    """

    directory = op.abspath(birrp_output_directory)

    if not op.isdir(directory):
        raise MTex.MTpyError_inputarguments('Directory {0} does not exist'.format(directory))

    stationname = stationname.upper()
    #locate file names

    # only for second stage coherences...:
    
    cohfilenames = sorted([ op.abspath(op.join(directory,i)) for i in fnmatch.filter(
                os.listdir(directory), '*%s*.[12]r.[12]c2'%stationname.upper()) ] )
    if len(cohfilenames) == 0:
        lo_files = [i.lower() for i in os.listdir(directory)]
        cohfilenames = sorted([ op.abspath(op.join(directory,i)) for i in fnmatch.filter(
                lo_files, '*%s*.[12]r.[12]c2'%stationname.lower()) ] ) 
    
    #print cohfilenames

    if len(cohfilenames) < 1:
        print 'No coherence files for station %s found in: %s'%(stationname, directory)
        raise MTex.MTpyError_file_handling()#'No coherence files for station %s found in: %s'%(stationname, directory))


    # if len(cohfilenames) > 3:
    #     print 'Too many coherence files for station %s found in: %s'%(stationname, directory)
    #     raise MTex.MTpyError_file_handling()#'Too many coherence files for station %s found in: %s'%(stationname, directory))
    try:
        for fn in cohfilenames:
            if fn.lower().endswith('1r.2c2'):
                break
        period,freq,coh1,zcoh1 = MTfh.read_2c2_file(fn)
        for fn in cohfilenames:
            if fn.lower().endswith('2r.2c2'):
                break
        period,freq,coh2,zcoh2 = MTfh.read_2c2_file(fn)

    except:
        print 'Cannot read coherence files for station %s found in: %s'%(stationname, directory)
        raise MTex.MTpyError_file_handling()#'Cannot read coherence files for station %s found in: %s'%(stationname, directory))


    twostage=False

    for fn in cohfilenames:
        if fn.lower().endswith('1r.1c2') or fn.lower().endswith('2r.1c2'):
            twostage = True
            break
    if twostage is True:
        coh3 = None
        coh4 = None
        zcoh3 = None
        zcoh4 = None
        try:
            for fn in cohfilenames:
                if fn.lower().endswith('1r.1c2'):
                    break
            period,freq,coh3,zcoh3 = MTfh.read_2c2_file(fn)
        except:
            coh3 = None

        try:
            for fn in cohfilenames:
                if fn.lower().endswith('2r.1c2'):
                    break
            period,freq,coh4,zcoh4 = MTfh.read_2c2_file(fn)
        except:
            coh4 = None



    fn = '%s.coh'%(stationname)
    out_fn = op.abspath(op.join(directory,fn))
    #out_fn = MTfh.make_unique_filename(out_fn)


    F_out =  open(out_fn,'w')
    if twostage is False:
        F_out.write('#period \t freq \t\t cohEx \t zcohEx \t\t cohEy \t zcohEy \n'.expandtabs(4))
    else:
        F_out.write('#period \t freq \t\t cohEx \t zcohEx \t\t cohEy \t zcohEy'\
                    ' \t\t cohBx \t zcohBx \t\t cohBy \t zcohBx\n'.expandtabs(4))


    for ff in range(len(period)):
        tmp_string =''
        tmp_string += '{0:.5f} \t {1:.5f}\t\t'.format(period[ff], freq[ff])

        try:
            c1 = float(coh1[ff])
            zc1 = float(zcoh1[ff])
        except :
            c1= 0.
            zc1= 0.
        tmp_string += '{0:.5f} \t {1:.5f}\t\t'.format(c1,zc1)
        
        try:
            c2 = float(coh2[ff])
            zc2 = float(zcoh2[ff])
        except :
            c2 = 0.
            zc2 = 0.
        tmp_string += '{0:.5f} \t {1:.5f}\t\t'.format(c2,zc2)

        if twostage is True:
            c3 = 0.
            zc3 = 0.
            c4 = 0.
            zc4 = 0.
            if coh3 is not None:
                try:
                    c3 = float(coh3[ff])
                    zc3 = float(zcoh3[ff])
                except:
                    pass
            tmp_string += '{0:.5f} \t {1:.5f}\t\t'.format(c3,zc3)
            if coh4 is not None:
                try:
                    c4 = float(coh4[ff])
                    zc4 = float(zcoh4[ff])
                except:
                    pass
            tmp_string += '{0:.5f} \t {1:.5f}\t\t'.format(c4,zc4)

        tmp_string += '\n'
        F_out.write(tmp_string.expandtabs(4))
                   
        #F_out.write(('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n'%(period[ff], freq[ff], c1, zc1, c2, zc2, c3, zc3)).expandtabs(4))
    F_out.close()

    return out_fn
