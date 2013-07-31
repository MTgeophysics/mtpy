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
import StringIO
import numpy as np
import re
import sys, os
import glob
import os.path as op
import subprocess
import time
import datetime 
import fnmatch
import math

import mtpy.utils.exceptions as MTex
import mtpy.utils.format as MTft
import mtpy.utils.filehandling as MTfh
import mtpy.utils.configfile as MTcf

import mtpy.utils.misc as MTmc

reload(MTcf)
reload(MTft)
reload(MTfh)

#=================================================================


def runbirrp2in2out_simple(birrp_exe, stationname, ts_directory, 
                           coherence_threshold = 0.5, output_dir = None):
    """
    Call BIRRP for 2 input and 2 output channels with the simplemost setup. 

    Provide stationname and directory containing the data folders. Data must 
    be in 1-column ascii format, including one header line for identification 
    of station, channel and sampling rate. For this function, the files must 
    have aligned time series!
    
    All parameters are automatically determined, the threshold for coherence 
    is set to 0.5 by default.
    
    If no output directory is specified, output files are put into a subfolder 
    of the source directory, named 'birrp_processed'.

    Additionally, a configuration file is created. It contains information 
    about the processing paramters for the station. Keys are generic for the 
    common parameters and named after BIRRP input keywords for the processing 
    parameters.

    """

    if not op.isfile(birrp_exe):
        raise MTex.MTpyError_inputarguments('birrp executable not found:'+\
                                                       '{0}'.format(birrp_exe))
    if not op.isdir(ts_directory):
        raise MTex.MTpyError_inputarguments('time series files directory not'+\
                                            'existing: {0}'.format(ts_directory))

    current_dir = op.abspath(os.curdir)

    wd = op.abspath(op.realpath(op.join(ts_directory,'birrp_processed')))
    if output_dir != None:
        output_dir = op.abspath(output_dir)
        if not op.isdir(output_dir) :
            try:
                os.makedirs(output_dir)
                wd = output_dir
            except:
                print 'Could not find or generate specified output '+\
                      'directory {0} - using default instead!'.format(output_dir)
        else:
            wd = output_dir 

    if not op.isdir(wd):
        try:
            os.makedirs(wd) 
        except:
            raise MTex.MTpyError_file_handling('cannot create working directory:%s'%(wd))
    
    print "Successfully generated output directory", wd

    os.chdir(wd)

    inputstring, birrp_stationdict, inputfilename = generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold)

    print "generated inputstring and configuration dictionary for station {0}".format(stationname)
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

    print 'starting Birrp processing...'
    #os.system("{0} < {1}".format(birrp_exe,inputfilename))

    birrpprocess = subprocess.Popen(birrp_exe, stdin=subprocess.PIPE, stdout=logfile,stderr=logfile)
    #instringhandler = StringIO.StringIO(inputstring)
    #birrpprocess = subprocess.Popen(birrp_exe, stdin=instringhandler, stdout=logfile,stderr=logfile)

    out,err = birrpprocess.communicate(inputstring)
    
    #sys.stdout = dummy1
    #sys.stderr = dummy2
    logfile.close()

    print 'logfile closed: {0}'.format(logfile.name)

    #generate a local configuration file, containing information about all BIRRP and station parameters
    #required for the header of the EDI file 
    print out
    print err
    print 'generating configuration file containing the applied processing parameters'
    station_config_file = '%s_birrpconfig.cfg'%(stationname)
    MTcf.write_dict_to_configfile(birrp_stationdict, station_config_file)
    print 'Wrote BIRRP and time series configurations to file: {0}'.format(op.abspath(station_config_file))

    #go back to initial directory
    os.chdir(current_dir)

    print '\n \t\tDONE !!!\n\n'



def generate_birrp_inputstring_simple(stationname, ts_directory, coherence_threshold, output_channels=2):

    if not output_channels in [2,3]:
        raise MTex.MTpyError_inputarguments( 'Output channels must be 2 or 3' )


    print 'setting basic input components,e.g. filenames, samplingrate,...'
    input_filename, length, sampling_rate, birrp_stationdict = set_birrp_input_file_simple(stationname, ts_directory, output_channels, op.join(ts_directory,'birrp_wd'))

    print 'calculate optimal time window bisection parameters'
    longest_section, number_of_bisections = get_optimal_window_bisection(length, sampling_rate)

    birrp_stationdict['max_window_length'] = longest_section
    birrp_stationdict['n_bisections'] = number_of_bisections
    birrp_stationdict['coherence_threshold'] = coherence_threshold

    #self referencing:
    birrp_stationdict['rr_station'] = birrp_stationdict['station']

    birrp_stationdict = MTmc.add_birrp_simple_parameters_to_dictionary(birrp_stationdict)


    if output_channels == 2:
        birrp_stationdict['nout'] = 2

        inputstring = '0\n2\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n%i\n4,3,4\n%s\n0\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,input_filename)

    elif output_channels == 3:
        birrp_stationdict['nout'] = 3
        birrp_stationdict['nz'] = 2


        inputstring = '0\n3\n2\n2\n-%f\n%i,%i\ny\n0,0.999\n%f\n2\n%s\n0\n1\n3\n2\n0\n0\n0\n0\n0\n0\n0\n4,1,2,3,4\n%s\n0\n\i\n4,3,4\n%s\n0\n0,90,0\n0,90,0\n0,90,0\n'%(sampling_rate,longest_section,number_of_bisections,coherence_threshold,stationname,input_filename,length,stationname)

    string_file = op.join(ts_directory,'birrp_wd','birrp_input_string.txt')
    string_file = MTfh.make_unique_filename(string_file)
    with open(string_file,'w') as F:
        F.write(inputstring)
        F.write('\n')
    

    return inputstring, birrp_stationdict,string_file



def set_birrp_input_file_simple(stationname, ts_directory, output_channels, w_directory = '.'):
    """
    File handling: collect longest possible input for BIRRP from different files for the given station name. Generate a new input file in the working directory and return the name of this file, as well as time series and processing properties in form of a dictionary. 

    Scan all files in the directory by their headers: if the stationname matches, add the data file to the list.
    Additionally read out channels and start-/end times. Find longest consecutive time available on all channels.
    Then generate nx4/5 array for n consecutive time steps. Array in order Ex, Ey, Bx, By (,Bz) saved into file 'birrp_input_data.txt'

    Output_channels determine the number of output channels: 2 for Ex, Ey - 3 for additional Bz

    Output:
    - filename for birrp input data_in
    - length of data (samples)
    - sampling_rate (in Hz)
    - configuration dictionary for processed station

    """


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
        if not op.isfile(fn):
            continue


        try:
            header = MTfh.read_ts_header(fn)
        except:
            continue

        stationname_read = header['station']

        if not stationname_read == stationname.upper():
            continue
        if not header['channel'] in channels:
            continue

        lo_files.append(fn)

                
        lo_channels.append(header['channel'])
        lo_sampling_rates.append(float(header['samplingrate']))
        lo_starttimes.append(float(header['t_min']))
        ta = np.arange(int(float(header['nsamples']))+1)/float(header['samplingrate']) + float(header['t_min'])
        endtime = ta[-1]
        lo_endtimes.append(endtime)

    if (len(lo_sampling_rates) == 0) or (len(lo_files) == 0):
        sys.exit( '\n\tERROR - no MTpy data files found in directory {0} !!\n'.format(ts_directory))
        

    #take the most common sampling rate, if there are more than one:
    from collections import Counter
    tmp_dummy1 = lo_sampling_rates
    tmp_dummy2 = Counter(tmp_dummy1)
    sampling_rate = tmp_dummy2.most_common(1)[0][0]


    if not len(set(lo_channels)) in [4,5]:
        sys.exit( 'Missing data files in directory {0} - not all channels found'.format(ts_directory))

    #get a list with all existing time windows of consecutive data for all the channels
    lo_time_windows = []
    #find sorting of the files by their start time:
    starttime_sorting = np.argsort(lo_starttimes)

    print 'start loop over all components...'
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


    print 'find longest common time window for all channels...'
    longest_common_time_window = MTmc.find_longest_common_time_window_from_list(lo_time_windows, sampling_rate)


    #data array to hold time series for longest possible time window for the files given 
    #order Ex, Ey, Bx, By (,Bz)
    totalmemory, freememory = MTmc.show_memory()
    print '\nTotal memory available: {0} MB - free: {1} MB'.format(totalmemory,freememory)

    data = np.zeros((longest_common_time_window[2],output_channels+2))

    print 'Size of data array: {0} MB\n'.format(np.round(data.nbytes/1024.**2,2))
    #define time axis for referencing time and position in the output array
    #correct by rounding for internal floating point errors
    #ta = np.array([ np.round( i , -int(np.log10(1./sampling_rate))) for i in  np.linspace(*longest_common_time_window)])
    #alternative : no rounding
    ta = np.linspace(*longest_common_time_window, endpoint=False)

    print 'data array ({0}) and time axis ({1}) initialised...start looping over channels...'.format(data.shape, len(ta))

    for idx_ch, ch in enumerate(channels):
        for st in starttime_sorting:
            ch_read = lo_channels[st]
            if not ch == ch_read:
                continue

            sampling_rate_read = lo_sampling_rates[st]
            if not sampling_rate_read == sampling_rate:
                continue

            #read in data
            print 'reading data from file {0}'.format(lo_files[st])
            data_in = np.loadtxt(lo_files[st])
            print len(data_in), sampling_rate , lo_starttimes[st]
            #define time axis for read in data
            ta_file = np.arange(len(data_in))/sampling_rate + lo_starttimes[st]
            #find overlap of overall time axis and the ta of current data set:
            for min_idx,dummy in enumerate(ta_file):
                if ta[0]<= dummy <= ta[-1]:
                    break
            for max_idx in range(-1,-(len(ta_file)+1), -1):
                if ta[0]<=ta_file[max_idx] <= ta[-1]:
                    break
            #include the last sample:
            max_idx = len(ta_file) + max_idx + 1
            print min_idx,max_idx
            overlap = ta_file[min_idx:max_idx]
            print ta_file[0],ta_file[-1], '   ', ta[0],ta[-1]

            #find starting index of overlap for current data file time axis
            idx_ta_file = min_idx#np.argmin(np.abs(ta_file - overlap[0]))

            #find starting index of overlap for overall time axis
            idx_overall_ta = np.argmin(np.abs(ta - overlap[0])) 

            #set data entries
            print 'fill data from file into temporary array'
            print idx_overall_ta, len(overlap), idx_ch, idx_ta_file
            print data[idx_overall_ta:idx_overall_ta+len(overlap), idx_ch].shape, data_in[idx_ta_file:idx_ta_file+len(overlap)].shape
            data[idx_overall_ta:idx_overall_ta+len(overlap), idx_ch] = data_in[idx_ta_file:idx_ta_file+len(overlap)]

            gc.collect()

    #define output file for storing the output data array to:
    w_directory = op.abspath(op.join(os.curdir, w_directory))
    if not op.isdir(w_directory):
        os.makedirs(w_directory)
        print 'created directory: {0}'.format(w_directory)

    print 'size of usable data arry: ',data.shape

    try:
        outfn = op.join(w_directory, 'birrp_input_data.txt') 
        outfn = MTfh.make_unique_filename(outfn)
        print 'save input data array to file: {0}'.format(outfn)
        np.savetxt(outfn, data)
    except:
        raise MTex.MTpyError_file_handling('Error - cannot write data to file:{0}'.format(outfn))

    print 'Done...'

    birrp_stationdict = {}
    birrp_stationdict['station'] =  stationname.upper()
    birrp_stationdict['n_output_channels'] = output_channels
    birrp_stationdict['sampling_rate'] = sampling_rate
    birrp_stationdict['n_samples'] =  len(data)
    birrp_stationdict['processing_window_start'] = longest_common_time_window[0]
    birrp_stationdict['recording_start'] =  lo_starttimes[0]

    gc.collect()

    return op.abspath(outfn), birrp_stationdict['n_samples'] , sampling_rate, birrp_stationdict


def get_optimal_window_bisection(length, sampling_rate):
    """

    input:
    - data series length in samples
    - sampling rate in Hz

    """

    # longest time window cannot exceed 1/30 of the total length in order to obtain statistically sound results (maybe correcting this to 1/100 later); value given in samples
    longest_window = int(length/30.)

    #shortest_window cannot be smaller than 2^4=16 times the sampling interval in order to suit Nyquist and other criteria; value given in samples
    shortest_window = 16 


    #find maximal number of bisections so that the current window is still longer than the shortest window allowed
    number_of_bisections = int(np.ceil(np.log(16./longest_window) / np.log(0.5) ))

    return longest_window, number_of_bisections



def write_script_file(processing_dict, save_path=None):
    """
    writeScriptfile(processingdict will write a script file for BIRRP using 
    info in processingdict which is a dictionary with keys:
    
    ================== ========================================================
    parameter          description
    ================== ======================================================== 
    station            station name
    fn_lst             list of file names to be processed
    rrfn_lst           list of remote reference file names
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
    jmode              input file mode (0=user defined; 1=start time 
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
                            os.path.join(os.path.dirname(fn_lst[0]),'BF')
        
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
    pdict = processing_dict
    
    try:
        fn_array = np.array(pdict['fn_lst'])
    except KeyError:
        raise KeyError('fn_lst --> Need to input a list of files to process')
    
    try:
        nds, ndf = fn_array.shape
    except ValueError:
        ndf = fn_array.shape[0]
        nds = 0
    if save_path is None:
        if nds == 0 or nds == 1:
            bfpath = os.path.join(os.path.dirname(pdict['fn_lst'][0]),
                                 'BF')
        else:
            bfpath = os.path.join(os.path.dirname(pdict['fn_lst'][0][0]),
                                  'BF')
    else:
        bfpath = save_path
    
    
    if nds == 0:
        npcs=1
    elif nds==1:
        nds=0
        npcs=1
        pdict['fn_lst'] = pdict['fn_lst'][0]
        try:
            pdict['rrfn_lst'] = pdict['rrfn_lst'][0]  
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
        
    elif ilev==1:
        print 'Writing Advanced mode'
        fid.write('{0:d} \n'.format(ilev))
        fid.write('{0:d} \n'.format(nout))
        fid.write('{0:d} \n'.format(ninp))
        fid.write('{0:d} \n'.format(nref))
        if nref>3:
            nrrlst=np.array([len(rrlst) 
                            for rrlst in pdict['rrfn_lst']])
            nr3=len(np.where(nrrlst==3)[0])
            nr2=len(np.where(nrrlst==2)[0])
            fid.write('{0:d},{1:d} \n'.format(nr3,nr2))
        fid.write('{0:d} \n'.format(nrr))
        #if remote referencing
        if int(nrr)==0:
            fid.write('{0:.3f} \n'.format(tbw))
            fid.write('{0:.3f} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshe))
            #parameters for bz component if ninp=3
            if nout==3:
                if c2threshe!=0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                else:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                if c2threshe1!=0.0 or c2threshe!=0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            else:
                if c2threshe!=0.0:
                    fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej=[prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        #if 2 stage processing
        elif int(nrr)==1:
            fid.write('{0:.5g} \n'.format(tbw))
            fid.write('{0:.5g} \n'.format(deltat))
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nfft,nsctinc,nsctmax))        
            fid.write('{0:d},{1:.2g},{2:d} \n'.format(nf1,nfinc,nfsect))
            fid.write('y \n')
            fid.write('{0:.2g} \n'.format(mfft))        
            fid.write('{0:.5g},{1:.5g},{2:.5g} \n'.format(uin,ainlin,ainuin))
            fid.write('{0:.3f} \n'.format(c2threshb))        
            fid.write('{0:.3f} \n'.format(c2threshe))
            if nout==3:
                if c2threshb!=0 or c2threshe!=0:
                    fid.write('{0:d} \n'.format(nz))
                    fid.write('{0:.3f} \n'.format(c2threshe1))
                elif c2threshb==0 and c2threshe==0:
                    fid.write('{0:d} \n'.format(0))
                    fid.write('{0:.3f} \n'.format(0))
            if c2threshb!=0.0 or c2threshe!=0.0:
                fid.write('{0:.6g},{1:.6g} \n'.format(perlo,perhi))
            fid.write(ofil+'\n')
            fid.write('{0:d} \n'.format(nlev))
            fid.write('{0:d} \n'.format(nprej))
            if nprej!=0:
                if type(prej) is not list:
                    prej=[prej]
                fid.writelines(['{0:.5g} \n'.format(nn) for nn in prej])
        
    fid.write('{0:d} \n'.format(npcs))    
    fid.write('{0:d} \n'.format(nar))    
    fid.write('{0:d} \n'.format(imode))    
    fid.write('{0:d} \n'.format(jmode))    
    #write in filenames 
    if npcs != 1:
        if jmode == 0:
            fid.write(str(nread[0])+'\n')
            #--> write filenames to process with other information for first
            #    time section
            for tt, tfile in enumerate(pdict['fn_lst'][0]):
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
                    if hx_cal is not None:
                        fid.write('-2\n')
                        fid.write(hz_cal+'\n')
                    else:
                        fid.write(str(nfil)+'\n')
                else:
                    fid.write(str(nfil)+'\n')
                fid.write(tfile+'\n')
                fid.write(str(nskip[0])+'\n')
                
            #--> write remote reference time series
            for rr, rfile in enumerate(pdict['rrfn_lst'][0]):
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
                for tfile in pdict['fn_lst'][nn]:
                    fid.write(tfile+'\n')
                    fid.write(str(nskip[0])+'\n')
                for rfile in pdict['rrfn_lst'][nn]:
                    fid.write(rfile+'\n')
                    fid.write(str(nskipr[nn])+'\n')
                    
        #--> if start and end time are give write in those
        elif jmode == 1:
            #write filenames
            for tt, tfile in enumerate(pdict['fn_lst'][0]):
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
                    if hx_cal is not None:
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
            for rr, rfile in enumerate(pdict['rrfn_lst'][0]):
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
                for tfile in pdict['fn_lst'][nn]:
                    fid.write(tfile+'\n')
                    fid.write(dstim+'\n')
                    fid.write(wstim+'\n')
                    fid.write(wetim+'\n')
                for rfile in pdict['rrfn_lst'][nn]:
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
                for tt, tfile in enumerate(pdict['fn_lst']):
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
                        if hx_cal is not None:
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
                for rr, rfile in enumerate(pdict['rrfn_lst']):
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
                for tfile in pdict['fn_lst'][0]:
                    fid.write(str(nfil)+'\n')
                    fid.write(tfile+'\n')
                    if type(nskip) is list:
                        fid.write(str(nskip[0])+'\n')
                    else:
                        fid.write(str(nskip)+'\n')
                for rfile in pdict['rrfn_lst'][0]:
                    fid.write(str(nfil)+'\n')
                    fid.write(rfile+'\n')
                    if type(nskipr) is list:
                        fid.write(str(nskipr[0])+'\n')
                    else:
                        fid.write(str(nskipr)+'\n')
                        
        elif jmode == 1:
            #write filenames
            if nds==0:
                for tt, tfile in enumerate(pdict['fn_lst']):
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
                        if hx_cal is not None:
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
                for rr, rfile in enumerate(pdict['rrfn_lst']):
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
                for tt, tfile in enumerate(pdict['fn_lst'][0]):
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
                        if hx_cal is not None:
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
                for rr, rfile in enumerate(pdict['rrfn_lst'][0]):
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
    j_filename_list = [i for i in os.listdir(input_dir) if op.basename(i).upper() == ('%s.j'%stationname).upper() ]
    try:
        j_filename = op.join(input_dir, j_filename_list[0])
    except:
        raise MTex.MTpyError_file_handling('j-file for station %s not found in directory %s'%(stationname, input_dir))
    
    if len(j_filename_list) > 1:
        raise MTex.MTpyError_file_handling('More than one j-file for station %s found in directory %s'%(stationname, input_dir))

    #Having now:
    # station_config_dict - contains information about station setup
    # birrp_config_dict - contains information about the processing (BIRRP parameters, selected time window, Rem.Ref.,..)
    # directory - contains BIRRP output files, coded by stationname

    # To be converted into .EDI
    # Dictionaries information goes into EDI header: HEAD and INFO section - check for other sections though
    # output EDI file is out_fn
      
    periods, Z_array, tipper_array = read_j_file(j_filename)


    HEAD = _set_edi_head(station_config_dict,birrp_config_dict)

    INFO = _set_edi_info(station_config_dict,birrp_config_dict)

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
                                        instr_response_file, out_dir = None):
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
    try:
        birrp_config_dict = MTcf.read_configfile(birrp_configfile)
    except:
        raise MTex.MTpyError_config_file( 'Config file with BIRRP processing parameters could not be read: %s' % (birrp_configfile) )

    #find the birrp-output j-file for the current station 
    j_filename_list = [i for i in os.listdir(input_dir) if op.basename(i).upper() == ('%s.j'%stationname).upper() ]
    try:
        j_filename = j_filename_list[0]
    except:
        raise MTex.MTpyError_file_handling('j-file for station %s not found in directory %s'%(stationname, input_dir))
    
    if len(j_filename_list) > 1:
        raise MTex.MTpyError_file_handling('More than one j-file for station %s found in directory %s'%(stationname, input_dir))

    #Having now:
    # station_config_dict - contains information about station setup
    # birrp_config_dict - contains information about the processing (BIRRP parameters, selected time window, Rem.Ref.,..)
    # directory - contains BIRRP output files, coded by stationname

    # To be converted into .EDI
    # Dictionaries information goes into EDI header: HEAD and INFO section - check for other sections though
    # output EDI file is out_fn
    
    periods, Z_array, tipper_array = read_j_file(j_filename)

    frequencies = 1./periods
    
    def correct_z_for_instrument_response(Z_array, instr_resp, frequencies):
        
        for idx_f, freq in enumerate(frequencies):
            if not (instr_resp[0,0] <= np.abs(freq) <= instr_resp[-1,0]):
                print 'no instrument response in this frequency range - array values set to zero here: ', freq
                correction_factor = 0.
                continue

            #find the appropriate frequencies ( since the current freq-value is most likely inbetween two values on the instr_freqs-axis) - get the respective value by linear interpolation

            #find the value closest to the current freq, assume it's lower
            closest_lower = np.abs(freq-instr_resp[:,0]).argmin()
            
            #if it coincides with the highest frequency/last entry:
            if closest_lower == len(instr_resp)-1:
                correction_factor = np.complex(instr_resp[-1,1],instr_resp[-1,2])
            # or the lowest
            elif closest_lower == 0:
                correction_factor = np.complex(instr_resp[0,1],instr_resp[0,2])
            else:
                #in case the closest frequency value is not lower but higher, take the freq value below as lower bound for the interval:        
                if instr_resp[closest_lower,0] > freq:
                    closest_lower -= 1
                
                #define the interval:
                instrfreq1 = instr_resp[closest_lower,0]
                instrfreq2 = instr_resp[closest_lower+1,0]
                instrfactor1 = np.complex( instr_resp[closest_lower,1]  , instr_resp[closest_lower,2] )
                instrfactor2 = np.complex( instr_resp[closest_lower+1,1], instr_resp[closest_lower+1,2])
                intervallength = instrfreq2 - instrfreq1
                weight = (freq-instrfreq1)/intervallength

                correction_factor = weight * instrfactor2 + (1-weight) * instrfactor1
                #lo_freqs.append([freq,instrfreq1,instrfreq2,weight,instrfactor1,instrfactor2,factor])
            
            #finally correct Z for the instrument influence by multiplying with the instrument response value: 
            for i in range(4):
                zentry = np.complex(Z_array[idx_f,0,i], Z_array[idx_f,1,i] )
                corrected = zentry * correction_factor
                Z_array[idx_f,0,i] = np.real(corrected)
                Z_array[idx_f,1,i] = np.imag(corrected)
                #correct the error: stretch by the ratio of amplitudes of original and correct Z value
                Z_array[idx_f,2,i] = Z_array[idx_f,2,i]/np.abs(zentry)*np.abs(corrected)


        return Z_array


    Z_array = correct_z_for_instrument_response(Z_array, instr_resp, frequencies)



    HEAD = _set_edi_head(station_config_dict,birrp_config_dict)
    INFO = _set_edi_info(station_config_dict,birrp_config_dict)
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
    #datastring += '>!****FREQUENCIES****!\n'
    datastring += '>FREQ nfreq=%i // %i\n'%(len(periods),len(periods))
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
            datastring += '>%s%s // %i\n'%(compstrings[Z_comp], Z_entries[entry], len(periods))
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
            datastring += '>%s%s // %i\n'%(compstrings[T_comp], T_entries[entry], len(periods))
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


def _set_edi_info(station_config_dict,birrp_config_dict):

    infostring = ''
    infostring += '>INFO\t max lines=1000\n'
    infostring += '\tStation parameters\n'

    for key in sorted(station_config_dict.iterkeys()):
        infostring += '\t\t{0}: {1}  \n'.format(str(key),
                                                str(station_config_dict[key]))   
    

    infostring += '\n'
    infostring += '\tProcessing parameters\n'

    for key in sorted(birrp_config_dict.iterkeys()):
        infostring += '\t\t{0}: {1}  \n'.format(str(key),
                                                str(birrp_config_dict[key]))   
    infostring += '\n'    
    if len(birrp_config_dict) == 0 :
        infostring += '\t\tunknown\n\n'


    infostring += '\n'

    return infostring.expandtabs(4)


def _set_edi_head(station_config_dict,birrp_config_dict):
    """
    set header string
    
    set date to format YYYY-MM-DD,HH:MM:SS
    """
    frmt = '%Y-%m-%d,%H:%M:%S'

    headstring = ''
    headstring += '>HEAD\n\n'
    if len((station_config_dict['station'].split())) != 1: 
        headstring += '\tdataid="%s"\n'%(station_config_dict['station'])
    else:
        headstring += '\tdataid={0}\n'.format(station_config_dict['station'])


    if station_config_dict.has_key('company'):
        acqby = station_config_dict.has_key('company')
    else:
        acqby = ''

    if len(acqby.split()) != 1:
        headstring += '\tacqby="%s"\n'%(acqby)
    else:
        headstring += '\tacqby={0}\n'.format(acqby)


    if len(birrp_config_dict) !=0 :
        sampling_rate = float(birrp_config_dict['sampling_rate'])
        try:
            n_samples = int(birrp_config_dict['n_samples'])
        except ValueError:
            n_samples = sum([int(ii) for ii in 
                             birrp_config_dict['n_samples'][1:-1].strip().split(',')])
        #new:
        try:
            acq_starttime = float(birrp_config_dict['processing_window_start'])
            dt_start = datetime.datetime.fromtimestamp(acq_starttime)
            acq_start = dt_start.combine(dt_start.date(),dt_start.time()).strftime(frmt+'.%f')

            dt_end =  acq_starttime + 1./sampling_rate*(n_samples) 
            acq_end = dt_end.combine(dt_end.date(),dt_end.time()).strftime(frmt+'.%f')
            
            headstring +='\tacqdate=%s \n'%(acq_start)
            headstring +='\tenddate=%s \n'%(acq_end)
        except KeyError:
            try:
                acq_start = station_config_dict.has_key('acq_date')
                headstring += '\tacqdate={0} \n'.format(acq_start)
            except KeyError:
                 headstring += '\tacqdate={0} \n'.format('1970-01-01')
                

    todaystring = datetime.datetime.now().strftime(frmt)
    headstring += '\tfiledate=%s\n'%(todaystring)


    network = ''
    if station_config_dict.has_key('network'):
        location = station_config_dict.has_key('network')
    headstring += '\tprospect="%s"\n'%(network)

    location = ''
    if station_config_dict.has_key('location'):
        location = station_config_dict.has_key('location')
    headstring += '\tloc="%s"\n'%(location)

    headstring += '\tlat=%.5f\n'%station_config_dict['latitude']
    headstring += '\tlong=%.5f\n'%station_config_dict['longitude']
    headstring += '\telev=%.1f\n'%station_config_dict['elevation']

    headstring += '\n'

    return headstring.expandtabs(4)


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
    #dmeasstring += '\treflat=%f\n'%station_config_dict['latitude']
    #dmeasstring += '\treflong=%f\n'%station_config_dict['longitude']
    #dmeasstring += '\trefelev=%.1f\n'%station_config_dict['elevation']
    
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
    """   

    j_fn = op.abspath(fn)
    if not op.isfile(j_fn):
        raise MTex.MTpyError_inputarguments('Cannot read j-file %s - file is not existing'%(j_fn))

    
    
    with open(j_fn,'r') as F_in:
        j_lines = F_in.readlines()
    
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

    return _check_j_file_content(periods, Z, tipper)

 

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
    lo_all_periods_raw.sort()
    lo_all_periods = np.array(lo_all_periods_raw)

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
    cohfilenames = [ op.abspath(op.join(directory,i)) for i in fnmatch.filter(
                os.listdir(directory), '*%s*.[123]r.2c2'%stationname.upper()) ] 
    
    if len(cohfilenames) < 1:
        print 'No coherence files for station %s found in: %s'%(stationname, directory)
        raise MTex.MTpyError_file_handling()#'No coherence files for station %s found in: %s'%(stationname, directory))


    if len(cohfilenames) > 3:
        print 'Too many coherence files for station %s found in: %s'%(stationname, directory)
        raise MTex.MTpyError_file_handling()#'Too many coherence files for station %s found in: %s'%(stationname, directory))

    try:
        period,freq,coh1,zcoh1 = MTfh.read_2c2_file(cohfilenames[0])
        period,freq,coh2,zcoh2 = MTfh.read_2c2_file(cohfilenames[1])

        if len(cohfilenames) == 3:

            period,freq,coh3,zcoh3 = MTfh.read_2c2_file(cohfilenames[2])
    except:
        print 'Cannot read coherence files for station %s found in: %s'%(stationname, directory)
        raise MTex.MTpyError_file_handling()#'Cannot read coherence files for station %s found in: %s'%(stationname, directory))

    fn = '%s.coh'%(stationname)
    out_fn = op.abspath(op.join(directory,fn))
    out_fn = MTfh.make_unique_filename(out_fn)


    F_out =  open(out_fn,'w')
    F_out.write('period \t freq \t coh1 \t zcoh1 \t coh2 \t zcoh2 \t coh3 \t zcoh3 \n'.expandtabs(4))


    for ff in range(len(period)):
        try:
            c1 = float(coh1[ff])
            zc1 = float(zcoh1[ff])
        except :
            c1= 0.
            zc1= 0.
        
        try:
            c2 = float(coh2[ff])
            zc2 = float(zcoh2[ff])
        except :
            c2 = 0.
            zc2 = 0.
        
        c3 = 0.
        zc3 = 0.

        if len(cohfilenames) == 3:
  
            try:
                c3 = float(coh3[ff])
                zc3 = float(zcoh3[ff])
            except:
                pass
    
                   
        F_out.write(('%f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n'%(period[ff], freq[ff], c1, zc1, c2, zc2, c3, zc3)).expandtabs(4))
    F_out.close()

    return out_fn
