#!/usr/bin/env python
"""
This is a convenience script for the removal of instrument response from a set of MTpy time series data files within a directory (non-recursive). The data files have to contain a MTpy style header line, which specifies station, channel, timestamps.

It needs the location of the directory and the location of the instrument response file. 
The latter has to consist of an array with three columns: frequencies, real, imaginary 

If no output folder is specified, a subfolder 'instr_resp_corrected' is set up within the input directory

"""

import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


import mtpy.utils.exceptions as EX
import mtpy.processing.calibration as PC
import mtpy.utils.filehandling as FH
import mtpy.processing.instrument as PI

reload(FH)
reload(EX)
reload(PC)
reload(PI)



def main():

    if len(sys.argv) < 3:
        raise EX.MTpyError_inputarguments('Need at least 2 arguments: <path to files> <response file> [<output dir>] [<channel(s)>] ')


    pathname_raw = sys.argv[1] 
    directory = op.abspath(op.realpath(pathname_raw))

    responsefilename_raw = sys.argv[2]
    responsefile = op.abspath(op.realpath(responsefilename_raw))


    if not op.isdir(directory):
        raise EX.MTpyError_inputarguments('Directory not existing: %s' % (directory))

    if not op.isfile(responsefile):
        raise EX.MTpyError_inputarguments('Response file not existing: %s' % (responsefile))
    
    #check, if response file is in proper shape (3 columns freq,re,im of real values):
    try:
        responsedata = np.loadtxt(responsefile)
        s = responsedata.shape
        if s[1] != 3:
            raise
        freq_min = responsedata[0,0]
        freq_max = responsedata[-1,0]

    except: 
        raise EX.MTpyError_inputarguments('Response file (%s) in wrong format - must be 3 columns: freq,real,imag' % (responsefile))

    #set up output directory: 
    try:
        outdir_raw = sys.argv[3]
        outdir = op.abspath(outdir_raw)
    except:
        outdir = op.join(directory,'instr_resp_corrected')

    try:
        if not op.isdir(outdir):
            os.makedirs(outdir)
    except:
        raise EX.MTpyError_inputarguments('Output directory cannot be generated: %s' % (outdir))

    #define channels to be considered for correction:
    try:
        lo_channels = list(set([i.upper() if len(i)==2 else 'B'+i.upper() for i in  sys.argv[4].split(',')]))
    except:
        print 'No channel list found - using BX, BY, HX, HY'
        lo_channels = ['BX', 'BY', 'HX', 'HY', 'BZ', 'HZ']


    #collect file names  within the folder 
    oldwd = os.getcwd()
    os.chdir(directory)
    lo_allfiles = glob.glob('*')
    lo_allfiles = [op.abspath(i)  for i in lo_allfiles if op.isfile(i)==True]
    os.chdir(oldwd)

    #generate list of list-of-files-for-each-channel:
    lo_lo_files_for_channels = [[] for i in lo_channels  ]

    #check the files for information about the determined channels:
    for fn in lo_allfiles:
        header_dict = FH.read_ts_header(fn)
        if len(header_dict.keys()) == 0 :
            continue

        ch = header_dict['channel'].upper()
        try:
            ch_idx = lo_channels.index(ch)
        except ValueError:
            continue

        # use the current file, if it contains a header line and contains signal from the requested channel:
        lo_lo_files_for_channels[ch_idx].append(fn)

    #if no files had header lines or did not contain data from the appropriate channel(s):
    if np.sum([len(i) for i in lo_lo_files_for_channels]) == 0:
        print 'channels: ', lo_channels, ' - directory: ',directory
        raise EX.MTpyError_inputarguments('No information for channels found in the directory - Check header lines!')

    #=============================================
    # start the instrument correction
        
    # looping over all requested channels:
    for ch in lo_channels:
        #skip, if no data are available for the current channel:
        if [len(i) for i in lo_lo_files_for_channels][lo_channels.index(ch)] == 0:
            continue

        #set up lists for the infos needed later, esp. for the file handling
        lo_files = lo_lo_files_for_channels[lo_channels.index(ch)]
        lo_t_mins = []
        lo_headers = []

        #read in header lines and sort files by increasing starttimes t_min
        for fn in lo_files:
            header_dict = FH.read_ts_header(fn)
            lo_t_mins.append(header_dict['t_min'])
            lo_headers.append(header_dict)

        #sort all the collected lists by t_min
        idxs = np.array(lo_t_mins).argsort()

        lo_t_mins = [lo_t_mins[i] for i in idxs]
        lo_files = [lo_files[i] for i in idxs]
        lo_headers = [lo_headers[i] for i in idxs]
           

        # finding consecutive, continuous time axes:
        lo_timeaxes = []
        ta_old = None

        for idx, header in enumerate(lo_headers):
            ta_cur = np.arange(int(header['nsamples']))/float(header['samplingrate']) + float(header['t_min'])

            #if there is no old ta:
            if ta_old == None:
                ta_old = ta_cur 
                continue

            # if gap between old and new ta is too big:
            if (ta_cur[0] - ta_old[-1]) > (2*1./float(header['samplingrate'])):
                lo_timeaxes.append(np.array(ta_old))
                ta_old = ta_cur 
                continue

            #find index of new ta which is closest to the end of old_ta - most commonly it's '0' !
            overlap = np.abs(ta_cur - ta_old[-1]).argmin()
            ta_cur = ta_cur[overlap:]
            ta_old = np.concatenate([ta_old,ta_cur])
        
        #append last active time axis ta:
        lo_timeaxes.append(np.array(ta_old))

        #determine maximal period from response file and existinng time axes. 
        #win = get_windowlength() = max([ i[-1]-i[0] for i in lo_timeaxes] ) 
        # the minimum of the maximal resolvable signal period and the longest continuous time axis:
        winmax = 1./freq_min
        #for debugging set large window size:
        winmax = 10e7
        #later on, if the TS is longer than 3 times this time window, we want to cut out subsections of the time series. These cuts shall consist of triplets of subwindows, each of which shall not be longer than this maximum period.

        #Now the data set has to be corrected/deconvolved by looping over the collected time axes:
        for ta in lo_timeaxes:
            print '\nhandling time axis: {0} - {1} ({2} samples) '.format(ta[0],ta[-1],len(ta))

            #if the time section is shorter than 3 times the maximum defined by the response function, read in the whole data set at once for this interval
            if (ta[-1] - ta[0]) < (3 * winmax) : 
                print 'time axis short enough - reading all at once'

                #collect the appropriate files in a list
                #after the MTpy preprocessing the start end end of the time series coincide with files start and endings, so no "half files" are involved. 
                cur_time = ta[0]
                data = []
                files = []
                headers = []
                starttimes = []

                while cur_time < ta[-1]:
                    for idx,header in enumerate(lo_headers):
                        ta_cur = np.arange(int(header['nsamples']))/float(header['samplingrate']) + float(header['t_min'])
                        if cur_time in ta_cur:
                            start_idx = np.where(ta_cur == cur_time)[0][0]
                            break
                    fn = lo_files[idx]
                    files.append(fn)
                    headers.append(header)
                    starttimes.append(float(header['t_min']))
                    cur_data = np.loadtxt(fn)

                    print 'current data section length: ',len(cur_data)
                    if ta_cur[-1] <= ta[-1]:
                        data.extend(cur_data[start_idx:].tolist())
                        cur_time = ta_cur[-1] + 1./float(header['samplingrate'])  
                    else:
                        end_idx = where(ta_cur == ta[-1])[0][0] 
                        data.extend(cur_data[start_idx:end_idx+1].tolist())
                        cur_time = ta[-1]
                    print 'current data length: ',len(data)

                #at this point, the data set should be set up for the given time axis
                corrected_timeseries = PI.correct_for_instrument_response(np.array(data),float(header['samplingrate']), responsedata)  

                print 'corrected TS starting at {0}, length {1}'.format(ta[0],len(corrected_timeseries))

                #now, save this TS back into the appropriate files, including headers
                for idx,fn in enumerate(files):

                    # output file name: use input file name and append '_true'
                    inbasename = op.basename(fn)
                    outbasename = ''.join([op.splitext(inbasename)[0]+'_true',op.splitext(inbasename)[1]])
                    outfn = op.join(outdir,outbasename)

                    outF = open(outfn,'w')
                    header = headers[idx]
                    unit = header['unit']
                    if unit[-6:].lower() != '(true)':
                        unit +='(true)'
                    header['unit'] = unit
                    headerline = FH.get_ts_header_string(header)
                    outF.write(headerline+'\n')
                    starttime = starttimes[idx]
                    length = int(float(header['nsamples']))
                    
                    startidx = (np.abs(starttime - ta)).argmin()
                    print startidx,length,len(corrected_timeseries),len(ta)
                    print 'handling file {0} - starttime {1}, - nsamples {2}'.format(outfn,starttime,length)
                    print outdir,outfn
                    data = corrected_timeseries[startidx:startidx+length]
                    np.savetxt(outF,data)
                    outF.close()


                #To do so, use the time axis and run over the input files again,determine the filenames. Use them, and put them (slightly modified, perhaps?) into the given output directory  
                #return corrected_timeseries

            else:
                #find partition into pieces of length 'winmax'. the remainder is equally split between start and end:

                #total time axis length:
                ta_length = ta[-1] - ta[0]
                
                #partition into winmax long windows 
                n_windows = int(ta_length/winmax)
                remainder = ta_length%winmax
                lo_windowstarts = [ta[0]]
                for i in range(n_windows+1):
                    t0 = ta[0] + remainder/2. + i * winmax
                    lo_windowstarts.append(t0)

                # loop over the winmax long sections:
                for idx_t0, t0 in enumerate(lo_windowstarts):
                    #for each step (except for he last one obviously), 3 consecutive parts are read in, concatenated and deconvolved. Then the central part is taken as 'true' data. 
                    #only for the first and the last bit (start and end pieces) are handled together with the respective following/preceding section

                    #a list of input file(s) containing the data:
                    lo_input_files = []
                    lo_output_files = []
                    #the data currently under processing:
                    data = []
                    #since this loop is working as a moving window, there is no need to re-read a section 3 times:
                    old_data = []

                    1. find the 3 data pieces in the files 
                    2. store the file name(s), perhaps including their respective time axes 
                    3. define 'data'
                    4. deconvolve 'data'
                    5. cut the middle section
                    6. find, from which files that came
                    7. write it to the equivalent output file(s)
                    8. for the first and last section: pre/append them to the middle section and write out the whole bit 
















            


if __name__=='__main__':
    main()