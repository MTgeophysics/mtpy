#!/usr/bin/env python
"""
This is a convenience script for the removal of instrument response from a set of MTpy time series data files within a directory (non-recursive). The data files have to contain a header line, which specifies station, channel, timestamps.

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
reload(C)



def main():

    if len(sys.argv) < 3:
        raise EX.MTpyError_inputarguments('Need at least 2 arguments: <path to files> <config file> [<output dir>] [<channel(s)>] ')


    pathname_raw = sys.argv[1] 
    directory = op.abspath(op.realpath(pathname_raw))

    responsefilename_raw = sys.argv[2]
    responsefile = op.abspath(op.realpath(configfilename_raw))


    if not op.isdir(directory):
        raise EX.MTpyError_inputarguments('Directory not existing: %s' % (directory))

    if not op.isfile(responsefile):
        raise EX.MTpyError_inputarguments('Response file not existing: %s' % (responsefile))
    
    #check, if response file is in proper shape (3 columns freq,re,im of real values):
    try:
        responsedata = np.loadtxt(responsefile)
        s = dd.shape
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
         lo_channels = list(set([i.upper() if len(i)==2 else 'B'+i.upper() for i in  sys.argv[3].split(',')]))
     except:
        print 'No channel list found - using BX, BY, HX, HY'
        lo_channels = ['BX', 'BY', 'HX', 'HY', 'BZ', 'HZ']


    #collect file names with information about the determined channels within the folder 
    oldwd = os.getcwd()
    os.chdir(directory)
    lo_allfiles = glob.glob('*')
    lo_allfiles = [op.abspath(i) for i in lo_allfiles]
    os.chdir(oldwd)

    #generate list of lists-of-files-for-each-channel:
    lo_lo_files_for_channels = [[] for i in lo_channels  ]

    for fn in lo_allfiles:
        header_dict = FH.read_header(fn)
        if len(header_dict.keys()) == 0 :
            continue

        ch = header_dict['channel'].upper()
        try:
            ch_idx = lo_channels.index(ch)
        except ValueError:
            continue
        #use file, if it contains a header line and contains signal from the requested channel:
        lo_lo_files_for_channels[ch_idx].append(fn)

    #if no files had header lines or did not contain data from the appropriate channel(s):
    if np.sum([len(i) for i in lo_lo_files_for_channels]) == 0:
        print 'No information for channels found in the directory - Check header lines!'
        print 'channels: ', lo_channels, ' - directory: ',directory
        sys.exit()

    #now dealing with the data of interest
    for ch in lo_channels:
        if [len(i) for i in lo_lo_files_for_channels][lo_channels.index(ch)] == 0:
            continue
        lo_files = lo_lo_files_for_channels[lo_channels.index(ch)]
        lo_t_mins = []
        
        lo_headers = []
        lo_data = []
        #sorting files by increasing t_min
        for fn in lo_files:
            header_dict = FH.read_header(fn)
            lo_t_mins.append(header_dict['t_min'])
            lo_headers.append(header_dict)
            lo_data.append(np.loadtxt(fn))

        idxs = np.array(lo_t_mins).argsort()
        #sort the collections by t_min
        lo_t_mins = [lo_t_mins[i] for i in idxs]
        lo_files = [lo_files[i] for i in idxs]
        lo_headers = [lo_headers[i] for i in idxs]
        lo_data = [lo_data[i] for i in idxs]
           
        # finding consecutive time axes:
        lo_timeaxes = []
        ta_old = None

        for idx, header in lo_headers:
            ta_cur = np.arange(int(header['nsamples']))/float(header['samplingrate']) + float(header['t_min'])

            #if there is no old ta:
            if ta_old == None:
                ta_old = ta_cur 
                continue

            # if gap between old and new ta is too big:
            if ta_cur[0] - ta_old[-1] > 2*1./float(header['samplingrate'])
                lo_timeaxes.append(ta_old)
                ta_old = ta_cur 
                continue

            #find index of new ta which is closest to the end of old_ta - most commonly it's 0 !
            overlap = np.abs(ta_cur - ta_old[-1]).argmin()
            ta_cur = ta_cur[overlap:]
            ta_old = np.concatenate([ta_old,ta_cur])
        
        #append last active ta:
        lo_timeaxes.append(ta_old)

        #determine maximal period from response file. 
        maxwindowlength = max([ i[-1]-i[0] for i in lo_timeaxes] ) 
        #later on, if the TS is longer than 4 times this time window, we want to cut out subsections of the time series. These cuts shall consist of triplets of subwindows, each of which shall not be longer than this maximum period.

        #Now prepare the data set to be corrected by looping over the time axes:

        for ta in lo_timeaxes:

            #if the time section is shorter than 4 times the maximum defined by the response function, read in the whole data set at once for this interval
            if (ta[-1] - ta[0]) < 4*1./freq_min : 
                #use full dataset belonging to this time axis
                cur_time = ta[0]
                data = []
                while cur_time < ta[-1]:
                    for idx,header in enumerate(lo_headers):
                        ta_cur = np.arange(int(header['nsamples']))/float(header['samplingrate']) + float(header['t_min'])
                        if cur_time in ta_cur:
                            start_idx = where(ta_cur == cur_time)[0][0]
                            break
                    if ta_cur[-1] <= ta[-1]:
                        cur_data = lo_data[idx]
                        data.extend(cur_data[start_idx:].tolist())
                        cur_time = ta_cur[-1] + 1./float(header['samplingrate'])  
                    else:
                        end_idx = where(ta_cur == ta[-1])[0][0] 
                        cur_data = lo_data[idx]
                        data.extend(cur_data[start_idx:end_idx+1].tolist())
                        cur_time = ta[-1]
                #at this point, the data set should be set up for the given time axis

                corrected_timeseries = PI.correct_for_instrument_resonse(np.array(data),float(header['samplingrate']), responsedata)   

                #now, save this TS back into the appropriate files, including headers
                #To do so, use the time axis and run over the input files again,determine the filenames. Use them, and put them (slightly modified) into the given output directory  







            


        #cut 







    #determine time window length

    #----------------------------------------------------------------------------
    #to be improved later - rather rely on header lines than filenames!! :

    #select files by suffix, since header information is not necessarily present
    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    oldwd = os.getcwd()
    os.chdir(directory)
    lo_allfiles = glob.glob('*.??')
    lo_allfiles = [op.abspath(i) for i in lo_allfiles]
    os.chdir(oldwd)

    lo_files = []

    for f in lo_allfiles:
        if f[-2:].lower() in components:
            lo_files.append(f)

    #check, if list of files is empty
    if len(lo_files) == 0:
        raise EX.MTpyError_inputarguments('Directory does not contain files to calibrate: %s' % (wd))
    #-------------------------------------------------------


    for f in lo_files:

        #find station 
        #try reading in a potentially existing header line
        try:
            F = open(f,'r')
            firstline = F.readline()
            F.close()
            firstlinesplit = firstline.strip().split()
            if firstlinesplit[0][0] == '#':
                #check for missing whitespace after commenting symbol #:
                if len(firstlinesplit[0]) == 1:
                    stationname = firstlinesplit[1].upper()
                    channel = firstlinesplit[2].lower()
                #otherwise take the rest of the first string as stationname    
                else:
                    stationname = firstlinesplit[0][1:].upper()
                    channel = firstlinesplit[1].lower()
            else:
                raise

        except:
            try:

                stationname = FH.EDL_get_stationname_fromfilename(f)
                channel = f[-2:].lower()
            except:
                print 'stationname or channel for file %s could not be determined - skipping file'%(f)
                continue
        

        #get configuration dictionary for this station
        try:
            stationdict = config_dir[stationname]
        except:
            print 'no entry for station %s found in configuration file %s skipping file'%(stationname, configfile )
            continue

        latitude = stationdict['latitude']
        longitude = stationdict['longitude']
        elevation = stationdict['elevation']

        field = channel[0]
        direction = channel[1]


        if field == 'e':
            #check North-South axis orientation
            if direction == 'x':
                angle = float(stationdict['e_xaxis_azimuth'])
                dipolelength = float(stationdict['e_xaxis_length'])
                if np.abs(180. - angle) < angleaccuracy:
                    #X-axis points southwards
                    dipolelength *= 1

                elif not np.abs(angle) < angleaccuracy:
                    print 'Configuration file error. E-field X-axis angle for station %s invalid: %f'%(stationname,angle)

            #check East-West axis orientation
            if direction == 'y':
                angle = float(stationdict['e_yaxis_azimuth'])
                dipolelength = float(stationdict['e_yaxis_length'])
                if np.abs(270. - angle) < angleaccuracy:
                    #X-axis points southwards
                    dipolelength *= 1
                elif not np.abs(90. - angle) < angleaccuracy:
                    print 'Configuration file error. E-field Y-axis angle for station %s invalid: %f'%(stationname,angle)

            logger = stationdict['e_logger_type']
            gain = stationdict['e_logger_gain']
            instrument = stationdict['e_instrument_type']
            instrument_amplification = stationdict['e_instrument_amplification']

        elif field == 'b':

            dipolelength = 1.
            logger = stationdict['b_logger_type']
            gain = stationdict['b_logger_gain']
            instrument = stationdict['b_instrument_type']
            instrument_amplification = stationdict['b_instrument_amplification']


        C.calibrate_file(f, outdir, instrument, logger, gain, dipolelength, stationname, channel, latitude, longitude, elevation,  offset = 0 )

if __name__=='__main__':
    main()