#!/usr/bin/env python

"""
This modules contains helper functions for file handling. 

The various functions deal with renaming, sorting, concatenation of time series, extraction of names and times from filenames, ....


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import re
import sys, os
import glob
import os.path as op
import glob
import calendar
import time


from mtpy.utils.exceptions import *
#=================================================================

#define uncertainty for differences between time steps
epsilon = 1e-9




#=================================================================



def get_sampling_interval_fromdatafile(filename, length = 3600):
    """ 
    Find sampling interval from data file.

    Provide data file (purely numerical content) and total data length in seconds (default 3600). 
    Data are read in by loadtxt-function, the lentgh of the data array yields the sampling interval. 
    Lines beginning with # are ignored.

    """

    fn = op.abspath(op.realpath(filename))
    dd = np.loadtxt(fn)
    sampling_interval = length/float(len(dd))

    return sampling_interval




def EDL_make_dayfiles(foldername, sampling , stationname = None):
    """

    Concatenate ascii time series to dayfiles (calendar day, UTC reference).

    Files in the directory have to be for one station only !! 

    If the time series are interrupted, a new file will be started after that point,
    where the index idx is increased by 1.
    If no stationname is given, the leading non-datetime characters in the first filename are used.


    Files are named as 'stationname_samplingrate_date_idx.channel'
    Stationname, sampling, and timestasmps for the first and the last 
    sample are written to a header line.

    Attention: Midnight must be at the start of a file, because only file starts are checked for a new day!!

    """


    wd = op.abspath(op.realpath(foldername))
    

    if not op.isdir(wd):
        raise MTpyError_inputarguments('Directory not existing: %s' % (wd))

    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz', 'tp']

    oldwd = os.getcwd()
    os.chdir(wd)
    lo_allfiles = glob.glob('*.??')
    os.chdir(oldwd)


    #check, if list of files is empty
    if len(lo_allfiles) == 0:
        raise MTpyError_inputarguments('Directory does not contain files to combine: %s' % (wd))

    #define subfolder for storing dayfiles
    outpath = op.join(wd,'dayfiles')

    #generate subfolder, if not existing
    if not op.exists(outpath):
        os.makedirs(outpath)

    #outer loop over all components
    for comp in components:

        #make list of files for thte current component
        lo_files = np.array([op.join(wd,i) for i in lo_allfiles if (i.lower()[-2:] == comp)])

        #make list of starting times for the respective files
        lo_starttimes = np.array([EDL_get_starttime_fromfilename(f) for f in lo_files])
        
        #sort the files by their starting times
        idx_chronologic = np.argsort(lo_starttimes)
        
        #obtain sorted lists of files and starting times
        lo_sorted_files = lo_files[idx_chronologic]
        lo_sorted_starttimes = lo_starttimes[idx_chronologic]

        #set stationname, either from arguments or from filename
        if not stationname:
            stationname = EDL_get_stationname_fromfilename(lo_sorted_files[0])


        #set counting variables - needed for handling of consecutive files

        sameday = 0
        fileopen = 0
        incomplete = 0
        fileindex = 0

        #loop over all (sorted) files for the current component
        for idx_f,f in enumerate(lo_sorted_files):

            print 'read in file %s' %(f)
            #starting time of current file
            file_start_time = lo_sorted_starttimes[idx_f]

            #get tuple with the starting time of the current file
            file_start = time.gmtime(file_start_time)
            
            #read in raw data
            data_in = np.loadtxt(f)

            no_samples = len(data_in)

            file_time_axis = (np.arange(no_samples)*sampling + file_start_time).tolist()


            #time of the last sample + 1x sampling-interval
            file_end_time =  file_time_axis[-1] + sampling
         



            #set the time as starting time for output file, if no output file is open already
            if fileopen == 0:
                outfile_starttime =  file_start_time
                outfile_timeaxis = file_time_axis

                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    outfile_data = data_in.tolist()
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    outfile_data = data_in[:,1].tolist()


                file_date = '%i%02i%02i'%(file_start[0], file_start[1], file_start[2]) 


                #define output filename
                new_fn = '%s_%.1fHz_%s_%i.%s'%(stationname, 1./sampling, file_date, fileindex, comp)
                new_file = op.abspath(op.join(outpath,new_fn))
                
                #open output file 
                F = open(new_file,'w')
                
                fileopen = 1


            
            else:
                #check, if the new file ends earlier than data in buffer.
                #if yes, just skip this file:
                if file_end_time < outfile_timeaxis[-1]:
                    continue 

                #if current file starts earlier than the endtime of data in buffer then delete ambiguous  parts of the buffer:
                elif (outfile_timeaxis[-1] - file_start_time) > epsilon:

                    #find point on the outfile time axis for the beginning of current file:
                    overlap_idx = np.argmin(np.abs(np.array(outfile_timeaxis) - file_start_time)) 

                    #re-define outfile time axis and data
                    outfile_timeaxis = np.delete(outfile_timeaxis, np.arange(len(outfile_timeaxis)- overlap_idx) + overlap_idx).tolist()

                    outfile_data = np.delete(outfile_data, np.arange(len(outfile_data) - overlap_idx) + overlap_idx).tolist()
                

                #append current file's time axis
                outfile_timeaxis.extend(file_time_axis)
                    
                #append current data                  
                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    outfile_data.extend(data_in.tolist())
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    outfile_data.extend(data_in[:,1].tolist())




            #-----------

            #check, if there is a next file:
            try:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
            except:
                incomplete = 1

            #if there is a next file, 
            # - check, if it's the same day
            # - check, if it continues at the end of the current one:
            if incomplete == 0:
                next_file_start_time = lo_sorted_starttimes[idx_f + 1]
                next_file_start = time.gmtime(next_file_start_time)
                
                if next_file_start[2] == file_start[2] :
                    #print 'sameday',file_start[:]
                    sameday = 1
                else:
                    incomplete = 1
                    sameday = 0
                    fileindex = 0
                    #print '\t NOT sameday', fileindex


                if next_file_start_time - file_end_time > epsilon: 
                    incomplete = 1

            if incomplete == 1 and sameday == 1 : 
                fileindex +=1        

           

            #check, if the file has to be closed and written now
            if incomplete == 1 :

                #define header info
                headerline = '# %s %i Hz first sample: %i - last sample: %i (epochs)\n'%(stationname, int(1./sampling), outfile_timeaxis[0],outfile_timeaxis[-1] )

                F.write(headerline)

                np.savetxt(F, np.array(outfile_data))

                F.close()
                print '\t wrote file %s'%(new_file)

                fileopen = 0
                incomplete = 0
    




def EDL_get_starttime_fromfilename(filename): 
    """ 
    Return starttime of data file in epoch seconds.

    Starting time is determined by the filename. This has to be of the form
    'somthing/*.stationname.ddmmyyHHMMSS.??'

    """     
    #clip parent paths and structure
    bn = op.basename(filename)
    parts_of_bn = bn.split('.')
    timestamp = parts_of_bn[-2]
    
    secs = int(float(timestamp[-2:]))
    mins = int(float(timestamp[-4:-2]))
    hours =int(float(timestamp[-6:-4]))
    day =  int(float(timestamp[-8:-6]))
    month =int( float(timestamp[-10:-8]))
    year = int(float(timestamp[-12:-10]))
    if year < 50:
        year += 2000
    else:
        year += 1900

    timetuple = (year, month, day, hours, mins, secs)

    epochtime = calendar.timegm(timetuple)

    return epochtime


def EDL_get_stationname_fromfilename(fn):

    bn = op.basename(fn)
    parts_of_bn = bn.split('.')
    stationtime = parts_of_bn[-2]

    stationname = stationtime[:-12].upper()

    if len(stationname) == 0:
        stationname = 'DUMMYSTATION'


    return stationname

