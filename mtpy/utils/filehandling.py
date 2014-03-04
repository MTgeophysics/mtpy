#!/usr/bin/env python

"""
Helper functions for file handling. 

The various functions deal with renaming, sorting, 
concatenation of time series, extraction of names and times from filenames,
reading configuration files, ....


@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys
import os
import os.path as op
import calendar
import time
import fnmatch
import shutil

import mtpy.utils.calculator as MTcc
import mtpy.processing.general as MTgn
import mtpy.utils.exceptions as MTex
import mtpy.utils.format as MTft
import mtpy.utils.configfile as MTcf

reload(MTgn)
reload(MTcc)
reload(MTex)
reload(MTcf)

#=================================================================

#define uncertainty for differences between time steps
epsilon = 1e-9

#=================================================================

lo_headerelements = ['station', 'channel','samplingrate','t_min',
                    'nsamples','unit','lat','lon','elev']

#=================================================================

def make_unique_filename(infn):

    fn = op.abspath(infn)
    outfn = fn
    i = 1
    while op.isfile(outfn):
        filebase = op.splitext(fn)[0]
        outfn = filebase +'_%i'%i+ op.splitext(fn)[1]
        i += 1

    return outfn


def get_sampling_interval_fromdatafile(filename, length = 3600):
    """ 
    Find sampling interval from data file.

    Provide data file (purely numerical content) and total 
    data length in seconds (default 3600). Data are read in by 
    the Numpy 'loadtxt'-function, the lentgh of the data array yields
    the sampling interval.  

    Lines beginning with # are ignored.

    """

    fn = op.abspath(op.realpath(filename))
    dd = np.loadtxt(fn)
    sampling_interval = length/float(len(dd))

    return sampling_interval




def EDL_make_dayfiles(inputdir, sampling , stationname = None, outputdir = None):
    """

    Concatenate ascii time series to dayfiles (calendar day, UTC reference).

    Data can be within a single directory or a list of directories. 
    However, the files in the directory(ies) 'inputdir' have to be for 
    one station only, and named with a 2 character suffix, defining the channel! 

    If the time series are interrupted/discontinuous at some point, a new file 
    will be started after that point, where the file index 'idx' is increased by 1.
    If no stationname is given, the leading non-datetime characters in the first 
    filename are used.


    Files are named as 'stationname_samplingrate_date_idx.channel'
    Stationname, channel, and sampling are written to a header line.

    Output data consists of a single column float data array. The data are 
    stored into one directory. If 'outputdir' is not specified, a subdirectory 
    'dayfiles' will be created witihn the current working directory. 

    Note: 
    Midnight cannot be in the middle of a file, because only file starts are 
    checked for a new day!!

    """
    try:
        if type(inputdir)==str:
            raise
        lo_foldernames = [i for i in inputdir]
    except TypeError:
        lo_foldernames = [inputdir]

    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    lo_allfiles = []
    pattern = '*.[ebEB][xyzXYZ]'
    if stationname is not None:
        pattern = '*{0}*.[ebEB][xyzXYZ]'.format(stationname)
    print '\nSearching for files with pattern: ',pattern

    for folder in lo_foldernames:
        wd = op.abspath(op.realpath(folder)) 
        if not op.isdir(wd):
            #print 'Directory not existing: %s' % (wd)
            lo_foldernames.remove(wd)
            continue    

        lo_dirfiles = [op.abspath(op.join(wd,i))  for i in os.listdir(wd) 
                        if fnmatch.fnmatch(i,pattern) is True]
        lo_allfiles.extend(lo_dirfiles)   

    #check, if list of files is empty
    if len(lo_allfiles) == 0:
        if stationname is not None:
            raise MTex.MTpyError_inputarguments('Directory(ies) do(es) not contain'\
            ' files to combine for station {0}:\n {1}'.format(stationname, inputdir))

        raise MTex.MTpyError_inputarguments('Directory does not contain files'\
                                            ' to combine:\n {0}'.format(inputdir))

    #define subfolder for storing dayfiles
    outpath = op.join(os.curdir,'dayfiles')    
    if outputdir is not None:
        try:
            outpath = op.abspath(op.join(os.curdir,outputdir))
            if not op.exists(outpath):
                try:
                    os.makedirs(outpath)
                except:
                    raise
            if not os.access(outpath, os.W_OK):
                raise
        except:
            print 'Cannot generate writable output directory {0} - using'\
                    ' generic location "dayfiles" instead'.format(outpath)
            outpath = op.join(wd,'dayfiles')    
            pass

    #generate subfolder, if not existing
    if not op.exists(outpath):
        try:
            os.makedirs(outpath)
        except:
            MTex.MTpyError_inputarguments('Cannot generate output'\
                                ' directory {0} '.format(outpath))

    #outer loop over all components
    for comp in components:

        #make list of files for the current component
        lo_files = np.array([op.join(wd,i) for i in lo_allfiles 
                            if (i.lower()[-2:] == comp)])

        #make list of starting times for the respective files
        lo_starttimes = np.array([EDL_get_starttime_fromfilename(f) 
                                    for f in lo_files])
        
        #sort the files by their starting times
        idx_chronologic = np.argsort(lo_starttimes)
        
        #obtain sorted lists of files and starting times
        lo_sorted_files = list(lo_files[idx_chronologic])
        lo_sorted_starttimes = list(lo_starttimes[idx_chronologic])
        idx = 0
        while idx < len(lo_sorted_starttimes):
            try:
                val = lo_sorted_starttimes[idx]
            except:
                break
            if val is None:
                dummy = lo_sorted_files.pop(idx)
                dummy = lo_sorted_starttimes.pop(idx)
            else:
                idx += 1


        #set stationname, either from arguments or from filename
        if stationname is None:
            stationname = EDL_get_stationname_fromfilename(lo_sorted_files[0]).upper()

        #set counting variables - needed for handling of consecutive files

        sameday = 0
        fileopen = 0
        incomplete = 0
        fileindex = 0

        #allocate a data array to fill
        # this is more memory efficient than extending lists!!
        #cater for potential rounding errors:
        if sampling < 1:
            max_n_data = 86400 * (int(1./sampling)+1)
        else:
            max_n_data = int(86400./sampling) + 1

        day_data = np.zeros(max_n_data,'float32')


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

            tmp_file_time_axis = np.arange(no_samples)*sampling+file_start_time
            #file_time_axis = (np.arange(no_samples)*sampling +
            #                 file_start_time).tolist()


            #time of the last sample + 1x sampling-interval
            #file_end_time =  file_time_axis[-1] + sampling
            file_end_time =  tmp_file_time_axis[-1] + sampling
         



            #set the time as starting time for output file, if no output file is open already
            if fileopen == 0:
                outfile_starttime =  file_start_time
                #outfile_timeaxis = file_time_axis
                old_time_axis = tmp_file_time_axis[:]
                
                arrayindex = 0

                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in                    
                    #outfile_data = data_in.tolist()
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in[:,1]
                    #outfile_data = data_in[:,1].tolist()
                
                arrayindex += len(data_in)


                file_date = '{0}{1:02}{2:02}'.format(file_start[0],
                                                 file_start[1], file_start[2]) 


                #define output filename
                new_fn = '{0}_1day_{1}_{2}.{3}'.format(stationname,
                                                 file_date, fileindex, comp)
                
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
                #elif (outfile_timeaxis[-1] - file_start_time) > epsilon:
                elif (old_time_axis[-1] - file_start_time) > epsilon:

                    #find point on the outfile time axis for the beginning of current file:
                    overlap_idx = np.argmin(np.abs(np.array(
                                            outfile_timeaxis) - file_start_time)) 

                    #set the array index back
                    arrayindex = overlap_idx
                   
                    #re-define outfile time axis and data
                    # outfile_timeaxis = np.delete(outfile_timeaxis,
                    #                              np.arange(len(outfile_timeaxis) - 
                    #                              overlap_idx) + 
                    #                              overlap_idx).tolist()

                    # outfile_data = np.delete(outfile_data, 
                    #                             np.arange(len(outfile_data) - 
                    #                             overlap_idx) + 
                    #                             overlap_idx).tolist()
                

                old_time_axis = tmp_file_time_axis[:]
                #append current file's time axis
                #outfile_timeaxis.extend(file_time_axis)
                    
                #append current data                  
                #if it's a single column of data
                if np.size(data_in.shape) == 1:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in                    
                    #outfile_data.extend(data_in.tolist())
                #otherwise assuming that the first column is time, so just take the second one
                else:
                    day_data[arrayindex:arrayindex+len(data_in)] = data_in[:,1]                    
                    #outfile_data.extend(data_in[:,1].tolist())


                arrayindex += len(data_in)


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
                headerline = '# {0} {1} {2:.1f} {3} {4} \n'.format(
                                    stationname, comp.lower(), 1./sampling, 
                                    outfile_timeaxis[0], len(outfile_timeaxis))

                F.write(headerline)

                #outfile_array = np.zeros((len(outfile_timeaxis),2))
                #outfile_array[:,0] = outfile_timeaxis
                #outfile_array[:,1] = outfile_data

                np.savetxt(F,day_data[:arrayindex])
                #np.savetxt(F, np.array(outfile_data))
                arrayindex = 0
                
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
    
    try:
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

    except:
        epochtime = None

    return epochtime


def EDL_get_stationname_fromfilename(filename):

    bn = op.basename(filename)
    parts_of_bn = bn.split('.')
    stationtime = parts_of_bn[-2]

    stationname = stationtime[:-12].upper()

    if len(stationname) == 0:
        stationname = 'DUMMYSTATION'


    return stationname


def read_data_header(fn_raw):
    """
    Deprecated!!!
    USE 
              read_ts_header

    INSTEAD


        Read the header line of MTpy TS data files.

    input
    -----
    MTpy TS data file name

    output
    -------
    list of header elements:
    stationname, channel, sampling rate, starttime first sample, 
    starttime last sample, unit, lat, lon, elevation

    """


    fn = op.abspath(op.realpath(fn_raw))

    if not op.isfile(fn):
        raise MTex.MTpyError_inputarguments('Not a file:%s'%fn)
    try:
        F = open(fn, 'r')
    except:
        raise MTex.MTpyError_inputarguments('File not readable:%s'%fn)

    firstline = F.readline().strip().split()
    if not firstline[0][0] == '#':
        raise MTex.MTpyError_ts_data('Time series data file does '
            'not have a proper header:%s'%fn)

    F.close()

    header_list = []

    idx_header = 0

    if len(firstline[0]) > 1:
        header_list.append(firstline[0][1:].upper())
    else:
        header_list.append(firstline[1].upper())
        idx_header += 1

    header_list.append( firstline[idx_header+1].lower() )
    header_list.append( float(firstline[idx_header+2]) )
    header_list.append( float(firstline[idx_header+3]) )
    header_list.append( int(float(firstline[idx_header+4])) )
    header_list.append( firstline[idx_header+5].lower() )
    header_list.append( float(firstline[idx_header+6]) )
    header_list.append( float(firstline[idx_header+7]) )
    header_list.append( float(firstline[idx_header+8]) )


    return header_list


def read_2c2_file(filename):
    """
    Read in BIRRP 2c2 coherence files and return 4 lists 
    containing [period],[freq],[coh],[zcoh]. Note if any of the coherences are 
    negative a value of 0 will be given to them.

    """

    period = []
    freq = []
    coh1 = []
    zcoh1 = []

    F_in = open(filename,'r')
    data_raw = F_in.readlines()
    
    for ii in range(len(data_raw)):

        coh_row = data_raw[ii].strip().split()
        
        try:
            period.append(float(coh_row[0]))
        except:
            period.append(0.)
        try:
            freq.append(  float(coh_row[1]))
        except:
            period.append(0.)
        try:
            coh1.append(  float(coh_row[2]))
        except:
            period.append(0.)
        try:
            zcoh1.append( float(coh_row[3]))
        except:
            period.append(0.)

    return period, freq, coh1, zcoh1

def validate_ts_file(tsfile):
    """
        Validate MTpy timeseries (TS) data file
        Return Boolean value True/False .

    """ 
    tsfile = op.abspath(tsfile)

    try:
        header = read_ts_header(tsfile)

        if header['station'] is None:
            #print 'header'
            raise
        if header['channel'] is None:
            #print 'channel'
            raise
        
        sr = float(header['samplingrate'])
        t0 = float(header['t_min'])
        ns = int(float(header['nsamples']))
        
        data = np.loadtxt(tsfile)
        
        if len(data) != ns:
            #print 'data length'
            raise
        if data.dtype not in [int, float]:
            #print 'data type'
            raise

    except:
        #print 'number'
        return False


    return True



def read_ts_header(tsfile):
    """ Read in the header line from MTpy timeseries data files.
        
        Return header as dictionary. Return empty dict, 
        if no header line was found.
    """

    header_dict = {}

    tsfile = op.abspath(tsfile)
    
    if not op.isfile(tsfile):
        raise MTex.MTpyError_inputarguments('Error - '
            'input file not existing: {0}'.format(tsfile))

    try:
        with open(tsfile,'r') as F:
            firstline =''
            #ignoring empty lines or lines with just the '#' character in it
            while len(firstline) == 0:
                firstline = F.readline().strip()
                if firstline == '#':
                    firstline = ''
        if firstline[0] != '#':
            raise
    except:
        raise MTex.MTpyError_ts_data('No header line found -'
            ' check file: {0}'.format(tsfile))
        

    firstline = firstline.replace('#','')
    headerlist = firstline.split()


    for i in range(len(headerlist)):
        header_dict[lo_headerelements[i]] = headerlist[i]
        #old header had tmax instead of n_samples:
        if ((i == 4) and float(headerlist[4])%1 != 0 
            and float(headerlist[i]) > float(headerlist[i-1])):
            header_dict[lo_headerelements[i]] = int(
                    (float(headerlist[i]) - float(headerlist[i-1])
                    )*float(headerlist[i-2]) )+1

    headerlements = ['samplingrate','t_min','nsamples','lat','lon','elev']

    for h in headerlements:                   
        try: 
            header_dict[h] = float(header_dict[h])
        except:
            pass

    return header_dict


def get_ts_header_string(header_dictionary):
    """
        Return a MTpy time series data file header string from a dictionary.

    """
    
    header_string = '# '
    for headerelement in lo_headerelements:
        if header_dictionary.has_key(headerelement):
            header_string += '{0} '.format(str(header_dictionary[headerelement]))
        else:
            header_string += '\t '   

    header_string += '\n'

    return header_string



def write_ts_file_from_tuple(outfile,ts_tuple, fmt='%.8e'):
    """
        Write an MTpy TS data file, where the content is provided as tuple:

        (station, channel,samplingrate,t_min,nsamples,unit,lat,lon,elev, data)

        todo:
        needs tuple-validation

    """

    
    header_dict = {}
    for i in range(len(ts_tuple) -1):
        if ts_tuple[i] is not None:
            header_dict[lo_headerelements[i]] = ts_tuple[i]

    header_string = get_ts_header_string(header_dict)
    data = ts_tuple[-1]

    outfilename = make_unique_filename(outfile)


    try:
        outF = open(outfilename,'w')
        outF.write(header_string)
        np.savetxt(outF, data, fmt=fmt)
        outF.close()
    except:
        raise MTex.MTpyError_inputarguments('ERROR - could not write content'
                            ' of TS tuple to file : {0}'.format(outfilename))

    return outfilename


def read_ts_file(mtdatafile):
    """
        Read an MTpy TS data file and provide the content as tuple:

        (station, channel,samplingrate,t_min,nsamples,unit,lat,lon,elev, data)
        If header information is incomplete, the tuple is filled up with 'None'

    """

    infile = op.abspath(mtdatafile)
    if not op.isfile(infile):
        raise MTex.MTpyError_inputarguments('ERROR - Data file not '
                                                'existing: {0}'.format(infile))

    header = read_ts_header(infile)
    if len(header) == 0 :
        raise MTex.MTpyError_inputarguments('ERROR - Data file not valid - '
                                        'header is missing : {0}'.format(infile))

    data = np.loadtxt(infile)
    if len(data) != int(float(header['nsamples'])):
        raise MTex.MTpyError_inputarguments('ERROR - Data file not valid '
                                    '- wrong number of samples in data ({1} '
                                    'instead of {2}): {0}'.format(
                                    infile,len(data) , int(float(
                                        header['nsamples']))) )

    lo_header_contents = []

    for i in lo_headerelements:
        if i in header:
            lo_header_contents.append(header[i])
        else:
            lo_header_contents.append(None)
 
    lo_header_contents.append(data)

    return tuple(lo_header_contents)


def reorient_files(lo_files, configfile, lo_stations = None, outdir = None):

    #read config file
    try:
        config_dict = MTcf.read_survey_configfile(configfile)
    except:
        raise MTex.MTpyError_config_file( 'Config file cannot be read:'
                                                    ' {0}'.format(configfile) )

    if lo_stations is not None:
        try:
            if type(lo_stations) == str:
                raise
            #check, if it's iterable:
            dummy = [i for i in lo_stations]
        except:
            raise MTex.MTpyError_inputarguments('ERROR - "lo_stations"'
                                                ' argument must be iterable!')
    
    #Do not require list of headers as input, as this function can be called directly rather than from a 'calibratefiles.py'-like script - so the list not necessarily exists in beforehand - 
    #collect header lines of files in list
    lo_headers = []
    lo_stationnames = []
    for file_idx, filename in enumerate(lo_files):
        header = read_ts_header(filename)
        station = header['station']
        if station.upper() not in [i.upper() for i in lo_stations]:
            #TODO: check, if this causes problems with the indices for the current loop:
            lo_files.remove(filename)
            continue
        lo_headers.append(header)
        lo_stationnames.append(station.upper())



    if len(lo_headers) == 0 :
        if lo_stations is not None:
            print 'ERROR - No files with header lines found for station(s)'\
                                                    ' {0}'.format(lo_stations)
        else:
            print 'ERROR - No files with header lines found'
        return 1

    lo_stationnames = list(set(lo_stationnames))

    # set up output directory 
    ori_outdir = op.abspath(op.join(os.curdir,'reoriented'))

    if outdir is not None:
        try:
            ori_outdir = op.abspath(op.join(os.curdir,outdir))
            if not op.isdir(ori_outdir):
                os.makedirs(ori_outdir)
        except:
            print 'Output directory cannot be generated: {0} - using generic'\
                                                ' location'.format(ori_outdir)
            ori_outdir = op.abspath(op.join(os.curdir,'reoriented'))
    try:
        if not op.isdir(ori_outdir):
            os.makedirs(ori_outdir)
    except:
        #this only comes up, if the generic location cannot be generated
        raise MTex.MTpyError_inputarguments('Generic directory cannot be'
                                        ' generated: {0}'.format(ori_outdir))

    #----------------------
    #start re-orientation
    #present: list of all files, list of all headers, list of all stations
    
    for sta_idx, sta in enumerate(lo_stationnames):
        #print sta
        try:
            stationconfig = config_dict[sta]
        except:
            print 'Warning - No config file entry for station {0} -'\
                                        ' no processing possible'.format(sta)
            continue
        
        declination = float(stationconfig.get('declination',0.))


        for sensor in ['e','b']:
            #TODO:
            # reduce this function to the re-orientation of files that have the same length for X and Y. Do the puzzlling for varying lengths later!!

            for idx_h_x, header_x in enumerate(lo_headers):
                #looking for one specific station
                if not header_x['station'].upper() == sta.upper():
                    continue
                #looking for the specific sensor type
                if not header_x['channel'].lower()[0] == sensor:
                    continue
                #looking for the X channel (the to-be-North)
                if not header_x['channel'].lower()[1] == 'x':
                    continue

                x_file = lo_files[idx_h_x]
                x_header_string = get_ts_header_string(header_x)

                t0 = float(header_x['t_min'])
                #print t0 
                #now look for the respective y-file and possible z-file - unfortunately by another loop over all headers:
                y_file = None
                z_file = None
                for idx_h_y, header_y in enumerate(lo_headers):
                    if (header_y['station'].upper() == sta.upper()) and \
                        (header_y['channel'].lower()[0] == sensor) and \
                        (float(header_y['t_min']) == float(header_x['t_min'] ) ):
                        if (header_y['channel'].lower()[1] == 'y') :
                            y_file = lo_files[idx_h_y]
                            y_header_string = get_ts_header_string(header_y)

                        elif   (header_y['channel'].lower()[1] == 'z') :
                            z_file = lo_files[idx_h_y]

                    else:
                        continue
                if y_file == None:
                    continue

                x_outfn = op.abspath(op.join(ori_outdir,op.basename(x_file)))
                y_outfn = op.abspath(op.join(ori_outdir,op.basename(y_file)))
                if z_file is not None:
                    z_outfn = op.abspath(op.join(ori_outdir,op.basename(z_file)))
                

                xdata = np.loadtxt(x_file)
                ydata = np.loadtxt(y_file)

                #declination is positive, if magnetic North is east of true North.
                # the measured angles are w.r.t. magnetic North, so the given 
                # azimuths do not include the declination 
                #-> thus the declination value is added to azimuths
                if sensor == 'e':
                    xangle = float(stationconfig.get(
                                    'e_xaxis_azimuth', 0.)) + declination
                    yangle = float(stationconfig.get(
                                    'e_yaxis_azimuth',90.)) + declination
                else:
                    xangle = float(stationconfig.get(
                                    'b_xaxis_azimuth', 0.)) + declination
                    yangle = float(stationconfig.get(
                                    'b_yaxis_azimuth',90.)) + declination                


                newx, newy =  MTcc.reorient_data2D(xdata, ydata, 
                            x_sensor_angle = xangle , y_sensor_angle = yangle)
                #print xdata.shape, ydata.shape, newx.shape, newy.shape 

                #continue
                outFx = open(x_outfn,'w')
                outFx.write(x_header_string)
                np.savetxt(outFx,newx)
                outFx.close()
                outFy = open(y_outfn,'w')
                outFy.write(y_header_string)
                np.savetxt(outFy,newy)
                outFy.close()
                written_files = [x_outfn,y_outfn]
                if z_file is not None:
                    shutil.copyfile(z_file, z_outfn)
                    written_files.append(z_outfn)
                print 'written files {0}'.format(written_files)


            
    return 0
