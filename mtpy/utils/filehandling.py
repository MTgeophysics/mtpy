#!/usr/bin/env python

"""
This module contains helper functions for file handling. 

The various functions deal with renaming, sorting, 
concatenation of time series, extraction of names and times from filenames,
reading configuration files, ....


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
import ConfigParser

from mtpy.utils.exceptions import *
import mtpy.utils.format as MTformat
#=================================================================

#define uncertainty for differences between time steps
epsilon = 1e-9


#=================================================================

def read_configfile(filename):
    """
        Read a general config file and return the content as dictionary.

        Config files without sections or only DEFAULT section -> return dictionary
        Config files with sections -> return nested dictionary (main level keys are section heads)
        Config files with sections as well as section-less entries -> return nested dictionary, which includes a 'DEFAULT' key
    """

    #generate config parser instance
    configobject = ConfigParser.ConfigParser()
    
    #check, if file is present
    if not op.isfile(filename):
        raise MTpyError_inputarguments( 'File does not exist: %s'%filename )

    # try to parse file - exit, if not a config file
    try:
        configobject.read(filename)
    except:
        raise MTpyError_inputarguments( 'File is not a proper configuration file: %s'%filename )


    #if 0:#len(configobject.sections()) != 0:
    config_dict = configobject._sections
        
    if len (config_dict.keys()) != 0:

        config_dict['DEFAULT'] = configobject.defaults()
    else:
        config_dict = configobject.defaults()


    return config_dict


def read_survey_configfile(filename):
    """
    Read in a survey configuration file and return a dictionary.

    Input config file must contain station names as section headers!

    The output dictionary keys are station names (capitalised), the values are (sub-)dictionaries.
    The configuration file must contain sections for all stations, each containing all mandatory keywords:
    
    - latitude (deg)
    - longitude (deg)
    - elevation (in meters)
    - sampling_interval (in seconds)
    - station_type (MT, E, B)

    Depending on the type of station the following entries are required.

    E-field recorded:

    - E_logger_type (edl/elogger)
    - E_logger_gain (factor/gain level)
    - E_instrument_type (electrodes)
    - E_instrument_amplification (applied amplification factor)
    - E_Xaxis_azimuth (degrees)
    - E_Xaxis_length (in meters)
    - E_Yaxis_azimuth (degrees)
    - E_Yaxis_length (in meters)

    B-field recorded:
  
    - B_logger_type (edl)
    - B_logger_gain (factor/gain level)
    - B_instrument_type (coil, fluxgate)
    - B_instrument_amplification (applied amplification factor)

    """

    list_of_required_keywords = ['latitude',
                                'longitude',
                                'elevation',
                                'sampling_interval',
                                'station_type'
                                ]

    list_of_efield_keywords = [ 'E_logger_type',
                                'E_logger_gain',
                                'E_instrument_type',
                                'E_instrument_amplification',
                                'E_Xaxis_azimuth',
                                'E_Xaxis_length',
                                'E_Yaxis_azimuth',
                                'E_Yaxis_length'
                                ]

    list_of_bfield_keywords = [ 'B_logger_type',
                                'B_logger_gain',
                                'B_instrument_type',
                                'B_instrument_amplification'
                              ]


    dict_of_allowed_values_efield = {'E_logger_type':['edl','elogger'] ,
                                    'E_logger_gain': ['low', 'verylow','high', 0.4, 1, 10, 11],
                                    'E_instrument_type':['electrodes'],
                                    'E_instrument_amplification':[1,10]
                                    }
    
    dict_of_allowed_values_bfield = {'B_logger_type':['edl'] ,
                                    'B_logger_gain': ['low', 'verylow','high', 0.4, 1, 10],
                                    'B_instrument_type':['fluxgate', 'coil'],
                                    'B_instrument_amplification':[1]
                                    }



    list_of_station_types = ['mt','e','b']


    error_counter = 0

    #generate config parser instance
    configobject = ConfigParser.ConfigParser()

    #check, if file is present
    if not op.isfile(filename):
        raise MTpyError_inputarguments( 'File does not exist: %s'%filename )


    # try to parse file - exit, if not a config file
    try:
        configobject.read(filename)
    except:
        raise MTpyError_inputarguments( 'File is not a proper configuration file: %s'%filename )

    #obtain dict of dicts containing the input file's sections (station names)
    #excludes DEFAULT section and key-value pairs without section header
    configobject_dict = configobject._sections

    #initialise the output dictionary
    config_dict = {}

    #loop over the sections (stations) of the config file
    for station in configobject_dict:
        #read in the sub-dictionary for the current station - bringing all keys to lowercase!
        temp_dict_in = dict((k.lower(),v) for k,v in configobject_dict[station].items())

        #initialise output sub-directory for current station 
        stationdict = {}

        #stationnames are uppercase in MTpy
        stationname = station.upper()


        #check for presence of all mandatory keywords for the current station
        #case insensitive - allow for short forms 'lat', 'lon', and 'ele'
        for req_keyword in list_of_required_keywords:
            if req_keyword.lower() in temp_dict_in.keys():
                stationdict[req_keyword.lower()] = temp_dict_in[req_keyword.lower()].lower()
            elif req_keyword in ['latitude', 'longitude', 'elevation']:
                if req_keyword[:3] in temp_dict_in.keys():
                    stationdict[req_keyword] = temp_dict_in[req_keyword[:3]]
            else:  
                print 'Station %s - keyword %s missing'%(stationname, req_keyword)
                error_counter += 1
                continue

        #check format of lat/lon - convert to degrees, if given in (deg,min,sec)-triple
        for coordinate in ['latitude', 'longitude', 'elevation']:
            value = stationdict[coordinate]
            try:
                new_value = MTformat._assert_position_format(coordinate,value)
            except:
                raise MTpyError_config_file('Error - wrong coordinate format for station %s'%(stationname))
            stationdict[coordinate] = new_value

        if not stationdict['station_type'] in list_of_station_types:
            raise MTpyError_config_file( 'Station type not valid' )


        if stationdict['station_type'] in ['mt','e']:
            #check for required electric field parameters
            for req_keyword in list_of_efield_keywords:
                if req_keyword.lower() in temp_dict_in.keys():
                    stationdict[req_keyword.lower()] = temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print 'Station %s - keyword %s missing'%(stationname, req_keyword)
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_efield)
            

        if stationdict['station_type'] in ['mt','b']:
            #check for required magnetic field parameters
            for req_keyword in list_of_bfield_keywords:
                if req_keyword.lower() in temp_dict_in.keys():
                    stationdict[req_keyword.lower()] = temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print 'Station %s - keyword %s missing'%(stationname, req_keyword)
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_bfield)
            

        #add the station's sub-dictionary to the config dictionary
        config_dict[stationname] = stationdict

    #re-loop for setting up correct remote reference station information :
    #if rem.ref. station key is present, its information must be contained in the config file!
    for station in config_dict.iterkeys():
        stationdict = config_dict[station]

        stationdict['rr_station'] = None
        stationdict['rr_station_latitude'] = None
        stationdict['rr_station_longitude'] = None
        stationdict['rr_station_elevation'] = None


        if stationdict.has_key('rr_station'):
            rem_station = stationdict['rr_station'] 
            try:
                #check, if values are contained in dict 
                float(stationdict['rr_station_latitude'] )
                float(stationdict['rr_station_longitude'])
                float(stationdict['rr_station_elevation'])
            except:
                try:
                    #check for shortened form
                    stationdict['rr_station_latitude']  = float(stationdict['rr_station_lat'] )
                    stationdict['rr_station_longitude'] = float(stationdict['rr_station_lon'] )
                    stationdict['rr_station_elevation'] = float(stationdict['rr_station_ele'] )                 

                except:
                    try:
                        #read from other config dict entry
                        stationdict['rr_station_latitude'] = config_dict[rem_station]['latitude']
                        stationdict['rr_station_longitude'] = config_dict[rem_station]['longitude']
                        stationdict['rr_station_elevation'] = config_dict[rem_station]['elevation']

                    except:
                        #if finally failed to read rr_station info, set rr_station back to None
                        stationdict['rr_station'] = None
                        stationdict['rr_station_latitude'] = None
                        stationdict['rr_station_longitude'] = None
                        stationdict['rr_station_elevation'] = None

        #check consistency of coordinates, if rr_station is present
        if stationdict['rr_station'] != None:
            try:
                stationdict['rr_station_latitude'] = MTformat._assert_position_format('latitude',stationdict['rr_station_latitude'])
                stationdict['rr_station_longitude'] = MTformat._assert_position_format('longitude',stationdict['rr_station_longitude'])
                stationdict['rr_station_elevation'] = MTformat._assert_position_format('elevation',stationdict['rr_station_elevation'])

            except:
                print 'Problem with remote reference station (%s) - remote reference (%s) coordinates invalid - remote reference set to None'%(station, stationdict['rr_station'] )
                stationdict['rr_station'] = None
                stationdict['rr_station_latitude'] = None
                stationdict['rr_station_longitude'] = None
                stationdict['rr_station_elevation'] = None


    if error_counter != 0:
        print 'Could not read all mandatory sections and options in config file - found %i errors - check configuration file before continuing!' %error_counter
    
    else:
        return config_dict

#=================================================================

def write_dict_to_configfile(dictionary, output_filename):
    """
    Write a dictionary into a configuration file.

    The dictionary can contain pure key-value pairs as well as a level-1 nested dictionary. In the first case, the entries are stored in a 'DEFAULT' section. In the latter case, the dictionary keys are taken as section heads and the sub-dictionaries key-value pairs fill up the respective section  

    """

    configobject = ConfigParser.ConfigParser()

    #check for nested dictionary - 
    #if the dict entry is a key-value pair, it's stored in a section with head 'DEFAULT' 
    #otherwise, the dict key is taken as section header
    for key,val in dictionary.items():
        try:
            for subkey, subval in val.items():
                sectionhead = key
                if not configobject.has_section(sectionhead):
                    configobject.add_section(sectionhead)
                configobject.set(sectionhead,subkey, subval)

        except:
            #if not configobject.has_section('DEFAULT'):
            #    configobject.add_section('')
            configobject.set('',key,val)


    with open(output_filename, 'w') as F:
        configobject.write(F)
    

#=================================================================
 
def _validate_dictionary(dict2validate,referencedict):

    for key, value in referencedict.items():
        #make everything to strings - easier to compare
        #in case of numbers, make to float first
        try:
            value2validate = str(float(dict2validate[key.lower()]))
        except:
            value2validate = str( dict2validate[key.lower()] ).lower()
        
        tmp = []
        for i in value:
            try: 
                tmp.append(str(float(i)))
            except:
                tmp.append(str(i))  
        value = tmp
        if not value2validate in value:
            raise MTpyError_config_file( 'Config file error -- key %s, value %s not valid'%(key, value2validate) )


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
    Stationname, channel,  and sampling are written to a header line.

    Output data consists of two column array: (time, data); timestamp given in epoch seconds.

    Attention: 
    Midnight cannot be in the middle of a file, because only file starts are checked for a new day!!

    """


    wd = op.abspath(op.realpath(foldername))
    

    if not op.isdir(wd):
        raise MTpyError_inputarguments('Directory not existing: %s' % (wd))

    #typical suffixes for EDL output file names
    components = ['ex', 'ey', 'bx', 'by', 'bz']

    oldwd = os.getcwd()
    os.chdir(wd)
    lo_allfiles = glob.glob('*.??')
    lo_allfiles = [op.abspath(i) for i in lo_allfiles]
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
                new_fn = '%s_1day_%s_%i.%s'%(stationname, file_date, fileindex, comp)
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
                headerline = '# %s %s %.1f %f %i \n'%(stationname, comp.lower(), 1./sampling, outfile_timeaxis[0], len(outfile_timeaxis) )

                F.write(headerline)

                #outfile_array = np.zeros((len(outfile_timeaxis),2))
                #outfile_array[:,0] = outfile_timeaxis
                #outfile_array[:,1] = outfile_data

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
        Read the header line of MTpy data files.

    input: 
    MTpy data file name

    output:
    list of header elements -
    stationname, channel, sampling rate, starttime first sample, starttime last sample, unit, lat, lon, elevation

    """


    fn = op.abspath(op.realpath(fn_raw))

    if not op.isfile(fn):
        raise MTpyError_inputarguments('Not a file:%s'%fn)
    try:
        F = open(fn, 'r')
    except:
        raise MTpyError_inputarguments('File not readable:%s'%fn)

    firstline = F.readline().strip().split()
    if not firstline[0][0] == '#':
        raise MTpyError_ts_data('Time series data file does not have a proper header:%s'%fn)

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

