#!/usr/bin/env python

"""
Helper functions for the handling of configuration files 
(survey.cfg  and BIRRP.cfg style).



@UofA, 2013
(LK)

"""

#=================================================================


import numpy as np
import sys
import os
import glob
import os.path as op
import glob
import calendar
import time
import ConfigParser
import fnmatch
import shutil

import mtpy.utils.calculator as MTcc
import mtpy.processing.general as MTgn
import mtpy.utils.exceptions as MTex
import mtpy.utils.format as MTft

reload(MTgn)
reload(MTcc)
reload(MTex)

#=================================================================

def read_configfile(filename):
    """
        Read a general config file and return the content as dictionary.

        Config files without sections or only DEFAULT section 
        -> return dictionary
        
        Config files with sections 
        -> return nested dictionary (main level keys are section heads)
        
        Config files with sections as well as section-less entries 
        -> return nested dictionary, which includes a 'DEFAULT' key
    """

    #generate config parser instance
    configobject = ConfigParser.ConfigParser()
    
    #check, if file is present
    if not op.isfile(filename):
        raise MTex.MTpyError_inputarguments( 'File does not exist: {0}'.format(filename))

    # try to parse file - exit, if not a config file
    try:
        configobject.read(filename)
    except:
        raise MTex.MTpyError_inputarguments( 'File is not a proper '
                        'configuration file: {0}'.format(filename) )


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

    The output dictionary keys are station names (capitalised), 
    the values are (sub-)dictionaries. The configuration file must contain
    sections for all stations, each containing all mandatory keywords:
    
    - latitude (deg)
    - longitude (deg)
    - elevation (in meters)
    - sampling_interval (in seconds)
    - station_type (MT, E, B)

    Not mandatory, but recommended
    - declination (in degrees, positive to East) - this is set to '0.0', if omitted

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
    - B_Xaxis_azimuth (degrees)
    - B_Yaxis_azimuth (degrees)

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
                                'B_instrument_amplification',
                                'B_Xaxis_azimuth',
                                'B_Yaxis_azimuth'
                              ]


    dict_of_allowed_values_efield = {'E_logger_type':['edl','elogger', 'zen'] ,
                                    'E_logger_gain': ['low', 'verylow','high', 
                                                      0.4, 1, 10, 11, 2, 4, 
                                                      8, 16, 32, 64],
                                    'E_instrument_type':['electrodes', 
                                                         'cu-cuso4 electrodes'],
                                    'E_instrument_amplification':[1,10]
                                    }
    
    dict_of_allowed_values_bfield = {'B_logger_type':['edl', 'zen'] ,
                                    'B_logger_gain': ['low', 'verylow','high',
                                                      0.4, 1, 10, 2, 4, 
                                                      8, 16, 32, 64],
                                    'B_instrument_type':['fluxgate', 'coil']
                                    }

    list_of_station_types = ['mt','e','b']


    error_counter = 0

    #generate config parser instance
    configobject = ConfigParser.ConfigParser()

    #check, if file is present
    if not op.isfile(filename):
        raise MTex.MTpyError_inputarguments( 'File does not'
                            ' exist: {0}'.format(filename) )


    # try to parse file - exit, if not a config file
    try:
        configobject.read(filename)
    except:
        raise MTex.MTpyError_inputarguments( 'File is not a '
                    'proper configuration file: {0}'.format(filename) )

    #obtain dict of dicts containing the input file's sections (station names)
    #excludes DEFAULT section and key-value pairs without section header
    configobject_dict = configobject._sections

    #initialise the output dictionary
    config_dict = {}

    #loop over the sections (stations) of the config file
    for station in configobject_dict:
        #read in the sub-dictionary for the current station - bringing all keys
        #to lowercase!
        temp_dict_in = dict((k.lower(),v) 
                            for k, v in configobject_dict[station].items())

        #initialise output sub-directory for current station 
        stationdict = temp_dict_in

        #stationnames are uppercase in MTpy
        stationname = station.upper()
        stationdict['station'] = stationname

        #check for presence of all mandatory keywords for the current station
        #case insensitive - allow for short forms 'lat', 'lon', and 'ele'
        for req_keyword in list_of_required_keywords:
            if req_keyword.lower() in temp_dict_in.keys():
                stationdict[req_keyword.lower()] = \
                                      temp_dict_in[req_keyword.lower()].lower()
            elif req_keyword in ['latitude', 'longitude', 'elevation']:
                if req_keyword[:3] in temp_dict_in.keys():
                    stationdict[req_keyword] = temp_dict_in[req_keyword[:3]]
            else:  
                print 'Station {0} - keyword {1} missing'.format(stationname,
                                                                 req_keyword)
                error_counter += 1
                continue

        #check format of lat/lon - convert to degrees, if given in 
        #(deg,min,sec)-triple
        for coordinate in ['latitude', 'longitude', 'elevation']:
            value = stationdict[coordinate]
            try:
                new_value = MTft._assert_position_format(coordinate,value)
            except:
                raise MTex.MTpyError_config_file('Error - wrong '
                        'coordinate format for station {0}'.format(stationname))
            stationdict[coordinate] = new_value

        if not stationdict['station_type'] in list_of_station_types:
            raise MTex.MTpyError_config_file( 'Station type not valid' )


        if stationdict['station_type'] in ['mt','e']:
            #check for required electric field parameters
            for req_keyword in list_of_efield_keywords:
                if req_keyword.lower() in temp_dict_in.keys():
                    stationdict[req_keyword.lower()] = \
                                      temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print 'Station {0} - keyword {1} missing'.format(stationname,
                                                                  req_keyword)
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_efield)
            

        if stationdict['station_type'] in ['mt','b']:
            #check for required magnetic field parameters
            for req_keyword in list_of_bfield_keywords:
                if req_keyword.lower() in temp_dict_in.keys():
                    stationdict[req_keyword.lower()] = \
                                     temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print 'Station {0} - keyword {1} missing'.format(stationname,
                                                                     req_keyword)
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_bfield)
            

        #add the station's sub-dictionary to the config dictionary
        config_dict[stationname] = stationdict

    #re-loop for setting up correct remote reference station information :
    #if rem.ref. station key is present, its information must be contained 
    #in the config file!
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
                    stationdict['rr_station_latitude']  = float(
                                            stationdict['rr_station_lat'] )
                    stationdict['rr_station_longitude'] = float(
                                            stationdict['rr_station_lon'] )
                    stationdict['rr_station_elevation'] = float(
                                            stationdict['rr_station_ele'] )                 

                except:
                    try:
                        #read from other config dict entry
                        stationdict['rr_station_latitude'] = \
                                          config_dict[rem_station]['latitude']
                        stationdict['rr_station_longitude'] = \
                                         config_dict[rem_station]['longitude']
                        stationdict['rr_station_elevation'] = \
                                         config_dict[rem_station]['elevation']

                    except:
                        #if finally failed to read rr_station info,\
                        #set rr_station back to None
                        stationdict['rr_station'] = None
                        stationdict['rr_station_latitude'] = None
                        stationdict['rr_station_longitude'] = None
                        stationdict['rr_station_elevation'] = None

        #check consistency of coordinates, if rr_station is present
        if stationdict['rr_station'] != None:
            try:
                stationdict['rr_station_latitude'] = \
                            MTft._assert_position_format(
                                'latitude',stationdict['rr_station_latitude'])
                stationdict['rr_station_longitude'] = \
                            MTft._assert_position_format(
                                'longitude',stationdict['rr_station_longitude'])
                stationdict['rr_station_elevation'] = \
                            MTft._assert_position_format(
                                'elevation',stationdict['rr_station_elevation'])

            except:
                print 'Problem with remote reference station ({0}) -'
                ' remote reference ({1}) coordinates invalid -'
                ' remote reference set to None'.format(station, 
                                                       stationdict['rr_station'])

                stationdict['rr_station'] = None
                stationdict['rr_station_latitude'] = None
                stationdict['rr_station_longitude'] = None
                stationdict['rr_station_elevation'] = None


    if error_counter != 0:
        print 'Could not read all mandatory sections and options'\
                ' in config file - found {0} errors - check configuration'\
                ' file before continuing!'.format(error_counter)
    
    else:
        return config_dict

#=================================================================

def write_dict_to_configfile(dictionary, output_filename):
    """
    Write a dictionary into a configuration file.

    The dictionary can contain pure key-value pairs as well as a 
    level-1 nested dictionary. In the first case, the entries are 
    stored in a 'DEFAULT' section. In the latter case, the dictionary 
    keys are taken as section heads and the sub-dictionaries key-value 
    pairs fill up the respective section  

    """

    configobject = ConfigParser.ConfigParser()

    #check for nested dictionary - 
    #if the dict entry is a key-value pair, it's stored in a section with head 'DEFAULT' 
    #otherwise, the dict key is taken as section header
    for key, val in sorted(dictionary.items()):
        try:
            for subkey, subval in sorted(val.items()):
                sectionhead = key
                if not configobject.has_section(sectionhead):
                    configobject.add_section(sectionhead)
                configobject.set(sectionhead, subkey, subval)

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
            raise MTex.MTpyError_config_file('Config file error --'
                ' key {0}, value {1} not valid'.format(key, value2validate) )

#==============================================================================
def read_survey_txt_file(survey_file, delimiter='\t'):
    """
    read survey file and return a dictionary of dictionaries where the first
    nested dictionary is keyed by the station name.  Each station dictionarly
    includes all the information input in the survey file with keywords 
    verbatim as the headers in survey file, all lower case.  
    
    *Must be included in survey file*
    ================= =========================================================
    key word           description
    ================= =========================================================
    station           station name
    lat(itude)        latitude (decimal degrees is best)
    long(itude)       longitude (decimal degrees is best)
    elev(ation)       elevation (in meters)
    ex/E_Xaxis_length dipole length in north direction (in meters)
    ey/E_Yaxis_length dipole length in east direction (in meters)
    E_Xaxis_azimuth   orientaion of Ex (degrees)
    E_Yaxis_azimuth   orientaion of Ey (degrees)


    sampling_interval sampling interval in seconds
    hx                coil number in north direction for calibration
    hy                coil number in east direction for calibration
    hz                coil number in vertical direction for calibration
    date              date of deployment
    notes             any notes that might help later
    station_type      type of data collected  (MT, E, B)
    declination       declination in degrees (N = 0 and East = 90)
    ================= =========================================================
     
    *Information on E-field data*:
    ========================== ================================================
    key word                    description
    ========================== ================================================
    E_logger_type               type of data logger used to record data 
    E_logger_gain               factor/gain level
    E_instrument_type           type of electrodes used
    E_instrument_amplification  applied amplification factor
    E_Xaxis_azimuth             orientaion of Ex (degrees)
    E_Xaxis_length              length of dipole for Ex (in meters)
    E_Yaxis_azimuth             orientaion of Ey (degrees)
    E_Yaxis_length              length of dipole for Ey (in meters)
    ========================== ================================================
   
   *Information on B-field data*:
    ========================== ================================================
    key word                    description
    ========================== ================================================
    B_logger_type              type of data logger used to record data
    B_logger_gain              factor/gain level
    B_instrument_type          type of magnetometer used (coil, fluxgate)
    B_instrument_amplification applied amplification factor
    B_Xaxis_azimuth            orientation of Bx (degrees)
    B_Yaxis_azimuth            orientation of By (degrees)
    ================= =========================================================
    
    Arguments:
    -----------
        **survey_file** : string (full path to file)
        
    Outputs:
    ---------
        **survey_lst** : list
                         list of dictionaries with key words the same as the
                         headers in survey file, all lower case
    """                  
        
    with open(survey_file, 'r') as sfid:
        slines = sfid.readlines()

    skeys = slines[0].rstrip()
    skeys = skeys.split(delimiter)
    survey_dict = {}
    
    for ss, sline in enumerate(slines[1:]):
        sstr = sline.rstrip()
        sstr = sstr.split(delimiter)
        if len(sstr)>1:
            sdict={}
            for kk, skey in enumerate(skeys):
                sstr[kk] = sstr[kk].replace('"','')
                sdict[skey.lower()] = sstr[kk]
        try:
            survey_dict[sdict['station']] = sdict
        except KeyError:
            try: 
                survey_dict[sdict['station_name']] = sdict
            except KeyError:
                survey_dict['MT{0:03}'.format(ss)] = sdict
        
    sfid.close()
    
    return survey_dict
    
#==============================================================================
def write_config_from_survey_txt_file(survey_file, save_path=None, 
                                      delimiter='\t'):
    """
    write a survey configuration file from a survey txt file 
    
    Arguments:
    ----------
        **survey_file** : string
                          full path to survey text file.  
                          See read_survey_txt_file for the assumed header 
                          information.
                          
        **save_path** : string
                        directory or full path to save file to.  
                        If save_path = None, then file saved as 
                        os.path.join(os.path.dirname(survey_file,
                                            os.path.basename(survey_file).cfg)
                                            
    Outputs:
    ---------
        **cfg_fn** : string
                    full path to saved config file
    """
    
    survey_dict = read_survey_txt_file(survey_file, delimiter=delimiter)
    
    #get the filename to save to
    if save_path is None:
        save_dir = os.path.dirname(survey_file)
        save_fn = os.path.splitext(os.path.basename(survey_file))[0]+'.cfg'
        save_path = os.path.join(save_dir, save_fn)
    elif os.path.isfile(save_path):
        pass
    elif os.path.isdir(save_path):
        save_fn = os.path.splitext(os.path.basename(survey_file))[0]+'.cfg'
        save_path = os.path.join(save_path, save_fn)
    
    #write the config file
    write_dict_to_configfile(survey_dict, save_path)
    
    return save_path






