#!/usr/bin/env python

"""
Helper functions for the handling of configuration files 
(survey.cfg  and BIRRP.cfg style).



@UofA, 2013
(LK)

"""

#=================================================================

import sys
import os
import os.path as op
import configparser
import copy
import io
import mtpy.utils.exceptions as MTex
import mtpy.utils.gis_tools as gis_tools
#=================================================================

list_of_required_keywords = ['latitude',
                            'longitude',
                            'elevation',
                            'sampling_interval',
                            'station_type'
                            ]
list_of_required_keywords_short = ['lat',
                                   'lon',
                                   'elev',
                                   'sampling',
                                   'type'
                                  ]

list_of_keyword_defaults_general = [0.,
                                    0.,
                                    0.,
                                    1.,
                                    'mt'
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

list_of_keyword_defaults_efield = ['edl',
                                    1,
                                    'electrodes',
                                    1.,
                                    0.,
                                    50.,
                                    90.,
                                    50.
                                    ]

list_of_bfield_keywords = [ 'B_logger_type',
                            'B_logger_gain',
                            'B_instrument_type',
                            'B_instrument_amplification',
                            'B_Xaxis_azimuth',
                            'B_Yaxis_azimuth'
                          ]


list_of_keyword_defaults_bfield = ['edl',
                                    1,
                                    'coil',
                                    1.,
                                    0.,
                                    90.
                                    ]


dict_of_allowed_values_efield = {'E_logger_type':['edl','elogger', 'zen','qel'] ,
                                'E_logger_gain': ['low', 'verylow','high', 
                                                  0.4, 1, 10, 11, 2, 4, 
                                                  8, 16, 32, 64],
                                'E_instrument_type':['electrodes','dipole', 
                                                     'cu-cuso4 electrodes',
                                                     'cuso4_electrodes',
                                                     'pbcl2_electrodes'],
                                'E_instrument_amplification':[1,10,11]
                                }

dict_of_allowed_values_bfield = {'B_logger_type':['edl', 'zen','qel_blogger'] ,
                                'B_logger_gain': ['low', 'verylow','high',
                                                  0.4, 1, 10, 2, 4, 
                                                  8, 16, 32, 64],
                                'B_instrument_type':['fluxgate', 'coil','coils']
                                }

list_of_station_types = ['mt','e','b','qe','qb']


#=================================================================

def read_configfile(filename):
    """
        Read a general config file and return the content as dictionary.

        Config files without sections or only DEFAULT section 
        -> return dictionary
        
        Config files with sections 
        -> return nested dictionary (main level keys are section heads)
        
        Config files with sections as well as section-less entries 
        -> return nested dictionary, which includes a top level 'DEFAULT' key
    """

    
    #check, if file is present
    if not op.isfile(filename):
        raise MTex.MTpyError_inputarguments( 'File does not exist: {0}'.format(filename))

    # try to parse file - exit, if not a config file
    try:
        #generate config parser instance
        configobject = configparser.SafeConfigParser()
        #do NOT ask, why it does not work with reading from filename directly...:
        with open(filename) as F:
            d = F.read()
        FH = io.StringIO(d)
        configobject.readfp(d)#filename)
    except:
        try:
            dummy_String = '[DEFAULT]\n' + open(filename, 'r').read()
            FH = io.StringIO(dummy_String)
            #generate config parser instance
            configobject = configparser.SafeConfigParser()
            configobject.readfp(FH)
        except:
            raise MTex.MTpyError_inputarguments( 'File is not a proper '
                                    'configuration file: {0}'.format(filename) )

    config_dict = configobject._sections      

    if len (list(config_dict.keys())) != 0:
        defaults = configobject.defaults()
        if len(list(defaults.keys())) != 0:
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
    - station_type (MT, (Q)E, (Q)B)

    Not mandatory, but recommended
    - declination (in degrees, positive to East) - this is set to '0.0', if omitted

    Depending on the type of station the following entries are required.

    E-field recorded:

    - E_logger_type ('edl'/'elogger'/'qel')
    - E_logger_gain (factor/gain-level)
    - E_instrument_type ('electrodes'/'dipole')
    - E_instrument_amplification (applied amplification factor)
    - E_Xaxis_azimuth (degrees)
    - E_Xaxis_length (in meters)
    - E_Yaxis_azimuth (degrees)
    - E_Yaxis_length (in meters)

    B-field recorded:
  
    - B_logger_type ('edl'/'qel_blogger')
    - B_logger_gain (factor/gain level)
    - B_instrument_type ('coil(s)', 'fluxgate')
    - B_instrument_amplification (applied amplification factor)
    - B_Xaxis_azimuth (degrees)
    - B_Yaxis_azimuth (degrees)


    A global section can be used to include parameters for all stations.
    The name of the section must be one of:

        global/main/default/general 

    """





    error_counter = 0

    #generate config parser instance
    configobject = configparser.ConfigParser()

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
        temp_dict_in = dict((k.lower(), v) 
                            for k, v in list(configobject_dict[station].items()))

        #initialise output sub-directory for current station 
        stationdict = temp_dict_in

        #stationnames are uppercase in MTpy
        stationname = station
        if stationname in ['GLOBAL','MAIN','DEFAULT','GENERAL']:
            stationname = 'GLOBAL'

        stationdict['station'] = stationname

        #add the station's sub-dictionary to the config dictionary
        config_dict[stationname] = stationdict


    # Check if a global section is present
    if 'GLOBAL' in config_dict: 
        globaldict = config_dict['GLOBAL']
    else:
        #set defaults for location
        globaldict={}
    # for i in ['latitude', 'longitude', 'elevation']:
    #     #skip if values are present
    #     if i in globaldict.keys() or i[:3] in globaldict.keys():
    #         continue
    #     #otherwise set defaults
    #     globaldict[i] = 0
        

    #remove other general sections to avoid redundancy
    for i in ['MAIN','DEFAULT','GENERAL']:
        if i in config_dict:
            dummy = config_dict.pop(i)

    # RE-loop to check for each station if required keywords are present,
    # if not if they can be pulled from the global section  
    
    #============================================================
    # local function definition
    def fromglobals(key,stationdict,globaldict):
        """
            Check if stationdict contains key. 
            If not search for key in global dict and add it to station dict.

            Return if global dict is not defined.
            Return True if key was present in either dictionary, False if not.

        """

        if key in list(stationdict.keys()):
            return True, stationdict.get(key)

        if globaldict is None or len(globaldict) == 0:
            return False, None


        if key in globaldict:
            stationdict[key] = globaldict[key]
            return True,globaldict.get(key)

        return False, None

    #============================================================


    for station in sorted(config_dict):
        #do not alter the global section
        if station == 'GLOBAL':
            continue
        
        stationdict =  config_dict[station]
        


        #check for presence of all mandatory keywords for the current station
        #case insensitive - allow for short forms 'sampling', 'lat', 'lon', and 'elev'
        for idx,req_keyword in enumerate(list_of_required_keywords):
            shortform = list_of_required_keywords_short[idx]

            try:
                found = False
                #import ipdb
                #ipdb.set_trace()
                if fromglobals(req_keyword,stationdict,globaldict)[0] is False:
                    #try short form instead
                    found,value = fromglobals(shortform,stationdict,globaldict)
                    #print shortform,value

                    if found is True:
                        stationdict[req_keyword] = value
  
                else:  
                    found = True

                if found is False:
                    print('Station {0} - keyword {1} missing'.format(stationname,
                                                                     req_keyword))
                    error_counter += 1
                    raise Exception

                if req_keyword in ['elevation','latitude', 'longitude']:
                    #check format of lat/lon - convert to degrees, if given in 
                    #(deg,min,sec)-triple#assert correct format
                    value = stationdict[req_keyword]
                    try:
                        if req_keyword in 'latitude':
                            new_value = gis_tools.assert_lat_value(value)
                        elif req_keyword in 'longitude':
                            new_value = gis_tools.assert_lon_value(value)
                        elif req_keyword in 'elevation':
                            new_value = gis_tools.assert_elevation_value(value)
                    except:
                        raise MTex.MTpyError_config_file('Error - wrong '
                                'coordinate format for station {0}'.format(stationname))
                    
                    stationdict[req_keyword] = new_value



            except:
                raise
                print('Missing information on station {0} in config file'\
                        ' - setting default (dummy) value'.format(station))
                stationdict[req_keyword] = list_of_keyword_defaults_general[idx]
            
            #to avoid duplicates remove the now obsolete short form from 
            #the station dictionary
            dummy = stationdict.pop(shortform,None)


        if not stationdict['station_type'] in list_of_station_types:
            raise MTex.MTpyError_config_file( 'Station type not valid' )


        if stationdict['station_type'] in ['mt','e']:
            #check for required electric field parameters - not done for QEL loggers yet
            for req_keyword in list_of_efield_keywords:
                if req_keyword.lower() in list(temp_dict_in.keys()):
                    stationdict[req_keyword.lower()] = \
                                      temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print('Station {0} - keyword {1} missing'.format(stationname,
                                                                  req_keyword))
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_efield)
            

        if stationdict['station_type'] in ['mt','b']:
            #check for required magnetic field parameters
            for req_keyword in list_of_bfield_keywords:
                if req_keyword.lower() in list(temp_dict_in.keys()):
                    stationdict[req_keyword.lower()] = \
                                     temp_dict_in[req_keyword.lower()].lower()
                else:  
                    print('Station {0} - keyword {1} missing'.format(stationname,
                                                                     req_keyword))
                    error_counter += 1
                    continue

            _validate_dictionary(stationdict,dict_of_allowed_values_bfield)
            


    #re-loop for setting up correct remote reference station information :
    #if rem.ref. station key is present, its information must be contained 
    #in the same config file!
    for station in config_dict.keys():
        stationdict = config_dict[station]
        if 'rr_station' not in stationdict:
            continue

        #stationdict['rr_station'] = None
        stationdict['rr_latitude'] = None
        stationdict['rr_longitude'] = None
        stationdict['rr_elevation'] = None


        rem_station = stationdict['rr_station'] 
        try:
            #check, if values are contained in dict 
            float(stationdict['rr_latitude'] )
            float(stationdict['rr_longitude'])
            float(stationdict['rr_elevation'])
        except:
            try:
                #check for shortened form
                stationdict['rr_latitude']  = float(stationdict['rr_lat'] )
                stationdict['rr_longitude'] = float(stationdict['rr_lon'] )
                stationdict['rr_elevation'] = float(stationdict['rr_elev'] )                 

            except:
                try:
                    #read from other config dict entry
                    stationdict['rr_latitude'] = \
                                      config_dict[rem_station]['latitude']
                    stationdict['rr_longitude'] = \
                                     config_dict[rem_station]['longitude']
                    stationdict['rr_elevation'] = \
                                     config_dict[rem_station]['elevation']

                except:
                    #if finally failed to read rr_station info,\
                    #set rr_station back to None
                    stationdict['rr_station'] = None
                    stationdict['rr_latitude'] = None
                    stationdict['rr_longitude'] = None
                    stationdict['rr_elevation'] = None

        #check consistency of coordinates, if rr_station is present
        if stationdict['rr_station'] != None:
            try:
                stationdict['rr_latitude'] = \
                            gis_tools.assert_lat_value(stationdict['rr_latitude'])
                stationdict['rr_longitude'] = \
                            gis_tools.assert_lon_value(stationdict['rr_longitude'])
                stationdict['rr_elevation'] = \
                            gis_tools.assert_elevation_value(stationdict['rr_elevation'])

            except:
                print('Problem with remote reference station ({0}) -')
                ' remote reference ({1}) coordinates invalid -'
                ' remote reference set to None'.format(station, 
                                                       stationdict['rr_station'])

                stationdict['rr_station'] = None
                stationdict['rr_latitude'] = None
                stationdict['rr_longitude'] = None
                stationdict['rr_elevation'] = None


    if error_counter != 0:
        print('Could not read all mandatory sections and options'\
                ' in config file - found {0} errors - check configuration'\
                ' file before continue!'.format(error_counter))
        answer = 5
        while not answer in ['y','n']:
            answer = input('\n\tDo you want to continue anyway? (y/n)')
            try:
                answer = answer.strip().lower()[0]
            except:
                continue
            if answer == 'n':
                sys.exit() 
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

    configobject = configparser.ConfigParser()

    #check for nested dictionary - 
    #if the dict entry is a key-value pair, it's stored in a section with head 'DEFAULT' 
    #otherwise, the dict key is taken as section header
    for key, val in sorted(dictionary.items()):
        try:
            for subkey, subval in sorted(val.items()):
                sectionhead = key
                if not configobject.has_section(sectionhead):
                    configobject.add_section(sectionhead)
                configobject.set(sectionhead, subkey, str(subval))

        except KeyError:
            #if not configobject.has_section('DEFAULT'):
            #    configobject.add_section('')
            configobject.set('',key, str(val))


    with open(output_filename, 'w') as F:
        configobject.write(F)
    

#=================================================================
 
def _validate_dictionary(dict2validate,referencedict):

    """Check, if there are lists of allowed entries for all 
    keys of the current dictionary. If yes, test, if the current 
    value is among the allowed ones. 
    """

    for key, value in list(dict2validate.items()):
        #make everything to strings - easier to compare
        #in case of numbers, make to float first

        try:
            allowed_vals = referencedict[key]
        except:
            try:
                key = key.lower()
                allowed_vals = referencedict[key]
            except:
                #no reference entry found - skip key
                continue

        tmp = []
        #allowed values must be given as a list (iterable)!!
        for i in allowed_vals:
            try: 
                tmp.append(str(float(i)))
            except:
                tmp.append(str(i))  
        tmp = [i.lower() for  i in tmp]
        allowed_vals = tmp

        #compare case-insensitive
        value = value.lower()

        if not value in allowed_vals:
            raise MTex.MTpyError_config_file('Config file error --'
                ' key {0}, value {1} not valid'.format(key, value) )

#==============================================================================
def read_survey_txt_file(survey_file, delimiter=None):
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
    
    slines = [i.replace('"','') for i in slines]

    
    skeys = slines[0].rstrip()
    if delimiter is not None:
        skeys = skeys.split(delimiter)
    else:
        skeys = skeys.split()

    skeys = [i.strip().replace(' ','_') for i in skeys]
       
    survey_dict = {}
    #print skeys, len(skeys)

    for ss, sline in enumerate(slines[1:]):

        sstr = sline.strip()
        if sstr[0]=='#':
            continue

        if delimiter is not None:
            sstr = sstr.split(delimiter)
        else:
            sstr = sstr.split()
        #print sstr
        #get rid of quotations
        sstr = [i.replace('"','') for i in sstr]  
        #get rid of spaces
        sstr = [i.replace(' ','_') for i in sstr]  
        #print sstr,len(sstr)
       
        if len(sstr) != len(skeys):
            print('cannot read line {0} - wrong number of entries - need {2}\
                                                    '.format(ss+2,len(skeys)))
            continue

        
        sdict={}
        #set default values for mandatory entries:
        sdict['E_Xaxis_azimuth'] = 0
        sdict['E_Yaxis_azimuth'] = 90
        sdict['B_Xaxis_azimuth'] = 0
        sdict['B_Yaxis_azimuth'] = 90
        sdict['station_type'] = 'MT'
        sdict['declination'] = 0.
        sdict['sampling_interval'] = 0
        sdict['E_instrument_amplification'] = 1.
        sdict['B_instrument_amplification'] = 1.
        sdict['E_logger_gain'] = 1.
        sdict['B_logger_gain'] = 1.
        sdict['B_instrument_type'] = 'coil'
        sdict['E_instrument_type'] = 'electrodes'
        sdict['E_logger_type'] = 'edl'
        sdict['B_logger_type'] = 'edl'


        #fill dictionary with given values
        for kk, skey in enumerate(skeys):
            #get rid of quotations
            skey.replace('"','')
            #get rid of blank spaces in keys
            skey.replace(' ','_')

            #do not include empty entries
            if len(sstr[kk])>0:
                #sstr[kk] = sstr[kk].replace('"','')
                sdict[skey.lower()] = sstr[kk]
        
        #print sorted(sdict)
        # print sdict['sampling_interval']
        #sys.exit()
        #assigne values to the standard keys
        for key in list(sdict.keys()):

            if key.lower() in ['ex','e_xaxis_length'] :
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['E_Xaxis_length'] = val
            if key.lower() in ['ey','e_yaxis_length'] :
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['E_Yaxis_length'] = val

            if key.lower() == 'station':
                sdict[key] = sdict[key].upper()

            if key.lower() in ['df', 'sampling_rate','sampling']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['sampling_interval'] = 1./float(val)

            if key.lower() in ['dt', 'sampling_interval']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['sampling_interval'] = float(val)

            if key.lower() == 'dlgain':
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['B_logger_gain'] = val
                sdict['E_logger_gain'] = val
                

            if key.lower() in ['b_logger_gain']:
                sdict['B_logger_gain'] = sdict[key]

            if key.lower() in ['e_logger_gain']:
                sdict['E_logger_gain'] = sdict[key]

            if key.lower() in ['e_logger_type']:
                sdict['E_logger_type'] = sdict[key]
            if key.lower() in ['b_logger_type']:
                sdict['B_logger_type'] = sdict[key]


            if key.lower() in ['egain', 'e_instrument_amplification']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['E_instrument_amplification'] = val
            
            if key.lower() in ['bgain', 'b_instrument_amplification']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['B_instrument_amplification'] = val


            if key.lower() in ['magtype','b_instrument_type']:
                if sdict[key].lower() in ['bb','coil','coils']: 
                    sdict['B_instrument_type'] = 'coil'
                if sdict[key].lower() in ['lp','fluxgate','fg']: 
                    sdict['B_instrument_type'] = 'fluxgate'
                sdict.pop(key)

            if key.lower() in ['e_instrument_type']:
                if sdict[key].lower() in ['electrode','electrodes']: 
                    sdict['E_instrument_type'] = 'electrodes'
                if (sdict[key].lower().find('lead') >= 0) or (sdict[key].lower().find('pb') >= 0): 
                    sdict['E_instrument_type'] = 'pbcl_electrodes'
                sdict.pop(key)

            if key.lower() in ['declination','decl']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['declination'] = val


            if key.lower() in ['lat','latitude']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['latitude'] = val

            if key.lower() in ['lon','long','longitude']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['longitude'] = val

            if key.lower() in ['ele','elev','elevation','height']:
                val = copy.copy(sdict[key])
                sdict.pop(key)
                sdict['elevation'] = val


        try:
            survey_dict[sdict['station']] = sdict
        except KeyError:
            try: 
                survey_dict[sdict['station_name']] = sdict
            except KeyError:
                survey_dict['MT{0:03}'.format(ss)] = sdict

    return survey_dict
    
#==============================================================================
def write_config_from_survey_txt_file(survey_file, save_name=None, 
                                      delimiter='\t'):
    """
    write a survey configuration file from a survey txt file 
    
    Arguments:
    ----------
        **survey_file** : string
                          full path to survey text file.  
                          See read_survey_txt_file for the assumed header 
                          information.
                          
        **save_name** : string
                        name to save file to.  
                        If save_name = None, then file saved as 
                        os.path.join(os.path.dirname(survey_file,
                                            os.path.basename(survey_file).cfg)
                                            
    Outputs:
    ---------
        **cfg_fn** : string
                    full path to saved config file
    """
    
    survey_dict = read_survey_txt_file(survey_file, delimiter=delimiter)
    
    #get the filename to save to
    if save_name is None:
        save_dir = os.path.dirname(survey_file)
        save_fn = os.path.splitext(os.path.basename(survey_file))[0]+'.cfg'
        save_name = os.path.join(save_dir, save_fn)
    elif os.path.isfile(save_name):
        pass
    elif os.path.isdir(save_name):
        save_fn = os.path.splitext(os.path.basename(survey_file))[0]+'.cfg'
        save_name = os.path.join(save_name, save_fn)

    if not save_name.lower().endswith('.cfg'):
        save_name += '.cfg'

    
    #write the config file
    write_dict_to_configfile(survey_dict, save_name)
    
    return save_name






