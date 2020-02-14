# -*- coding: utf-8 -*-
"""
===============
NIMS
===============

    * deals with reading in NIMS DATA.BIN files
    
Created on Thu Oct 31 10:03:20 2019

@author: jpeacock
"""

# =============================================================================
# Imports
# =============================================================================
import os
import numpy as np
import struct
import datetime
import dateutil

import pandas as pd
import logging

from mtpy.core import ts

from matplotlib import pyplot as plt

### setup logger
logging.basicConfig(filename='ReadNIMSData.log', 
                    filemode='w', 
                    level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')


# =============================================================================
# Exceptions
# =============================================================================
class NIMSError(Exception):
    pass

class GPSError(Exception):
    pass

class ResponseError(Exception):
    pass

# =============================================================================
# class objects
# =============================================================================
class GPS(object):
    """
    class to parse GPS stamp from the NIMS
    
    Depending on the type of Stamp different attributes will be filled.
    
    GPRMC has full date and time information and declination
    GPGGA has elevation data
    
    .. note:: GPGGA date is set to 1980-01-01 so that the time can be estimated.
              Should use GPRMC for accurate date/time information.  
    """
    
    def __init__(self, gps_string, index=0):
        
        self.gps_string = gps_string
        self.index = index
        self._time = None
        self._date = '010180'
        self._latitude = None
        self._latitude_hemisphere = None
        self._longitude = None
        self._longitude_hemisphere = None
        self._declination = None
        self._declination_hemisphere = None
        self._elevation = None
        self.valid = False
        self.elevation_units = 'meters'
        
        self.type_dict = {'gprmc':{0:'type', 
                                   1:'time', 
                                   2:'fix',
                                   3:'latitude',
                                   4:'latitude_hemisphere',
                                   5:'longitude',
                                   6:'longitude_hemisphere',
                                   7:'skip',
                                   8:'skip',
                                   9:'date',
                                   10:'declination',
                                   11:'declination_hemisphere',
                                   'length':[12],
                                   'type':0, 
                                   'time':1, 
                                   'fix':2,
                                   'latitude':3,
                                   'latitude_hemisphere':4,
                                   'longitude':5,
                                   'longitude_hemisphere':6,
                                   'date':9,
                                   'declination':10},
                          'gpgga':{0:'type', 
                                   1:'time', 
                                   2:'latitude',
                                   3:'latitude_hemisphere',
                                   4:'longitude',
                                   5:'longitude_hemisphere',
                                   6:'var_01',
                                   7:'var_02',
                                   8:'var_03',
                                   9:'elevation',
                                   10:'elevation_units',
                                   11:'elevation_error',
                                   12:'elevation_error_units',
                                   13:'null_01',
                                   14:'null_02',
                                   'length':[14,15],
                                   'type':0, 
                                   'time':1, 
                                   'latitude':2,
                                   'latitude_hemisphere':3,
                                   'longitude':4,
                                   'longitude_hemisphere':5,
                                   'elevation':9,
                                   'elevation_units':10,
                                   'elevation_error':11,
                                   'elevation_error_units':12}}
        self.parse_gps_string(self.gps_string)
        
    def validate_gps_string(self, gps_string):
        """
        make sure the string is valid, remove any binary numbers and find
        the end of the string as '*'
        
        :param string gps_string: raw GPS string to be validated
        
        :returns: validated string or None if there is something wrong
        """
        for replace_str in [b'\xd9', b'\xc7', b'\xcc']:
            gps_string = gps_string.replace(replace_str, b'')
            
        ### sometimes the end is set with a zero for some reason
        gps_string = gps_string.replace(b'\x00', b'*')
        
        if gps_string.find(b'*') < 0:
            logging.error('GPSError: No end to stamp {0}'.format(gps_string))
        else:
            try:
                gps_string = gps_string[0:gps_string.find(b'*')].decode()
                return gps_string
            except UnicodeDecodeError:
                logging.error('GPSError: stamp not correct format, {0}'.format(gps_string))
                return None
        
    def parse_gps_string(self, gps_string):
        """
        Parse a raw gps string from the NIMS and set appropriate attributes.
        GPS string will first be validated, then parsed. 
        
        :param string gps_string: raw GPS string to be parsed
        """
        gps_string = self.validate_gps_string(gps_string)
        if gps_string is None:
            self.valid = False
            return
        
        if isinstance(gps_string, bytes):
            gps_list = gps_string.strip().split(b',')
            gps_list = [value.decode() for value in gps_list]
        else:
            gps_list = gps_string.strip().split(',')
        
        ### validate the gps list to make sure it is usable
        gps_list, error_list = self.validate_gps_list(gps_list)
        if len(error_list) > 0:
            for error in error_list:
                logging.error('GPSError:' + error)
        if gps_list is None:
            return

        attr_dict = self.type_dict[gps_list[0].lower()]
            
        for index, value in enumerate(gps_list):
            setattr(self, '_'+attr_dict[index], value)
            
        if None not in gps_list:
            self.valid = True
            self.gps_string = gps_string
                  
    def validate_gps_list(self, gps_list):
        """
        check to make sure the gps stamp is the correct format
        """
        error_list = []
        try:
            gps_list = self._validate_gps_type(gps_list)
        except GPSError as error:
            error_list.append(error.args[0])
            return None, error_list
        
        ### get the string type
        g_type = gps_list[0].lower()
        
        ### first check the length, if it is not the proper length then
        ### return, cause you never know if everything else is correct
        try:
            self._validate_list_length(gps_list)
        except GPSError as error:
            error_list.append(error.args[0])    
            return None, error_list
        
        try:
            gps_list[self.type_dict[g_type]['time']] = \
                  self._validate_time(gps_list[self.type_dict[g_type]['time']])
        except GPSError as error:
            error_list.append(error.args[0])
            gps_list[self.type_dict[g_type]['time']] = None
            
        try:
            gps_list[self.type_dict[g_type]['latitude']] =\
                self._validate_latitude(gps_list[self.type_dict[g_type]['latitude']],
                                        gps_list[self.type_dict[g_type]['latitude_hemisphere']])
        except GPSError as error:
            error_list.append(error.args[0])
            gps_list[self.type_dict[g_type]['latitude']] = None
        
        try:
            gps_list[self.type_dict[g_type]['longitude']] =\
                self._validate_longitude(gps_list[self.type_dict[g_type]['longitude']],
                                         gps_list[self.type_dict[g_type]['longitude_hemisphere']])
        except GPSError as error:
            error_list.append(error.args[0])
            gps_list[self.type_dict[g_type]['longitude']] = None
        
        if g_type == 'gprmc':
            try:
                gps_list[self.type_dict['gprmc']['date']] =\
                    self._validate_date(gps_list[self.type_dict['gprmc']['date']])
            except GPSError as error:
                error_list.append(error.args[0])
                gps_list[self.type_dict[g_type]['date']] = None
            
        elif g_type == 'gpgga':
            try:
                gps_list[self.type_dict['gpgga']['elevation']] =\
                    self._validate_elevation(gps_list[self.type_dict['gpgga']['elevation']])
            except GPSError as error:
                error_list.append(error.args[0])
                gps_list[self.type_dict['gpgga']['elevation']] = None
            
        return gps_list, error_list
    
    def _validate_gps_type(self, gps_list):
        """Validate gps type should be gpgga or gprmc"""
        gps_type = gps_list[0].lower()
        if 'gpg' in gps_type:
            if len(gps_type) > 5:
                gps_list = ['GPGGA', gps_type[-6:]] + gps_list[1:]
            elif len(gps_type) < 5:
                gps_list[0] = 'GPGGA'
        elif 'gpr' in gps_type:
            if len(gps_type) > 5:
                gps_list = ['GPRMC', gps_type[-6:]] + gps_list[1:]
            elif len(gps_type) < 5:
                gps_list[0] = 'GPRMC'
        
        gps_type = gps_list[0].lower()
        if gps_type not in ['gpgga', 'gprmc']:
            raise GPSError('GPS String type not correct.  '+\
                              'Expect GPGGA or GPRMC, got {0}'.format(gps_type.upper()))
            
        return gps_list
    
    def _validate_list_length(self, gps_list):
        """validate gps list length based on type of string"""
        
        gps_list_type = gps_list[0].lower()
        expected_len = self.type_dict[gps_list_type]['length']
        if len(gps_list) not in expected_len:
            raise GPSError('GPS string not correct length for {0}.  '.format(gps_list_type.upper())+\
                           'Expected {0}, got {1} \n{2}'.format(expected_len, 
                                                              len(gps_list),
                                                              ','.join(gps_list)))
            
    def _validate_time(self, time_str):
        """ validate time string, should be 6 characters long and an int """
        if len(time_str) != 6:
            raise GPSError('Lenght of time string {0} not correct.  '.format(time_str)+\
                           'Expected 6 got {0}'.format(len(time_str)))
        try:
            int(time_str)
        except ValueError:
            raise GPSError('Could not convert time string {0}'.format(time_str))
        
        return time_str
    
    def _validate_date(self, date_str):
        """ validate date string, should be 6 characters long and an int """
        if len(date_str) != 6:
            raise GPSError('Length of date string not correct {0}.  '.format(date_str)+\
                           'Expected 6 got {0}'.format(len(date_str)))
        try:
            int(date_str)
        except ValueError:
            raise GPSError('Could not convert date string {0}'.format(date_str))
        
        return date_str
    
    def _validate_latitude(self, latitude_str, hemisphere_str):
        """validate latitude, should have hemisphere string with it"""
        
        if len(latitude_str) < 8:
            raise GPSError('Latitude string should be larger than 7 characters.  '+\
                           'Got {0}'.format(len(latitude_str)))
        if len(hemisphere_str) != 1:
            raise GPSError('Latitude hemisphere should be 1 character.  '+\
                           'Got {0}'.format(len(hemisphere_str)))
        if hemisphere_str.lower() not in ['n', 's']:
            raise GPSError('Latitude hemisphere {0} not understood'.format(hemisphere_str.upper()))
        try:
            float(latitude_str)
        except ValueError:
            raise GPSError('Could not convert latitude string {0}'.format(latitude_str))
            
        return latitude_str
    
    def _validate_longitude(self, longitude_str, hemisphere_str):
        """validate longitude, should have hemisphere string with it"""
        
        if len(longitude_str) < 8:
            raise GPSError('Longitude string should be larger than 7 characters.  '+\
                           'Got {0}'.format(len(longitude_str)))
        if len(hemisphere_str) != 1:
            raise GPSError('Longitude hemisphere should be 1 character.  '+\
                           'Got {0}'.format(len(hemisphere_str)))
        if hemisphere_str.lower() not in ['e', 'w']:
            raise GPSError('Longitude hemisphere {0} not understood'.format(hemisphere_str.upper()))
        try:
            float(longitude_str)
        except ValueError:
            raise GPSError('Could not convert longitude string {0}'.format(longitude_str))
            
        return longitude_str
    
    def _validate_elevation(self, elevation_str):
        """validate elevation, check for converstion to float"""
        try:
            float(elevation_str)
        except ValueError:
            raise GPSError('Elevation could not be converted {0}'.format(elevation_str))
            
        return elevation_str
                      
    @property
    def latitude(self):
        """
        Latitude in decimal degrees, WGS84
        """
        if self._latitude is not None and self._latitude_hemisphere is not None:
            index = len(self._latitude) - 7
            lat = float(self._latitude[0:index]) + float(self._latitude[index:])/60
            if 's' in self._latitude_hemisphere.lower():
                lat *= -1
            return lat
        else:
            return 0.0
    
    @property
    def longitude(self):
        """
        Latitude in decimal degrees, WGS84
        """
        if self._longitude is not None and self._longitude_hemisphere is not None:
            index = len(self._longitude) - 7
            lon = float(self._longitude[0:index]) + float(self._longitude[index:])/60
            if 'w' in self._longitude_hemisphere.lower():
                lon *= -1
            return lon
        else:
            return 0.0
        
    @property
    def elevation(self):
        """
        elevation in meters
        """
        if self._elevation is not None:
            try:
                return float(self._elevation)
            except ValueError:
                logging.error('GPSError: Could not get elevation GPS string'+\
                              'not complete {0}'.format(self.gps_string))
        else:
            return 0.0
        
    @property
    def time_stamp(self):
        """
        return a datetime object of the time stamp
        """
        if self._time is None:
            return None
        if self._date is None:
            self._date = '010180'
        try:
            return dateutil.parser.parse('{0} {1}'.format(self._date, self._time),
                                         dayfirst=True)
        except ValueError:
            logging.error('GPSError: bad date string {0}'.format(self.gps_string))
            return None
        
    @property
    def declination(self):
        """
        geomagnetic declination in degrees from north
        """
        if self._declination is None or self._declination_hemisphere is None:
            return None
        
        dec = float(self._declination)
        if 'w' in self._declination_hemisphere.lower():
            dec *= -1
        return dec

        
    @property
    def gps_type(self):
        """GPRMC or GPGGA"""
        return self._type
        
    @property
    def fix(self):
        """
        GPS fixed
        """
        if hasattr(self, '_fix'):
            return self._fix
        else:
            return None

class NIMSHeader(object):
    """
    class to hold the NIMS header information.  
    
    A typical header looks like
    
    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    >>>user field>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    SITE NAME: Budwieser Spring
    STATE/PROVINCE: CA
    COUNTRY: USA
    >>> The following code in double quotes is REQUIRED to start the NIMS <<
    >>> The next 3 lines contain values required for processing <<<<<<<<<<<<
    >>> The lines after that are optional <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    "300b"  <-- 2CHAR EXPERIMENT CODE + 3 CHAR SITE CODE + RUN LETTER
    1105-3; 1305-3  <-- SYSTEM BOX I.D.; MAG HEAD ID (if different)
    106  0 <-- N-S Ex WIRE LENGTH (m); HEADING (deg E mag N)
    109  90 <-- E-W Ey WIRE LENGTH (m); HEADING (deg E mag N)
    1         <-- N ELECTRODE ID
    3          <-- E ELECTRODE ID
    2          <-- S ELECTRODE ID
    4          <-- W ELECTRODE ID
    Cu          <-- GROUND ELECTRODE INFO
    GPS INFO: 01/10/19 16:16:42 1616.7000 3443.6088 115.7350 W 946.6
    OPERATOR: KP
    COMMENT: N/S CRS: .95/.96 DCV: 3.5 ACV:1
    E/W CRS: .85/.86 DCV: 1.5 ACV: 1
    Redeployed site for run b b/c possible animal disturbance
    
    """
    
    def __init__(self, fn=None):
        self.fn = fn
        self._max_header_length = 1000
        self.header_dict = None
        self.site_name = None
        self.state_province = None
        self.country = None
        self.box_id = None
        self.mag_id = None
        self.ex_length = None
        self.ex_azimuth = None
        self.ey_length = None
        self.ey_azimuth = None
        self.n_electrode_id = None
        self.s_electrode_id = None
        self.e_electrode_id = None
        self.w_electrode_id = None
        self.ground_electrode_info = None
        self.header_gps_stamp = None
        self.header_gps_latitude = None
        self.header_gps_longitude = None
        self.header_gps_elevation = None
        self.operator = None
        self.comments = None
        self.run_id = None
        self.data_start_seek = 0
        
        
    def read_header(self, fn=None):
        """
        read header information
        """
        if fn is not None:
            self.fn = fn
            
        if not os.path.exists(self.fn):
            raise NIMSError('Could not find file {0}'.format(self.fn))
            
        print('Reading NIMS file {0}'.format(self.fn))
        logging.info('='*72)
        logging.info('Reading NIMS file {0}'.format(self.fn))
            
        ### load in the entire file, its not too big
        with open(self.fn, 'rb') as fid:
            header_str = fid.read(self._max_header_length)
            header_list = header_str.split(b'\r')
            
        self.header_dict = {}
        last_index = len(header_list)
        last_line = header_list[-1]
        for ii, line in enumerate(header_list[0:-1]):
            if ii == last_index:
                break
            if b'comments' in line.lower():
                last_line = header_list[ii+1]
                last_index = ii + 1
                
            line = line.decode()
            if line.find('>') == 0:
                continue
            elif line.find(':') > 0:
                key, value = line.split(':', 1)
                self.header_dict[key.strip().lower()] = value.strip()
            elif line.find('<--') > 0:
                value, key = line.split('<--')
                self.header_dict[key.strip().lower()] = value.strip()
        ### sometimes there are some spaces before the data starts
        if last_line.count(b' ') > 0:
            if last_line[0:1] == b' ':
                last_line = last_line.strip()
            else:
                last_line = last_line.split()[1].strip()
        data_start_byte = last_line[0:1]
        ### sometimes there are rogue $ around
        if data_start_byte in [b'$', b'g']:
            data_start_byte = last_line[1:2]
        self.data_start_seek = header_str.find(data_start_byte)
        
        self.parse_header_dict()
    
    def parse_header_dict(self, header_dict=None):
        """
        parse the header dictionary into something useful
        """
        if header_dict is not None:
            self.header_dict = header_dict
            
        assert isinstance(self.header_dict, dict)
        
        for key, value in self.header_dict.items():
            if 'wire' in key:
                if key.find('n') == 0:
                    self.ex_length = float(value.split()[0])
                    self.ex_azimuth = float(value.split()[1])
                elif key.find('e') == 0:
                    self.ey_length = float(value.split()[0])
                    self.ey_azimuth = float(value.split()[1])
            elif 'system' in key:
                self.box_id = value.split(';')[0].strip()
                self.mag_id = value.split(';')[1].strip()
            elif 'gps' in key:
                gps_list = value.split()
                self.header_gps_stamp = dateutil.parser.parse(' '.join(gps_list[0:2]),
                                                              dayfirst=True)
                self.header_gps_latitude = self._get_latitude(gps_list[2], gps_list[3])
                self.header_gps_longitude = self._get_longitude(gps_list[4], gps_list[5])
                self.header_gps_elevation = float(gps_list[6])
            elif 'run' in key:
                self.run_id = value.replace('"', '')
            else:
                setattr(self, key.replace(' ', '_').replace('/','_'), value)
    
    def _get_latitude(self, latitude, hemisphere):
        if not isinstance(latitude, float):
            latitude = float(latitude)
        if hemisphere.lower() == 'n':
            return latitude
        if hemisphere.lower() == 's':
            return -1*latitude

    def _get_longitude(self, longitude, hemisphere):
        if not isinstance(longitude, float):
            longitude = float(longitude)
        if hemisphere.lower() == 'e':
            return longitude
        if hemisphere.lower() == 'w':
            return -1*longitude
        
class NIMS(NIMSHeader):
    """
    NIMS Class will read in a NIMS DATA.BIN file.
    
    A fast way to read the binary files are to first read in the GPS strings, 
    the third byte in each block as a character and parse that into valid 
    GPS stamps.
    
    Then read in the entire data set as unsigned 8 bit integers and reshape
    the data to be n seconds x block size.  Then parse that array into the 
    status information and data.
    
    I only have a limited amount of .BIN files to test so this will likely 
    break if there are issues such as data gaps.  
    
    .. todo:: deal with timing issues, right now a warning is sent to the user
              need to figure out a way to find where the gap is and adjust
              accordingly.
              
    .. warning:: 
        Currently Only 8 Hz data is supported 
    """
    
    def __init__(self, fn=None):

        super().__init__(fn)
        self.block_size = 131
        self.block_sequence = [1, self.block_size]
        self.sampling_rate = 8 ### samples/second
        self.e_conversion_factor = 2.44141221047903e-06
        self.h_conversion_factor = 0.01
        self.t_conversion_factor = 70
        self.t_offset = 18048
        self._int_max = 8388608
        self._int_factor = 16777216
        self._block_dict = {'soh':0,
                            'block_len':1,
                            'status':2,
                            'gps':3,
                            'sequence':4,
                            'elec_temp':(5, 6),
                            'box_temp':(7, 8),
                            'logic':81,
                            'end':130}
        self.info_array = None
        self.stamps = None
        self.ts = None
        self.gaps = None
        self.duplicate_list = None
        
        self.indices = self._make_index_values()
        
        if self.fn is not None:
            self.read_nims()
            
    @property
    def latitude(self):
        """
        median latitude value from all the GPS stamps in decimal degrees
        WGS84
        
        Only get from the GPRMC stamp as they should be duplicates
        """
        if self.stamps is not None:
            latitude = np.zeros(len(self.stamps))
            for ii, stamp in enumerate(self.stamps):
                latitude[ii] = stamp[1][0].latitude
            return np.median(latitude[np.nonzero(latitude)])
        else:
            return None
    
    @property
    def longitude(self):
        """
        median longitude value from all the GPS stamps in decimal degrees
        WGS84
        
        Only get from the first stamp within the sets
        """
        if self.stamps is not None:
            longitude = np.zeros(len(self.stamps))
            for ii, stamp in enumerate(self.stamps):
                longitude[ii] = stamp[1][0].longitude
            return np.median(longitude[np.nonzero(longitude)])
        else:
            return None
    
    @property
    def elevation(self):
        """
        median elevation value from all the GPS stamps in decimal degrees
        WGS84
        
        Only get from the first stamp within the sets
        """
        if self.stamps is not None:
            elevation = np.zeros(len(self.stamps))
            for ii, stamp in enumerate(self.stamps):
                if len(stamp[1]) == 1:
                    elev = stamp[1][0].elevation
                if len(stamp[1]) == 2:
                    elev = stamp[1][1].elevation
                if elev is None:
                    continue
                elevation[ii] = elev 
            return np.median(elevation[np.nonzero(elevation)])
        else:
            return None
    
    @property
    def start_time(self):
        """
        start time is the first good GPS time stamp minus the seconds to the
        beginning of the time series.
        """
        if self.stamps is not None:
            return self.ts.index[0]
        else:
            return None
        
    @property
    def end_time(self):
        """
        start time is the first good GPS time stamp minus the seconds to the
        beginning of the time series.
        """
        if self.stamps is not None:
            return self.ts.index[-1]
        else:
            return None
        
    @property
    def hx(self):
        """HX"""
        if self.ts is not None:
            ts_obj = ts.MTTS()
            ts_obj.fn = self.fn
            ts_obj.station = self.run_id
            ts_obj.lat = self.latitude
            ts_obj.lon = self.longitude
            ts_obj.elev = self.elevation
            ts_obj.azimuth = 0
            ts_obj.component = 'hx'
            ts_obj.data_logger = self.box_id
            ts_obj.instrument_id = self.mag_id
            ts_obj.channel_number = 1
            ts_obj.sampling_rate = self.sampling_rate
            ts_obj.ts = pd.DataFrame({'data':self.ts.hx})
            return ts_obj
        else:
            return None
    
    @property
    def hy(self):
        """HY"""
        if self.ts is not None:
            ts_obj = ts.MTTS()
            ts_obj.fn = self.fn
            ts_obj.station = self.run_id
            ts_obj.lat = self.latitude
            ts_obj.lon = self.longitude
            ts_obj.elev = self.elevation
            ts_obj.azimuth = 90
            ts_obj.component = 'hy'
            ts_obj.data_logger = self.box_id
            ts_obj.instrument_id = self.mag_id
            ts_obj.channel_number = 2
            ts_obj.sampling_rate = self.sampling_rate
            ts_obj.ts = pd.DataFrame({'data':self.ts.hy})
            return ts_obj
        else:
            return None
        
    @property
    def hz(self):
        """HZ"""
        if self.ts is not None:
            ts_obj = ts.MTTS()
            ts_obj.fn = self.fn
            ts_obj.station = self.run_id
            ts_obj.lat = self.latitude
            ts_obj.lon = self.longitude
            ts_obj.elev = self.elevation
            ts_obj.azimuth = 90
            ts_obj.component = 'hz'
            ts_obj.data_logger = self.box_id
            ts_obj.instrument_id = self.mag_id
            ts_obj.channel_number = 3
            ts_obj.sampling_rate = self.sampling_rate
            ts_obj.ts = pd.DataFrame({'data':self.ts.hz})
            return ts_obj
        else:
            return None
        
    @property
    def ex(self):
        """EX"""
        if self.ts is not None:
            ts_obj = ts.MTTS()
            ts_obj.fn = self.fn
            ts_obj.station = self.run_id
            ts_obj.lat = self.latitude
            ts_obj.lon = self.longitude
            ts_obj.elev = self.elevation
            ts_obj.azimuth = self.ex_azimuth
            ts_obj.dipole_length = self.ex_length
            ts_obj.component = 'ex'
            ts_obj.data_logger = self.box_id
            ts_obj.instrument_id = 1
            ts_obj.channel_number = 4
            ts_obj.sampling_rate = self.sampling_rate
            ts_obj.ts = pd.DataFrame({'data':self.ts.ex})
            return ts_obj
        else:
            return None
    @property
    def ey(self):
        """EY"""
        if self.ts is not None:
            ts_obj = ts.MTTS()
            ts_obj.fn = self.fn
            ts_obj.station = self.run_id
            ts_obj.lat = self.latitude
            ts_obj.lon = self.longitude
            ts_obj.elev = self.elevation
            ts_obj.azimuth = self.ey_azimuth
            ts_obj.dipole_length = self.ey_length
            ts_obj.component = 'ey'
            ts_obj.data_logger = self.box_id
            ts_obj.instrument_id = 1
            ts_obj.channel_number = 5
            ts_obj.sampling_rate = self.sampling_rate
            ts_obj.ts = pd.DataFrame({'data':self.ts.ey})
            return ts_obj
        else:
            return None        
        
    def _make_index_values(self):
        """
        Index values for the channels recorded
        """
        ### make an array of index values for magnetics and electrics
        indices = np.zeros((8,5), dtype=np.int)
        for kk in range(8):
            ### magnetic blocks
            for ii in range(3):
                indices[kk, ii] = 9 + (kk) * 9 + (ii) * 3
            ### electric blocks
            for ii in range(2):
                indices[kk, 3+ii] = 82 + (kk) * 6 + (ii) * 3
        return indices
                
    def _get_gps_string_list(self, nims_string):
        """
        get the gps strings from the raw string output by the NIMS.  This will
        take the 3rd value in each block, concatenate into a long string and
        then make a list by splitting by '$'.  The index values of where the
        '$' are found are also calculated.
        
        :param str nims_string: raw binary string output by NIMS
        
        :returns: list of index values associated with the location of the '$'
        
        :returns: list of possible raw GPS strings
        
        .. note:: This assumes that there are an even amount of data blocks.  
                  Might be a bad assumption          
        """
        ### get index values of $ and gps_strings
        index_values = []
        gps_str_list = []
        for ii in range(int(len(nims_string)/self.block_size)):
            index = ii*self.block_size+3
            g_char = struct.unpack('c', 
                                   nims_string[index:index+1])[0]
            if g_char == b'$':
                index_values.append((index-3)/self.block_size)
            gps_str_list.append(g_char)
        gps_raw_stamp_list = b''.join(gps_str_list).split(b'$')
        return index_values, gps_raw_stamp_list
    
    def get_stamps(self, nims_string):
        """
        get a list of valid GPS strings and match synchronous GPRMC with GPGGA
        stamps if possible.
        
        :param str nims_string: raw GPS string output by NIMS
        """
        ### read in GPS strings into a list to be parsed later
        index_list, gps_raw_stamp_list = self._get_gps_string_list(nims_string)
        
        gps_stamp_list = []
        ### not we are skipping the first entry, it tends to be not 
        ### complete anyway
        for index, raw_stamp in zip(index_list, gps_raw_stamp_list[1:]):
            gps_obj = GPS(raw_stamp, index)
            if gps_obj.valid:
                gps_stamp_list.append(gps_obj)
            
        return self._gps_match_gprmc_gpgga_strings(gps_stamp_list)
            
    def _gps_match_gprmc_gpgga_strings(self, gps_obj_list):
        """
        match GPRMC and GPGGA strings together into a list
        
        [[GPRMC, GPGGA], ...]
        
        :param list gps_obj_list: list of GPS objects
        
        :returns: list of matched GPRMC and GPGGA stamps 
        """
        ### match up the GPRMC and GPGGA together
        gps_match_list = []
        for ii in range(0, len(gps_obj_list)-1, 2):
            if gps_obj_list[ii] is None:
                continue
            time_stamp = gps_obj_list[ii].time_stamp
            match_list = [gps_obj_list[ii]]
            try:
                if gps_obj_list[ii+1].time_stamp.time() == time_stamp.time():
                    match_list.append(gps_obj_list[ii+1])
            except AttributeError:
                pass
            gps_match_list.append(match_list)
        return gps_match_list
        
        
    def _get_gps_stamp_indices_from_status(self, status_array):
        """
        get the index location of the stamps from the status array assuming 
        that 0 indicates GPS lock.
        
        :param array status_array: an array of status values from data blocks
        
        :returns: array of index values where GPS lock was acquired ignoring
                  sequential locks.   
        """
        
        index_values = np.where(status_array == 0)[0]
        status_index = np.zeros_like(index_values)
        for ii in range(index_values.size):
            if index_values[ii] - index_values[ii-1] == 1:
                continue
            else:
                status_index[ii] = index_values[ii]
        status_index = status_index[np.nonzero(status_index)]
        
        return status_index
    
    def match_staus_with_gps_stamps(self, status_array, gps_list):
        """
        Match the index values from the status array with the index values of 
        the GPS stamps.
        
        :param array status_array: array of status values from each data block
        :param list gps_list: list of valid GPS stamps [[GPRMC, GPGGA], ...]
        
        .. note:: I think there is a 2 second gap between the lock and the 
                  first stamp character.
        """
        
        stamp_indices = self._get_gps_stamp_indices_from_status(status_array)
        gps_stamps = []
        for index in stamp_indices:
            stamp_find = False
            for ii, stamps in enumerate(gps_list):
                index_diff = stamps[0].index - index
                ### check the index value, should be 2 or 74, if it is off by
                ### a value left or right apply a correction.
                if index_diff == 1 or index_diff == 73:
                    index += 1
                    stamps[0].index += 1
                elif index_diff == 2 or index_diff == 74:
                    index = index
                elif index_diff == 3 or index_diff == 75:
                    index -= 1
                    stamps[0].index -= 1
                if stamps[0].gps_type in ['GPRMC', 'gprmc']:
                    if index_diff in [1, 2, 3]:
                        gps_stamps.append((index, stamps))
                        stamp_find = True
                        del gps_list[ii]
                        break
                elif stamps[0].gps_type in ['GPGGA', 'gpgga']:
                    if index_diff in [73, 74, 75]:
                        gps_stamps.append((index, stamps))
                        stamp_find = True
                        del gps_list[ii]
                        break
            if not stamp_find:
                logging.warning('No good GPS stamp at {0} seconds'.format(index))

        return gps_stamps
    
    def find_sequence(self, data_array, block_sequence=None):
        """
        find a sequence in a given array
        
        :param array data_array: array of the data with shape [n, m]
                                 where n is the number of seconds recorded
                                 m is the block length for a given sampling
                                 rate.
        :param list block_sequence: sequence pattern to locate
                                    *default* is [1, 131] the start of a 
                                    data block.
                                    
        :returns: array of index locations where the sequence is found.
        """
        if block_sequence is not None:
            self.block_sequence = block_sequence
        
        n_data = data_array.size
        n_sequence = len(self.block_sequence)
        
        slices = [np.s_[ii:n_data-n_sequence+1+ii] for ii in range(n_sequence)]
        
        sequence_search =[data_array[slices[ii]] == self.block_sequence[ii]
                          for ii in range(n_sequence)][0]
        find_index = np.where(sequence_search == True)[0]
        
        return find_index
    
    def unwrap_sequence(self, sequence):
        """
        unwrap the sequence to be sequential numbers instead of modulated by
        256.  sets the first number to 0
        """
        count = 0
        unwrapped = np.zeros_like(sequence)
        for ii, seq in enumerate(sequence):
            unwrapped[ii] = seq + count * 256
            if seq == 255:
                count += 1
                
        unwrapped -= unwrapped[0]
        
        return unwrapped
    
    def _locate_duplicate_blocks(self, sequence):
        """
        locate the sequence number where the duplicates exist
        """
        
        duplicates = np.where(np.abs(np.diff(sequence)) == 0)[0]
        if len(duplicates) == 0:
            return None
        duplicate_list = []
        for dup in duplicates:
            dup_dict = {}
            dup_dict['sequence_index'] = dup
            dup_dict['ts_index_0'] = dup * self.sampling_rate
            dup_dict['ts_index_1'] = dup * self.sampling_rate + self.sampling_rate
            dup_dict['ts_index_2'] = (dup + 1) * self.sampling_rate
            dup_dict['ts_index_3'] = (dup + 1) * self.sampling_rate + self.sampling_rate
            duplicate_list.append(dup_dict)
        return duplicate_list
    
    def _check_duplicate_blocks(self, block_01, block_02, info_01, info_02):
        """
        make sure the blocks are truly duplicates
        """
        if np.array_equal(block_01, block_02):
            if np.array_equal(info_01, info_02):
                return True
            else:
                return False
        else:
            return False
    
    def remove_duplicates(self, info_array, data_array):
        """
        remove duplicate blocks, removing the first duplicate as suggested by
        Paul and Anna. Checks to make sure that the mag data are identical for 
        the duplicate blocks.  Removes the blocks from the information and
        data arrays and returns the reduced arrays.  This should sync up the
        timing of GPS stamps and index values.
        
        :param np.array info_array: structured array of block information
        :param np.array data_array: structured array of the data
        
        :returns: reduced information array
        :returns: reduced data array
        :returns: index of duplicates in raw data
        """
        ### locate 
        duplicate_test_list = self._locate_duplicate_blocks(self.info_array['sequence'])
        if duplicate_test_list is None:
            return info_array, data_array, None
        
        duplicate_list = []
        for d in duplicate_test_list:
            if self._check_duplicate_blocks(data_array[d['ts_index_0']:d['ts_index_1']],
                                            data_array[d['ts_index_2']:d['ts_index_3']],
                                            info_array[d['sequence_index']],
                                            info_array[d['sequence_index'] + 1]):
                duplicate_list.append(d)
        
        print('    Deleting {0} duplicate blocks'.format(len(duplicate_list)))
        ### get the index of the blocks to be removed, namely the 1st duplicate
        ### block
        remove_sequence_index = [d['sequence_index'] for d in duplicate_list]
        remove_data_index = np.array([np.arange(d['ts_index_0'], d['ts_index_1'], 1)
                                    for d in duplicate_list]).flatten() 
        ### remove the data
        return_info_array = np.delete(info_array, remove_sequence_index)
        return_data_array = np.delete(data_array, remove_data_index)
        
        ### set sequence to be monotonic
        return_info_array['sequence'][:] = np.arange(return_info_array.shape[0])
        
        return return_info_array, return_data_array, duplicate_list
        
    def read_nims(self, fn=None):
        """
        Read NIMS DATA.BIN file.
        
        1. Read in the header information and stores those as attributes
           with the same names as in the header file.
        
        2. Locate the beginning of the data blocks by looking for the 
           first [1, 131, ...] combo.  Anything before that is cut out.
        
        3. Make sure the data is a multiple of the block length, if the
           data is longer the extra bits are cut off.
        
        4. Read in the GPS data (3rd byte of each block) as characters.
           Parses those into valid GPS stamps with appropriate index locations
           of where the '$' was found.
          
        5. Read in the data as unsigned 8-bit integers and reshape the array
           into [N, data_block_length].  Parse this array into the status
           information and the data.
           
        6. Remove duplicate blocks, by removing the first of the duplicates
           as suggested by Anna and Paul.  
           
        7. Match the GPS locks from the status with valid GPS stamps.
                
        8. Check to make sure that there is the correct number of seconds
           between the first and last GPS stamp.  If there is not a warning
           message will appear. 
        
        .. note:: The data and information array returned have the duplicates
                  removed and the sequence reset to be monotonic.
        
        :param str fn: full path to DATA.BIN file
        
        """
        if fn is not None:
            self.fn = fn

        st = datetime.datetime.now()
        ### read in header information and get the location of end of header
        self.read_header(self.fn)
        
        ### load in the entire file, its not too big, start from the 
        ### end of the header information.
        with open(self.fn, 'rb') as fid:
            fid.seek(self.data_start_seek)
            data_str = fid.read()
            
        ### read in full string as unsigned integers
        data = np.frombuffer(data_str, dtype=np.uint8)
        
        ### need to make sure that the data starts with a full block
        find_first = self.find_sequence(data[0:self.block_size*5])[0]
        data = data[find_first:]
        
        ### get GPS stamps from the binary string first
        self.gps_list = self.get_stamps(data_str[find_first:])
        
        ### check the size of the data, should have an equal amount of blocks
        if (data.size % self.block_size) != 0:
            logging.warning('odd number of bytes {0}, not even blocks'.format(data.size)+\
                            'cutting down the data by {0}'.format(data.size % self.block_size))
            end_data = (data.size - (data.size % self.block_size))
            data = data[0:end_data]
            
        data = data.reshape((int(data.size/self.block_size), 
                             self.block_size))

        ### need to parse the data
        ### first get the status information
        self.info_array = np.zeros(data.shape[0],
                                   dtype=[('soh', np.int),
                                          ('block_len', np.int),
                                          ('status', np.int),
                                          ('gps', np.int),
                                          ('sequence', np.int),
                                          ('elec_temp', np.float),
                                          ('box_temp', np.float),
                                          ('logic', np.int),
                                          ('end', np.int)])    
        
        for key, index in self._block_dict.items():
            if 'temp' in key:
                value = ((data[:, index[0]] * 256 + data[:, index[1]]) - \
                         self.t_offset)/self.t_conversion_factor
            else:
                value = data[:, index]
            self.info_array[key][:] = value
            
        ### unwrap sequence
        self.info_array['sequence'] = self.unwrap_sequence(self.info_array['sequence'])
         
        ### get data
        data_array = np.zeros(data.shape[0]*self.sampling_rate,
                              dtype=[('hx', np.float),
                                     ('hy', np.float), 
                                     ('hz', np.float),
                                     ('ex', np.float),
                                     ('ey', np.float)])
        
        ### fill the data
        for cc, comp in enumerate(['hx', 'hy', 'hz', 'ex', 'ey']):
            channel_arr = np.zeros((data.shape[0], 8), dtype=np.float)
            for kk in range(self.sampling_rate):
                index = self.indices[kk, cc]
                value = (data[:, index]*256 + data[:, index+1]) * np.array([256]) + \
                        data[:, index+2]
                value[np.where(value > self._int_max)] -= self._int_factor
                channel_arr[:, kk] = value
            data_array[comp][:] = channel_arr.flatten()
            
        ### clean things up
        ### I guess that the E channels are opposite phase?
        for comp in ['ex', 'ey']:
            data_array[comp] *= -1
            
        ### remove duplicates 
        self.info_array, data_array, self.duplicate_list = self.remove_duplicates(self.info_array,
                                                                                  data_array)
        ### get GPS stamps with index values
        self.stamps = self.match_staus_with_gps_stamps(self.info_array['status'],
                                                       self.gps_list)
        ### align data 
        self.ts = self.align_data(data_array, self.stamps) 
        et = datetime.datetime.now()
        
        print('--> Took {0:.2f} seconds'.format((et-st).total_seconds()))

    def _get_first_gps_stamp(self, stamps):
        """
        get the first GPRMC stamp
        """ 
        for stamp in stamps:
            if stamp[1][0].gps_type in ['gprmc', 'GPRMC']:
                return stamp
        return None
    
    def _get_last_gps_stamp(self, stamps):
        """
        get the last gprmc stamp
        """
        for stamp in stamps[::-1]:
            if stamp[1][0].gps_type in ['gprmc', 'GPRMC']:
                return stamp
        return None
    
    def _locate_timing_gaps(self, stamps):
        """
        locate timing gaps in the data by comparing the stamp index with the 
        GPS time stamp.  The number of points and seconds should be the same
        
        :param list stamps: list of GPS stamps [[status_index, [GPRMC, GPGGA]]]
        
        :returns: list of gap index values
        """
        stamp_01 = self._get_first_gps_stamp(stamps)[1][0]
        diff_arr = np.zeros(len(stamps))
        diff_arr[0] = -666
        for ii, stamp in enumerate(stamps[1:], 1):
            stamp = stamp[1][0]
            if stamp._date == '010180':
                diff_arr[ii] = -666
                continue
            time_diff = (stamp.time_stamp - stamp_01.time_stamp).total_seconds()
            index_diff = stamp.index - stamp_01.index
            
            diff_arr[ii] = index_diff - time_diff
        
        gap_max = int(diff_arr.max())
        gap_beginning = []
        if gap_max > 0:
            print('    Check times:')
            for ii in range(1, gap_max+1, 1):
                step_index = np.where(diff_arr == ii)[0][0]
                gap_beginning.append(step_index)
                print('{0}{1} is off from start time by {2} seconds'.format(
                      ' '*4, 
                      stamps[step_index][1][0].time_stamp.isoformat(),
                      ii))

        return gap_beginning
           
    def check_timing(self, stamps):
        """
        make sure that there are the correct number of seconds in between
        the first and last GPS GPRMC stamps
        
        :param list stamps: list of GPS stamps [[status_index, [GPRMC, GPGGA]]]
        
        :returns: [ True | False ] if data is valid or not.
        :returns: gap index locations
        
        .. note:: There is currently no solution to fix the gap or to 
                  locate where the gap occurs.  Still trying to figure out if
                  there is an acctual data gap or there is something wrong
                  with location with in the file of the stamps.
        """
        gaps = None
        first_stamp = self._get_first_gps_stamp(stamps)[1][0]
        last_stamp = self._get_last_gps_stamp(stamps)[1][0]
        
        time_diff = last_stamp.time_stamp - first_stamp.time_stamp
        index_diff = last_stamp.index - first_stamp.index
        
        difference = index_diff - time_diff.total_seconds() 
        if difference != 0:
            
            gaps = self._locate_timing_gaps(stamps)
            if len(gaps) > 0:
                print('-'*50)
                print('Timing might be off by {0} seconds'.format(difference))
                print('-'*50)
            
            return False, gaps
        else:
            return True, gaps 
                   
    def align_data(self, data_array, stamps):
        """
        Need to match up the first good GPS stamp with the data
        
        Do this by using the first GPS stamp and assuming that the time from
        the first time stamp to the start is the index value.
        
        put the data into a pandas data frame that is indexed by time
        
        :param array data_array: structure array with columns for each 
                                 component [hx, hy, hz, ex, ey]
        :param list stamps: list of GPS stamps [[status_index, [GPRMC, GPGGA]]]
        
        :returns: pandas DataFrame with colums of components and indexed by 
                  time initialized by the start time.
        
        .. note:: There is currently no solution to fix the gap or to 
                  locate where the gap occurs.  Just a message of where the 
                  gap may occur.
        """
        ### check timing first to make sure there is no drift
        timing_valid, self.gaps = self.check_timing(stamps)
        
        ### first GPS stamp within the data is at a given index that is 
        ### assumed to be the number of seconds from the start of the run.
        ### therefore make the start time the first GPS stamp time minus
        ### the index value for that stamp.
        ### need to be sure that the first GPS stamp has a date, need GPRMC
        first_stamp = self._get_first_gps_stamp(stamps)
        first_index = first_stamp[0]
        start_time = first_stamp[1][0].time_stamp - \
                            datetime.timedelta(seconds=int(first_index))

        dt_index = self.make_dt_index(start_time.isoformat(),
                                      self.sampling_rate,
                                      n_samples=data_array.shape[0])
        
        return pd.DataFrame(data_array, index=dt_index)
        
    def calibrate_data(self, ts):
        """
        Apply calibrations to data
        
        .. note:: this needs work, would not use this now.
        """
        
        ts[['hx', 'hy', 'hz']] *= self.h_conversion_factor
        ts[['ex', 'ey']] *= self.e_conversion_factor
        ts['ex'] /= self.ex_length/1000.
        ts['ey'] /= self.ey_length/1000.
        
        return ts
        
    def make_dt_index(self, start_time, sampling_rate, stop_time=None,
                      n_samples=None):
        """
        make time index array

        .. note:: date-time format should be YYYY-M-DDThh:mm:ss.ms UTC

        :param start_time: start time
        :type start_time: string

        :param end_time: end time
        :type end_time: string

        :param sampling_rate: sampling_rate in samples/second
        :type sampling_rate: float
        """

        # set the index to be UTC time
        dt_freq = '{0:.0f}N'.format(1./(sampling_rate)*1E9)
        if stop_time is not None:
            dt_index = pd.date_range(start=start_time,
                                     end=stop_time,
                                     freq=dt_freq,
                                     closed='left',
                                     tz='UTC')
        elif n_samples is not None:
            dt_index = pd.date_range(start=start_time,
                                     periods=n_samples,
                                     freq=dt_freq,
                                     tz='UTC')
        else:
            raise ValueError('Need to input either stop_time or n_samples')

        return dt_index
    
    def plot_time_series(self, fig_num=1, order=['hx', 'hy', 'hz', 'ex', 'ey']):
        """
        plot time series
        """
        
        fig = plt.figure(fig_num)
        ax_list = []
        n = len(order)
        for ii, comp in enumerate(order, 1):
            if ii == 1:
                ax = fig.add_subplot(n, 1, ii)
            else:
                ax = fig.add_subplot(n, 1, ii, sharex=ax_list[0])
            l1, = ax.plot(getattr(self, comp).ts.data)
            ax_list.append(ax)
            ax.set_ylabel(comp.upper())
            
        return ax_list
                 
class Response(object):
    """
    class for instrument response functions.
    
    """
    
    def __init__(self, system_id=None, **kwargs):
        self.system_id = system_id
        self.hardware = 'PC'
        self.instrument_type = 'backbone'
        self.sampling_rate = 8
        self.e_conversion_factor = 2.44141221047903e-06
        self.h_conversion_factor = 0.01
        
        self.time_delays_dict = {'hp200':{'hx':-0.0055,
                                          'hy':-0.0145, 
                                          'hz':-0.0235,
                                          'ex':0.1525,
                                          'ey':0.0275},
                                 1: {'hx':-0.1920,
                                     'hy':-0.2010,
                                     'hz':-0.2100,
                                     'ex':-0.2850,
                                     'ey':-0.2850},
                                 8: {'hx':0.2455,
                                     'hy':0.2365,
                                     'hz':0.2275,
                                     'ex':0.1525,
                                     'ey':0.1525}}
        self.mag_low_pass = {'name': '3 pole butterworth',
                             'type': 'poles-zeros',
                             'parameters':{'zeros': [0, 3, 1984.31],
                                           'poles': [complex(-6.28319, 10.8825),
                                                     complex(-6.28319, 10.8825),
                                                     complex(-12.5664, 0)]}}
        self.electric_low_pass = {'name': '5 pole butterworth',
                                  'type': 'poles-zeros',
                                  'parameters':{'zeros':[0, 5, 313384],
                                                'poles':[complex(-3.88301,11.9519),
                                                         complex(-3.88301,-11.9519),
                                                         complex(-10.1662,7.38651),
                                                         complex(-10.1662,-7.38651),
                                                         complex(-12.5664,0.0)]}}
        self.electric_high_pass_pc = {'name': '1 pole butterworth',
                                      'type': 'poles-zeros',
                                      'parameters':{'zeros': [1, 1, 1],
                                                    'poles': [complex(0.0, 0.0),
                                                              complex(-3.333333E-05, 0.0)]},
                                      't0':2 * np.pi * 30000}
        self.electric_high_pass_hp = {'name': '1 pole butterworth',
                                      'type': 'poles-zeros',
                                      'parameters':{'zeros': [1, 1, 1],
                                                    'poles': [complex(0.0, 0.0),
                                                              complex(-1.66667E-04, 0.0)]},
                                      't0':2 * np.pi * 6000}
        
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def get_electric_high_pass(self, hardware='pc'):
        """
        get the electric high pass filter based on the hardware
        """
        
        self.hardware = hardware
        if 'pc' in hardware.lower():
            return self.electric_high_pass_pc
        elif 'hp' in hardware.lower():
            return self.electric_high_pass_hp
        else:
            raise ResponseError('Hardware value {0} not understood'.format(self.hardware))
            
    def _get_dt_filter(self, channel, sampling_rate):
        """
        get the DT filter based on channel ans sampling rate
        """
        dt_filter = {'type':'dt',
                     'name':'time_offset',
                     'parameters':{'offset':self.time_delays_dict[sampling_rate][channel]}}
        return dt_filter
    
    def _get_mag_filter(self, channel):
        """
        get mag filter, seems to be the same no matter what
        """
        filter_list = [self.mag_low_pass]
        filter_list.append(self._get_dt_filter(channel, self.sampling_rate))
        
        return_dict = {'channel_id':channel,
                       'gain':1,
                       'conversion_factor':self.h_conversion_factor,
                       'units':'nT',
                       'filters': filter_list}
        return return_dict
    
    def _get_electric_filter(self, channel):
        """
        Get electric filter
        """
        filter_list = []
        if self.instrument_type in ['backbone']:
            filter_list.append(self.get_electric_high_pass(self.hardware))
        filter_list.append(self.electric_low_pass)
        filter_list.append(self._get_dt_filter(channel, self.sampling_rate))
            
        
        return_dict = {'channel_id':channel,
                       'gain':1,
                       'conversion_factor':self.e_conversion_factor,
                       'units':'nT',
                       'filters': filter_list}
        return return_dict
    
    @property
    def hx_filter(self):
        """HX filter"""
        
        return self._get_mag_filter('hx')
    
    @property
    def hy_filter(self):
        """HY Filter"""
        return self._get_mag_filter('hy')
    
    @property
    def hz_filter(self):
        return self._get_mag_filter('hz')
    
    @property
    def ex_filter(self):
        return self._get_electric_filter('ex')
    
    @property
    def ey_filter(self):
        return self._get_electric_filter('ey')
    
    
        
    
