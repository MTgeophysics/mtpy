# -*- coding: utf-8 -*-
"""
====================
ZenTools
====================

    * Tools for reading and writing files for Zen and processing software
    * Tools for copying data from SD cards
    * Tools for copying schedules to SD cards
    
    
Created on Tue Jun 11 10:53:23 2013

@author: jpeacock-pr
"""

#==============================================================================

import numpy as np
import scipy.signal as sps
import time
import datetime
import os
import struct
import string
import win32api
import shutil
from collections import Counter
import mtpy.utils.filehandling as mtfh
import mtpy.processing.birrp as birrp
import mtpy.utils.configfile as mtcfg
import mtpy.utils.exceptions as mtex
import mtpy.utils.configfile as mtcf
import matplotlib.pyplot as plt
import mtpy.imaging.plotspectrogram as plotspectrogram
import mtpy.imaging.plotnresponses as plotnresponses
import mtpy.imaging.plotresponse as plotresponse
from cStringIO import StringIO
import sys
import mtpy.processing.filter as mtfilt

try:
    import mtpy.utils.mseed as mtmseed
except ImportError:
    print ('Can not convert data to mini seed format need to install Obspy, '
           'good luck! You can find information on Obspy at '
           'https://github.com/obspy/obspy/wiki') 
    
#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'
datetime_sec = '%Y-%m-%d %H:%M:%S'
#==============================================================================
# 
#==============================================================================

class Z3D_Header(object):
    """
    class for z3d header.  This will read in the header information of a 
    Z3D file and make each metadata entry an attirbute
    
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
                 
        **fid** : file object
                  ie. open(Z3D_file, 'rb')
                  
    ======================== ==================================================
    Attributes               Definition    
    ======================== ==================================================
    _header_len              lenght of header in bits (512)
    ad_gain                  gain of channel
    ad_rate                  sampling rate in Hz 
    alt                      altitude of the station (not reliable)
    attenchannelsmask        not sure
    box_number               ZEN box number
    box_serial               ZEN box serial number  
    channel                  channel number of the file 
    channelserial            serial number of the channel board
    duty                     duty cycle of the transmitter
    fpga_buildnum            build number of one of the boards
    gpsweek                  GPS week
    header_str               full header string
    lat                      latitude of station
    logterminal              not sure
    long                     longitude of the station
    main_hex_buildnum        build number of the ZEN box in hexidecimal
    numsats                  number of gps satelites
    period                   period of the transmitter 
    tx_duty                  transmitter duty cycle
    tx_freq                  transmitter frequency
    version                  version of the firmware
    ======================== ==================================================
    
    
    ======================== ==================================================
    Methods                  Description
    ======================== ==================================================
    convert_values           convert the read in header metadata to 
                             appropriate units and data types.
    read_header              read in the header data from the given file
    ======================== ==================================================
    
    Example
    --------------
        >>> import mtpy.usgs.zen as zen
        >>> z3d_fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3d_Header()
        >>> header_obj.read_header()
    
    """
    
    def __init__(self, fn=None, fid=None, **kwargs):
        self.fn = fn
        self.fid = fid
        
        self.header_str = None
        self._header_len = 512
        
        self.ad_gain = None
        self.ad_rate = None
        self.alt = None
        self.attenchannelsmask = None
        self.box_number = None
        self.box_serial = None
        self.channel = None
        self.channelserial = None
        self.duty = None
        self.fpga_buildnum = None
        self.gpsweek = 1740
        self.lat = None
        self.logterminal = None
        self.long = None
        self.main_hex_buildnum = None
        self.numsats = None
        self.period = None
        self.tx_duty = None
        self.tx_freq = None
        self.version = None 
        
        for key in kwargs:
            setattr(self, key, kwargs[key])
    
    def read_header(self, fn=None, fid=None):
        """
        read in the header string
        """
        if fn is not None:
            self.fn = fn
            
        if fid is not None:
            self.fid = fid

        if self.fn is None and self.fid is None:
            print 'no file to read'
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(0)
                self.header_str = self.fid.read(self._header_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = file(self.fn, 'rb')
                self.header_str = self.fid.read(self._header_len)
            else:
                self.fid.seek(0)
                self.header_str = self.fid.read(self._header_len)

        header_list = self.header_str.split('\n')
        for h_str in header_list:
            if h_str.find('=') > 0:
                h_list = h_str.split('=')
                h_key = h_list[0].strip().lower()
                h_key = h_key.replace(' ', '_').replace('/', '').replace('.', '_')
                h_value = self.convert_value(h_key, h_list[1].strip())
                setattr(self, h_key, h_value)

    def convert_value(self, key_string, value_string):
        """
        convert the value to the appropriate units given the key
        """
        
        try:
            return_value = float(value_string)
        except ValueError:
            return_value = value_string
        
        if key_string.lower() == 'lat' or key_string.lower() == 'long':
            return_value = np.rad2deg(float(value_string))
            
        return return_value
            
#==============================================================================
# meta data 
#==============================================================================
class Z3D_Schedule_metadata(object):
    """
    class object for metadata of Z3d file.  This will read in the schedule
    information of a Z3D file and make each metadata entry an attirbute.
    The attributes are left in capitalization of the Z3D file.
    
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
                 
        **fid** : file object
                  ie. open(Z3D_file, 'rb')
                  
    ======================== ==================================================
    Attributes               Definition    
    ======================== ==================================================
    AutoGain                 Auto gain for the channel  
    Comment                  Any comments for the schedule
    Date                     Date of when the schedule action was started
                             YYYY-MM-DD  
    Duty                     Duty cycle of the transmitter 
    FFTStacks                FFT stacks from the transmitter
    Filename                 Name of the file that the ZEN gives it
    Gain                     Gain of the channel
    Log                      Log the data [ Y | N ]
    NewFile                  Create a new file [ Y | N ]
    Period                   Period of the transmitter
    RadioOn                  Turn on the radio [ Y | N ]
    SR                       Sampling Rate in Hz
    SamplesPerAcq            Samples per aquisition for transmitter 
    Sleep                    Set the box to sleep [ Y | N ]
    Sync                     Sync with GPS [ Y | N ]
    Time                     Time the schedule action started
                             HH:MM:SS (GPS time)
    _header_len              length of header in bits (512)
    _schedule_metadata_len   length of schedule metadata in bits (512)
    fid                      file object of the file
    fn                       file name to read in
    meta_string              string of the schedule
    ======================== ==================================================
    
    
    ======================== ==================================================
    Methods                  Description
    ======================== ==================================================
    read_schedule_metadata   read in the schedule information from the given
                             file
    ======================== ==================================================
    
    Example
    --------------
        >>> import mtpy.usgs.zen as zen
        >>> z3d_fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3d_Schedule_metadata()
        >>> header_obj.read_schedule_metadata()
    """
    def __init__(self, fn=None, fid=None, **kwargs):
        self.fn = fn
        self.fid = None
        self.meta_string = None

        self._schedule_metadata_len = 512
        self._header_len = 512 
        
        self.AutoGain = None
        self.Comment = None
        self.Date = None
        self.Duty = None
        self.FFTStacks  = None
        self.Filename = None
        self.Gain = None
        self.Log = None
        self.NewFile = None
        self.Period = None
        self.RadioOn = None
        self.SR = None
        self.SamplesPerAcq = None
        self.Sleep = None
        self.Sync = None
        self.Time = None
        
        for key in kwargs:
            setattr(self, key, kwargs[key])

            
    def read_schedule_metadata(self, fn=None, fid=None):
        """
        read meta data string
        """
        if fn is not None:
            self.fn = fn
            
        if fid is not None:
            self.fid = fid

        if self.fn is None and self.fid is None:
            print 'no file to read'
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = file(self.fn, 'rb')
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)
            else:
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)
 
        meta_list = self.meta_string.split('\n')
        for m_str in meta_list:
            if m_str.find('=') > 0:
                m_list = m_str.split('=')
                m_key = m_list[0].split('.')[1].strip()
                m_key = m_key.replace('/', '')
                m_value = m_list[1].strip()
                setattr(self, m_key, m_value)

#==============================================================================
#  Meta data class    
#==============================================================================
class Z3D_Metadata(object):
    """
    class object for metadata of Z3d file.  This will read in the metadata
    information of a Z3D file and make each metadata entry an attirbute.
    The attributes are left in capitalization of the Z3D file.
    
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
                 
        **fid** : file object
                  ie. open(Z3D_file, 'rb')
                  
    ======================== ==================================================
    Attributes               Definition    
    ======================== ==================================================
    _header_length           length of header in bits (512)  
    _metadata_length         length of metadata blocks (512) 
    _schedule_metadata_len   length of schedule meta data (512)
    board_cal                board calibration np.ndarray()
    cal_ant                  antenna calibration
    cal_board                board calibration
    cal_ver                  calibration version
    ch_azimuth               channel azimuth
    ch_cmp                   channel component
    ch_length                channel length (or # of coil)
    ch_number                channel number on the ZEN board
    ch_xyz1                  channel xyz location (not sure)
    ch_xyz2                  channel xyz location (not sure)
    coil_cal                 coil calibration np.ndarray (freq, amp, phase)
    fid                      file object
    find_metadata            boolean of finding metadata
    fn                       full path to Z3D file
    gdp_operator             operater of the survey
    gdp_progver              program version
    job_by                   job preformed by  
    job_for                  job for 
    job_name                 job name
    job_number               job number
    m_tell                   location in the file where the last metadata 
                             block was found.
    rx_aspace                electrode spacing
    rx_sspace                not sure
    rx_xazimuth              x azimuth of electrode
    rx_xyz0                  not sure
    rx_yazimuth              y azimuth of electrode
    survey_type              type of survey 
    unit_length              length units (m)
    ======================== ==================================================
    
    
    ======================== ==================================================
    Methods                  Description
    ======================== ==================================================
    read_metadata            read in the metadata information from the given
                             file
    ======================== ==================================================
    
    Example
    --------------
        >>> import mtpy.usgs.zen as zen
        >>> z3d_fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3d_Metadata()
        >>> header_obj.read_metadata()
    
    
    """

    def __init__(self, fn=None, fid=None, **kwargs):
        self.fn = fn
        self.fid = None
        self.find_metadata = True
        self.board_cal = None
        self.coil_cal = None
        self._metadata_length = 512
        self._header_length = 512
        self._schedule_metadata_len = 512
        self.m_tell = 0
        
        self.cal_ant = None
        self.cal_board = None
        self.cal_ver = None
        self.ch_azimuth = None
        self.ch_cmp = None
        self.ch_length = None
        self.ch_number = None
        self.ch_xyz1 = None
        self.ch_xyz2 = None
        self.gdp_operator = None
        self.gdp_progver = None
        self.job_by = None
        self.job_for = None
        self.job_name = None
        self.job_number = None
        self.rx_aspace = None
        self.rx_sspace = None
        self.rx_xazimuth = None
        self.rx_xyz0 = None
        self.rx_yazimuth = None
        self.survey_type = None
        self.unit_length = None
        
        for key in kwargs:
            setattr(self, key, kwargs[key])

    def read_metadata(self, fn=None, fid=None):
        """
        read meta data
        """
        if fn is not None:
            self.fn = fn
            
        if fid is not None:
            self.fid = fid

        if self.fn is None and self.fid is None:
            print 'no file to read'
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(self._header_length+self._schedule_metadata_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = file(self.fn, 'rb')
                self.fid.seek(self._header_length+self._schedule_metadata_len)
            else:
                self.fid.seek(self._header_length+self._schedule_metadata_len)
        
        # read in calibration and meta data
        self.find_metadata = True
        self.board_cal = []
        self.coil_cal = []
        self.count = 0
        while self.find_metadata == True:
            test_str = self.fid.read(self._metadata_length)
            if test_str.lower().find('metadata record') > 0:
                self.count += 1
                cal_find = False
                test_str = test_str.strip().split('\n')[1] 
                if test_str.count('|') > 1:
                    for t_str in test_str.split('|'):
                        if t_str.find('=') == -1:
                            pass
                        else:
                            t_list = t_str.split('=')
                            t_key = t_list[0].strip().replace('.', '_')
                            t_value = t_list[1].strip()
                            setattr(self, t_key.lower(), t_value)
                elif test_str.lower().find('cal.brd') >= 0:
                    t_list = test_str.split(',')
                    t_key = t_list[0].strip().replace('.', '_')
                    setattr(self, t_key.lower(), t_list[1])
                    for t_str in t_list[2:]:
                        t_str = t_str.replace('\x00', '').replace('|', '')
                        self.board_cal.append([float(tt.strip()) 
                                           for tt in t_str.strip().split(':')])
                # some times the coil calibration does not start on its own line
                # so need to parse the line up and I'm not sure what the calibration
                # version is for so I have named it odd
                elif test_str.lower().find('cal.ant') >= 0:
                    # check to see if the coil calibration exists  
                    cal_find = True
                    if test_str.find('|') > 0:
                        odd_str = test_str.split('|')[0]
                        odd_list = odd_str.split(',')
                        odd_key = odd_list[0].strip().replace('.', '_')
                        setattr(self, odd_key.lower(), odd_list[1].strip())
                        
                        #this may be for a specific case so should test this
                        test_str = test_str.split('|')[1]
                        test_list = test_str.split(',')
                        if test_list[0].lower().find('cal.ant') >= 0:
                            m_list = test_list[0].split('=')
                            m_key = m_list[0].strip().replace('.', '_')
                            setattr(self, m_key.lower(), m_list[1].strip())
                        else:
                            for t_str in test_list[1:]:
                                self.coil_cal.append([float(tt.strip()) 
                                                 for tt in t_str.split(':')])
                elif cal_find:
                    t_list = test_str.split(',')
                    for t_str in t_list:
                        if t_str.find('\x00') >= 0:
                            pass
                        else:
                            self.coil_cal.append([float(tt.strip()) 
                                                for tt in t_str.strip().split(':')])
            else:
                self.find_metadata = False
                # need to go back to where the meta data was found wo
                # we don't skip a gps time stamp
                self.m_tell = self.fid.tell()-self._metadata_length
                    
        # make coil calibration and board calibration structured arrays
        if len(self.coil_cal) > 0:
            self.coil_cal = np.core.records.fromrecords(self.coil_cal, 
                                           names='frequency, amplitude, phase',
                                           formats='f8, f8, f8')
        if len(self.board_cal) > 0:  
            self.board_cal = np.core.records.fromrecords(self.board_cal, 
                                   names='frequency, rate, amplitude, phase',
                                   formats='f8, f8, f8, f8')
        
        
#==============================================================================
# 
#==============================================================================

class Zen3D(object):
    """
    Deals with the raw Z3D files output by zen.
    
    Arguments
    -----------
        **fn** : string
                 full path to .Z3D file to be read in
                 
    ======================== ================================ =================
    Attributes               Description                      Default Value
    ======================== ================================ =================
    _block_len               length of data block to read in  65536
                             as chunks faster reading
    _counts_to_mv_conversion conversion factor to convert     9.53674316406e-10
                             counts to mv                   
    _gps_bytes               number of bytes for a gps stamp  16
    _gps_dtype               data type for a gps stamp        see below
    _gps_epoch               starting date of GPS time
                             format is a tuple                (1980, 1, 6, 0, 
                                                               0, 0, -1, -1, 0)
 
    _gps_f0                  first gps flag in raw binary     
    _gps_f1                  second gps flag in raw binary     
    _gps_flag_0              first gps flag as an int32       2147483647       
    _gps_flag_1              second gps flag as an int32      -2147483648  
    _gps_stamp_length        bit length of gps stamp          64 
    _leap_seconds            leap seconds, difference         16
                             between UTC time and GPS 
                             time.  GPS time is ahead 
                             by this much
    _week_len                week length in seconds           604800
    df                       sampling rate of the data        256
    fn                       Z3D file name                    None
    gps_flag                 full gps flag                    _gps_f0+_gps_f1 
    gps_stamps               np.ndarray of gps stamps         None
    header                   Z3D_Header object                Z3D_Header                
    metadata                 Z3D_Metadata                     Z3D_Metadata 
    schedule                 Z3D_Schedule_metadata            Z3D_Schedule
    time_series              np.ndarra(len_data)              None
    units                    units in which the data is in    counts           
    zen_schedule             time when zen was set to         None
                             run
    ======================== ================================ =================
    
    * gps_dtype is formated as np.dtype([('flag0', np.int32),
                                        ('flag1', np.int32),
                                        ('time', np.int32),
                                        ('lat', np.float64),
                                        ('lon', np.float64),
                                        ('num_sat', np.int32),
                                        ('gps_sens', np.int32),
                                        ('temperature', np.float32),
                                        ('voltage', np.float32),
                                        ('num_fpga', np.int32),
                                        ('num_adc', np.int32),
                                        ('pps_count', np.int32),
                                        ('dac_tune', np.int32),
                                        ('block_len', np.int32)])    
    
    ============================ ==============================================
    Methods                       Description
    ============================ ==============================================
    apply_addaptive_notch_filter apply a notch filter to the data, usually
                                 to remove 60 Hz noise and harmonics 
        
    get_gps_time                 converts the gps counts to relative epoch 
                                 seconds according to gps week.      
    get_UTC_date_time            converts gps seconds into the actual date and
                                 time in UTC.  Note this is different than GPS
                                 time which is how the zen is scheduled, so 
                                 the time will be off by the current amount of
                                 leap seconds.
    plot_timeseries              make a generic plot of the time series
    plot_spectra                 plot a the spectra in loglog scales.
    plot_spectrogram             plot the spectragram of the data.
    read_z3d                      read 3D file making sure all the time stamps
                                 are correctly spaced.  Returned time series 
                                 starts at  the first stamp which has the 
                                 correct amount of data points between it and 
                                 the next time stamp.  Note there are usually
                                 a few seconds at the end and maybe beginning 
                                 that aren't correct because the internal 
                                 computer is busy switchin sampling rate.
    read_header                  read just the header data from the Z3D file
    read_metadata                read just the metadata from the Z3D file
    read_schedule                read just the schedule info from the Z3D file
    validate_gps_time            make sure each time stamp is 1 second apart
    validate_time_blocks         make sure that the size of each time block
                                 between stamps is equal to the sampling rate
    write_ascii_mt_file          write an mtpy ascii file of the data
    ============================ ==============================================
      
      
    Example
    ----------------
        >>> import mtpy.usgs.zen as zen
        >>> zt = zen.Zen3D(r"/home/mt/mt00/mt00_20150522_080000_256_EX.Z3D")
        >>> zt.read_z3d()
        >>> ------- Reading /home/mt/mt00/mt00_20150522_080000_256_EX.Z3D -----
            --> Reading data took: 0.322 seconds
            Scheduled time was 2015-05-22,08:00:16 (GPS time)
            1st good stamp was 2015-05-22,08:00:18 (GPS time)
            difference of 2.00 seconds
            found 6418 GPS time stamps
            found 1642752 data points
        >>> zt.plot_time_series()

    """

    def __init__(self, fn=None, **kwargs):
        self.fn = fn

        self.header = Z3D_Header(fn)
        self.schedule = Z3D_Schedule_metadata(fn)
        self.metadata = Z3D_Metadata(fn)
        
        self._gps_stamp_length = kwargs.pop('stamp_len', 64)
        self._gps_bytes = self._gps_stamp_length/4
        
        self.gps_stamps = None
        
        self._gps_flag_0 = np.int32(2147483647)
        self._gps_flag_1 = np.int32(-2147483648)
        self._gps_f0 = self._gps_flag_0.tostring()
        self._gps_f1 = self._gps_flag_1.tostring()
        self.gps_flag = self._gps_f0+self._gps_f1

        self._gps_dtype = np.dtype([('flag0', np.int32),
                                     ('flag1', np.int32),
                                     ('time', np.int32),
                                     ('lat', np.float64),
                                     ('lon', np.float64),
                                     ('num_sat', np.int32),
                                     ('gps_sens', np.int32),
                                     ('temperature', np.float32),
                                     ('voltage', np.float32),
                                     ('num_fpga', np.int32),
                                     ('num_adc', np.int32),
                                     ('pps_count', np.int32),
                                     ('dac_tune', np.int32),
                                     ('block_len', np.int32)])
                                     
        self._week_len = 604800
        self._gps_epoch = (1980, 1, 6, 0, 0, 0, -1, -1, 0)
        self._leap_seconds = 16
        self._block_len = 2**16
        self.zen_schedule = None
        self._counts_to_mv_conversion = 9.5367431640625e-10
        self.units = 'counts'
        self.df = None
        
        self.time_series = None
        
    #====================================== 
    def read_header(self, fn=None, fid=None):
        """
        read header information from Z3D file
        
        Arguments
        ---------------
            **fn** : string
                     full path to Z3D file to read
                     
            **fid** : file object
                      if the file is open give the file id object
                      
        Outputs:
        ----------
            * fills the Zen3ZD.header object's attributes
            
        Example with just a file name
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_header(fn)
            
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_fid = open(fn, 'rb')
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_header(fid=z3d_fid)
            
        """
        
        if fn is not None:
            self.fn = fn
            
        self.header.read_header(fn=self.fn, fid=fid)
        self.df = self.header.ad_rate
        
    #====================================== 
    def read_schedule(self, fn=None, fid=None):
        """
        read schedule information from Z3D file
        
        Arguments
        ---------------
            **fn** : string
                     full path to Z3D file to read
                     
            **fid** : file object
                      if the file is open give the file id object
                      
        Outputs:
        ----------
            * fills the Zen3ZD.schedule object's attributes
            
        Example with just a file name
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_schedule(fn)
            
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_fid = open(fn, 'rb')
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_schedule(fid=z3d_fid)
        """
        
        if fn is not None:
            self.fn = fn
            
        self.schedule.read_schedule_metadata(fn=self.fn, fid=fid)
        # set the zen schedule time
        self.zen_schedule = '{0},{1}'.format(self.schedule.Date, 
                                             self.schedule.Time)
    
    #======================================     
    def read_metadata(self, fn=None, fid=None):
        """
        read header information from Z3D file
        
        Arguments
        ---------------
            **fn** : string
                     full path to Z3D file to read
                     
            **fid** : file object
                      if the file is open give the file id object
                      
        Outputs:
        ----------
            * fills the Zen3ZD.metadata object's attributes
            
        Example with just a file name
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_metadata(fn)
            
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_fid = open(fn, 'rb')
            >>> z3d_obj = zen.Zen3D()
            >>> z3d_obj.read_metadata(fid=z3d_fid)
        """
        
        if fn is not None:
            self.fn = fn
            
        self.metadata.read_metadata(fn=self.fn, fid=fid)
    
    #======================================    
    def read_z3d(self):
        """
        read in z3d file and populate attributes accordingly
        
        read in the entire file as if everything but header and metadata are
        np.int32, then extract the gps stamps and convert accordingly
        
        Checks to make sure gps time stamps are 1 second apart and incrementing
        as well as checking the number of data points between stamps is the
        same as the sampling rate.  
        
        Converts gps_stamps['time'] to seconds relative to header.gps_week 
        
        We skip the first two gps stamps because there is something wrong with 
        the data there due to some type of buffering.  
        
        Therefore the first GPS time is when the time series starts, so you
        will notice that gps_stamps[0]['block_len'] = 0, this is because there
        is nothing previous to this time stamp and so the 'block_len' measures
        backwards from the corresponding time index.
        
        
        """
        
        print '------- Reading {0} ---------'.format(self.fn)
        st = time.time()
        
        #get the file size to get an estimate of how many data points there are
        file_size = os.path.getsize(self.fn)
        
        # using the with statement works in Python versions 2.7 or higher
        # the added benefit of the with statement is that it will close the
        # file object upon reading completion.
        with open(self.fn, 'rb') as file_id:
        
            self.read_header(fid=file_id)
            self.read_schedule(fid=file_id)
            self.read_metadata(fid=file_id)
            
            # move the read value to where the end of the metadata is
            file_id.seek(self.metadata.m_tell)
            
            # initalize a data array filled with zeros, everything goes into
            # this array then we parse later
            data = np.zeros((file_size-512*(2+self.metadata.count))/4, 
                             dtype=np.int32)
            # go over a while loop until the data cound exceed the file size
            data_count = 0
            while data_count+self.metadata.m_tell/4 < data.size:
                test_str = np.fromstring(file_id.read(self._block_len), 
                                         dtype=np.int32)
                data[data_count:data_count+len(test_str)] = test_str
                data_count += test_str.size

        # find the gps stamps
        gps_stamp_find = np.where(data==self._gps_flag_0)[0]
        
        # skip the first two stamps and trim data
        data = data[gps_stamp_find[3]:]
        gps_stamp_find = np.where(data==self._gps_flag_0)[0]
        
        self.gps_stamps = np.zeros(len(gps_stamp_find), dtype=self._gps_dtype)
        
        for ii, gps_find in enumerate(gps_stamp_find):
            if data[gps_find+1] == self._gps_flag_1:
                gps_str = struct.pack('<'+'i'*self._gps_bytes,
                                      *data[gps_find:gps_find+self._gps_bytes])
                self.gps_stamps[ii] = np.fromstring(gps_str, 
                                                   dtype=self._gps_dtype)
                if ii > 0:
                    self.gps_stamps[ii]['block_len'] = gps_find-\
                                           gps_stamp_find[ii-1]-self._gps_bytes 
                elif ii == 0:
                    self.gps_stamps[ii]['block_len'] = 0
                data[gps_find:gps_find+self._gps_bytes] = 0

        # trim the data after taking out the gps stamps
        self.time_series = data[np.nonzero(data)]
        
        # time it
        et = time.time()
        print '--> Reading data took: {0:.3f} seconds'.format(et-st)
        
        self.validate_time_blocks()
        self.convert_gps_time()
        self.check_start_time()
        
        print '    found {0} GPS time stamps'.format(self.gps_stamps.shape[0])
        print '    found {0} data points'.format(self.time_series.size)
    
    #=======================================    
    def read_z3d_slow(self):
        """
        read Z3D file out put by Zen, this a slow method but if you want to 
        be sure the data is read in correctly, this method is the most 
        correct way.  It will be deprecated as soon as I field test the 
        read_z3d method
        
        """        
        print '------- Reading {0} ---------'.format(self.fn)
        st = time.time()
        
        file_id = file(self.fn, 'rb')
        
        self.read_header(fid=file_id)
        self.read_schedule(fid=file_id)
        self.read_metadata(fid=file_id)
        
        file_id.seek(self.metadata.m_tell)
        f_str = file_id.read()
        
        
       # find number of gps stamps from the data string
        num_gps_stamps = int(f_str.count(self.gps_flag))
        #num_data = int(num_gps_stamps*header.ad_rate)
        #df = int(header.ad_rate)
        
        # make and empty array for gps time stamps in the appropriate format
        self.gps_stamps = np.zeros(num_gps_stamps, dtype=self._gps_dtype)
        
        #--> find first time stamp
        find_0 = f_str.find(self.gps_flag)
        
        gps_stamp_str = f_str[find_0:find_0+self._gps_stamp_length]
        gps_stamp = np.fromstring(gps_stamp_str, dtype=self._gps_dtype)
        gps_stamp['block_len'] = 0
        self.gps_stamps[0] = gps_stamp
        
        # make the input string start from the end of the first GPS stamp
        f_str = f_str[find_0+self._gps_stamp_length:]
        
        # make an empty list to append each time block, this is relatively quick
        # though might get caught up in large data files.  But this is the safest
        # to avoid misplacing the time series data into a predefined array
        # aslo insures sequential time series data
        ts_string_list = []
        
        # this loop starts from the end of the first time stamp and searches for the
        # next time stamp.  It will convert the time stamp and append the time series
        # bintary string between the time stamps to ts_string_list to insure
        # the time series is sequential.
        # boolean for testing whether a gps time stamp is found or not
        stamp_count = 1
        while stamp_count < num_gps_stamps:
            # get the index of the next gps stamp
            stamp_index = f_str.find(self.gps_flag)
        
            # get the gps string and convert it according to the format
            gps_stamp_str = f_str[stamp_index:stamp_index+self._gps_stamp_length]
            gps_stamp = np.fromstring(gps_stamp_str, dtype=self._gps_dtype)
            
            # get length between time stamp and put it in the empty slot of the 
            # gps time stamp
            gps_stamp['block_len'] = len(f_str[0:stamp_index])/4
            self.gps_stamps[stamp_count] = gps_stamp
            
            # append the time series binary string to the list, this is faster 
            # than string concatenation
            ts_string_list.append(f_str[0:stamp_index])
        
            # remove the found time stamp and time series data and start the 
            # data binary string from the end of the found time stamp.
            f_str = f_str[stamp_index+self._gps_stamp_length:]
            stamp_count += 1
        
        # convert the time series data into int32 format
        self.time_series = np.fromstring(''.join(ts_string_list), 
                                         dtype=np.int32)
        
        # time it
        et = time.time()
        print '--> Reading data took: {0:.3f} seconds'.format(et-st)
        
        self.trim_data()
        self.validate_time_blocks()
        self.convert_gps_time()
        self.check_start_time()
        
        print '    found {0} GPS time stamps'.format(self.gps_stamps.shape[0])
        print '    found {0} data points'.format(self.time_series.size)
        
    #=================================================
    def trim_data(self):
        """
        apparently need to skip the first 3 seconds of data because of 
        something to do with the SD buffer
        
        This method will be deprecated after field testing
        """ 
        # the block length is the number of data points before the time stamp
        # therefore the first block length is 0.  The indexing in python 
        # goes to the last index - 1 so we need to put in 3
        ts_skip = self.gps_stamps['block_len'][0:3].sum()
        self.gps_stamps = self.gps_stamps[2:]
        self.gps_stamps[0]['block_len'] = 0
        self.time_series = self.time_series[ts_skip:]
        
    #=================================================
    def check_start_time(self):
        """
        check to make sure the scheduled start time is similar to 
        the first good gps stamp
        """
        
        # make sure the time is in gps time
        zen_start = self.get_UTC_date_time(self.header.gpsweek,
                                           self.gps_stamps['time'][0]+\
                                                            self._leap_seconds)
        # set the zen schedule to the first gps stamp
        self.zen_schedule = zen_start
        zen_time = time.strptime(zen_start, datetime_fmt)
        
        # calculate the scheduled start time
        s_start = '{0},{1}'.format(self.schedule.Date, self.schedule.Time)
        schedule_time = time.strptime(s_start, datetime_fmt)
        
        # reset the data and time in the schedule meta data so there is no
        # confusion on when the time series starts
        self.schedule.Date = zen_start.split(',')[0]
        self.schedule.Time = zen_start.split(',')[1]
        
        # estimate the time difference between the two                                               
        time_diff = time.mktime(zen_time)-time.mktime(schedule_time)
        print '    Scheduled time was {0} (GPS time)'.format(s_start)
        print '    1st good stamp was {0} (GPS time)'.format(zen_start)
        print '    difference of {0:.2f} seconds'.format(time_diff)
        
    #==================================================
    def validate_gps_time(self):
        """
        make sure each time stamp is 1 second apart
        """
        
        t_diff = np.zeros_like(self.gps_stamps['time'])
        
        for ii in range(len(t_diff)-1):
            t_diff[ii] = self.gps_stamps['time'][ii]-self.gps_stamps['time'][ii+1]
        
        bad_times = np.where(abs(t_diff) > 0.5)[0]
        if len(bad_times) > 0:
            print '-'*50
            for bb in bad_times:
                print 'bad time at index {0} > 0.5 s'.format(bb) 
    
    #================================================== 
    def validate_time_blocks(self):
        """
        validate gps time stamps and make sure each block is the proper length
        """
        # first check if the gps stamp blocks are of the correct length
        bad_blocks = np.where(self.gps_stamps['block_len'][1:] != 
                                self.header.ad_rate)[0]
                                
        if len(bad_blocks) > 0:
            if bad_blocks.max() < 5:
                ts_skip = self.gps_stamps['block_len'][0:bad_blocks[-1]+1].sum()
                self.gps_stamps = self.gps_stamps[bad_blocks[-1]:]
                self.time_series = self.time_series[ts_skip:]
                
                print '{0}Skipped the first {1} seconds'.format(' '*4,
                                                                bad_blocks[-1])
                print '{0}Skipped first {1} poins in time series'.format(' '*4,
                                                                      ts_skip)
            
    #================================================== 
    def convert_gps_time(self):
        """
        convert gps time integer to relative seconds from gps_week
        """
        # need to convert gps_time to type float from int
        dt = self._gps_dtype.descr
        dt[2] = ('time', np.float32)
        self.gps_stamps = self.gps_stamps.astype(np.dtype(dt))
        
        # convert to seconds
        # these are seconds relative to the gps week
        time_conv = self.gps_stamps['time'].copy()/1024.
        time_ms = (time_conv-np.floor(time_conv))*1.024
        time_conv = np.floor(time_conv)+time_ms
        
        self.gps_stamps['time'][:] = time_conv 
            
    #==================================================    
    def convert_counts(self):
        """
        convert the time series from counts to millivolts

        """
        
        return self.time_series*self._counts_to_mv_conversion
    
    #==================================================    
    def convert_mV(self):
        """
        convert millivolts to counts assuming no other scaling has been applied
        
        """
        
        return self.time_series/self._counts_to_mv_conversion
        
    #==================================================        
    def get_gps_time(self, gps_int, gps_week=0):
        """
        from the gps integer get the time in seconds.
        
        Arguments
        -------------
            **gps_int**: int
                         integer from the gps time stamp line
                         
            **gps_week**: int
                          relative gps week, if the number of seconds is 
                          larger than a week then a week is subtracted from 
                          the seconds and computed from gps_week += 1
                          
        Returns
        ---------
            **gps_time**: int
                          number of seconds from the beginning of the relative
                          gps week.
        
        """
            
        gps_seconds = gps_int/1024.
        
        gps_ms = (gps_seconds-np.floor(gps_int/1024.))*(1.024)
        
        cc = 0
        if gps_seconds > self._week_len:
            gps_week += 1
            cc = gps_week*self._week_len
            gps_seconds -= self._week_len
        
        gps_time = np.floor(gps_seconds)+gps_ms+cc
        
        return gps_time, gps_week
    
    #==================================================
    def get_UTC_date_time(self, gps_week, gps_time):
        """
        get the actual date and time of measurement as UTC. 
        
        Note that GPS time is curently off by 16 seconds from actual UTC time.
        
        Arguments
        -------------
            **gps_week**: int
                          integer value of gps_week that the data was collected
            
            **gps_time**: int
                          number of seconds from beginning of gps_week
            
            **leap_seconds**: int
                              number of seconds gps time is off from UTC time.
                              It is currently off by 16 seconds.
                              
        Returns
        ------------
            **date_time**: YYYY-MM-DD,HH:MM:SS
                           formated date and time from gps seconds.
        
        
        """
        # need to check to see if the time in seconds is more than a gps week
        # if it is add 1 to the gps week and reduce the gps time by a week
        if gps_time > self._week_len:
            gps_week += 1
            gps_time -= self._week_len
            
        mseconds = gps_time % 1
        
        #make epoch in seconds, mktime computes local time, need to subtract
        #time zone to get UTC
        epoch_seconds = time.mktime(self._gps_epoch)-time.timezone
        
        #gps time is 14 seconds ahead of GTC time, but I think that the zen
        #receiver accounts for that so we will leave leap seconds to be 0        
        gps_seconds = epoch_seconds+(gps_week*self._week_len)+gps_time-\
                                                        self._leap_seconds

        #compute date and time from seconds
        (year, month, day, hour, minutes, seconds, dow, jday, dls) = \
                                                    time.gmtime(gps_seconds)
        
        date_time = time.strftime(datetime_fmt ,(year,
                                                 month, 
                                                 day, 
                                                 hour, 
                                                 minutes, 
                                                 int(seconds+mseconds), 
                                                 0, 0, 0))
        return date_time
        
    #==================================================    
    def apply_adaptive_notch_filter(self, notch_dict):
        """
        apply notch filter to the data that finds the peak around each 
        frequency.
        
        see mtpy.processing.filter.adaptive_notch_filter
        
        Arguments
        -------------
            **notch_dict** : dictionary
                             dictionary of filter parameters.
                             if an empty dictionary is input the filter looks
                             for 60 Hz and harmonics to filter out.
        
        
        """
        
        try:
            self.time_series
        except AttributeError:
            self.read_3d()
        
        notches = notch_dict.pop('notches', list(np.arange(60, 2048, 60)))
        notchradius = notch_dict.pop('notchradius', 0.5)
        freqrad = notch_dict.pop('freqrad', 0.5)
        rp = notch_dict.pop('rp', 0.1)
        kwargs = {'df':self.df, 'notches':notches, 'notchradius':notchradius,
                  'freqrad':freqrad, 'rp':rp}
                  
        self.time_series, self.filt_list = \
                    mtfilt.adaptive_notch_filter(self.time_series, **kwargs) 
        
    #==================================================
    def write_ascii_mt_file(self, save_fn=None, save_station='mb', fmt='%.8e',
                            ex=100., ey=100., notch_dict=None):
        """
        write an mtpy time series data file
        
        Arguments
        -------------
            **save_fn** : full path to save file, if None file is saved as:
                          station_YYYYMMDD_hhmmss_df.component
                          
                          ex. mt01_20130206_120000_256.HX
                          
            **save_station** : string
                               prefix string to add to station number as only
                               integers can be input into metadata of the zen
                               boxes.  ex. mb001
            
            **fmt** : string format
                      format of data numbers output to ascii file.
                      *default* is '%.8e' for 8 significan figures in 
                      scientific notation.
            **ex** : float
                     scaling parameter of ex line, the length of the dipole
                     be careful to not scale when creating an .edi file
                     *default* is 1
                     
            **ey** : float
                     scaling parameter of ey line, the length of the dipole
                     be careful to not scale when creating an .edi file
                     *default* is 1
            **notch_dict** : dictionary
                             dictionary of notch filter parameters
                             *default* is None
                             if an empty dictionary is input then the 
                             filter looks for 60 Hz and harmonics to filter
                      
        Output
        -------------
            **fn_mt_ascii** : full path to saved file
            
        Example
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> z3d_obj = zen.Zen3D(fn)
            >>> asc_fn = z3d.write_ascii_mt_file(save_station='mt', 
                                                 notch_dict={})
        
        """
        if self.time_series is None:
            self.read_3d()
            
        time_series = self.convert_counts()
        if save_fn is None:
            svfn_directory = os.path.join(os.path.dirname(self.fn), 'TS')
            if not os.path.exists(svfn_directory):
                os.mkdir(svfn_directory)
                
            svfn_date = ''.join(self.schedule.Date.split('-'))
            svfn_time = ''.join(self.schedule.Time.split(':'))
            svfn_station = save_station+self.metadata.rx_xyz0.split(':')[0]
            save_fn = os.path.join(svfn_directory, 
                                   '{0}_{1}_{2}_{3}.{4}'.format(svfn_station,
                                                   svfn_date,
                                                   svfn_time,
                                                   int(self.df),
                                                   self.metadata.ch_cmp.upper()))
        #calibrate electric channels 
        if self.metadata.ch_cmp == 'ex':
            time_series /= ex
        elif self.metadata.ch_cmp == 'ey':
            time_series /= ey

        #apply notch filter if desired
        if notch_dict is not None:
            self.apply_adaptive_notch_filter(notch_dict)
            print 'Filtered notches: '
            for nfilt in self.filt_list:
                if type(nfilt[0]) != str:
                    print '{0}{1:.2f} Hz'.format(' '*4, nfilt[0])
                                                   
        header_tuple = (save_station+self.metadata.rx_xyz0.split(':')[0], 
                        self.metadata.ch_cmp.lower(), 
                        self.df,
                        time.mktime(time.strptime(self.zen_schedule,
                                                  datetime_fmt )), 
                        time_series.shape[0], 
                        'mV', 
                        '{0:.3f}'.format(np.median(np.rad2deg(self.gps_stamps['lat']))), 
                        '{0:.3f}'.format(np.median(np.rad2deg(self.gps_stamps['lon']))), 
                        0.0, 
                        time_series)                
        self.fn_mt_ascii = mtfh.write_ts_file_from_tuple(save_fn, header_tuple,
                                                         fmt=fmt)
        
        print 'Wrote mtpy timeseries file to {0}'.format(self.fn_mt_ascii)
    
    #==================================================                                                           
    def plot_time_series(self, fig_num=1):
        """
        plots the time series
        """                                                               
        
        time_series = self.convert_counts()
        fig = plt.figure(fig_num )
        ax = fig.add_subplot(1,1,1)
        ax.plot(time_series)
        
        #ax.xaxis.set_minor_locator(MultipleLocator(self.df))
        #ax.xaxis.set_major_locator(MultipleLocator(self.df*15))
        #ax.xaxis.set_ticklabels([self.date_time[ii] 
        #                        for ii in range(0,len(self.date_time), 15)])
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude (mV)')
        plt.show()
        
        self.convert_mV()
        return fig, ax
    
    #==================================================    
    def plot_spectrogram(self, time_window=2**8, time_step=2**6, s_window=11,
                         frequency_window=1, n_freq_bins=2**9, sigma_L=None):
        """
        plot the spectrogram of the data using the S-method
        
        Arguments:
        -----------
            **s_window** : int (should be odd)
                    length of window for S-method calculation, higher numbers tend 
                    toward WVD
                    
            **time_window** : int (should be power of 2)
                     window length for each time step
                     *default* is 2**8 = 256
                     
            **frequency_window** : int (should be odd)
                     length of smoothing window along frequency plane
            
            **time_step** : int 
                        number of sample between short windows
                        *default* is 2**7 = 128
                     
            **sigmaL** : float
                     full width half max of gaussian window for L
            
            
            **n_freq_bins** : int 
                            (should be power of 2 and equal or larger than nh)
                            number of frequency bins
                            
        Returns:
        ---------
            **ptf** : mtpy.imaging.plotspectrogram.PlotTF object
            
        """
        
        time_series = self.convert_counts()
        
        kwargs = {'nh':time_window, 'tstep':time_step, 'L':s_window, 
                  'ng':frequency_window, 'df':self.df, 'nfbins':n_freq_bins,
                  'sigmaL': sigma_L}
        ptf = plotspectrogram.PlotTF(time_series, **kwargs)
        
        return ptf
        
    #==================================================
    def plot_spectra(self, fig_num=2):
        """
        plot the spectra of time series
        """
        if self.time_series is None:
            self.read_3d()
            
        time_series = self.convert_counts()
            
        spect = np.fft.fft(mtfilt.zero_pad(time_series))
        plot_freq = np.fft.fftfreq(spect.shape[0], 1./self.df)
        
        fig = plt.figure(fig_num, [4,4], dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.loglog(plot_freq, abs(spect)**2, lw=.5)
        ax.grid(which='both', lw=.25)
        
        ax.set_xlabel('Frequency (Hz)')
        #ax.set_xlim(1./plot_freq.max(), 1./plot_freq.min())
        ax.set_ylabel('Amplitude')
        
        plt.show()
        
        return fig, ax
        
#==============================================================================
#  for older Z3d files
#==============================================================================
class Zen3D_old(object):
    """
    Deal with the raw data output from the Zen box as Z3D files, which is in
    a formatted binary file with GPS time stamps every second.  Each time
    stamp contains information about position, GPS lock, temperature and a 
    few other things.  The stamp is 32 bytes long.  
    
    The read_file makes sure that there is the correct number of points 
    between time stamps, and the time stamp is correct.  The data read from 
    the file begins on the first coherent time stamp and ends on the last one.
    
    Usually the first coherent time stamp is a few seconds after scheduled 
    start time.  
    
    Arguments:
    ------------
        **fn**: string
                full path to .Z3D file to be manipulated
                
                
    ====================== ====================================================
    Methods                 Description
    ====================== ====================================================
    read_3d                read 3D file making sure all the time stamps are 
                           correctly spaced.  Returned time series starts at 
                           the first stamp which has the correct amount of data
                           points between it and the next time stamp.  Note
                           there are usually a few seconds at the end and maybe 
                           beginning that aren't correct because the internal 
                           computer is busy switchin sampling rate.
                     
    get_gps_stamp_location locates the gps stamp location
        
    get_gps_time           converts the gps counts to relative epoch seconds
                           according to gps week.      
    get_date_time          converts gps seconds into the actual date and time
                           in UTC.  Note this is different than GPS time which
                           is how the zen is scheduled, so the time will be 
                           off by the current amount of leap seconds.
    ====================== ====================================================
        
    =================== =======================================================
    Attributes           Description
    =================== =======================================================
    ch_adcard_sn        serial number of a/d card in channel
    ch_cmp              MT component of channel
    ch_length           distance between electrodes for channel, 
                        doesn't matter for magnetic channels
    ch_number           number of channel
    date_time           np.ndarray of date,time of gps stamps
    df                  sampling rate 
    fn                  full path to file name read in
    gps_diff            difference between gps time stamps
    gps_list             list of gps stamps
    gps_time            np.ndarray of gps times from time stamps
    gps_week            gps week
    header_dict         dictionary of header parameters
    log_lines           list of information to write into a log file later
    meta_dict           dictionary of meta data parameters
    rx_stn              name of station
    start_time          starting time and date of first time stamp with 
                         correct number of samples
    temperature         np.ndarray of temperature measurements at each time 
                        stamp
    time_series         np.ndarray of time series data in counts
    tx_id               name of transmitter if used
    units               [ 'counts' | 'mv' ] units of time series *default* is
                        counts. Plotting will convert to mV. 
    verbose             [ True | False ] for printing information to the  
                        interpreter
                         
    _data_type          np.dtype to convert binary formatted string
    _data_types         list of data types in binary formatted string
    _gps_epoch          gps_epoch in time.gmtime format. 
    _gps_stamp          string of gps_stamp 
    _header_len         length of header string in bytes. (512)
    _meta_len           length of meta data in bytes. (512)
    _raw_data           data in binary format
    _seconds_diff       difference in seconds from start time to look for 
                         gps stamp. *default* is 5
    _stamp_len          length of gps time stamp in bits
    _stamp_list          list of gps time stamp variables
    _week_len           length of a gps week in seconds
    =================== =======================================================
    """
    
    def __init__(self, fn=None, **kwargs):
        
        self.fn = fn
        self._header_len = kwargs.pop('header_len', 512)
        self._meta_len = kwargs.pop('meta_len', 512)
        self._stamp_len = kwargs.pop('stamp_len', 36)
        self._gps_stamp = kwargs.pop('gps_stamp', '\xff\xff\xff\xff')
        
        self._stamp_list = ['gps', 'time', 'lat', 'lon', 'status', 
                           'gps_accuracy', 'temperature']
                           
        self._data_types = [np.int32, np.int32, np.float64, np.float64, 
                            np.uint32, np.int32, np.float32]
                            
        self._data_type = np.dtype([(st, dt) for st, dt in 
                                     zip(self._stamp_list, self._data_types)])
                                     
        self._week_len = 604800
        self._gps_epoch = (1980, 1, 6, 0, 0, 0, -1, -1, 0)
        self._leap_seconds = 16
        
        #seconds different between scheduling time and actual collection time
        self._seconds_diff = 5 
        
        self.log_lines = []
        self.verbose = True
        self._skip_sample_tolerance = 5
        self.sample_diff_list = []
        self.counts_to_mv_conversion = 9.5367431640625e-10
        self.units = 'counts'
        self.gps_week = 1740
        self.time_series = None
        self.date_time = None
        
        self.header_dict = None
        self.df = None
        self.gain = None
        self.gps_week = None
        self.schedule_date = None
        self.schedule_time = None
        self.start_dt = None
        self.start_time = None
        self.start_date = None
        self.ch_adcard_sn = None
        
        self.meta_dict = None  
        self.ch_number = None
        self.ch_cmp = None
        self.ch_length = None
        self.rx_stn = None
        self.tx_id = None
        
        self.gps_diff = None
        self.gps_time = None
        self.gps_list = None
        self.temperature = None
        self.lat = None
        self.lon = None
        
    #==================================================    
    def read_header(self, header_string):
        """
        read header information and fill attribute:
        
            * header_dict   --> dictionary of head information
            * df            --> sampling frequency in Hz
            * gain          --> gain from within Zen box for that channel
            * gps_week      --> current gps week
            * schedule_time --> schedule start time
            * schedule_date --> schedule start date
            * start_dt      --> schedule start date and time
            * ch_adcard_sn  --> a/d card serial number
            
        **Note:** there are different versions of the header keywords from 
                  different generations of the Zen firmware.
        """
        
        #----read in header information----------------------------------------
        header_list = header_string.replace('\n', ',').split(',')
        
        header_dict = {}
        for hh in header_list:
            if hh != '' and hh.find('builddate') == -1:
                hkv = hh.split(':')
                if len(hkv) == 2:
                    if hkv[0].lower() == 'period' or \
                        hkv[0].lower() == 'duty':
                        try:
                            header_dict[hkv[0].strip().lower()] +=\
                                                                hkv[1].strip()
                        except KeyError:
                            header_dict[hkv[0].strip().lower()] =\
                                                                hkv[1].strip()
                    else:
                        header_dict[hkv[0].strip().lower()] = hkv[1].strip()
                elif len(hkv) == 3:
                    header_dict['start_time'] = hh.strip()
                else:
                    pass
            elif hh == '':
                pass
            else:
                hline = hh.split(';')
                for ll in hline:
                    if ll.find('builddate') > 0:
                        hlist = ll.split('&')
                        for kk in hlist:
                            klist = kk.split(':')
                            header_dict[klist[0].strip().lower()] = klist[1].strip()
                    else:
                        hlist = ll.split(':')
                        try:
                            header_dict[hlist[0].strip().lower()] = hlist[1].strip()
                        except IndexError:
                            pass
        #make attributes that will be useful latter
        self.header_dict = header_dict
        self.df = float(header_dict['a/d rate'])
        self.gain = float(header_dict['a/d gain'])
        self.gps_week = int(header_dict['gpsweek'])
        try:
            self.schedule_date = header_dict['schedule for this file']
        except KeyError:
            self.schedule_date = header_dict['schedule']
        self.schedule_time = header_dict['start_time']
        
        #get the start date/time in UTC time
        self.start_dt = self.compute_schedule_start(self.schedule_date, 
                                                      self.schedule_time)
        self.start_time = self.schedule_time
        self.start_date = self.schedule_date
                                            
        
        #--> get serial number of a/d board
        try:
            self.ch_adcard_sn = header_dict['serial']
        except KeyError:
            self.ch_adcard_sn = header_dict['brd339 serial']
    
    #==================================================        
    def read_metadata(self, meta_data_string):
        """
        read in meta data and make important information attributes
        
        Fills attributes:
        
            * meta_dict      --> dictionary of metadata 
            * ch_number      --> channel number
            * ch_cmp         --> channel component
            * ch_length      --> length of dipole
            * rx_stn         --> station name (can only be an integer)
            * tx.id          --> name of transmitter if used
            
        """
        meta_list = meta_data_string.replace('\n','|').split('|') 
        meta_dict = {}
        for mm in meta_list:
            mlist = mm.split(',')
            if len(mlist) == 2:
                meta_dict[mlist[0].strip().lower()] = mlist[1].strip().lower()
            else:
                pass
        self.meta_dict = meta_dict  
        self.ch_number = meta_dict['ch.number']
        self.ch_cmp = meta_dict['ch.cmp'].replace('b','h')
        self.ch_length = meta_dict['ch.varasp']
        self.rx_stn = meta_dict['rx.stn']
        self.tx_id = meta_dict['tx.id']
    
    #==================================================    
    def get_info(self):
        """
        read header and meta data
        """
        
        #beginning index of data blocks
        ds = self._header_len+self._meta_len
        
        #read in as a binary file.
        rfid = open(self.fn, 'rb')
        raw_data = rfid.read(ds+4)
        self._raw_data = raw_data
        rfid.close()
        
        if len(raw_data) < ds:
            print 'Data file is not complete cannot read header information'
            return

        try:
            self.log_lines[0] != '-'*72+'\n'
        except IndexError:
            self.log_lines.append('-'*72+'\n')
            self.log_lines.append('--> Reading File: {0}\n'.format(self.fn))
        
        #----read in header information----------------------------------------
        header_string = raw_data[0:self._header_len]
        self.read_header(header_string)
        
        print('-'*40)
        print('   ad card sn     =  {0}'.format(self.ch_adcard_sn))
        print('   sampling rate  =  {0:.0f}'.format(self.df))
        print('   gain           =  {0:.1f}'.format(self.gain))
        print('   gps_week       =  {0:.0f}'.format(self.gps_week))
        print('   schedule date  =  {0}'.format(self.schedule_date))
        print('   schedule time  =  {0}'.format(self.schedule_time))
        
        #---read in meta raw_data----------------------------------------------
        meta_string = raw_data[self._header_len-1:ds]
        self.read_metadata(meta_string)

        print('   channel no     =  {0}'.format(self.ch_number))
        print('   channel comp   =  {0}'.format(self.ch_cmp))
        print('   channel len    =  {0}'.format(self.ch_length))
        print('   rx station     =  {0}'.format(self.rx_stn))
        print('   tx id          =  {0}'.format(self.tx_id))
        print('-'*40)
 
    #==================================================
    def read_3d(self):
        """
        read in the time series and gps time stamps.
        
        Makes sure that the number of samples between each time stamp is
        the sampling rate.  If it is not an error is raised if the difference
        is more than _skip_sample_tolerance.  
        
        Creates a time series that starts at the time where the first gps
        time stamp has the correct number of points, and stops where the first
        incorrect number of points occurs.  A corresponding time,date array
        is created.

        """
        #read in as a binary file.
        raw_data = open(self.fn, 'rb').read()
        self._raw_data = raw_data
        
        try:
            self.log_lines[0] != '-'*72+'\n'
        except IndexError:
            self.log_lines.append('-'*72+'\n')
            self.log_lines.append('--> Reading File: {0}\n'.format(self.fn))
        
        #number of bytes in the file
        num_bytes = len(raw_data)
        
        #beginning index of data blocks
        ds = self._header_len+self._meta_len
        
        #----read in header information----------------------------------------
        header_string = raw_data[0:self._header_len]
        self.read_header(header_string)
        
        #---read in meta raw_data----------------------------------------------
        meta_string = raw_data[self._header_len-1:ds]
        self.read_metadata(meta_string)
        
                
        #---read in gps raw_data-----------------------------------------------
        #sampling rate times 4 bytes for 32 bit measurement
        df = int(self.df)      
        dt = df*4
        
        #length of data block plus gps stamp
        block_len = self._stamp_len+dt
        
        #number of data blocks
        num_blocks = int(np.ceil(num_bytes/float(block_len)))
        
        #get position of gps stamps
        gps_list = np.zeros(num_blocks, dtype=np.int)
        
        gps_dict = dict([(key, np.zeros(num_blocks, dtype=dtp)) 
                          for key, dtp in zip(self._stamp_list, 
                                              self._data_types)])
        #make the time array floats instead of ints so can get the decimal 
        #place if it isn't 0.
        gps_dict['time'] = gps_dict['time'].astype(np.float32)
        
        #get gps information from the data
        #get first time stamp that matches the starting time
        s1 = 0
        gps_list[0] = self.get_gps_stamp_location()
        gps_info = np.fromstring(raw_data[gps_list[0]:gps_list[0]+self._stamp_len], 
                                 dtype=self._data_type)
        gps_info['time'] = gps_info['time'].astype(np.float32)
        gps_info['time'] = self.get_gps_time(gps_info['time'])[0]
        start_test = self.get_date_time(self.gps_week, gps_info['time'])
        
        #--> test to make sure the first time corresponds to the scheduled 
        #start time
        time_stop = 0
        while start_test != self.start_dt and s1 <= self._seconds_diff and \
                time_stop <= self._seconds_diff:
            s1 += 1
            gps_list[0] = self.get_gps_stamp_location(gps_list[0]+7)
            gps_info = np.fromstring(raw_data[gps_list[0]:gps_list[0]+\
                                                self._stamp_len], 
                                     dtype=self._data_type)
                                     
            gps_info['time'] = gps_info['time'].astype(np.float32)
            gps_info['time'], gps_dweek = self.get_gps_time(gps_info['time'])
            
            start_test = self.get_date_time(self.gps_week+gps_dweek, 
                                            gps_info['time'])
            if s1 == self._seconds_diff:
                s1 = 0
                self.start_dt = self.start_dt[:-2]+\
                                 '{0:02}'.format(int(self.start_dt[-2:])+1)
                gps_list[0] = self.get_gps_stamp_location()
                time_stop += 1  
       
        #----Raise an error if the first gps stamp is more than allowed time
        #    difference.
        if time_stop >= self._seconds_diff:
            print ('GPS start time is more than '+\
                           '{0} '.format(self._seconds_diff)+\
                           'seconds different than scheduled start time of '+\
                           '{0}. \n '.format(self.start_dt)+\
                           'Estimated start time is {0} +/- {1} sec'.format(
                           start_test, self._seconds_diff))
                     
        #put the information into the correct arrays via dictionary                         
        for jj, key in enumerate(self._stamp_list):
            gps_dict[key][0] = gps_info[0][jj]
  
        #find the next time stamp
        for ii in range(s1,num_blocks-1):
            sfind = self.get_gps_stamp_location(gps_list[ii-1]+7)
            #make sure it isn't the same time stamp as before
            if sfind != gps_list[ii-1] and sfind != -1:
                gps_info, gps_index, gps_week = self.get_gps_stamp(sfind)
                gps_list[ii] = gps_index
                
                if gps_info is not None:
                    for jj, key in enumerate(self._stamp_list):
                        gps_dict[key][ii] = gps_info[0][jj]
        
        #get only the values that are non zero
        gps_dict['time'] = gps_dict['time'][np.nonzero(gps_dict['time'])] 

        num_samples = len(gps_dict['time'])
        
        #calculate the difference between time stamps
        gps_diff = np.array([gps_dict['time'][ii+1]-gps_dict['time'][ii] 
                             for ii in range(num_samples-1)])
        
        #check for any spots where gps was not locked or mised a sampling interval
        bad_lock = np.where(gps_diff[np.nonzero(gps_diff)] != 1.0)[0]
        
        if len(bad_lock) > 0:
            for bb in bad_lock:
                if gps_diff[bb] > 5:
                    self.log_lines.append(' '*4+\
                                      'point {0:^15},'.format(gps_list[bb])+\
                                      'gps diff {0:^15}\n'.format(gps_diff[bb]))
            
            self.log_lines.append(' '*4+'*'*52+'\n')

        #need to be sure that the number of data points between time stamps is 
        #equal to the sampling rate, if it is not then remove that interval.  
        #Most likely it is at the beginning or end of time series.
        dsamples = np.array([(gps_list[nn+1]-gps_list[nn]-self._stamp_len-df*4)/4 
                              for nn in range(num_samples)])
        
        bad_interval = np.where(abs(dsamples)>self._skip_sample_tolerance)[0]
        bmin = 0
        bmax = num_samples
        if len(bad_interval) > 0:        
            #need to locate the bad interval numbers
            for bb in bad_interval:
                if bb <= 10:
                    bmin = bb+1
                if bb > num_samples-10:
                    bmax = bb
        
            gps_list = gps_list[bmin:bmax]
            
        num_samples = len(gps_list)
        if self.verbose:
            print 'Found {0} gps time stamps, '.format(num_samples)+\
                  'with equal intervals of {0} samples'.format(int(self.df))
              
        self.log_lines.append(' '*4+\
                            'Found {0} gps time stamps, '.format(num_samples)+\
                  'with equal intervals of {0} samples\n'.format(int(self.df)))
        
        #read in data
        data_array = np.zeros((num_samples+1)*df, dtype=np.float32)
        for ll, kk in enumerate(gps_list[0:-1]):
            pdiff = ((gps_list[ll+1]-(kk+self._stamp_len))-(df*4))/4
            self.sample_diff_list.append(pdiff)
            dblock = raw_data[kk+self._stamp_len:gps_list[ll+1]] 
            try:
                data_array[ll*df:(ll+1)*df+pdiff] = np.fromstring(dblock, 
                                                                dtype=np.int32)
            except ValueError:
                print 'samples between time step {0} is off by {1} samples'.format(ll,
                                                                   abs(pdiff))
                          
        if sum(self.sample_diff_list) != 0:
            if self.verbose:
                print 'time series is off by {0} seconds'.format(
                                           float(sum(self.sample_diff_list))/df)
                self.log_lines.append('time series is off by {0} seconds'.format(
                                          float(sum(self.sample_diff_list))/df))
                                           
        #get only the non-zero data bits, this is dangerous if there is 
        #actually an exact 0 in the data, but rarely happens 
        self.time_series = data_array[np.nonzero(data_array)]
        
        #need to cut all the data arrays to have the same length and corresponding 
        #data points
        for key in gps_dict.keys():
            gps_dict[key] = gps_dict[key][bmin:bmax]
        
        #make attributes of imporant information
        self.gps_diff = gps_diff[bmin:bmax]
        self.gps_time = gps_dict['time']
        self.gps_list = gps_list
        self.temperature = gps_dict['temperature']
        self.lat = gps_dict['lat']
        self.lon = gps_dict['lon']

        self.date_time = np.zeros_like(gps_dict['time'], dtype='|S24')

        for gg, gtime in enumerate(gps_dict['time']):
            self.date_time[gg]= self.get_date_time(self.gps_week, gtime)
        
        try:
            self.start_dt = self.date_time[0]
            self.start_date = self.date_time[0].split(',')[0]
            self.start_time = self.date_time[0].split(',')[1]
            if self.verbose:
                print 'Starting time of time series is '+\
                        '{0} UTC'.format(self.date_time[0])
            self.log_lines.append(' '*4+'Starting time of time series is '+\
                                  '{0} UTC\n'.format(self.date_time[0]))
        except IndexError:
            print 'No quality data was collected'
            self.log_lines.append(' '*4+'No quality data was collected\n')
            self.start_dt = None
            self.start_date = None
            self.start_time = None
            
        if self.units == 'mv':
            self.time_series = self.convert_counts()
    
    #==================================================    
    def convert_counts(self):
        """
        convert the time series from counts to millivolts

        """
        
        return self.time_series*self.counts_to_mv_conversion
    
    #==================================================    
    def convert_mV(self):
        """
        convert millivolts to counts assuming no other scaling has been applied
        
        """
        
        return self.time_series/self.counts_to_mv_conversion
    
    #==================================================    
    def compute_schedule_start(self, start_date, start_time, 
                               leap_seconds=None):
        """
        compute the GMT time for scheduling from start time of the gps 
        according to the leap seconds.
        
        Arguments:
        -----------
            **start_date**: YYYY-MM-DD
                            schedule start date
                            
            **start_time**: hh:mm:ss
                            time of schedule start on a 24 hour basis
            
            **leap_seconds**: int
                              number of seconds that GPS is off from UTC time.
                              as of 2013 GPS is ahead by 16 seconds.
                              
        Returns:
        --------
            **ndate_time**: YYYY-MM-DD,hh:mm:ss
                            calibrated date and time in UTC time.
        
        """                                
        month_dict = {1:31, 2:28, 3:31, 4:30, 5:31, 6:30, 7:31, 8:31, 9:30, 
                      10:31, 11:30, 12:31}
        if leap_seconds is not None:
            self._leap_seconds = leap_seconds
            
        year, month, day = start_date.split('-')
        
        hour, minutes, seconds = start_time.split(':')
        
        new_year = int(year)
        new_month = int(month)
        new_day = int(day)
        new_hour = int(hour)
        new_minutes = int(minutes)
        new_seconds = int(seconds)-self._leap_seconds
       
        if new_seconds < 0:
            new_seconds = (int(seconds)-self._leap_seconds)%60
            new_minutes = int(minutes)-1
            if new_minutes < 0:
                new_minutes = (int(minutes)-1)%60
                new_hour = int(hour)-1
                if new_hour < 0:
                    new_hour = (int(hour)-1)%24
                    new_day = int(day)-1
                    if new_day <= 0:
                        new_day = (int(day)-1)%30
                        new_month = int(month)-1
                        if new_month <= 0:
                            new_month = (12-new_month)
                        new_day = month_dict[new_month]-int(day)+1
                        print 'need to check date, have not implemented '+\
                              'leap years yet'
                     
                              
        ndate_time = time.strftime(datetime_fmt ,
                                   (new_year, 
                                    new_month, 
                                    new_day, 
                                    new_hour, 
                                    new_minutes, 
                                    new_seconds, 0, 0, 0))
                                    
        return ndate_time
    
    #==================================================
    def get_gps_stamp_location(self, start_index=None):
        """
        get the location in the data file where there is a gps stamp.  Makes
        sure that the time stamp is what it should be.
        
        Arguments:
        -----------
            **start_index**: int
                             starting point to look for the time stamp within
                             the file.
                             
        Returns:
        ---------
            **gps_index**: int
                           the index in the file where the start of the 
                           time stamp is.
        
        """
        
        gps_index = self._raw_data.find(self._gps_stamp, start_index)
        if self._raw_data[gps_index+4] == '\xff':
            gps_index += 1
            if self._raw_data[gps_index+4] == '\xff':
                gps_index += 1
                if self._raw_data[gps_index+4] == '\xff':
                    gps_index += 1
                    if self._raw_data[gps_index+4] == '\xff':
                        gps_index += 1
                        
        return gps_index
    
    #==================================================    
    def get_gps_stamp(self, gps_index):
        """
        get the gps stamp data
        
        """
        #get numbers from binary format
        try:
            
            gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                     dtype=self._data_type)
            while gps_info['time'] < 0:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'time ', gps_index
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
            
            while gps_info['status'] < 0:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'status ', gps_index                
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
            
            while abs(gps_info['temperature']) > 80:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'temperature ', gps_index    
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
                
            while abs(gps_info['lat']) > np.pi:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'lat ', gps_index
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
                
            #convert lat and lon into decimal degrees
            gps_info['lat'] = self.get_degrees(gps_info['lat'])
            gps_info['lon'] = self.get_degrees(gps_info['lon'])
            gps_info['time'] = gps_info['time'].astype(np.float32)
            gps_info['time'], gps_week = self.get_gps_time(gps_info['time'])
            
            if gps_info == []:
                print gps_index
                raise ZenGPSError('Something is fucked')
            if gps_index == -1:
                print gps_info
                raise ZenGPSError('Something is fucked')
                
            return gps_info, gps_index, gps_week 
            
        except ValueError:
            print 'Ran into end of file, gps stamp not complete.'+\
                  ' Only {0} points.'.format(len(self._raw_data[gps_index:]))
            return None, gps_index, 0
    
    #==================================================        
    def get_gps_time(self, gps_int, gps_week=0):
        """
        from the gps integer get the time in seconds.
        
        Arguments:
        ----------
            **gps_int**: int
                         integer from the gps time stamp line
                         
            **gps_week**: int
                          relative gps week, if the number of seconds is 
                          larger than a week then a week is subtracted from 
                          the seconds and computed from gps_week += 1
                          
        Returns:
        ---------
            **gps_time**: int
                          number of seconds from the beginning of the relative
                          gps week.
        
        """
            
        gps_seconds = gps_int/1024.
        
        gps_ms = (gps_seconds-np.floor(gps_int/1024.))*(1.024)
        
        cc = 0
        if gps_seconds > self._week_len:
            gps_week += 1
            cc = gps_week*self._week_len
            gps_seconds -= self._week_len
        
        gps_time = np.floor(gps_seconds)+gps_ms+cc
        
        return gps_time, gps_week
    
    #==================================================
    def get_date_time(self, gps_week, gps_time):
        """
        get the actual date and time of measurement as UTC. 
        
        Note that GPS time is curently off by 16 seconds from actual UTC time.
        
        Arguments:
        ----------
            **gps_week**: int
                          integer value of gps_week that the data was collected
            
            **gps_time**: int
                          number of seconds from beginning of gps_week
            
            **leap_seconds**: int
                              number of seconds gps time is off from UTC time.
                              It is currently off by 16 seconds.
                              
        Returns:
        --------
            **date_time**: YYYY-MM-DD,HH:MM:SS
                           formated date and time from gps seconds.
        
        
        """
        
        mseconds = gps_time % 1
        
        #make epoch in seconds, mktime computes local time, need to subtract
        #time zone to get UTC
        epoch_seconds = time.mktime(self._gps_epoch)-time.timezone
        
        #gps time is 14 seconds ahead of GTC time, but I think that the zen
        #receiver accounts for that so we will leave leap seconds to be 0        
        gps_seconds = epoch_seconds+(gps_week*self._week_len)+gps_time-\
                                                        self._leap_seconds

        #compute date and time from seconds
        (year, month, day, hour, minutes, seconds, dow, jday, dls) = \
                                                    time.gmtime(gps_seconds)
        
        date_time = time.strftime(datetime_fmt ,(year,
                                                 month, 
                                                 day, 
                                                 hour, 
                                                 minutes, 
                                                 int(seconds+mseconds), 
                                                 0, 0, 0))
        return date_time
    
    #==================================================    
    def get_degrees(self, radian_value):
        """
        convert lat or lon into decimal degrees
        
        """
        
        degrees = radian_value*180/np.pi
        
        return degrees
    
    #==================================================    
    def apply_adaptive_notch_filter(self, notch_dict):
        """
        apply notch filter to the data that finds the peak around each 
        frequency.
        
        see mtpy.processing.filter.adaptive_notch_filter
        
        
        """
        
        try:
            self.time_series
        except AttributeError:
            self.read_3d()
        
        notches = notch_dict.pop('notches', list(np.arange(60,2048,60)))
        notchradius = notch_dict.pop('notchradius', 0.5)
        freqrad = notch_dict.pop('freqrad', 0.5)
        rp = notch_dict.pop('rp', 0.1)
        kwargs = {'df':self.df, 'notches':notches, 'notchradius':notchradius,
                  'freqrad':freqrad, 'rp':rp}
                  
        self.time_series, self.filt_list = \
                    mtfilt.adaptive_notch_filter(self.time_series, **kwargs) 
        
    #==================================================
    def write_ascii_mt_file(self, save_fn=None, save_station='mb', fmt='%.8e',
                            ex=1, ey=1, notch_dict=None):
        """
        write an mtpy time series data file
        
        Arguments:
        -----------
            **save_fn** : full path to save file, if None file is saved as:
                          station_YYYYMMDD_hhmmss_df.component
                          
                          ex. mt01_20130206_120000_256.HX
                          
            **save_station** : string
                               prefix string to add to station number as only
                               integers can be input into metadata of the zen
                               boxes.  ex. mb001
            
            **fmt** : string format
                      format of data numbers output to ascii file.
                      *default* is '%.8e' for 8 significan figures in 
                      scientific notation.
                      
        Output:
        --------
            **fn_mt_ascii** : full path to saved file
        
        """
        try:
            self.start_date
        except AttributeError:
            self.read_3d()
            
        time_series = self.convert_counts()
        if save_fn is None:
            svfn_directory = os.path.join(os.path.dirname(self.fn), 'TS')
            if not os.path.exists(svfn_directory):
                os.mkdir(svfn_directory)
                
            svfn_date = ''.join(self.start_date.split('-'))
            svfn_time = ''.join(self.start_time.split(':'))
            svfn_station = save_station+self.rx_stn
            save_fn = os.path.join(svfn_directory, 
                                   '{0}_{1}_{2}_{3}.{4}'.format(svfn_station,
                                                   svfn_date,
                                                   svfn_time,
                                                   int(self.df),
                                                   self.ch_cmp.upper()))
        #calibrate electric channels 
        if self.ch_cmp == 'ex':
            time_series /= ex
        elif self.ch_cmp == 'ey':
            time_series /= ey

        #apply notch filter if desired
        if notch_dict is not None:
            self.apply_adaptive_notch_filter(notch_dict)
            print 'Filtered notches: '
            for nfilt in self.filt_list:
                if type(nfilt[0]) != str:
                    print '{0}{1:.2f} Hz'.format(' '*4, nfilt[0])
                                                   
        header_tuple = (save_station+self.rx_stn, 
                        self.ch_cmp, 
                        self.df,
                        time.mktime(time.strptime(self.start_dt,
                                                  datetime_fmt )), 
                        time_series.shape[0], 
                        'mV', 
                        np.median(self.lat), 
                        np.median(self.lon), 
                        0.0, 
                        time_series)
                        
        self.fn_mt_ascii = mtfh.write_ts_file_from_tuple(save_fn, header_tuple,
                                                         fmt=fmt)
        
        print 'Wrote mtpy timeseries file to {0}'.format(self.fn_mt_ascii)
   
    #==================================================     
    def write_mseed_mt_file(self, save_fn=None, save_station='mb', 
                            location='Mono Basin', network='USGS'):
        """
        write a miniseed file, note need to have Obspy installed.  This proves
        to be difficult under windows
        
        Arguments:
        ----------
            **save_fn** : full path to file to save data to
            
            **save_station** : string 
                               prefix to add onto station name
            
            **location** : string
                           description of where the data was collected
            
            **network** : string
                          network or company that collected the data
                          
        Outputs:
        --------
            **save_fn** : string
                          full path to file where data was saved.
        """
        
        try:
            self.start_date
        except AttributeError:
            self.read_3d()
            
        time_series = self.convert_counts()
        
        svfn_date = ''.join(self.start_date.split('-'))
        svfn_time = ''.join(self.start_time.split(':'))
        svfn_station = save_station+self.rx_stn
        svfn_chn = self.ch_cmp.upper()
        delta_t = 1./self.df
        t0 = self.start_dt.replace(',','T')
        if save_fn is None:
            save_fn = os.path.join(os.path.dirname(self.fn), 
                               '{0}_{1}_{2}_{3}_{4}.mseed'.format(svfn_station,
                                                   svfn_date,
                                                   svfn_time,
                                                   int(self.df),
                                                   svfn_chn))
                                            
                                                   
        self.fn_mt_mseed = mtmseed.writefile_obspy_singletrace(save_fn, 
                                                               svfn_station,
                                                               svfn_chn,
                                                               network,
                                                               location,
                                                               delta_t,
                                                               t0,
                                                               time_series)
        return save_fn
    
    #==================================================                                                           
    def plot_time_series(self, fig_num=1):
        """
        plots the time series
        """                                                               
        
        time_series = self.convert_counts()
        fig = plt.figure(fig_num, dpi=300)
        ax = fig.add_subplot(1,1,1)
        ax.plot(time_series)
        
        #ax.xaxis.set_minor_locator(MultipleLocator(self.df))
        #ax.xaxis.set_major_locator(MultipleLocator(self.df*15))
        #ax.xaxis.set_ticklabels([self.date_time[ii] 
        #                        for ii in range(0,len(self.date_time), 15)])
        
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Amplitude (mV)')
        plt.show()
        
        self.convert_mV()
        return fig, ax
    
    #==================================================    
    def plot_spectrogram(self, time_window=2**8, time_step=2**6, s_window=11,
                         frequency_window=1, n_freq_bins=2**9, sigma_L=None):
        """
        plot the spectrogram of the data using the S-method
        
        Arguments:
        -----------
            **s_window** : int (should be odd)
                    length of window for S-method calculation, higher numbers tend 
                    toward WVD
                    
            **time_window** : int (should be power of 2)
                     window length for each time step
                     *default* is 2**8 = 256
                     
            **frequency_window** : int (should be odd)
                     length of smoothing window along frequency plane
            
            **time_step** : int 
                        number of sample between short windows
                        *default* is 2**7 = 128
                     
            **sigmaL** : float
                     full width half max of gaussian window for L
            
            
            **n_freq_bins** : int 
                            (should be power of 2 and equal or larger than nh)
                            number of frequency bins
                            
        Returns:
        ---------
            **ptf** : mtpy.imaging.plotspectrogram.PlotTF object
            
        """
        
        time_series = self.convert_counts()
        
        kwargs = {'nh':time_window, 'tstep':time_step, 'L':s_window, 
                  'ng':frequency_window, 'df':self.df, 'nfbins':n_freq_bins,
                  'sigmaL': sigma_L}
        ptf = plotspectrogram.PlotTF(time_series, **kwargs)
        
        return ptf
        
    #==================================================
    def plot_spectra(self, fig_num=2):
        """
        plot the spectra of time series
        """
        if self.time_series is None:
            self.read_3d()
            
        time_series = self.convert_counts()
            
        spect = np.fft.fft(mtfilt.zero_pad(time_series))
        plot_freq = np.fft.fftfreq(spect.shape[0], 1./self.df)
        
        fig = plt.figure(fig_num, [4,4], dpi=200)
        ax = fig.add_subplot(1,1,1)
        ax.loglog(plot_freq, abs(spect)**2, lw=.5)
        ax.grid(which='both', lw=.25)
        
        ax.set_xlabel('Frequency (Hz)')
        #ax.set_xlim(1./plot_freq.max(), 1./plot_freq.min())
        ax.set_ylabel('Amplitude')
        
        plt.show()
        
        return fig, ax
        
#==============================================================================
# Cache files
#==============================================================================
class Cache_Metadata(object):
    def __init__(self, fn=None, **kwargs):
        self.fn = fn
        self.ch_adcardsn = None
        self.ch_azimuth = None
        self.ch_cmp = None
        self.ch_cres = None
        self.ch_factor = None
        self.ch_gain = None
        self.ch_gainfactor = None
        self.ch_gdpslot = None
        self.ch_length = None
        self.ch_lowpass = None
        self.ch_number = None
        self.ch_numon = None
        self.data_version = None
        self.gdp_cardtype = None
        self.gdp_date = None
        self.gdp_operator = None
        self.gdp_progver = None
        self.gdp_time = None
        self.gdp_type = None
        self.gps_alt = None
        self.gps_lat = None
        self.gps_lon = None
        self.gps_numsat = None
        self.gps_sec = None
        self.gps_utmzone = None
        self.gps_week = None
        self.header_type = None
        self.job_by = None
        self.job_for = None
        self.job_name = None 
        self.job_number = None
        self.line_name = None
        self.rx_aspace = None
        self.rx_sspace = None
        self.rx_utm0 = None
        self.rx_utm1 = None
        self.rx_utm2 = None
        self.rx_xyz0 = None
        self.rx_xyz1 = None
        self.rx_xyz2 = None
        self.survey_acqmethod = None
        self.survey_type = None
        self.ts_adfreq = None
        self.ts_npnt = None
        self.unit_length = None
        self.station_number = None
        

        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def read_meta_string(self, meta_string=None):
        """
        read in a meta from the raw string
        """
        
        if meta_string is not None:
            self._meta_string = meta_string
            
            meta_list = self._meta_string.split('\n')
            for m_str in meta_list:
                line_list = m_str.strip().split(',')
                l_key = line_list[0].replace('.', '_').lower()
                l_value = line_list[1:]
                if len(l_value) == 1:
                    try:
                        l_value = float(l_value[0])
                    except ValueError:
                        l_value = l_value[0]
                setattr(self, l_key, l_value)
            self._get_station_number()
    
    def _get_station_number(self):
        """
        get station name from metadata from all versions of .cac files
        """
        
        try: 
            self.station_number = str(int(self.rx_stn))
        except AttributeError:
            try:
                self.station_number = self.rx_xyz0.split(':')[0]
            except AttributeError:
                print ('Could not find station number in rx.stn or rx.xyz0'
                        ' setting station_number to 000')
                
class Board_Calibration(object):
    """
    deal with baord calibration 
    """
    
    def __init__(self, board_cal_str=None, **kwargs):
        self.board_cal_str = board_cal_str
        
        self.cal_sys = {}
        self.cal_ant = {}
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def read_board_cal_str(self, board_cal_str=None):
        """
        read board calibration data
        """
        
        if board_cal_str is not None:
            self.board_cal_str = board_cal_str
            
            
        cal_list = self.board_cal_str.split('\n')
        for c_str in cal_list:
            c_list = c_str.split(',')
            c_key = c_list[0].replace('.', '_').lower()
            if len(c_list) == 2:
               c_value = c_list[1]
               setattr(self, c_key, c_value)
            elif len(c_list) > 2:
                c_key2 = c_list[1]
                c_arr = np.zeros(len(c_list[2:]), 
                                 dtype=[('frequency', np.float),
                                        ('amplitude', np.float),
                                        ('phase', np.float)])
                for ii, cc in enumerate(c_list[2:]):
                    c_arr[ii] = np.array([float(kk) for kk in cc.split(':')])
        
                self.__dict__[c_key][c_key2] = c_arr
                
class Cache(object):
    """
    deal with Zonge .cac files
    """
    def __init__(self, fn=None, **kwargs):
        self.fn = fn
        
        self.metadata = None
        self.time_series = None
        self.other = None
        self.calibration = None
        
        self._flag_len = 10
        self._len_bytes = 4
        self._flag_dtype = [('length', np.int32),
                            ('flag', np.int32),
                            ('type', np.int16)]
              

        self._type_dict = {4:'navigation',
                           514:'metadata',
                           768:'calibration',
                           16:'time_series',
                           15:'other',
                           640:'status'}
                           
        self._f_tell = 0
                
    def _read_file_block(self, file_id):
        """
        read a cache block
        """
        file_pointer = np.fromstring(file_id.read(self._flag_len),
                                     dtype=self._flag_dtype)
        f_str = file_id.read(file_pointer['length']-2)
        end_len = np.fromstring(file_id.read(self._len_bytes),
                                dtype=np.int32)
        
        if self._validate_block_len(file_pointer, end_len) is True:
            self._f_tell = file_id.tell()
            return file_pointer, f_str
        
    def _validate_block_len(self, file_pointer, end_length):
        """
        validate that the block lengths as defined at the beginning and 
        the end are the same
        """
        
        try:
            assert file_pointer['length'] == end_length
            return True
        except AssertionError:
            raise ValueError('File pointer length {0} != end length {1}'.format(
                             file_pointer['length'], end_length))
                             
    def read_cache_metadata(self, fn=None):
        """
        read .cac file
        """
        if fn is not None:
            self.fn = fn
        
        f_pointer = True
        with open(self.fn, 'rb') as fid:
            while f_pointer:            
                # read in first pointer            
                f_pointer, f_str = self._read_file_block(fid)
                
                # if the data type is the meta data
                if int(f_pointer['type']) == 514:
                    meta_obj = Cache_Metadata()
                    meta_obj.read_meta_string(f_str)

                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, meta_obj)
                    print 'Read in metadata'
                    return
                            
        
    def read_cache_file(self, fn=None):
        """
        read .cac file
        """
        if fn is not None:
            self.fn = fn
        
        f_pointer = True
        with open(self.fn, 'rb') as fid:
            while f_pointer:            
                # read in first pointer            
                f_pointer, f_str = self._read_file_block(fid)
                
                # if the data type is the meta data
                if int(f_pointer['type']) == 514:
                    meta_obj = Cache_Metadata()
                    meta_obj.read_meta_string(f_str)

                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, meta_obj)
                    print 'Read in metadata'
                    continue
                
                # if the data type is calibration
                elif int(f_pointer['type']) == 768:
                    cal_obj = Board_Calibration(f_str)
                    cal_obj.read_board_cal_str()
                    
                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, cal_obj)
                    print 'Read in calibration'
                    continue
                    
                # if the data type is time series
                elif int(f_pointer['type']) == 16:
                    ts_arr = np.fromstring(f_str, dtype=np.int32)
                    ts_arr = np.resize(ts_arr, (int(self.metadata.ts_npnt), 
                                                len(self.metadata.ch_cmp)))
                    
                    
                    ts = np.zeros(1, 
                                  dtype=[(cc.lower(), np.int32, 
                                          (int(self.metadata.ts_npnt),)) for 
                                          cc in self.metadata.ch_cmp])
                    
                    for ii, cc in enumerate(self.metadata.ch_cmp):
                        ts[cc.lower()][:] = ts_arr[:, ii]
                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, ts)
                    print 'Read in time series,  # points = {0}'.format(
                                                        self.metadata.ts_npnt)
                    return
                # if the data type is time series
                elif int(f_pointer['type']) == 15:
                    ts = np.fromstring(f_str, dtype=np.int32)
                    

                
                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, ts)
                    print 'Read in other'
                    continue
                
#==============================================================================
# 
#==============================================================================
class ZenCache(object):
    """
    deals with cache files or combined time series files.
    
    This will combine all coincident files into a .cac file for use in the
    Zonge processing software.  It will start at the first coherent time
    stamp and go to the longest coherent time stamp, such that each channel 
    will have the same start time and end time and same number of data points.
    
    ================== ========================================================
     Attributes         Description
    ================== ======================================================== 
    cal_data            list of calibrations, as is from file
    fn_list              list of filenames merged together
    log_lines           list of information to put into a log file
    meta_data           dictionary of meta data key words and values
    nav_data            list of navigation data, as is from file
    save_fn             file to save merged file to
    ts                  np.ndarray(len(ts), num_channels) of time series
    verbose             [ True | False ] True prints information to console
    zt_list              list of class: Zen3D objects
    _ch_factor          scaling factor for the channels, got this from Zonge
    _ch_gain            gain on channel, not sure of the format
    _ch_lowpass_dict    dictionary of values for lowpass filter, not sure how
                        they get the values
    _data_type          np.dtype of data type for cache block
    _flag               flag for new data block
    _nav_len            length of navigation information in bytes
    _stamp_len          length of new block stamp in bytes
    _type_dict          dictionary of cache block types, from Zonge.
    ================== ========================================================
    
    Methods:
    ---------
        * *check_sampling_rate* : makes sure sampling rate is the same for all
                                  files being merged.
        
        * *check_time_series* : makes sure all time series start at the same
                                time and have the same length.
                                
        * *write_cache_file* : writes a cache file for given filenames.
        
        * *read_cache* : reads in a cache file.
        
    :Example: ::
    
        >>> import ZenTools as zen
        >>> zc = zen.ZenCache()
        >>> # read a cache file
        >>> zc.read_cache(fn=r"/home/MT/mt01_20130601_190001_4096.cac")
        >>> # write a cache file
        >>> import os
        >>> file_path = r"/home/MT/Station1"
        >>> fn_list = [os.path.join(file_path, fn) 
        >>> ...       for fn in os.listdir(file_path)
        >>> ...       if fn.find('.Z3D')>0]
        >>> zc.write_cache_file(fn_list, r"/home/MT/Station1", station='s1')
        >>> Saved File to: /home/MT/Station1/Merged/s1_20130601_190001_4096.cac
        
    """   
    
    def __init__(self):
        
        self.fn_list = None
        self.zt_list = None
        self.save_fn = None
        self._ch_factor = '9.5367431640625e-10'
        self._ch_gain = '01-0'
        self._ch_lowpass_dict = {'256':'112', 
                           '1024':'576',
                           '4096':'1792'}
        self._flag = -1
        self._ts_dtype = np.int32
        
        self._type_dict = {'nav' : 4,
                          'meta' :  514,
                          'cal' : 768,
                          'ts' : 16}
                          
        self._data_type = np.dtype([('len',np.int32), 
                                    ('flag', np.int32), 
                                    ('input_type', np.int16)])
        self._stamp_len = 10
        self._nav_len = 43
        
        self.nav_data = None
        self.cal_data = None
        self.ts = None
        self.verbose = True
        self.log_lines = []
        self.chn_order = ['hx','hy','hz','ex','ey']
        
        self.meta_data = {'SURVEY.ACQMETHOD' : ',timeseries',
                           'SURVEY.TYPE' : ',',
                           'LENGTH.UNITS' : ',m',
                           'DATA.DATE0' : '',
                           'DATA.TIME0' : '',
                           'TS.ADFREQ' : '',
                           'TS.NPNT': '',              
                           'CH.NUNOM' : ',',
                           'CH.FACTOR' : '',
                           'CH.GAIN' : '',
                           'CH.NUMBER' : '',
                           'CH.CMP' : '',
                           'CH.LENGTH' : '',
                           'CH.EXTGAIN' : '',
                           'CH.NOTCH' : '',
                           'CH.HIGHPASS' : '',
                           'CH.LOWPASS' : '',
                           'CH.ADCARDSN' : '',
                           'CH.STATUS' : ',',
                           'CH.SP' : ',',
                           'CH.GDPSLOT' : ',',
                           'RX.STN' : '',
                           'RX.AZIMUTH' : ',',
                           'LINE.NAME' : ',',
                           'LINE.NUMBER' : ',',
                           'LINE.DIRECTION' : ',',
                           'LINE.SPREAD' : ',',
                           'JOB.NAME' : ',',
                           'JOB.FOR' : ',',
                           'JOB.BY' : ',',
                           'JOB.NUMBER' : ',',
                           'GDP.File' : ',',
                           'GDP.SN' : ',',
                           'GDP.TCARDSN' : ',',
                           'GDP.NUMCARD' : ',',
                           'GDP.ADCARDSN' : ',',
                           'GDP.ADCARDSND' : ',',
                           'GDP.CARDTYPE' : ',',
                           'GDP.BAT' : ',',
                           'GDP.TEMP' : ',',
                           'GDP.HUMID' : ',',
                           'TS.NCYCLE' : ',',
                           'TS.NWAVEFORM' : ',',
                           'TS.DECFAC' : ',',
                           'TX.SN,NONE' : ',',
                           'TX.STN' : ',',
                           'TX.FREQ' : ',',
                           'TX.DUTY' : ',',
                           'TX.AMP' : ',',
                           'TX.SHUNT' : ','}
    
    #==================================================
    def check_sampling_rate(self, zt_list):
        """
        check to make sure the sampling rate is the same for all channels
        
        Arguments:
        -----------
            **zt_list** : list of Zen3D instances
        
        Outputs:
        --------
            **None** : raises an error if sampling rates are not all the same
            
        """
        
        nz = len(zt_list)
        
        df_list = np.zeros(nz)
        for ii, zt in enumerate(zt_list):
            df_list[ii] = zt.df
            
        tf_array = np.zeros((nz, nz))
        
        for jj in range(nz):
            tf_array[jj] = np.in1d(df_list, [df_list[jj]])
        
        false_test = np.where(tf_array==False)
        
        if len(false_test[0]) != 0:
            raise IOError('Sampling rates are not the same for all channels '+\
                          'Check file(s)'+zt_list[false_test[0]])
        
    #==================================================
    def check_time_series(self, zt_list, decimate=1):
        """
        check to make sure timeseries line up with eachother.
        
        """
        
        n_fn = len(zt_list)
        
        #test start time
        #st_list = np.array([int(zt.date_time[0][-2:]) for zt in zt_list])
        st_list = np.array([int(zt.gps_time[0]) for zt in zt_list])
        time_max = max(st_list)
        #time_max = np.where(st_list==st_list.max())[0]
        
        #get the number of seconds each time series is off by
        skip_dict = {}
        for ii, zt in enumerate(list(zt_list)):
            try:
                skip_dict[ii] = np.where(zt.gps_time==time_max)[0][0]
            except IndexError:
                zt_list.remove(zt_list[ii])
                print '***SKIPPING {0} '.format(zt.fn)
                print '   because it does not contain correct gps time'
                print '   {0} --> {1}'.format(time_max, 
                                             zt.get_date_time(zt.gps_week, 
                                                             time_max))    
        
        #change data by amount needed        
        for ii, zt in zip(skip_dict.keys(), zt_list):
            if skip_dict[ii] != 0:
                skip_points = skip_dict[ii]*zt.df
                print 'Skipping {0} points for {1}'.format(skip_points,
                                                            zt.ch_cmp)
                zt.time_series = zt.time_series[skip_points:]
                zt.gps_diff = zt.gps_diff[skip_dict[ii]:]
                zt.gps_list = zt.gps_list[skip_dict[ii]:]
                zt.date_time = zt.date_time[skip_dict[ii]:]
                zt.gps_time = zt.gps_time[skip_dict[ii]:]
            
        #test length of time series
        ts_len_list = np.array([len(zt.time_series) for zt in zt_list])
        
        #get the smallest number of points in the time series
        ts_min = ts_len_list.min()
        
        #make a time series array for easy access
        ts_min /= decimate
            
        ts_array = np.zeros((ts_min, n_fn))
        
        #trim the time series if needed
        for ii, zt in enumerate(zt_list):
            if decimate > 1:
                zt.time_series = sps.resample(zt.time_series, 
                                              zt.time_series.shape[0]/decimate,
                                              window='hanning')
            if len(zt.time_series) != ts_min:
                ts_trim = zt.time_series[:ts_min]
            else:
                ts_trim = zt.time_series
            zt.time_series = ts_trim
            
            ts_array[:, ii] = ts_trim
            
            if self.verbose:
                print 'TS length for channel {0} '.format(zt.ch_number)+\
                      '({0}) '.format(zt.ch_cmp)+\
                      '= {0}'.format(len(ts_trim))
                print '    T0 = {0}\n'.format(zt.date_time[0])
            self.log_lines.append(' '*4+\
                                  'TS length for channel {0} '.format(zt.ch_number)+\
                                  '({0}) '.format(zt.ch_cmp)+\
                                  '= {0}'.format(len(ts_trim)))
            self.log_lines.append(', T0 = {0}\n'.format(zt.date_time[0]))
            
        if decimate is not 1:
            ts_array = sps.resample(ts_array, ts_min/decimate, 
                                    window='hanning')
            ts_min = ts_array.shape[0]
            
        
        return ts_array, ts_min
    
    #==================================================    
    def write_cache_file(self, fn_list, save_fn, station='ZEN', decimate=1):
        """
        write a cache file from given filenames
        
        """
        #sort the files so they are in order
        fn_sort_list = []
        for cs in self.chn_order:
            for fn in fn_list:
                if cs in fn.lower():
                    fn_sort_list.append(fn)

        fn_list = fn_sort_list
        print fn_list
            
        n_fn = len(fn_list)
        self.zt_list = []
        for fn in fn_list:
            zt1 = Zen3D(fn=fn)
            zt1.verbose = self.verbose
            try:
                zt1.read_3d()
            except ZenGPSError:
                zt1._seconds_diff = 59
                zt1.read_3d()
            self.zt_list.append(zt1)
        
            #fill in meta data from the time series file
            self.meta_data['DATA.DATE0'] = ','+zt1.date_time[0].split(',')[0]
            self.meta_data['DATA.TIME0'] = ','+zt1.date_time[0].split(',')[1]
            self.meta_data['TS.ADFREQ'] = ',{0}'.format(int(zt1.df))
            self.meta_data['CH.FACTOR'] += ','+self._ch_factor 
            self.meta_data['CH.GAIN'] += ','+self._ch_gain
            self.meta_data['CH.CMP'] += ','+zt1.ch_cmp.upper()
            self.meta_data['CH.LENGTH'] += ','+zt1.ch_length
            self.meta_data['CH.EXTGAIN'] += ',1'
            self.meta_data['CH.NOTCH'] += ',NONE'
            self.meta_data['CH.HIGHPASS'] += ',NONE'
            self.meta_data['CH.LOWPASS'] += ','+\
                                       self._ch_lowpass_dict[str(int(zt1.df))]
            self.meta_data['CH.ADCARDSN'] += ','+zt1.ch_adcard_sn
            self.meta_data['CH.NUMBER'] += ',{0}'.format(zt1.ch_number)
            self.meta_data['RX.STN'] += ','+zt1.rx_stn
            
        #make sure all files have the same sampling rate
        self.check_sampling_rate(self.zt_list)
        
        #make sure the length of time series is the same for all channels
        self.ts, ts_len = self.check_time_series(self.zt_list,
                                                 decimate=decimate)
        
        self.meta_data['TS.NPNT'] = ',{0}'.format(ts_len)
        
        #get the file name to save to 
        if save_fn[-4:] == '.cac':
            self.save_fn = save_fn
        elif save_fn[-4] == '.':
            raise ZenInputFileError('File extension needs to be .cac, not'+\
                                    save_fn[-4:])
        else:
            general_fn = station+'_'+\
                         self.meta_data['DATA.DATE0'][1:].replace('-','')+\
                         '_'+self.meta_data['DATA.TIME0'][1:].replace(':','')+\
                         '_'+self.meta_data['TS.ADFREQ'][1:]+'.cac'
            
            if os.path.basename(save_fn) != 'Merged':             
                save_fn = os.path.join(save_fn, 'Merged')
                if not os.path.exists(save_fn):
                    os.mkdir(save_fn)
            self.save_fn = os.path.join(save_fn, general_fn)
                
                
            
        cfid = file(self.save_fn, 'wb+')
        #--> write navigation records first        
        cfid.write(struct.pack('<i', self._nav_len))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['nav']))
        for nd in range(self._nav_len-2):
            cfid.write(struct.pack('<b', 0))
        cfid.write(struct.pack('<i', self._nav_len))
        
        #--> write meta data
        meta_str = ''.join([key+self.meta_data[key]+'\n' 
                             for key in np.sort(self.meta_data.keys())])
        
        meta_len = len(meta_str)
        
        cfid.write(struct.pack('<i', meta_len+2))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['meta']))
        cfid.write(meta_str)
        cfid.write(struct.pack('<i', meta_len+2))
        
        #--> write calibrations
        cal_data1 = 'HEADER.TYPE,Calibrate\nCAL.VER,019\nCAL.SYS,0000,'+\
                   ''.join([' 0.000000: '+'0.000000      0.000000,'*3]*27)
        cal_data2 = '\nCAL.SYS,0000,'+\
                    ''.join([' 0.000000: '+'0.000000      0.000000,'*3]*27)
                    
        cal_data = cal_data1+(cal_data2*(n_fn-1))
        cal_len = len(cal_data)
        
        cfid.write(struct.pack('<i', cal_len+2))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['cal']))
        cfid.write(cal_data[:-1]+'\n')
        cfid.write(struct.pack('<i', cal_len+2))
        
        #--> write data
        
        ts_block_len = int(ts_len)*n_fn*4+2
        
        #--> Need to scale the time series into counts cause that is apparently
        #    what MTFT24 expects
        self.ts = self.ts.astype(np.int32)
        
        #--> make sure none of the data is above the allowed level
        self.ts[np.where(self.ts>2.14e9)] = 2.14e9
        self.ts[np.where(self.ts<-2.14e9)] = -2.14e9
        
        #--> write time series block
        cfid.write(struct.pack('<i', ts_block_len))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['ts']))
        
        #--> need to pack the data as signed integers
        for zz in range(ts_len):
            cfid.write(struct.pack('<'+'i'*n_fn, *self.ts[zz]))
                                
        cfid.write(struct.pack('<i', ts_block_len))
         
        
        cfid.close()
        
        if self.verbose:
            print 'Saved File to: ', self.save_fn
        self.log_lines.append('='*72+'\n')
        self.log_lines.append('Saved File to: \n')
        self.log_lines.append(' '*4+'{0}\n'.format(self.save_fn))
        self.log_lines.append('='*72+'\n')
    
    #==================================================    
    def rewrite_cache_file(self):
        """
        rewrite a cache file if parameters changed
        
        assuming data that was read in is in counts.
        
        """
        self.save_fn_rw = mtfh.make_unique_filename(self.save_fn)
        
        cfid = file(self.save_fn_rw, 'wb+')
        
        n_fn = self.ts.shape[1]
        
        #--> write navigation records first        
        cfid.write(struct.pack('<i', self._nav_len))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['nav']))
        for nd in range(self._nav_len-2):
            cfid.write(struct.pack('<b', 0))
        cfid.write(struct.pack('<i', self._nav_len))
        
        #--> write meta data
        meta_str = ''.join([key+','+','.join(self.meta_data[key])+'\n' 
                             for key in np.sort(self.meta_data.keys())
                             if key != ''])
        
        meta_len = len(meta_str)
        
        cfid.write(struct.pack('<i', meta_len+2))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['meta']))
        cfid.write(meta_str)
        cfid.write(struct.pack('<i', meta_len+2))
        
        #--> write calibrations
        cal_data1 = 'HEADER.TYPE,Calibrate\nCAL.VER,019\nCAL.SYS,0000,'+\
                   ''.join([' 0.000000: '+'0.000000      0.000000,'*3]*1)
        cal_data2 = '\nCAL.SYS,0000,'+\
                    ''.join([' 0.000000: '+'0.000000      0.000000,'*3]*1)
                    
        cal_data = cal_data1+(cal_data2*(self.ts.shape[1]-1))
        cal_len = len(cal_data)
        
        cfid.write(struct.pack('<i', cal_len+2))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['cal']))
        cfid.write(cal_data[:-1]+'\n')
        cfid.write(struct.pack('<i', cal_len+2))
        
        #--> write data
        ts_block_len = self.ts.shape[0]*n_fn*4+2
        
        #--> make sure none of the data is above the allowed level
        self.ts[np.where(self.ts>2.14e9)] = 2.14e9
        self.ts[np.where(self.ts<-2.14e9)] = -2.14e9
        
        #--> write time series block
        cfid.write(struct.pack('<i', ts_block_len))
        cfid.write(struct.pack('<i', self._flag))
        cfid.write(struct.pack('<h', self._type_dict['ts']))
        for zz in range(self.ts.shape[0]):
            cfid.write(struct.pack('<'+'i'*n_fn, *self.ts[zz]))
                                
        cfid.write(struct.pack('<i', ts_block_len))
                 
        cfid.close()
        
        print 'Rewrote {0}\n to {1}'.format(self.save_fn, self.save_fn_rw)        
    
    #==================================================    
    def read_cache_metadata(self, cache_fn):
        """
        read only the meta data from the cache file
        """
        
        self.save_fn = cache_fn
        #open cache file to read in as a binary file
        cfid = file(cache_fn, 'rb')
        
        #read into a long string
        cdata = cfid.read(1050)
        
        #--> read navigation data
        nav_block = np.fromstring(cdata[0:self._stamp_len], 
                                  dtype=self._data_type)
        
        #get starting and ending indices for navigation block
        ii = int(self._stamp_len)
        jj = self._stamp_len+nav_block['len']-2
        self.nav_data = np.fromstring(cdata[ii:jj], dtype=np.int8)
        
        #get indicies for length of block
        ii = int(jj)
        jj = ii+4
        nav_len_check = np.fromstring(cdata[ii:jj], np.int32)
        if nav_len_check != nav_block['len']:
            if self.verbose:
                print 'Index for second navigation length is {0}'.format(ii)
            raise CacheNavigationError('Navigation length in data block are'
                                       'not equal: {0} != {1}'.format(
                                       nav_block['len'], nav_len_check))
        
        #--> read meta data
        ii = int(jj)
        jj = ii+self._stamp_len
        
        meta_block = np.fromstring(cdata[ii:jj], dtype=self._data_type)
        ii = int(jj)
        jj = ii+meta_block['len']-2
        self.meta_data = {}
        meta_list = cdata[ii:jj].split('\n')
        for mm in meta_list:
            mfind = mm.find(',')
            self.meta_data[mm[0:mfind]] = [ms.strip() for ms in 
                                            mm[mfind+1:].split(',')]
        
        #get index for second length test
        ii = int(jj)
        jj = ii+4
        meta_len_check = np.fromstring(cdata[ii:jj], dtype=np.int32)
        if meta_len_check != meta_block['len']:
            if self.verbose:
                print 'Index for second meta length is {0}'.format(ii)
            raise CacheMetaDataError('Meta length in data blocks are not '
                                     'equal: {0} != {1}'.format(
                                     meta_block['len'], meta_len_check))
        cfid.close()
        
    #==================================================
    def read_cache(self, cache_fn):
        """
        read a cache file
        
        """
        
        self.save_fn = cache_fn
        #open cache file to read in as a binary file
        cfid = file(cache_fn, 'rb')
        
        #read into a long string
        cdata = cfid.read()
        
        #--> read navigation data
        nav_block = np.fromstring(cdata[0:self._stamp_len], 
                                  dtype=self._data_type)
        
        #get starting and ending indices for navigation block
        ii = int(self._stamp_len)
        jj = self._stamp_len+nav_block['len']-2
        self.nav_data = np.fromstring(cdata[ii:jj], dtype=np.int8)
        
        #get indicies for length of block
        ii = int(jj)
        jj = ii+4
        nav_len_check = np.fromstring(cdata[ii:jj], np.int32)
        if nav_len_check != nav_block['len']:
            if self.verbose:
                print 'Index for second navigation length is {0}'.format(ii)
            raise CacheNavigationError('Navigation length in data block are'
                                       'not equal: {0} != {1}'.format(
                                       nav_block['len'], nav_len_check))
        
        #--> read meta data
        ii = int(jj)
        jj = ii+self._stamp_len
        
        meta_block = np.fromstring(cdata[ii:jj], dtype=self._data_type)
        ii = int(jj)
        jj = ii+meta_block['len']-2
        self.meta_data = {}
        meta_list = cdata[ii:jj].split('\n')
        for mm in meta_list:
            mfind = mm.find(',')
            self.meta_data[mm[0:mfind]] = mm[mfind+1:].split(',')
        
        #get index for second length test
        ii = int(jj)
        jj = ii+4
        meta_len_check = np.fromstring(cdata[ii:jj], dtype=np.int32)
        if meta_len_check != meta_block['len']:
            if self.verbose:
                print 'Index for second meta length is {0}'.format(ii)
            raise CacheMetaDataError('Meta length in data blocks are not'
                                     'equal: {0} != {1}'.format(
                                     meta_block['len'], meta_len_check))
        
        #--> read calibrations
        ii = int(jj)
        jj = ii+self._stamp_len
        cal_block = np.fromstring(cdata[ii:jj], dtype=self._data_type)
        
        ii = int(jj)
        jj = ii+cal_block['len']-2
        self.cal_data = cdata[ii:jj]
                
        
        ii = int(jj)
        jj = ii+4
        cal_len_check = np.fromstring(cdata[ii:jj], dtype=np.int32)
        if cal_len_check != cal_block['len']:
            if self.verbose:
                print 'Index for second cal length is {0}'.format(ii)
            raise CacheCalibrationError('Cal length in data blocks are not'
                                        'equal: {0} != {1}'.format(
                                        cal_block['len'], cal_len_check))
        
        #--> read data
        ii = int(jj)
        jj = ii+self._stamp_len
        
        ts_block = np.fromstring(cdata[ii:jj], dtype=self._data_type)
        
        #get time series data
        ii = int(jj)
        jj = ii+ts_block['len']-2
        self.ts = np.fromstring(cdata[ii:jj], dtype=self._ts_dtype)
        #resize time series to be length of each channel
        num_chn = len(self.meta_data['ch.cmp'.upper()])
        if self.ts.shape[0]%num_chn != 0:
            print 'Trimming TS by {0} points'.format(self.ts.shape[0]%num_chn)
        self.ts = np.resize(self.ts, (int(self.ts.shape[0]/num_chn), num_chn))
        
        ii = int(jj)
        jj = ii+4
        ts_len_check = np.fromstring(cdata[ii:jj], dtype=np.int32)
        if ts_len_check != ts_block['len']:
            if self.verbose:
                print 'Index for second ts length is {0}'.format(ii)
            raise CacheTimeSeriesError('ts length in data blocks are not'
                                       'equal: {0} != {1}'.format(
                                       ts_block['len'], ts_len_check))

#==============================================================================
# read and write a zen schedule 
#==============================================================================
class ZenSchedule(object):
    """
    deals with reading, writing and copying schedule

    Creates a repeating schedule based on the master_schedule.  It will
    then change the first scheduling action to coincide with the master 
    schedule, such that all deployed boxes will have the same schedule.
    
    :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> zs = zen.ZenSchedule()
        >>> zs.write_schedule('MT01', dt_offset='2013-06-23,04:00:00')
        
    ====================== ====================================================
    Attributes              Description    
    ====================== ====================================================
    ch_cmp_dict            dictionary for channel components with keys being
                           the channel number and values being the channel
                           label
    ch_num_dict            dictionary for channel components whith keys 
                           being channel label and values being channel number
    df_list                 sequential list of sampling rates to repeat in 
                           schedule
    df_time_list            sequential list of time intervals to measure for
                           each corresponding sampling rate
    dt_format              date and time format. *default* is 
                           YYY-MM-DD,hh:mm:ss
    dt_offset              start date and time of schedule in dt_format
    gain_dict              dictionary of gain values for channel number 
    initial_dt             initial date, or dummy zero date for scheduling
    light_dict             dictionary of light color values for schedule
    master_schedule        the schedule that all data loggers should schedule
                           at.  Will taylor the schedule to match the master
                           schedule according to dt_offset
    meta_dict              dictionary for meta data
    meta_keys              keys for meta data dictionary
    sa_keys                keys for schedule actions
    sa_list                 list of schedule actions including time and df
    sr_dict                dictionary of sampling rate values
    verbose                [ True | False ] True to print information to 
                           console 
    ====================== ====================================================

    
    """
    
    def __init__(self):
        
        self.verbose = True
        self.sr_dict = {'256':'0', '512':'1', '1024':'2', '2048':'3', 
                        '4096':'4'}
        self.gain_dict = dict([(mm, 2**mm) for mm in range(7)])
        self.sa_keys = ['date', 'time', 'resync_yn', 'log_yn', 'tx_duty', 
                        'tx_period', 'sr', 'gain', 'nf_yn']
        self.sa_list = []
        self.ch_cmp_dict = {'1':'hx', '2':'hy', '3':'hz', '4':'ex', '5':'ey',
                            '6':'hz'}
        self.ch_num_dict = dict([(self.ch_cmp_dict[key], key) 
                                 for key in self.ch_cmp_dict])
        
        self.meta_keys = ['TX.ID', 'RX.STN', 'Ch.Cmp', 'Ch.Number', 
                          'Ch.varAsp']
        self.meta_dict = {'TX.ID':'none', 'RX.STN':'01', 'Ch.Cmp':'HX',
                          'Ch.Number':'1', 'Ch.varAsp':50}
        self.light_dict = {'YellowLight':0,
                          'BlueLight':1,
                          'RedLight':0,
                          'GreenLight':1}
                          
        self.dt_format = datetime_fmt 
        self.initial_dt = '2000-01-01,00:00:00'
        self.dt_offset = time.strftime(datetime_fmt ,time.gmtime())
        self.df_list = (4096, 1024, 256)
        self.df_time_list = ('00:05:00','00:15:00','05:40:00')
        self.master_schedule = self.make_schedule(self.df_list, 
                                                  self.df_time_list,
                                                  repeat=21)
                                                  
    #==================================================
    def read_schedule(self, fn):
        """
        read zen schedule file
        
        """
        
        sfid = file(fn, 'r')
        lines = sfid.readlines()
        
        for line in lines:
            if line.find('scheduleaction') == 0:
                line_list = line.strip().split(' ')[1].split(',')
                sa_dict = {}
                for ii, key in enumerate(self.sa_keys):
                    sa_dict[key] = line_list[ii]
                self.sa_list.append(sa_dict)
                
            elif line.find('metadata'.upper()) == 0:
                line_list = line.strip().split(' ')[1].split('|')
                for md in line_list[:-1]:
                    md_list = md.strip().split(',')
                    self.meta_dict[md_list[0]] = md_list[1]
                    
            elif line.find('offset') == 0:
                line_str = line.strip().split(' ')
                self.offset = line_str[1]
                
            elif line.find('Light') > 0:
                line_list = line.strip().split(' ')
                try:
                    self.light_dict[line_list[0]]
                    self.light_dict[line_list[0]] = line_list[1]
                except KeyError:
                    pass
    
    #==================================================            
    def add_time(self, date_time, add_minutes=0, add_seconds=0, add_hours=0,
                 add_days=0):
        """
        add time to a time string
        
        assuming date_time is in the format  YYYY-MM-DD,HH:MM:SS
        
        """    
        
        fulldate = datetime.datetime.strptime(date_time, self.dt_format)
                                     
        fulldate = fulldate + datetime.timedelta(days=add_days,
                                                 hours=add_hours,
                                                 minutes=add_minutes,
                                                 seconds=add_seconds)
        return fulldate
    
    #==================================================
    def make_schedule(self, df_list, df_length_list, repeat=5, t1_dict=None):
        """
        make a repeated schedule given list of sampling frequencies and
        duration for each.
        
        Arguments:
        -----------
            **df_list** : list   
                         list of sampling frequencies in Hz, note needs to be
                         powers of 2 starting at 256
            **df_length_list** : list
                                list of durations in hh:mm:ss format
            **repeat** : int
                         number of times to repeat the sequence
            
            **t1_dict** : dictionary
                          dictionary returned from get_schedule_offset
                          
        Returns:
        --------
            **time_list**: list of dictionaries with keys:
                            * 'dt' --> date and time of schedule event
                            * 'df' --> sampling rate for that event
        """
    
        df_list = np.array(df_list)
        df_length_list = np.array(df_length_list)
        ndf = len(df_list)
        
    
        if t1_dict is not None:
            time_list = [{'dt':self.initial_dt,'df':t1_dict['df']}]
    
            kk = np.where(np.array(df_list)==t1_dict['df'])[0]-ndf+1
            df_list = np.append(df_list[kk:], df_list[:kk])
            df_length_list = np.append(df_length_list[kk:], df_length_list[:kk])
            time_list.append(dict([('dt',t1_dict['dt']), ('df',df_list[0])]))
            ii = 1
        else:
            time_list = [{'dt':self.initial_dt,'df':df_list[0]}]
            ii = 0
            
        for rr in range(1,repeat+1):
            for df, df_length, jj in zip(df_list, df_length_list, range(ndf)):
                dtime = time.strptime(df_length, '%H:%M:%S')
                ndt = self.add_time(time_list[ii]['dt'], 
                                    add_hours=dtime.tm_hour,
                                    add_minutes=dtime.tm_min,
                                    add_seconds=dtime.tm_sec)
                time_list.append({'dt':ndt.strftime(self.dt_format),
                                 'df':df_list[jj-ndf+1]})
                ii += 1
                
        for nn, ns in enumerate(time_list):
            sdate, stime = ns['dt'].split(',')
            ns['date'] = sdate
            ns['time'] = stime
            ns['log_yn'] = 'Y'
            ns['nf_yn'] = 'Y'
            ns['sr'] = self.sr_dict[str(ns['df'])]
            ns['tx_duty'] = '0'
            ns['tx_period'] = '0'
            ns['resync_yn'] = 'Y'
            ns['gain'] = '0'
            
        return time_list
    
    #==================================================    
    def get_schedule_offset(self, time_offset, schedule_time_list):
        """
        gets the offset in time from master schedule list and time_offset so 
        that all schedules will record at the same time according to master
        schedule list schedule_time_list
        
        Attributes:
        -----------
            **time_offset** : hh:mm:ss
                              the time offset given to the zen reciever
                              
            **schedule_time_list** : list
                                    list of actual schedule times returned 
                                    from make_schedule
                                    
        Returns:
        --------
            **s1** : dictionary
                     dictionary with keys:
                         * 'dt' --> date and time of offset from next schedule
                                    event from schedule_time_list
                         * 'df' --> sampling rate of that event
        """
        
        dt_offset = '{0},{1}'.format('2000-01-01', time_offset)
        t0 = time.mktime(time.strptime('2000-01-01,00:00:00', self.dt_format))
        
        for ii, tt in enumerate(schedule_time_list):
            ssec = time.mktime(time.strptime(tt['dt'], self.dt_format))
            osec = time.mktime(time.strptime(dt_offset, self.dt_format))
            
            if ssec > osec:
                sdiff = time.localtime(t0+(ssec-osec))
                t1 = self.add_time('2000-01-01,00:00:00', 
                                   add_hours=sdiff.tm_hour,
                                   add_minutes=sdiff.tm_min,
                                   add_seconds=sdiff.tm_sec)
                s1 = {'dt':t1.strftime(self.dt_format), 
                      'df':schedule_time_list[ii-1]['df']}
                return s1
    
    #==================================================            
    def write_schedule(self, station, clear_schedule=True, 
                       clear_metadata=True, varaspace=100, 
                       savename=0, dt_offset=None, 
                       df_list=None, 
                       df_time_list=None, 
                       repeat=8, gain=0):
        """
        write a zen schedule file
        
        **Note**: for the older boxes use 'Zeus3Ini.cfg' for the savename
        
        Arguments:
        ----------
            **station** : int
                          station name must be an integer for the Zen, can
                          be changed later
            
            **clear_schedule** : [ True | False ]
                                 write the line clearschedule in .cfg file
                                 
            **clear_metadata** : [ True | False ]
                                 write the line metadata clear in .cfg file
                                 
            **varaspace** : electrode spacing in meters, can be changed later
            
            **savename** : [ 0 | 1 | 2 | string]
                           * 0 --> saves as zenini.cfg
                           * 1 --> saves as Zeus2Ini.cfg
                           * 2 --> saves as ZEN.cfg
                           * string --> saves as the string, note the zen
                                        boxes look for either 0 or 1, so this
                                        option is useless
                                        
            **dt_offset** : YYYY-MM-DD,hh:mm:ss
                            date and time off offset to start the scheduling.
                            if this is none then current time on computer is
                            used. **In UTC Time**
                            
                            **Note**: this will shift the starting point to 
                                      match the master schedule, so that all
                                      stations have the same schedule.
                                      
            **df_list** : list
                         list of sampling rates in Hz
            
            **df_time_list** : list
                              list of time intervals corresponding to df_list
                              in hh:mm:ss format
            **repeat** : int
                         number of time to repeat the cycle of df_list
            
            **gain** : int
                       gain on instrument, 2 raised to this number.
                       
        Returns:
        --------
            * writes .cfg files to any connected SD card according to channel
              number and ch_num_dict
                    
        """
        
        if dt_offset is not None:
            self.dt_offset = dt_offset
        s1_dict = self.get_schedule_offset(self.dt_offset.split(',')[1],
                                           self.master_schedule)

        if df_list is not None:
            self.df_list = df_list
        if df_time_list is not None:
            self.df_time_list = df_time_list
        
        self.master_schedule =  self.make_schedule(self.df_list, 
                                                  self.df_time_list,
                                                  repeat=repeat*3)
                                                  
        self.sa_list = self.make_schedule(self.df_list,
                                          self.df_time_list,
                                          t1_dict=s1_dict, repeat=repeat)

        drive_names = get_drive_names()
        self.meta_dict['RX.STN'] = station
        self.meta_dict['Ch.varAsp'] = '{0}'.format(varaspace)
        
        if savename == 0:
            save_name = 'zenini.cfg'
        elif savename == 1:
            save_name = 'Zeus3Ini.cfg'
        elif savename == 2:
            save_name = 'ZEN.cfg'
            sfid = file(os.path.normpath(os.path.join('c:\\MT', save_name)),
                        'w')
            for sa_dict in self.sa_list:
                new_time = self.add_time(self.dt_offset,
                                         add_hours=int(sa_dict['time'][0:2]),
                                         add_minutes=int(sa_dict['time'][3:5]),
                                         add_seconds=int(sa_dict['time'][6:]))
                sa_line = ','.join([new_time.strftime(self.dt_format),
                                    sa_dict['resync_yn'], 
                                    sa_dict['log_yn'],
                                    '2047',
                                    '1999999999',
                                    sa_dict['sr'],
                                    '0','0','0','y','n','n','n'])
                sfid.write('scheduleaction '.upper()+sa_line[:-1]+'\n')
            meta_line = ''.join(['{0},{1}|'.format(key,self.meta_dict[key]) 
                                 for key in self.meta_keys])
            sfid.write('METADATA '+meta_line+'\n')
            for lkey in self.light_dict.keys():
                sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
            sfid.close()
            #print 'Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
            #                                       self.ch_cmp_dict[dname[-1]])
            
            for dd in drive_names.keys():
                dname = drive_names[dd]
                sfid = file(os.path.normpath(os.path.join(dd+':\\', save_name)),
                            'w')
                for sa_dict in self.sa_list:
                    new_time = self.add_time(self.dt_offset,
                                             add_hours=int(sa_dict['time'][0:2]),
                                             add_minutes=int(sa_dict['time'][3:5]),
                                             add_seconds=int(sa_dict['time'][6:]))
                    sa_line = ','.join([new_time.strftime(self.dt_format),
                                        sa_dict['resync_yn'], 
                                        sa_dict['log_yn'],
                                        '2047',
                                        '1999999999',
                                        sa_dict['sr'],
                                        '0','0','0','y','n','n','n'])
                    sfid.write('scheduleaction '.upper()+sa_line[:-1]+'\n')
                
                self.meta_dict['Ch.Cmp'] = self.ch_cmp_dict[dname[-1]]
                self.meta_dict['Ch.Number'] = dname[-1]
                meta_line = ''.join(['{0},{1}|'.format(key,self.meta_dict[key]) 
                                     for key in self.meta_keys])
                sfid.write('METADATA '+meta_line+'\n')
                for lkey in self.light_dict.keys():
                    sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
                sfid.close()
                print 'Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
                                                   self.ch_cmp_dict[dname[-1]])
            return
        else:
            save_name = savename
         
        for dd in drive_names.keys():
            dname = drive_names[dd]
            sfid = file(os.path.normpath(os.path.join(dd+':\\', save_name)),
                        'w')
            if clear_schedule:
                sfid.write('clearschedule\n')
            if clear_metadata:
                sfid.write('metadata clear\n')
            for sa_dict in self.sa_list:
                if gain != 0:
                    sa_dict['gain'] = gain
                sa_line = ''.join([sa_dict[key]+',' for key in self.sa_keys])
                sfid.write('scheduleaction '+sa_line[:-1]+'\n')
            sfid.write('offsetschedule {0}\n'.format(self.dt_offset))
            
            self.meta_dict['Ch.Cmp'] = self.ch_cmp_dict[dname[-1]]
            self.meta_dict['Ch.Number'] = dname[-1]
            meta_line = ''.join(['{0},{1}|'.format(key,self.meta_dict[key]) 
                                 for key in self.meta_keys])
            sfid.write('METADATA '+meta_line+'\n')
            for lkey in self.light_dict.keys():
                sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
            sfid.close()
            print 'Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
                                                   self.ch_cmp_dict[dname[-1]])
                                                   
    def write_schedule_for_gui(self, zen_start=None, df_list=None, 
                               df_time_list=None, repeat=8, gain=0,
                               save_path=None, 
                               schedule_fn='zen_schedule.MTsch'):
        
        """
        write a zen schedule file
        
        **Note**: for the older boxes use 'Zeus3Ini.cfg' for the savename
        
        Arguments:
        ----------
                                        
            **zen_start** : hh:mm:ss
                            start time you want the zen to start collecting
                            data. 
                            if this is none then current time on computer is
                            used. **In UTC Time**
                            
                            **Note**: this will shift the starting point to 
                                      match the master schedule, so that all
                                      stations have the same schedule.
                                      
            **df_list** : list
                         list of sampling rates in Hz
            
            **df_time_list** : list
                              list of time intervals corresponding to df_list
                              in hh:mm:ss format
            **repeat** : int
                         number of time to repeat the cycle of df_list
            
            **gain** : int
                       gain on instrument, 2 raised to this number.
                       
        Returns:
        --------
            * writes a schedule file to input into the ZenAcq Gui
        """
    
        if df_list is not None:
            self.df_list = df_list
            
        if df_time_list is not None:
            self.df_time_list = df_time_list
            
        if save_path is None:
            save_path = os.getcwd()
            
        
        # make a master schedule first        
        self.master_schedule =  self.make_schedule(self.df_list, 
                                                  self.df_time_list,
                                                  repeat=repeat*3)
        # estimate the first off set time
        t_offset_dict = self.get_schedule_offset(zen_start,
                                                 self.master_schedule)
        
        # make the schedule with the offset of the first schedule action                                          
        self.sa_list = self.make_schedule(self.df_list,
                                          self.df_time_list,
                                          t1_dict=t_offset_dict,
                                          repeat=repeat)
        
        # make a list of lines to write to a file for ZenAcq
        zacq_list = []
        for ii, ss in enumerate(self.sa_list[:-1]):
            t0 = self._convert_time_to_seconds(ss['time'])
            t1 = self._convert_time_to_seconds(self.sa_list[ii+1]['time'])
            if ss['date'] != self.sa_list[ii+1]['date']:
                t1 += 24*3600                
            
            t_diff = t1-t0
            zacq_list.append('$schline{0:.0f} = {1:.0f},{2:.0f},{3:.0f}\n'.format(
                                ii+1, 
                                t_diff,
                                int(self.sr_dict[str(ss['df'])]),
                                1))
        
        fn = os.path.join(save_path, schedule_fn)
        fid = file(fn, 'w')
        fid.writelines(zacq_list[0:16])
        fid.close()
        
        print 'Wrote schedule file to {0}'.format(fn)
        print '+--------------------------------------+'
        print '|   SET ZEN START TIME TO: {0}    |'.format(zen_start)
        print '+--------------------------------------+'

    def _convert_time_to_seconds(self, time_string):
        """
        convert a time string given as hh:mm:ss into seconds
        """
        t_list = [float(tt) for tt in time_string.split(':')]
        t_seconds = t_list[0]*3600+t_list[1]*60+t_list[2]
        
        return t_seconds
                                                   
                                                  
#==============================================================================
# interface with birrp                        
#==============================================================================
class ZenBIRRP():
    """
    class to deal with Birrp from Zen outputs
    
    survey file is .cfg file
    
    Need to create a processing file which has information on how to 
    process the data.  See read_processing_file and BIRRP documentation
    for details.
    
    The program will run BIRRP from python and convert the outputs into 
    an .edi file.  
    
    Arguments:
    ------------
        **station_path** : string   
                           full path to station data that will be processed
                           
        **station** : string
                      name of station to be processes.  
                      *default* is os.path.basename(station_path) 
                      
        **birrp_exe** : string
                        full path to BIRRP executable.
                        
        **calibration_path** : string
                               full path to calibration file directory
                               In this directory should be the calibration
                               files for the coils named by the coil number.
                               You need to make these files from the 
                               amtant.cal, basically it into individual files
                               which are seperated by commas (csv) files.
                               ex: Ant2344_cal.csv
                           
        **processing_fn** : string
                            full path to processing file, see BIRRP 
                            documentation and mtpy.zen.read_processing_file 
                            for more details on the structure and key words.
       
       **survey_config_fn** : string
                               full path to survey configuration file.
                               This file contains all the important information
                               on how the data was collected. For more see
                               mtpy.utils.configfile.read_survey_configfile
      
        **df** : float
                 sampling rate in Hz of the data being processed.  

        **rr_path** : string
                      full path to remote reference data.

        **rr_station** : string
                         name of remote reference station
                         *default* is os.path.basename(rr_path)
                    
    ======================== ==================================================
      Attributes               Description    
    ======================== ==================================================
    birrp_config_fn          configuration file written once BIRRP runs for
                             convenience if you want to rember what you did.  
    birrp_exe                full path to the BIRRP executable 
    birrp_dict               dictionary of birrp parameters, *default* is None 
    calibration_path         full path to where calibration files exist
    calibration_list         list of coils numbers used in the measurement
    calibration_dict         dictionary of calibration values with keys
                             as coil numbers and values as calbration values.
    df                       sampling frequency (Hz))
    output_path              path to put BIRRP output files
                             *default* is station_path/BF
    processing_dict          dictionary of porcessing information from 
                             processin_fn  
    processing_fn            full path to processing file.  This contains the
                             the information BIRRP needs to process the 
                             station.  For more details on what key words and
                             values would be useful see BIRRP documentation and
                             mtpy.zen.read_processing_file  
    rr_path                  full path to remote reference station /home/mt/rr
    rr_station               name of remote reference station
    rr_survey_dict           dictionary of survey parameters from 
                             survey_config_fn for remote reference 
    script_file              full path to script file used to process BIRRP
    station                  name of station to process
    station_path             full path to station directory ex. /home/mt/mt01
    survey_config_fn         full path to survey cofiguration file which 
                             contains all the important information about how 
                             the data was collected. For more details see on
                             what key words and values to put in see
                             mtpy.utils.configfile.read_survey_configfile
    survey_dict              dictionary with information about survey 
                             parameters from survey_config_fn
    ======================== ==================================================
    
    ======================== ==================================================
     Methods                  Description
    ======================== ==================================================
    get_birrp_parameters     gets the birrp parameters from processing_fn
    get_calibrations         reads in the files in calibration_path and gets
                             data for coil numbers in calibration_list
    get_survey_parameters    get the survey info from survey_config_fn
    set_remote_reference_path     set the remote refernce station and get 
                             survey information
    get_fn_list               get filenames of data files to process
    run_birrp                writes a script file, run's BIRRP from Python
                             and then converts the outputs of BIRRP to .edi
    write_edi_file           writes and .edi file from the outputs of BIRRP
    write_script_file        writes a script file to control how BIRRP 
                             processes the data    
    ======================== ==================================================
    
    
    :Example: ::
        
        >>> import mtpy.usgs.zen as zen
        >>> zen_bp = zen.ZenBIRRP(r"/home/mt/mt01")
        >>> zen.processing_fn = r"/home/mt/mt01/processing.txt"
        >>> zen.survey_config_fn = r"/home/mt/survey.cfg"
        >>> zen.df = 256
        >>> zen_bp.birrp_exe = r"/home/bin/birrp.exe"
        >>> zen.calibration_list = ['2234', '2244', '2254']
        >>> zen.calibration_path = r"/home/zonge/ant_calibrations"
        >>> zen.rr_path = r"/home/mt/rr01"
        >>> zen.run_birrp()
        
    """              

    def __init__(self, station_path, **kwargs):
        
        self.station_path = station_path
        self.rr_path = kwargs.pop('rr_path', None)
        self.survey_config_fn = kwargs.pop('survey_config_fn', None)
        self.processing_fn = kwargs.pop('processing_fn', None)
        self.calibration_path = kwargs.pop('calibration_path', 
                                         r"d:\Peacock\MTData\Ant_calibrations")
        self.calibration_list = ['2254', '2264', '2274', '2284', '2294',
                                '2304', '2314', '2324', '2334', '2344']
        self.birrp_dict = kwargs.pop('birrp_dict', None)
        self.station = kwargs.pop('station', 
                        os.path.basename(os.path.dirname(self.station_path)))
        self.rr_station = None
        self.rr_survey_dict = None
        self.df = kwargs.pop('df', 256)
        self.processing_dict = kwargs.pop('processing_dict', None)
        self.survey_dict = kwargs.pop('survey_dict', None)
        self.birrp_exe = r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts_big.exe"
        self.script_file = None
        self.output_path = None
        self.birrp_config_fn = None
        
        self.calibration_dict = {}
        
    def write_processing_fn(self, station_path=None, **kwargs):
        """
        write a processing station file from the data files
        """
        pass 
        
        
    def get_calibrations(self):
        """
        get coil calibrations
        """
        for cal_fn in os.listdir(self.calibration_path):
            for cal_num in self.calibration_list:
                if cal_num in cal_fn:
                    self.calibration_dict[cal_num] = \
                                    os.path.join(self.calibration_path, cal_fn)
                    break
            
        
    def get_birrp_parameters(self, processing_fn=None):
        """
        get parameters to put into birrp from file
        
        """
        if processing_fn is not None:
            self.processing_fn = processing_fn
        if self.processing_fn is None:
            raise IOError('Need to input a processing file')
            
        processing_list = read_processing_fn(self.processing_fn)
        for pdict in processing_list:
            if pdict['station'] == self.station and \
                                        float(pdict['df']) == self.df:
                return pdict

    
    def get_survey_parameters(self, survey_config_fn=None, rr_station=None):
        """
        get survey parameters from file
        
        """
        if survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
        if self.survey_config_fn is None:
            raise IOError('Need to input a survey config file')
            
        survey_dict_list = mtcf.read_survey_configfile(self.survey_config_fn)
        
        try:
            self.survey_dict = survey_dict_list[self.station.upper()]
        except KeyError:
            print 'Did not find station information in {0}'.format(
                                                        self.survey_config_fn)
        
        if self.rr_station is not None:                                                
            try:
                self.rr_survey_dict = survey_dict_list[self.rr_station.upper()]
            except KeyError:
                print 'Did not find remote station information in {0}'.format(
                                                      self.survey_config_fn)
                    
    def set_remote_reference_path(self, rr_station, rr_path=None):
        """
        set remote reference station and find survey information and filenames
        """ 
        self.rr_station = rr_station
        
        if rr_path is not None:
            self.rr_path = rr_path
            return
        
        # look for path if none is given
        if self.rr_station is not self.station:
            rr_path = self.station_path
            kk = 0
            while os.path.basename(rr_path) != self.station and kk < 5:
                rr_path = os.path.dirname(rr_path)
                kk += 1
            self.rr_path = os.path.join(os.path.dirname(rr_path), 
                                        self.rr_station, 'TS')
            if not os.path.exists(self.rr_path):
                raise IOError('Need to input rrpath, could not find it')
        else:
            self.rr_path = self.station_path
              
        self.get_survey_parameters()
                    
                
    def get_fn_list(self, df=None, start_dt=None, end_dt=None, ncomps=5):
        """
        get the file name list to process
        
        """
        if df is not None:
            self.df = df
            
        comp_dict = dict([(cc, ii) 
                          for ii, cc in enumerate(['ex','ey','hz','hx','hy'])])
        rrcomp_dict = dict([(cc, ii) 
                          for ii, cc in enumerate(['hx','hy'])])
                              
        if start_dt is not None:
            start_seconds = time.mktime(time.strptime(start_dt, datetime_fmt))
        else:
            start_seconds = 0
        if end_dt is not None:
            end_seconds = time.mktime(time.strptime(end_dt, datetime_fmt))
        else:
            end_seconds = 10E11
            
        if self.rr_path is None:
            self.rr_path = self.station_path
            
        fn_list = []
        ii = 0
        for fn in os.listdir(self.station_path):
            try:
                if np.remainder(ii, ncomps) == 0:
                    tarr = np.zeros(ncomps, dtype=[('fn','|S100'),
                                                   ('npts',np.int),
                                                   ('start_dt','|S19'),
                                                   ('end_dt','|S19')])
                header_dict = \
                        mtfh.read_ts_header(os.path.join(self.station_path,fn))
                if header_dict['t_min'] >= start_seconds and \
                   header_dict['t_min'] <= end_seconds and \
                   header_dict['samplingrate'] == float(self.df):
                       
                    kk = comp_dict[header_dict['channel'].lower()]
                    tarr[kk]['fn'] = os.path.join(self.station_path,fn)
                    tarr[kk]['npts'] = int(header_dict['nsamples'])
                    ts_start_dt = time.strftime(datetime_fmt.replace(',',' '), 
                                                time.localtime(header_dict['t_min']))
                    tarr[kk]['start_dt'] = ts_start_dt
                    ts_end_seconds = header_dict['t_min']+\
                         float(header_dict['nsamples']/header_dict['samplingrate'])
                    tarr[kk]['end_dt'] = time.strftime(datetime_fmt.replace(',',' '),
                                                time.localtime(ts_end_seconds))
                    
                    
                    ii += 1
                if ii == ncomps:
                    fn_list.append(tarr)
                    ii = 0
            except mtex.MTpyError_ts_data:
                pass
            except mtex.MTpyError_inputarguments:
                pass
        
        #get remote reference time series
        rrfn_list = []
        ii = 0
        for fn in os.listdir(self.rr_path):
            try:
                if np.remainder(ii, 2) == 0:
                    tarr = np.zeros(2, dtype=[('fn','|S100'),
                                               ('npts',np.int),
                                               ('start_dt','|S19'),
                                               ('end_dt','|S19')])
                header_dict = \
                        mtfh.read_ts_header(os.path.join(self.rr_path,fn))
                        
                if header_dict['t_min'] >= start_seconds and \
                   header_dict['t_min'] <= end_seconds and \
                   header_dict['samplingrate'] == float(self.df):
                    
                    try:
                        kk = rrcomp_dict[header_dict['channel'].lower()]
                        tarr[kk]['fn'] = os.path.join(self.rr_path,fn)
                        tarr[kk]['npts'] = int(header_dict['nsamples'])
                        ts_start_dt = time.strftime(datetime_fmt.replace(',',' '), 
                                                    time.localtime(header_dict['t_min']))
                        tarr[kk]['start_dt'] = ts_start_dt
                        ts_end_seconds = header_dict['t_min']+\
                             float(header_dict['nsamples']/header_dict['samplingrate'])
                        tarr[kk]['end_dt'] = time.strftime(datetime_fmt.replace(',',' '),
                                                    time.localtime(ts_end_seconds))
                        
                        
                        ii += 1
                    except KeyError:
                        pass
                if ii == 2:
                    rrfn_list.append(tarr)
                    ii = 0
            except mtex.MTpyError_ts_data:
                print 'MTpyError_ts_data'
            except mtex.MTpyError_inputarguments:
                print 'MTpyError_inputarguments'
            
        if len(fn_list) > 3:
            fn_list = fn_list[0:3]
            
        if len(rrfn_list) > 3:
            rrfn_list = rrfn_list[0:3]
        
        return fn_list, rrfn_list
                                    
    def write_script_file(self, df=None, processing_fn=None,
                          processing_dict=None, start_dt=None, end_dt=None,
                          ncomps=5, jmode=0, survey_config_fn=None):
        """
        write a script file to guide birrp
        
        """       
        self.get_calibrations()
        
        if df is not None:
            self.df = df
        
        #--> get survey parameters
        if survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
            self.get_survey_parameters()
        elif self.survey_dict is None:
            self.get_survey_parameters()
            
        #--> get processing dictionary 
        if processing_fn is not None:
            self.processing_fn = processing_fn
        
        self.processing_dict = self.get_birrp_parameters()
        
        if processing_dict is not None:
            self.processing_dict = processing_dict
            
        if self.processing_fn is None and self.processing_dict is None:
            raise IOError('Need to input a processing file')
        
        #--> set jmode (how files are read in) as points
        try:
            self.processing_dict['jmode']
        except KeyError:
            self.processing_dict['jmode'] = jmode
            
        #make sure that deltat is set to sampling rate    
        self.processing_dict['deltat'] = -self.df
        
        #get start and end date and time if available
        try:
            start_dt = self.processing_dict['start_dt']
        except KeyError:
            pass
        try:
            end_dt = self.processing_dict['stop_dt']
        except KeyError:
            pass
        
        try:
            self.set_remote_reference_path(self.processing_dict['rrstation'])
        except KeyError:
            self.set_remote_reference_path(self.station)
        
        #get list of files to process from the station folder
        fn_list, rrfn_list = self.get_fn_list(self.df, 
                                           start_dt=start_dt, 
                                           end_dt=end_dt, 
                                           ncomps=ncomps)
        
        
        self.processing_dict['fn_list'] = [fnlist['fn'] for fnlist in fn_list]
        self.processing_dict['rrfn_list'] = [rrfnlist['fn'] 
                             for rrfnlist in rrfn_list]
        
        #need to skip the header string                         
        try:
            self.processing_dict['nskip']
        except KeyError:
            self.processing_dict['nskip'] = 1
        try:
            self.processing_dict['nskipr']
        except KeyError:
            self.processing_dict['nskipr'] = 1
            
        
        #if jmode == 0 for number of points
        if self.processing_dict['jmode'] == 0:
            self.processing_dict['nread'] = [fnlist['npts'].min() 
                                  for fnlist in fn_list]
        
        #if jmode == 1 for entering start and end times
        elif self.processing_dict['jmode'] == 1:
            self.processing_dict['dstim'] = [fnlist['start_dt'] 
                                             for fnlist in fn_list]
            self.processing_dict['wstim'] = [fnlist['start_dt'] 
                                             for fnlist in fn_list]
            self.processing_dict['wetim'] = [fnlist['end_dt'] 
                                             for fnlist in fn_list]
                                                 
        #get calibration files
        #--> HX                                                 
        try:
            self.processing_dict['hx_cal'] = \
                                self.calibration_dict[self.survey_dict['hx']]
        except KeyError:
            print 'Did not find HX calibration in {0}'.format(
                                                    self.survey_config_fn)
            self.processing_dict['hx_cal'] = self.calibration_dict['2284'] 
            print 'Setting calibration coil number to 2284 as default.'            
            
        #--> HY
        try:
            self.processing_dict['hy_cal'] = \
                                self.calibration_dict[self.survey_dict['hy']]
        except KeyError:
            print 'Did not find HZ calibration in {0}'.format(
                                                    self.survey_config_fn)
            self.processing_dict['hy_cal'] = self.calibration_dict['2284']
            print 'Setting calibration coil number to 2284 as default.' 
            
        #--> HZ
        try:
            self.processing_dict['hz_cal'] = \
                                self.calibration_dict[self.survey_dict['hz']]
        except KeyError:
            print 'Did not find HZ calibration in {0}'.format(
                                                    self.survey_config_fn)
            self.processing_dict['hz_cal'] = self.calibration_dict['2284']
            print 'Setting calibration coil number to 2284 as default.' 
            
        if self.rr_survey_dict is not None:
            try:
                self.processing_dict['rrhx_cal'] = \
                                self.calibration_dict[self.rr_survey_dict['hx']]
            except KeyError:
                print 'Did not find RRHX calibration in {0}'.format(
                                                    self.survey_config_fn)
                self.processing_dict['rrhx_cal'] = \
                                            self.calibration_dict['2284']
                print 'Setting calibration coil number to 2284 as default.' 
                
            try:
                self.processing_dict['rrhy_cal'] = \
                                self.calibration_dict[self.rr_survey_dict['hy']]
            except KeyError:
                print 'Did not find RRHY calibration in {0}'.format(
                                                        self.survey_config_fn)
                self.processing_dict['rrhy_cal'] = \
                                            self.calibration_dict['2284']
                print 'Setting calibration coil number to 2284 as default.' 
        
        #set the save path to include the sampling rate
        self.output_path = os.path.join(os.path.dirname(
                                        self.processing_dict['fn_list'][0][0]), 
                                        'BF_{0}'.format(self.df))
        #write script file using mtpy.processing.birrp    
        script_file, birrp_dict = birrp.write_script_file(dict(self.processing_dict),
                                                    save_path=self.output_path)
        
        cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(
                                                             script_file[:-7]))
                                                             
        mtcf.write_dict_to_configfile(birrp_dict, cfg_fn)
        print 'Wrote BIRRP config file for edi file to {0}'.format(cfg_fn)
        
        self.birrp_config_fn = cfg_fn
        self.script_file = script_file
        self.birrp_dict = birrp_dict
        
    def run_birrp(self, script_file=None, birrp_exe=None):
        """
        run birrp given the specified files
        
        """
        if script_file is not None:
            self.script_file = script_file
        if self.script_file is None:
            self.write_script_file()
        
        if birrp_exe is not None:
            self.birrp_exe = birrp_exe
            
        birrp.run(self.birrp_exe, self.script_file)
        
        self.edi_fn = self.write_edi_file(self.output_path, 
                                          self.survey_config_fn,
                                          self.birrp_config_fn)
        
    def write_edi_file(self, birrp_output_path=None, survey_config_fn=None, 
                       birrp_config_fn=None):
        """
        write an edi file from outputs of birrp
        """
        
        if birrp_output_path is not None and self.output_path is None:
            self.run_birrp()
        elif birrp_output_path is not None:
            self.output_path = birrp_output_path
            
        if survey_config_fn is None and self.survey_config_fn is None:
            self.get_survey_parameters()
        elif survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
        
        if self.birrp_config_fn is None and birrp_config_fn is None:
            self.write_script_file()
        elif birrp_config_fn is not None:
            self.birrp_config_fn = birrp_config_fn
           
            
        edi_fn = birrp.convert2edi(self.station, 
                                   self.output_path, 
                                   self.survey_config_fn, 
                                   self.birrp_config_fn)
        
        return edi_fn
                                                  
#==============================================================================
#  Error instances for Zen
#==============================================================================
class ZenGPSError(Exception):
    """
    error for gps timing
    """
    pass

class ZenSamplingRateError(Exception):
    """
    error for different sampling rates
    """
    pass

class ZenInputFileError(Exception):
    """
    error for input files
    """
    pass

class CacheNavigationError(Exception):
    """
    error for navigation block in cache file
    """
    pass

class CacheMetaDataError(Exception):
    """
    error for meta data block in cache file
    """
    pass

class CacheCalibrationError(Exception):
    """
    error for calibration block in cache file
    """
    pass

class CacheTimeSeriesError(Exception):
    """
    error for time series block in cache file
    """
    pass

#==============================================================================
# make a class to go from Z3d to .edi
#==============================================================================
class BIRRP_processing(object):
    """
    configuration file for birrp processing
    """
    
    def __init__(self, **kwargs):
        self.jmode = 0
        self.nskip = 1
        self.nskipr = 1
        self.calibration_path = kwargs.pop('calibration_path', 
                                         r"d:\Peacock\MTData\Ant_calibrations")
        self.calibration_list = ['2254', '2264', '2274', '2284', '2294',
                                '2304', '2314', '2324', '2334', '2344']
                                
        self.mcomps = 5
        self.elecori = "EX,EY"
        self.tbw = 2
        self.ainuin = .9999
        self.magtype = 'bb'
        self.nfft = 2**18
        self.nsctmax = 14
        self.ilev = 0
        self.nar = 5
        self.nrr = 0
        self.c2thresb = 0.45        
        
    def get_calibrations(self, calibration_path=None):
        """
        get coil calibrations
        """
        if calibration_path is not None:
            self.calibration_path = calibration_path
            
        calibration_dict = {}
        for cal_fn in os.listdir(self.calibration_path):
            for cal_num in self.calibration_list:
                if cal_num in cal_fn:
                    calibration_dict[cal_num] = \
                                    os.path.join(self.calibration_path, cal_fn)
                    break
        return calibration_dict
        
    def get_processing_dict(self, fn_birrp_list, hx=2284, hy=2284, hz=2284):
        """
        from fn_birrp_arr make a processing dictionary to input into writing
        a birrp script file
        
        fn_birrp_list = fn_birrp_arr[df]
        """
        comp_dict = {4:{'EX':0, 'EY':1, 'HX':2, 'HY':3},
                     5:{'EX':0, 'EY':1, 'HZ':2, 'HX':3, 'HY':4}}
        rr_comp_dict = {'HX':0, 'HY':1}
        
        self.fn_list = [fn_list['fn'] for fn_list in fn_birrp_list]
        # need to sort the fn list so that the files are in the correct
        # order for input and output as defined by birrp
        for ii, f_list in enumerate(self.fn_list):
            sort_list = list(f_list)
            num_comps = len(f_list)
            for fn in f_list:
                key = fn[-2:]
                sort_list[comp_dict[num_comps][key.upper()]] = fn
            self.fn_list[ii] = sort_list
            
        # get remote reference file names, same as input, just hx and hy
        self.rrfn_list = []
        for fn_list in fn_birrp_list:
            rr_list = [1, 2]
            for fn in fn_list['fn']:
                key = fn[-2:].upper()
                if key == 'HX' or key == 'HY':
                   rr_list[rr_comp_dict[key]] = fn
            self.rrfn_list.append(rr_list)

        self.nread = [fn_list['npts'].min() for fn_list in fn_birrp_list] 
        self.mcomps = len(fn_birrp_list[0])
        
        if self.mcomps == 5:
            self.magori = "HZ,HX,HY"
            
        elif self.mcomps == 4:
            self.magori = "HX,HY"
        else:
            raise IOError('Number of components is {0}'.format(self.mcomps))
        
        # get calibrations for coil responses        
        cal_dict = self.get_calibrations()
        #get calibration files
        #--> HX                                                 
        try:
            self.hx_cal = cal_dict[str(hx)]
            self.rrhx_cal = cal_dict[str(hx)]
        except KeyError:
            print 'Did not find HX calibration for {0}'.format(hx)
            self.hx_cal = cal_dict['2284'] 
            self.rrhx_cal = cal_dict['2284'] 
            print 'Setting calibration coil number to 2284 as default.'            
        #--> HY                                                 
        try:
            self.hy_cal = cal_dict[str(hy)]
            self.rrhy_cal = cal_dict[str(hy)]
        except KeyError:
            print 'Did not find HX calibration for {0}'.format(hy)
            self.hy_cal = cal_dict['2284'] 
            self.rrhy_cal = cal_dict['2284'] 
            print 'Setting calibration coil number to 2284 as default.'            
        #--> HZ                                                 
        try:
            self.hz_cal = cal_dict[str(hz)]
        except KeyError:
            print 'Did not find HX calibration for {0}'.format(hz)
            self.hz_cal = cal_dict['2284'] 
            print 'Setting calibration coil number to 2284 as default.'            
        
        return self.__dict__
        
        
class Survey_Config(object):
    """
    survey config class
    """
    def __init__(self, **kwargs):
        self.b_instrument_amplification = 1
        self.b_instrument_type = 'coil'
        self.b_logger_gain = 1
        self.b_logger_type = 'zen'
        self.b_xaxis_azimuth = 0
        self.b_yaxis_azimuth = 90
        self.box = 24
        self.date = '01/01/00'
        self.e_instrument_amplification = 1
        self.e_instrument_type = 'Cu-CuSO4 electrodes'
        self.e_logger_gain = 1
        self.e_logger_type = 'zen'
        self.e_xaxis_azimuth = 0
        self.e_xaxis_length = 100
        self.e_yaxis_azimuth = 90
        self.e_yaxis_length = 100
        self.elevation = 2113.2
        self.hx = 2324
        self.hy = 2314
        self.hz = 2334
        self.lat = 37.8861
        self.location = 'Earth'
        self.lon = -119.05417
        self.network = 'USGS'
        self.notes = 'Generic config file'
        self.sampling_interval = 'all'
        self.station = 'mb000'
        self.station_type = 'mt'
        self.save_path = None
        
        for key in kwargs:
            setattr(self, key, kwargs[key])
    
    def write_survey_config_file(self, save_path=None):
        """
        write a survey config file to save path
        """
        
        if save_path is not None:
            self.save_path = save_path
        fn = os.path.join(self.save_path, '{0}.cfg'.format(self.station))
        mtcfg.write_dict_to_configfile({self.station:self.__dict__}, fn)
        
        print 'Wrote survey config file to {0}'.format(fn)
        
        return fn


class Z3D_to_edi(object):
    """
    go from z3d files to .edi
    """
    
    def __init__(self, station_dir=None, **kwargs):
        
        self.station_dir = station_dir
        #ZenBIRRP.__init__(self, self.station_dir)
        self.survey_config = Survey_Config(save_path=self.station_dir)
        self.survey_config_fn = None
        self.birrp_config_fn = None
        self.birrp_exe = r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts.exe"
        self.coil_cal_path = r"c:\MT\Ant_calibrations"
        self.num_comp = 5

        
    def make_survey_config_file(self, survey_config_dict=None):
        """
        make a survey configuration file from the data
        """
        
        self.survey_config_fn = self.survey_config.write_survey_config_file()
        
    def get_schedules_fn_from_dir(self, station_ts_dir):
        """
        get the birrp fn list from a directory of TS files
        """
        
        self.station_dir = station_ts_dir
        
        fn_arr = np.zeros(len(os.listdir(station_ts_dir)),
                          dtype=[('fn','|S100'),
                                 ('npts',np.int),
                                 ('start_dt','|S19'),
                                 ('end_dt','|S19'),
                                 ('df', np.float)])
        fn_count = 0
        for fn in os.listdir(station_ts_dir):
            fn = os.path.join(station_ts_dir, fn)
            try:
                header_dict = mtfh.read_ts_header(fn)
                fn_arr[fn_count]['fn'] = fn
                fn_arr[fn_count]['npts'] = header_dict['nsamples']
                fn_arr[fn_count]['df'] = header_dict['samplingrate']
                start_sec = header_dict['t_min']
                num_sec = float(header_dict['nsamples'])/\
                                                   header_dict['samplingrate']
                fn_arr[fn_count]['start_dt'] = time.strftime(datetime_fmt, 
                                                time.localtime(start_sec)) 
                fn_arr[fn_count]['end_dt'] = time.strftime(datetime_fmt, 
                                                time.localtime(start_sec+\
                                                num_sec))
                fn_count += 1
            except mtex.MTpyError_ts_data:
                print '  Skipped {0}'.format(fn)
            except mtex.MTpyError_inputarguments:
                print '  Skipped {0}'.format(fn)
        
        # be sure to trim the array
        fn_arr = fn_arr[np.nonzero(fn_arr['npts'])]
        
        return self.get_schedules_fn(fn_arr)
            
        
    def get_schedules_fn(self, fn_arr):
        """
        seperate out the different schedule blocks and frequencies so the
        can be processed
        
        Returns
        ---------
            **schedule_fn_dict** : dictionary
                                   keys are sampling rates and values are
                                   lists of file names for each schedule
                                   block up to 3 blocks
        """
        # get the sampling rates used
        s_keys = set(fn_arr['df'])
        
        # make a dictionary with keys as the sampling rates 
        s_dict = dict([(skey, []) for skey in s_keys])
        
        # loop over the sampling rates and find the schedule blocks
        for df in s_keys:
            # find startind dates for sampling rate
            s_dates = set(fn_arr['start_dt'][np.where(fn_arr['df']==df)])
            for sdate in s_dates:
                s_fn_arr = fn_arr[np.where(fn_arr['start_dt']==sdate)]
                s_fn_birrp_arr = np.zeros(len(s_fn_arr), 
                                          dtype=[('fn','|S100'),
                                                 ('npts',np.int),
                                                 ('start_dt','|S19'),
                                                 ('end_dt','|S19')])
                s_fn_birrp_arr['fn'] = s_fn_arr['fn']
                s_fn_birrp_arr['npts'][:] = s_fn_arr['npts'].min()
                s_fn_birrp_arr['start_dt'][:] = sdate
                start_seconds = time.mktime(time.strptime(sdate, 
                                                          datetime_fmt))
 
                end_seconds = start_seconds+s_fn_arr['npts'].min()/float(df)
                s_fn_birrp_arr['end_dt'][:] = time.strftime(datetime_sec,
                                                time.localtime(end_seconds))  
                s_dict[df].append(s_fn_birrp_arr)
        
        return s_dict
        
    def make_mtpy_ascii_files(self, station_dir=None, fmt='%.8', 
                              station_name='mb', notch_dict={},
                              df_list=None, max_blocks=3, ex=100., ey=100.,): 
        """
        makes mtpy_mt files from .Z3D files
        
        Arguments:
        -----------
            **dirpath** : full path to .Z3D files
            
            **station_name** : prefix for station names
            
            **fmt** : format of data numbers for mt_files
            
        Outputs:
        --------
            **fn_arr** : np.ndarray(file, length, df, start_dt)
            
        :Example: ::
        
            >>> import mtpy.usgs.zen as zen
            >>> fn_list = zen.copy_from_sd('mt01')
            >>> mtpy_fn = zen.make_mtpy_files(fn_list, station_name='mt')
        """
        
        if station_dir is not None:
            self.station_dir = station_dir
            
        fn_list = [os.path.join(self.station_dir, fn) 
                    for fn in os.listdir(self.station_dir) 
                    if fn[-4:] == '.Z3D']
        if len(fn_list) == 0:
            raise IOError('Could not find any .Z3D files in {0}'.format(
                            self.station_dir))
                            
        # make an array that has all the information about each file
        fn_arr = np.zeros(len(fn_list), 
                          dtype=[('station','|S6'), 
                                 ('npts',np.int), 
                                 ('df', np.int),
                                 ('start_dt', '|S22'), 
                                 ('comp','|S2'),
                                 ('fn','|S100')])
        fn_lines = []
        z3d_count = 0             
        for ii, fn in enumerate(fn_list):
            if z3d_count > len(df_list)*self.num_comp*max_blocks-1:
                   break
            if df_list is not None:
               zd = Zen3D(fn)
               zd.read_header()
               
               
               if zd.header.ad_rate in df_list:
                   zd.read_z3d()
                   z3d_count += 1

               else:
                   continue
            else:
                zd = Zen3D(fn)
                zd.read_z3d()
            
            if zd.metadata.ch_cmp.lower() == 'hx':
                self.survey_config.hx = zd.metadata.ch_number
            if zd.metadata.ch_cmp.lower() == 'hy':
                self.survey_config.hy = zd.metadata.ch_number
            if zd.metadata.ch_cmp.lower() == 'hz':
                self.survey_config.hz = zd.metadata.ch_number
            if zd.metadata.ch_cmp.lower() == 'ex':
                self.survey_config.e_xaxis_length = zd.metadata.ch_length
            if zd.metadata.ch_cmp.lower() == 'ey':
                self.survey_config.e_yaxis_length = zd.metadata.ch_length

            # get station configuration from the first Z3D file            
            if ii == 0:
                self.survey_config.lat = zd.header.lat
                self.survey_config.lon = zd.header.long
                self.survey_config.date = zd.schedule.Date.replace('-','/')
                self.survey_config.box = int(zd.header.box_number)
            
            #write mtpy mt file
            zd.write_ascii_mt_file(notch_dict=notch_dict, ex=ex, ey=ey)
            
            #create lines to write to a log file                       
            station = zd.metadata.rx_xyz0.split(':')[0]
            fn_arr[ii]['station'] = '{0}{1}'.format(station_name, station)
            fn_arr[ii]['npts'] = zd.time_series.shape[0]
            fn_arr[ii]['df'] = zd.df
            fn_arr[ii]['start_dt'] = zd.zen_schedule
            fn_arr[ii]['comp'] = zd.metadata.ch_cmp.lower()
            fn_arr[ii]['fn'] = zd.fn_mt_ascii
            fn_lines.append(''.join(['--> station: {0}{1}\n'.format(station_name, station),
                                     '    ts_len = {0}\n'.format(zd.time_series.shape[0]),
                                     '    df = {0}\n'.format(zd.df),
                                     '    start_dt = {0}\n'.format(zd.zen_schedule),
                                     '    comp = {0}\n'.format(zd.metadata.ch_cmp),
                                     '    fn = {0}\n'.format(zd.fn)]))
                                     
        self.station_dir = os.path.join(self.station_dir, 'TS')
        self.survey_config.save_path = self.station_dir
        # write survey configuration file
        self.survey_config.write_survey_config_file()
        
            
        return fn_arr[np.nonzero(fn_arr['npts'])], fn_lines
        
    def write_script_files(self, fn_birrp_dict, save_path=None):
        """
        write a script file from a generic processing dictionary
        """
        
        if save_path is None:
            save_path = os.path.join(self.station_dir, 'BF')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
        s_keys = fn_birrp_dict.keys()
        script_fn_list = []
        for skey in s_keys:
            bf_path = os.path.join(save_path, '{0:.0f}'.format(skey))
            fn_birrp_arr = fn_birrp_dict[skey]
            pro_obj = BIRRP_processing()
            pro_obj.calibration_path = self.coil_cal_path
            pro_obj.station = self.survey_config.station
            pro_obj.deltat = -float(skey)
            pro_dict = pro_obj.get_processing_dict(fn_birrp_arr, 
                                                   hx=self.survey_config.hx,
                                                   hy=self.survey_config.hy,
                                                   hz=self.survey_config.hz)
            
            #write script file using mtpy.processing.birrp    
            script_fn, birrp_dict = birrp.write_script_file(pro_dict,
                                                        save_path=bf_path)
                                                        
            script_fn_list.append(script_fn)
            
            cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(
                                                                 script_fn[:-7]))
                                                                 
            mtcfg.write_dict_to_configfile(birrp_dict, cfg_fn)
            print 'Wrote BIRRP config file for edi file to {0}'.format(cfg_fn)
    
            self.birrp_config_fn = cfg_fn
        
        return script_fn_list   
        
    def run_birrp(self, script_fn_list=None, birrp_exe=None):
        """
        run birrp given the specified files
        
        """
        
        if script_fn_list is None:
            raise IOError('Need to input a script file or list of script files')
        
        if birrp_exe is not None:
            self.birrp_exe = birrp_exe
            
            
        if type(script_fn_list) is list:
            self.edi_fn = []
            for script_fn in script_fn_list:
                birrp.run(self.birrp_exe, script_fn)
                
                output_path = os.path.dirname(script_fn)
                self.edi_fn.append(self.write_edi_file(output_path, 
                                      survey_config_fn=self.survey_config_fn,
                                      birrp_config_fn=self.birrp_config_fn))
        elif type(script_fn_list) is str:
            birrp.run(self.birrp_exe, script_fn_list)
            
            output_path = os.path.dirname(script_fn)
            self.edi_fn = self.write_edi_file(output_path, 
                                      survey_config_fn=self.survey_config_fn,
                                      birrp_config_fn=self.birrp_config_fn)
        
    def write_edi_file(self, birrp_output_path, survey_config_fn=None, 
                       birrp_config_fn=None):
        """
        write an edi file from outputs of birrp
        """
        if self.survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn        
        
        if self.survey_config_fn is None:
            ts_find = birrp_output_path.find('TS')
            if ts_find > 0:
                ts_dir = birrp_output_path[0:ts_find+2]
                for fn in os.listdir(ts_dir):
                    if fn[-4:] == '.cfg':
                        self.survey_config_fn = os.path.join(ts_dir, fn)
        
        edi_fn = birrp.convert2edi(self.survey_config.station, 
                                   birrp_output_path, 
                                   self.survey_config_fn, 
                                   self.birrp_config_fn)
        
        return edi_fn
        
    def plot_responses(self, edi_fn_list=None):
        """
        
        plot all the edi files that were created.
        """
        if edi_fn_list is not None:
            self.edi_fn = edi_fn_list
            

        
        if type(self.edi_fn) is list:
            # check file lengths to make sure there are no zero length files
            for edi_fn in self.edi_fn:
                fn_size = os.path.getsize(edi_fn)
                if fn_size < 3000:
                    self.edi_fn.remove(edi_fn)
                if len(self.edi_fn) == 0:
                    raise ValueError('No good .edi files where produced')
            resp_plot = plotnresponses.PlotMultipleResponses(fn_list=self.edi_fn,
                                                         plot_style='compare',
                                                         plot_tipper='yri')
        elif type(self.edi_fn) is str:
            if os.path.getsize(self.edi_fn) < 3000:
                raise ValueError('No good .edi files where produced')
            resp_plot = plotresponse.PlotResponse(fn=self.edi_fn,
                                                  plot_tipper='yri')
                                                         
        return resp_plot
 
    def process_data(self, df_list=None, max_blocks=2, num_comp=5):
        """
        from the input station directory, convert files to ascii, run through
        BIRRP, convert to .edi files and plot
        """
        
        st = time.time()
        self.num_comp = num_comp
        
        if df_list is not None:
            if type(df_list) is float or type(df_list) is int or\
               type(df_list) is str:
               df_list = [df_list] 
        
        # make files into mtpy files
        z3d_fn_list, log_lines = self.make_mtpy_ascii_files(df_list=df_list,
                                                            max_blocks=max_blocks)
        
        # get all information from mtpy files
        schedule_dict = self.get_schedules_fn(z3d_fn_list)
            
        # write script files for birrp
        sfn_list = self.write_script_files(schedule_dict)
        
        # run birrp
        self.run_birrp(sfn_list)
        
        # plot the output
        r_plot = self.plot_responses()
        
        et = time.time()
        print 'took {0} seconds'.format(et-st)
        
        return r_plot


#==============================================================================
# read processing file
#==============================================================================
def read_processing_fn(processing_fn, delimiter='\t'):
    """
    Read in the information from processing file and output
    as a list of dictionaries.
    
    can include:
    
    ================== ========================================================
    parameter          description
    ================== ======================================================== 
    station            station name
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
    ================== ========================================================
    
    ..see also::                
                
        => see BIRRP Manual for more details on the parameters
        => see A. D. Chave and D. J. Thomson [1989,2003,2004] for more
            information on Bounded influence and robust processing.
            
    Arguments:
    -----------
        **processing_fn** : string (full path to file)
                              tab delimited text file with appropriate 
                              information.
                              
        **station_path** : directory path to where station folders are
        
    Outputs:
    ----------
        **slist** : list of dictionaries with key words related to headers of 
                   txt file
                   
    :Example File: ::
    
    station	df	start_dt	stop	rrstation	rrstart	rrstop	mcomps\
    magori	elecori	rrmagori	tbw	ainuin	magtype	nfft	nsctmax\
    ilev	nar	nrr	c2thresb	nsctinc	nf1	nfinc	nfsect	ainlin\
    declination	thetae
    mb037	256	2013-06-28,00:00:00	2013-06-28,18:00:00	mbrr	\
    2013-06-28,00:00:00	2013-06-28,18:00:00	5	"HZ,HX,HY"	"EX,EY"\
    "HX,HY"	2	0.9999	bb	262144	14	1	5	0	0.45	2\
    3	1	3	0.0001	-13.367	0,90,180
    
    """
    
    pfid = open(processing_fn, 'r')
    plines = pfid.readlines()
    pkeys = plines[0].rstrip()
    pkeys = pkeys.split('\t')
    plist=[]
    for pline in plines[1:]:
        pstr = pline.rstrip()
        pstr = pstr.split(delimiter)
        if len(pstr)>1:
            pdict={}
            for kk, pkey in enumerate(pkeys):
                pstr[kk] = pstr[kk].replace('"','')
                pdict[pkey.lower()] = pstr[kk]
        plist.append(pdict)
        
    pfid.close()
    
    return plist
#==============================================================================
# get the external drives for SD cards
#==============================================================================
def get_drives():
    """
    get a list of logical drives detected on the machine
    
    Note this only works for windows.
    
    Outputs:
    ----------
        **drives** : list of drives as letters
        
    :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> zen.get_drives()
    
    """
    drives = []
    bitmask = win32api.GetLogicalDrives()
    for letter in string.uppercase:
        if bitmask & 1:
            drives.append(letter)
        bitmask >>= 1

    return drives
   
#==============================================================================
# get the names of the drives which should correspond to channels
#==============================================================================
def get_drive_names():
    """
    get a list of drive names detected assuming the cards are names by box 
    and channel.
    
    Outputs:
    ----------
        **drive_dict** : dictionary
                         keys are the drive letters and values are the 
                         drive names
    :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> zen.get_drives_names()
    """
    
    drives = get_drives()
    
    drive_dict = {}
    for drive in drives:
        try:
            drive_name = win32api.GetVolumeInformation(drive+':\\')[0]
            if drive_name.find('CH') > 0:
                drive_dict[drive] = drive_name
        except:
            pass
    
    if drives == {}:
        print 'No external drives detected, check the connections.'
        return None
    else:
        return drive_dict

#==============================================================================
# copy files from SD cards   
#==============================================================================
def copy_from_sd(station, save_path=r"d:\Peacock\MTData", 
                 channel_dict={'1':'HX', '2':'HY', '3':'HZ',
                               '4':'EX', '5':'EY', '6':'HZ'},
                 copy_date=None, copy_type='all'):
    """
    copy files from sd cards into a common folder (save_path)
    
    do not put an underscore in station, causes problems at the moment
    
    Arguments:
    -----------
        **station** : string
                      full name of station from which data is being saved
        
        **save_path** : string
                       full path to save data to
                       
        **channel_dict** : dictionary
                           keys are the channel numbers as strings and the
                           values are the component that corresponds to that 
                           channel, values are placed in upper case in the 
                           code
                           
        **copy_date** : YYYY-MM-DD
                        date to copy from depending on copy_type
                        
        **copy_type** : [ 'all' | 'before' | 'after' | 'on' ]
                        * 'all' --> copy all files on the SD card
                        * 'before' --> copy files before and on this date
                        * 'after' --> copy files on and after this date
                        * 'on' --> copy files on this date only
                        
    Outputs:
    -----------
        **fn_list** : list
                     list of filenames copied to save_path
                     
    :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> fn_list = zen.copy_from_sd('mt01', save_path=r"/home/mt/survey_1")
    
    """
    
    drive_names = get_drive_names()
    if drive_names is None:
        raise IOError('No drives to copy from.')
    save_path = os.path.join(save_path,station)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    log_fid = file(os.path.join(save_path,'copy_from_sd.log'),'w')
    
    
    st_test = time.ctime()
    fn_list = []
    for key in drive_names.keys():
        dr = r"{0}:\\".format(key)
        print '='*25+drive_names[key]+'='*25
        log_fid.write('='*25+drive_names[key]+'='*25+'\n')
        for fn in os.listdir(dr):
            full_path_fn = os.path.normpath(os.path.join(dr, fn))
            if fn[-4:] == '.cfg':
                shutil.copy(full_path_fn, os.path.join(save_path, fn))
                    
            try:
                file_size = os.stat(full_path_fn)[6]
                if file_size >= 1600L and fn.find('.cfg') == -1:
                    zt = Zen3D(fn=full_path_fn)
                    #zt.get_info()
                    zt.read_header()
                    zt.read_schedule()
                    zt.read_metadata()
                    schedule_date = '{0}'.format(zt.schedule.Date)
                    
                    if zt.metadata.rx_xyz0.find(station[2:]) >= 0:
                        fn_find = True
                        if copy_date is not None:
                            cp_date = int(''.join(copy_date.split('-')))
                            
                            fn_find = False
                            
                            zt_date = int(''.join(schedule_date.split('-')))
                            if copy_type == 'before':
                                if zt_date <= cp_date:
                                    fn_find = True
                            elif copy_type == 'after':
                                if zt_date >= cp_date:
                                    fn_find = True
                            elif copy_type == 'on':
                                if zt_date == cp_date:
                                    fn_find = True
                                                                
                        if fn_find:
                            channel = zt.metadata.ch_cmp.upper()
                            st = zt.schedule.Time.replace(':','')
                            sd = zt.schedule.Date.replace('-','')
                            sv_fn = '{0}_{1}_{2}_{3}_{4}.Z3D'.format(station, 
                                                                     sd, 
                                                                     st,
                                                                     int(zt.df),
                                                                     channel)
                                                                 
                            full_path_sv = os.path.join(save_path, sv_fn)
                            fn_list.append(full_path_sv)
                            
                            shutil.copy(full_path_fn, full_path_sv)
                            print 'copied {0} to {1}\n'.format(full_path_fn, 
                                                             full_path_sv)
                                                             
                            #log_fid.writelines(zt.log_lines)
                                                             
                            log_fid.write('copied {0} to \n'.format(full_path_fn)+\
                                          '       {0}\n'.format(full_path_sv))
                        else:
                            pass
#                            print '+++ SKIPPED {0}+++\n'.format(zt.fn)
#                            log_fid.write(' '*4+\
#                                          '+++ SKIPPED {0}+++\n'.format(zt.fn))
                        
                    else:
                        pass
#                        print '{0} '.format(full_path_fn)+\
#                               'not copied due to bad data.'
#                               
#                        log_fid.write(' '*4+'***{0} '.format(full_path_fn)+\
#                                      'not copied due to bad data.\n\n')
            except WindowsError:
                print 'Faulty file at {0}'.format(full_path_fn)
                log_fid.write('---Faulty file at {0}\n\n'.format(full_path_fn))
    log_fid.close()
    
    et_test = time.ctime()
    
    print 'Started at: {0}'.format(st_test)
    print 'Ended at: {0}'.format(et_test)
    return fn_list
 
#==============================================================================
# merge files into cache files for each sample block   
#==============================================================================
def merge_3d_files(fn_list, save_path=None, verbose=False, 
                   calibration_fn=r"c:\MT\amtant.cal"):
    """
    merge .Z3D files into cache files.  Looks through the file list and 
    Combines files with the same start time and sampling rate into a 
    cache file.  The calibration file is copied to the merged path for 
    later use with mtft24.exe processing code.
    
    Arguments:
    ----------
        **fn_list** : list
                     list of files to be merged
                     
        **save_path** : directory to save cach files to
        
        **verbose** : [ True | False ]
                      * True --> prints out information about the merging
                      * False--> surpresses print statements
        
        **calibration_fn** : string
                             full path to calibration file for ANT's
                             
    Outputs:
    --------
        **merged_fn_list** : nested list of files that were merged together
        
        A log file is written to save_path\station_merged_log.log that contains
        information about the files that were merged together.
        
     :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> fn_list = zen.copy_from_sd('mt01', save_path=r"/home/mt/survey_1")
        >>> zen.merge_3d_files(fn_list, calibration_fn=r"/home/mt/amtant.cal")
    
    """
    
    start_time = time.ctime()
    merge_list = np.array([[fn]+\
                          os.path.basename(fn)[:-4].split('_')
                          for fn in fn_list if fn[-4:]=='.Z3D'])
                              
    merge_list = np.array([merge_list[:,0], 
                          merge_list[:,1],  
                          np.core.defchararray.add(merge_list[:,2],
                                                   merge_list[:,3]),
                          merge_list[:,4],
                          merge_list[:,5]])
    merge_list = merge_list.T
                              
    time_counts = Counter(merge_list[:,2])
    time_list = time_counts.keys()
    
    log_lines = []
  
    merged_fn_list = []
    for tt in time_list:
        log_lines.append('+'*72+'\n')
        log_lines.append('Files Being Merged: \n')
        cache_fn_list = merge_list[np.where(merge_list==tt)[0],0].tolist()
        
        for cfn in cache_fn_list:
            log_lines.append(' '*4+cfn+'\n')
        if save_path is None:
            save_path = os.path.dirname(cache_fn_list[0])
            station_name = merge_list[np.where(merge_list==tt)[0][0],1]
        else:
            save_path = save_path
            station_name = 'ZEN'
            
        zc = ZenCache()
        zc.verbose = verbose
        zc.write_cache_file(cache_fn_list, save_path, station=station_name)
            
        for zt in zc.zt_list:
            log_lines.append(zt.log_lines)
        merged_fn_list.append(zc.save_fn)
        log_lines.append('\n---> Merged Time Series Lengths and Start Time \n')
        log_lines.append(zc.log_lines)
        log_lines.append('\n')
    
    end_time = time.ctime()
    
    #copy the calibration file into the merged folder for mtft24
    try:
        copy_cal_fn = os.path.join(save_path, 'Merged',
                                 os.path.basename(calibration_fn))
    except:
        copy_cal_fn = os.path.join(save_path, os.path.basename(calibration_fn))
        
    shutil.copy(calibration_fn, copy_cal_fn)
    print 'copied {0} to {1}'.format(calibration_fn, copy_cal_fn)
    
    print 'Start time: {0}'.format(start_time)
    print 'End time:   {0}'.format(end_time)
    
    if os.path.basename(save_path) != 'Merged':
        log_fid = file(os.path.join(save_path, 'Merged', 
                                    station_name+'_Merged.log'), 'w')
    else:
        log_fid = file(os.path.join(save_path, station_name+'_Merged.log'),
                       'w')
    for line in log_lines:
        log_fid.writelines(line)
    log_fid.close()
        
    return merged_fn_list
    
#==============================================================================
# delete files from sd cards    
#==============================================================================
def delete_files_from_sd(delete_date=None, delete_type=None, 
                         delete_folder=r"d:\Peacock\MTData\Deleted",
                         verbose=True):
    """
    delete files from sd card, if delete_date is not None, anything on this 
    date and before will be deleted.  Deletes just .Z3D files, leaves 
    zenini.cfg
    
    Agruments:
    -----------
        **delete_date** : YYYY-MM-DD
                         date to delete files from 
                         
        **delete_type** : [ 'all' | 'before' | 'after' | 'on' ]
                          * 'all' --> delete all files on sd card
                          * 'before' --> delete files on and before delete_date
                          * 'after' --> delete files on and after delete_date
                          * 'on' --> delete files on delete_date
                          
        **delete_folder** : string
                            full path to a folder where files will be moved to
                            just in case.  If None, files will be deleted 
                            for ever.
                            
    Returns:
    ---------
        **delete_fn_list** : list
                            list of deleted files.
                            
     :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> # Delete all files before given date, forever. 
        >>> zen.delete_files_from_sd(delete_date='2004/04/20', 
                                     delete_type='before',
                                     delete_folder=None)
        >>> # Delete all files into a folder just in case 
        >>> zen.delete_files_from_sd(delete_type='all',
                                     delete_folder=r"/home/mt/deleted_files")
    
    """
    
    drive_names = get_drive_names()
    if drive_names is None:
        raise IOError('No drives to copy from.')

    log_lines = []
    if delete_folder is not None:
        if not os.path.exists(delete_folder):
            os.mkdir(delete_folder)
        log_fid = file(os.path.join(delete_folder,'Log_file.log'),'w')
    
    if delete_date is not None:
        delete_date = int(delete_date.replace('-',''))
    
    delete_fn_list = []
    for key in drive_names.keys():
        dr = r"{0}:\\".format(key)
        log_lines.append('='*25+drive_names[key]+'='*25+'\n')
        for fn in os.listdir(dr):
            if fn[-4:].lower() == '.Z3D'.lower():
                full_path_fn = os.path.normpath(os.path.join(dr, fn))
                zt = Zen3D(full_path_fn)
                #zt.get_info()
                if delete_type == 'all' or delete_date is None:
                    if delete_folder is None:
                        os.remove(full_path_fn)
                        delete_fn_list.append(full_path_fn)
                        log_lines.append('Deleted {0}'.format(full_path_fn))
                    else:
                        shutil.move(full_path_fn, 
                                    os.path.join(delete_folder,
                                    os.path.basename(full_path_fn)))
                        delete_fn_list.append(full_path_fn)
                        log_lines.append('Moved {0} '.format(full_path_fn)+
                                         'to {0}'.format(delete_folder))
                else:
                    zt_date = int(zt.schedule_date.replace('-',''))
                   
                    if delete_type == 'before':
                        if zt_date <= delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
                    elif delete_type == 'after':
                        if zt_date >= delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
                    elif delete_type == 'on':
                        if zt_date == delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_list.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
    if delete_folder is not None:
        log_fid = file(os.path.join(delete_folder, 'Delete_log.log'), 'w')
        log_fid.writelines(log_lines)
        log_fid.close()
    if verbose:
        for lline in log_lines:
            print lline
    
    return delete_fn_list
    
#==============================================================================
# copy and merge Z3D files from SD cards          
#==============================================================================
def copy_and_merge(station, z3d_save_path=None, merge_save_path=None, 
                   channel_dict={'1':'HX', '2':'HY', '3':'HZ','4':'EX', 
                                 '5':'EY', '6':'HZ'},
                   copy_date=None, copy_type='all'):
    """
    copy files from sd card then merge them together and run mtft24.exe
    
    Arguments:
    ----------
        **station** : string
                      full station name
                      
        **z3d_save_path** : string
                          full path to save .Z3D files
                          
        **merge_save_path** : string
                             full path to save merged cache files.  If None
                             saved to z3d_save_path\Merged
                             
        
        **channel_dict** : dictionary
                           keys are the channel numbers as strings and the
                           values are the component that corresponds to that 
                           channel, values are placed in upper case in the 
                           code
                           
        **copy_date** : YYYY-MM-DD
                        date to copy from depending on copy_type
                        
        **copy_type** : [ 'all' | 'before' | 'after' | 'on' ]
                        * 'all' --> copy all files on the SD card
                        * 'before' --> copy files before and on this date
                        * 'after' --> copy files on and after this date
                        * 'on' --> copy files on this date only
                        
    Returns:
    ------------
        **mfn_list** : list
                      list of merged file names
                      
    :Example: ::
    
        >>> import mpty.usgs.zen as zen
        >>> mfn_list = zen.copy_and_merge('mt01', z3d_save_path=r"/home/mt")
        >>> #copy only after a certain date
        >>> mfn_list = zen.copy_and_merge('mt01', z3d_save_path=r"/home/mt",\
                                          copy_date='2014/04/20', \
                                          copy_type='after')
    
    """
    
    #--> copy files from sd cards
    cpkwargs = {}
    cpkwargs['channel_dict'] = channel_dict
    cpkwargs['copy_date'] = copy_date
    cpkwargs['copy_type'] = copy_type
    if z3d_save_path != None:
        cpkwargs['save_path'] = z3d_save_path
    
    fn_list = copy_from_sd(station, **cpkwargs)
    
    #--> merge files into cache files
    mfn_list = merge_3d_files(fn_list, save_path=merge_save_path)
    
    return mfn_list
    
#==============================================================================
#   Make mtpy_mt files  
#==============================================================================
    
def make_mtpy_mt_files(fn_list, station_name='mb', fmt='%.8e', 
                       ex=1, ey=1, notch_dict=None, ey_skip=False):
    """
    makes mtpy_mt files from .Z3D files
    
    Arguments:
    -----------
        **dirpath** : full path to .Z3D files
        
        **station_name** : prefix for station names
        
        **fmt** : format of data numbers for mt_files
        
    Outputs:
    --------
        **fn_arr** : np.ndarray(file, length, df, start_dt)
        
    :Example: ::
    
        >>> import mtpy.usgs.zen as zen
        >>> fn_list = zen.copy_from_sd('mt01')
        >>> mtpy_fn = zen.make_mtpy_files(fn_list, station_name='mt')
    """
    
    fn_arr = np.zeros(len(fn_list), 
                      dtype=[('station','|S6'), ('len',np.int), ('df', np.int),
                             ('start_dt', '|S22'), ('comp','|S2'),
                             ('fn','|S100')])
    fn_lines = []
                             
    for ii, fn in enumerate(fn_list):
        zd = Zen3D(fn)
        
        #read in Z3D data
        try:
            zd.read_3d()
        except ZenGPSError:
            try:
                zd._seconds_diff = 59
                zd.read_3d()
            except ZenGPSError:
                pass
        if ey_skip and zd.ch_cmp == 'ey':
            pass
        else:
            #write mtpy mt file
            zd.write_ascii_mt_file(save_station=station_name, 
                                   fmt=fmt, 
                                   ex=ex, 
                                   ey=ey, 
                                   notch_dict=notch_dict)
            
            #create lines to write to a log file                       
            fn_arr[ii]['station'] = '{0}{1}'.format(station_name, zd.rx_stn)
            fn_arr[ii]['len'] = zd.time_series.shape[0]
            fn_arr[ii]['df'] = zd.df
            fn_arr[ii]['start_dt'] = zd.start_dt
            fn_arr[ii]['comp'] = zd.ch_cmp
            fn_arr[ii]['fn'] = zd.fn
            fn_lines.append(''.join(['--> station: {0}{1}\n'.format(station_name, 
                                                                     zd.rx_stn),
                                     '    ts_len = {0}\n'.format(zd.time_series.shape[0]),
                                     '    df = {0}\n'.format(zd.df),
                                     '    start_dt = {0}\n'.format(zd.start_dt),
                                     '    comp = {0}\n'.format(zd.ch_cmp),
                                     '    fn = {0}\n'.format(zd.fn)]))
        
    return fn_arr, fn_lines
    
#==============================================================================
# make time series loop
#==============================================================================
def make_mtpy_ts_loop(station_path, station_list, survey_file=None, 
                      station_name='mb', fmt='%.8e', notch_dict=None,
                      ey_skip=False):
    """
    loop over station folder to write mtpy time series
    
    Arguments:
    ----------
        **station_path** : directory of station folders
        
        **station_list** : list of stations to process
        
        **survey_file** : string
                          full path to survey_config file created by 
                          mtpy.utils.configfile
                          
        **station_name** : string
                           prefix to append to station name from Z3D files
                           
        **fmt** : string format of how the numbers are formated in new file
        
        **notch_dict** : dictionary
                         if the data has noise at single frequencies, such
                         as power line noise input a dictionary with keys:
                         
                        * df --> float sampling frequency in Hz
                 
                        * notches --> list of frequencies (Hz) to filter
                      
                        * notchradius --> float radius of the notch in 
                                          frequency domain (Hz)
        
                        * freqrad --> float radius to searching for peak about 
                                      notch from notches
                                  
                        * rp --> float ripple of Chebyshev type 1 filter, 
                                 lower numbers means less ripples
                             
                        * dbstop_limit --> float (in decibels) limits the 
                                           difference between the peak at the 
                                           notch and surrounding spectra. 
                                           Any difference above dbstop_limit
                                           will be filtered, anything
                                           less will not
                         
    """
    
    log_fid = file(os.path.join(station_path, 'TS_log.log'), 'a')
    
    if survey_file is not None:
        survey_dict = mtcf.read_survey_configfile(survey_file)
    
    for station in station_list:
        spath = os.path.join(station_path, station)
        if survey_file is not None:
            try:
                sdict = survey_dict[station.upper()]
                ex = float(sdict['e_xaxis_length'])
                ey = float(sdict['e_yaxis_length'])
            except KeyError:
                ex = 1.
                ey = 1.
        else:
            ex = 1.
            ey = 1.
        
        log_fid.write('-'*72+'\n')
        fn_list = [os.path.join(spath, fn) for fn in os.listdir(spath)
                  if fn[-3:]=='Z3D']
        sfn_arr, sfn_lines = make_mtpy_mt_files(fn_list, 
                                                station_name=station_name,
                                                fmt=fmt, 
                                                ex=ex, 
                                                ey=ey,
                                                notch_dict=notch_dict,
                                                ey_skip=ey_skip)
        log_fid.writelines(sfn_lines)
    log_fid.close()
    
#==============================================================================

# this should capture all the print statements from zen
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


#==============================================================================
def compute_mt_response(survey_dir, station='mt000', copy_date=None, 
                        birrp_exe=r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts.exe", 
                        ant_calibrations=r"c:\MT\Ant_calibrations",
                        process_df_list=[256],
                        num_comp=5):
    """
    This code will down load Z3D files from a Zen that is in SD Mode, 
    convert the Z3D files to ascii format, then process them for each
    sampling rate using Alan Chave's BIRRP code.  The outputs are then
    converted to .edi files and plotted.
    
    You need 2 things to run this code:
        * mtpy --> a Python package for MT and can be found at
    	           https://github.com/geophysics/mtpy
        * BIRRP executable --> you can get this from Alan Chave at WHOI 
                               if you are using it for non-commercial projects.
                               
    ..note:: This code is quite specific to my setup, so let me know what
             doesn't work for you so I can generalize it.
    
    Arguments
    ----------------
        **survey_dir** : string
                         full path to the directory where you are storing 
                         the station data.  ie. /home/MT/Survey_00
                         
        **station** : string
                      name of the station you are down loading.
                      *default* is 'mt000'
                      
        **copy_date** : string
                        copy all files on and after this date
                        format is YYYY-MM-DD
                        *default* is None, which copies all files on the SD
                        cards.
                                  
        **birrp_exe** : string
                        full path to the BIRRP executable on your machine
                        *default* is the location on my machine
        
        **ant_calibrations** : string
                               full path to a folder that contains the coil
                               calibration data.  These must be in seperate
                               .csv files for each coil named by corresponding
                               coil name. If you're coil is 2884, then you
                               need a calibration file named Ant2884_cal.csv
                               in which the data is freq,real,imaginary 
                               
        **process_df_list** : list
                              list of sampling rates to process
                              
                               
    Returns
    -----------------
        **rp_plot** : mtpy.imaging.plotnresponses object
                      ploting object of data, if you want to change how the
                      output is plot change the attributes of rp_plot
                      
    Outputs
    -----------------
        **copy_from_sd.log** : file
                               contains information on how files were copied
                               from the SD cards.
                               
        **processing.log** : file
                             a log file of how the program ran
        
        **survey_dir/station/TS** : directory
                                   contains the time series data in .ascii 
                                   format
                                   
        **survey_dir/station/TS/BF** : directory
                                       contains the processing results from
                                       BIRRP for each sampling rate in the
                                       data in subsequent directories
                             
        **survey_dir/station/TS/station.cfg** : file
                                                configuration file of the
                                                station parameters
                                                
        
                             
    Example
    ------------------------
        >>> import zen_processing_data as zpd
        >>> zpd.compute_mt_response(r"/home/mt/survey_00", 
                                    station='mt010',
                                    copy_date='2015-05-22',
                                    birrp_exe=r"/home/bin/birrp52.exe",
                                    ant_calibrations=r"/home/ant_calibrations",
                                    process_df_list=[1024, 256])
    """
                        
    station_dir = os.path.join(survey_dir, station)
               
    st = time.time()
    #--> Copy data from files
    try:
        if copy_date is None:
            copy_from_sd(station, save_path=survey_dir)
        else:
            copy_from_sd(station, save_path=survey_dir, 
                             copy_date=copy_date, copy_type='after')
    except IOError:
        print 'No files copied from SD cards'
        print 'Looking in  {0} for Z3D files'.format(station_dir)
    
    #--> process data
     
    with Capturing() as output:
        z2edi = Z3D_to_edi(station_dir)
        z2edi.birrp_exe = birrp_exe
        z2edi.coil_cal_path = ant_calibrations
        try:
            rp = z2edi.process_data(df_list=process_df_list, num_comp=num_comp)
        except mtex.MTpyError_inputarguments:
            print '==> Data not good!! Did not produce a proper .edi file' 
            et = time.time()
            print '--> took {0} seconds'.format(et-st)
            rp = None
    
    #--> write log file
    log_fid = open(os.path.join(station_dir, 'Processing.log'), 'w')
    log_fid.write('\n'.join(output))
    log_fid.close()
    
    return rp
    

#==============================================================================
def rename_cac_files(station_dir, station='mb'):
    """
    rename and move .cac files to something more useful
    """
    fn_list = [os.path.join(station_dir, fn) for fn in os.listdir(station_dir)
                if fn[-4:].lower() == '.cac']
                    
    if len(fn_list) == 0:
        raise IOError('Could not find any .cac files')
        
    save_path = os.path.join(station_dir, 'Merged')
    if not os.path.exists(save_path) :
        os.mkdir(save_path)
    
    for fn in fn_list:
        cac_obj = Cache(fn)
        cac_obj.read_cache_metadata()
        station_name = '{0}{1}'.format(station, 
                                       cac_obj.metadata.station_number)
        station_date = cac_obj.metadata.gdp_date.replace('-', '')
        station_time = cac_obj.metadata.gdp_time.replace(':', '')
        new_fn = '{0}_{1}_{2}_{3:.0f}.cac'.format(station_name,
                                                  station_date, 
                                                  station_time,
                                                  cac_obj.metadata.ts_adfreq)
        new_fn = os.path.join(save_path, new_fn)
        shutil.move(fn, new_fn)
        print 'moved {0} to {1}'.format(fn, new_fn)
