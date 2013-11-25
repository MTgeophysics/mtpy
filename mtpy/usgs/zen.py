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
import mtpy.utils.exceptions as mtex
import mtpy.utils.configfile as mtcf
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import mtpy.imaging.plotspectrogram as plotspectrogram

try:
    import mtpy.utils.mseed as mtmseed
except ImportError:
    pass

try:
    import mtpy.processing.filter as mtfilt
except ImportError:
    print 'no scipy.signal, its fucked thanks obspy'
    
#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'

class Zen3D(object):
    """
    Deal with the raw data output from the Zen box as Z3D files.
    
    Arguments:
    ----------
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
    gps_lst             list of gps stamps
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
    _stamp_lst          list of gps time stamp variables
    _week_len           length of a gps week in seconds
    =================== =======================================================
    """
    
    def __init__(self, fn=None, **kwargs):
        
        self.fn = fn
        self._header_len = kwargs.pop('header_len', 512)
        self._meta_len = kwargs.pop('meta_len', 512)
        self._stamp_len = kwargs.pop('stamp_len', 36)
        self._gps_stamp = kwargs.pop('gps_stamp', '\xff\xff\xff\xff')
        
        self._stamp_lst = ['gps', 'time', 'lat', 'lon', 'status', 
                           'gps_accuracy', 'temperature']
                           
        self._data_types = [np.int32, np.int32, np.float64, np.float64, 
                            np.uint32, np.int32, np.float32]
                            
        self._data_type = np.dtype([(st, dt) for st, dt in 
                                     zip(self._stamp_lst, self._data_types)])
                                     
        self._week_len = 604800
        self._gps_epoch = (1980, 1, 6, 0, 0, 0, -1, -1, 0)
        self._leap_seconds = 16
        
        #seconds different between scheduling time and actual collection time
        self._seconds_diff = 5 
        
        self.log_lines = []
        self.verbose = True
        self._skip_sample_tolerance = 5
        self.sample_diff_lst = []
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
        self.gps_lst = None
        self.temperature = None
        self.lat = None
        self.lon = None
        
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
        header_lst = header_string.replace('\n', ',').split(',')
        
        header_dict = {}
        for hh in header_lst:
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
                        hlst = ll.split('&')
                        for kk in hlst:
                            klst = kk.split(':')
                            header_dict[klst[0].strip().lower()] = klst[1].strip()
                    else:
                        hlst = ll.split(':')
                        try:
                            header_dict[hlst[0].strip().lower()] = hlst[1].strip()
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
        meta_lst = meta_data_string.replace('\n','|').split('|') 
        meta_dict = {}
        for mm in meta_lst:
            mlst = mm.split(',')
            if len(mlst) == 2:
                meta_dict[mlst[0].strip().lower()] = mlst[1].strip().lower()
            else:
                pass
        self.meta_dict = meta_dict  
        self.ch_number = meta_dict['ch.number']
        self.ch_cmp = meta_dict['ch.cmp'].replace('b','h')
        self.ch_length = meta_dict['ch.varasp']
        self.rx_stn = meta_dict['rx.stn']
        self.tx_id = meta_dict['tx.id']
        
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
        
        #---read in meta raw_data----------------------------------------------
        meta_string = raw_data[self._header_len-1:ds]
        self.read_metadata(meta_string)
        
        
        
    
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
        gps_lst = np.zeros(num_blocks, dtype=np.int)
        
        gps_dict = dict([(key, np.zeros(num_blocks, dtype=dtp)) 
                          for key, dtp in zip(self._stamp_lst, 
                                              self._data_types)])
        #make the time array floats instead of ints so can get the decimal 
        #place if it isn't 0.
        gps_dict['time'] = gps_dict['time'].astype(np.float32)
        
        #get gps information from the data
        #get first time stamp that matches the starting time
        s1 = 0
        gps_lst[0] = self.get_gps_stamp_location()
        gps_info = np.fromstring(raw_data[gps_lst[0]:gps_lst[0]+self._stamp_len], 
                                 dtype=self._data_type)
        gps_info['time'] = gps_info['time'].astype(np.float32)
        gps_info['time'] = self.get_gps_time(gps_info['time'])
        start_test = self.get_date_time(self.gps_week, gps_info['time'])
        
        #--> test to make sure the first time corresponds to the scheduled 
        #start time
        time_stop = 0
        while start_test != self.start_dt and s1 <= self._seconds_diff and \
                time_stop <= self._seconds_diff:
            s1 += 1
            gps_lst[0] = self.get_gps_stamp_location(gps_lst[0]+7)
            gps_info = np.fromstring(raw_data[gps_lst[0]:gps_lst[0]+\
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
                gps_lst[0] = self.get_gps_stamp_location()
                time_stop += 1  
       
       #----Raise an error if the first gps stamp is more than allowed time
        #    difference.
        if time_stop >= self._seconds_diff:
#            raise ZenGPSError('GPS start time is more than '+\
#                           '{0} '.format(self._seconds_diff)+\
#                           'seconds different than scheduled start time of '+\
#                           '{0}. \n '.format(self.start_dt)+\
#                           'Estimated start time is {0} +/- {1} sec'.format(
#                           start_test, self._seconds_diff))
            print ('GPS start time is more than '+\
                           '{0} '.format(self._seconds_diff)+\
                           'seconds different than scheduled start time of '+\
                           '{0}. \n '.format(self.start_dt)+\
                           'Estimated start time is {0} +/- {1} sec'.format(
                           start_test, self._seconds_diff))
                     
        #put the information into the correct arrays via dictionary                         
        for jj, key in enumerate(self._stamp_lst):
            gps_dict[key][0] = gps_info[0][jj]
  
        #find the next time stamp
        for ii in range(s1,num_blocks-1):
            sfind = self.get_gps_stamp_location(gps_lst[ii-1]+7)
            #make sure it isn't the same time stamp as before
            if sfind != gps_lst[ii-1] and sfind != -1:
                gps_info, gps_index, gps_week = self.get_gps_stamp(sfind)
                gps_lst[ii] = gps_index
                
                if gps_info is not None:
                    for jj, key in enumerate(self._stamp_lst):
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
                                      'point {0:^15},'.format(gps_lst[bb])+\
                                      'gps diff {0:^15}\n'.format(gps_diff[bb]))
            
            self.log_lines.append(' '*4+'*'*52+'\n')

        #need to be sure that the number of data points between time stamps is 
        #equal to the sampling rate, if it is not then remove that interval.  
        #Most likely it is at the beginning or end of time series.
        dsamples = np.array([(gps_lst[nn+1]-gps_lst[nn]-self._stamp_len-df*4)/4 
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
        
            gps_lst = gps_lst[bmin:bmax]
            
        num_samples = len(gps_lst)
        if self.verbose:
            print 'Found {0} gps time stamps, '.format(num_samples)+\
                  'with equal intervals of {0} samples'.format(int(self.df))
              
        self.log_lines.append(' '*4+\
                            'Found {0} gps time stamps, '.format(num_samples)+\
                  'with equal intervals of {0} samples\n'.format(int(self.df)))
        
        #read in data
        data_array = np.zeros((num_samples+1)*df, dtype=np.float32)
        for ll, kk in enumerate(gps_lst[0:-1]):
            pdiff = ((gps_lst[ll+1]-(kk+self._stamp_len))-(df*4))/4
            self.sample_diff_lst.append(pdiff)
            dblock = raw_data[kk+self._stamp_len:gps_lst[ll+1]] 
            try:
                data_array[ll*df:(ll+1)*df+pdiff] = np.fromstring(dblock, 
                                                                dtype=np.int32)
            except ValueError:
                print 'samples between time step {0} is off by {1} samples'.format(ll,
                                                                   abs(pdiff))
                          
        if sum(self.sample_diff_lst) != 0:
            if self.verbose:
                print 'time series is off by {0} seconds'.format(
                                           float(sum(self.sample_diff_lst))/df)
                self.log_lines.append('time series is off by {0} seconds'.format(
                                          float(sum(self.sample_diff_lst))/df))
                                           
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
        self.gps_lst = gps_lst
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
        
    def convert_counts(self):
        """
        convert the time series from counts to millivolts

        """
        
        return self.time_series*self.counts_to_mv_conversion
        
    def convert_mV(self):
        """
        convert millivolts to counts assuming no other scaling has been applied
        
        """
        
        return self.time_series/self.counts_to_mv_conversion
        
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
                print 'time', gps_index
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
            
            while gps_info['status'] < 0:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'status', gps_index                
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
            
            while abs(gps_info['temperature']) > 80:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'temperature', gps_index    
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
                
            while abs(gps_info['lat']) > np.pi:
                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
                print 'lat', gps_index
                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
                                         dtype=self._data_type)
                
#            while np.log10(abs(gps_info['lat'])) < -3:
#                gps_index = self.get_gps_stamp_location(start_index=gps_index+7)
#                print 'lat_log', gps_index    
#                gps_info = np.fromstring(self._raw_data[gps_index:gps_index+self._stamp_len], 
#                                         dtype=self._data_type)

            
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
#            if gps_info:
#                return -1, -1, -1
                
            return gps_info, gps_index, gps_week 
            
        except ValueError:
            print 'Ran into end of file, gps stamp not complete.'+\
                  ' Only {0} points.'.format(len(self._raw_data[gps_index:]))
            return None, gps_index, 0
            
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
        
    def get_degrees(self, radian_value):
        """
        convert lat or lon into decimal degrees
        
        """
        
        degrees = radian_value*180/np.pi
        
        return degrees
        
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
                  
        self.time_series, self.filt_lst = \
                    mtfilt.adaptive_notch_filter(self.time_series, **kwargs) 
        
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
            for nfilt in self.filt_lst:
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
        

class ZenCache(object):
    """
    deals with cache files or combined time series files.
    
    ================== ========================================================
     Attributes         Description
    ================== ======================================================== 
    cal_data            list of calibrations, as is from file
    fn_lst              list of filenames merged together
    log_lines           list of information to put into a log file
    meta_data           dictionary of meta data key words and values
    nav_data            list of navigation data, as is from file
    save_fn             file to save merged file to
    ts                  np.ndarray(len(ts), num_channels) of time series
    verbose             [ True | False ] True prints information to console
    zt_lst              list of class: Zen3D objects
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
        >>> fn_lst = [os.path.join(file_path, fn) 
        >>> ...       for fn in os.listdir(file_path)
        >>> ...       if fn.find('.Z3D')>0]
        >>> zc.write_cache_file(fn_lst, r"/home/MT/Station1", station='s1')
        >>> Saved File to: /home/MT/Station1/Merged/s1_20130601_190001_4096.cac
        
    """   
    
    def __init__(self):
        
        self.fn_lst = None
        self.zt_lst = None
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
    
    def check_sampling_rate(self, zt_lst):
        """
        check to make sure the sampling rate is the same for all channels
        
        Arguments:
        -----------
            **zt_lst** : list of Zen3D instances
        
        Outputs:
        --------
            **None** : raises an error if sampling rates are not all the same
            
        """
        
        nz = len(zt_lst)
        
        df_lst = np.zeros(nz)
        for ii, zt in enumerate(zt_lst):
            df_lst[ii] = zt.df
            
        tf_array = np.zeros((nz, nz))
        
        for jj in range(nz):
            tf_array[jj] = np.in1d(df_lst, [df_lst[jj]])
        
        false_test = np.where(tf_array==False)
        
        if len(false_test[0]) != 0:
            raise IOError('Sampling rates are not the same for all channels '+\
                          'Check file(s)'+zt_lst[false_test[0]])
        
    
    def check_time_series(self, zt_lst, decimate=1):
        """
        check to make sure timeseries line up with eachother.
        
        """
        
        n_fn = len(zt_lst)
        
        #test start time
        #st_lst = np.array([int(zt.date_time[0][-2:]) for zt in zt_lst])
        st_lst = np.array([int(zt.gps_time[0]) for zt in zt_lst])
        time_max = max(st_lst)
        #time_max = np.where(st_lst==st_lst.max())[0]
        
        #get the number of seconds each time series is off by
        skip_dict = {}
        for ii, zt in enumerate(list(zt_lst)):
            try:
                skip_dict[ii] = np.where(zt.gps_time==time_max)[0][0]
            except IndexError:
                zt_lst.remove(zt_lst[ii])
                print '***SKIPPING {0} '.format(zt.fn)
                print '   because it does not contain correct gps time'
                print '   {0} --> {1}'.format(time_max, 
                                             zt.get_date_time(zt.gps_week, 
                                                             time_max))
#        skip_dict = dict([(ii, np.where(zt.gps_time==time_max)[0][0]) 
#                          for ii, zt in enumerate(zt_lst)])
#        if len(time_max) != n_fn:
#            for ii in range(n_fn):
#                skip_dict[ii] = st_lst.max()-st_lst[ii]       
        
        #change data by amount needed        
        for ii, zt in zip(skip_dict.keys(), zt_lst):
            if skip_dict[ii] != 0:
                skip_points = skip_dict[ii]*zt.df
                print 'Skipping {0} points for {1}'.format(skip_points,
                                                            zt.ch_cmp)
                zt.time_series = zt.time_series[skip_points:]
                zt.gps_diff = zt.gps_diff[skip_dict[ii]:]
                zt.gps_lst = zt.gps_lst[skip_dict[ii]:]
                zt.date_time = zt.date_time[skip_dict[ii]:]
                zt.gps_time = zt.gps_time[skip_dict[ii]:]
            
        #test length of time series
        ts_len_lst = np.array([len(zt.time_series) for zt in zt_lst])
        
        #get the smallest number of points in the time series
        ts_min = ts_len_lst.min()
        
        #make a time series array for easy access
        ts_min /= decimate
            
        ts_array = np.zeros((ts_min, n_fn))
        
        #trim the time series if needed
        for ii, zt in enumerate(zt_lst):
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
        
    def write_cache_file(self, fn_lst, save_fn, station='ZEN', decimate=1):
        """
        write a cache file from given filenames
        
        """
        #sort the files so they are in order
        fn_sort_lst = []
        for cs in self.chn_order:
            for fn in fn_lst:
                if cs in fn.lower():
                    fn_sort_lst.append(fn)

        fn_lst = fn_sort_lst
        print fn_lst
            
        n_fn = len(fn_lst)
        self.zt_lst = []
        for fn in fn_lst:
            zt1 = Zen3D(fn=fn)
            zt1.verbose = self.verbose
            try:
                zt1.read_3d()
            except ZenGPSError:
                zt1._seconds_diff = 59
                zt1.read_3d()
            self.zt_lst.append(zt1)
        
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
        self.check_sampling_rate(self.zt_lst)
        
        #make sure the length of time series is the same for all channels
        self.ts, ts_len = self.check_time_series(self.zt_lst,
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
        meta_lst = cdata[ii:jj].split('\n')
        for mm in meta_lst:
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
        meta_lst = cdata[ii:jj].split('\n')
        for mm in meta_lst:
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
    
    :Example: ::
    
        >>> import myp.usgs.zen as zen
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
    df_lst                 sequential list of sampling rates to repeat in 
                           schedule
    df_time_lst            sequential list of time intervals to measure for
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
    sa_lst                 list of schedule actions including time and df
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
        self.sa_lst = []
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
        self.df_lst = (4096, 1024, 256)
        self.df_time_lst = ('00:05:00','00:15:00','05:40:00')
        self.master_schedule = self.make_schedule(self.df_lst, 
                                                  self.df_time_lst,
                                                  repeat=21)
                                                  

    def read_schedule(self, fn):
        """
        read zen schedule file
        
        """
        
        sfid = file(fn, 'r')
        lines = sfid.readlines()
        
        for line in lines:
            if line.find('scheduleaction') == 0:
                line_lst = line.strip().split(' ')[1].split(',')
                sa_dict = {}
                for ii, key in enumerate(self.sa_keys):
                    sa_dict[key] = line_lst[ii]
                self.sa_lst.append(sa_dict)
                
            elif line.find('metadata'.upper()) == 0:
                line_lst = line.strip().split(' ')[1].split('|')
                for md in line_lst[:-1]:
                    md_lst = md.strip().split(',')
                    self.meta_dict[md_lst[0]] = md_lst[1]
                    
            elif line.find('offset') == 0:
                line_str = line.strip().split(' ')
                self.offset = line_str[1]
                
            elif line.find('Light') > 0:
                line_lst = line.strip().split(' ')
                try:
                    self.light_dict[line_lst[0]]
                    self.light_dict[line_lst[0]] = line_lst[1]
                except KeyError:
                    pass
                
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
        
    def make_schedule(self, df_lst, df_length_lst, repeat=5, t1_dict=None):
        """
        make a repeated schedule given list of sampling frequencies and
        duration for each.
        
        Arguments:
        -----------
            **df_lst** : list   
                         list of sampling frequencies in Hz, note needs to be
                         powers of 2 starting at 256
            **df_length_lst** : list
                                list of durations in hh:mm:ss format
            **repeat** : int
                         number of times to repeat the sequence
            
            **t1_dict** : dictionary
                          dictionary returned from get_schedule_offset
                          
        Returns:
        --------
            **time_lst**: list of dictionaries with keys:
                            * 'dt' --> date and time of schedule event
                            * 'df' --> sampling rate for that event
        """
    
        df_lst = np.array(df_lst)
        df_length_lst = np.array(df_length_lst)
        ndf = len(df_lst)
        
    
        if t1_dict is not None:
            time_lst = [{'dt':self.initial_dt,'df':t1_dict['df']}]
    
            kk = np.where(np.array(df_lst)==t1_dict['df'])[0]-ndf+1
            df_lst = np.append(df_lst[kk:], df_lst[:kk])
            df_length_lst = np.append(df_length_lst[kk:], df_length_lst[:kk])
            print df_lst, df_length_lst
            time_lst.append(dict([('dt',t1_dict['dt']), ('df',df_lst[0])]))
            ii = 1
        else:
            time_lst = [{'dt':self.initial_dt,'df':df_lst[0]}]
            ii = 0
            
        for rr in range(1,repeat+1):
            for df, df_length, jj in zip(df_lst, df_length_lst, range(ndf)):
                dtime = time.strptime(df_length, '%H:%M:%S')
                ndt = self.add_time(time_lst[ii]['dt'], 
                                    add_hours=dtime.tm_hour,
                                    add_minutes=dtime.tm_min,
                                    add_seconds=dtime.tm_sec)
                time_lst.append({'dt':ndt.strftime(self.dt_format),
                                 'df':df_lst[jj-ndf+1]})
                ii += 1
                
        for nn, ns in enumerate(time_lst):
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
            
        return time_lst
        
    def get_schedule_offset(self, time_offset, schedule_time_lst):
        """
        gets the offset in time from master schedule list and time_offset so 
        that all schedules will record at the same time according to master
        schedule list schedule_time_lst
        
        Attributes:
        -----------
            **time_offset** : hh:mm:ss
                              the time offset given to the zen reciever
                              
            **schedule_time_lst** : list
                                    list of actual schedule times returned 
                                    from make_schedule
                                    
        Returns:
        --------
            **s1** : dictionary
                     dictionary with keys:
                         * 'dt' --> date and time of offset from next schedule
                                    event from schedule_time_lst
                         * 'df' --> sampling rate of that event
        """
        
        dt_offset = '{0},{1}'.format('2000-01-01', time_offset)
        t0 = time.mktime(time.strptime('2000-01-01,00:00:00', self.dt_format))
        
        for ii, tt in enumerate(schedule_time_lst):
            ssec = time.mktime(time.strptime(tt['dt'], self.dt_format))
            osec = time.mktime(time.strptime(dt_offset, self.dt_format))
            
            if ssec > osec:
                sdiff = time.localtime(t0+(ssec-osec))
                t1 = self.add_time('2000-01-01,00:00:00', 
                                   add_hours=sdiff.tm_hour,
                                   add_minutes=sdiff.tm_min,
                                   add_seconds=sdiff.tm_sec)
                s1 = {'dt':t1.strftime(self.dt_format), 
                      'df':schedule_time_lst[ii-1]['df']}
                return s1
                
    def write_schedule(self, station, clear_schedule=True, 
                       clear_metadata=True, varaspace=100, 
                       savename=0, dt_offset=None, 
                       df_lst=None, 
                       df_time_lst=None, 
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
                                      
            **df_lst** : list
                         list of sampling rates in Hz
            
            **df_time_lst** : list
                              list of time intervals corresponding to df_lst
                              in hh:mm:ss format
            **repeat** : int
                         number of time to repeat the cycle of df_lst
            
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

        if df_lst is not None:
            self.df_lst = df_lst
        if df_time_lst is not None:
            self.df_time_lst = df_time_lst
        
        self.master_schedule =  self.make_schedule(self.df_lst, 
                                                  self.df_time_lst,
                                                  repeat=repeat*3)
                                                  
        self.sa_lst = self.make_schedule(self.df_lst,
                                          self.df_time_lst,
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
            for sa_dict in self.sa_lst:
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
                for sa_dict in self.sa_lst:
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
            for sa_dict in self.sa_lst:
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
                                                   
                                                  
#==============================================================================
# interface with birrp                        
#==============================================================================
class ZenBIRRP():
    """
    class to deal with Birrp from Zen outputs
    
    survey file is .cfg file
    
    """              

    def __init__(self, station_path, **kwargs):
        
        self.station_path = station_path
        self.rr_path = kwargs.pop('rr_path', self.station_path)
        self.survey_config_fn = kwargs.pop('survey_config_fn', None)
        self.processing_file = kwargs.pop('processing_file', None)
        self.calibration_path = kwargs.pop('calibration_path', 
                                         r"d:\Peacock\MTData\Ant_calibrations")
        self.calibration_lst = ['2254', '2264', '2274', '2284', '2294',
                                '2304', '2314', '2324', '2334', '2344']
        self.birrp_dict = kwargs.pop('birrp_dict', None)
        self.station = kwargs.pop('station', 
                        os.path.basename(os.path.dirname(self.station_path)))
        self.rrstation = None
        self.rrsurvey_dict = None
        self.df = kwargs.pop('df', 256)
        self.processing_dict = kwargs.pop('processing_dict', None)
        self.survey_dict = kwargs.pop('survey_dict', None)
        self.birrp_exe = r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts_big.exe"
        self.script_file = None
        self.output_path = None
        self.birrp_config_fn = None
        self.calibration_dict = {}
        for cal_fn in os.listdir(self.calibration_path):
            for cal_num in self.calibration_lst:
                if cal_num in cal_fn:
                    self.calibration_dict[cal_num] = \
                                    os.path.join(self.calibration_path, cal_fn)
                    break
            
        
    def get_birrp_parameters(self, processing_file=None):
        """
        get parameters to put into birrp from file
        
        """
        if processing_file is not None:
            self.processing_file = processing_file
        if self.processing_file is None:
            raise IOError('Need to input a processing file')
            
        processing_lst = read_processing_file(self.processing_file)
        for pdict in processing_lst:
            if pdict['station'] == self.station and \
                                        float(pdict['df']) == self.df:
                return pdict

    
    def get_survey_parameters(self, survey_config_fn=None, rrstation=None):
        """
        get survey parameters from file
        
        """
        if survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
        if self.survey_config_fn is None:
            raise IOError('Need to input a survey config file')
            
        survey_dict_lst = mtcf.read_survey_configfile(self.survey_config_fn)
        
        try:
            self.survey_dict = survey_dict_lst[self.station.upper()]
        except KeyError:
            print 'Did not find station information in {0}'.format(
                                                        self.survey_config_fn)
        
        if self.rrstation is not None:                                                
            try:
                self.rrsurvey_dict = survey_dict_lst[self.rrstation.upper()]
            except KeyError:
                print 'Did not find remote station information in {0}'.format(
                                                      self.survey_config_fn)
                    
    def set_remote_reference(self, rrstation, rrpath=None):
        """
        set remote reference station and find survey information and filenames
        """ 
        self.rrstation = rrstation
        if rrpath is not None:
            self.rr_path = rrpath
        elif not os.path.exists(self.rr_path):
            rr_path = self.station_path
            kk = 0
            while os.path.basename(rr_path) != self.station and kk < 5:
                rr_path = os.path.dirname(rr_path)
                kk += 1
            self.rr_path = os.path.join(os.path.dirname(rr_path), 
                                        self.rrstation, 'TS')
            if not os.path.exists(self.rr_path):
                raise IOError('Need to input rrpath, could not find it')
            
        self.get_survey_parameters()
                    
                
    def get_fn_lst(self, df=None, start_dt=None, end_dt=None, ncomps=5):
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
            
        fn_lst = []
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
                    fn_lst.append(tarr)
                    ii = 0
            except mtex.MTpyError_ts_data:
                pass
            except mtex.MTpyError_inputarguments:
                pass
        
        #get remote reference time series
        rrfn_lst = []
        ii = 0
        for fn in os.listdir(self.rr_path):
            try:
                if np.remainder(ii, 2) == 0:
                    tarr = np.zeros(2, dtype=[('fn','|S100'),
                                               ('npts',np.int),
                                               ('start_dt','|S19'),
                                               ('end_dt','|S19')])
                header_dict = \
                        mtfh.read_ts_header(os.path.join(self.station_path,fn))
                if header_dict['t_min'] >= start_seconds and \
                   header_dict['t_min'] <= end_seconds and \
                   header_dict['samplingrate'] == float(self.df):
                    
                    try:
                        kk = rrcomp_dict[header_dict['channel'].lower()]
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
                    except KeyError:
                        pass
                if ii == 2:
                    rrfn_lst.append(tarr)
                    ii = 0
            except mtex.MTpyError_ts_data:
                pass
            except mtex.MTpyError_inputarguments:
                pass
        
        return fn_lst, rrfn_lst
                                    
    def write_script_file(self, df=None, processing_file=None,
                          processing_dict=None, start_dt=None, end_dt=None,
                          ncomps=5, jmode=0, survey_config_fn=None):
        """
        write a script file to guide birrp
        
        """       
        
        if df is not None:
            self.df = df
        
        #--> get survey parameters
        if survey_config_fn is not None:
            self.survey_config_fn = survey_config_fn
            self.get_survey_parameters()
        elif self.survey_dict is None:
            self.get_survey_parameters()
            
        #--> get processing dictionary 
        if processing_file is not None:
            self.processing_file = processing_file
        
        self.processing_dict = self.get_birrp_parameters()
        
        if processing_dict is not None:
            self.processing_dict = processing_dict
            
        if self.processing_file is None and self.processing_dict is None:
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
            self.set_remote_reference(self.processing_dict['rrstation'])
        except KeyError:
            self.set_remote_reference(self.station)
        
        #get list of files to process from the station folder
        fn_lst, rrfn_lst = self.get_fn_lst(self.df, 
                                           start_dt=start_dt, 
                                           end_dt=end_dt, 
                                           ncomps=ncomps)
        
        
        self.processing_dict['fn_lst'] = [fnlst['fn'] for fnlst in fn_lst]
        self.processing_dict['rrfn_lst'] = [rrfnlst['fn'] 
                             for rrfnlst in rrfn_lst]
        
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
            self.processing_dict['nread'] = [fnlst['npts'].min() 
                                  for fnlst in fn_lst]
        
        #if jmode == 1 for entering start and end times
        elif self.processing_dict['jmode'] == 1:
            self.processing_dict['dstim'] = [fnlst['start_dt'] 
                                             for fnlst in fn_lst]
            self.processing_dict['wstim'] = [fnlst['start_dt'] 
                                             for fnlst in fn_lst]
            self.processing_dict['wetim'] = [fnlst['end_dt'] 
                                             for fnlst in fn_lst]
                                                 
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
            
        if self.rrsurvey_dict is not None:
            try:
                self.processing_dict['rrhx_cal'] = \
                                self.calibration_dict[self.rrsurvey_dict['hx']]
            except KeyError:
                print 'Did not find RRHX calibration in {0}'.format(
                                                    self.survey_config_fn)
                self.processing_dict['rrhx_cal'] = \
                                            self.calibration_dict['2284']
                print 'Setting calibration coil number to 2284 as default.' 
                
            try:
                self.processing_dict['rrhy_cal'] = \
                                self.calibration_dict[self.rrsurvey_dict['hy']]
            except KeyError:
                print 'Did not find RRHY calibration in {0}'.format(
                                                        self.survey_config_fn)
                self.processing_dict['rrhy_cal'] = \
                                            self.calibration_dict['2284']
                print 'Setting calibration coil number to 2284 as default.' 
        
        #set the save path to include the sampling rate
        self.output_path = os.path.join(os.path.dirname(
                                        self.processing_dict['fn_lst'][0][0]), 
                                        'BF_{0}'.format(self.df))
        #write script file using mtpy.processing.birrp    
        script_file, birrp_dict = birrp.write_script_file(dict(self.processing_dict),
                                                    save_path=self.output_path)
        
        cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(
                                                             script_file[:-7]))
        print cfg_fn        
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
# read processing file
#==============================================================================
def read_processing_file(processing_file, delimiter='\t'):
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
        **processing_file** : string (full path to file)
                              tab delimited text file with appropriate 
                              information.
                              
        **station_path** : directory path to where station folders are
        
    Outputs:
    --------
        **slst** : list of dictionaries with key words related to headers of 
                   txt file
    
    """
    
    pfid = open(processing_file, 'r')
    plines = pfid.readlines()
    pkeys = plines[0].rstrip()
    pkeys = pkeys.split('\t')
    plst=[]
    for pline in plines[1:]:
        pstr = pline.rstrip()
        pstr = pstr.split(delimiter)
        if len(pstr)>1:
            pdict={}
            for kk, pkey in enumerate(pkeys):
                pstr[kk] = pstr[kk].replace('"','')
                pdict[pkey.lower()] = pstr[kk]
        plst.append(pdict)
        
    pfid.close()
    
    return plst
#==============================================================================
# get the external drives for SD cards
#==============================================================================
def get_drives():
    """
    get a list of logical drives detected on the machine
    
    Note this only works for windows.
    
    Outputs:
    --------
        **drives** : list of drives as letters
    
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
    --------
        **drive_dict** : dictionary
                         keys are the drive letters and values are the 
                         drive names
    
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
def copy_from_sd(station, savepath=r"d:\Peacock\MTData", 
                 channel_dict={'1':'HX', '2':'HY', '3':'HZ',
                               '4':'EX', '5':'EY', '6':'HZ'},
                 copy_date=None, copy_type='all'):
    """
    copy files from sd cards into a common folder
    
    do not put an underscore in station, causes problems at the moment
    
    Arguments:
    -----------
        **station** : string
                      full name of station from which data is being saved
        
        **savepath** : string
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
    --------
        **fn_lst** : list
                     list of filenames copied to savepath
    
    """
    
    drive_names = get_drive_names()
    if drive_names is None:
        raise IOError('No drives to copy from.')
    save_path = os.path.join(savepath,station)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    log_fid = file(os.path.join(save_path,'Log_file.log'),'w')
    
    
    st_test = time.ctime()
    fn_lst = []
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
                    zt.get_info()
                    
                    if zt.schedule_date is not None:
                        fn_find = True
                        if copy_date is not None:
                            cp_date = int(''.join(copy_date.split('-')))
                            
                            fn_find = False
                            
                            zt_date = int(''.join(zt.schedule_date.split('-')))
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
                            channel = channel_dict[drive_names[key][-1]]
                            st = zt.schedule_time.replace(':','')
                            sd = zt.schedule_date.replace('-','')
                            sv_fn = '{0}_{1}_{2}_{3}_{4}.Z3D'.format(station, 
                                                                     sd, 
                                                                     st,
                                                                     int(zt.df),
                                                                     channel)
                                                                 
                            full_path_sv = os.path.join(save_path, sv_fn)
                            fn_lst.append(full_path_sv)
                            
                            shutil.copy(full_path_fn, full_path_sv)
                            print 'copied {0} to {1}\n'.format(full_path_fn, 
                                                             full_path_sv)
                                                             
                            log_fid.writelines(zt.log_lines)
                                                             
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
    return fn_lst
 
#==============================================================================
# merge files into cache files for each sample block   
#==============================================================================
def merge_3d_files(fn_lst, savepath=None, verbose=False, 
                   calibration_fn=r"c:\MT\amtant.cal"):
    """
    merge .Z3D files into cache files.  Looks through the file list and 
    Combines files with the same start time and sampling rate into a 
    cache file.  The calibration file is copied to the merged path for 
    later use with mtft24.exe processing code.
    
    Arguments:
    ----------
        **fn_lst** : list
                     list of files to be merged
                     
        **savepath** : directory to save cach files to
        
        **verbose** : [ True | False ]
                      * True --> prints out information about the merging
                      * False--> surpresses print statements
        
        **calibration_fn** : string
                             full path to calibration file for ANT's
                             
    Outputs:
    --------
        **merged_fn_lst** : nested list of files that were merged together
        
        A log file is written to savepath\station_merged_log.log that contains
        information about the files that were merged together.
    
    """
    
    start_time = time.ctime()
    merge_lst = np.array([[fn]+\
                          os.path.basename(fn)[:-4].split('_')
                          for fn in fn_lst if fn[-4:]=='.Z3D'])
                              
    merge_lst = np.array([merge_lst[:,0], 
                          merge_lst[:,1],  
                          np.core.defchararray.add(merge_lst[:,2],
                                                   merge_lst[:,3]),
                          merge_lst[:,4],
                          merge_lst[:,5]])
    merge_lst = merge_lst.T
                              
    time_counts = Counter(merge_lst[:,2])
    time_lst = time_counts.keys()
    
    log_lines = []
  
    merged_fn_lst = []
    for tt in time_lst:
        log_lines.append('+'*72+'\n')
        log_lines.append('Files Being Merged: \n')
        cache_fn_lst = merge_lst[np.where(merge_lst==tt)[0],0].tolist()
        
        for cfn in cache_fn_lst:
            log_lines.append(' '*4+cfn+'\n')
        if savepath is None:
            save_path = os.path.dirname(cache_fn_lst[0])
            station_name = merge_lst[np.where(merge_lst==tt)[0][0],1]
        else:
            save_path = savepath
            station_name = 'ZEN'
            
        zc = ZenCache()
        zc.verbose = verbose
        zc.write_cache_file(cache_fn_lst, save_path, station=station_name)
            
        for zt in zc.zt_lst:
            log_lines.append(zt.log_lines)
        merged_fn_lst.append(zc.save_fn)
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
        
    return merged_fn_lst
    
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
        **delete_fn_lst** : list
                            list of deleted files.
    
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
    
    delete_fn_lst = []
    for key in drive_names.keys():
        dr = r"{0}:\\".format(key)
        log_lines.append('='*25+drive_names[key]+'='*25+'\n')
        for fn in os.listdir(dr):
            if fn[-4:].lower() == '.Z3D'.lower():
                full_path_fn = os.path.normpath(os.path.join(dr, fn))
                zt = Zen3D(full_path_fn)
                zt.get_info()
                if delete_type == 'all' or delete_date is None:
                    if delete_folder is None:
                        os.remove(full_path_fn)
                        delete_fn_lst.append(full_path_fn)
                        log_lines.append('Deleted {0}'.format(full_path_fn))
                    else:
                        shutil.move(full_path_fn, 
                                    os.path.join(delete_folder,
                                    os.path.basename(full_path_fn)))
                        delete_fn_lst.append(full_path_fn)
                        log_lines.append('Moved {0} '.format(full_path_fn)+
                                         'to {0}'.format(delete_folder))
                else:
                    zt_date = int(zt.schedule_date.replace('-',''))
                   
                    if delete_type == 'before':
                        if zt_date <= delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
                    elif delete_type == 'after':
                        if zt_date >= delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
                    elif delete_type == 'on':
                        if zt_date == delete_date:
                            if delete_folder is None:
                                os.remove(full_path_fn)
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Deleted {0}\n'.format(full_path_fn))
                            else:
                                shutil.move(full_path_fn, 
                                            os.path.join(delete_folder,
                                            os.path.basename(full_path_fn)))
                                delete_fn_lst.append(full_path_fn)
                                log_lines.append('Moved {0} '.format(full_path_fn)+
                                                 'to {0}\n'.format(delete_folder))
    if delete_folder is not None:
        log_fid = file(os.path.join(delete_folder, 'Delete_log.log'), 'w')
        log_fid.writelines(log_lines)
        log_fid.close()
    if verbose:
        for lline in log_lines:
            print lline
    
    return delete_fn_lst
    
#==============================================================================
# copy and merge Z3D files from SD cards          
#==============================================================================
def copy_and_merge(station, z3d_savepath=None, merge_savepath=None, 
                   channel_dict={'1':'HX', '2':'HY', '3':'HZ','4':'EX', 
                                 '5':'EY', '6':'HZ'},
                   copy_date=None, copy_type='all'):
    """
    copy files from sd card then merge them together and run mtft24.exe
    
    Arguments:
    ----------
        **station** : string
                      full station name
                      
        **z3d_savepath** : string
                          full path to save .Z3D files
                          
        **merge_savepath** : string
                             full path to save merged cache files.  If None
                             saved to z3d_savepath\Merged
                             
        
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
        
    
    """
    
    #--> copy files from sd cards
    cpkwargs = {}
    cpkwargs['channel_dict'] = channel_dict
    cpkwargs['copy_date'] = copy_date
    cpkwargs['copy_type'] = copy_type
    if z3d_savepath != None:
        cpkwargs['savepath'] = z3d_savepath
    
    fn_lst = copy_from_sd(station, **cpkwargs)
    
    #--> merge files into cache files
    mfn_lst = merge_3d_files(fn_lst, savepath=merge_savepath)
    
    return mfn_lst
    
#==============================================================================
#   Make mtpy_mt files  
#==============================================================================
    
def make_mtpy_mt_files(fn_lst, station_name='mb', fmt='%.8e', 
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
    """
    
    fn_arr = np.zeros(len(fn_lst), 
                      dtype=[('station','|S6'), ('len',np.int), ('df', np.int),
                             ('start_dt', '|S22'), ('comp','|S2'),
                             ('fn','|S100')])
    fn_lines = []
                             
    for ii, fn in enumerate(fn_lst):
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
def make_mtpy_ts_loop(station_path, station_lst, survey_file=None, 
                      station_name='mb', fmt='%.8e', notch_dict=None,
                      ey_skip=False):
    """
    loop over station folder to write mtpy time series
    
    Arguments:
    ----------
        **station_path** : directory of station folders
        
        **station_lst** : list of stations to process
        
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
    
    for station in station_lst:
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
        fn_lst = [os.path.join(spath, fn) for fn in os.listdir(spath)
                  if fn[-3:]=='Z3D']
        sfn_arr, sfn_lines = make_mtpy_mt_files(fn_lst, 
                                                station_name=station_name,
                                                fmt=fmt, 
                                                ex=ex, 
                                                ey=ey,
                                                notch_dict=notch_dict,
                                                ey_skip=ey_skip)
        log_fid.writelines(sfn_lines)
    log_fid.close()
