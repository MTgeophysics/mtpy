# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:37:13 2017

@author: jpeacock
"""

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

# ==============================================================================
import time
import datetime
import os
import struct
import string
import shutil
from cStringIO import StringIO
import sys
import numpy as np
import scipy.signal as sps
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
import mtpy.processing.filter as mtfilt
import mtpy.core.edi as mtedi
import mtpy.core.ts as mtts

try:
    import win32api
except ImportError:
    print "Cannot find win32api, will not be able to detect drive names"

try:
    import mtpy.utils.mseed as mtmseed
except ImportError:
    print (
        "Can not convert data to mini seed format need to install Obspy, "
        "good luck! You can find information on Obspy at "
        "https://github.com/obspy/obspy/wiki"
    )

# ==============================================================================
datetime_fmt = "%Y-%m-%d,%H:%M:%S"
datetime_sec = "%Y-%m-%d %H:%M:%S"
# ==============================================================================
#
# ==============================================================================
# ==============================================================================
#  for older Z3d files
# ==============================================================================
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
        self._header_len = kwargs.pop("header_len", 512)
        self._meta_len = kwargs.pop("meta_len", 512)
        self._stamp_len = kwargs.pop("stamp_len", 36)
        self._gps_stamp = kwargs.pop("gps_stamp", "\xff\xff\xff\xff")

        self._stamp_list = [
            "gps",
            "time",
            "lat",
            "lon",
            "status",
            "gps_accuracy",
            "temperature",
        ]

        self._data_types = [
            np.int32,
            np.int32,
            np.float64,
            np.float64,
            np.uint32,
            np.int32,
            np.float32,
        ]

        self._data_type = np.dtype(
            [(st, dt) for st, dt in zip(self._stamp_list, self._data_types)]
        )

        self._week_len = 604800
        self._gps_epoch = (1980, 1, 6, 0, 0, 0, -1, -1, 0)
        self._leap_seconds = 16

        # seconds different between scheduling time and actual collection time
        self._seconds_diff = 5

        self.log_lines = []
        self.verbose = True
        self._skip_sample_tolerance = 5
        self.sample_diff_list = []
        self.counts_to_mv_conversion = 9.5367431640625e-10
        self.units = "counts"
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

    # ==================================================
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

        # ----read in header information----------------------------------------
        header_list = header_string.replace("\n", ",").split(",")

        header_dict = {}
        for hh in header_list:
            if hh != "" and hh.find("builddate") == -1:
                hkv = hh.split(":")
                if len(hkv) == 2:
                    if hkv[0].lower() == "period" or hkv[0].lower() == "duty":
                        try:
                            header_dict[hkv[0].strip().lower()] += hkv[1].strip()
                        except KeyError:
                            header_dict[hkv[0].strip().lower()] = hkv[1].strip()
                    else:
                        header_dict[hkv[0].strip().lower()] = hkv[1].strip()
                elif len(hkv) == 3:
                    header_dict["start_time"] = hh.strip()
                else:
                    pass
            elif hh == "":
                pass
            else:
                hline = hh.split(";")
                for ll in hline:
                    if ll.find("builddate") > 0:
                        hlist = ll.split("&")
                        for kk in hlist:
                            klist = kk.split(":")
                            header_dict[klist[0].strip().lower()] = klist[1].strip()
                    else:
                        hlist = ll.split(":")
                        try:
                            header_dict[hlist[0].strip().lower()] = hlist[1].strip()
                        except IndexError:
                            pass
        # make attributes that will be useful latter
        self.header_dict = header_dict
        self.df = float(header_dict["a/d rate"])
        self.gain = float(header_dict["a/d gain"])
        self.gps_week = int(header_dict["gpsweek"])
        try:
            self.schedule_date = header_dict["schedule for this file"]
        except KeyError:
            self.schedule_date = header_dict["schedule"]
        self.schedule_time = header_dict["start_time"]

        # get the start date/time in UTC time
        self.start_dt = self.compute_schedule_start(
            self.schedule_date, self.schedule_time
        )
        self.start_time = self.schedule_time
        self.start_date = self.schedule_date

        # --> get serial number of a/d board
        try:
            self.ch_adcard_sn = header_dict["serial"]
        except KeyError:
            self.ch_adcard_sn = header_dict["brd339 serial"]

    # ==================================================
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
        meta_list = meta_data_string.replace("\n", "|").split("|")
        meta_dict = {}
        for mm in meta_list:
            mlist = mm.split(",")
            if len(mlist) == 2:
                meta_dict[mlist[0].strip().lower()] = mlist[1].strip().lower()
            else:
                pass
        self.meta_dict = meta_dict
        self.ch_number = meta_dict["ch.number"]
        self.ch_cmp = meta_dict["ch.cmp"].replace("b", "h")
        self.ch_length = meta_dict["ch.varasp"]
        self.rx_stn = meta_dict["rx.stn"]
        self.tx_id = meta_dict["tx.id"]

    # ==================================================
    def get_info(self):
        """
        read header and meta data
        """

        # beginning index of data blocks
        ds = self._header_len + self._meta_len

        # read in as a binary file.
        rfid = open(self.fn, "rb")
        raw_data = rfid.read(ds + 4)
        self._raw_data = raw_data
        rfid.close()

        if len(raw_data) < ds:
            print "Data file is not complete cannot read header information"
            return

        try:
            self.log_lines[0] != "-" * 72 + "\n"
        except IndexError:
            self.log_lines.append("-" * 72 + "\n")
            self.log_lines.append("--> Reading File: {0}\n".format(self.fn))

        # ----read in header information----------------------------------------
        header_string = raw_data[0 : self._header_len]
        self.read_header(header_string)

        print ("-" * 40)
        print ("   ad card sn     =  {0}".format(self.ch_adcard_sn))
        print ("   sampling rate  =  {0:.0f}".format(self.df))
        print ("   gain           =  {0:.1f}".format(self.gain))
        print ("   gps_week       =  {0:.0f}".format(self.gps_week))
        print ("   schedule date  =  {0}".format(self.schedule_date))
        print ("   schedule time  =  {0}".format(self.schedule_time))

        # ---read in meta raw_data----------------------------------------------
        meta_string = raw_data[self._header_len - 1 : ds]
        self.read_metadata(meta_string)

        print ("   channel no     =  {0}".format(self.ch_number))
        print ("   channel comp   =  {0}".format(self.ch_cmp))
        print ("   channel len    =  {0}".format(self.ch_length))
        print ("   rx station     =  {0}".format(self.rx_stn))
        print ("   tx id          =  {0}".format(self.tx_id))
        print ("-" * 40)

    # ==================================================
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
        # read in as a binary file.
        raw_data = open(self.fn, "rb").read()
        self._raw_data = raw_data

        try:
            self.log_lines[0] != "-" * 72 + "\n"
        except IndexError:
            self.log_lines.append("-" * 72 + "\n")
            self.log_lines.append("--> Reading File: {0}\n".format(self.fn))

        # number of bytes in the file
        num_bytes = len(raw_data)

        # beginning index of data blocks
        ds = self._header_len + self._meta_len

        # ----read in header information----------------------------------------
        header_string = raw_data[0 : self._header_len]
        self.read_header(header_string)

        # ---read in meta raw_data----------------------------------------------
        meta_string = raw_data[self._header_len - 1 : ds]
        self.read_metadata(meta_string)

        # ---read in gps raw_data-----------------------------------------------
        # sampling rate times 4 bytes for 32 bit measurement
        df = int(self.df)
        dt = df * 4

        # length of data block plus gps stamp
        block_len = self._stamp_len + dt

        # number of data blocks
        num_blocks = int(np.ceil(num_bytes / float(block_len)))

        # get position of gps stamps
        gps_list = np.zeros(num_blocks, dtype=np.int)

        gps_dict = dict(
            [
                (key, np.zeros(num_blocks, dtype=dtp))
                for key, dtp in zip(self._stamp_list, self._data_types)
            ]
        )
        # make the time array floats instead of ints so can get the decimal
        # place if it isn't 0.
        gps_dict["time"] = gps_dict["time"].astype(np.float32)

        # get gps information from the data
        # get first time stamp that matches the starting time
        s1 = 0
        gps_list[0] = self.get_gps_stamp_location()
        gps_info = np.fromstring(
            raw_data[gps_list[0] : gps_list[0] + self._stamp_len], dtype=self._data_type
        )
        gps_info["time"] = gps_info["time"].astype(np.float32)
        gps_info["time"] = self.get_gps_time(gps_info["time"])[0]
        start_test = self.get_date_time(self.gps_week, gps_info["time"])

        # --> test to make sure the first time corresponds to the scheduled
        # start time
        time_stop = 0
        while (
            start_test != self.start_dt
            and s1 <= self._seconds_diff
            and time_stop <= self._seconds_diff
        ):
            s1 += 1
            gps_list[0] = self.get_gps_stamp_location(gps_list[0] + 7)
            gps_info = np.fromstring(
                raw_data[gps_list[0] : gps_list[0] + self._stamp_len],
                dtype=self._data_type,
            )

            gps_info["time"] = gps_info["time"].astype(np.float32)
            gps_info["time"], gps_dweek = self.get_gps_time(gps_info["time"])

            start_test = self.get_date_time(self.gps_week + gps_dweek, gps_info["time"])
            if s1 == self._seconds_diff:
                s1 = 0
                self.start_dt = self.start_dt[:-2] + "{0:02}".format(
                    int(self.start_dt[-2:]) + 1
                )
                gps_list[0] = self.get_gps_stamp_location()
                time_stop += 1

        # ----Raise an error if the first gps stamp is more than allowed time
        #    difference.
        if time_stop >= self._seconds_diff:
            print (
                "GPS start time is more than "
                + "{0} ".format(self._seconds_diff)
                + "seconds different than scheduled start time of "
                + "{0}. \n ".format(self.start_dt)
                + "Estimated start time is {0} +/- {1} sec".format(
                    start_test, self._seconds_diff
                )
            )

        # put the information into the correct arrays via dictionary
        for jj, key in enumerate(self._stamp_list):
            gps_dict[key][0] = gps_info[0][jj]

        # find the next time stamp
        for ii in range(s1, num_blocks - 1):
            sfind = self.get_gps_stamp_location(gps_list[ii - 1] + 7)
            # make sure it isn't the same time stamp as before
            if sfind != gps_list[ii - 1] and sfind != -1:
                gps_info, gps_index, gps_week = self.get_gps_stamp(sfind)
                gps_list[ii] = gps_index

                if gps_info is not None:
                    for jj, key in enumerate(self._stamp_list):
                        gps_dict[key][ii] = gps_info[0][jj]

        # get only the values that are non zero
        gps_dict["time"] = gps_dict["time"][np.nonzero(gps_dict["time"])]

        num_samples = len(gps_dict["time"])

        # calculate the difference between time stamps
        gps_diff = np.array(
            [
                gps_dict["time"][ii + 1] - gps_dict["time"][ii]
                for ii in range(num_samples - 1)
            ]
        )

        # check for any spots where gps was not locked or mised a sampling interval
        bad_lock = np.where(gps_diff[np.nonzero(gps_diff)] != 1.0)[0]

        if len(bad_lock) > 0:
            for bb in bad_lock:
                if gps_diff[bb] > 5:
                    self.log_lines.append(
                        " " * 4
                        + "point {0:^15},".format(gps_list[bb])
                        + "gps diff {0:^15}\n".format(gps_diff[bb])
                    )

            self.log_lines.append(" " * 4 + "*" * 52 + "\n")

        # need to be sure that the number of data points between time stamps is
        # equal to the sampling rate, if it is not then remove that interval.
        # Most likely it is at the beginning or end of time series.
        dsamples = np.array(
            [
                (gps_list[nn + 1] - gps_list[nn] - self._stamp_len - df * 4) / 4
                for nn in range(num_samples)
            ]
        )

        bad_interval = np.where(abs(dsamples) > self._skip_sample_tolerance)[0]
        bmin = 0
        bmax = num_samples
        if len(bad_interval) > 0:
            # need to locate the bad interval numbers
            for bb in bad_interval:
                if bb <= 10:
                    bmin = bb + 1
                if bb > num_samples - 10:
                    bmax = bb

            gps_list = gps_list[bmin:bmax]

        num_samples = len(gps_list)
        if self.verbose:
            print "Found {0} gps time stamps, ".format(
                num_samples
            ) + "with equal intervals of {0} samples".format(int(self.df))

        self.log_lines.append(
            " " * 4
            + "Found {0} gps time stamps, ".format(num_samples)
            + "with equal intervals of {0} samples\n".format(int(self.df))
        )

        # read in data
        data_array = np.zeros((num_samples + 1) * df, dtype=np.float32)
        for ll, kk in enumerate(gps_list[0:-1]):
            pdiff = ((gps_list[ll + 1] - (kk + self._stamp_len)) - (df * 4)) / 4
            self.sample_diff_list.append(pdiff)
            dblock = raw_data[kk + self._stamp_len : gps_list[ll + 1]]
            try:
                data_array[ll * df : (ll + 1) * df + pdiff] = np.fromstring(
                    dblock, dtype=np.int32
                )
            except ValueError:
                print "samples between time step {0} is off by {1} samples".format(
                    ll, abs(pdiff)
                )

        if sum(self.sample_diff_list) != 0:
            if self.verbose:
                print "time series is off by {0} seconds".format(
                    float(sum(self.sample_diff_list)) / df
                )
                self.log_lines.append(
                    "time series is off by {0} seconds".format(
                        float(sum(self.sample_diff_list)) / df
                    )
                )

        # get only the non-zero data bits, this is dangerous if there is
        # actually an exact 0 in the data, but rarely happens
        self.time_series = data_array[np.nonzero(data_array)]

        # need to cut all the data arrays to have the same length and corresponding
        # data points
        for key in gps_dict.keys():
            gps_dict[key] = gps_dict[key][bmin:bmax]

        # make attributes of imporant information
        self.gps_diff = gps_diff[bmin:bmax]
        self.gps_time = gps_dict["time"]
        self.gps_list = gps_list
        self.temperature = gps_dict["temperature"]
        self.lat = gps_dict["lat"]
        self.lon = gps_dict["lon"]

        self.date_time = np.zeros_like(gps_dict["time"], dtype="|S24")

        for gg, gtime in enumerate(gps_dict["time"]):
            self.date_time[gg] = self.get_date_time(self.gps_week, gtime)

        try:
            self.start_dt = self.date_time[0]
            self.start_date = self.date_time[0].split(",")[0]
            self.start_time = self.date_time[0].split(",")[1]
            if self.verbose:
                print "Starting time of time series is " + "{0} UTC".format(
                    self.date_time[0]
                )
            self.log_lines.append(
                " " * 4
                + "Starting time of time series is "
                + "{0} UTC\n".format(self.date_time[0])
            )
        except IndexError:
            print "No quality data was collected"
            self.log_lines.append(" " * 4 + "No quality data was collected\n")
            self.start_dt = None
            self.start_date = None
            self.start_time = None

        if self.units == "mv":
            self.time_series = self.convert_counts()

    # ==================================================
    def convert_counts(self):
        """
        convert the time series from counts to millivolts

        """

        return self.time_series * self.counts_to_mv_conversion

    # ==================================================
    def convert_mV(self):
        """
        convert millivolts to counts assuming no other scaling has been applied
        
        """

        return self.time_series / self.counts_to_mv_conversion

    # ==================================================
    def compute_schedule_start(self, start_date, start_time, leap_seconds=None):
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
        month_dict = {
            1: 31,
            2: 28,
            3: 31,
            4: 30,
            5: 31,
            6: 30,
            7: 31,
            8: 31,
            9: 30,
            10: 31,
            11: 30,
            12: 31,
        }
        if leap_seconds is not None:
            self._leap_seconds = leap_seconds

        year, month, day = start_date.split("-")

        hour, minutes, seconds = start_time.split(":")

        new_year = int(year)
        new_month = int(month)
        new_day = int(day)
        new_hour = int(hour)
        new_minutes = int(minutes)
        new_seconds = int(seconds) - self._leap_seconds

        if new_seconds < 0:
            new_seconds = (int(seconds) - self._leap_seconds) % 60
            new_minutes = int(minutes) - 1
            if new_minutes < 0:
                new_minutes = (int(minutes) - 1) % 60
                new_hour = int(hour) - 1
                if new_hour < 0:
                    new_hour = (int(hour) - 1) % 24
                    new_day = int(day) - 1
                    if new_day <= 0:
                        new_day = (int(day) - 1) % 30
                        new_month = int(month) - 1
                        if new_month <= 0:
                            new_month = 12 - new_month
                        new_day = month_dict[new_month] - int(day) + 1
                        print "need to check date, have not implemented " + "leap years yet"

        ndate_time = time.strftime(
            datetime_fmt,
            (new_year, new_month, new_day, new_hour, new_minutes, new_seconds, 0, 0, 0),
        )

        return ndate_time

    # ==================================================
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
        if self._raw_data[gps_index + 4] == "\xff":
            gps_index += 1
            if self._raw_data[gps_index + 4] == "\xff":
                gps_index += 1
                if self._raw_data[gps_index + 4] == "\xff":
                    gps_index += 1
                    if self._raw_data[gps_index + 4] == "\xff":
                        gps_index += 1

        return gps_index

    # ==================================================
    def get_gps_stamp(self, gps_index):
        """
        get the gps stamp data
        
        """
        # get numbers from binary format
        try:

            gps_info = np.fromstring(
                self._raw_data[gps_index : gps_index + self._stamp_len],
                dtype=self._data_type,
            )
            while gps_info["time"] < 0:
                gps_index = self.get_gps_stamp_location(start_index=gps_index + 7)
                print "time ", gps_index
                gps_info = np.fromstring(
                    self._raw_data[gps_index : gps_index + self._stamp_len],
                    dtype=self._data_type,
                )

            while gps_info["status"] < 0:
                gps_index = self.get_gps_stamp_location(start_index=gps_index + 7)
                print "status ", gps_index
                gps_info = np.fromstring(
                    self._raw_data[gps_index : gps_index + self._stamp_len],
                    dtype=self._data_type,
                )

            while abs(gps_info["temperature"]) > 80:
                gps_index = self.get_gps_stamp_location(start_index=gps_index + 7)
                print "temperature ", gps_index
                gps_info = np.fromstring(
                    self._raw_data[gps_index : gps_index + self._stamp_len],
                    dtype=self._data_type,
                )

            while abs(gps_info["lat"]) > np.pi:
                gps_index = self.get_gps_stamp_location(start_index=gps_index + 7)
                print "lat ", gps_index
                gps_info = np.fromstring(
                    self._raw_data[gps_index : gps_index + self._stamp_len],
                    dtype=self._data_type,
                )

            # convert lat and lon into decimal degrees
            gps_info["lat"] = self.get_degrees(gps_info["lat"])
            gps_info["lon"] = self.get_degrees(gps_info["lon"])
            gps_info["time"] = gps_info["time"].astype(np.float32)
            gps_info["time"], gps_week = self.get_gps_time(gps_info["time"])

            if gps_info == []:
                print gps_index
                raise ZenGPSError("Something is fucked")
            if gps_index == -1:
                print gps_info
                raise ZenGPSError("Something is fucked")

            return gps_info, gps_index, gps_week

        except ValueError:
            print "Ran into end of file, gps stamp not complete." + " Only {0} points.".format(
                len(self._raw_data[gps_index:])
            )
            return None, gps_index, 0

    # ==================================================
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

        gps_seconds = gps_int / 1024.0

        gps_ms = (gps_seconds - np.floor(gps_int / 1024.0)) * (1.024)

        cc = 0
        if gps_seconds > self._week_len:
            gps_week += 1
            cc = gps_week * self._week_len
            gps_seconds -= self._week_len

        gps_time = np.floor(gps_seconds) + gps_ms + cc

        return gps_time, gps_week

    # ==================================================
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

        # make epoch in seconds, mktime computes local time, need to subtract
        # time zone to get UTC
        epoch_seconds = time.mktime(self._gps_epoch) - time.timezone

        # gps time is 14 seconds ahead of GTC time, but I think that the zen
        # receiver accounts for that so we will leave leap seconds to be 0
        gps_seconds = (
            epoch_seconds + (gps_week * self._week_len) + gps_time - self._leap_seconds
        )

        # compute date and time from seconds
        (year, month, day, hour, minutes, seconds, dow, jday, dls) = time.gmtime(
            gps_seconds
        )

        date_time = time.strftime(
            datetime_fmt,
            (year, month, day, hour, minutes, int(seconds + mseconds), 0, 0, 0),
        )
        return date_time

    # ==================================================
    def get_degrees(self, radian_value):
        """
        convert lat or lon into decimal degrees
        
        """

        degrees = radian_value * 180 / np.pi

        return degrees

    # ==================================================
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

        notches = notch_dict.pop("notches", list(np.arange(60, 2048, 60)))
        notchradius = notch_dict.pop("notchradius", 0.5)
        freqrad = notch_dict.pop("freqrad", 0.5)
        rp = notch_dict.pop("rp", 0.1)
        kwargs = {
            "df": self.df,
            "notches": notches,
            "notchradius": notchradius,
            "freqrad": freqrad,
            "rp": rp,
        }

        self.time_series, self.filt_list = mtfilt.adaptive_notch_filter(
            self.time_series, **kwargs
        )

    # ==================================================
    def write_ascii_mt_file(
        self, save_fn=None, save_station="mb", fmt="%.8e", ex=1, ey=1, notch_dict=None
    ):
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
            svfn_directory = os.path.join(os.path.dirname(self.fn), "TS")
            if not os.path.exists(svfn_directory):
                os.mkdir(svfn_directory)

            svfn_date = "".join(self.start_date.split("-"))
            svfn_time = "".join(self.start_time.split(":"))
            svfn_station = save_station + self.rx_stn
            save_fn = os.path.join(
                svfn_directory,
                "{0}_{1}_{2}_{3}.{4}".format(
                    svfn_station,
                    svfn_date,
                    svfn_time,
                    int(self.df),
                    self.ch_cmp.upper(),
                ),
            )
        # calibrate electric channels
        if self.ch_cmp == "ex":
            time_series /= ex
        elif self.ch_cmp == "ey":
            time_series /= ey

        # apply notch filter if desired
        if notch_dict is not None:
            self.apply_adaptive_notch_filter(notch_dict)
            print "Filtered notches: "
            for nfilt in self.filt_list:
                if type(nfilt[0]) != str:
                    print "{0}{1:.2f} Hz".format(" " * 4, nfilt[0])

        header_tuple = (
            save_station + self.rx_stn,
            self.ch_cmp,
            self.df,
            time.mktime(time.strptime(self.start_dt, datetime_fmt)),
            time_series.shape[0],
            "mV",
            np.median(self.lat),
            np.median(self.lon),
            0.0,
            time_series,
        )

        self.fn_mt_ascii = mtfh.write_ts_file_from_tuple(save_fn, header_tuple, fmt=fmt)

        print "Wrote mtpy timeseries file to {0}".format(self.fn_mt_ascii)

    # ==================================================
    def write_mseed_mt_file(
        self, save_fn=None, save_station="mb", location="Mono Basin", network="USGS"
    ):
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

        svfn_date = "".join(self.start_date.split("-"))
        svfn_time = "".join(self.start_time.split(":"))
        svfn_station = save_station + self.rx_stn
        svfn_chn = self.ch_cmp.upper()
        delta_t = 1.0 / self.df
        t0 = self.start_dt.replace(",", "T")
        if save_fn is None:
            save_fn = os.path.join(
                os.path.dirname(self.fn),
                "{0}_{1}_{2}_{3}_{4}.mseed".format(
                    svfn_station, svfn_date, svfn_time, int(self.df), svfn_chn
                ),
            )

        self.fn_mt_mseed = mtmseed.writefile_obspy_singletrace(
            save_fn, svfn_station, svfn_chn, network, location, delta_t, t0, time_series
        )
        return save_fn

    # ==================================================
    def plot_time_series(self, fig_num=1):
        """
        plots the time series
        """

        time_series = self.convert_counts()
        fig = plt.figure(fig_num, dpi=300)
        ax = fig.add_subplot(1, 1, 1)
        ax.plot(time_series)

        # ax.xaxis.set_minor_locator(MultipleLocator(self.df))
        # ax.xaxis.set_major_locator(MultipleLocator(self.df*15))
        # ax.xaxis.set_ticklabels([self.date_time[ii]
        #                        for ii in range(0,len(self.date_time), 15)])

        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Amplitude (mV)")
        plt.show()

        self.convert_mV()
        return fig, ax

    # ==================================================
    def plot_spectrogram(
        self,
        time_window=2 ** 8,
        time_step=2 ** 6,
        s_window=11,
        frequency_window=1,
        n_freq_bins=2 ** 9,
        sigma_L=None,
    ):
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

        kwargs = {
            "nh": time_window,
            "tstep": time_step,
            "L": s_window,
            "ng": frequency_window,
            "df": self.df,
            "nfbins": n_freq_bins,
            "sigmaL": sigma_L,
        }
        ptf = plotspectrogram.PlotTF(time_series, **kwargs)

        return ptf

    # ==================================================
    def plot_spectra(self, fig_num=2):
        """
        plot the spectra of time series
        """
        if self.time_series is None:
            self.read_3d()

        time_series = self.convert_counts()

        spect = np.fft.fft(mtfilt.zero_pad(time_series))
        plot_freq = np.fft.fftfreq(spect.shape[0], 1.0 / self.df)

        fig = plt.figure(fig_num, [4, 4], dpi=200)
        ax = fig.add_subplot(1, 1, 1)
        ax.loglog(plot_freq, abs(spect) ** 2, lw=0.5)
        ax.grid(which="both", lw=0.25)

        ax.set_xlabel("Frequency (Hz)")
        # ax.set_xlim(1./plot_freq.max(), 1./plot_freq.min())
        ax.set_ylabel("Amplitude")

        plt.show()

        return fig, ax
