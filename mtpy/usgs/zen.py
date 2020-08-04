  # -*- coding: utf-8 -*-
"""
====================
Zen
====================
    * Tools for reading and writing files for Zen and processing software
    * Tools for copying data from SD cards
    * Tools for copying schedules to SD cards
Created on Tue Jun 11 10:53:23 2013
@author: jpeacock-pr
"""

#==============================================================================
#from __future__ import unicode_literals

import time
import datetime
import dateutil.parser
import os
import struct
import string
import shutil
import numpy as np

import mtpy.imaging.plotspectrogram as plotspectrogram
import mtpy.core.ts as mtts

try:
    import win32api
except ImportError:
    print("WARNING: Cannot find win32api, will not be able to detect"+
          " drive names")

#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'
datetime_sec = '%Y-%m-%d %H:%M:%S.%f'
#==============================================================================
# =============================================================================
#  Get leap seconds
# =============================================================================
def calculate_leap_seconds(year, month, day):
    """
    get the leap seconds for the given year to convert GPS time to UTC time
    .. note:: GPS time started in 1980
    .. note:: GPS time is leap seconds ahead of UTC time, therefore you
              should subtract leap seconds from GPS time to get UTC time.
    =========================== ===============================================
    Date Range                  Leap Seconds
    =========================== ===============================================
    1981-07-01 - 1982-07-01     1
    1982-07-01 - 1983-07-01     2
    1983-07-01 - 1985-07-01     3
    1985-07-01 - 1988-01-01     4
    1988-01-01 - 1990-01-01     5
    1990-01-01 - 1991-01-01     6
    1991-01-01 - 1992-07-01     7
    1992-07-01 - 1993-07-01     8
    1993-07-01 - 1994-07-01     9
    1994-07-01 - 1996-01-01     10
    1996-01-01 - 1997-07-01     11
    1997-07-01 - 1999-01-01     12
    1999-01-01 - 2006-01-01     13
    2006-01-01 - 2009-01-01     14
    2009-01-01 - 2012-07-01     15
    2012-07-01 - 2015-07-01     16
    2015-07-01 - 2017-01-01        17
    2017-01-01 - ????-??-??        18
    =========================== ===============================================
    """
    leap_second_dict = {0: {'min':datetime.date(1980, 1, 1),
                            'max':datetime.date(1981, 7, 1)},
                        1: {'min':datetime.date(1981, 7, 1),
                            'max':datetime.date(1982, 7, 1)},
                        2: {'min':datetime.date(1982, 7, 1),
                            'max':datetime.date(1983, 7, 1)},
                        3: {'min':datetime.date(1983, 7, 1),
                            'max':datetime.date(1985, 7, 1)},
                        4: {'min':datetime.date(1985, 7, 1),
                            'max':datetime.date(1988, 1, 1)},
                        5: {'min':datetime.date(1988, 1, 1),
                            'max':datetime.date(1990, 1, 1)},
                        6: {'min':datetime.date(1990, 1, 1),
                            'max':datetime.date(1991, 1, 1)},
                        7: {'min':datetime.date(1991, 1, 1),
                            'max':datetime.date(1992, 7, 1)},
                        8: {'min':datetime.date(1992, 7, 1),
                            'max':datetime.date(1993, 7, 1)},
                        9: {'min':datetime.date(1993, 7, 1),
                            'max':datetime.date(1994, 7, 1)},
                        10: {'min':datetime.date(1994, 7, 1),
                             'max':datetime.date(1996, 1, 1)},
                        11: {'min':datetime.date(1996, 1, 1),
                             'max':datetime.date(1997, 7, 1)},
                        12: {'min':datetime.date(1997, 7, 1),
                             'max':datetime.date(1999, 1, 1)},
                        13: {'min':datetime.date(1999, 1, 1),
                             'max':datetime.date(2006, 1, 1)},
                        14: {'min':datetime.date(2006, 1, 1),
                             'max':datetime.date(2009, 1, 1)},
                        15: {'min':datetime.date(2009, 1, 1),
                             'max':datetime.date(2012, 6, 30)},
                        16: {'min':datetime.date(2012, 6, 30),
                             'max':datetime.date(2015, 6, 30)},
                        17: {'min':datetime.date(2015, 6, 30),
                             'max':datetime.date(2016, 12, 31)},
                        18: {'min':datetime.date(2016, 12, 31),
                             'max':datetime.date(2020, 12, 1)}}

    year = int(year)
    month = int(month)
    day = int(day)

    # make the date a datetime object, easier to test
    given_date = datetime.date(year, month, day)

    # made an executive decision that the date can be equal to the min, but
    # no the max, otherwise get an error.
    for leap_key in sorted(leap_second_dict.keys()):
        if given_date < leap_second_dict[leap_key]['max'] and\
           given_date >= leap_second_dict[leap_key]['min']:
            return int(leap_key)

    return None
#==============================================================================
class Z3DHeader(object):
    """
    class for z3d header.  This will read in the header information of a
    Z3D file and make each metadata entry an attirbute
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
        **fid** : file object
                  ie. open(Z3Dfile, 'rb')
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
        >>> Z3Dfn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3DHeader()
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
        self.old_version = False

        for key in kwargs:
            setattr(self, key, kwargs[key])

    @property
    def data_logger(self):
        return 'ZEN{0:03}'.format(int(self.box_number))

    def read_header(self, fn=None, fid=None):
        """
        read in the header string
        """
        if fn is not None:
            self.fn = fn

        if fid is not None:
            self.fid = fid

        if self.fn is None and self.fid is None:
            print('WARNING: No file to read')
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(0)
                self.header_str = self.fid.read(self._header_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = open(self.fn, 'rb')
                self.header_str = self.fid.read(self._header_len)
            else:
                self.fid.seek(0)
                self.header_str = self.fid.read(self._header_len)

        header_list = self.header_str.split(b'\n')
        for h_str in header_list:
            h_str = h_str.decode()
            if h_str.find('=') > 0:
                h_list = h_str.split('=')
                h_key = h_list[0].strip().lower()
                h_key = h_key.replace(' ', '_').replace('/', '').replace('.', '_')
                h_value = self.convert_value(h_key, h_list[1].strip())
                setattr(self, h_key, h_value)
            elif len(h_str) == 0:
                continue
            # need to adjust for older versions of z3d files
            elif h_str.count(',') > 1:
                self.old_version = True
                if h_str.find('Schedule') >= 0:
                    h_str = h_str.replace(',', 'T', 1)
                for hh in h_str.split(','):
                    if hh.find(';') > 0:
                        m_key, m_value = hh.split(';')[1].split(':')

                    elif len(hh.split(':', 1)) == 2:
                        m_key, m_value = hh.split(':', 1)
                    else:
                         print(hh)

                    m_key = m_key.strip().lower().replace(' ', '_').replace('/', '').replace('.', '_')
                    m_value = self.convert_value(m_key, m_value.strip())
                    setattr(self, m_key, m_value)

    def convert_value(self, key_string, value_string):
        """
        convert the value to the appropriate units given the key
        """

        try:
            return_value = float(value_string)
        except ValueError:
            return_value = value_string

        if key_string.lower() in ['lat', 'lon', 'long']:
            return_value = np.rad2deg(float(value_string))
            if 'lat' in key_string.lower():
                if abs(return_value) > 90:
                    return_value = 0.0
            elif 'lon' in key_string.lower():
                if abs(return_value) > 180:
                    return_value = 0.0

        return return_value

#==============================================================================
# meta data
#==============================================================================
class Z3DSchedule(object):
    """
    class object for metadata of Z3d file.  This will read in the schedule
    information of a Z3D file and make each metadata entry an attirbute.
    The attributes are left in capitalization of the Z3D file.
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
        **fid** : file object
                  ie. open(Z3Dfile, 'rb')
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
    read_schedule   read in the schedule information from the given
                             file
    ======================== ==================================================
    Example
    --------------
        >>> import mtpy.usgs.zen as zen
        >>> Z3Dfn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3DSchedule()
        >>> header_obj.read_schedule()
    """
    def __init__(self, fn=None, fid=None, **kwargs):
        self.fn = fn
        self.fid = fid
        self.meta_string = None

        self._schedule_metadata_len = 512
        self._header_len = 512

        self.AutoGain = None
        self.Comment = None
        self.Date = None
        self.Duty = None
        self.FFTStacks = None
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
        self.datetime = None

        for key in kwargs:
            setattr(self, key, kwargs[key])


    def read_schedule(self, fn=None, fid=None):
        """
        read meta data string
        """
        if fn is not None:
            self.fn = fn

        if fid is not None:
            self.fid = fid

        if self.fn is None and self.fid is None:
            print('WARNING: No file to read')
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = open(self.fn, 'rb')
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)
            else:
                self.fid.seek(self._header_len)
                self.meta_string = self.fid.read(self._header_len)

        meta_list = self.meta_string.split(b'\n')
        for m_str in meta_list:
            m_str = m_str.decode()
            if m_str.find('=') > 0:
                m_list = m_str.split('=')
                m_key = m_list[0].split('.')[1].strip()
                m_key = m_key.replace('/', '')
                m_value = m_list[1].strip()
                setattr(self, m_key, m_value)

        # the first good GPS stamp is on the 3rd, so need to add 2 seconds
        try:
            self.Time = '{0}{1:02}'.format(self.Time[0:6],
                                           int(self.Time[6:])+2)
        except TypeError:
            return

        self.datetime = datetime.datetime.strptime('{0},{1}'.format(self.Date,
                                                                    self.Time),
                                                   datetime_fmt)

#==============================================================================
#  Meta data class
#==============================================================================
class Z3DMetadata(object):
    """
    class object for metadata of Z3d file.  This will read in the metadata
    information of a Z3D file and make each metadata entry an attirbute.
    The attributes are left in capitalization of the Z3D file.
    Arguments
    ------------
        **fn** : string
                 full path to Z3D file
        **fid** : file object
                  ie. open(Z3Dfile, 'rb')
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
        >>> Z3Dfn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
        >>> header_obj = zen.Z3DMetadata()
        >>> header_obj.read_metadata()
    """

    def __init__(self, fn=None, fid=None, **kwargs):
        self.fn = fn
        self.fid = fid
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
        self.line_name = None
        self.survey_type = None
        self.unit_length = None
        self.station = None
        self.count = 0
        self.notes = None

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
            print('WARNING: No file to read')
        elif self.fn is None:
            if self.fid is not None:
                self.fid.seek(self._header_length+self._schedule_metadata_len)
        elif self.fn is not None:
            if self.fid is None:
                self.fid = open(self.fn, 'rb')
                self.fid.seek(self._header_length+self._schedule_metadata_len)
            else:
                self.fid.seek(self._header_length+self._schedule_metadata_len)

        # read in calibration and meta data
        self.find_metadata = True
        self.board_cal = []
        self.coil_cal = []
        self.count = 0
        m_list = []
        cal_find = False
        while self.find_metadata == True:
            try:
                test_str = self.fid.read(self._metadata_length).decode().lower()
            except UnicodeDecodeError:
                self.find_metadata = False
                break
            
            if 'metadata' in test_str:
                self.count += 1
                test_str = test_str.strip().split('record')[1].strip()
                if test_str.count('|') > 1:
                    for t_str in test_str.split('|'):
                        # get metadata name and value
                        if t_str.find('=') == -1 and \
                           t_str.lower().find('line.name') == -1:
                            # get metadata for older versions of z3d files
                            if len(t_str.split(',')) == 2:
                                t_list = t_str.lower().split(',')
                                t_key = t_list[0].strip().replace('.', '_')
                                if t_key == 'ch_varasp':
                                    t_key = 'ch_length'
                                t_value = t_list[1].strip()
                                setattr(self, t_key, t_value)
                            if t_str.count(' ') > 1:
                                self.notes = t_str
                        # get metadata for just the line that has line name
                        # because for some reason that is still comma separated
                        elif t_str.lower().find('line.name') >= 0:
                            t_list = t_str.split(',')
                            t_key = t_list[0].strip().replace('.', '_')
                            t_value = t_list[1].strip()
                            setattr(self, t_key.lower(), t_value)
                        # get metadata for newer z3d files
                        else:
                            t_list = t_str.split('=')
                            t_key = t_list[0].strip().replace('.', '_')
                            t_value = t_list[1].strip()
                            setattr(self, t_key.lower(), t_value)

                elif 'cal.brd' in test_str:
                    t_list = test_str.split(',')
                    t_key = t_list[0].strip().replace('.', '_')
                    setattr(self, t_key.lower(), t_list[1])
                    for t_str in t_list[2:]:
                        t_str = t_str.replace('\x00', '').replace('|', '')
                        try:
                            self.board_cal.append([float(tt.strip())
                                                   for tt in t_str.strip().split(':')])
                        except ValueError:
                            self.board_cal.append([tt.strip()
                                                   for tt in t_str.strip().split(':')])
                # some times the coil calibration does not start on its own line
                # so need to parse the line up and I'm not sure what the calibration
                # version is for so I have named it odd
                elif 'cal.ant' in test_str:
                    # check to see if the coil calibration exists
                    cal_find = True
                    test_list = test_str.split(',')
                    coil_num = test_list[1].split('|')[1]
                    coil_key, coil_value = coil_num.split('=')
                    setattr(self, coil_key.replace('.', '_').lower(),
                            coil_value.strip())
                    for t_str in test_list[2:]:
                        if '\x00' in t_str:
                            break
                        self.coil_cal.append([float(tt.strip())
                                              for tt in t_str.split(':')])
                elif cal_find and self.count > 3:
                    t_list = test_str.split(',')
                    for t_str in t_list:
                        if '\x00' in t_str:
                            break
                        else:
                            self.coil_cal.append([float(tt.strip())
                                                  for tt in t_str.strip().split(':')])

            else:
                self.find_metadata = False
                # need to go back to where the meta data was found so
                # we don't skip a gps time stamp
                self.m_tell = self.fid.tell() - self._metadata_length

        # make coil calibration and board calibration structured arrays
        if len(self.coil_cal) > 0:
            self.coil_cal = np.core.records.fromrecords(self.coil_cal,
                                                        names='frequency, amplitude, phase')
        if len(self.board_cal) > 0:
            try:
                self.board_cal = np.core.records.fromrecords(self.board_cal,
                                                             names='frequency, rate, amplitude, phase')
            except ValueError:
                self.board_cal = None

        try:
            self.station = '{0}{1}'.format(self.line_name,
                                           self.rx_xyz0.split(':')[0])
        except AttributeError:
            if hasattr(self, 'rx_stn'):
                self.station = f"{self.rx_stn}"
            elif hasattr(self, 'ch_stn'):
                self.station = f"{self.ch_stn}"
            else:
                self.station = None
                print("WARNING: Need to input station name")

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
    header                   Z3DHeader object                Z3DHeader
    metadata                 Z3DMetadata                     Z3DMetadata
    schedule                 Z3DSchedule            Z3DSchedule
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

        self.header = Z3DHeader(fn)
        self.schedule = Z3DSchedule(fn)
        self.metadata = Z3DMetadata(fn)

        self._gps_stamp_length = kwargs.pop('stamp_len', 64)
        self._gps_bytes = self._gps_stamp_length/4

        self.gps_stamps = None

        self._gps_flag_0 = np.int32(2147483647)
        self._gps_flag_1 = np.int32(-2147483648)
        self._gps_f0 = self._gps_flag_0.tostring()
        self._gps_f1 = self._gps_flag_1.tostring()
        self.gps_flag = self._gps_f0 + self._gps_f1

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
        self._leap_seconds = 18
        self._block_len = 2**16
        # the number in the cac files is for volts, we want mV
        self._counts_to_mv_conversion = 9.5367431640625e-10 * 1E3
        self.num_sec_to_skip = 3

        self.units = 'counts'
        self.df = None

        self.ts_obj = mtts.MTTS()

    @property
    def station(self):
        """
        station name
        """
        return self.metadata.station

    @station.setter
    def station(self, station):
        """
        station name
        """
        self.metadata.station = station

    @property
    def dipole_len(self):
        """
        dipole length
        """
        if self.metadata.ch_length is not None:
            return self.metadata.ch_length
        elif hasattr(self.metadata, 'ch_offset_xyz1'):
            # only ex and ey have xyz2
            if hasattr(self.metadata, 'ch_offset_xyz2'):
                x1, y1, z1 = [float(offset) for offset in
                              self.metadata.ch_offset_xyz1.split(':')]
                x2, y2, z2 = [float(offset) for offset in
                              self.metadata.ch_offset_xyz2.split(':')]
                length = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                return np.round(length, 2)
            else:
                return 0
        elif self.metadata.ch_xyz1 is not None:
            x1, y1 = [float(d) for d in self.metadata.ch_xyz1.split(':')]
            x2, y2 = [float(d) for d in self.metadata.ch_xyz2.split(':')]
            length = np.sqrt((x2-x1)**2+(y2-y1)**2)*100.
            return np.round(length, 2)

    @property
    def azimuth(self):
        """
        azimuth of instrument setup
        """
        if self.metadata.ch_azimuth is not None:
            return float(self.metadata.ch_azimuth)
        elif self.metadata.rx_xazimuth is not None:
            return float(self.metadata.rx_xazimuth)
        else:
            return None

    @property
    def component(self):
        """
        channel
        """
        return self.metadata.ch_cmp.lower()

    @property
    def lat(self):
        """
        latitude in decimal degrees
        """
        return self.header.lat

    @property
    def lon(self):
        """
        longitude in decimal degrees
        """
        return self.header.long

    @property
    def elev(self):
        """
        elevation in meters
        """
        return self.header.alt

    @property
    def df(self):
        """
        sampling rate
        """
        return self.header.ad_rate

    @df.setter
    def df(self, sampling_rate):
        """
        sampling rate
        """
        if sampling_rate is not None:
            self.header.ad_rate = float(sampling_rate)

    @property
    def zen_schedule(self):
        """
        zen schedule data and time
        """

        if self.header.old_version is True:
            dt_str = self.header.schedule.replace('T', ',')
            self.schedule.Date = dt_str.split(',')[0]
            self.schedule.Time = dt_str.split(',')[1]
            # the first good GPS stamp is on the 3rd, so need to add 2 seconds
            self.schedule.Time = '{0}{1:02}'.format(self.schedule.Time[0:6],
                                                    int(self.schedule.Time[6:])+2)

            self.schedule.datetime = datetime.datetime.strptime('{0},{1}'.format(self.schedule.Date,
                                                                                 self.schedule.Time),
                                                                datetime_fmt)

        return self.schedule.datetime

    @zen_schedule.setter
    def zen_schedule(self, schedule_dt):
        """
        on setting set schedule datetime
        """
        if type(schedule_dt) is not datetime.datetime:
            raise TypeError('New schedule datetime must be type datetime.datetime')
        self.schedule.datetime = schedule_dt

        # set the leap seconds
        self._leap_seconds = calculate_leap_seconds(self.schedule.datetime.year,
                                                    self.schedule.datetime.month,
                                                    self.schedule.datetime.day)

    @property
    def coil_num(self):
        """
        coil number
        """
        if self.metadata.cal_ant is not None:
            return self.metadata.cal_ant
        elif self.metadata.ch_number is not None:
            return self.metadata.ch_number
        else:
            return None

    def _get_gps_stamp_type(self, old_version=False):
        """
        get the correct stamp type.
        Older versions the stamp length was 36 bits
        New versions have a 64 bit stamp
        """

        if old_version is True:
            self._gps_dtype = np.dtype([('gps', np.int32),
                                        ('time', np.int32),
                                        ('lat', np.float64),
                                        ('lon', np.float64),
                                        ('block_len', np.int32),
                                        ('gps_accuracy', np.int32),
                                        ('temperature', np.float32)])
            self._gps_stamp_length = 36
            self._gps_bytes = self._gps_stamp_length/4
            self._gps_flag_0 = -1
            self._block_len = int(self._gps_stamp_length+self.df*4)
            self.gps_flag = self._gps_f0

        else:
            return

    #======================================
    def _read_header(self, fn=None, fid=None):
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
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_header(fn)
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> Z3Dfid = open(fn, 'rb')
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_header(fid=Z3Dfid)
        """

        if fn is not None:
            self.fn = fn

        self.header.read_header(fn=self.fn, fid=fid)
        if self.header.old_version:
            if self.header.box_number is None:
                self.header.box_number = '6666'

    #======================================
    def _read_schedule(self, fn=None, fid=None):
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
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_schedule(fn)
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> Z3Dfid = open(fn, 'rb')
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_schedule(fid=Z3Dfid)
        """

        if fn is not None:
            self.fn = fn

        self.schedule.read_schedule(fn=self.fn, fid=fid)
        if self.header.old_version:
            dt_str = self.header.schedule.replace('T', ',')
            self.schedule.Date = dt_str.split(',')[0]
            self.schedule.Time = dt_str.split(',')[1]
            year, month, day = [int(dd) for dd in self.schedule.Date.split('-')]
            hour, minute, second = [int(dd) for dd in self.schedule.Time.split(':')]
            self.schedule.datetime = datetime.datetime(year, month, day,
                                                       hour, minute, second)

        # set the leap seconds
        self._leap_seconds = calculate_leap_seconds(self.schedule.datetime.year,
                                                    self.schedule.datetime.month,
                                                    self.schedule.datetime.day)

    #======================================
    def _read_metadata(self, fn=None, fid=None):
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
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_metadata(fn)
        Example with file object
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> Z3Dfid = open(fn, 'rb')
            >>> Z3Dobj = zen.Zen3D()
            >>> Z3Dobj.read_metadata(fid=Z3Dfid)
        """

        if fn is not None:
            self.fn = fn

        if self.header.old_version:
            self.metadata._schedule_metadata_len = 0
        self.metadata.read_metadata(fn=self.fn, fid=fid)

    #=====================================
    def read_all_info(self):
        """
        Read header, schedule, and metadata
        """
        with open(self.fn, 'rb') as file_id:

            self._read_header(fid=file_id)
            self._read_schedule(fid=file_id)
            self._read_metadata(fid=file_id)

    #======================================
    def read_z3d(self, Z3Dfn=None):
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
        if Z3Dfn is not None:
            self.fn = Z3Dfn

        #print(u'------- Reading {0} ---------'.format(self.fn))
        st = time.time()

        #get the file size to get an estimate of how many data points there are
        file_size = os.path.getsize(self.fn)

        # using the with statement works in Python versions 2.7 or higher
        # the added benefit of the with statement is that it will close the
        # file object upon reading completion.
        with open(self.fn, 'rb') as file_id:

            self._read_header(fid=file_id)
            self._read_schedule(fid=file_id)
            self._read_metadata(fid=file_id)

            if self.header.old_version is True:
                self._get_gps_stamp_type(True)

            # move the read value to where the end of the metadata is
            file_id.seek(self.metadata.m_tell)

            # initalize a data array filled with zeros, everything goes into
            # this array then we parse later
            data = np.zeros(int((file_size-512*(1+self.metadata.count))/4 + 8*self.df),
                            dtype=np.int32)
            # go over a while loop until the data cound exceed the file size
            data_count = 0
            while True:
                # need to make sure the last block read is a multiple of 32 bit
                read_len = min([self._block_len,
                                int(32*((file_size-file_id.tell())//32))])
                test_str = np.fromstring(file_id.read(read_len),
                                         dtype=np.int32)
                if len(test_str) == 0:
                    break
                data[data_count:data_count+len(test_str)] = test_str
                data_count += test_str.size

        self.raw_data = data.copy()
        # find the gps stamps
        gps_stamp_find = self.get_gps_stamp_index(data, self.header.old_version)

        # skip the first two stamps and trim data
        try:
            data = data[gps_stamp_find[self.num_sec_to_skip]:]
        except IndexError:
            raise ZenGPSError("Data is bad, cannot open file {0}".format(self.fn))

        # find gps stamps of the trimmed data
        gps_stamp_find = self.get_gps_stamp_index(data, self.header.old_version)

        self.gps_stamps = np.zeros(len(gps_stamp_find), dtype=self._gps_dtype)

        for ii, gps_find in enumerate(gps_stamp_find):
            try:
                data[gps_find+1]
            except IndexError:
                pass
                print('***Failed gps stamp***')
                print('    stamp {0} out of {1}'.format(ii+1,
                                                        len(gps_stamp_find)))
                break

            if self.header.old_version is True or data[gps_find+1] == self._gps_flag_1:
                gps_str = struct.pack('<'+'i'*int(self._gps_bytes),
                                      *data[int(gps_find):int(gps_find+self._gps_bytes)])
                self.gps_stamps[ii] = np.fromstring(gps_str,
                                                    dtype=self._gps_dtype)
                if ii > 0:
                    self.gps_stamps[ii]['block_len'] = gps_find-\
                                           gps_stamp_find[ii-1]-self._gps_bytes
                elif ii == 0:
                    self.gps_stamps[ii]['block_len'] = 0
                data[int(gps_find):int(gps_find+self._gps_bytes)] = 0

        # fill the time series object
        self._fill_ts_obj(data[np.nonzero(data)])

        print('    found {0} GPS time stamps'.format(self.gps_stamps.shape[0]))
        print('    found {0} data points'.format(self.ts_obj.ts.data.size))

        # time it
        et = time.time()
        print('INFO: --> Reading data took: {0:.3f} seconds'.format(et-st))

    #=================================================
    def _fill_ts_obj(self, ts_data):
        """
        fill time series object
        """
        # fill the time series object
        self.ts_obj = mtts.MTTS()
        self.ts_obj.ts = ts_data

        # convert data to mV
        self.convert_counts_to_mv()
        self.ts_obj.ts = self.ts_obj.ts.astype(np.float32)

        self.validate_time_blocks()
        self.convert_gps_time()
        self.zen_schedule = self.check_start_time()

        # fill time series object metadata
        self.ts_obj.station = self.station
        self.ts_obj.sampling_rate = float(self.df)
        self.ts_obj.start_time_utc = self.zen_schedule.isoformat()
        self.ts_obj.component = self.component
        self.ts_obj.coordinate_system = 'geomagnetic'
        try:
            self.ts_obj.dipole_length = float(self.dipole_len)
        except TypeError:
            self.ts_obj.dipole_length = -666
        try:
            self.ts_obj.azimuth = float(self.azimuth)
        except TypeError:
            self.ts_obj.azimuth = -666
        self.ts_obj.units = 'mV'
        self.ts_obj.lat = self.lat
        self.ts_obj.lon = self.lon
        self.ts_obj.datum = 'WGS84'
        self.ts_obj.data_logger = self.header.data_logger
        self.ts_obj.elev = self.elev
        self.ts_obj.instrument_id = self.coil_num
        self.ts_obj.calibration_fn = None
        self.ts_obj.declination = 0.0
        self.ts_obj.conversion = self._counts_to_mv_conversion
        self.ts_obj.gain = self.header.ad_gain
        self.ts_obj.channel_number = int(self.header.channel)
        self.ts_obj.fn = os.path.basename(self.fn)

    # =================================================
    def get_gps_stamp_index(self, ts_data, old_version=False):
        """
        locate the time stamps in a given time series.

        Looks for gps_flag_0 first, if the file is newer, then makes sure the
        next value is gps_flag_1

        :returns: list of gps stamps indicies
        """

        # find the gps stamps
        gps_stamp_find = np.where(ts_data == self._gps_flag_0)[0]

        if old_version is False:
            gps_stamp_find = [gps_find for gps_find in gps_stamp_find
                              if ts_data[gps_find+1] == self._gps_flag_1]

        return gps_stamp_find


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
        zen_start_utc = self.get_UTC_date_time(self.header.gpsweek,
                                               self.gps_stamps['time'][0])

        # calculate the scheduled start time
        s_start = '{0},{1}'.format(self.schedule.Date, self.schedule.Time)
        schedule_time = datetime.datetime.strptime(s_start, datetime_fmt)

        # reset the data and time in the schedule meta data so there is no
        # confusion on when the time series starts
        self.schedule.Date = zen_start_utc.strftime(u'%Y-%m-%d')
        self.schedule.Time = zen_start_utc.strftime(u'%H:%M:%S')

        # estimate the time difference between the two
        time_diff = self.zen_schedule - schedule_time
        print('    Scheduled time was {0} (GPS time)'.format(s_start))
        print('    1st good stamp was {0} (UTC time)'.format(zen_start_utc.isoformat()))
        print('    difference of {0:.2f} seconds'.format(time_diff.total_seconds()))

        return zen_start_utc

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
            print('-'*50)
            for bb in bad_times:
                print('WARNING: bad time at index {0} > 0.5 s'.format(bb))

    #===================================================
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

                print('WARNING: Skipped the first {0} seconds'.format(
                    bad_blocks[-1]))
                print('WARNING: Skipped first {0} poins in time series'.format(
                                                                      ts_skip))

    #==================================================
    def convert_gps_time(self):
        """
        convert gps time integer to relative seconds from gps_week
        """
        # need to convert gps_time to type float from int
        dt = self._gps_dtype.descr
        if self.header.old_version is True:
            dt[1] = ('time', np.float32)
        else:
            dt[2] = ('time', np.float32)
        self.gps_stamps = self.gps_stamps.astype(np.dtype(dt))

        # convert to seconds
        # these are seconds relative to the gps week
        time_conv = self.gps_stamps['time'].copy()/1024.
        time_ms = (time_conv-np.floor(time_conv))*1.024
        time_conv = np.floor(time_conv)+time_ms

        self.gps_stamps['time'][:] = time_conv

    #==================================================
    def convert_counts_to_mv(self):
        """
        convert the time series from counts to millivolts
        """

        self.ts_obj.ts.data *= self._counts_to_mv_conversion

    #==================================================
    def convert_mv_to_counts(self):
        """
        convert millivolts to counts assuming no other scaling has been applied
        """

        self.ts_obj.ts.data /= self._counts_to_mv_conversion
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
        .. note:: GPS time is curently ahead by 18 (after 2016) seconds from
        UTC time.
        Arguments
        -------------
            **gps_week**: int
                          integer value of gps_week that the data was collected
            **gps_time**: int
                          number of seconds from beginning of gps_week
            **leap_seconds**: int
                              number of seconds gps time is off from UTC time.
                              It is currently off by 18 seconds (after 2016).
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

        #mseconds = gps_time % 1

        #make epoch in seconds, mktime computes local time, need to subtract
        #time zone to get UTC
        epoch_seconds = time.mktime(self._gps_epoch)-time.timezone

        #gps time is 18 seconds ahead of GTC time
        utc_seconds = epoch_seconds+(gps_week*self._week_len)+gps_time-\
                                                        self._leap_seconds

        # compute date and time from seconds and return a datetime object
        # easier to manipulate later
        date_time = datetime.datetime.utcfromtimestamp(utc_seconds)

        return date_time

    #==================================================
    def apply_adaptive_notch_filter(self, notch_dict={'notches':np.arange(60, 1860, 60),
                                                      'notch_radius':0.5,
                                                      'freq_rad':0.5,
                                                      'rp':0.1}):
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
            notch_dict['notches']
        except KeyError:
            return
        try:
            self.ts_obj.ts.data
        except AttributeError:
            self.read_z3d()

        self.ts_obj.ts = self.ts_obj.apply_addaptive_notch_filter(**notch_dict)

    #==================================================
    def write_ascii_mt_file(self, save_fn=None, fmt='%.8e', notch_dict=None,
                            dec=1):
        """
        write an mtpy time series data file
        Arguments
        -------------
            **save_fn** : full path to save file, if None file is saved as:
                          station_YYYYMMDD_hhmmss_df.component
                          ex. mt01_20130206_120000_256.HX
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
            **dec** : int
                      decimation factor
                      *default* is 1
        Output
        -------------
            **fn_mt_ascii** : full path to saved file
        Example
        ------------
            >>> import mtpy.usgs.zen as zen
            >>> fn = r"/home/mt/mt01/mt01_20150522_080000_256_EX.Z3D"
            >>> Z3Dobj = zen.Zen3D(fn)
            >>> asc_fn = Z3Dobj.write_ascii_mt_file(save_station='mt', notch_dict={})
        """
        if self.station is None:
            self.read_all_info()

        if dec > 1:
            print('INFO: Decimating data by factor of {0}'.format(dec))
            self.df = self.df/dec

        # make a new file name to save to that includes the meta information
        if save_fn is None:
            svfn_directory = os.path.join(os.path.dirname(self.fn), 'TS')
            if not os.path.exists(svfn_directory):
                os.mkdir(svfn_directory)

            svfn_date = ''.join(self.schedule.Date.split('-'))
            svfn_time = ''.join(self.schedule.Time.split(':'))
            self.fn_mt_ascii = os.path.join(svfn_directory,
                                            '{0}_{1}_{2}_{3}.{4}'.format(self.station,
                                                                         svfn_date,
                                                                         svfn_time,
                                                                         int(self.df),
                                                                         self.metadata.ch_cmp.upper()))
        else:
            self.fn_mt_ascii = save_fn
        # if the file already exists skip it
        if os.path.isfile(self.fn_mt_ascii) == True:
            print('\t************')
            print('\tmtpy file already exists for {0} --> {1}'.format(self.fn,
                                                                    self.fn_mt_ascii))
            print('\tskipping')
            print('\t************')
            # if there is a decimation factor need to read in the time
            # series data to get the length.
            c = self.ts_obj.read_ascii_header(self.fn_mt_ascii)
            self.zen_schedule = dateutil.parser.parse(self.ts_obj.start_time_utc)

            return

        # read in time series data if haven't yet.
        if len(self.ts_obj.ts) <=1:
            self.read_z3d()

        # decimate the data.  try resample at first, see how that goes
        # make the attribute time series equal to the decimated data.
        if dec > 1:
            self.ts_obj.decimate(dec)

        # apply notch filter if desired
        if notch_dict is not None:
            self.apply_adaptive_notch_filter(notch_dict)

        # convert counts to mV and scale accordingly
        # self.convert_counts() #--> data is already converted to mV
        # calibrate electric channels should be in mV/km
        if self.component in ['ex', 'ey']:
            e_scale = float(self.dipole_len)
            self.ts_obj.ts.data /= e_scale/1000.
            print('INFO: Using scales {0} = {1} m'.format(self.metadata.ch_cmp.upper(),
                                                    e_scale))
            self.ts_obj.units = 'mV/km'

        self.ts_obj.write_ascii_file(fn_ascii=self.fn_mt_ascii)

        print('INFO: Wrote mtpy timeseries file to {0}'.format(self.fn_mt_ascii))

    #==================================================
    def plot_time_series(self, fig_num=1):
        """
        plots the time series
        """

        self.ts_obj.ts.plot(x_compat=True)

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

        kwargs = {'nh':time_window, 'tstep':time_step, 'L':s_window,
                  'ng':frequency_window, 'df':self.df, 'nfbins':n_freq_bins,
                  'sigmaL': sigma_L}
        ptf = plotspectrogram.PlotTF(self.ts_obj.ts.data.to_numpy(), **kwargs)

        return ptf

    #==================================================
    def plot_spectra(self, fig_num=2):
        """
        plot the spectra of time series
        """
        self.ts_obj.plot_spectra(fig_num=fig_num)


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
        self.df_list = (4096, 256)
        self.df_time_list = ('00:10:00','07:50:00')
        self.master_schedule = self.make_schedule(self.df_list,
                                                  self.df_time_list,
                                                  repeat=16)
        self._resync_pause = 20

    #==================================================
    def read_schedule(self, fn):
        """
        read zen schedule file
        """

        sfid = open(fn, 'r')
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
            time_list = [{'dt':self.initial_dt, 'df':t1_dict['df']}]

            kk = np.where(np.array(df_list) == t1_dict['df'])[0][0]-ndf+1
            df_list = np.append(df_list[kk:], df_list[:kk])
            df_length_list = np.append(df_length_list[kk:], df_length_list[:kk])
            time_list.append(dict([('dt', t1_dict['dt']), ('df', df_list[0])]))
            ii = 1
        else:
            time_list = [{'dt':self.initial_dt, 'df':df_list[0]}]
            ii = 0

        for rr in range(1, repeat+1):
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

        self.master_schedule = self.make_schedule(self.df_list,
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
            sfid = open(os.path.normpath(os.path.join('c:\\MT', save_name)),
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
                                    '0', '0', '0', 'y', 'n', 'n', 'n'])
                sfid.write('scheduleaction '.upper()+sa_line[:-1]+'\n')
            meta_line = ''.join(['{0},{1}|'.format(key,self.meta_dict[key])
                                 for key in self.meta_keys])
            sfid.write('METADATA '+meta_line+'\n')
            for lkey in list(self.light_dict.keys()):
                sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
            sfid.close()
            # print('Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
            #                                       self.ch_cmp_dict[dname[-1]]))

            for dd in list(drive_names.keys()):
                dname = drive_names[dd]
                sfid = open(os.path.normpath(os.path.join(dd+':\\', save_name)),
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
                                        '0', '0', '0', 'y', 'n', 'n', 'n'])
                    sfid.write('scheduleaction '.upper()+sa_line[:-1]+'\n')

                self.meta_dict['Ch.Cmp'] = self.ch_cmp_dict[dname[-1]]
                self.meta_dict['Ch.Number'] = dname[-1]
                meta_line = ''.join(['{0},{1}|'.format(key,self.meta_dict[key])
                                     for key in self.meta_keys])
                sfid.write('METADATA '+meta_line+'\n')
                for lkey in list(self.light_dict.keys()):
                    sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
                sfid.close()

                print('Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
                                                   self.ch_cmp_dict[dname[-1]]))
            return
        else:
            save_name = savename

        for dd in list(drive_names.keys()):
            dname = drive_names[dd]
            sfid = open(os.path.normpath(os.path.join(dd+':\\', save_name)),
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
            for lkey in list(self.light_dict.keys()):
                sfid.write('{0} {1}\n'.format(lkey, self.light_dict[lkey]))
            sfid.close()

            print('Wrote {0}:\{1} to {2} as {3}'.format(dd, save_name, dname,
                                                   self.ch_cmp_dict[dname[-1]]))

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
        self.master_schedule = self.make_schedule(self.df_list,
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

            # subtract 10 seconds for transition between schedule items.
            t_diff = t1-t0-self._resync_pause
            zacq_list.append('$schline{0:.0f} = {1:.0f},{2:.0f},{3:.0f}\n'.format(
                                ii+1,
                                t_diff,
                                int(self.sr_dict[str(ss['df'])]),
                                1))

        fn = os.path.join(save_path, schedule_fn)
        fid = open(fn, 'w')
        fid.writelines(zacq_list[0:16])
        fid.close()

        print('Wrote schedule file to {0}'.format(fn))
        print('+--------------------------------------+')
        print('|   SET ZEN START TIME TO: {0}    |'.format(zen_start))
        print('+--------------------------------------+')

    def _convert_time_to_seconds(self, time_string):
        """
        convert a time string given as hh:mm:ss into seconds
        """
        t_list = [float(tt) for tt in time_string.split(':')]
        t_seconds = t_list[0] * 3600 +t_list[1] * 60 + t_list[2]

        return t_seconds

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
    for letter in string.ascii_uppercase:
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
        print('No external drives detected, check the connections.')
        return None
    return drive_dict

def split_station(station):
    """
    split station name into name and number
    """

    for ii, ss in enumerate(station):
        try:
            int(ss)
            find = ii
            break
        except ValueError:
            continue

    name = station[0:find]
    number = station[find:]

    return (name, number)

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
    s_name, s_int = split_station(station)
    drive_names = get_drive_names()
    if drive_names is None:
        raise IOError('No drives to copy from.')
    save_path = os.path.join(save_path, station)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    log_fid = open(os.path.join(save_path, 'copy_from_sd.log'), 'w')
    
    # make a datetime object from copy date
    if copy_date is not None:
        c_date = dateutil.parser.parse(copy_date)
        
    st_test = time.ctime()
    fn_list = []
    for key in list(drive_names.keys()):
        dr = r"{0}:\\".format(key)
        print('='*25+drive_names[key]+'='*25)
        log_fid.write('='*25+drive_names[key]+'='*25+'\n')

        for fn in os.listdir(dr):
            if 'cal' in fn.lower():
                continue
            full_path_fn = os.path.normpath(os.path.join(dr, fn))
            if fn[-4:] == '.cfg':
                shutil.copy(full_path_fn, os.path.join(save_path, fn))
            # test for copy date
            if copy_date is not None:
                file_date = datetime.datetime.fromtimestamp(
                                    os.path.getmtime(full_path_fn))
                if copy_type == 'after':
                    if file_date < c_date:
                        continue
                elif copy_type == 'before':
                    if file_date > c_date:
                        continue
                elif copy_type == 'on':
                    if file_date.date() != c_date.date():
                        continue
                        
            try:
                file_size = os.stat(full_path_fn)[6]
                if file_size >= 1600 and fn.find('.cfg') == -1:
                    zt = Zen3D(fn=full_path_fn)
                    zt.read_all_info()

                    if zt.metadata.station.find(s_int) >= 0:
                        channel = zt.metadata.ch_cmp.upper()
                        st = zt.schedule.Time.replace(':', '')
                        sd = zt.schedule.Date.replace('-', '')
                        sv_fn = '{0}_{1}_{2}_{3}_{4}.Z3D'.format(station,
                                                                 sd,
                                                                 st,
                                                                 int(zt.df),
                                                                 channel)

                        full_path_sv = os.path.join(save_path, sv_fn)
                        fn_list.append(full_path_sv)

                        shutil.copy(full_path_fn, full_path_sv)

                        print('copied {0} to {1}\n'.format(full_path_fn, 
                                                           full_path_sv))
                                                         
                        #log_fid.writelines(zt.log_lines)

                        log_fid.write('copied {0} to \n'.format(full_path_fn)+\
                                      '       {0}\n'.format(full_path_sv))

                else:
                    log_fid.write('+++ Skipped {0} because file to small {1}'.format(full_path_fn, 
                                                                                      file_size))
            except WindowsError:
                print('Faulty file at {0}'.format(full_path_fn))
                log_fid.write('---Faulty file at {0}\n\n'.format(full_path_fn))
    log_fid.close()

    et_test = time.ctime()
    print('Started copying at: {0}'.format(st_test))
    print('Ended copying at: {0}'.format(et_test))

    return fn_list

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
        log_fid = open(os.path.join(delete_folder, 'Log_file.log'), 'w')

    if delete_date is not None:
        delete_date = int(delete_date.replace('-', ''))

    delete_fn_list = []
    for key, value in drive_names.items():
        dr = r"{0}:\\".format(key)
        log_lines.append('='*25+value+'='*25+'\n')
        for fn in os.listdir(dr):
            if fn[-4:].lower() == '.Z3D'.lower():
                full_path_fn = os.path.normpath(os.path.join(dr, fn))
                zt = Zen3D(full_path_fn)
                zt.read_all_info()
                zt_date = int(zt.schedule.Date.replace('-', ''))
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
                    #zt_date = int(zt.schedule_date.replace('-',''))

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
        log_fid = open(os.path.join(delete_folder, 'Delete_log.log'), 'w')
        log_fid.writelines(log_lines)
        log_fid.close()
    if verbose:
        for lline in log_lines:
            print(lline)

    return delete_fn_list

#==============================================================================
#   Make mtpy_mt files
#==============================================================================
def make_mtpy_mt_files(fn_list, station_name='mb', fmt='%.8e',
                       notch_dict=None):
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
                      dtype=[('station', '|S6'), ('len', np.int), ('df', np.int),
                             ('start_dt', '|S22'), ('comp', '|S2'),
                             ('fn', '|S100')])
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
        else:
            #write mtpy mt file
            zd.write_ascii_mt_file(save_station=station_name,
                                   fmt=fmt,
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

