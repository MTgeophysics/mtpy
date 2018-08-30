# -*- coding: utf-8 -*-
"""
USGS Archive
==============
    
    * Collect z3d files into logical scheduled blocks
    * Merge Z3D files into USGS ascii format
    * Collect metadata information
    * make .csv, .xml, .shp files.

Created on Tue Aug 29 16:38:28 2017

@author: jpeacock
"""
#==============================================================================

import os
import time
import datetime
from collections import Counter

import zipfile
import gzip
import urllib2 as url
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

import numpy as np
import scipy.signal as sps
import pandas as pd

import mtpy.usgs.zen as zen
import mtpy.usgs.zonge as zonge
import mtpy.utils.gis_tools as gis_tools
import mtpy.utils.configfile as mtcfg

# for writing shape file
import geopandas as gpd
from shapely.geometry import Point


# =============================================================================
# 
# =============================================================================
def get_nm_elev(lat, lon):
    """
    Get national map elevation for a given lat and lon.
    
    Queries the national map website for the elevation value.
    
    :param lat: latitude in decimal degrees
    :type lat: float
    
    :param lon: longitude in decimal degrees
    :type lon: float
    
    :return: elevation (meters)
    :rtype: float
    
    :Example: ::
        
        >>> import mtpy.usgs.usgs_archive as archive
        >>> archive.get_nm_elev(35.467, -115.3355)
        >>> 809.12
        
    .. note:: Needs an internet connection to work.
        
    """
    nm_url = r"https://nationalmap.gov/epqs/pqs.php?x={0:.5f}&y={1:.5f}&units=Meters&output=xml"

    # call the url and get the response
    response = url.urlopen(nm_url.format(lon, lat))
    
    # read the xml response and convert to a float
    info = ET.ElementTree(ET.fromstring(response.read()))
    info = info.getroot()
    for elev in info.iter('Elevation'):
        nm_elev = float(elev.text) 
    return nm_elev

# =============================================================================
# Collect Z3d files
# =============================================================================
class Z3DCollection(object):
    """
    Collects .z3d files into useful arrays and lists
    
    ================= ============================= ===========================
    Attribute         Description                   Default
    ================= ============================= ===========================
    chn_order         list of the order of channels [hx, ex, hy, ey, hz]
    meta_notes        extraction of notes from      None
                      the .z3d files 
    leap_seconds      number of leap seconds for    16 [2016]
                      a given year
    ================= ============================= ===========================
    
    ===================== =====================================================
    Methods               Description
    ===================== =====================================================
    get_time_blocks       Get a list of files for each schedule action
    check_sampling_rate   Check the sampling rate a given time block
    check_time_series     Get information for a given time block
    merge_ts              Merge a given schedule block making sure that they
                          line up in time.  
    get_chn_order         Get the appropriate channels, in case some are 
                          missing 
    ===================== =====================================================
    
    :Example: ::
        
        >>> import mtpy.usgs.usgs_archive as archive
        >>> z3d_path = r"/Data/Station_00"
        >>> zc = archive.Z3DCollection()
        >>> fn_list = zc.get_time_blocks(z3d_path)
        
    """   
    
    def __init__(self):

        self.chn_order = ['hx','ex','hy','ey','hz']
        self.meta_notes = None
        self.leap_seconds = 16
        self.verbose = True
        
    def get_time_blocks(self, z3d_dir):
        """
        Organize z3d files into blocks based on start time and sampling rate
        in the file name.
        
        .. note:: This assumes the z3d file is named
                  * station_date_time_samplingrate_chn.z3d
                  
        :param z3d_dir: full path to z3d files
        :type z3d_dir: string
        
        :returns: nested list of files for each time block, sorted by time
        
        :Example: ::
            
            >>> import mtpy.usgs.usgs_archive as archive
            >>> zc = archive.Z3DCollection()
            >>> fn_list = zc.get_time_blocks(r"/home/mt_data/station_01")
         
        """
        
        fn_list = os.listdir(z3d_dir)
        merge_list = np.array([[fn]+\
                              os.path.basename(fn)[:-4].split('_')
                              for fn in fn_list if fn.lower().endswith('.z3d')])
                                  
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
        
        return_list = []
        for tt in sorted(time_list):
            block_list = []
            log_lines.append('+'*72+'\n')
            log_lines.append('Files Being Merged: \n')
            time_fn_list = merge_list[np.where(merge_list == tt)[0], 0].tolist()
            for cfn in time_fn_list:
                log_lines.append(' '*4+cfn+'\n')
                block_list.append(os.path.join(z3d_dir, cfn))

            return_list.append(block_list)
            log_lines.append('\n---> Merged Time Series Lengths and Start Time \n')
            log_lines.append('\n')
            
        with open(os.path.join(z3d_dir, 'z3d_merge.log'), 'w') as fid:
            fid.writelines(log_lines)
        
        return return_list
   
    #==================================================
    def check_sampling_rate(self, time_array):
        """
        Check to make sure the sampling rate is the same for all channels
        
        :param time_array: structured array with key df for sampling rate
                           made by check_time_series 
        :type time_array: np.ndarray
        
        :returns: sampling rate
        :rtype: float
        """
        
        nz = len(time_array)
        
        df_array = time_array['df']
            
        tf_array = np.zeros((nz, nz))
        
        for jj in range(nz):
            tf_array[jj] = np.in1d(df_array, [df_array[jj]])
        
        false_test = np.where(tf_array == False)
        
        if len(false_test[0]) != 0:
            raise IOError('Sampling rates are not the same for all channels '+\
                          'Check file(s)'+time_array[false_test[0]]['fn'])
            
        return df_array.mean()
        
    #==================================================
    def check_time_series(self, fn_list):
        """
        Check to make sure timeseries line up with eachother.
        
        :param fn_list: list of files to merge
        :type fn_list: list
        
        :return: time_array
        :rtype: np.ndarray including all important information for merging
        """
        # get number of channels
        n_fn = len(fn_list)
        
        # make an empty array to put things into
        t_arr = np.zeros(n_fn, 
                         dtype=[('comp', 'S3'),
                                ('start', np.int64),
                                ('stop', np.int64),
                                ('fn', 'S140'),
                                ('df', np.float32),
                                ('lat', np.float32),
                                ('lon', np.float32),
                                ('elev', np.float32), 
                                ('ch_azm', np.float32),
                                ('ch_length', np.float32),
                                ('ch_num', np.int32),
                                ('ch_box', 'S6'),
                                ('n_samples', np.int32),
                                ('t_diff', np.int32),
                                ('std', np.float32)])
    
        t_arr['ch_num'] = np.arange(1, n_fn+1)
            
        print('-'*50)
        for ii, fn in enumerate(fn_list):
            z3d_obj = zen.Zen3D(fn)
            z3d_obj._leap_seconds = self.leap_seconds
            try:
                z3d_obj.read_z3d()
            except zen.ZenGPSError:
                print('xxxxxx BAD FILE: Skipping {0} xxxx'.format(fn))
                continue
            
            # convert the time index into an integer
            dt_index = z3d_obj.ts_obj.ts.data.index.astype(np.int64)/10.**9
            
            # extract useful data
            t_arr[ii]['comp'] = z3d_obj.metadata.ch_cmp.lower()
            t_arr[ii]['start'] = dt_index[0]
            t_arr[ii]['stop'] = dt_index[-1]
            t_arr[ii]['fn'] = fn
            t_arr[ii]['df'] = z3d_obj.df
            t_arr[ii]['lat'] = z3d_obj.header.lat
            t_arr[ii]['lon'] = z3d_obj.header.long
            t_arr[ii]['elev'] = z3d_obj.header.alt
            t_arr[ii]['ch_azm'] = z3d_obj.metadata.ch_azimuth
            if 'e' in t_arr[ii]['comp']:
                t_arr[ii]['ch_length'] = z3d_obj.metadata.ch_length
            if 'h' in t_arr[ii]['comp']:
                t_arr[ii]['ch_num'] = int(float(z3d_obj.metadata.ch_number))
            t_arr[ii]['ch_box'] = int(z3d_obj.header.box_number)
            t_arr[ii]['n_samples'] = z3d_obj.ts_obj.ts.shape[0]
            t_arr[ii]['t_diff'] = int((dt_index[-1]-dt_index[0])*z3d_obj.df)-\
                                      z3d_obj.ts_obj.ts.shape[0]
            t_arr[ii]['std'] = z3d_obj.ts_obj.ts.std()
            try:
                self.meta_notes = z3d_obj.metadata.notes.replace('\r', ' ').replace('\x00', '').rstrip()
            except AttributeError:
                pass

        # cut the array to only those channels with data
        t_arr = t_arr[np.nonzero(t_arr['start'])]
        
        return t_arr
    
    def merge_ts(self, fn_list, decimate=1):
        """
        Merge z3d's based on a mutual start and stop time
        
        :param fn_list: list of Z3D files to merge together (full paths)
        :type fn_list: list
        
        :return: pandas database of merged time series
        :rtype: pandas database
        
        :return: time_array that includes all important information
        :rtype: np.ndarray
        
        """
        meta_arr = self.check_time_series(fn_list)
        df = self.check_sampling_rate(meta_arr)
        
        # get the start and stop times that correlates with all time series
        start = meta_arr['start'].max()
        stop = meta_arr['stop'].min()
        
        # figure out the max length of the array, getting the time difference into
        # seconds and then multiplying by the sampling rate
        max_ts_len = int((stop-start)*df)
        ts_len = min([meta_arr['n_samples'].max(), max_ts_len])
        
        if decimate > 1:
            ts_len /= decimate
        
        print(ts_len, meta_arr.size)
        ts_db = pd.DataFrame(np.zeros((ts_len, meta_arr.size)),
                             columns=list(meta_arr['comp']),
                             dtype=np.float32)
        
        for ii, m_arr in enumerate(meta_arr):
            z3d_obj = zen.Zen3D(m_arr['fn'])
            z3d_obj._leap_seconds = self.leap_seconds
            z3d_obj.read_z3d()
            
            dt_index = z3d_obj.ts_obj.ts.data.index.astype(np.int64)/10**9
            index_0 = np.where(dt_index == start)[0][0]
            #index_1 = np.where(dt_index == stop)[0][0]
            index_1 = min([ts_len-index_0, z3d_obj.ts_obj.ts.shape[0]-index_0])
            t_diff = ts_len-(index_1-index_0)
            meta_arr[ii]['t_diff'] = t_diff

            if t_diff != 0:
                if self.verbose:
                    print '{0} off by {1} points --> {2} sec'.format(z3d_obj.ts_obj.fn,
                                                                     t_diff,
                                                                     t_diff/z3d_obj.ts_obj.sampling_rate)
            if decimate > 1:
                 ts_db[:, ii] = sps.resample(z3d_obj.ts_obj.ts.data[index_0:index_1],
                                              ts_len, 
                                              window='hanning')

            else:
                ts_db[m_arr['comp']][0:ts_len-t_diff] = z3d_obj.ts_obj.ts.data[index_0:index_1]

        # reorder the columns        
        ts_db = ts_db[self.get_chn_order]
        
        # return the pandas database and the metadata array 
        return ts_db, meta_arr 
    
    def get_chn_order(self, chn_list):
        """
        Get the order of the array channels according to the components.
        
        .. note:: If you want to change the channel order, change the 
                  parameter Z3DCollection.chn_order
        
        :param chn_list: list of channels in data
        :type chn_list: list
        
        :return: channel order list 
        """
        
        if len(chn_list) == 5:
            return self.chn_order
        else:
            chn_order = []
            for chn_00 in self.chn_order:
                for chn_01 in chn_list:
                    if chn_00.lower() == chn_01.lower():
                        chn_order.append(chn_00.lower())
                        continue
            
            return chn_order

#==============================================================================
# Need a dummy utc time zone for the date time format
#==============================================================================
class UTC(datetime.tzinfo):
    def utcoffset(self, df):
        return datetime.timedelta(hours=0)
    def dst(self, df):
        return datetime.timedelta(0)
    def tzname(self, df):
        return "UTC"

# =============================================================================
#  Metadata for usgs ascii file
# =============================================================================
class Metadata(object):
    """
    Container for all the important metadata in a USGS ascii file.
    
    ========================= =================================================
    Attributes                Description
    ========================= =================================================
    SurveyID                  Survey name
    SiteID                    Site name
    RunID                     Run number
    SiteLatitude              Site latitude in decimal degrees WGS84
    SiteLongitude             Site longitude in decimal degrees WGS84
    SiteElevation             Site elevation according to national map meters
    AcqStartTime              Start time of station YYYY-MM-DDThh:mm:ss UTC
    AcqStopTime               Stop time of station YYYY-MM-DDThh:mm:ss UTC
    AcqSmpFreq                Sampling rate samples/second
    AcqNumSmp                 Number of samples 
    Nchan                     Number of channels
    CoordinateSystem          [ Geographic North | Geomagnetic North ]
    ChnSettings               Channel settings, see below
    MissingDataFlag           Missing data value
    ========================= =================================================
    
    *ChnSettings*
    ========================= =================================================
    Keys                      Description
    ========================= =================================================
    ChnNum                    SiteID+channel number
    ChnID                     Component [ ex | ey | hx | hy | hz ]
    InstrumentID              Data logger + sensor number
    Azimuth                   Setup angle of componet in degrees relative to 
                              CoordinateSystem
    Dipole_Length             Dipole length in meters 
    ========================= =================================================                

    
    """
    def __init__(self, fn=None, **kwargs):
        self.fn = fn
        self.SurveyID = None
        self.RunID = None
        self._station = None
        self._latitude = None
        self._longitude = None
        self._elevation = None
        self._start_time = None
        self._stop_time = None
        self._sampling_rate = None
        self._n_samples = None
        self.channel_dict = None
        self.MissingDataFlag = np.NaN
        self._chn_num = None
        self.CoordinateSystem = None
        self._time_fmt = '%Y-%m-%dT%H:%M:%S %Z'
        self._metadata_len = 30
        self.declination = 0.0
        
        self._key_list = ['SurveyID',
                          'SiteID',
                          'RunID',
                          'SiteLatitude',
                          'SiteLongitude',
                          'SiteElevation',
                          'AcqStartTime',
                          'AcqStopTime',
                          'AcqSmpFreq',
                          'AcqNumSmp',
                          'Nchan',
                          'CoordinateSystem',
                          'ChnSettings',
                          'MissingDataFlag',
                          'DataSet']
        
        self._chn_settings = ['ChnNum',
                              'ChnID',
                              'InstrumentID',
                              'Azimuth',
                              'Dipole_Length']
        self._chn_fmt = {'ChnNum':'<8',
                         'ChnID':'<6',
                         'InstrumentID':'<12',
                         'Azimuth':'>7.1f',
                         'Dipole_Length':'>14.1f'}

        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    @property
    def SiteID(self):
        return self._station
    @SiteID.setter
    def SiteID(self, station):
        self._station = station
        
    @property
    def SiteLatitude(self):
        return self._latitude
        #return gis_tools.convert_position_float2str(self._latitude)
    
    @SiteLatitude.setter
    def SiteLatitude(self, lat):
        self._latitude = gis_tools.assert_lat_value(lat)
        
    @property
    def SiteLongitude(self):
        return self._longitude
        #return gis_tools.convert_position_float2str(self._longitude)
    
    @SiteLongitude.setter
    def SiteLongitude(self, lon):
        self._longitude = gis_tools.assert_lon_value(lon)
        
    @property
    def SiteElevation(self):
        """
        get elevation from national map
        """  
        # the url for national map elevation query
        nm_url = r"https://nationalmap.gov/epqs/pqs.php?x={0:.5f}&y={1:.5f}&units=Meters&output=xml"

        # call the url and get the response
        try:
            response = url.urlopen(nm_url.format(self._longitude, self._latitude))
        except url.HTTPError:
            print nm_url.format(self._longitude, self._latitude)
            return -666
        
        # read the xml response and convert to a float
        info = ET.ElementTree(ET.fromstring(response.read()))
        info = info.getroot()
        for elev in info.iter('Elevation'):
            nm_elev = float(elev.text) 
        return nm_elev 
        
    @property
    def AcqStartTime(self):
        return self._start_time.strftime(self._time_fmt)
    
    @AcqStartTime.setter
    def AcqStartTime(self, time_string):
        if type(time_string) in [int, np.int64]:
            dt = datetime.datetime.utcfromtimestamp(time_string)
        elif type(time_string) in [str]:
            dt = datetime.datetime.strptime(time_string, self._time_fmt)
        self._start_time = datetime.datetime(dt.year, dt.month, dt.day,
                                             dt.hour, dt.minute, dt.second,
                                             dt.microsecond, tzinfo=UTC())
        
    @property
    def AcqStopTime(self):
        return self._stop_time.strftime(self._time_fmt)
    
    @AcqStopTime.setter
    def AcqStopTime(self, time_string):
        if type(time_string) in [int, np.int64]:
            dt = datetime.datetime.utcfromtimestamp(time_string)
        elif type(time_string) in [str]:
            dt = datetime.datetime.strptime(time_string, self._time_fmt)
        self._stop_time = datetime.datetime(dt.year, dt.month, dt.day,
                                            dt.hour, dt.minute, dt.second,
                                            dt.microsecond, tzinfo=UTC())
    
    @property
    def Nchan(self):
        return self._chn_num
    
    @Nchan.setter
    def Nchan(self, n_channel):
        try:
            self._chn_num = int(n_channel)
        except ValueError:
            print("{0} is not a number, setting Nchan to 0".format(n_channel))
            
    @property
    def AcqSmpFreq(self):
        return self._sampling_rate
    @AcqSmpFreq.setter
    def AcqSmpFreq(self, df):
        self._sampling_rate = float(df)

    @property
    def AcqNumSmp(self):
        return self._n_samples

    @AcqNumSmp.setter
    def AcqNumSmp(self, n_samples):
        self._n_samples = int(n_samples)  

    def read_metadata(self, fn=None, meta_lines=None):
        """
        Read in a meta from the raw string or file.  Populate all metadata
        as attributes.
        
        :param fn: full path to USGS ascii file
        :type fn: string
        
        :param meta_lines: lines of metadata to read
        :type meta_lines: list
        """
        chn_find = False
        comp = 0
        self.channel_dict = {}
        if fn is not None:
            self.fn = fn
        if self.fn is not None:
            with open(self.fn, 'r') as fid:
                meta_lines = [fid.readline() for ii in range(self._metadata_len)]
        for ii, line in enumerate(meta_lines):
            if line.find(':') > 0:
                key, value = line.strip().split(':', 1)
                value = value.strip()
                if len(value) < 1 and key == 'DataSet':
                    chn_find = False
                    # return the line that the data starts on that way can
                    # read in as a numpy object or pandas
                    return ii+1
                elif len(value) < 1:
                    chn_find = True
                setattr(self, key, value)
            elif 'coordinate' in line:
                self.CoordinateSystem = ' '.join(line.strip().split()[-2:])
            else:
                if chn_find is True:
                    if 'chnnum' in line.lower():
                        ch_key = line.strip().split()
                    else:
                        line_list = line.strip().split()
                        if len(line_list) == 5:
                            comp += 1
                            self.channel_dict[comp] = {}
                            for key, value in zip(ch_key, line_list):
                                if key.lower() in ['azimuth', 'dipole_length']:
                                    value = float(value)
                                self.channel_dict[comp][key] = value
                        else:
                            print('Not sure what line this is')
                            
    def write_metadata(self, chn_list=['Ex', 'Ey', 'Hx', 'Hy', 'Hz']):
        """
        Write out metadata in the format of USGS ascii.
        
        :return: list of metadate lines.
        
        .. note:: meant to use '\n'.join(lines) to write out in a file.
        """
        
        lines = []
        for key in self._key_list:
            if key in ['ChnSettings']:
                lines.append('{0}:'.format(key))
                lines.append(' '.join(self._chn_settings))
                for chn_key in chn_list:
                    chn_line = []
                    try:
                        for comp_key in self._chn_settings:
                            chn_line.append('{0:{1}}'.format(self.channel_dict[chn_key][comp_key],
                                            self._chn_fmt[comp_key]))
                        lines.append(''.join(chn_line))
                    except KeyError:
                        pass
            elif key in ['DataSet']:
                lines.append('{0}:'.format(key))
                return lines
            else:
                if key in ['SiteLatitude', 'SiteLongitude']:
                    lines.append('{0}: {1:.5f}'.format(key, getattr(self, key)))
                else:
                    lines.append('{0}: {1}'.format(key, getattr(self, key)))
        
        return lines


# =============================================================================
# Class for the asc file
# =============================================================================
class USGSasc(Metadata):
    """
    Read and write USGS ascii formatted time series. 
    
    =================== =======================================================
    Attributes          Description
    =================== =======================================================
    ts                  Pandas dataframe holding the time series data
    fn                  Full path to .asc file
    station_dir         Full path to station directory
    meta_notes          Notes of how the station was collected
    =================== =======================================================
    
    ============================== ============================================
    Methods                        Description
    ============================== ============================================
    get_z3d_db                     Get Pandas dataframe for schedule block
    locate_mtft24_cfg_fn           Look for a mtft24.cfg file in station_dir
    get_metadata_from_mtft24       Get metadata from mtft24.cfg file
    get_metadata_from_survey_csv   Get metadata from survey csv file
    fill_metadata                  Fill Metadata container from a meta_array
    read_asc_file                  Read in USGS ascii file
    convert_electrics              Convert electric channels to mV/km
    write_asc_file                 Write an USGS ascii file
    write_station_info_metadata    Write metadata to a .cfg file
    ============================== ============================================
    
    :Example: ::
        
        >>> zc = Z3DCollection()
        >>> fn_list = zc.get_time_blocks(z3d_path)
        >>> zm = USGSasc()
        >>> zm.SurveyID = 'iMUSH'
        >>> zm.get_z3d_db(fn_list[0])
        >>> zm.read_mtft24_cfg()
        >>> zm.CoordinateSystem = 'Geomagnetic North'
        >>> zm.SurveyID = 'MT'
        >>> zm.write_asc_file(str_fmt='%15.7e')
        >>> zm.write_station_info_metadata()
    """
    
    def __init__(self, **kwargs):
        Metadata.__init__(self)
        self.ts = None
        self.fn = None
        self.station_dir = os.getcwd()
        self.meta_notes = None
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def get_z3d_db(self, fn_list):
        """
        Merge time series from Z3D files into a pandas database
        
        :param fn_list: list of Z3D files to merge (full paths)
        :type fn_list: list
        
        Fills ts attribute as pandas database.
        """
        zc_obj = Z3DCollection()
        self.ts, meta_arr = zc_obj.merge_ts(fn_list)
        self.fill_metadata(meta_arr)
        self.meta_notes = zc_obj.meta_notes
        
    def locate_mtft24_cfg_fn(self):
        """
        Try to automatically detect mtft24 file
        
        :return: path to mtft24.cfg file
        """
        for fn in os.listdir(self.station_dir):
            if 'mtft24' in fn and fn.endswith('cfg'):
                return os.path.join(self.station_dir, fn)
            
        return None
            
    def get_metadata_from_mtft24_cfg(self, mtft24_cfg_fn=None):
        """
        Read in a MTFT24 configuration file and fill in meta data
        
        :param mtft24_cfg_fn: full path to mtft24.cfg
        :type: string
        """
        if mtft24_cfg_fn is None:
            mtft24_cfg_fn = self.locate_mtft24_cfg_fn()
            
        zm = zonge.ZongeMTFT()
        try:
            zm.read_cfg(mtft24_cfg_fn)
        except TypeError:
            print('*** No MTFT24 file for {0} ***'.format(self.SiteID))
            return False
        
        # need to update channel dict
        # figure out channel order first
        chn_order = dict([(cc.lower(), ii) for ii, cc in enumerate(zm.Chn_Cmp)])
        for chn in self.channel_dict.keys():
            index = chn_order[chn.lower()]
            # get coil number
            if 'h' in chn.lower():
                stem = self.channel_dict[chn]['InstrumentID'].split('-')[0]
                try:
                    inst = self.channel_dict[chn]['InstrumentID'].split('-')[1]
                    if inst != zm.Chn_ID[index]:
                        self.channel_dict[chn]['InstrumentID'] = '{0}-{1}'.format(stem, 
                                                                                 zm.Chn_ID[index])
                except IndexError:
                    self.channel_dict[chn]['InstrumentID'] = '{0}-{1}'.format(stem, 
                                                                             zm.Chn_ID[index])
            
            # get azimuth direction
            if 'geographic' in self.CoordinateSystem.lower() and chn.lower() != 'hz':
                azm = float(zm.Chn_Azimuth[index])+self.declination
            else:
                azm = float(zm.Chn_Azimuth[index])
            if np.isnan(self.channel_dict[chn]['Azimuth']):
                self.channel_dict[chn]['Azimuth'] = azm
            elif int(np.nan_to_num(self.channel_dict[chn]['Azimuth'])) != int(azm):
                self.channel_dict[chn]['Azimuth'] = azm
                
            # get dipole length
            if float(np.nan_to_num(self.channel_dict[chn]['Dipole_Length'])) != float(zm.Chn_Length[index]):
                if 'e' in chn.lower(): 
                    self.channel_dict[chn]['Dipole_Length'] = float(zm.Chn_Length[index])
        return True
        
    def get_metadata_from_survey_csv(self, survey_fn):
        """
        get station information from a survey .csv file
        
        :param survey_fn: full path to survey summary .csv file
        :type survey_fn: string
        
        """
        
        s_cfg = USGScfg()
        try:
            station_db = s_cfg.get_station_info_from_csv(survey_fn, self.SiteID)
        except ValueError:
            print('Could not find information for {0}'.format(self.SiteID))
            return False
        
        # fill metadata
        for chn in self.channel_dict.keys():
            if 'h' in chn.lower():
                #stem = self.channel_dict[chn]['InstrumentID'].split('-', 1)[0]
                stem = station_db.zen_num
                h_attr = '{0}_{1}'.format(chn.lower(), 'id')
                h_id = getattr(station_db, h_attr)
                self.channel_dict[chn]['InstrumentID'] = '{0}-{1}'.format(stem, 
                                                                          h_id)
            elif 'e' in chn.lower():
                self.channel_dict[chn]['InstrumentID'] = station_db.zen_num
                e_attr = '{0}_{1}'.format(chn.lower(), 'len')
                e_len = getattr(station_db, e_attr)
                self.channel_dict[chn]['Dipole_Length'] = e_len
        
            azm_attr = '{0}_{1}'.format(chn.lower(), 'azm')
            azm_value = getattr(station_db, azm_attr)
            if 'geographic' in self.CoordinateSystem:
                azm_value += self.declination
            self.channel_dict[chn]['Azimuth'] = azm_value
            self.channel_dict[chn]['ChnNum'] = '{0}{1}'.format(self.channel_dict[chn]['ChnNum'][:-1],
                                                               getattr(station_db, 
                                                                       '{0}_num'.format(chn.lower())))
            
            
        # get location
        self.SiteLatitude = float(station_db.lat)
        self.SiteLongitude = float(station_db.lon)
            
        return True
        
    def fill_metadata(self, meta_arr):
        """
        Fill in metadata from time array made by 
        Z3DCollection.check_time_series.
        
        :param meta_arr: structured array of metadata for the Z3D files to be
                         combined.
        :type meta_arr: np.ndarray
        """
        try:
            self.AcqNumSmp = self.ts.shape[0]
        except AttributeError:
            pass
        self.AcqSmpFreq = meta_arr['df'].mean()
        self.AcqStartTime = meta_arr['start'].max()
        self.AcqStopTime = meta_arr['stop'].min()
        try:
            self.Nchan = self.ts.shape[1]
        except AttributeError:
            self.Nchan = meta_arr.shape[0]
        self.RunID = 1
        self.SiteLatitude = np.median(meta_arr['lat'])
        self.SiteLongitude = np.median(meta_arr['lon'])
        self.SiteID = os.path.basename(meta_arr['fn'][0]).split('_')[0]
        self.station_dir = os.path.dirname(meta_arr['fn'][0])
        
        # if geographic coordinates add in declination
        if 'geographic' in self.CoordinateSystem.lower():
            meta_arr['ch_azm'][np.where(meta_arr['comp'] != 'hz')] += self.declination

        # fill channel dictionary with appropriate values
        self.channel_dict = dict([(comp.capitalize(),
                                   {'ChnNum':'{0}{1}'.format(self.SiteID, ii+1),
                                    'ChnID':meta_arr['comp'][ii].capitalize(),
                                    'InstrumentID':meta_arr['ch_box'][ii],
                                    'Azimuth':meta_arr['ch_azm'][ii],
                                    'Dipole_Length':meta_arr['ch_length'][ii],
                                    'n_samples':meta_arr['n_samples'][ii],
                                    'n_diff':meta_arr['t_diff'][ii],
                                    'std':meta_arr['std'][ii],
                                    'start':meta_arr['start'][ii]})
                                   for ii, comp in enumerate(meta_arr['comp'])])
        for ii, comp in enumerate(meta_arr['comp']):
            if 'h' in comp.lower():
                self.channel_dict[comp.capitalize()]['InstrumentID'] += '-{0}'.format(meta_arr['ch_num'])

    def read_asc_file(self, fn=None):
        """
        Read in a USGS ascii file and fill attributes accordingly.
        
        :param fn: full path to .asc file to be read in
        :type fn: string
        """
        if fn is not None:
            self.fn = fn
        st = datetime.datetime.now()    
        data_line = self.read_metadata()
        self.ts = pd.read_csv(self.fn,
                              delim_whitespace=True,
                              skiprows=data_line,
                              dtype=np.float32)
        et = datetime.datetime.now()
        read_time = et-st
        print('Reading took {0}'.format(read_time.total_seconds()))
        
    def convert_electrics(self):
        """
        Convert electric fields into mV/km
        """
        
        try:
            self.ts.ex /= (self.channel_dict['Ex']['Dipole_Length']/1000.)
        except AttributeError:
            print('No EX')
            
                
        try:
            self.ts.ey /= (self.channel_dict['Ey']['Dipole_Length']/1000.)
        except AttributeError:
            print('No EY')
            
    def _make_file_name(self, save_path=None, compression=True, 
                        compress_type='zip'):
        """
        get the file name to save to
        
        :param save_path: full path to directory to save file to
        :type save_path: string
        
        :param compression: compress file
        :type compression: [ True | False ]
        
        :return: save_fn
        :rtype: string
        
        """
        # make the file name to save to
        if save_path is not None:
            save_fn = os.path.join(save_path, 
                                   '{0}_{1}T{2}_{3:.0f}.asc'.format(self.SiteID,
                                    self._start_time.strftime('%Y-%m-%d'),
                                    self._start_time.strftime('%H%M%S'),
                                    self.AcqSmpFreq))
        else:
            save_fn = os.path.join(self.station_dir, 
                                   '{0}_{1}T{2}_{3:.0f}.asc'.format(self.SiteID,
                                    self._start_time.strftime('%Y-%m-%d'),
                                    self._start_time.strftime('%H%M%S'),
                                    self.AcqSmpFreq))
            
        if compression:
            if compress_type == 'zip':
                save_fn = save_fn + '.zip'
            elif compress_type == 'gzip':
                save_fn = save_fn + '.gz'
            
        return save_fn
        
    def write_asc_file(self, save_fn=None, chunk_size=1024, str_fmt='%15.7e', 
                       full=True, compress=False, save_dir=None, 
                       compress_type='zip'):
        """
        Write an ascii file in the USGS ascii format.
        
        :param save_fn: full path to file name to save the merged ascii to
        :type save_fn: string
        
        :param chunck_size: chunck size to write file in blocks, larger numbers
                            are typically slower. 
        :type chunck_size: int
        
        :param str_fmt: format of the data as written
        :type str_fmt: string
        
        :param full: write out the complete file, mostly for testing.
        :type full: boolean [ True | False ]
        
        :param compress: compress file 
        :type compress: boolean [ True | False ]
        
        :param compress_type: compress file using zip or gzip
        :type compress_type: boolean [ zip | gzip ]
        """
        # get the filename to save to
        save_fn = self._make_file_name(save_path=save_dir, 
                                       compression=compress,
                                       compress_type=compress_type)
        # get the number of characters in the desired string
        s_num = int(str_fmt[1:str_fmt.find('.')])
        
        # convert electric fields into mV/km
        self.convert_electrics()
        
        print('==> {0}'.format(save_fn))
        print('START --> {0}'.format(time.ctime()))
        st = datetime.datetime.now()
        
        # write meta data first
        # sort channel information same as columns
        meta_lines = self.write_metadata(chn_list=[c.capitalize() for c in self.ts.columns])
        if compress == True and compress_type == 'gzip':
            with gzip.open(save_fn, 'wb') as fid:
                h_line = [''.join(['{0:>{1}}'.format(c.capitalize(), s_num) 
                          for c in self.ts.columns])]
                fid.write('\n'.join(meta_lines+h_line) + '\n')
                
                # write out data
                if full is False:
                    out = np.array(self.ts[0:chunk_size])
                    out[np.where(out == 0)] = float(self.MissingDataFlag)
                    out = np.char.mod(str_fmt, out)
                    lines = '\n'.join([''.join(out[ii, :]) for ii in range(out.shape[0])])
                    fid.write(lines+'\n')
                    print('END --> {0}'.format(time.ctime()))
                    et = datetime.datetime.now()
                    write_time = et-st
                    print('Writing took: {0} seconds'.format(write_time.total_seconds()))
                    return 
                
                for chunk in range(int(self.ts.shape[0]/chunk_size)):
                    out = np.array(self.ts[chunk*chunk_size:(chunk+1)*chunk_size])
                    out[np.where(out == 0)] = float(self.MissingDataFlag)
                    out = np.char.mod(str_fmt, out)
                    lines = '\n'.join([''.join(out[ii, :]) for ii in range(out.shape[0])])
                    fid.write(lines+'\n')

        else:
            if compress == True and compress_type == 'zip':
                print('ZIPPING')
                save_fn = save_fn[0:-4]
                zip_file = True
                print(zip_file)
            with open(save_fn, 'w') as fid:
                h_line = [''.join(['{0:>{1}}'.format(c.capitalize(), s_num) 
                          for c in self.ts.columns])]
                fid.write('\n'.join(meta_lines+h_line) + '\n')
                
                # write out data
                if full is False:
                    out = np.array(self.ts[0:chunk_size])
                    out[np.where(out == 0)] = float(self.MissingDataFlag)
                    out = np.char.mod(str_fmt, out)
                    lines = '\n'.join([''.join(out[ii, :]) for ii in range(out.shape[0])])
                    fid.write(lines+'\n')
                    print('END --> {0}'.format(time.ctime()))
                    et = datetime.datetime.now()
                    write_time = et-st
                    print('Writing took: {0} seconds'.format(write_time.total_seconds()))
                    return 
                
                for chunk in range(int(self.ts.shape[0]/chunk_size)):
                    out = np.array(self.ts[chunk*chunk_size:(chunk+1)*chunk_size])
                    out[np.where(out == 0)] = float(self.MissingDataFlag)
                    out = np.char.mod(str_fmt, out)
                    lines = '\n'.join([''.join(out[ii, :]) for ii in range(out.shape[0])])
                    fid.write(lines+'\n')
                    
        # for some fucking reason, all interal variables don't exist anymore
        # and if you try to do the zipping nothing happens, so have to do
        # it externally.  WTF
        print('END -->   {0}'.format(time.ctime()))
        et = datetime.datetime.now()
        write_time = et-st
        print('Writing took: {0} seconds'.format(write_time.total_seconds()))
        
    def write_station_info_metadata(self, save_dir=None, mtft_bool=False):
        """
        write out station info that can later be put into a data base
        
        the data we need is
            - site name
            - site id number
            - lat
            - lon
            - national map elevation
            - hx azimuth
            - ex azimuth
            - hy azimuth
            - hz azimuth
            - ex length
            - ey length
            - start date
            - end date
            - instrument type (lp, bb)
            - number of channels
            
        """
        if save_dir is not None:
            save_fn = os.path.join(save_dir, 
                                   '{0}_{1}T{2}_{3:.0f}.cfg'.format(self.SiteID,
                                    self._start_time.strftime('%Y-%m-%d'),
                                    self._start_time.strftime('%H%M%S'),
                                    self.AcqSmpFreq))
        else:
            save_fn = os.path.join(self.station_dir, 
                                       '{0}_{1}T{2}_{3:.0f}.cfg'.format(self.SiteID,
                                        self._start_time.strftime('%Y-%m-%d'),
                                        self._start_time.strftime('%H%M%S'),
                                        self.AcqSmpFreq))
        meta_dict = {}
        key = '{0}_{1}T{2}_{3:.0f}'.format(self.SiteID,
                                    self._start_time.strftime('%Y-%m-%d'),
                                    self._start_time.strftime('%H%M%S'),
                                    self.AcqSmpFreq)
        meta_dict[key] = {}
        meta_dict[key]['site'] = self.SiteID
        meta_dict[key]['lat'] = self._latitude
        meta_dict[key]['lon'] = self._longitude
        meta_dict[key]['elev'] = self.SiteElevation
        meta_dict[key]['mtft_file'] = mtft_bool
        try:
            meta_dict[key]['hx_azm'] = self.channel_dict['Hx']['Azimuth']
            meta_dict[key]['hx_id'] = self.channel_dict['Hx']['InstrumentID'].split('-')[1]
            meta_dict[key]['hx_nsamples'] = self.channel_dict['Hx']['n_samples']
            meta_dict[key]['hx_ndiff'] = self.channel_dict['Hx']['n_diff']
            meta_dict[key]['hx_std'] = self.channel_dict['Hx']['std']
            meta_dict[key]['hx_start'] = self.channel_dict['Hx']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hx']['InstrumentID'].split('-')[0]
            meta_dict[key]['hx_num'] = self.channel_dict['Hx']['ChnNum'][-1]
        except KeyError:
            meta_dict[key]['hx_azm'] = None
            meta_dict[key]['hx_id'] = None
            meta_dict[key]['hx_nsamples'] = None
            meta_dict[key]['hx_ndiff'] = None
            meta_dict[key]['hx_std'] = None
            meta_dict[key]['hx_start'] = None
            meta_dict[key]['hx_num'] = None
            
        try:
            meta_dict[key]['hy_azm'] = self.channel_dict['Hy']['Azimuth']
            meta_dict[key]['hy_id'] = self.channel_dict['Hy']['InstrumentID'].split('-')[1]
            meta_dict[key]['hy_nsamples'] = self.channel_dict['Hy']['n_samples']
            meta_dict[key]['hy_ndiff'] = self.channel_dict['Hy']['n_diff']
            meta_dict[key]['hy_std'] = self.channel_dict['Hy']['std']
            meta_dict[key]['hy_start'] = self.channel_dict['Hy']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hy']['InstrumentID'].split('-')[0]
            meta_dict[key]['hy_num'] = self.channel_dict['Hy']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['hy_azm'] = None
            meta_dict[key]['hy_id'] = None
            meta_dict[key]['hy_nsamples'] = None
            meta_dict[key]['hy_ndiff'] = None
            meta_dict[key]['hy_std'] = None
            meta_dict[key]['hy_start'] = None
            meta_dict[key]['hy_num'] = None
        try:
            meta_dict[key]['hz_azm'] = self.channel_dict['Hz']['Azimuth']
            meta_dict[key]['hz_id'] = self.channel_dict['Hz']['InstrumentID'].split('-')[1]
            meta_dict[key]['hz_nsamples'] = self.channel_dict['Hz']['n_samples']
            meta_dict[key]['hz_ndiff'] = self.channel_dict['Hz']['n_diff']
            meta_dict[key]['hz_std'] = self.channel_dict['Hz']['std']
            meta_dict[key]['hz_start'] = self.channel_dict['Hz']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hz']['InstrumentID'].split('-')[0]
            meta_dict[key]['hz_num'] = self.channel_dict['Hz']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['hz_azm'] = None
            meta_dict[key]['hz_id'] = None
            meta_dict[key]['hz_nsamples'] = None
            meta_dict[key]['hz_ndiff'] = None
            meta_dict[key]['hz_std'] = None
            meta_dict[key]['hz_start'] = None
            meta_dict[key]['hz_num'] = None
        
        try:
            meta_dict[key]['ex_azm'] = self.channel_dict['Ex']['Azimuth']
            meta_dict[key]['ex_id'] = self.channel_dict['Ex']['InstrumentID']
            meta_dict[key]['ex_len'] = self.channel_dict['Ex']['Dipole_Length']
            meta_dict[key]['ex_nsamples'] = self.channel_dict['Ex']['n_samples']
            meta_dict[key]['ex_ndiff'] = self.channel_dict['Ex']['n_diff']
            meta_dict[key]['ex_std'] = self.channel_dict['Ex']['std']
            meta_dict[key]['ex_start'] = self.channel_dict['Ex']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Ex']['InstrumentID']
            meta_dict[key]['ex_num'] = self.channel_dict['Ex']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['ex_azm'] = None
            meta_dict[key]['ex_id'] = None
            meta_dict[key]['ex_len'] = None
            meta_dict[key]['ex_nsamples'] = None
            meta_dict[key]['ex_ndiff'] = None
            meta_dict[key]['ex_std'] = None
            meta_dict[key]['ex_start'] = None
            meta_dict[key]['ex_num'] = None
        try:
            meta_dict[key]['ey_azm'] = self.channel_dict['Ey']['Azimuth']
            meta_dict[key]['ey_id'] = self.channel_dict['Ey']['InstrumentID']
            meta_dict[key]['ey_len'] = self.channel_dict['Ey']['Dipole_Length']
            meta_dict[key]['ey_nsamples'] = self.channel_dict['Ey']['n_samples']
            meta_dict[key]['ey_ndiff'] = self.channel_dict['Ey']['n_diff']
            meta_dict[key]['ey_std'] = self.channel_dict['Ey']['std']
            meta_dict[key]['ey_start'] = self.channel_dict['Ey']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Ey']['InstrumentID']
            meta_dict[key]['ey_num'] = self.channel_dict['Ey']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['ey_azm'] = None
            meta_dict[key]['ey_id'] = None
            meta_dict[key]['ey_len'] = None
            meta_dict[key]['ey_nsamples'] = None
            meta_dict[key]['ey_ndiff'] = None
            meta_dict[key]['ey_std'] = None
            meta_dict[key]['ey_start'] = None
            meta_dict[key]['ey_num'] = None
        
        meta_dict[key]['start_date'] = self.AcqStartTime
        meta_dict[key]['stop_date'] = self.AcqStopTime
        meta_dict[key]['sampling_rate'] = self.AcqSmpFreq
        meta_dict[key]['n_samples'] = self.AcqNumSmp
        meta_dict[key]['n_chan'] = self.Nchan
        
        
        if meta_dict[key]['zen_num'] in [24, 25, 26, 46, '24', '25', '26', '46',
                                        'ZEN24', 'ZEN25', 'ZEN26', 'ZEN46']:
            meta_dict[key]['collected_by'] = 'USGS'
        else:
            meta_dict[key]['collected_by'] = 'OSU'
        
        # in the old OSU z3d files there are notes in the metadata section
        # pass those on
        meta_dict[key]['notes'] = self.meta_notes
            
        mtcfg.write_dict_to_configfile(meta_dict, save_fn)
        
        return save_fn
        
# =============================================================================
# Functions to help analyze config files
# =============================================================================
class USGScfg(object):
    """
    container to deal with configuration files needed for USGS archiving
    
    =========================== ===============================================
    Attributes                  Description
    =========================== ===============================================
    db                          Pandas dataframe
    std_tol                     Tolerance for std 
    note_names                  look for names to put in notes 
    =========================== ===============================================
    
    =========================== ===============================================
    Methods                     Description
    =========================== ===============================================
    combine_run_cfg             Get information from all .cfg files in  
                                a station directory and put in a Pandas 
                                dataframe
    summarize_runs              Summarize a run database
    make_station_db             Make a station database
    make_location_db            Make a location database
    combine_all_station_info    Get all information from each station and put
                                in a Pandas dataframe
    check_data                  Check to make sure all data are in right format
    check_std                   Check for erroneous standard deviation in TS
    write_shp_file              Make a shapefile with important information
    read_survey_csv             Read survey summary csv file
    get_station_info_from_csv   Get station information from survey summar csv
    =========================== ===============================================
    """
    
    def __init__(self, **kwargs):
        self.db = None
        self.std_tol = .005
        self.note_names = ['_start', '_nsamples']
        
    def combine_run_cfg(self, cfg_dir, write=True):
        """
        Combine all the cfg files for each run into a master spreadsheet
        
        :param cfg_dir: directory to cfg files for a station
                        /home/mtdata/station
        :type cfg_dir: string
        
        :param write: write a csv file summarizing the runs
        :type write: boolean [ True | False ]
        
        .. note:: the station name is assumed to be the same as the folder name
        
        """
        station = os.path.basename(cfg_dir)
            
        cfg_fn_list = sorted([os.path.join(cfg_dir, fn) for fn in os.listdir(cfg_dir)
                              if 'mt' not in fn and 'runs' not in fn 
                              and fn.endswith('.cfg')])
        
        count = 0
        for cfg_fn in cfg_fn_list:
            cfg_dict = mtcfg.read_configfile(cfg_fn)
            if count == 0:
                cfg_db = self.check_db(pd.DataFrame([cfg_dict[cfg_dict.keys()[0]]]))
                count += 1
            else:
                cfg_db = cfg_db.append(self.check_db(pd.DataFrame([cfg_dict[cfg_dict.keys()[0]]])),
                                       ignore_index=True)
                count += 1
        try:
            cfg_db = cfg_db.replace('None', '0')
        except UnboundLocalError:
            return None, None
        
        cfg_db.lat = cfg_db.lat.astype(np.float)
        cfg_db.lon = cfg_db.lon.astype(np.float)
        
        if write:
            csv_fn = os.path.join(cfg_dir, '{0}_runs.csv'.format(station))
            cfg_db.to_csv(csv_fn, index=False)
            return cfg_db, csv_fn
        else:
            return cfg_db, None
        
    def summarize_runs(self, run_db):
        """
        summarize the runs into a single row database
        
        :param run_db: run dataframe
        :type run_db: Pandas Dataframe
        
        :return: station_db
        :rtype: Pandas Dataframe
        """
        print(run_db.site.iloc[0])
        station_dict = pd.compat.OrderedDict()
        station_dict['site_name'] = run_db.site.iloc[0]
        station_dict['siteID'] = run_db.site.iloc[0]
        station_dict['lat'] = run_db.lat.astype(np.float).mean()
        station_dict['lon'] = run_db.lon.astype(np.float).mean()
        station_dict['nm_elev'] = get_nm_elev(station_dict['lat'],
                                              station_dict['lon'])
        
        station_dict['hx_azm'] = run_db.hx_azm.astype(np.float).median()
        station_dict['hy_azm'] = run_db.hy_azm.astype(np.float).median()
        station_dict['hz_azm'] = run_db.hz_azm.astype(np.float).median()
        
        try:
            station_dict['hx_id'] = run_db.hx_id.astype(np.float).median()
            station_dict['hy_id'] = run_db.hy_id.astype(np.float).median()
            station_dict['hz_id'] = run_db.hz_id.astype(np.float).median()
        except ValueError:
            station_dict['hx_id'] = 2254
            station_dict['hy_id'] = 2264
            station_dict['hz_id'] = 2274
        
        station_dict['ex_len'] = run_db.ex_len.astype(np.float).median()
        station_dict['ey_len'] = run_db.ey_len.astype(np.float).median()
        
        station_dict['ex_azm'] = run_db.ex_azm.astype(np.float).median()
        station_dict['ey_azm'] = run_db.ey_azm.astype(np.float).median()
        
        station_dict['ex_num'] = run_db.ex_num.median()
        station_dict['ey_num'] = run_db.ey_num.median()
        station_dict['hx_num'] = run_db.hx_num.median()
        station_dict['hy_num'] = run_db.hy_num.median()
        station_dict['hz_num'] = run_db.hz_num.median()
        
        station_dict['n_chan'] = run_db.n_chan.max()
        
        station_dict['sampling_rate'] = run_db.sampling_rate.astype(np.float).median()
        
        station_dict['zen_num'] = run_db.zen_num.iloc[0]
        
        station_dict['collected_by'] = run_db.collected_by.iloc[0]
        station_dict['notes'] = ''.join([run_db.notes.iloc[ii] for ii in range(len(run_db))])
        station_dict['mtft_file'] = run_db.mtft_file.iloc[0]
        
        station_dict['start_date'] = run_db.start_date.min()
        station_dict['stop_date'] = run_db.stop_date.max()
        
        station_dict['type'] = 'wb'
        station_dict['quality'] = 5
        
        return pd.DataFrame([station_dict])
        
    def make_station_db(self, cfg_db, station):
        """
        Following Danny's instructions make a file with the following 
        information
        # make a single file that summarizes
            * (1) site name
            * (2) siteID
            * (3) lat
            * (4) lon
            * (5) national map elevation
            * (6) Hx azimuth
            * (7) Ex azimuth
            * (8) start date [yyyymmdd]
            * (9) Ex dipole length [m]
            * (10) Ey dipole length [m]
            * (11) wideband channels
            * (12) long period channel
            
        :param cfg_db: summarized dataframe for one station
        :type cfg_db: pandas datafram
        
        :param station: station name
        :type station: string
        
        :return: station database
        :rtype: pandas daraframe
        """

        station_dict = pd.compat.OrderedDict()
        station_dict['site_name'] = station
        station_dict['siteID'] = station
        station_dict['lat'] = cfg_db.lat.astype(np.float).mean()
        station_dict['lon'] = cfg_db.lon.astype(np.float).mean()
        station_dict['nm_elev'] = get_nm_elev(station_dict['lat'],
                                              station_dict['lon'])
        station_dict['hx_azm'] = cfg_db.hx_azm.astype(np.float).median()
        station_dict['ex_azm'] = cfg_db.ex_azm.astype(np.float).median()
        station_dict['start_date'] = cfg_db.start_date.min().split('T')[0].replace('-', '')
        station_dict['ex_len'] = cfg_db.ex_len.astype(np.float).median()
        station_dict['ey_len'] = cfg_db.ey_len.astype(np.float).median()
        station_dict['wb'] = cfg_db.n_chan.astype(np.int).median()
        station_dict['lp'] = 0
        
        s_db = pd.DataFrame([station_dict])
        return s_db
    
    def make_location_db(self, cfg_db, station):
        """
        make a file for Danny including:
                * (1) site name
                * (2) lat
                * (3) lon
                * (4) national map elevation
                * (5) start date
                * (6) end date
                * (7) instrument type [W = wideband, L = long period]
                * (8) quality factor from 1-5 [5 is best]
                
        :param cfg_db: summarized dataframe for one station
        :type cfg_db: pandas datafram
        
        :param station: station name
        :type station: string
        
        :return: station database
        :rtype: pandas daraframe
        """

        loc_dict = pd.compat.OrderedDict()
        loc_dict['Stn_name'] = station
        loc_dict['Lat_WGS84'] = np.round(cfg_db.lat.astype(np.float).mean(), 5)
        loc_dict['Lon_WGS84'] = np.round(cfg_db.lon.astype(np.float).mean(), 5)
        loc_dict['Z_NAVD88'] = np.round(get_nm_elev(loc_dict['Lat_WGS84'],
                                          loc_dict['Lon_WGS84']), 2)
        loc_dict['Start_date'] = cfg_db.start_date.min().split('T')[0].replace('-', '')
        loc_dict['End_date'] = cfg_db.stop_date.max().split('T')[0].replace('-', '')
        loc_dict['Data_type'] = 'W'
        loc_dict['Qual_fac'] = 5
        
        l_db = pd.DataFrame([loc_dict])
        return l_db
    
    def combine_all_station_info(self, survey_dir, skip_stations=None, write=True):
        """
        A convinience function to:
            * combine all cfg files for each run into a single spreadsheet
            * combine all runs to make a station and location file
            * combine all stations into a spreadsheet for Danny
            
        :param survey_dir: directory where all the station folders are
        :type survey_dir: string
        
        :param skip_stations: list of stations to skip
        :type skip_stations: list
        
        """
        if type(skip_stations) is not list:
            skip_stations = [skip_stations]
            
        #s_fn = os.path.join(survey_dir, 'usgs_station_info.csv') 
        l_fn = os.path.join(survey_dir, 'usgs_location_info.csv')
        
        s_count = 0
        for station in os.listdir(survey_dir):
            cfg_dir = os.path.join(survey_dir, station)
            if not os.path.isdir(cfg_dir):
                continue
            if station in skip_stations:
                continue

            # get the database and write a csv file            
            cfg_db, csv_fn = self.combine_run_cfg(cfg_dir)
            if cfg_db is None:
                print('*** No Information for {0} ***'.format(station))
                continue
            
            s_db = self.summarize_runs(cfg_db)
            # get station and location information
            if s_count == 0:
                survey_db = s_db
#                s_db = self.make_station_db(cfg_db, station)
                l_db = self.make_location_db(cfg_db, station)
                s_count += 1
            else:
                survey_db = survey_db.append(s_db, ignore_index=True)
#                s_db = s_db.append(self.make_station_db(cfg_db, station))
                l_db = l_db.append(self.make_location_db(cfg_db, station))
                s_count += 1 
        
#        s_db.to_csv(s_fn, index=False)
        l_db.to_csv(l_fn, index=False)
        
        survey_db.lat = survey_db.lat.astype(np.float)
        survey_db.lon = survey_db.lon.astype(np.float)
#        
#        return csv_fn, s_fn, l_fn
        if write:
            csv_fn = os.path.join(survey_dir, 'survey_summary.csv')
            survey_db.to_csv(csv_fn, index=False)
            return survey_db, csv_fn, l_fn
        else:
            return survey_db, None, l_fn

    def check_data(self, database, name):
        """
        check the columns with base name to make sure all values are equal
        """
        if database.notes[0] in [None, 'none', 'None']:
            database.notes[0] = ''
        elif database.notes[0] == '':
            pass
        elif database.notes[0][-1] != ';':
            database.notes[0] += ';'
        # check start times
        column_bool = database.columns.str.contains(name)
        name_labels = database.columns[column_bool].tolist()
        #value_arr = np.array(database[database.columns[column_bool]],
        #                     dtype=np.int).flatten()
        value_arr = np.array(database[database.columns[column_bool]]).flatten()
        value_arr[value_arr == 'None'] = '0'
        value_arr = value_arr.astype(np.int)
        test = np.array([value_arr[0] == value_arr[ii] for ii in range(value_arr.size)])
        if not (test == True).sum() == test.size:
            if test[0] == True and (test[1:] == False).sum() == test[1:].size:
                database.notes[0] += ' {0} is off;'.format(name_labels[0])
            else:
                for comp, test in zip(name_labels, test):
                    if test == False:
                        database.notes[0] += ' {0} is off;'.format(comp)
        return database
    
    def check_std(self, database):
        """
        check standard deviation for bad data
        """
        if database.notes[0] in [None, 'none', 'None']:
            database.notes[0] = ''
        elif database.notes[0] == '':
            pass
        elif database.notes[0][-1] != ';':
            database.notes[0] += ';'
            
        column_bool = database.columns.str.contains('_std')
        name_labels = database.columns[column_bool].tolist()
        value_arr = np.array(database[database.columns[column_bool]]).flatten()
        value_arr[value_arr == 'None'] = 'nan'
        value_arr = value_arr.astype(np.float)
        for name, value in zip(name_labels, value_arr):
            if value > self.std_tol:
                database.notes[0] += ' {0} is large;'.format(name)
        
        return database
    
    def check_db(self, database):
        """
        convinience function to check all values in database
        """
        for key in self.note_names:
            database = self.check_data(database, key)
        database = self.check_std(database)
        
        return database
    
    def write_shp_file(self, survey_csv_fn, save_path=None):
        """
        write a shape file with important information
        
        :param survey_csv_fn: full path to survey_summary.csv
        :type survey_csf_fn: string
        
        :param save_path: directory to save shape file to
        :type save_path: string
        
        :return: full path to shape files
        :rtype: string
        """
        if save_path is not None:
            save_fn = save_path
        else:
            save_fn = os.path.join(os.path.dirname(survey_csv_fn),
                                   'survey_sites.shp')
    
        survey_db = pd.read_csv(survey_csv_fn)
        geometry = [Point(x, y) for x, y in zip(survey_db.lon, survey_db.lat)]
        crs = {'init':'epsg:4326'}
        survey_db = survey_db.drop(['lat', 'lon'], axis=1)
        survey_db = survey_db.rename(columns={'collected_by':'operator', 
                                              'nm_elev':'elev',
                                              'zen_num':'instr_id'})
        
        # list of columns to take from the database
        col_list = ['siteID', 
                    'elev', 
                    'hx_azm',
                    'hy_azm',
                    'hz_azm', 
                    'hx_id', 
                    'hy_id', 
                    'hz_id', 
                    'ex_len',
                    'ey_len',
                    'ex_azm',
                    'ey_azm', 
                    'n_chan', 
                    'instr_id',
                    'operator',
                    'type',
                    'quality',
                    'start_date',
                    'stop_date']
        
        survey_db = survey_db[col_list]
        
        geo_db = gpd.GeoDataFrame(survey_db, 
                                  crs=crs,
                                  geometry=geometry)
        
        geo_db.to_file(save_fn)
        
        print('*** Wrote survey shapefile to {0}'.format(save_fn))
        return save_fn
    
    def read_survey_csv(self, survey_csv):
        """
        Read in a survey .csv file that will overwrite existing metadata
        parameters.
        
        :param survey_csv: full path to survey_summary.csv file
        :type survey_csv: string
        
        :return: survey summary database
        :rtype: pandas dataframe
        """
        db = pd.read_csv(survey_csv, 
                         dtype={'lat':np.float,
                                'lon':np.float})
        for key in ['hx_id', 'hy_id', 'hz_id']:
            db[key] = db[key].fillna(0)
            db[key] = db[key].astype(np.int)
            
        return db
    
    def get_station_info_from_csv(self, survey_csv, station):
        """
        get station information from a survey .csv file
        
        :param survey_csv: full path to survey_summary.csv file
        :type survey_csv: string
        
        :param station: station name
        :type station: string
        
        :return: station database
        :rtype: pandas dataframe
        
        .. note:: station must be verbatim for whats in summary.
        """
        
        db = self.read_survey_csv(survey_csv)
        try:
            station_index = db.index[db.siteID == station].tolist()[0]
        except IndexError:
            raise ValueError('Could not find {0}, check name'.format(station))
        
        return db.iloc[station_index]

# =============================================================================
# XML data
# =============================================================================
class Person(object):
    """
    Data type for the submitter
    """

    def __init__(self):
        self.name = None
        self.org = None
        self.address = None
        self.city = None
        self.state = None
        self.postal = None
        self.phone = None
        self.email = None
        self.country = None
        self.position = None
        self.science_center = None
        self.program = None
        self.liability = None
        self.funding_source = None
        
class Citation(object):
    """
    Data type for a citation
    """
    
    def __init__(self):
        self._author = None
        self.title = None
        self.journal = None
        self.date = None
        self.issue = None
        self.volume = None
        self.doi_url = None
        self._orcid = []
        
    @property
    def author(self):
        return self._author
    
    @author.setter
    def author(self, author):
        if type(author) is not list:
            author = [author]
        self._author = author
        
    @property
    def orcid(self):
        return self._orcid
    
    @orcid.setter
    def orcid(self, orcid):
        if type(orcid) is not list:
            orcid = [orcid]
        self._orcid = orcid
        
class Survey(object):
    """
    data type to hold survey information
    """
    def __init__(self):
        self.begin_date = None
        self.end_date = None
        self._east = None
        self._west = None
        self._north = None
        self._south = None
        self._elev_min = None
        self._elev_max = None
        
    @property
    def east(self):
        return self._east
    @east.setter
    def east(self, east):
        self._east = float(east)
    
    @property
    def west(self):
        return self._west
    @west.setter
    def west(self, west):
        self._west = float(west)
        
    @property
    def north(self):
        return self._north
    @north.setter
    def north(self, north):
        self._north = float(north)
        
    @property
    def south(self):
        return self._south
    @south.setter
    def south(self, south):
        self._south = float(south)
    
    @property
    def elev_min(self):
        return self._elev_min
    @elev_min.setter
    def elev_min(self, elev_min):
        self._elev_min = float(elev_min)
        
    @property
    def elev_max(self):
        return self._elev_max
    @elev_max.setter
    def elev_max(self, elev_max):
        self._elev_max = float(elev_max)
        
class Processing(object):
    """
    Data type for processing steps
    """
    
    def __init__(self):
        self.step_01 = None
        self.date_01 = None
class Attachment(object):
    """
    Data type for attachments
    """
    
    def __int__(self):
        self.fn = None
        self.description = None

class XMLMetadata(object):
    """
    Container for important information to put in the metadata xml file
    
    The assumed workflow is to make a text file with key = value pairs for the 
    important information.
    
        ###==================================================###
        ### Metadata Configuration File for Science Base XML ###
        ###==================================================###
        
        ### General Information
        authors = [Jared R. Peacock, Kevin Denton, Dave Ponce]
        title = Magnetotelluric data from Mountain Pass, CA
        doi_url = https://doi.org/10.5066/F7610XTR
        
        ### Key words to associate with dataset
        keywords_general = [Magnetotellurics, MT, Time Series, Impedance, 
                            Apparent Resistivity, Phase, Tipper]
        
        ### Key words that are synonymous with USGS
        keywords_thesaurus = [Earth magnetic field, Geophysics, 
                              Electromagnetic surveying]
        
    Then convert that to an XML format.
    
    :Example: ::
        
        >>> import mtpy.usgs.science_base as sb 
        >>> sb_xml = sb.SurveyMetadata()
        >>> sb_xml.read_config_file(r"/home/data/survey_config.txt")
        >>> sb_xml.write_xml_file(r"/home/data/science_base.xml")
    """
    
    def __init__(self, **kwargs):

        self.usgs_str = 'U.S. Geological Survey'
        self.xml_fn = None
        self.authors = None
        self.title = None
        self.doi_url = None
        self.file_types = 'ASCII data, shapefile, image'
        self.journal_citation = Citation()
        self.purpose = None
        self.abstract = None
        self.supplement_info = None
        self.survey = Survey()
        self.submitter = Person()
        
        self.sc_dict = {'GGGSC': 'Geology, Geophysics, and Geochemistry Science Center',
                        'GMEGSC': 'Geology, Minerals, Energy, and Geophysics Science Center'}
        self.program_dict = {'MRP': 'Mineral Resources Program',
                             'VHP': 'Volcano Hazards Program'}
        
        self.keywords_general = None
        self.keywords_thesaurus = None
        self.locations = None
        self.temporal = None
        
        self.use_constraint = None
        self.science_base = Person() 
        self.complete_warning = None
        self.horizontal_accuracy = None
        self.vertical_accuracy = None
        self.processing = Processing()
        self.guide = Attachment()
        self.dictionary = Attachment()
        self.shapefile = Attachment()
        
        self.udom = 'Station identifier of MT sounding used to distinguish '+\
                    'between the soundings associated with this survey.'
        self.lat_def = 'Latitude coordinate of station referenced to the '+\
                       'World Geodetic Service Datum of 1984 (WGS84).'
        self.lon_def = 'Longitude coordinate of station, referenced to the '+\
                       'World Geodetic Service Datum of 1984 (WGS84)'
        self.elev_def = 'Elevation, referenced to the North American '+\
                        'Vertical Datum of 1988 (NAVD 88)'
        
        self.metadata = ET.Element('metadata')
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
    def _get_date(self, date_time_str):
        """
        get the date from date_time_string
        """
        
        try:
            date, time = date_time_str.split('T', 1)
        except ValueError:
            date = date_time_str
        date = date.replace('-', '').replace("'", '')
        return date
    
    def _get_time(self, date_time_str):
        """
        get time string 
        """
        try:
            date, time = date_time_str.split('T', 1)
        except ValueError:
            time = '00:00:00 UTC'
        time, zone = time.split()
        time = time.replace(':', '')
        time += '00Z'

        return time         
        
    def read_config_file(self, config_fn):
        """
        Read in configuration file
        """
        # read in the configuration file
        with open(config_fn, 'r') as fid:
            lines = fid.readlines()
            
        for line in lines:
            # skip comment lines
            if line.find('#') == 0 or len(line.strip()) < 2:
                continue
            # make a key = value pair
            key, value = [item.strip() for item in line.split('=', 1)]
            if value == 'usgs_str':
                value = self.usgs_str
            if value.find('[') >= 0 and value.find(']') >= 0 and value.find('<') != 0:
                value = value.replace('[', '').replace(']', '')
                value = [v.strip() for v in value.split(',')]
            
            # if there is a dot, meaning an object with an attribute separate 
            if key.find('.') > 0:
                obj, obj_attr = key.split('.')
                setattr(getattr(self, obj), obj_attr, value)
            else:
                setattr(self, key, value)
                
    def _set_id_info(self):
        """
        set the ID information
        """
        add_info_str = 'Additional information about Originator: '
        idinfo = ET.SubElement(self.metadata, 'idinfo')
        
        citation = ET.SubElement(idinfo, 'citation')
        citeinfo = ET.SubElement(citation, 'citeinfo')
        for author in self.authors:
            ET.SubElement(citeinfo, 'origin').text = author
        ET.SubElement(citeinfo, 'pubdate').text = datetime.datetime.now().strftime('%Y')
        ET.SubElement(citeinfo, 'title').text = self.title
        ET.SubElement(citeinfo, 'geoform').text = self.file_types
        
        pubinfo = ET.SubElement(citeinfo, 'pubinfo')
        ET.SubElement(pubinfo, 'pubplace').text = 'Denver, CO'
        ET.SubElement(pubinfo, 'publish').text = self.usgs_str 
        ET.SubElement(citeinfo, 'onlink').text = self.doi_url
        # journal publication
        if self.journal_citation:
            journal = ET.SubElement(citeinfo, 'lworkcit')
            jciteinfo = ET.SubElement(journal, 'citeinfo')
            for author in self.journal_citation.author:
                ET.SubElement(jciteinfo, 'origin').text = author
            ET.SubElement(jciteinfo, 'pubdate').text = self._get_date(self.journal_citation.date)
            ET.SubElement(jciteinfo, 'title').text = self.journal_citation.title
            ET.SubElement(jciteinfo, 'geoform').text = 'Publication'
            serinfo = ET.SubElement(jciteinfo, 'serinfo')
            ET.SubElement(serinfo, 'sername').text = self.journal_citation.journal
            ET.SubElement(serinfo, 'issue').text = self.journal_citation.volume
            
            jpubinfo = ET.SubElement(jciteinfo, 'pubinfo')
            ET.SubElement(jpubinfo, 'pubplace').text = 'Denver, CO'
            ET.SubElement(jpubinfo, 'publish').text = self.usgs_str
            # set orcid id #s if given
            if len(self.journal_citation.orcid) > 0:
                orcid_str = ', '.join(['{0}, https://orcid.org/{1}'.format(author, orcnum) 
                                      for author, orcnum in zip(self.journal_citation.author,
                                                                self.journal_citation.orcid)
                                      if orcnum not in [None, 'none', 'None']])
                ET.SubElement(jciteinfo, 'othercit').text = add_info_str+orcid_str
            ET.SubElement(jciteinfo, 'onlink').text = self.journal_citation.doi_url
            
        # description
        description = ET.SubElement(idinfo, 'descript')
        ET.SubElement(description, 'abstract').text = self.abstract
        ET.SubElement(description, 'purpose').text = self.purpose
        ET.SubElement(description, 'supplinf').text = self.supplement_info
        
        # dates covered
        time_period = ET.SubElement(idinfo, 'timeperd')
        time_info = ET.SubElement(time_period, 'timeinfo')
        dates = ET.SubElement(time_info, 'rngdates')
        # start and stop date and time
        ET.SubElement(dates, 'begdate').text = self._get_date(self.survey.begin_date)
        ET.SubElement(dates, 'begtime').text = self._get_time(self.survey.begin_date)
        ET.SubElement(dates, 'enddate').text = self._get_date(self.survey.end_date)
        ET.SubElement(dates, 'endtime').text = self._get_time(self.survey.end_date)
        ET.SubElement(time_period, 'current').text = 'ground condition'

        
        # status
        status = ET.SubElement(idinfo, 'status')
        ET.SubElement(status, 'progress').text = 'Complete'
        ET.SubElement(status, 'update').text = 'As needed'
         
        # extent
        extent = ET.SubElement(idinfo, 'spdom')
        bounding = ET.SubElement(extent, 'bounding')
        for name in ['westbc', 'eastbc', 'northbc', 'southbc']:
            ET.SubElement(bounding, name).text = '{0:.5f}'.format(getattr(self.survey, 
                                                                         name[:-2]))
            
        ### keywords
        keywords = ET.SubElement(idinfo, 'keywords')
        t1 = ET.SubElement(keywords, 'theme')
        ET.SubElement(t1, 'themekt').text = 'None'
        ET.SubElement(t1, 'themekey').text = self.sc_dict[self.submitter.science_center]
        ET.SubElement(t1, 'themekey').text = self.submitter.science_center
        ET.SubElement(t1, 'themekey').text = self.program_dict[self.submitter.program]
        ET.SubElement(t1, 'themekey').text = self.submitter.program
        for kw in self.keywords_general:
            ET.SubElement(t1, 'themekey').text = kw     
        
        # categories
        t2 = ET.SubElement(keywords, 'theme')
        ET.SubElement(t2, 'themekt').text = 'ISO 19115 Topic Categories'
        ET.SubElement(t2, 'themekey').text = 'GeoscientificInformation'
        
        # USGS thesaurus
        t3 = ET.SubElement(keywords, 'theme')
        ET.SubElement(t3, 'themekt').text = 'USGS Thesaurus'
        for kw in self.keywords_thesaurus:
            ET.SubElement(t3, 'themekey').text = kw
        
        # places
        place = ET.SubElement(keywords, 'place')
        ET.SubElement(place, 'placekt').text = 'Geographic Names Information System (GNIS)'
        for loc in self.locations:
            ET.SubElement(place, 'placekey').text = loc
            
        # time periods
        temporal = ET.SubElement(keywords, 'temporal')
        ET.SubElement(temporal, 'tempkt').text = 'None'
        for temp in self.temporal:
            ET.SubElement(temporal, 'tempkey').text = temp
            
        ## constraints
        ET.SubElement(idinfo, 'accconst').text = 'None'
        ET.SubElement(idinfo, 'useconst').text = self.use_constraint
        
        ### contact information
        ptcontact = ET.SubElement(idinfo, 'ptcontac')
        contact_info = ET.SubElement(ptcontact, 'cntinfo')
        c_perp = ET.SubElement(contact_info, 'cntperp')
        ET.SubElement(c_perp, 'cntper').text = self.submitter.name
        ET.SubElement(c_perp, 'cntorg').text = self.submitter.org
        c_address = ET.SubElement(contact_info, 'cntaddr')
        ET.SubElement(c_address, 'addrtype').text = 'Mailing and physical'
        for key in ['address', 'city', 'state', 'postal', 'country']:
            ET.SubElement(c_address, key).text = getattr(self.submitter, key)
        
        ET.SubElement(contact_info, 'cntvoice').text = self.submitter.phone
        ET.SubElement(contact_info, 'cntemail').text = self.submitter.email
        # funding source
        try:
            ET.SubElement(idinfo, 'datacred').text = self.program_dict[self.submitter.funding_source]
        except KeyError:
            ET.SubElement(idinfo, 'datacred').text = self.submitter.funding_source
            
        
    def _set_data_quality(self):
        """
        Set data quality section
        """
        data_quality = ET.SubElement(self.metadata, 'dataqual')

        accuracy = ET.SubElement(data_quality, 'attracc')
        ET.SubElement(accuracy, 'attraccr').text = 'No formal attribute accuracy '+\
                                                   'tests were conducted.'
        ET.SubElement(data_quality, 'logic').text = 'No formal logical accuracy tests'+\
                                                    ' were conducted.'
        ET.SubElement(data_quality, 'complete').text = self.complete_warning
        
        # accuracy
        position_acc = ET.SubElement(data_quality, 'posacc')
        h_acc = ET.SubElement(position_acc, 'horizpa')
        ET.SubElement(h_acc, 'horizpar').text = self.horizontal_accuracy
        v_acc = ET.SubElement(position_acc, 'vertacc')
        ET.SubElement(v_acc, 'vertaccr').text = self.vertical_accuracy
        
        # lineage
        lineage = ET.SubElement(data_quality, 'lineage')
        step_num = len(self.processing.__dict__.keys())/2
        for step in range(1, step_num+1, 1):
            processing_step = ET.SubElement(lineage, 'procstep')
            ET.SubElement(processing_step, 'procdesc').text = getattr(self.processing,
                                                                         'step_{0:02}'.format(step))
            if step == 1:
                ET.SubElement(processing_step, 'procdate').text = self._get_date(self.survey.begin_date)
            else:
                ET.SubElement(processing_step, 'procdate').text = self._get_date(getattr(self.processing,
                                                                             'date_{0:02}'.format(step)))
        
    def _set_spational_info(self):
        """
        set spatial information
        """
        spref = ET.SubElement(self.metadata, 'spref')
        
        horizontal_sys = ET.SubElement(spref, 'horizsys')
        h_geographic = ET.SubElement(horizontal_sys, 'geograph')
        ET.SubElement(h_geographic, 'latres').text = '0.0197305745'
        ET.SubElement(h_geographic, 'longres').text = '0.0273088247'
        ET.SubElement(h_geographic, 'geogunit').text = 'Decimal seconds'
        
        h_geodetic = ET.SubElement(horizontal_sys, 'geodetic')
        ET.SubElement(h_geodetic, 'horizdn').text = 'D_WGS_1984'
        ET.SubElement(h_geodetic, 'ellips').text = 'WGS_1984'
        ET.SubElement(h_geodetic, 'semiaxis').text = '6378137.0'
        ET.SubElement(h_geodetic, 'denflat').text = '298.257223563'
        
    def _set_extra_info(self, station=False):
        """
        set the extra info, still not sure what that stands for
        """
        eainfo = ET.SubElement(self.metadata, 'eainfo')
        
        if station is False:
            detailed = ET.SubElement(eainfo, 'detailed')    
            entry_type = ET.SubElement(detailed, 'enttyp')
            ET.SubElement(entry_type, 'enttypl').text = self.shapefile.fn
            ET.SubElement(entry_type, 'enttypd').text = self.shapefile.description
            ET.SubElement(entry_type, 'enttypds').text = self.usgs_str
            
            entry_attr = ET.SubElement(detailed, 'attr')
            ET.SubElement(entry_attr, 'attrlabl').text = 'Station'
            ET.SubElement(entry_attr, 'attrdef').text = 'Individual station name within MT survey.'
            ET.SubElement(entry_attr, 'attrdefs').text = self.usgs_str
            entry_attr_dom = ET.SubElement(entry_attr, 'attrdomv')
            ET.SubElement(entry_attr_dom, 'udom').text = self.udom 
            
            lat_attr = ET.SubElement(detailed, 'attr')
            ET.SubElement(lat_attr, 'attrlabl').text = 'Lat_WGS84'
            ET.SubElement(lat_attr, 'attrdef').text = self.lat_def
            ET.SubElement(lat_attr, 'attrdefs').text = self.usgs_str
            lat_dom = ET.SubElement(lat_attr, 'attrdomv')
            lat_rdom = ET.SubElement(lat_dom, 'rdom')
            ET.SubElement(lat_rdom, 'rdommin').text = '{0:.5f}'.format(self.survey.south)
            ET.SubElement(lat_rdom, 'rdommax').text = '{0:.5f}'.format(self.survey.north)
            ET.SubElement(lat_rdom, 'attrunit').text = 'Decimal degrees'
            
            lon_attr = ET.SubElement(detailed, 'attr')
            ET.SubElement(lon_attr, 'attrlabl').text = 'Lon_WGS84'
            ET.SubElement(lon_attr, 'attrdef').text = self.lon_def
            ET.SubElement(lon_attr, 'attrdefs').text = self.usgs_str
            lon_dom = ET.SubElement(lon_attr, 'attrdomv')
            lon_rdom = ET.SubElement(lon_dom, 'rdom')
            ET.SubElement(lon_rdom, 'rdommin').text = '{0:.5f}'.format(self.survey.west)
            ET.SubElement(lon_rdom, 'rdommax').text = '{0:.5f}'.format(self.survey.east)
            ET.SubElement(lon_rdom, 'attrunit').text = 'Decimal degrees'
            
            elev_attr = ET.SubElement(detailed, 'attr')
            ET.SubElement(elev_attr, 'attrlabl').text = 'Elev_NAVD88'
            ET.SubElement(elev_attr, 'attrdef').text = self.elev_def
            ET.SubElement(elev_attr, 'attrdefs').text = self.usgs_str
            elev_dom = ET.SubElement(elev_attr, 'attrdomv')
            elev_rdom = ET.SubElement(elev_dom, 'rdom')
            ET.SubElement(elev_rdom, 'rdommin').text = '{0:.1f}'.format(self.survey.elev_min)
            ET.SubElement(elev_rdom, 'rdommax').text = '{0:.1f}'.format(self.survey.elev_max)
            ET.SubElement(elev_rdom, 'attrunit').text = 'Meters'
        
        overview = ET.SubElement(eainfo, 'overview')
        ET.SubElement(overview, 'eaover').text = self.guide.fn
        ET.SubElement(overview, 'eadetcit').text = self.guide.description
        
        overview_02 = ET.SubElement(eainfo, 'overview')
        ET.SubElement(overview_02, 'eaover').text = self.dictionary.fn
        ET.SubElement(overview_02, 'eadetcit').text = self.dictionary.description
        
        
            
    def _set_distribution_info(self):
        """
        write distribution information block
        """

        distinfo = ET.SubElement(self.metadata, 'distinfo')
        
        distribute = ET.SubElement(distinfo, 'distrib')
        center_info = ET.SubElement(distribute, 'cntinfo')
        center_perp = ET.SubElement(center_info, 'cntperp')
        ET.SubElement(center_perp, 'cntper').text = self.science_base.name
        ET.SubElement(center_perp, 'cntorg').text = self.science_base.org
        center_address = ET.SubElement(center_info, 'cntaddr')
        ET.SubElement(center_address, 'addrtype').text = 'Mailing and physical'
        for key in ['address', 'city', 'state', 'postal', 'country']:
            ET.SubElement(center_address, key).text = getattr(self.science_base,
                                                              key)
        ET.SubElement(center_info, 'cntvoice').text = self.science_base.phone
        ET.SubElement(center_info, 'cntemail').text = self.science_base.email
        ET.SubElement(distinfo, 'distliab').text = self.science_base.liability
        
    def _set_meta_info(self):
        """
        set the metadata info section
        """
        meta_info = ET.SubElement(self.metadata, 'metainfo')
        
        ET.SubElement(meta_info, 'metd').text = datetime.datetime.now().strftime('%Y%m%d')
        meta_center = ET.SubElement(meta_info, 'metc')
        
        ### contact information
        meta_contact = ET.SubElement(meta_center, 'cntinfo')
        meta_perp = ET.SubElement(meta_contact, 'cntperp')
        ET.SubElement(meta_contact, 'cntpos').text = self.submitter.position
        ET.SubElement(meta_perp, 'cntper').text = self.submitter.name
        ET.SubElement(meta_perp, 'cntorg').text = self.submitter.org
        meta_address = ET.SubElement(meta_contact, 'cntaddr')
        ET.SubElement(meta_address, 'addrtype').text = 'Mailing and physical'
        for key in ['address', 'city', 'state', 'postal', 'country']:
            ET.SubElement(meta_address, key).text = getattr(self.submitter, key)
        
        ET.SubElement(meta_contact, 'cntvoice').text = self.submitter.phone
        ET.SubElement(meta_contact, 'cntemail').text = self.submitter.email
        
        ET.SubElement(meta_info, 'metstdn').text = 'Content Standard for Digital '+\
                                                        'Geospatial Metadata'
        ET.SubElement(meta_info, 'metstdv').text = 'FGDC-STD-001-1998'
#
    def write_xml_file(self, xml_fn, write_station=False):
        """
        write xml file
        """
        self.xml_fn = xml_fn
        
        self._set_id_info()
        self._set_data_quality()
        self._set_spational_info()
        self._set_extra_info(write_station)
        self._set_distribution_info()
        self._set_meta_info()
        
        # write the xml in a readable format
        xml_str = ET.tostring(self.metadata)
        xml_str = minidom.parseString(xml_str).toprettyxml(indent="    ", 
                                                           encoding='UTF-8')
        with open(self.xml_fn, 'w') as fid:
            fid.write(xml_str)
            
        return self.xml_fn