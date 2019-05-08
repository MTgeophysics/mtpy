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

import h5py
import gzip
import urllib2 as url
import xml.etree.ElementTree as ET

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

# science base
import sciencebasepy as sb

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
    try:
        response = url.urlopen(nm_url.format(lon, lat))
    except url.HTTPError:
        print('GET_ELEVATION_ERROR: Could not connect to internet')
        return -666

    # read the xml response and convert to a float
    info = ET.ElementTree(ET.fromstring(response.read()))
    info = info.getroot()
    for elev in info.iter('Elevation'):
        nm_elev = float(elev.text)
    return nm_elev

# =============================================================================
# schedule
# =============================================================================
class ScheduleDB(object):
    """
    Container for a single schedule item
    """

    def __init__(self, time_series_database, meta_db=None):

        self.ts_db = time_series_database
        self.meta_db = meta_db

    @property
    def start_time(self):
        """
        Start time in UTC string format
        """
        return '{0} UTC'.format(self.ts_db.index[0].isoformat())

    @property
    def stop_time(self):
        """
        Stop time in UTC string format
        """
        return '{0} UTC'.format(self.ts_db.index[-1].isoformat())
    
    @property
    def start_seconds(self):
        """
        Start time in epoch seconds
        """
        return self.ts_db.index[0].to_datetime64().astype(np.int64)/1e9
    
    @property
    def stop_seconds(self):
        """
        sopt time in epoch seconds
        """
        return self.ts_db.index[-1].to_datetime64().astype(np.int64)/1e9

    @property
    def n_chan(self):
        """
        number of channels
        """
        return self.ts_db.shape[1]

    @property
    def sampling_rate(self):
        """
        sampling rate
        """
        return np.round(1.0e9/self.ts_db.index[0].freq.nanos, decimals=1)

    @property
    def n_samples(self):
        """
        number of samples
        """
        return self.ts_db.shape[0]
    
    def write_metadata_csv(self, csv_dir):
        """
        write metadata to a csv file
        """
        
        csv_fn = self._make_csv_fn(csv_dir)
        self.meta_db.to_csv(csv_fn)
            
    def _make_csv_fn(self, csv_dir):
        if not isinstance(self.meta_db, pd.Series):
            raise ValueError('meta_db is not a Pandas Series, {0}'.format(type(self.meta_db)))
        csv_fn = '{0}_{1}_{2}_{3}.csv'.format(self.meta_db.station,
                                              self.ts_db.index[0].strftime('%Y%m%d'),
                                              self.ts_db.index[1].strftime('%H%M%S'),
                                              int(self.sampling_rate))
        
        return os.path.join(csv_dir, csv_fn)
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
        self.verbose = True
        self._pd_dt_fmt = '%Y-%m-%d %H:%M:%S.%f'
        self._meta_dtype = np.dtype([('comp', 'S3'),
                                     ('start', np.int64),
                                     ('stop', np.int64),
                                     ('fn', 'S140'),
                                     ('sampling_rate', np.float32),
                                     ('latitude', np.float32),
                                     ('longitude', np.float32),
                                     ('elevation', np.float32),
                                     ('ch_azimuth', np.float32),
                                     ('ch_length', np.float32),
                                     ('ch_num', np.int32),
                                     ('ch_sensor', 'S10'),
                                     ('n_samples', np.int32),
                                     ('t_diff', np.int32),
                                     ('standard_deviation', np.float32),
                                     ('station', 'S12')])

    def _empty_meta_arr(self):
        """
        Create an empty pandas Series
        """
        dtype_list = [('station', 'S10'),
                      ('latitude', np.float),
                      ('longitude', np.float),
                      ('elevation', np.float),
                      ('start', np.int64),
                      ('stop', np.int64),
                      ('sampling_rate', np.float),
                      ('n_chan', np.int),
                      ('n_samples', np.int),
                      ('instrument_id', 'S10'),
                      ('collected_by', 'S30'),
                      ('notes', 'S200'),
                      ('comp', 'S24')]
        
        for cc in ['ex', 'ey', 'hx', 'hy', 'hz']:
            for name, n_dtype in self._meta_dtype.fields.items():
                if name in ['station', 'latitude', 'longitude', 'elevation',
                            'sampling_rate', 'comp']:
                    continue
                elif 'ch' in name:
                    m_name = name.replace('ch', cc)
                else:
                    m_name = '{0}_{1}'.format(cc, name)
                dtype_list.append((m_name, n_dtype[0]))
        ### make an empy data frame, for now just 1 set.
        df = pd.DataFrame(np.zeros(1, dtype=dtype_list))
        
        ### return a pandas series, easier to access than dataframe
        return df.iloc[0]

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
    
    def merge_z3d(self, fn_list):
        """
        Merge a block of z3d files and fill metadata.
        
        :param fn_list: list of z3d files from same schedule action
        :type fn_list: list of strings
        
        :returns: ScheduleDB object that contains metadata and TS dataframes
        :rtype: ScheduleDB
        
        """
        # length of file list
        n_fn = len(fn_list)
        
        ### get empty series to fill
        meta_db = self._empty_meta_arr()
        
        ### need to get some statistics from the files, sometimes a file can
        ### be corrupt so we can make some of these lists
        lat = np.zeros(n_fn)
        lon = np.zeros(n_fn)
        elev = np.zeros(n_fn)
        station = np.zeros(n_fn, dtype='S12')
        sampling_rate = np.zeros(n_fn)
        zen_num = np.zeros(n_fn, dtype='S4')
        start = []
        stop = []
        n_samples = []
        ts_list = []
        
        print('-'*50)
        for ii, fn in enumerate(fn_list):
            z3d_obj = zen.Zen3D(fn)
            try:
                z3d_obj.read_z3d()
            except zen.ZenGPSError:
                print('xxxxxx BAD FILE: Skipping {0} xxxx'.format(fn))
                continue
            # get the components from the file
            comp = z3d_obj.metadata.ch_cmp.lower()
            # convert the time index into an integer
            dt_index = z3d_obj.ts_obj.ts.data.index.astype(np.int64)/10.**9

            # extract useful data that will be for the full station
            sampling_rate[ii] = z3d_obj.df
            lat[ii] = z3d_obj.header.lat
            lon[ii] = z3d_obj.header.long
            elev[ii] = z3d_obj.header.alt
            station[ii] = z3d_obj.station
            zen_num[ii] = int(z3d_obj.header.box_number)
            
            #### get channel setups
            meta_db['comp'] += '{} '.format(comp)
            meta_db['{0}_{1}'.format(comp, 'start')] = dt_index[0]
            meta_db['{0}_{1}'.format(comp, 'stop')] = dt_index[-1]
            start.append(dt_index[0])
            stop.append(dt_index[-1])
            meta_db['{0}_{1}'.format(comp, 'fn')] = fn
            meta_db['{0}_{1}'.format(comp, 'azmimuth')] = z3d_obj.metadata.ch_azimuth
            if 'e' in comp:
                meta_db['{0}_{1}'.format(comp, 'length')] = z3d_obj.metadata.ch_length
            ### get sensor number
            elif 'h' in comp:
                meta_db['{0}_{1}'.format(comp, 'sensor')] = int(z3d_obj.metadata.ch_number.split('.')[0])
            meta_db['{0}_{1}'.format(comp, 'num')] = ii+1
            meta_db['{0}_{1}'.format(comp,'n_samples')] = z3d_obj.ts_obj.ts.shape[0]
            n_samples.append(z3d_obj.ts_obj.ts.shape[0])
            meta_db['{0}_{1}'.format(comp,'t_diff')] = int((dt_index[-1]-dt_index[0])*z3d_obj.df)-\
                                      z3d_obj.ts_obj.ts.shape[0]
            # give deviation in percent
            meta_db['{0}_{1}'.format(comp,'standard_deviation')] = \
                                100*abs(z3d_obj.ts_obj.ts.std()[0]/\
                                        z3d_obj.ts_obj.ts.median()[0])
            try:
                meta_db['notes'] = z3d_obj.metadata.notes.replace('\r', ' ').replace('\x00', '').rstrip()
            except AttributeError:
                pass
            
            ts_list.append(z3d_obj.ts_obj)

        ### fill in meta data for the station
        meta_db.latitude = self._median_value(lat)
        meta_db.longitude = self._median_value(lon)
        meta_db.elevation = get_nm_elev(meta_db.latitude,
                                        meta_db.longitude)
        meta_db.station = self._median_value(station)
        meta_db.instrument_id = 'ZEN{0}'.format(self._median_value(zen_num))
        meta_db.sampling_rate = self._median_value(sampling_rate)
        
        ### merge time series into a single data base
        sch_obj = self.merge_ts_list(ts_list)
        meta_db.start = sch_obj.start_seconds
        meta_db.stop = sch_obj.stop_seconds
        meta_db.n_chan = sch_obj.n_chan
        meta_db.n_samples = sch_obj.n_samples
        sch_obj.meta_db = meta_db
        
        return sch_obj
    
    def merge_ts_list(self, ts_list, decimate=1):
        """
        Merge time series from a list of TS objects.
        
        Looks for the maximum start time and the minimum stop time to align
        the time series.  Indexed by UTC time.
        
        :param ts_list: list of mtpy.core.ts.TS objects from z3d files
        :type ts_list: list
        
        :param decimate: factor to decimate the data by
        :type decimate: int
        
        :returns: merged time series
        :rtype: pandas DataFrame indexed by UTC time
        """
        comp_list = [ts_obj.component.lower() for ts_obj in ts_list]
        df = ts_list[0].sampling_rate
        dt_index_list = [ts_obj.ts.data.index.astype(np.int64)/10.**9
                         for ts_obj in ts_list]
        
        # get start and stop times
        start = max([dt[0] for dt in dt_index_list])
        stop = min([dt[-1] for dt in dt_index_list])
        
        ### make start time in UTC
        dt_struct = datetime.datetime.utcfromtimestamp(start)
        start_time_utc = datetime.datetime.strftime(dt_struct, self._pd_dt_fmt)

        # figure out the max length of the array, getting the time difference into
        # seconds and then multiplying by the sampling rate
        max_ts_len = int((stop-start)*df)
        ts_len = min([ts_obj.ts.size for ts_obj in ts_list]+[max_ts_len])
        if decimate > 1:
            ts_len /= decimate

        ### make an empty pandas dataframe to put data into, seems like the
        ### fastes way so far.
        ts_db = pd.DataFrame(np.zeros((ts_len, len(comp_list))),
                             columns=comp_list,
                             dtype=np.float32)
        for ts_obj in ts_list:
            comp = ts_obj.component.lower()
            dt_index = ts_obj.ts.data.index.astype(np.int64)/10**9
            index_0 = np.where(dt_index == start)[0][0]
            index_1 = min([ts_len-index_0, ts_obj.ts.shape[0]-index_0])

            ### check to see what the time difference is, should be 0,
            ### but sometimes not, then need to account for that.
            t_diff = ts_len-(index_1-index_0)
            if decimate > 1:
                 ts_db[comp][0:(ts_len-t_diff)/decimate] = \
                                 sps.resample(ts_obj.ts.data[index_0:index_1],
                                              ts_len,
                                              window='hanning')

            else:
                ts_db[comp][0:ts_len-t_diff] = ts_obj.ts.data[index_0:index_1]

        # reorder the columns
        ts_db = ts_db[self.get_chn_order(comp_list)]

        # set the index to be UTC time
        dt_freq = '{0:.0f}N'.format(1./(df)*1E9)
        dt_index = pd.date_range(start=start_time_utc,
                                 periods=ts_db.shape[0],
                                 freq=dt_freq)
        ts_db.index = dt_index
        
        return ScheduleDB(ts_db)

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
        t_arr = np.zeros(n_fn, dtype=self._meta_dtype)
    
        t_arr['ch_num'] = np.arange(1, n_fn+1)

        print('-'*50)
        for ii, fn in enumerate(fn_list):
            z3d_obj = zen.Zen3D(fn)
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
            t_arr[ii]['sampling_rate'] = z3d_obj.df
            t_arr[ii]['latitude'] = z3d_obj.header.lat
            t_arr[ii]['longitude'] = z3d_obj.header.long
            t_arr[ii]['elevation'] = z3d_obj.header.alt
            t_arr[ii]['ch_azimuth'] = z3d_obj.metadata.ch_azimuth
            if 'e' in t_arr[ii]['comp']:
                t_arr[ii]['ch_length'] = z3d_obj.metadata.ch_length
            ### get sensor number
            if 'h' in t_arr[ii]['comp']:
                t_arr[ii]['ch_num'] = int(float(z3d_obj.metadata.ch_number))
            t_arr[ii]['ch_box'] = int(z3d_obj.header.box_number)
            t_arr[ii]['n_samples'] = z3d_obj.ts_obj.ts.shape[0]
            t_arr[ii]['t_diff'] = int((dt_index[-1]-dt_index[0])*z3d_obj.df)-\
                                      z3d_obj.ts_obj.ts.shape[0]
            t_arr[ii]['standard_deviation'] = z3d_obj.ts_obj.ts.std()
            t_arr[ii]['station'] = z3d_obj.station
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

        dt_struct = datetime.datetime.utcfromtimestamp(start)
        start_time_utc = datetime.datetime.strftime(dt_struct, self._pd_dt_fmt)

        # figure out the max length of the array, getting the time difference into
        # seconds and then multiplying by the sampling rate
        max_ts_len = int((stop-start)*df)
        ts_len = min([meta_arr['n_samples'].max(), max_ts_len])

        if decimate > 1:
            ts_len /= decimate

        ### make an empty pandas dataframe to put data into, seems like the
        ### fastes way so far.
        ts_db = pd.DataFrame(np.zeros((ts_len, meta_arr.size)),
                             columns=list(meta_arr['comp']),
                             dtype=np.float32)
        ### loop over each time series and find the earliest time all TS start
        ### and latest time all TS end.
        for ii, m_arr in enumerate(meta_arr):
            z3d_obj = zen.Zen3D(m_arr['fn'])
            z3d_obj.read_z3d()

            dt_index = z3d_obj.ts_obj.ts.data.index.astype(np.int64)/10**9
            index_0 = np.where(dt_index == start)[0][0]
            index_1 = min([ts_len-index_0, z3d_obj.ts_obj.ts.shape[0]-index_0])

            ### check to see what the time difference is, should be 0,
            ### but sometimes not, then need to account for that.
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
        ts_db = ts_db[self.get_chn_order(self.chn_order)]

        # set the index to be UTC time
        dt_freq = '{0:.0f}N'.format(1./(df)*1E9)
        dt_index = pd.date_range(start=start_time_utc,
                                 periods=ts_db.shape[0],
                                 freq=dt_freq)
        ts_db.index = dt_index

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
        
    def _median_value(self, value_array):
        """
        get the median value from a metadata entry
        """
        try:
            return np.median(value_array[np.nonzero(value_array)])
        except TypeError:
            return list(set(value_array))[0]
        
    def convert_metadata_to_db(self, metadata_arr):
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
        meta_dict = {}

        meta_dict = {}
        meta_dict['station'] = self._median_value(metadata_arr['station'])
        meta_dict['latitude'] = self._median_value(metadata_arr['latitude'])
        meta_dict['longitude'] = self._median_value(metadata_arr['longitude'])
        meta_dict['elevation'] = get_nm_elev(meta_dict['latitude'],
                                             meta_dict['longitude'])
        meta_dict['start'] = metadata_arr['start'].max()
        meta_dict['stop'] = metadata_arr['stop'].min()
        meta_dict['sampling_rate'] = self._median_value(metadata_arr['sampling_rate'])
        meta_dict['n_chan'] = len(metadata_arr['comp'])
        meta_dict['n_samples'] = metadata_arr['n_samples'].min()
        meta_dict['instrument_id'] = self._median_value(metadata_arr['ch_box'])
        
        for cc in ['ex', 'ey', 'hx', 'hy', 'hz']:
            try:
                c_find = np.where(metadata_arr['comp'] == cc)[0][0]
                find = True
            except IndexError:
                find = False
            for name in self._meta_dtype.names:
                if name in ['station', 'latitude', 'longitude', 'elevation']:
                    continue
                elif 'ch' in name:
                    m_name = name.replace('ch', cc)
                else:
                    m_name = '{0}_{1}'.format(cc, name)
                if find:
                    meta_dict[m_name] = metadata_arr[c_find][name]
                else:
                    meta_dict[m_name] = None

        if meta_dict['instrument_id'] in [24, 25, 26, 46, '24', '25', '26', '46',
                                          'ZEN24', 'ZEN25', 'ZEN26', 'ZEN46']:
            meta_dict['collected_by'] = 'USGS'
        else:
            meta_dict['collected_by'] = 'OSU'

        # in the old OSU z3d files there are notes in the metadata section
        # pass those on
        meta_dict['notes'] = ''

        return pd.Series(meta_dict)

# =============================================================================
#  Metadata for usgs ascii file
# =============================================================================
class AsciiMetadata(object):
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
class USGSasc(AsciiMetadata):
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
        AsciiMetadata.__init__(self)
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
                                                               int(getattr(station_db,
                                                                       '{0}_num'.format(chn.lower()))))


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
            meta_arr['ch_azimuth'][np.where(meta_arr['comp'] != 'hz')] += self.declination

        # fill channel dictionary with appropriate values
        self.channel_dict = dict([(comp.capitalize(),
                                   {'ChnNum':'{0}{1}'.format(self.SiteID, ii+1),
                                    'ChnID':meta_arr['comp'][ii].capitalize(),
                                    'InstrumentID':meta_arr['ch_box'][ii],
                                    'Azimuth':meta_arr['ch_azimuth'][ii],
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
            meta_dict[key]['hx_azimuth'] = self.channel_dict['Hx']['Azimuth']
            meta_dict[key]['hx_id'] = self.channel_dict['Hx']['InstrumentID'].split('-')[1]
            meta_dict[key]['hx_nsamples'] = self.channel_dict['Hx']['n_samples']
            meta_dict[key]['hx_ndiff'] = self.channel_dict['Hx']['n_diff']
            meta_dict[key]['hx_std'] = self.channel_dict['Hx']['std']
            meta_dict[key]['hx_start'] = self.channel_dict['Hx']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hx']['InstrumentID'].split('-')[0]
            meta_dict[key]['hx_num'] = self.channel_dict['Hx']['ChnNum'][-1]
        except KeyError:
            meta_dict[key]['hx_azimuth'] = None
            meta_dict[key]['hx_id'] = None
            meta_dict[key]['hx_nsamples'] = None
            meta_dict[key]['hx_ndiff'] = None
            meta_dict[key]['hx_std'] = None
            meta_dict[key]['hx_start'] = None
            meta_dict[key]['hx_num'] = None

        try:
            meta_dict[key]['hy_azimuth'] = self.channel_dict['Hy']['Azimuth']
            meta_dict[key]['hy_id'] = self.channel_dict['Hy']['InstrumentID'].split('-')[1]
            meta_dict[key]['hy_nsamples'] = self.channel_dict['Hy']['n_samples']
            meta_dict[key]['hy_ndiff'] = self.channel_dict['Hy']['n_diff']
            meta_dict[key]['hy_std'] = self.channel_dict['Hy']['std']
            meta_dict[key]['hy_start'] = self.channel_dict['Hy']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hy']['InstrumentID'].split('-')[0]
            meta_dict[key]['hy_num'] = self.channel_dict['Hy']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['hy_azimuth'] = None
            meta_dict[key]['hy_id'] = None
            meta_dict[key]['hy_nsamples'] = None
            meta_dict[key]['hy_ndiff'] = None
            meta_dict[key]['hy_std'] = None
            meta_dict[key]['hy_start'] = None
            meta_dict[key]['hy_num'] = None
        try:
            meta_dict[key]['hz_azimuth'] = self.channel_dict['Hz']['Azimuth']
            meta_dict[key]['hz_id'] = self.channel_dict['Hz']['InstrumentID'].split('-')[1]
            meta_dict[key]['hz_nsamples'] = self.channel_dict['Hz']['n_samples']
            meta_dict[key]['hz_ndiff'] = self.channel_dict['Hz']['n_diff']
            meta_dict[key]['hz_std'] = self.channel_dict['Hz']['std']
            meta_dict[key]['hz_start'] = self.channel_dict['Hz']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Hz']['InstrumentID'].split('-')[0]
            meta_dict[key]['hz_num'] = self.channel_dict['Hz']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['hz_azimuth'] = None
            meta_dict[key]['hz_id'] = None
            meta_dict[key]['hz_nsamples'] = None
            meta_dict[key]['hz_ndiff'] = None
            meta_dict[key]['hz_std'] = None
            meta_dict[key]['hz_start'] = None
            meta_dict[key]['hz_num'] = None

        try:
            meta_dict[key]['ex_azimuth'] = self.channel_dict['Ex']['Azimuth']
            meta_dict[key]['ex_id'] = self.channel_dict['Ex']['InstrumentID']
            meta_dict[key]['ex_len'] = self.channel_dict['Ex']['Dipole_Length']
            meta_dict[key]['ex_nsamples'] = self.channel_dict['Ex']['n_samples']
            meta_dict[key]['ex_ndiff'] = self.channel_dict['Ex']['n_diff']
            meta_dict[key]['ex_std'] = self.channel_dict['Ex']['std']
            meta_dict[key]['ex_start'] = self.channel_dict['Ex']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Ex']['InstrumentID']
            meta_dict[key]['ex_num'] = self.channel_dict['Ex']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['ex_azimuth'] = None
            meta_dict[key]['ex_id'] = None
            meta_dict[key]['ex_len'] = None
            meta_dict[key]['ex_nsamples'] = None
            meta_dict[key]['ex_ndiff'] = None
            meta_dict[key]['ex_std'] = None
            meta_dict[key]['ex_start'] = None
            meta_dict[key]['ex_num'] = None
        try:
            meta_dict[key]['ey_azimuth'] = self.channel_dict['Ey']['Azimuth']
            meta_dict[key]['ey_id'] = self.channel_dict['Ey']['InstrumentID']
            meta_dict[key]['ey_len'] = self.channel_dict['Ey']['Dipole_Length']
            meta_dict[key]['ey_nsamples'] = self.channel_dict['Ey']['n_samples']
            meta_dict[key]['ey_ndiff'] = self.channel_dict['Ey']['n_diff']
            meta_dict[key]['ey_std'] = self.channel_dict['Ey']['std']
            meta_dict[key]['ey_start'] = self.channel_dict['Ey']['start']
            meta_dict[key]['zen_num'] = self.channel_dict['Ey']['InstrumentID']
            meta_dict[key]['ey_num'] = self.channel_dict['Ey']['ChnNum'][-1:]
        except KeyError:
            meta_dict[key]['ey_azimuth'] = None
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
#  HDF5 object
# =============================================================================
class USGSHDF5(object):
    """
    Container for HDF5 time series format to store in Science Base.

    """

    def __init__(self, hdf5_fn=None, **kwargs):
        self.hdf5_fn = hdf5_fn
        self.datum = 'WGS84'
        self.coordinate_system = 'Geomagnetic North'
        self.instrument_id = 'mt01'
        self.station = 'mt01'
        self.units = 'mV'
        self.declination = 0.0
        self._attr_list = ['station',
                           'latitude',
                           'longitude',
                           'elevation',
                           'declination',
                           'start',
                           'stop',
                           'datum',
                           'coordinate_system',
                           'units',
                           'instrument_id',
                           'ex_azimuth',
                           'ex_length',
                           'ex_sensor',
                           'ex_num',
                           'ey_azimuth',
                           'ey_length',
                           'ey_sensor',
                           'ey_num',
                           'hx_azimuth',
                           'hx_sensor',
                           'hx_num',
                           'hy_azimuth',
                           'hy_sensor',
                           'hy_num',
                           'hz_azimuth',
                           'hz_sensor',
                           'hz_num']
                
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def update_metadata(self, schedule_obj, csv_fn):
        """
        Update metadata extracted from Z3D files with data from a csv file
        """
        ### get the station data base
        try:
            station_db = get_station_info_from_csv(csv_fn, self.station)
        except ValueError:
            print('Could not find information for {0}'.format(self.station))
            return False
        
        ### loop over the columns assuming they have the same keys as the 
        ### metadata series
        for index in station_db.index:
            try:
                schedule_obj.meta_db[index] = station_db[index]
            except AttributeError:
                continue
        return schedule_obj
        
    def write_hdf5(self, z3d_dir, save_dir=None, hdf5_fn=None, csv_fn=None,
                   compress=True, station=None):
        """
        Write an hdf5 file to archive in science base.
        
        :param z3d_dir: full path to directory of station z3d files
        :type z3d_dir: string
        
        :param hdf5_fn: full path to save the hdf5 file
        :type hdf5_fn: string
        
        :param csv_fn: full path to station csv file to overwrite metadata
        :type csv_fn: string
        
        :param compress: boolean to compress the file
        :type compress: boolean [ True | False ]
        
        :param station: station name if different from the folder name
        :type station: string
        
        :returns: full path to saved hdf5 file
        :rtype: string
        
        .. note:: If hdf5_fn is None, the saved file name will be z3d_dir.hdf5
        
        .. note:: If you want to overwrite metadata you can input a csv file
                  that is formatted with specific headings. 
                  
                  ===================== =======================================
                  name                  description
                  ===================== =======================================
                  station               station name 
                  latitude              latitude of station (decimal degrees)
                  longitude             longitude of station (decimal degrees)
                  hx_azimuth            azimuth of HX (degrees from north=0)
                  hy_azimuth            azimuth of HY (degrees from north=0)
                  hz_azimuth            azimuth of HZ (degrees from horizon=0)
                  ex_azimuth            azimuth of EX (degrees from north=0)
                  ey_azimuth            azimuth of EY (degrees from north=0)
                  hx_sensor             instrument id number for HX
                  hy_sensor             instrument id number for HY
                  hz_sensor             instrument id number for HZ
                  ex_sensor             instrument id number for EX
                  ey_sensor             instrument id number for EY
                  ex_length             dipole length (m) for EX
                  ey_length             dipole length (m) for EX
                  ex_num                channel number of EX
                  ey_num                channel number of EX
                  hx_num                channel number of EX
                  hy_num                channel number of EX
                  hz_num                channel number of EX 
                  instrument_id         instrument id 
                  ===================== =======================================
        """
        if station is not None:
            self.station = station
        else:
            self.station = os.path.basename(z3d_dir)
        if save_dir is not None:
            archive_dir = save_dir
        else:
            archive_dir = z3d_dir
            
        if hdf5_fn is not None:
            self.hdf5_fn = hdf5_fn
        else:
            self.hdf5_fn = os.path.join(archive_dir,
                                        '{0}.hdf5'.format(self.station))
  
        # need to over write existing files
        if os.path.exists(self.hdf5_fn):
            os.remove(self.hdf5_fn)

        st = datetime.datetime.now()

        ### get the file names for each block of z3d files
        zc = Z3DCollection()
        fn_list = zc.get_time_blocks(z3d_dir)

        ### Use with so that it will close if something goes amiss
        with h5py.File(self.hdf5_fn, 'w') as h5_obj:
            lat_list = []
            lon_list = []
            instr_id_list = []
            start_list = []
            stop_list = []
            for ii, fn_block in enumerate(fn_list, 1):
                sch_obj = zc.merge_z3d(fn_block)
                if csv_fn is not None:
                    sch_obj = self.update_metadata(sch_obj, csv_fn)
    
                # get lat and lon for statistics
                lat_list.append(sch_obj.meta_db.latitude)
                lon_list.append(sch_obj.meta_db.longitude)
                instr_id_list.append(sch_obj.meta_db.instrument_id)
    
                for key in self._attr_list:
                    try:
                        h5_obj.attrs[key] = sch_obj.meta_db[key]
                    except TypeError:
                        h5_obj.attrs[key] = 'None'
                    except KeyError:
                        try:
                            h5_obj.attrs[key] = getattr(self, key)
                        except AttributeError:
                            print('\txxx No information for {0}'.format(key))
    
                ### create group for schedule action
                schedule = h5_obj.create_group('schedule_{0:02}'.format(ii))
                ### add metadata
                schedule.attrs['start_time'] = sch_obj.start_time
                schedule.attrs['stop_time'] = sch_obj.stop_time
                schedule.attrs['n_samples'] = sch_obj.n_samples
                schedule.attrs['n_channels'] = sch_obj.n_chan
                schedule.attrs['sampling_rate'] = sch_obj.sampling_rate
                
                ### want to get the earliest and latest times
                start_list.append(sch_obj.start_time)
                stop_list.append(sch_obj.stop_time)
    
                ### add datasets for each channel
                for comp in sch_obj.ts_db.columns:
                    if compress:
                        d_set = schedule.create_dataset(comp, data=sch_obj.ts_db[comp],
                                                        compression='gzip',
                                                        compression_opts=4)
                    else:
                        d_set = schedule.create_dataset(comp, data=sch_obj.ts_db[comp])
                    ### might be good to have some notes, will make space for it
                    d_set.attrs['notes'] = ''
    
                sch_obj.write_metadata_csv(archive_dir)
            ### calculate the lat and lon
            station_lat = np.median(np.array(lat_list, dtype=np.float))
            station_lon = np.median(np.array(lon_list, dtype=np.float))
    
            ### set main attributes
            h5_obj.attrs['station'] = self.station
            h5_obj.attrs['coordinate_system'] = self.coordinate_system 
            h5_obj.attrs['datum'] = self.datum
            h5_obj.attrs['latitude'] = station_lat
            h5_obj.attrs['longitude'] = station_lon
            h5_obj.attrs['elevation'] = get_nm_elev(station_lat, station_lon)
            h5_obj.attrs['instrument_id'] = list(set(instr_id_list))[0]
            h5_obj.attrs['units'] = self.units
            h5_obj.attrs['start'] = sorted(start_list)[0]
            h5_obj.attrs['stop'] = sorted(stop_list)[-1]

        run_df, run_csv = combine_station_runs(archive_dir)
        
        et = datetime.datetime.now()
        t_diff = et-st
        print('Took --> {0:.2f} seconds'.format(t_diff.total_seconds()))

        return self.hdf5_fn
    
    def read_hdf5(self, hdf5_fn):
        """
        Read in an HDF5 time series file and make the attributes easy to access
        
        :param hdf5_fn: full path to hdf5 file
        :type hdf5_fn: string
        
        :returns: h5py.File object with read/write privelages
        :rtype: h5py.File
        
        :Example: ::
            
            >>> import mtpy.usgs.usgs_archive as archive
            >>> h5_obj = archive.USGSHDF5(r"/home/mt/station.hdf5")
            >>> h5_obj.attrs['latitude'] = 46.85930
            >>> h5_obj.close()
        
        """ 
        return h5py.File(hdf5_fn, 'r+')

# =============================================================================
# Functions to analyze csv files
# =============================================================================
def read_pd_series(csv_fn):
    """
    read a pandas series and turn it into a dataframe
    
    :param csv_fn: full path to schedule csv
    :type csv_fn: string
    
    :returns: pandas dataframe 
    :rtype: pandas.DataFrame
    """
    series = pd.read_csv(csv_fn, index_col=0, header=None, squeeze=True)
    
    return pd.DataFrame(dict([(k, [v]) for k, v in zip(series.index,
                                                       series.values)]))

def combine_station_runs(csv_dir):
    """
    combine all scheduled runs into a single data frame
    
    :param csv_dir: full path the station csv files
    :type csv_dir: string
    
    """
    station = os.path.basename(csv_dir)

    csv_fn_list = sorted([os.path.join(csv_dir, fn) for fn in os.listdir(csv_dir)
                          if 'runs' not in fn and fn.endswith('.csv')])

    count = 0
    for ii, csv_fn in enumerate(csv_fn_list):
        if ii == 0:
            run_df = read_pd_series(csv_fn)

        else:
            run_df = run_df.append(read_pd_series(csv_fn), ignore_index=True)
            count += 1
            
    ### replace any None with 0, makes it easier
    try:
        run_df = run_df.replace('None', '0')
    except UnboundLocalError:
        return None, None

    ### make lat and lon floats
    run_df.latitude = run_df.latitude.astype(np.float)
    run_df.longitude = run_df.longitude.astype(np.float)

    ### write combined csv file
    csv_fn = os.path.join(csv_dir, '{0}_runs.csv'.format(station))
    run_df.to_csv(csv_fn, index=False)
    return run_df, csv_fn

def summarize_station_runs(run_df):
    """
    summarize all runs into a single row dataframe to be appended to survey df
    
    :param run_df: combined run dataframe for a single station
    :type run_df: pd.DataFrame
    
    :returns: single row data frame with summarized information
    :rtype: pd.DataFrame
    """
    station_dict = pd.compat.OrderedDict() 
    for col in run_df.columns:
        if col == 'start':
            value = run_df['start'].max()
        elif col == 'stop':
            value = run_df['stop'].min()
        else:
            try:
                value = run_df[col].median()
            except TypeError:
                value = list(set(run_df[col]))[0]
        station_dict[col] = value
        
    return pd.DataFrame(station_dict)

def combine_survey_csv(survey_dir, skip_stations=None):
    """
    Combine all stations into a single data frame
    
    :param survey_dir: full path to survey directory
    :type survey_dir: string
    
    :param skip_stations: list of stations to skip
    :type skip_stations: list
    
    :returns: data frame with all information summarized
    :rtype: pandas.DataFrame
    
    :returns: full path to csv file
    :rtype: string
    """
    
    if not isinstance(skip_stations, list):
        skip_stations = [skip_stations]
        
    count = 0
    for station in os.listdir(survey_dir):
        station_dir = os.path.join(survey_dir, station)
        if not os.path.isdir(station_dir):
            continue
        if station in skip_stations:
            continue
        
        # get the database and write a csv file
        run_df, run_fn = combine_station_runs(station_dir)
        if run_df is None:
            print('*** No Information for {0} ***'.format(station))
            continue
        if count == 0:
            survey_df = summarize_station_runs(run_df)
            count += 1
        else:
            survey_df = survey_df.append(summarize_station_runs(run_df))
            count += 1
            
    survey_df.latitude = survey_df.latitude.astype(np.float)
    survey_df.longitude = survey_df.longitude.astype(np.float)
    
    csv_fn = os.path.join(survey_dir, 'survey_summary.csv')
    survey_df.to_csv(csv_fn, index=False)
    
    return survey_df, csv_fn

def read_survey_csv(survey_csv):
    """
    Read in a survey .csv file that will overwrite existing metadata
    parameters.

    :param survey_csv: full path to survey_summary.csv file
    :type survey_csv: string

    :return: survey summary database
    :rtype: pandas dataframe
    """
    db = pd.read_csv(survey_csv,
                     dtype={'latitude':np.float,
                            'longitude':np.float})
    for key in ['hx_sensor', 'hy_sensor', 'hz_sensor']:
        db[key] = db[key].fillna(0)
        db[key] = db[key].astype(np.int)

    return db

def get_station_info_from_csv(survey_csv, station):
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

    db = read_survey_csv(survey_csv)
    try:
        station_index = db.index[db.station == station].tolist()[0]
    except IndexError:
        raise ValueError('Could not find {0}, check name'.format(station))

    return db.iloc[station_index]

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

        station_dict['hx_azimuth'] = run_db.hx_azimuth.astype(np.float).median()
        station_dict['hy_azimuth'] = run_db.hy_azimuth.astype(np.float).median()
        station_dict['hz_azimuth'] = run_db.hz_azimuth.astype(np.float).median()

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

        station_dict['ex_azimuth'] = run_db.ex_azimuth.astype(np.float).median()
        station_dict['ey_azimuth'] = run_db.ey_azimuth.astype(np.float).median()

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
        station_dict['hx_azimuth'] = cfg_db.hx_azimuth.astype(np.float).median()
        station_dict['ex_azimuth'] = cfg_db.ex_azimuth.astype(np.float).median()
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
                    'hx_azimuth',
                    'hy_azimuth',
                    'hz_azimuth',
                    'hx_id',
                    'hy_id',
                    'hz_id',
                    'ex_len',
                    'ey_len',
                    'ex_azimuth',
                    'ey_azimuth',
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
# Science Base Functions
# =============================================================================
def sb_locate_child_item(sb_session, station, sb_page_id):
    """
    See if there is a child item already for the given station.  If there is
    not an existing child item returns False.

    :param sb_session: sciencebase session object
    :type sb_session: sciencebasepy.SbSession

    :param station: station to archive
    :type station: string

    :param sb_page_id: page id of the sciencebase database to download to
    :type sb_page_id: string

    :returns: page id for the station child item
    :rtype: string or False
    """
    for item_id in sb_session.get_child_ids(sb_page_id):
        ### for some reason there is a child item that doesn't play nice
        ### so have to skip it
        try:
            item_title = sb_session.get_item(item_id, {'fields':'title'})['title']
        except:
            continue
        if station in item_title:
            return item_id

    return False

def sb_sort_fn_list(fn_list):
    """
    sort the file name list to xml, edi, png

    :param fn_list: list of files to sort
    :type fn_list: list

    :returns: sorted list ordered by xml, edi, png, zip files
    """

    fn_list_sort = [None, None, None]
    index_dict = {'xml':0, 'edi':1, 'png':2}

    for ext in ['xml', 'edi', 'png']:
        for fn in fn_list:
            if fn.endswith(ext):
                fn_list_sort[index_dict[ext]] = fn
                fn_list.remove(fn)
                break
    fn_list_sort += sorted(fn_list)

    # check to make sure all the files are there
    if fn_list_sort[0] is None:
        print('\t\t!! No .xml file found !!')
    if fn_list_sort[1] is None:
        print('\t\t!! No .edi file found !!')
    if fn_list_sort[2] is None:
        print('\t\t!! No .png file found !!')

    # get rid of any Nones in the list in case there aren't all the files
    fn_list_sort[:] = (value for value in fn_list_sort if value is not None)

    return fn_list_sort

def sb_session_login(sb_session, sb_username, sb_password=None):
    """
    login in to sb session using the input credentials.  Checks to see if
    you are already logged in.  If no password is given, the password will be
    requested through the command prompt.

    .. note:: iPython shells will echo your password.  Use a Python command
              shell to not have your password echoed.

    :param sb_session: sciencebase session object
    :type sb_session: sciencebasepy.SbSession

    :param sb_username: sciencebase username, typically your full USGS email
    :type sb_username: string

    :param sb_password: AD password
    :type sb_password: string

    :returns: logged in sciencebasepy.SbSession
    """

    if not sb_session.is_logged_in():
        if sb_password is None:
            sb_session.loginc(sb_username)
        else:
            sb_session.login(sb_username, sb_password)
        time.sleep(5)

    return sb_session

def sb_get_fn_list(archive_dir):
    """
    Get the list of files to archive looking for .zip, .edi, .png within the
    archive directory.  Sorts in the order of xml, edi, png, zip

    :param archive_dir: full path to the directory to be archived
    :type archive_dir: string

    :returns: list of files to archive ordered by xml, edi, png, zip

    """

    fn_list = [os.path.join(archive_dir, fn)
               for fn in os.listdir(archive_dir)
               if fn.endswith('.zip') or fn.endswith('.xml') or
               fn.endswith('.edi') or fn.endswith('.png')]

    return sb_sort_fn_list(fn_list)


def sb_upload_data(sb_page_id, archive_station_dir, sb_username,
                   sb_password=None):
    """
    Upload a given archive station directory to a new child item of the given
    sciencebase page.

    .. note:: iPython shells will echo your password.  Use a Python command
              shell to not have your password echoed.


    :param sb_page_id: page id of the sciencebase database to download to
    :type sb_page_id: string

    :param archive_station_dir: full path to the station directory to archive
    :type archive_station_dir: string

    :param sb_username: sciencebase username, typically your full USGS email
    :type sb_username: string

    :param sb_password: AD password
    :type sb_password: string

    :returns: child item created on the sciencebase page
    :rtype: dictionary

    :Example: ::
        >>> import mtpy.usgs.usgs_archive as archive
        >>> sb_page = '521e213451bd21451n'
        >>> child_item = archive.sb_upload_data(sb_page,
                                                r"/home/mt/archive_station",
                                                'jdoe@usgs.gov')
    """
    ### initialize a session
    session = sb.SbSession()

    ### login to session, note if you run this in a console your password will
    ### be visible, otherwise run from a command line > python sciencebase_upload.py
    sb_session_login(session, sb_username, sb_password)

    station = os.path.basename(archive_station_dir)

    ### File to upload
    upload_fn_list = sb_get_fn_list(archive_station_dir)

    ### check if child item is already created
    child_id = sb_locate_child_item(session, station, sb_page_id)
    ## it is faster to remove the child item and replace it all
    if child_id:
        session.delete_item(session.get_item(child_id))
        sb_action = 'Updated'

    else:
        sb_action = 'Created'

    ### make a new child item
    new_child_dict = {'title':'station {0}'.format(station),
                      'parentId':sb_page_id,
                      'summary': 'Magnetotelluric data'}
    new_child = session.create_item(new_child_dict)

    # sort list so that xml, edi, png, zip files
    # upload data
    try:
        item = session.upload_files_and_upsert_item(new_child, upload_fn_list)
    except:
        sb_session_login(session, sb_username, sb_password)
        # if you want to keep the order as read on the sb page,
        # need to reverse the order cause it sorts by upload date.
        for fn in upload_fn_list[::-1]:
            try:
                item = session.upload_file_to_item(new_child, fn)
            except:
                print('\t +++ Could not upload {0} +++'.format(fn))
                continue

    print('==> {0} child for {1}'.format(sb_action, station))

    session.logout()

    return item