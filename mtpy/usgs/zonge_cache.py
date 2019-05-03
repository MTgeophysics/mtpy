# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:38:28 2017

@author: jpeacock
"""
#==============================================================================

import os
import struct
import shutil
import time
from collections import Counter

import numpy as np
import scipy.signal as sps

import mtpy.usgs.zen as zen
import mtpy.utils.filehandling as mtfh
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
        

        for key in list(kwargs.keys()):
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
        
        for key in list(kwargs.keys()):
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
                    print('Read in metadata')
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
                    print('Read in metadata')
                    continue
                
                # if the data type is calibration
                elif int(f_pointer['type']) == 768:
                    cal_obj = Board_Calibration(f_str)
                    cal_obj.read_board_cal_str()
                    
                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, cal_obj)
                    print('Read in calibration')
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
                    print('Read in time series,  # points = {0}'.format(
                                                        self.metadata.ts_npnt))
                    return
                # if the data type is time series
                elif int(f_pointer['type']) == 15:
                    ts = np.fromstring(f_str, dtype=np.int32)
                    

                
                    key = self._type_dict[int(f_pointer['type'])]        
                    setattr(self, key, ts)
                    print('Read in other')
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
                print('***SKIPPING {0} '.format(zt.fn))
                print('   because it does not contain correct gps time')
                print('   {0} --> {1}'.format(time_max, 
                                             zt.get_date_time(zt.gps_week, 
                                                             time_max)))    
        
        #change data by amount needed        
        for ii, zt in zip(list(skip_dict.keys()), zt_list):
            if skip_dict[ii] != 0:
                skip_points = skip_dict[ii]*zt.df
                print('Skipping {0} points for {1}'.format(skip_points,
                                                            zt.ch_cmp))
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
                print('TS length for channel {0} '.format(zt.ch_number)+\
                      '({0}) '.format(zt.ch_cmp)+\
                      '= {0}'.format(len(ts_trim)))
                print('    T0 = {0}\n'.format(zt.date_time[0]))
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
        print(fn_list)
            
        n_fn = len(fn_list)
        self.zt_list = []
        for fn in fn_list:
            zt1 = zen.Zen3D(fn=fn)
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
                             for key in np.sort(list(self.meta_data.keys()))])
        
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
            print('Saved File to: ', self.save_fn)
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
                             for key in np.sort(list(self.meta_data.keys()))
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
        
        print('Rewrote {0}\n to {1}'.format(self.save_fn, self.save_fn_rw))        
    
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
                print('Index for second navigation length is {0}'.format(ii))
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
                print('Index for second meta length is {0}'.format(ii))
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
                print('Index for second navigation length is {0}'.format(ii))
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
                print('Index for second meta length is {0}'.format(ii))
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
                print('Index for second cal length is {0}'.format(ii))
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
            print('Trimming TS by {0} points'.format(self.ts.shape[0]%num_chn))
        self.ts = np.resize(self.ts, (int(self.ts.shape[0]/num_chn), num_chn))
        
        ii = int(jj)
        jj = ii+4
        ts_len_check = np.fromstring(cdata[ii:jj], dtype=np.int32)
        if ts_len_check != ts_block['len']:
            if self.verbose:
                print('Index for second ts length is {0}'.format(ii))
            raise CacheTimeSeriesError('ts length in data blocks are not'
                                       'equal: {0} != {1}'.format(
                                       ts_block['len'], ts_len_check))

    
#==============================================================================
# 
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
        print('moved {0} to {1}'.format(fn, new_fn))
        
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
    
    fn_list = zen.copy_from_sd(station, **cpkwargs)
    
    #--> merge files into cache files
    mfn_list = merge_3d_files(fn_list, save_path=merge_save_path)
    
    return mfn_list

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
    time_list = list(time_counts.keys())
    
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
    print('copied {0} to {1}'.format(calibration_fn, copy_cal_fn))
    
    print('Start time: {0}'.format(start_time))
    print('End time:   {0}'.format(end_time))
    
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
