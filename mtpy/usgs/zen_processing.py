# -*- coding: utf-8 -*-
"""
ZEN PROCESSING TOOLS
=======================

    * Interface with BIRRP
    * TODO: interface with EMTF
    * Prepare files for Zonge Processing codes

Created on Fri Sep 16 14:29:43 2016

@author: jpeacock

"""
#==============================================================================
import numpy as np
import time
import datetime
import os
import sys
#from io import StringIO
from io import BytesIO


import mtpy.utils.filehandling as mtfh
import mtpy.processing.birrp as birrp
import mtpy.utils.configfile as mtcfg
import mtpy.utils.exceptions as mtex
import mtpy.imaging.plotnresponses as plotnresponses
import mtpy.imaging.plotresponse as plotresponse
import mtpy.usgs.zen as zen
import mtpy.core.edi as mtedi
import mtpy.core.ts as mtts

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator

#==============================================================================


#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'
datetime_sec = '%Y-%m-%d %H:%M:%S'
#==============================================================================

#==============================================================================
# make a class to deal with birrp inputs
#==============================================================================
class BIRRP_processing(birrp.BIRRP_Parameters):
    """
    configuration file for birrp processing
    
    Takes a list of a structured np.ndarrays that includes filenames and 
    metadata and creates a dictionary that can be input into 
    birrp.write_script_file.
    
    Methods
    ----------
        *get_calibrations* : get the coil calibrations from files that are
                             setup as frequency,amplitude,phase
                             
        *get_processing_dict* : create a processing dictionary with all the 
                                appropriate key words for 
    
    """
    
    def __init__(self, **kwargs):

        self.deltat = 256
        super(BIRRP_processing, self).__init__(**kwargs)

        self.calibration_path = kwargs.pop('calibration_path', 
                                         r"d:\Peacock\MTData\Ant_calibrations")
        self.calibration_list = ['2254', '2264', '2274', '2284', '2294',
                                '2304', '2314', '2324', '2334', '2344',
                                '2844', '2854']
        self._max_nread = 16000000
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
                                
 
        
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
        
    def get_processing_dict(self, fn_birrp_list, hx=2284, hy=2284, hz=2284,
                            **kwargs):
        """
        from fn_birrp_arr make a processing dictionary to input into writing
        a birrp script file
        
        fn_birrp_list = fn_birrp_arr[df]
        """
        for key in kwargs:
            setattr(self, key, kwargs[key])
            
        self.fn_list = []
        self.rrfn_list = []
        # need to sort the fn list so that the files are in the correct
        # order for input and output as defined by birrp
        self.nread = []
        self.nskip = []
        self.nskipr = []
        self.nread = []
        #print fn_birrp_list
        for block_arr in fn_birrp_list:
            s_list = np.zeros(len(block_arr), dtype='|S100')
            r_list = np.zeros(2, dtype='|S100')
            #print '-'*50+'\n'+str(len(block_arr))
            if len(block_arr) == 5 or len(block_arr) == 7:
                comp_dict = {'ex':0, 'ey':1, 'hz':2, 'hx':3, 'hy':4, 
                             'rrhx':0, 'rrhy':1}
                self.nout = 3
            elif len(block_arr) == 4 or len(block_arr) == 6:
                comp_dict = {'ex':0, 'ey':1, 'hx':2, 'hy':3, 
                             'rrhx':0, 'rrhy':1}
                self.nout = 2
            # get the time to start, number of points to read for each
            # segment.
            start_dt_list = sorted(list(set(block_arr['start_dt'])))
            end_dt_list = sorted(list(set(block_arr['end_dt'])))
            # if there is only one date time for the block, easy
            if len(start_dt_list) == 1:
                # skip the header of each file
                nskip = 0
                nskipr = 0
                start_dt = time.mktime(time.strptime(start_dt_list[0],
                                                     datetime_fmt)) 
                                                
            else:
                s_dt_list = np.array([time.mktime(time.strptime(sdt, datetime_fmt))
                                      for sdt in start_dt_list])
                
                # get the latest starting time
                start_dt = s_dt_list.max()
            
            # get the ending time
            if len(end_dt_list) == 1:
                end_dt = time.mktime(time.strptime(end_dt_list[0],
                                                   datetime_fmt)) 
            else:
                e_dt_list = np.array([time.mktime(time.strptime(edt, datetime_fmt))
                                      for edt in end_dt_list])
                                          
                # get the earliest ending time
                end_dt = e_dt_list.min()
            
            # calculate the number of data points to skip for the station
            # data and the remote reference data
            nskip = 0
            nskipr = 0
            for fn_arr in block_arr:
                # figure out which component the file is
                comp = fn_arr['comp'].lower()
                st_time = time.mktime(time.strptime(fn_arr['start_dt'],
                                                    datetime_fmt))

                if comp.find('rr') == -1:
                    s_list[comp_dict[comp]] = fn_arr['fn']
                    
                    # fill remote references just in case there are none
                    if comp in ['hx', 'hy']:
                        r_list[comp_dict['rr{0}'.format(comp)]] = fn_arr['fn']
                    
                    # calculate the number of points to skip given the 
                    # start time
                    nskip_ii = (start_dt-st_time)*abs(self.deltat) 
                    if nskip_ii > nskip:
                        nskip = nskip_ii
                else:
                    r_list[comp_dict[comp]] = fn_arr['fn']
                    nskipr_ii = (start_dt-st_time)*abs(self.deltat)
                    if nskipr_ii > nskipr:
                        nskipr = nskipr_ii
            
            nread_ii = min(self._max_nread, 
                           int((end_dt-start_dt)*abs(self.deltat)))
                           
            # need to check the length of the block.  If it isn't big enough
            # then we need to skip it
            if np.floor(np.log2(nread_ii-nskipr)) > np.log2(self.nfft):
                # need to add a 1 in to skip the header line
                # but for 16 we need to skip the first 3 lines because of the
                # filtering process gives bogus numbers at the beginning
                if self.deltat == -16:
                    self.nskip.append(3+int(nskip))
                    self.nskipr.append(3+int(nskipr))
                    self.nread.append(nread_ii-6)
                else:
                    self.nskip.append(1+int(nskip))
                    self.nskipr.append(1+int(nskipr))
                    self.nread.append(nread_ii)
    
                # append file names appropriately removing any empty names
                self.fn_list.append(list(s_list[np.where(s_list != '')]))
                if len(r_list[np.where(r_list != '')]) > 0:
                    self.rrfn_list.append(list(r_list[np.where(r_list != '')]))
                    
            else:
                print(('Not enough points {0}'.format(nread_ii)))
                print(('skipping time block {0} for sampling rate {1}'.format(
                      time.strftime(datetime_fmt, time.localtime(start_dt)),
                      -self.deltat)))

        # need to check for the number of points to be read in, there is 
        # a memory max, for this computer the max is self._max_nread        
        if sum(self.nread) > self._max_nread:
            self.nread[-1] = self._max_nread-sum(self.nread[0:-1])
            print('processing {0} points'.format(sum(self.nread)))
        
        try:
            self.mcomps = len(self.fn_list[0])
        except IndexError:
            raise IOError('''Time blocks are not long enough to process, 
                             consider changing BIRRP parameters''')
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
            print('Did not find HX calibration for {0}'.format(hx))
            self.hx_cal = cal_dict['2284'] 
            self.rrhx_cal = cal_dict['2284'] 
            print('Setting calibration coil number to 2284 as default.')            
        #--> HY                                                 
        try:
            self.hy_cal = cal_dict[str(hy)]
            self.rrhy_cal = cal_dict[str(hy)]
        except KeyError:
            print('Did not find HY calibration for {0}'.format(hy))
            self.hy_cal = cal_dict['2284'] 
            self.rrhy_cal = cal_dict['2284'] 
            print('Setting calibration coil number to 2284 as default.')            
        #--> HZ                                                 
        try:
            self.hz_cal = cal_dict[str(hz)]
        except KeyError:
            print('Did not find HZ calibration for {0}'.format(hz))
            self.hz_cal = cal_dict['2284'] 
            print('Setting calibration coil number to 2284 as default.')
            
        return self.__dict__
                        
    def read_config_file(self, birrp_config_fn):
        """
        read in a configuration file and fill in the appropriate parameters
        """
        
        birrp_dict = mtcfg.read_configfile(birrp_config_fn)
        
        for birrp_key in list(birrp_dict.keys()):
            setattr(self, birrp_key, birrp_dict[birrp_key])
        
#==============================================================================
# Survey configuration file
#==============================================================================
class Survey_Config(object):
    """
    survey config class
    
    will setup a survey configuration file that has the form of:
    
    [station]
        b_instrument_amplification = 1
        b_instrument_type = coil
        b_logger_gain = 1
        b_logger_type = zen
        b_xaxis_azimuth = 0
        b_yaxis_azimuth = 90
        box = 26
        date = 2015/06/09
        e_instrument_amplification = 1
        e_instrument_type = Ag-Agcl electrodes
        e_logger_gain = 1
        e_logger_type = zen
        e_xaxis_azimuth = 0
        e_xaxis_length = 100
        e_yaxis_azimuth = 90
        e_yaxis_length = 100
        elevation = 2113.2
        hx = 2274
        hy = 2284
        hz = 2254
        lat = 37.7074236995
        location = Earth
        lon = -118.999542099
        network = USGS
        notes = Generic config file
        rr_box = 25
        rr_date = 2015/06/09
        rr_hx = 2334
        rr_hy = 2324
        rr_lat = 37.6909139779
        rr_lon = -119.028707542
        rr_station = 302
        sampling_interval = all
        save_path = \home\mtdata\survey_01\mt_01
        station = 300
        station_type = mt
    """
    def __init__(self, **kwargs):
        self.b_instrument_amplification = 1
        self.b_instrument_type = 'induction coil'
        self.b_logger_gain = 1
        self.b_logger_type = 'zen'
        self.b_xaxis_azimuth = 0
        self.b_yaxis_azimuth = 90
        self.box = 24
        self.date = '01/01/00'
        self.e_instrument_amplification = 1
        self.e_instrument_type = 'Ag-Agcl electrodes'
        self.e_logger_gain = 1
        self.e_logger_type = 'zen'
        self.e_xaxis_azimuth = 0
        self.e_xaxis_length = 100
        self.e_yaxis_azimuth = 90
        self.e_yaxis_length = 100
        self.elevation = 0.0
        self.hx = 2324
        self.hy = 2314
        self.hz = 2334
        self.lat = 0.0
        self.location = 'Earth'
        self.lon = 0.0
        self.network = 'USGS'
        self.notes = 'Generic config file'
        self.sampling_interval = 'all'
        self.station = 'mb000'
        self.station_type = 'mt'
        self.save_path = None
        
        self.rr_lat = None
        self.rr_lon = None
        self.rr_station = None
        self.rr_date = None
        self.rr_box = None
#        self.survey_config.rr_lat = zd.header.lat

        
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
        
        print('Wrote survey config file to {0}'.format(fn))
        
        return fn

#==============================================================================
# Z3D files to EDI using BIRRP
#==============================================================================
class Z3D_to_edi(object):
    """
    go from z3d files to .edi using BIRRP as the processing code
    
    Arguments
    ---------------
    
        **station_dir** : string
                          full path to station directory where the Z3D files
                          to be processed are.
                          
        **rr_station_dir** : string
                            full path to remote reference directory where 
                            Z3D files exist.
                            
        **birrp_exe** : string
                        full path to birrp executable
                        
        **coil_cal_path** : string
        
    Methods
    ----------
    
        **make_survey_config_file : make a survey configuration file from
                                    a dictionary of information
        
        **get_z3d_fn_blocks** : get z3d files in blocks of sampling rate and 
                                date
                                
        **make_mtpy_ascii_files** : make mtpy ascii files from given blocks 
                                    of file names
                                    
        **get_schedules_fn_from_dir** : get mtpy ascii file names in blocks
                                        by sampling rate and date from a TS 
                                        directory.
                                        
        **get_schedules_fn** : get mtpy ascii file names in blocks
                               by sampling rate and date from an array.
                               
        **write_script_files** : write birrp script files for each sampling
                                 rate schedule block
                                 
        **run_birrp** : run birrp from a script file and write an .edi file
        
        **write_edi_file** : write edi file from birrp outputs and given 
                            survey parameters.
                            
        **plot_responses** : plots all edi files output from each sampling rate
        
        **process_data** : a convinience function to go from Z3D files to 
                           edi files.
                           
        **combine_edi_files** : combine all edi files from each sampling rate
                                given a frequency range for each.
                                
    Example
    ------------
    
        >>> import mtpy.usgs.zen_processing as zp
        >>> zp_obj = zp.Z3D_to_edi()
        >>> zp_obj.station_dir = r"/home/data/mt01"
        >>> zp_obj.rr_station_dir = r"/home/data/mt02"
        >>> zp_obj.birrp_exe = r"/home/bin/birrp52"
        >>> zp_obj.coil_cal_path = r"/home/data/ant_calibration"
        >>> plot_obj, comb_edi = zp_obj.process_data(df_list=[4096, 256, 16])
                                        
        
    """
    
    def __init__(self, station_dir=None, **kwargs):
        
        self.station_dir = station_dir
        self.rr_station_dir = kwargs.pop('rr_station_dir', None)
        self.survey_config = Survey_Config(save_path=self.station_dir)
        self.survey_config_fn = None
        self.birrp_config_fn = None
        self.birrp_exe = r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts.exe"
        self.coil_cal_path = r"d:\Peacock\MTData\Ant_calibrations\rsp_cal"
        self._coil_calibration_list = ['2254', '2264', '2274', '2284', '2294',
                                       '2304', '2314', '2324', '2334', '2344',
                                       '2844', '2854']
        self.num_comp = 5
        self.df_list = [4096, 256, 16]
        self.max_blocks = 3
        
        # data types for different aspects of getting information
        self._ts_fn_dtype = np.dtype([('station','S6'), 
                                      ('npts', np.int), 
                                      ('df', np.int),
                                      ('start_dt', 'S22'),
                                      ('end_dt', 'S22'),     
                                      ('comp', 'S2'),
                                      ('fn', 'S100'),
                                      ('calibration_fn', 'S100')])
                                   
        self._birrp_fn_dtype = np.dtype([('fn', 'S100'),
                                         ('nread', np.int),
                                         ('nskip', np.int),
                                         ('comp', 'S2'),
                                         ('calibration_fn', 'S100'),
                                         ('rr', np.bool),
                                         ('rr_num', np.int),
                                         ('start_dt', 'S22'),
                                         ('end_dt', 'S22')])
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])

        
    def make_survey_config_file(self, survey_config_dict=None):
        """
        make a survey configuration file from the data
        """
        
        self.survey_config_fn = self.survey_config.write_survey_config_file()
    
    def get_z3d_fn_blocks(self, station_dir, remote=False,
                          df_list=None, max_blocks=None):
        """
        get z3d file names in an array of blocks
        """
            
        if df_list is not None:
            self.df_list = df_list
        
        if max_blocks is not None:
            self.max_blocks = max_blocks

        fn_block_dict = dict([(df, {}) for df in self.df_list])
        fn_count = 0
        for fn in os.listdir(station_dir):
            if fn.lower().endswith('.z3d'):
            
                z3d_fn = os.path.join(station_dir, fn)
                z3d_obj = zen.Zen3D(z3d_fn)
                z3d_obj.read_all_info()
                if remote is True:
                    if not z3d_obj.metadata.ch_cmp.lower() in ['hx', 'hy']:
                        continue
                if z3d_obj.df in self.df_list:
                    fn_count += 1
                    try:
                        fn_block_dict[z3d_obj.df][z3d_obj.zen_schedule].append(z3d_fn)
                    except KeyError:
                        fn_block_dict[z3d_obj.df][z3d_obj.zen_schedule] = [z3d_fn]


        if fn_count == 0:
            raise ValueError('No Z3D files found for in {0}'.format(station_dir))
        else:
            print('Found {0} Z3D files in {1}'.format(fn_count, station_dir))
        
        # check for maximum number of blocks
        for df_key in list(fn_block_dict.keys()):
            date_dict = fn_block_dict[df_key]
            dates = sorted(date_dict.keys())
            if len(dates) == 0:
                print('No Z3D files found for {0} in {1}'.format(str(df_key),
                                                                 station_dir))

            if len(dates) > self.max_blocks:
                for pop_date in dates[-(len(dates)-self.max_blocks):]:
                    fn_block_dict[df_key].pop(pop_date)
        
        return fn_block_dict
        
    def get_calibrations(self, coil_cal_path=None):
        """
        get coil calibrations
        """
        if coil_cal_path is not None:
            self.coil_cal_path = coil_cal_path
            
        if self.coil_cal_path is None or not os.path.isdir(self.coil_cal_path):
            print('could not find calibration path {0}'.format(self.coil_cal_path))
            self.calibration_dict = dict([(cc, '') for cc in 
                                          self._coil_calibration_list])
            return 
            
        calibration_dict = {}
        for cal_fn in os.listdir(self.coil_cal_path):
            for cal_num in self._coil_calibration_list:
                if cal_num in cal_fn:
                    calibration_dict[cal_num] = \
                                    os.path.join(self.coil_cal_path, cal_fn)
                    break
                
        self.calibration_dict = calibration_dict
        
    def make_mtpy_ascii_files(self, station_dir=None, rr_station_dir=None, 
                              notch_dict={4096:{}, 256:None, 16:None}, 
                              df_list=None, max_blocks=None,
                              ex=50., ey=50., coil_cal_path=None): 
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
            
        if rr_station_dir is not None:
            self.rr_station_dir = rr_station_dir
        
        # make rr_station_dir a list so we can loop over it later             
        if type(self.rr_station_dir) is str:
                self.rr_station_dir = [self.rr_station_dir]
                
        if df_list is not None:
            self.df_list = df_list
            
        if max_blocks is not None:
            self.max_blocks = max_blocks
            
        if coil_cal_path is not None:
            self.coil_cal_path = coil_cal_path
            
        # get coil calibrations
        self.get_calibrations() 
        
        #-------------------------------------------------
        # make station z3d's into mtpy ts
        fn_dict = self.get_z3d_fn_blocks(self.station_dir)
        
        # make an empty array to put things
        n_files, n_comps = self._get_num_files(fn_dict)
        fn_arr = self._make_empty_fn_arr(n_files)
        
        # number of components to process        
        self.num_comp = n_comps
        print(' --> processing {0} components'.format(self.num_comp))
                    
        jj = 0  # index for fn_array            
        for df_key in self.df_list:
            try:
                df_notch_dict = notch_dict[df_key]
            except KeyError:
                df_notch_dict = None
                
            for date_key in list(fn_dict[df_key].keys()):
                for fn in fn_dict[df_key][date_key]:
                    fn_arr[jj] = self._convert_z3d_to_mt_ts(fn, df_notch_dict)
                    
                    # account for decimation, need to make a new instance
                    # so as to include the original file as well
                    if fn_arr[jj]['df'] == 256 and 16 in self.df_list:
                        jj += 1
                        fn_arr[jj] = self._convert_z3d_to_mt_ts(fn, 
                                                                df_notch_dict,
                                                                16)  
                    jj += 1
        # get only the non zero entries
        fn_arr = fn_arr[np.nonzero(fn_arr['npts'])]
        
        #--------------------------------------------------
        # Remote Reference z3d's into mtpy ts 
        if rr_station_dir is not None:
            self.rr_station_dir = rr_station_dir
            
        if self.rr_station_dir is not None:
            print('*** Tranforming remote reference Z3D to mtpy format ***')
            
            # get the maximum number of remote reference time series
            # multiply by 3 just to be save
            num_ref = (len(fn_arr)/self.num_comp)*3*len(self.rr_station_dir)
            rr_fn_arr = self._make_empty_fn_arr(num_ref)
            rr = 0
            for rr_dir in self.rr_station_dir:
                rr_fn_dict = self.get_z3d_fn_blocks(rr_dir, remote=True)
           
                for df_key in self.df_list:
                    try:
                        df_notch_dict = notch_dict[df_key]
                    except KeyError:
                        df_notch_dict = None
                        
                    for date_key in list(rr_fn_dict[df_key].keys()):
                        for fn in rr_fn_dict[df_key][date_key]:
                            rr_fn_arr[rr] = self._convert_z3d_to_mt_ts(fn,
                                                                       df_notch_dict,
                                                                       remote=True)
                    
                            # account for decimation, need to make a new instance
                            # so as to include the original file as well
                            if rr_fn_arr[rr]['df'] == 256 and 16 in self.df_list:
                                rr += 1
                                rr_fn_arr[rr] = self._convert_z3d_to_mt_ts(fn, 
                                                                           df_notch_dict,
                                                                           dec=16,
                                                                           remote=True)
                            rr += 1

            # change the time series directory to reflect where the time series
            # are.
            self.rr_station_dir = [os.path.join(rr_dir, 'TS') for rr_dir in 
                                   self.rr_station_dir] 
            
            # get only the non-empty time series
            rr_fn_arr = rr_fn_arr[np.nonzero(rr_fn_arr['npts'])]
        else:
            rr_fn_arr = None
        # change directories to be associated with where the mtpy ts are
        self.station_dir = os.path.join(self.station_dir, 'TS')
        self.survey_config.save_path = self.station_dir
        
        # write survey configuration file
        self.survey_config_fn = self.survey_config.write_survey_config_file()
        
        
        log_lines = []
        for f_arr in fn_arr:
            log_lines += [self._make_config_lines(f_arr)]
        if rr_fn_arr is not None:
            for rr_f_arr in rr_fn_arr:
                log_lines += [self._make_config_lines(rr_f_arr)]
 
        return fn_arr, rr_fn_arr, log_lines
        
    def _get_num_files(self, fn_dict):
        """
        get number of files from fn dict
        """
        num_files = 0
        num_comp = []
        for df_key in list(fn_dict.keys()): 
            for date_key in list(fn_dict[df_key].keys()):
                # in case there is decimation, then we need to double the
                # number of 256 files
                if df_key == 256 and 16 in self.df_list:
                    num_files += len(fn_dict[df_key][date_key])*2
                else:
                    num_files += len(fn_dict[df_key][date_key])
                num_comp.append(len(fn_dict[df_key][date_key]))
   
        num_comp = int(np.mean(np.array(num_comp)))
        return num_files, num_comp
       
    def _make_empty_fn_arr(self, num_files):
        """
        make an empty array with the proper dtype for TS files.
       
        """
        # make an array that has all the information about each file
        fn_arr = np.zeros(num_files, dtype=self._ts_fn_dtype)
                                 
        return fn_arr
        
    def _convert_z3d_to_mt_ts(self, fn, notch_dict=None, dec=1, remote=False):
        """
        convert a z3d file to mtpy ts format and return important information
        
        fn - filename
        notch_dict - dictionary of notches to filter
        dec - decimation level
        """  
        
        return_fn_arr = np.zeros(1, dtype=self._ts_fn_dtype)
        zd = zen.Zen3D(fn)
        zd.read_all_info()
        
        if dec != 1:
            zd.read_z3d()

        #write mtpy mt file
        zd.write_ascii_mt_file(notch_dict=notch_dict, dec=dec)
        
        #create lines to write to a log file                       
        station_num = zd.metadata.rx_xyz0.split(':')[0]
        station_label = zd.metadata.line_name
        station = '{0}{1}'.format(station_label, station_num) 
        
        return_fn_arr['station'] = station
        return_fn_arr['npts'] = zd.ts_obj.n_samples
        return_fn_arr['df'] = zd.df
        return_fn_arr['start_dt'] = zd.zen_schedule
        return_fn_arr['comp'] = zd.metadata.ch_cmp.lower()
        return_fn_arr['fn'] = zd.fn_mt_ascii

        # add metadata to survey configuration file
        if zd.metadata.ch_cmp.lower() in ['hx', 'hy', 'hz']: 
            
            comp = zd.metadata.ch_cmp.lower()
            chn_num = zd.metadata.ch_number
            try:
                cal_fn = self.calibration_dict[chn_num]
            except KeyError:
                print('Did not find calibration for {0}, number {1}'.format(comp, 
                                                                        chn_num)) 
                cal_fn = None
                
            if remote == False:
                setattr(self.survey_config, comp, chn_num)
                setattr(self.survey_config, comp+'_cal_fn', cal_fn)
                
            elif remote == True:
                if hasattr(self.survey_config, 'rr_{0}_00'.format(comp)) == False:
                    setattr(self.survey_config, 
                            'rr_{0}_00'.format(comp),
                            chn_num)
                    setattr(self.survey_config,
                            'rr_{0}_00_cal_fn'.format(comp),
                            cal_fn)
                else:
                    count = 1
#                    while hasattr(self.survey_config, 'rr_{0}_{1:02}'.format(comp, count)):
#                        count += 1
#                        print count 
                    
                    setattr(self.survey_config, 
                            'rr_{0}_{1:02}'.format(comp, count),
                            chn_num)
                    setattr(self.survey_config,
                            'rr_{0}_{1:02}_cal_fn'.format(comp, count),
                            cal_fn)
                    
        elif zd.metadata.ch_cmp.lower() == 'ex':
            self.survey_config.e_xaxis_length = zd.metadata.ch_length
        elif zd.metadata.ch_cmp.lower() == 'ey':
            self.survey_config.e_yaxis_length = zd.metadata.ch_length

        
        if remote == False:
            if self.survey_config.lat is None or self.survey_config.lat == 0.0:
                self.survey_config.lat = zd.header.lat
                self.survey_config.lon = zd.header.long
                self.survey_config.date = zd.schedule.Date
                self.survey_config.box = int(zd.header.box_number)
                self.survey_config.station = station
                self.survey_config.location = zd.metadata.job_name
                self.survey_config.network = zd.metadata.job_by
                self.survey_config.elevation = zd.header.alt
        else:
            if self.survey_config.rr_lat is not None:
                if self.survey_config.rr_lat != zd.header.lat:
                    self.survey_config.rr_lat_01 = zd.header.lat
               
                if self.survey_config.rr_lon != zd.header.long:
                    self.survey_config.rr_lon_01 = zd.header.long
                
                if self.survey_config.rr_date != zd.schedule.Date:
                    self.survey_config.rr_date_01 = zd.schedule.Date
    
                if self.survey_config.rr_station != station:
                    self.survey_config.rr_station_01 = station
                
                if self.survey_config.rr_box != int(zd.header.box_number):
                    self.survey_config.rr_box_01 = int(zd.header.box_number)
            else:
                self.survey_config.rr_lat = zd.header.lat
                self.survey_config.rr_lon = zd.header.long
                self.survey_config.rr_date = zd.schedule.Date
                self.survey_config.rr_box = int(zd.header.box_number)
                self.survey_config.rr_station = station
                self.survey_config.rr_elevation = zd.header.alt

        return return_fn_arr
        
    def _make_config_lines(self, fn_arr):
        
        cfg_list = ['--> station: {0}\n'.format(fn_arr['station']),
                    '    ts_len = {0}\n'.format(fn_arr['npts']),
                    '    df = {0}\n'.format(fn_arr['df']),
                    '    start_dt = {0}\n'.format(fn_arr['start_dt']),
                    '    comp = {0}\n'.format(fn_arr['comp']),
                    '    fn = {0}\n'.format(fn_arr['fn'])]
        if fn_arr['df'] < 256:
            cfg_list.append('    decimated to {0:.0f} samples/s'.format(fn_arr['df']))
                            
        return ''.join(cfg_list)
                   
        
    def get_schedules_fn_from_dir(self, station_ts_dir=None, rr_ts_dir=None,
                                  df_list=None, max_blocks=None,
                                  use_blocks_dict={4096:'all', 256:'all', 16:'all'}):
        """
        get the birrp fn list from a directory of TS files
        """
        
        if station_ts_dir is not None:
            self.station_dir = station_ts_dir
        
        if rr_ts_dir is not None:
            self.rr_station_dir = rr_ts_dir
            
        if type(self.rr_station_dir) is str:
            self.rr_station_dir = [self.rr_station_dir]
            
        if df_list is not None:
            self.df_list = df_list
            
        if max_blocks is not None:
            self.max_blocks = max_blocks
            
        if not os.path.isdir(self.station_dir):
            print('{0} is not a valid directory, check path.'.format(self.station_dir)) 
            return None
            
        fn_arr = np.zeros(len(os.listdir(self.station_dir)),
                          dtype=self._ts_fn_dtype)
        fn_count = 0
        for fn in os.listdir(self.station_dir):
            fn = os.path.join(self.station_dir, fn)
            f_arr, count = self._make_ts_arr_entry(fn)
            if f_arr is None:
                continue
            fn_arr[fn_count] = f_arr
            fn_count += count

        # be sure to trim the array
        fn_arr = fn_arr[np.nonzero(fn_arr['npts'])]   

        #--------------------------------------------------------        
        # get remote reference filenames
        if self.rr_station_dir is not None:
            n_rr = (len(fn_arr)/self.num_comp)*10*len(self.rr_station_dir)
            rr_fn_arr = np.zeros(n_rr, dtype=self._ts_fn_dtype) 
        
            rr_count = 0
        
            for rr_dir in self.rr_station_dir:
                for rr_fn in os.listdir(rr_dir):
                    rr_fn = os.path.join(rr_dir, rr_fn)
                    rr_arr, count = self._make_ts_arr_entry(rr_fn, remote=True)
                    if rr_arr is not None:
                        rr_fn_arr[rr_count] = rr_arr
                    rr_count += count     
            # be sure to trim the array
            rr_fn_arr = rr_fn_arr[np.nonzero(rr_fn_arr['npts'])]

        else:
            rr_fn_arr = None
            
        
        return self.get_schedules_fn(fn_arr, rr_fn_arr,
                                     use_blocks_dict=use_blocks_dict)
#        return (fn_arr, rr_fn_arr)
            
    
    def _make_ts_arr_entry(self, fn, remote=False):
        """
        make ts_arr entries
        """
        fn_ext = os.path.splitext(fn)[1][1:].lower()
        if fn_ext not in ['ex', 'ey', 'hx', 'hy', 'hz']:
            return None, 0
            
        return_fn_arr = np.zeros(1, dtype=self._ts_fn_dtype)
        
        ts_obj = mtts.MT_TS()
        try:
            ts_obj.read_ascii_header(fn)

            if remote == True:
                if ts_obj.component.lower() not in ['hx', 'hy']:
                    return None, 0
                    
            if ts_obj.sampling_rate in self.df_list:
                return_fn_arr['fn'] = fn
                return_fn_arr['npts'] = ts_obj.n_samples
                return_fn_arr['df'] = ts_obj.sampling_rate
                return_fn_arr['comp'] = ts_obj.component.lower()
                    
                start_sec = ts_obj.start_time_epoch_sec
                
                num_sec = float(ts_obj.n_samples)/ts_obj.sampling_rate
                return_fn_arr['start_dt'] = time.strftime(datetime_fmt, 
                                                time.localtime(start_sec)) 
                return_fn_arr['end_dt'] = time.strftime(datetime_fmt, 
                                                time.localtime(start_sec+\
                                                num_sec))
                count = 1
            else:
                count = 0
        except mtts.MT_TS_Error:
            print('  Skipped {0}'.format(fn))
            count = 0
            
        return return_fn_arr, count
        
    def get_schedules_fn(self, fn_arr, rr_fn_arr=None, df_list=None, 
                         max_blocks=None, 
                         use_blocks_dict={4096:'all', 256:'all', 16:'all'}):
        """
        seperate out the different schedule blocks and frequencies so the
        can be processed
        
        Returns
        ---------
            **schedule_fn_dict** : dictionary
                                   keys are sampling rates and values are
                                   lists of file names for each schedule
                                   block up to max blocks
        """
        
        if df_list is not None:
            self.df_list = df_list
            
        if max_blocks is not None:
            self.max_blocks = max_blocks
        
        # get the sampling rates used
        s_keys = set(self.df_list)
        
        # make a dictionary with keys as the sampling rates 
        s_dict = dict([(skey, []) for skey in s_keys])
        
        # loop over the sampling rates and find the schedule blocks
        for df in s_keys:
            # find startind dates for sampling rate
            s_dates = sorted(list(set(fn_arr['start_dt'][np.where(fn_arr['df']==df)])))
            use_blocks = use_blocks_dict[df]            
            if use_blocks == 'all':
                date_list = s_dates[0:self.max_blocks]
            elif type(use_blocks) is list:
                date_list = [s_dates[ii] for ii in use_blocks]
                if len(date_list) > self.max_blocks:
                    date_list = date_list[0:self.max_blocks]
                
            else:
                raise ValueError('Do not understand use_blocks type {0}'.format(type(use_blocks)))
                
            for sdate in date_list:
                s_fn_arr = fn_arr[np.where((fn_arr['start_dt'] == sdate) &
                                            (fn_arr['df'] == df))]
                s_fn_birrp_arr = self._fill_birrp_fn_arr(s_fn_arr) 
                
                # get remote reference information if input
                if rr_fn_arr is not None:
                    # make an empty array to append to
                    rr_birrp_fn_arr = np.zeros(0, dtype=self._birrp_fn_dtype)

                    # find elements where df is the same                    
                    df_find = np.where(rr_fn_arr['df'] == df)
                    
                    # make an array of just the sampling rates                     
                    df_rr_arr = rr_fn_arr[df_find]

                    # first check to see if any of the dates match
                    dt_find = np.where(df_rr_arr['start_dt'] == sdate)
                    
                    # if there are dates that match, fill them into rr array
                    rr_station_find = []
                    if len(dt_find[0]) >= 1:
                        dt_rr_arr = self._fill_birrp_fn_arr(df_rr_arr[dt_find],
                                                            remote=True)

                        rr_birrp_fn_arr = np.append(rr_birrp_fn_arr,
                                                    dt_rr_arr)
                        # need to make a list of remote reference station
                        # already found so there are no duplicates
                        rr_station_find += list(set(['{0}_{1}'.format(os.path.basename(f['fn'])[0:4],
                                                                      f['comp']) 
                                                                      for f in 
                                                                      dt_rr_arr]))

                    # find the where dates do not match
                    dt_not_find = np.where(df_rr_arr['start_dt'] != sdate)
                    
                    # loop over each remote reference file to check for time
                    # that might be close enough, but skip any that have the 
                    # same station name
                    for rr_arr in df_rr_arr[dt_not_find]:
                        # check to see if this station already has been added
                        # as a remote reference
                        f_station = os.path.basename(rr_arr['fn'])[0:4]
                        f_find = '{0}_{1}'.format(f_station, rr_arr['comp'])
                        if f_find in rr_station_find:
                            continue
                        
                        # find the next closest date
                        # estimate the time difference between station date
                        # and remote reference date
                        rr_sec = time.mktime(time.strptime(rr_arr['start_dt'], 
                                                           datetime_fmt))
                        s_sec = time.mktime(time.strptime(sdate, 
                                                          datetime_fmt))
                        # estimate time difference                            
                        t_diff = s_sec-rr_sec
                        
                        # look for the correct time difference
                        dt_arr = None
                        if df == 256:
                            # if the difference is more than 5 hours pass
                            if abs(t_diff) > 5*3600:
                                continue
                            else:
                                print('Using rr {0} for TS starting on {1}'.format(f_station,
                                                                                   rr_arr['start_dt']))
                                print('For station TS starting on      {0}'.format(sdate))
                        
                                dt_arr = self._fill_birrp_fn_arr(rr_arr,
                                                                 remote=True)
                                rr_station_find.append(f_find)
                        elif df == 4096:
                            # if the difference is more than 3 minutes
                            if abs(t_diff) > 4*60:
                                continue
                            else:
                                print('Using rr {0} for TS starting on {1}'.format(f_station,
                                                                                   rr_arr['start_dt']))
                                print('For station TS starting on      {0}'.format(sdate))
                        
                                dt_arr = self._fill_birrp_fn_arr(rr_arr,
                                                                 remote=True)
                                rr_station_find.append(f_find)
                    
                        elif df == 16:
                            # if the difference is more than 3 minutes
                            if abs(t_diff) > 5*3600:
                                continue
                            else:
                                print('Using rr {0} for TS starting on {1}'.format(f_station,
                                                                                   rr_arr['start_dt']))
                                print('For station TS starting on      {0}'.format(sdate))
                        
                                dt_arr = self._fill_birrp_fn_arr(rr_arr,
                                                                 remote=True)
                                rr_station_find.append(f_find)
                        #--------------------------------------------------
                        # add in nskip values
                        # when t_diff is positive then skip values in 
                        # the remote reference file
                        n_skip = abs(rr_arr['npts']-s_fn_arr['npts'].min())
                        if t_diff > 0 and dt_arr is not None:
                            dt_arr['nskip'] = n_skip
                        elif t_diff < 0:
                            #need to test if nskip is already there
                            if s_fn_birrp_arr['nskip'][0] != 1:
                                if n_skip > s_fn_birrp_arr['nskip'][0]:
                                    s_fn_birrp_arr['nskip'][:] = n_skip
                                else:
                                    pass

                            else:
                                s_fn_birrp_arr['nskip'][:] = n_skip
                          
                        # if there was a remote referenc channel found 
                        # append it to the array
                        if dt_arr is not None:
                            rr_birrp_fn_arr = np.append(rr_birrp_fn_arr,
                                                        dt_arr)
                                                        
                    # need to fill in what number remote references are
                    rr_index = 1
                    rr_count = 0
                    for rr_b_arr in rr_birrp_fn_arr:
                        # for now have a default calibration file
                        # nearly all the calibration files are identical
                        # this will work fine for now, till I can figure 
                        # out how to fill it in automatically, need to 
                        # change the format and metadate of mtpy TS
                        # for now get it from the calibration file
                        try:
                            rr_b_arr['calibration_fn'] = getattr(self.survey_config,
                                                                 'rr_{0}_{1:02}_cal_fn'.format(rr_b_arr['comp'], 
                                                                 rr_index-1))
                        except AttributeError:
                            print('Could not find calibration for {0}'.format(rr_b_arr['fn']))
                            
                        rr_b_arr['rr_num'] = rr_index
                        rr_count += 1
                        if rr_count%2 == 0 and rr_count != 0:
                            rr_index += 1
                        
                        # need to make sure there is n_skip for other
                        # remote references
                        rr_min_read = rr_birrp_fn_arr['nread'].min()
                        for rr_b_arr in rr_birrp_fn_arr:
                            rr_n_skip = abs(rr_b_arr['nread']-rr_min_read)
                            if rr_b_arr['nread'] != rr_min_read:
                                rr_b_arr['nskip'] = rr_n_skip
                                print(rr_b_arr['fn'], rr_n_skip)
                        
                    # append the remote reference data to station data                                                                        
                    s_fn_birrp_arr = np.append(s_fn_birrp_arr, rr_birrp_fn_arr)
                
                # fill in the dictionary accordingly
                s_dict[df].append(s_fn_birrp_arr[np.nonzero(s_fn_birrp_arr['nread'])])
        
        # need to check for maximum number of points
                
        
        # return the station dictionary        
        return s_dict
        
    def _fill_birrp_fn_arr(self, fn_arr, remote=False):
        """
        fill a birrp_fn_arr with a ts_fn_arr
        """
        # test whether fn_arr has just one entry
        if fn_arr.shape == ():
            nc = 1
        # if there are more than 1 entry get the set
        else:
            nc = fn_arr.shape[0]
        
        # make an empty array to fill
        s_fn_birrp_arr = np.zeros(nc, dtype=self._birrp_fn_dtype)
                                  
        # fill array with data
        s_fn_birrp_arr['fn'][:] = fn_arr['fn']
        s_fn_birrp_arr['nread'][:] = fn_arr['npts'].min()
        s_fn_birrp_arr['nskip'][:] = 22
        s_fn_birrp_arr['start_dt'][:] = fn_arr['start_dt']
        s_fn_birrp_arr['comp'][:] = fn_arr['comp']
        
        # be sure to fill in calibration file for station mags
        for sfb_arr in s_fn_birrp_arr:
            if sfb_arr['comp'] in ['hx', 'hy', 'hz'] and remote == False:
                try:
                    sfb_arr['calibration_fn'] = getattr(self.survey_config, 
                                                        '{0}_cal_fn'.format(sfb_arr['comp']))
                except AttributeError:
                    sfb_arr['calibration_fn'] = None
        
        if remote == True:
            s_fn_birrp_arr['rr'] = True
        
    
        # compute start and end times
        try:
            start_dt = fn_arr[0]['start_dt']
        except TypeError:
            start_dt = fn_arr['start_dt']
            
        start_seconds = time.mktime(time.strptime(start_dt, 
                                                  datetime_fmt))
 
        end_seconds = start_seconds+fn_arr['npts'].min()/float(fn_arr['df'].mean())

        s_fn_birrp_arr['end_dt'][:] = time.strftime(datetime_fmt,
                                        time.localtime(end_seconds))
                                        
        return s_fn_birrp_arr
        
    def write_script_files(self, birrp_arr_dict, save_path=None, 
                           birrp_params_dict={}, **kwargs):
        """
        write script files in the new method
        """
        
        # make save path
        if save_path is None:
            save_path = os.path.join(self.station_dir, 'BF')
        if not os.path.exists(save_path):
            os.mkdir(save_path)  
            
        script_fn_list = []            
        # loop through by keys, which should be sampling rates
        df_keys = list(birrp_arr_dict.keys())
        for df_key in df_keys:
            # make a path unique to the sampling rate
            bf_path = os.path.join(save_path, '{0:.0f}'.format(df_key))
            if not os.path.exists(bf_path):
                os.mkdir(bf_path)
            
            # get the fn_array, and make sure that it is a ndarray type
            birrp_fn_arr = np.array(birrp_arr_dict[df_key])
            
            # get station name
            ex_find = np.where(birrp_fn_arr[0]['comp']) == 'ex'
            station = os.path.splitext(os.path.basename(birrp_fn_arr[0][ex_find]['fn']))[0]
            station = station.split('_')[0]
            
            # add parameters to birrp_params_dict 
            birrp_params_dict['deltat'] = -1*df_key
            birrp_params_dict['ofil'] = os.path.join(bf_path, station)
            if df_key == 16:
                birrp_params_dict['nfft'] = 2**16
                birrp_params_dict['nsctmax'] = 11
            
            # make a script object passing on the desired birrp parameters
            birrp_script_obj = birrp.ScriptFile(**birrp_params_dict)
            birrp_script_obj.fn_arr = birrp_fn_arr
            # write script file
            birrp_script_obj.write_script_file(script_fn=os.path.join(bf_path,
                                                                      station+'.script'))
            # write a birrp parameter configuration file
            birrp_script_obj.write_config_file(os.path.join(bf_path, 
                                                            station))            
            script_fn_list.append(birrp_script_obj.script_fn)
        
        return script_fn_list
        
    def write_script_files_old(self, fn_birrp_dict, save_path=None,
                               birrp_params_dict={}, **kwargs):
        """
        write a script file from a generic processing dictionary
        """
        
        if save_path is None:
            save_path = os.path.join(self.station_dir, 'BF')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
        s_keys = list(fn_birrp_dict.keys())
        script_fn_list = []
        for df in s_keys:
            bf_path = os.path.join(save_path, '{0:.0f}'.format(df))
            fn_birrp_arr = fn_birrp_dict[df]
            pro_obj = BIRRP_processing()
            if df == 16:
                pro_obj.nfft = 2**16
                pro_obj.nsctmax = 11
            pro_obj.calibration_path = self.coil_cal_path
            pro_obj.station = self.survey_config.station
            pro_obj.deltat = -float(df)
            
            # for advanced processing
            if self.rr_station_dir is not None:
                pro_obj.ilev = 1
                pro_obj.tbw = 3
                pro_obj.nar = 9
                pro_obj.c2threshb = .35
                pro_obj.c2threshe = .45
                pro_obj.c2threshe1 = .45
                pro_obj.nf1 = 4
                pro_obj.nfinc = 2
                pro_obj.nsctinc = 2
                pro_obj.nfsect = 2
                pro_obj.ainlin = .0001
            for b_key in list(birrp_params_dict.keys()):
                setattr(pro_obj, b_key, birrp_params_dict[b_key])
            
            pro_dict = pro_obj.get_processing_dict(fn_birrp_arr, 
                                                   hx=self.survey_config.hx,
                                                   hy=self.survey_config.hy,
                                                   hz=self.survey_config.hz,
                                                   **kwargs)
                                                   
            
            #write script file using mtpy.processing.birrp    
            script_fn, birrp_dict = birrp.write_script_file(pro_dict,
                                                            save_path=bf_path)
                                                        
            script_fn_list.append(script_fn)
            
            cfg_fn = mtfh.make_unique_filename('{0}_birrp_params.cfg'.format(
                                                                 script_fn[:-7]))
                                                                 
            mtcfg.write_dict_to_configfile(birrp_dict, cfg_fn)
            print('Wrote BIRRP config file for edi file to {0}'.format(cfg_fn))
    
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
                        
        print(self.survey_config_fn)
        
        j2edi_obj = birrp.J_To_Edi(station=self.survey_config.station,
                                   survey_config_fn=self.survey_config_fn,
                                   birrp_dir=birrp_output_path,
                                   birrp_config_fn=self.birrp_config_fn)
                    
        edi_fn = j2edi_obj.write_edi_file()

        
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
 
    def process_data(self, df_list=None, max_blocks=3,
                     notch_dict={}, 
                     sr_dict={4096:(1000., 4),
                              1024:(3.99, 1.),
                              256:(3.99, .126), 
                              16:(.125, .0001)},
                     use_blocks_dict={4096:'all', 256:'all', 16:'all'},
                     birrp_param_dict={},
                     **kwargs):
        """
        process_data is a convinience function that will process Z3D files
        and output an .edi file.  The workflow is to convert Z3D files to 
        MTpy ascii format, apply a notch filter if desired, decimate if
        desired.  Do the same for a remote reference if specified.  Make a 
        BIRRP script file for each sampling rate, process data, take outputs 
        and write an edi file and plot resutl.  The whole process will take 
        a few minutes so be patient.  
        
        Arguments
        ---------------
        
            **df_list** : list of sampling rates
                          list of sampling rates (int) to process, options are
                          - 4096
                          - 1024
                          - 256
                          - 16
                          
            **max_blocks** : int
                            maximum number of blocks to process, this cannot
                            exceed the number of blocks that BIRRP was 
                            compiled with.  *default* is 3
                            
            **notch_dict** : dict(sampling_rate: notches to filter)
                            dictionary of notches to filter for each sampling
                            rate.  Keys are sampling rates (int) and values
                            are dictionary of notches to filter
                            ex. {4096:{60, 120, 180}} to filter out 60 Hz noise
                            *default* is {} which filters out 60 Hz noise and
                            harmonics.
                            
            **sr_dict** : dict(sampling_rate: (max_freq, max_freq))
                          dictionary of min and max frequencies to use for 
                          each sampling rate when making an .edi file.  The 
                          defaults usually work well, but check the plot
                          to see if there is a better frequency range for
                          each sampling rate.  
                          *default* is 
                              {4096:(1000., 4),
                              1024:(3.99, 1.),
                              256:(3.99, .126), 
                              16:(.125, .0001)}
                              
            **use_blocks_dict** : dict(sampling_rate: list of blocks to use)
                                  dictionary with sampling rate as keys and 
                                  list of which blocks to use as values
                                  ex. {4096: [0, 2]} to skip the second block
                                  use value 'all' to use all blocks
                                  *default* is 
                                      {4096:'all', 256:'all', 16:'all'}
                                      
            **birrp_param_dict** : dict(birrp_param: value)
                                   dictionary of birrp parameters where the 
                                   key is the birrp parameter and value is
                                   what you want that parameter to be.  See
                                   mtpy.processing.birrp.write_script_file
                                   or the BIRRP manual for details on birrp 
                                   parameters
                                   ex. {'nar':5}
                                   *default* is {}

        Returns
        -------------

            **plot_obj** : plot_response object

            **comb_edi_fn** : string
                              full path to combined edi file                          
            
        """
        if df_list is not None:
            self.df_list = df_list
            
        if max_blocks is not None:
            self.max_blocks = max_blocks
        
        # get start time
        st = time.time()
        
        if self.df_list is not None:
            if type(self.df_list) is float or type(self.df_list) is int or\
               type(self.df_list) is str:
               self.df_list = [self.df_list] 

        
        # make files into mtpy files
        z3d_fn_list, rr_fn_list, log_lines = self.make_mtpy_ascii_files(notch_dict=notch_dict)
        
        # get all information from mtpy files
        schedule_dict = self.get_schedules_fn_from_dir(use_blocks_dict=use_blocks_dict)
            
        # write script files for birrp
        sfn_list = self.write_script_files(schedule_dict,
                                           birrp_params_dict=birrp_param_dict,
                                           **kwargs)
        
        # run birrp
        self.run_birrp(sfn_list)
        
        # combine edi files
        comb_edi_fn = self.combine_edi_files(self.edi_fn, sr_dict)
        
        if comb_edi_fn is not None:
            self.edi_fn.append(comb_edi_fn)
        
        # plot the output
        r_plot = self.plot_responses()
        
        et = time.time()
        print('--> Processing took {0:02.0f}:{1:02.0f} minutes'.format((et-st)//60, (et-st)%60))
        
        return r_plot, comb_edi_fn
        
    def combine_edi_files(self, edi_fn_list, 
                          sr_dict={4096:(1000., 4),
                                   1024:(3.99, 1.),
                                   256:(3.99, .126), 
                                   16:(.125, .0001)}):
        """
        combine the different edi files that are computed for each sampling 
        rate.
        
        for now just a simple cutoff
        """
        
        if type(edi_fn_list) is str:
            print('Only one edi file, skipping combining')
            return
            
        if type(edi_fn_list) is list and len(edi_fn_list) == 1:
            print('Only one edi file, skipping combining')
            return
            
        data_arr = np.zeros(100, 
                            dtype=[('freq', np.float),
                                   ('z', (np.complex, (2, 2))),
                                   ('z_err', (np.float, (2, 2))), 
                                   ('tipper', (np.complex, (2, 2))),
                                   ('tipper_err', (np.float, (2, 2)))])
                                   
        count = 0
        for edi_fn in edi_fn_list:
            # get sampling rate from file name 
            fn_list = edi_fn[edi_fn.find('BF'):].split(os.path.sep)
            for ss in fn_list:
                try:
                    sr_key = int(ss)
                    break
                except ValueError:
                    pass
            if sr_key in list(sr_dict.keys()):
                try:
                    edi_obj = mtedi.Edi(edi_fn)
                    # locate frequency range
                    f_index = np.where((edi_obj.Z.freq >= sr_dict[sr_key][1]) & 
                                       (edi_obj.Z.freq <= sr_dict[sr_key][0]))
                                       
                                       
                    data_arr['freq'][count:count+len(f_index[0])] = edi_obj.Z.freq[f_index]
                    data_arr['z'][count:count+len(f_index[0])] = edi_obj.Z.z[f_index]
                    data_arr['z_err'][count:count+len(f_index[0])] = edi_obj.Z.z_err[f_index]
                    if edi_obj.Tipper.tipper is not None:                    
                        data_arr['tipper'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper[f_index]
                        data_arr['tipper_err'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper_err[f_index]
        
                    count += len(f_index[0])
                except IndexError:
                    print('Something went wrong with processing {0}'.format(edi_fn))
                
            else:
                print('{0} was not in combining dictionary'.format(sr_key))
                
        # now replace
        data_arr = data_arr[np.nonzero(data_arr['freq'])]
        sort_index = np.argsort(data_arr['freq'])
        
        # check to see if the sorted indexes are descending or ascending,
        # make sure that frequency is descending
        if data_arr['freq'][0] > data_arr['freq'][1]:
            sort_index = sort_index[::-1]
            
        data_arr = data_arr[sort_index]
        new_z = mtedi.MTz.Z(data_arr['z'],
                            data_arr['z_err'],
                            data_arr['freq'])
        
        # check for all zeros in tipper, meaning there is only 
        # one unique value                    
        if np.unique(data_arr['tipper']).size > 1:
            new_t = mtedi.MTz.Tipper(data_arr['tipper'], 
                                     data_arr['tipper_err'],
                                     data_arr['freq'])
                         
        else:
            new_t = mtedi.MTz.Tipper()
            
        edi_obj = mtedi.Edi(edi_fn_list[0])
        edi_obj.Z = new_z
        edi_obj.Tipper = new_t
        edi_obj.Data_sect.nfreq = new_z.z.shape[0]

        n_edi_fn = os.path.join(self.station_dir, 
                                '{0}_comb.edi'.format(os.path.basename(self.station_dir)))        
        n_edi_fn = edi_obj.write_edi_file(new_edi_fn=n_edi_fn)
        
        return n_edi_fn
                
def get_remote_reference_schedule(survey_path, plot=True):
    """
    Get a detailed list of which stations recorded at the same time, just from
    the raw Z3D files.
    
    Arguments
    ----------------
        **survey_path** : string
                          full path to the survey directory
                          
        **plot** : [ True | False]
                   True to plot the schedules as line bars. 
                   *default* is True
                   
    Outputs
    -----------------
        **station_path/Remote_Reference_List.txt** : file with a list of
                                                     dates with station
                                                     and sampling rate
                                                     
        **plot** if plot is True.
        
    :Example: ..
    
        >>> import mtpy.usgs.zen_processing as zp
        >>> zp.get_remote_reference_schedule(r"\home\MT_Data\Survey_01")
    """
    date_dict = {}
    rr_dict = {}
    for station in os.listdir(survey_path):
        station_path = os.path.join(survey_path, station)
        if os.path.isdir(station_path) is True:
            fn_list = [os.path.join(station_path, fn) 
                       for fn in os.listdir(station_path)
                       if fn.lower().endswith('ex.z3d')]
            
            if len(fn_list) == 0:
                print('No Z3D files found in folder: {0} '.format(station))
                continue
    
            station_date_arr = np.zeros(len(fn_list), dtype=[('df', np.float),
                                                             ('start_dt', '|S20')])
            for f_index, fn in enumerate(fn_list):
                zd = zen.Zen3D(fn)
                zd.read_all_info()
                station_date_arr[f_index]['df'] = zd.df
                station_date_arr[f_index]['start_dt'] = zd.zen_schedule
                try:
                    rr_dict[zd.zen_schedule]    
                except KeyError:
                    rr_dict[zd.zen_schedule] = []
                    
                rr_dict[zd.zen_schedule].append((station, zd.df))
                    
            if len(np.nonzero(station_date_arr)[0]) != 0:
                date_dict[station] = station_date_arr[np.nonzero(station_date_arr['df'])]
    
    #---------------------------------------------
    # print out in a useful way
    
    lines = []
    for key in sorted(rr_dict.keys()):
        k_list = key.split(',')    
        k_date = k_list[0]
        k_time = k_list[1]
        lines.append('Date: {0}, Time: {1}'.format(k_date, k_time))
        lines.append('-'*60)
        for k_tuple in rr_dict[key]:
            lines.append('\tStation: {0}, Sampling Rate: {1:.0f}'.format(k_tuple[0],
                         k_tuple[1]))
                         
        lines.append('='*60)
    with open(os.path.join(survey_path, 'Remote_Reference_List.txt'), 'w') as fid:
        fid.write('\n'.join(lines))
    
    
    
    #-------------------------------------
    
    if plot is True:
        plt.rcParams['font.size'] = 14
        
        datetime_display = '%m-%d, %H:%M:%S'
        
        df_dict = {4096:(.7, .1, 0), 1024:(.5, .5, 0), 256:(0, .2, .8)}                       
        # plot the results is a compeling graph
        fig = plt.figure(2, [12, 10])
        ax = fig.add_subplot(1, 1, 1)
        y_labels = ['', '']
        
        for k_index, key in enumerate(sorted(date_dict.keys())):
            y_labels.append(key)
            for ii, k_arr in enumerate(date_dict[key]):
                x_date_0 = datetime.datetime.strptime(k_arr['start_dt'], datetime_fmt)
                try:
                    x_date_1 = datetime.datetime.strptime(date_dict[key]['start_dt'][ii+1],
                                                          datetime_fmt)
                except IndexError:
                    x_date_1 = x_date_0
                
                y_values = [k_index, k_index]
                
                l1, = ax.plot([x_date_0, x_date_1], y_values,
                              lw=8, color=df_dict[k_arr['df']])
        
        fig.autofmt_xdate(rotation=60)
        ax.xaxis.set_major_formatter(mdates.DateFormatter(datetime_display))
#        ax.xaxis.set_major_locator(MultipleLocator((1)))
#        ax.xaxis.set_minor_locator(MultipleLocator((.25)))
        ax.xaxis.set_tick_params(width=2, size=5)
       
        ax.yaxis.set_major_locator(MultipleLocator(1))
        ax.yaxis.set_ticklabels(y_labels)
        ax.yaxis.set_tick_params(width=2, size=5)
        ax.set_ylim(-1, len(y_labels)-2)
        
        ax.grid(which='major', linestyle='--', color=(.7, .7, .7)) 
        ax.set_axisbelow(True)
        
        l_4096 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[4096])
        l_1024 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[1024])
        l_256 = plt.Line2D([0, 1], [0, 0], lw=8, color=df_dict[256])
      
        fig.legend([l_4096, l_1024, l_256], ['4096', '1024', '256'], ncol=3,
                   loc='upper center')
                             
        plt.show()   

#==============================================================================

# this should capture all the print statements from zen
class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
#        sys.stdout = self._stringio = StringIO()
        sys.stdout = self._stringio = BytesIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout


#==============================================================================
def compute_mt_response(survey_dir, station='mt000', copy_date=None, 
                        birrp_exe=r"c:\MinGW32-xy\Peacock\birrp52\birrp52_3pcs6e9pts.exe", 
                        ant_calibrations=r"c:\MT\Ant_calibrations",
                        process_df_list=[256],
                        max_blocks=2,
                        notch_dict={256:None}):
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
            zen.copy_from_sd(station, save_path=survey_dir)
        else:
            zen.copy_from_sd(station, save_path=survey_dir, 
                             copy_date=copy_date, copy_type='after')
    except IOError:
        print('No files copied from SD cards')
        print('Looking in  {0} for Z3D files'.format(station_dir))
    
    #--> process data
     
    with Capturing() as output:
        z2edi = Z3D_to_edi(station_dir)
        z2edi.birrp_exe = birrp_exe
        z2edi.coil_cal_path = ant_calibrations
        try:
            rp, comb_edi = z2edi.process_data(df_list=process_df_list,  
                                              max_blocks=max_blocks,
                                              notch_dict=notch_dict)
        except mtex.MTpyError_inputarguments:
            print('==> Data not good!! Did not produce a proper .edi file') 
            et = time.time()
            print('--> took {0} seconds'.format(et-st))
            rp = None
    
    #--> write log file
    log_fid = open(os.path.join(station_dir, 'Processing.log'), 'w')
    log_fid.write('\n'.join(output))
    log_fid.close()
    
    return rp         
        
        