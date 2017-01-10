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

import mtpy.utils.filehandling as mtfh
import mtpy.processing.new_birrp as birrp
import mtpy.utils.configfile as mtcfg
import mtpy.utils.exceptions as mtex
import mtpy.imaging.plotnresponses as plotnresponses
import mtpy.imaging.plotresponse as plotresponse
import mtpy.usgs.zen as zen
import mtpy.core.edi as mtedi

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
        
        self.fn_list = []
        self.rrfn_list = []
        # need to sort the fn list so that the files are in the correct
        # order for input and output as defined by birrp
        self.nread = []
        self.nskip = []
        self.nskipr = []
        self.nread = []
        for block_arr in fn_birrp_list:
            s_list = np.zeros(len(block_arr), dtype='|S100')
            r_list = np.zeros(2, dtype='|S100')
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
                    self.nread.append(nread_ii-3)
                else:
                    self.nskip.append(1+int(nskip))
                    self.nskipr.append(1+int(nskipr))
                    self.nread.append(nread_ii)
    
                # append file names appropriately removing any empty names
                self.fn_list.append(list(s_list[np.where(s_list != '')]))
                if len(r_list[np.where(r_list != '')]) > 0:
                    self.rrfn_list.append(list(r_list[np.where(r_list != '')]))
                    
            else:
                print('Not enough points {0}'.format(nread_ii))
                print('skipping time block {0} for sampling rate {1}'.format(
                      time.strftime(datetime_fmt, time.localtime(start_dt)),
                      -self.deltat))

        # need to check for the number of points to be read in, there is 
        # a memory max, for this computer the max is self._max_nread        
        if sum(self.nread) > self._max_nread:
            self.nread[-1] = self._max_nread-sum(self.nread[0:-1])
            print 'processing {0} points'.format(sum(self.nread))
        
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
                        
    def read_config_file(self, birrp_config_fn):
        """
        read in a configuration file and fill in the appropriate parameters
        """
        
        birrp_dict = mtcfg.read_configfile(birrp_config_fn)
        
        for birrp_key in birrp_dict.keys():
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
        self.b_instrument_type = 'coil'
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
        self.coil_cal_path = r"d:\Peacock\MTData\Ant_calibrations"
        self.num_comp = 5
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])

        
    def make_survey_config_file(self, survey_config_dict=None):
        """
        make a survey configuration file from the data
        """
        
        self.survey_config_fn = self.survey_config.write_survey_config_file()
    
    def get_z3d_fn_blocks(self, station_dir=None, remote=False,
                          df_list=[4096, 1024, 256], max_blocks=3):
        """
        get z3d file names in an array of blocks
        """

        fn_block_dict = dict([(df, {}) for df in df_list])
        fn_count = 0
        for fn in os.listdir(station_dir):
            if fn.lower().endswith('.z3d'):
            
                z3d_fn = os.path.join(station_dir, fn)
                z3d_obj = zen.Zen3D(z3d_fn)
                z3d_obj.read_all_info()
                if remote is True:
                    if not z3d_obj.metadata.ch_cmp.lower() in ['hx', 'hy']:
                        continue
                if z3d_obj.df in df_list:
                    fn_count += 1
                    try:
                        fn_block_dict[z3d_obj.df][z3d_obj.zen_schedule].append(z3d_fn)
                    except KeyError:
                        fn_block_dict[z3d_obj.df][z3d_obj.zen_schedule] = [z3d_fn]


        if fn_count == 0:
            raise ValueError('No Z3D files found for in {0}'.format(station_dir))
        else:
            print 'Found {0} Z3D files in {1}'.format(fn_count, station_dir)
        
        # check for maximum number of blocks
        for df_key in fn_block_dict.keys():
            date_dict = fn_block_dict[df_key]
            dates = sorted(date_dict.keys())
            if len(dates) == 0:
                print 'No Z3D files found for {0} in {1}'.format(str(df_key),
                                                                station_dir)

            if len(dates) > max_blocks:
                for pop_date in dates[-(len(dates)-max_blocks):]:
                    fn_block_dict[df_key].pop(pop_date)
        
        return fn_block_dict
        
    def make_mtpy_ascii_files(self, station_dir=None, rr_station_dir=None, 
                              fmt='%.8', station_name='mb', 
                              notch_dict={4096:{}, 256:None, 1024:None}, 
                              df_list=[4096, 1024, 256], max_blocks=3,
                              ex=100., ey=100.,): 
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
            
        fn_dict = self.get_z3d_fn_blocks(self.station_dir,
                                         df_list=df_list,
                                         max_blocks=max_blocks)
        
        num_files = 0
        num_comp = []
        for df_key in fn_dict.keys(): 
            for date_key in fn_dict[df_key].keys():
                # in case there is decimation, then we need to double the
                # number of 256 files
                if df_key == 256 and 16 in df_list:
                    num_files += len(fn_dict[df_key][date_key])*2
                else:
                    num_files += len(fn_dict[df_key][date_key])
                num_comp.append(len(fn_dict[df_key][date_key]))
                
        self.num_comp = int(np.mean(np.array(num_comp)))
        print ' --> processing {0} components'.format(self.num_comp)
                    
        # make an array that has all the information about each file
        fn_arr = np.zeros(num_files, 
                          dtype=[('station','|S6'), 
                                 ('npts',np.int), 
                                 ('df', np.int),
                                 ('start_dt', '|S22'), 
                                 ('comp','|S2'),
                                 ('fn','|S100')])
        fn_lines = []
        jj = 0  # index for fn_array            
        for df_key in df_list: 
            for date_key in fn_dict[df_key].keys():
                for fn in fn_dict[df_key][date_key]:
                    zd = zen.Zen3D(fn)
                    zd.read_all_info()
               
                    # account for decimation, need to make a new instance
                    # so as to include the original file as well
                    if zd.header.ad_rate == 256 and 16 in df_list:
                        zdec = zen.Zen3D(fn)
                        zdec.read_all_info()

                        ex = float(zdec.metadata.ch_length)
                        ey = float(zdec.metadata.ch_length)
                        #write mtpy mt file
                        zdec.write_ascii_mt_file(notch_dict=notch_dict[df_key], 
                                                 ex=ex, ey=ey, dec=16)
                        
                        #create lines to write to a log file                       
                        station = zdec.metadata.rx_xyz0.split(':')[0]
                        fn_arr[jj]['station'] = '{0}{1}'.format(station_name, 
                                                                station)
                        fn_arr[jj]['npts'] = zdec.time_series_len
                        fn_arr[jj]['df'] = zdec.df
                        fn_arr[jj]['start_dt'] = zdec.zen_schedule
                        fn_arr[jj]['comp'] = zdec.metadata.ch_cmp.lower()
                        fn_arr[jj]['fn'] = zdec.fn_mt_ascii
                        fn_lines.append(''.join(['--> station: {0}{1}\n'.format(station_name, station),
                                                 '    ts_len = {0}\n'.format(zdec.time_series_len),
                                                 '    df = {0}\n'.format(zdec.df),
                                                 '    start_dt = {0}\n'.format(zdec.zen_schedule),
                                                 '    comp = {0}\n'.format(zdec.metadata.ch_cmp),
                                                 '    fn = {0}\n'.format(zdec.fn),
                                                 '    decimated to 16 samples/s']))
                   
        
                        jj += 1
                
                
                
                    # get metadata information about the channel to put into 
                    # survey configuration file
                    if zd.metadata.ch_cmp.lower() == 'hx':
                        self.survey_config.hx = zd.metadata.ch_number
                    elif zd.metadata.ch_cmp.lower() == 'hy':
                        self.survey_config.hy = zd.metadata.ch_number
                    elif zd.metadata.ch_cmp.lower() == 'hz':
                        self.survey_config.hz = zd.metadata.ch_number
                    elif zd.metadata.ch_cmp.lower() == 'ex':
                        self.survey_config.e_xaxis_length = zd.metadata.ch_length
                        ex = float(zd.metadata.ch_length)
                    elif zd.metadata.ch_cmp.lower() == 'ey':
                        self.survey_config.e_yaxis_length = zd.metadata.ch_length
                        ey = float(zd.metadata.ch_length)
                        
                    #write mtpy mt file
                    zd.write_ascii_mt_file(notch_dict=notch_dict[df_key],
                                           ex=ex, ey=ey)
        
                    # get station configuration from the first Z3D file            
                    if jj == 1:
                        self.survey_config.lat = zd.header.lat
                        self.survey_config.lon = zd.header.long
                        self.survey_config.date = zd.schedule.Date.replace('-','/')
                        self.survey_config.box = int(zd.header.box_number)
                        self.survey_config.station = zd.metadata.rx_xyz0.split(':')[0]
     
                    #create lines to write to a log file                       
                    station = zd.metadata.rx_xyz0.split(':')[0]
                    fn_arr[jj]['station'] = '{0}{1}'.format(station_name, station)
                    fn_arr[jj]['npts'] = zd.time_series_len
                    fn_arr[jj]['df'] = zd.df
                    fn_arr[jj]['start_dt'] = zd.zen_schedule
                    fn_arr[jj]['comp'] = zd.metadata.ch_cmp.lower()
                    fn_arr[jj]['fn'] = zd.fn_mt_ascii
                    fn_lines.append(''.join(['--> station: {0}{1}\n'.format(station_name, station),
                                             '    ts_len = {0}\n'.format(zd.time_series_len),
                                             '    df = {0}\n'.format(zd.df),
                                             '    start_dt = {0}\n'.format(zd.zen_schedule),
                                             '    comp = {0}\n'.format(zd.metadata.ch_cmp),
                                             '    fn = {0}\n'.format(zd.fn)]))
                    jj += 1               
        #-------------------------------
        # Remote Reference
        #-------------------------------          
        if rr_station_dir is not None:
            self.rr_station_dir = rr_station_dir
            
        if self.rr_station_dir is not None:
            print 'Tranforming remote reference Z3D to mtpy format'
            rr_fn_dict = self.get_z3d_fn_blocks(self.rr_station_dir,
                                                df_list=df_list,
                                                max_blocks=max_blocks,
                                                remote=True)
            
            num_files = 0
            for df_key in rr_fn_dict.keys(): 
                for date_key in rr_fn_dict[df_key].keys():
                    # in case there is decimation, then we need to double the
                    # number of 256 files
                    if df_key == 256 and 16 in df_list:
                        num_files += len(rr_fn_dict[df_key][date_key])*2
                    else:
                        num_files += len(rr_fn_dict[df_key][date_key])
       
            # make an array that has all the information about each file
            rr_fn_arr = np.zeros(num_files, 
                                  dtype=[('station','|S6'), 
                                         ('npts',np.int), 
                                         ('df', np.int),
                                         ('start_dt', '|S22'), 
                                         ('comp','|S2'),
                                         ('fn','|S100')])

            jj = 0  # index for fn_array            
            for df_key in df_list: 
                for date_key in rr_fn_dict[df_key].keys():
                    for fn in rr_fn_dict[df_key][date_key]:
                        zd = zen.Zen3D(fn)
                        zd.read_all_info()
                    
                        # account for decimation, need to make a new instance
                        # so as to include the original file as well
                        if zd.header.ad_rate == 256 and 16 in df_list:
                            zdec = zen.Zen3D(fn)
                            zdec.read_all_info()
                            
                            #write mtpy mt file
                            zdec.write_ascii_mt_file(notch_dict=notch_dict[df_key],
                                                     dec=16)
                            
                            #create lines to write to a log file                       
                            rr_station = zdec.metadata.rx_xyz0.split(':')[0]
                            rr_fn_arr[jj]['station'] = '{0}{1}'.format(station_name, 
                                                                       rr_station)
                            rr_fn_arr[jj]['npts'] = zdec.time_series_len
                            rr_fn_arr[jj]['df'] = zdec.df
                            rr_fn_arr[jj]['start_dt'] = zdec.zen_schedule
                            rr_fn_arr[jj]['comp'] = 'rr{0}'.format(zdec.metadata.ch_cmp.lower())
                            rr_fn_arr[jj]['fn'] = zdec.fn_mt_ascii
                            fn_lines.append(''.join(['--> rr_station: {0}{1}\n'.format(station_name, rr_station),
                                                     '    ts_len = {0}\n'.format(zdec.time_series_len),
                                                     '    df = {0}\n'.format(zdec.df),
                                                     '    start_dt = {0}\n'.format(zdec.zen_schedule),
                                                     '    comp = rr{0}\n'.format(zdec.metadata.ch_cmp),
                                                     '    fn = {0}\n'.format(zdec.fn),
                                                     '    decimated to 16 samples/s']))
            
                            jj += 1
            
                    
                        #write mtpy mt file
                        zd.write_ascii_mt_file(notch_dict=notch_dict[df_key])
                        
                        #create lines to write to a log file                       
                        rr_station = zd.metadata.rx_xyz0.split(':')[0]
                        rr_fn_arr[jj]['station'] = '{0}{1}'.format(station_name, rr_station)
                        rr_fn_arr[jj]['npts'] = zd.time_series_len
                        rr_fn_arr[jj]['df'] = zd.df
                        rr_fn_arr[jj]['start_dt'] = zd.zen_schedule
                        rr_fn_arr[jj]['comp'] = 'rr{0}'.format(zd.metadata.ch_cmp.lower())
                        rr_fn_arr[jj]['fn'] = zd.fn_mt_ascii
                        fn_lines.append(''.join(['--> station: {0}{1}\n'.format(station_name, rr_station),
                                                 '    ts_len = {0}\n'.format(zd.time_series_len),
                                                 '    df = {0}\n'.format(zd.df),
                                                 '    start_dt = {0}\n'.format(zd.zen_schedule),
                                                 '    comp = rr{0}\n'.format(zd.metadata.ch_cmp),
                                                 '    fn = {0}\n'.format(zd.fn)]))
                
                        # get metadata information about the channel to put into 
                        # survey configuration file
                        if zd.metadata.ch_cmp.lower() == 'hx':
                            self.survey_config.rr_hx = zd.metadata.ch_number
                        elif zd.metadata.ch_cmp.lower() == 'hy':
                            self.survey_config.rr_hy = zd.metadata.ch_number
                        
                        # input remote reference information into survey file
                        if jj == 1:
                            self.survey_config.rr_station = rr_station
                            self.survey_config.rr_latitude = zd.header.lat
                            self.survey_config.rr_longitude = zd.header.long
                            self.survey_config.rr_date = zd.schedule.Date.replace('-','/')
                            self.survey_config.rr_box = int(zd.header.box_number)
                
            # change the time series directory to reflect where the time series
            # are.
            self.rr_station_dir = os.path.join(self.rr_station_dir, 'TS')
            
            # get only the non-empty time series
            rr_fn_arr = rr_fn_arr[np.nonzero(rr_fn_arr['npts'])]
        else:
            rr_fn_arr = None
        # change directories to be associated with where the mtpy ts are
        self.station_dir = os.path.join(self.station_dir, 'TS')
        self.survey_config.save_path = self.station_dir
        
        # write survey configuration file
        self.survey_config.write_survey_config_file()
 
        return fn_arr[np.nonzero(fn_arr['npts'])], rr_fn_arr, fn_lines
        
    def get_schedules_fn_from_dir(self, station_ts_dir=None, rr_ts_dir=None,
                                  df_list=[4096, 256, 16], max_blocks=4,
                                  use_blocks_dict={4096:'all', 256:'all', 16:'all'}):
        """
        get the birrp fn list from a directory of TS files
        """
        
        if station_ts_dir is not None:
            self.station_dir = station_ts_dir
        
        if rr_ts_dir is not None:
            self.rr_station_dir = rr_ts_dir
            
        if not os.path.isdir(self.station_dir):
            print '{0} is not a valid directory, check path.'.format(self.station_dir) 
            return None
            
        fn_arr = np.zeros(len(os.listdir(self.station_dir)),
                          dtype=[('fn','|S100'),
                                 ('npts', np.int),
                                 ('start_dt', '|S19'),
                                 ('end_dt', '|S19'),
                                 ('df', np.float),
                                 ('comp', '|S4')])
        fn_count = 0
        for fn in os.listdir(self.station_dir):
            fn = os.path.join(self.station_dir, fn)
            try:
                header_dict = mtfh.read_ts_header(fn)
                if header_dict['samplingrate'] in df_list:
                    fn_arr[fn_count]['fn'] = fn
                    fn_arr[fn_count]['npts'] = header_dict['nsamples']
                    fn_arr[fn_count]['df'] = header_dict['samplingrate']
                    fn_arr[fn_count]['comp'] = header_dict['channel']
                    start_sec = header_dict['t_min']
                    num_sec = float(header_dict['nsamples'])/\
                                                       header_dict['samplingrate']
                    fn_arr[fn_count]['start_dt'] = time.strftime(datetime_fmt, 
                                                    time.localtime(start_sec)) 
                    fn_arr[fn_count]['end_dt'] = time.strftime(datetime_fmt, 
                                                    time.localtime(start_sec+\
                                                    num_sec))
                    fn_count += 1
                else:
                    pass
            except mtex.MTpyError_ts_data:
                print '  Skipped {0}'.format(fn)
            except mtex.MTpyError_inputarguments:
                print '  Skipped {0}'.format(fn)
        # be sure to trim the array
        fn_arr = fn_arr[np.nonzero(fn_arr['npts'])]   
        
        # get remote reference filenames
        if self.rr_station_dir is not None:
            if not os.path.isdir(self.rr_station_dir):
                print '{0} is not a valid directory, check path.'.format(self.rr_station_dir) 
                return fn_arr[np.nonzero(fn_arr['npts'])]
            rr_fn_arr = np.zeros(len(os.listdir(self.rr_station_dir)),
                                  dtype=[('fn','|S100'),
                                         ('npts', np.int),
                                         ('start_dt', '|S19'),
                                         ('end_dt', '|S19'),
                                         ('df', np.float),
                                         ('comp', '|S4')])  
            rr_fn_count = 0
            for rr_fn in os.listdir(self.rr_station_dir):
                rr_fn = os.path.join(self.rr_station_dir, rr_fn)
                try:
                    # only get the hx and hy channels                    
                    header_dict = mtfh.read_ts_header(rr_fn)                    
                    if header_dict['channel'].lower() in ['hx', 'hy'] and\
                       header_dict['samplingrate'] in df_list:
                        rr_fn_arr[rr_fn_count]['fn'] = rr_fn
                        rr_fn_arr[rr_fn_count]['npts'] = header_dict['nsamples']
                        rr_fn_arr[rr_fn_count]['df'] = header_dict['samplingrate']
                        rr_fn_arr[rr_fn_count]['comp'] = 'rr{0}'.format(header_dict['channel'])
                        start_sec = header_dict['t_min']
                        num_sec = float(header_dict['nsamples'])/\
                                        header_dict['samplingrate']
                        rr_fn_arr[rr_fn_count]['start_dt'] = time.strftime(datetime_fmt, 
                                                        time.localtime(start_sec)) 
                        rr_fn_arr[rr_fn_count]['end_dt'] = time.strftime(datetime_fmt, 
                                                        time.localtime(start_sec+\
                                                        num_sec))
                        rr_fn_count += 1
                except mtex.MTpyError_ts_data:
                    print '  Skipped {0}'.format(rr_fn)
                except mtex.MTpyError_inputarguments:
                    print '  Skipped {0}'.format(rr_fn) 
            rr_fn_arr = rr_fn_arr[np.nonzero(rr_fn_arr['npts'])]
        else:
            rr_fn_arr = None
        
        
        return self.get_schedules_fn(fn_arr, rr_fn_arr, df_list,
                                     max_blocks=max_blocks,
                                     use_blocks_dict=use_blocks_dict)
#        return (fn_arr, rr_fn_arr)
            
        
    def get_schedules_fn(self, fn_arr, rr_fn_arr=None, df_list=[4096, 256, 16],
                         max_blocks=4, 
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
        # get the sampling rates used
        s_keys = set(df_list)
        
        # make a dictionary with keys as the sampling rates 
        s_dict = dict([(skey, []) for skey in s_keys])
        
        # loop over the sampling rates and find the schedule blocks
        for df in s_keys:
            # find startind dates for sampling rate
            s_dates = sorted(list(set(fn_arr['start_dt'][np.where(fn_arr['df']==df)])))
            use_blocks = use_blocks_dict[df]            
            if use_blocks == 'all':
                date_list = s_dates[0:max_blocks]
            elif type(use_blocks) is list:
                date_list = [s_dates[ii] for ii in use_blocks]
                if len(date_list) > max_blocks:
                    date_list = date_list[0:max_blocks]
                
            else:
                raise ValueError('Do not understand use_blocks type {0}'.format(type(use_blocks)))
                
            for sdate in date_list:
                s_fn_arr = fn_arr[np.where((fn_arr['start_dt']==sdate) &
                                            (fn_arr['df']==df))]
                num_comp = len(s_fn_arr) 
                s_fn_birrp_arr = np.zeros(num_comp+2, 
                                          dtype=[('fn','|S100'),
                                                 ('npts',np.int),
                                                 ('start_dt','|S19'),
                                                 ('end_dt','|S19'),
                                                 ('comp', '|S4')])
                s_fn_birrp_arr['fn'][0:num_comp] = s_fn_arr['fn']
                s_fn_birrp_arr['npts'][0:num_comp] = s_fn_arr['npts'].min()
                s_fn_birrp_arr['start_dt'][0:num_comp] = sdate
                s_fn_birrp_arr['comp'][0:num_comp] = s_fn_arr['comp']
                start_seconds = time.mktime(time.strptime(sdate, 
                                                          datetime_fmt))
 
                end_seconds = start_seconds+s_fn_arr['npts'].min()/float(df)
                s_fn_birrp_arr['end_dt'][0:num_comp] = time.strftime(datetime_fmt,
                                                time.localtime(end_seconds))  
                
                
                # get remote reference information if input
                if rr_fn_arr is not None:
                    # need to find the date closest to the station date
                    if sdate in set(rr_fn_arr['start_dt']):
                        s_rr_arr = rr_fn_arr[np.where((rr_fn_arr['start_dt']==sdate) &
                                            (rr_fn_arr['df']==df))]
                        s_fn_birrp_arr['fn'][num_comp:] = s_rr_arr['fn']
                        s_fn_birrp_arr['npts'][num_comp:] = s_rr_arr['npts'].min()
                        s_fn_birrp_arr['start_dt'][num_comp:] = sdate
                        s_fn_birrp_arr['comp'][num_comp:] = s_rr_arr['comp']
                        start_seconds = time.mktime(time.strptime(sdate, 
                                                          datetime_fmt))
                        end_seconds = start_seconds+s_rr_arr['npts'].min()/float(df)
                        s_fn_birrp_arr['end_dt'][num_comp:] = time.strftime(datetime_fmt,
                                                        time.localtime(end_seconds))
                    # locate the closest time for remote referencing
                    else:
                        rr_df_arr = rr_fn_arr[np.where(rr_fn_arr['df']==df)]
                        rr_dates = sorted(list(set(rr_df_arr['start_dt'])))
                        t_diff = np.zeros(len(rr_dates))
                        # estimate the time difference between station date
                        # and remote reference date
                        for ii, rr_date in enumerate(rr_dates):
                            rr_sec = time.mktime(time.strptime(rr_date, 
                                                          datetime_fmt))
                            s_sec = time.mktime(time.strptime(sdate, 
                                                          datetime_fmt))
                                                          
                            # this might be dangerous, but taking the 
                            # absolute value of the time difference and 
                            # assuming the smallest one is the closest                                                          
                            t_diff[ii] = abs(s_sec-rr_sec)
                        
                        # find the index where the minimum is
                        try:
                            rr_date = rr_dates[np.where(t_diff==t_diff.min())[0]]
                            print 'Using remote reference time series starting on {0}'.format(rr_date)
                            print 'For station time series starting on            {0}'.format(sdate)
                            
                            s_rr_arr = rr_fn_arr[np.where((rr_fn_arr['start_dt']==rr_date) &
                                                (rr_fn_arr['df']==df))]
                            s_fn_birrp_arr['fn'][num_comp:] = s_rr_arr['fn']
                            s_fn_birrp_arr['npts'][num_comp:] = s_rr_arr['npts'].min()
                            s_fn_birrp_arr['start_dt'][num_comp:] = rr_date
                            s_fn_birrp_arr['comp'][num_comp:] = s_rr_arr['comp']
                            start_seconds = time.mktime(time.strptime(rr_date, 
                                                              datetime_fmt))
                            end_seconds = start_seconds+s_rr_arr['npts'].min()/float(df)
                            s_fn_birrp_arr['end_dt'][num_comp:] = time.strftime(datetime_fmt,
                                                            time.localtime(end_seconds))
                        except ValueError:
                            pass
                        
                        
                        
                
                s_dict[df].append(s_fn_birrp_arr[np.nonzero(s_fn_birrp_arr['npts'])])
        return s_dict
        
    
        
    def write_script_files(self, fn_birrp_dict, save_path=None,
                           birrp_params_dict={}):
        """
        write a script file from a generic processing dictionary
        """
        
        if save_path is None:
            save_path = os.path.join(self.station_dir, 'BF')
        if not os.path.exists(save_path):
            os.mkdir(save_path)
            
        s_keys = fn_birrp_dict.keys()
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
            for b_key in birrp_params_dict.keys():
                setattr(pro_obj, b_key, birrp_params_dict[b_key])
            
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
                     birrp_param_dict={}):
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
        
        st = time.time()
        
        if df_list is not None:
            if type(df_list) is float or type(df_list) is int or\
               type(df_list) is str:
               df_list = [df_list] 
        
        # make files into mtpy files
        z3d_fn_list, rr_fn_list, log_lines = self.make_mtpy_ascii_files(df_list=df_list,
                                                            max_blocks=max_blocks,
                                                            notch_dict=notch_dict)
        
        # get all information from mtpy files
        schedule_dict = self.get_schedules_fn_from_dir(df_list=df_list,
                                                       max_blocks=max_blocks,
                                                       use_blocks_dict=use_blocks_dict)
            
        # write script files for birrp
        sfn_list = self.write_script_files(schedule_dict,
                                           birrp_params_dict=birrp_param_dict)
        
        # run birrp
        self.run_birrp(sfn_list)
        
        # combine edi files
        comb_edi_fn = self.combine_edi_files(self.edi_fn, sr_dict)
        
        if comb_edi_fn is not None:
            self.edi_fn.append(comb_edi_fn)
        
        # plot the output
        r_plot = self.plot_responses()
        
        et = time.time()
        print '--> Processing took {0:02.0f}:{1:02.0f} minutes'.format((et-st)//60, (et-st)%60)
        
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
            print 'Only one edi file, skipping combining'
            return
            
        if type(edi_fn_list) is list and len(edi_fn_list) == 1:
            print 'Only one edi file, skipping combining'
            return
            
        data_arr = np.zeros(100, 
                            dtype=[('freq', np.float),
                                   ('z', (np.complex, (2, 2))),
                                   ('z_err', (np.float, (2, 2))), 
                                   ('tipper', (np.complex, (2, 2))),
                                   ('tipper_err', (np.float, (2, 2)))])
                                   
        count = 0
        for edi_fn in edi_fn_list:
            edi_obj = mtedi.Edi(edi_fn)
            # get sampling rate from directory path
            for fkey in sorted(sr_dict.keys(), reverse=True):
                if str(fkey) in edi_fn:
                    # locate frequency range
                    f_index = np.where((edi_obj.Z.freq >= sr_dict[fkey][1]) & 
                                       (edi_obj.Z.freq <= sr_dict[fkey][0]))
                                       
                                       
                    data_arr['freq'][count:count+len(f_index[0])] = edi_obj.Z.freq[f_index]
                    data_arr['z'][count:count+len(f_index[0])] = edi_obj.Z.z[f_index]
                    data_arr['z_err'][count:count+len(f_index[0])] = edi_obj.Z.z_err[f_index]
                    if edi_obj.Tipper.tipper is not None:                    
                        data_arr['tipper'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper[f_index]
                        data_arr['tipper_err'][count:count+len(f_index[0])] = edi_obj.Tipper.tipper_err[f_index]
        
                    count += len(f_index[0])
                    
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
                print 'No Z3D files found in folder: {0} '.format(station)
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
        ax.xaxis.set_major_locator(MultipleLocator((1)))
        ax.xaxis.set_minor_locator(MultipleLocator((.25)))
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
        
        