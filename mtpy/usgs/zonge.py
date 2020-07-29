# -*- coding: utf-8 -*-
"""
====================
zonge
====================
    * Tools for interfacing with MTFT24
    * Tools for interfacing with MTEdit
    
    
Created on Tue Jul 11 10:53:23 2013
@author: jpeacock-pr
"""

#==============================================================================

import numpy as np
import time
import os
import shutil

import mtpy.core.z as mtz
import mtpy.imaging.plotresponse as plotresponse
import mtpy.utils.gis_tools as gis_tools
import mtpy.utils.configfile as mtcf
import mtpy.core.edi as mtedi
import mtpy.usgs.zen as zen

#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'

#==============================================================================
# class for  mtft24
#==============================================================================
class ZongeMTFT():
    """
    Reads and writes config files for MTFT24 version 1.10
    
    The important thing to have is the survey configuration file.  This is a
    configuration file that has all the important information about the
    survey in it.  And entry in this file should look something like:
    
    [MB093]
    battery = Li
    b_xaxis_azimuth = 0
    b_yaxis_azimuth = 90
    b_instrument_amplification = 1.0
    b_instrument_type = coil
    b_logger_gain = 1.0
    b_logger_type = zen
    date = 2013-11-05
    data_logger = 25
    e_xaxis_azimuth = 0
    e_yaxis_azimuth = 90
    e_instrument_amplification = 1.0
    e_instrument_type = electrodes
    e_logger_gain = 1.0
    e_logger_type = zen
    e_xaxis_length = 107
    e_yaxis_length = 104
    declination = 0.0
    elevation = 2324
    latitude = 37.82125
    location = Mono Basin
    longitude = -119.04723
    hx = 2304
    hy = 2344
    hz = 2314
    network = USGS
    notes = near an intrusion, HX power off on pickup, HY loose connection
    sampling_interval = 0
    station = MB093
    station_type = MT
    
    
    ========================= =================================================
    Attributes                Description
    ========================= =================================================
     Ant_FrqMax               Frequency max of antenna
     Ant_FrqMin               Frequency min of antenna
     Chn_Cmp                  Channel components (string)
     Chn_Cmp_lst              List of channel components
     Chn_Gain                 Channel gain
     Chn_ID                   Channel ID (number))
     Chn_Length               length of dipole for that channel
     Chn_dict                 dictionary of channel information
     MTFT_BandFrq             Frequency band to process
     MTFT_BandFrqMax          frequency band max to process
     MTFT_BandFrqMin          frequency band min to process
     MTFT_DeTrend             remove DC components from data 
     MTFT_Despike             remove spikes from the data
     MTFT_MHAFreq             Not sure
     MTFT_NDecFlt             Number of decade filter to apply to data
     MTFT_NPWCoef             number of prewhitening coefficients to apply
     MTFT_NotchFlt            apply a notch filter
     MTFT_NotchFrq            notch filter frequencies
     MTFT_NotchWidth          notch width
     MTFT_PWFilter            prewhitening filter
     MTFT_SpikeDev            spike deviation 
     MTFT_SpikePnt            number of points per spike
     MTFT_StackFlt            stack data
     MTFT_StackFrq            stack frequencies
     MTFT_StackTaper          taper of stack filter
     MTFT_SysCal              system calibraions 
     MTFT_T0OffsetMax         offset of time 0
     MTFT_TSPlot_ChnRange     time series plot channel range
     MTFT_TSPlot_PntRange     time series plot number of points 
     MTFT_Version             version of mtft
     MTFT_WindowLength        length of decimation window
     MTFT_WindowOverlap       amount of overlap between decimation windows
     MTFT_WindowTaper         taper on decimation window
     Remote_Component         remote reference components
     Remote_Path              path to remote reference data
     Remote_Rotation          rotation on remote reference data
     Rx_HPR                   rotation of data
     Setup_ID                 setup ID, in case there were different setups
     Setup_Number             setup number
     Setup_Use                use setup in processing
     TS_FrqBand               time series frequency band
     TS_Number                time series number
     TS_T0Error               time series time zero error
     TS_T0Offset              time series time zero offset  
     Unit_Length              length units
     cache_path               path to .cac files
     log_lines                list of lines to write to a log file 
     meta_dict                dictionary of meta data 
     meta_keys                keys of meta data
     new_remote_path          new remote referenc path
     num_comp                 number of components
     rr_tdiff_dict            difference in remote reference and data
     setup_keys               key words for setup
     setup_lst                list of setup values
     sort_ts_lst              sorted time series list
     ts_info_keys             keys for time series information
     ts_info_lst              list of time series information
     value_lst                values of meta data
     verbose                  [ True | False ] to write things to screen
    ========================= =================================================
    
    ========================== ================================================
     methods                   description
    ========================== ================================================
     compute_number_of_setups  compute the number of setups from the data
     get_rr_ts                 get remote reference time series 
     get_survey_info           get survey information
     get_ts_info_lst           get time series information
     make_value_dict           make values for time series
     read_cfg                  read configuration file
     set_remote_reference_info set important remote reference data
     set_values                set values of meta_dict
     write_mtft_cfg           write MTFT configuration file
    ========================== ================================================
    
    :Example: ::
        
        >>> import mtpy.usgs.zonge as zonge
        >>> mtft = zonge.ZongeMTFT()
        >>> mtft.write_mtft_cfg(r"/home/mt/mt01/Merged", \
                                'mt01', \
                                rrstation='rr01',\
                                survey_file=r"/home/mt/survey/cfg")
    
    """
    def __init__(self):
        
        #--> standard MTFT meta data
        self.MTFT_Version = '1.12v'
        self.MTFT_MHAFreq = 3
        self.MTFT_WindowTaper = '4 Pi Prolate'
        self.MTFT_WindowLength = 64
        self.MTFT_WindowOverlap = 48
        self.MTFT_NDecFlt = 5
        self.MTFT_PWFilter = 'Auto-Regression'
        self.MTFT_NPWCoef = 5
        self.MTFT_DeTrend = 'Yes'
        self.MTFT_Despike = 'Yes'
        self.MTFT_SpikePnt = 1
        self.MTFT_SpikeDev = 4
        self.MTFT_NotchFlt = 'Yes'
        self.MTFT_NotchFrq =  list(np.arange(60, 600, 60))
        self.MTFT_NotchWidth = list(range(1,len(self.MTFT_NotchFrq)+1))
        self.MTFT_StackFlt = 'No'
        self.MTFT_StackTaper = 'Yes'
        self.MTFT_StackFrq = [60,1808,4960]
        self.MTFT_SysCal = 'Yes'
        self.MTFT_BandFrq = [32, 256, 512, 1024, 2048, 4096, 32768] 
        self.MTFT_BandFrqMin = [7.31000E-4, 7.31000E-4, 7.31000E-4, 0.25, 0.25, 
                                0.25, 1]
        self.MTFT_BandFrqMax = [10, 80, 160, 320, 640, 1280, 10240]
        self.MTFT_TSPlot_PntRange = 4096
        self.MTFT_TSPlot_ChnRange = '1000'+',1000'*17 
        
        #--> time series meta data
        self.TS_Number = 3
        self.TS_FrqBand = [32, 256, 512]
        self.TS_T0Offset = [0, 0, 0]
        self.TS_T0Error = [0, 0, 0]

        #--> setup parameters
        self.Setup_Number = 1
        self.Setup_ID = 1
        self.Setup_Use = 'Yes'
        self.setup_lst = []

        #--> survey parameters
        self.Unit_Length = 'm'
        self.Chn_Cmp = ['Hx', 'Hy', 'Hz', 'Ex', 'Ey']
        self.Chn_ID = ['2314', '2324', '2334', '1', '1']
        self.Chn_Gain = [1, 1, 1, 1, 1]
        self.Chn_Length = [100]*5
        self.Chn_Azimuth = [0, 270, 90, 0, 270]
        self.Chn_dict = dict([(chkey, [cid, cg, cl]) for chkey, cid, cg, cl in 
                               zip(self.Chn_Cmp, self.Chn_ID, self.Chn_Gain,
                                   self.Chn_Length)])
                                   
        self.Chn_Cmp_lst = []
        self.num_comp = len(self.Chn_Cmp)
        self.Ant_FrqMin = 7.31E-4
        self.Ant_FrqMax = 10240
        self.Rx_HPR = [0, 0, 180]
        self.Remote_Component = 'Hx,Hy'
        self.Remote_HPR = [0, 0, 0]
        self.Remote_Rotation = 0
        self.Remote_Path = ''
        self.cache_path = None
        self.new_remote_path = ''
        
        self.verbose = True
        self.log_lines = []
        
        #info dict
        self.ts_info_keys = ['File#', 'Setup', 'SkipWgt', 'LocalFile', 
                             'RemoteFile', 'LocalBlock', 'RemoteBlock', 
                             'LocalByte', 'RemoteByte', 'Date', 'Time0', 
                             'T0Offset', 'ADFrequency', 'NLocalPnt',
                             'NRemotePnt', 'ChnGain1', 'ChnGain2', 'ChnGain3',
                             'ChnGain4', 'ChnGain5']
        self.ts_info_lst = []
                             
        self.meta_keys = ['MTFT.Version', 
                          'MTFT.MHAFreq',
                          'MTFT.WindowTaper',
                          'MTFT.WindowLength',
                          'MTFT.WindowOverlap', 
                          'MTFT.NDecFlt',
                          'MTFT.PWFilter',
                          'MTFT.NPWCoef',
                          'MTFT.DeTrend', 
                          'MTFT.Despike',
                          'MTFT.SpikePnt',
                          'MTFT.SpikeDev',
                          'MTFT.NotchFlt',
                          'MTFT.NotchFrq',
                          'MTFT.NotchWidth',
                          'MTFT.StackFlt',
                          'MTFT.StackTaper',
                          'MTFT.StackFrq',
                          'MTFT.SysCal',
                          'MTFT.BandFrq',
                          'MTFT.BandFrqMin',
                          'MTFT.BandFrqMax',
                          'MTFT.TSPlot.PntRange',
                          'MTFT.TSPlot.ChnRange',
                          'Setup.Number',
                          'TS.Number',
                          'TS.FrqBand',
                          'TS.T0Offset',
                          'TS.T0Error',
                          'setup_lst']
                          
        self.setup_keys = ['Setup.ID',
                           'Setup.Use',
                           'Unit.Length',
                           'Chn.Cmp',
                           'Chn.ID',
                           'Chn.Length',
                           'Chn.Azimuth',
                           'Chn.Gain',
                           'Ant.FrqMin',
                           'Ant.FrqMax',
                           'Rx.HPR',
                           'Remote.HPR',
                           'Remote.Rotation',
                           'Remote.Path']
                          
        self.value_lst = []
        self.meta_dict = None
        self.make_value_dict()
        
        self.rr_tdiff_dict = {'256':'060000', '1024':'002000', '4096':'000500'}
        
     
    def make_value_dict(self):
        """
        make value dictionary with all the important information
        """
        self.value_lst = [self.__dict__[key.replace('.', '_')] 
                          for key in self.meta_keys]
                          
        self.meta_dict = dict([(mkey, mvalue) for mkey, mvalue in
                               zip(self.meta_keys, self.value_lst)])
                               
    def set_values(self):
        """
        from values in meta dict set attribute values
        
        """
        for key in list(self.meta_dict.keys()):
            setattr(self, key.replace('.', '_'), self.meta_dict[key])
    
    def sort_ts_lst(self):
        """
        sort the time series list such that all the same sampling rates are
        in sequential order, this needs to be done to get reasonable 
        coefficients out of mtft
        """
        
        if self.ts_info_lst == []:
            return
        
        new_ts_lst = sorted(self.ts_info_lst, key=lambda k: k['ADFrequency'])
        
        for ii, new_ts in enumerate(new_ts_lst, 1):
            new_ts['File#'] = ii
            for tkey in self.ts_info_keys:
                if not type(new_ts[tkey]) is str:
                    new_ts[tkey] = str(new_ts[tkey])
                    
        self.ts_info_lst = new_ts_lst
        
    def get_rr_ts(self, ts_info_lst, remote_path=None):
        """
        get remote reference time series such that it has the same starting
        time and number of points as the collected time series
        
        Arguments:
        ----------
            **ts_info_lst** : list of dictionaries that relate to the time 
                              series to be processed.
            
            **remote_path** : directory path of merged cache files to use as
                              a remote reference.
        
        """

        if remote_path is not None:
            self.Remote_Path = remote_path
            
        if self.Remote_Path is None or self.Remote_Path == '':
            return 
            
        self.new_remote_path = os.path.join(self.cache_path, 'RR')
        if not os.path.exists(self.new_remote_path):
            os.mkdir(self.new_remote_path)
            
            
        new_ts_info_lst = []  
        rrfnlst = [rrfn for rrfn in os.listdir(self.Remote_Path) 
                   if rrfn.find('.cac')>0]
                       
        for ts_dict in ts_info_lst:
            local_zc = zen.ZenCache()
            local_zc.read_cache_metadata(os.path.join(self.cache_path,
                                                      ts_dict['LocalFile']))
            try:
                local_start_date = local_zc.meta_data['DATA.DATE0'][0]
                local_start_time = \
                        local_zc.meta_data['DATA.TIME0'][0].replace(':','')
            except KeyError:
                lsd_lst = \
                        local_zc.meta_data['DATE0'][0].split('/')
                local_start_date = '20{0}-{1}-{2}'.format(lsd_lst[2],
                                                          lsd_lst[0],
                                                          lsd_lst[1]) 
                        
                local_start_time = \
                        local_zc.meta_data['TIMEO'][0].replace(':','') 
                                                                 
            local_df = local_zc.meta_data['TS.ADFREQ'][0]
            local_npts = int(local_zc.meta_data['TS.NPNT'][0])
            tdiff = self.rr_tdiff_dict[local_df]
            
            print('='*60)
            print(ts_dict['LocalFile'], local_start_date, local_start_time)
            self.log_lines.append('='*60+'\n')
            self.log_lines.append('{0} {1} {2} \n'.format(ts_dict['LocalFile'],
                                  local_start_date, local_start_time)) 
            
            rrfind = False
            #look backwards because if a new file was already created it will
            #be found before the original file
            for rrfn in rrfnlst[::-1]:
                remote_zc = zen.ZenCache()
                remote_zc.read_cache_metadata(os.path.join(self.Remote_Path,
                                                           rrfn))
                
                try:
                    remote_start_date = remote_zc.meta_data['DATA.DATE0'][0]
                    remote_start_time = \
                        remote_zc.meta_data['DATA.TIME0'][0].replace(':','')
                except KeyError:
                    remote_start_date = \
                        remote_zc.meta_data['DATE0'][0].replace('/','-')
                    remote_start_time = \
                        remote_zc.meta_data['TIMEO'][0].replace(':','')
                remote_df = remote_zc.meta_data['TS.ADFREQ'][0]
                remote_npts = int(remote_zc.meta_data['TS.NPNT'][0])

                if local_start_date == remote_start_date and \
                   local_df == remote_df and \
                   local_start_time[0:2] == remote_start_time[0:2]:
                    print(rrfn, remote_start_date, remote_start_time)
                    self.log_lines.append('{0} {1} {2}\n'.format(rrfn, 
                                          remote_start_date, 
                                          remote_start_time))
                                          
                    if local_start_time == remote_start_time:
                        if local_npts == remote_npts:
                            ts_dict['RemoteFile'] = rrfn
                            ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                            ts_dict['RemoteByte'] = ts_dict['LocalByte']
                            ts_dict['NRemotePnt'] = ts_dict['NLocalPnt']
                            for ii in range(self.num_comp+1,
                                            self.num_comp+3):
                                ts_dict['ChnGain{0}'.format(ii)] = '1'
                            new_ts_info_lst.append(ts_dict)
                            
                            #copy remote referenc data to local directory
                            shutil.copy(os.path.join(self.Remote_Path, rrfn),
                                        os.path.join(self.new_remote_path, 
                                                     rrfn))
                            rrfind = True
                            break
                        
                        #if time series is longer than remote reference
                        elif remote_npts < local_npts:
                            print('{0} local_npts > remote_npts {0}'.format('*'*4))
                            self.log_lines.append('{0} local_npts > remote_npts {0}\n'.format('*'*4))
                            #read in cache file
                            local_zc.read_cache(os.path.join(self.cache_path,
                                                  ts_dict['LocalFile']))
                            #resize local ts accordingly
                            local_zc.ts = np.resize(local_zc.ts, 
                                                    (remote_npts,
                                                     local_zc.ts.shape[1]))

                            #reset some meta data 
                            local_zc.meta_data['TS.NPNT'] = \
                                            [str(local_zc.ts.shape[0])]
                            
                            print('Resized Local TS in {0} to {1}'.format(
                               os.path.join(self.cache_path, 
                                            ts_dict['LocalFile']),
                               local_zc.ts.shape))
                            self.log_lines.append('Resized Local TS in {0} to {1}\n'.format(
                               os.path.join(self.cache_path, 
                                            ts_dict['LocalFile']),
                               local_zc.ts.shape))
                            #rewrite the cache file
                            local_zc.rewrite_cache_file()

                            #reset some of the important parameters
                            ts_dict['LocalFile'] = \
                                    os.path.basename(local_zc.save_fn_rw)
                            ts_dict['RemoteFile'] = rrfn
                            ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                            ts_dict['RemoteByte'] = ts_dict['LocalByte']
                            ts_dict['NLocalPnt'] = local_zc.ts.shape[0]
                            ts_dict['NRemotePnt'] = local_zc.ts.shape[0]
                            for ii in range(self.num_comp+1,
                                            self.num_comp+3):
                                ts_dict['ChnGain{0}'.format(ii)] = '1'
                            new_ts_info_lst.append(ts_dict)
                            
                            #copy remote referenc data to local directory
                            shutil.copy(os.path.join(self.Remote_Path, rrfn),
                                        os.path.join(self.new_remote_path, 
                                                     rrfn))
                            break
                                                        
                        #if remote reference is longer than time series
                        elif remote_npts > local_npts:
                            print('{0} local_npts < remote_npts {0}'.format('*'*4))
                            self.log_lines.append('{0} local_npts < remote_npts {0}\n'.format('*'*4))
                            
                            remote_zc.read_cache(os.path.join(self.Remote_Path,
                                                              rrfn))
                            #resize remote ts accordingly
                            remote_zc.ts = np.resize(remote_zc.ts, 
                                                      (local_npts,
                                                       remote_zc.ts.shape[1]))
                            #reset some meta data 
                            remote_zc.meta_data['TS.NPNT'] = \
                                            [str(remote_zc.ts.shape[0])]
                            
                            print('Resized Remote TS in {0} to {1}'.format(
                                    os.path.join(self.Remote_Path, rrfn),
                                    remote_zc.ts.shape))
                            self.log_lines.append('Resized Remote TS in {0} to {1}\n'.format(
                                    os.path.join(self.Remote_Path, rrfn),
                                    remote_zc.ts.shape))
                                    
                            #rewrite the remote cache file 
                            remote_zc.rewrite_cache_file()

                            #reset some of the important parameters
                            ts_dict['RemoteFile'] = \
                                    os.path.basename(remote_zc.save_fn_rw)
                            ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                            ts_dict['RemoteByte'] = ts_dict['LocalByte']
                            ts_dict['NLocalPnt'] = local_npts
                            ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                            for ii in range(self.num_comp+1,
                                            self.num_comp+3):
                                ts_dict['ChnGain{0}'.format(ii)] = '1'
                            new_ts_info_lst.append(ts_dict)
                            
                            #copy remote referenc data to local directory
                            shutil.move(remote_zc.save_fn_rw,
                                        os.path.join(self.new_remote_path, 
                                     os.path.basename(remote_zc.save_fn_rw)))
                            
                            rrfind = True
                            break
                            
                    #if the starting time is different
                    elif abs(int(local_start_time)-int(remote_start_time)) < \
                                                    int(tdiff):
                        
                        local_hour = int(local_start_time[0:2])
                        local_minute = int(local_start_time[2:4])
                        local_second = int(local_start_time[4:])
                        
                        rr_hour = int(remote_start_time[0:2])
                        rr_minute = int(remote_start_time[2:4])
                        rr_second = int(remote_start_time[4:])
                        
                        hour_diff = (rr_hour-local_hour)*3600
                        minute_diff = (rr_minute-local_minute)*60
                        second_diff = rr_second-local_second
                        
                        time_diff = hour_diff+minute_diff+second_diff
                        skip_points = int(local_df)*abs(time_diff)
                        
                        remote_zc.read_cache(os.path.join(self.Remote_Path,
                                                          rrfn))
                        
                        #remote start time is later than local
                        if time_diff > 0:
                            
                            print(('Time difference is {0} seconds'.format(
                                                                time_diff)))
                            self.log_lines.append('Time difference is {0} seconds\n'.format(
                                                                time_diff))
                            print('Skipping {0} points in {1}'.format(
                                                skip_points,
                                                os.path.join(self.cache_path, 
                                                        ts_dict['LocalFile'])))
                            self.log_lines.append('Skipping {0} points in {1}\n'.format(
                                                skip_points,
                                                os.path.join(self.cache_path, 
                                                        ts_dict['LocalFile'])))
                                                
                            local_zc.read_cache(os.path.join(self.cache_path, 
                                                        ts_dict['LocalFile']))
                            
                            #resize local ts
                            local_zc.ts = local_zc.ts[skip_points:, :]
                            local_zc.meta_data['DATA.TIME0'] = \
                                    ['{0}:{1}:{2}'.format(
                                            local_hour+int(hour_diff/3600.),
                                            local_minute+int(minute_diff/60.),
                                            local_second+int(second_diff))]
                            
                            #if for some reason after reshaping the remote
                            #the local time series is still larger, cull
                            #the local to match the remote so mtft doesn't
                            #get angry
                            if remote_zc.ts.shape[0] < local_zc.ts.shape[0]:
                                print('{0} local_npts > remote_npts {0}'.format('*'*4))
                                self.log_lines.append('{0} local_npts > remote_npts {0}\n'.format('*'*4))                                
                                #read in cache file
                                local_zc.read_cache(os.path.join(self.cache_path,
                                                      ts_dict['LocalFile']))
                                #resize local ts accordingly
                                local_zc.ts = np.resize(local_zc.ts, 
                                                        (remote_zc.ts.shape[0],
                                                         local_zc.ts.shape[1]))

                                #reset some meta data 
                                local_zc.meta_data['TS.NPNT'] = \
                                                [str(local_zc.ts.shape[0])]
                                
                                
                                print('Resized Local TS in {0} to {1}'.format(
                                   os.path.join(self.cache_path, 
                                                ts_dict['LocalFile']),
                                   local_zc.ts.shape))
                                self.log_lines.append('Resized Local TS in {0} to {1}\n'.format(
                                   os.path.join(self.cache_path, 
                                                ts_dict['LocalFile']),
                                   local_zc.ts.shape))
                                #rewrite the cache file
                                local_zc.rewrite_cache_file()

                                #reset some of the important parameters
                                ts_dict['LocalFile'] = \
                                        os.path.basename(local_zc.save_fn_rw)
                                ts_dict['RemoteFile'] = rrfn
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = local_zc.ts.shape[0]
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                
                                #copy remote to local directory
                                shutil.copy(os.path.join(self.Remote_Path, 
                                                         rrfn),
                                            os.path.join(self.new_remote_path, 
                                                         rrfn))
                                rrfind = True
                                break
                                
                            #reshape local file if number of points is larger
                            #than the remote reference.
                            elif remote_zc.ts.shape[0] > local_zc.ts.shape[0]:
                                print('{0} local_npts < remote_npts {0}'.format('*'*4))
                                self.log_lines.append('{0} local_npts < remote_npts {0}\n'.format('*'*4))                                
                                #reset local meta data 
                                local_zc.meta_data['TS.NPNT'] = \
                                                [str(local_zc.ts.shape[0])]
                                
                                #rewrite the local cache file
                                local_zc.rewrite_cache_file()
                            
                                #resize remote ts accordingly
                                remote_zc.ts = np.resize(remote_zc.ts, 
                                                          (local_zc.ts.shape[0],
                                                           remote_zc.ts.shape[1]))
                                #reset some meta data 
                                remote_zc.meta_data['TS.NPNT'] = \
                                                [str(remote_zc.ts.shape[0])]
                                
                                print('Resized Remote TS in {0} to {1}'.format(
                                        os.path.join(self.Remote_Path, rrfn),
                                        remote_zc.ts.shape))
                                self.log_lines.append('Resized Remote TS in {0} to {1}\n'.format(
                                        os.path.join(self.Remote_Path, rrfn),
                                        remote_zc.ts.shape))
                                        
                                #rewrite the remote cache file 
                                remote_zc.rewrite_cache_file()

                                
                                #reset some of the important parameters
                                ts_dict['LocalFile'] = \
                                        os.path.basename(local_zc.save_fn_rw)
                                ts_dict['RemoteFile'] = \
                                        os.path.basename(remote_zc.save_fn_rw)
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = local_zc.ts.shape[0]
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                
                                #copy remote referenc data to local directory
                                shutil.move(remote_zc.save_fn_rw,
                                        os.path.join(self.new_remote_path, 
                                       os.path.basename(remote_zc.save_fn_rw)))
                                rrfind = True
                                break
                            
                            elif remote_zc.ts.shape[0] == local_npts:
                                #reset local meta data 
                                local_zc.meta_data['TS.NPNT'] = \
                                                [str(local_zc.ts.shape[0])]
                                
                                #rewrite the local cache file
                                local_zc.rewrite_cache_file()
                                
                                #reset some of the important parameters
                                ts_dict['LocalFile'] = \
                                        os.path.basename(local_zc.save_fn_rw)
                                ts_dict['RemoteFile'] = rrfn
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = local_zc.ts.shape[0]
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                
                                #copy remote reference to local directory
                                shutil.copy(os.path.join(self.Remote_Path, 
                                                         rrfn),
                                            os.path.join(self.new_remote_path, 
                                                         rrfn))
                                rrfind = True
                                break
                    
                        #local start time is later than remote start time                        
                        elif time_diff < 0:
                            
                            print(('Time difference is {0} seconds'.format(
                                                                    time_diff)))
                            self.log_lines.append('Time difference is {0} seconds\n'.format(
                                                                    time_diff))
                            print('Skipping {0} points in {1}'.format(
                                                skip_points,
                                                os.path.join(self.Remote_Path,
                                                             rrfn)))
                            self.log_lines.append('Skipping {0} points in {1}\n'.format(
                                                skip_points,
                                                os.path.join(self.Remote_Path,
                                                             rrfn)))
                                
                            
                            #resize remote reference
                            new_rr_ts = remote_zc.ts[skip_points:, :]
                            
                            remote_zc.ts = new_rr_ts
                            remote_zc.meta_data['DATA.TIME0'] = \
                                    ['{0}:{1}:{2}'.format(
                                            rr_hour-int(hour_diff/3600.),
                                            rr_minute-int(minute_diff/60.),
                                            rr_second-int(second_diff))]
                            
                            #if for some reason after reshaping the remote
                            #the local time series is still larger, cull
                            #the local to match the remote so mtft doesn't
                            #get angry
                            if remote_zc.ts.shape[0] < local_npts:
                                print('{0} local_npts > remote_npts {0}'.format('*'*4))
                                self.log_lines.append('{0} local_npts > remote_npts {0}\n'.format('*'*4))                                
                                #reset remote meta data 
                                remote_zc.meta_data['TS.NPNT'] = \
                                                [str(remote_zc.ts.shape[0])]
                                
                                #rewrite the remote cache file 
                                remote_zc.rewrite_cache_file()
                                
                                #read in cache file
                                local_zc.read_cache(os.path.join(self.cache_path,
                                                      ts_dict['LocalFile']))
                                #resize local ts accordingly
                                local_zc.ts = np.resize(local_zc.ts, 
                                                        (remote_zc.ts.shape[0],
                                                         local_zc.ts.shape[1]))

                                #reset some meta data 
                                local_zc.meta_data['TS.NPNT'] = \
                                                [str(local_zc.ts.shape[0])]
                               
                                print('Resized Local TS in {0} to {1}'.format(
                                   os.path.join(self.cache_path, 
                                                ts_dict['LocalFile']),
                                   local_zc.ts.shape))
                                self.log_lines.append('Resized Local TS in {0} to {1}\n'.format(
                                   os.path.join(self.cache_path, 
                                                ts_dict['LocalFile']),
                                   local_zc.ts.shape))
                                   
                                #rewrite the cache file
                                local_zc.rewrite_cache_file()

                                #reset some of the important parameters
                                ts_dict['LocalFile'] = \
                                        os.path.basename(local_zc.save_fn_rw)
                                ts_dict['RemoteFile'] = \
                                        os.path.basename(remote_zc.save_fn_rw)
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = local_zc.ts.shape[0]
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                
                                #copy remote referenc data to local directory
                                shutil.move(remote_zc.save_fn_rw,
                                        os.path.join(self.new_remote_path, 
                                       os.path.basename(remote_zc.save_fn_rw)))
                                rrfind = True
                                break
                                
                            #reshape local file if number of points is larger
                            #than the remote reference.
                            elif remote_zc.ts.shape[0] > local_npts:
                                print('{0} local_npts < remote_npts {0}'.format('*'*4))
                                self.log_lines.append('{0} local_npts < remote_npts {0}\n'.format('*'*4))                                
                                #resize remote ts accordingly
                                remote_zc.ts = np.resize(remote_zc.ts, 
                                                          (local_npts,
                                                           remote_zc.ts.shape[1]))
                                #reset some meta data 
                                remote_zc.meta_data['TS.NPNT'] = \
                                                [str(remote_zc.ts.shape[0])]
                                
                                print('Resized Remote TS in {0} to {1}'.format(
                                        os.path.join(self.Remote_Path, rrfn),
                                        remote_zc.ts.shape))
                                self.log_lines.append('Resized Remote TS in {0} to {1}\n'.format(
                                        os.path.join(self.Remote_Path, rrfn),
                                        remote_zc.ts.shape))
                                        
                                #rewrite the remote cache file 
                                remote_zc.rewrite_cache_file()
                                
                                #reset some of the important parameters
                                ts_dict['RemoteFile'] = \
                                        os.path.basename(remote_zc.save_fn_rw)
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = local_npts
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                #copy remote referenc data to local directory
                                shutil.move(remote_zc.save_fn_rw,
                                        os.path.join(self.new_remote_path, 
                                       os.path.basename(remote_zc.save_fn_rw)))
                                rrfind = True
                                break
                            
                            elif remote_zc.ts.shape[0] == local_npts:
                                #reset local meta data 
                                remote_zc.meta_data['TS.NPNT'] = \
                                                [str(remote_zc.ts.shape[0])]
                                
                                #rewrite the local cache file
                                remote_zc.rewrite_cache_file()
                                
                                #reset some of the important parameters
                                ts_dict['RemoteFile'] = rrfn
                                ts_dict['RemoteBlock'] = ts_dict['LocalBlock']
                                ts_dict['RemoteByte'] = ts_dict['LocalByte']
                                ts_dict['NLocalPnt'] = remote_zc.ts.shape[0]
                                ts_dict['NRemotePnt'] = remote_zc.ts.shape[0]
                                for ii in range(self.num_comp+1,
                                                self.num_comp+3):
                                    ts_dict['ChnGain{0}'.format(ii)] = '1'
                                new_ts_info_lst.append(ts_dict)
                                
                                #copy remote to local directory
                                shutil.move(remote_zc.save_fn_rw,
                                        os.path.join(self.new_remote_path, 
                                       os.path.basename(remote_zc.save_fn_rw)))
                                rrfind = True
                                break
                        
            if rrfind == False:
                print(('Did not find remote reference time series '
                       'for {0}'.format(ts_dict['LocalFile'])))
                self.log_lines.append('Did not find remote reference time series '
                       'for {0}\n'.format(ts_dict['LocalFile']))
                       
                new_ts_info_lst.append(ts_dict)
                        
        self.ts_info_lst = new_ts_info_lst
        
    def get_ts_info_lst(self, cache_path):
        """
        get information about time series and put it into dictionaries with
        keys according to header line in .cfg file
        
        Arguments:
        -----------
            **cache_path** : directory to cache files that are to be processed
        
        """

        #--> get .cac files and read meta data
        cc = 0
        
        self.cache_path = cache_path
        
        if len(self.ts_info_lst) == 0:
            for cfn in os.listdir(cache_path):
                if cfn[-4:] == '.cac' and cfn.find('$') == -1:
                    zc = zen.ZenCache()
                    zc.read_cache_metadata(os.path.join(cache_path, cfn))
                    self.Chn_Cmp_lst.append([md.capitalize() 
                                    for md in zc.meta_data['CH.CMP']
                                    if md.capitalize() in self.Chn_Cmp])
                    
                    info_dict = dict([(key, []) for key in self.ts_info_keys])
                    #put metadata information into info dict
                    info_dict['File#'] = cc+1
                    info_dict['Setup'] = 1
                    info_dict['SkipWgt'] = 1
                    info_dict['LocalFile'] = cfn
                    info_dict['RemoteFile'] = ''
                    info_dict['LocalBlock'] = cc
                    info_dict['RemoteBlock'] = ''
                    info_dict['LocalByte'] = 65
                    info_dict['RemoteByte'] = ''
                    try:
                        info_dict['Date'] = zc.meta_data['DATA.DATE0'][0]
                    except KeyError:
                        info_dict['Date'] = zc.meta_data['DATE0'][0]
                    if info_dict['Date'].find('/') >= 0:
                        dlst = info_dict['Date'].split('/')
                        info_dict['Date'] = '20{0}-{1}-{2}'.format(dlst[2], 
                                                                   dlst[0],
                                                                   dlst[1])
                    info_dict['Time0'] = '0'
                    #try:
                    #    info_dict['Time0'] = zc.meta_data['Data.Time0'][0]
                    #except KeyError:
                    #    info_dict['Time0'] = zc.meta_data['TIMEO'][0]
                    info_dict['T0Offset'] = '0'
                    info_dict['ADFrequency'] = int(zc.meta_data['TS.ADFREQ'][0])
                    info_dict['NLocalPnt'] = int(zc.meta_data['TS.NPNT'][0])
                    info_dict['NRemotePnt'] = ''
                    for ii in range(1,self.num_comp+1):
                        info_dict['ChnGain{0}'.format(ii)] = '1'               
                    
                    cc += 1
                    self.ts_info_lst.append(info_dict)
    
    def compute_number_of_setups(self):
        """
        get number of setups and set all the necessary values
        
        **Need to match setups with time series info**
        
        """
        self.setup_lst = []
        len_lst = []
        ii = 1
        for cc in self.Chn_Cmp_lst:
            comp_num = len(cc)
            if comp_num not in len_lst:
                len_lst.append(comp_num)
                setup_dict = {}
                setup_dict['Setup.ID'] = ii
                setup_dict['Setup.Use'] = 'Yes'
                setup_dict['Unit.Length'] = 'm'
                setup_dict['Chn.Cmp'] = cc
                setup_dict['Chn.ID'] = list(range(1,comp_num+1))
                setup_dict['Chn.Length'] = [100]*len(cc)
                setup_dict['Chn.Gain'] = [1]*len(cc)
                setup_dict['Ant.FrqMin'] = self.Ant_FrqMin
                setup_dict['Ant.FrqMax'] = self.Ant_FrqMax
                setup_dict['Rx.HPR'] = self.Rx_HPR
                setup_dict['Remote.Component'] = self.Remote_Component
                setup_dict['Remote.Rotation'] = self.Remote_Rotation
                setup_dict['Remote.Path'] = self.Remote_Path
                setup_dict['chn_dict'] = dict([(chkey, [cid, cg, cl]) 
                                                for chkey, cid, cg, cl in 
                                                zip(setup_dict['Chn.Cmp'],
                                                    setup_dict['Chn.ID'],
                                                    setup_dict['Chn.Gain'],
                                                    setup_dict['Chn.Length'])])
                ts_key_skip = len(self.ts_info_keys)-(5-len(cc)) 
                setup_dict['ts_info_keys'] = self.ts_info_keys[:ts_key_skip]
                self.setup_lst.append(setup_dict)
                ii += 1
        self.Setup_Number = len(self.setup_lst)
                
    def set_remote_reference_info(self, remote_path):
        """
        set the remote reference information in ts_info_lst
        
        Arguments:
        ----------
            **remote_path** : directory of remote reference cache files
            
        """
        
        if remote_path is None:
            return
        
        self.Remote_Path = remote_path
        for setup_dict in self.setup_lst:
            setup_dict['Chn.Cmp'] += ['Hxr', 'Hyr']
            setup_dict['Chn.ID'] += ['2284', '2274']
            setup_dict['Chn.Length'] += [100]*2
            setup_dict['Chn.Gain'] += [1]*2
            setup_dict['Ant.FrqMin'] = self.Ant_FrqMin
            setup_dict['Ant.FrqMax'] = self.Ant_FrqMax
            setup_dict['Rx.HPR'] = self.Rx_HPR
            setup_dict['Remote.Component'] = self.Remote_Component
            setup_dict['Remote.Rotation'] = self.Remote_Rotation
            setup_dict['Remote.Path'] = self.new_remote_path+os.path.sep
            setup_dict['chn_dict'] = dict([(chkey, [cid, cg, cl]) 
                                            for chkey, cid, cg, cl in 
                                            zip(setup_dict['Chn.Cmp'],
                                                setup_dict['Chn.ID'],
                                                setup_dict['Chn.Gain'],
                                                setup_dict['Chn.Length'])])
            num_comp = len(setup_dict['Chn.Cmp'])
            setup_dict['ts_info_keys'] += ['ChnGain{0}'.format(ii) 
                                           for ii in range(num_comp-1, 
                                                           num_comp+1)]
            
                                                           
    def get_survey_info(self, survey_file, station_name, rr_station_name=None):
        """
        extract information from survey file
        
        Arguments:
        ----------
            **survey_file** : string
                              full path to survey config file created by 
                              mtpy.utils.configfile
                              
            **station_name** : string
                               full station name
            
            **rr_station_name** : string
                                  full station name of remote reference
        
        """

        if survey_file is None:
            return
            
        #--> get information from survey file
        for setup_dict in self.setup_lst:
            sdict = mtcf.read_survey_configfile(survey_file)
            
            try:
                survey_dict = sdict[station_name.upper()]
                try:
                    setup_dict['chn_dict']['Hx'][0] = survey_dict['hx']
                except KeyError:
                    print('No hx data')
                    self.log_lines.append('No hx data from survey file.\n')
                try:
                    setup_dict['chn_dict']['Hy'][0] = survey_dict['hy']
                except KeyError:
                    print('No hy data')
                    self.log_lines.append('No hy data from survey file.\n')
                
                try:
                    if survey_dict['hz'].find('*') >= 0:
                        setup_dict['chn_dict']['Hz'][0] = '3'
                    else:
                        setup_dict['chn_dict']['Hz'][0] = survey_dict['hz']
                except KeyError:
                    print('No hz data')
                    self.log_lines.append('No hz data from survey file.\n')
                
                try:
                    setup_dict['chn_dict']['Ex'][2] = \
                                                 survey_dict['e_xaxis_length']
                except KeyError:
                    print('No ex data')
                    self.log_lines.append('No ex data from survey file.\n')
                try:
                    setup_dict['chn_dict']['Ey'][2] = \
                                                 survey_dict['e_yaxis_length']
                except KeyError:
                    print('No ey data')
                    self.log_lines.append('No ey data from survey file.\n')
                    
            except KeyError:
                print(('Could not find survey information from ' 
                       '{0} for {1}'.format(survey_file, station_name)))
                self.log_lines.append('Could not find survey information from ' 
                       '{0} for {1}\n'.format(survey_file, station_name))
            
            if rr_station_name is not None:
                try:
                    survey_dict = sdict[rr_station_name.upper()]
                    try:
                        setup_dict['chn_dict']['Hxr'][0] = survey_dict['hx']
                    except KeyError:
                        print('No hxr data')
                        self.log_lines.append('No hxr data from survey file.\n')
                        
                    try:
                        setup_dict['chn_dict']['Hyr'][0] = survey_dict['hy']
                    except KeyError:
                        print('No hyr data')
                        self.log_lines.append('No hyr data from survey file.\n')
               
                except KeyError:
                    print(('Could not find survey information from ' 
                           '{0} for {1}'.format(survey_file, rr_station_name)))
                    self.log_lines.append('Could not find survey information from ' 
                           '{0} for {1}\n'.format(survey_file, rr_station_name))
                           
    def write_mtft_cfg(self, cache_path, station, rrstation=None,
                       remote_path=None, survey_file=None, save_path=None):
        """
        write a config file for mtft24 from the cache files in cache_path
        
        Arguments:
        ----------
            **cache_path** : string
                             directory to cache files to be processed
            
            **station** : string
                          full name of station to be processed
            
            **rrstation** : string
                            full name of remote reference station
                            
            **remote_path** : string
                              directory to remote reference cache files
                              
            **survey_file** : string
                             full path to survey config file, written in the
                             format of mtpy.utils.configfile
                             
            **save_path** : string
                            path to save mtft24.cfg file, if none saved to
                            cache_path\mtft24.cfg
        
        """
        
        if save_path is None:
            save_path = os.path.join(cache_path, 'mtft24.cfg')

        self.cache_path = cache_path
        
        #--> get information about the time series
        self.get_ts_info_lst(self.cache_path)
        
        #--> get number of components
        self.compute_number_of_setups()
        print('Number of setups = {0}'.format(self.Setup_Number))
        
        #--> get remote reference information if needed
        if remote_path is not None:
            if rrstation is None:
                rrstation = os.path.basename(os.path.dirname(
                                             os.path.dirname(remote_path)))

        self.get_rr_ts(self.ts_info_lst, remote_path=remote_path)
        
        self.set_remote_reference_info(remote_path)
        
        #--> sort the time series such that each section with the same sampling
        #    rate is in sequential order
        self.sort_ts_lst()
        self.TS_Number = len(self.ts_info_lst)
        
        #--> fill in data from survey file
        self.get_survey_info(survey_file, station, rr_station_name=rrstation)
        
        #make a dictionary of all the values to write file
        if self.new_remote_path is not '' or self.new_remote_path is not None:
            self.new_remote_path += os.path.sep
        
        #--> set a dictionary with all attributes
        self.make_value_dict()
       
        #--> write mtft24.cfg file
        cfid = open(save_path, 'w')
        cfid.write('\n')
        #---- write processing parameters ----
        for ii, mkey in enumerate(self.meta_keys[:-1]):

            if type(self.meta_dict[mkey]) is list:
                cfid.write('${0}={1}\n'.format(mkey, ','.join(['{0}'.format(mm) 
                                             for mm in self.meta_dict[mkey]])))
            else:
                cfid.write('${0}={1}\n'.format(mkey, self.meta_dict[mkey]))
            
            #blanks line before setup and ts number
            if ii == 24:
               cfid.write('\n')
        cfid.write('\n')
        #---- write setup parameters ----
        for setup_dict in self.setup_lst:
            #set channel information                  
            setup_dict['Chn.Gain'] = [setup_dict['chn_dict'][ckey][1] 
                                      for ckey in setup_dict['Chn.Cmp']]
            setup_dict['Chn.ID'] = [setup_dict['chn_dict'][ckey][0] 
                                      for ckey in setup_dict['Chn.Cmp']]
            setup_dict['Chn.Length'] = [setup_dict['chn_dict'][ckey][2] 
                                        for ckey in setup_dict['Chn.Cmp']]
            
            for ii, mkey in enumerate(self.setup_keys):
                #write setups
                if type(setup_dict[mkey]) is list:
                    cfid.write('${0}={1}\n'.format(mkey,
                               ','.join(['{0}'.format(mm) 
                               for mm in setup_dict[mkey]])))
                else:
                    cfid.write('${0}={1}\n'.format(mkey, setup_dict[mkey]))
            cfid.write('\n')
                                    

        #---- write time series information ----
        ts_key_len = np.array([len(sd['ts_info_keys']) 
                                for sd in self.setup_lst])
        ts_key_find = np.where(ts_key_len==ts_key_len.max())[0][0]
        self.ts_info_keys = self.setup_lst[ts_key_find]['ts_info_keys']

        cfid.write(','.join(self.ts_info_keys)+'\n')
        for cfn in self.ts_info_lst:
            try:
                cfid.write(','.join([cfn[ikey] 
                                     for ikey in self.ts_info_keys])+'\n')   
            except KeyError:
                pass                            
        
        cfid.close()
        print('Wrote config file to {0}'.format(save_path))
        self.log_lines.append('Wrote config file to {0}\n'.format(save_path))
        
        #write log file
        lfid = open(os.path.join(os.path.dirname(save_path), 'MTFTcfg.log'),'w')
        lfid.writelines(self.log_lines)
        lfid.close()
        
    def read_cfg(self, cfg_fn):
        """
        read a mtft24.cfg file
        
        Arguments:
        ----------
            **cfg_fn** : full path to mtft24.cfg file to read.
        """
        
        if not os.path.isfile(cfg_fn):
            raise IOError('{0} does not exist'.format(cfg_fn))
            
        cfid = open(cfg_fn, 'r')
        clines = cfid.readlines()
        info_lst = []
        setup_lst = []
        self.meta_dict = {}
        for cline in clines:
            if cline[0] == '$':
                clst = cline[1:].strip().split('=')
                key = clst[0]
                value = clst[1].split(',')
                
                if key.find('MTFT') == 0 or \
                   key.find('Setup.Number') == 0 or\
                   key.find('TS') == 0:
                    self.meta_dict[key] = value
                elif clst[0].find('Setup.ID') == 0:
                    setup_lst.append({key:value})
                    ss = int(clst[1])-1
                    setup_lst[ss][key] = value
                else:
                    setup_lst[ss][key] = value
                    if key.find('Chn') == 0:
                        self.meta_dict[key] = value
                
            elif cline.find('.cac') > 0:
                info_dict = {}
                clst = cline.strip().split(',')
                for ckey, cvalue in zip(self.ts_info_keys, clst):
                    info_dict[ckey] = cvalue
                info_lst.append(info_dict)
        self.ts_info_lst = info_lst
        
        cfid.close()
        
        self.meta_dict['setup_lst'] = setup_lst
        
        self.set_values()
        self.make_value_dict()
        
#==============================================================================
# Deal with mtedit outputs  
#==============================================================================
class ZongeMTEdit():
    """
    deal with input and output config files for mtedi.  
    
    This is not used as much, but works if you need it
    """
    
    def __init__(self):
        self.meta_keys = ['MTEdit:Version', 'Auto.PhaseFlip', 
                          'PhaseSlope.Smooth', 'PhaseSlope.toZMag',
                          'DPlus.Use', 'AutoSkip.onDPlus', 'AutoSkip.DPlusDev']
                          
        self.mtedit_version = '3.10d applied on {0}'.format(time.ctime())
        self.phase_flip = 'No'
        self.phaseslope_smooth = 'Minimal'
        self.phaseslope_tozmag = 'Yes'
        self.dplus_use = 'Yes'
        self.autoskip_ondplus = 'No'
        self.autoskip_dplusdev = 500.0
        
        self.meta_dict = None
        self.meta_lst = None
        
        self.cfg_fn = None
        
        self.param_header = ['Frequency  ', 'AResXYmin', 'AResXYmax', 
                              'ZPhzXYmin', 'ZPhzXYmax', 'AResYXmin', 
                              'AResYXmax', 'ZPhzYXmin', 'ZPhzYXmax', 
                              'CoherXYmin', 'CoherYXmin', 'CoherXYmax', 
                              'CoherYXmax', 'ExMin', 'ExMax', 'EyMin', 
                              'EyMax', 'HxMin', 'HxMax', 'HyMin', 'HyMax', 
                              'NFC/Stack']
                              
        self.freq_lst = [7.32420000e-04, 9.76560000e-04, 1.22070000e-03, 
                         1.46480000e-03, 1.95310000e-03, 2.44140000e-03,  
                         2.92970000e-03, 3.90620000e-03, 4.88280000e-03, 
                         5.85940000e-03, 7.81250000e-03, 9.76560000e-03, 
                         1.17190000e-02, 1.56250000e-02, 1.95310000e-02, 
                         2.34380000e-02, 3.12500000e-02, 3.90620000e-02, 
                         4.68750000e-02, 6.25000000e-02, 7.81250000e-02, 
                         9.37500000e-02, 1.25000000e-01, 1.56200000e-01, 
                         1.87500000e-01, 2.50000000e-01, 3.12500000e-01, 
                         3.75000000e-01, 5.00000000e-01, 6.25000000e-01, 
                         7.50000000e-01, 1.00000000e+00, 1.25000000e+00, 
                         1.50000000e+00, 2.00000000e+00, 2.50000000e+00, 
                         3.00000000e+00, 4.00000000e+00, 5.00000000e+00, 
                         6.00000000e+00, 8.00000000e+00, 1.00000000e+01, 
                         1.20000000e+01, 1.60000000e+01, 2.00000000e+01, 
                         2.40000000e+01, 3.20000000e+01, 4.00000000e+01, 
                         4.80000000e+01, 6.40000000e+01, 8.00000000e+01, 
                         9.60000000e+01, 1.28000000e+02, 1.60000000e+02, 
                         1.92000000e+02, 2.56000000e+02, 3.20000000e+02, 
                         3.84000000e+02, 5.12000000e+02, 6.40000000e+02, 
                         7.68000000e+02, 1.02400000e+03, 1.28000000e+03, 
                         1.53600000e+03, 2.04800000e+03, 2.56000000e+03, 
                         3.07200000e+03, 4.09600000e+03, 5.12000000e+03, 
                         6.14400000e+03, 8.19200000e+03, 1.02400000e+04, 
                         0.00000000e+00]
                         
        self.num_freq = len(self.freq_lst)
        
        #--> default parameters for MTEdit                 
        self.AResXYmin = [1.0e-2]*self.num_freq
        self.AResXYmax = [1.0e6]*self.num_freq
        self.ZPhzXYmin = [-3150.]*self.num_freq
        self.ZPhzXYmax = [3150.]*self.num_freq
        
        self.AResYXmin = [1.0e-2]*self.num_freq
        self.AResYXmax = [1.0e6]*self.num_freq
        self.ZPhzYXmin = [-3150.]*self.num_freq
        self.ZPhzYXmax = [3150.]*self.num_freq
        
        self.CoherXYmin = [0.6]*self.num_freq
        self.CoherYXmin = [0.6]*self.num_freq
        self.CoherXYmax = [0.999]*self.num_freq
        self.CoherYXmax = [0.999]*self.num_freq
        
        self.ExMin = [0]*self.num_freq
        self.ExMax = [1.0e6]*self.num_freq
        
        self.EyMin = [0]*self.num_freq
        self.EyMax = [1.0e6]*self.num_freq
        
        self.HxMin = [0]*self.num_freq
        self.HxMax = [1.0e6]*self.num_freq
        
        self.HyMin = [0]*self.num_freq
        self.HyMax = [1.0e6]*self.num_freq
        
        self.NFCStack = [8]*self.num_freq
        
        self.param_dict = None
        self.param_lst = None
        
        self.string_fmt_lst = ['.4e', '.4e', '.4e', '.1f', '.1f', '.4e', '.4e',
                               '.1f', '.1f', '.3f', '.3f', '.3f', '.3f', '.1g',
                               '.4e', '.1g', '.4e', '.1g', '.4e', '.1g', '.4e',
                               '.0f']
        
    def make_meta_dict(self):
        """
        make meta data dictionary
        """
        if not self.meta_lst:
            self.make_meta_lst()
            
        self.meta_dict = dict([(mkey, mvalue) for mkey, mvalue in 
                                zip(self.meta_keys, self.meta_lst)])
                                
    def make_meta_lst(self):
        """
        make metadata list
        """
        
        self.meta_lst = [self.mtedit_version,
                         self.phase_flip,
                         self.phaseslope_smooth,
                         self.phaseslope_tozmag,
                         self.dplus_use,
                         self.autoskip_ondplus,
                         self.autoskip_dplusdev]
        
    def make_param_dict(self):
        """
        make a parameter dictionary
        """
        if not self.param_lst:
            self.make_param_lst()
            
        self.param_dict = dict([(mkey, mvalue) for mkey, mvalue in
                                 zip(self.param_header, self.param_lst)])
        
    def make_param_lst(self):
        """
        make a list of parameters
        """
        
        self.param_lst = [self.freq_lst,
                          self.AResXYmin, 
                          self.AResXYmax, 
                          self.ZPhzXYmin, 
                          self.ZPhzXYmax,
                          self.AResYXmin, 
                          self.AResYXmax, 
                          self.ZPhzYXmin, 
                          self.ZPhzYXmax,
                          self.CoherXYmin,
                          self.CoherYXmin,
                          self.CoherXYmax,
                          self.CoherYXmax,
                          self.ExMin,
                          self.ExMax,
                          self.EyMin,
                          self.EyMax,
                          self.HxMin,
                          self.HxMax,
                          self.HyMin,
                          self.HyMax,
                          self.NFCStack]
                          

    def read_config(self, cfg_fn):
        """
        read a MTEdit.cfg file
        """

        if not os.path.isfile(cfg_fn):
            raise IOError('{0} does not exist'.format(cfg_fn))
            
        self.cfg_fn = cfg_fn
        self.meta_dict = {}
        self.param_dict = {}
        
        cfid = open(cfg_fn, 'r')
        clines = cfid.readlines()
        for ii, cline in enumerate(clines):
            #--> get metadata 
            if cline[0] == '$':
                clst = cline[1:].strip().split('=')
                self.meta_dict[clst[0]] = clst[1]
            
            #--> get filter parameters header            
            elif cline.find('Frequency') == 0:
                pkeys = [cc.strip() for cc in cline.strip().split(',')]
                nparams = len(clines)-ii
                self.param_dict = dict([(pkey, np.zeros(nparams))
                                           for pkey in pkeys])
                jj = 0
            
            #--> get filter parameters as a function of frequency
            else:
                if len(cline) > 3:
                    clst = [cc.strip() for cc in cline.strip().split(',')]
                    for nn, pkey in enumerate(pkeys):
                        self.param_dict[pkey][jj] = float(clst[nn])
                    jj += 1
        self.num_freq = len(self.param_dict['Frequency'])
                    
    def write_config(self, save_path, mtedit_params_dict=None):
        """
        write a mtedit.cfg file        
        
        """

        if os.path.isdir(save_path) == True:
            save_path = os.path.join(save_path, 'mtedit.cfg')
        
        self.cfg_fn = save_path
        if not self.meta_dict:
            self.make_meta_dict()
            
        if not self.param_dict:
            self.make_param_dict()
        
        #--- write file ---
        cfid = open(self.cfg_fn, 'w')
        
        #--> write metadata
        for mkey in self.meta_keys:
            cfid.write('${0}={1}\n'.format(mkey, self.meta_dict[mkey]))
        
        #--> write parameter header
        for header in self.param_header[:-1]:
            cfid.write('{0:>11},'.format(header))
        cfid.write('{0:>11}\n'.format(self.param_header[-1]))
        
        #--> write filter parameters
        for ii in range(self.num_freq):
            for jj, pkey in enumerate(self.param_header[:-1]):
                cfid.write('{0:>11},'.format('{0:{1}}'.format(
                                self.param_dict[pkey][ii], 
                                self.string_fmt_lst[jj])))
            cfid.write('{0:>11}\n'.format('{0:{1}}'.format(
                                self.param_dict[self.param_header[-1]][ii], 
                                self.string_fmt_lst[-1])))
    
        cfid.close()
        print('Wrote mtedit config file to {0}'.format(self.cfg_fn))
        
    
#==============================================================================
# deal with avg files output from mtedit
#==============================================================================    
class ZongeMTAvg():
    """
    deal with avg files output from mtedit and makes an .edi file.
    
    
    =============================== ===========================================
    Attributes                       Description     
    =============================== ===========================================
     MTEdit3Auto_PhaseFlip          [ yes | no ] flip phase automatically
     MTEdit3DPlus_Use               [ yes | no ] use D+ smoothing
     MTEdit3PhaseSlope_Smooth       [ yes | no ] smooth data using phase
     MTEdit3PhaseSlope_toMag        [ yes | no ] use phase to predict mag
     MTEdit3Version                 version of mtedit
     Rx_GdpStn                      station name
     Rx_HPR                         station rotation (N, E, Z)
     Rx_Length                      dipole lenghts
     Survey_Array                   survey array
     Survey_Type                    survey type (MT)
     Tipper                         mtpy.core.z.Tipper object
     Tx_Type                        Transmitter type
     Unit_Length                    units of length (m) 
     Z                              mtpy.core.z.Z object
     avg_dict                       dictionary of all meta data for MTAvg 
     comp                           components
     comp_dict                      dictionary of components
     comp_flag                      component flag
     comp_index                     index of component
     comp_lst_tip                   list of tipper information
     comp_lst_z                     list of z information
     freq_dict                      dictionary of frequencies
     freq_dict_x                    dictionary of frequencies in x direction
     freq_dict_y                    dictionary of frequencies in y direction
     header_dict                    dictionary of header information
     info_dtype                     numpy.dtype for information 
     info_keys                      keys for information
     info_type                      keys type
     nfreq                          number of frequencies
     nfreq_tipper                   number of frequencies for tipper
     z_coordinate                   coordinate of z
    =============================== ===========================================
    
    =============================== ===========================================
    Methods                         Description
    =============================== ===========================================
    convert2complex                 convert res/phase to Z          
    fill_Tipper                     fill tipper data in to Tipper             
    fill_Z                          fill z data to Z
    read_avg_file                   read in .avg file output by MTEdit 
    write_edi                       write .edi from .avg file   
    =============================== ===========================================
    
    
    :Example: ::
        
        >>> import mtpy.usgs.zonge as zonge
        >>> zm = zonge.ZongeMTAvg(r"/home/mt01/Merged"\
                                  'mt01', \
                                  survey_cfg_file=r"/home/mt/survey.cfg",\
                                  mtft_cfg_file=r"/home/mt/mt01/Merged/mtft24.cfg"\,
                                  mtedit_cfg_file=r"/home/bin/mtedit.cfg",\
                                  copy_path=r"/home/mt/edi_files")
    """                     

    def __init__(self):
        
        
        self.Survey_Type = 'NSAMT'
        self.Survey_Array = 'Tensor'
        self.Tx_Type = 'Natural'
        self.MTEdit3Version = '3.001 applied on 2010-11-19'
        self.MTEdit3Auto_PhaseFlip = 'No'
        self.MTEdit3PhaseSlope_Smooth = 'Moderate'
        self.MTEdit3PhaseSlope_toMag = 'No'
        self.MTEdit3DPlus_Use = 'No'
        self.Rx_GdpStn = 4
        self.Rx_Length = 100
        self.Rx_HPR = [90, 0, 0]
        self.GPS_Lat = 0.0
        self.GPS_Lon = 0.0
        self.Unit_Length = 'm'
        self.header_dict = {'Survey.Type':self.Survey_Type,
                            'Survey.Array':self.Survey_Array,
                            'Tx.Type':self.Tx_Type,
                            'MTEdit:Version':self.MTEdit3Version,
                            'MTEdit:Auto.PhaseFlip':self.MTEdit3Auto_PhaseFlip,
                            'MTEdit:PhaseSlope.Smooth':self.MTEdit3PhaseSlope_Smooth,
                            'MTEdit:PhaseSlope.toZmag':self.MTEdit3PhaseSlope_toMag,
                            'MTEdit:DPlus.Use':self.MTEdit3DPlus_Use,
                            'Rx.GdpStn':self.Rx_GdpStn,
                            'Rx.Length':self.Rx_Length,
                            'Rx.HPR':self.Rx_HPR,
                            'GPS.Lat':self.GPS_Lat,
                            'GPS.Lon':self.GPS_Lon,
                            'Unit.Length':self.Unit_Length}
                            
        self.info_keys = ['Skp', 'Freq', 'E.mag', 'B.mag', 'Z.mag', 'Z.phz',
                          'ARes.mag', 'ARes.%err', 'Z.perr', 'Coher', 
                          'FC.NUse', 'FC.NTry']
        self.info_type = [np.int, np.float, np.float, np.float, np.float, 
                          np.float, np.float, np.float, np.float, np.float, 
                          np.int, np.int]
        self.info_dtype = np.dtype([(kk.lower(), tt) 
                                    for kk, tt in zip(self.info_keys, 
                                                      self.info_type)])
                          
        self.Z = mtz.Z()
        self.Tipper = mtz.Tipper()
        self.comp_lst_z = ['zxx','zxy','zyx','zyy']
        self.comp_lst_tip = ['tzx','tzy']
        self.comp_index = {'zxx':(0,0), 'zxy':(0,1), 'zyx':(1,0), 'zyy':(1,1),
                           'tzx':(0,0), 'tzy':(0,1)}
        self.comp_flag = {'zxx':False, 'zxy':False, 'zyx':False, 'zyy':False,
                          'tzx':False, 'tzy':False}
        self.comp_dict = None
        self.comp = None
        self.nfreq = None
        self.nfreq_tipper = None
        self.freq_dict = None
        self.freq_dict_x = None
        self.freq_dict_y = None
        self.avg_dict = {'ex':'4', 'ey':'5'}
        self.z_coordinate = 'down'

        
    def read_avg_file(self, avg_fn):
        """
        read in average file        
        """
        
        if not os.path.isfile(avg_fn):
            raise IOError('{0} does not exist, check file'.format(avg_fn))
        
        self.comp = os.path.basename(avg_fn)[0]
        with open(avg_fn, 'r') as fid:
            alines = fid.readlines()
        self.comp_flag = {'zxx':False, 'zxy':False, 'zyx':False, 'zyy':False,
                          'tzx':False, 'tzy':False}
        
        if not self.comp_dict:
            # check to see if all 4 components are in the .avg file
            if len(alines) > 140:  
                self.comp_dict = dict([(ckey, np.zeros(int(len(alines)/4), 
                                                       dtype=self.info_dtype))
                                        for ckey in list(self.comp_flag.keys())])
            # if there are only 2
            else:
                self.comp_dict = dict([(ckey, np.zeros(int(len(alines)/2), 
                                                       dtype=self.info_dtype))
                                        for ckey in list(self.comp_flag.keys())])
        self.comp_lst_z = []
        self.comp_lst_tip = []
        ii = 0                        
        for aline in alines:
            if aline.find('=') > 0 and aline.find('$') == 0:
                alst = [aa.strip() for aa in aline.strip().split('=')]
                if alst[1].lower() in list(self.comp_flag.keys()):
                    akey = alst[1].lower()
                    self.comp_flag[akey] = True
                    if akey[0] == 'z':
                        self.comp_lst_z.append(akey)
                    elif akey[0] == 't':
                        self.comp_lst_tip.append(akey)
                    ii = 0
                else:
                    akey = alst[0][1:].replace('.', '_')
                    if akey.lower().find('length'):
                        alst[1] = alst[1][0:-1]
                    try:
                        self.__dict__[akey] = float(alst[1])
                    except ValueError:
                        self.__dict__[akey] = alst[1]
                    #self.header_dict[alst[0][1:]] = al            print(aline)st[1]
            elif aline[0] == 'S':
                pass
            # read the data line.
            elif len(aline) > 2:
                aline = aline.replace('*', '0.50')
                alst = [aa.strip() for aa in aline.strip().split(',')]
                for cc, ckey in enumerate(self.info_keys):
                    self.comp_dict[akey][ii][ckey.lower()] = alst[cc]
                ii += 1
     
        self.fill_Z()
        self.fill_Tipper()
        
        print('Read file {0}'.format(avg_fn))
        
    def convert2complex(self, zmag, zphase):
        """
        outputs of mtedit are magnitude and phase of z, convert to real and
        imaginary parts, phase is in milliradians
        
        """
        
        if type(zmag) is np.ndarray:
            assert len(zmag) == len(zphase)
        
        if self.z_coordinate == 'up':
            zreal = zmag*np.cos((zphase/1000)%np.pi)
            zimag = zmag*np.sin((zphase/1000)%np.pi)
        else:
            zreal = zmag*np.cos((zphase/1000))
            zimag = zmag*np.sin((zphase/1000))
        
        return zreal, zimag
        
    def _match_freq(self, freq_list1, freq_list2):
        """
        fill the frequency dictionary where keys are freqeuency and 
        values are index of where that frequency should be in the array of z
        and tipper
        """
#        
#        if set(freq_list1).issubset(freq_list2) == True:
#            return dict([(freq, ff) for ff, freq in enumerate(freq_list1)])
#        else:
        comb_freq_list = list(set(freq_list1).intersection(freq_list2))+\
                      list(set(freq_list1).symmetric_difference(freq_list2))
        comb_freq_list.sort()
        return dict([(freq, ff) for ff, freq in enumerate(comb_freq_list)])
        
    def fill_Z(self):
        """
        create Z array with data
        """
        flst = np.array([len(np.nonzero(self.comp_dict[comp]['freq'])[0])
                         for comp in self.comp_lst_z])
        
        nz = flst.max()
        freq = self.comp_dict[self.comp_lst_z[np.where(flst==nz)[0][0]]]['freq']
        freq = freq[np.nonzero(freq)]

        if self.nfreq:
            self.freq_dict_y = dict([(ff, nn) for nn, ff in enumerate(freq)])
            #get new frequency dictionary to match index values
            new_freq_dict = self._match_freq(sorted(self.freq_dict_x.keys()), 
                                             freq)
            
            new_nz = len(list(new_freq_dict.keys()))
            self.freq_dict = new_freq_dict
            #fill z according to index values
            new_Z = mtz.Z()
            new_Z.freq = sorted(new_freq_dict.keys())
            new_Z.z = np.zeros((new_nz, 2, 2), dtype='complex')
            new_Z.z_err = np.ones((new_nz, 2, 2))
            nzx, nzy, nzz = self.Z.z.shape
            
            
            
            #need to fill the new array with the old values, but they
            # need to be stored in the correct position
            clst = ['zxx', 'zxy', 'zyx', 'zyy']
            for cc in self.comp_lst_z:
                clst.remove(cc)
            for ikey in clst:
                for kk, zz in enumerate(self.Z.z):
                    ii, jj = self.comp_index[ikey]
                    if zz[ii, jj].real != 0.0:
                        #index for new Z array
                        ll = self.freq_dict[self.comp_dict[ikey]['freq'][kk]]
                        
                        #index for old Z array
                        try:
                            mm = self.freq_dict_x[self.comp_dict[ikey]['freq'][kk]]

                            new_Z.z[ll] = self.Z.z[mm]
                            new_Z.z_err[ll] = self.Z.z_err[mm]
                        except KeyError:
                            pass
                        
            #fill z with values from comp_dict
            for ikey in self.comp_lst_z:
                ii, jj = self.comp_index[ikey]

                zr, zi = self.convert2complex(self.comp_dict[ikey]['z.mag'][:nz].copy(),
                                              self.comp_dict[ikey]['z.phz'][:nz].copy())
                for kk, zzr, zzi in zip(list(range(len(zr))), zr, zi):
                    ll = self.freq_dict[self.comp_dict[ikey]['freq'][kk]]
                    if ikey.find('yx') > 0 and self.z_coordinate == 'up':
                        new_Z.z[ll, ii, jj] = -1*(zzr+zzi*1j)
                    else:
                        new_Z.z[ll, ii, jj] = zzr+zzi*1j
                    new_Z.z_err[ll,ii, jj] = \
                                self.comp_dict[ikey]['ares.%err'][kk]*.005
                
            
            self.Z = new_Z
        
        #fill for the first time
        else:
            self.nfreq = nz
            self.freq_dict_x = dict([(ff, nn) for nn, ff in enumerate(freq)])
            #fill z with values
            z = np.zeros((nz, 2, 2), dtype='complex')
            z_err = np.ones((nz, 2, 2))
            
            for ikey in self.comp_lst_z:
                ii, jj = self.comp_index[ikey]
                    
                zr, zi = self.convert2complex(self.comp_dict[ikey]['z.mag'][:nz].copy(),
                                              self.comp_dict[ikey]['z.phz'][:nz].copy())
                
                if ikey.find('yx') > 0 and self.z_coordinate == 'up':
                    z[:, ii, jj] = -1*(zr+zi*1j)
                else:
                    z[:, ii, jj] = zr+zi*1j

                z_err[:,ii, jj] = self.comp_dict[ikey]['ares.%err'][:nz]*.005 

            self.Z.freq = freq
            self.Z.z = z
            self.Z.z_err = z_err
            
            
        self.Z.z = np.nan_to_num(self.Z.z)
        self.Z.z_err = np.nan_to_num(self.Z.z_err)
                
                
    def fill_Tipper(self):
        """
        fill tipper values
        """
        
        if self.comp_flag['tzy'] == False and self.comp_flag['tzx'] == False:
            print('No Tipper found')
            return
            
        flst = np.array([len(np.nonzero(self.comp_dict[comp]['freq'])[0])
                         for comp in self.comp_lst_tip])
        nz = flst.max()
        freq = self.comp_dict[self.comp_lst_tip[np.where(flst==nz)[0][0]]]['freq']
        freq = freq[np.nonzero(freq)]
        if self.nfreq_tipper and self.Tipper.tipper is not None:
            #get new frequency dictionary to match index values
            new_freq_dict = self._match_freq(sorted(self.freq_dict.keys()), 
                                             freq)
            
            new_nz = len(list(new_freq_dict.keys()))
            #fill z according to index values
            new_Tipper = mtz.Tipper()
            new_Tipper.tipper = np.zeros((new_nz, 1, 2), dtype='complex')
            new_Tipper.tipper_err = np.ones((new_nz, 1, 2))
            
            self.freq_dict = new_freq_dict
            
            #need to fill the new array with the old values, but they
            # need to be stored in the correct position
            for ikey in ['tzx', 'tzy']:
                for kk, tt in enumerate(self.Tipper.tipper):
                    ii, jj = self.comp_index[ikey]
                    if tt[ii, jj].real != 0.0:
                        #index for new tipper array
                        ll = self.freq_dict[self.comp_dict[ikey]['freq'][kk]]
                        
                        #index for old tipper array
                        try:
                            mm = self.freq_dict_x[self.comp_dict[ikey]['freq'][kk]]

                            new_Tipper.tipper[ll] = self.Tipper.tipper[mm]
                            new_Tipper.tipper_err[ll] = self.Tipper.tipper_err[mm]
                        except KeyError:
                            pass

                        
            #fill z with values from comp_dict
            for ikey in self.comp_lst_tip:
                ii, jj = self.comp_index[ikey]

                tr, ti = self.convert2complex(self.comp_dict[ikey]['z.mag'][:nz],
                                              self.comp_dict[ikey]['z.phz'][:nz])
                for kk, tzr, tzi in zip(list(range(len(tr))), tr, ti):
                    ll = self.freq_dict[self.comp_dict[ikey]['freq'][kk]]
                    
                    if self.z_coordinate == 'up':
                        new_Tipper.tipper[ll, ii, jj] = -1*(tzr+tzi*1j)
                    else:
                        new_Tipper.tipper[ll, ii, jj] = tzr+tzi*1j
                    #error estimation
                    new_Tipper.tipper_err[ll,ii, jj] += \
                                self.comp_dict[ikey]['ares.%err'][kk]*\
                                                .05*np.sqrt(tzr**2+tzi**2)
                
            new_Tipper.freq = sorted(self.freq_dict.keys())
            self.Tipper = new_Tipper
           
        else:
            self.nfreq_tipper = nz
            self.freq_dict_x = dict([(ff, nn) for nn, ff in enumerate(freq)])
            #fill z with values
            tipper = np.zeros((nz, 1, 2), dtype='complex')
            tipper_err = np.ones((nz, 1, 2))
            
            for ikey in self.comp_lst_tip:
                ii, jj = self.comp_index[ikey]
                    
                tzr, tzi = self.convert2complex(self.comp_dict[ikey]['z.mag'][:nz],
                                              self.comp_dict[ikey]['z.phz'][:nz])
                
                if self.z_coordinate == 'up':
                    tipper[:, ii, jj] = -1*(tzr+tzi*1j)
                else:
                    tipper[:, ii, jj] = tzr+tzi*1j
                tipper_err[:, ii, jj] = self.comp_dict[ikey]['ares.%err'][:nz]*\
                                                     .05*np.sqrt(tzr**2+tzi**2)
                    
            self.Tipper.freq = sorted(self.freq_dict_x.keys())
            self.Tipper.tipper = tipper
            self.Tipper.tipper_err = tipper_err
            
            
        self.Tipper.tipper = np.nan_to_num(self.Tipper.tipper)
        self.Tipper.tipper_err = np.nan_to_num(self.Tipper.tipper_err)
        
    def write_edi(self, avg_fn, station, survey_dict=None, 
                  survey_cfg_file=None,  mtft_cfg_file=None, 
                  mtedit_cfg_file=r"c:\MinGW32-xy\Peacock\zen\bin\mtedit.cfg", 
                  save_path=None, rrstation=None, 
                  copy_path=r"d:\Peacock\MTData\EDI_Files", avg_ext='.avg'):
        """
        write an edi file from the .avg files
        
        Arguments:
        ----------
            **avg_fn** : string
                         full path to avg file name
            
            **survey_dict** : dictionary
                              dictionary containing the survey parameters
                              such as lat, lon, elevation, date, etc.
                              
            **survey_cfg_file** : string (full path to survey file)
                              file contains all the important information 
                              about the setup of the station, input file if
                              survey_dict is None.  This is created by 
                              mtpy.configfile
                              
            **mtft_cfg_file** : string (full path to mtft24.cfg file)
                               this file contains information on how the
                               Fourier coefficients were calculated
                               
            **mtedit_cfg_file** : string (full path to MTEdit.cfg file)
                                  this file contains information on how 
                                  the transfer functions were estimated
            
            **save_path** : string (full path or directory to where .edi file 
                                    will be saved)
                                    
        Outputs:
        ---------
            **edi_fn** : string (full path to .edi file)
                      
                      
                      
        """
        
        if save_path is None:
            save_dir = os.path.dirname(avg_fn)
            save_path = os.path.join(save_dir, station+'.edi')
        
        #create an mtedi instance
        self.edi = mtedi.Edi()
        self.edi.Z = self.Z
        self.edi.Tipper = self.Tipper

        #read in avg file
        if os.path.isfile(avg_fn) == True:
            self.read_avg_file(avg_fn)
            self.edi.Z = self.Z
            self.edi.Tipper = self.Tipper
        else:
            raise NameError('Could not find {0}'.format(avg_fn))
        
        #read in survey file
        survey_dict = None
        if survey_cfg_file is not None:
            sdict = mtcf.read_survey_configfile(survey_cfg_file)
            
        try:
            survey_dict = sdict[station.upper()]
        except KeyError:
            if survey_dict is not None:
                try:
                    survey_dict['station']
                except KeyError:
                    try:
                        survey_dict['station_name']
                    except KeyError:
                        raise KeyError('Could not find station information in'
                                       ', check inputs')
            else:
                raise KeyError('Could not find {0} in survey file'.format(
                                station.upper()))
                                 
        #get remote reference information if desired
        if rrstation:
            try:
                rrsurvey_dict = sdict[rrstation.upper()]
                survey_dict['rr_station'] = rrsurvey_dict['station']
                survey_dict['rr_station_elevation'] = rrsurvey_dict['elevation']
                survey_dict['rr_station_latitude'] = gis_tools.assert_lat_value(
                                               rrsurvey_dict.pop('latitude',0.0))
                survey_dict['rr_station_longitude'] = gis_tools.assert_lon_value(
                                               rrsurvey_dict.pop('longitude',0.0))
            except KeyError:
                print('Could not find station information for remote reference')
        else:
            rrsurvey_dict = None
            
        #read in mtft24.cfg file
        if mtft_cfg_file is None:
            try:
                mtft_cfg_file = os.path.join(save_dir, 'mtft24.cfg')
                zmtft = ZongeMTFT()
                zmtft.read_cfg(mtft_cfg_file)
                mtft_dict = zmtft.meta_dict
            except:
                mtft_dict = None
        else:
            zmtft = ZongeMTFT()
            zmtft.read_cfg(mtft_cfg_file)
            mtft_dict = zmtft.meta_dict
            
        #read in mtedit.cfg file
        if mtedit_cfg_file:
            zmtedit = ZongeMTEdit()
            zmtedit.read_config(mtedit_cfg_file)
            mtedit_dict = zmtedit.meta_dict
        else:
            mtedit_dict = None
            
        #----------------HEAD BLOCK------------------
        #from survey dict get information

        #--> data id
        try:
            self.edi.Header.dataid = survey_dict['station']
        except KeyError:
            self.edi.Header.dataid = station
            
        #--> acquired by
        self.edi.Header.acqby = survey_dict.pop('network','USGS')
        
        #--> file by
        self.edi.Header.fileby = survey_dict.pop('network','MTpy')
        
        #--> acquired date
        self.edi.Header.acqdate = survey_dict.pop('date', 
                                    time.strftime('%Y-%m-%d',time.localtime()))
        
        #--> prospect
        self.edi.Header.loc = survey_dict.pop('location', 'Earth')
        
        #--> latitude
        self.edi.Header.lat = survey_dict.pop('latitude',0.0)
        
        #--> longitude
        self.edi.Header.lon = survey_dict.pop('longitude',0.0)
        
        #--> elevation
        self.edi.Header.elev = survey_dict.pop('elevation', 0)
       
        #-----------------INFO BLOCK---------------------------
        self.edi.Info.info_list = []
        self.edi.Info.info_list.append('MAX LINES: 999')
        
        #--> put the rest of the survey parameters in the info block
        for skey in sorted(survey_dict.keys()):
            self.edi.Info.info_list.append('{0}: {1}'.format(skey, 
                                                           survey_dict[skey]))
        
        #--> put parameters about how fourier coefficients were found
        if mtft_dict is not None:
            for mkey in sorted(mtft_dict.keys()):
                if mkey == 'setup_lst' or \
                   mkey.lower() == 'mtft.tsplot.chnrange':
                    pass
                else:
                    self.edi.Info.info_list.append('{0}: {1}'.format(mkey,
                                                           mtft_dict[mkey]))
        
        #--> put parameters about how transfer function was found
        if mtedit_dict is not None:
            for mkey in list(mtedit_dict.keys()):
                self.edi.Info.info_list.append('{0}: {1}'.format(mkey,
                                                           mtedit_dict[mkey]))
        
        #----------------DEFINE MEASUREMENT BLOCK------------------
        self.edi.Define_measurement.maxchan = 5
        self.edi.Define_measurement.maxrun = 999
        self.edi.Define_measurement.maxmeas = 99999

        try:
            self.edi.Define_measurement.units = mtedit_dict['unit.length']
        except (TypeError, KeyError):
            self.edi.Define_measurement.units = 'm'
            
        self.edi.Define_measurement.reftype = 'cartesian'
        self.edi.Define_measurement.reflat = self.edi.Header.lat
        self.edi.Define_measurement.reflon = self.edi.Header.lon
        self.edi.Define_measurement.refelev = self.edi.Header.elev

        
        #------------------HMEAS_EMEAS BLOCK--------------------------
        if mtft_dict:
            chn_lst = mtft_dict['setup_lst'][0]['Chn.Cmp']
            chn_id = mtft_dict['setup_lst'][0]['Chn.ID']
            chn_len_lst = mtft_dict['setup_lst'][0]['Chn.Length']
            
        else:
            chn_lst = ['hx', 'hy', 'hz', 'ex', 'ey']
            chn_id = [1, 2, 3, 4, 5]
            chn_len_lst = [100]*5
            
        chn_id_dict = dict([(comp.lower(), (comp.lower(), cid, clen)) 
                            for comp, cid, clen in zip(chn_lst, chn_id, 
                                                       chn_len_lst)])
                                                       
        
        #--> hx component                
        try:
            hxazm = survey_dict['b_xaxis_azimuth']
        except KeyError:
            hxazm = 0
        try:
            hdict = {'id': chn_id_dict['hx'][1], 
                     'chtype': '{0}'.format(chn_id_dict['hx'][0].upper()), 
                     'x':0, 
                     'y':0, 
                     'azm':hxazm,
                     'acqchan':'{0}'.format(chn_id_dict['hx'][0].upper())}
        except KeyError:
            hdict = {'id': 1, 
                     'chtype': '{0}'.format('hx'), 
                     'x':0, 
                     'y':0, 
                     'azm':hxazm,
                     'acqchan':'hx'}
        self.edi.Define_measurement.meas_hx = mtedi.HMeasurement(**hdict)
            
        #--> hy component
        try:
            hyazm = survey_dict['b_yaxis_azimuth']
        except KeyError:
            hyazm = 90
        try:
            hdict = {'id': chn_id_dict['hy'][1], 
                     'chtype': '{0}'.format(chn_id_dict['hy'][0].upper()), 
                     'x':0, 
                     'y':0, 
                     'azm':hyazm,
                     'acqchan':'{0}'.format(chn_id_dict['hy'][0].upper())}

        except KeyError:
            hdict = {'id': 2, 
                     'chtype': 'hy', 
                     'x':0, 
                     'y':0, 
                     'azm':hyazm,
                     'acqchan':'hy'}
        self.edi.Define_measurement.meas_hy = mtedi.HMeasurement(**hdict)
        
        #--> hz component
        try:
            hdict = {'id': chn_id_dict['hz'][1], 
                     'chtype': '{0}'.format(chn_id_dict['hz'][0].upper()), 
                     'x':0, 
                     'y':0, 
                     'azm':0,
                     'acqchan':'{0}'.format(chn_id_dict['hz'][0].upper())}

        except KeyError:
            hdict = {'id': 3, 
                     'chtype': 'hz', 
                     'x':0, 
                     'y':0, 
                     'azm':0}
        self.edi.Define_measurement.meas_hz = mtedi.HMeasurement(**hdict)
        
        #--> ex component
        try:
            edict = {'id':chn_id_dict['ex'][1], 
                     'chtype':'{0}'.format(chn_id_dict['ex'][0].upper()), 
                     'x':0, 
                     'y':0, 
                     'x2':chn_id_dict['ex'][2],
                     'y2':0}
        except KeyError:
            edict = {'id':4, 
                     'chtype':'ex', 
                     'x':0, 
                     'Y':0, 
                     'x2':100,
                     'y2':0}
        self.edi.Define_measurement.meas_ex = mtedi.EMeasurement(**edict)
                           
        #--> ey component
        try:
            edict = {'id':chn_id_dict['ey'][1], 
                     'chtype':'{0}'.format(chn_id_dict['ey'][0].upper()), 
                     'x':0, 
                     'y':0, 
                     'x2':0,
                     'y2':chn_id_dict['ey'][2]}
        except KeyError:
            edict = {'id':5, 
                     'chtype':'ey', 
                     'x':0, 
                     'Y':0, 
                     'x2':0,
                     'y2':100}
        self.edi.Define_measurement.meas_ey = mtedi.EMeasurement(**edict)
                           
        #--> remote reference 
        if rrsurvey_dict:
            hxid = rrsurvey_dict.pop('hx', 6)
            hyid = rrsurvey_dict.pop('hy', 7)
            hxazm = rrsurvey_dict.pop('b_xaxis_azimuth', 0)
            hyazm = rrsurvey_dict.pop('b_xaxis_azimuth', 90)
        else:
            hxid =  6
            hyid =  7
            hxazm = 0
            hyazm = 90
                
        #--> rhx component
        hdict = {'id': hxid, 
                 'chtype': 'rhx', 
                 'x':0, 
                 'y':0, 
                 'azm':hxazm,
                 'acqchan':'rhx'}
        self.edi.Define_measurement.meas_rhx = mtedi.HMeasurement(**hdict)

        #--> rhy component
        hdict = {'id': hyid, 
                 'chtype': 'rhy', 
                 'x':0, 
                 'y':0, 
                 'azm':hyazm,
                 'acqchan':'rhy'}
        self.edi.Define_measurement.meas_rhy = mtedi.HMeasurement(**hdict)
        
        #----------------------MTSECT-----------------------------------------
        self.edi.Data_sect.nfreq = len(self.Z.freq)
        self.edi.Data_sect.sectid = station        
        self.edi.Data_sect.nchan = len(chn_lst)
        for chn, chnid in zip(chn_lst, chn_id):
            setattr(self.edi.Data_sect, chn, chnid)
        
        #----------------------ZROT BLOCK--------------------------------------
        self.edi.zrot = np.zeros(len(self.edi.Z.z))
        
        #----------------------FREQUENCY BLOCK---------------------------------
        self.edi.freq = self.Z.freq
        
            
        #============ WRITE EDI FILE ==========================================
        edi_fn = self.edi.write_edi_file(new_edi_fn=save_path)
        
        print('Wrote .edi file to {0}'.format(edi_fn))
        
        if copy_path is not None:
            copy_edi_fn = os.path.join(copy_path, os.path.basename(edi_fn))
            if not os.path.exists(copy_path):
                os.mkdir(copy_path)
            shutil.copy(edi_fn, copy_edi_fn)
            print('Copied {0} to {1}'.format(edi_fn, copy_edi_fn))
        
        return edi_fn
        
    def plot_mt_response(self, avg_fn, **kwargs):
        """
        plot an mtv file
        """
        
        if os.path.isfile(avg_fn) is False:
            raise IOError('Could not find {0}, check path'.format(avg_fn))
        
        self.read_avg_file(avg_fn)
        
        plot_resp = plotresponse.PlotResponse(z_object=self.Z,
                                              tipper_object=self.Tipper, 
                                              plot_tipper='yri',
                                              **kwargs)
        
        return plot_resp
        
        
    def write_edi_from_avg(self, avg_fn, station, survey_dict=None, 
                  survey_cfg_file=None,  mtft_cfg_file=None, 
                  mtedit_cfg_file=r"c:\MinGW32-xy\Peacock\zen\bin\mtedit.cfg", 
                  save_path=None, rrstation=None, 
                  copy_path=r"d:\Peacock\MTData\EDI_Files", avg_ext='.avg'):
        """
        write an edi file from the .avg files
        
        Arguments:
        ----------
            **fnx** : string (full path to electric north file)
                      file for Zxx, Zxy
                      
            **fny** : string (full path to electric east file)
                      file for Zyx, Zyy
            
            **survey_dict** : dictionary
                              dictionary containing the survey parameters
                              such as lat, lon, elevation, date, etc.
                              
            **survey_cfg_file** : string (full path to survey file)
                              file contains all the important information 
                              about the setup of the station, input file if
                              survey_dict is None.  This is created by 
                              mtpy.configfile
                              
            **mtft_cfg_file** : string (full path to mtft24.cfg file)
                               this file contains information on how the
                               Fourier coefficients were calculated
                               
            **mtedit_cfg_file** : string (full path to MTEdit.cfg file)
                                  this file contains information on how 
                                  the transfer functions were estimated
            
            **save_path** : string (full path or directory to where .edi file 
                                    will be saved)
                                    
        Outputs:
        ---------
            **edi_fn** : string (full path to .edi file)
                      
                      
                      
        """
        
        if save_path is None:
            save_dir = os.path.dirname(avg_fn)
            save_path = os.path.join(save_dir, station+'.edi')
        
        #create an mtedi instance
        self.edi = mtedi.Edi()
        self.edi.Z = self.Z
        self.edi.Tipper = self.Tipper
        
        if os.path.isfile(avg_fn) == True:
            self.read_avg_file(avg_fn)
            self.edi.Z = self.Z
            self.edi.Tipper = self.Tipper

        #read in survey file
        survey_dict = {}
        survey_dict['latitude'] = gis_tools.assert_lat_value(self.GPS_Lat)
        survey_dict['longitude'] = gis_tools.assert_lon_value( self.GPS_Lon)
        survey_dict['elevation'] = gis_tools.assert_elevation_value(self.Rx_Length)
        survey_dict['station'] = station
        if survey_cfg_file is not None:
            sdict = mtcf.read_survey_configfile(survey_cfg_file)
            
            try:
                survey_dict = sdict[station.upper()]
            except KeyError:
                if survey_dict is not None:
                    try:
                        survey_dict['station']
                    except KeyError:
                        try:
                            survey_dict['station_name']
                        except KeyError:
                            print('Could not find station information in'
                                           ', check inputs')
                                 
        #get remote reference information if desired
        if rrstation:
            try:
                rrsurvey_dict = sdict[rrstation.upper()]
                survey_dict['rr_station'] = rrsurvey_dict['station']
                survey_dict['rr_station_elevation'] = rrsurvey_dict['elevation']
                survey_dict['rr_station_latitude'] = gis_tools.assert_lat_value(
                                               rrsurvey_dict.pop('latitude',0.0))
                survey_dict['rr_station_longitude'] = gis_tools.assert_lon_value(
                                               rrsurvey_dict.pop('longitude',0.0))
            except KeyError:
                print('Could not find station information for remote reference')
        else:
            rrsurvey_dict = None
            
        #read in mtft24.cfg file
        if mtft_cfg_file is None:
            try:
                mtft_cfg_file = os.path.join(save_dir, 'mtft24.cfg')
                zmtft = ZongeMTFT()
                zmtft.read_cfg(mtft_cfg_file)
                mtft_dict = zmtft.meta_dict
            except:
                mtft_dict = None
        else:
            zmtft = ZongeMTFT()
            zmtft.read_cfg(mtft_cfg_file)
            mtft_dict = zmtft.meta_dict
            
        #read in mtedit.cfg file
        if mtedit_cfg_file:
            zmtedit = ZongeMTEdit()
            zmtedit.read_config(mtedit_cfg_file)
            mtedit_dict = zmtedit.meta_dict
        else:
            mtedit_dict = None
            
        #----------------HEAD BLOCK------------------
        #from survey dict get information
        head_dict = {}

        #--> data id
        try:
            head_dict['dataid'] = survey_dict['station']
        except KeyError:
            head_dict['dataid'] = station
            
        #--> acquired by
        head_dict['acqby'] = survey_dict.pop('network','')
        
        #--> file by
        head_dict['fileby'] = survey_dict.pop('network','')
        
        #--> acquired date
        head_dict['acqdate'] = survey_dict.pop('date', 
                                    time.strftime('%Y-%m-%d',time.localtime()))
        
        #--> prospect
        head_dict['loc'] = survey_dict.pop('location', '')
        
        #--> latitude
        head_dict['lat'] = gis_tools.assert_lat_value(survey_dict.pop('latitude',
                                                                      0.0))
        
        #--> longitude
        head_dict['long'] = gis_tools.assert_lon_value(survey_dict.pop('longitude',
                                                                       0.0))
        
        #--> elevation
        head_dict['elev'] = survey_dict.pop('elevation', 0.0)
        
        #--> set header dict as attribute of edi
        self.edi.head = head_dict
       
       #-----------------INFO BLOCK---------------------------
        info_dict = {}
        info_dict['max lines'] = 1000
        
        #--> put the rest of the survey parameters in the info block
        for skey in list(survey_dict.keys()):
            info_dict[skey] = survey_dict[skey]
        
        #--> put parameters about how fourier coefficients were found
        if mtft_dict:
            for mkey in list(mtft_dict.keys()):
                if mkey == 'setup_lst' or \
                   mkey.lower() == 'mtft.tsplot.chnrange':
                    pass
                else:
                    info_dict[mkey] = mtft_dict[mkey]
        
        #--> put parameters about how transfer function was found
        if mtedit_dict:
            for mkey in list(mtedit_dict.keys()):
                info_dict[mkey] = mtedit_dict[mkey]
        
        #--> set info dict as attribute of edi
        self.edi.info_dict = info_dict
                
        #----------------DEFINE MEASUREMENT BLOCK------------------
        definemeas_dict = {}
        
        definemeas_dict['maxchan'] = 5
        definemeas_dict['maxrun'] = 999
        definemeas_dict['maxmeas'] = 99999
        try:
            definemeas_dict['units'] = mtedit_dict['unit.length']
        except (TypeError, KeyError):
            definemeas_dict['units'] = 'm'
        definemeas_dict['reftypy'] = 'cartesian'
        definemeas_dict['reflat'] = head_dict['lat']
        definemeas_dict['reflon'] = head_dict['long']
        definemeas_dict['refelev'] = head_dict['elev']
        
        #--> set definemeas as attribure of edi
        self.edi.definemeas = definemeas_dict
        
        #------------------HMEAS_EMEAS BLOCK--------------------------
        hemeas_lst = []
        if mtft_dict:
            chn_lst = mtft_dict['setup_lst'][0]['Chn.Cmp']
            chn_id = mtft_dict['setup_lst'][0]['Chn.ID']
            chn_len_lst = mtft_dict['setup_lst'][0]['Chn.Length']
            
        else:
            chn_lst = ['hx', 'hy', 'hz', 'ex', 'ey']
            chn_id = [1, 2, 3, 4, 5]
            chn_len_lst = [100]*5
            
        chn_id_dict = dict([(comp.lower(), (comp.lower(), cid, clen)) 
                            for comp, cid, clen in zip(chn_lst, chn_id, 
                                                       chn_len_lst)])
        
        #--> hx component                
        try:
            hxazm = survey_dict['b_xaxis_azimuth']
        except KeyError:
            hxazm = 0
        try:
            hemeas_lst.append(['HMEAS', 
                               'ID={0}'.format(chn_id_dict['hx'][1]), 
                               'CHTYPE={0}'.format(chn_id_dict['hx'][0].upper()), 
                               'X=0', 
                               'Y=0', 
                               'AZM={0}'.format(hxazm),
                               ''])
        except KeyError:
            hemeas_lst.append(['HMEAS', 
                               'ID={0}'.format(1), 
                               'CHTYPE={0}'.format('HX'), 
                               'X=0', 
                               'Y=0', 
                               'AZM={0}'.format(hxazm),
                               ''])
            
        #--> hy component
        try:
            hyazm = survey_dict['b_yaxis_azimuth']
        except KeyError:
            hyazm = 90
        try:
            hemeas_lst.append(['HMEAS', 
                               'ID={0}'.format(chn_id_dict['hy'][1]), 
                               'CHTYPE={0}'.format(chn_id_dict['hy'][0].upper()), 
                               'X=0', 
                               'Y=0', 
                               'AZM={0}'.format(hxazm),
                               ''])
        except KeyError:
            hemeas_lst.append(['HMEAS', 
                               'ID={0}'.format(1), 
                               'CHTYPE={0}'.format('HY'), 
                               'X=0', 
                               'Y=0', 
                               'AZM={0}'.format(hxazm),
                               ''])
        #--> ex component
        try:
            hemeas_lst.append(['EMEAS', 
                               'ID={0}'.format(chn_id_dict['ex'][1]), 
                               'CHTYPE={0}'.format(chn_id_dict['ex'][0].upper()), 
                               'X=0', 
                               'Y=0', 
                               'X2={0}'.format(chn_id_dict['ex'][2]),
                               'Y2=0'])
        except KeyError:
            hemeas_lst.append(['EMEAS', 
                               'ID={0}'.format(1), 
                               'CHTYPE={0}'.format('EX'), 
                               'X=0', 
                               'Y=0', 
                               'X2={0}'.format(100),
                               'Y2=0'])
                           
        #--> ey component
        try:
            hemeas_lst.append(['EMEAS', 
                               'ID={0}'.format(chn_id_dict['ey'][1]), 
                               'CHTYPE={0}'.format(chn_id_dict['ey'][0].upper()), 
                               'X=0', 
                               'Y=0', 
                               'X2=0',
                               'Y2={0}'.format(chn_id_dict['ey'][2])])
        except KeyError:
            hemeas_lst.append(['EMEAS', 
                               'ID={0}'.format(1), 
                               'CHTYPE={0}'.format('EY'), 
                               'X=0', 
                               'Y=0', 
                               'X2=0',
                               'Y2={0}'.format(100)])
                           
        #--> remote reference 
        if rrsurvey_dict:
            hxid = rrsurvey_dict.pop('hx', 6)
            hyid = rrsurvey_dict.pop('hy', 7)
            hxazm = rrsurvey_dict.pop('b_xaxis_azimuth', 0)
            hyazm = rrsurvey_dict.pop('b_xaxis_azimuth', 90)
        else:
            hxid =  6
            hyid =  7
            hxazm = 0
            hyazm = 90
                
        #--> rhx component
        hemeas_lst.append(['HMEAS', 
                           'ID={0}'.format(hxid), 
                           'CHTYPE={0}'.format('rhx'.upper()), 
                           'X=0', 
                           'Y=0', 
                           'AZM={0}'.format(hxazm),
                           ''])
        #--> rhy component
        hemeas_lst.append(['HMEAS', 
                           'ID={0}'.format(hyid), 
                           'CHTYPE={0}'.format('rhy'.upper()), 
                           'X=0', 
                           'Y=0', 
                           'AZM={0}'.format(hyazm),
                           ''])
        hmstring_lst = []
        for hm in hemeas_lst:
            hmstring_lst.append(' '.join(hm))
        #--> set hemeas as attribute of edi
        self.edi.hmeas_emeas = hmstring_lst
        
        #----------------------MTSECT-----------------------------------------
        mtsect_dict = {}
        mtsect_dict['sectid'] = station
        mtsect_dict['nfreq'] = len(self.Z.freq)
        for chn, chnid in zip(chn_lst, chn_id):
            mtsect_dict[chn] = chnid
        
        #--> set mtsect as attribure of edi
        self.edi.mtsect = mtsect_dict
        
        #----------------------ZROT BLOCK--------------------------------------
        self.edi.zrot = np.zeros(len(self.edi.Z.z))
        
        #----------------------FREQUENCY BLOCK---------------------------------
        self.edi.freq = self.Z.freq
        
            
        #============ WRITE EDI FILE ==========================================
        edi_fn = self.edi.write_edi_file(save_path)
        
        print('Wrote .edi file to {0}'.format(edi_fn))
        
        if copy_path is not None:
            copy_edi_fn = os.path.join(copy_path, os.path.basename(edi_fn))
            if not os.path.exists(copy_path):
                os.mkdir(copy_path)
            shutil.copy(edi_fn, copy_edi_fn)
            print('Copied {0} to {1}'.format(edi_fn, copy_edi_fn))
        
        return edi_fn
        
            
            
        
        
            
            
        
        
        
        

        
                    
        
        
        
        
        
        
        
    

    
    
    

    
    
