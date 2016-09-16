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
import mtpy.core.edi as mtedi

#==============================================================================


#==============================================================================
datetime_fmt = '%Y-%m-%d,%H:%M:%S'
datetime_sec = '%Y-%m-%d %H:%M:%S'
#==============================================================================

#==============================================================================
# make a class to deal with birrp inputs
#==============================================================================
class BIRRP_processing(object):
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
        self._max_nread = 16000000 
        self.deltat = 256
        
        for key in kwargs:
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
        
    def get_processing_dict(self, fn_birrp_list, hx=2284, hy=2284, hz=2284):
        """
        from fn_birrp_arr make a processing dictionary to input into writing
        a birrp script file
        
        fn_birrp_list = fn_birrp_arr[df]
        """
        comp_dict = {'ex':0, 'ey':1, 'hz':2, 'hx':3, 'hy':4, 
                     'rrhx':0, 'rrhy':1}
        
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
            print st_time, nread_ii, self.nfft*2
                           
            if nread_ii > self.nfft*2:
                # need to add a 1 in to skip the header line            
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
                      time.mktime(start_dt, datetime_fmt), self.deltat))

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
        
