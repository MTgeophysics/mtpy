# -*- coding: utf-8 -*-
"""
Created on Thu Jul 06 14:24:18 2017

@author: jpeacock
"""

#==============================================================================
# Imports
#==============================================================================
import os
import time

import tables
import mtpy.usgs.zen as zen

#==============================================================================

fn = r"d:\Peacock\MTData\Umatilla\hf05\hf05_20170517_193018_256_EX.Z3D" 
#==============================================================================
class MT_TS(object):
    """
    MT time series object that will read/write data in different formats
    including hdf5, txt, miniseed.
    
    Metadata
    -----------
    
        ==================== ==================================================
        Name                 Description        
        ==================== ==================================================
        azimuth              clockwise angle from coordinate system N (deg)
        calibration_fn       file name for calibration data
        component            component name [ 'ex' | 'ey' | 'hx' | 'hy' | 'hz']
        coordinate_system    [ geographic | geomagnetic ]
        datum                datum of geographic location ex. WGS84
        declination          geomagnetic declination (deg)
        dipole_length        length of dipole (m)
        data_logger          data logger type
        instrument_id        ID number of instrument for calibration
        lat                  latitude of station in decimal degrees
        lon                  longitude of station in decimal degrees
        n_samples            number of samples in time series
        sampling_rate        sampling rate in samples/second
        start_time_epoch_sec start time in epoch seconds
        start_time_utc       start time in UTC
        station              station name
        units                units of time series
        ==================== ==================================================
    
    .. note:: Currently only supports hdf5 and text files
    """
    
    def __init__(self, **kwargs):
        
        self.station = 'mt00'
        self.sampling_rate = 1
        self.start_time_epoch_sec = 0.0
        self.start_time_utc = '1980-01-01,00:00:00'
        self.n_samples = 0
        self.component = None
        self.coordinate_system = 'geomagnetic'
        self.dipole_length = 0
        self.azimuth = 0
        self.units = 'mV'
        self.lat = 0.0
        self.lon = 0.0
        self.datum = 'WGS84'
        self.data_logger = 'Zonge Zen'
        self.instrument_id = None
        self.calibration_fn = None
        self.declination = 0.0
        self.ts = None 
        self.fn_hdf5 = None
        self.fn_txt = None
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
        
    def write_hdf5(self, fn_hdf5, ts=None):
        """
        Write a hdf5 time series file
        
        Arguments
        ---------------
            **fn_hdf5** : string
                          full path to save hdf5 file
                          
            **ts** : numpy.ndarray (n_samples)
                     time series data 
                     
        """
        
        
        self.fn_hdf5 = fn_hdf5
        if ts is not None:
            self.ts = ts
            
        if self.ts is None:
            raise MT_TS_Error('Time series is None')
        
        # check to see if a file already exists, it doesn't like it when
        # one already does
        if os.path.exists(self.fn_hdf5) is True:
            print('{0} already exists'.format(self.fn_hdf5))
            delete_bool = input('Delete old file (True/False)')
            if delete_bool == True:
                os.remove(self.fn_hdf5)
            else:
                print('Did not overwrite {0}'.format(self.fn_hdf5))
                return 
        
        # create the hdf5 object
        hdf5_obj = tables.open_file(self.fn_hdf5, mode='w', title='Test')
        
        # make a time series array        
        ts_arr = hdf5_obj.create_array('/', 'time_series', self.ts)
        
        # add attributes
        ts_arr.attrs.station = self.station
        ts_arr.attrs.sampling_rate = self.sampling_rate
        ts_arr.attrs.start_time_epoch_sec = self.start_time_epoch_sec 
        ts_arr.attrs.start_time_utc = self.start_time_utc
        ts_arr.attrs.n_samples = self.n_samples
        ts_arr.attrs.component = self.component
        ts_arr.attrs.coordinate_system = self.coordinate_system
        ts_arr.attrs.dipole_length = self.dipole_length
        ts_arr.attrs.azimuth = self.azimuth
        ts_arr.attrs.units = self.units
        ts_arr.attrs.lat = self.lat
        ts_arr.attrs.lon = self.long
        ts_arr.attrs.datum = self.datum
        ts_arr.attrs.data_logger = self.data_logger
        ts_arr.attrs.instrument_id = self.instrument_id
        ts_arr.attrs.calibration_fn = self.calibration_fn
        ts_arr.attrs.declination = self.declination
        
        z1_h5.close()
        
    def read_hdf5(self, fn_hdf5):
        """
        read in an hdf5 file and fill attributes accordingly
        """
        self.fn_hdf5 = fn_hdf5        
        
        if not os.path.isfile(self.fn_hdf5):
            raise MT_TS_Error('Could not find {0}'.format(self.fn_hdf5))

        hdf5_obj = tables.open_file(self.fn_hdf5, 'r')
        
        self.ts = hdf5_obj.get_node('/', 'time_series')
        
        self.station = self.ts.attrs.station
        self.sampling_rate = self.ts.attrs.sampling_rate
        self.start_time_epoch_sec = self.ts.attrs.start_time_epoch_sec
        self.start_time_utc = self.ts.attrs.start_time_utc
        self.n_samples = self.ts.attrs.n_samples
        self.component = self.ts.attrs.component
        self.coordinate_system = self.ts.attrs.coordinate_system
        self.dipole_length = self.ts.attrs.dipole_length
        self.azimuth = self.ts.attrs.azimuth
        self.units = self.ts.attrs.units
        self.lat = self.ts.attrs.lat
        self.lon = self.ts.attrs.lon
        self.datum = self.ts.attrs.datum
        self.data_logger = self.ts.attrs.data_logger
        self.instrument_id = self.ts.attrs.instrument_id
        self.calibration_fn = self.ts.attrs.calibration_fn
        self.declination = self.ts.attrs.declination
        
        print 'Read in {0}'.format(self.fn_hdf5)
        
        
                
#==============================================================================
# Error classes
#==============================================================================
class MT_TS_Error(Exception):
    pass        
    
z1 = zen.Zen3D(fn)
z1.read_z3d()
z1.station = '{0}{1}'.format(z1.metadata.line_name, z1.metadata.rx_xyz0[0:2])

h5_fn = fn[0:-4]+'.h5'
pd_h5_fn = fn[0:-4]+'_pd.h5' 

if not os.path.exists(h5_fn):

    z1_h5 = tables.open_file(h5_fn, mode='w', title='Test')
    ts_arr = z1_h5.create_array('/', 'time_series', z1.convert_counts())
    
    ts_arr.attrs.station = z1.station
    ts_arr.attrs.sampling_rate = int(z1.df)
    ts_arr.attrs.start_time_epoch_sec = time.mktime(time.strptime(z1.zen_schedule, 
                                                                  zen.datetime_fmt))
    ts_arr.attrs.start_time_utc = z1.zen_schedule
    ts_arr.attrs.n_samples = int(z1.time_series.size)
    ts_arr.attrs.component = z1.metadata.ch_cmp
    ts_arr.attrs.coordinate_system = 'geomagnetic'
    ts_arr.attrs.dipole_length = float(z1.metadata.ch_length)
    ts_arr.attrs.azimuth = float(z1.metadata.ch_azimuth)
    ts_arr.attrs.units = 'mV'
    ts_arr.attrs.lat = z1.header.lat
    ts_arr.attrs.lon = z1.header.long
    ts_arr.attrs.datum = 'WGS84'
    ts_arr.attrs.data_logger = 'Zonge Zen'
    ts_arr.attrs.instrument_num = None
    ts_arr.attrs.calibration_fn = None
    ts_arr.attrs.declination = 3.6
    
    z1_h5.close()
    
#if not os.path.exists(pd_h5_fn):
#
#    cols = pd.MultiIndex.from_product([[z1.station], 
#                                       [int(z1.df)]],
#                                      names=['station', 
#                                             'sampling_rate'])
#
#    z1_df = pd.DataFrame(data=z1.convert_counts(), columns=cols)
#    
#    z1_store = pd.HDFStore(pd_h5_fn, 'w')
#    z1_store.put('time_series', z1_df)
#        
##    z1_store.get_storer('time_series').attrs.station = z1.station
##    z1_store.get_storer('time_series').attrs.sampling_rate = int(z1.df)
##    z1_store.get_storer('time_series').attrs.start_time_epoch_sec = time.mktime(time.strptime(z1.zen_schedule, 
##                                                                  zen.datetime_fmt))
##    z1_store.get_storer('time_series').attrs.start_time_utc = z1.zen_schedule
##    z1_store.get_storer('time_series').attrs.n_samples = int(z1.time_series.size)
##    z1_store.get_storer('time_series').attrs.component = z1.metadata.ch_cmp
##    z1_store.get_storer('time_series').attrs.coordinate_system = 'geomagnetic'
##    z1_store.get_storer('time_series').attrs.dipole_length = float(z1.metadata.ch_length)
##    z1_store.get_storer('time_series').attrs.azimuth = float(z1.metadata.ch_azimuth)
##    z1_store.get_storer('time_series').attrs.units = 'mV'
##    z1_store.get_storer('time_series').attrs.lat = z1.header.lat
##    z1_store.get_storer('time_series').attrs.lon = z1.header.long
##    z1_store.get_storer('time_series').attrs.datum = 'WGS84'
##    z1_store.get_storer('time_series').attrs.instrument = 'Zonge Zen'
##    z1_store.get_storer('time_series').attrs.calibration_fn = None
##    z1_store.get_storer('time_series').attrs.declination = 3.6
#    
#    z1_store.close()


#store = pd.HDFStore(h5_fn, 'r')


