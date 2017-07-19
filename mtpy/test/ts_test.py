# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 12:04:39 2017

@author: jpeacock
"""

import mtpy.core.ts as mtts
reload(mtts)


fn = r"d:\Peacock\MTData\Umatilla\hf05\hf05_20170517_230518_256_EX.Z3D"

def read_z3d(fn):
    import mtpy.usgs.zen as zen
    ## TEST Writing
    z1 = zen.Zen3D(fn)
    z1.read_z3d()
    z1.station = '{0}{1}'.format(z1.metadata.line_name, z1.metadata.rx_xyz0[0:2])
    
    #h5_fn = r"d:\Peacock\MTData\Umatilla\hf05\hf05_20170517_193018_256_EX.h5"
    ts_obj = mtts.MT_TS()
    
    ts_obj.ts = z1.convert_counts()
    ts_obj.station = z1.station
    ts_obj.sampling_rate = int(z1.df)
    ts_obj.start_time_utc = z1.zen_schedule
    ts_obj.n_samples = int(z1.time_series.size)
    ts_obj.component = z1.metadata.ch_cmp
    ts_obj.coordinate_system = 'geomagnetic'
    ts_obj.dipole_length = float(z1.metadata.ch_length)
    ts_obj.azimuth = float(z1.metadata.ch_azimuth)
    ts_obj.units = 'mV'
    ts_obj.lat = z1.header.lat
    ts_obj.lon = z1.header.long
    ts_obj.datum = 'WGS84'
    ts_obj.data_logger = 'Zonge Zen'
    ts_obj.instrument_num = None
    ts_obj.calibration_fn = None
    ts_obj.declination = 3.6
    
    return ts_obj

def make_hdf5_from_z3d(z3d_fn):
    
    ts_obj = read_z3d(z3d_fn)
    
    h5_fn = z3d_fn[0:-4]+'.h5'
    ts_obj.write_hdf5(h5_fn)
    
    return h5_fn
    
def make_txt_from_hdf5(hdf5_fn, chunk=4096):
    ts_obj = mtts.MT_TS()
    ts_obj.read_hdf5(hdf5_fn)
    ts_obj.write_ascii_file(chunk_size=chunk)
    
    return ts_obj.fn_ascii
    
def read_txt(txt_fn):
    ts_obj = mtts.MT_TS()
    ts_obj.read_ascii(txt_fn)
    
    return ts_obj
    
#==============================================================================
# try a test
#==============================================================================
ts_obj = read_z3d(fn)

#h5_fn = make_hdf5_from_z3d(fn)
#txt_fn = make_txt_from_hdf5(h5_fn, chunk=8192)
#ts_obj = read_txt(txt_fn)


