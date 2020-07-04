# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:34:23 2017
@author: jrpeacock
"""

#==============================================================================
# Imports
#==============================================================================
import numpy as np
import mtpy.core.z as mtz
import mtpy.utils.gis_tools as gis_tools
#==============================================================================
class ZMMError(Exception):
    pass

class Channel(object):
    """
    class to hold channel information
    """
    def __init__(self, channel_dict=None):
        self.number = None
        self.azimuth = None
        self.tilt = None
        self.dl = None
        self.channel = None
        
        if channel_dict is not None:
            self.from_dict(channel_dict)
    
    @property
    def index(self):
        if self.number is not None:
            return self.number - 1
        else:
            return None
        
    def from_dict(self, channel_dict):
        """
        fill attributes from a dictionary
        """
        
        for key, value in channel_dict.items():
            if key in ['azm', 'azimuth']:
                self.azimuth = value
            elif key in ['chn_num', 'number']:
                self.number = value
            elif key in ['tilt']:
                self.tilt = value
            elif key in ['dl']:
                self.dl = value
            elif key in ['channel']:
                self.channel = value
            

class ZMMHeader(object):
    """
    Container for Header of an Egbert file
    """
    
    def __init__(self, z_fn=None, **kwargs):
        
        self.description = None
        self.processing_type = None
        self.station = None
        self._lat = None
        self._lon = None
        self.declination = None
        self.num_channels = None
        self.num_freq = None
        self._header_count = 0
        self._component_dict = None
        self.ex = None
        self.ey = None
        self.hx = None
        self.hy = None
        self.hz = None
        
    @property
    def lat(self):
        return self._lat
    
    @lat.setter
    def lat(self, lat):
        self._lat = gis_tools.assert_lat_value(lat)
        
    @property
    def lon(self):
        return self._lon
    
    @lon.setter
    def lon(self, lon):
        self._lon = gis_tools.assert_lon_value(lon)
        
    def read_header(self, z_fn=None):
        """
        read header information
        """
        
        if z_fn is not None:
            self.z_fn = z_fn
            
        with open(self.z_fn, 'r') as fid:
            line = fid.readline()
 
            self._header_count = 0           
            header_list = []
            while 'period' not in line:
                header_list.append(line)
                self._header_count += 1
                
                line = fid.readline() 
                
        self.description = ''
        self.station = header_list[3].lower().strip()
        for ii, line in enumerate(header_list):
            if line.find('**') >= 0:
                self.description += line.replace('*', '').strip()
            elif ii == 2:
                self.processing_type = line.lower().strip()
            elif 'station' in line:
                self.station = line.split(':')[1].strip()
            elif 'coordinate' in line:
                line_list = line.strip().split()
                self.lat = line_list[1]
                try:
                    self.lon = line_list[2]
                except ValueError:
                    self.lon = float(line_list[2])%180
                    
                self.declination = float(line_list[-1])
            elif 'number' in line:
                line_list = line.strip().split()
                self.num_channels = int(line_list[3])
                self.num_freq = int(line_list[-1])
            elif 'orientations' in line:
                pass
            elif line.strip()[-2:].lower() in ['ex', 'ey', 'hx', 'hy', 'hz']:
                line_list = line.strip().split()
                comp = line_list[-1].lower()
                channel_dict = {'channel':comp}
                channel_dict['chn_num'] = int(line_list[0]) % self.num_channels
                channel_dict['azm'] = float(line_list[1])
                channel_dict['tilt'] = float(line_list[2])
                channel_dict['dl'] = line_list[3]
                if channel_dict['chn_num'] == 0:
                    channel_dict['chn_num'] = self.num_channels
                setattr(self, comp, Channel(channel_dict))
    

class ZMM(ZMMHeader):
    """
    Container for Egberts zrr format.
    
    """
    
    def __init__(self, z_fn=None, **kwargs):
        
#        EgbertHeader.__init__(self, **kwargs)
        super(ZMM, self).__init__()
        
        self.z_fn = z_fn
        self._header_count = 0
        self.Z = mtz.Z()
        self.Tipper = mtz.Tipper()
        self.transfer_functions = None
        self.sigma_e = None
        self.sigma_s = None
        self.period = None
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    @property
    def frequency(self):
        if self.period is None:
            return None
        return 1./self.period
            
    def initialize_arrays(self):
        """
        make initial arrays based on number of frequencies and channels
        """
        if self.num_freq is None:
            return
        
        self.period = np.zeros(self.num_freq)
        self.transfer_functions = np.zeros((self.num_freq, 
                                            self.num_channels - 2, 2),
                                           dtype=np.complex64)

        # residual covariance -- square matrix with dimension as number of
        # predicted channels
        self.sigma_e = np.zeros((self.num_freq, 
                                 self.num_channels- 2, self.num_channels - 2),
                                dtype=np.complex64)

        # inverse coherent signal power -- square matrix, with dimension as the
        #    number of predictor channels
        # since EMTF and this code assume N predictors is 2,
        #    this dimension is hard-coded
        self.sigma_s = np.zeros((self.num_freq, 2, 2), dtype=np.complex64)
                  
    def read_zmm_file(self, z_fn=None):
        """
        Read in Egbert zrr file
        """
        if z_fn is not None:
            self.z_fn = z_fn
        
        self.read_header()
        self.initialize_arrays() 
        
        ### read each data block and fill the appropriate array
        for ii, period_block in enumerate(self._get_period_blocks()):
            data_block = self._read_period_block(period_block)
            self.period[ii] = data_block['period']
            
            self._fill_tf_array_from_block(data_block['tf'], ii)
            self._fill_sig_array_from_block(data_block['sig'], ii)
            self._fill_res_array_from_block(data_block['res'], ii)
            
        ### make Z and Tipper
        self.Z = self.calculate_impedance()
         
        
        try:
            self.Tipper = self.calculate_tippers()
        except ZMMError:
            self.Tipper = mtz.Tipper()
            print('*** No HZ found cannot calculate induction vectors. ***')
        
    def _get_period_blocks(self):
        """
        split file into period blocks
        """
        
        with open(self.z_fn, 'r') as fid:
            fn_str = fid.read()
        
        period_strings = fn_str.lower().split('period')
        period_blocks = []
        for per in period_strings:
            period_blocks.append(per.split('\n'))
        
        return period_blocks[1:]
        
    def _read_period_block(self, period_block):
        """
        read block:
            period :      0.01587    decimation level   1    freq. band from   46 to   80
            number of data point  951173 sampling freq.   0.004 Hz
             Transfer Functions
              0.1474E+00 -0.2049E-01  0.1618E+02  0.1107E+02
             -0.1639E+02 -0.1100E+02  0.5559E-01  0.1249E-01
             Inverse Coherent Signal Power Matrix
              0.2426E+03 -0.2980E-06
              0.9004E+02 -0.2567E+01  0.1114E+03  0.1192E-06
             Residual Covaraince
              0.8051E-05  0.0000E+00
             -0.2231E-05 -0.2863E-06  0.8866E-05  0.0000E+00
        """
        
        period = float(period_block[0].strip().split(':')[1].split()[0].strip())
        
        data_dict = {'period':period, 'tf':[], 'sig':[], 'res':[]}
        key = 'tf'
        for line in period_block[2:]:
            if 'transfer' in line.lower():
                key = 'tf'
                continue
            elif 'signal' in line.lower():
                key = 'sig'
                continue
            elif 'residual' in line.lower():
                key = 'res'
                continue
            
            line_list = [float(xx) for xx in line.strip().split()]
            values = [complex(line_list[ii], line_list[ii+1])
                      for ii in range(0, len(line_list), 2)]
            data_dict[key].append(values)
            # for ii in range(0, len(line_list), 2):
            #     data_dict[key].append(complex(line_list[ii], line_list[ii+1]))
        
        return data_dict
    
    def _flatten_list(self, x_list):
        """
        flatten = lambda l: [item for sublist in l for item in sublist]

        Returns
        -------
        None.

        """
        
        flat_list = [item for sublist in x_list for item in sublist]
        
        return flat_list
    
    def _fill_tf_array_from_block(self, tf_block, index):
        """
        fill tf arrays from data blocks
        """
        tf_block = self._flatten_list(tf_block)
        for kk, jj in enumerate(range(0, len(tf_block), 2)):
            self.transfer_functions[index, kk, 0] = tf_block[jj]
            self.transfer_functions[index, kk, 1] = tf_block[jj+1]
            
    def _fill_sig_array_from_block(self, sig_block, index):
        """
        fill signal array
        """
        sig_block = self._flatten_list(sig_block)
        self.sigma_s[index, 0, 0] = sig_block[0]
        self.sigma_s[index, 1, 0] = sig_block[1]
        self.sigma_s[index, 0, 1] = sig_block[1]
        self.sigma_s[index, 1, 1] = sig_block[2]
        
    def _fill_res_array_from_block(self, res_block, index):
        """
        fill residual covariance array
        """
        for jj in range(self.num_channels - 2):
            values = res_block[jj]
            for kk in range(jj + 1):
                if jj == kk:
                    self.sigma_e[index, jj, kk] = values[kk]
                else:
                    self.sigma_e[index, jj, kk] = values[kk]
                    self.sigma_e[index, kk, jj] = values[kk].conjugate()
    
    def calculate_impedance(self, angle=0.):
        """
        calculate the impedances from the transfer functions
        """
        
        # check to see if there are actually electric fields in the TFs
        if not hasattr(self, 'ex') or not hasattr(self, 'ey'):
            raise ZMMError("Cannot return apparent resistivity and phase "
                             "data because these TFs do not contain electric "
                             "fields as a predicted channel.")

        # transform the TFs first...
        # build transformation matrix for predictor channels
        #    (horizontal magnetic fields)
        hx_index = self.hx.index
        hy_index = self.hy.index
        u = np.eye(2, 2)
        u[hx_index, hx_index] = np.cos(np.deg2rad(self.hx.azimuth - angle))
        u[hx_index, hy_index] = np.sin(np.deg2rad(self.hx.azimuth - angle))
        u[hy_index, hx_index] = np.cos(np.deg2rad(self.hy.azimuth - angle))
        u[hy_index, hy_index] = np.sin(np.deg2rad(self.hy.azimuth - angle))
        u = np.linalg.inv(u)

        # build transformation matrix for predicted channels (electric fields)
        ex_index = self.ex.index
        ey_index = self.ey.index
        v = np.eye(self.transfer_functions.shape[1],
                   self.transfer_functions.shape[1])
        v[ex_index-2, ex_index-2] = np.cos(np.deg2rad(self.ex.azimuth - angle))
        v[ey_index-2, ex_index-2] = np.sin(np.deg2rad(self.ex.azimuth - angle))
        v[ex_index-2, ey_index-2] = np.cos(np.deg2rad(self.ey.azimuth - angle))
        v[ey_index-2, ey_index-2] = np.sin(np.deg2rad(self.ey.azimuth - angle))

        # matrix multiplication...
        rotated_transfer_functions = np.matmul(v, 
                                               np.matmul(self.transfer_functions,
                                                         u.T))
        rotated_sigma_s = np.matmul(u, np.matmul(self.sigma_s, u.T))
        rotated_sigma_e = np.matmul(v, np.matmul(self.sigma_e, v.T))

        # now pull out the impedance tensor
        z = np.zeros((self.num_freq, 2, 2), dtype=np.complex64)
        z[:, 0, 0] = rotated_transfer_functions[:, ex_index-2, hx_index]   # Zxx
        z[:, 0, 1] = rotated_transfer_functions[:, ex_index-2, hy_index]   # Zxy
        z[:, 1, 0] = rotated_transfer_functions[:, ey_index-2, hx_index]   # Zyx
        z[:, 1, 1] = rotated_transfer_functions[:, ey_index-2, hy_index]   # Zyy

        # and the variance information
        var = np.zeros((self.num_freq, 2, 2))
        var[:, 0, 0] = np.real(rotated_sigma_e[:, ex_index-2, ex_index-2] *
                               rotated_sigma_s[:, hx_index, hx_index])
        var[:, 0, 1] = np.real(rotated_sigma_e[:, ex_index-2, ex_index-2] *
                               rotated_sigma_s[:, hy_index, hy_index])
        var[:, 1, 0] = np.real(rotated_sigma_e[:, ey_index-2, ey_index-2] *
                               rotated_sigma_s[:, hx_index, hx_index])
        var[:, 1, 1] = np.real(rotated_sigma_e[:, ey_index-2, ey_index-2] *
                               rotated_sigma_s[:, hy_index, hy_index])

        error = np.sqrt(var)

        z_object = mtz.Z(z, error, self.frequency)
        return z_object
    
    def calculate_tippers(self, angle=0.):
        """
        calculate induction vectors
        """
        
        # check to see if there is a vertical magnetic field in the TFs
        if self.hz is None:
            raise ZMMError("Cannot return tipper data because the TFs do not "
                             "contain the vertical magnetic field as a "
                             "predicted channel.")

        # transform the TFs first...
        # build transformation matrix for predictor channels
        #    (horizontal magnetic fields)
        hx_index = self.hx.index
        hy_index = self.hy.index
        u = np.eye(2, 2)
        u[hx_index, hx_index] = np.cos(np.deg2rad(self.hx.azimuth - angle))
        u[hx_index, hy_index] = np.sin(np.deg2rad(self.hx.azimuth - angle))
        u[hy_index, hx_index] = np.cos(np.deg2rad(self.hy.azimuth - angle))
        u[hy_index, hy_index] = np.sin(np.deg2rad(self.hy.azimuth - angle))
        u = np.linalg.inv(u)

        # don't need to transform predicated channels (assuming no tilt in Hz)
        hz_index = self.hz.index
        v = np.eye(self.transfer_functions.shape[1],
                   self.transfer_functions.shape[1])

        # matrix multiplication...
        rotated_transfer_functions = \
            np.matmul(v, np.matmul(self.transfer_functions, u.T))
        rotated_sigma_s = np.matmul(u, np.matmul(self.sigma_s, u.T))
        rotated_sigma_e = np.matmul(v, np.matmul(self.sigma_e, v.T))

        # now pull out tipper information
        tipper = np.zeros((self.num_freq, 2), dtype=np.complex64)
        tipper[:, 0] = rotated_transfer_functions[:, hz_index-2, hx_index]   # Tx
        tipper[:, 1] = rotated_transfer_functions[:, hz_index-2, hy_index]   # Ty

        # and the variance/error information
        var = np.zeros((self.num_freq, 2))
        var[:, 0] = np.real(rotated_sigma_e[:, hz_index-2, hz_index-2] *
                            rotated_sigma_s[:, hx_index, hx_index])   # Tx
        var[:, 1] = np.real(rotated_sigma_e[:, hz_index-2, hz_index-2] *
                            rotated_sigma_s[:, hy_index, hy_index])   # Ty
        error = np.sqrt(var)
        
        tipper = tipper.reshape((self.num_freq, 1, 2))
        error = error.reshape((self.num_freq, 1, 2))
        
        tipper_obj = mtz.Tipper(tipper, error, self.frequency)

        return tipper_obj

