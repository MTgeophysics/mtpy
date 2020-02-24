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
class EgbertHeader(object):
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
        self._lon = -180 + gis_tools.assert_lon_value(lon)
        
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
        self.component_dict = {}
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
                self.component_dict[comp] = {}
                self.component_dict[comp]['chn_num'] = int(line_list[0])
                self.component_dict[comp]['azm'] = float(line_list[1])
                self.component_dict[comp]['tilt'] = float(line_list[2])
                self.component_dict[comp]['dl'] = line_list[3]
    

class EgbertZ(EgbertHeader):
    """
    Container for Egberts zrr format.
    
    """
    
    def __init__(self, z_fn=None, **kwargs):
        
#        EgbertHeader.__init__(self, **kwargs)
        super(EgbertZ, self).__init__()
        
        self.z_fn = z_fn
        self._header_count = 0
        self.Z = None
        self.Tipper = None
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    def read_egbert_file(self, z_fn=None):
        """
        Read in Egbert zrr file
        """
        if z_fn is not None:
            self.z_fn = z_fn
        
        self.read_header()
        
        period_list = []
        z_list = []
        for period_block in self._get_period_blocks():
            z_dict = self._read_period_block(period_block)
            period_list.append(z_dict['period'])
            z_list.append(z_dict['z'])
        
        # make the lists arrays
        z_list = np.array(z_list)
        period_list = np.array(period_list)
            
        # make z objects
        self.Z = mtz.Z(z_array=z_list,
                       z_err_array=np.zeros_like(z_list, dtype=np.float),
                       freq=1./period_list)
        
        
#        self.mt_obj = mt.MT()
#        self.mt_obj.lat = self.lat
#        self.mt_obj.lon = self.lon
#        self.mt_obj.station = self.station
#        self.mt_obj.Z = self.Z
#        self.mt_obj.Tipper = mtz.Tipper(np.zeros((z_list.shape[0], 1, 2), 
#                                                 dtype=np.complex),
#                                        np.zeros((z_list.shape[0], 1, 2), 
#                                                 dtype=np.float),
#                                        1./np.array(period_list))
#        self.mt_obj.Notes.info_dict = {'notes':'processed with EMTF'}
        
        
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
        z_dict = {}
        
        p_list = period_block[0].strip().split(':')
        z_dict['period'] = float(p_list[1].split()[0].strip())
        
        data_dict = {'tf':{}, 'sig':{}, 'res':{}}
        key = 'tf'
        # for line in period_block[2:]:
        #     if line in 
        
        
        zx_list = [float(xx) for xx in period_block[3].strip().split()]
        zy_list = [float(yy) for yy in period_block[4].strip().split()]
        z_arr = np.zeros((2, 2), dtype=np.complex)
        z_arr[0, 0] = zx_list[0]+1j*zx_list[1]
        z_arr[0, 1] = zx_list[2]+1j*zx_list[3]
        z_arr[1, 0] = zy_list[0]+1j*zy_list[1]
        z_arr[1, 1] = zy_list[2]+1j*zy_list[3]
        
        z_dict['z'] = z_arr
        
        # get errors
        
        return z_dict
    
#    def zmm_to_edi(self, edi_fn=None):
#        """
#        write a simple edi file from zmm
#        """
#        if edi_fn is None:
#            edi_dir = os.path.dirname(self.z_fn)
#            edi_basename = '{0}.edi'.format(self.station)
#        else:
#            edi_dir = os.path.dirname(edi_fn)
#            edi_basename = os.path.join(edi_fn)
#        self.mt_obj.write_mt_file(save_dir=edi_dir,
#                                  fn_basename=edi_basename)
