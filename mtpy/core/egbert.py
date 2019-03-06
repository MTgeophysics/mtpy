# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:34:23 2017

@author: jrpeacock
"""

#==============================================================================
# Imports
#==============================================================================
import os

import numpy as np
import mtpy.core.z as z
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
        
        for key in list(kwargs.keys()):
            setattr(self, key, kwargs[key])
            
    def read_egbert_file(self, z_fn=None):
        """
        Read in Egbert zrr file
        """
        if z_fn is not None:
            self.z_fn = z_fn
        
        self.read_header()
            