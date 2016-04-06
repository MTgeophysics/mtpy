# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:03:31 2015

@author: jpeacock
"""

import os
import numpy as np
import time
import datetime

import mtpy.utils.format as MTft
import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex
import mtpy.utils.filehandling as MTfh
import mtpy.core.z as MTz

tab = ' '*4

class Edi(object):
    """
    
    """

    def __init__(self, edi_fn=None, **kwargs):
        
        self.edi_fn = edi_fn
        self.header = Header(edi_fn=self.edi_fn)
        self.info = Information(edi_fn=self.edi_fn)

#==============================================================================
#  Header object        
#==============================================================================
class Header(object):
    """
    Header object
    """
    
    def __init__(self, edi_fn=None, **kwargs):
        self.edi_fn = edi_fn
        self.dataid = None
        self.acqby = None
        self.fileby = None
        self.acqdate = None
        self.units = None
        self.filedate = datetime.datetime.utcnow().strftime(
                                                    '%Y/%m/%d %H:%M:%S UTC')       
        self.loc = None
        self.lat = None
        self.lon = None
        self.elev = None
        self.empty = 1E32
        self.progvers = None
        self.progdate = None
        self.phoenix_edi = False
        
        self.header_list = None
        
        self._header_keys = ['acqby', 
                             'acqdate',
                             'dataid',
                             'elev',
                             'fileby',
                             'lat',
                             'loc',
                             'lon',
                             'filedate',
                             'empty',
                             'progdate',
                             'progvers']
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
        if self.edi_fn is not None:
            self.read_header()
            
    def get_header_list(self):
        """
        Get the header information from the .edi file in the form of a list,
        where each item is a line in the header section.
        """
       
        if self.edi_fn == None:
            print 'No edi file to read'
            return
        if os.path.isfile(self.edi_fn) == False:
            print 'Could not find {0}, check path'.format(self.edi_fn)
            
        self.header_list = []
        head_find = False
        count = 0
        with open(self.edi_fn, 'r') as fid:
            for line in fid:
                if line.find('>') == 0:
                    count += 1
                    if line.lower().find('head') > 0:
                        head_find = True
                    else:
                        head_find = False
                    if count == 2 and head_find == False:
                        break
                elif count == 1 and line.find('>') != 0 and head_find == True:
                    line = line.strip()
                    # skip any blank lines
                    if len(line) > 2:
                        self.header_list.append(line.strip())
                        
        self.header_list = self._validate_header_list(self.header_list)
    
    def read_header(self, header_list=None):
        """
        read a header information
        """

        if header_list is not None:
            self.header_list = self._validate_header_list(header_list)
            
        if self.header_list is None and self.edi_fn is None:
            print 'Nothing to read. header_list and edi_fn are None'
            
        if self.header_list is None and self.edi_fn is not None:
            self.get_header_list()
        
        for h_line in self.header_list:
            h_list = h_line.split('=')
            key = h_list[0].lower()
            value = h_list[1].replace('"', '')
            
            if key in 'latitude':
                key = 'lat'
                value = MTft._assert_position_format(key, value)
            
            elif key in 'longitude':
                key = 'lon'
                value = MTft._assert_position_format(key, value)
                
            elif key in 'elevation':
                key = 'elev'
                try:
                    value = float(value)
                except ValueError:
                    value = 0.0
                    print 'No elevation data'
                    
            elif key in ['country', 'state', 'loc', 'location', 'prospect']:
                key = 'loc'
                try:
                    if self.__dict__[key] is not None:
                        value = '{0}, {1}'.format(self.__dict__[key], value)
                except KeyError:
                    pass
            # test if its a phoenix formated .edi file
            elif key in ['progvers']:
                if value.lower().find('mt-editor') != -1:
                    self.phoenix_edi = True
                    
            elif key in ['fileby']:
                if value == '':
                    value = 'mtpy'
                
            setattr(self, key, value)
            
    def write_header(self, header_list=None):
        """
        write header information to a list of lines
        """

        if header_list is not None:
            self.header_list = self._validate_header_list(header_list)
            
        if self.header_list is None and self.edi_fn is None:
            print 'Nothing to read. header_list and edi_fn are None'
            
        if self.header_list is None and self.edi_fn is not None:
            self.get_header_list()
            
        header_lines = ['>HEAD\n\n']
        for key in sorted(self._header_keys):
            value = self.__dict__[key]
            if key in ['progdate', 'progvers']:
                if value is None:
                    value = 'mtpy'
            elif key in ['lat', 'lon']:
                value = MTft.convert_dms_tuple2string(
                                        MTft.convert_degrees2dms_tuple(value))
            if key in ['elev']:
                try:
                    value = '{0:.3f}'.format(value)
                except ValueError:
                    value = '0.000'
                    
            if key in ['filedate']:
                value = datetime.datetime.utcnow().strftime(
                                                    '%Y/%m/%d %H:%M:%S UTC')
                                                    
            header_lines.append('{0}{1}={2}\n'.format(tab, key.upper(), value))
        header_lines.append('\n')
        return header_lines
        
    def _validate_header_list(self, header_list):
        """
        make sure the input header list is valid
        
        returns a validated header list
        """
        
        if header_list is None:
            print 'No header information to read'
            return None
            
        new_header_list = []
        for h_line in header_list:
            h_line = h_line.strip().replace('"', '')
            if len(h_line) > 1:
                h_list = h_line.split('=')
                if len(h_list) == 2:
                    key = h_list[0]
                    value = h_list[1]
                    new_header_list.append('{0}={1}'.format(key, value))
        
        return new_header_list
        
#==============================================================================
# Info object
#==============================================================================
class Information(object):
    """
    Contain, read, and write info section of .edi file
    
    not much to really do here, but just keep it in the same format that it is
    read in as, except if it is in phoenix format then split the two paragraphs
    up
    """
    
    def __init__(self, edi_fn=None):
        self.edi_fn = edi_fn
        self.info_list = None
        
        if self.edi_fn is not None:
            self.read_info()
            
    def get_info_list(self):
        """
        get a list of lines from the info section
        """
        
        if self.edi_fn is None:
            print 'no edi file input, check edi_fn attribute'            
            return
        if os.path.isfile(self.edi_fn) is False:
            print 'Could not find {0}, check path'.format(self.edi_fn)
            return
            
        self.info_list = []
        info_find = False
        phoenix_file = False
        phoenix_list_02 = []
        count = 0
        with open(self.edi_fn, 'r') as fid:
            for line in fid:
                if line.find('>') == 0:
                    count += 1
                    if line.lower().find('info') > 0:
                        info_find = True
                    else:
                        info_find = False
                    if count > 2 and info_find == False:
                        break
                elif count > 1 and line.find('>') != 0 and info_find == True:
                    if line.lower().find('run information') >= 0:
                        phoenix_file = True
                    if phoenix_file == True and len(line) > 40:
                        self.info_list.append(line[0:37].strip())
                        phoenix_list_02.append(line[38:].strip())
                    else:
                        if len(line.strip()) > 1:
                            self.info_list.append(line.strip())
                        
        self.info_list += phoenix_list_02
        # validate the information list
        self.info_list = self._validate_info_list(self.info_list)
        
    def read_info(self, info_list=None):
        """
        read information section of the .edi file
        """
        
        if info_list is not None:
            self.info_list = self._validate_info_list(info_list)
            
        if self.edi_fn is not None and self.info_list is None:
            self.get_info_list()
            
        if self.info_list is None:
            print "Could not read information"
            return
            
    def write_info(self, info_list=None):
        """
        
        """
        
        if info_list is not None:
            self.info_list = self._validate_info_list(info_list)
        
            
        info_lines = ['>INFO\n\n']
        for line in self.info_list:
            info_lines.append('{0}{1}\n'.format(tab, line))
        
        return info_lines
            
         
    def _validate_info_list(self, info_list):
        """
        check to make sure the info list input is valid, really just checking
        for Phoenix format where they put two columns in the file and remove 
        any blank lines and the >info line
        """
        
        new_info_list = []
        for line in info_list:
            # get rid of empty lines
            lt = str(line).strip()
            if len(lt) > 1:
                if line.find('>') == 0:
                    pass
                else:
                    new_info_list.append(line.strip())
                    
        return new_info_list
        
            
#==============================================================================
#  Define measurement class       
#==============================================================================
class DefineMeasurement(object):
    """
    hold information about the measurement
    """
    
    def __init__(self, edi_fn=None):
        self.edi_fn = edi_fn
        self.measurement_list = None
        
        self.maxchan = None
        self.maxmeas = 99999
        self.maxrun = 999
        self.refelev = None
        self.reflat = None
        self.reflon = None
        self.reftype = None
        self.units = 'm'
        
        if self.edi_fn is not None:
            self.read_define_measurement()
        
    def get_measurement_lists(self):
        """
        get measurement list including measurement setup
        """
        if self.edi_fn is None:
            print 'No edi file input, check edi_fn attribute'
            return 
            
        if os.path.isfile(self.edi_fn) is False:
            print 'Could not find {0}, check path'.format(self.edi_fn)
            
        self.measurement_list = []
        meas_find = False
        count = 0
        with open(fn, 'r') as fid:
            for line in fid:
                if line.find('>=') == 0:
                    count += 1
                    if line.lower().find('definemeas') > 0:
                        meas_find = True
                    else:
                        meas_find = False
                    if count == 2 and meas_find == False:
                        break
                elif count == 1 and line.find('>') != 0 and meas_find == True:
                    line = line.strip()
                    if len(line) > 2:
                        self.measurement_list.append(line.strip())
                
                # look for the >XMEAS parts
                elif count == 1 and line.find('>') == 0 and meas_find == True:
                    if line.find('!') > 0:
                        pass
                    else:
                        line_list = line.strip().split()
                        m_dict = {}
                        for ll in line_list[1:]:
                            ll_list = ll.split('=')
                            key = ll_list[0].lower()
                            value = ll_list[1]
                            m_dict[key] = value
                        self.measurement_list.append(m_dict)
                        
    def read_define_measurement(self, measurement_list=None):
        """
        read the define measurment section of the edi file
        
        should be a list with lines for:
            - maxchan 
            - maxmeas
            - maxrun
            - refelev
            - reflat
            - reflon
            - reftype
            - units
            - dictionaries for >XMEAS with keys:
                - id
                - chtype
                - x
                - y
                - axm
                -acqchn
        
        """
        
        if measurement_list is not None:
            self.measurement_list = measurement_list
            
        if self.measurement_list is None and self.edi_fn is not None:
            self.get_measurement_lists()
            
        if self.measurement_list is None and self.edi_fn is None:
            print 'Nothing to read, check edi_fn or measurement_list attributes'
            return
       
        m_count = 1    
        for line in self.measurement_list:
            if type(line) is str:
                line_list = line.split('=')
                key = line_list[0].lower()
                value = line_list[1]
                if key in 'reflatitude':
                    key = 'reflat'
                    value = MTft._assert_position_format(key, value)
                elif key in 'reflongitude':
                    key = 'reflon'
                    value = MTft._assert_position_format(key, value)
                elif key in 'refelevation':
                    key = 'refelev'
                    try:
                        value = float(value)
                    except ValueError:
                        value = 0.0
                elif key in 'maxchannels':
                    key = 'maxchan'
                    value = int(value)
                elif key in 'maxmeasurements':
                    key = 'maxmeas'
                    value = int(value)
                setattr(self, key, value)
        
            elif type(line) is dict:
                try:
                    key = 'meas_{0:02.0f}'.format(float(line['id']))
                except KeyError:
                    key = 'meas_{0:02}'.format(m_count)
                value = Measurement(**line)
                
                setattr(self, key, value)
       
class Measurement(object):
    """
    class to put the MT measurements in
    """
    
    def __init__(self, **kwargs):
        
        self._kw_list = ['id', 'chtype', 'x', 'y', 'azm', 'acqchan']
        for key in self._kw_list:
            setattr(self, key, None)
        
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
        
        
#==============================================================================
# Test            
#==============================================================================
#fn = r"d:\Peacock\MTData\EDI_Files\mb018.edi"
fn = r"c:\Users\jpeacock\Documents\SaudiArabia\edi_files_fixed_lon\101_rr.edi"
#h = Header(edi_fn=fn)
#h.read_header()

m = DefineMeasurement(edi_fn=fn)

#print ''.join(f.write_info())


#for key in sorted(h.__dict__.keys()):
#    print key, h.__dict__[key]
    