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
        self.define_measurement = DefineMeasurement(edi_fn=self.edi_fn)
        self.data_sect = DataSection(edi_fn=self.edi_fn)
        self.Z = MTz.Z()
        self.Tipper = MTz.Tipper()
        
        if self.edi_fn is not None:
            self.read_data()
        
        
    def read_data(self):
        """
        read either impedance or spectra data
        """
        
        if self.edi_fn is None:
            raise MTex.MTpyError_EDI('No edi file input, check edi_fn')
        if os.path.isfile(self.edi_fn) is False:
            raise MTex.MTpyError_EDI('No edi file input, check edi_fn')
            
        with open(self.edi_fn, 'r') as fid:
            lines = fid.readlines()[self.data_sect.line_num+2:]
        
        if self.data_sect.data_type == 'spectra':
            self._read_spectra(lines)
        
        elif self.data_sect.data_type == 'z':
            self._read_mt(lines)
            
    def _read_mt(self, data_lines):
        """
        read in impedance and tipper data if its there
        """
        data_dict = {}
        data_find = False
        for line in data_lines:
            if line.find('>') == 0 and line.find('!') == -1:
                line_list = line[1:].strip().split()
                key = line_list[0].lower()
                if key[0] == 'z' or key[0] == 't' or key == 'freq':
                    data_find = True
                    data_dict[key] = []
                else:
                    data_find = False
                
        
            elif data_find == True and line.find('>') == -1 and line.find('!') == -1:
                d_lines = line.strip().split()
                for ii, dd in enumerate(d_lines):
                    # check for empty values and set them to 0, check for any
                    # other characters sometimes there are ****** for a null
                    # component
                    try:
                        d_lines[ii] = float(dd)
                        if d_lines[ii] == 1.0e32:
                            d_lines[ii] = 0.0
                    except ValueError:
                        d_lines[ii] = 0.0
                data_dict[key] += d_lines
        
        ## fill useful arrays
        freq_arr = np.array(data_dict['freq'], dtype=np.float)
        
        ## fill impedance tensor
        self.Z.freq = freq_arr.copy()
        self.Z.z = np.zeros((self.data_sect.nfreq, 2, 2), dtype=np.complex)
        self.Z.zerr = np.zeros((self.data_sect.nfreq, 2, 2), dtype=np.float)
        
        self.Z.z[:, 0, 0] = np.array(data_dict['zxxr'])+\
                             np.array(data_dict['zxxi'])*1j
        self.Z.z[:, 0, 1] = np.array(data_dict['zxyr'])+\
                            np.array(data_dict['zxyi'])*1j
        self.Z.z[:, 1, 0] = np.array(data_dict['zyxr'])+\
                            np.array(data_dict['zyxi'])*1j
        self.Z.z[:, 1, 1] = np.array(data_dict['zyyr'])+\
                            np.array(data_dict['zyyi'])*1j
        
        self.Z.zerr[:, 0, 0] = np.array(data_dict['zxx.var'])
        self.Z.zerr[:, 0, 1] = np.array(data_dict['zxy.var'])
        self.Z.zerr[:, 1, 0] = np.array(data_dict['zyx.var'])
        self.Z.zerr[:, 1, 1] = np.array(data_dict['zyy.var'])

        
        ## fill tipper data if there it exists
        self.Tipper.tipper = np.zeros((self.data_sect.nfreq, 1, 2), 
                                      dtype=np.complex) 
        self.Tipper.tippererr = np.zeros((self.data_sect.nfreq, 1, 2),
                                         dtype=np.float) 
        self.Tipper.freq = freq_arr.copy()

        if 'txr.exp' in data_dict.keys():
            self.Tipper.tipper[:, 0, 0] = np.array(data_dict['txr.exp'])+\
                                            np.array(data_dict['txi.exp'])*1j
            self.Tipper.tipper[:, 0, 1] = np.array(data_dict['tyr.exp'])+\
                                            np.array(data_dict['tyi.exp'])*1j
            
            self.Tipper.tippererr[:, 0, 0] = np.array(data_dict['txvar.exp'])    
            self.Tipper.tippererr[:, 0, 1] = np.array(data_dict['tyvar.exp'])
              
        else:
            print 'Could not find any Tipper data.'
            
    def _read_spectra(self, data_lines):
        """
        read in spectra data
        """
        
        data_dict = {}
        
        
                        

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
                    if getattr(self, key) is not None:
                        value = '{0}, {1}'.format(getattr(self, key), value)
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
            self.read_header(header_list)
            
        if self.header_list is None and self.edi_fn is not None:
            self.get_header_list()
            
        header_lines = ['>HEAD\n\n']
        for key in sorted(self._header_keys):
            value = getattr(self, key)
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
        self.maxmeas = 7
        self.maxrun = 999
        self.refelev = None
        self.reflat = None
        self.reflon = None
        self.reftype = 'cartesian'
        self.units = 'm'
        
        self._define_meas_keys = ['maxchan',
                                  'maxrun',
                                  'maxmeas',
                                  'reflat',
                                  'reflon',
                                  'refelev',
                                  'reftype',
                                  'units']
                                  
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
                    value = MTft._assert_position_format('lat', value)
                elif key in 'reflongitude':
                    key = 'reflon'
                    value = MTft._assert_position_format('lon', value)
                elif key in 'refelevation':
                    key = 'refelev'
                    value = MTft._assert_position_format('elev', value)
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
                if line['chtype'].lower().find('h') >= 0:
                    value = HMeasurement(**line)
                elif line['chtype'].lower().find('e') >= 0:
                    value = EMeasurement(**line)
                setattr(self, key, value)
                
    def write_define_measurement(self, measurement_list=None):
        """
        write the define measurement block as a list of strings
        """
        
        if measurement_list is not None:
            self.read_define_measurement(measurement_list=measurement_list)
            
        measurement_lines = ['>=DEFINEMEAS\n\n']
        for key in self._define_meas_keys:
            value = getattr(self, key)
            if key == 'reflat' or key == 'reflon':
                value = MTft.convert_dms_tuple2string(
                                        MTft.convert_degrees2dms_tuple(value))
            elif key == 'refelev':
                value = '{0:.3f}'.format(value)
            
            measurement_lines.append('{0}{1}={2}\n'.format(tab,
                                                           key.upper(),
                                                           value))
        measurement_lines.append('\n')
                                                           
        ## need to write the >XMEAS type
        m_key_list = [kk for kk in self.__dict__.keys() if kk.find('meas_')==0]
        if len(m_key_list) == 0:
            print 'No XMEAS information.'
        else:
            for key in sorted(m_key_list):
                m_obj = getattr(self, key)
                if m_obj.chtype.lower().find('h') >= 0:
                    head = 'hmeas'
                elif m_obj.chtype.lower().find('e') >= 0:
                    head = 'emeas'                
                else:
                    head = None
                
                m_list = ['>{0}'.format(head.upper())]
                for mkey, mfmt in zip(m_obj._kw_list, m_obj._fmt_list):
                    print mkey, mfmt, getattr(m_obj, mkey)
                    m_list.append(' {0}={1:{2}}'.format(mkey.upper(),
                                                        getattr(m_obj, mkey),
                                                        mfmt))
                m_list.append('\n')
                measurement_lines.append(''.join(m_list))
        
        return measurement_lines
            
#==============================================================================
# magnetic measurements
#==============================================================================
class HMeasurement(object):
    """
    class to put the MT measurements in
    """
    
    def __init__(self, **kwargs):
        
        self._kw_list = ['id', 'chtype', 'x', 'y', 'azm', 'acqchan']
        self._fmt_list = ['<4.4g','<3', '<4.1f', '<4.1f', '<4.1f', '<4']
        for key in self._kw_list:
            setattr(self, key, None)
        
        for key in kwargs.keys():
            try:
                setattr(self, key, float(kwargs[key]))
            except ValueError:
                setattr(self, key, kwargs[key])

#==============================================================================
# electric measurements            
#==============================================================================
class EMeasurement(object):
    """
    class to put the MT measurements in
    """
    
    def __init__(self, **kwargs):
        
        self._kw_list = ['id', 'chtype', 'x', 'y', 'x2', 'y2', 'acqchan']
        self._fmt_list = ['<4.4g', '<3', '<4.1f', '<4.1f', '<4.1f', '<4.1f',
                          '<4']
        for key in self._kw_list:
            setattr(self, key, None)
        
        for key in kwargs.keys():
            try:
                setattr(self, key, float(kwargs[key]))
            except ValueError:
                setattr(self, key, kwargs[key])
        
        
#==============================================================================
# data section        
#==============================================================================
class DataSection(object):
    """
    read the data section (XSECT part)
    """
    def __init__(self, edi_fn=None):
        self.edi_fn = edi_fn
        
        self.data_type = 'z'
        self.line_num = 0
        self.data_sect_list = None
        
        self._kw_list = ['ex', 
                         'ey',
                         'hx',
                         'hy',
                         'hz',
                         'nfreq',
                         'sectid',
                         'nchan',
                         'maxblks']
                         
        for key in self._kw_list:
            setattr(self, key, None)
            
        if self.edi_fn is not None:
            self.read_data_sect()
        
        
    def get_data_sect(self):
        """
        read in the data of the file, will detect if reading spectra or 
        impedance.
        """
        
        if self.edi_fn is None:
            raise MTex.MTpyError_EDI('No edi file to read. Check edi_fn')
            
        if os.path.isfile(self.edi_fn) is False:
            raise MTex.MTpyError_EDI('Could not find {0}. Check path'.format(self.edi_fn))
        
        self.data_sect_list = []
        data_sect_find = False
        count = 0
        with open(self.edi_fn) as fid:
            for ii, line in enumerate(fid):
                if line.find('>=') == 0:
                    count += 1
                    if line.lower().find('sect') > 0:
                        data_sect_find = True
                        self.line_num = ii
                        if line.lower().find('spect') > 0:
                            self.data_type = 'spectra'
                        elif line.lower().find('mt') > 0:
                            self.data_type = 'z'
                    else:
                        data_sect_find = False
                    if count > 2 and data_sect_find == False:
                        break
                elif count == 2 and line.find('>') != 0 and \
                    data_sect_find == True:
                    if len(line.strip()) > 2:
                        self.data_sect_list.append(line.strip())
                        
    def read_data_sect(self, data_sect_list=None):
        """
        read data section
        """
        
        if data_sect_list is not None:
            self.data_sect_list = data_sect_list
            
        if self.edi_fn is not None and self.data_sect_list is None:
            self.get_data_sect()
            
        for d_line in self.data_sect_list:
            d_list = d_line.split('=')
            if len(d_list) > 1:
                key = d_list[0].lower()
                try:
                    value = int(d_list[1].strip())
                except ValueError:
                    value = d_list[1].strip().replace('"', '')
            
                setattr(self, key, value)
                
    def write_data_sect(self, data_sect_list=None):
        """
        write a data section
        """
        
        if data_sect_list is not None:
            self.read_data_sect(data_sect_list)
            
        if self.data_type == 'spectra':
            data_sect_lines = ['\n>=spectrasect\n'.upper()]
            
        if self.data_type == 'z':
            data_sect_lines = ['\n>=mtsect\n'.upper()]
        for key in self._kw_list:
            data_sect_lines.append('{0}{1}={2}\n'.format(tab, 
                                                         key, 
                                                         getattr(self, key)))
        
        return data_sect_lines
                        

#==============================================================================
# Test            
#==============================================================================
#fn = r"d:\Peacock\MTData\EDI_Files\mb018.edi"
fn = r"c:\Users\jpeacock\Documents\ShanesBugs\Jess\EDI_files\104A.edi"
#fn = r"c:\Users\jpeacock\Documents\SaudiArabia\edi_files_fixed_lon\101_rr.edi"
#h = Header(edi_fn=fn)
#h.read_header()

#m = DefineMeasurement(edi_fn=fn)
#
#print ''.join(m.write_define_measurement())

edi_obj = Edi(edi_fn=fn)

#d = DataSection(edi_fn=fn)
#for key in sorted(h.__dict__.keys()):
#    print key, h.__dict__[key]
    