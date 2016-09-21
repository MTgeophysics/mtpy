# -*- coding: utf-8 -*-
"""
BIRRP
===========
    * deals with inputs and outputs from BIRRP

Created on Tue Sep 20 14:33:20 2016

@author: jrpeacock
"""

#==============================================================================
import numpy as np
import os

import mtpy.core.z as mtz

#==============================================================================
class BIRRP_Parameters(object):
    """
    class to hold and produce the appropriate parameters given the input 
    parameters.
    """
    
    def __init__(self, ilev=0, **kwargs):
        # set the input level
        self.ilev = ilev

        # get default parameters        
        self._get_parameters()
        
        # set any attributes that are given
        for key in kwargs.keys():
            setattr(self, key, kwargs[key])
            
        self._validate_parametrs()
        
    def _get_parameters(self):
        """
        get appropriate parameters
        """
        
        if self.ilev not in [0, 1]:
            raise ValueError('Input level, ilev, must be 0 for basic or 1 d'+ 
                             'for advanced, not {0}'.format(self.ilev))
        
        if self.ilev == 0:
            self.nout = 3
            self.ninp = 2
            self.tbw = 2.0
            self.nfft = 2**18
            self.nsctmax = 14
            self.uin = 0
            self.ainuin = .9999
            self.c2threshe = 0
            self.nz = 0
            self.c2threshe1 = 0
            self.ofil = 'mt'
            self.nlev = 0
            self.nar = 5
            self.imode = 0
            self.jmode = 0
            self.nfil = 0
            self.thetae = [0, 90, 0]
            self.thetab = [0, 90, 0]
            self.thetaf = [0, 90, 0]
            
        elif self.ilev == 1:
            self.nout = 3
            self.ninp = 2
            self.nref = 2
            self.nrr = 1
            self.tbw = 2
            self.nfft = 2**18
            self.nsctinc = 2
            self.nsctmax = np.floor(np.log2(self.nfft))-4
            self.nf1 = self.tbw+2
            self.nfinc = self.tbw
            self.nfsect = 2
            self.uin = 0
            self.ainlin = 0.0001
            self.ainuin = 0.9999
            if self.nrr == 1:
                self.c2threshb = 0.0
                self.c2threshe = 0.0
                if self.c2threshe == 0 and self.c2threshb == 0:
                    self.nz = 0
                else:
                    self.nz = 0
                    self.perlo = 1000
                    self.perhi = .0001
            elif self.nrr == 0:
                self.c2threshb = 0
                self.c2threshe = 0
            self.nprej = 0
            self.prej = None
            self.c2threshe1 = 0
            self.ofil = 'mt'
            self.nlev = 0
            self.nar = 5
            self.imode = 0
            self.jmode = 0
            self.nfil = 0
            self.thetae = [0, 90, 0]
            self.thetab = [0, 90, 0]
            self.thetaf = [0, 90, 0]
            
    def _validate_parametrs(self):
        """
        check to make sure the parameters are legit.
        """
        
        # be sure the 
        if self.ninp not in [1, 2, 3]:
            print 'Number of inputs {0} not allowed.'.format(self.ninp)
            self.ninp = 2
            print '  --> setting ninp to {0}'.format(self.ninp)            
        
        if self.nout not in [2, 3]:
            print 'Number of outputs {0} not allowed.'.format(self.nout)
            self.nout = 2
            print '  --> setting nout to {0}'.format(self.nout)
        
        if self.tbw < 0 or self.tbw > 4:
            print 'Total bandwidth of slepian window {0} not allowed.'.format(self.tbw)
            self.tbw = 2
            print '  --> setting tbw to {0}'.format(self.tbw)
            
        if np.remainder(np.log2(self.nfft), 2) != 0:
            print 'Window length nfft should be a power of 2 not {0}'.format(self.nfft)
            self.nfft = 2**np.floor(np.log2(self.nfft))
            print '  -- > setting nfft to {0}, (2**{1:.0f})'.format(self.nfft,
                                                              np.log2(self.nfft))
                                                              
        if np.log2(self.nfft)-self.nsctmax < 4:
            print 'Maximum number of windows {0} is too high'.format(self.nsctmax)
            self.nsctmax = np.log2(self.nfft)-4
            print '  --> setting nsctmax to {0}'.format(self.nsctmax)
            
        if self.uin != 0:
            print 'You\'re playing with fire if uin is not 0.'
            self.uin = 0
            print '  --> setting uin to 0, if you don\'t want that change it back'
        
        if self.imode not in [0, 1, 2, 3]:
            raise BIRRP_Parameter_Error('Invalid number for time series mode,'
                                       'imode, {0}, should be 0, 1, 2, or 3'.format(self.imode))
        
        if self.jmode not in [0, 1]:
            raise BIRRP_Parameter_Error('Invalid number for time mode,'
                                       'imode, {0}, should be 0, or 1'.format(self.imode))
        if self.ilev == 1:
            if self.nsctinc != 2:
                print '!WARNING! Check decimation increment nsctinc, should be 2 not {0}'.format(self.nsctinc)
                
            if self.nfsect != 2:
                print 'Will get an error from BIRRP if nfsect is not 2.'
                print 'number of frequencies per section is {0}'.format(self.nfsect)
                self.nfsect = 2
                print '  --> setting nfsect to 2'
                
            if self.nf1 != self.tbw+2:
                print '!WARNING! First frequency should be around tbw+2.'
                print 'nf1 currently set to {0}'.format(self.nf1)
                
            if self.nfinc != self.tbw:
                print '!WARNING! sequence of frequencies per window should be around tbw.'
                print 'nfinc currently set to {0}'.format(self.nfinc) 
                
            if self.nprej != 0:
                if self.prej is None or type(self.prej) is not list:
                    raise BIRRP_Parameter_Error('Need to input a prejudice list if nprej != 0'+
                                     '\nInput as a list of frequencies' )
                                     
            if self.nrr not in [0, 1]:
                print('!WARNING! Value for picking remote reference or '+ 
                     'two stage processing, nrr, '+ 
                     'should be 0 or 1 not {0}'.format(self.nrr))
                self.nrr = 0
                print '  --> setting nrr to {0}'.format(self.nrr)
                
                
            if self.c2threshe != 0 or self.c2threshb != 0:
                if not self.perhi:
                    raise BIRRP_Parameter_Error('Need to input a high period (s) threshold as perhi')
        
                if not self.perlo:
                    raise BIRRP_Parameter_Error('Need to input a low period (s) threshold as perlo')
         
        if len(self.thetae) != 3:
            print 'Electric rotation angles not input properly {0}'.format(self.thetae)
            print 'input as north, east, orthogonal rotation'            
            self.thetae = [0, 90, 0]
            print '  --> setting thetae to {0}'.format(self.thetae)
        
        if len(self.thetab) != 3:
            print 'Magnetic rotation angles not input properly {0}'.format(self.thetab)
            print 'input as north, east, orthogonal rotation'            
            self.thetab = [0, 90, 0]
            print '  --> setting thetab to {0}'.format(self.thetab)
            
                
        if len(self.thetaf) != 3:
            print 'Fiedl rotation angles not input properly {0}'.format(self.thetaf)
            print 'input as north, east, orthogonal rotation'            
            self.thetaf = [0, 90, 0]
            print '  --> setting thetaf to {0}'.format(self.thetaf)

#==============================================================================
# Error classes            
#==============================================================================
class BIRRP_Parameter_Error(Exception):
    pass

#==============================================================================
# Class to read j_file
#==============================================================================
class JFile(object):
    """
    be able to read and write a j-file
    """
    
    def __init__(self, j_fn=None):
        self.j_fn = j_fn
        self.header_dict = None
        self.metadata_dict = None
        self._j_lines = None
        self.Z = None
        self.Tipper = None
        
    def _get_j_lines(self, j_fn=None):
        """
        read in the j_file as a list of lines, put the lines in attribute
        _j_lines
        """
        if j_fn is not None:
            self.j_fn = j_fn
            
        if os.path.isfile(os.path.abspath(self.j_fn)) is False:
            raise IOError('Could not find {0}, check path'.format(self.j_fn))
            
        with open(self.j_fn, 'r') as fid:
            self._j_lines = fid.readlines()
        
    def read_header(self, j_lines=None, j_fn=None):
        """
        Parsing the header lines of a j-file to extract processing information.
    
        Input:
        - j-file as list of lines (output of readlines())
    
        Output:
        - Dictionary with all parameters found

        """
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
            
        header_lines = [j_line for j_line in self._j_lines if '#' in j_line]
        header_dict = {'title':header_lines[0][1:].strip()}
        
        fn_count = 0
        theta_count = 0
        # put the information into a dictionary 
        for h_line in header_lines[1:]:
            # replace '=' with a ' ' to be sure that when split is called there is a
            # split, especially with filenames
            h_list = h_line[1:].strip().replace('=', ' ').split()
            # skip if there is only one element in the list
            if len(h_list) == 1:
                continue
            # get the key and value for each parameter in the given line
            for h_index in range(0, len(h_list), 2):
                h_key = h_list[h_index]
                # if its the file name, make the dictionary value be a list so that 
                # we can append nread and nskip to it, and make the name unique by
                # adding a counter on the end
                if h_key == 'filnam':
                    h_key = '{0}_{1:02}'.format(h_key, fn_count)
                    fn_count += 1
                    h_value = [h_list[h_index+1]]
                    header_dict[h_key] = h_value
                    continue
                elif h_key == 'nskip' or h_key == 'nread':
                    h_key = 'filnam_{0:02}'.format(fn_count-1)
                    h_value = int(h_list[h_index+1])
                    header_dict[h_key].append(h_value)
                    
                # if its the line of angles, put them all in a list with a unique key
                elif h_key == 'theta1':
                    h_key = '{0}_{1:02}'.format(h_key, theta_count)
                    theta_count += 1
                    h_value = float(h_list[h_index+1])
                    header_dict[h_key] = [h_value]
                elif h_key == 'theta2' or h_key == 'phi':
                    h_key = '{0}_{1:02}'.format('theta1', theta_count-1)
                    h_value = float(h_list[h_index+1])
                    header_dict[h_key].append(h_value)
                    
                else:
                    try:
                        h_value = float(h_list[h_index+1])
                    except ValueError:
                        h_value = h_list[h_index+1]
                    
                    header_dict[h_key] = h_value
            
        self.header_dict = header_dict
        
    def read_metadata(self, j_lines=None, j_fn=None):
        """
        read in the metadata of the station, or information of station 
        logistics like: lat, lon, elevation
        
        Not really needed for a birrp output since all values are nan's
        """
        
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
        
        metadata_lines = [j_line for j_line in self._j_lines if '>' in j_line]
    
        metadata_dict = {}
        for m_line in metadata_lines:
            m_list = m_line.strip().split('=')
            m_key = m_list[0][1:].strip().lower()
            try:
                m_value = float(m_list[0].strip())
            except ValueError:
                m_value = 0.0
                
            metadata_dict[m_key] = m_value
            
        self.metadata_dict = metadata_dict
        
    def read_j_file(self, j_lines=None, j_fn=None):
        """
        read_j_file will read in a *.j file output by BIRRP (better than reading lots of *.<k>r<l>.rf files)
    
        Input:
        j-filename
    
        Output: 4-tuple
        - periods : N-array
        - Z_array : 2-tuple - values and errors
        - tipper_array : 2-tuple - values and errors
        - processing_dict : parsed processing parameters from j-file header
    
        """   
    
        # read data
        z_index_dict = {'zxx':(0, 0),
                        'zxy':(0, 1),
                        'zyx':(1, 0),
                        'zyy':(1, 1)}
        t_index_dict = {'tzx':(0, 0),
                        'tzy':(0, 1)}
                        
        if j_lines is not None:
            self._j_lines = j_lines
            
        if j_fn is not None:
            self.j_fn = j_fn
            
        if self._j_lines is None:
            self._get_j_lines()
            
        self.header_dict = self.read_header()
        self.metadata_dict = self.read_metadata()
        
        data_lines = [j_line for j_line in self._j_lines 
                      if not '>' in j_line and not '#' in j_line][1:]
                
        # sometimes birrp outputs some missing periods, so the best way to deal with 
        # this that I could come up with was to get things into dictionaries with 
        # key words that are the period values, then fill in Z and T from there
        # leaving any missing values as 0
        
        # make empty dictionary that have keys as the component 
        z_dict = dict([(z_key, {}) for z_key in z_index_dict.keys()])
        t_dict = dict([(t_key, {}) for t_key in t_index_dict.keys()])
        for d_line in data_lines:
            # check to see if we are at the beginning of a component block, if so 
            # set the dictionary key to that value
            if 'z' in d_line.lower():
                d_key = d_line.strip().split()[0].lower()
            
            # if we are at the number of periods line, skip it
            elif len(d_line.strip().split()) == 1:
                continue
            # get the numbers into the correct dictionary with a key as period and
            # for now we will leave the numbers as a list, which we will parse later
            else:
                # split the line up into each number
                d_list = d_line.strip().split()
                
                # make a copy of the list to be sure we don't rewrite any values,
                # not sure if this is necessary at the moment
                d_value_list = list(d_list)
                for d_index, d_value in enumerate(d_list):
                    # check to see if the column number can be converted into a float
                    # if it can't, then it will be set to 0, which is assumed to be
                    # a masked number when writing to an .edi file
                    try:
                        d_value_list[d_index] = float(d_value)
                    except ValueError:
                        d_value_list[d_index] = 0.0
                
                # put the numbers in the correct dictionary as:
                # key = period, value = [real, imaginary, error]
                if d_key in z_index_dict.keys():
                    z_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                elif d_key in t_index_dict.keys():
                    t_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                    
        # now we need to get the set of periods for all components
        all_periods = sorted(list(set(np.array([z_dict[z_key].keys() for z_key in z_index_dict.keys()]+\
                                    [t_dict[t_key].keys() for t_key in t_index_dict.keys()]).flatten())))
        
        num_per = len(all_periods)
        
        # fill arrays using the period key from all_periods
        z_arr = np.zeros((num_per, 2, 2), dtype=np.complex)
        z_err_arr = np.zeros((num_per, 2, 2), dtype=np.float)
        
        t_arr = np.zeros((num_per, 1, 2), dtype=np.complex)
        t_err_arr = np.zeros((num_per, 1, 2), dtype=np.float)
        
        for p_index, per in enumerate(all_periods):
            for z_key in sorted(z_index_dict.keys()):
                kk = z_index_dict[z_key][0]
                ll = z_index_dict[z_key][1]
                z_value = z_dict[z_key][per][0]+1j*z_dict[z_key][per][1]
                z_arr[p_index, kk, ll] = z_value
                z_err_arr[p_index, kk, ll] = z_dict[z_key][per][2]
            for t_key in sorted(t_index_dict.keys()):
                kk = t_index_dict[t_key][0]
                ll = t_index_dict[t_key][1]
                t_value = t_dict[t_key][per][0]+1j*t_dict[t_key][per][1]
                t_arr[p_index, kk, ll] = t_value
                t_err_arr[p_index, kk, ll] = t_dict[t_key][per][2]
        
        # put the results into mtpy objects
        freq = 1./np.array(all_periods)    
        self.Z = mtz.Z(z_arr, z_err_arr, freq)
        self.Tipper = mtz.Tipper(t_arr, t_err_arr, freq)    
