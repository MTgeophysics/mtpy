# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 15:47:52 2017

@author: jrpeacock
"""

#==============================================================================
import numpy as np
import os
from datetime import datetime

import mtpy.core.z as mtz
import mtpy.utils.configfile as mtcfg
import mtpy.utils.filehandling as mtfh
import mtpy.utils.exceptions as mtex
import mtpy.core.edi as mtedi

#==============================================================================
# Class to read j_file
#==============================================================================
class JFile(object):
    """
    be able to read and write a j-file
    """
    
    def __init__(self, j_fn=None):
        self.j_lines = j_fn
        self.header_dict = None
        self.metadata_dict = None
        self.Z = None
        self.Tipper = None
        
        if self.j_fn is not None:
            self.read_j_file(self.j_fn)

        
    def _get_j_lines(self):
        """
        read in the j_file as a list of lines, put the lines in attribute
        _j_lines
        """
        if self.j_fn is None:
            print 'j_fn is None'
            return
            
        if os.path.isfile(os.path.abspath(self.j_fn)) is False:
            raise IOError('Could not find {0}, check path'.format(self.j_fn))
        
        self._validate_j_file()
        
        with open(self.j_fn, 'r') as fid:
            self._j_lines = fid.readlines()
        print 'read in {0}'.format(self.j_fn)
        
    def _validate_j_file(self):
        """
        change the lat, lon, elev lines to something machine readable,
        if they are not.
        """
        
        # need to remove any weird characters in lat, lon, elev
        with open(self.j_fn, 'r') as fid:
            j_str = fid.read()
            
        # change lat
        j_str = self._rewrite_line('latitude', j_str)
        
        # change lon
        j_str = self._rewrite_line('longitude', j_str)
        
        # change elev
        j_str = self._rewrite_line('elevation', j_str)
        
        with open(self.j_fn, 'w') as fid:
            fid.write(j_str)
            
        print 'rewrote j-file {0} to make lat, lon, elev float values'.format(self.j_fn)
        
    def _get_str_value(self, string):

        value = string.split('=')[1].strip()
        try:
            value = float(value)
        except ValueError:
            value = 0.0
            
        return value
        
    def _rewrite_line(self, variable, file_str):
        variable_str = '>'+variable.upper()
        index_begin = file_str.find(variable_str)
        index_end = index_begin+file_str[index_begin:].find('\n')
        
        value = self._get_str_value(file_str[index_begin:index_end])
        print 'Changed {0} to {1}'.format(variable.upper(), value)
       
        new_line = '{0} = {1:<.2f}'.format(variable_str, value)
        file_str = file_str[0:index_begin]+new_line+file_str[index_end:]        
         
        return file_str
        
    def read_header(self):
        """
        Parsing the header lines of a j-file to extract processing information.
    
        Input:
        - j-file as list of lines (output of readlines())
    
        Output:
        - Dictionary with all parameters found

        """
            
        if self._j_lines is None:
            print "specify a file with jfile.j_fn = path/to/j/file"
            
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
        
        if self._j_lines is None:
            print "specify a file with jfile.j_fn = path/to/j/file"
        
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
        
    def read_j_file(self):
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
                        
        self._get_j_lines()

        self.read_header()
        self.read_metadata()       
        
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
            elif len(d_line.strip().split()) == 1 and 'r' not in d_line.lower():
                continue
            elif 'r' in d_line.lower():
                break
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
                        d_value = float(d_value)
                        # need to check for masked points represented by
                        # birrp as -999, apparently
                        if d_value == -999 or np.isnan(d_value):
                            d_value_list[d_index] = 0.0
                        else:
                            d_value_list[d_index] = d_value
                    except ValueError:
                        d_value_list[d_index] = 0.0
                
                # put the numbers in the correct dictionary as:
                # key = period, value = [real, imaginary, error]
                if d_key in z_index_dict.keys():
                    z_dict[d_key][d_value_list[0]] = d_value_list[1:4]
                elif d_key in t_index_dict.keys():
                    t_dict[d_key][d_value_list[0]] = d_value_list[1:4]
        
        # --> now we need to get the set of periods for all components  
        # check to see if there is any tipper data output          

        all_periods = []            
        for z_key in z_index_dict.keys():
            for f_key in z_dict[z_key].keys():
                all_periods.append(f_key)
        
        if len(t_dict['tzx'].keys()) == 0:
            print 'Could not find any Tipper data in {0}'.format(self.j_fn)
            find_tipper = False
  
        else:
            for t_key in t_index_dict.keys():
                for f_key in t_dict[t_key].keys():
                    all_periods.append(f_key)
            find_tipper = True
        
        all_periods = np.array(sorted(list(set(all_periods))))
        all_periods = all_periods[np.nonzero(all_periods)]
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
                try:
                    z_value = z_dict[z_key][per][0]+1j*z_dict[z_key][per][1]
                    z_arr[p_index, kk, ll] = z_value
                    z_err_arr[p_index, kk, ll] = z_dict[z_key][per][2]
                except KeyError:
                    print 'No value found for period {0:.4g}'.format(per)
                    print 'For component {0}'.format(z_key)
            if find_tipper is True:
                for t_key in sorted(t_index_dict.keys()):
                    kk = t_index_dict[t_key][0]
                    ll = t_index_dict[t_key][1]
                    try:
                        t_value = t_dict[t_key][per][0]+1j*t_dict[t_key][per][1]
                        t_arr[p_index, kk, ll] = t_value
                        t_err_arr[p_index, kk, ll] = t_dict[t_key][per][2]
                    except KeyError:
                        print 'No value found for period {0:.4g}'.format(per)
                        print 'For component {0}'.format(t_key)
        
        # put the results into mtpy objects
        freq = 1./all_periods    
        self.Z = mtz.Z(z_arr, z_err_arr, freq)
        self.Tipper = mtz.Tipper(t_arr, t_err_arr, freq)    
