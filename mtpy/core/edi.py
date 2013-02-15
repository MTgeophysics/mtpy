#!/usr/bin/env python

"""
mtpy/mtpy/core/edi.py

Contains classes and functions for handling EDI files. 
 
    Class:
    Edi - an edi object, containing all information from or for an edi file. Sections of edi files are given as respective attributes, their information are stored as dictionaries.

        Methods:
        - readfile()
        - writefile()
        - validate()
        - z2resphase()
        - rotate()
        - set/get_head()
        - set/get_info()
        - set/get_z()
        - set/get_tipper()
        - set/get_definemeas()
        - set/get_mtsect()
        - set/get_datacomponent()
        - set/get_frequencies()
        - set/get_edi_dict()
        - set/get_zrot()

    Functions:
    - read_edifile()
    - write_edifile()
    - combine_edifiles()
    - validate_edifile()
    - rotate_edifile()




@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np
import os
import sys
import os.path as op

import mtpy.utils.exceptions as MTexceptions
reload(MTexceptions)


#=================================================================



class Edi(object):
    """
        Edi class - generates an edi object.

        Methods  include reading and writing from and to edi-files, combination of edi-files, as well as 'get' and 'set' for all edi file sections

    """



    def __init__(self, fn = None):
    
            self.filename = fn
            if fn != None:
                if op.isfile(op.abspath(fn)):
                    self.filename = op.abspath(fn)
                else:
                    self.filename = None

            self.raw_filestring = None
            self.edi_dict = {}
            self.head = {}
            self.info_string = None
            self.info_dict = {}
            self.definemeas_dict = {}
            self.hmeas_emeas = None
            self.mtsect = {}
            self.freq = None
            self.n_freqs = 0.
            self.zrot = None
            self.data = {}
            self.z = {}
            self.tipper = None
            self.rho = None
            self.phase = None
            self.frequencies = None
            self.periods = None
            self.data = {}

    def readfile(self, fn):
        
        infile = op.abspath(fn)


        if not op.isfile(infile):
            raise MTexceptions.MTpyError_edi_file('File is not existing: %s'%infile)

        with open(infile,'r') as F:
            edistring = F.read()

        if not self._validate_edifile_string(edistring):
            raise MTexceptions.MTpyError_edi_file('%s is no proper edi file'%infile)

        self.filename = infile

        if 1:
            self._read_head(edistring)
        #except:
        #    raise MTexceptions.MTpyError_edi_file('Could not read HEAD section: %s'%fn)

        try:
            self._read_info(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read INFO section:%s'%fn)

        try:
            self._read_definemeas_dict(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read DEFINEMEAS section:%s'%fn)

        try:
            self._read_hmeas_emeas(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read HMEAS/EMEAS sub-section:%s'%fn)

        try:
            self._read_mtsect(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read MTSECT section:%s'%fn)

        try:
            self._read_freq(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read FREQ section:%s'%fn)

        try:
            self._read_z(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read Z section:%s'%fn)

        try:
            self._read_tipper(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read Tipper section:%s'%fn)

        try:
            self._read_zrot(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read Zrot section:%s'%fn)



    def _read_head(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'HEAD')
        except:
            raise

        head_dict = {}
        t1 = temp_string.strip().split('\n')
        t2 = [i.strip() for i in t1 if '=' in i]
        for j in t2:
            k = j.split('=')
            key = str(k[0]).lower()
            value = k[1].replace('"','')           
            head_dict[key] = value

        self.head = head_dict

    def _read_info(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'INFO')
        except:
            raise

        self.info_string = temp_string.strip()

        info_dict = {}

        t1 = temp_string.strip().split('\n')
        t2 = [i.strip() for i in t1 if '=' in i or ':' in i]

        for tmp_str in t2:
            #fill dictionary 
            #ignore lines with no information after '='' or ':'

            if '=' in tmp_str:
                t3 = tmp_str.split('=')
                key = str(t3[0]).lower()
                value = t3[1].replace('"','')           
                if not len(value) == 0:
                    info_dict[key] = value
            
            elif ':' in tmp_str:
                #consider potential ':' characters in coordinates!
                t3 = tmp_str.split(':')
                key = str(t3[0]).lower()
                value = t3[1:]
                value = [i.strip().replace('"','') for i in value]
                print value
                if len(value) > 1:
                    value = ':'.join(value)
                else:
                    value = value[0]

                if not len(value) == 0:
                    info_dict[key] = value

        self.info_dict = info_dict



    def _read_definemeas_dict(self, edistring):

        pass



    def _read_hmeas_emeas(self, edistring):

        pass


    def _read_mtsect(self, edistring):

        pass


    def _read_freq(self, edistring):

        pass


    def _read_z(self, edistring):

        pass

    def _read_tipper(self, edistring):

        pass

    def _read_zrot(self, edistring):

        pass



    def writefile(self,fn):
        pass

    def _validate_edifile_string(self, edistring):
        """
            Read the file as string and check, if blocks 'HEAD, INFO, =DEFINEMEAS, =MTSECT, FREQ, Z, END' are present.

            Within the blocks look for mandatory entries:
            HEAD: 'DATAID'
            INFO: None
            DEFINEMEAS: subblocks 'HMEAS, EMEAS' 
                        ('REFLAT, REFLONG, REFELEV' have to be present for measured data though)
            MTSECT: 'NFREQ'
            FREQ: non empty list
            Z: all components xx, yy, xy, yx ; real, imag and var ; each containing a non-empty list

        """
        isvalid = False
        found = 1

        #adding 1 to position of find to correct for possible occurrence at position 0 )
        found *= np.sign(edistring.upper().find('>HEAD') + 1 )
        found *= np.sign(edistring.upper().find('DATAID') + 1 )
        found *= np.sign(edistring.upper().find('>HMEAS') + 1 )
        found *= np.sign(edistring.upper().find('>EMEAS') + 1 )
        found *= np.sign(edistring.upper().find('NFREQ') + 1 )
        found *= np.sign(edistring.upper().find('>FREQ') + 1 )
        found *= np.sign(edistring.upper().find('>END') + 1 )
        found *= np.sign(edistring.upper().find('>=DEFINEMEAS') + 1 )
        found *= np.sign(edistring.upper().find('>=MTSECT') + 1 )


        if found < 1 :
            print 'Could not find all mandatory sections for a valid EDI file!\n (Most basic version must contain: "HEAD, INFO, =DEFINEMEAS, =MTSECT, FREQ, Z, END") '
            return False


        compstrings = ['ZXX','ZXY','ZYX','ZYY']
        Z_entries = ['R','I','.VAR']
        
        for comp in compstrings:
            for zentry in Z_entries:
                searchstring = '>'+comp+zentry
                z_comp_start_idx = edistring.upper().find(searchstring)
                found *= np.sign(z_comp_start_idx + 1 )
                #checking for non empty value list:
                next_block_start = edistring.upper().find('>',z_comp_start_idx+1)
                string_dummy_1 = edistring[z_comp_start_idx:next_block_start]
                lo_string_dummy_1 = string_dummy_1.strip().split()
                n_numbers = 0 
                for i in lo_string_dummy_1:
                    try:
                        n = float(i)
                        n_numbers +=1
                    except:
                        continue

                if n_numbers == 0:
                    print  MTexceptions.MTpyError_edi_file('Error in %s block: no values found'%(comp+zentry))

                    found *= 0 


        #checking for non empty frequency list:
        freq_start_idx = edistring.upper().find('>FREQ')
        next_block_start = edistring.upper().find('>',freq_start_idx + 1)
        string_dummy_2 = edistring[freq_start_idx:next_block_start]
        lo_string_dummy_2 = string_dummy_1.strip().split()
        #check, if there are actually one/some valid numbers:
        n_numbers = 0 
        for i in lo_string_dummy_2:
            try:
                n = float(i)
                n_numbers +=1
            except:
                continue

        if n_numbers == 0:
            print  MTexceptions.MTpyError_edi_file('Error in FREQ block: no frequencies found')

            found *= 0 


        if found > 0: isvalid = True

        return isvalid    


    def z2resphase(self):
        pass

    def rotate(self,angle):
        pass
        

    def get_head():
        pass
        

    def set_head():
        pass
        
      

    def get_info():
        pass
        
    
    def set_info():
        pass
        

    
    def get_z():
        pass
        

    def set_z():
        pass
        


    def get_tipper():
        pass
        
    
    def set_tipper():
        pass
        


    def get_definemeas():
        pass
        
    
    def set_definemeas():
        pass
        


    def get_mtsect():
        pass
        
    
    def set_mtsect():
        pass
        


    def get_datacomponent():
        pass
        
    
    def set_datacomponent():
        pass
        


    def get_frequencies():
        pass
        
   
    def set_frequencies():
        pass
        


    def get_edi_dict():
        pass
        
    
    def set_edi_dict():
        pass
        


    def get_zrot():
        pass
        
    
    def set_zrot():
        pass
        



def read_edifile(fn):
    pass

    return edi_object


def write_edifile(out_fn = None):
    pass

    return out_filename


def combine_edifiles(fn1, fn2, out_fn = None):
    pass

    return out_filename


def validate_edifile(fn):
    is_edi  = False

    return is_edi

def rotate_edifile(fn, out_fn = None):
    pass

    return out_filename


def _cut_sectionstring(edistring,sectionhead):

    start_idx = edistring.upper().find('>'+sectionhead.upper())
    if start_idx == -1:
        start_idx = edistring.upper().find('>='+sectionhead.upper())
        if start_idx == -1:
            raise 
        #correct for the = character
        start_idx += 1
    #start cut behind the section keyword
    start_idx += (1+len(sectionhead))


    next_block_start = edistring.upper().find('>', start_idx + 1)
    cutstring = edistring[start_idx:next_block_start]

    if len(cutstring) == 0 :
        raise


    return cutstring



