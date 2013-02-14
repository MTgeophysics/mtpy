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
    
        if fn != None:
            if op.isfile(op.abspath(fn)):
                self.filename = op.abspath(fn)
            else:
                self.filename = None

            self.raw_filestring = None
            self.edi_dict = {}
            self.head = {}
            self.info = {}
            self.definemeas = {}
            self.mtsect = {}
            self.freq = {}
            self.n_freqs = 0.
            self.zrot = None
            self.data = {}
            self.z = None
            self.tipper = None
            self.rho = None
            self.phase = None


    def readfile(self, fn):
        
        infile = op.abspath(fn)


        if not op.isfile(infile):
            raise MTexceptions.MTpyError_edi_file('File is not existing: %s'%infile)

        with open(infile,'r') as F:
            edistring = F.read()

        if not self._validate_edifile_string(edistring):
            raise MTexceptions.MTpyError_edi_file('%s is no proper edi file'%infile)




    def writefile(self,fn):
        pass

    def _validate_edifile_string(self, edistring):
        """
            Read the file as string and check, if blocks 'HEAD, INFO, DEFINEMEAS, MTSECT, FREQ, Z' are present.

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

        if found < 1 :
            print 'Could not find all mandatory sections for a valid EDI file!'
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




