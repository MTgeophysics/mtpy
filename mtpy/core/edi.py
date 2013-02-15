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
import math

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
            self.definemeas = {}
            self.hmeas_emeas = None
            self.mtsect = {}
            self.freq = None
            self.n_freqs = 0.
            self.zrot = None
            self.data = {}
            self.z_dict = {}
            self.z = None
            self.zerr = None
            self.tipper = None
            self.tippererr = None
            self.tipper_dict = {}
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
        self.raw_filestring = edistring


        try:
            self._read_head(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read HEAD section: %s'%infile)

        try:
            self._read_info(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read INFO section: %s'%infile)

        try:
            self._read_definemeas(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read DEFINEMEAS section: %s'%infile)

        try:
            self._read_hmeas_emeas(edistring)
        except:
            print 'Could not read HMEAS/EMEAS sub-section: %s'%infile

        try:
            self._read_mtsect(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read MTSECT section: %s'%infile)

        try:
            self._read_freq(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read FREQ section: %s'%infile)

        try:
            self._read_z(edistring)
        except:
            raise MTexceptions.MTpyError_edi_file('Could not read Z section: %s'%infile)

        try:
            self._read_tipper(edistring)
        except:
            print 'Could not read Tipper section: %s'%infile

        try:
            self._read_zrot(edistring)
        except:
            print 'Could not read Zrot section: %s'%infile

        #collect all data information in one dictionary
        data_dict = {}

        data_dict['z'] = self.z
        data_dict['tipper'] = self.tipper
        data_dict['zrot'] = self.zrot
        data_dict['frequencies'] = self.frequencies
        data_dict['zerr'] = self.zerr
        data_dict['tippererr'] = self.tippererr
        
        self.data = data_dict


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
                if len(value) > 1:
                    value = ':'.join(value)
                else:
                    value = value[0]

                if not len(value) == 0:
                    info_dict[key] = value

        self.info_dict = info_dict



    def _read_definemeas(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'DEFINEMEAS')
        except:
            raise

        d_dict = {} 

        t1 = temp_string.strip().split('\n')

        for tmp_str in t1:
            if '=' in tmp_str:
                k = tmp_str.strip().split('=')
                key = k[0].lower()
                value = k[1].replace('"','')
                if len(value) != 0:
                    d_dict[key] = value
         
        if len(d_dict.keys()) == 0:
            raise


        self.definemeas = d_dict



    def _read_hmeas_emeas(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'HMEAS_EMEAS')
        except:
            raise

        t1 = temp_string.strip().split('\n')
        lo_hmeas_emeas = []
        for j in t1:
            j = j.replace('>','')
            lo_j = j.split()
            lo_hmeas_emeas.append(tuple(lo_j))

        self.hmeas_emeas = lo_hmeas_emeas


    def _read_mtsect(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'MTSECT')
        except:
            raise
        m_dict = {}

        t1 = temp_string.strip().split('\n')

        for tmp_str in t1:
            if '=' in tmp_str:
                k = tmp_str.strip().split('=')
                key = k[0].lower()
                value = k[1].replace('"','')
                if len(value) != 0:
                    m_dict[key] = value
         
        if len(m_dict.keys()) == 0:
            raise


        self.mtsect = m_dict


    def _read_freq(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'FREQ')
        except:
            raise
        
        lo_freqs = []

        t1 = temp_string.strip().split('\n')[1:]

        for j in t1:
            lo_j = j.strip().split()
            for k in lo_j:
                try:
                    lo_freqs.append(float(k))
                except:
                    pass


        self.n_freqs = len(lo_freqs)
        self.frequencies = lo_freqs
        self.freq = lo_freqs
        self.periods = list(1./np.array(lo_freqs))


    def _read_z(self, edistring):
        """
        Read in impedances information. 
        Store it as dictionary and complex array (incl. Zvar values in 'zerr' array)

        """

        compstrings = ['ZXX','ZXY','ZYX','ZYY']
        Z_entries = ['R','I','.VAR']

        z_array = np.zeros((self.n_freqs,2,2),dtype=np.complex)
        zerr_array = np.zeros((self.n_freqs,2,2),dtype=np.float)
        z_dict = {}


        for idx_comp,comp in enumerate(compstrings):
            for idx_zentry,zentry in enumerate(Z_entries):
                sectionhead = comp + zentry
                try:
                    temp_string = _cut_sectionstring(edistring,sectionhead)
                except:
                    pass
  
                lo_z_vals = []

                #check, if correct number of entries are given in the block
                t0 = temp_string.strip().split('\n')[0]
                n_dummy = int(float(t0.split('//')[1].strip()))
                if not n_dummy == self.n_freqs:
                    raise


                t1 = temp_string.strip().split('\n')[1:]
                for j in t1:
                    lo_j = j.strip().split()
                    for k in lo_j:
                        try:
                            lo_z_vals.append(float(k))
                        except:
                            pass

                z_dict[sectionhead] = lo_z_vals

        self.z_dict = z_dict

        for idx_freq  in range( self.n_freqs):
            z_array[idx_freq,0,0] = np.complex(self.z_dict['ZXXR'][idx_freq], self.z_dict['ZXXI'][idx_freq])
            zerr_array[idx_freq,0,0] = self.z_dict['ZXX.VAR'][idx_freq]

            z_array[idx_freq,0,1] = np.complex(self.z_dict['ZXYR'][idx_freq], self.z_dict['ZXYI'][idx_freq])
            zerr_array[idx_freq,0,1] = self.z_dict['ZXY.VAR'][idx_freq]

            z_array[idx_freq,1,0] = np.complex(self.z_dict['ZYXR'][idx_freq], self.z_dict['ZYXI'][idx_freq])
            zerr_array[idx_freq,1,0] = self.z_dict['ZYX.VAR'][idx_freq]

            z_array[idx_freq,1,1] = np.complex(self.z_dict['ZYYR'][idx_freq], self.z_dict['ZYYI'][idx_freq])
            zerr_array[idx_freq,1,1] = self.z_dict['ZYY.VAR'][idx_freq]

        self.z = z_array
        self.zerr = zerr_array


    def _read_tipper(self, edistring):

        compstrings = ['TX','TY']
        T_entries = ['R','I','VAR']
    
        tipper_array = np.zeros((self.n_freqs,1,2),dtype=np.complex)
        tippererr_array = np.zeros((self.n_freqs,1,2),dtype=np.float)
        t_dict = {}


        for idx_comp,comp in enumerate(compstrings):
            for idx_tentry,tentry in enumerate(T_entries):

                try:
                    sectionhead = comp + tentry + '.EXP'
                    temp_string = _cut_sectionstring(edistring,sectionhead)
                except:
                    try:
                        sectionhead = comp + tentry
                        temp_string = _cut_sectionstring(edistring,sectionhead)
                    except:
                        pass
  
                lo_t_vals = []
                
                #check, if correct number of entries are given in the block
                t0 = temp_string.strip().split('\n')[0]
                n_dummy = int(float(t0.split('//')[1].strip()))
                if not n_dummy == self.n_freqs:
                    raise

                t1 = temp_string.strip().split('\n')[1:]
                for j in t1:
                    lo_j = j.strip().split()
                    for k in lo_j:
                        try:
                            lo_t_vals.append(float(k))
                        except:
                            pass

                t_dict[comp + tentry] = lo_t_vals

        self.tipper_dict = t_dict

        for idx_freq  in range( self.n_freqs):
            tipper_array[idx_freq,0,0] = np.complex(self.tipper_dict['TXR'][idx_freq], self.tipper_dict['TXI'][idx_freq])
            tippererr_array[idx_freq,0,0] = self.tipper_dict['TXVAR'][idx_freq]

            tipper_array[idx_freq,0,1] = np.complex(self.tipper_dict['TYR'][idx_freq], self.tipper_dict['TYI'][idx_freq])
            tippererr_array[idx_freq,0,1] = self.tipper_dict['TYVAR'][idx_freq]



        self.tipper = tipper_array
        self.tippererr = tippererr_array



    def _read_zrot(self, edistring):

        try:
            temp_string = _cut_sectionstring(edistring,'ZROT')
        except:
            lo_angles = list( np.zeros((self.n_freqs)) )            
            self.zrot = lo_angles
            return


        lo_angles = []

        t1 = temp_string.strip().split('\n')[1:]

        for j in t1:
            lo_j = j.strip().split()
            for k in lo_j:
                try:
                    lo_angles.append(float(k))
                except:
                    pass


        if len(lo_angles) != self.n_freqs:
            raise

        self.zrot = lo_angles




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
        """
            Rotate the Z and tipper information in the Edi object. Change the rotation angles in Zrot respectively.

            Rotation angle must be given in degrees. All angles are referenced to geographic North, positive in clockwise direction. (Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.

        """
        n_freqs = self.n_freqs
        lo_original_angles = self.zrot
        z = self.z
        tipper = self.tipper
        zerr = self.zerr
        angle = angle%360

        #1. rotation positive in clockwise direction
        #2. orientation of new X axis X' given by rotation angle
        #3. express contents of Z/tipper (points P) in this new system (points P')
        #4. rotation for points calculated as P' = ([cos , sin ],[-sin, cos]) * P <=> P' = R * P
        #5. => B' = R * B and E' = R * E
        # (Rt is the inverse rotation matrix)
        #6. E = Z * B => Rt * E' = Z * Rt * B' => E' = (R*Z*Rt) * B' => Z' = (R*Z*Rt)  

        #7. Bz = T * B => Bz = T * Rt * B' => T' = (T * Rt)

        # Rotation of the uncertainties:
        # a) rotate Z into Z'
        # b) use propagation of errors on Z' to obtain the rotated Z'err
        # That is NOT the same as the rotated error matrix Zerr (although the result is similar)

        z_rot = z.copy()
        zerr_rot = zerr.copy()

        for idx_freq in range(self.n_freqs):

            phi = math.radians(angle)

            cphi = np.cos(phi)
            sphi = np.sin(phi)

            z_orig = z[idx_freq,:,:]

            z_rot[idx_freq,0,0] = cphi**2 * z_orig[0,0] + cphi*sphi*(z_orig[0,1]+z_orig[1,0]) + sphi**2 * z_orig[1,1]
            z_rot[idx_freq,0,1] = cphi**2 * z_orig[0,1] + cphi*sphi*(z_orig[1,1]-z_orig[0,0]) - sphi**2 * z_orig[1,0]
            z_rot[idx_freq,1,0] = cphi**2 * z_orig[1,0] + cphi*sphi*(z_orig[1,1]-z_orig[0,0]) - sphi**2 * z_orig[0,1]
            z_rot[idx_freq,1,1] = sphi**2 * z_orig[0,0] - cphi*sphi*(z_orig[0,1]+z_orig[1,0]) + cphi**2 * z_orig[1,1]

            zerr_orig = zerr[idx_freq,:,:]

            zerr_rot[idx_freq,0,0] = np.sqrt( (cphi**2 * zerr_orig[0,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[0,1])**2 + (zerr_orig[1,0])**2) + (sphi**2 * zerr_orig[1,1])**2)

            zerr_rot[idx_freq,0,1] = np.sqrt( (cphi**2 * zerr_orig[0,1])**2 + cphi**2 * sphi**2 * ( (zerr_orig[1,1])**2 + (zerr_orig[0,0])**2) + (sphi**2 * zerr_orig[1,0])**2) 

            zerr_rot[idx_freq,1,0] = np.sqrt( (cphi**2 * zerr_orig[1,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[1,1])**2 + (zerr_orig[0,0])**2) + (sphi**2 * zerr_orig[0,1])**2) 

            zerr_rot[idx_freq,1,1] = np.sqrt( (sphi**2 * zerr_orig[0,0])**2 + cphi**2 * sphi**2 * ( (zerr_orig[0,1])**2 + (zerr_orig[1,0])**2) + (cphi**2 * zerr_orig[1,1])**2) 


        #todo: check !!!!
        #todo : tipper rotation


        self.z = z_rot
        self.zerr = zerr_rot


        self.zrot = list( np.zeros((len(self.zrot))) + angle) 
        

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
        




    def get_zrot():
        pass
        
    
    def set_zrot():
        pass



#end of Edi Class
#=========================


def read_edifile(fn):

    edi_object = Edi()

    edi_object.readfile(fn) 
   

    return edi_object


def write_edifile(out_fn = None):
    pass

    return out_filename


def combine_edifiles(fn1, fn2, out_fn = None):
    
    edi_object1 = Edi()
    edi_object1.readfile(fn1)
    edi_object2 = Edi()
    edi_object2.readfile(fn2)



    return out_filename


def validate_edifile(fn):

    edi_object = Edi()

    try:
        edi_object.readfile(fn) 
        return True
    except:
        return False


def rotate_edifile(fn, out_fn = None):
    pass

    return out_filename


def _cut_sectionstring(edistring,sectionhead):

    #in this case, several blocks have to be handled together, therefore, a simple cut to the next block start does not work:
    if sectionhead.upper() == 'HMEAS_EMEAS':
        #required for finding HMEAS and EMEAS at once:
        import re

        lo_start_idxs = [m.start() for m in re.finditer('>[HE]MEAS', edistring) ]
        if len(lo_start_idxs) == 0 :
            del re
            raise

        start_idx = lo_start_idxs[0]

        end_idx = edistring[(lo_start_idxs[-1]+1):].upper().find('>') + lo_start_idxs[-1]

        hmeas_emeas_string = edistring[start_idx:end_idx]

        if len(hmeas_emeas_string) == 0:
            del re
            raise
        
        del re
        return hmeas_emeas_string
 


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



