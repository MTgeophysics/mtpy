#!/usr/bin/env python

"""
mtpy/mtpy/core/edi.py

Contains classes and functions for handling EDI files.

    Class:
    "Edi" contains all information from or for an EDI file. Sections of EDI files are given as respective attributes, section-keys and values are stored in dictionaries.

        Methods:

        - readfile
        - edi_dict
        - data_dict
        - z_dict
        - tipper_dict
        - periods
        - frequencies
        - n_freqs
        - _read_head
        - _read_info
        - _read_definemeas
        - _read_hmeas_emeas
        - _read_mtsect
        - _read_freq
        - _read_z
        - _read_tipper
        - _read_zrot
        - writefile
        - z2resphase
        - rotate
        - rho
        - phi
        - set_rho_phi
        - set_head
        - set_info_dict
        - set_info_string
        - set_z
        - set_zerr
        - set_tipper
        - set_tippererr
        - set_definemeas
        - set_mtsect
        - get_datacomponent
        - set_frequencies
        - set_zrot


    Functions:

    - read_edifile
    - write_edifile
    - combine_edifiles
    - validate_edifile
    - rotate_edifile
    - _generate_edifile_string
    - _cut_sectionstring
    - _validate_edifile_string

@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np
import os,sys
import os.path as op
import math, cmath
import time, calendar
import copy
#required for finding HMEAS and EMEAS at once:
import re

import mtpy.utils.format as MTformat
import mtpy.utils.calculator as MTc
import mtpy.utils.exceptions as MTexceptions

#reload(MTexceptions)
reload(MTformat)
reload(MTc)


#=================================================================

class Edi(object):
    """
        Edi class - generates an edi-object.

        Methods  include reading and writing from and to edi-files, rotations/combinations of edi-files, as well as 'get' and 'set' for all edi file sections

        Errors are given as standard deviations (sqrt(VAR))


    """

    def __init__(self, fn = None):

        """
            Initialise an instance of the Edi class.

            Optional input:
            fn : filename to connect with the object - not read in automatically!

            Initialise the attributes with None/empty dictionary
        """

        self.filename = fn
        if fn != None:
            if op.isfile(op.abspath(fn)):
                self.filename = op.abspath(fn)
            else:
                self.filename = None

        self.in_filestring = None
        self._head = {}
        self._info_string = None
        self._info_dict = {}
        self._definemeas = {}
        self._hmeas_emeas = None
        self._mtsect = {}
        self._freq = None
        self.zrot = None
        self.z = None
        self.zerr = None
        self.tipper = None
        self.tippererr = None


    def readfile(self, fn, Bscale = 1, datatype = 'z'):
        """
            Read in an EDI file.

            Returns an exception, if the file is invalid (following MTpy standards).


            Bscale is the factor that needs to be applied to the data in the file to end up with km/s as unit for Z.
            E.g., if the Z-data are given in E/H with E measured in V/m and H in A/m, the factor is 1e-3 / mu_0
            (If the data are in basic units "mV/m" and "A/m", just use the letter 'H' instead )

            datatype determines the way data  are provided. Default is 'z', so the full impedance tensor is expected to be present in the file. Other possibilities are 'rhophi' and 'spectra' - they exclude the reading of a potentially present Z information.
            TODO: 'spectra' - not implemented yet


        """
        infile = op.abspath(fn)

        #define the scaling factor for obtaining Z in Ohm
        if Bscale in ['h','H']:
            Bscale =  1./MTc.mu0
        try:
            Bscale = float(Bscale)
        except:
            raise MTexceptions.MTpyError_edi_file('ERROR - B field scaling factor not understood ')

        try:
            datatype = datatype.lower()
            if not datatype in ['z' , 'rhophi', 'spectra']:
                raise
        except:
            raise MTexceptions.MTpyError_edi_file('ERROR - datatype not understood')

        #check for existence
        if not op.isfile(infile):
            raise MTexceptions.MTpyError_edi_file('File is not existing: %s'%infile)

        with open(infile,'r') as F:
            edistring = F.read()

        #validate edi file string following MTpy standard
        if not _validate_edifile_string(edistring):
            raise MTexceptions.MTpyError_edi_file('%s is no proper EDI file'%infile)

        self.filename = infile
        self.in_filestring = edistring

        #read out the mandatory EDI file sections from the raw string
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
            print 'Could not read MTSECT section: %s'%infile

        try:
            self._read_freq(edistring)
        except:
            print 'Could not read FREQ section: %s'%infile

        if datatype == 'z':
            try:
                self._read_z(edistring, Bscale)
            except:
                raise MTexceptions.MTpyError_edi_file('Could not read Z section: %s'%infile)

        elif datatype == 'rhophi':
            try:
                self._read_rhophi(edistring)
            except:
                raise MTexceptions.MTpyError_edi_file('Could not read RhoPhi section: %s'%infile)
            #rotation is optional
            try:
                self._read_rhorot(edistring)
            except:
                self.zrot = list(np.zeros((self.n_freqs())))
                print 'Could not read Rhorot section: %s'%infile


        elif datatype == 'spectra':
            if 1:
                self._read_spectra(edistring)
            # except:
            #     raise MTexceptions.MTpyError_edi_file('Could not read Spectra section: %s'%infile)


        #Tipper is optional
        try:
            self._read_tipper(edistring)
        except:
            self.tipper = None
            self.tippererr = None
            print 'Could not read Tipper section: %s'%infile

        #rotation is optional
        try:
            self._read_zrot(edistring)
        except:
            self.zrot = None
            print 'Could not read Zrot section: %s'%infile


    def edi_dict(self):
        """
            Collect sections of the EDI file and return them as a dictionary.

        """

        edi_dict = {}

        edi_dict['HEAD'] = self.head
        edi_dict['INFO'] = self.info_dict
        edi_dict['DEFINEMEAS'] = self.definemeas
        edi_dict['HMEAS_EMEAS'] = self.hmeas_emeas
        edi_dict['MTSECT'] = self.mtsect
        edi_dict['FREQ'] = self.freq

        #update the dictionary information from the z and tipper arrays (those may have changed by rotation):
        edi_dict['Z'] = self.z_dict()
        edi_dict['TIPPER'] = self.tipper_dict()
        edi_dict['ZROT'] = self.zrot


        return  edi_dict


    def data_dict(self):
        """
            Return collected raw data information in one dictionary:
            Z, Tipper, Zrot, frequencies, Zerror, TipperError

        """
        data_dict = {}

        data_dict['z'] = self.z
        data_dict['tipper'] = self.tipper
        data_dict['zrot'] = self.zrot
        data_dict['frequencies'] = self.freq
        data_dict['zerr'] = self.zerr
        data_dict['tippererr'] = self.tippererr

        return data_dict

    def z_dict(self):
        """
            Return the content of the Z and Zerror arrays in a dictionary.
        """

        new_z_dict = {}
        z_array = self.z
        zerr_array = self.zerr
        compstrings = ['ZXX','ZXY','ZYX','ZYY']
        Z_entries = ['R','I','.VAR']

        for idx_comp,comp in enumerate(compstrings):
            for idx_zentry,zentry in enumerate(Z_entries):
                section = comp + zentry
                if idx_zentry == 0:
                    new_z_dict[section] = list(np.real(z_array[:,idx_comp/2, idx_comp%2]))
                elif idx_zentry == 1:
                    new_z_dict[section] = list(np.imag(z_array[:,idx_comp/2, idx_comp%2]))
                elif idx_zentry == 2:
                    #squaring the errors (stddev) to get VAR values
                    new_z_dict[section] = list( (zerr_array[:,idx_comp/2, idx_comp%2])**2 )

        return new_z_dict

    def tipper_dict(self):
        """
            Return the content of the Tipper and TipperError arrays in a dictionary.
        """

        new_t_dict = {}
        t_array = self.tipper
        if  t_array is None:
            return None
        terr_array = self.tippererr
        compstrings = ['TX','TY']
        T_entries = ['R','I','VAR']

        for idx_comp,comp in enumerate(compstrings):
            for idx_tentry,tentry in enumerate(T_entries):
                section = comp + tentry
                if idx_tentry == 0:
                    new_t_dict[section] = list(np.real(t_array[:,idx_comp/2, idx_comp%2]))
                elif idx_tentry == 1:
                    new_t_dict[section] = list(np.imag(t_array[:,idx_comp/2, idx_comp%2]))
                elif idx_tentry == 2:
                    #square errors (stddev) to get VAR values
                    new_t_dict[section] = list( (terr_array[:,idx_comp/2, idx_comp%2])**2)


        return new_t_dict

    def _get_periods(self):
        """
            Return a list of periods (output values in seconds).
        """

        return list( 1./np.array(self.freq) )
    
    def _set_periods(self, list_of_periods):
        """
            Set frequencies by a list of periods (values in seconds).
        """
        if len(list_of_periods) is not len(self.z):
            print 'length of periods list not correct (%i instead of %i)'%(len(list_of_periods), len(self.z))
            return
        self.freq = list(1./np.array(list_of_periods) )

    periods = property(_get_periods, _set_periods, doc='List of periods (values in seconds)')    

    # def frequencies(self):
    #     """
    #         Return a list of the frequencies (output values in Hertz).
    #     """

    #     return list(self.freq)

    def n_freqs(self):
        """
            Return the number of frequencies/length of the Z array.
        """

        return len(self.freq)


    def _getlat(self): return self.head['lat']
    def _setlat(self, value): self.head['lat'] = MTformat._assert_position_format('lat',value)
    lat = property(_getlat, _setlat, doc='Location latitude in degrees')
   

    def _getlong(self): return self.head['long']
    def _setlong(self, value): self.head['long'] = MTformat._assert_position_format('long',value)
    long = property(_getlong, _setlong, doc='Location longitude in degrees')
 
    def _getelev(self): return self.head['elev']
    def _setelev(self, value): self.head['elev'] = MTformat._assert_position_format('elev',value)
    elev = property(_getelev, _setelev, doc='Location elevation in meters')
 


    def _read_head(self, edistring):
        """
            Read in the HEAD  section from the raw edi-string.
        """

        try:
            temp_string = _cut_sectionstring(edistring,'HEAD')
        except:
            raise

        head_dict = {}
        t1 = temp_string.strip().split('\n')
        t2 = [i.strip() for i in t1 if '=' in i]
        for j in t2:
            k = j.split('=')
            key = str(k[0]).lower().strip()
            value = k[1].replace('"','')
            if key in ['lat','long','lon','latitude','longitude','ele','elev','elevation']:
                value = MTformat._assert_position_format(key,value)

            head_dict[key] = value

        self._head = head_dict

    def _read_info(self, edistring):
        """
            Read in the INFO  section from the raw edi-string.
        """

        try:
            temp_string = _cut_sectionstring(edistring,'INFO')
        except:
            raise

        self._info_string = temp_string.strip()
        

        t1 = temp_string.strip().split('\n')
        t2 = [i.strip() for i in t1 if '=' in i or ':' in i]

        info_dict = {}

        for tmp_str in t2:
            #fill dictionary
            #ignore lines with no information after '='' or ':'

            if '=' in tmp_str:
                t3 = tmp_str.split('=')
                key = str(t3[0]).lower().strip()
                value = t3[1].replace('"','').strip()
                if not len(value) == 0:
                    info_dict[key] = value

            elif ':' in tmp_str:
                #consider potential ':' characters in coordinates!
                t3 = tmp_str.split(':')
                key = str(t3[0]).lower().strip()
                value = t3[1:]
                value = [i.strip().replace('"','').strip() for i in value]
                if len(value) > 1:
                    value = ':'.join(value)
                else:
                    value = value[0]

                if not len(value) == 0:
                    info_dict[key] = value

        if 'station' not in info_dict:
            station = self.head['dataid']
            if len(station) == 0:
                station = op.splitext(op.basename(self.filename))[0]
            info_dict['station'] = station

        info_dict['Z_unit'] = 'km/s'

        self._info_dict = info_dict



    def _read_definemeas(self, edistring):
        """
            Read in the DEFINEMEAS  section from the raw edi-string.
        """

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
                value = k[1].replace('"','').strip()
                if len(value) != 0:
                    if 'lat' in key:
                        value = MTformat._assert_position_format('lat',value)
                    if 'lon' in key:
                        value = MTformat._assert_position_format('lon',value)


                    d_dict[key] = value

        if len(d_dict.keys()) == 0:
            raise


        self._definemeas = d_dict



    def _read_hmeas_emeas(self, edistring):
        """
            Read in the HMEAS/EMEAS  section from the raw edi-string.
        """
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

        self._hmeas_emeas = lo_hmeas_emeas

        

    def _read_mtsect(self, edistring):
        """
            Read in the MTSECT  section from the raw edi-string.
        """

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
                value = k[1].replace('"','').strip()
                if len(value) != 0:
                    m_dict[key] = value

        if len(m_dict.keys()) == 0:
            raise


        self._mtsect = m_dict


    def _read_freq(self, edistring):
        """
            Read in the FREQ  section from the raw edi-string.
        """

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
                    passs

        self._freq = lo_freqs


    def _read_z(self, edistring, Bscale):
        """
        Read in impedances information from a raw EDI-string.
        Store it as attribute (complex array).

        """

        compstrings = ['ZXX','ZXY','ZYX','ZYY']
        Z_entries = ['R','I','.VAR']

        z_array = np.zeros((self.n_freqs(),2,2),dtype=np.complex)
        zerr_array = np.zeros((self.n_freqs(),2,2),dtype=np.float)
        z_dict = {}

        for idx_comp,comp in enumerate(compstrings):
            for idx_zentry,zentry in enumerate(Z_entries):
                sectionhead = comp + zentry
                try:
                    temp_string = _cut_sectionstring(edistring,sectionhead)
                except:
                    continue

                lo_z_vals = []

                #check, if correct number of entries are given in the block
                t0 = temp_string.strip().split('\n')[0]
                n_dummy = int(float(t0.split('//')[1].strip()))
                if not n_dummy == self.n_freqs():
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


        for idx_freq  in range( self.n_freqs()):
            try:
                z_array[idx_freq,0,0] = np.complex(z_dict['ZXXR'][idx_freq], z_dict['ZXXI'][idx_freq])
            except:
                pass
            try:
                z_array[idx_freq,0,1] = np.complex(z_dict['ZXYR'][idx_freq], z_dict['ZXYI'][idx_freq])
            except:
                pass
            try:
                z_array[idx_freq,1,0] = np.complex(z_dict['ZYXR'][idx_freq], z_dict['ZYXI'][idx_freq])
            except:
                pass
            try:
                z_array[idx_freq,1,1] = np.complex(z_dict['ZYYR'][idx_freq], z_dict['ZYYI'][idx_freq])
            except:
                pass


            for idx_comp,comp in enumerate(compstrings):
                sectionhead = comp + '.VAR'
                if sectionhead in z_dict:
                    zerr_array[idx_freq, idx_comp/2, idx_comp%2] = z_dict[sectionhead][idx_freq]


        self.z = Bscale * z_array

        #errors are stddev, not VAR :
        self.zerr = np.sqrt(Bscale * zerr_array)


    def _read_tipper(self, edistring):
        """
        Read in Tipper information from a raw EDI-string.
        Store it as attribute (complex array).

        """

        compstrings = ['TX','TY']
        T_entries = ['R','I','VAR']

        tipper_array = np.zeros((self.n_freqs(),1,2),dtype=np.complex)
        tippererr_array = np.zeros((self.n_freqs(),1,2),dtype=np.float)
        t_dict = {}


        for idx_comp,comp in enumerate(compstrings):
            for idx_tentry,tentry in enumerate(T_entries):
                temp_string = None
                try:
                    sectionhead = comp + tentry + '.EXP'
                    temp_string = _cut_sectionstring(edistring,sectionhead)
                except:
                    try:
                        sectionhead = comp + tentry
                        temp_string = _cut_sectionstring(edistring,sectionhead)
                    except:
                        # if tipper is given with sectionhead "TX.VAR"
                        if (idx_tentry == 2) and (temp_string is None):
                            try:
                                sectionhead = comp + '.' + tentry
                                temp_string = _cut_sectionstring(edistring,sectionhead)
                            except:
                                pass
                        pass

                lo_t_vals = []

                #check, if correct number of entries are given in the block
                t0 = temp_string.strip().split('\n')[0]
                n_dummy = int(float(t0.split('//')[1].strip()))

                if not n_dummy == self.n_freqs():
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


        for idx_freq  in range( self.n_freqs()):
            tipper_array[idx_freq,0,0] = np.complex(t_dict['TXR'][idx_freq], t_dict['TXI'][idx_freq])
            tippererr_array[idx_freq,0,0] = t_dict['TXVAR'][idx_freq]
            tipper_array[idx_freq,0,1] = np.complex(t_dict['TYR'][idx_freq], t_dict['TYI'][idx_freq])
            tippererr_array[idx_freq,0,1] = t_dict['TYVAR'][idx_freq]


        self.tipper = tipper_array
        #errors are stddev, not VAR :
        self.tippererr = np.sqrt(tippererr_array)

    def _read_rhophi(self, edistring):
        """
            Read in RhoPhi information from a raw EDI-string.
            Convert the information into Z and Zerr.
            Store this as attribute (complex array).

        """
        #using the loop over all  components. For each component check, if Rho and Phi are given, raise exception if not! Then convert the polar RhoPhi representation into the cartesian Z. Rho is assumed to be in Ohm m, Phi in degrees. Z will be in km/s.
        z_array = np.zeros((self.n_freqs(),2,2), 'complex')
        zerr_array = np.zeros((self.n_freqs(),2,2))

        rhophistrings = ['RHO','PHS']
        compstrings = ['XX','XY','YX','YY']
        entries = ['','.ERR']

        rhophi_dict = {}
        for rp in rhophistrings:
            for comp in compstrings:
                for entry in entries:
                    sectionhead = rp + comp + entry
                    try:
                        temp_string = _cut_sectionstring(edistring,sectionhead)
                    except:
                        continue

                    lo_vals = []
                    #check, if correct number of entries are given in the block
                    t0 = temp_string.strip().split('\n')[0]
                    n_dummy = int(float(t0.split('//')[1].strip()))
                    if not n_dummy == self.n_freqs():
                        raise

                    t1 = temp_string.strip().split('\n')[1:]
                    for j in t1:
                        lo_j = j.strip().split()
                        for k in lo_j:
                            try:
                                lo_vals.append(float(k))
                            except:
                                pass

                    rhophi_dict[sectionhead] = lo_vals
        #print rhophi_dict
        #sys.exit()

        for idx_freq  in range( self.n_freqs()):
            r = np.zeros((2,2))
            phi = np.zeros((2,2))
            rerr = np.zeros((2,2))
            phierr = np.zeros((2,2))
            zerr = np.zeros((2,2))

            for idx_c, comp in enumerate(compstrings):
                #convert rho to amplitude:
                try:
                    r[idx_c/2,idx_c%2] = np.sqrt(rhophi_dict['RHO'+comp][idx_freq] * 5 * self.freq[idx_freq] )
                except:
                    pass
                try:
                    phi[idx_c/2,idx_c%2] = rhophi_dict['PHS'+comp][idx_freq]
                except:
                    pass
                try:
                    try:
                        #check for small amplitude close to zero
                        if r[idx_c/2,idx_c%2] == 0:
                            raise
                        f1 = np.abs(np.sqrt(2.5*self.freq[idx_freq]/r[idx_c/2,idx_c%2])*(rhophi_dict['RHO'+comp + '.ERR'][idx_freq]))
                    except:
                        f1 = np.sqrt(rhophi_dict['RHO'+comp+ '.ERR'][idx_freq] * 5 * self.freq[idx_freq] )
                    rerr[idx_c/2,idx_c%2] = f1
                except:
                    pass
                try:
                    phierr[idx_c/2,idx_c%2] = rhophi_dict['PHS'+comp + '.ERR'][idx_freq]
                except:
                    pass
                zerr[idx_c/2,idx_c%2] = max( MTc.propagate_error_polar2rect( r[idx_c/2,idx_c%2], rerr[idx_c/2,idx_c%2], \
                                                                            phi[idx_c/2,idx_c%2], phierr[idx_c/2,idx_c%2]))

            z_array[idx_freq] = MTc.rhophi2z(r, phi)
            zerr_array[idx_freq] = zerr


        self.z = z_array
        self.zerr = zerr_array


    def _read_rhorot(self, edistring):
        """
            Read in the (optional) RhoRot  section from the raw edi-string for data file containing data in  RhoPhi style. Angles are stored in the ZROT attribute. 
        """

        try:
            temp_string = _cut_sectionstring(edistring,'RHOROT')
        except:
            lo_angles = list( np.zeros((self.n_freqs())) )
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


        if len(lo_angles) != self.n_freqs():
            raise


    def _read_spectra(self,edistring):
        """
            Read in Spectra information from a raw EDI-string.
            Convert the information into Z and Zerr.
            Store this as attribute (complex array).

        """

        #identify and cut spectrasect part:
        specset_string = _cut_sectionstring(edistring,'SPECTRASECT')
        dummy4 = specset_string.upper().find('NCHAN')
        n_chan = int(float(specset_string[dummy4:].strip().split('=')[1].split()[0]))
        id_list = specset_string.split('//')[1].split('\n')[1].strip().split()

        dummy5 = specset_string.upper().find('NFREQ')
        n_freq = int(float(specset_string[dummy5:].strip().split('=')[1].split()[0]))
        
        lo_spectra_strings = []
        tmp_string = copy.copy(edistring)

        #read in all SPECTRA subsections into a list 
        while True:
            try:
                dummy3 = tmp_string.find('>SPECTRA')
                #check, if SPECTRA subsection exists
                if dummy3 <0 :
                    raise               
                # cut the respective sub string
                tmp_cut_string = _cut_sectionstring(tmp_string,'SPECTRA')
                #append to the list
                lo_spectra_strings.append(tmp_cut_string)
                # reduce the input string by the keyword 'SPECTRA', so the subsequent subsection will be read in the next loop-run instead
                tmp_string = tmp_string[:dummy3] + tmp_string[dummy3+8:]
            #exceptions for breaking the while, called, if no more SPECTRA subsections can be found
            except:
                break

        #assert that the list of read in SPECTRA subsection is not empty:
        if len(lo_spectra_strings) == 0:
            raise MTexceptions.MTpyError_EDI('ERROR - EDI file does not contain readable SPECTRA sections!')

        z_array = np.zeros((len(lo_spectra_strings),2,2), 'complex')
        zerr_array = np.zeros((len(lo_spectra_strings),2,2))  

        lo_freqs = []

        id_channel_dict = _build_id_channel_dict(self.hmeas_emeas)
        channellist = [id_channel_dict[i] for i in id_list]
        for j in ['HX', 'HY', 'EX', 'EY'] :
            if j not in channellist:
                raise MTexceptions.MTpyError_edi_file('Mandatory channel data for %s not given !'%j)

        for s_idx, spectra in enumerate(lo_spectra_strings):
            firstline = spectra.split('\n')[0]
            freq = float(_find_key_value('FREQ','=',firstline))
            lo_freqs.append(freq)
            datalist = []
            for innerline in spectra.split('\n')[1:]:
                datalist.extend(innerline.strip().split())
            data = np.array([float(i) for i in datalist]).reshape(7,7)


        sys.exit()




        pass
        # z_array = np.zeros((self.n_freqs(),2,2), 'complex')
        # zerr_array = np.zeros((self.n_freqs(),2,2))
 
        # self.z = z_array
        # self.zerr = zerr_array


    def _read_zrot(self, edistring):
        """
            Read in the (optional) Zrot  section from the raw edi-string.
        """

        try:
            temp_string = _cut_sectionstring(edistring,'ZROT')
        except:
            lo_angles = list( np.zeros((self.n_freqs())) )
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


        if len(lo_angles) != self.n_freqs():
            raise

        self.zrot = lo_angles


    def writefile(self, *fn):
        """
            Write out the edi object into an EDI file.
        """

        if len(fn) == 0 :
            fn = None
        else:
            fn = fn[0]

        self.info_dict['edifile_generated_with'] = 'MTpy'

        outstring, stationname = _generate_edifile_string(self.edi_dict())

        if not _validate_edifile_string(outstring):
            #return outstring
            raise MTexceptions.MTpyError_edi_file('Cannot write EDI file...output string is invalid')


        if fn != None:
            try:
                outfilename = op.abspath(fn)
                if not outfilename.lower().endswith('.edi'):
                    outfilename += '.edi'
            except:
                fn = None

        if fn == None:
            outfilename = op.abspath(stationname.upper()+'.edi')

        if op.isfile(outfilename):
            newfile = outfilename

            i = 0
            while op.isfile(newfile):
                i += 1
                newfile = outfilename+'_%i'%i

            outfilename = newfile

        try:
            with open(outfilename , 'w') as F:
                F.write(outstring)
        except:
            raise MTexceptions.MTpyError_edi_file('Cannot write EDI file: %s'%(outfilename))

        return outfilename




    def rotate(self,angle):
        """
            Rotate the Z and tipper information in the Edi object. Change the rotation angles in Zrot respectively.

            Rotation angle must be given in degrees. All angles are referenced to geographic North, positive in clockwise direction. (Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.

            Updates the attributes "z, zerr, zrot, tipper, tippererr".

        """

        angle = angle%360
        zerr_rot = None
        tipper_rot = None
        tippererr_rot = None

        z_rot = copy.copy(self.z)
        if self.zerr is not None:
            zerr_rot = copy.copy(self.zerr)
        if self.tipper is not None:
            tipper_rot = copy.copy(self.tipper)
        if self.tippererr is not None:
            tippererr_rot = copy.copy(self.tippererr)

        for idx_freq in range(self.n_freqs()):

            if self.zerr is not None:
                z_rot[idx_freq], zerr_rot[idx_freq] = MTc.rotatematrix_incl_errors(self.z[idx_freq,:,:], angle, self.zerr[idx_freq,:,:])
            else:
                z_rot[idx_freq], zerr_rot = MTc.rotatematrix_incl_errors(self.z[idx_freq,:,:], angle)


            if self.tipper is not None:

                if self.tippererr is not None:
                    tipper_rot[idx_freq], tippererr_rot[idx_freq] = MTc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], angle,self.tippererr[idx_freq,:,:] )
                else:
                    tipper_rot[idx_freq], tippererr_rot = MTc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], angle)



        self.z = z_rot
        if zerr_rot is not None:
            self.zerr = zerr_rot
        if tipper_rot is not None:
            self.tipper = tipper_rot
        if tippererr_rot is not None:
            self.tippererr = tippererr_rot

        self.zrot = list( (np.array(self.zrot) + angle)%360)


    def _get_rho_phi(self):
        """
            Return values for resistivity (rho - in Ohm m) and phase (phi - in degrees).

            Output is a 4-tuple of arrays:
            (Rho, Phi, RhoError, PhiError)
        """

        if self.z is None:
            print 'Z array is None - cannot calculate rho/phi'
            return
        rhoerr = None
        phierr = None
        if self.zerr is not None:
            rhoerr = np.zeros(self.zerr.shape)
            phierr = np.zeros(self.zerr.shape)

        rho = np.zeros(self.z.shape)
        phi = np.zeros(self.z.shape)


        for idx_f in range(len(self.z)):
            for i in range(2):
                for j in range(2):

                    rho[idx_f,i,j] = np.abs(self.z[idx_f,i,j])**2 /self.freq[idx_f] *0.2
                    phi[idx_f,i,j] = math.degrees(cmath.phase(self.z[idx_f,i,j]))

                    if self.zerr is not None:
                        r_err, phi_err = MTc.propagate_error_rect2polar( np.real(self.z[idx_f,i,j]), self.zerr[idx_f,i,j], np.imag(self.z[idx_f,i,j]), self.zerr[idx_f,i,j])
                        rhoerr[idx_f,i,j] = 0.4 * np.abs(self.z[idx_f,i,j])/self.freq[idx_f] * r_err
                        phierr[idx_f,i,j] = phi_err

        return rho, phi, rhoerr, phierr



    def _set_rho_phi(self, rho_array, phi_array):
        """
            Set values for resistivity (rho - in Ohm m) and phase (phi - in degrees).

            Updates the attributes "z, zerr".

        """

        if self.z is not None:
            z_new = copy.copy(self.z)

            if self.z.shape != rho_array.shape:
                print 'Error - shape of "rho" array does not match shape of Z array: %s ; %s'%(str(rho_array.shape),str(self.z.shape))
                return

            if self.z.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of Z array: %s ; %s'%(str(phi_array.shape),str(self.z.shape))
                return
        else:
            z_new = p.zeros(rho_array.shape,'complex')
            if rho_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of "rho" array: %s ; %s'%(str(phi_array.shape),str(rho_array.shape))
                return

        if (self.freq is None) or (len(self.freq) != len(rho_array)) :
            raise MTexceptions.MTpyError_EDI('ERROR - cannot set rho without proper frequency information - proper "freq" attribute must be defined ')

        #assert real array:
        if np.linalg.norm(np.imag(rho_array )) != 0 :
            print 'Error - array "rho" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(z_new)):
            for i in range(2):
                for j in range(2):
                    abs_z = np.sqrt(5 * self.freq[idx_f] * rho_array[idx_f,i,j])
                    z_new[idx_f,i,j] = cmath.rect( abs_z, math.radians(phi_array[idx_f,i,j] ))

        self.z = z_new

    rho_phi = property(_get_rho_phi,_set_rho_phi,doc='Values for resistivity (rho - in Ohm m) and phase (phi - in degrees). Updates the attributes "z, zerr"')



    def _set_head(self, head_dict):
        """
            Set the attribute 'head'.

            Input:
            HEAD section dictionary

            No test for consistency!

        """

        self._head = head_dict
    def _get_head(self): return self._head
    head = property(_get_head, _set_head, doc='HEAD attribute of EDI file')



    def _set_info_dict(self,info_dict):
        """
            Set the attribute 'info_dict'.

            Input:
            INFO section dictionary

            No test for consistency!

        """

        self._info_dict = info_dict
    def _get_info_dict(self): return self._info_dict
    info_dict = property(_get_info_dict, _set_info_dict, doc='INFO section dictionary')


    def _set_info_string(self,info_string):
        """
            Set the attribute 'info_string'.

            Input:
            INFO section string

            No test for consistency!

        """

        self._info_string = info_string
    def _get_info_string(self): return self._info_string
    info_string = property(_get_info_string, _set_info_string, doc='INFO section string')


    def set_z(self, z_array):
        """
            Set the attribute 'z'.

            Input:
            Z array

            Test for shape, but no test for consistency!

        """


        if (self.z is not None) and (self.z.shape != z_array.shape):
            print 'Error - shape of "z" array does not match shape of existing Z array: %s ; %s'%(str(z_array.shape),str(self.z.shape))
            return

        self.z = z_array


    def set_zerr(self, zerr_array):
        """
            Set the attribute 'zerr'.

            Input:
            Zerror array

            Test for shape, but no test for consistency!

        """

        if (self.zerr is not None) and (self.zerr.shape != zerr_array.shape):
            print 'Error - shape of "zerr" array does not match shape of Zerr array: %s ; %s'%(str(zerr_array.shape),str(self.zerr.shape))
            return

        self.zerr = zerr_array


    def set_tipper(self, tipper_array):
        """
            Set the attribute 'tipper'.

            Input:
            tipper array

            Test for shape, but no test for consistency!

        """

        if (self.tipper is not None) and (self.tipper.shape != tipper_array.shape):
            print 'Error - shape of "tipper" array does not match shape of tipper-array: %s ; %s'%(str(tipper_array.shape),str(self.tipper.shape))
            return

        self.tipper = tipper_array


    def set_tippererr(self, tippererr_array):
        """
            Set the attribute 'tippererr'.

            Input:
            TipperError array

            Test for shape, but no test for consistency!

        """


        if (self.tippererr is not None) and (self.tippererr.shape != tippererr_array.shape):
            print 'Error - shape of "tippererr" array does not match shape of tippererr array: %s ; %s'%(str(tippererr_array.shape),str(self.tippererr.shape))
            return

        self.tippererr = tippererr_array


    def _set_definemeas(self,definemeas_dict):
        """
            Set the attribute 'definemeas'.

            Input:
            DEFINEMEAS section dictionary

            No test for consistency!

        """
        self._definemeas = definemeas_dict
    def _get_definemeas(self): return self._definemeas
    definemeas = property(_get_definemeas, _set_definemeas, doc='DEFINEMEAS section dictionary')

    def _set_hmeas_emeas(self,hmeas_emeas_list):
        """
            Set the attribute 'hmeas_emeas'.

            Input:
            hmeas_emeas section list of 7-tuples

            No test for consistency!

        """
        self._hmeas_emeas = hmeas_emeas_list
    def _get_hmeas_emeas(self): return self._hmeas_emeas
    hmeas_emeas = property(_get_hmeas_emeas, _set_hmeas_emeas, doc='hmeas_emeas section list of 7-tuples')


    def _set_mtsect(self, mtsect_dict):
        """
            Set the attribute 'mtsect'.

            Input:
            MTSECT section dictionary

            No test for consistency!

        """

        self._mtsect = mtsect_dict
    def _get_mtsect(self): return self._mtsect
    mtsect = property(_get_mtsect, _set_mtsect, doc='MTSECT section dictionary')



    def get_datacomponent(self, componentname):
        """
            Return a specific data component.

            Input:
            specification of the data component (Z or Tipper components)
        """

        data_dict = self.data_dict()
        if componentname.lower() in data_dict:
            return data_dict[componentname.lower()]

        compstrings = ['ZXX','ZXY','ZYX','ZYY']
        Z_entries = ['R','I','.VAR']
        for idx_comp,comp in enumerate(compstrings):
            for idx_zentry,zentry in enumerate(Z_entries):
                section = comp + zentry
                if section.lower() == componentname.lower():
                    return self.z_dict()[section]

        compstrings = ['TX','TY']
        T_entries = ['R','I','VAR']

        for idx_comp,comp in enumerate(compstrings):
            for idx_tentry,tentry in enumerate(T_entries):
                section = comp + tentry
                if section.lower() == componentname.lower():
                    return self.tipper_dict()[section]

        print 'unknown data component: %s'%componentname.lower()
        return


    def set_datacomponent(self, componentname, value):
        """
            Set a specific data component.

            Input:
            specification of the data component (Z or Tipper components)
            new value

            No test for consistency!
        """
        pass


    def _set_frequencies(self, lo_frequencies):
        """
            Set the list of frequencies.

            Input:
            list of frequencies

            No test for consistency!
        """

        if len(lo_frequencies) is not len(self.z):
            print 'length of frequency list not correct (%i instead of %i)'%(len(lo_frequencies), len(self.z))
            return

        self._freq = lo_frequencies
    def _get_frequencies(self): return self._freq
    freq = property(_get_frequencies, _set_frequencies, doc='list of frequencies')


    def set_zrot(self, angle):
        """
            Set the list of rotation angles.

            Input:
            single angle or list of angles (in degrees)

            No test for consistency!
        """
        if np.iterable(angle):
            if len(angle) is not len(self.z):
                print 'length of angle list not correct (%i instead of %i)'%(len(angle), len(self.z))
                return
            try:
                angle = [float(i%360) for i in angle]
            except:
                raise MTexceptions.MTpyError_edi_file('list of angles contains non-numercal values')
        else:
            try:
                angle = [float(angle%360) for i in self.z]
            except:
                raise MTexceptions.MTpyError_edi_file('Angles is a non-numercal value')


        self.zrot = angle


#end of Edi Class
#=========================


def read_edifile(fn, Bscale = 1):
    """
        Read in an EDI file.

        Return an instance of the Edi class.
    """

    if Bscale in ['h','H']:
        Bscale =  1./MTc.mu0
    try:
        Bscale = float(Bscale)
    except:
        raise MTexceptions.MTpyError_edi_file('ERROR - B field scaling factor not understood ')


    edi_object = Edi()

    edi_object.readfile(fn, Bscale)


    return edi_object


def write_edifile(edi_object, out_fn = None):
    """
        Write an EDI file from an instance of the Edi class.

        optional input:
        EDI file name
    """

    if not isinstance(z_object, MTedi.Edi):
        raise MTexceptions.MTpyError_EDI('Input argument is not an instance of the Edi class')

    if out_fn is not None:
        dirname = op.dirname(op.abspath(op.join('.',out_fn)))
        fn = op.basename(op.abspath(op.join('.',out_fn)))
        if not op.isdir(dirname):
            try:
                os.makedirs(dirname)
                out_fn = op.join(dirname,fn)
            except:
                out_fn = None
        else:
            out_fn = op.join(dirname,fn)

    outfilename = None
    try:
        outfilename = edi_object.writefile(out_fn)
    except:
        print 'Cannot write EDI file...output string invalid!'

    return outfilename


def combine_edifiles(fn1, fn2,  merge_frequency=None, out_fn = None, allow_gaps = True, Bscale1 = 1, Bscale2 = 1 ):
    """
        Combine two EDI files.

        Inputs:
        - name of EDI file 1
        - name of EDI file 2

        optional input:
        - merge_frequency : frequency in Hz, on which to merge the files - default is the middle of the overlap
        - out_fn : output EDI file name
        - allow_gaps : allow merging EDI files whose frequency ranges does not overlap
        - Bscale1/2 : scaling factor to bring the Z values to unit km/s


        Outputs:
        - instance of Edi class, containing merged information
        - full path of the output EDI file
    """



    if Bscale1 in ['h','H']:
        Bscale1 =  1./MTc.mu0
    try:
        Bscale1 = float(Bscale1)
    except:
        raise MTexceptions.MTpyError_edi_file('ERROR - B field scaling factor 1 not understood ')

    if Bscale2 in ['h','H']:
        Bscale2 =  1./MTc.mu0
    try:
        Bscale2 = float(Bscale2)
    except:
        raise MTexceptions.MTpyError_edi_file('ERROR - B field scaling factor 2 not understood ')


    #edi objects:
    eo1 = Edi()
    eo1.readfile(fn1, Bscale1)
    eo2 = Edi()
    eo2.readfile(fn2, Bscale2)
    #edi object merged
    eom = Edi()

    #check frequency lists
    lo_freqs1 = eo1.frequencies()
    lo_freqs2 = eo2.frequencies()


    lo_eos = []

    #check for overlap of the frequency regimes:
    if (not min(lo_freqs1) > max(lo_freqs2)) and (not max(lo_freqs1) > min(lo_freqs2)):
        if allow_gaps is False:
            raise MTexceptions.MTpyError_edi_file('Cannot merge files %s and %s - frequency ranges do not overlap and "allow_gaps" is set to False')


    #determine, which is the low frequency part
    lo_eos = [eo1, eo2]

    if min(lo_freqs1) <= min(lo_freqs2):
        if max(lo_freqs1) >= max(lo_freqs2):
            print 'Frequency range of file %s fully contained in range of file %s => no merging of files!'%(fn2, fn1)
            return

    if min(lo_freqs1) >= min(lo_freqs2):
        if max(lo_freqs1) <= max(lo_freqs2):
            print 'Frequency range of file %s fully contained in range of file %s => no merging of files!'%(fn1, fn2)
            return
        else:
            lo_eos = [eo2, eo1]

    #find sorting indices for obtaining strictly increasing frequencies:
    inc_freq_idxs_lower = np.array(lo_eos[0].frequencies()).argsort()
    inc_freq_idxs_upper = np.array(lo_eos[1].frequencies()).argsort()

    #determine overlap in frequencies
    upper_bound = max(lo_eos[0].frequencies())
    lower_bound = min(lo_eos[1].frequencies())

    overlap_mid_freq = 0.5*(upper_bound + lower_bound)

    if merge_frequency is not None:
        try:
            merge_frequency = float(merge_frequency)
        except:
            print 'could not read "merge frequency" argument (float expected)...taking mean of frequency overlap instead: %f Hz'%overlap_mid_freq
            merge_frequency = overlap_mid_freq
    else:
        merge_frequency = overlap_mid_freq


    #find indices for all frequencies from the frequency lists, which are below(lower part) or above (upper part) of the merge frequency - use sorted frequency lists !:

    lower_idxs = list(np.where( np.array(lo_eos[0].frequencies())[inc_freq_idxs_lower] <= merge_frequency)[0])
    upper_idxs = list(np.where( np.array(lo_eos[1].frequencies())[inc_freq_idxs_upper]  > merge_frequency)[0])


    #total of frequencies in new edi object
    n_total_freqs = len(lower_idxs) + len(upper_idxs)

    #------------
    # fill data fields

    eom.z = np.zeros((eom.n_freqs(),2,2),dtype=np.complex)
    eom.zerr = np.zeros((eom.n_freqs(),2,2),dtype=np.float)

    #check, if tipper exists for both files:
    if (eo1.tipper  is not None ) and (eo2.tipper  is not None ):
        eom.tipper = np.zeros((eom.n_freqs(),1,2),dtype=np.complex)
        eom.tippererr = np.zeros((eom.n_freqs(),1,2),dtype=np.float)

    freq_idx = 0
    zrot = []
    lo_freqs = []
    #first read out z, zerr (and tipper) of lower freq. edi object:
    in_z_lower = eo1.z[inc_freq_idxs_lower]
    in_zerr_lower = eo1.zerr[inc_freq_idxs_lower]
    if eom.tipper is not None:
        in_t_lower = eo1.tipper[inc_freq_idxs_lower]
        in_terr_lower = eo1.tippererr[inc_freq_idxs_lower]

    for li in lower_idxs:
        lo_freqs.append(np.array(lo_eos[0].frequencies())[inc_freq_idxs_lower][li])
        eom.z[freq_idx,:,:] = in_z_lower[li,:,:]
        eom.zerr[freq_idx,:,:] = in_zerr_lower[li,:,:]
        if eom.tipper is not None:
            eom.tipper[freq_idx,:,:] =  in_t_lower[li,:,:]
            eom.tippererr[freq_idx,:,:] =  in_terr_lower[li,:,:]
        try:
            zrot.append(eo1.zrot[freq_idx])
        except:
            zrot.append(0.)

        freq_idx += 1

    #then read upper freq. edi object:
    in_z_upper = eo2.z[inc_freq_idxs_upper]
    in_zerr_upper = eo2.zerr[inc_freq_idxs_upper]
    if eom.tipper is not None:
        in_t_upper = eo2.tipper[inc_freq_idxs_upper]
        in_terr_upper = eo2.tippererr[inc_freq_idxs_upper]

    for ui in upper_idxs:
        lo_freqs.append(np.array(lo_eos[1].frequencies())[inc_freq_idxs_upper][ui])
        eom.z[freq_idx,:,:] = in_z_upper[ui,:,:]
        eom.zerr[freq_idx,:,:] = in_zerr_upper[ui,:,:]

        if eom.tipper is not None:
            eom.tipper[freq_idx,:,:] =  in_t_upper[ui,:,:]
            eom.tippererr[freq_idx,:,:] =  in_terr_upper[ui,:,:]
        try:
            zrot.append(eo1.zrot[freq_idx])
        except:
            zrot.append(0.)

        freq_idx += 1

    eom.zrot = zrot
    eom.freq = lo_freqs


    #------------
    # fill header information

    #I) HEAD
    head1 = dict((k.lower(),v) for k,v in eo1.head.items())
    head2 = dict((k.lower(),v) for k,v in eo2.head.items())

    so_headsections = set(head1.keys() + head2.keys())

    head_dict = {}
    for element in so_headsections:

        if (element in head1) and (element not in head2):
            head_dict[element] = str(head1[element])
            continue
        if  (element in head2) and (element not in head1):
            head_dict[element] = str(head2[element])
            continue
        if head1[element] == head2[element]:
            head_dict[element] = str(head2[element])
            continue

        if element in ['lat','long','lon','latitude','longitude', 'elevation','elev','ele']:
            try:
                head_dict[element] = 0.5 * (float(head1[element]) + float(head2[element]))
            except:
                raise MTexceptions.MTpyError_edi_file('Cannot merge files: wrong format of "%s" coordinate'%element)
            continue

        if element == 'dataid':
            head_dict[element] = head1[element]+'+'+head2[element]
            continue

        if 'date' in element:
            dateformat1 = '%d/%m/%y'
            dateformat2 = '%d.%m.%y'
            dateformat3 = '%d.%m.%Y'
            date1 = None
            try:
                date1 = calendar.timegm(time.strptime(head1[element], dateformat1))
                date2 = calendar.timegm(time.strptime(head2[element], dateformat1))
            except:
                pass
            try:
                date1 = calendar.timegm(time.strptime(head1[element], dateformat2))
                date2 = calendar.timegm(time.strptime(head2[element], dateformat2))
            except:
                pass
            try:
                date1 = calendar.timegm(time.strptime(head1[element], dateformat3))
                date2 = calendar.timegm(time.strptime(head2[element], dateformat3))
            except:
                pass
            if date1 is None:
                raise MTexceptions.MTpyError_edi_file('Cannot merge file, because data format is not understood: %s=%s|%s'%(element,head1[element],head2[element]))


            if element in ['acqdate']:
                date = min(date1, date2 )

            elif element in ['enddate']:
                date = max(date1, date2  )

            elif element in ['filedate']:
                date = calendar.timegm(time.gmtime())

            datetuple = time.gmtime(date)
            head_dict[element] = '%02i/%02i/%2i'%(datetuple[2],datetuple[1],datetuple[0])

    eom.head = head_dict

    #II) INFO
    info1 = dict((k.lower(),v) for k,v in eo1.info_dict.items())
    info2 = dict((k.lower(),v) for k,v in eo2.info_dict.items())

    so_infosections = set(info1.keys() + info2.keys())

    info_dict = {}
    info_dict['merge_frequency'] = merge_frequency

    for element in so_infosections:

        if (element in info1) and (element not in info2):
            info_dict[element] = str(info1[element])
            continue
        if  (element in info2) and (element not in info1):
            info_dict[element] = str(info2[element])
            continue
        if info1[element] == info2[element]:
            info_dict[element] = str(info2[element])
            continue

        if element in ['lat','long','lon','latitude','longitude', 'elevation','elev','ele']:
            try:
                info_dict[element] = 0.5 * (float(info1[element]) + float(info2[element]))
            except:
                raise MTexceptions.MTpyError_edi_file('Cannot merge files: wrong format of "%s" coordinate'%element)
            continue

        if element == 'dataid':
            info_dict[element] = info1[element]+'_merged_with_'+info2[element]
            continue

        if 'date' in element:
            dateformat1 = '%d/%m/%y'
            dateformat2 = '%d.%m.%y'
            dateformat3 = '%d.%m.%Y'
            date1 = None
            try:
                date1 = calendar.timegm(time.strptime(info1[element], dateformat1))
                date2 = calendar.timegm(time.strptime(info2[element], dateformat1))
            except:
                pass
            try:
                date1 = calendar.timegm(time.strptime(info1[element], dateformat2))
                date2 = calendar.timegm(time.strptime(info2[element], dateformat2))
            except:
                pass
            try:
                date1 = calendar.timegm(time.strptime(info1[element], dateformat3))
                date2 = calendar.timegm(time.strptime(info2[element], dateformat3))
            except:
                pass
            if date1 is None:
                raise MTexceptions.MTpyError_edi_file('Cannot merge file, because data format is not understood: %s=%s|%s'%(element,info1[element],info2[element]))


            if element in ['acqdate']:
                date = min(date1, date2 )

            elif element in ['enddate']:
                date = max(date1, date2  )

            elif element in ['filedate']:
                date = calendar.timegm(time.gmtime())

            #arbitrarily choisen to take information from low frequency file:
            else:
                date = date1

            datetuple = time.gmtime(date)
            info_dict[element] = '%02i/%02i/%2i'%(datetuple[2],datetuple[1],datetuple[0])
            continue

        if element == 'station':
            info_dict['station'] = info1[element] + '+' + info2[element]


    eom.info_dict = info_dict

    #III) DEFINEMEAS
    dmeas1 = dict((k.lower(),v) for k,v in eo1.definemeas.items())
    dmeas2 = dict((k.lower(),v) for k,v in eo2.definemeas.items())

    so_dmeassections = set(dmeas1.keys() + dmeas2.keys())

    dmeas_dict = {}

    for element in so_dmeassections:

        if element == 'refloc':
            dmeas_dict[element] = ''
            continue

        if (element in dmeas1) and (element not in dmeas2):
            dmeas_dict[element] = str(dmeas1[element])
            continue
        if  (element in dmeas2) and (element not in dmeas1):
            dmeas_dict[element] = str(dmeas2[element])
            continue
        if dmeas1[element] == dmeas2[element]:
            dmeas_dict[element] = str(dmeas2[element])
            continue

        if 'lat' in element or 'lon' in element or 'elev' in element:
            try:
                dmeas_dict[element] = 0.5 * (float(dmeas1[element]) + float(dmeas2[element]))
            except:
                raise MTexceptions.MTpyError_edi_file('Cannot merge files: wrong format of "%s" coordinate'%element)
            continue



    eom.definemeas = dmeas_dict

    #take hmeas/dmeas section directly from file1:

    eom.hmeas_emeas = eo1.hmeas_emeas

    #IV) MTSECT

    msec1 = dict((k.lower(),v) for k,v in eo1.mtsect.items())
    msec2 = dict((k.lower(),v) for k,v in eo2.mtsect.items())

    so_msecsections = set(msec1.keys() + msec2.keys())

    msec_dict = {}

    for element in so_msecsections:
        #completely unimportant, kept just for the sake of the format:
        if element in ['ex','ey','hx','hy','hz','bx','by','bz']:
            msec_dict[element] = msec1[element]
        if element == 'nfreq':
            msec_dict[element] = eom.n_freqs()
        if element == 'sectid':
            msec_dict[element] = msec1[element]+'+'+msec2[element]


    eom.mtsect = msec_dict


    if out_fn is not None:
        dirname = op.dirname(op.abspath(op.join('.',out_fn)))
        fn = op.basename(op.abspath(op.join('.',out_fn)))
        if not op.isdir(dirname):
            try:
                os.makedirs(dirname)
                out_fn = op.join(dirname,fn)
            except:
                out_fn = None
        else:
            out_fn = op.join(dirname,fn)


        eom.writefile(out_fn)


    return eom, out_fn


def validate_edifile(fn):
    """
        Validate an EDI file following MTpy standard.

        Return boolean result.
    """

    edi_object = Edi()

    try:
        edi_object.readfile(fn)
        return True
    except:
        return False


def rotate_edifile(fn, angle, out_fn = None):
    """
        Rotate data contents (Z and Tipper) of an EDI file and write it to a new EDI file.
        (Use a script with consecutive renaming of the file for in place rotation. MTpy does not overwrite.)

        Input:
        - angle/list of angles for the rotation

        optional input:
        - name of output EDI file

        Output:
        - full path to the new (rotated) EDI file
    """

    ediobject = Edi()

    ediobject.readfile(fn)

    ediobject.rotate(angle)

    if out_fn is not None:
        dirname = op.dirname(op.abspath(op.join('.',out_fn)))
        fn = op.basename(op.abspath(op.join('.',out_fn)))
        if not op.isdir(dirname):
            try:
                os.makedirs(dirname)
                out_fn = op.join(dirname,fn)
            except:
                out_fn = None
        else:
            out_fn = op.join(dirname,fn)


    ediobject.writefile(out_fn)


    return out_fn



def _generate_edifile_string(edidict):
    """
    Generate a string to write out to an EDI file.

    Reading in information from an edi file dictionary. Using the standard sections:
    HEAD, INFO, DEFINEMEAS, HMEAS_EMEAS, MTSECT, ZROT, FREQ, Z, TIPPER

    Can be extended later on...

    """
    # define section heads explicitely instead of iteration over the dictionary for getting the correct order!
    lo_sectionheads = ['HEAD', 'INFO', 'DEFINEMEAS', 'HMEAS_EMEAS', 'MTSECT', 'ZROT', 'FREQ', 'Z', 'TIPPER']

    edistring = ''
    stationname = None
    ZROTflag = 0

    if len(edidict.keys()) == 0:
        raise MTexceptions.MTpyError_edi_file('Cannot generate string from empty EDI dictionary. Fill dict or read in file first')


    for sectionhead in lo_sectionheads:

        if sectionhead == 'HEAD':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "HEAD" missing!')
            edistring += '>HEAD\n'
            head_dict = edidict['HEAD']
            for k in  sorted(head_dict.iterkeys()):
                v = str(head_dict[k])
                if len(v) == 0 or len(v.split()) > 1:
                    edistring += '\t%s=""\n'%(k.upper())
                else:
                    try:
                        v = v.upper()
                    except:
                        pass
                    edistring += '\t%s=%s\n'%(k.upper(),v)

        if sectionhead == 'INFO':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "INFO" missing!')
            info_dict = edidict['INFO']
            info_dict = dict((k.lower(),v) for k,v in info_dict.items())

            if 'max lines' in info_dict:
                edistring += '>INFO  MAX LINES=%i\n'%(int(float(info_dict.pop('max lines'))))
            else:
                edistring += '>INFO \n'

            for k in sorted(info_dict.iterkeys()):
                v = str(info_dict[k])
                if len(v) == 0  or len(v.split()) > 1:
                    edistring += '\t%s: ""\n'%(k)
                else:
                    edistring += '\t%s: %s\n'%(k,v)

                #get station name (to be returned aside with the edistring, allowing for proper naming of output file)
                if k == 'station':
                    stationname = v.upper()


        if sectionhead == 'DEFINEMEAS':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "DEFINEMEAS" missing!')
            defm_dict = edidict['DEFINEMEAS']
            defm_dict = dict((k.upper(),v) for k,v in defm_dict.items())

            edistring += '>=DEFINEMEAS \n'

            for k in sorted(defm_dict.iterkeys()):
                v = str(defm_dict[k])
                if len(v) == 0  or len(v.split()) > 1:
                    edistring += '\t%s=""\n'%(k)
                else:
                    edistring += '\t%s=%s\n'%(k,v)


        if sectionhead == 'HMEAS_EMEAS':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required subsection "HMEAS_EMEAS" missing!')
            lo_hemeas = edidict['HMEAS_EMEAS']

            for hemeas in lo_hemeas:
                edistring += ('>'+' '.join(hemeas)+'\n').upper()


        if sectionhead == 'MTSECT':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "MTSECT" missing!')
            mtsct_dict = edidict['MTSECT']
            mtsct_dict = dict((k.upper(),v) for k,v in mtsct_dict.items())

            edistring += '>=MTSECT \n'

            for k in sorted(mtsct_dict.iterkeys()):
                v = str(mtsct_dict[k])
                if len(v) == 0 or len(v.split()) > 1:
                    edistring += '\t%s=""\n'%(k)
                else:
                    edistring += '\t%s=%s\n'%(k,v)


        if sectionhead == 'FREQ':
            if not sectionhead in edidict:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "FREQ" missing!')
            lo_freqs = edidict['FREQ']

            edistring+= '>FREQ // %i\n'%(len(lo_freqs))

            for i,freq in enumerate(lo_freqs):
                edistring += '\t%E'%(freq)
                if (i+1)%5 == 0 and (i != len(lo_freqs) - 1) and i > 0:
                    edistring += '\n'

        if sectionhead == 'ZROT':

            try:
                lo_rots = edidict['ZROT']
            except:
                continue

            edistring+= '>ZROT // %i\n'%(len(lo_rots))

            for i,angle in enumerate(lo_rots):
                edistring += '\t%E'%(angle)
                if (i+1)%5 == 0 and (i != len(lo_rots) - 1) and i > 0:
                    edistring += '\n'

            ZROTflag = 1

        if sectionhead == 'Z':

            compstrings = ['ZXX','ZXY','ZYX','ZYY']
            Z_entries = ['R','I','.VAR']

            try:
                z_dict = edidict['Z']
            except:
                raise MTexceptions.MTpyError_edi_file('Cannot write file - required section "Z" missing!')


            for idx_comp,comp in enumerate(compstrings):
                for idx_zentry,zentry in enumerate(Z_entries):
                    section = comp + zentry
                    if not section in z_dict:
                        raise MTexceptions.MTpyError_edi_file('Cannot write file - required subsection "%s" missing!'%(section))
                    lo_vals = z_dict[section]

                    if ZROTflag == 1:
                        edistring += '>%s ROT=ZROT // %i\n'%(section,len(lo_freqs))
                    else:
                        edistring += '>%s // %i\n'%(section,len(lo_freqs))

                    for i,val in enumerate(lo_vals):
                        edistring += '\t%E'%(float(val))
                        if (i+1)%5 == 0 and (i != len(lo_vals) - 1) and i > 0:
                            edistring += '\n'
                    edistring += '\n'


        if sectionhead == 'TIPPER':

            compstrings = ['TX','TY']
            T_entries = ['R','I','VAR']
            Tout_entries = ['R.EXP','I.EXP','VAR.EXP']

            try:
                t_dict = edidict['TIPPER']
                if t_dict == None:
                    continue
            except:
                continue

            for idx_comp,comp in enumerate(compstrings):
                for idx_tentry,tentry in enumerate(T_entries):
                    section = comp + tentry
                    outsection = comp + Tout_entries[idx_tentry]
                    if not section in t_dict:
                        raise MTexceptions.MTpyError_edi_file('Cannot write file - required subsection "%s" missing!'%(section))
                    lo_vals = t_dict[section]

                    if ZROTflag == 1:
                        edistring += '>%s ROT=ZROT // %i\n'%(outsection,len(lo_freqs))
                    else:
                        edistring += '>%s // %i\n'%(outsection,len(lo_freqs))

                    for i,val in enumerate(lo_vals):
                        edistring += '\t%E'%(float(val))
                        if (i+1)%5 == 0 and (i != len(lo_vals) - 1) and i > 0:
                            edistring += '\n'

                    edistring += '\n'



        edistring += '\n'


    edistring += '>END\n'


    return edistring.expandtabs(4), stationname



def _cut_sectionstring(edistring,sectionhead):
    """
        Cut an edi-string for the specified section.

        Input:
        - name of the section

        Output:
        - string : part of the raw edi-string containing starting at the head of the section and ends at beginnig of the next section.
    """

    #in this case, several blocks have to be handled together, therefore, a simple cut to the next block start does not work:
    if sectionhead.upper() == 'HMEAS_EMEAS':

        lo_start_idxs = [m.start() for m in re.finditer('>[HE]MEAS', edistring) ]
        if len(lo_start_idxs) == 0 :
            raise

        start_idx = lo_start_idxs[0]

        end_idx = edistring[(lo_start_idxs[-1]+1):].upper().find('>') + lo_start_idxs[-1]

        hmeas_emeas_string = edistring[start_idx:end_idx]

        if len(hmeas_emeas_string) == 0:
            raise

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


def _validate_edifile_string(edistring):
    """
        Read the file as string and check, if blocks 'HEAD, INFO, =DEFINEMEAS, =MTSECT, FREQ, (Z,) END' are present. If 'Z' is missing, check for 'spectra' or 'rho'/'phs'!

        Within the blocks look for mandatory entries:
        HEAD: 'DATAID'
        INFO: None
        DEFINEMEAS: subblocks 'HMEAS, EMEAS'
                    ('REFLAT, REFLONG, REFELEV' have to be present for measured data though)
        MTSECT: 'NFREQ'
        FREQ: non empty list

        Z: at least one component xx, yy, xy, yx ; real, imag and var ; containing a non-empty list
        Otherwise check for presence of 'RHO'/'PHS' OR 'spectra'


    """
    isvalid = False
    found = 1

    #adding 1 to position of find to correct for possible occurrence at position 0 )
    found *= np.sign(edistring.upper().find('>HEAD') + 1 )
    found *= np.sign(edistring.upper().find('DATAID') + 1 )
    found *= np.sign(edistring.upper().find('>HMEAS') + 1 )
    found *= np.sign(edistring.upper().find('>EMEAS') + 1 )
    found *= np.sign(edistring.upper().find('NFREQ') + 1 )
    found *= np.sign(edistring.upper().find('>END') + 1 )
    found *= np.sign(edistring.upper().find('>=DEFINEMEAS') + 1 )

    #allow spectral information as alternative:
    if np.sign(edistring.upper().find('>FREQ') + 1 ) == 0:
        if np.sign(edistring.upper().find('>SPECTRA') + 1 ) == 0 :
            found *= 0
    if np.sign(edistring.upper().find('>=MTSECT') + 1 ) == 0:
        if np.sign(edistring.upper().find('>=SPECTRASECT') + 1 ) == 0:
            found *= 0


    if found < 1 :
        print 'Could not find all mandatory sections for a valid EDI file!\n (Most basic version must contain: "HEAD, INFO, =DEFINEMEAS, =MTSECT or =SPECTRASECT, FREQ or SPECTRA, (Z,) END") '
        return False




    #checking for non empty frequency list:
    freq_start_idx = edistring.upper().find('>FREQ')
    next_block_start = edistring.upper().find('>',freq_start_idx + 1)
    string_dummy_2 = edistring[freq_start_idx:next_block_start]
    lo_string_dummy_2 = string_dummy_2.strip().split()
    #check, if there are actually one/some valid numbers:
    n_numbers = 0
    for i in lo_string_dummy_2:
        try:
            n = float(i)
            n_numbers +=1
        except:
            continue

    if n_numbers == 0:
        print  MTexceptions.MTpyError_edi_file('Problem in FREQ block: no frequencies found...checking for spectra instead')
        #found *= 0
    #Check for data entry following priority:
    # 1. Z
    z_found = 0
    rhophi_found = 0
    spectra_found = 0

    compstrings = ['ZXX','ZXY','ZYX','ZYY']
    Z_entries = ['R','I','.VAR']

    for comp in compstrings:
        n_entries = 0

        for zentry in Z_entries:
            searchstring = '>'+comp+zentry
            z_comp_start_idx = edistring.upper().find(searchstring)
            if z_comp_start_idx < 0:
                continue
            #found *= np.sign(z_comp_start_idx + 1 )
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
                continue

            if zentry in ['R','I']:
                n_entries += 1
        if n_entries > 1:
            z_found += 1


    # If no Z entry is found continue searching for RhoPhase information
    # 2. RHO,PHS
    if z_found == 0:

        rhophistrings = ['RHO','PHS']
        compstrings = ['XX','XY','YX','YY']

        for comp in compstrings:
            n_entries = 0
            for rp in rhophistrings:
                sectionhead = rp + comp
                try:
                    temp_string = _cut_sectionstring(edistring,sectionhead)
                    lo_vals = []
                    t0 = temp_string.strip().split('\n')[0]
                    n_dummy = int(float(t0.split('//')[1].strip()))
                    t1 = temp_string.strip().split('\n')[1:]
                    for j in t1:
                        lo_j = j.strip().split()
                        for k in lo_j:
                            try:
                                lo_vals.append(float(k))
                            except:
                                pass
                    if len(lo_vals) == 0:
                        raise
                except:
                    continue

                n_entries += 1
            if n_entries > 1 :
                rhophi_found += 1

    # If neither Z nor RHO/PHS  entries are found continue searching for spectra information
    # 3. spectra
    if z_found == 0 and rhophi_found == 0:

        spectrasect = _cut_sectionstring(edistring, '=SPECTRASECT')
        if len(spectrasect) == 0 :
            found *=0
        dummy4 = spectrasect.upper().find('NCHAN')
        n_chan = int(float(spectrasect[dummy4:].strip().split('=')[1].split()[0]))

        if n_chan not in [4,5,6,7]:
            found *= 0 

        dummy5 = spectrasect.upper().find('NFREQ')
        n_freq = int(float(spectrasect[dummy5:].strip().split('=')[1].split()[0]))
       
        firstspectrum = _cut_sectionstring(edistring, 'SPECTRA')
        if len(firstspectrum) == 0 :
            found *=0

        no_values = int(float(firstspectrum.split('\n')[0].strip().split('//')[1]))

        if firstspectrum.upper().find('FREQ') <0 :
            found *= 0
        
        if not n_chan**2 == no_values:
            found *= 0

        lo_valuelines = firstspectrum.split('\n')[1:]
        dummy6 = ''
        for i in lo_valuelines:
            dummy6 += i

        if not len(dummy6.split()) == no_values:
            found *= 0

        if not edistring.upper().count('>SPECTRA') ==  n_freq:
            found *= 0
        if found > 0:
            print 'Found spectra data'
            spectra_found = 1
 
    if z_found == 0 and rhophi_found == 0 and spectra_found == 0 :
        print 'ERROR - no data found in terms of "Z" or "RHO/PHS" or "SPECTRA" - reading of multiple stations is not supported (yet)!'
        found *= 0

    if found > 0: isvalid = True

    return isvalid



def _build_id_channel_dict(lo_hmeas_emeas):

    id_dict = {}

    for line in lo_hmeas_emeas:
        if len(''.join(line).strip()) == 0:
            continue

        channel = _find_key_value('CHTYPE','=',' '.join(line),valuelength=2)
        ID = _find_key_value('ID','=',' '.join(line))

        id_dict[ID] = channel

    return id_dict

def _find_key_value(key, separator, instring, valuelength=None):

    line = instring.strip().split()

    #loop over list/line elements
    for idx, element in enumerate(line):
        #if keyword is not found in entry:
        if element.upper().find(key.upper()) < 0:
            continue
        #else check, if the separator is present in the same element (equiv. to no spacing)
        if element.upper().find(separator) >= 0:

            #, if the splitting worked out 
            if len(element.split(separator)) == 2 :
                #if all fine until now, read in the part after the separator as value
                value = element.split(separator)[1].upper()
                
                #if the separator was at the end of the element, read the next element as value
                if  len(element.split(separator)[1]) == 0 :
                    value = line[idx+1]
        
        #else, the separator is in the next element
        else:
            #check, if the next line is entirely defined by separator -> value must be one later
            if line[idx+1] == separator:
                value = line[idx+2]
            #else, cut off the separator from the value
            else:
                value = line[idx+1].split(separator)[1]
                #check for correct length of value, if specified 
        if valuelength is not None:
            if len(value) != valuelength :
                continue



    return value

