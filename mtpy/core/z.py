#!/usr/bin/env python

"""
mtpy/mtpy/core/z.py

Contains classes and functions for handling impedance tensors (Z). 
 
    Class:
    "Z" contains information about an impedance tensor Z. 

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
import time, calendar 

import mtpy.core.edi as MTedi 
reload (MTedi)

import mtpy.utils.format as MTformat
reload(MTformat)

import mtpy.utils.exceptions as MTexceptions
reload(MTexceptions)


#=================================================================


#------------------------
class Z(object):
    """
        Z class - generates an impedance tensor (Z-) object.

        Methods  include reading and writing from and to edi-objects, rotations/combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        'get' and 'set' for all edi file sections
        
        Z is a complex array of the form (n_frequencies, 2, 2), 
        with indices in the following order: 
            Zxx: (0,0) - Zxy: (0,1) - Zyx: (1,0) - Zyy: (1,1)   


    """

    def __init__(self, z_array = None, zerr_array = None, edi_object = None):
    

        self.z = None
        self.zerr = None
              
        try:
            if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                if z_array.dtype in ['complex', 'float']:
                    self.z = z_array
        except:
            pass

        try:
            if len(zerr_array.shape) == 3 and zerr_array.shape[1:3] == (2,2):
                if zerr_array.dtype in ['float']:
                    self.zerr = zerr_array
        except:
            pass
            

        self.edi_object = None
        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            self.z = edi_object.z
            if  edi_object.zerr  is not None:
                self.zerr = edi_object.zerr
        try:
            if len(self.z) != len(zerr):
                self.zerr = None
        except:
            pass

        self.frequencies = None
        self.rotation_angle = 0.
        self.static_shift = None
        self.distortion = None
        self.invariants = None
        self.no_distortion = None
        self.no_ss = None
        self.no_ss_no_distortion = None
        self.ellipticity = None

        #intialise variables, if Z has been generated non-empty
        self.update()


        
    def set_z(self, z_array):
        pass

    def set_zerr(self, zerr_array):
        pass


    def real(self):
        return np.real(self.z)

        
    def set_real(self, real_array):
        pass


    def imag(self):
        return np.imag(self.z)

        
    def set_imag(self, imag_array):
        pass

    def rho(self):
        
        
        
    def set_rho(self, rho_array):
        pass

    def phi(self):
        pass
        
    def set_phi(self, phi_array):
        pass



    def inverse(self):
        
        inverse = self.z.copy()
        for mat in range(len(inverse)):
            inverse[mat,:,:] = np.array( (np.matrix(self.z[mat,:,:])).I )

        return inverse




    def rotate(self, alpha):
        pass

    def remove_ss(self, rho_x = 1., rho_y = 1.):
        pass


    def _remove_distortion(self):
        pass


#     def update(self):

#         #0. check, if information available in z, edi_object and/or rho/phi arrays:
#         if self.z == None and self.edi_object == None and (self.rho == None or self.phi == None):
#             print 'No information to update Z object from -- skipping update'
#             return

#         if self.edi_object is not None:
#             if not isinstance(edi_object,MTedi.Edi):
#                 print '"edi_object" is not a valid instance of the Edi class'
#                 wrong_edi_object = self.edi_object
#                 self.edi_object == None
#                 return wrong_edi_object
        
#             try:
#                 new_z = edi_object.z
#                 if  len(new_z.shape) != 3 or new_z.shape[1:3] != (2,2):
#                     raise
#                 if new_z.dtype not in ['complex', 'float']:
#                     raise
#                 if  edi_object.zerr is not None:
#                     new_zerr = edi_object.zerr
#                     if  len(new_zerr.shape) != 3 or new_zerr.shape[1:3] != (1,2):
#                         raise
#                     if new_zerr.dtype not in ['float']:
#                         raise
#                     if len(new_zerr) != len(new_z):
#                         raise

#                 self.z = new_z
#                 self.zerr = new_zerr

#             except:
#                 print '"edi_object" is not a valid instance of the Edi class'
#                 wrong_edi_object = self.edi_object
#                 self.edi_object == None
#                 return wrong_edi_object

# ...



#         #1. Use rho/phi data to set up Z array, if not existing
#         if self.z is None:


        
#         #2. otherwise start directly with Z information:
#         try:
#             if len(z_array.shape) != 3 or z_array.shape[1:3] != (2,2):
#                 raise
#             if z_array.dtype not in ['complex', 'float']:
#                 raise
#         except:
#             try:
#                 print 'Z array is in wrong shape: %s'%(str(self.z.shape))
#             except:
#                 print 'Error - shape of Z array cannot be determined - check, if Z is array!'

    



#         pass

#------------------------

class Tipper(object):
    """
        Tipper class - generates a Tipper-object.



    """

    def __init__(self, tipper_array = None, tippererr_array = None, edi_object = None):
    

        self.tipper = None        
        self.tippererr = None
        try:
            if len(tipper_array.shape) == 3 and tipper_array.shape[1:3] == (1,2):
                if tipper_array.dtype in ['complex', 'float']:
                    self.tipper = tipper_array
        except:
            pass
        try:
            if len(tippererr_array.shape) == 3 and tippererr_array.shape[1:3] == (1,2):
                if tippererr_array.dtype in ['float']:
                    self.tippererr = tippererr_array
        except:
            pass


        self.edi_object = None
        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            if edi_object.tipper is not None:
                self.tipper = edi_object.tipper
            if edi_object.tippererr is not None:
                self.tippererr = edi_object.tippererr

        try:
            if len(self.tipper) != len(tippererr):
                self.tippererr = None
        except:
            pass

        self.frequencies = None

        self.real = None
        self.imag = None
        self.amplitude = None
        self.phase = None
        self.rotation_angle = 0.


    def rotate(self, alpha):
        pass


#------------------------


def rotate_z(z_array, alpha, zerr_array = None):

    z_object = Z(z_array=z_array, zerr_array=zerr_array)

    z_object.rotate(alpha)

    return z_object.z, z_object.zerr


def remove_distortion(z_array, zerr_array = None ):

    z_object = Z(z_array=z_array, zerr_array=zerr_array)

    return z_object.no_distortion, z.zerr


def remove_static_shift(z_array, zerr_array = None):

    z_object = Z(z_array=z_array, zerr_array=zerr_array)

    return z_object.no_ss, z_object.zerr_array


def z2rhophi(z_array):
    pass



def rotate_tipper(tipper_array, alpha, tippererr_array = None):

    tipper_object = Tipper(tipper_array=tipper_array, tippererr_array=tippererr_array)

    tipper_object.rotate(alpha)

    return tipper_object.tipper, tipper_object.tippererr


