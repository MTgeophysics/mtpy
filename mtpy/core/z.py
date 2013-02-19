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
import math, cmath
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
            
        self.frequencies = None
        self.edi_object = None

        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            self.frequencies = edi_object.frequencies
            self.z = edi_object.z
            if  edi_object.zerr  is not None:
                self.zerr = edi_object.zerr
        try:
            if len(self.z) != len(zerr):
                self.zerr = None
        except:
            pass

        self.rotation_angle = 0.
        self.static_shift = None
        self.distortion = None
        self.invariants = None
        self.no_distortion = None
        self.no_ss = None
        self.no_ss_no_distortion = None
        self.ellipticity = None

    def set_edi_object(self, edi_object):

        if not isinstance(edi_object,MTedi.Edi):
            print 'Object is not a valid Edi instance - Z object not updated'
            return


        self.edi_object = edi_object
        self.frequencies = edi_object.frequencies
        self.z = edi_object.z
        self.zerr = edi_object.zerr
       
    def set_z(self, z_array):

        z_orig = self.z 

        if z_orig.shape != z_array.shape:
            print 'Error - shape of "z" array does not match shape of Z array: %s ; %s'%(str(z_array.shape),str(z_orig.shape))
            return

        self.z = z_array


    def set_zerr(self, zerr_array):

        zerr_orig = self.zerr 

        if zerr_orig.shape != zerr_array.shape:
            print 'Error - shape of "zerr" array does not match shape of Zerr array: %s ; %s'%(str(zerr_array.shape),str(zerr_orig.shape))
            return

        self.zerr = zerr_array


    def real(self):

        return np.real(self.z)

        
    def set_real(self, real_array):
        
        z_orig = self.z 

        if z_orig.shape != real_array.shape:
            print 'shape of "real" array does not match shape of Z array: %s ; %s'%(str(real_array.shape),str(z_orig.shape))
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return

        i = np.complex(0,1)

        z_new = real_array + i * self.imag() 

        self.z = z_new


    def imag(self):
        
        return np.imag(self.z)

        
    def set_imag(self, imag_array):

        z_orig = self.z 

        if z_orig.shape != imag_array.shape:
            print 'Error - shape of "imag" array does not match shape of Z array: %s ; %s'%(str(imag_array.shape),str(z_orig.shape))
            return
        
        #assert real array:
        if np.linalg.norm(np.imag(imag_array )) != 0 :
            print 'Error - array "imag" is not real valued !'
            return

        i = np.complex(0,1)

        z_new = self.real() + i * imag_array 

        self.z = z_new


    def rho(self):
        rho = np.zeros(self.z.shape)

        for idx_f in range(len(rho)):                         
            rho[idx_f,:,:] = np.abs(self.z[idx_f,:,:])

        return rho

        
    def phi(self):
        phi = np.zeros(self.z.shape)
        
        for idx_f in range(len(phi)):
            for i in range(2):
                for j in range(2):
                    phi[idx_f,i,j] = math.degrees(cmath.phase(self.z[idx_f,i,j]))

        return phi


    def set_rho_phi(self, rho_array, phi_array):

        z_new = self.z.copy() 

        if self.z.shape != rho_array.shape:
            print 'Error - shape of "rho" array does not match shape of Z array: %s ; %s'%(str(rho_array.shape),str(self.z.shape))
            return

        if self.z.shape != phi_array.shape:
            print 'Error - shape of "phi" array does not match shape of Z array: %s ; %s'%(str(phi_array.shape),str(self.z.shape))
            return
        
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
                    z_new[idx_f,i,j] = cmath.rect( rho_array[idx_f,i,j], math.radians(phi_array[idx_f,i,j] ))

        self.z = z_new


    def inverse(self):
        
        inverse = self.z.copy()
        for idx_f in range(len(inverse)):
            inverse[idx_f,:,:] = np.array( (np.matrix(self.z[idx_f,:,:])).I )

        return inverse


    def rotate(self, alpha):
        
        degreeangle = alpha%360
        angle = math.radians(degreeangle)
        
        cphi = np.cos(angle)
        sphi = np.sin(angle)

        z = self.z
        zerr = self.zerr

        z_new = z.copy()
        zerr_new = zerr.copy()

        for idx_freq in range(len(z)):

            z_new[idx_freq,0,0] = cphi**2 * z[idx_freq,0,0] + cphi*sphi*(z[idx_freq,0,1]+z[idx_freq,1,0]) + sphi**2 * z[idx_freq,1,1]
            z_new[idx_freq,0,1] = cphi**2 * z[idx_freq,0,1] + cphi*sphi*(z[idx_freq,1,1]-z[idx_freq,0,0]) - sphi**2 * z[idx_freq,1,0]
            z_new[idx_freq,1,0] = cphi**2 * z[idx_freq,1,0] + cphi*sphi*(z[idx_freq,1,1]-z[idx_freq,0,0]) - sphi**2 * z[idx_freq,0,1]
            z_new[idx_freq,1,1] = cphi**2 * z[idx_freq,1,1] + cphi*sphi*(-z[idx_freq,0,1]-z[idx_freq,1,0]) + sphi**2 *z[ idx_freq,0,0]

            zerr_new[idx_freq,0,0] = cphi**2 * zerr[idx_freq,0,0] + np.abs(cphi * sphi) * (zerr[idx_freq,0,1] + zerr[idx_freq,1,0]) + sphi**2 * zerr[idx_freq,1,1]

            zerr_new[idx_freq,0,1] = cphi**2 * zerr[idx_freq,0,1] + np.abs(cphi * sphi) * (zerr[idx_freq,1,1] + zerr[idx_freq,0,0]) + sphi**2 * zerr[idx_freq,1,0] 

            zerr_new[idx_freq,1,0] = cphi**2 * zerr[idx_freq,1,0] + np.abs(cphi * sphi) * (zerr[idx_freq,1,1] + zerr[idx_freq,0,0]) + sphi**2 * zerr[idx_freq,0,1] 

            zerr_new[idx_freq,1,1] = cphi**2 * zerr[idx_freq,1,1] + np.abs(cphi * sphi) * (zerr[idx_freq,0,1] + zerr[idx_freq,1,0]) + sphi**2 * zerr[idx_freq,0,0] 

        self.z = z_new
        self.zerr = zerr_new
        self.rotation_angle = (self.rotation_angle + degreeangle)%360


    def remove_ss(self, rho_x = 1., rho_y = 1.):
        
        return z_corrected, static_shift



    def remove_distortion(self):
        

        return z_corrected, distortion



    def remove_ss_and_distortion(self, rho_x = 1., rho_y = 1.):


        return z_corrected, static_shift, distortion
 


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

        self.frequencies = None
        self.edi_object = None

        if isinstance(edi_object,MTedi.Edi):
            self.edi_object = edi_object
            self.frequencies = edi_object.frequencies
            if edi_object.tipper is not None:
                self.tipper = edi_object.tipper
            if edi_object.tippererr is not None:
                self.tippererr = edi_object.tippererr

        try:
            if len(self.tipper) != len(tippererr):
                self.tippererr = None
        except:
            pass


        self.rotation_angle = 0.


    def set_edi_object(self, edi_object):

        if not isinstance(edi_object,MTedi.Edi):
            print 'Object is not a valid Edi instance - Tipper object not updated'
            return


        self.edi_object = edi_object
        self.frequencies = edi_object.frequencies
        
        try:
            if edi_object.tipper is None or edi_object.tippererr:
                raise
            
            tipper_new = edi_object.tipper
            tippererr_new = edi_object.tippererr

            if len(tipper_new) != len(tippererr_new):
                raise
            self.tipper = tipper_new
            self.tippererr = tippererr_new

        except:
            print 'Edi object does not contain correct tipper information - Tipper object not updated'

        
    def set_tipper(self, tipper_array):

        tipper_orig = self.tipper  

        if (tipper_orig is not None) and (tipper_orig.shape != tipper_array.shape):
            print 'Error - shape of "tipper" array does not match shape of tipper-array: %s ; %s'%(str(tipper_array.shape),str(tipper_orig.shape))
            return

        self.tipper = tipper_array


    def set_tippererr(self, tippererr_array):


        if (self.tippererr is not None) and (self.tippererr.shape != tippererr_array.shape):
            print 'Error - shape of "tippererr" array does not match shape of tippererr array: %s ; %s'%(str(tippererr_array.shape),str(self.tippererr.shape))
            return

        self.tippererr = tippererr_array


    def real(self):
        if self.tipper is None:
            print 'tipper array is None - cannot calculate real'
            return

        return np.real(self.tipper)

        
    def set_real(self, real_array):
        

        if (self.tipper is not None ) and (self.tipper.shape != real_array.shape):
            print 'shape of "real" array does not match shape of tipper array: %s ; %s'%(str(real_array.shape),str(self.tipper.shape))
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return


        i = np.complex(0,1)

        if self.tipper is not None:
            tipper_new = real_array + i * self.imag() 
        else:
            tipper_new = real_array

        self.tipper = tipper_new


    def imag(self):
        if self.tipper is None:
            print 'tipper array is None - cannot calculate imag'
            return

        return np.imag(self.tipper)

        
    def set_imag(self, imag_array):


        if (self.tipper is not None) and (self.tipper.shape != imag_array.shape):
            print 'Error - shape of "imag" array does not match shape of tipper array: %s ; %s'%(str(imag_array.shape),str(self.tipper.shape))
            return
        
        #assert real array:
        if np.linalg.norm(np.imag(imag_array )) != 0 :
            print 'Error - array "imag" is not real valued !'
            return

        i = np.complex(0,1)
        if self.tipper is not None:
            tipper_new = self.real() + i * imag_array 
        else:
            tipper_new = i * imag_array 


        self.tipper = tipper_new


    def rho(self):
        
        if self.tipper is None:
            print 'tipper array is None - cannot calculate rho'
            return

        rho = np.zeros(self.tipper.shape)

        for idx_f in range(len(rho)):                         
            rho[idx_f,:,:] = np.abs(self.tipper[idx_f,:,:])

        return rho

        
    def phi(self):
        if self.tipper is None:
            print 'tipper array is None - cannot calculate phi'
            return

        phi = np.zeros(self.tipper.shape)
        
        for idx_f in range(len(phi)):
                for j in range(2):
                    phi[idx_f,0,j] = math.degrees(cmath.phase(self.tipper[idx_f,0,j]))

        return phi


    def set_rho_phi(self, rho_array, phi_array):

        if self.tipper is not None: 

            tipper_new = self.tipper.copy() 

            if self.tipper.shape != rho_array.shape:
                print 'Error - shape of "rho" array does not match shape of tipper array: %s ; %s'%(str(rho_array.shape),str(self.tipper.shape))
                return

            if self.tipper.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of tipper array: %s ; %s'%(str(phi_array.shape),str(self.tipper.shape))
                return
        else:

            tipper_new = np.zeros(rho_array.shape,'complex')

            if rho_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of "rho" array: %s ; %s'%(str(phi_array.shape),str(rho_array.shape))
                return
       
        #assert real array:
        if np.linalg.norm(np.imag(rho_array )) != 0 :
            print 'Error - array "rho" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(rho_array)):
                for j in range(2):
                    tipper_new[idx_f,0,j] = cmath.rect( rho_array[idx_f,0,j], math.radians(phi_array[idx_f,0,j] ))

        self.tipper = tipper_new


    def rotate(self, alpha):

        if self.tipper is None :
            print 'tipper array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length 1 or same as len(tipper):
        if iterable(alpha) == 0:
            try:
                degreeangle = alpha%360
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            self.rotation_angle = degreeangle
            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.tipper]
        else:
            try:
                lo_angles = [ i%360 for i in alpha]
            except:
                print '"Angles" must be valid numbers (in degrees)'
                return
            
            self.rotation_angle = lo_angles

        if len(lo_angles) != len(tipper):
            print 'Wrong number Number of "angles" - need %i '%(len(tipper))
            self.rotation_angle = 0.
            return

        tipper_rot = self.tipper.copy()
        if self.tippererr is not None:
            tippererr_rot = self.tippererr.copy()

        for idx_freq in range(len(tipper_rot)):

            degreeangle = lo_angles[idx_freq]

            angle = math.radians(degreeangle)
            
            cphi = np.cos(angle)
            sphi = np.sin(angle)


            t_orig = self.tipper[idx_freq,:,:]

            tipper_rot[idx_freq,0,0] =  cphi * t_orig[0,0] + sphi * t_orig[0,1]
            tipper_rot[idx_freq,0,1] = -sphi * t_orig[0,0] + cphi * t_orig[0,1]

            if self.tippererr is not None:
                #absolute error propagation
                terr_orig = self.tippererr[idx_freq,:,:]

                tippererr_rot[idx_freq,0,0] = np.abs(cphi * terr_orig[0,0])  + np.abs(sphi * terr_orig[0,1])
                tippererr_rot[idx_freq,0,1] = np.abs(-sphi * terr_orig[0,0]) + np.abs(cphi * terr_orig[0,1])
 
        self.tipper = tipper_rot
        if self.tippererr is not None:
            self.tippererr = tippererr_rot


#------------------------


def rotate_z(z_array, alpha, zerr_array = None):

    z_object = _read_z_array(z_array, zerr_array)

    z_object.rotate(alpha)

    return z_object.z, z_object.zerr



def remove_distortion(z_array, zerr_array = None ):

    z_object = _read_z_array(z_array, zerr_array)
    
    z_object.remove_distortion()

    return z_object.no_distortion, z.zerr


def remove_static_shift(z_array, zerr_array = None):


    z_object = _read_z_array(z_array, zerr_array)

    z_object.remove_static_shift()

    return z_object.no_ss, z_object.zerr_array


def z2rhophi(z_array):
    
    z_object = _read_z_array(z_array)

    return z_object.rho(), z_object.phi()


def rotate_tipper(tipper_array, alpha, tippererr_array = None):

    tipper_object = _read_tipper_array(tipper_array, tippererr_array)

    tipper_object.rotate(alpha)

    return tipper_object.tipper, tipper_object.tippererr



def _read_z_array(z_array, zerr_array = None):

    try:
        z_object = Z( z_array=z_array, zerr_array=zerr_array )
    except:
        raise MTexceptions.MTpyError_Z('Cannot generate Z instance - check z-array dimensions/type: (N,2,2)/complex ; %s'%(str(z_array.shape)))

    return z_object


def _read_tipper_array(tipper_array, tippererr_array = None):

    try:
        tipper_object = tipper( tipper_array=tipper_array, tippererr_array=tippererr_array )
    except:
        raise MTexceptions.MTpyError_tipper('Cannot generate Tipper instance - check dimensions/type: (N,1,2)/complex ; %s'%(str(tipper_array.shape)))

    return tipper_object
