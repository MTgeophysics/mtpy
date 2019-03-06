# -*- coding: utf-8 -*-
"""
Created on Wed May 08 09:40:42 2013

Interpreted from matlab code written by Stephan Thiel 2005

@author: jpeacock
"""

import numpy as np


class Zinvariants:
    """
    calculates invariants from Weaver et al. [2000, 2003].  At the moment it 
    does not calculate the error for each invariant, only the strike.
    
    Arguments
    ----------
        **z_object** : type mtpy.core.z
                       needs to have attributes:
                           *z --> np.array((nf, 2, 2), dtype='complex')
                           
                           *z_err --> np.array((nf, 2, 2), dtype='real')
                           
                           *freq --> np.array(nf)
                           
        **z** : complex np.array(nf,2,2)
                impedance tensor array
                
        **z_err** : real np.array(nf,2,2)
                impedance tensor error array
                
        **freq** : np.array(nf)
                          array of freq cooresponding to the impedance 
                          tensor elements.
                          
    Attributes
    -----------
        **inv1**       : real off diaganol part normalizing factor
        
        **inv2**       : imaginary off diaganol normalizing factor
        
        **inv3**       : real anisotropy factor (range from [0,1])
        
        **inv4**       : imaginary anisotropy factor (range from [0,1])
        
        **inv5**       : suggests electric field twist
        
        **inv6**       : suggests in phase small scale distortion
        
        **inv7**       : suggests 3D structure
        
        **strike**     : strike angle (deg) assuming positive clockwise 0=N
        
        **strike_err** : strike angle error (deg)
        
        **q**          : dependent variable suggesting dimensionality
        
        
    Further reading
    ----------------
        Weaver, J. T., Agarwal, A. K., Lilley, F. E. M., 2000,
           Characterization of the magnetotelluric tensor in terms of its 
           invariants, Geophysical Journal International, 141, 321--336.
           
        Weaver, J. T., Agarwal, A. K., Lilley, F. E. M., 2003,
            The relationship between the magnetotelluric tensor invariants and 
            the phase tensor of Caldwell, Bibby and Brown, 
            presented at 3D Electromagnetics III, ASEG, paper 43.
           
        Lilley, F. E. M, 1998, Magnetotelluric tensor dcomposition: 1: Theory
            for a basic procedure, Geophysics, 63, 1885--1897.
           
        Lilley, F. E. M, 1998, Magnetotelluric tensor dcomposition: 2: Examples
            of a basic procedure, Geophysics, 63, 1898--1907.
           
        Szarka, L. and Menvielle, M., 1997, Analysis of rotational invariants 
            of the magnetotelluric impedance tensor, Geophysical Journal 
            International, 129, 133--142.
    
    """

    def __init__(self, z_object=None, z_array=None, z_err_array=None,
                 freq=None, rot_z=0):

        # --> read in z_object
        if z_object is not None:
            if z_object.freq is None:
                raise AttributeError('z_object needs to have attrtibute' + \
                                     'freq filled')

            # --> make the z_object an attribute
            self.z = z_object.z
            self.freq = z_object.freq
            self.z_err = z_object.z_err

        if freq is not None:
            self.freq = freq

        # --> if an array is input read it in and make it a z_object
        if z_array is not None:
            self.z = z_array.copy()

            assert len(freq) == len(self.z), \
                'length of freq is not the same as z'
        # --> if an array is input read it in and make it a z_object
        if z_err_array is not None:
            self.z_err = z_err_array.copy()

            assert len(freq) == len(self.z), \
                'length of freq is not the same as z'

        if self.freq is None:
            raise AttributeError('z_object needs to have attrtibute' + \
                                 'freq filled')

        # --> rotate data if desired
        self.rotate(rot_z)

        # compute the invariants
        self.compute_invariants()

    def compute_invariants(self):
        """
        Computes the invariants according to Weaver et al., [2000, 2003]
        
        Mostly used to plot Mohr's circles
        
        In a 1D case: rho = mu (inv1**2+inv2**2)/w & phi = arctan(inv2/inv1)
        
        Sets the invariants as attributes:
            **inv1**       : real off diaganol part normalizing factor
            
            **inv2**       : imaginary off diaganol normalizing factor
            
            **inv3**       : real anisotropy factor (range from [0,1])
            
            **inv4**       : imaginary anisotropy factor (range from [0,1])
            
            **inv5**       : suggests electric field twist
            
            **inv6**       : suggests in phase small scale distortion
            
            **inv7**       : suggests 3D structure
            
            **strike**     : strike angle (deg) assuming positive clockwise 0=N
            
            **strike_err** : strike angle error (deg)
            
            **q**          : dependent variable suggesting dimensionality
            
        """
        print("computing invariants")
        # get the length of z to initialize some empty arrays           
        nz = self.z.shape[0]

        # set some empty arrays to put stuff into
        self.inv1 = np.zeros(nz)
        self.inv2 = np.zeros(nz)
        self.inv3 = np.zeros(nz)
        self.inv4 = np.zeros(nz)
        self.inv5 = np.zeros(nz)
        self.inv6 = np.zeros(nz)
        self.inv7 = np.zeros(nz)
        self.q = np.zeros(nz)
        self.strike = np.zeros(nz)
        self.strike_err = np.zeros(nz)

        c_tf = np.all(self.z == 0.0)
        if c_tf == True:
            return
        # loop over each freq
        for ii in range(nz):
            # compute the mathematical invariants
            x1 = .5 * (self.z[ii, 0, 0].real + self.z[ii, 1, 1].real)  # trace
            x2 = .5 * (self.z[ii, 0, 1].real + self.z[ii, 1, 0].real)
            x3 = .5 * (self.z[ii, 0, 0].real - self.z[ii, 1, 1].real)
            x4 = .5 * (self.z[ii, 0, 1].real - self.z[ii, 1, 0].real)  # berd
            e1 = .5 * (self.z[ii, 0, 0].imag + self.z[ii, 1, 1].imag)  # trace
            e2 = .5 * (self.z[ii, 0, 1].imag + self.z[ii, 1, 0].imag)
            e3 = .5 * (self.z[ii, 0, 0].imag - self.z[ii, 1, 1].imag)
            e4 = .5 * (self.z[ii, 0, 1].imag - self.z[ii, 1, 0].imag)  # berd
            ex = x1 * e1 - x2 * e2 - x3 * e3 + x4 * e4

            if ex == 0.0:
                print('Could not compute invariants for {0:5e} Hz'.format(
                       self.freq[ii]))
                self.inv1[ii] = np.nan
                self.inv2[ii] = np.nan
                self.inv3[ii] = np.nan
                self.inv4[ii] = np.nan
                self.inv5[ii] = np.nan
                self.inv6[ii] = np.nan
                self.inv7[ii] = np.nan
                self.q[ii] = np.nan
                self.strike[ii] = np.nan
                self.strike_err[ii] = np.nan
            else:

                d12 = (x1 * e2 - x2 * e1) / ex
                d34 = (x3 * e4 - x4 * e3) / ex
                d13 = (x1 * e3 - x3 * e1) / ex
                d24 = (x2 * e4 - x4 * e2) / ex
                d41 = (x4 * e1 - x1 * e4) / ex
                d23 = (x2 * e3 - x3 * e2) / ex

                inv1 = np.sqrt(x4 ** 2 + x1 ** 2)
                inv2 = np.sqrt(e4 ** 2 + e1 ** 2)
                inv3 = np.sqrt(x2 ** 2 + x3 ** 2) / inv1
                inv4 = np.sqrt(e2 ** 2 + e3 ** 2) / inv2

                s41 = (x4 * e1 + x1 * e4) / ex

                inv5 = s41 * ex / (inv1 * inv2)
                inv6 = d41 * ex / (inv1 * inv2)

                q = np.sqrt((d12 - d34) ** 2 + (d13 + d24) ** 2)

                inv7 = (d41 - d23) / q

                # if abs(inv7)>1.0:
                #     print("debug value inv7=", inv7)
                strikeang = .5 * np.arctan2(d12 - d34, d13 + d24) * (180 / np.pi)
                strikeangerr = abs(.5 * np.arcsin(inv7)) * (180 / np.pi)

                self.inv1[ii] = inv1
                self.inv2[ii] = inv2
                self.inv3[ii] = inv3
                self.inv4[ii] = inv4
                self.inv5[ii] = inv5
                self.inv6[ii] = inv6
                self.inv7[ii] = inv7
                self.q[ii] = q
                self.strike[ii] = strikeang
                self.strike_err[ii] = strikeangerr

    def rotate(self, rot_z):
        """
        Rotates the impedance tensor by the angle rot_z clockwise positive
        assuming 0 is North
        
        """

        self.rot_z = rot_z
        # rotate the data
        # self._Z.rotate(self.rot_z)

        # --> update the invariants (though they should be rotationally
        # invariant except for the strike angle)
        self.compute_invariants()

    def set_z(self, z_array):
        """
        set the z array.

        If the shape changes or the freq are changed need to input 
        those as well.        
        """

        self.z = z_array

        # --> update the invariants
        self.compute_invariants()

    def set_z_err(self, z_err_array):
        """
        set the z_err array.

        If the shape changes or the freq are changed need to input 
        those as well.        
        """

        self.z_err = z_err_array

        # --> update the invariants
        self.compute_invariants()

    def set_freq(self, freq):
        """
        set the freq array, needs to be the same length at z
        """
        
        self.freq = freq

    def __str__(self):
        return "Computes the invariants of the impedance tensor according " + \
               "Weaver et al., [2000, 2003]."
