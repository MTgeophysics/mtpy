#!/usr/bin/env python

"""
mtpy/mtpy/analysis/pt.py

Contains classes and functions for handling Phase Tensor analysis of given impedance tensors (Z). 
 
    Class:
    "PhaseTensor" contains information about a Phase tensor PT. 

        Methods:



    Functions:


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

import mtpy.core.z as MTz 
reload (MTz)

import mtpy.utils.format as MTformat
reload(MTformat)

import mtpy.utils.exceptions as MTexceptions
reload(MTexceptions)


#=================================================================

class PhaseTensor(object):
    """
        PhaseTensor class - generates a Phase Tensor (PT) object.

        Methods  include reading and writing from and to edi-objects, rotations/combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        
        PT is a complex array of the form (n_frequencies, 2, 2), 
        with indices in the following order: 
            PTxx: (0,0) - PTxy: (0,1) - PTyx: (1,0) - PTyy: (1,1)   

        All internal methods are based on (Caldwell et al.,2004) , where they use the canonical cartesian 2D reference (x1, x2). However, all components, coordinates, and angles for in- and outputs are given in the geographical reference frame:
            x-axis = North ; y-axis = East (; z-axis = Down) 

    """

    def __init__(self, pt_array = None, pterr_array = None, z_array = None, zerr_array = None, z_object = None, edi_object = None):
    

        self.pt = None
        self.pterr = None  

        #A) check for direct import of a provided pt array 
        try:
            if len(pt_array.shape) == 3 and pt_array.shape[1:3] == (2,2):
                if pt_array.dtype in ['complex', 'float']:
                    self.pt = pt_array
        except:
            pass

        try:
            if len(pt_array.shape) == 2 and pt_array.shape == (2,2):
                if pt_array.dtype in ['complex', 'float']:
                    self.pt = np.zeros((1,2,2),'complex')
                    self.pt[0] = pt_array            
        except:
            pass

        #B) check, if z array is given as well
        #1. check, if valid Edi object is given
        if isinstance(edi_object,MTedi.Edi):

            try:
                z_array = edi_object.z
            except:
                pass

        #2. otherwise check, if valid Z object is given 
        elif isinstance(z_object,MTedi.Z):
            
            try:
                z_array = z_object.z
            except:
                pass

        #3. if provided PT array was invalid, try to use the Z array 
        if self.pt is None:
            try:
                if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                    if z_array.dtype in ['complex', 'float']:
                        try:
                            self.pt = np.zeros((len(z_array),2,2),'complex')
                            for idx_f in range(len(z_array)):
                                self.pt[idx_f] = z2pt(z_array[idx_f])

            except:
                pass

            try:
                if len(z_array.shape) == 2 and z_array.shape == (2,2):
                    if z_array.dtype in ['complex', 'float']:
                        self.pt = np.zeros((1,2,2),'complex')
                        self.pt[0] = z2pt(z_array)       
            except:
                pass
           


        self.frequencies = None
        self.rotation_angle = 0.


    def set_pt(self, pt_array):

        pass

    def set_pterr(self, pterr_array):
        pass


    def read_edi_file(self,fn):

        pass

    def read_edi(self,edi_object):

        pass

    def read_z(self,z_object):
        pass

    def read_z_array(self,z_array, zerr_array = None):
        pass


    def alpha(self):

        pass

    def beta(self):
        pass

    def invariants(self):

        pass

    def skew(self):
        pass

    def phimin(self):

        pass

    def phimax(self):
        pass

    def rotate(self,angle):

        pass






def z2pt(z_array, zerr_array = None):
    pass


    return pt_array, pterr_array

def z_object2pt(z_object):

    return pt_array, pterr_array


def edi_object2pt(edi_object):

    return pt_array, pterr_array
   