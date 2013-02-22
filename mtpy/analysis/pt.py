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
import mtpy.core.z as MTz 
import mtpy.utils.format as MTformat
import mtpy.utils.exceptions as MTexceptions

reload(MTexceptions)
reload(MTedi)
reload(MTz)
reload(MTformat)


#=================================================================

class PhaseTensor(object):
    """
        PhaseTensor class - generates a Phase Tensor (PT) object.

        Methods  include reading and writing from and to edi-objects, rotations/combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        
        PT is a complex array of the form (n_frequencies, 2, 2), 
        with indices in the following order: 
            PTxx: (0,0) - PTxy: (0,1) - PTyx: (1,0) - PTyy: (1,1)   

        All internal methods are based on (Caldwell et al.,2004), in which they use the canonical cartesian 2D reference (x1, x2). However, all components, coordinates, and angles for in- and outputs are given in the geographical reference frame:
            x-axis = North ; y-axis = East (; z-axis = Down) 
            Therefore, all results from using those methods are consisting (angles are referenced from North rather than x1)

  

    """

    def __init__(self, pt_array = None, pterr_array = None, z_array = None, zerr_array = None, z_object = None, edi_object = None):
    
        self.pt = None
        self.pterr = None  

        #A) check for direct import of a provided pt array 
        try:
            if len(pt_array.shape) == 3 and pt_array.shape[1:3] == (2,2):
                if pt_array.dtype in ['float']:
                    self.pt = np.zeros_like(pt_array)
                    for idx_f in range(pt_array.shape[0]):
                        self.pt[idx_F] = pt_array[idx_F]

        except:
            pass

        try:
            if len(pt_array.shape) == 2 and pt_array.shape == (2,2):
                if pt_array.dtype in ['float']:
                    self.pt = np.zeros((1,2,2))
                    self.pt[0] = pt_array
        except:
            pass

        try:
            if len(pterr_array.shape) == 3 and ptrr_array.shape[1:3] == (2,2):
                if pterr_array.dtype in ['float']:
                    if self.pt is not None:
                        if self.pt.shape == pterr_array.shape:
                            for idx_f in range(pterr_array.shape[0]):
                                self.pterr[idx_f] = pterr_array[idx_f]
        except:
            pass

        try:
            if len(pterr_array.shape) == 2 and pterr_array.shape == (2,2):
                if pterr_array.dtype in ['float']:
                    if self.pt is not None:
                        if self.pt.shape[-2:] == pterr_array.shape:
                            self.pterr = np.zeros((1,2,2))
                            self.pterr[0] =  pt_array           
        except:
            pass

        #B) check, if z array is given as well
        #1. check, if valid Edi object is given
        if isinstance(edi_object,MTedi.Edi):

            try:
                z_array = edi_object.z
            except:
                pass

            try:
                zerr_array = edi_object.zerr_array
            except:
                pass

        #2. otherwise check, if valid Z object is given 
        elif isinstance(z_object,MTz.Z):
            
            try:
                z_array = z_object.z
            except:
                pass
            try:
                zerr_array = edi_object.zerr_array
            except:
                pass

        if z_array is not None and zerr_array is not None:
            if z_array.shape != zerr_array.shape:
                zerr_array = None

        #3. if provided PT array was invalid, try to use the Z array 
        if self.pt is None:
            try:
                if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                    if z_array.dtype in ['complex', 'float']:
                        try:
                            self.pt = np.zeros((len(z_array),2,2))
                            if zerr_array is not None:
                                    if zerr_array.shape == z_array.shape:
                                        self.pterr = np.zeros_like(self.pt)
                                        for idx_f in range(len(z_array)):
                                            self.pt[idx_f], self.pterr[idx_f] = z2pt(z_array[idx_f], zerr_array[idx_f] )

                            else:
                                for idx_f in range(len(z_array)):
                                    self.pt[idx_f] = z2pt( z_array[idx_f])[0]
                                    
                        except:
                            self.pt = None

            except:
                pass

            try:
                if len(z_array.shape) == 2 and z_array.shape == (2,2):
                    if z_array.dtype in ['complex', 'float']:
                        if zerr_array is not None:
                            self.pt = np.zeros((1,2,2))
                            self.pterr = np.zeros((1,2,2))
                            self.pt[0], self.pterr[0] = z2pt( z_array, zerr_array )    
                        else:
                            self.pt = np.zeros((1,2,2))
                            self.pt[0] = z2pt(z_array)  
            except:
                pass
           
        self.frequencies = None

        self.ptrot = None
        if self.pt is not None:
            self.ptrot = [ 0. for i in self.pt ]


    def set_pt(self, pt_array):

        pass

    def set_pterr(self, pterr_array):
        pass

    def set_frequencies(self,lo_frequencies):
        pass

    def read_edi_file(self,fn):

        pass

    def read_edi(self,edi_object):

        pass

    def read_z(self,z_object):
        pass

    def read_z_array(self,z_array, zerr_array = None):
        pass


    def invariants(self):

        inv_dict = {}
        inv_dict['trace'] =  np.array( [np.trace(i) for i in self.pt])
        inv_dict['skew'] = self.skew()[0]
        inv_dict['det'] = np.array( [np.linalg.det(i) for i in self.pt])

        inv_dict['phimax'] = self.phimax()[0] 
        inv_dict['phimin'] = self.phimin()[0] 
        inv_dict['beta'] = self.beta()[0] 

        return inv_dict

    def alpha(self):

        alpha = 0.5 * np.arctan2( self.pt[:,0,1] +self.pt[:,1,0]  , self.pt[:,0,0] - self.pt[:,1,1] )  
        
        alphaerr = np.zeros_like(alpha)


        return alpha, alphaerr
        

    def beta(self):
        
        beta = 0.5 * np.arctan2( inv_dict['skew'], inv_dict['trace'])  


        return beta, betaerr



    def skew(self):
        
        skew =  np.array( [ i[0,1] - i[1,0] for i in self.pt ] )


        return skew, skewerr

    def phimin(self):

        det = np.array( [np.linalg.det(i) for i in self.pt])

        phimin = np.zeros_like(det)

        for i in range(len(P.pt)):
            s = 1.
            if det[i] < 0 :
                s = -1.

        phinmin[i] = s * (0.5 * np.sqrt( self.trace()[i]**2 + self.skew()[i]**2 ) - np.sqrt( 0.25* self.trace()[i]**2 + 0.25*self.skew()[i]**2 - np.abs(  det[i] )) )

 
        return phimin, phiminerr

    def phimax(self):
        
        det = np.array( [np.linalg.det(i) for i in self.pt])

        phimax = np.zeros_like(det)


        phinmax[i] = 0.5 * np.sqrt( self.trace()[i]**2 + self.skew()[i]**2 ) +  np.sqrt( 0.25* self.trace()[i]**2 + 0.25*self.skew()[i]**2 - np.abs(  det[i] ))

        return phimax, phimaxerr


    def rotate(self,angle):

        pass



def z2pt(z_array, zerr_array = None):
    """
    
    """
    try:
        if not  len(z_array.shape) in [2,3]:
            raise
        if not z_array.shape[-2:] == (2,2):
            raise
        if not z_array.dtype in ['complex', 'float']:
            raise
    except:
        raise MTexceptions.MTpyError_PT('Error - incorrect z array: %s;%s instead of (N,2,2);complex'%(str(z_array.shape), str(z_array.dtype)))
    

    if zerr_array is not None:
        try:
            if not  len(zerr_array.shape) in [2,3]:
                raise
            if not zerr_array.shape[-2:] == (2,2):
                raise
            if not zerr_array.dtype in ['float']:
                raise
        except:
            raise MTexceptions.MTpyError_PT('Error - incorrect z-err-array: %s;%s instead of (N,2,2);real'%(str(zerr_array.shape), str(zerr_array.dtype)))

        if not z_array.shape == zerr_array.shape:
            raise MTexceptions.MTpyError_PT('Error - z-array and z-err-array have different shape: %s;%s'%(str(z_array.shape), str(zerr_array.shape)))
    #for a single matrix as input:
    if len(z_array.shape) == 2:
        pt_array = np.zeros((2,2))

        realz = np.real(z_array)
        imagz = np.imag(z_array)
        detreal = np.linalg.det(realz)
        if detreal == 0 :
            raise MTexceptions.MTpyError_PT('Error - z-array contains a singular matrix, thus it cannot be converted into a PT!' )


        pt_array[0,0] =  realz[1,1] * imagz[0,0] - realz[0,1] * imagz[1,0] 
        pt_array[0,1] =  realz[1,1] * imagz[0,1] - realz[0,1] * imagz[1,1] 
        pt_array[1,0] =  realz[0,0] * imagz[1,0] - realz[1,0] * imagz[0,0] 
        pt_array[1,1] =  realz[0,0] * imagz[1,1] - realz[1,0] * imagz[0,1] 

        pt_array /= detreal

        if zerr_array is None:
            return pt_array, None

        pterr_array = np.zeros_like(pt_array)
        pterr_array[0,0] = 1/detreal * (np.abs( -pt_array[0,0] * realz[1,1] * zerr_array[0,0]) + \
                                        np.abs(  pt_array[0,0] * realz[0,1] * zerr_array[1,0]) + \
                                        np.abs(  (imagz[0,0] - pt_array[0,0] * realz[0,0] ) * zerr_array[1,1]) +\
                                        np.abs(  (-imagz[1,0]+ pt_array[0,0] * realz[1,0] ) * zerr_array[0,1]) + \
                                        np.abs(  realz[1,1] * zerr_array[0,0]) + np.abs( realz[0,1] * zerr_array[1,0]) )

        pterr_array[0,1] = 1/detreal * (np.abs( -pt_array[0,1] * realz[1,1] * zerr_array[0,0]) + \
                                        np.abs(  pt_array[0,1] * realz[0,1] * zerr_array[1,0]) + \
                                        np.abs(  (imagz[0,1] - pt_array[0,1] * realz[0,0] ) * zerr_array[1,1]) +\
                                        np.abs(  (-imagz[1,1]+ pt_array[0,1] * realz[1,0] ) * zerr_array[0,1]) + \
                                        np.abs(  realz[1,1] * zerr_array[0,1]) + np.abs( realz[0,1] * zerr_array[1,1]) )

        pterr_array[1,0] = 1/detreal * (np.abs(  (imagz[1,0] - pt_array[1,0] * realz[1,1] ) * zerr_array[0,0]) +\
                                        np.abs( pt_array[1,0] * realz[1,0] * zerr_array[0,1]) + \
                                        np.abs(  (-imagz[0,0] + pt_array[1,0] * realz[0,1] ) * zerr_array[1,0]) + \
                                        np.abs( -pt_array[1,0] * realz[0,0] * zerr_array[1,1]) + \
                                        np.abs(  realz[0,0] * zerr_array[1,0]) + np.abs( -realz[1,0] * zerr_array[0,0]) )

        pterr_array[1,1] = 1/detreal * (np.abs(  (imagz[1,1] - pt_array[1,1] * realz[1,1] ) * zerr_array[0,0]) +\
                                        np.abs( pt_array[1,1] * realz[1,0] * zerr_array[0,1]) + \
                                        np.abs(  (-imagz[0,1] + pt_array[1,1] * realz[0,1] ) * zerr_array[1,0]) + \
                                        np.abs( -pt_array[1,1] * realz[0,0] * zerr_array[1,1]) + \
                                        np.abs(  realz[0,0] * zerr_array[1,1]) + np.abs( -realz[1,0] * zerr_array[0,1]) )

        return pt_array, pterr_array


    #else:
    pt_array = np.zeros((z_array.shape[0],2,2))

    for idx_f in range(len(z_array)):       
        
        realz = np.real(z_array[idx_f])
        imagz = np.imag(z_array[idx_f])

        detreal = np.linalg.det(realz)
        if detreal == 0 :
            raise MTexceptions.MTpyError_PT('Error - z-array contains a singular matrix, thus it cannot be converted into a PT!' )

        pt_array[idx_f,0,0] =  realz[1,1] * immagz[0,0] - realz[0,1] * immagz[1,0] 
        pt_array[idx_f,0,1] =  realz[1,1] * immagz[0,1] - realz[0,1] * immagz[1,1] 
        pt_array[idx_f,1,0] =  realz[0,0] * immagz[1,0] - realz[1,0] * immagz[0,0] 
        pt_array[idx_f,1,1] =  realz[0,0] * immagz[1,1] - realz[1,0] * immagz[0,1] 

        pt_array /= detreal

        if zerr_array is None:
            return pt_array, pterr_array

        pterr_array = np.zeros_like(pt_array)
        pterr_array[idx_f,0,0] = 1/detreal * (np.abs( -pt_array[idx_f,0,0] * realz[1,1] * zerr_array[0,0]) + \
                                        np.abs(  pt_array[idx_f,0,0] * realz[0,1] * zerr_array[1,0]) + \
                                        np.abs(  (imagz[0,0] - pt_array[idx_f,0,0] * realz[0,0] ) * zerr_array[1,1]) +\
                                        np.abs(  (-imagz[1,0]+ pt_array[idx_f,0,0] * realz[1,0] ) * zerr_array[0,1]) + \
                                        np.abs(  realz[1,1] * zerr_array[0,0]) + np.abs( realz[0,1] * zerr_array[1,0]) )

        pterr_array[idx_f,0,1] = 1/detreal * (np.abs( -pt_array[idx_f,0,1] * realz[1,1] * zerr_array[0,0]) + \
                                        np.abs(  pt_array[idx_f,0,1] * realz[0,1] * zerr_array[1,0]) + \
                                        np.abs(  (imagz[0,1] - pt_array[idx_f,0,1] * realz[0,0] ) * zerr_array[1,1]) +\
                                        np.abs(  (-imagz[1,1]+ pt_array[idx_f,0,1] * realz[1,0] ) * zerr_array[0,1]) + \
                                        np.abs(  realz[1,1] * zerr_array[0,1]) + np.abs( realz[0,1] * zerr_array[1,1]) )

        pterr_array[idx_f,1,0] = 1/detreal * (np.abs(  (imagz[1,0] - pt_array[idx_f,1,0] * realz[1,1] ) * zerr_array[0,0]) +\
                                        np.abs( pt_array[idx_f,1,0] * realz[1,0] * zerr_array[0,1]) + \
                                        np.abs(  (-imagz[0,0] + pt_array[idx_f,1,0] * realz[0,1] ) * zerr_array[1,0]) + \
                                        np.abs( -pt_array[idx_f,1,0] * realz[0,0] * zerr_array[1,1]) + \
                                        np.abs(  realz[0,0] * zerr_array[1,0]) + np.abs( -realz[1,0] * zerr_array[0,0]) )

        pterr_array[idx_f,1,1] = 1/detreal * (np.abs(  (imagz[1,1] - pt_array[idx_f,1,1] * realz[1,1] ) * zerr_array[0,0]) +\
                                        np.abs( pt_array[idx_f,1,1] * realz[1,0] * zerr_array[0,1]) + \
                                        np.abs(  (-imagz[0,1] + pt_array[idx_f,1,1] * realz[0,1] ) * zerr_array[1,0]) + \
                                        np.abs( -pt_array[idx_f,1,1] * realz[0,0] * zerr_array[1,1]) + \
                                        np.abs(  realz[0,0] * zerr_array[1,1]) + np.abs( -realz[1,0] * zerr_array[0,1]) )

    return pt_array, pterr_array



def z_object2pt(z_object):
    pass
    return pt_array, pterr_array


def edi_object2pt(edi_object):
    pass
    return pt_array, pterr_array



