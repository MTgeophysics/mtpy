#!/usr/bin/env python

"""
mtpy/mtpy/analysis/pt.py

Contains classes and functions for handling Phase Tensor analysis of given impedance tensors (Z). 
 
    Class:
    "PhaseTensor" contains information about a Phase tensor PT. 

        Methods:

        - set_pt
        - set_pterr
        - set_frequencies
        - read_edi_file
        - read_edi
        - read_z
        - read_z_array
        - invariants
        - trace
        - alpha
        - beta
        - skew
        - det
        - _pi1
        - _pi2
        - phimin
        - phimax
        - rotate
        - only1d
        - only2d

    Class:
    "ResidualPhaseTensor" contains information about a REsidual Phase tensor ResPT.

        Methods:

        - read_pt_objects
        - read_pts
        - set_rpt
        - set_rpterr


    Functions:

    - z2pt
    - z_object2pt
    - edi_object2pt
    - edi_file2pt


@UofA, 2013
(LK)

"""

#=================================================================
import numpy as np
import copy

import mtpy.core.edi as MTedi 
import mtpy.core.z as MTz 
import mtpy.utils.exceptions as MTex
import mtpy.utils.calculator as MTcc

#reload(MTex)
#reload(MTedi)
#reload(MTz)
#reload(MTcc)


#=================================================================

class PhaseTensor(object):
    """
        PhaseTensor class - generates a Phase Tensor (PT) object.

        Methods  include reading of Z(arrays and instances), rotations/combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        
        PT is a complex array of the form (n_frequencies, 2, 2), 
        with indices in the following order: 
            PTxx: (0,0) - PTxy: (0,1) - PTyx: (1,0) - PTyy: (1,1)   

        All internal methods are based on (Caldwell et al.,2004) and (Bibby et al.,2005), in which they use the canonical cartesian 2D reference (x1, x2). However, all components, coordinates, and angles for in- and outputs are given in the geographical reference frame:
            x-axis = North ; y-axis = East (; z-axis = Down) 
            
            Therefore, all results from using those methods are consistent (angles are referenced from North rather than x1).

    """

    def __init__(self, pt_array = None, pterr_array = None, z_array = None, zerr_array = None, z_object = None):
        """
            Initialise an instance of the PhaseTensor class.

            Optional input:
            pt_array : Numpy array containing Phase-Tensor values
            pterr_array : Numpy array containing Phase-Tensor-error values
            z_array : Numpy array containing Z values
            zerr_array : Numpy array containing Z-error values (NOT variance, but stddev!)
            z_object: MTpy core.z Z class instance

            (Initialise attributes with None)
        """
    
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
            self.pt = None
            pass

        try:
            if len(pt_array.shape) == 2 and pt_array.shape == (2,2):
                if pt_array.dtype in ['float']:
                    self.pt = np.zeros((1,2,2))
                    self.pt[0] = pt_array
        except:
            self.pt = None
            pass

        try:
            if len(pterr_array.shape) == 3 and ptrr_array.shape[1:3] == (2,2):
                if pterr_array.dtype in ['float']:
                    if self.pt is not None:
                        if self.pt.shape == pterr_array.shape:
                            for idx_f in range(pterr_array.shape[0]):
                                self.pterr[idx_f] = pterr_array[idx_f]
        except:
            self.pterr = None  
            pass

        try:
            if len(pterr_array.shape) == 2 and pterr_array.shape == (2,2):
                if pterr_array.dtype in ['float']:
                    if self.pt is not None:
                        if self.pt.shape[-2:] == pterr_array.shape:
                            self.pterr = np.zeros((1,2,2))
                            self.pterr[0] =  pt_array           
        except:
            self.pterr = None  
            pass

 
        #B) check, if z array is given in some form as well

        # check, if valid Z object is given 
        if isinstance(z_object, MTz.Z):
            
            try:
                z_array = z_object.z

                try:
                    zerr_array = z_object.zerr
                    if z_array.shape != zerr_array.shape:
                        zerr_array = None
                except:
                    pass
            except:
                z_array = None
                zerr_array = None
                pass

        


        #3. if provided PT array was invalid, try to use the Z array 
        if self.pt is None:
            #nulling the error array:
            self.pterr = None
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
                            self.pterr = None                        
            except:
                self.pt = None
                self.pterr = None                        
                

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
           
        self._rotation_angle = None

        if self.pt is not None:
            self._rotation_angle = np.zeros(len(self.pt))

        self._frequencies = None



    def set_pt(self, pt_array):
        """
            Set the attribute 'pt'.

            Input:
            Phase-Tensor array

            Test for shape, but no test for consistency!

        """         

        if not len(pt_array.shape) in [2,3]:
            raise MTex.MTpyError_PT('ERROR - I cannot set new pt array! Invalid dimensions')
        if not pt_array.shape[-2:] == (2,2):
            raise MTex.MTpyError_PT('ERROR - I cannot set new pt array! Invalid dimensions')
        try:
            if not pt_array.dtype in ['float']:
                raise
        except:
            raise MTex.MTpyError_PT('ERROR - I cannot set new pt array! Invalid data type (float expected)')

        if len(pt_array.shape) == 3:
            self.pt = pt_array
        else:
            self.pt = np.zeros((1,pt_array.shape[0],pt_array.shape[1])) 
            self.pt[0] = pt_array

        #testing existing atributes for consistent shapes:
        try:
            if shape(self.pt) != shape(self.pterr):
                raise
        except:
            print 'Shape of new PT array and existing PTerror do not match - setting PTerror to "None"'
            self.pterr = None
        try:
            if len(self.pt) != len(self.freq):
                raise
        except:
            print 'Shape of new PT array and existing "Frequencies" do not match - setting frequencies to "None"'
            self.freq = None
        try:
            if len(self.pt) != len(self.rotation_angle):
                raise
        except:
            print 'Shape of new PT array and existing "Rotation_angle" do not match - setting rotation_angle to "None"'
            self.rotation_angle = None


    def set_pterr(self, pterr_array):
        """
            Set the attribute 'pterr'.

            Input:
            Phase-Tensor-error array

            Test for shape, but no test for consistency!

        """         

        if not len(pterr_array.shape) in [2,3]:
            raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! Invalid dimensions')
        if not pterr_array.shape[-2:] == (2,2):
            raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! Invalid dimensions')
        try:
            if not pterr_array.dtype in ['float']:
                raise
        except:
            raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! Invalid data type (float expected)')

        if self.pt is not None:
            if self.pt.shape != pterr_array.shape:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! Invalid dimensions')
 

        if len(pterr_array.shape) == 3:
            self.pterr = pterr_array
        else:
            self.pterr = np.zeros((1,pt_array.shape[0],pt_array.shape[1])) 
            self.pterr[0] = pterr_array



    def _set_frequencies(self, lo_frequencies):
        """
            Set array of frequencies.

            Input:
            list/array of frequencies (iterable)

            No test for consistency!
        """ 

        if (self.pt is not None):
            if lo_frequencies is not None:
                if (len(lo_frequencies) is not len(self.pt)):
                    print 'length of frequency list not correct (%i instead of %i)'%(len(lo_frequencies), len(self.pt))
                    return
        try:
            self._frequencies = np.array(lo_frequencies)
        except:
            self._frequencies = None

    def _get_frequencies(self): return self._frequencies
    
    freq = property(_get_frequencies, _set_frequencies,doc="")



    def _read_edi_file(self,fn):
        """
            Read in EDI file and convert information into PhaseTensor object attributes.
        """

        e = MTedi.Edi()
        self.read_edi(e)


    def _read_edi(self,edi_object):
        """
            Read in EDI object and convert information into PhaseTensor object attributes.
        """

        z = edi_object.z
        zerr = edi_object.zerr
        self.pt = np.zeros((len(z),2,2))
 
        if zerr is not None:
            self.pterr = np.zeros_like(self.pt) 

            for idx_f in range(len(z)):
                self.pt[idx_f], self.pterr[idx_f] = z2pt(z[idx_f], zerr[idx_f] )

        else:
            for idx_f in range(len(z)):
                self.pt[idx_f] = z2pt( z[idx_f])[0]


        pt.freq = edi_object.freq
        pt.rotation_angle = edi_object.zrot


    def read_z(self,z_object):
        """
            Read in Z object and convert information into PhaseTensor object attributes.
        """
        
        z = z_object.z
        zerr = z_object.zerr
        self.pt = np.zeros((len(z),2,2))
 
        if zerr is not None:
            self.pterr = np.zeros_like(self.pt) 

            for idx_f in range(len(z)):
                self.pt[idx_f], self.pterr[idx_f] = z2pt(z[idx_f], zerr[idx_f] )

        else:
            for idx_f in range(len(z)):
                self.pt[idx_f] = z2pt( z[idx_f])[0]

        self.freq = z_object.freq
        self.rotation_angle = z_object.rotation_angle


    def read_z_array(self,z_array, zerr_array = None):
        """
            Read in Z array (optional Z-error) and convert information into PhaseTensor object attributes.
        """

        z = z_array
        zerr = zerr_array 

        self.pt = np.zeros((len(z),2,2))
 
        if zerr is not None:
            self.pterr = np.zeros_like(self.pt) 

            for idx_f in range(len(z)):
                self.pt[idx_f], self.pterr[idx_f] = z2pt(z[idx_f], zerr[idx_f] )

        else:
            for idx_f in range(len(z)):
                self.pt[idx_f] = z2pt( z[idx_f])[0]

        if self.freq is not None:
            if len(self.freq) != len(z_array):
                self.freq = None


    def _get_invariants(self):
        """
            Return a dictionary of PT-invariants.

            Contains:
            trace, skew, det, phimax, phimin, beta
        """
        if self.pt is None:
            return None

        inv_dict = {}
        inv_dict['trace'] = self.trace[0] 
        inv_dict['skew'] = self.skew[0]
        inv_dict['det'] = self.det[0] 

        inv_dict['phimax'] = self.phimax[0] 
        inv_dict['phimin'] = self.phimin[0] 
        inv_dict['beta'] = self.beta[0] 

        return inv_dict

    invariants = property(_get_invariants, doc="")

    def _get_trace(self):
        """
            Return the trace of PT (incl. uncertainties).

            Output:
            - Trace(PT) - Numpy array
            - Error of Trace(PT) - Numpy array

        """
        if self.pt is None:
            return None, None
        
        tr = np.array( [np.trace(i) for i in self.pt])

        tr_err = None
        if self.pterr is not None:
            tr_err = np.zeros_like(tr)
            tr_err[:] = self.pterr[:,0,0] + self.pterr[:,1,1]


        return tr, tr_err

    trace = property(_get_trace, doc= "")


    def _get_alpha(self):
        """
            Return the principal axis angle (strike) of PT in degrees (incl. uncertainties).

            Output:
            - Alpha - Numpy array
            - Error of Alpha - Numpy array

        """

        if self.pt is None:
            return None, None

        alpha = np.degrees(0.5 * np.arctan2( self.pt[:,0,1] + self.pt[:,1,0]  , self.pt[:,0,0] - self.pt[:,1,1] )  )
        
        alphaerr = None
        if self.pterr is not None:
            alphaerr = np.zeros_like(alpha)
            y = self.pt[:,0,1] + self.pt[:,1,0]
            yerr = np.sqrt( self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2  )
            x = self.pt[:,0,0] - self.pt[:,1,1] 
            xerr = np.sqrt( self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2  )

            alphaerr[:] = 0.5 / ( x**2 + y**2) * np.sqrt( y**2 * xerr**2 + x**2 * yerr**2 )

        return alpha, alphaerr
        
    alpha = property(_get_alpha, doc = "")

    def _get_beta(self):
        """
            Return the 3D-dimensionality angle Beta of PT in degrees (incl. uncertainties).

            Output:
            - Beta - Numpy array
            - Error of Beta - Numpy array

        """

        if self.pt is None:
            return None, None
        
        beta = np.degrees(0.5 * np.arctan2( self.skew[0], self.trace[0])  )
        betaerr = None

        if self.pterr is not None:
            betaerr = np.zeros_like(beta)

            y = self.skew[0]
            yerr = self.skew[1]
            x = self.trace[0]
            xerr = self.trace[1]

            betaerr[:] = 0.5 / ( x**2 + y**2) * np.sqrt( y**2 * xerr**2 + x**2 * yerr**2 )

        return beta, betaerr

    beta = property(_get_beta, doc = "")

    def _get_skew(self):
        """
            Return the skew of PT (incl. uncertainties).

            Output:
            - Skew(PT) - Numpy array
            - Error of Skew(PT) - Numpy array

        """
        if self.pt is None:
            return None, None
       
        skew =  np.array( [ i[0,1] - i[1,0] for i in self.pt ] )
        
        skewerr = None
        if self.pterr is not None:
            skewerr = np.zeros_like(skew)
            skewerr[:] = self.pterr[:,0,1] + self.pterr[:,1,0]

        return skew, skewerr

    skew = property(_get_skew, doc = "")

    def _get_det(self):
        """
            Return the determinant of PT (incl. uncertainties).

            Output:
            - Det(PT) - Numpy array
            - Error of Det(PT) - Numpy array

        """
        if self.pt is None:
            return None, None

        det_phi = np.array( [np.linalg.det(i) for i in self.pt])
        
        det_phi_err = None
        if self.pterr is not None:
            det_phi_err = np.zeros_like(det_phi)
            det_phi_err[:] = np.abs(self.pt[:,1,1] * self.pterr[:,0,0]) + np.abs(self.pt[:,0,0] * self.pterr[:,1,1]) + np.abs(self.pt[:,0,1] * self.pterr[:,1,0]) + np.abs(self.pt[:,1,0] * self.pterr[:,0,1])

        return det_phi, det_phi_err

    det = property(_get_det, doc = "")

    def _pi1(self):
        """
            Return Pi1 (incl. uncertainties).
            
            Pi1 is calculated according to Bibby et al. 2005: Pi1 = 0.5 * sqrt( PT[0,0] - PT[1,1] )**2 + (PT[0,1] + PT[1,0] )**2 )

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """
        #after bibby et al. 2005

        pi1 = 0.5 * np.sqrt( (self.pt[:,0,0] - self.pt[:,1,1] )**2 + (self.pt[:,0,1] + self.pt[:,1,0] )**2 )
        pi1err = None

        if self.pterr is not None:
            pi1err = 1./ pi1 * np.sqrt( (self.pt[:,0,0] - self.pt[:,1,1] )**2 * (self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2)  +\
                                       (self.pt[:,0,1] + self.pt[:,1,0] )**2 * (self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2) )

        return pi1, pi1err

    def _pi2(self):
        """
            Return Pi1 (incl. uncertainties).
            
            Pi1 is calculated according to Bibby et al. 2005: Pi1 = 0.5 * sqrt( PT[0,0] + PT[1,1] )**2 + (PT[0,1] - PT[1,0] )**2 )

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """
        #after bibby et al. 2005

        pi2 = 0.5 * np.sqrt( (self.pt[:,0,0] + self.pt[:,1,1] )**2 + (self.pt[:,0,1] - self.pt[:,1,0] )**2 )
        pi2err = None

        if self.pterr is not None:
            pi2err = 1./ pi2 * np.sqrt( (self.pt[:,0,0] + self.pt[:,1,1] )**2 * (self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2)  +\
                                       (self.pt[:,0,1] - self.pt[:,1,0] )**2 * (self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2) )


        return pi2, pi2err
   


    def _get_phimin(self):
        """
            Return the angle Phi_min of PT (incl. uncertainties).
            
            Phi_min is calculated according to Bibby et al. 2005: Phi_min = Pi2 - Pi1 

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """


        #following caldwell et al 2004:

        #det = np.array( [np.linalg.det(i) for i in self.pt])

        #phimin = np.zeros_like(det)

        # for i in range(len(self.pt)):
        #     s = 1.
        #     if det[i] < 0 :
        #         s = -1.
       
        #     phimin[i] = s * (0.5 * np.sqrt( self.trace()[0][i]**2 + self.skew()[0][i]**2 ) - np.sqrt( 0.25* self.trace()[0][i]**2 + 0.25*self.skew()[0][i]**2 - np.abs(  det[i] )) )

        #following bibby et al 2005:

        if self.pt is None:
            return None, None

        phimin = self._pi2()[0] - self._pi1()[0]

        phiminerr = None
        if self.pterr is not None:
            phiminerr = np.sqrt( self._pi2()[1]**2 + self._pi1()[1]**2 )
 
        return np.degrees(np.arctan(phimin)), np.degrees(np.arctan(phiminerr))

    phimin = property(_get_phimin, doc ="")

    def _get_phimax(self):
        """
            Return the angle Phi_max of PT (incl. uncertainties).
            
            Phi_max is calculated according to Bibby et al. 2005: Phi_max = Pi2 + Pi1 

            Output:
            - Phi_max - Numpy array
            - Error of Phi_max - Numpy array

        """
        
        #following caldwell et al 2004:

        # det = np.array( [np.linalg.det(i) for i in self.pt])

        # phimax = np.zeros_like(det)

        # for i in range(len(phimax)):
        #     phimax[i] = 0.5 * np.sqrt( self.trace()[0][i]**2 + self.skew()[0][i]**2 ) +  np.sqrt( 0.25* self.trace()[0][i]**2 + 0.25*self.skew()[0][i]**2 - np.abs(  det[i] ))

        #following bibby et al 2005:

        if self.pt is None:
            return None, None

        phimax = self._pi2()[0] + self._pi1()[0]

        phimaxerr = None
        if self.pterr is not None:
            phimaxerr = np.sqrt( self._pi2()[1]**2 + self._pi1()[1]**2 )
 
        return np.degrees(np.arctan(phimax)), np.degrees(np.arctan(phimaxerr))

    phimax = property(_get_phimax, doc = "")

    def rotate(self,alpha):
        """
            Rotate PT array. Change the rotation angles attribute respectively.

            Rotation angle must be given in degrees. All angles are referenced to geographic North, positive in clockwise direction. (Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.


        """


        if self.pt is None :
            print 'pt-array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length 1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha%360)
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.pt]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha%360)
                except:
                    print '"Angle" must be a valid number (in degrees)'
                    return
                #make an n long list of identical angles
                lo_angles = [degreeangle for i in self.pt]
            else:                    
                try:
                    lo_angles = [ float(i%360) for i in alpha]
                except:
                    print '"Angles" must be valid numbers (in degrees)'
                    return
            
        self.rotation_angle = list( (np.array(lo_angles) + np.array(self.rotation_angle))%360)

        if len(lo_angles) != len(self.pt):
            print 'Wrong number Number of "angles" - need %i '%(len(self.pt))
            self.rotation_angle = 0.
            return
        
        pt_rot = copy.copy(self.pt)
        pterr_rot = copy.copy(self.pterr)
       
        for idx_freq in range(len(self.pt)):
                    
            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.pterr is not None:
                pt_rot[idx_freq], pterr_rot[idx_freq] = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle, self.pterr[idx_freq,:,:])
            else:
                pt_rot[idx_freq], pterr_rot = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle)

        
        self.pt = pt_rot
        self.pterr = pterr_rot


    def _get_only1d(self):
        """
            Return PT in 1D form.

            If PT is not 1D per se, the diagonal elements are set to zero, the off-diagonal elements keep their signs, but their absolute is set to the mean of the original PT off-diagonal absolutes.
        """

        if self.pt is None:
            return None

        pt1d = copy.copy(self.pt)

        for i in range(len(pt1d)):
            pt1d[i,0,1] = 0
            pt1d[i,1,0] = 0
            
            mean1d = 0.5* (pt1d[i,0,0]+pt1d[i,1,1])
            pt1d[i,0,0] = mean1d
            pt1d[i,1,1] = mean1d

        return pt1d

    only1d = property(_get_only1d, doc = "")

    def _get_only2d(self):
        """
            Return PT in 2D form.

            If PT is not 2D per se, the diagonal elements are set to zero.
        """
        if self.pt is None:
            return None

        pt2d = copy.copy(self.pt)

        for i in range(len(pt2d)):
            pt2d[i,0,1] = 0
            pt2d[i,1,0] = 0
            
            pt2d[i,0,0] = self.phimax[0][i]
            pt2d[i,1,1] = self.phimin[0][i]
            
        return pt2d

    only2d = property(_get_only2d, doc = "")


class ResidualPhaseTensor(PhaseTensor):
    """
        PhaseTensor class - generates a Phase Tensor (PT) object DeltaPhi

        DeltaPhi = 1 - (Phi2.I*Phi1 + Phi*Phi2.I)/2

    """

    def __init__(self, pt_object1 = None, pt_object2 = None):
        """
            Initialise an instance of the ResidualPhaseTensor class.

            Optional input:
            pt_object1 : instance of the PhaseTensor class
            pt_object2 : instance of the PhaseTensor class

            Initialise the attributes with None
        """

        self.rpt = None
        self.rpterr = None
        self._pt1 = None  
        self._pt2 = None  
        self._pt1err = None  
        self._pt2err = None  

        if pt_object1 is not None or  pt_object2 is not None:
            if not (( isinstance(pt_object1,PhaseTensor) and isinstance(pt_object2,PhaseTensor))):
                raise MTex.MTpyError_PT('ERROR - arguments must be instances of the PhaseTensor class')
            
            self.read_pt_objects(pt_object1,pt_object2)


    def read_pt_objects(self, pt_o1, pt_o2):
        """
            Read in two instance of the MTpy PhaseTensor class.

            Update attributes:
            rpt, rpterr, _pt1, _pt2, _pt1err, _pt2err  

        """

        if not ( (isinstance(pt_o1, PhaseTensor)) and (isinstance(pt_o2, PhaseTensor)) ):
            raise MTex.MTpyError_PT('ERROR - both arguments must be instances of the PhaseTensor class')

        pt1 = pt_object1.pt
        pt2 = pt_object2.pt

        if pt1 is not None and pt2 is not None:
            try:
                if (not np.dtype(pt1) in ['float']) or (not np.dtype(pt2) in ['float']):
                    raise
                if not pt1.shape == pt2.shape:
                    raise
                if (not len(pt1.shape) in [2,3] ) :
                    raise

                if len(pt1.shape) == 3:
                    self.rpt = np.zeros((len(pt1),2,2))

                    for idx in range(len(pt1)):
                        self.rpt[idx] = np.eye(2) - 0.5 * np.array( np.dot( np.matrix(pt2[idx]).I, np.matrix(pt1[idx]) ) + np.dot( np.matrix(pt1[idx]), np.matrix(pt2[idx]).I ) ) 
                    self._pt1 = pt1  
                    self._pt2 = pt2  

                else:
                    self.rpt = np.zeros((1,2,2))
                    self.rpt[0] = np.eye(2) - 0.5 * np.array( np.dot( np.matrix(pt2).I, np.matrix(pt1) ) + np.dot( np.matrix(pt1), np.matrix(pt2).I ) ) 
                    
                    self._pt1 =  np.zeros((1,2,2))  
                    self._pt1[0] = pt1 
                    self._pt2 =  np.zeros((1,2,2))  
                    self._pt2[0] = pt2 

            except:
                raise MTex.MTpyError_PT('ERROR - both PhaseTensor objects must contain PT arrays of the same shape')

        else:
            print  'Could not determine ResPT - both PhaseTensor objects must contain PT arrays of the same shape'


        pt1err = pt_object1.pterr
        pt2err = pt_object2.pterr

        if pt1err is not None and pt2err is not None:
            try:
                if (not np.dtype(pt1err) in ['float']) or (not np.dtype(pt2err) in ['float']):
                    raise
                if not pt1err.shape == pt2err.shape:
                    raise
                if (not len(pt1err.shape) in [2,3] ):
                    raise

                if self.rpterr.shape != self.rpt.shape:
                    raise

                if len(pt1err.shape) == 3:
                    self.rpt = np.zeros((len(pt1),2,2))

                    for idx in range(len(pt1err)):
                        matrix1 = pt1[idx]
                        matrix1err = pt1err[idx]                        
                        matrix2, matrix2err = invertmatrix_incl_errors(pt2[idx], inmatrix_err = pt2err[idx])

                        summand1,err1 = multiplymatrices_incl_errors(matrix2, matrix1, inmatrix1_err = matrix2err,inmatrix2_err =  matrix1err)
                        summand2,err2 = multiplymatrices_incl_errors(matrix1, matrix2, inmatrix1_err = matrix1err,inmatrix2_err =  matrix2err)

                        self.rpterr[idx] = np.sqrt( 0.25 * err1**2 + 0.25 * err2**2 )

                    self._pterr1 = pt1err  
                    self._pterr2 = pt2err  

                else:
                    self.rpt = np.zeros((1,2,2))
                    self.rpt[0] = np.eye(2) - 0.5 * np.array( np.dot( np.matrix(pt2).I, np.matrix(pt1) ) + np.dot( np.matrix(pt1), np.matrix(pt2).I ) ) 
            
                    self._pt1err =  np.zeros((1,2,2))  
                    self._pt1err[0] = pt1err
                    self._pt2err =  np.zeros((1,2,2))  
                    self._pt2err[0] = pt2err 

            except:
                raise MTex.MTpyError_PT('ERROR - both PhaseTensor objects must contain PT-error arrays of the same shape')
                
        else:
            print  'Could not determine ResPT uncertainties - both PhaseTensor objects must contain PT-error arrays of the same shape'
 

    def read_pts(self, pt1, pt2, pt1err = None, pt2err = None):
        """
            Read two PT arrays and calculate the ResPT array (incl. uncertainties).


            Input:
            - 2x PT array

            Optional:
            - 2x PTerror array
            
        """

        try:
            if pt1.shape != pt2.shape:
                raise

        except:
            raise MTex.MTpyError_PT('ERROR - could not build ResPT array from given PT arrays - check shapes! ')
        #TODO - check arrays here:


        pt_o1 = PhaseTensor(pt_array = pt1, pterr_array = pt1err)
        pt_o2 = PhaseTensor(pt_array = pt2, pterr_array = pt2err)

        self.read_pt_objects(pt_o1,pt_o2)
        

    def set_rpt(self, rpt_array):
        """
            Set the attribute 'rpt'.

            Input:
            ResPT array

            Test for shape, but no test for consistency!

        """ 
        if (self.rpt is not None) and (self.rpt.shape != rpt_array.shape):
            print 'Error - shape of "ResPT" array does not match shape of existing rpt array: %s ; %s'%(str(rpt_array.shape),str(self.rpt.shape))
            return

        self.rpt = rpt_array

    def set_rpterr(self, rpterr_array):
        """
            Set the attribute 'rpterr'.

            Input:
            ResPT-error array

            Test for shape, but no test for consistency!

        """ 
        if (self.rpterr is not None) and (self.rpterr.shape != rpterr_array.shape):
            print 'Error - shape of "ResPT-error" array does not match shape of existing rpterr array: %s ; %s'%(str(rpterr_array.shape),str(self.rpterr.shape))
            return

        self.rpterr = rpterr_array


#=======================================================================

def z2pt(z_array, zerr_array = None):
    """
        Calculate Phase Tensor from Z array (incl. uncertainties)

        Input:
        - Z : 2x2 complex valued Numpy array

        Optional:
        - Z-error : 2x2 real valued Numpy array

        Return:
        - PT : 2x2 real valued Numpy array
        - PT-error : 2x2 real valued Numpy array

    """
    if z_array is not None:
        try:
            if not  len(z_array.shape) in [2,3]:
                raise
            if not z_array.shape[-2:] == (2,2):
                raise
            if not z_array.dtype in ['complex', 'float']:
                raise
        except:
            raise MTex.MTpyError_PT('Error - incorrect z array: %s;%s instead of (N,2,2);complex'%(str(z_array.shape), str(z_array.dtype)))    


    if zerr_array is not None:
        try:
            if not  len(zerr_array.shape) in [2,3]:
                raise
            if not zerr_array.shape[-2:] == (2,2):
                raise
            if not zerr_array.dtype in ['float']:
                raise
        except:
            raise MTex.MTpyError_PT('Error - incorrect z-err-array: %s;%s instead of (N,2,2);real'%(str(zerr_array.shape), str(zerr_array.dtype)))

        if not z_array.shape == zerr_array.shape:
            raise MTex.MTpyError_PT('Error - z-array and z-err-array have different shape: %s;%s'%(str(z_array.shape), str(zerr_array.shape)))




    #for a single matrix as input:
    if len(z_array.shape) == 2:
        pt_array = np.zeros((2,2))

        realz = np.real(z_array)
        imagz = np.imag(z_array)
        detreal = np.linalg.det(realz)

        if detreal == 0 :
            if np.linalg.norm(realz) == 0 and np.linalg.norm(imagz) == 0:
                pterr_array = np.zeros_like(pt_array)
                if zerr_array is None:
                    pterr_array = None                
                return pt_array, pterr_array

            else:
                raise MTex.MTpyError_PT('Error - z-array contains a singular matrix, thus it cannot be converted into a PT!' )



        pt_array[0,0] =  realz[1,1] * imagz[0,0] - realz[0,1] * imagz[1,0] 
        pt_array[0,1] =  realz[1,1] * imagz[0,1] - realz[0,1] * imagz[1,1] 
        pt_array[1,0] =  realz[0,0] * imagz[1,0] - realz[1,0] * imagz[0,0] 
        pt_array[1,1] =  realz[0,0] * imagz[1,1] - realz[1,0] * imagz[0,1] 

        pt_array /= detreal

        if zerr_array is None:
            return pt_array, None

        pterr_array = np.zeros_like(pt_array)

        #Z entries are independent -> use Gaussian error propagation (squared sums/2-norm)

        pterr_array[0,0] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs( -pt_array[0,0] * realz[1,1] * zerr_array[0,0])**2,
                                                                np.abs(  pt_array[0,0] * realz[0,1] * zerr_array[1,0])**2,
                                                                np.abs(  ( (imagz[0,0] * realz[1,0] - realz[0,0] * imagz[1,0]) / np.abs(detreal) * realz[0,0] ) * zerr_array[0,1])**2, 
                                                                np.abs(  ( (imagz[1,0] * realz[0,0] - realz[1,0] * imagz[1,1]) / np.abs(detreal) * realz[0,1] ) * zerr_array[1,1])**2,
                                                                np.abs(  realz[1,1] * zerr_array[0,0])**2,
                                                                np.abs( realz[0,1] * zerr_array[1,0])**2 ]))


        pterr_array[0,1] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs( -pt_array[0,1] * realz[1,1] * zerr_array[0,0])**2,
                                                                np.abs(  pt_array[0,1] * realz[0,1] * zerr_array[1,0])**2,
                                                                np.abs(  ( (imagz[0,1] * realz[1,0] - realz[0,0] * imagz[1,1]) / np.abs(detreal) * realz[1,1] ) * zerr_array[0,1])**2, 
                                                                np.abs(  ( (imagz[1,1] * realz[0,0] - realz[0,1] * imagz[1,0]) / np.abs(detreal) * realz[0,1] ) * zerr_array[1,1])**2,
                                                                np.abs(  realz[1,1] * zerr_array[0,1])**2,
                                                                np.abs( realz[0,1] * zerr_array[1,1])**2 ]))

        pterr_array[1,0] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs(  pt_array[1,0] * realz[1,0] * zerr_array[0,1])**2,
                                                                np.abs( -pt_array[1,0] * realz[0,0] * zerr_array[1,1])**2,
                                                                np.abs(  ( (imagz[0,0] * realz[1,1] - realz[0,1] * imagz[1,1]) / np.abs(detreal) * realz[1,0] ) * zerr_array[0,0])**2, 
                                                                np.abs(  ( (imagz[1,0] * realz[0,1] - realz[1,1] * imagz[0,0]) / np.abs(detreal) * realz[0,0] ) * zerr_array[0,1])**2,
                                                                np.abs(  realz[1,0] * zerr_array[0,0])**2,
                                                                np.abs( realz[0,0] * zerr_array[1,0])**2 ]))


        pterr_array[1,1] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs(  pt_array[1,1] * realz[1,0] * zerr_array[0,1])**2,
                                                                np.abs( -pt_array[1,1] * realz[0,0] * zerr_array[1,1])**2,
                                                                np.abs(  ( (imagz[0,1] * realz[1,1] - realz[0,1] * imagz[1,1]) / np.abs(detreal) * realz[1,0] ) * zerr_array[0,0])**2, 
                                                                np.abs(  ( (imagz[1,1] * realz[0,1] - realz[1,1] * imagz[0,1]) / np.abs(detreal) * realz[0,0] ) * zerr_array[0,1])**2,
                                                                np.abs( - realz[1,0] * zerr_array[0,1])**2,
                                                                np.abs( realz[0,0] * zerr_array[1,1])**2 ]))


        return pt_array, pterr_array


    #else:
    pt_array = np.zeros((z_array.shape[0],2,2))

    for idx_f in range(len(z_array)):       
        
        realz = np.real(z_array[idx_f])
        imagz = np.imag(z_array[idx_f])

        detreal = np.linalg.det(realz)
        if detreal == 0 :
            raise MTex.MTpyError_PT('Error - z-array contains a singular matrix, thus it cannot be converted into a PT!' )

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
    """
        Calculate Phase Tensor from Z object (incl. uncertainties)

        Input:
        - Z-object : instance of the MTpy Z class


        Return:
        - PT : nx2x2 real valued Numpy array
        - PT-error : nx2x2 real valued Numpy array

    """

    if not isinstance(z_object, MTz.Z):
        raise MTex.MTpyError_Z('Input argument is not an instance of the Z class')

    p = PhaseTensor(z_object = z_object)

    pt_array = p.pt
    pterr_array = p.pterr

    return pt_array, pterr_array


def edi_object2pt(edi_object):
    """
        Calculate Phase Tensor from Edi object (incl. uncertainties)

        Input:
        - Edi-object : instance of the MTpy Edi class


        Return:
        - PT : nx2x2 real valued Numpy array
        - PT-error : nx2x2 real valued Numpy array

    """


    if not isinstance(z_object, MTedi.Edi):
        raise MTex.MTpyError_EDI('Input argument is not an instance of the Edi class')
    p = PhaseTensor(edi_object = edi_object)

    pt_array = p.pt
    
    pterr_array = p.pterr

    return pt_array, pterr_array


def edi_file2pt(filename):
    """
        Calculate Phase Tensor from Edi-file (incl. uncertainties)

        Input:
        - Edi-file : full path to the Edi-file


        Return:
        - PT : nx2x2 real valued Numpy array
        - PT-error : nx2x2 real valued Numpy array

    """

    e = MTedi.Edi()
    e.readfile(filename)

    p = PhaseTensor(edi_object = e)

    pt_array = p.pt
    
    pterr_array = p.pterr

    return pt_array, pterr_array



