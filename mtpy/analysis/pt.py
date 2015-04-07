#!/usr/bin/env python

"""
mtpy/mtpy/analysis/pt.py

Contains classes and functions for handling Phase Tensor analysis of given impedance tensors (Z). 
 
    Class:
    "PhaseTensor" contains information about a Phase tensor PT. 

        Methods:

        - set_pt
        - set_pterr
        - set_freq
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
import mtpy.utils.exceptions as MTex
import mtpy.utils.calculator as MTcc


#=================================================================

class PhaseTensor(object):
    """
        PhaseTensor class - generates a Phase Tensor (PT) object.

        Methods  include reading and writing from and to edi-objects, rotations
        combinations of Z instances, as well as 
        calculation of invariants, inverse, amplitude/phase,...

        
        PT is a complex array of the form (n_freq, 2, 2), 
        with indices in the following order: 
            PTxx: (0,0) - PTxy: (0,1) - PTyx: (1,0) - PTyy: (1,1)   

        All internal methods are based on (Caldwell et al.,2004) and 
		(Bibby et al.,2005), in which they use the canonical cartesian 2D 
            reference (x1, x2). However, all components, coordinates, 
            and angles for in- and outputs are given in the geographical 
            reference frame:
                x-axis = North ; y-axis = East (; z-axis = Down) 
            
            Therefore, all results from using those methods are consistent 
			(angles are referenced from North rather than x1).

    """

    def __init__(self, pt_array = None, pterr_array = None, z_array = None, 
			zerr_array = None, z_object = None, freq=None, pt_rot=0.0):
        """
            Initialise an instance of the PhaseTensor class.

            Optional input:
            pt_array : Numpy array containing Phase-Tensor values
            pterr_array : Numpy array containing Phase-Tensor-error values
            z_array : Numpy array containing Z values
            zerr_array : Numpy array containing Z-error values (NOT variance, but stddev!)
            z_object: MTpy core.z Z class instance
            freq : numpy array containing freq values

            (Initialise attributes with None)
        """

        self._pt = pt_array
        self._pterr = pterr_array
        self._z = z_array
        self._z_err = zerr_array
        self._freq = freq
        self.rotation_angle = pt_rot
        
        #if a z object is input be sure to set the z and z_err so that the
        #pt will be calculated
        # print type(z_object)==type(MTz.Z()),isinstance(z_object, MTz.Z)
        # if isinstance(z_object, MTz.Z):
        if z_object is not None:
            try:
                self.set_z_object(z_object)
            except:
                print '\tWarning - could not digest provided Z-Object'

        elif z_array is not None:

            try:
                self._set_z(z_array)
            except:
                self._z = None
                self._z_err = None
                print 'Can not calculate pt from z==None'

            if zerr_array is not None:
                
                try:
                    self._set_z_err(zerr_array)
                    if z_array.shape != zerr_array.shape:
                        self._set_z_err(None)
                except:
                    pass

        
        if self._freq is None:
            print 'Should input a freq array to know which index of the'+\
                  ' PT array corresponds to which freq.'




    #==========================================================================
    #  define get/set functions and properties
    #==========================================================================
    #---phase tensor array----------------------------------------    
    def _set_pt(self, pt_array):
        """
            Set the attribute 'pt'.

            Input:
            Phase-Tensor array

            Test for shape, but no test for consistency!

        """         
        self._pt = pt_array
        
        #check for dimensions
        if pt_array is not None:
            #--> large array
            if not len(pt_array.shape) in [2,3]:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!'+\
                      ' Invalid dimensions')
            
            #--> single matrix
            if not pt_array.shape[-2:] == (2,2):
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!'+\
                      ' Invalid dimensions')
           
           #--> make sure values are floats
            try:
                if not pt_array.dtype in ['float']:
                    raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!'+\
                          'Invalid dimensions')
            except:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!'+\
                      'Invalid data type (float expected)')
            
            if len(pt_array.shape) == 3:
                self._pt = pt_array
            else:
                self._pt = np.zeros((1,pt_array.shape[0],pt_array.shape[1])) 
                self._pt[0] = pt_array

            #testing existing atributes for consistent shapes:
            try:
                if np.shape(self.pt) != np.shape(self.pterr):
                    raise MTex.MTpyError_inputarguments('pt and pterr are not'+\
                                                        ' the same shape')
            except:
                print 'Shape of new PT array and existing PTerror do not match'+\
                      '- setting PTerror to "None"'
                self._pterr = None
            try:
                if len(self.pt) != len(self.freq):
                    raise MTex.MTpyError_inputarguments('pt and freq are'+\
                                                        'not the same shape')
            except:
                print 'Shape of new PT array and existing "freq" do not'+\
                      'match - setting freq to "None"'
                self._freq = None
            try:
                if len(self.pt) != len(self.rotation_angle):
                    raise MTex.MTpyError_inputarguments('pt and rotation angles'+\
                                                        'are not the same shape')
            except:
                print 'Shape of new PT array and existing "Rotation_angle" do '+\
                       'not match - setting rotation_angle to "None"'
                self.rotation_angle = None
                
        else:
            pass
            
    def _get_pt(self):
        return self._pt
        
    pt = property(_get_pt, _set_pt, doc="Phase tensor array")
    
    #---phase tensor Error-----------------------------------------------------
    def _set_pterr(self, pterr_array):
        """
            Set the attribute 'pterr'.

            Input:
            Phase-Tensor-error array

            Test for shape, but no test for consistency!

        """         
        self._pterr = pterr_array
        
        #check dimensions
        if pterr_array is not None:
            if not len(pterr_array.shape) in [2,3]:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! '+\
                      'Invalid dimensions')
            if not pterr_array.shape[-2:] == (2,2):
                raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! '+\
                      'Invalid dimensions')
            try:
                if not pterr_array.dtype in ['float']:
                    raise ('ERROR - I cannot set new pterr array! '+\
                      'Invalid type')
            except:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pterr array! '+\
                      'Invalid data type (float expected)')
    
            if self.pt is not None:
                if self.pt.shape != pterr_array.shape:
                    raise MTex.MTpyError_PT('ERROR - I cannot set new pterr '+\
                          'array! Invalid dimensions')
     
    
            if len(pterr_array.shape) == 3:
                self.pterr = pterr_array
            else:
                self.pterr = np.zeros((1,self.pt.shape[0],self.pt.shape[1])) 
                self.pterr[0] = pterr_array
        
        else:
            pass
            
    def _get_pterr(self):
        return self._pterr
        
    pterr = property(_get_pterr, _set_pterr, 
                      doc='Phase tensor error array, must be same shape as pt')

    #---freq------------------------------------------------------------
    def _set_freq(self, lo_freq):
        """
            Set array of freq.

            Input:
            list/array of freq (iterable)

            No test for consistency!
        """ 

        if (self._pt is not None):
            if lo_freq is not None:
                if (len(lo_freq) is not len(self._pt)):
                    print 'length of freq list not correct'+\
                          '(%i instead of %i)'%(len(lo_freq), 
                                                    len(self._pt))
                    return
        try:
            self._freq = np.array(lo_freq)
        except:
            self._freq = None

    def _get_freq(self): return self._freq
    
    freq = property(_get_freq, _set_freq,doc="freq array")

    #---z_object---------------------------------------------------------------
    
    def set_z_object(self, z_object):
        """
            Read in Z object and convert information into PhaseTensor object 
            attributes.
        """
        
        self._z = z_object.z
        self._z_err = z_object.zerr
        self._freq = z_object.freq
        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pterr = np.zeros_like(self._z, dtype=np.float)
 
        if self._z_err is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f], self._pterr[idx_f] = z2pt(self._z[idx_f], 
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print 'Singular Matrix at {0:.5g} Hz'.format(
                                                       self._freq[idx_f])
                    except AttributeError:
                        print 'Computed singular matrix'
                        print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)
            
        #--> if there is not error to the impedance tensor
        else:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                except MTex.MTpyError_PT:
                    try:
                        print 'Singular Matrix at {0:.5g}'.format(
                                                   self._freq[idx_f])
                    except AttributeError:
                        print 'Computed singular matrix'
                        print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)

        self.rotation_angle = z_object.rotation_angle
        
    # def _get_z_object(self):
    #     z_object = MTz.Z(z_array=self._z, zerr_array=self._z_err)
    #     z_object.freq = self._freq
    #     z_object.rotation_angle = self.rotation_angle
        
    #     return z_object
        
    # _z_object = property(_get_z_object, _set_z_object, 
    #                     doc="class mtpy.core.z.Z")


    #---z array---------------------------------------------------------------
    def _set_z(self, z_array):
        """
            Set  Z array as PhaseTensor object attribute.
        """

        self._z = z_array
        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pterr = np.zeros_like(self._z, dtype=np.float)

        if self._z_err is not None and self._z is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f], self._pterr[idx_f] = z2pt(self._z[idx_f], 
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print 'Singular Matrix at {0:.5g} Hz'.format(
                                                       self._freq[idx_f])
                    except AttributeError:
                        print 'Computed singular matrix'
                        print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)
            
        #--> if there is not error to the impedance tensor
        elif self._z is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                except MTex.MTpyError_PT:
                    try:
                        print 'Singular Matrix at {0:.5g}'.format(
                                                   self._freq[idx_f])
                    except AttributeError:
                        print 'Computed singular matrix'
                        print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)
                

                
    # def _get_z(self):
    #     return self._z
        
    # z = property(_get_z, _set_z, 
    #              doc="impedance tensor numpy.array((nf, 2, 2))")
                 
    #---Z Error array---------------------------------------------------------------
    def _set_z_err(self,z_err_array):
        """
            Set  Z-error array as PhaseTensor object attribute.
        """

        self._z_err = z_err_array
        if self._z.shape!=self._z_err.shape:
            print 'z and z_err are not the not the same shape, setting '+\
                  'z_err to None'

        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pterr = np.zeros_like(self._z, dtype=np.float)
 
        if self._z_err is not None:
            for idx_f in range(len(self._z)):
                try:
                    self.pt[idx_f], self.pterr[idx_f] = z2pt(self._z[idx_f], 
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print 'Singular Matrix at {0:.5g} Hz'.format(
                                                       self._freq[idx_f])
                    except AttributeError:
                        print 'Computed singular matrix'
                        print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)
            
            #--> if there is not error to the impedance tensor
            else:
                for idx_f in range(len(self.z)):
                    try:
                        self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                    except MTex.MTpyError_PT:
                        try:
                            print 'Singular Matrix at {0:.5g}'.format(
                                                       self._freq[idx_f])
                        except AttributeError:
                            print 'Computed singular matrix'
                            print '  --> pt[{0}]=np.zeros((2,2))'.format(idx_f)

    # def _get_z_err(self):
    #     return self._z_err
        
    # z_err = property(_get_z_err, _set_z_err, 
    #                  doc="impedance tensor numpy.array((nf, 2, 2))")
        
                

    #==========================================================================
    #  define get methods for read only properties
    #==========================================================================
    #---invariants-------------------------------------------------------------
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

    #---trace-------------------------------------------------------------
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


        return [tr, tr_err]

    trace = property(_get_trace, doc= "")

    #---alpha-------------------------------------------------------------
    def _get_alpha(self):
        """
            Return the principal axis angle (strike) of PT in degrees 
			(incl. uncertainties).

            Output:
            - Alpha - Numpy array
            - Error of Alpha - Numpy array

        """

        if self.pt is None:
            return None, None

        alpha = np.degrees(0.5 * np.arctan2( self.pt[:,0,1] + self.pt[:,1,0],
                                            self.pt[:,0,0] - self.pt[:,1,1]))
        
        alphaerr = None
        if self.pterr is not None:
            alphaerr = np.zeros_like(alpha)
            y = self.pt[:,0,1] + self.pt[:,1,0]
            yerr = np.sqrt( self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2  )
            x = self.pt[:,0,0] - self.pt[:,1,1] 
            xerr = np.sqrt( self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2  )

            alphaerr[:] = 0.5 / ( x**2 + y**2) * np.sqrt(y**2 * xerr**2 + \
                                                         x**2 * yerr**2 )

        return [alpha, alphaerr]
        
    alpha = property(_get_alpha, doc = "")

    #---beta-------------------------------------------------------------
    def _get_beta(self):
        """
            Return the 3D-dimensionality angle Beta of PT in degrees 
            (incl. uncertainties).

            Output:
            - Beta - Numpy array
            - Error of Beta - Numpy array

        """

        if self.pt is None:
            return None, None
        
        beta = np.degrees(0.5 * np.arctan2( self.pt[:,0,1] - self.pt[:,1,0],
                                            self.pt[:,0,0] + self.pt[:,1,1]))

        betaerr = None

        if self.pterr is not None:
            betaerr = np.zeros_like(beta)

            y = self.pt[:,0,1] - self.pt[:,1,0]
            yerr = np.sqrt( self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2  )
            x = self.pt[:,0,0] + self.pt[:,1,1] 
            xerr = np.sqrt( self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2  )

            betaerr[:] = 0.5 / ( x**2 + y**2) * np.sqrt( y**2 * xerr**2 +\
                                                         x**2 * yerr**2 )

        return [beta, betaerr]

    beta = property(_get_beta, doc="")

    #---skew-------------------------------------------------------------
    def _get_skew(self):
        """
            Return the skew of PT (incl. uncertainties).

            Output:
            - Skew(PT) - Numpy array
            - Error of Skew(PT) - Numpy array

        """
        if self.pt is None:
            return None, None
       
        skew = np.array([i[0,1] - i[1,0] for i in self.pt])
        
        skewerr = None
        if self.pterr is not None:
            skewerr = np.zeros_like(skew)
            skewerr[:] = self.pterr[:,0,1] + self.pterr[:,1,0]

        return [skew, skewerr]

    skew = property(_get_skew, doc="Skew angle in degrees")

    #---azimuth (strike angle)-------------------------------------------------
    def _get_azimuth(self):
        """
        Returns the azimuth angle related to geoelectric strike in degrees
        including uncertainties
        
        Returns:
        --------
            **azimuth(pt)** : numpy.array(nf)
                              azimuth angles in degrees assuming North is 0
                              and angle is positive clockwise
                              
            **azimuth_err** : numpy.array(nf)
                              azimuth angle errors in degrees
                              
        """

        if self.pt is None:
            return None, None
            
        az = self.alpha[0] - self.beta[0]
        
        if self.pterr is not None:
            az_err = np.sqrt(self.alpha[1]+self.beta[1])
        else:
            az_err = None
            
        return [az, az_err]
        
    azimuth = property(_get_azimuth, 
                       doc="Azimuth angle (deg) related to geoelectric strike")
                       
    #---ellipticity----------------------------------------------------
    def _get_ellipticity(self):
        """
        Returns the ellipticity of the phase tensor, related to dimesionality
        
        Returns:
        --------
            **ellipticity** : np.array(nf)
                              ellipticity values
                              
            **ellipticity_err** : np.array(nf)
                                  ellipticity errors
                                  
        """
        
        if self.pt is None:
            return None, None
            
        ellip = (self.phimax[0]-self.phimin[0])/(self.phimax[0]+self.phimin[0])
        
        if self.pterr is not None:
            ellip_err = ellip * np.sqrt(self.phimax[1]+self.phimin[1])*\
                        np.sqrt((1/(self.phimax[0]-self.phimin[0]))**2+\
                        (1/(self.phimax[0]+self.phimin[0]))**2)
        else:
            ellip_err = None
            
        return [ellip, ellip_err]
        
    ellipticity = property(_get_ellipticity,
                           doc="Ellipticity of phase tensor related to "+\
                               "dimensionality")
    #---det-------------------------------------------------------------
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
            det_phi_err[:] = np.abs(self.pt[:,1,1] * self.pterr[:,0,0]) +\
                             np.abs(self.pt[:,0,0] * self.pterr[:,1,1]) +\
                             np.abs(self.pt[:,0,1] * self.pterr[:,1,0]) +\
                             np.abs(self.pt[:,1,0] * self.pterr[:,0,1])

        return [det_phi, det_phi_err]

    det = property(_get_det, doc = "")

    #---principle component 1----------------------------------------------
    def _pi1(self):
        """
            Return Pi1 (incl. uncertainties).
            
            Pi1 is calculated according to Bibby et al. 2005: 
			Pi1 = 0.5 * sqrt(PT[0,0]-PT[1,1])**2 + (PT[0,1]+PT[1,0])**2)

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """
        #after bibby et al. 2005

        pi1 = 0.5 * np.sqrt((self.pt[:,0,0] - self.pt[:,1,1])**2 +\
                             (self.pt[:,0,1] + self.pt[:,1,0])**2)
        pi1err = None

        if self.pterr is not None:
            pi1err = 1./ pi1 * np.sqrt( (self.pt[:,0,0] - self.pt[:,1,1] )**2 *\
                              (self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2)  +\
                              (self.pt[:,0,1] + self.pt[:,1,0] )**2 * \
                              (self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2) )
        return pi1, pi1err
        
    #---principle component 2----------------------------------------------
    def _pi2(self):
        """
            Return Pi1 (incl. uncertainties).
            
            Pi1 is calculated according to Bibby et al. 2005: 
			Pi1 = 0.5 * sqrt(PT[0,0]+PT[1,1])**2 + (PT[0,1]-PT[1,0])**2)

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """
        #after bibby et al. 2005

        pi2 = 0.5 * np.sqrt((self.pt[:,0,0] + self.pt[:,1,1])**2 +\
                              (self.pt[:,0,1] - self.pt[:,1,0])**2)
        pi2err = None

        if self.pterr is not None:
            pi2err = 1./ pi2 * np.sqrt( (self.pt[:,0,0] + self.pt[:,1,1] )**2 *
                    (self.pterr[:,0,0]**2 + self.pterr[:,1,1]**2)  +\
                    (self.pt[:,0,1] - self.pt[:,1,0] )**2 *\
                    (self.pterr[:,0,1]**2 + self.pterr[:,1,0]**2) )


        return pi2, pi2err
   

    #---phimin----------------------------------------------
    def _get_phimin(self):
        """
            Return the angle Phi_min of PT (incl. uncertainties).
            
            Phi_min is calculated according to Bibby et al. 2005: 
                Phi_min = Pi2 - Pi1 

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """

        if self.pt is None:
            return None, None

        phimin = self._pi2()[0] - self._pi1()[0]

        phiminerr = None
        if self.pterr is not None:
            phiminerr = np.sqrt(self._pi2()[1]**2+self._pi1()[1]**2)
 
            return [np.degrees(np.arctan(phimin)), np.degrees(np.arctan(phiminerr))]
        else:
            return [np.degrees(np.arctan(phimin)),None]

    phimin = property(_get_phimin, doc =" Minimum phase in degrees")

    #---phimax----------------------------------------------
    def _get_phimax(self):
        """
            Return the angle Phi_max of PT (incl. uncertainties).
            
            Phi_max is calculated according to Bibby et al. 2005: Phi_max = Pi2 + Pi1 

            Output:
            - Phi_max - Numpy array
            - Error of Phi_max - Numpy array

        """
        

        if self.pt is None:
            return None, None

        phimax = self._pi2()[0] + self._pi1()[0]

        phimaxerr = None
        if self.pterr is not None:
            phimaxerr = np.sqrt(self._pi2()[1]**2+self._pi1()[1]**2)
 
            return [np.degrees(np.arctan(phimax)), np.degrees(np.arctan(phimaxerr))]
        else:
            return [np.degrees(np.arctan(phimax)),None]

    phimax = property(_get_phimax, doc = "Maximum phase in degrees")

    
    def rotate(self,alpha):
        """
            Rotate PT array. Change the rotation angles attribute respectively.

            Rotation angle must be given in degrees. All angles are referenced to 
			geographic North, positive in clockwise direction. 
			(Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.


        """


        if self._pt is None :
            print 'pt-array is "None" - I cannot rotate that'
            return
            
        if np.iterable(self.rotation_angle) == 0:
            self.rotation_angle = np.array([self.rotation_angle 
                                             for ii in self.pt])

        #check for iterable list/set of angles - if so, it must have length 1 
        #or same as len(pt):
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
            
        self.rotation_angle = list((np.array(lo_angles) + \
                                    np.array(self.rotation_angle))%360)

        if len(lo_angles) != len(self._pt):
            print 'Wrong number Number of "angles" - need %i '%(len(self._pt))
            self.rotation_angle = 0.
            return
        
        pt_rot = copy.copy(self._pt)
        pterr_rot = copy.copy(self._pterr)
       
        for idx_freq in range(len(self._pt)):
                    
            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.pterr is not None:
                pt_rot[idx_freq], pterr_rot[idx_freq] = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle, self.pterr[idx_freq,:,:])
            else:
                pt_rot[idx_freq], pterr_rot = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle)

        
        #--> set the rotated tensors as the current attributes
        self._pt = pt_rot
        self._pterr = pterr_rot

    #---only 1d----------------------------------------------
    def _get_only1d(self):
        """
            Return PT in 1D form.

            If PT is not 1D per se, the diagonal elements are set to zero, 
            the off-diagonal elements keep their signs, but their absolute is 
            set to the mean of the original PT off-diagonal absolutes.
        """

        if self._pt is None:
            return None

        pt1d = copy.copy(self._pt)

        for i in range(len(pt1d)):
            pt1d[i,0,1] = 0
            pt1d[i,1,0] = 0
            
            mean1d = 0.5* (pt1d[i,0,0]+pt1d[i,1,1])
            pt1d[i,0,0] = mean1d
            pt1d[i,1,1] = mean1d

        return pt1d

    only1d = property(_get_only1d, doc = "")

    #---only 2d----------------------------------------------
    def _get_only2d(self):
        """
            Return PT in 2D form.

            If PT is not 2D per se, the diagonal elements are set to zero.
        """
        if self._pt is None:
            return None

        pt2d = copy.copy(self._pt)

        for i in range(len(pt2d)):
            pt2d[i,0,1] = 0
            pt2d[i,1,0] = 0
            
            pt2d[i,0,0] = self.phimax[0][i]
            pt2d[i,1,1] = self.phimin[0][i]
            
        return pt2d

    only2d = property(_get_only2d, doc = "")


class ResidualPhaseTensor():
    """
        PhaseTensor class - generates a Phase Tensor (PT) object DeltaPhi

        DeltaPhi = 1 - Phi1^-1*Phi2

    """

    def __init__(self, pt_object1 = None, pt_object2 = None):
        """
            Initialise an instance of the ResidualPhaseTensor class.

            Optional input:
            pt_object1 : instance of the PhaseTensor class
            pt_object2 : instance of the PhaseTensor class

            Initialise the attributes with None
        """

        self.residual_pt = None
        self.rpt = None
        self.rpterr = None
        self._pt1 = None  
        self._pt2 = None  
        self._pt1err = None  
        self._pt2err = None 
        self.freq = None

        if pt_object1 is not None or  pt_object2 is not None:
            if not ((isinstance(pt_object1,PhaseTensor) and\
                     isinstance(pt_object2,PhaseTensor))):
                raise MTex.MTpyError_PT('ERROR - arguments must be instances '
                                        'of the PhaseTensor class')
            
        self.compute_residual_pt(pt_object1, pt_object2)


    def compute_residual_pt(self, pt_o1, pt_o2):
        """
            Read in two instance of the MTpy PhaseTensor class.

            Update attributes:
            rpt, rpterr, _pt1, _pt2, _pt1err, _pt2err  

        """

        if not ((isinstance(pt_o1, PhaseTensor)) and \
               (isinstance(pt_o2, PhaseTensor)) ):
            raise MTex.MTpyError_PT('ERROR - both arguments must be instances'
                                    'of the PhaseTensor class')

        pt1 = pt_o1.pt
        pt2 = pt_o2.pt
        self.freq = pt_o1.freq

        #--> compute residual phase tensor
        if pt1 is not None and pt2 is not None:
            try:
                if pt1.dtype not in [float,int]:
                    raise
                if pt2.dtype not in [float,int]:
                    raise
                if not pt1.shape == pt2.shape:
                    raise
                if (not len(pt1.shape) in [2,3]) :
                    raise

                if len(pt1.shape) == 3:
                    self.rpt = np.zeros_like(pt1)

                    for idx in range(len(pt1)):
                        self.rpt[idx] = np.eye(2)-np.dot(np.matrix(pt1[idx]).I,
                                                         np.matrix(pt2[idx]))
                 
                    self._pt1 = pt1  
                    self._pt2 = pt2  

                else:
                    self.rpt = np.zeros((1,2,2)) 
                    self.rpt[0] = np.eye(2)-np.dot(np.matrix(pt1).I, 
                                                   np.matrix(pt2))
                    
                    self._pt1 =  np.zeros((1,2,2))  
                    self._pt1[0] = pt1 
                    self._pt2 =  np.zeros((1,2,2))  
                    self._pt2[0] = pt2 

            except:
                raise MTex.MTpyError_PT('ERROR - both PhaseTensor objects must'
                                  ' contain valid PT arrays of the same shape')

        else:
            print  ('Could not determine ResPT - both PhaseTensor objects must'
                   'contain PT arrays of the same shape')

        
        #--> compute residual error
        pt1err = pt_o1.pterr
        pt2err = pt_o2.pterr

        if pt1err is not None and pt2err is not None:
            self.rpterr = np.zeros(self.rpt.shape)
            try:
                if (pt1err.dtype not in [float,int]) or \
                    (pt2err.dtype not in [float,int]):
                    raise
                if not pt1err.shape == pt2err.shape:
                    raise
                if (not len(pt1err.shape) in [2,3] ):
                    raise
                if self.rpterr is not None:
                    if self.rpterr.shape != pt1err.shape:
                        raise

                if len(pt1err.shape) == 3:
                    self.rpterr = np.zeros((len(pt1),2,2))

                    for idx in range(len(pt1err)):
                        matrix1 = pt1[idx]
                        matrix1err = pt1err[idx]                        
                        matrix2, matrix2err = MTcc.invertmatrix_incl_errors(
                                        pt2[idx], inmatrix_err = pt2err[idx])

                        summand1,err1 = MTcc.multiplymatrices_incl_errors(
                                            matrix2, matrix1, 
                                            inmatrix1_err = matrix2err,
                                            inmatrix2_err =  matrix1err)
                        summand2,err2 = MTcc.multiplymatrices_incl_errors(
                                            matrix1, matrix2, 
                                            inmatrix1_err = matrix1err,
                                            inmatrix2_err =  matrix2err)

                        self.rpterr[idx] = np.sqrt(0.25*err1**2 +0.25*err2**2)

                    self._pterr1 = pt1err  
                    self._pterr2 = pt2err  

                else:
                    self.rpterr = np.zeros((1,2,2))
                    self.rpterr[0] = np.eye(2) - 0.5 * np.array( 
                                    np.dot( np.matrix(pt2).I, np.matrix(pt1) ) 
                                    + np.dot( np.matrix(pt1), np.matrix(pt2).I))
                    matrix1 = pt1
                    matrix1err = pt1err                        
                    matrix2, matrix2err = MTcc.invertmatrix_incl_errors(
                                                   pt2, inmatrix_err = pt2err)

                    summand1,err1 = MTcc.multiplymatrices_incl_errors(
                                        matrix2, matrix1, 
                                        inmatrix1_err = matrix2err,
                                        inmatrix2_err =  matrix1err)
                    summand2,err2 = MTcc.multiplymatrices_incl_errors(
                                        matrix1, matrix2, 
                                        inmatrix1_err = matrix1err,
                                        inmatrix2_err =  matrix2err)

                    self.rpterr = np.sqrt(0.25*err1**2 +0.25*err2**2)
            
                    self._pt1err =  np.zeros((1,2,2))  
                    self._pt1err[0] = pt1err
                    self._pt2err =  np.zeros((1,2,2))  
                    self._pt2err[0] = pt2err 

            except:
                raise MTex.MTpyError_PT('ERROR - both PhaseTensor objects must'
                                   'contain PT-error arrays of the same shape')
                
        else:
            print  ('Could not determine Residual PT uncertainties - both'
                    ' PhaseTensor objects must contain PT-error arrays of the'
                    'same shape')
 
        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt, 
                                       pterr_array=self.rpterr,
                                       freq=self.freq)
                                       
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

        self.compute_residual_pt(pt_o1,pt_o2)
        

    def set_rpt(self, rpt_array):
        """
            Set the attribute 'rpt' (ResidualPhaseTensor array).

            Input:
            ResPT array

            Test for shape, but no test for consistency!

        """ 
        if (self.rpt is not None) and (self.rpt.shape != rpt_array.shape):
            print 'Error - shape of "ResPT" array does not match shape of existing rpt array: %s ; %s'%(str(rpt_array.shape),str(self.rpt.shape))
            return

        self.rpt = rpt_array
        
        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt, 
                                       pterr_array=self.rpterr,
                                       freq=self.freq)
        

    def set_rpterr(self, rpterr_array):
        """
            Set the attribute 'rpterr' (ResidualPhaseTensor-error array).

            Input:
            ResPT-error array

            Test for shape, but no test for consistency!

        """ 
        if (self.rpterr is not None) and (self.rpterr.shape != rpterr_array.shape):
            print 'Error - shape of "ResPT-error" array does not match shape of existing rpterr array: %s ; %s'%(str(rpterr_array.shape),str(self.rpterr.shape))
            return

        self.rpterr = rpterr_array
        
        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt, 
                                       pterr_array=self.rpterr,
                                       freq=self.freq)


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
        pterr_array[0,0] = 1/np.abs(detreal) * np.sqrt(np.sum([np.abs(-pt_array[0,0] * realz[1,1] * zerr_array[0,0])**2,
                                                                np.abs( pt_array[0,0] * realz[0,1] * zerr_array[1,0])**2,
                                                                np.abs(((imagz[0,0] * realz[1,0] - realz[0,0] * imagz[1,0]) / np.abs(detreal) * realz[0,0] ) * zerr_array[0,1])**2, 
                                                                np.abs(((imagz[1,0] * realz[0,0] - realz[1,0] * imagz[1,1]) / np.abs(detreal) * realz[0,1] ) * zerr_array[1,1])**2,
                                                                np.abs(realz[1,1] * zerr_array[0,0])**2,
                                                                np.abs(realz[0,1] * zerr_array[1,0])**2 ]))


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
            raise MTex.MTpyError_Z('Warning - z-array no. {0} contains a singular matrix,'\
            ' thus it cannot be converted into a PT!'.format(idx_f))            

        pt_array[idx_f,0,0] =  realz[1,1] * imagz[0,0] - realz[0,1] * imagz[1,0] 
        pt_array[idx_f,0,1] =  realz[1,1] * imagz[0,1] - realz[0,1] * imagz[1,1] 
        pt_array[idx_f,1,0] =  realz[0,0] * imagz[1,0] - realz[1,0] * imagz[0,0] 
        pt_array[idx_f,1,1] =  realz[0,0] * imagz[1,1] - realz[1,0] * imagz[0,1] 

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
        - PT object
    """
    #     - PT : nx2x2 real valued Numpy array
    #     - PT-error : nx2x2 real valued Numpy array

    # """

    try:
        p = PhaseTensor(z_object = z_object)
    except:
        raise MTex.MTpyError_Z('Input argument is not a valid instance of the Z class')


    # pt_array = p.pt
    # pterr_array = p.pterr

    # return pt_array, pterr_array
    return p

def _edi_object2pt(edi_object):
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
        - PT object

    """
        # Return:
        # - PT : nx2x2 real valued Numpy array
        # - PT-error : nx2x2 real valued Numpy array

    #"""

    e = MTedi.Edi()
    e.readfile(filename)

    p = PhaseTensor(z_object = e.Z)

    # pt_array = p.pt
    
    # pterr_array = p.pterr

    # return pt_array, pterr_array
    
    return p


