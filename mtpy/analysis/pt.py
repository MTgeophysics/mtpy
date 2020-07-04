#!/usr/bin/env python

"""
======================
Phase Tensor
======================

Following Caldwell et al, 2004

Residual Phase Tensor following Heise et al., [2008]

@UofA, 2013
(LK)

Revised by Peacock, 2016

"""

import copy

import numpy as np

import mtpy.core.edi as MTedi
import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex


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
  
    ====================== ====================================================
    Attributes             Description    
    ====================== ====================================================
    freq                   array of frequencies associated with elements of 
                           impedance tensor.
    pt                     phase tensor array
    pt_err                 phase tensor error
    z                      impedance tensor
    z_err                  impedance error
    rotation_angle         rotation angle in degrees  
    ====================== ====================================================

    """

    def __init__(self, pt_array=None, pt_err_array=None, z_array=None,
                 z_err_array=None, z_object=None, freq=None, pt_rot=0.0):

        self._pt = pt_array
        self._pt_err = pt_err_array
        self._z = z_array
        self._z_err = z_err_array
        self._freq = freq
        self.rotation_angle = pt_rot

        # if a z object is input be sure to set the z and z_err so that the
        # pt will be calculated
        # print type(z_object)==type(MTz.Z()),isinstance(z_object, MTz.Z)
        # if isinstance(z_object, MTz.Z):
        if z_object is not None:
            try:
                self.set_z_object(z_object)
            except:
                print('\tWarning - could not digest provided Z-Object')

        elif z_array is not None:

            try:
                self._set_z(z_array)
            except:
                self._z = None
                self._z_err = None
                print('Can not calculate pt from z==None')

            if z_err_array is not None:

                try:
                    self._set_z_err(z_err_array)
                    if z_array.shape != z_err_array.shape:
                        self._set_z_err(None)
                except:
                    pass

        if self._freq is None:
            print('Should input a freq array to know which index of the' + \
                  ' PT array corresponds to which freq.')

    # ==========================================================================
    #  define get/set functions and properties
    # ==========================================================================
    # ---phase tensor array----------------------------------------
    def _set_pt(self, pt_array):
        """
            Set the attribute 'pt'.

            Input:
            Phase-Tensor array

            Test for shape, but no test for consistency!

        """
        self._pt = pt_array

        # check for dimensions
        if pt_array is not None:
            # --> large array
            if not len(pt_array.shape) in [2, 3]:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!' + \
                                        ' Invalid dimensions')

            # --> single matrix
            if not pt_array.shape[-2:] == (2, 2):
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!' + \
                                        ' Invalid dimensions')

                # --> make sure values are floats
            try:
                if not pt_array.dtype in ['float']:
                    raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!' + \
                                            'Invalid dimensions')
            except:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt array!' + \
                                        'Invalid data type (float expected)')

            if len(pt_array.shape) == 3:
                self._pt = pt_array
            else:
                self._pt = np.zeros((1, pt_array.shape[0], pt_array.shape[1]))
                self._pt[0] = pt_array

            # testing existing atributes for consistent shapes:
            try:
                if np.shape(self.pt) != np.shape(self.pt_err):
                    raise MTex.MTpyError_inputarguments('pt and pt_err are not'+\
                                                        ' the same shape')
            except:
                print('Shape of new PT array and existing pt_error do not match'+\
                      '- setting pt_error to "None"')
                self._pt_err = None
            try:
                if len(self.pt) != len(self.freq):
                    raise MTex.MTpyError_inputarguments('pt and freq are' + \
                                                        'not the same shape')
            except:
                print('Shape of new PT array and existing "freq" do not' + \
                      'match - setting freq to "None"')
                self._freq = None
            try:
                if len(self.pt) != len(self.rotation_angle):
                    raise MTex.MTpyError_inputarguments('pt and rotation angles' + \
                                                        'are not the same shape')
            except:
                print('Shape of new PT array and existing "Rotation_angle" do ' + \
                      'not match - setting rotation_angle to "None"')
                self.rotation_angle = None

        else:
            pass

    def _get_pt(self):
        return self._pt

    pt = property(_get_pt, _set_pt, doc="Phase tensor array")

    # ---phase tensor Error-----------------------------------------------------
    def _set_pt_err(self, pt_err_array):
        """
            Set the attribute 'pt_err'.

            Input:
            Phase-Tensor-error array

            Test for shape, but no test for consistency!

        """
        self._pt_err = pt_err_array

        # check dimensions
        if pt_err_array is not None:
            if not len(pt_err_array.shape) in [2, 3]:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt_err array! '+\
                      'Invalid dimensions')
            if not pt_err_array.shape[-2:] == (2, 2):
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt_err array! '+\
                      'Invalid dimensions')
            try:
                if not pt_err_array.dtype in ['float']:
                    raise 'ERROR - I cannot set new pt_err array! '
            except:
                raise MTex.MTpyError_PT('ERROR - I cannot set new pt_err array! '+\
                      'Invalid data type (float expected)')

            if self.pt is not None:
                if self.pt.shape != pt_err_array.shape:
                    raise MTex.MTpyError_PT('ERROR - I cannot set new pt_err '+\
                          'array! Invalid dimensions')


            if len(pt_err_array.shape) == 3:
                self.pt_err = pt_err_array
            else:
                self.pt_err = np.zeros((1, self.pt.shape[0], self.pt.shape[1]))
                self.pt_err[0] = pt_err_array

        else:
            pass

    def _get_pt_err(self):
        return self._pt_err

    pt_err = property(_get_pt_err, _set_pt_err,
                      doc='Phase tensor error array, must be same shape as pt')

    # ---freq------------------------------------------------------------
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
                    print('length of freq list not correct' + \
                          '(%i instead of %i)' % (len(lo_freq),
                                                  len(self._pt)))
                    return
        try:
            self._freq = np.array(lo_freq)
        except:
            self._freq = None

    def _get_freq(self):
        return self._freq

    freq = property(_get_freq, _set_freq, doc="freq array")

    # ---z_object---------------------------------------------------------------

    def set_z_object(self, z_object):
        """
            Read in Z object and convert information into PhaseTensor object 
            attributes.
        """

        self._z = z_object.z
        self._z_err = z_object.z_err
        self._freq = z_object.freq
        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pt_err = np.zeros_like(self._z, dtype=np.float)

        if self._z_err is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f], self._pt_err[idx_f] = z2pt(self._z[idx_f],
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print('Singular Matrix at {0:.5g} Hz'.format(
                            self._freq[idx_f]))
                    except AttributeError:
                        print('Computed singular matrix')
                        print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

        # --> if there is not error to the impedance tensor
        else:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                except MTex.MTpyError_PT:
                    try:
                        print('Singular Matrix at {0:.5g}'.format(
                            self._freq[idx_f]))
                    except AttributeError:
                        print('Computed singular matrix')
                        print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

        self.rotation_angle = z_object.rotation_angle

    # def _get_z_object(self):
    #     z_object = MTz.Z(z_array=self._z, z_err_array=self._z_err)
    #     z_object.freq = self._freq
    #     z_object.rotation_angle = self.rotation_angle

    #     return z_object

    # _z_object = property(_get_z_object, _set_z_object, 
    #                     doc="class mtpy.core.z.Z")


    # ---z array---------------------------------------------------------------
    def _set_z(self, z_array):
        """
            Set  Z array as PhaseTensor object attribute.
        """

        self._z = z_array
        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pt_err = np.zeros_like(self._z, dtype=np.float)

        if self._z_err is not None and self._z is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f], self._pt_err[idx_f] = z2pt(self._z[idx_f],
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print('Singular Matrix at {0:.5g} Hz'.format(
                            self._freq[idx_f]))
                    except AttributeError:
                        print('Computed singular matrix')
                        print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

        # --> if there is not error to the impedance tensor
        elif self._z is not None:
            for idx_f in range(len(self._z)):
                try:
                    self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                except MTex.MTpyError_PT:
                    try:
                        print('Singular Matrix at {0:.5g}'.format(
                            self._freq[idx_f]))
                    except AttributeError:
                        print('Computed singular matrix')
                        print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

    # def _get_z(self):
    #     return self._z

    # z = property(_get_z, _set_z, 
    #              doc="impedance tensor numpy.array((nf, 2, 2))")

    # ---Z Error array---------------------------------------------------------------
    def _set_z_err(self, z_err_array):
        """
            Set  Z-error array as PhaseTensor object attribute.
        """

        self._z_err = z_err_array
        if self._z.shape != self._z_err.shape:
            print('z and z_err are not the not the same shape, setting ' + \
                  'z_err to None')

        self._pt = np.zeros_like(self._z, dtype=np.float)
        self._pt_err = np.zeros_like(self._z, dtype=np.float)

        if self._z_err is not None:
            for idx_f in range(len(self._z)):
                try:
                    self.pt[idx_f], self.pt_err[idx_f] = z2pt(self._z[idx_f],
                                                             self._z_err[idx_f])
                except MTex.MTpyError_PT:
                    try:
                        print('Singular Matrix at {0:.5g} Hz'.format(
                            self._freq[idx_f]))
                    except AttributeError:
                        print('Computed singular matrix')
                        print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

            # --> if there is not error to the impedance tensor
            else:
                for idx_f in range(len(self.z)):
                    try:
                        self._pt[idx_f] = z2pt(self._z[idx_f])[0]
                    except MTex.MTpyError_PT:
                        try:
                            print('Singular Matrix at {0:.5g}'.format(
                                self._freq[idx_f]))
                        except AttributeError:
                            print('Computed singular matrix')
                            print('  --> pt[{0}]=np.zeros((2,2))'.format(idx_f))

    # def _get_z_err(self):
    #     return self._z_err

    # z_err = property(_get_z_err, _set_z_err, 
    #                  doc="impedance tensor numpy.array((nf, 2, 2))")



    # ==========================================================================
    #  define get methods for read only properties
    #==========================================================================
    #---invariants-------------------------------------------------------------
    @property
    def invariants(self):
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

    #---trace-------------------------------------------------------------
    @property
    def trace(self):
        """
            Return the trace of PT (incl. uncertainties).

            Output:
            - Trace(PT) - Numpy array
            - Error of Trace(PT) - Numpy array

        """
        if self.pt is None:
            return None

        return np.array([np.trace(i) for i in self.pt])

    @property
    def trace_err(self):
        tr_err = None
        if self.pt_err is not None:
            tr_err = np.zeros_like(self.trace)
            tr_err[:] = self.pt_err[:,0,0] + self.pt_err[:,1,1]
        return tr_err

    #---alpha-------------------------------------------------------------
    @property
    def alpha(self):
        """
            Return the principal axis angle (strike) of PT in degrees 
			(incl. uncertainties).

            Output:
            - Alpha - Numpy array
            - Error of Alpha - Numpy array

        """
        if self.pt is None:
            return None

        return np.degrees(0.5 * np.arctan2( self.pt[:,0,1] + self.pt[:,1,0],
                                            self.pt[:,0,0] - self.pt[:,1,1]))

    @property
    def alpha_err(self):
        alpha_err = None
        if self.pt_err is not None:
            alphaerr = np.zeros_like(self.alpha)
            y = self.pt[:,0,1] + self.pt[:,1,0]
            yerr = np.sqrt( self.pt_err[:,0,1]**2 + self.pt_err[:,1,0]**2  )
            x = self.pt[:,0,0] - self.pt[:,1,1]
            xerr = np.sqrt( self.pt_err[:,0,0]**2 + self.pt_err[:,1,1]**2  )

            alphaerr[:] = 0.5 / (x ** 2 + y ** 2) * np.sqrt(y ** 2 * xerr ** 2 + \
                                                            x ** 2 * yerr ** 2)

        return alpha_err

    #---beta-------------------------------------------------------------
    @property
    def beta(self):
        """
            Return the 3D-dimensionality angle Beta of PT in degrees 
            (incl. uncertainties).

            Output:
            - Beta - Numpy array
            - Error of Beta - Numpy array

        """

        if self.pt is None:
            return None

        return np.degrees(0.5 * np.arctan2( self.pt[:,0,1] - self.pt[:,1,0],
                                            self.pt[:,0,0] + self.pt[:,1,1]))
    @property
    def beta_err(self):
        betaerr = None

        if self.pt_err is not None:
            beta_err = np.zeros_like(self.beta)

            y = self.pt[:, 0, 1] - self.pt[:, 1, 0]
            yerr = np.sqrt(self.pt_err[:, 0, 1] ** 2 + self.pt_err[:, 1, 0] ** 2)
            x = self.pt[:, 0, 0] + self.pt[:, 1, 1]
            xerr = np.sqrt(self.pt_err[:, 0, 0] ** 2 + self.pt_err[:, 1, 1] ** 2)

            beta_err[:] = 0.5 / ( x**2 + y**2) * np.sqrt( y**2 * xerr**2 +\
                                                          x**2 * yerr**2 )

        return betaerr

    #---skew-------------------------------------------------------------
    @property
    def skew(self):
        """
            Return the skew of PT (incl. uncertainties).

            Output:
            - Skew(PT) - Numpy array
            - Error of Skew(PT) - Numpy array

        """
        if self.pt is None:
            return None
       
        return np.array([i[0,1] - i[1,0] for i in self.pt])

    @property
    def skew_err(self):
        skew_err = None
        if self.pt_err is not None:
            skew_err = np.zeros_like(self.skew)
            skew_err[:] = self.pt_err[:,0,1] + self.pt_err[:,1,0]

        return skew_err

    #---azimuth (strike angle)-------------------------------------------------
    @property
    def azimuth(self):
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
            return None
            
        return self.alpha - self.beta

    @property
    def azimuth_err(self):
        if self.pt_err is not None:
            az_err = np.sqrt(self.alpha+self.beta)
        else:
            az_err = None

        return az_err

    #---ellipticity----------------------------------------------------
    @property
    def ellipticity(self):
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
            return None

        result = None
        with np.errstate(divide='ignore', invalid='ignore'):
            result = (self.phimax-self.phimin)/(self.phimax+self.phimin)
        return result

    @property
    def ellipticity_err(self):
        if self.pt_err is not None:
            ellip_err = self.ellipticity * np.sqrt(self.phimax_err+self.phimin_err)*\
                        np.sqrt((1/(self.phimax-self.phimin))**2+\
                        (1/(self.phimax+self.phimin))**2)
        else:
            ellip_err = None

        return ellip_err

    #---det-------------------------------------------------------------
    @property
    def det(self):
        """
            Return the determinant of PT (incl. uncertainties).

            Output:
            - Det(PT) - Numpy array
            - Error of Det(PT) - Numpy array

        """
        if self.pt is None:
            return None

        return np.array([np.linalg.det(pt_arr) for pt_arr in self.pt])

    @property
    def det_err(self):
        det_phi_err = None
        if self.pt_err is not None:
            det_phi_err = np.zeros_like(self.det)
            det_phi_err[:] = np.abs(self.pt[:,1,1] * self.pt_err[:,0,0]) +\
                             np.abs(self.pt[:,0,0] * self.pt_err[:,1,1]) +\
                             np.abs(self.pt[:,0,1] * self.pt_err[:,1,0]) +\
                             np.abs(self.pt[:,1,0] * self.pt_err[:,0,1])
        return det_phi_err

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
        # after bibby et al. 2005

        pi1 = 0.5 * np.sqrt((self.pt[:, 0, 0] - self.pt[:, 1, 1]) ** 2 + \
                            (self.pt[:, 0, 1] + self.pt[:, 1, 0]) ** 2)
        pi1err = None

        if self.pt_err is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                pi1err = 1./ pi1 * np.sqrt((self.pt[:,0,0] - self.pt[:,1,1])**2*\
                                  (self.pt_err[:,0,0]**2 + self.pt_err[:,1,1]**2)+\
                                  (self.pt[:,0,1] + self.pt[:,1,0])**2 *\
                                  (self.pt_err[:,0,1]**2 + self.pt_err[:,1,0]**2))
        return pi1, pi1err

    # ---principle component 2----------------------------------------------
    def _pi2(self):
        """
            Return Pi1 (incl. uncertainties).
            
            Pi1 is calculated according to Bibby et al. 2005: 
			Pi1 = 0.5 * sqrt(PT[0,0]+PT[1,1])**2 + (PT[0,1]-PT[1,0])**2)

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """
        # after bibby et al. 2005

        pi2 = 0.5 * np.sqrt((self.pt[:, 0, 0] + self.pt[:, 1, 1]) ** 2 + \
                            (self.pt[:, 0, 1] - self.pt[:, 1, 0]) ** 2)
        pi2err = None

        if self.pt_err is not None:
            with np.errstate(divide='ignore', invalid='ignore'):
                pi2err = 1./ pi2 * np.sqrt( (self.pt[:,0,0] + self.pt[:,1,1] )**2*\
                            (self.pt_err[:,0,0]**2 + self.pt_err[:,1,1]**2) +\
                            (self.pt[:,0,1] - self.pt[:,1,0])**2*\
                            (self.pt_err[:,0,1]**2 + self.pt_err[:,1,0]**2))

        return pi2, pi2err

    #---phimin----------------------------------------------
    @property
    def phimin(self):
        """
            Return the angle Phi_min of PT (incl. uncertainties).
            
            Phi_min is calculated according to Bibby et al. 2005: 
                Phi_min = Pi2 - Pi1 

            Output:
            - Phi_min - Numpy array
            - Error of Phi_min - Numpy array

        """

        if self.pt is None:
            return None
        
#        return self._pi2()[0] - self._pi1()[0]
        return np.degrees(np.arctan(self._pi2()[0] - self._pi1()[0]))

    @property
    def phimin_err(self):
        phiminerr = None
        if self.pt_err is not None:
            phiminerr = np.sqrt(self._pi2()[1]**2+self._pi1()[1]**2)
            return np.degrees(np.arctan(phiminerr))
        else:
            return None

    #---phimax----------------------------------------------
    @property
    def phimax(self):
        """
            Return the angle Phi_max of PT (incl. uncertainties).
            
            Phi_max is calculated according to Bibby et al. 2005: Phi_max = Pi2 + Pi1 

            Output:
            - Phi_max - Numpy array
            - Error of Phi_max - Numpy array

        """

        if self.pt is None:
            return None

#        return self._pi2()[0] + self._pi1()[0]
        return np.degrees(np.arctan(self._pi2()[0] + self._pi1()[0]))

    @property
    def phimax_err(self):
        phimaxerr = None
        if self.pt_err is not None:
            phimaxerr = np.sqrt(self._pi2()[1]**2+self._pi1()[1]**2)
 
            return np.degrees(np.arctan(phimaxerr))
        else:
            return None

    def rotate(self, alpha):
        """
            Rotate PT array. Change the rotation angles attribute respectively.

            Rotation angle must be given in degrees. All angles are referenced to 
			geographic North, positive in clockwise direction. 
			(Mathematically negative!)

            In non-rotated state, X refs to North and Y to East direction.


        """

        if self._pt is None :
            print('pt-array is "None" - I cannot rotate that')
            return

        if np.iterable(self.rotation_angle) == 0:
            self.rotation_angle = np.array([self.rotation_angle
                                            for ii in self.pt])

        # check for iterable list/set of angles - if so, it must have length 1
        # or same as len(pt):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha % 360)
            except:
                print('"Angle" must be a valid number (in degrees)')
                return

            # make an n long list of identical angles
            lo_angles = [degreeangle for i in self.pt]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha % 360)
                except:
                    print('"Angle" must be a valid number (in degrees)')
                    return
                # make an n long list of identical angles
                lo_angles = [degreeangle for i in self.pt]
            else:
                try:
                    lo_angles = [float(i % 360) for i in alpha]
                except:
                    print('"Angles" must be valid numbers (in degrees)')
                    return

        self.rotation_angle = list((np.array(lo_angles) + \
                                    np.array(self.rotation_angle)) % 360)

        if len(lo_angles) != len(self._pt):
            print('Wrong number Number of "angles" - need %i ' % (len(self._pt)))
            self.rotation_angle = 0.
            return

        pt_rot = copy.copy(self._pt)
        pt_err_rot = copy.copy(self._pt_err)

        for idx_freq in range(len(self._pt)):

            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.pt_err is not None:
                pt_rot[idx_freq], pt_err_rot[idx_freq] = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle, self.pt_err[idx_freq,:,:])
            else:
                pt_rot[idx_freq], pt_err_rot = MTcc.rotatematrix_incl_errors(self.pt[idx_freq,:,:], angle)

        # --> set the rotated tensors as the current attributes
        self._pt = pt_rot
        self._pt_err = pt_err_rot

    # ---only 1d----------------------------------------------
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
            pt1d[i, 0, 1] = 0
            pt1d[i, 1, 0] = 0

            mean1d = 0.5 * (pt1d[i, 0, 0] + pt1d[i, 1, 1])
            pt1d[i, 0, 0] = mean1d
            pt1d[i, 1, 1] = mean1d

        return pt1d

    only1d = property(_get_only1d, doc="")

    # ---only 2d----------------------------------------------
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

            pt2d[i,0,0] = self.phimax[i]
            pt2d[i,1,1] = self.phimin[i]

        return pt2d

    only2d = property(_get_only2d, doc="")


class ResidualPhaseTensor():
    """
        PhaseTensor class - generates a Phase Tensor (PT) object DeltaPhi

        DeltaPhi = 1 - Phi1^-1*Phi2

    """

    def __init__(self, pt_object1=None, pt_object2=None, residualtype='heise'):
        """
            Initialise an instance of the ResidualPhaseTensor class.

            Optional input:
            pt_object1 : instance of the PhaseTensor class
            pt_object2 : instance of the PhaseTensor class

            Initialise the attributes with None
        """

        self.residual_pt = None
        self.rpt = None
        self.rpt_err = None
        self._pt1 = None
        self._pt2 = None
        self._pt1err = None
        self._pt2err = None
        self.freq = None
        self.residualtype = residualtype

        if pt_object1 is not None or  pt_object2 is not None:
            if not ((isinstance(pt_object1, PhaseTensor) and\
                     isinstance(pt_object2, PhaseTensor))):
                print(type(pt_object1), type(pt_object2))
                raise MTex.MTpyError_PT('ERROR - arguments must be instances '
                                        'of the PhaseTensor class')

        self.compute_residual_pt(pt_object1, pt_object2)

    def compute_residual_pt(self, pt_o1, pt_o2):
        """
            Read in two instance of the MTpy PhaseTensor class.

            Update attributes:
            rpt, rpt_err, _pt1, _pt2, _pt1err, _pt2err

        """

        if not ((isinstance(pt_o1, PhaseTensor)) and \
                        (isinstance(pt_o2, PhaseTensor))):
            raise MTex.MTpyError_PT('ERROR - both arguments must be instances'
                                    'of the PhaseTensor class')

        pt1 = pt_o1.pt
        pt2 = pt_o2.pt
        self.freq = pt_o1.freq

        # --> compute residual phase tensor
        if pt1 is not None and pt2 is not None:
                if pt1.dtype not in [float, int]:
                    raise ValueError
                if pt2.dtype not in [float, int]:
                    raise ValueError
                if not pt1.shape == pt2.shape:
                    raise MTex.MTpyError_PT('PT arrays not the same shape')
                if (not len(pt1.shape) in [2, 3]):
                    raise MTex.MTpyError_PT('PT array is not a valid shape')
                if self.residualtype == 'heise':
                    if len(pt1.shape) == 3:
                        self.rpt = np.zeros_like(pt1)
        
                        for idx in range(len(pt1)):
                            try:
#                                self.rpt[idx] = np.eye(2) - np.dot(np.matrix(pt1[idx]).I,
#                                                                   np.matrix(pt2[idx]))
                                self.rpt[idx] = np.eye(2) - 0.5*(np.dot(np.matrix(pt1[idx]).I,np.matrix(pt2[idx]))+\
                                                                 np.dot(np.matrix(pt2[idx]),np.matrix(pt1[idx]).I))
                            except np.linalg.LinAlgError:
                                #print 'Singular matrix at index {0}, frequency {1:.5g}'.format(idx, self.freq[idx])
                                #print 'Setting residual PT to zeros. '
                                self.rpt[idx] = np.zeros((2, 2))
        
                        self._pt1 = pt1
                        self._pt2 = pt2
        
                    else:
                        self.rpt = np.zeros((1,2,2))
                        try:
#                            self.rpt[0] = np.eye(2)-np.dot(np.matrix(pt1).I,
#                                                             np.matrix(pt2))
                            self.rpt[idx] = np.eye(2) - 0.5*(np.dot(np.matrix(pt2[idx]).I,np.matrix(pt1[idx]))+\
                                                             np.dot(np.matrix(pt1[idx]),np.matrix(pt2[idx]).I))
                            
                        except np.linalg.LinAlgError:
                            #print 'Singular matrix at frequency {0:.5g}'.format(self.freq)
                            #print 'Setting residual PT to zeros. '
                            pass
        
                        self._pt1 =  np.zeros((1,2,2))
                        self._pt1[0] = pt1
                        self._pt2 =  np.zeros((1,2,2))
                        self._pt2[0] = pt2
                elif self.residualtype == 'booker':
                    self.rpt = pt1 - pt2

        else:
            print  ('Could not determine ResPT - both PhaseTensor objects must'
                    'contain PT arrays of the same shape')


        #--> compute residual error
        pt1err = pt_o1.pt_err
        pt2err = pt_o2.pt_err

        if pt1err is not None and pt2err is not None:
            self.rpt_err = np.zeros(self.rpt.shape)
            try:
                if (pt1err.dtype not in [float,int]) or \
                    (pt2err.dtype not in [float,int]):
                    raise MTex.MTpyError_value
                if not pt1err.shape == pt2err.shape:
                    raise MTex.MTpyError_value
                if (not len(pt1err.shape) in [2,3] ):
                    raise MTex.MTpyError_value
                if self.rpt_err is not None:
                    if self.rpt_err.shape != pt1err.shape:
                        raise MTex.MTpyError_value
                if self.residualtype == 'heise':
                    if len(pt1err.shape) == 3:
                        self.rpt_err = np.zeros((len(pt1),2,2))
    
                        for idx in range(len(pt1err)):
                            matrix1 = pt1[idx]
                            matrix1err = pt1err[idx]
                            try:
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
                                self.rpt_err[idx] = np.sqrt(0.25*err1**2 +0.25*err2**2)
                            except MTex.MTpyError_inputarguments:
                                self.rpt_err[idx] = 1e10
    
    
                        self._pt_err1 = pt1err
                        self._pt_err2 = pt2err
    
                    else:
                        self.rpt_err = np.zeros((1, 2, 2))
                        try:
                            self.rpt_err[0] = np.eye(2) - 0.5 * np.array(
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
    
                            self.rpt_err = np.sqrt(0.25*err1**2 +0.25*err2**2)
                        except MTex.MTpyError_inputarguments:
                            self.rpt_err[idx] = 1e10
                
                        self._pt1err =  np.zeros((1,2,2))  
                        self._pt1err[0] = pt1err
                        self._pt2err = np.zeros((1, 2, 2))
                        self._pt2err[0] = pt2err
                elif self.residualtype == 'booker':
                    self.rpt_err = pt1err + pt2err

            except MTex.MTpyError_value:
                raise MTex.MTpyError_PT('ERROR - both PhaseTensor objects must'
                                        'contain PT-error arrays of the same shape')

        else:
            print  ('Could not determine Residual PT uncertainties - both'
                    ' PhaseTensor objects must contain PT-error arrays of the'
                    'same shape')

        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt,
                                       pt_err_array=self.rpt_err,
                                       freq=self.freq)

    def read_pts(self, pt1, pt2, pt1err=None, pt2err=None):
        """
            Read two PT arrays and calculate the ResPT array (incl. uncertainties).


            Input:
            - 2x PT array

            Optional:
            - 2x pt_error array
            
        """

        try:
            if pt1.shape != pt2.shape:
                raise

        except:
            raise MTex.MTpyError_PT('ERROR - could not build ResPT array from given PT arrays - check shapes! ')
        # TODO - check arrays here:


        pt_o1 = PhaseTensor(pt_array = pt1, pt_err_array = pt1err)
        pt_o2 = PhaseTensor(pt_array = pt2, pt_err_array = pt2err)

        self.compute_residual_pt(pt_o1, pt_o2)

    def set_rpt(self, rpt_array):
        """
            Set the attribute 'rpt' (ResidualPhaseTensor array).

            Input:
            ResPT array

            Test for shape, but no test for consistency!

        """
        if (self.rpt is not None) and (self.rpt.shape != rpt_array.shape):
            print('Error - shape of "ResPT" array does not match shape of existing rpt array: %s ; %s' % (
                str(rpt_array.shape), str(self.rpt.shape)))
            return

        self.rpt = rpt_array

        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt,
                                       pt_err_array=self.rpt_err,
                                       freq=self.freq)

    def set_rpt_err(self, rpt_err_array):
        """
            Set the attribute 'rpt_err' (ResidualPhaseTensor-error array).

            Input:
            ResPT-error array

            Test for shape, but no test for consistency!

        """
        if (self.rpt_err is not None) and (self.rpt_err.shape != rpt_err_array.shape):
            print('Error - shape of "ResPT-error" array does not match shape of existing rpt_err array: %s ; %s'%(str(rpt_err_array.shape),str(self.rpt_err.shape)))
            return

        self.rpt_err = rpt_err_array

        #--> make a pt object that is the residual phase tensor
        self.residual_pt = PhaseTensor(pt_array=self.rpt,
                                       pt_err_array=self.rpt_err,
                                       freq=self.freq)


# =======================================================================

def z2pt(z_array, z_err_array=None):
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
            if not len(z_array.shape) in [2, 3]:
                raise
            if not z_array.shape[-2:] == (2, 2):
                raise
            if not z_array.dtype in ['complex', 'float']:
                raise
        except:
            raise MTex.MTpyError_PT('Error - incorrect z array: %s;%s instead of (N,2,2);complex' % (
                str(z_array.shape), str(z_array.dtype)))

    if z_err_array is not None:
        try:
            if not len(z_err_array.shape) in [2, 3]:
                raise
            if not z_err_array.shape[-2:] == (2, 2):
                raise
            if not z_err_array.dtype in ['float']:
                raise
        except:
            raise MTex.MTpyError_PT('Error - incorrect z-err-array: %s;%s instead of (N,2,2);real' % (
                str(z_err_array.shape), str(z_err_array.dtype)))

        if not z_array.shape == z_err_array.shape:
            raise MTex.MTpyError_PT('Error - z-array and z-err-array have different shape: %s;%s' % (
                str(z_array.shape), str(z_err_array.shape)))

    # for a single matrix as input:
    if len(z_array.shape) == 2:

        pt_array = np.zeros((2, 2))

        realz = np.real(z_array)
        imagz = np.imag(z_array)
        detreal = np.linalg.det(realz)

        if detreal == 0:
            if np.linalg.norm(realz) == 0 and np.linalg.norm(imagz) == 0:
                pt_err_array = np.zeros_like(pt_array)
                if z_err_array is None:
                    pt_err_array = None
                return pt_array, pt_err_array

            else:
                raise MTex.MTpyError_PT(
                    'Error - z-array contains a singular matrix, thus it cannot be converted into a PT!')

        pt_array[0, 0] = realz[1, 1] * imagz[0, 0] - realz[0, 1] * imagz[1, 0]
        pt_array[0, 1] = realz[1, 1] * imagz[0, 1] - realz[0, 1] * imagz[1, 1]
        pt_array[1, 0] = realz[0, 0] * imagz[1, 0] - realz[1, 0] * imagz[0, 0]
        pt_array[1, 1] = realz[0, 0] * imagz[1, 1] - realz[1, 0] * imagz[0, 1]

        pt_array /= detreal

        if z_err_array is None:
            return pt_array, None

        pt_err_array = np.zeros_like(pt_array)

        #Z entries are independent -> use Gaussian error propagation (squared sums/2-norm)
        pt_err_array[0,0] = 1/np.abs(detreal) * np.sqrt(np.sum([np.abs(-pt_array[0,0] * realz[1,1] * z_err_array[0,0])**2,
                                                                np.abs( pt_array[0,0] * realz[0,1] * z_err_array[1,0])**2,
                                                                np.abs(((imagz[0,0] * realz[1,0] - realz[0,0] * imagz[1,0]) / np.abs(detreal) * realz[0,0] ) * z_err_array[0,1])**2,
                                                                np.abs(((imagz[1,0] * realz[0,0] - realz[1,0] * imagz[1,1]) / np.abs(detreal) * realz[0,1] ) * z_err_array[1,1])**2,
                                                                np.abs(realz[1,1] * z_err_array[0,0])**2,
                                                                np.abs(realz[0,1] * z_err_array[1,0])**2 ]))


        pt_err_array[0,1] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs( -pt_array[0,1] * realz[1,1] * z_err_array[0,0])**2,
                                                                np.abs(  pt_array[0,1] * realz[0,1] * z_err_array[1,0])**2,
                                                                np.abs(  ( (imagz[0,1] * realz[1,0] - realz[0,0] * imagz[1,1]) / np.abs(detreal) * realz[1,1] ) * z_err_array[0,1])**2,
                                                                np.abs(  ( (imagz[1,1] * realz[0,0] - realz[0,1] * imagz[1,0]) / np.abs(detreal) * realz[0,1] ) * z_err_array[1,1])**2,
                                                                np.abs(  realz[1,1] * z_err_array[0,1])**2,
                                                                np.abs( realz[0,1] * z_err_array[1,1])**2 ]))

        pt_err_array[1,0] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs(  pt_array[1,0] * realz[1,0] * z_err_array[0,1])**2,
                                                                np.abs( -pt_array[1,0] * realz[0,0] * z_err_array[1,1])**2,
                                                                np.abs(  ( (imagz[0,0] * realz[1,1] - realz[0,1] * imagz[1,1]) / np.abs(detreal) * realz[1,0] ) * z_err_array[0,0])**2,
                                                                np.abs(  ( (imagz[1,0] * realz[0,1] - realz[1,1] * imagz[0,0]) / np.abs(detreal) * realz[0,0] ) * z_err_array[0,1])**2,
                                                                np.abs(  realz[1,0] * z_err_array[0,0])**2,
                                                                np.abs( realz[0,0] * z_err_array[1,0])**2 ]))


        pt_err_array[1,1] = 1/np.abs(detreal) * np.sqrt( np.sum([np.abs(  pt_array[1,1] * realz[1,0] * z_err_array[0,1])**2,
                                                                np.abs( -pt_array[1,1] * realz[0,0] * z_err_array[1,1])**2,
                                                                np.abs(  ( (imagz[0,1] * realz[1,1] - realz[0,1] * imagz[1,1]) / np.abs(detreal) * realz[1,0] ) * z_err_array[0,0])**2,
                                                                np.abs(  ( (imagz[1,1] * realz[0,1] - realz[1,1] * imagz[0,1]) / np.abs(detreal) * realz[0,0] ) * z_err_array[0,1])**2,
                                                                np.abs( - realz[1,0] * z_err_array[0,1])**2,
                                                                np.abs( realz[0,0] * z_err_array[1,1])**2 ]))


        return pt_array, pt_err_array

    # else:
    pt_array = np.zeros((z_array.shape[0], 2, 2))

    for idx_f in range(len(z_array)):

        realz = np.real(z_array[idx_f])
        imagz = np.imag(z_array[idx_f])

        detreal = np.linalg.det(realz)
        if detreal == 0:
            raise MTex.MTpyError_Z('Warning - z-array no. {0} contains a singular matrix,' \
                                   ' thus it cannot be converted into a PT!'.format(idx_f))

        pt_array[idx_f, 0, 0] = realz[1, 1] * imagz[0, 0] - realz[0, 1] * imagz[1, 0]
        pt_array[idx_f, 0, 1] = realz[1, 1] * imagz[0, 1] - realz[0, 1] * imagz[1, 1]
        pt_array[idx_f, 1, 0] = realz[0, 0] * imagz[1, 0] - realz[1, 0] * imagz[0, 0]
        pt_array[idx_f, 1, 1] = realz[0, 0] * imagz[1, 1] - realz[1, 0] * imagz[0, 1]

        pt_array /= detreal

        if z_err_array is None:
            return pt_array, pt_err_array

        pt_err_array = np.zeros_like(pt_array)
        pt_err_array[idx_f,0,0] = 1/detreal * (np.abs( -pt_array[idx_f,0,0] * realz[1,1] * z_err_array[0,0]) + \
                                        np.abs(  pt_array[idx_f,0,0] * realz[0,1] * z_err_array[1,0]) + \
                                        np.abs(  (imagz[0,0] - pt_array[idx_f,0,0] * realz[0,0] ) * z_err_array[1,1]) +\
                                        np.abs(  (-imagz[1,0]+ pt_array[idx_f,0,0] * realz[1,0] ) * z_err_array[0,1]) + \
                                        np.abs(  realz[1,1] * z_err_array[0,0]) + np.abs( realz[0,1] * z_err_array[1,0]) )

        pt_err_array[idx_f,0,1] = 1/detreal * (np.abs( -pt_array[idx_f,0,1] * realz[1,1] * z_err_array[0,0]) + \
                                        np.abs(  pt_array[idx_f,0,1] * realz[0,1] * z_err_array[1,0]) + \
                                        np.abs(  (imagz[0,1] - pt_array[idx_f,0,1] * realz[0,0] ) * z_err_array[1,1]) +\
                                        np.abs(  (-imagz[1,1]+ pt_array[idx_f,0,1] * realz[1,0] ) * z_err_array[0,1]) + \
                                        np.abs(  realz[1,1] * z_err_array[0,1]) + np.abs( realz[0,1] * z_err_array[1,1]) )

        pt_err_array[idx_f,1,0] = 1/detreal * (np.abs(  (imagz[1,0] - pt_array[idx_f,1,0] * realz[1,1] ) * z_err_array[0,0]) +\
                                        np.abs( pt_array[idx_f,1,0] * realz[1,0] * z_err_array[0,1]) + \
                                        np.abs(  (-imagz[0,0] + pt_array[idx_f,1,0] * realz[0,1] ) * z_err_array[1,0]) + \
                                        np.abs( -pt_array[idx_f,1,0] * realz[0,0] * z_err_array[1,1]) + \
                                        np.abs(  realz[0,0] * z_err_array[1,0]) + np.abs( -realz[1,0] * z_err_array[0,0]) )

        pt_err_array[idx_f,1,1] = 1/detreal * (np.abs(  (imagz[1,1] - pt_array[idx_f,1,1] * realz[1,1] ) * z_err_array[0,0]) +\
                                        np.abs( pt_array[idx_f,1,1] * realz[1,0] * z_err_array[0,1]) + \
                                        np.abs(  (-imagz[0,1] + pt_array[idx_f,1,1] * realz[0,1] ) * z_err_array[1,0]) + \
                                        np.abs( -pt_array[idx_f,1,1] * realz[0,0] * z_err_array[1,1]) + \
                                        np.abs(  realz[0,0] * z_err_array[1,1]) + np.abs( -realz[1,0] * z_err_array[0,1]) )

    return pt_array, pt_err_array


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
        p = PhaseTensor(z_object=z_object)
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

    if not isinstance(edi_object, MTedi.Edi):
        raise MTex.MTpyError_EDI('Input argument is not an instance of the Edi class')
    p = PhaseTensor(edi_object=edi_object)

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

    # """

    e = MTedi.Edi()
    e.readfile(filename)

    p = PhaseTensor(z_object=e.Z)

    # pt_array = p.pt

    # pterr_array = p.pterr

    # return pt_array, pterr_array

    return p
