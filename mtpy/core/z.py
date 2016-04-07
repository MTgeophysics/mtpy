#!/usr/bin/env python

"""
=============
z module
=============

Classes
---------
    * Z --> deals with impedance tesnsor.
    * Tipper --> deals with Tipper matrix.
    
        
LK, JP 2013
"""

#=================================================================
import numpy as np
import math, cmath
import copy

import mtpy.utils.calculator as MTcc
import mtpy.utils.exceptions as MTex

reload(MTcc)

#=================================================================


#------------------------
class Z(object):
    """
    Z class - generates an impedance tensor (Z-) object.

    Methods  include rotations/combinations of Z instances, as well as 
    calculation of invariants, inverse, amplitude/phase,...
    Calculation of invariants, inverse, amplitude/phase,...

    
    Z is a complex array of the form (n_freq, 2, 2), 
    with indices in the following order: 
    
        - Zxx: (0,0) 
        - Zxy: (0,1)
        - Zyx: (1,0) 
        - Zyy: (1,1)   

    All errors are given as standard deviations (sqrt(VAR))
    
    Arguments
    
        **z_array** : numpy.ndarray(n_freq, 2, 2)
                    array containing complex impedance values
        
        **zerr_array** : numpy.ndarray(n_freq, 2, 2)
                       array containing error values (standard deviation) 
                       of impedance tensor elements 
        **freq** : np.ndarray(n_freq)
                 array of frequency values corresponding to impedance tensor
                 elements.
        
    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    freq             array of frequencies corresponding to elements of z 
    rotation_angle   angle of which data is rotated by
    z                impedance tensor
    zerr             estimated errors of impedance tensor
    resistivity      apparent resisitivity estimated from z in Ohm-m
    resistivity_err  apparent resisitivity error
    phase            impedance phase (deg)
    phase_err        error in impedance phase
    =============== ===========================================================
        
    =================== =======================================================
    Methods             Description
    =================== =======================================================
    det                  calculates determinant of z with errors
    invariants           calculates the invariants of z
    inverse              calculates the inverse of z
    no_distortion        removes distortion given a distortion matrix
    no_ss                removes static shift by assumin Z = S * Z_0 
    no_ss_no_distortion  not implemented yet
    norm                 calculates the norm of Z
    only1d               zeros diagonal components and computes 
                         the absolute valued mean of the off-diagonal 
                         components.
    only2d               zeros diagonal components 
    res_phase            computes resistivity and phase 
    rotate               rotates z positive clockwise, angle assumes
                         North is 0.
    set_res_phase        recalculates z and zerr, needs attribute freq
    skew                 calculates the invariant skew (off diagonal trace)
    trace                calculates the trace of z
    =================== =======================================================
    
    Example
	-----------
    
        >>> import mtpy.core.z as mtz
        >>> import numpy as np
        >>> z_test = np.array([[0+0j, 1+1j], [-1-1j, 0+0j]])
        >>> z_object = mtz.Z(z_array=z_test, freq=[1])
        >>> z_object.rotate(45)
        >>> z_object.res_phase


    """

    def __init__(self, z_array=None, zerr_array=None, freq=None):
        """
        Initialise an instance of the Z class.

        Arguments
    
            **z_array** : numpy.ndarray(n_freq, 2, 2)
                        array containing complex impedance values
            
            **zerr_array** : numpy.ndarray(n_freq, 2, 2)
                           array containing error values (standard deviation) 
                           of impedance tensor elements 
            
            **freq** : np.ndarray(n_freq)
                     array of frequency values corresponding to impedance
                     tensor elements.

        Initialises the attributes with None
        """    

        self._z = z_array
        self._zerr = zerr_array

        self._freq = freq
        if z_array is not None:
            if len(z_array.shape) == 2 and z_array.shape == (2,2):
                if z_array.dtype in ['complex', 'float','int']:
                    self._z = np.zeros((1, 2, 2), 'complex')
                    self._z[0] = z_array            

        if zerr_array is not None:
            if len(zerr_array.shape) == 2 and zerr_array.shape == (2,2):
                if zerr_array.dtype in ['complex', 'float','int']:
                    self._zerr = np.zeros((1, 2, 2), 'complex')
                    self._zerr[0] = zerr_array            

        self.rotation_angle = 0.
        if self._z is not None:
            self.rotation_angle = np.zeros((len(self._z)))
        
        #make attributes for resistivity and phase
        self._resistivity= None
        self._resistivity_err = None
        
        self._phase = None
        self._phase_err = None
        if self._freq is not None:
            if self._z is not None:
                self._compute_res_phase()

    #---frequency-------------------------------------------------------------
    def _set_freq(self, lo_freq):
        """
        Set the array of freq.

        Arguments
		-------------
            
            **lo_freq** : list or array of frequnecies (Hz)

        No test for consistency!
        """

        if not np.iterable(lo_freq):
            lo_freq= np.array([lo_freq])

        if self.z is not None:
            if len(self.z.shape) == 3:
                if len(lo_freq) is not len(self.z):
                    print ('length of freq list/array not correct'
                           '({0} instead of {1})'.format(len(lo_freq), 
                                                         len(self.z)))
                    return
         
        self._freq = np.array(lo_freq)
        
        #for consistency recalculate resistivity and phase
        if self._z is not None:
            try:
                self._compute_res_phase()
            except IndexError:
                print 'Need to input frequency array'

    def _get_freq(self):
            if self._freq is None:
                return None
            else:
                return np.array(self._freq)
		
    freq = property(_get_freq, _set_freq, doc='array of frequencies in Hz')
    
    #----impedance tensor -----------------------------------------------------
    def _set_z(self, z_array):
        """
        Set the attribute 'z'.

        Arguments
		-------------
            
            **z_array** : np.ndarray(nfreq, 2, 2) 
                        complex impedance tensor array

        Test for shape, but no test for consistency!

        Nulling the rotation_angle

        """         

        try:
            if len(z_array.shape) == 3 and z_array.shape[1:3] == (2,2):
                if z_array.dtype in ['complex', 'float', 'int']:
                    self._z = z_array
        except:
            try:
                if len(z_array.shape) == 2 and z_array.shape == (2,2):
                    if z_array.dtype in ['complex', 'float','int']:
                        self._z = np.zeros((1,2,2),'complex')
                        self._z[0] = z_array            
            except:
                print ('provided Z array does not have correct dimensions'
                       '- Z unchanged')

        
        if type(self.rotation_angle) is float:
            self.rotation_angle = np.array([self.rotation_angle 
                                             for ii in self._z])
                                                 
        #for consistency recalculate resistivity and phase
        if self._z is not None:
            try:
                self._compute_res_phase()
            except IndexError:
                print 'Need to input frequency array'

        
    def _get_z(self):
        return self._z
    
    z = property(_get_z, _set_z, doc="impedance tensor")

    #----impedance error-----------------------------------------------------
    def _set_zerr(self, zerr_array):
        """
        Set the attribute zerr

        Arguments
		------------
            
            **zerr_array** : np.ndarray(nfreq, 2, 2) 
                           error of impedance tensor array as standard deviation

        Test for shape, but no test for consistency!

        """ 
        if zerr_array.shape!=self.z.shape:
            print 'zerr_array shape {0} is not same shape as z {1}'.format(
                                                             zerr_array.shape,
                                                             self.z.shape) 
        self._zerr = zerr_array
        
        #for consistency recalculate resistivity and phase
        if self._zerr is not None and self._z is not None:
            try:
                self._compute_res_phase()
            except IndexError:
                print 'Need to input frequency array'
        
    def _get_zerr(self):
        return self._zerr
    
    zerr = property(_get_zerr, _set_zerr, doc='impedance tensor error')


    #---real part of impedance tensor-----------------------------------------
    def _get_real(self):
        """
        Return the real part of Z.
        """ 

        if self.z is None:
            print 'z array is None - cannot calculate real'
            return

        return np.real(self.z)

        
    def _set_real(self, real_array):
        """
        Set the real part of 'z'.

        Arguments
		-------------
            
            **real_array** : np.ndarray(nfreq, 2, 2) 
                          real part of impedance tensor array

        Test for shape, but no test for consistency!

        """         

        if (self.z is not None) and (self.z.shape != real_array.shape):
            print 'shape of "real" array does not match shape of'+\
                  'Z array: {0} ; {1}'.format(real_array.shape, self.z.shape)
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return

        i = np.complex(0,1)

        if self.z is not None:
            z_new = real_array + i * self.imag() 
        else:
            z_new = real_array 
                 

        self.z = z_new
        
        #for consistency recalculate resistivity and phase
        self._compute_res_phase()

    #real = property(_get_real, _set_real, doc='Real part of Z')
    #---imaginary part of impedance tensor------------------------------------
    def _get_imag(self):
        """
        Return the imaginary part of Z.

        """

        if self.z is None:
            print 'z array is None - cannot calculate imag'
            return


        return np.imag(self.z)

        
    def _set_imag(self, imag_array):
        """
        Set the imaginary part of 'z'.

        Arguments
		-------------
            
            **imag_array** : np.ndarray(nfreq, 2, 2) 
                           imaginary part of impedance tensor array

        Test for shape, but no test for consistency!

        """         


        if (self.z is not None) and (self.z.shape != imag_array.shape):
            print 'Error - shape of "imag" array does not match shape of'+\
                  'Z array: {0} ; {1}'.format(imag_array.shape, self.z.shape)
            return
        
        #assert real array:
        if np.linalg.norm(np.imag(imag_array )) != 0 :
            print 'Error - array "imag" is not real valued !'
            return

        i = np.complex(0,1)
        
        if self.z is not None:
            z_new = self.real() + i * imag_array 
        else:
            z_new = i * imag_array 

        self.z = z_new
        
        #for consistency recalculate resistivity and phase
        self._compute_res_phase()

    #imag = property(_get_imag, _set_imag, doc='Imaginary part of Z ')

    #-----resistivity and phase------------------------------------------------
    def _compute_res_phase(self):
        """
        Sets attributes
			* resistivity
	 		* phase
			* resistivity_err
			* phase_err*
        
        values for resistivity are in in Ohm m and phase in degrees.

        """ 
        if self.freq is None:
            print 'Need to input frequency list'
            return
            
        if self.z is None:
            print 'Z array is None - cannot calculate Res/Phase'
            return 
        
        

        self._resistivity_err = None
        self._phase_err = None
        if self.zerr is not None:
            self._resistivity_err = np.zeros_like(self.zerr)
            self._phase_err = np.zeros_like(self.zerr)

        self._resistivity = np.zeros_like(self.z, dtype='float')
        self._phase = np.zeros_like(self.z, dtype='float')

        #calculate resistivity and phase
        for idx_f in range(len(self.z)): 
            for i in range(2):                        
                for j in range(2):
                    self._resistivity[idx_f,i,j] = np.abs(self.z[idx_f,i,j])**2/\
                                                  self.freq[idx_f]*0.2
                    self._phase[idx_f,i,j] = math.degrees(cmath.phase(
                                                    self.z[idx_f,i,j]))
                
                    if self.zerr is not None:
                        
                        r_err, phi_err = MTcc.zerror2r_phi_error(
                                                 np.real(self.z[idx_f,i,j]), 
                                                 self.zerr[idx_f,i,j], 
                                                 np.imag(self.z[idx_f,i,j]), 
                                                 self.zerr[idx_f,i,j])

                        

                        self._resistivity_err[idx_f,i,j] = \
                                               0.4*np.abs(self.z[idx_f,i,j])/\
                                               self.freq[idx_f]*r_err
                        self._phase_err[idx_f,i,j] = phi_err

    
    def _get_resistivity(self): return self._resistivity
    def _get_resistivity_err(self): return self._resistivity_err
    def _get_phase(self): return self._phase
    def _get_phase_err(self): return self._phase_err
    def _set_resistivity(self, *kwargs): print "cannot be set individually - use method 'set_res_phase' !"
    def _set_resistivity_err(self, *kwargs): print "cannot be set individually - use method 'set_res_phase' !"
    def _set_phase(self, *kwargs): print "cannot be set individually - use method 'set_res_phase' !"
    def _set_phase_err(self, *kwargs): print "cannot be set individually - use method 'set_res_phase' !"



    resistivity = property(_get_resistivity, _set_resistivity, doc='Resistivity array')
    resistivity_err = property(_get_resistivity_err,_set_resistivity_err, doc='Resistivity error array')

    phase = property(_get_phase,_set_phase, doc='Phase array')
    phase_err = property(_get_phase_err,_set_phase_err, doc='Phase error array')

    def set_res_phase(self, res_array, phase_array, reserr_array = None, 
                      phaseerr_array = None):
        """
        Set values for resistivity (res - in Ohm m) and phase 
        (phase - in degrees), including error propagation.

        Updates the attributes
			* z
			* zerr

        """ 
        if self.z is not None: 
            z_new = copy.copy(self.z) 

            if self.z.shape != res_array.shape:
                print 'Error - shape of "res" array does not match shape'+\
                       'of Z array: {0} ; {1}'.format(res_array.shape,
                                                      self.z.shape)
                return

            if self.z.shape != phase_array.shape:
                print 'Error - shape of "phase" array does not match shape'+\
                      'of Z array: {0} ; {1}'.format(phase_array.shape,
                                                     self.z.shape)
                return
        else:
            z_new = np.zeros(res_array.shape,'complex')
			
            if res_array.shape != phase_array.shape:
                print 'Error - shape of "phase" array does not match shape'+\
                      'of "res" array: {0} ; {1}'.format(phase_array.shape,
                                                         res_array.shape)
                return


        if (self.freq is None) or (len(self.freq) != len(res_array)) :
            raise MTex.MTpyError_EDI('ERROR - cannot set res without correct'+\
                                     'freq information - proper "freq" '+\
                                     'attribute must be defined')

            
        #assert real array:
        if np.linalg.norm(np.imag(res_array )) != 0 :
            raise MTex.MTpyError_inputarguments( 'Error - array "res" is not'+\
                                                 'real valued !')
            
        if np.linalg.norm(np.imag(phase_array )) != 0 :
            raise MTex.MTpyError_inputarguments( 'Error - array "phase" is'+\
                                                 'not real valued !')
            
        for idx_f in range(len(z_new)):
            for i in range(2):
                for j in range(2):
                    abs_z = np.sqrt(5 * self.freq[idx_f] *\
                                   res_array[idx_f,i,j])
                    z_new[idx_f,i,j] = cmath.rect(abs_z,
                                            np.radians(phase_array[idx_f,i,j]))

        self.z = z_new
        
        #---------------------------
        # error propagation:
        if reserr_array is None or  phaseerr_array is None:
            return

        if self.zerr is not None: 
            zerr_new = copy.copy(self.zerr) 

            try:
                if self.zerr.shape != reserr_array.shape:
                    print 'Error - shape of "reserr" array does not match'+\
                          'shape of Zerr array: {0} ; {1}'.format(
                                                         reserr_array.shape,
                                                         self.zerr.shape)
                    return

                if self.zerr.shape != phaseerr_array.shape:
                    print 'Error - shape of "phase" array does not match'+\
                          'shape of Zerr array: {0} ; {1}'.format(
                                                           phase_array.shape,
                                                           self.z.shape)
                    return
            except:
                print 'Error - "phaseerr" or "reserr" is/are not array(s)'+\
                                '- Zerr not set'
                self.zerr = None
                return 

        else:
            zerr_new = np.zeros(reserr_array.shape,'float')
            try:
                if reserr_array.shape != phaseerr_array.shape:
                    print 'Error - shape of "phase" array does not match'+\
                          'shape of Zerr array: {0} ; {1}'.format(
                                                         reserr_array.shape,
                                                         self.zerr.shape)
                    return
            except:
                print 'Error - "phaseerr" or "reserr" is/are not array(s) -'+\
                      ' Zerr not set'
                return 
               
        for idx_f in range(len(zerr_new)):
            for i in range(2):
                for j in range(2):
                    abs_z = np.sqrt(5 * self.freq[idx_f] * \
                                    res_array[idx_f,i,j])
                    rel_error_res = reserr_array[idx_f,i,j]/\
                                                 res_array[idx_f,i,j]
                    #relative error varies by a factor of 0.5, which is the
                    #exponent in the relation between them:
                    abs_z_error = 0.5 * abs_z * rel_error_res

                    zerr_new[idx_f,i,j] = max(MTcc.propagate_error_polar2rect(
                                                        abs_z, 
                                                        abs_z_error, 
                                                        phase_array[idx_f,i,j], 
                                                        phaseerr_array[idx_f,i,j]))

        self.zerr = zerr_new
        
        #for consistency recalculate resistivity and phase
        self._compute_res_phase()




    def _get_inverse(self):
        """
            Return the inverse of Z.

            (no error propagtaion included yet)

        """

        if self.z is None :
            print 'z array is "None" - I cannot invert that'
            return
        
        inverse = copy.copy(self.z)
        for idx_f in range(len(inverse)):
            try:
                inverse[idx_f,:,:] = np.array( (np.matrix(self.z[idx_f,:,:])).I )
            except:
                raise MTex.MTpyError_Z('The {0}ith impedance'.format(idx_f+1)+\
                                        'tensor cannot be inverted')

        return inverse

    inverse = property(_get_inverse, doc='Inverse of Z')

    def rotate(self, alpha):
        """
        Rotate Z array by angle alpha. 

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North, positive in clockwise direction. 
        (Mathematically negative!)

        In non-rotated state, X refs to North and Y to East direction.

        Updates the attributes
            - *z*
            - *zerr*
            - *zrot*
            - *resistivity*
            - *phase*
            - *resistivity_err*
            - *phase_err*

        """

        if self.z is None :
            print 'Z array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length
        #1 or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha%360)
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.z]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha%360)
                except:
                    print '"Angle" must be a valid number (in degrees)'
                    return
                #make an n long list of identical angles
                lo_angles = [degreeangle for i in self.z]
            else:                    
                try:
                    lo_angles = [float(i%360) for i in alpha]
                except:
                    print '"Angles" must be valid numbers (in degrees)'
                    return
            
        self.rotation_angle = np.array([(oldangle + lo_angles[i])%360 
					  for i,oldangle in enumerate(self.rotation_angle)])

        if len(lo_angles) != len(self.z):
            print 'Wrong number Number of "angles" - I need %i '%(len(self.z))
            #self.rotation_angle = 0.
            return

        z_rot = copy.copy(self.z)
        zerr_rot = copy.copy(self.zerr)

        for idx_freq in range(len(self.z)):
                    
            angle = lo_angles[idx_freq]
            if np.isnan(angle):
                angle = 0.

            if self.zerr is not None:
                z_rot[idx_freq], zerr_rot[idx_freq] = \
                        MTcc.rotatematrix_incl_errors(self.z[idx_freq,:,:], 
                                                      angle, 
                                                      self.zerr[idx_freq,:,:])
            else:
                z_rot[idx_freq], zerr_rot = \
                            MTcc.rotatematrix_incl_errors(self.z[idx_freq,:,:],
                                                          angle)


        self.z = z_rot
        if self.zerr is not None:
            self.zerr = zerr_rot
    
        #for consistency recalculate resistivity and phase
        self._compute_res_phase()
        
    def no_ss(self, reduce_res_factor_x = 1., reduce_res_factor_y = 1.):
        """
        Remove the static shift by providing the respective correction factors
        for the resistivity in the x and y components.
        (Factors can be determined by using the "Analysis" module for the 
        impedance tensor)

        Assume the original observed tensor Z is built by a static shift S 
        and an unperturbated "correct" Z0 :
             
             * Z = S * Z0
            
        therefore the correct Z will be :
            * Z0 = S^(-1) * Z
            
        Arguments
        ------------
		
            **reduce_res_factor_x** : float or iterable list or array
                                    static shift factor to be applied to x
                                    components (ie z[:, 0, 1]).  This is 
                                    assumed to be in resistivity scale
            
            **reduce_res_factor_y** : float or iterable list or array
                                    static shift factor to be applied to y
                                    components (ie z[:, 1, 0]).  This is 
                                    assumed to be in resistivity scale

        Returns
		--------------
        
            **S** : np.ndarray ((2, 2))
                    static shift matrix, 
            
            **Z0**: corrected Z   (over all freq)

        .. note:: The factors are in resistivity scale, so the
                  entries of  the matrix "S" need to be given by their
                  square-roots! 

        """
        
        #check for iterable list/set of reduce_res_factor_x - if so, it must
        #have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_x) == 0:
            try:
                x_factor = float(reduce_res_factor_x)
            except:
                print '"reduce_res_factor_x" must be a valid numbers'
                return

            lo_x_factors = [x_factor for i in self.z]
        else:
            if len(reduce_res_factor_x) == 1:
                try:
                    x_factor = float(reduce_res_factor_x)
                except:
                    print '"reduce_res_factor_x" must be a valid numbers'
                    return
                lo_x_factors = [x_factor for i in self.z]
            else:                    
                try:
                    lo_x_factors = [x_factor for i in reduce_res_factor_x]
                except:
                    print '"reduce_res_factor_x" must be valid numbers'
                    return
            
        if len(lo_x_factors) != len(self.z):
            print 'Wrong number Number of "reduce_res_factor_x"'+\
                  '- need {0}'.format(len(self.z))
            return
  
        #check for iterable list/set of reduce_res_factor_y - if so, 
        #it must have length 1 or same as len(z):
        if np.iterable(reduce_res_factor_y) == 0:
            try:
                y_factor = float(reduce_res_factor_y)
            except:
                print '"reduce_res_factor_y" must be a valid numbers'
                return

            lo_y_factors = [y_factor for i in self.z]
        else:
            if len(reduce_res_factor_y) == 1:
                try:
                    y_factor = float(reduce_res_factor_y)
                except:
                    print '"reduce_res_factor_y" must be a valid numbers'
                    return
                lo_y_factors = [y_factor for i in self.z]
            else:                    
                try:
                    lo_y_factors = [y_factor for i in reduce_res_factor_y]
                except:
                    print '"reduce_res_factor_y" must be valid numbers'
                    return
            
        if len(lo_y_factors) != len(self.z):
            print 'Wrong number Number of "reduce_res_factor_y"'+\
                  '- need {0} '.format(len(self.z))
            return
  

        z_corrected = copy.copy(self.z)
        static_shift = np.zeros((len(self.z),2,2))

        for idx_f in range(len(self.z)):
            #correct for x-direction
            z_corrected[idx_f, 0, :] = self.z[idx_f, 0, :]/\
                                         np.sqrt(lo_x_factors[idx_f])
            #correct for y-direction
            z_corrected[idx_f, 1, :] = self.z[idx_f, 1, :]/\
                                         np.sqrt(lo_y_factors[idx_f])
            #make static shift array
            static_shift[idx_f, 0, 0] = np.sqrt(lo_x_factors[idx_f])
            static_shift[idx_f, 1, 1] = np.sqrt(lo_y_factors[idx_f])

        return  static_shift, z_corrected



    def no_distortion(self, distortion_tensor, distortion_err_tensor=None):
        """
        Remove distortion D form an observed impedance tensor Z to obtain
        the uperturbed "correct" Z0:
        Z = D * Z0

        Propagation of errors/uncertainties included
		
		Arguments
		------------
			**distortion_tensor** : np.ndarray(2, 2, dtype=real)
			                      real distortion tensor as a 2x2
			
			**distortion_err_tensor** : np.ndarray(2, 2, dtype=real), 
									  default is None
			
		Returns
		-----------
			**distortion_tensor** :  np.ndarray(2, 2, dtype='real')
			                       input distortion tensor
			**z_corrected** : np.ndarray(num_freq, 2, 2, dtype='complex')
			                impedance tensor with distorion removed
							
			**z_corrected_err** : np.ndarray(num_freq, 2, 2, dtype='complex')
							    impedance tensor error after distortion is removed
								
		Example
		----------
			>>> import mtpy.core.z as mtz
			>>> distortion = np.array([[1.2, .5],[.35, 2.1]])
			>>> d, new_z, new_z_err = z_obj.no_distortion(distortion)
        """

        if distortion_err_tensor is None:
            distortion_err_tensor = np.zeros_like(distortion_tensor)
        #for all freq, calculate D.Inverse, then obtain Z0 = D.I * Z
        try:
            if not ( len(distortion_tensor.shape) in [2,3] ) and \
                   (len(distortion_err_tensor.shape) in [2,3]):
                raise
            if len(distortion_tensor.shape) == 3 or \
                len(distortion_err_tensor.shape) == 3:
                print 'Distortion is not time-dependent - take only first'+\
                      'of given distortion tensors'
                try:
                    distortion_tensor = distortion_tensor[0]
                    distortion_err_tensor = distortion_err_tensor[0]
                except:
                    raise

            if not (distortion_tensor.shape == (2,2) ) and \
                    (distortion_err_tensor.shape == (2,2) ):
                raise

            distortion_tensor = np.matrix(np.real(distortion_tensor))

        except:
            raise MTex.MTpyError_Z('The array provided is not a proper'+\
                                   'distortion tensor')

        try: 
            DI = distortion_tensor.I
        except:
            raise MTex.MTpyError_Z('The provided distortion tensor is'+\
                                   'singular - I cannot invert that!')

        #propagation of errors (using 1-norm) - step 1 - inversion of D:
        DI_err = np.zeros_like(distortion_err_tensor)

        #todo :include error on  determinant!!
        #D_det = np.linalg.det(distortion_tensor)

        dummy, DI_err = MTcc.invertmatrix_incl_errors(distortion_tensor,
                                                      distortion_err_tensor)

        #propagation of errors - step 2 - product of D.inverse and Z;
        #D.I * Z, making it 4 summands for each component:
        z_corrected = np.zeros_like(self.z)
        z_corrected_err = np.zeros_like(self.zerr)
 
        for idx_f in range(len(self.z)):
            z_corrected[idx_f] = np.array(np.dot(DI, np.matrix(self.z[idx_f])))           
            for i in range(2):
                for j in range(2):
                    z_corrected_err[idx_f,i,j] = np.sum(np.abs(
                                                        np.array([DI_err[i,0]*\
                                                        self.z[idx_f,0,j],\
                                                        DI[i,0]*\
                                                        self.zerr[idx_f,0,j],\
                                                        DI_err[i,1]*\
                                                        self.z[idx_f,1,j],\
                                                        DI[i,1]*\
                                                        self.zerr[idx_f,1,j]]))) 

        return distortion_tensor , z_corrected, z_corrected_err



    def no_ss_no_distortion(self, rho_x = 1., rho_y = 1.):
        """
            Not implemented yet!!
        """

        pass


    def _get_only1d(self):
        """
        Return Z in 1D form.

        If Z is not 1D per se, the diagonal elements are set to zero, 
        the off-diagonal elements keep their signs, but their absolute 
        is set to the mean of the original Z off-diagonal absolutes.
        """

        z1d = copy.copy(self.z)

        for i in range(len(z1d)):
            z1d[i,0,0] = 0
            z1d[i,1,1] = 0
            sign01 = np.sign(z1d[i,0,1])
            sign10 = np.sign(z1d[i,1,0])
            mean1d = 0.5* (z1d[i,1,0]+z1d[i,0,1])
            z1d[i,0,1] = sign01 * mean1d
            z1d[i,1,0] = sign10 * mean1d

        return z1d

    only1d = property(_get_only1d, 
                      doc=""" Return Z in 1D form. If Z is not 1D per se, 
                              the diagonal elements are set to zero, the 
                              off-diagonal elements keep their signs, but 
                              their absolute is set to the mean of the 
                              original Z off-diagonal absolutes.""")

    def _get_only2d(self):
        """
        Return Z in 2D form.

        If Z is not 2D per se, the diagonal elements are set to zero.
        """

        z2d = copy.copy(self.z)

        for i in range(len(z2d)):
            z2d[i,0,0] = 0
            z2d[i,1,1] = 0
            
        return z2d
    
    only2d = property(_get_only2d, 
                      doc="""Return Z in 2D form. If Z is not 2D per se,
                             the diagonal elements are set to zero. """)


    def _get_trace(self):
        """
        Return the trace of Z (incl. uncertainties).

        Returns
	    -------------
            **tr** : np.ndarray(nfreq, 2, 2) 
                    Trace(z)
            **tr_err** : np.ndarray(nfreq, 2, 2)
                       Error of Trace(z)

        """

        tr = np.array( [np.trace(i) for i in self.z])

        tr_err = None
        if self.zerr is not None:
            tr_err = np.zeros_like(tr)
            tr_err[:] = self.zerr[:,0,0] + self.zerr[:,1,1]


        return tr, tr_err

    trace = property(_get_trace, doc='Trace of Z, incl. error')

    def _get_skew(self):
        """
        Return the skew of Z (incl. uncertainties).

        Returns
	    -----------
            **skew**: np.ndarray(nfreq, 2, 2) 
                    skew(z)
            **skew_err** : np.ndarray(nfreq, 2, 2)
                         Error of skew(z)

        """
        
        skew =  np.array( [ i[0,1] - i[1,0] for i in self.z ] )
        
        skewerr = None
        if self.zerr is not None:
            skewerr = np.zeros_like(skew)
            skewerr[:] = self.zerr[:,0,1] + self.zerr[:,1,0]

        return skew, skewerr
    skew = property(_get_skew, doc='Skew of Z, incl. error')

    def _get_det(self):
        """
        Return the determinant of Z (incl. uncertainties).

        Returns
		----------
            **det_Z** : np.ndarray(nfreq) 
                      det(z)
            **det_Z_err** : np.ndarray(nfreq)
                          Error of det(z)

        """

        det_Z = np.array( [np.linalg.det(i) for i in self.z])
        
        det_Z_err = None
        if self.zerr is not None:
            det_Z_err = np.zeros_like(det_Z)
            det_Z_err[:] = np.abs(self.z[:,1,1] * self.zerr[:,0,0]) +\
                           np.abs(self.z[:,0,0] * self.zerr[:,1,1]) +\
                           np.abs(self.z[:,0,1] * self.zerr[:,1,0]) +\
                           np.abs(self.z[:,1,0] * self.zerr[:,0,1])

        return det_Z, det_Z_err
    det = property(_get_det, doc='Determinant of Z, incl. error')


    def _get_norm(self):
        """
        Return the 2-/Frobenius-norm of Z (NO uncertainties yet).
        
        Returns
		---------
            **znorm** : np.ndarray(nfreq) 
                      norm(z)
            **znormerr** : np.ndarray(nfreq)
                         Error of norm(z)

        """

        znorm = np.array( [np.linalg.norm(i) for i in self.z ])            
        znormerr = None

        if self.zerr is not None:
            znormerr = np.zeros_like(znorm)
            for idx,z_tmp in enumerate(self.z):
                value = znorm[idx]
                error_matrix = self.zerr[idx]
                radicand = 0.
                for i in range(2):
                    for j in range(2):
                        radicand += (error_matrix[i,j]*np.real(z_tmp[i,j]))**2
                        radicand += (error_matrix[i,j]*np.imag(z_tmp[i,j]))**2
                       
                znormerr[idx] = 1./value*np.sqrt(radicand)


        return znorm, znormerr

    norm = property(_get_norm, doc='Norm of Z, incl. error')

    def _get_invariants(self):
        """
        Return a dictionary of Z-invariants.

        Contains
		-----------
			* z1
			* det
			* det_real
			* det_imag
			* trace
			* skew
			* norm
			* lambda_plus/minus,
			* sigma_plus/minus
        """


        invariants_dict = {}

        z1 = (self.z[:,0,1] - self.z[:,1,0])/2.
        invariants_dict['z1'] = z1 

        invariants_dict['det'] = self.det[0]
        
        det_real = np.array([np.linalg.det(i) for i in np.real(self.z)])
        invariants_dict['det_real'] = det_real
        
        det_imag = np.array([np.linalg.det(i) for i in np.imag(self.z)])
        invariants_dict['det_imag'] = det_imag

        invariants_dict['trace'] = self.trace[0]
        
        invariants_dict['skew'] = self.skew[0]
        
        invariants_dict['norm'] = self.norm[0]
        
        lambda_plus = np.array([z1[i] + np.sqrt(z1[i] * z1[i] -\
                                self.det[0][i]) for i in range(len(z1))])
        invariants_dict['lambda_plus'] = lambda_plus
        
        lambda_minus = np.array([z1[i] - np.sqrt(z1[i] * z1[i] -\
                                 self.det[0][i]) for i in range(len(z1))])
        invariants_dict['lambda_minus'] = lambda_minus
        
        sigma_plus = np.array([0.5*self.norm[0][i]**2 + \
                              np.sqrt( 0.25*self.norm[0][i]**4 + \
                              np.abs(self.det[0][i])**2) 
                              for i in range(len(self.norm[0]))])
                              
        invariants_dict['sigma_plus'] = sigma_plus
        
        sigma_minus = np.array([0.5*self.norm[0][i]**2 - \
                               np.sqrt(0.25*self.norm[0][i]**4 + \
                               np.abs(self.det[0][i])**2) 
                               for i in range(len(self.norm[0]))])
        invariants_dict['sigma_minus'] = sigma_minus

        return invariants_dict
        
    invariants = property(_get_invariants, 
                          doc="""Dictionary, containing the invariants of
                                 Z: z1, det, det_real, det_imag, trace, 
                                 skew, norm, lambda_plus/minus, 
                                 sigma_plus/minus""")

#======================================================================
#                               TIPPER
#======================================================================

class Tipper(object):
    """
    Tipper class --> generates a Tipper-object.

    Errors are given as standard deviations (sqrt(VAR))
        
    Arguments
	-----------

        **tipper_array** : np.ndarray((nf, 1, 2), dtype='complex')
                          tipper array in the shape of [Tx, Ty]
                          *default* is None
                           
        **tippererr_array** : np.ndarray((nf, 1, 2))
                             array of estimated tipper errors
                               in the shape of [Tx, Ty].
                               Must be the same shape as tipper_array.
                               *default* is None
                               
        **freq** : np.ndarray(nf)
                   array of frequencies corresponding to the tipper elements.
                   Must be same length as tipper_array.
                   *default* is None
                   
    =============== ===========================================================
    Attributes      Description
    =============== ===========================================================
    freq            array of frequencies corresponding to elements of z 
    rotation_angle  angle of which data is rotated by
        
    tipper          tipper array
    tippererr       tipper error array
    =============== ===========================================================
        
    =============== ===========================================================
    Methods         Description
    =============== ===========================================================
    mag_direction   computes magnitude and direction of real and imaginary
                    induction arrows.
    amp_phase       computes amplitude and phase of Tx and Ty.
    rotate          rotates the data by the given angle
    =============== ===========================================================
                           
        

    """

    def __init__(self, tipper_array=None, tippererr_array=None, 
                 freq=None):
        """
        Initialise an instance of the Tipper class.
        
        Arguments
		--------------

            **tipper_array** : np.ndarray((nf, 1, 2), dtype='complex')
                               tipper array in the shape of [Tx, Ty]
                               *default* is None
                               
            **tippererr_array** : np.ndarray((nf, 1, 2))
                                   array of estimated tipper errors
                                   in the shape of [Tx, Ty].
                                   Must be the same shape as tipper_array.
                                   *default* is None
                                   
            **freq** : np.ndarray(nf)
                       array of frequencies corresponding to the tipper 
                       elements.
                       Must be same length as tipper_array.
                       *default* is None
                          
        """    

        self._tipper = tipper_array        
        self._tippererr = tippererr_array
        self._freq = freq

        self.rotation_angle = 0.
        if self.tipper is not None:
            self.rotation_angle = np.zeros((len(self.tipper)))
            
        self.amplitude = None
        self.amplitude_err = None
        self._phase = None
        self._phase_err = None
        
        self.mag_real = None
        self.mag_imag = None
        self.angle_real = None
        self.angle_imag = None


    #==========================================================================
    # Define get/set and properties
    #==========================================================================
    #----freq----------------------------------------------------------    
    def _set_freq(self, lo_freq):
        """
        Set the array of freq.

        Arguments
		-----------
            
            **lo_freq** : list or array of frequnecies (Hz)

        No test for consistency!
        """

        if len(lo_freq) is not len(self.tipper):
            print 'length of freq list/array not correct'+\
                  ' (%i instead of %i)'%(len(lo_freq), len(self.tipper))
            return

        self._freq = np.array(lo_freq)
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase() 

    def _get_freq(self): 
        if self._freq is not None:
            self._freq = np.array(self._freq)
        return self._freq
        
    freq = property(_get_freq, _set_freq, 
                           doc='array of freq')

    #---tipper--------------------------------------------------------------  
    def _set_tipper(self, tipper_array):
        """
        Set the attribute *tipper*

        Arguments
		-------------
            
            **tipper_array** : np.ndarray((nf, 1, 2), dtype='complex')
                               tipper array in the shape of [Tx, Ty]
                               *default* is None

        Test for shape, but no test for consistency!

        """         
        #make sure the array is of required shape
        try:
            if len(tipper_array.shape)==3 and tipper_array.shape[1:3]==(1,2):
                if tipper_array.dtype in ['complex', 'float','int']:
                    self._tipper = tipper_array
        except:
            pass

        #check to see if the new tipper array is the same shape as the old
        if (self._tipper!=None) and (self._tipper.shape!=tipper_array.shape):
            print 'Error - shape of "tipper" array does not match shape of '+\
                  'tipper-array: %s ; %s'%(str(tipper_array.shape),
                                           str(self.tipper.shape))
            return

        self._tipper = tipper_array
        
        #neeed to set the rotation angle such that it is an array
        if self.rotation_angle is float:
            self.rotation_angle = np.array([self.rotation_angle 
                                            for ii in self._tipper])
       
       #for consistency recalculate mag and angle
        self._compute_mag_direction()
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase() 
    
    def _get_tipper(self):
        return self._tipper
        
    tipper = property(_get_tipper, _set_tipper, doc="Tipper array")
    
    #----tipper error---------------------------------------------------------
    def _set_tippererr(self, tippererr_array):
        """
        Set the attribute *tippererr*.

        Arguments
		--------------
            **tippererr_array** : np.ndarray((nf, 1, 2))
                                   array of estimated tipper errors
                                   in the shape of [Tx, Ty].
                                   Must be the same shape as tipper_array.
                                   *default* is None

        Test for shape, but no test for consistency!

        """         

        #make sure the input array is of required shape
        try:
            if len(tippererr_array.shape)==3 and \
                                        tippererr_array.shape[1:3]==(1,2):
                if tippererr_array.dtype in ['float','int']:
                    self._tippererr = tippererr_array
        except:
            pass
        
        #make sure the error array is the same shape as tipper
        try:
            if len(self.tipper) != len(self._tippererr):
                self._tippererr = None
        except:
            pass

        
        if (self.tippererr!=None) and \
                            (self._tippererr.shape!=tippererr_array.shape):
            print 'Error - shape of "tippererr" array does not match shape '+\
                  'of tippererr array: %s ; %s'%(str(tippererr_array.shape),
                                                 str(self._tippererr.shape))
            return

        self._tippererr = tippererr_array
        
        #for consistency recalculate mag and angle
        self._compute_mag_direction()
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase()
        
    def _get_tippererr(self):
        return self._tippererr
        
    tippererr = property(_get_tippererr, _set_tippererr,
                          doc="Estimated Tipper errors")
                          
    #----real part---------------------------------------------------------
    def _get_real(self):
        """
        Return the real part of the Tipper.

        """ 
        if self.tipper is None:
            print 'tipper array is None - cannot calculate real'
            return

        return np.real(self.tipper)

    def _set_real(self, real_array):
        """
        Set the real part of 'tipper'.

        Arguments
		--------------
        
            **tipper_array** : np.ndarray((nf, 1, 2)) real part
                               tipper array in the shape of [Tx, Ty]
                               *default* is None

        Test for shape, but no test for consistency!

        """         
        

        if (self.tipper is not None) and (self.tipper.shape!=real_array.shape):
            print 'shape of "real" array does not match shape of tipper '+\
                  'array: %s ; %s'%(str(real_array.shape),
                                    str(self.tipper.shape))
            return

        #assert real array:
        if np.linalg.norm(np.imag(real_array )) != 0 :
            print 'Error - array "real" is not real valued !'
            return

        if self.tipper is not None:
            tipper_new = real_array + 1j* self.imag() 
        else:
            tipper_new = real_array

        self.tipper = tipper_new
        
        #for consistency recalculate mag and angle
        self._compute_mag_direction()
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase()

    _real = property(_get_real, _set_real, doc='Real part of the Tipper')

    #---imaginary part------------------------------------------------------
    def _get_imag(self):
        """
        Return the imaginary part of the Tipper.

        """ 

        if self.tipper is None:
            print 'tipper array is None - cannot calculate imag'
            return

        return np.imag(self.tipper)

        
    def _set_imag(self, imag_array):
        """
        Set the imaginary part of 'tipper'.

        Arguments
		--------------
        
            **tipper_array** : np.ndarray((nf, 1, 2)) imaginary part
                               tipper array in the shape of [Tx, Ty]
                               *default* is None

        Test for shape, but no test for consistency!

        """         

        if (self.tipper is not None) and (self.tipper.shape!=imag_array.shape):
            print 'shape of "real" array does not match shape of tipper '+\
                  'array: %s ; %s'%(str(imag_array.shape),
                                    str(self.tipper.shape))
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
        
        #for consistency recalculate mag and angle
        self._compute_mag_direction()
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase()

    _imag = property(_get_imag, _set_imag, doc='Imaginary part of the Tipper')

    #----amplitude and phase
    def _compute_amp_phase(self):
        """
        Sets attributes
			* *amplitude*
			* *phase*
			* *amplitude_err*
			* *phase_err*
        
        values for resistivity are in in Ohm m and phase in degrees.
        """ 
 
        if self.tipper is None:
            #print 'tipper array is None - cannot calculate rho/phi'
            return None

        self.amplitude_err = None
        self._phase_err = None
        if self.tippererr is not None:
            self.amplitude_err = np.zeros(self.tippererr.shape)
            self._phase_err = np.zeros(self.tippererr.shape)

        self.amplitude = np.zeros(self.tipper.shape)
        self._phase = np.zeros(self.tipper.shape)


        for idx_f in range(len(self.tipper)):                         
            for j in range(2):
                self.amplitude[idx_f,0,j] = np.abs(self.tipper[idx_f,0,j])
                self._phase[idx_f,0,j] = math.degrees(cmath.phase(
                                                      self.tipper[idx_f,0,j]))
                
                if self.tippererr is not None:
                    r_err, phi_err = MTcc.propagate_error_rect2polar(
                                            np.real(self.tipper[idx_f,0,j]), 
                                            self.tippererr[idx_f,0,j], 
                                            np.imag(self.tipper[idx_f,0,j]), 
                                            self.tippererr[idx_f,0,j])
                                            
                    self.amplitude_err[idx_f,0,j] = r_err
                    self._phase_err[idx_f,0,j] = phi_err

    def set_amp_phase(self, r_array, phi_array):
        """
        Set values for amplitude(r) and argument (phi - in degrees).

        Updates the attributes 
			* tipper
			* tippererr

        """ 

        if self.tipper is not None: 
                
            tipper_new = copy.copy(self.tipper) 

            if self.tipper.shape != r_array.shape:
                print 'Error - shape of "r" array does not match shape of '+\
                      'tipper array: %s ; %s'%(str(r_array.shape),
                                               str(self.tipper.shape))
                return

            if self.tipper.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape of '+\
                      'tipper array: %s ; %s'%(str(phi_array.shape),
                                               str(self.tipper.shape))
                return
        else:

            tipper_new = np.zeros(r_array.shape,'complex')

            if r_array.shape != phi_array.shape:
                print 'Error - shape of "phi" array does not match shape '+\
                       'of "r" array: %s ; %s'%(str(phi_array.shape),
                                                str(r_array.shape))
                return
       
        #assert real array:
        if np.linalg.norm(np.imag(r_array )) != 0 :
            print 'Error - array "r" is not real valued !'
            return
        if np.linalg.norm(np.imag(phi_array )) != 0 :
            print 'Error - array "phi" is not real valued !'
            return

        for idx_f in range(len(r_array)):
                for j in range(2):
                    tipper_new[idx_f,0,j] = cmath.rect(r_array[idx_f,0,j], 
                                            math.radians(phi_array[idx_f,0,j]))

        self.tipper = tipper_new
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase()
                       
    #----magnitude and direction----------------------------------------------
    def _compute_mag_direction(self):
        """
        computes the magnitude and direction of the real and imaginary 
        induction vectors.  
        
        Returns
		------------

            **mag_real** : np.array(nf)
                           magnitude of the real induction vector
            
            **ang_real** : np.array(nf)
                           angle (deg) of the real induction vector assuming 
                           that North is 0 and angle is positive clockwise
            
            **mag_imag** : np.array(nf)
                           magnitude of the imaginary induction vector
                           
            **ang_imag** : np.array(nf)
                           angle (deg) of the imaginary induction vector 
                           assuming that North is 0 and angle is positive 
                           clockwise
                           
                
        """
        

        if self.tipper is None:
            return None
        self.mag_real = np.sqrt(self.tipper[:,0,0].real**2 + \
                                self.tipper[:,0,1].real**2)
        self.mag_imag = np.sqrt(self.tipper[:,0,0].imag**2 + 
                                self.tipper[:,0,1].imag**2)
        #get the angle, need to make both parts negative to get it into the
        #parkinson convention where the arrows point towards the conductor
    
        self.angle_real = np.rad2deg(np.arctan2(-self.tipper[:,0,1].real,
                                              -self.tipper[:,0,0].real))
                                       
        self.angle_imag = np.rad2deg(np.arctan2(-self.tipper[:,0,1].imag,
                                                -self.tipper[:,0,0].imag))
        
        ## estimate error: THIS MAYBE A HACK                                        
        if self.tippererr is not None:
            self.mag_err = np.sqrt(self.tippererr[:, 0, 0]**2+ \
                                   self.tippererr[:, 0, 1]**2)
            self.angle_err = np.rad2deg(np.arctan2(self.tippererr[:, 0, 0],
                                                   self.tippererr[:, 0, 1]))%45
        else:
            self.mag_err = None
            self.angle_err = None
        
    def set_mag_direction(self, mag_real, ang_real, mag_imag, ang_imag):
        """
        computes the tipper from the magnitude and direction of the real
        and imaginary components.
        
        Updates tipper
        
        No error propagation yet
        """
        
        self.tipper[:,0,0].real = np.sqrt((mag_real**2*np.arctan(ang_real)**2)/\
                                          (1-np.arctan(ang_real)**2))
                                       
        self.tipper[:,0,1].real = np.sqrt(mag_real**2/\
                                          (1-np.arctan(ang_real)**2))
                                       
        self.tipper[:,0,0].imag = np.sqrt((mag_imag**2*np.arctan(ang_imag)**2)/\
                                       (1-np.arctan(ang_imag)**2))
                                       
        self.tipper[:,0,1].imag = np.sqrt(mag_imag**2/\
                                         (1-np.arctan(ang_imag)**2))
        #for consistency recalculate mag and angle
        self._compute_mag_direction()
                             
    #----rotate---------------------------------------------------------------
    def rotate(self, alpha):
        """
        Rotate  Tipper array.

        Rotation angle must be given in degrees. All angles are referenced
        to geographic North=0, positive in clockwise direction. 
        (Mathematically negative!)

        In non-rotated state, 'X' refs to North and 'Y' to East direction.

        Updates the attributes
			* *tipper*
			* *tippererr*
			* *rotation_angle*

        """

        if self.tipper is None :
            print 'tipper array is "None" - I cannot rotate that'
            return

        #check for iterable list/set of angles - if so, it must have length 1
        #or same as len(tipper):
        if np.iterable(alpha) == 0:
            try:
                degreeangle = float(alpha%360)
            except:
                print '"Angle" must be a valid number (in degrees)'
                return

            #make an n long list of identical angles
            lo_angles = [degreeangle for i in self.tipper]
        else:
            if len(alpha) == 1:
                try:
                    degreeangle = float(alpha%360)
                except:
                    print '"Angle" must be a valid number (in degrees)'
                    return
                #make an n long list of identical angles
                lo_angles = [degreeangle for i in self.tipper]
            else:                    
                try:
                    lo_angles = [ float(i%360) for i in alpha]
                except:
                    print '"Angles" must be valid numbers (in degrees)'
                    return
           
        self.rotation_angle = np.array([(oldangle + lo_angles[i])%360 
                              for i,oldangle in enumerate(self.rotation_angle)] )

        if len(lo_angles) != len(self.tipper):
            print 'Wrong number Number of "angles" - need %i '%(len(self.tipper))
            self.rotation_angle = 0.
            return

        tipper_rot = copy.copy(self.tipper)
        tippererr_rot = copy.copy(self.tippererr)

        for idx_freq in range(len(tipper_rot)):
            angle = lo_angles[idx_freq]

            if self.tippererr is not None:
                tipper_rot[idx_freq], tippererr_rot[idx_freq] =  \
                    MTcc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], 
                                                  angle,
                                                  self.tippererr[idx_freq,:,:] )
            else:
                tipper_rot[idx_freq], tippererr_rot = \
                    MTcc.rotatevector_incl_errors(self.tipper[idx_freq,:,:], 
                                                  angle)


 
        self.tipper = tipper_rot
        self.tippererr = tippererr_rot
        
        #for consistency recalculate mag and angle
        self._compute_mag_direction()
        
        #for consistency recalculate amplitude and phase
        self._compute_amp_phase()


#------------------------


def rotate_z(z_array, alpha, zerr_array = None):
    """
	Rotate a Z array assuming N==0 and E==90

	Arguments
	------------
		**z_array** : np.ndarray(1, 2, 2) or (2, 2) 
					  impedance tensor
		**alpha** : float
					rotation angle in degrees clockwise 
		**z_err_array** : np.ndarray(z_array.shape)
						  impedance tensor error
						  needs to be the same shape as
						  z_array							  

	Returns
	-----------
		**z_rot** : np.ndarray(z_array.shape)
					rotated impedance tensor
		**z_rot_err** : np.ndarray(z_array.shape)
						rotated impedance tensor error

	Example
	-------------
		>>> import mtpy.core.z as mtz
		>>> z_array = np.array([[0-0j, 1-1j], [0-0j, -1+1j]])
		>>> rot_z = mtz.rotate_z(z_array, 45)
    """

    z_object = _read_z_array(z_array, zerr_array)

    z_object.rotate(alpha)

    return z_object.z, z_object.zerr



def remove_distortion(z_array, distortion_tensor, 
                      distortion_err_tensor = None, zerr_array = None):
    """
	Remove the distortion from a given Z array.

	Arguments
	-------------
		**z_array** : np.ndarray(num_freq, 2, 2) 
					  impedance tensor
		**distortion_tensor** : np.ndarray(2 ,2) 
								real distorion tensor
		**distortion_err_tensor : np.ndarray(2, 2)
								error in distortion tensor
		**zerr_array** : np.ndarray(z_array.shape)
						 error in impedance tensor, 
						 same shape as z_array

	Returns
	-----------
		**z_corrected** : np.ndarray(z_array.shape)
						  impedance tensor with distorion
						  removed
		**z_corrected_err** : np.ndarray(z_array.shape)
							  error in impedance after
							  distortion is removed
		**z_array** : np.ndarray(z_array.shape)
					  input z_array

    """
    
    z_object = _read_z_array(z_array, zerr_array)
    
    distortion, z_corrected, z_corrected_err = \
               z_object.no_distortion(distortion_tensor, distortion_err_tensor)

    return  z_corrected, z_corrected_err, z_array


def remove_ss(z_array, zerr_array = None, res_x = 1., res_y = 1.):
    """
	Remove the static shift from a given Z array.

	Arguments
	-------------
		**z_array** : np.ndarray(num_freq, 2, 2) 
					  impedance tensor
					  
		**zerr_array** : np.ndarray(z_array.shape)
						 error in impedance tensor, 
						 same shape as z_array
		**res_x** : float
					amount of static shift to remove from 
					Zxx and Zxy 
		**res_y** : float
					amount of static shift to remove from 
					Zyx and Zyy 

	Returns
	-----------
		**z_corrected** : np.ndarray(z_array.shape)
						  impedance tensor with static shift
						  removed
		**static_shift** : np.ndarray(2, 2)
						   static shift removed
		**z_array** : np.ndarray(z_array.shape)
					  input z_array

    """
 

    z_object = _read_z_array(z_array, zerr_array)

    z_corrected, static_shift = z_object.no_ss()

    return z_corrected, static_shift, z_array


def remove_ss_and_distortion(z_array, zerr_array = None, res_x = 1.,
                             res_y = 1.):
    """
        Not implemented yet !!

    """
    pass



def z2resphi(z_array, periods, zerr_array = None):
    """
	Return the resistivity/phase information for Z 
	(Z must be given in units "km/s" = "mu m / nT"!!).

	Arguments
	------------
		**z_array** : np.ndarray(num_periods, 2, 2) 
					  impedance tensor
		**periods** : np.ndarray(num_periods)
					  array of periods corresponding to 
					  the elements in z_array
					  
		**zerr_array** : np.ndarray(z_array.shape)
						 error in impedance tensor, 
						 same shape as z_array

	Returns
	-----------
		**resistivity** : np.ndarray(z_array.shape)
						  apparent resistivity in Ohm-m	
		**resistivity_err** : np.ndarray(z_array.shape)
						  apparent resistivity error in Ohm-m
		**phase** : np.ndarray(z_array.shape)
					impedance phase in degrees
		**phase_err** : np.ndarray(z_array.shape)
						impedance phase error in degrees
							
	Example
	----------
		>>> res, res_err, phase, phase_err = mtz.z2resphi(z_array, 
														  periods)
    """

    z_object = _read_z_array(z_array,zerr_array)
    z_object.freq= np.array(1./periods)

    return (z_object.resistivity, z_object.resistivity_err, 
            z_object.phase, z_object.phase_err)
                

def rotate_tipper(tipper_array, alpha, tippererr_array = None):
    """
	Rotate a Tipper array

		Arguments
	------------
		**tipper_array** : np.ndarray(num_freq, 2, 2) 
					       tipper array
		**alpha** : float
					rotation angle in degrees clockwise 
		**tippererr_array** : np.ndarray(tipper_array.shape)
						  tipper error
						  needs to be the same shape as
						  tipper_array							  

	Returns
	-----------
		**tipper_rot** : np.ndarray(z_array.shape)
					rotated impedance tensor
		**tipper_rot_err** : np.ndarray(z_array.shape)
						rotated impedance tensor error

	Example
	-------------
		>>> import mtpy.core.z as mtz
		>>> rot_tip = mtz.rotate_tipper(tipper_array, 45)

    """

    tipper_object = _read_tipper_array(tipper_array, tippererr_array)

    tipper_object.rotate(alpha)

    return tipper_object.tipper, tipper_object.tippererr


def tipper2rhophi(tipper_array, tippererr_array = None):
    """
	Return values for amplitude (rho) and argument (phi - in degrees).

	Output is a 4-tuple of arrays:
	(Rho, Phi, RhoError, PhiError)
    """ 
     
    tipper_object = _read_tipper_array(tipper_array, tippererr_array )

    return tipper_object.rho_phi()


def _read_z_array(z_array, zerr_array = None):
    """
	Read a Z array and return an instance of the Z class.


	Arguments
	------------
		**z_array** : np.ndarray (num_freq, 2, 2, dtype='complex')
					  impedance tensor
		**zerr_array** : np.ndarray(z_array.shape)
		                 impedance tensor error, same shape as 
						 z_array

	Returns
	-----------
		**z_object** : mtpy.core.z.Z object
    """


    try:


        z_object = Z(z_array=z_array, zerr_array=zerr_array)
    except:
        raise MTex.MTpyError_Z('Cannot generate Z instance - check z-array'+\
                               'dimensions/type: (N,2,2)/complex ; '+\
                               '{0}'.format(z_array.shape))

    return z_object


def _read_tipper_array(tipper_array, tippererr_array = None):
    """
        Read a Tipper array and return an instance of the Tipper class.


	Arguments
	------------
		**tipper_array** : np.ndarray (num_freq, 2, 2, dtype='complex')
					       tipper tensor
		**tippererr_array** : np.ndarray(z_array.shape)
		                      tipper error, same shape as 
							  tipper_array

	Returns
	-----------
		**tipper_object** : mtpy.core.z.Tipper object

    """

    try:
        tipper_object = Tipper( tipper_array=tipper_array, 
                               tippererr_array=tippererr_array )
    except:
        raise MTex.MTpyError_tipper('Cannot generate Tipper instance -'+\
                                    'check dimensions/type: (N,1,2)/complex'+
                                    '; {0}'.format(tipper_array.shape))

    return tipper_object


def correct4sensor_orientation(Z_prime, Bx=0, By=90, Ex=0, Ey=90, 
                               Z_prime_error = None):
    """
    Correct a Z-array for wrong orientation of the sensors.

    Assume, E' is measured by sensors orientated with the angles
        E'x: a
        E'y: b
    
    Assume, B' is measured by sensors orientated with the angles
        B'x: c
        B'y: d
     
    With those data, one obtained the impedance tensor Z':
        E' = Z' * B'

    Now we define change-of-basis matrices T,U so that
        E = T * E'
        B = U * B'

    =>   T contains the expression of the E'-basis in terms of E 
    (the standard basis) 
    and  U contains the expression of the B'-basis in terms of B 
    (the standard basis)
    The respective expressions for E'x-basis vector and E'y-basis 
    vector are the columns of T.
    The respective expressions for B'x-basis vector and B'y-basis
    vector are the columns of U.

    We obtain the impedance tensor in default coordinates as: 

    E' = Z' * B' => T^(-1) * E = Z' * U^(-1) * B 
                 => E = T * Z' * U^(-1) * B
                 => Z = T * Z' * U^(-1)


    Arguments
	---------------
		**Z_prime** : np.ndarray(num_freq, 2, 2, dtype='complex')
					  impedance tensor to be adjusted

		**Bx** : float (angle in degrees)
		         orientation of Bx relative to geographic north (0)
				 *default* is 0		
		**By** : float (angle in degrees)
		         orientation of By relative to geographic north (0)
				 *default* is 90
		**Ex** : float (angle in degrees)
		         orientation of Ex relative to geographic north (0)
				 *default* is 0		
		**Ey** : float (angle in degrees)
		         orientation of Ey relative to geographic north (0)
				 *default* is 90
				 
		Z_prime_error : np.ndarray(Z_prime.shape)
		                impedance tensor error (std)
						*default* is None

	Returns
	-------------
		**Z** : np.ndarray(Z_prime.shape, dtype='complex')
		        adjusted impedance tensor
				
		**Z_err** : np.ndarray(Z_prime.shape, dtype='real')
		            impedance tensor standard deviation in 
					default orientation

	
    """
    try:
        if len(Z_prime.shape) != 2:
            raise
        if Z_prime.shape != (2,2):
            raise
        
        if Z_prime.dtype not in ['complex', 'float', 'int']:
            raise

        Z_prime = np.matrix(Z_prime)

    except:
        raise MTex.MTpyError_inputarguments('ERROR - Z array not valid!'+\
                                            'Must be 2x2 complex array')

    if Z_prime_error is not None:
        try:
            if len(Z_prime_error.shape) != 2:
                raise
            if Z_prime_error.shape != (2,2):
                raise
        
            if Z_prime_error.dtype not in [ 'float', 'int']:
                raise

        except:
            raise MTex.MTpyError_inputarguments('ERROR - Z-error array not'+\
                                               'valid! Must be 2x2 real array')


    T = np.matrix(np.zeros((2,2)))
    U = np.matrix(np.zeros((2,2)))

    dummy1 = cmath.rect(1,math.radians(Ex))


    T[0,0] = np.real(dummy1)
    T[1,0] = np.imag(dummy1)
    dummy2 = cmath.rect(1,math.radians(Ey))
    T[0,1] = np.real(dummy2)
    T[1,1] = np.imag(dummy2)

    dummy3 = cmath.rect(1,math.radians(Bx))
    U[0,0] = np.real(dummy3)
    U[1,0] = np.imag(dummy3)
    dummy4 = cmath.rect(1,math.radians(By))
    U[0,1] = np.real(dummy4)
    U[1,1] = np.imag(dummy4)

    try:
        Z = np.array(np.dot(T,np.dot(Z_prime, U.I)))
    except:
        raise MTex.MTpyError_inputarguments("ERROR - Given angles do not"+\
                           "define basis for 2 dimensions - cannot convert Z'")

    Zerr = copy.copy(Z_prime_error)

    #TODO: calculate error propagation

    return Z, Zerr
